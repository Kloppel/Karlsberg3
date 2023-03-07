from kbp2 import file_parser

__author__ = 'Tim Meyer'

import os
import cPickle as pickle

# from file_parser import Simple_struct_parser
from charmm import Charmm_manager
import shutil, re
import kbp_tools

class MD_kbp2_project(object):
    """ A class that contains and manages several md_kbp2_job objects and all general parameters.
    """

    def __init__(self, path, structure_filename, template=None):
        """Constructor of the md_kbp2_jobs object. If path points to an existing directory the function expects that
        it contains a already existing job and tries to read previous results.

        @param path: str: A path to a folder where the job data will be stored. If the folder exists previous results
        are read. If the folder does not exist it will created. Can be None, if the object is only used as a template.
        @param structure_filename: str: Path to a PDB or CRD file containing the structure to work with. Can be None if
        the the object is only used as a template.
        @param template: md_kbp2_jobs: most parameters are copied from this object. Not allowed if 'path' points to an
        existing folder. The following parameters are not copied:
        self.structure
        self.project_folder
        self.jobs
        @return: None
        """
        if (path is not None) and (path[-1] != '/'):
            path += '/'

        # Path to the main job folder. Sub folders will be created for each job.
        self.project_folder = path

        # Path to a PDB or CRD file containing the structure to work with.
        self.structure_filename = structure_filename

        # A 'Simple_struct_parser' object containing the structure.
        # self.structure = None


        # A list that contains 'md_kbp2_task' objects.
        self.tasks = []

        # CHARMM topology and parameter files.
        self.top = None
        self.par = None

        # List of residues that are of special interest. If set, two additional MDs are used: pH5 and pH7 with the
        # protonation beeing completely dependend of the initial pKa calulation.
        self.roi = None

        # Path to a folder containing the template external_scripts for preparation and starting of the molecular dynamic
        # simulation.
        self.md_template_folder = None

        # A representation of the structure as Charmm_manger object.
        # self.__set_structure will initialize the structure (if provided) below and the first call of self.add_task
        # will execute self.charmm_struct.check_structures(). So all modelling instructions have to be set until then.
        self.charmm_struct = None

        if path is not None:
            if os.path.exists(path):
                # if template is not None:
                #     error = "The 'template' keyword is not allowed if the job folder does exist!"
                #     raise AssertionError(error)
                pass
                # Do something?
            else:
                os.mkdir(path)

            if template is not None:
                assert(isinstance(template, MD_kbp2_project))

                self.top = list(template.top)
                self.par = list(template.par)
                self.md_template_folder = template.md_template_folder
                self.roi = template.roi

        if structure_filename is not None:
            self.__set_structure(structure_filename)

        # Results from a Karlsberg+ calculation for the initial structure, stored as a kbp_results.KbpResult() object.
        self.kbp_result = None

        ### file names ###
        # Contains (names, titr_residue_dicts)
        # with: names: List of str: <task name>
        #       titr_residue_dicts: dict created by charmm.Charmm_manager.get_titr_residue_dict()
        # self.pickle_prot_tasks = 'prot_tasks.pickle'
        # !!! not used anymore !!!

    def __set_structure(self, structure):
        """Set the structure to work with. The property self.charmm_struct is initialized as self.charmm_struct object.
        The first call of self.add_task will execute self.charmm_struct.check_structures(). So all modelling
        instructions have to be set until then.

        @param structure: Can be a str containing the path to a PDB or CRD file or a Simple_struct_parser object.
        @return: None
        """
        structure_name = self.structure_filename.strip('/').split('/')[-1]
        structure_name = structure_name.rsplit('.', 1)[0]

        self.charmm_out_prefix = structure_name + '_m'

        self.charmm_struct = Charmm_manager(top=self.top, par=self.par)
        self.charmm_struct.add_structure(structure)
        self.charmm_struct.charmm_out_prefix = self.charmm_out_prefix

    def write_pickle(self, filename):
        """Stores the class as a pickle.
        @param filename: Name of the file to write the pickle into.
        @return: None
        """
        f = open(filename, 'w')
        pickle.dump(self.__dict__, f)
        f.close()

    def restore_from_pickle(self, filename):
        f = open(filename, 'r')
        self.__dict__ = pickle.load(f)
        f.close()

    def add_task(self, name, titr_residue_dict):
        # Create task folder
        project = self
        new_task = MD_kbp2_task(project, name, titr_residue_dict)

        if not new_task.modelling_competed:
            if not self.charmm_struct.structures_checked:
                self.charmm_struct.check_structures(quiet=True)

            new_task.modell_structure()

            # Store the protonation pattern in the task folder
            filename = self.project_folder + name + '/prot_task.pickle'
            f = open(filename, 'w')
            pickle.dump((name, titr_residue_dict), f)
            f.close()

            # pickle_charmm_struct_filename = self.project_folder + name + '/modelling/charmm_struct.pkl'
            #
            # f = open(pickle_charmm_struct_filename, 'wb')
            # pickle.dump(self.charmm_struct, f, protocol=2)
            # f.close()

        # if new_task.modelling_prefix is None:
        new_task.modelling_prefix = self.charmm_struct.charmm_out_prefix

        self.tasks.append(new_task)

        if os.path.exists(new_task.task_folder + 'md'):
            return

        ################################################################################################################
        ### The following section checks if there is an job started by the old script that can be partially recycled. ##
        ################################################################################################################
        # # print("Looking for previous jobs that can be used for task: %s" % name)
        # # (<path to folder>, <filter suffix>) -> only folder with th given suffix are considered.
        # cache_folders = []
        # # cache_folders.append(('/public/local/scratch/md_cache/', ''))
        # cache_folders.append(('/scratch/scratch/tmeyer/md_pka/runs/general/', 'large_longer'))
        # cache_folders.append(('/scratch/scratch/tmeyer/md_pka/runs/snase_nocal/', ''))
        #
        # usable_kbp_folders = ['kbp_e4_noMin_new2bb_wi']
        #
        # folder_found = False
        # for cache_folder, filter_suffix in cache_folders:
        #     folder_list = os.listdir(cache_folder)
        #     for folder in folder_list:
        #         ### Check whether this is a useful job folder ###
        #         folder += '/'
        #         if filter_suffix and (folder.find(filter_suffix) == -1):
        #             continue
        #         md_folder = cache_folder + folder + 'md/'
        #         if not os.path.exists(md_folder):
        #             continue
        #         kbp_folders = []
        #         for kbp_folder in usable_kbp_folders:
        #             full_kbp_folder = cache_folder + folder + kbp_folder + '/'
        #             if os.path.exists(full_kbp_folder):
        #                 kbp_folders.append(kbp_folder)
        #         # if not kbp_folders:
        #         #     continue
        #
        #         ### Check whether this is job matches the current task ###
        #         job_pickle = cache_folder + folder + 'job_dump.pickle'
        #         f = open(job_pickle, 'r')
        #         job_status = pickle.load(f)
        #         f.close()
        #         job_status = job_status['pdb']
        #         if self.structure_filename != job_status:
        #             continue
        #
        #         #### Check if the protonation states matches ###
        #         # Get the protonation state of the old job
        #         residue_list = self.charmm_struct.get_titr_residue_list(titr_residue_dict)
        #         residue_list_kbp = []
        #         for residue_descr in residue_list:
        #             residue_descr_kbp = kbp_tools.get_kbp_residue(residue_descr)
        #             residue_list_kbp.append(residue_descr_kbp)
        #
        #         reread = False
        #         prot_state_pickel_filename = cache_folder + folder + 'prot_state.pickle'
        #         if os.path.exists(prot_state_pickel_filename):
        #             # Read the protonation state from a previous run of this script.
        #             f = open(prot_state_pickel_filename, 'r')
        #             titr_residue_dict_md = pickle.load(f)
        #             f.close()
        #
        #             # Make sure that the residues in titr_residue_dict_md object loaded from the folder, matches the
        #             # residues in residue_list. This is important if residue_list has changed since the last time the
        #             # folder was checked. (e.g. because new titratable residues have been introduced)
        #             residue_descr_list = sorted(residue_list)
        #             residue_descr_list_md = []
        #             for residue_tuple in titr_residue_dict_md.keys():
        #                 resname, resid, segname = residue_tuple
        #                 residue_descr = '%s-%i_%s' % (resname, resid, segname)
        #                 residue_descr_list_md.append(residue_descr)
        #             residue_descr_list_md = sorted(residue_descr_list_md)
        #             if not residue_descr_list == residue_descr_list_md:
        #                 reread = True
        #         else:
        #             reread = True
        #         if reread:
        #             # Guess the state and write it
        #             titratable_yaml = \
        #                 kbp_tools.parse_titratable_yaml('/scratch/scratch/tmeyer/kbplus2/titratable.yaml')
        #
        #             prot_state_yaml = kbp_tools.determine_md_protonation_pattern(md_folder, residue_list_kbp,
        #                                                                          titratable_yaml)
        #
        #             titr_residue_dict_md = {}
        #             for residue_tuple, state in titr_residue_dict.iteritems():
        #                 titr_residue_dict_md[tuple(residue_tuple)] = state
        #             self.charmm_struct.mod_titr_residue_dict_from_titrable_yaml(titr_residue_dict_md, prot_state_yaml,
        #                                                                         titratable_yaml)
        #
        #             # Write the prot_state
        #             f = open(prot_state_pickel_filename, 'w')
        #             pickle.dump(titr_residue_dict_md, f)
        #             f.close()
        #
        #         # diff_str = '# Task: %s\n# Folder: %s\n  <residue> -> state md / state here\n' % (name, folder)
        #         prot_state_differs = 0
        #         for key in titr_residue_dict.keys():
        #             state = titr_residue_dict[key]
        #             state_md = titr_residue_dict_md[key]
        #             if state is None:
        #                 state = 0
        #             if state != state_md:
        #                 prot_state_differs += 1
        #                 # diff_str += "  %s -> %i / %i\n" % (key, state_md, state)
        #         # if prot_state_differs < 10:
        #         #     print diff_str
        #
        #         if prot_state_differs == 0:
        #             print("Found an useful old job for task %s -> %s" % (name,  folder))
        #
        #             # Copy useful files
        #             shutil.copytree(md_folder, new_task.task_folder + 'md')
        #             for kbp_folder in kbp_folders:
        #                 found_kbp_folder = cache_folder + folder + kbp_folder + '/'
        #                 shutil.copytree(found_kbp_folder, new_task.task_folder + kbp_folder)
        #
        #             # The folder name is used as prefix for the files produced by the MD simulations.
        #             md_prefix = folder.strip('/')
        #             new_task.md_prefix = md_prefix
        #
        #             folder_found = True
        #             break
        #     if folder_found:
        #         break
        #
        # if not folder_found:
        #     print("No useful old job for task %s found." % name)



    def __get_std_prot_tasks(self):
        """Generates and returns a standard protocol in form of a list of titr_residue_dict objects created by
        charmm_struct.get_titr_residue_dict().
        If self.kb_results is set, it is used to adjust the protocol.

        @return: list of titr_residue_dict objects created by charmm_struct.get_titr_residue_dict()
        """
        titr_residue_dicts = []
        names = []
        titr_residue_dict_template = self.charmm_struct.get_titr_residue_dict()
        mod_titr_residue_dict = self.charmm_struct.mod_titr_residue_dict
        # Creat shortcut for the function
        copy_titr_residue_dict = self.charmm_struct.copy_titr_residue_dict

        # Charged HIS
        name = 'ph5'
        titr_residue_dict = copy_titr_residue_dict(titr_residue_dict_template)
        mod_titr_residue_dict(titr_residue_dict, 'HIS', charge=1)
        if self.kbp_result is not None:
            # self.__adjust_prot_to_ph(titr_residue_dict, 5)
            self.__adjust_prot_to_ph(titr_residue_dict, 5, restrict=('ASP', None))
            self.__adjust_prot_to_ph(titr_residue_dict, 5, restrict=('GLU', None))
            self.__adjust_prot_to_ph(titr_residue_dict, 5, restrict=('NTE', None))
            self.__adjust_prot_to_ph(titr_residue_dict, 5, restrict=('CTE', None))
            self.__adjust_prot_to_ph(titr_residue_dict, 5, restrict=('CYS', None))
            self.__adjust_prot_to_ph(titr_residue_dict, 5, restrict=('HIS', None))
        # Check that there are no double entries in titr_residue_dicts
        if titr_residue_dict not in titr_residue_dicts:
            titr_residue_dicts.append(titr_residue_dict)
            names.append(name)
        else:
            print("Skipping %s since it offers no new protonation state." % name)

        # Neutral HIS
        name = 'ph7'
        titr_residue_dict = copy_titr_residue_dict(titr_residue_dict_template)

        mod_titr_residue_dict(titr_residue_dict, 'HIS', charge=0)
        if self.kbp_result is not None:
            # self.__adjust_prot_to_ph(titr_residue_dict, 7)
            # self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('ARG', None))
            self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('LYS', None))
            self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('TYR', None))
            self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('NTE', None))
            self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('CTE', None))
            self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('CYS', None))
            self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('HIS', None))
        # Check that there are no double entries in titr_residue_dicts
        if titr_residue_dict not in titr_residue_dicts:
            titr_residue_dicts.append(titr_residue_dict)
            names.append(name)
        else:
            print("Skipping %s since it offers no new protonation state." % name)

        # Neutral ASP, GLU, charged HIS
        name = 'ph-10'
        titr_residue_dict = copy_titr_residue_dict(titr_residue_dict_template)
        mod_titr_residue_dict(titr_residue_dict, 'HIS', charge=1)
        mod_titr_residue_dict(titr_residue_dict, 'ASP', charge=0)
        mod_titr_residue_dict(titr_residue_dict, 'GLU', charge=0)
        mod_titr_residue_dict(titr_residue_dict, 'CTE', charge=0)
        # Check that there are no double entries in titr_residue_dicts
        if titr_residue_dict not in titr_residue_dicts:
            titr_residue_dicts.append(titr_residue_dict)
            names.append(name)
        else:
            print("Skipping %s since it offers no new protonation state." % name)

        # Neutral LYS, HIS, charged TYR
        name = 'ph11'
        titr_residue_dict = copy_titr_residue_dict(titr_residue_dict_template)
        mod_titr_residue_dict(titr_residue_dict, 'HIS', charge=0)
        mod_titr_residue_dict(titr_residue_dict, 'LYS', charge=0)
        mod_titr_residue_dict(titr_residue_dict, 'NTE', charge=0)
        mod_titr_residue_dict(titr_residue_dict, 'TYR', charge=-1)
        mod_titr_residue_dict(titr_residue_dict, 'CYS', charge=-1)
        if self.kbp_result is not None:
            self.__adjust_prot_to_ph(titr_residue_dict, 11, restrict=('HIS', 0))
        # Check that there are no double entries in titr_residue_dicts
        if titr_residue_dict not in titr_residue_dicts:
            titr_residue_dicts.append(titr_residue_dict)
            names.append(name)
        else:
            print("Skipping %s since it offers no new protonation state." % name)

        # Neutral LYS, ARG, HIS, charged TYR
        # name = 'ph20'
        # titr_residue_dict = copy_titr_residue_dict(titr_residue_dict_template)
        # mod_titr_residue_dict(titr_residue_dict, 'HIS', charge=0)
        # mod_titr_residue_dict(titr_residue_dict, 'LYS', charge=0)
        # mod_titr_residue_dict(titr_residue_dict, 'ARG', charge=0)
        # mod_titr_residue_dict(titr_residue_dict, 'NTE', charge=0)
        # mod_titr_residue_dict(titr_residue_dict, 'TYR', charge=-1)
        # if self.kbp_result is not None:
        #     self.__adjust_prot_to_ph(titr_residue_dict, 20, restrict=('HIS', 0))
        # # Check that there are no double entries in titr_residue_dicts
        # if titr_residue_dict not in titr_residue_dicts:
        #     titr_residue_dicts.append(titr_residue_dict)
        #     names.append(name)
        # else:
        #     print("Skipping %s since it offers no new protonation state." % name)



        # # Charged HIS
        # name = 'ph5_roifree'
        # titr_residue_dict = copy_titr_residue_dict(titr_residue_dict_template)
        # mod_titr_residue_dict(titr_residue_dict, 'HIS', charge=1)
        # if self.kbp_result is not None:
        #     # self.__adjust_prot_to_ph(titr_residue_dict, 5)
        #     self.__adjust_prot_to_ph(titr_residue_dict, 5, restrict=('ASP', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 5, restrict=('GLU', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 5, restrict=('NTE', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 5, restrict=('CTE', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 5, restrict=('CYS', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 5, restrict=('HIS', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 5, residue_list=self.roi)
        # # Check that there are no double entries in titr_residue_dicts
        # if titr_residue_dict not in titr_residue_dicts:
        #     titr_residue_dicts.append(titr_residue_dict)
        #     names.append(name)
        # else:
        #     print("Skipping %s since it offers no new protonation state." % name)
        #
        # # Neutral HIS
        # name = 'ph7_roifree'
        # titr_residue_dict = copy_titr_residue_dict(titr_residue_dict_template)
        #
        # mod_titr_residue_dict(titr_residue_dict, 'HIS', charge=0)
        # if self.kbp_result is not None:
        #     # self.__adjust_prot_to_ph(titr_residue_dict, 7)
        #     # self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('ARG', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('LYS', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('TYR', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('NTE', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('CTE', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('CYS', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 7, restrict=('HIS', None))
        #     self.__adjust_prot_to_ph(titr_residue_dict, 7, residue_list=self.roi)
        # # Check that there are no double entries in titr_residue_dicts
        # if titr_residue_dict not in titr_residue_dicts:
        #     titr_residue_dicts.append(titr_residue_dict)
        #     names.append(name)
        # else:
        #     print("Skipping %s since it offers no new protonation state." % name)




        # Check that residues in salt bridges are not made neutral.
        # Todo: Move this to a function
        # get the modelled structure from The pre MD Karlsberg+ calculation.
        kbp_folder = self.kbp_result.folder
        structure_filename = kbp_tools.parse_input_structure_filename(kbp_folder)
        ssp = file_parser.MDAnalysis_ssp(structure_filename)
        (res_in_sb, res_not_in_sb, salt_bridges) = ssp.get_residues_in_sb()
        # salt_bridges = []

        s = file_parser.Simple_struct_parser()
        s.read_pdb(structure_filename)
        s.sasa(estimate=True)

        buried_sb_residues = []
        for sb in salt_bridges:
            sb_residues = []
            for (res, residue) in sb:
                resname, resid, segname = re.split(r'[-_]', residue)
                resid = int(resid)
                sasa = s.struct[segname][resid].sasa
                if sasa < 1.0:
                    sb_residues.append(residue)
                else:
                    # If one residue if a salt bridge is solvent exposed, no residue of that salt bridge is excluded.
                    sb_residues = []
                    break
            for residue in sb_residues:
                print("Residue %s will not be made neutral, since the sasa values of all salt bridge partners are "
                      "smaller that 20.0 A^2." % residue)
                buried_sb_residues.append(residue)

        if buried_sb_residues:
            for residue_descr in buried_sb_residues:
                resname, resid, segname = re.split(r'[-_]', residue_descr)
                residue_tuple = (resname, int(resid), segname)
                restype = self.charmm_struct.titr_restypes[resname]
                if restype in ['ARG', 'LYS', 'HIS']:
                    charge = 1
                elif restype in ['ASP', 'GLU']:
                    charge = -1
                else:
                    error = "Residue %s is in a buried salt bridge, but no rule is defined to handle this situation"
                    raise AssertionError(error)

                for titr_residue_dict in titr_residue_dicts:
                    self.charmm_struct.mod_titr_residue_dict_single(titr_residue_dict, residue_tuple, charge=charge)

        return names, titr_residue_dicts

    def __adjust_prot_to_ph(self, titr_residue_dict, ph, restrict=None, residue_list=None):
        """Requires self.kbp_result to be set.

        @param titr_residue_dict: dict
        @param ph: int
        @param residue_list: list of residues described by residue tuple: (resname, resid, segname). Only these
        residues are modified.
        @return:
        """
        # if restrict is not None:
        #     restrict_restype, restrict_charge = restrict

        ph_index = self.kbp_result.descr.ph_values.index(ph)
        prot_state_yaml = {}
        for residue_descr in titr_residue_dict.keys():
            # residue_tuple = tuple(re.split(r'[-_]', residue_descr)[:2])

            resname, resid, segname = residue_descr
            resname_kbp = kbp_tools.get_kbp_resname(resname)
            residue_descr_kb = (resname_kbp, resid, segname)
            if not self.kbp_result.descr.has_residue(residue_descr_kb):
                continue
            kbp_index = self.kbp_result.descr.get_residue_index(residue_descr_kb)
            nr_of_states = len(self.kbp_result.occs[kbp_index])

            best_state = None
            best_occ = -1
            for state in range(nr_of_states):
                occ = self.kbp_result.occs[kbp_index][state][ph_index]
                if occ > best_occ:
                    best_state = state
                    best_occ = occ

            kbp_residue_descr = self.kbp_result.descr.sorted_residue_list[kbp_index]
            prot_state_yaml[kbp_residue_descr] = best_state

        # Copy dict to keep track of changes
        titr_residue_dict_old = {}
        for residue_tuple, state in titr_residue_dict.iteritems():
            titr_residue_dict_old[tuple(residue_tuple)] = state

        titratable_yaml = self.kbp_result.descr.titratable_residues
        self.charmm_struct.mod_titr_residue_dict_from_titrable_yaml(titr_residue_dict, prot_state_yaml,
                                                                    titratable_yaml, restrict=restrict,
                                                                    residue_list=residue_list)
        for residue_tuple in titr_residue_dict.keys():
            old_state = titr_residue_dict_old[residue_tuple]
            new_state = titr_residue_dict[residue_tuple]
            if (old_state is not None) and (old_state != new_state):
                kbp_residue_descr = kbp_tools.get_kbp_residue(residue_tuple)
                resname = residue_tuple[0]
                if resname != 'HIS':
                    print "Residue %9s has been switched from state %i to %i for pH %.1f" \
                          % (kbp_residue_descr, old_state, new_state, ph)


    def __write_prot_tasks(self, names, titr_residue_dicts):
        # filename = self.project_folder + self.pickle_prot_tasks
        filename = self.project_folder + 'prot_tasks.pickle'
        # changed for PROPKA pKas
        # filename = self.project_folder + 'prot_tasks_propka.pickle'
        f = open(filename, 'w')
        pickle.dump((names, titr_residue_dicts), f)
        f.close()


    def __read_prot_tasks(self):
        # Find latest iteration
        files_in_project_folder = os.listdir(self.project_folder)
        max_iteration = 0

        reg = re.compile(r'^prot_tasks_propka_?i?(\d*).pickle$')
        # changed for PROPKA pKas
        # reg = re.compile(r'^prot_tasks_propka_?i?(\d*).pickle$')
        last_iteration = 0
        for file in files_in_project_folder:
            reg_m = reg.match(file)
            if reg_m is not None:
                # if len(reg_m.groups()) > 0:
                iteration = reg_m.groups()[0]
                if iteration:
                    iteration = int(iteration)
                    if iteration > last_iteration:
                        last_iteration = iteration
        if max_iteration != -1:
            if last_iteration > max_iteration:
                last_iteration = max_iteration

        # Read protonation tasks
        if last_iteration != 0:
            pickle_prot_tasks = 'prot_tasks_i%i.pickle' % last_iteration
            # changed for PROPKA pKas
            # pickle_prot_tasks = 'prot_tasks_propka_i%i.pickle' % last_iteration
            # print("\nReading protonation from %s" % pickle_prot_tasks)
        else:
            pickle_prot_tasks = 'prot_tasks.pickle'
            # changed for PROPKA pKas
            # pickle_prot_tasks = 'prot_tasks_propka.pickle'
            # print("\nReading protonation from %s" % pickle_prot_tasks)

        # filename = self.project_folder + self.pickle_prot_tasks
        filename = self.project_folder + pickle_prot_tasks
        if os.path.exists(filename):
            f = open(filename, 'r')
            names, titr_residue_dicts = pickle.load(f)
            f.close()
        else:
            names, titr_residue_dicts = None, None

        return names, titr_residue_dicts

    def is_written_prot_available(self):
        names, titr_residue_dicts = self.__read_prot_tasks()
        if titr_residue_dicts is None:
            return False
        else:
            return True

    def setup_tasks(self):
        # Try to read a previous stored protocol for protonation tasks.
        names, titr_residue_dicts = self.__read_prot_tasks()
        if titr_residue_dicts is None:
            # If reading was not successful, use the standard protocol.
            names, titr_residue_dicts = self.__get_std_prot_tasks()
            names = self.check_task_names(names, titr_residue_dicts)
            self.__write_prot_tasks(names, titr_residue_dicts)

        # Create the tasks.
        for name, titr_residue_dict in zip(names, titr_residue_dicts):
            self.add_task(name, titr_residue_dict)

    def check_task_names(self, names, titr_residue_dicts):
        # Search the project folder for exiting tasks and change the name of the new tasks, if they have been calculated
        # already.
        new_names = [None] * len(names)
        new_tasks_already_calculated = 0
        folders = os.listdir(self.project_folder)
        for folder in folders:
            prot_pickel_filename = self.project_folder + folder + '/prot_task.pickle'
            if os.path.exists(prot_pickel_filename):
                f = open(prot_pickel_filename, 'r')
                name, found_titr_residue_dict = pickle.load(f)
                f.close()

                # Check if this task matches one of the new tasks
                found_residue_tuples = sorted(found_titr_residue_dict.keys())
                for i, titr_residue_dict in enumerate(titr_residue_dicts):
                    residue_tuples = sorted(titr_residue_dict.keys())
                    if residue_tuples != found_residue_tuples:
                        continue
                    for residue_tuple, state in titr_residue_dict.iteritems():
                        found_state = found_titr_residue_dict[residue_tuple]
                        if state is None:
                            state = 0
                        if found_state is None:
                            found_state = 0
                        if state != found_state:
                            break
                    else:
                        # Match found! Use the existing task name
                        new_names[i] = name
                        print("Renaming %s into %s" % (names[i], new_names[i]))

        # If no matching task has been found make sure, that the name has not been used yet.
        for i, new_name in enumerate(new_names):
            if new_name is None:
                name = names[i]
                suffix_nr = 0
                while True:
                    if suffix_nr == 0:
                        unused_name = name
                    else:
                        unused_name = name + '_' + str(suffix_nr)

                    task_folder = self.project_folder + unused_name
                    if not os.path.exists(task_folder):
                        break
                    suffix_nr += 1
                new_names[i] = unused_name
                if suffix_nr > 0:
                    print("Renaming %s into %s" % (names[i], new_names[i]))


        # # Check if a task with that name does exist. If it does compare the protonation states and use another name if
        # # they mismatch.
        # new_names = []
        # for name, titr_residue_dict in zip(names, titr_residue_dicts):
        #     suffix_nr = 0
        #     while True:
        #         if suffix_nr == 0:
        #             new_name = name
        #         else:
        #             new_name = name + '_' + str(suffix_nr)
        #
        #         task_folder = self.project_folder + new_name
        #         if os.path.exists(task_folder):
        #             filename = task_folder + '/prot_task.pickle'
        #             f = open(filename)
        #             exisiting_name, exisiting_titr_residue_dict = pickle.load(f)
        #             f.close()
        #             # print task_folder
        #             # print "to existing: " + str(exisiting_titr_residue_dict)
        #
        #             match = True
        #             for residue_tuple, existing_state in exisiting_titr_residue_dict.iteritems():
        #                 if existing_state is None:
        #                     existing_state = 0
        #                 if residue_tuple in titr_residue_dict:
        #                     state = titr_residue_dict[residue_tuple]
        #                     if state is None:
        #                         state = 0
        #                     if existing_state != state:
        #                         match = False
        #                         break
        #                 else:
        #                     match = False
        #                     break
        #
        #             if not match:
        #                 # print "Mismatch in folder %s and residue %s" % (task_folder, residue_tuple)
        #                 suffix_nr += 1
        #             else:
        #                 break
        #         else:
        #             break
        #     new_names.append(new_name)
        return new_names

class MD_kbp2_task(object):
    """A class that contains all parameters and results for a pka calculation based on the trajectory of a molecular
    dynamic simulation.
    """
    def __init__(self, project, name, titr_residue_dict):
        """Constructor of MD_kbp2_task.
        """
        # Link to the parent project.
        self.parent_project = project

        # Job name, a folder with this name will be created.
        self.taskname = name

        # A dictionary defining the protonation state of the structure. For details about the definition see:
        # charmm.get_titr_residue_dict()
        self.titr_residue_dict = self.parent_project.charmm_struct.copy_titr_residue_dict(titr_residue_dict)

        # Folder of the task.
        self.task_folder = self.parent_project.project_folder + name + '/'

        self.modelling_folder = self.task_folder + 'modelling/'


        # Prefix of the files created by modelling.
        self.modelling_prefix = ''

        # This prefix is used for the files produced by the MD simulation.
        self.md_prefix = None

        # Check if the folder already contains a modelled structure.
        # Checks: 1) Check if self.modelling_folder exists
        #            -> if it does: Don't do any modelling, just submit the job
        #            -> if it doesn't: Check if self.task_folder exists
        #               -> if it does: raise error since the folder could contain a crashed or incompatible task.
        #               -> if it doesn't: create the task folder
        if os.path.exists(self.modelling_folder):
            # Todo: check if the modelling was successful?
            modelling_competed = True
        else:
            modelling_competed = False
            if os.path.exists(self.task_folder):
                error = "Task folder %s exists, but it does not contain a modelled structure." % self.task_folder
                raise AssertionError(error)
            else:
                os.mkdir(self.task_folder)

        self.modelling_competed = modelling_competed


    def modell_structure(self):
        os.mkdir(self.modelling_folder)

        charmm_struct = self.parent_project.charmm_struct
        charmm_struct.apply_titr_residue_dict(self.titr_residue_dict)
        charmm_struct.workdir = self.modelling_folder
        charmm_struct.run_charmm()

        # Todo: Die if modelling failed.
        return


