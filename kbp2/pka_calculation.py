# coding=utf-8

__author__ = 'Tim Meyer & Jovan Dragelj'

__version__ = '0.2'


import kbp2
import os
import shutil
import re
import pickle 
from multiprocessing import Pool

class PkaCalcSettings():
    """
    Main class with all settings and input information for pKa computations.
    There are several options that will be explained below, close to the variable/object.
    Suggestion: basic usage is better explained in an example script.
    """

    def __init__(self, pdb_filename = '', workdir = None):


        ###########################
        ### Required parameters ###
        ###########################

        # PDB filepath: should be ste with set_structure
        self.pdb_filename = ''

        # PDB code, will be set automatically from the PDB file
        self.pdb_name = ''

        # working directory: can be set with a set_workdir
        self.workdir = workdir

        ###########################
        ### Optional parameters ###
        ###########################

        # a list of topology files
        self.top = []

        # a list of parameter files
        self.par = []

        # Karlsberg protocols: PACs and defined modelling procedures
        self.protocol = []
        self.protocol.append((-10, 'open_sb_acids'))
        self.protocol.append((  7, 'h_min'))
        self.protocol.append(( 20, 'open_sb_bases'))

        # log file: automatically set
        self.log_file = []

        # titratable yaml: should be set with set_yaml
        dir = os.path.dirname(__file__)
        self.titratable_yaml = dir + '/additional_files/titratable.yaml'

        # modellign decision from charmm.py
        self.modelling_decisions = []

        # if there is more than one PAC, PACs will be parallelized
        self.processes = 3

        # prinout control
        self.quiet_mode = False

        # how to submitt
        self.qsub_parameter = ''

        # some important files will be copied as a final output, suggested to be kept with a value True
        self.copy_files = True

        # options for deletion of folders upon completion:
        # 'all' - all folders ,
        # 'PACs' - PAC folders,
        # 'keep' - keep everything
        self.remove_folders = 'PACs'

        # binary for CHARMM software
        self.charmm_bin = "charmm"

        # definitions of dielectric contstants for energy computations, should not be changed!
        self.sdie = 80
        self.pdie = 4
        self.init_die80 = True

        # For titration of frames of MDs in order to combine them to calculate pka values this option should be used!
        # This is needed for conformational energy computation in order to properly complete the procedure in the paper
        # Meyer et al. This also requires the usage of a special yaml file. See examples
        self.md_evaluation_mode = False

        # filling internal cavities with water, details in the paper Meyer et al
        # list with three parameters for example: 'cavity 0.8 0.2 0.0'
        self.cavity_par = []

        # if there is only one PAC, which happens often (usually the PAC ph=7, 'h_min')
        self.force_conf_ene_calc = False

        # the choice of software for protonation energy computations
        # note: mFES should be implemented in future
        self.prot_energy_method = 'tapbs'
        # binary for tapbs software
        self.tapbs_bin = '/home/users/j/jdragelj/bin/tapbs_1.3_cav_enere'

        # the choice of software for confromational energy computations
        # note: mFES should be implemented in future
        self.conf_energy_method = 'apbs'
        self.apbs_bin    = "/home/users/j/jdragelj/bin/apbs"
        # apbs resolution for conformational energy
        self.apbs_res = 0.3
        # obsolete parameter as it is better to use CHARMM energy for coulomb component
        self.coulomb_bin = "/home/users/j/jdragelj/bin/coulomb"

        # mFES construction site - not ready!
        self.mfes_settings = {}
        self.pka_cycle0 = {}

        # if PACs ph=-10 and ph=20 are used salt-bridges will be opened.
        # The inromation about salt-bridge pairs will be stored here
        self.salt_bridges = []
        # cutoff value used to define a salt-bridge, as the distance in Ansgtroms between charged groups
        self.sb_cutoff =  4

        # Should any minimization be done when the structure is being prepared for pka computation?
        self.init_modelling_min = True
        # initial_protonation of titratable residues, should be dine via set_initial_protonation
        self.initial_protonation = {}
        # preoptimization feature from older version of Karlsberg
        self.preopt = {}
        self.preopt = { 'carb_oxi_relax' : True, \
                        'init_die_4' :  True \
                        }

        #These residues will be the only ones titrated
        self.select_residues = []

        #These residues will be excluded from titration
        self.exclude_residues = []

        # needed patches for the proper modelling, see examples
        self.patches = None
        # pacthes to be used after AUTOGENERATE ANGLES DIHEDRALS, see examples
        self.patches_no_autogen = None

        # an option to add something specific for the user or testing!
        self.tmp_settings = {}


    def set_initial_protonation(self, protonation_vector):
        """
        :param protonation_vector: dictionary containing residues and intial prtonation states according to the yamle file
        numbered from 0. example: protonation_vector = {'EPP-1_A':2, 'LYS-4_B':0}
        :return:
        """
        for residue, state in protonation_vector.iteritems():
            if type(residue) is str:
                entries = re.split(r'[-_]', residue)
                resname, resid, segname = entries
                resname = kbp2.kbp_tools.get_kbp_resname(resname)
                residue_tuple = (resname, int(resid), segname)
            else:
                residue_tuple = residue
            self.initial_protonation.update({residue_tuple:state})

    def set_yaml(self, titratable_yaml):
        """
        :param titratable_yaml: filepath to the yaml file
        :return:
        """
        self.titratable_yaml = titratable_yaml

    def set_selected_residues(self, residue_list):
        """
        :param residue_list: a list of residues that will be titrated, example: ['EPP-1_A', 'LYS-4_B']
        :return:
        """

        for residue_desc in residue_list:
            if type(residue_desc) == str:
                entries = re.split(r'[-_]', residue_desc)
                resname, resid, segname = entries
                resname = kbp2.kbp_tools.get_kbp_resname(resname)
                residue_tuple = (resname, int(resid), segname)
            else:
                (resname, resid, segname) = residue_desc
                resname = kbp2.kbp_tools.get_kbp_resname(resname)
                residue_tuple = (resname, resid, segname)

            self.select_residues.append(residue_tuple)

    def set_excluded_residues(self, residue_list):
        """
        :param residue_list: a list of residues that will not be titrated, example: ['DPP-2_A', 'TYR-8_B']
        :return:
        """

        for residue_desc in residue_list:
            if type(residue_desc) == str:
                entries = re.split(r'[-_]', residue_desc)
                resname, resid, segname = entries
                resname = kbp2.kbp_tools.get_kbp_resname(resname)
                residue_tuple = (resname, int(resid), segname)
            else:
                (resname, resid, segname) = residue_desc
                resname = kbp2.kbp_tools.get_kbp_resname(resname)
                residue_tuple = (resname, resid, segname)

            self.exclude_residues.append(residue_tuple)

    def set_structure(self, pdb_filename):
        """
        :param pdb_filename: filepath to the pdb file
        :return:
        """

        if not os.path.exists(pdb_filename):
            error = "This file does not exist: %s" % pdb_filename
            raise AssertionError(error)
        else:
            self.pdb_filename = pdb_filename

        end = self.pdb_filename.rfind('.')
        start = self.pdb_filename.rfind('/')
        if start == -1:
            start += 1
        self.pdb_name = self.pdb_filename[start+1:end]

    def set_workdir(self, workdir):
        """
        :param workdir: work directory filepath
        :return:
        """

        if workdir[-1] != '/':
            workdir += '/'
        self.workdir = workdir

    def check(self):

        files_to_check = []
        files_to_check += self.top
        files_to_check += self.par
        files_to_check.append(self.titratable_yaml)
        files_to_check.append(self.pdb_filename)

        for file in files_to_check:
            if not os.path.exists(file):
                error = "File %s does not exist!" % file
                raise AssertionError(error)

        if self.select_residues and self.exclude_residues:
            error = 'Selecting residues and excluding residues is not allowed at the same time.'
            raise AssertionError(error)


    def copy(self):

        new = PkaCalcSettings(pdb_filename=self.pdb_filename, workdir=self.workdir)

        new.pdb_filename = self.pdb_filename

        new.pdb_name = self.pdb_name

        new.workdir = self.workdir

        new.top = []
        new.top = list(self.top)

        new.par = []
        new.par = list(self.par)

        new.protocol = []
        new.protocol = list(self.protocol)

        new.quiet_mode = self.quiet_mode

        new.qsub_parameter = self.qsub_parameter

        new.titratable_yaml = self.titratable_yaml

        new.modelling_decisions = list(self.modelling_decisions)

        new.prot_energy_method = self.prot_energy_method

        new.conf_energy_method = self.conf_energy_method

        new.charmm_bin = self.charmm_bin

        new.cavity_par = []
        new.cavity_par = list(self.cavity_par)

        new.pka_cycle0 = {}
        if self.pka_cycle0:
            for residue, marks in self.pka_cycle0.iteritems():
                new.pka_cycle0[residue] = list(marks)
        new.mfes_settings = dict(self.mfes_settings)

        new.apbs_bin = self.apbs_bin

        new.coulomb_bin = self.coulomb_bin

        new.apbs_res = self.apbs_res

        new.processes = self.processes

        new.init_die80 = self.init_die80

        new.sdie = self.sdie

        new.pdie = self.pdie

        new.tapbs_bin = self.tapbs_bin

        new.copy_files = self.copy_files

        new.remove_folders = self.remove_folders

        new.log_file = self.log_file

        new.salt_bridges = []
        new.salt_bridges = list(self.salt_bridges)

        new.sb_cutoff = self.sb_cutoff

        new.preopt = {}
        new.preopt = dict(self.preopt)

        new.initial_protonation = {}
        for residue_desc, state in self.initial_protonation.iteritems():
            new.initial_protonation.update({residue_desc:state})

        if self.select_residues:
            for residue in self.select_residues:
                new.select_residues.append(residue)
        else:
            new.select_residues = list(self.select_residues)

        if self.exclude_residues:
            for residue in self.exclude_residues:
                new.exclude_residues.append(residue)
        else:
            new.exclude_residues = list(self.exclude_residues)

        if self.patches is not None:
             if type(self.patches) is list:
                 new.patches = []
                 for patch in self.patches:
                     for name, residues in patch.iteritems():
                         new_patch = {name: list(residues)}
                         new.patches.append(new_patch)
             else:
                 new.patches = self.patches
        else:
                 new.patches = self.patches

        if self.patches_no_autogen is not None:
             if type(self.patches_no_autogen) is list:
                 new.patches_no_autogen = []
                 for patch in self.patches_no_autogen:
                     for name, residues in patch.iteritems():
                         new_patch = {name: list(residues)}
                         new.patches_no_autogen.append(new_patch)
             else:
                 new.patches_no_autogen = self.patches_no_autogen
        else:
                 new.patches_no_autogen = self.patches_no_autogen

        new.init_modelling_min = self.init_modelling_min

        new.force_conf_ene_calc = self.force_conf_ene_calc

        new.md_evaluation_mode = self.md_evaluation_mode

        new.tmp_settings = {}
        for setting_name, setting in self.tmp_settings.iteritems():
            new.tmp_settings[setting_name] = setting

        return new

    def log_settings(self):

        filename = self.workdir + 'log.dat'
        log_file = open(filename, 'w')
        log_file.write('Log file for pka calculation for the structure %s \n' % self.pdb_name)
        log_file.write('\n')
        if not self.quiet_mode:
            print "Starting calculations of pKa values for protein: %s" % self.pdb_name
        log_file.write('# Setup: \n')

        files = ['pdb_file', 'workdir', 'topology', 'parameters', 'yaml', 'protocol', 'preopt', 'salt_bridge_cutoff']

        files_for_log = {}
        files_for_log['pdb_file'] = self.pdb_filename
        files_for_log['workdir'] = self.workdir
        files_for_log['topology'] = self.top
        files_for_log['parameters'] = self.par
        files_for_log['yaml'] = self.titratable_yaml
        files_for_log['protocol'] = self.protocol
        files_for_log['preopt'] = self.preopt
        files_for_log['salt_bridge_cutoff'] = self.sb_cutoff
        if self.cavity_par:
            files_for_log['cavity_parameters'] = self.cavity_par
            files.append('cavity_parameters')

        for file in files:
            data = file
            info = files_for_log[data]
            log_file.write('-%s:' % data)
            if type(info) is list:
                for item in info:
                    log_file.write('\n %s' % str(item))
            elif type(info) is dict:
                for key, value in info.iteritems():
                    log_file.write('\n %s: %s' % (str(key), str(value)))
            else:
                log_file.write('%s' % info)
            log_file.write('\n')
        log_file.write('\n')
        log_file.close()

    def set_log_file(self, log_folder=None, filename='log.dat'):

       if log_folder is None:
           log_folder = self.workdir

       if log_folder[-1] != '/':
           log_folder += '/'

       self.log_file = log_folder + filename

    def log(self, info, quiet=True):

        ''' Creates a short log of what has been done'''

        filename = self.log_file

        # check if this is needed
        if not os.path.exists(filename):
            log_file = open(filename, 'w')
        else:
            log_file = open(filename, 'a')
        if str(info) == '\n':
            log_file.write('%s' % str(info))
        else:
            log_file.write('%s \n' % str(info))

        if not quiet:
            quiet = self.quiet_mode
            if not quiet:
                print info


        log_file.close()


def prot_energy_karlsberg_run(iteration, kbp2_settings, work_folder, jobname, charmm_struct, residues_to_titrate, titratable_residues, ph = [-10, 20, 0.5]):
    """
    Get protonation from tapbs / karlsberg(ph) run
    """

    if not os.path.exists(work_folder):
        os.mkdir(work_folder)

    method = kbp2_settings.prot_energy_method

    prot_energy_folder = work_folder + method + '/'
    if not os.path.exists(prot_energy_folder):
        charmm_ssp = charmm_struct.structure
        if not charmm_ssp.par_read:
            for par in charmm_struct.par:
                charmm_ssp.read_par(par)

        info_to_log = "Starting prot_energy job to determine protonation"
        kbp2_settings.log(info_to_log)

        cavity_par = ''
        if kbp2_settings.cavity_par:
            cavity_par = 'cavity %.2f %.2f %.1f' % (kbp2_settings.cavity_par[0], kbp2_settings.cavity_par[1], kbp2_settings.cavity_par[2])

        if kbp2_settings.init_die80 and (iteration == 0):
            pdie = 80
            sdie = 80
        else:
            pdie = kbp2_settings.pdie
            sdie = kbp2_settings.sdie
        if method == 'tapbs':
            if kbp2_settings.md_evaluation_mode:
                tapbs_bin = kbp2_settings.tapbs_bin
                kbp2.tapbs.run_tapbs(prot_energy_folder, jobname, residues_to_titrate, charmm_ssp, titratable_residues, tapbs_bin, pdie=pdie, sdie=sdie, cavity_par=cavity_par, md_evaluation_mode=True)
            else:
                tapbs_bin = kbp2_settings.tapbs_bin
                kbp2.tapbs.run_tapbs(prot_energy_folder, jobname, residues_to_titrate, charmm_ssp, titratable_residues, tapbs_bin, pdie=pdie, sdie=sdie, cavity_par=cavity_par)
        elif method == 'mfes':
            # needs to be implemented!
            mfes_settings = kbp2_settings.mfes_settings
            pka_cycle0 = kbp2_settings.pka_cycle0
            kbp2.mfes.calc_prot_ener(prot_energy_folder, jobname, residues_to_titrate, charmm_ssp, titratable_residues, mfes_settings, pka_cycle0, pdie=pdie, sdie=sdie)
    else:
        info_to_log = 'prot_energy folder exists, results will be used from it'
        kbp2_settings.log(info_to_log)



    pkint_filename = prot_energy_folder + jobname + '.pkint'
    g_filename = prot_energy_folder + jobname + '.g'

    info_to_log = "Starting Karlsberg"
    kbp2_settings.log(info_to_log)
    kb_workdir = work_folder + 'karlsberg/'
    if not os.path.exists(kb_workdir):
        if type(ph) is list:
            ph_range = ph
        elif type(ph) is int:
            ph_range = [ph, ph, 0.5]

        kbp2.karlsberg.run_karlsberg(pkint_filename, g_filename, folder=kb_workdir, ph_range = ph_range)

    kbp_results = kbp2.karlsberg.parse_karlsberg_results(kb_workdir)
    kbp_results.descr.titratable_residues = titratable_residues
    kbp_results.find_pkas()

    # Store g and pkint file, it is needed to combine the PACs.
    pkaint, g = kbp2.kbp_tools.parse_g_pkint(pkint_filename, g_filename, kbp_results.descr.residue_list_ext)

    kbp_results.pkaint = pkaint
    kbp_results.g = g

    return kbp_results


def  conf_energy_calc(kbp2_settings, conf_ene_folder, charmm_pac, prot_pattern=None):

    '''Runs apbs, returns total conf energy, used internally!!! '''

    conf_energy_method = kbp2_settings.conf_energy_method

    if not os.path.exists(conf_ene_folder):
        os.mkdir(conf_ene_folder)

    pdb_name = kbp2_settings.pdb_name
    jobname = pdb_name

    info_to_log = "For this calculation coulomb energy from CHARMM will be used!"
    kbp2_settings.log(info_to_log)

    if prot_pattern is not None:
        info_to_log = "Protonation pattern found"
        kbp2_settings.log(info_to_log)

        conf_ener_modelling_folder = conf_ene_folder + 'modelling/'
        if not os.path.exists(conf_ener_modelling_folder):
            os.mkdir(conf_ener_modelling_folder)

        charmm_pac.backup_settings()
        charmm_pac.apply_titr_residue_dict(prot_pattern)

        #md_evaluation_mode requires special_charmm_settings!
        if kbp2_settings.md_evaluation_mode:
            charmm_commands = 'BOMBlev -4\n'
            charmm_pac.add_charmm_command(charmm_commands, 'read_structure_sequence')
            charmm_commands = 'BOMBlev 0\n'
            if charmm_pac.charmm_instructions['patches_no_autogen']:
               charmm_pac.add_charmm_command(charmm_commands, 'patches_no_autogen')
            elif charmm_pac.charmm_instructions['patches']:
               charmm_pac.add_charmm_command(charmm_commands, 'autogen')

        if not charmm_pac.charmm_instructions['minimize']:
            charmm_pac.charmm_instructions['minimize'] = []

        charmm_pac.workdir = conf_ener_modelling_folder

        if not charmm_pac.check_structures(quiet=True):
            charmm_pac.check_structures(quiet=True)

        charmm_pac.charmm_instructions['do_minimize'] = False

        # # copy files needed for mfes
        # if kbp2_settings.conf_energy_method == 'mfes':
        #     raise AssertionError("mFES must be checked for new method of calculation!")
        #     mfes_input_files = "/scratch/scratch/jdragelj/mfes_files/"
        #     needed_folders = os.listdir(mfes_input_files)
        #     for filename in needed_folders:
        #         shutil.copy2(mfes_input_files +'/'+filename, conf_ener_modelling_folder+filename)
        #     charmm_pac.add_charmm_command('mfes sele all .and. .not. resname tip3 end', adj_task='hbuild')

        charmm_pac.add_charmm_command('ENERgy EPS 4.0 E14Fac 0.0', adj_task='hbuild')

        charmm_pac.run_charmm()

        charmm_pac.update_charges_with_modelled_structure()
        structure = charmm_pac.structure

        charmm_pac.restore_settings()
        energies = kbp2.charmm.get_charmm_energy(conf_ener_modelling_folder)

    else:
        print('WARNING: CHARMM electrostatic term will not be calculated')
        # if kbp2_settings.conf_energy_method == 'mfes':
        #     print('mFES calculation will not proceed!')
        #     raise AssertionError("mFES calculation will not be performed!")
        structure = charmm_pac.structure

    if not structure.par_read:
        for par in charmm_pac.par:
            structure.read_par(par)
    run_folder = conf_ene_folder + 'run/'

    ####################
    ####### APBS #######
    ####################
    if  conf_energy_method == 'apbs':
        apbs_bin    = kbp2_settings.apbs_bin
        coulomb_bin = kbp2_settings.coulomb_bin
        qsub_parameter = kbp2_settings.qsub_parameter
        target_res = kbp2_settings.apbs_res

        if not structure.par_read:
            for par in charmm_pac.par:
                structure.read_par(par)

        kbp2.apbs.start_apbs_job(jobname, run_folder, structure, apbs_bin, coulomb_bin, target_res=target_res,
                                 verbose=False, qsub_parameter=qsub_parameter, cavity_par=kbp2_settings.cavity_par)

        results = kbp2.apbs.read_result(run_folder)
        if len(results) == 1:
            kbp2_settings.log("No output files found. Job may not be finished yet, come back later.")
        else:
            if None in results:
                kbp2_settings.log("Job not finished yet or either the solvation or coulomb calculation failed!")
            else:
                kbp2_settings.log("Everything is fine! Here are the results:\n")
                solvation_energy, coulomb_energy = results
                print "coulomb changed - CHARMM for MDs, APBS for crystal"
                if kbp2_settings.md_evaluation_mode:
                    coulomb_energy = energies['ELEC'][0]

    ####################

    # ####################
    # ####### mFES #######
    # ####################
    # elif kbp2_settings.conf_energy_method == 'mfes':
    #         mfes_filename = conf_ener_modelling_folder +  "result.out"
    #         results_file = open(mfes_filename, 'r')
    #         for line in results_file:
    #             line = line.split(" ")
    #             solvation_energy = float(line[0])
    #             coulomb_energy = energies['ELEC'][0]
    # ####################

    else:
        error = 'Unknown conformation energy method!'
        raise AssertionError(error)

    ############################
    ### CONFORMATION ENERGY ####
    ############################
    conformation_energy = coulomb_energy + solvation_energy
    return conformation_energy


def initial_modelling(kbp2_settings, initial_modelling_folder, titratable_residues):
    """
    Used internally, preparation of the initial structure for KBP2 computation.
    """

    info_to_log = '# Starting inital modelling'
    kbp2_settings.log(info_to_log, quiet=False)

    pdb_filename = kbp2_settings.pdb_filename
    workdir = kbp2_settings.workdir
    top = kbp2_settings.top
    par = kbp2_settings.par
    decisions = kbp2_settings.modelling_decisions
    patches = kbp2_settings.patches
    patches_no_autogen = kbp2_settings.patches_no_autogen

    if workdir[-1] != '/':
        workdir += '/'
    if not os.path.exists(workdir):
        os.mkdir(workdir)

    if not os.path.exists(initial_modelling_folder):
        os.mkdir(initial_modelling_folder)

    charmm_struct = kbp2.charmm.Charmm_manager(workdir=initial_modelling_folder, charmm_bin=kbp2_settings.charmm_bin, top=top, par=par)
    charmm_ssp = kbp2.file_parser.Simple_struct_parser()
    charmm_ssp.read_pdb(pdb_filename)

    ####################
    ##### RENAMING #####
    ####################

    for residue in charmm_ssp.struct.iter_residues():
        rename_residue = True
        resname = residue.resname
        resid = residue.resid
        segname = residue.segname
        # this will return HSP, DPP, EPP and all else in regular form
        kbp_resname = kbp2.kbp_tools.get_kbp_resname(resname)
        # this residue tuple will always be kbp compatible regarldes of residue
        kbp_residue_tuple = (kbp_resname, resid, segname)

        # set_initial_prot_histdines = False
        # if kbp2_settings.exclude_residues:
        #     if kbp_residue_tuple not in kbp2_settings.exclude_residues:
        #         set_initial_prot_histdines = True
        # if kbp2_settings.select_residues:
        #     if kbp_residue_tuple in kbp2_settings.select_residues:
        #         set_initial_prot_histdines = True
        # if set_initial_prot_histdines:
        #         if resname == 'HSD':
        #             if not kbp2_settings.md_evaluation_mode:
        #                 kbp2_settings.initial_protonation.update({kbp_residue_tuple:0})
        #             else:
        #                 kbp2_settings.initial_protonation.update({kbp_residue_tuple:1})
        #         if resname == 'HSP':
        #             if not kbp2_settings.md_evaluation_mode:
        #                 kbp2_settings.initial_protonation.update({kbp_residue_tuple:1})
        #             else:
        #                 kbp2_settings.initial_protonation.update({kbp_residue_tuple:2})
        #         if resname == 'HSE':
        #             if not kbp2_settings.md_evaluation_mode:
        #                 kbp2_settings.initial_protonation.update({kbp_residue_tuple:2})
        #             else:
        #                 kbp2_settings.initial_protonation.update({kbp_residue_tuple:3})

        if kbp2_settings.exclude_residues:
            if kbp_residue_tuple in kbp2_settings.exclude_residues:
                rename_residue = False
        if kbp2_settings.select_residues:
            if kbp_residue_tuple not in kbp2_settings.select_residues:
                rename_residue = False
        if rename_residue:
            residue.rename(kbp_resname)

    charmm_struct.add_structure(charmm_ssp)

    if decisions is not None:
        for decision_name, decision in decisions:
            charmm_struct.add_decision(decision_name, decision)

    if patches is not None:
        for patch in patches:
            patch_name = patch.keys()[0]
            residues = patch[patch_name]
            charmm_struct.add_patch(patch_name, residues)

    if patches_no_autogen is not None:
        for patch in patches_no_autogen:
            patch_name = patch.keys()[0]
            residues = patch[patch_name]
            charmm_struct.add_patch(patch_name, residues, no_autogen=True)

    #################################################
    ###########  Karlsberg compatibility ############
    #################################################

    template_state = {'charge' : 0,
                      'patch' : None,
                      'external_patches' : None,
                      'rename' : None,
                      'special' : None}
    titr_residues = {}

    for restype in titratable_residues:
        titr_residues[restype] = []
        for state in titratable_residues[restype]:
            new_state = dict(template_state)
            if 'patch' in state:
                new_state['patch'] = state['patch']
            if 'rename' in state:
                new_state['rename'] = state['rename']
            if 'external_patches' in state:
                new_state['external_patches'] = []
                for i, external_patch in enumerate(state['external_patches']):
                    external_patch_name, residues_in_patches = external_patch
                    new_residue_list = []
                    for residue_tuple in residues_in_patches:
                        new_residue_list.append(residue_tuple)
                    new_state['external_patches'].append((external_patch_name, new_residue_list))
            titr_residues[restype].append(new_state)

    charmm_struct.set_titr_residues(titr_residues)

    if kbp2_settings.initial_protonation:
        for residue_desc, state in  kbp2_settings.initial_protonation.iteritems():
            residue_tuple = residue_desc
            charmm_struct.set_prot_residue(residue_tuple, state=state)

    charmm_struct.check_structures(quiet=True)
    charmm_struct.backup_settings()

    preopt = kbp2_settings.preopt
    for task, value in preopt.iteritems():
       if value:
           if task == 'carb_oxi_relax':
               kbp2_settings.log(" -Preoptimization: Minimization of carobxylic oxigens will be done", quiet=True)
               charmm_struct.charmm_instructions['minimize_selections'] = ['type OE*', 'type OD*']
           if task == 'init_die_4':
               kbp2_settings.log(" -Preoptimization: Dielectric constant will be set to 4 ", quiet=True)
               charmm_struct.add_charmm_command('NBONDS EPS 4', adj_task='hbuild')

    if not kbp2_settings.init_modelling_min:
        charmm_struct.charmm_instructions['do_minimize'] = False

    #special feature for jovan start
    if kbp2_settings.tmp_settings:
        if 'bomblev_flip_charge' in kbp2_settings.tmp_settings.keys():
            if kbp2_settings.tmp_settings['bomblev_flip_charge']:
                charmm_struct.add_charmm_command('bomb -4', adj_task='patch_disu')
                charmm_struct.add_charmm_command('bomb 0', adj_task='patches_no_autogen')
    # special feature for jovan end

    charmm_struct.run_charmm()
    charmm_struct.restore_settings()
    charmm_struct.structure = charmm_struct.get_modelled_structure()

    initial_pdb_filename = workdir + charmm_struct.title + '_init.pdb'
    charmm_struct.structure.write_pdb(initial_pdb_filename)
    shutil.copy2(initial_modelling_folder+'modelling_decision.dat', workdir+'modelling_decision.dat')

    return charmm_struct

def pac_modelling(kbp2_settings, modelling_type, charmm_pac, protonation, iteration_folder, iteration):
    """
    :param modelling_type: definition of the modellign that shiuld be done for this PAC/protocol
    This allows implementation of new protcols in protocols.py module
    (for example, another protocol for histidines at ph=5 etc.)
    """

    info_to_log = ("# Starting PAC modelling")
    kbp2_settings.log(info_to_log)

    if not os.path.exists(iteration_folder):
        os.mkdir(iteration_folder)

    # give charmm new protonation pattern
    charmm_pac.apply_titr_residue_dict(protonation)

    if not charmm_pac.structures_checked:
        charmm_pac.check_structures(quiet=True)

    charmm_pac.workdir = iteration_folder

    methodToCall = getattr(kbp2.protocols, modelling_type)
    methodToCall(charmm_pac, kbp2_settings, iteration)


def get_protonation(titr_residue_dict, kbp_results, ph):
    """
    After the titration, results are used to determine protonation of all residues, usually for the purposes of finding
    the best protonation pattern. Used internally.
    """
    new_titr_residue_dict = kbp2.charmm.Charmm_manager.copy_titr_residue_dict(titr_residue_dict)

    for (residue_num, residue_kbp, resname_kbp, nr_of_states) in kbp_results.descr.residue_list_ext:
        ph_nr = kbp_results.descr.ph_values.index(ph)

        best_occ = 0.0
        best_state = None

        for state in range(nr_of_states):
            occ = kbp_results.occs[residue_num][state, ph_nr]
            if occ > best_occ:
                best_occ = occ
                best_state = state

        resname, resid, segname = re.split(r'[-_]', residue_kbp)
        residue_tuple = (resname, int(resid), segname)
        new_titr_residue_dict[residue_tuple] = best_state

    return new_titr_residue_dict

def get_his_protonation(kbp_results, ph):
    """
    Utility function that will determine exact histidine protonation according to the titration results
    :param kbp_results: kbp2.kbp_results.FrameKbpResults() object created in the calc_pkas.py
    :param ph: ph value of interest (likely used during the titration)
    :return: a printout of histidine states
    """

    for (residue_num, residue_kbp, resname_kbp, nr_of_states) in kbp_results.descr.residue_list_ext:

        ph_nr = kbp_results.descr.ph_values.index(ph)
        for state in range(nr_of_states):
            # histidine protonation determination
            if resname_kbp == 'HSP':
                print 'HSP', 'resid', residue_num, 'state', state, 'occ:', kbp_results.occs[residue_num][state, ph_nr]
        if resname_kbp == 'HSP':
            print '-----------------------'

def pac_run(kbp2_settings, charmm_structure, pac_folder, ph, modelling_type, residues_to_titrate, titratable_residues, current_protonation):
    """
    Internal function that involves modelling and titration of a specific PAC/protocol
    """
    # Copy of charmm object
    charmm_pac = charmm_structure.copy()
    kbp2_settings_pac = kbp2_settings.copy()

    filename = '%i_%s_pac.log' % (ph, modelling_type)
    kbp2_settings_pac.set_log_file(pac_folder, filename)
    pdb_name = kbp2_settings_pac.pdb_name

    if pac_folder[-1] != '/':
        pac_folder += '/'

    pickle_results_kbp2 = pac_folder + 'kbp2_res.pkl'
    pickle_results_charmm_pac = pac_folder + 'charmm_pac.pkl'

    if not os.path.exists(pac_folder):
        os.mkdir(pac_folder)
    else:
        if os.path.exists(pickle_results_kbp2) and os.path.exists(pickle_results_charmm_pac):
            info_to_log =  'This pac %i_%s was already done for %s' % (ph, modelling_type, pdb_name)
            kbp2_settings_pac.set_log_file()
            kbp2_settings_pac.log(info_to_log)
            # results are added to the list and will be included in combined results
            kbp_results = pickle.load( open( pickle_results_kbp2, "rb" ) )
            charmm_pac = pickle.load( open( pickle_results_charmm_pac, "rb" ) )
            return kbp_results, charmm_pac, ph
        else:
            kbp2_settings_pac.set_log_file(pac_folder, filename)
            error = 'Something went wrong, folder exists but there are not results!!!'
            kbp2_settings_pac.log(error)
            raise AssertionError(error)


    ##################################
    ### Initial energy calculation ###
    ##################################

    jobname = 'init_calc'
    elec_calc_folder = pac_folder + jobname + '/'

    kbp2_settings_pac.log('# Initial energy calculation for this PAC!')
    kbp_results = prot_energy_karlsberg_run(0, kbp2_settings_pac, elec_calc_folder, jobname, charmm_pac, residues_to_titrate,
                                            titratable_residues, ph=[-10, 20, 0.5])
    kbp2_settings_pac.log('\n')

    ############################
    ### Self consistent loop ###
    ############################

    previous_protonation = kbp2.charmm.Charmm_manager.copy_titr_residue_dict(current_protonation)
    next_to_previous_protonation = None

    info_to_log = "# Starting self-consistent-cycle!"
    kbp2_settings_pac.log(info_to_log)
    iteration = 0
    while True:

        iteration += 1
        kbp2_settings_pac.log('\n')
        info_to_log = "***Starting iteration %i***" % iteration
        kbp2_settings_pac.log(info_to_log)


        # Find the most likely protonation state
        current_protonation = get_protonation(previous_protonation, kbp_results, ph)
        kbp2_settings_pac.log('previous_protonation')
        kbp2_settings_pac.log(previous_protonation)
        kbp2_settings_pac.log('current_protonation')
        kbp2_settings_pac.log(current_protonation)

        # Check why is the iteration started
        changed_prot_resids = []
        for residue, state in previous_protonation.iteritems():
            if current_protonation[residue] != state:
                changed_prot_resids.append(residue)
            else:
                continue
        if changed_prot_resids:
            kbp2_settings_pac.log('Residues that have changed the protonation state are:')
            kbp2_settings_pac.log(changed_prot_resids)


        # Compare with previous protonation
        stop_iteration = False
        if current_protonation == previous_protonation:
            stop_iteration = True
        elif (next_to_previous_protonation is not None) and (current_protonation == next_to_previous_protonation):
            stop_iteration = True
        if stop_iteration and iteration != 1 :
            kbp2_settings_pac.log('\n')
            kbp2_settings_pac.set_log_file()
            kbp2_settings_pac.set_log_file(pac_folder, filename)
            kbp2_settings_pac.log("Pac %i_%s completed!" % (ph, modelling_type))

            if kbp2_settings_pac.copy_files:
                workdir = kbp2_settings_pac.workdir
                folder = iteration_folder + 'tapbs/' + jobname
                shutil.copy2(folder + '.pqr', workdir + jobname + '_last.pdb')
                shutil.copy2(folder + '.reference.pqr', workdir + jobname + '_last_reference.pdb')
                shutil.copy2(pac_folder + filename, workdir + filename)
            break

        next_to_previous_protonation = previous_protonation
        previous_protonation = current_protonation

        # Set separate folder for each iteration
        iteration_folder = pac_folder + 'iter' + str(iteration) + '/'

        #####################
        ### PAC modelling ###
        #####################
        pac_modelling(kbp2_settings_pac, modelling_type, charmm_pac, previous_protonation, iteration_folder, iteration)

        ##########################
        ### Energy calculation ###
        ##########################

        #TAPBS AND KARSLBERG#
        jobname = str(ph) + '_' + str(iteration)

        kbp_results = prot_energy_karlsberg_run(iteration, kbp2_settings_pac, iteration_folder, jobname, charmm_pac, residues_to_titrate,
                                          titratable_residues, ph=[-10, 20, 0.5])
        kbp2_settings_pac.log('\n')

    ##########################
    ### kbp_results pickle ###
    ##########################
    f = open(pickle_results_kbp2, 'wb')
    pickle.dump(kbp_results, f)
    f.close()

    f = open(pickle_results_charmm_pac, 'wb')
    pickle.dump(charmm_pac, f)
    f.close()

    return kbp_results, charmm_pac


def calc_pkas(kbp2_settings, titratable_residues):
    """
    Runs pka computation on a PDB structure.
    :param kbp2_settings: kbp2.pka_calculation.PkaCalcSettings() object with titration settings
    :param titratable_residues: definitions of states for all titratable residues by type, obtained by doing:
           kbp2.kbp_tools.parse_titratable_yaml(titratable_yaml)
    :return: results.dat file with pka values
    """

    kbp2_settings.check()

    workdir = kbp2_settings.workdir
    if workdir is None:
        error = "Workdir is not provided."
        raise AssertionError(error)
    if not os.path.exists(workdir):
        os.mkdir(kbp2_settings.workdir)

    kbp2_settings.log_settings()
    kbp2_settings.set_log_file()

    folders_to_remove = []

    protocols_unfinished = False

    protocol = kbp2_settings.protocol

    for ph, modelling_type in protocol:

        pac_subfolder = 'pac_ph_' + str(ph) + '_' + modelling_type
        pac_folder = workdir + pac_subfolder + '/'
        pickle_results_kbp2 = pac_folder + 'kbp2_res.pkl'

        if os.path.exists(pickle_results_kbp2):
            info_to_log =  'Pac %i_%s was already done for %s' % (ph, modelling_type, modelling_type)
            kbp2_settings.set_log_file()
            kbp2_settings.log(info_to_log)
        else:
            protocols_unfinished = True


    #Check if the job was already done
    pickle_results = workdir + 'results.pkl'
    if os.path.exists(pickle_results) and not protocols_unfinished:
        kbp2_settings.log("*********** CONGRATULATIONS ****************")
        kbp2_settings.log("This calculation has been completed in total!", quiet=False)

        numerical_results = workdir + 'results.dat'
        if os.path.exists(workdir + 'results.dat'):
            kbp2_settings.log("Results can be read from data file!", quiet=False)
            kbp2_settings.log(("Results are -> %s" % numerical_results), quiet=False)
        else:
            kbp2_settings.log("Results can be read from pickle file!", quiet=False)
            kbp2_settings.log(("Results are -> %s" % pickle_results), quiet=False)

            combined_results_done = pickle.load( open( pickle_results, "rb" ) )
            write_pka_results(combined_results_done, kbp2_settings)
        return
    else:
        if os.path.exists(pickle_results):
            os.remove(pickle_results)

    ###################################
    ### Initial structure modelling ###
    ###################################

    initial_modelling_folder = workdir + 'initial_modelling/'
    if kbp2_settings.remove_folders == 'all':
        folders_to_remove.append(initial_modelling_folder)
    charmm_structure = initial_modelling(kbp2_settings, initial_modelling_folder, titratable_residues)

    kbp2_settings.log("Structure preparation (initial modelling) is done!", quiet=False)
    kbp2_settings.log('\n')

    current_protonation = charmm_structure.get_titr_residue_dict()
    # current_protonation = charmm_structure.get_titr_residue_dict(all=True)

    residues_not_to_titrate = []

    if kbp2_settings.select_residues:
        for res_tuple, state in current_protonation.iteritems():
            if res_tuple not in kbp2_settings.select_residues:
                residues_not_to_titrate.append(res_tuple)
    if kbp2_settings.exclude_residues:
        for res_tuple, state in current_protonation.iteritems():
            if res_tuple in kbp2_settings.exclude_residues:
                residues_not_to_titrate.append(res_tuple)

    for res_tuple in residues_not_to_titrate:
        current_protonation.pop(res_tuple)

    residues_to_titrate = charmm_structure.get_titr_residue_list(current_protonation)

    residues_to_titrate_tuples = []
    for residue_to_titrate in residues_to_titrate:
        entries = re.split(r'[-_]', residue_to_titrate)
        resname, resid, segname = entries
        residue_tuple = (resname, int(resid), segname)
        residues_to_titrate_tuples.append(residue_tuple)
    residues_to_titrate = residues_to_titrate_tuples

    # Define the protonation pattern used for conformational energy calculations
    ref_protonation = {}
    for residue_tuple, state in current_protonation.iteritems():
        if residue_tuple in residues_to_titrate:
            ref_protonation[residue_tuple] = 0


    #################################################
    ### ph adopted conformations (PAC) generation ###
    #################################################

    pac_kbp_results = []


    if len(protocol) == 1:
        kbp2_settings.processes = 1

    processes = kbp2_settings.processes
    if processes != 1:
        pool = Pool(processes=processes)

    pac_results = []
    for ph, modelling_type in protocol:

        info_to_log = ("# Pac run for %i %s" % (ph, modelling_type))
        kbp2_settings.log(info_to_log, quiet=False)

        pac_subfolder = 'pac_ph_' + str(ph) + '_' + modelling_type
        pac_folder = workdir + pac_subfolder + '/'
        if kbp2_settings.remove_folders in ['all', 'PACs']:
            folders_to_remove.append(pac_folder)

        pickle_results_kbp2 = pac_folder + 'kbp2_res.pkl'
        if os.path.exists(pickle_results_kbp2):
            continue

        # PAC run
        if processes != 1:
            pac_result = pool.apply_async(pac_run, (kbp2_settings, charmm_structure, pac_folder, ph, modelling_type,
                                                    residues_to_titrate, titratable_residues, current_protonation))
        else:
            pac_result = pac_run(kbp2_settings, charmm_structure, pac_folder, ph, modelling_type, residues_to_titrate,
                                 titratable_residues, current_protonation)

        pac_results.append((pac_folder, pac_result, ph))

    for pac_folder, pac_result, ph in pac_results:

        if processes != 1:
            (kbp_results, charmm_pac) = pac_result.get()
        else:
            (kbp_results, charmm_pac) = pac_result


        #Conformational energy calculation
        if len(protocol) == 1:
           calc_conf_energy = False
        else:
            calc_conf_energy = True

        # even if there is one protocol this flag allows conformational energy computations, eg. frame weighting
        if kbp2_settings.force_conf_ene_calc:
            calc_conf_energy = True

        if calc_conf_energy:
            print "# Starting conformational energy calculation for pH %i" % ph
            conf_energy_method = kbp2_settings.conf_energy_method
            conf_folder = pac_folder + conf_energy_method +'/'
            kbp2_settings.set_log_file(pac_folder, filename= 'conf_energy.dat')

            conformation_energy = conf_energy_calc(kbp2_settings, conf_folder, charmm_pac, prot_pattern=ref_protonation)
            kbp_results.conf_energy = conformation_energy

            info_to_log = "Conformational energy calculation done in %s" % (conf_folder)
            kbp2_settings.log(info_to_log)

        else:
            info_to_log = "Conformational energy will not be calculated."
            kbp2_settings.log(info_to_log)
            kbp_results.conf_energy = 0.0


        pac_kbp_results.append(kbp_results)


    ########################
    ### Combine the PACs ###
    ########################

    combined_results = kbp2.kbp_results.FrameKbpResults()
    for kbp_results in pac_kbp_results:
        combined_results.add_task(kbp_results)

    combined_results.combine_frames_karlsberg(cpus=processes)

    kbp2_settings.set_log_file()
    kbp2_settings.log("Final results are combined, the pKa computation is done!", quiet=False)

    ##################################
    ### Pickling of important data ###
    ##################################

    f = open(pickle_results, 'wb')
    pickle.dump(combined_results, f)

    f.close()

    write_pka_results(combined_results, kbp2_settings)

    ######################
    ### Folder removal ###
    ######################

    if kbp2_settings.remove_folders != 'keep':
        for folder in folders_to_remove:
            shutil.rmtree(folder)

def write_pka_results(results, kbp2_settings):
    """
    :param results: kbp2.kbp_results.FrameKbpResults() object created in the calc_pkas.py
    :param kbp2_settings: kbp2.pka_calculation.PkaCalcSettings() object with titration settings
    :return: 'results.dat' file with residue pkas
    """

    folder = kbp2_settings.workdir
    final_pkas = results.combined_results.pkas

    filename = folder + 'results.dat'
    results_file = open(filename, 'w')

    sorted_data = {}
    segname_list = []
    residue_list = []

    for residue_descr, pka in final_pkas.iteritems():
        residue_list.append(residue_descr)
        resname, resid, segname = re.split(r'[-_]', residue_descr)

        if segname not in sorted_data.keys():
            sorted_data[segname] = []
            segname_list.append(segname)
        if resid not in sorted_data[segname]:
            sorted_data[segname].append(int(resid))

    segname_list = sorted(segname_list)

    for segname in segname_list:
        sorted_data[segname] = sorted(sorted_data[segname])
    for segname in segname_list:
        for resid in sorted_data[segname]:
            for residue_descr, pka in final_pkas.iteritems():
                if '-' + str(resid) + '_' + segname in residue_descr:
                    if residue_descr in residue_list:
                        if not kbp2_settings.quiet_mode:
                            print residue_descr, pka
                        residue_list.remove(residue_descr)
                        entries = re.split(r'[-_]', residue_descr)
                        resname, resid, segname = entries
                        resname = kbp2.kbp_tools.get_real_resname(resname)
                        residue_descr_real = '%s-%s_%s' % (resname, resid, segname)
                        s = " %13s: %0.2f " % (residue_descr_real, pka)
                        results_file.write(s)
                        results_file.write("\n")

def get_atom_charges(resname, states, top):

    """ states is a list of 'template_state' dicts

        input example:
        template_state = {'patch': 'patch-name',
                          'pka': 6.7,
                          'name': 'D',
                          'atoms': {}
                          'rename': 'new_name'
                          'external_patches':[('external-patch-1', [(<resname1>, <resid1>, <segnam1e>)]),
                                              ('external-patch-2', [(<resname2>, <resid2>, <segname2>)])],
                          'external_atoms':[]
                          }

        output example:
        template_state = {'patch': 'patch-name',
                          'pka': 6.7,
                          'name': 'D',
                          'atoms': {atom_name : charge}
                          'rename': 'new_name'
                          'external_patches':[('external-patch-1', [(<resname1>, <resid1>, <segnam1e>)]),
                                              ('external-patch-2', [(<resname2>, <resid2>, <segname2>)])],
                          'external_atoms':[[{atom_name : charge}],[{atom_name : charge},{atom_name : charge}]]
                          }

        @ in all states number of external patches must be the same
        @ in all states order of the same external patch type (regaring the same residue(s)) must be the same
    """
    charmm_struct = kbp2.charmm.Charmm_manager(top=top)

    titratable_definition = {}
    titratable_definition[resname] = []
    for state in states:
        atoms = {}
        state.update({'atoms':atoms})
        if 'external_patches' in state.keys():
            external_atoms = []
            state.update({'external_atoms':external_atoms})
            for i, patch_tuple in enumerate(state['external_patches']):
                patch_name, residue_list = patch_tuple
                list_of_residue_tuples = []
                for residue in residue_list:
                    if type(residue) == str:
                        entries = re.split(r'[-_]', residue)
                        resname, resid, segname = entries
                        residue_tuple = (resname, int(resid), segname)
                        list_of_residue_tuples.append(residue_tuple)
                del state['external_patches'][i]
                state['external_patches'].append((patch_name, list_of_residue_tuples))


        titratable_definition[resname].append(state)
    all_atoms = []

    nr_of_patch_atoms = 0
    for i,state in enumerate(titratable_definition.values()[0]):

        if 'rename' in state.keys():
            if ('patch' in state.keys()) or ('external_patches' in state.keys()):
                error = 'Renaming cannot be used with patches or external patches'
                raise AssertionError(error)

        #check if there are PATCHES
        patch_or_external = False
        if 'patch' in state.keys():
            patch_or_external = True
            patch_name = state['patch']
            atom_list = charmm_struct.top_content['patches'][patch_name]['atoms']
            if i != 1:
                if nr_of_patch_atoms == len(atom_list):
                    nr_of_patch_atoms = len(atom_list)
                else:
                    error = 'Number of atoms in the patches for this titratable residue are not the same'
                    raise AssertionError(error)
            for atom in atom_list:
                charge = atom['charge']
                name = atom['name']
                all_atoms.append(name)
                titratable_definition[resname][i]['atoms'][name]= float(charge)

        #check if there are EXTERNAL PATCHES
        if 'external_patches' in state.keys():
            patch_or_external = True
            for k, external_patch in enumerate(state['external_patches']):
                patch_name, patch_residues = external_patch

                titratable_definition[resname][i]['external_atoms'].append([])
                for j in range(len(patch_residues)):
                    titratable_definition[resname][i]['external_atoms'][k].append({})

                atom_list = charmm_struct.top_content['patches'][patch_name]['atoms']
                for atom in atom_list:
                    if atom['name'][0].isdigit():
                        number = ''
                        for letter in atom['name']:
                            if letter.isdigit():
                                number += letter
                            else:
                                break
                        name = atom['name'][len(number):]
                        charge = atom['charge']
                        titratable_definition[resname][i]['external_atoms'][k][int(number)-1].update({name:float(charge)})

                    else:
                        name = atom['name']
                        charge = atom['charge']
                        titratable_definition[resname][i]['external_atoms'][k][0].update({name:float(charge)})

        #if there are not patches or external patches, RESIDUE atoms are taken
        if not patch_or_external:
            #check is the is RENAME
            if 'rename' in state.keys():
                residue_to_look = state['rename']
            else:
                residue_to_look = resname

            atom_list = charmm_struct.top_content['residues'][residue_to_look]['atoms']
            for atom in atom_list:
                charge = atom['charge']
                name = atom['name']
                all_atoms.append(name)
                titratable_definition[resname][i]['atoms'][name]=float(charge)

    # number of atoms check -> external pacthes
    for i in range(len((titratable_definition[resname]))):
        if i == 0:
            continue
        for k in range(len(titratable_definition[resname][i]['external_atoms'])):
            for j in range(len((titratable_definition[resname][i]['external_atoms'][k]))):
                count = len(titratable_definition[resname][i]['external_atoms'][k][j])
                count_ref = len(titratable_definition[resname][0]['external_atoms'][k][j])
                if count != count_ref:
                    error = 'Number of atoms in the external patches for this titratable residue are not the same'
                    raise AssertionError(error)


    # number of atoms check -> patches & intersection of atom lists is made
    relevant_atoms = []
    for atom in all_atoms:
        is_everywhere = True
        is_changing = False
        charge = None
        for state in titratable_definition[resname]:
            if atom not in state['atoms']:
                is_everywhere = False
                if charge is not None:
                    charge = state['atoms'][atom]
            else:
                if charge != state['atoms'][atom]:
                    is_changing = True
        if is_everywhere and is_changing:
            relevant_atoms.append(atom)

    for i,state in enumerate(titratable_definition[resname]):
        new_atoms = {}
        for atom, charge in state['atoms'].iteritems():
            if atom in relevant_atoms:
                new_atoms[atom] = charge
        titratable_definition[resname][i]['atoms'] = dict(new_atoms)

    return titratable_definition

