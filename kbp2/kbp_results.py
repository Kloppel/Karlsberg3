# -*- coding: utf-8 -*-

import os
import re
import shutil
import tempfile

import storable

import numpy as np

import kbp_tools
import file_parser

import kbp2.cython.combine_md_tools


class KbpJobDescr(object):
    """Object that contains information about a karlsberg+ job. Does not contain results.
    """
    def __init__(self):
        # List of titrated residues, ordered by segment name and resid.
        self.sorted_residue_list = None

        # List of titrated residues, based on 'sorted_residue_list' but with additional information to allow efficient
        # iteration over the resides. Each entry contains a tuple:
        # (residue_num, residue_kbp, resname_kbp, nr_of_states)
        # with: residue_num : Index of the residue in 'sorted_residue_list'
        #       residue_kbp : Residue description within the framework of Karlsberg+ (e.g. 'ARG-12_A')
        #       resname_kbp : Residue name within the framework of Karlsberg+ (e.g. ARG or EPP)
        self.residue_list_ext = None

        # Sorted list of strings containing pH values used for the titration.
        self.ph_values = None

        # Contains information about the type of resides that where titratable in the Karlsberg+ calculation and their
        # properties. Format:
        # {<Karlsberg+ residue name> : [ state : {'name'  : str <state name e.g. 'R', 'P', 'D', 'PP'>
        #                                         'pka'   : float <energy of the residue at pH 0 in units of pKa.>
        #                                         'patch' : str <name of the CHARMM patch that creates the state.>
        #                                         'atoms' : dict {<atom name> : <charge>}
        #                                       }
        #                             ]
        #         }
        self.titratable_residues = None

        # Alternative representation of self.sorted_residue_list in residue tuples: (<resid>, <segname>)
        # Is not unique for terminal residues. In that case the reside if preferred over the termini.
        # Is set the first time self.get_residue_index is called.
        self.__sorted_residue_list_tuples = None

    def get_residue_index(self, residue):
        """
        @param residue: Can be (ASP, 12, 'A') or 'ARG-12_A'
        @return:
        """
        # Todo: Meaningful error when Residue does not exist.

        if type(residue) is str:
            return self.sorted_residue_list.index(residue)
        elif type(residue) is tuple:
            # Create self.__sorted_residue_list_tuples if required.
            if self.__sorted_residue_list_tuples is None:
                self.__sorted_residue_list_tuples = []
                for residue_descr in self.sorted_residue_list:
                    resname, resid, segname = re.split(r'[-_]', residue_descr)
                    residue_tuple = (resname, int(resid), segname)

                    self.__sorted_residue_list_tuples.append(residue_tuple)

            return self.__sorted_residue_list_tuples.index(residue)

    def has_residue(self, residue):
        """
        @param residue: Can be (ASP, 12, 'A') [may not be unique,due to termini -> the residue is preferred] or
        'ARG-12_A'
        @return:
        """

        if type(residue) is str:
            try:
                self.sorted_residue_list.index(residue)
                return 1
            except ValueError:
                return 0
        elif type(residue) is tuple:
            # Create self.__sorted_residue_list_tuples if required.
            if self.__sorted_residue_list_tuples is None:
                self.__sorted_residue_list_tuples = []
                for residue_descr in self.sorted_residue_list:
                    resname, resid, segname = re.split(r'[-_]', residue_descr)
                    residue_tuple = (resname, int(resid), segname)

                    self.__sorted_residue_list_tuples.append(residue_tuple)
            try:
                self.__sorted_residue_list_tuples.index(residue)
                return 1
            except ValueError:
                return 0

    @staticmethod
    def sort_residue_list(residue_list):
        """
        @param residue_list: list of residue descriptions (list of str) e.g. ['ASP-4_A']
        @return: a sorted copy of the residue list
        """
        residue_list_copy = list(residue_list)
        # Sort by resname.
        sorted_residue_list = sorted(residue_list_copy)
        # Sort by segid.
        sorted_residue_list = sorted(sorted_residue_list, key=lambda residue: int(re.split('[_-]', residue)[1]))
        # Sorty by segname.
        sorted_residue_list = sorted(sorted_residue_list, key=lambda residue: re.split('[_-]', residue)[2])

        return sorted_residue_list

    def set_residue_list(self, residue_list):
        """Sortes the given residue list and stores it.

        @param residue_list: list of residue descriptions (list of str) e.g. ['ASP-4_A']
        @return: a sorted copy of the residue list
        """

        sorted_residue_list = self.sort_residue_list(residue_list)
        self.sorted_residue_list = sorted_residue_list

        return sorted_residue_list

class KbpResult(object):
    """Container to store the results of a Karlsberg+ calculation.
    """

    def __init__(self, descr=None):
        """Constructor
        Parameters: descr : Object of type KbpJobDescr. The object is used directly, not copied!
        """

        # Contains general properties of the Karlsberg+ job.
        if descr is None:
            self.descr = KbpJobDescr()
        else:
            assert(isinstance(descr, KbpJobDescr))
            self.descr = descr

        # Path to the folder, the results have been read from.
        self.folder = None

        ### Karlsberg+ results ###
        # Intrinsic pKa values. (These are the pKa values for each residue, if all other residues would be in the
        # reference ('R') state. Format:
        # {str <Karlsberg+ residue description e.g. ARG-12_A> : numpy array [state : float <intrinsic pKa>]
        # }
        self.pkaint = None

        # Version 2
        # Numpy matrix that contains the interaction energy between all residues for all possible states.
        # Format:
        # fake 4-dimensional numpy array (array of arrays): [residue_num][residue2_num][state][state2]
        # with: residue_num  : Number of the first residue. This is the index of this residue in 'sorted_residue_list'.
        #       residue2_num : Number of the second residue. This is the index of this residue in 'sorted_residue_list'.
        #       state  : State of first residue. This is the index of the state in 'titratable_residues'.
        #       state2 : State of second residue. This is the index of the state in 'titratable_residues'.
        # Requirements:
        #   residue_num < residue2_num
        #   state < state2
        # All other entries are None.
        self.g = None

        # Numpy matrix that contains the occupancy of all states of all residues for each pH.
        # Format:
        # [residue_num: 2-dimensional numpy array: [state, ph_number]]
        # with: residue_num  : Number of the first residue. This is the index of this residue in 'sorted_residue_list'.
        #       state        : State of first residue. This is the index of the state in 'titratable_residues'.
        #       ph_number    : Number of the pH value. This is the index of the pH in 'ph_values'.
        self.occs = None

        # Contains the occupancy of each 'pH adopted conformation' (PAC).
        # Set if the results stored in the object are result of the combination of several PACs.
        ### old ### Format: list: [int <ph_number> : float <occupancy of the PAC>]

        # 2d numpy array
        # Format: conformer_occs[int <ph_number>, int <conformation number>] e.g.: [[0, 0.5, 1.0], [1.0, 0.5, 0.0]]
        # with: pH number : Number of the pH value. This is the index of the pH in 'ph_values'.
        self.conformer_occs = None

        ### Results from titration ###
        # List of float: pKa values extracted from the titration curves. The residue description for each entry can be
        # obtained from 'sorted_residue_list'.
        self.pkas = None

        # Titration curve for each residue, that shows the probability for being deprotonated. If the deprotonated
        # state is represented by more than one state in the Karlsberg calculation the corresponding curves have been
        # merged. The residue description for each entry can be obtained from 'sorted_residue_list'.
        # Format:
        # [residue: [numpy array: <occupancy>]]
        # with: residue: The residue description for each entry can be obtained from 'sorted_residue_list'.
        #       occupancy: Probability of the residue to be deprotonated. Each entry corresponds to one pH value. The
        #                  corresponding pH values can be found in 'ph_values'.
        self.deprot_curves = None

        # List containing the titration curves for each residue.
        # Format:
        # [residue: [state: [numpy array: <occupancy>]]]
        # with: residue: The residue description for each entry can be obtained from 'sorted_residue_list'.
        #       state: State of the residue. This is the index of the state in 'titratable_residues'.
        #       occupancy: Probability of this state. Each entry corresponds to one pH value. The corresponding pH
        #                  values can be found in 'ph_values'.
        # self.titr_curves = None

        # Dictionary containing conformational energy of the reference state. Each entry is a float or numpy.array.
        self.conf_energies = {}

        # Any other information about the calculation.
        self.misc_info = {}

        # The energy required to switch the protonation from the reference state ('R' states) to the protonation state
        # according to self.occs. The results is pH dependent. The x-axis for this data are the ph values stored in
        # self.descr.ph_values
        # Format: numpy array
        self.prot_energy = None

        # The energy required to switch the protonation from the reference state ('R' states) to the protonation state
        # according to self.occs for a single residue. The results is pH dependent. The x-axis for this data are the ph
        # values stored in self.descr.ph_values. The energy includes all interaction energies to all other residues.
        # That means that the sum of the entries for individual residues is NOT equal to self.prot_energy, since
        # interaction energies would be counted several time.
        # Format: {str: <residue description> : numpy array: <protonation energy>}]
        self.residue_prot_energy = None

        # The conformational energy for the reference state (state 'R).
        self.conf_energy = None

        # A numpy array containing the Solvent accessible surface area (SASA) for all titratable residues. The residue
        # description for each entry can be obtained from 'sorted_residue_list'.
        self.residue_sasa = None


    def read_kbp_results(self, folder_name, detailed=False, allow_unfinished_job=False):
        """Reads Karlsberg+ results from a folder: with detailed=True also .pkint and .g
         @return: 1 if successful, 0 otherwise.
        """
        if folder_name[-1] != '/':
            folder_name += '/'

        # read occupancy file
        fd = os.listdir( folder_name )
        occ_fn = None
        pkint_fn = None
        g_fn = None
        for entry in fd:
            if entry.split('.')[-1] == 'occ':
                if occ_fn is not None:
                    error = "More than one .occ files found in %s" % folder_name
                    raise AssertionError(error)
                occ_fn = entry
            elif entry.split('.')[-1] =='pkint':
                if pkint_fn is not None and detailed:
                    error = "More than one .pkint files found in %s" % folder_name
                    raise AssertionError(error)
                pkint_fn = entry
            elif entry.split('.')[-1] == 'g':
                if g_fn is not None and detailed:
                    error = "More than one .g files found in %s" % folder_name
                    raise AssertionError(error)
                g_fn = entry
        if occ_fn is None:
            if allow_unfinished_job:
                return 0
            else:
                error = "ERROR in 'read_folder': No .occ file found in " + folder_name
                raise AssertionError(error)

        occ_file = folder_name + '/' + occ_fn
        so = storable.retrieve(occ_file)

        # Parse titratable.yaml file.
        self.descr.titratable_residues = kbp_tools.parse_titratable_residues(folder_name)

        residue_list = so['bias'].keys()
        # # Sort by resname.
        # sorted_residue_list = sorted(residue_list)
        # # Sort by segid.
        # sorted_residue_list = sorted(sorted_residue_list, key=lambda residue: int(re.split('[_-]', residue)[1]))
        # # Sorty by segname.
        # sorted_residue_list = sorted(sorted_residue_list, key=lambda residue: re.split('[_-]', residue)[2])
        # self.descr.sorted_residue_list = sorted_residue_list
        sorted_residue_list = self.descr.set_residue_list(residue_list)


        ph_values = so['bias'][sorted_residue_list[0]]['300'].keys()
        ph_values = [float(x) for x in ph_values]
        ph_values.sort()
        self.descr.ph_values = ph_values

        self.occs = []
        for residue in sorted_residue_list:
            residue_entry = so['bias'][residue]['300']
            nr_of_states = len(residue_entry.values()[0]['0'])
            occs = np.zeros((nr_of_states, len(ph_values)), dtype=np.float)
            for ph_s, entry in residue_entry.iteritems():
                index = ph_values.index(float(ph_s))
                for state in range(nr_of_states):
                    occ = float(entry['0'][state])
                    occs[state, index] = occ
            self.occs.append(occs)
        self.occs = np.array(self.occs, dtype=object)

        # Generate residue_list_ext list.
        residue_list_ext = []
        for residue_num, residue_kbp in enumerate(sorted_residue_list):
            resname_kbp = residue_kbp[0:residue_kbp.find('-')]
            nr_of_states = len(self.descr.titratable_residues[resname_kbp])
            residue_list_ext.append( (residue_num, residue_kbp, resname_kbp, nr_of_states) )
        self.descr.residue_list_ext = residue_list_ext
        self.descr.sorted_residue_list = sorted_residue_list

        # Calculate pKas
        self.find_pkas()

        self.folder = folder_name

        if detailed:
            # Parse .g and .pkint files
            pkint, g = kbp_tools.parse_g_pkint(folder_name + pkint_fn, folder_name + g_fn, residue_list_ext)
            self.pkaint = pkint
            self.g = g

        return 1

    def create_deprot_curves(self):
        """
        Generates the titration curve for each residue, that shows the probability
        for being deprotonated. If the deprotonated state is represented by more
        than one state in the karlsberg calculation the corresponding curves are
        merged. The result is stored in self.deprot_curves
        """
        self.deprot_curves = []
        for (residue_num, residue_kbp, resname_kbp, nr_of_states) in self.descr.residue_list_ext:
            state_names = []
            for state in self.descr.titratable_residues[resname_kbp]:
                state_names.append(state['name'])
                if state['name'] not in ['R', '0', 'D', 'P']:
                    error = "Residue %s has state of type %s. Support for this type of residue is not " \
                            + "yet implemented." % state['name']
                    raise AssertionError(error)

            is_acid = is_base = False
            if 'D' in state_names:
                is_acid = True
            if 'P' in state_names:
                is_base = True
            if is_acid and is_base:
                error = "Residue %s has at least one 'P' and 'D' state. Support for this type of residue is not "\
                        + "yet implemented."  % resname_kbp
                AssertionError(error)
            if not is_acid and not is_base:
                error = "Residue %s has no 'P' or 'D' state. Support for this type of residue is not " \
                        + "yet implemented." % resname_kbp
                raise AssertionError(error)

            residue_occs = self.occs[residue_num]

            curve = np.zeros(len(self.descr.ph_values))
            for state in range(nr_of_states):
                state_occ = residue_occs[state]
                state_name = state_names[state]

                if is_base:
                    # Add '0' and 'R'
                    if state_name in ['R', '0']:
                        curve = np.add(curve, state_occ)
                elif is_acid:
                    # Add 'D'
                    if state_name in ['D']:
                        curve = np.add(curve, state_occ)
                else:
                    raise AssertionError("Residue %s is neither acid nor base!" % residue_kbp)

                if state_name == 'R':
                    sum_charge = 0.0
                    for charge in self.descr.titratable_residues[resname_kbp][state]['atoms'].itervalues():
                        sum_charge += abs(charge)
                    if sum_charge < 0.01 < np.sum(state):
                        error = "The reference state for residue %s with zero charges has a " % residue_kbp\
                                + "non vanishing occupancy!"
                        raise AssertionError(error)

            self.deprot_curves.append(curve)

    def find_pkas2(self):
        """
        old

        Extracts pKas from titration curves in self.deprot_curves with linear interpolation. self.create_deprot_curves
        is used to set the deprotonation curves if necessary.

        Creates a dictionary with one entry for each Residue. Each entry contains the pka.The result is stored in
        self.pkas
        """
        if self.deprot_curves is None:
            self.create_deprot_curves()

        ph_values = self.descr.ph_values

        self.pkas = {}
        for (residue_num, residue_kbp, resname_kbp, nr_of_states) in self.descr.residue_list_ext:
            deprot_curve = self.deprot_curves[residue_num]
            x1 = 0
            x2 = 0
            y1 = 0
            y2 = 0
            pka_found = False
            at_first_point = False
            for i, y in enumerate(deprot_curve):
                if y < 0.5:
                    x1 = ph_values[i]
                    y1 = y
                else:
                    x2 = ph_values[i]
                    y2 = y

                    pka_found = True
                    if i == 0:
                        at_first_point = True
                    break
            if pka_found:
                pka = (x2 - x1) / (y2 - y1) * (0.5 - y1) + x1
            else:
                if at_first_point:
                    #pka = '<-10'
                    pka = -10
                else:
                    #pka = '>20'
                    pka = 20

            self.pkas[residue_kbp] = pka

    def find_pkas(self):
        """
        Extracts pKas from titration curves in self.deprot_curves with linear interpolation. self.create_deprot_curves
        is used to set the deprotonation curves if necessary.

        Creates a dictionary with one entry for each Residue. Each entry contains the pka.The result is stored in
        self.pkas
        """
        if self.deprot_curves is None:
            self.create_deprot_curves()

        ph_values = self.descr.ph_values

        self.pkas = {}

        for (residue_num, residue_kbp, resname_kbp, nr_of_states) in self.descr.residue_list_ext:
            deprot_curve = self.deprot_curves[residue_num]
            x1 = 0
            x2 = 0
            y1 = 0
            y2 = 0
            pka_found = False
            at_first_point = False
            for i, y in enumerate(deprot_curve):
                if y < 0.5:
                    x1 = ph_values[i]
                    y1 = y
                    pka_found = False
                elif not pka_found:
                    x2 = ph_values[i]
                    y2 = y

                    pka_found = True
                    if i == 0:
                        at_first_point = True
                    pka = (x2 - x1) / (y2 - y1) * (0.5 - y1) + x1

            # if not pka_found:
            #     if at_first_point:
            #         #pka = '<-10'
            #         pka = -10
            #     else:
            #         #pka = '>20'
            #         pka = 20

            if not pka_found:
                if not at_first_point:
                    #pka = '>20'
                    pka = 20
            else:
                if at_first_point:
                    pka = -10

            self.pkas[residue_kbp] = pka

    def calc_protonation_energy(self):
        # kbp2.cython.combine_md_tools.calc_protonation_energy(self)
        kbp2.cython.analyse_md_pkas_tools.calc_protonation_energy(self)

    def calc_residue_sasa(self):
        folder = self.folder
        if folder[-1] != '/':
            folder += '/'
        out_prefix = kbp_tools.parse_sfc_prefix(folder)
        pqr_filename = folder + out_prefix + '.pqr'

        pqr = file_parser.Simple_struct_parser()
        pqr.read_pdb(pqr_filename, is_pqr=True)

        # Remove zero charge atoms
        for i in range(len(pqr.atoms)-1, -1, -1):
            atom_charge = pqr.atoms[i]['charge']
            if atom_charge == 0.0:
                pqr.atoms.pop(i)
        # Can be skipped, since the 'sasa' function and all the function it calls are iterating directly over the atom
        # list.
        # pqr.create_struct()

        # pqr.sasa(estimate=False)
        pqr.sasa(estimate=True)

        residue_sasa = np.zeros(len(self.descr.sorted_residue_list))
        for i, residue_descr in enumerate(self.descr.sorted_residue_list):
            resname, resid, segname = re.split(r'[-_]', residue_descr)
            resid = int(resid)
            residue = pqr.struct[segname][resid]

            for atom in residue.iter_atoms():
                residue_sasa[i] += atom['sasa']

        self.residue_sasa = residue_sasa

        total_sasa = 0.0
        for atom in pqr.atoms:
            total_sasa += atom['sasa']

        self.conf_energies['sasa'] = total_sasa

    def karlsberg_titration(self, folder=None, keep_files=False, cpus=1):
        if folder is not None:
            if folder[-1] != '/':
                folder += '/'

            if os.path.exists(folder):
                error = "Folder %s does exist." % folder
                raise AssertionError(error)

        pkints = [self.pkaint]
        g_matrices = [self.g]

        titratable_residues = self.descr.titratable_residues
        residue_list = self.descr.sorted_residue_list
        folder = kbp2.karlsberg.run_karlsberg_parallel(cpus, pkints, g_matrices, folder=folder,
                                     titratable_residues=titratable_residues, residue_list=residue_list)

        print("Karlsberg job finished in: %s" % folder)

        kbp2.karlsberg.parse_karlsberg_results(folder, kbp_result=self)
        self.find_pkas()

        if not keep_files:
            shutil.rmtree(folder)


class FrameKbpResults(object):
    """Object to store the results of Karlsberg+ calculations that are based on MD snapshots.
    """

    def __init__(self, descr=None):
        """ Constructor
        Parameters: descr : Object of type KbpJobDescr. The object is used directly, not copied! If not specified, the
                            KbpJobDescr object is taken from the first job added that has a KbpJobDescr object.
        """

        # Contains general properties of the Karlsberg+ job.
        if descr is not None:
            assert(isinstance(descr, KbpJobDescr))
        self.descr = descr

        # Contains 'kbp_results' objects, each representing the results of the Karlsberg+ calculation for a MD snapshot.
        # Format: {<frame number> : kbp_results <Karlsberg+ results>}
        self.kbp_results = {}

        # Float the translates frame number into nanoseconds. (optional)
        self.timestep = None

        # Specifies the frame range that is supposed to be analysed. Some functions use this variable, if set. (optional)
        # Format: [<first frame number>, <last frame number>]
        # -1 can be used for <last frame number> to select the last frame.
        self.frame_range = [0, -1]

        # Averaged values of the corresponding entries in the stored KbpResults. Set by 'average_over_frames'.
        # It is a KbpResult object.
        self.avg_results = None

        # Combined results of all frames. Set by 'combine_frames'.
        # It is a KbpResult object.
        self.combined_results = None

        # The ph dependent weight for each frame. The x-axis for the data are the ph values stored in 'descr.ph_values'.
        # Set by 'combine_frames'.
        # 2D numpy array: [ph, conformer]
        self.weights = None

    def __str__(self):
        # return self.__unicode__()
        msg = '<FrameKbpResults object containing %i frames>' % self.get_nr_of_frames()
        return msg

    def __iter__(self):
        frame_nrs = sorted(self.kbp_results.keys())
        for frame_nr in frame_nrs:
            yield frame_nr, self.kbp_results[frame_nr]

    def get_nr_of_frames(self):
        return len(self.kbp_results)

    def get_task_names(self):
        task_names = []
        for frame_nr, kbp_result in self:
            folder = kbp_result.folder
            task_name = folder.strip('/').split('/')[-1]
            task_names.append(task_name)
        return task_names

    def add_task(self, kbp_result, frame_nr=None):
        """ Adds a KbpResults object. If the property kbp_results.desc is None it is set to the self.descr. If
        self.descr is None it is set to kbp_results.desc
        """
        if frame_nr is None:
            frame_nr = len(self.kbp_results)

        assert(not self.kbp_results.has_key(frame_nr))
        assert(isinstance(kbp_result, KbpResult))

        if kbp_result.descr is None:
            kbp_result.descr = self.descr
        elif self.descr is None:
            self.descr = kbp_result.descr

        self.kbp_results[frame_nr] = kbp_result

    def set_task_conf_energies(self, name, frame_nr, energy):
        if type(energy) == np.ndarray:
            self.kbp_results[frame_nr].conf_energies[name] = energy.copy()
        else:
            self.kbp_results[frame_nr].conf_energies[name] = energy

    def get_task_conf_energies(self, name, frame_nr):
        return self.kbp_results[frame_nr].conf_energies[name]

    def set_misc_info_energies(self, name, frame_nr, value):
        self.kbp_results[frame_nr].misc_info[name] = value


    def is_in_range(self, frame_nr, frame_range=None):
        """Checks if a given frame_nr is within a interval defined by frame_range. it returns the result of:
        frame_range[0] < frame_nr < frame_range[1]
        If frame_range[1] has the value -1, no upper boundary is considered.
        @param frame_nr: int
        @param frame_range: [<first frame number>, <last first frame number or -1>] if not specified, it is taken from
                            self.frame_range
        @return: bool
        """
        if frame_range is None:
            frame_range = self.frame_range

        if frame_range is not None:
            x_min = frame_range[0]
            x_max = frame_range[1]
            if frame_nr < x_min:
                return False
            if x_max != -1:
                if frame_nr > x_max:
                    return False
        return True

    def average_over_frames(self, frame_range=None):
        """Averages the following properties of the stored frames_kbp_results object over the frames: pkaint, g and occs
        The frame selection of self.frame_range is used if given as parameter or set previously. The results are stored
        as a new KbpResult object in self.avg_results.

        @return: None
        """
        if frame_range is not None:
            self.frame_range = list(frame_range)

        nr_of_residues = len(self.descr.sorted_residue_list)
        nr_of_ph_values = len(self.descr.ph_values)
        max_states = 0
        for residue_kbp in self.descr.sorted_residue_list:
            resname_kbp = residue_kbp[0:residue_kbp.find('-')]
            # Skip the R state.
            nr_of_states = len(self.descr.titratable_residues[resname_kbp])
            if max_states < nr_of_states:
                max_states = nr_of_states

        first_conf_energies_dict = self.kbp_results.values()[0].conf_energies
        conf_energies_keys = first_conf_energies_dict.keys()
        conf_energies_avg_dict = {}
        for key in conf_energies_keys:
            first_entry = first_conf_energies_dict[key]
            if type(first_entry) == np.ndarray:
                conf_energies_avg_dict[key] = first_entry * 0.0
            else:
                conf_energies_avg_dict[key] = 0.0

        # g_avg = np.zeros((nr_of_residues, nr_of_residues, max_states, max_states), dtype=np.float)
        # Since the structure of self.g is only consisting of numpy arrays, this operation is allowed and creates a
        # copy of the structure with zeros.

        g_avg = self.kbp_results.values()[0].g * 0.0

        occs_avg = np.zeros(nr_of_residues, dtype=object)
        pkaint_avg = {}
        for (residue_num, residue_kbp, resname_kbp, nr_of_states) in self.descr.residue_list_ext:
            pkaint_avg[residue_kbp] = 0.0
            occs_avg[residue_num] = np.zeros((nr_of_states, len(self.descr.ph_values)), dtype=np.float)
        conf_energy_avg = 0.0
        residue_sasa_avg = np.zeros(nr_of_residues, dtype=np.float)

        nr_of_frames = 0
        for frame_nr, kbp_result in self.kbp_results.iteritems():

            if not self.is_in_range(frame_nr):
                continue

            nr_of_frames += 1

            g_avg += kbp_result.g
            for (residue_num, residue_kbp, resname_kbp, nr_of_states) in self.descr.residue_list_ext:
                pkaint_avg[residue_kbp] += kbp_result.pkaint[residue_kbp]
                # occs_avg[residue_num] += kbp_result.occs[residue_num]
                occs_avg[residue_num] = occs_avg[residue_num] + kbp_result.occs[residue_num]
            for key in conf_energies_keys:
                conf_energies_avg_dict[key] += kbp_result.conf_energies[key]
            if (kbp_result.conf_energy is not None) and (conf_energy_avg is not None):
                conf_energy_avg += kbp_result.conf_energy
            else:
                conf_energy_avg = None
            if (kbp_result.residue_sasa is not None) and (residue_sasa_avg is not None):
                residue_sasa_avg += kbp_result.residue_sasa
            else:
                residue_sasa_avg = None

        g_avg /= nr_of_frames
        for (residue_num, residue_kbp, resname_kbp, nr_of_states) in self.descr.residue_list_ext:
            pkaint_avg[residue_kbp] /= nr_of_frames
            occs_avg[residue_num] /= nr_of_frames
        for key in conf_energies_keys:
            conf_energies_avg_dict[key] /= nr_of_frames
        if conf_energy_avg is not None:
            conf_energy_avg /= nr_of_frames
        if residue_sasa_avg is not None:
            residue_sasa_avg /= nr_of_frames


        self.avg_results = KbpResult(self.descr)
        self.avg_results.g = g_avg
        self.avg_results.pkaint = pkaint_avg
        self.avg_results.occs = occs_avg
        self.avg_results.conf_energies = conf_energies_avg_dict
        self.avg_results.conf_energy = conf_energy_avg
        self.avg_results.residue_sasa = residue_sasa_avg



    def combine_frames(self):
        """ Combines the results from all stored frames by weighting each frame with conformational and protonation
        energy.
        @return: None
        """

        R = 8.3144621              # J/(K*mol)
        # kb = 1.3806488e-23       # J/K
        RT = R * 300.0 / 1000.0    # kJ/mol

        ############################
        ### Extract the energies ###
        ############################
        prot_energies = []
        conf_ernergies = []
        for frame_nr, kbp_result in self:
            prot_energies.append(kbp_result.prot_energy)
            conf_ernergies.append(kbp_result.conf_energy)
        # Set first frame as reference.
        # prot_energies = [e - prot_energies[0] for e in prot_energies]
        min_value = min(conf_ernergies)
        conf_ernergies = [e - min_value for e in conf_ernergies]

        ph_values = self.descr.ph_values
        nr_of_frames = self.get_nr_of_frames()

        weights = np.zeros((len(ph_values), nr_of_frames), np.float)
        for ph_number, ph in enumerate(ph_values):
            e_tot = np.zeros(nr_of_frames, np.float)

            ######################################
            ### Start with Karlsberg+ energies ###
            ######################################
            for i in range(nr_of_frames):
                e_tot[i] += prot_energies[i][ph_number]

            ###############################################
            ### Add conformational electrostatic energy ###
            ###############################################
            for i in range(nr_of_frames):
                e_tot[i] += conf_ernergies[i]

            e_offset = min(e_tot)
            e_tot = [x - e_offset for x in e_tot]

            #################################
            ### Energies -> Probabilities ###
            #################################
            z = 0
            for e in e_tot:
                z += np.exp(-e / RT)
            p_tot = []
            for e in e_tot:
                p_tot.append(np.exp(-e / RT) / z)
            weights[ph_number] = p_tot

        nr_of_residues = len(self.descr.sorted_residue_list)
        deprot_curves_weighted = np.zeros((nr_of_residues, len(ph_values)), dtype=np.float)
        occs_weighted = self.kbp_results.values()[0].occs * 0.0
        for frame_nr, kbp_result in self:
            weight = weights[:, frame_nr]
            deprot_curves_weighted += kbp_result.deprot_curves * weight
            for (residue_num, residue_kbp, resname_kbp, nr_of_states) in self.descr.residue_list_ext:
                occs_weighted[residue_num] += kbp_result.occs[residue_num] * np.vstack([weight] * nr_of_states)


        self.combined_results = KbpResult(descr=self.descr)
        self.combined_results.deprot_curves = deprot_curves_weighted
        self.combined_results.occs = occs_weighted
        self.combined_results.find_pkas()
        self.weights = weights


    def combine_frames_karlsberg(self, folder=None, keep_files=False, cpus=1, ph_range=[-10, 20, 0.5], quiet=True):

        if folder is not None:
            if folder[-1] != '/':
                folder += '/'

            if os.path.exists(folder):
                error = "Folder %s does exist." % folder
                raise AssertionError(error)


        # Get conformational energies
        conf_energies = []
        pkints = []
        g_matrices = []
        for frame_nr, kbp_result in self:
            conf_energies.append(kbp_result.conf_energy)
            pkints.append(kbp_result.pkaint)
            g_matrices.append(kbp_result.g)

        if len(conf_energies) == 1:
            conf_energies = [0.0]
        else:
            for conf_energy in conf_energies:
                assert(conf_energy is not None)

        titratable_residues = self.descr.titratable_residues
        residue_list = self.descr.sorted_residue_list

        if cpus > 1:
            folder = kbp2.karlsberg.run_karlsberg_parallel(cpus, pkints, g_matrices, conf_energies, folder=folder,
                                     titratable_residues=titratable_residues, residue_list=residue_list,
                                     ph_range=ph_range, seed=123456)
        else:
            folder = kbp2.karlsberg.run_karlsberg(pkints, g_matrices, conf_energies, folder=folder,
                                     titratable_residues=titratable_residues, residue_list=residue_list,
                                     ph_range=ph_range, seed=123456)



        if not quiet:
            print("Karlsberg job finished in: %s" % folder)
        # else:
        #     print("Skipping calculation, since folder %s exisits" % folder)
        self.combined_results = KbpResult(descr=self.descr)

        # The list ph values may change at this point. So generate a new pH list and compare it with the existing ph
        # list and overwrite it if necessary.
        ph_values = [ph_range[0]]
        while ph_values[-1] < ph_range[1]:
            new_value = ph_values[-1]+ph_range[2]
            new_value = float('%.1f' % new_value)
            ph_values.append(new_value)
        if not self.combined_results.descr.ph_values == ph_values:
            self.combined_results.descr.ph_values = ph_values

        kbp2.karlsberg.parse_karlsberg_results(folder, kbp_result=self.combined_results)
        self.combined_results.find_pkas()
        self.weights = self.combined_results.conformer_occs


        if not keep_files:
            shutil.rmtree(folder)