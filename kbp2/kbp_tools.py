# -*- coding: utf-8 -*-

import re, os
import numpy as np

def parse_titratable_residues(kbp_folder):
    if kbp_folder[-1] != '/':
        kbp_folder += '/'

    # Assume that the titatabe.yaml is the same for all MDs. So look into the first Karlsberg+ job in the first MD.
    # first_bbp_job_folder = kbp_folder
    # subfolders = sorted(os.listdir(first_bbp_job_folder))
    # for subfolder in subfolders:
    #     if subfolder.find("frame") == 0:
    #         first_bbp_job_folder += subfolder + '/'
    #         break
    # else:
    #     raise AssertionError("No Karlsberg+ job found.")

    # Find the name of the used titratable yml file.
    titratable_file = None
    include_path = None
    config_filename = kbp_folder + "titrate.yml"
    f = open(config_filename)
    for line in f:
        if line.find("titratable:") > -1:
            titratable_file = line.split()[1]
        if line.find("include:") > -1:
            include_path = line.split()[1]
            if include_path[-1] != '/':
                include_path += '/'
    if titratable_file[0] != '/':
        titratable_file = include_path + titratable_file
    f.close()

    titratable_residues = parse_titratable_yaml(titratable_file)

    return titratable_residues

def parse_input_structure_filename(kbp_folder):
    if kbp_folder[-1] != '/':
        kbp_folder += '/'

    # Find the name of the input structure.
    filename = ''
    config_filename = kbp_folder + "titrate.yml"
    f = open(config_filename)
    for line in f:
        if line.find("pdb:") > -1:
            filename = line.split()[1]
            break
    if filename[0] != '/':
        filename = kbp_folder + filename
    f.close()

    return filename

def parse_sfc_prefix(kbp_folder):
    if kbp_folder[-1] != '/':
        kbp_folder += '/'

    # Find the name of the input structure.
    filename = ''
    config_filename = kbp_folder + "titrate.yml"
    f = open(config_filename)
    # To look for:
    #sfc:
    # - name: c_pH7_frame0
    sfc_reached = False
    for line in f:
        if line.find("sfc:") > -1:
            sfc_reached = True
        elif sfc_reached:
            prefix = line.split()[-1]
            break
    f.close()

    return prefix


def parse_titratable_yaml(filename):
    # Parse the content of the titarable yml file.
    f = open(filename)
    last_residue = None
    last_residue_name = None
    last_state = None
    titratable_residues = {}
    for line in f:
        if line.strip(' ')[0] == '#':
            continue
        if line[0].strip('\n') not in [' ', '-', '']:
            # New residue found.
            if last_residue is not None:
                # Store the last residue.
                last_residue.append(last_state)
                last_state = None
                titratable_residues[last_residue_name] = last_residue
            last_residue = []
            last_residue_name = line.strip(' :\n')
        elif line.find('- sym:') > -1:
            # New state found.
            if last_state is not None:
                last_residue.append(last_state)
            last_state = {}
            last_state['name'] = line.split()[2]
            last_state['atoms'] = {}
        elif line.find('pka:') > -1:
            last_state['pka'] = float(line.split()[1])
        elif line.find('patch:') > -1:
            last_state['patch'] = line.split()[1].strip('"')
        else:
            reg = re.compile(r'^\s+(\w+):\s+([-.\d]+)')
            reg_m = reg.match(line)
            if reg_m is not None:
                name = reg_m.groups()[0]
                charge = float(reg_m.groups()[1])
                last_state['atoms'][name] = charge
        # Store the last residue.
    last_residue.append(last_state)
    titratable_residues[last_residue_name] = last_residue
    f.close()

    return titratable_residues

def guess_prot_state(resname_kbp, res, titratable_residues):
    # Compare charges in structure with charges in titratable yml file.
    states_abs_diff = []
    for state in titratable_residues[resname_kbp]:
        state_abs_diff = 0.0
        for atom_name in state['atoms'].keys():
            # The residue could be deprotonated in the PDB file.
            if atom_name not in res.keys():
                state_abs_diff += state['atoms'][atom_name]
            else:
                state_abs_diff += (state['atoms'][atom_name] - res[atom_name]['charge'])

        states_abs_diff.append(abs(state_abs_diff))

    # Make a decision.
    best_state = states_abs_diff.index( min(states_abs_diff) )
    best_state_value = states_abs_diff[best_state]
    if best_state_value > 0.001:
        raise AssertionError("No matching protonation state found for residue %s-%i_%s. "\
                             % (resname_kbp, res.resid, res.segname)\
                             + "Missmatch to next best state (%i) is: %f"\
                             % (best_state, best_state_value))
    else:
        return best_state

def determine_md_protonation_pattern(md_folder, residue_list, titratable_residues):
    # Guess the protonation state for each residue in the MD.

    from file_parser import Simple_struct_parser

    # What are the standart amino acids?
    protein_aa = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'SEC', 'GLY', 'PRO',
                      'ALA', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR', 'VAL', 'HSP', 'HSD', 'HSE']

    prot_pattern = {}

    pdb = Simple_struct_parser()
    pdb.read_pdb(md_folder + "protein_in_water.pdb")
    pdb.read_xplor_psf(md_folder + "protein_in_water.xplor.psf")
    for seg in pdb.struct.iter_segments():
        first_res = None
        last_res  = None
        for res in seg.iter_residues():
            resname = res.resname

            resname_kbp = get_kbp_resname(resname)

            if titratable_residues.has_key(resname_kbp):
                residue_kbp = "%s-%i_%s" %(resname_kbp, res.resid, res.segname)
                if residue_kbp in residue_list:
                    best_state = guess_prot_state(resname_kbp, res, titratable_residues)
                    prot_pattern[residue_kbp] = best_state
                    # print "%s -> %i" % (residue_kbp, best_state)

            if first_res is None:
                first_res = res
            last_res = res

        # Add N- and C_terminus.
        if first_res is not None and first_res.resname in protein_aa:
            nter = "NTE-%i_%s" % (first_res.resid, first_res.chainid)
            if nter in residue_list:
                best_state = guess_prot_state('NTE', first_res, titratable_residues)
                prot_pattern[nter] = best_state

        if last_res is not None and last_res.resname in protein_aa:
            cter = "CTE-%i_%s" % (last_res.resid, last_res.chainid)
            if cter in residue_list:
                best_state = guess_prot_state('CTE', last_res, titratable_residues)
                prot_pattern[cter] = best_state

    return prot_pattern

def get_real_resname(resname):
    real_name = {}
    # base
    real_name['ARG'] = 'ARG'
    real_name['HSP'] = 'HIS'
    real_name['HIS'] = 'HIS'
    real_name['LYS'] = 'LYS'
    # acid
    real_name['DPP'] = 'ASP'
    real_name['ASP'] = 'ASP'
    real_name['EPP'] = 'GLU'
    real_name['GLU'] = 'GLU'
    real_name['CYS'] = 'CYS'
    real_name['TYR'] = 'TYR'
    real_name['PRA'] = 'PRA'
    real_name['PRD'] = 'PRD'
    # termini
    real_name['CTE'] = 'CTE'
    real_name['NTE'] = 'NTE'
    # sugar
    real_name['SIA'] = 'SIA'

    # return real_name[resname]

    if real_name.has_key(resname):
        return real_name[resname]
    else:
        return resname

def get_kbp_resname(resname):
    kbp_name = {}
    # base
    kbp_name['HIS'] = 'HSP'
    kbp_name['HSP'] = 'HSP'
    kbp_name['HSD'] = 'HSP'
    kbp_name['HSE'] = 'HSP'
    # base
    kbp_name['ASP'] = 'DPP'
    kbp_name['GLU'] = 'EPP'

    if kbp_name.has_key(resname):
        return kbp_name[resname]
    else:
        return resname

def get_real_residue(residue):
    (resname, resid_chain) = residue.split("-")
    (resid, chain) = resid_chain.split("_")

    real_resname = get_real_resname(resname)
    real_residue = "%s-%s_%s" % (real_resname, resid, chain)

    return real_residue

def get_kbp_residue(residue):
    if type(residue) is str:
        (resname, resid, segname) = re.split(r'[-_]', residue)
    elif type(residue) is tuple:
        (resname, resid, segname) = residue

    kbp_resname = get_kbp_resname(resname)
    kbp_residue = "%s-%s_%s" % (kbp_resname, resid, segname)

    return kbp_residue

def parse_g_pkint(pkint_filename, g_filename, residue_list_ext=None):
    """
    Parses a .g and .pkint file created by Karlsberg+.

    Returns 'pkint, g_matrix_np, residue_list_ext' if residue_list_ext is None and
            'return pkint, g_matrix_np' otherwise
    """
    pkint = {}

    # The order of the residues according to the .pkint file (Is needed to interpret the numbers in the .g file.)
    pkaint_residue_list = []

    # Parse pkint file.
    f = open(pkint_filename)
    for line in f:
        # 0.000000e+00 R 2.506675e+01 D 3.403473e+00 0 EPP-35_A
        entries = line.split()
        if len(entries) < 5:
            continue
        residue = entries[-1]

        new_residue = []
        for state_nr in range(0, int(len(entries)-1), 2):
            state_pkint = float(entries[state_nr  ])
            state_name =        entries[state_nr+1]
            new_residue.append(state_pkint)

        pkaint_residue_list.append(residue)
        pkint[residue] = np.array(new_residue, dtype=np.float32)
    f.close()


    # Parse g file.
    g_matrix = {}
    f = open(g_filename)
    #    1 2   21 2    1.449279e-05
    for line in f:
        entries = line.split()
        if len(entries) != 5:
            continue
        # State 0 - x and x - 0 interactions are zero by definition.
        # Residue and state numbering start at 1 in the .g file.
        residue1_nr = int(entries[0]) - 1
        residue1 = pkaint_residue_list[residue1_nr]
        state_1  = int(entries[1]) - 1

        residue2_nr = int(entries[2]) - 1
        residue2 = pkaint_residue_list[residue2_nr]
        state_2  = int(entries[3]) - 1

        interaction = float(entries[4]) * 1388.4269

        # 1 e^2/A  = 331.842 kcal/mol
        # 1 e^2/A  = 331.842 * 4.184 kJ/mol = 1388.4269 kJ/mol
        #  From the MEAD documentation:
        #  331.842  (kcal/mole)/(protons squared/Angstrom)
        #            Conversion factor for going from energy units
        #            of (charge units squared)/(length unit) to
        #            to whatever energy units you want for your output.

        # if residue2_nr > residue1_nr:
        if not residue1 in g_matrix:
            g_matrix[residue1] = {}
        if not residue2 in g_matrix[residue1]:
            g_matrix[residue1, residue2] = {}
        if not state_1 in g_matrix[residue1, residue2]:
            g_matrix[residue1, residue2, state_1] = {}
        g_matrix[residue1, residue2, state_1, state_2] = interaction

    f.close()

    # If not specified, create one
    if residue_list_ext is None:
        residue_list_ext = []
        for residue_num, residue_kbp in enumerate(pkaint_residue_list):
            resname_kbp = re.split(r'[-_]', residue_kbp)[0]
            nr_of_states = len(pkint[residue_kbp])
            residue_list_ext.append((residue_num, residue_kbp, resname_kbp, nr_of_states))
        residue_list_ext_created = True
    else:
        residue_list_ext_created = False

    # Adjust shape and order of the g matrix to residue_list_ext
    g_matrix_np = []
    for (residue_num, residue_kbp, resname_kbp, nr_of_states) in residue_list_ext:
        array1 = []
        for state1 in range(1, nr_of_states):
            array2 = []
            for (residue2_num, residue2_kbp, resname2_kbp, nr_of_states2) in residue_list_ext:
                if residue2_num <= residue_num:
                    array2.append(np.array([]))
                else:
                    array3 = np.zeros(nr_of_states2 - 1, dtype=np.float32)
                    for state2 in range(1, nr_of_states2):
                        first_entry = g_matrix[residue_kbp, residue2_kbp, state1, state2]
                        second_entry = g_matrix[residue2_kbp, residue_kbp, state2, state1]
                        # assert(abs(first_entry - second_entry) < 0.001 * abs(first_entry))
                        assert(abs(first_entry - second_entry) < 0.001)
                        array3[state2 - 1] = (first_entry + second_entry) / 2.0
                        # array3[state2 - 1] = g_matrix[residue_kbp, residue2_kbp, state1, state2]
                    array2.append(np.array(array3))
            array1.append(np.array(array2))
        g_matrix_np.append(np.array(array1))

    # g_matrix_np = np.array(g_matrix_np)

    g_matrix_np1 = np.array(g_matrix_np)
    g_matrix_np = g_matrix_np1

    # g_array_np:
    # State entries are shifted by -1.

    if residue_list_ext_created:
        return pkint, g_matrix_np, residue_list_ext
    else:
        return pkint, g_matrix_np

# def translate_prot_state(prot_state, titratable_yaml, titr_residues):
#     """ Converts a protonation pattern referring to states defined by kbp_tools.determine_md_protonation_pattern into a
#      titr_residue_dict object as created by charmm.Charm_manager.get_titr_residue_dict.
#     @param prot_state: States defined by kbp_tools.determine_md_protonation_pattern.
#     @param titratable_yaml: state definition created by charmm.parse_titratable_yaml
#     @param titr_residues: state definition created by charmm.get_std_titr_residues
#     @return: titr_residue_dict: object as created by charmm.Charm_manager.get_titr_residue_dict
#     """
#
#     for residue_descr, state_yaml in prot_state.iteritems():
#         resname, resid, segname = re.split(r'[-_]', residue_descr)
#         residue_tuple = (resid, segname)
#         # print titratable_yaml
#         if  titratable_yaml[resname][state_yaml].has_key('patch'):
#             patch = titratable_yaml[resname][state_yaml]['patch']
#         else:
#             patch = None

def parse_kb_energy(kbp_folder):
    files = os.listdir(kbp_folder)
    kb_energies = {}
    reg = re.compile(r'^kb_frame\d+.o\d+$')
    for filename in files:
        reg_m = reg.match(filename)
        if reg_m is not None:
            f = open(kbp_folder + '/' + filename)
            last_E = None
            reg_line = re.compile(r'^.+the average energy is\s+([-.\deE]+)\s+kJ/mol\.$')
            reg_line_ph = re.compile(r'^.+kb_frame\d+_pH([-.\d]+)_0mV_300K$')
            for line in f:
                reg_line_m = reg_line.match(line)
                if reg_line_m is not None:
                    last_E = float(reg_line_m.groups()[0])

                reg_line_ph_m = reg_line_ph.match(line)
                if reg_line_ph_m is not None:
                    ph = float(reg_line_ph_m.groups()[0])
                    if last_E is None:
                        print('last_E is not set. There is something wrong with the output file!')
                        break
                    if kb_energies.has_key(ph):
                        print "### Warning: Energy in folder %s is defined more than once!" \
                              % kbp_folder
                    kb_energies[ph] = last_E
            f.close()
    return kb_energies

def parse_gernots_exp_pkas(filename, verbose=False):
    f = open(filename, 'r')

    pkas = {}
    if verbose:
        print("Parsing experimental pKa values from: " + filename)
        print("The following lines have been skiped:")
        print('---')
    for line in f:
        line = line.rstrip()
        # ASP-18      1.92   	 1.93 	 1.92 	2.93    3.54	  2.66(0.1)
        reg = re.compile(
            r'^\s*(\w+)-(\d+)_?(\w?)\s+[-.\d]+\s+[-.\d]+\s+[-.\d]+\s+[-.\d]+\s+[-.\d]*\s+([-.\d]+)\(?[.\d]*\)?(\s+#\s+.+)?$')
        reg_m = reg.match(line)
        if reg_m is not None:
            resname = reg_m.groups()[0]
            resid = int(reg_m.groups()[1])
            chain = reg_m.groups()[2]
            pka = float(reg_m.groups()[3])
            if resname == "HSP" or resname == "HSD":
                resname = "HIS"
            if chain == '':
                chain = "A"
            residue = "%s-%i_%s" % (resname, resid, chain)

            if pkas.has_key(residue):
                if not pkas[residue] == pka:
                    raise AssertionError("Different pka values found for the same residue in: %s" % filename)
            else:
                pkas[residue] = pka
        else:
            if verbose:
                if line != '':
                    print(line)
    if verbose:
        print('---')
    f.close()
    return pkas


