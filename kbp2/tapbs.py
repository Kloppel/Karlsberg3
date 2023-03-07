# -*- coding: utf-8 -*-

import re
import subprocess
import sys
import os
import numpy as np
import file_parser
import kbp_tools
import inspect
import charmm


def get_grid_information(structure, residue_list):

    '''This function returns grid sizes for Tapbs calculation AND appropriate resolution in case #3 calculations are needed

        fine_dime ->  dimensions of a grid with the best resolution (calculation  #1 protein & model)
        middle_dime -> dimensions of a grid with the middle resolution (calculation  #2 protein & model)
        coarse_resolution -> calculated resolution for the grid with middle_dime (calculation  #3 protein)
    '''

    # Find the largest dimension that fits all titratable residues.
    # We add 6 Angstroms to the size to get also parts of the environment
    # at this resolution. In each direction the titratable has a sharp
    # view at least within a radius of 3 Angstroms.
    max_residue_size = get_max_residue_size(structure, residue_list)
    max_residue_size = max_residue_size + 6
    protein_size = get_protein_size(structure)

    fine_resolution = 0.25
    middle_resolution = 1.0

    #### 1. FINE DIMENSION

    fine_dime = size2dime(max_residue_size, fine_resolution)

    #### 2. MIDDLE DIMENSION

    fine_size = dime2size(fine_dime, fine_resolution)
    #Todo: ADD COMMENT
    middle_dime = size2dime(protein_size + fine_size, middle_resolution)

    #### 3. COARSE DIMENSION

    middle_size = dime2size(middle_dime, middle_resolution)

    coarse_resolution = (4 * protein_size) / middle_size
    coarse_resolution_max = max(coarse_resolution)
    coarse_resolution = (coarse_resolution_max, coarse_resolution_max, coarse_resolution_max)

    return fine_dime, middle_dime, coarse_resolution


def dime2size(dime, target_res):

    new_size = dime - 1
    new_size = new_size * target_res

    return new_size


def size2dime(size, target_res):

    # dime = np.array([0,0,0])
    # for i in range(3):
    #     for d in [65, 97, 129, 161, 193, 225, 257, 289, 321, 353, 385, 417, 449]:
    #         res = size[i] / float(d)
    #         if res <= target_res:
    #             dime[i] = d
    #             break
    #     else:
    #         dime[i] = d
    #         res = size[i] / float(d)
    #         print("WARNING (tapbs): Resolution in dimension %i is larger than desired: %f" % (i, res))

    dime = np.array([0,0,0])
    for i in range(3):
        for d in [65, 97, 129, 161, 193, 225, 257, 289, 321, 353, 385, 417, 449]:
            dim = (size[i]/target_res)+1
            if dim <= d:
                dime[i] = d
                break
        else:
            dime[i] = d
            res = size[i] / float(d)
            raise AssertionError("WARNING (tapbs): Resolution in dimension %i is larger than desired: %f" % (i, res))

    return dime

def get_max_residue_size(structure, residue_list):

    max_size = [-np.inf, -np.inf, -np.inf]

    for resname, resid, segname in residue_list:
        residue = structure.struct[segname][resid]
        for atom in residue.iter_atoms():
            if not atom['vdw']:
                atom['vdw'] = structure.vdw_table[atom['type']]
        (min_coord, max_coord) = residue.get_size()
        min_coord = np.array(min_coord)
        max_coord = np.array(max_coord)
        size = max_coord - min_coord
        max_size = np.max([size, max_size], axis=0)

    return max_size

def get_protein_size(structure):

    (min_coord_p, max_coord_p) = structure.get_size()
    min_coord_p = np.array(min_coord_p)
    max_coord_p = np.array(max_coord_p)

    protein_size = max_coord_p - min_coord_p

    return protein_size


def termini_detection(residue_tuple, residue_list):

    resname_main, resid_main, segname_main = residue_tuple

    is_termini = False

    if resname_main in ['NTE', 'CTE']:
        return is_termini

    for resname, resid, segname in residue_list:
        if resid == resid_main and segname == segname_main and (resname in ['NTE', 'CTE']):
            is_termini = True
            break
        else:
            continue

    return is_termini


def run_tapbs(folder, jobname, residue_list, structure, titratable_residues, tapbs_bin, pdie=4, sdie=80, cavity_par='', quiet=True, md_evaluation_mode = False):
    """
    residue_list can be a string 'ASP-34_A' or tuple ('ASP', 34, 'A')
    """

    if folder[-1] != '/':
        folder += '/'

    if os.path.exists(folder):
        error = "Folder %s does exist." % folder
        raise AssertionError(error)

    os.mkdir(folder)

    if type(residue_list[0]) is str:
        residue_list_descr = residue_list
        residue_list = []
        for residue_descr in residue_list_descr:
            resname, resid, segname = re.split(r'[-_]', residue_descr)
            resid = int(resid)
            residue_list.append((resname, resid, segname))


    structure = structure.copy()
    titratable_residues_mod = {}


    for resname in titratable_residues:

        # check if the residue should really be titrated
        res_to_titrate = False
        for res_tuple in residue_list:
            if resname == res_tuple[0]:
                res_to_titrate = True
                break
        if not res_to_titrate:
            continue

        if not titratable_residues[resname][0].has_key('external_patches'):
            continue
        res_def = titratable_residues[resname]
        titratable_residues_mod_entry = []

        for nr_of_states in range(len(res_def)):
            titratable_residues_mod_entry.append({})

        for i, state in enumerate(res_def):
            for res_tuple in residue_list:
                if res_tuple[0] == resname:
                    main_residue_tuple = res_tuple
                    break
                else:
                    continue
                    # Todo: Residue is not titrated. asserttion error

            titratable_residues_mod_entry[i]['pka'] = titratable_residues[resname][i]['pka']
            # rename is excluded because it is not allowed to have patches and/or external_patches and rename
            titratable_residues_mod_entry[i]['name'] = titratable_residues[resname][i]['name']
            main_atom_list = titratable_residues[resname][i]['atoms']
            titratable_residues_mod_entry[i]['atoms'] = {}
            main_patch_set = False
            for name, charge in main_atom_list.iteritems():
                titratable_residues_mod_entry[i]['atoms'][name] = charge
                main_patch_set = True

            resname_main, resid_main, segname_main = main_residue_tuple
            residue_main = structure.struct[segname_main][resid_main]

            # everything needs to merge one way or another if there is an external patch
            for k, external_patch in enumerate(titratable_residues[resname][i]['external_patches']):
                external_patch_name, external_patch_residues = external_patch
                for m, residue_tuple in enumerate(external_patch_residues):
                    is_main_residue = False
                    if residue_tuple == main_residue_tuple:
                        is_main_residue = True
                        if main_patch_set:
                            continue

                    resname_src, resid_src, segname_src = residue_tuple
                    patch_atoms = titratable_residues[resname][i]['external_atoms'][k][m]

                    for atom_name, charge in patch_atoms.iteritems():
                        mod_atom_name = atom_name
                        # Find a new atom name if necessary
                        if mod_atom_name in residue_main.keys() and not is_main_residue:
                            if len(mod_atom_name) >= 4:
                                mod_atom_name = 'A'
                            counter = 1
                            while True:
                                new_atom_name = mod_atom_name + str(counter)
                                if len(new_atom_name) >= 4:
                                    mod_atom_name = 'A'
                                    counter = 1
                                    new_atom_name = mod_atom_name
                                if new_atom_name not in residue_main.keys():
                                    break
                                else:
                                    counter += 1
                        else:
                            new_atom_name = atom_name

                        # Atom is completely new - just add it
                        titratable_residues_mod_entry[i]['atoms'].update({new_atom_name: charge})
                        if not is_main_residue:
                            # Copy atom in structure
                            atom =  structure.struct[segname_src][resid_src][atom_name]
                            atom['segname'] = segname_main
                            atom['resid'] = resid_main
                            atom['resname'] = resname_main


        titratable_residues_mod[resname] = titratable_residues_mod_entry


    structure.create_struct()

    # Write .st files
    # Write .sites file
    # Write .pqr file
    # Write input script

    sites_file_str = ''
    resnames = []
    is_termini = False

    if md_evaluation_mode:
        # @todo: hotfixed solution when not true termini is titratable
        first_residues = []
        last_residues = []
        for segment in structure.struct.iter_segments():
            first_residue = None
            last_residue = None
            for residue in segment.iter_residues():
                if charmm.is_std_aa(residue.resname):
                    if first_residue is None:
                        first_residue = residue
                    else:
                        # break
                        last_residue = residue
            if first_residue is not None:
                residue_tuple_first = (first_residue.resname, first_residue.resid, first_residue.segname)
                first_residues.append(residue_tuple_first)
            if last_residue is not None:
                residue_tuple_last = (last_residue.resname, last_residue.resid, last_residue.segname)
                last_residues.append(residue_tuple_last)

    for resname, resid, segname in residue_list:

        st_filename = resname + '.st'

        if md_evaluation_mode:
            residue_tuple = (resname, resid, segname)
            is_termini = termini_detection(residue_tuple, residue_list)

            #@todo: part of the hotfix above
            if not is_termini:
                if residue_tuple in first_residues:
                    is_termini = True
                elif residue_tuple in last_residues:
                    is_termini = True

            if is_termini:
                st_filename = resname[0:2] + 'T.st'

        if (resname not in resnames) or is_termini:

            resnames.append(resname)
            # Create .st file for this residue
            st_file_str = ''
            if resname in titratable_residues_mod.keys():
                res_def = titratable_residues_mod[resname]
            else:
                res_def = titratable_residues[resname]

            for state_def in res_def:
                st_file_str += "%.2f pK %s\n" % (state_def['pka'], state_def['name'])

                for atom_nr, atom_name in enumerate(state_def['atoms']):

                    if is_termini:
                        if atom_name not in ['N', 'HN', 'CA', 'C', 'O']:
                            charge = state_def['atoms'][atom_name]
                            st_file_str += "ATOM   %4i %4s %3s A   1    9999.9999999.9999999.999%6.3f99.999      A\n" \
                                % (atom_nr, atom_name, resname, charge)
                    else:
                        charge = state_def['atoms'][atom_name]
                        st_file_str += "ATOM   %4i %4s %3s A   1    9999.9999999.9999999.999%6.3f99.999      A\n" \
                                % (atom_nr, atom_name, resname, charge)

            f = open(folder + st_filename, 'w')
            f.write(st_file_str)
            f.close()
        sites_file_str += "%s %i %s %s\n" % (segname, resid, resname, st_filename)
    sites_filename = jobname + '.sites'
    f = open(folder + sites_filename, 'w')
    f.write(sites_file_str)
    f.close()

    # Write pqr file
    pqr_filename = jobname + '.pqr'
    structure.write_pqr(folder + pqr_filename, kb_style=True)
    # structure.write_pqr(folder + pqr_filename)

    ####################################
    #### information about the grid ####
    ####################################

    (fine_dime, middle_dime, coarse_resolution) = get_grid_information(structure, residue_list)

    if 'cavity' in cavity_par:
        f = open(folder + 'cavities.conf','w')
        f.write(cavity_par)
        f.close()

    input_script = ''
    input_script += 'pqr ' + pqr_filename + '\n'
    input_script += 'sites ' + sites_filename + '\n'
    input_script += 'output ' + jobname + '\n'
    input_script += '\n\n'


    #########
    input_script += """

%s                                 """ % cavity_par + """
#writemap

bcfl sdh
pdie %s                                 """ % str(pdie) + """
sdie %s                                 """ % str(sdie) + """
temperature 300
ion -1.0 0.1 2.0
ion 1.0 0.1 2.0

srfm mol
srad 1.4
sdens 3
chgm spl2

dimension_of_protein_grid %i %i %i     """ % tuple(middle_dime) + """
spacing_of_protein_grid %.3f %.3f %.3f """ % tuple(coarse_resolution) + """
center_of_protein_grid oncent

dimension_of_protein_grid %i %i %i     """ % tuple(middle_dime) + """
spacing_of_protein_grid 1 1 1
center_of_protein_grid oncent

dimension_of_protein_grid %i %i %i     """ % tuple(fine_dime) + """
spacing_of_protein_grid 0.25 0.25 0.25
center_of_protein_grid onsite

dimension_of_model_grid %i %i %i       """ % tuple(fine_dime) + """
spacing_of_model_grid 1 1 1
center_of_model_grid onsite

dimension_of_model_grid %i %i %i       """ % tuple(fine_dime) + """
spacing_of_model_grid 0.25 0.25 0.25
center_of_model_grid onsite


errtol 1E-6
# this works most of the time (i haven't encountered any problem so far, but note that the APBS default is 255)
itmax 4
presmooth 2
postsmooth 2
iinfo 0
    """

    input_filename = jobname + '.tapbs.in'
    output_filename  = jobname + '.tapbs.out'
    f = open(folder + input_filename, 'w')
    f.write(input_script)
    f.close()


    shell = subprocess.Popen('bash\n',\
            stdin=subprocess.PIPE,\
            stdout=subprocess.PIPE,\
            stderr=subprocess.PIPE,\
            shell=True\
            )


    shell.stdin.write('cd ' + folder + '\n')
    shell.stdin.write("%s < %s > %s\n" % (tapbs_bin, input_filename, output_filename))
    shell.stdin.write('exit\n')


    if not quiet:
        print 'Waiting for TAPBS to finish..'
    while True:
        nextline = shell.stdout.readline()
        if not quiet:
            print nextline
        if not nextline:
            break

# #check again if it is ok
# def run_local(run_folder, quiet=True):
#     if run_folder[-1] != '/':
#         run_folder += '/'
#
#     shell = subprocess.Popen('csh\n',\
#             stdin=subprocess.PIPE,\
#             stdout=subprocess.PIPE,\
#             stderr=subprocess.PIPE,\
#             shell=True\
#             )
#
#     shell.stdin.write('%srun.sh\n' % run_folder )
#     shell.stdin.write('exit\n')
#
#     if not quiet:
#         print 'Waiting for APBS to finish..'
#     while True:
#         nextline = shell.stdout.readline()
#         if not nextline:
#             break


if __name__ == '__main__':
    tapbs_bin = '/scratch/scratch/tmeyer/kbplus2/tapbs_1.3_cav_enere'

    jobname = 'first_run'
    # folder = '/scratch/scratch/tmeyer/tmp/tapbs_test/'
    folder = '/user/jdragelj/python/karlsbergplus/tests/tapbs/'
    psf = '/scratch/scratch/tmeyer/md_pka/runs/general2/2lzt/ph7/modelling/2lzt_m.xplor.psf'
    crd = '/scratch/scratch/tmeyer/md_pka/runs/general2/2lzt/ph7/modelling/2lzt_m.crd'
    par = '/scratch/scratch/tmeyer/karlsbergplus/par.inp'
    titratable_yaml = '/scratch/scratch/tmeyer/kbplus2/titratable.yaml'


    structure = file_parser.Simple_struct_parser()
    structure.read_crd(crd)
    structure.read_xplor_psf(psf)
    structure.read_par(par)
    structure.struct['A'][15].rename('HSP')
    structure = structure.copy(segname='A')

    titratable_residues = kbp_tools.parse_titratable_yaml(titratable_yaml)

    residues_to_titrate = ['LYS-1_A', 'ARG-5_A', 'NTE-1_A', 'ASP-18_A']
    # , 'HSP-15_A'


    # Run
    run_tapbs(folder, jobname, residues_to_titrate, structure, titratable_residues, tapbs_bin)

    # Read the results
    pkint_filename = folder + jobname + '.pkint'
    g_filename = folder + jobname + '.g'

    pkint, g, residue_list_ext = kbp_tools.parse_g_pkint(pkint_filename, g_filename)

