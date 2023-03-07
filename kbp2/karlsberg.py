from collections import defaultdict
from multiprocessing import Process

import os
import re
import subprocess
import tempfile
import kbp2
import numpy as np

def write_pkint_g_files(pkint, g, basename, titratable_residues, residue_list=None):
    if residue_list is None:
        residue_list = pkint.keys()
        sorted_residue_list = kbp2.kbp_results.KbpJobDescr.sort_residue_list(residue_list)
    else:
        sorted_residue_list = residue_list


    f_pkint = open(basename + '.pkint', 'w')
    # 0.000000e+00 R -5.380468e+01 P LYS-1_A
    for residue in sorted_residue_list:
        resname = residue.split('-')[0]
        line = ''
        for state, energy in enumerate(pkint[residue]):
            # line += '%.6f ' % energy
            line += '{:6e} '.format(float(energy))
            line += titratable_residues[resname][state]['name'] + ' '
        line += residue + '\n'
        f_pkint.write(line)
    f_pkint.close()

    f_g = open(basename + '.g', 'w')
    #   1 2    1 2    0.000000e+00
    for res1, residue1 in enumerate(sorted_residue_list):
        resname1 = residue1.split('-')[0]
        for state1 in range(len(titratable_residues[resname1])):
            for res2, residue2 in enumerate(sorted_residue_list):
                resname2 = residue2.split('-')[0]
                for state2 in range(len(titratable_residues[resname2])):
                    if 0 in [state1, state2]:
                        continue
                    if (res1 == res2):
                        energy = 0.0
                    elif res1 < res2:
                        energy = g[res1][state1-1][res2][state2-1]
                    else:
                        energy = g[res2][state2-1][res1][state1-1]

                    # Convert energy from kJ/mol in e^2/A
                    energy /= 1388.4269

                    line = '{:4}{:2} {:4}{:2}    {:6e}\n'.format(res1+1, state1+1, res2+1, state2+1, float(energy))
                    f_g.write(line)

    f_g.close()


def run_karlsberg_parallel(cpus, pkint, g, conf_energies=[0.0], folder=None, overwrite_folder=False,  karlsberg_bin=None, \
                         ph_range=[-10, 20, 0.5], seed=123456, titratable_residues=None, residue_list=None):



    # Split ph_range into chunks
    ph_min, ph_max, ph_step = ph_range
    ph_list = []
    ph = ph_min
    while ph <= ph_max:
        ph_list.append(float('%.2f' % ph))
        ph += ph_step

    chunk_size = int(np.ceil(len(ph_list) / float(cpus)))

    ph_chunk_ranges = []
    for i in xrange(0, len(ph_list), chunk_size):
        ph_section = ph_list[i:i+chunk_size]
        ph_chunk_range = [ph_section[0], ph_section[-1], ph_step]
        ph_chunk_ranges.append(ph_chunk_range)


    if folder is not None:
        if os.path.exists(folder):
            if not overwrite_folder:
                error = "Folder %s does exist." % folder
                raise AssertionError(error)
        else:
            os.mkdir(folder)
    else:
        folder = tempfile.mkdtemp()

    if folder[-1] != '/':
        folder += '/'

    # Folder checks are done already, no need to repeat them in 'run_karlsberg'
    overwrite_folder = True


    if type(pkint) == str:
        pkint = [pkint]
    if type(g) == str:
        g = [g]

    if type(pkint[0]) != str:
        pkint_filenames = []
        g_filenames = []
        # Write pkint and g files
        for i, pkint_entry, g_entry in zip(range(len(pkint)), pkint, g):
            if titratable_residues is None:
                error = "To write pkint and g files, the parameter titratable_residues must be specified."
                raise(AssertionError(error))
            basename = folder + 'structure_' + str(i)
            write_pkint_g_files(pkint_entry, g_entry, basename, titratable_residues, residue_list=residue_list)
            pkint_filenames.append(basename + '.pkint')
            g_filenames.append(basename + '.g')
    else:
        pkint_filenames = pkint
        g_filenames = g

    # Stores the running processes.
    processes = []

    for core, ph_chunk_range in enumerate(ph_chunk_ranges):
        args = (pkint_filenames, g_filenames, conf_energies, folder, overwrite_folder, karlsberg_bin, ph_chunk_range,
                seed, titratable_residues, residue_list, '_c' + str(core))
        p = Process(target=run_karlsberg, args=args)
        p.start()
        processes.append(p)
        if cpus == 1:
            p.join()

    if cpus > 1:
        for p in processes:
            p.join()

    # Sequential version:
    # for ph_chunk_range in ph_chunk_ranges:
    #     run_karlsberg(pkint, g, conf_energies=conf_energies, folder=folder, overwrite_folder=overwrite_folder,
    #                   karlsberg_bin=karlsberg_bin, ph_range=ph_chunk_range, seed=seed,
    #                   titratable_residues=titratable_residues, residue_list=residue_list)

    return folder

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def run_karlsberg(pkint, g, conf_energies=[0.0], folder=None, overwrite_folder=False,  karlsberg_bin=None,
                  ph_range=[-10, 20, 0.5], seed=123456, titratable_residues=None, residue_list=None, script_suffix=''):
    """
    pkint: string or list of strings to .pkint files created by tapbs
    g: string or list of strings to .g files created by tapbs
    If the two parameters above are not strings or list of strings it is assumed, that the are list of objects that
    contain information about intrinsic pkas (pkaint) and interactions (g) in the format specified in
    kbp_results.KbpResult.g and kbp_results.KbpResult.pkint. In that case titratable_residues must be set.

    titratable_residues: dict containing information about the titratable residues in a format as specified for
    kbp_results.KbpJobDescr.titratable_residues.

    returns: the path to the folder, that contains the results.
    """


    if karlsberg_bin is None:
        karlsberg_bin = "../kb2plus_package/binaries_karlsberg/bin/karlsberg2.x86_64_fixed"

    if folder is not None:
        if os.path.exists(folder):
            if not overwrite_folder:
                error = "Folder %s does exist." % folder
                raise AssertionError(error)
        else:
            os.mkdir(folder)
    else:
        folder = tempfile.mkdtemp()

    if folder[-1] != '/':
        folder += '/'

    if type(pkint) == str:
        pkint = [pkint]
    if type(g) == str:
        g = [g]

    if type(pkint[0]) != str:
        pkint_filenames = []
        g_filenames = []
        # Write pkint and g files
        for i, pkint_entry, g_entry in zip(range(len(pkint)), pkint, g):
            if titratable_residues is None:
                error = "To write pkint and g files, the parameter titratable_residues must be specified."
                raise(AssertionError(error))
            basename = folder + 'structure_' + str(i)
            write_pkint_g_files(pkint_entry, g_entry, basename, titratable_residues, residue_list=residue_list)
            pkint_filenames.append(basename + '.pkint')
            g_filenames.append(basename + '.g')
    else:
        pkint_filenames = pkint
        g_filenames = g

    conf_energies = [x - conf_energies[0] for x in conf_energies]

    if not len(pkint_filenames) == len(g_filenames) == len(conf_energies):
        error = "Number of specified pkint, g files and conformational energies are not the same."
        raise(AssertionError(error))

    ph_start = ph_range[0]
    ph_end = ph_range[1]
    ph_incr = ph_range[2]


    conformations = ''
    for conf in range(len(pkint_filenames)):
        conformations += "conformation %s %s %.2f kJ/mol\n" \
                         % (pkint_filenames[conf], g_filenames[conf], conf_energies[conf])


    config_script = conformations + """

output kb

full_scans 2500
reduced_scans 10000
reduced_set_tolerance 0.1
conformation_moves_per_scan 1
tempering_moves_per_scan 1

correlation_limit 0.1
max_correlation_time 100

pH_start   %.1f """ % ph_start + """
pH_end     %.1f """ % ph_end + """
pH_incr    %.1f """ % ph_incr + """
redox_start 0
redox_end   0
redox_incr 25


temperature 300

min_int_pairs 2.5 pK
min_int_triples 5.0 pK

seed %i """ % seed



#     config_script = conformations + """
#
# output kb
#
# full_scans 2500
# reduced_scans 10000
# reduced_set_tolerance 0.1
# conformation_moves_per_scan 10
# tempering_moves_per_scan 1
#
# correlation_limit 0.1
# max_correlation_time 100
#
# pH_start   %.1f """ % ph_start + """
# pH_end     %.1f """ % ph_end + """
# pH_incr    %.1f """ % ph_incr + """
# redox_start 0
# redox_end   0
# redox_incr 25
#
#
# min_int_pairs 2.5 pK
# min_int_triples 5.0 pK
#
# temperature 300
# temperature 373
# temperature 465
# temperature 579
# temperature 722
# temperature 900
#
#
# seed %i """ % seed

#     config_script = conformations + """
#
# output kb
#
# full_scans 5000
# reduced_scans 25000
# reduced_set_tolerance 0.001
# conformation_moves_per_scan 10
# tempering_moves_per_scan 1
#
# correlation_limit 0.1
# max_correlation_time 100
#
# pH_start   %.1f """ % ph_start + """
# pH_end     %.1f """ % ph_end + """
# pH_incr    %.1f """ % ph_incr + """
# redox_start 0
# redox_end   0
# redox_incr 25
#
#
# min_int_pairs 2.5 pK
# min_int_triples 5.0 pK
#
# temperature 300
# temperature 373
# temperature 465
# temperature 579
# temperature 722
# temperature 900
#
#
# seed %i """ % seed


    in_out_basename = 'karlsberg' + script_suffix

    f = open(folder + in_out_basename + '.in', 'w')
    f.write(config_script)
    f.close()

    commands = "cd %s\n" % folder
    commands += karlsberg_bin + " < %s.in > %s.out\n" % (in_out_basename, in_out_basename)
    commands += "exit\n"
    process = subprocess.Popen("tcsh", \
                               shell=True, \
                               stdin=subprocess.PIPE, \
                               stdout=subprocess.PIPE, \
                               stderr=subprocess.PIPE, \
        )
    process.stdin.write(commands)

    kb_output = ''
    while True:
        next_line = process.stdout.readline()

        if not next_line:
            break

    return folder

def parse_karlsberg_results(folder, kbp_result=None, output_prefix='kb'):
    """

    @param folder: Folder to work
    @param kbp_result: KbpResult object to store the results, if not specified a new object is created.
    @return: A KbpResult object.
    The following fields of the KbpResult object are set:
    - kbp_result.descr.sorted_residue_list (if required)
    - kbp_result.descr.residue_list_ext (if required)
    - kbp_result.descr.ph_values (if required)
    - kbp_result.occs
    - kbp_result.conformer_occs
    """

    result_files = os.listdir(folder)

    occupancies = defaultdict(dict)

    reg = re.compile(r'^%s_pH([-.\d]+)_[0.]+mV_300K$' % output_prefix)
    for result_file in result_files:
        reg_m = reg.match(result_file)
        if reg_m is not None:
            ph = float(reg_m.groups()[0])

            f = open(folder + result_file)
            for line in f:
                entries = line.split()

                nr_of_states = (len(entries) - 2) / 3
                if len(entries) > 2:
                    residue = entries[1]
                    states_occ = []
                    for i in range(nr_of_states):
                        occ = entries[2 + i*3]
                        states_occ.append(occ)
                    occupancies[residue][ph] = states_occ

            f.close()


    conf_occupancies = defaultdict(dict)

    reg = re.compile(r'^%s_pH([-.\d]+)_[0.]+mV_300K_conf(\d+)$' % output_prefix)
    for result_file in result_files:
        reg_m = reg.match(result_file)
        if reg_m is not None:
            ph = float(reg_m.groups()[0])
            conf_nr = int(reg_m.groups()[1])

            f = open(folder + result_file)
            line = f.readline()
            f.close()

            entries = line.split()
            conf_occupancies[conf_nr][ph] = float(entries[0])

    conf_occupancies = dict(conf_occupancies)


    if kbp_result is None:
        kbp_result = kbp2.kbp_results.KbpResult()

        # Set residue list
        residue_list = occupancies.keys()
        sorted_residue_list = kbp_result.descr.set_residue_list(residue_list)

        # Generate residue_list_ext list.
        first_ph = occupancies[sorted_residue_list[0]].keys()[0]
        residue_list_ext = []
        for residue_num, residue_kbp in enumerate(sorted_residue_list):
            resname_kbp = residue_kbp[0:residue_kbp.find('-')]
            nr_of_states = len(occupancies[residue_kbp][first_ph])
            residue_list_ext.append( (residue_num, residue_kbp, resname_kbp, nr_of_states) )
        kbp_result.descr.residue_list_ext = residue_list_ext

        # Set ph_values
        kbp_result.descr.ph_values = sorted(occupancies.values()[0].keys())
    else:
        # Check that residue list in kbp_results matches the one found in the output file.
        residue_list = occupancies.keys()
        sorted_residue_list = kbp_result.descr.sort_residue_list(residue_list)
        if sorted_residue_list != kbp_result.descr.sorted_residue_list:
            error = "The residue list in the provided KbpResult object does not match the residues in the " \
                    "Karlsberg output files."
            raise(AssertionError(error))

        # Check that the ph values in kbp_results matches the one found in the output file.
        if kbp_result.descr.ph_values != sorted(occupancies.values()[0].keys()):
            error = "The ph list in the provided KbpResult object does not match the ph values in the " \
                    "Karlsberg output files."
            raise(AssertionError(error))

    # Set occupancies
    ph_values = kbp_result.descr.ph_values
    kbp_result.occs = []
    for (residue_num, residue_kbp, resname_kbp, nr_of_states) in kbp_result.descr.residue_list_ext:
        occs = np.zeros((nr_of_states, len(ph_values)), dtype=np.float)
        for state in range(nr_of_states):
            for ph_nr, ph in enumerate(ph_values):
                occs[state, ph_nr] = occupancies[residue_kbp][ph][state]

        kbp_result.occs.append(occs)
    kbp_result.occs = np.array(kbp_result.occs, dtype=object)

    kbp_result.conformer_occs = np.zeros((len(ph_values), len(conf_occupancies)), dtype=np.float)
    for conf in range(len(conf_occupancies)):
        for ph_nr, ph in enumerate(ph_values):
            # The numbering of conformations in karlsberg output files starts with one.
            kbp_result.conformer_occs[ph_nr, conf] = conf_occupancies[conf+1][ph]


    return kbp_result



if __name__ == '__main__':
    pkint_filename = '/scratch/scratch/tmeyer/tmp/tapbs_test/first_run.pkint'
    g_filename = '/scratch/scratch/tmeyer/tmp/tapbs_test/first_run.g'

    workdir = '/scratch/scratch/tmeyer/tmp/kb_test/'

    titratable_yaml = '/scratch/scratch/tmeyer/kbplus2/titratable.yaml'
    titratable_residues = kbp2.kbp_tools.parse_titratable_yaml(titratable_yaml)

    #########################################
    ### Test one: pkint and g files given ###
    #########################################
    run_karlsberg(pkint_filename, g_filename, folder=workdir)
    kbp_results = parse_karlsberg_results(workdir)

    kbp_results.descr.titratable_residues = titratable_residues

    kbp_results.find_pkas()

    print kbp_results.pkas
    # print kbp_results.conformer_occs


    # ####################################################
    # ### Test two: pkint and g files given as objects ###
    # ####################################################
    workdir = '/scratch/scratch/tmeyer/tmp/kb_test2/'
    pkint, g_matrix_np, residue_list_ext = kbp2.kbp_tools.parse_g_pkint(pkint_filename, g_filename)
    residue_list = [x[1] for x in residue_list_ext]
    run_karlsberg([pkint], [g_matrix_np], folder=workdir, titratable_residues=titratable_residues,
                  residue_list=residue_list)

    kbp_results = parse_karlsberg_results(workdir)

    kbp_results.descr.titratable_residues = titratable_residues

    kbp_results.find_pkas()

    print kbp_results.pkas
    # print kbp_results.conformer_occs


    #################################################
    ### Test three: Like test one but in parallel ###
    #################################################
    run_karlsberg_parallel(8, pkint_filename, g_filename)
    kbp_results = parse_karlsberg_results(workdir)
    kbp_results.descr.titratable_residues = titratable_residues

    kbp_results.find_pkas()

    print kbp_results.pkas
    # print kbp_results.conformer_occs
