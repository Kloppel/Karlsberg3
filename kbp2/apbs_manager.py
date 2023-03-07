# -*- coding: utf-8 -*-

import os, shutil, re
import subprocess
import stat

import file_parser


def prepare_job(run_folder, jobname, structure, apbs_bin, coulomb_bin, target_res=0.3, verbose=True, pqr_filename=None,\
                water_folder=None, conc=0.1, queues=None, dolly_run_folder=None, cavity_parameter=0.9):
    """
    Preapares a apbs job. That includes a config file for APBS and a shell script to submit.

    Parameters:
    run_folder:     A folder to place the config script, the shell script and the structure.
    jobname:        The name used for the pqr file and for the workdir on the dollies.
    structure:      A file_parser.Simple_struct_parser() object for the structure.
    apbs_bin:       Path to the APBS binary.
    coulomb_bin     Path to the coulomb binary. (Part of the APBS package)
    target_res:     The minimum resolution. The grid size is chosen to get a equal or larger resolution. Also the maximum
                    grid size is 289!
    verbose:        Print information about the chosen grid size in the shell.
    pqr_filename:   If set, this name will be used for the created pqr file.
    """

    if dolly_run_folder is None:
        dolly_run_folder = '/public/scratch/tmeyer/kb/apbs/'
    if dolly_run_folder[-1] != '/':
        dolly_run_folder += '/'

    if run_folder[-1] != '/':
        run_folder += '/'
    if not os.path.isdir(run_folder):
        os.mkdir(run_folder)

    if pqr_filename is None:
        apbs_pqr = jobname + '.pqr'
    else:
        apbs_pqr = pqr_filename
    structure.write_pqr(run_folder + apbs_pqr)


    ### Determine grid size ###
    (min, max) = structure.get_size()
    size = [max[dim] - min[dim] for dim in range(3)]
    size = [x + 5.0 for x in size]
    # max_dim = 289
    dime = [0,0,0]
    for i in range(3):
        for d in [65, 97, 129, 161, 193, 225, 257, 289, 321, 353, 385, 417, 449]:
            res = size[i] / float(d)
            if res <= target_res:
                dime[i] = d
                break
        else:
            dime[i] = d
            res = size[i] / float(d)
            print("WARNING: Resolution in dimension %i is larger than desired: %f" % (i, res))
    mem = 200. / 1e9
    clen = [0,0,0]
    for i in range(3):
        mem *= dime[i]
        clen[i] = size[i] * 4.
    if verbose:
        print("A grid size of %i/%i/%i will be used." % tuple(dime))
        print("About %.1f GB of memory required per job." % mem)

    ### Generate APBS input script ###
    input_script = """
    read
        mol pqr %s """ % apbs_pqr + """
    end

    elec name water
       mg-auto
       dime %i %i %i           """ % tuple(dime) + """
       cglen %5.1f %5.1f %5.1f """ % tuple(clen) + """
       fglen %5.1f %5.1f %5.1f """ % tuple(size) + """
       cgcent mol 1
       fgcent mol 1
       mol 1
       npbe
       bcfl sdh
       pdie 4.0
       sdie 80.0
       ion charge 1  conc %.3f radius 2.0 """ % conc + """
       ion charge -1 conc %.3f radius 2.0 """ % conc + """
       srfm mol
       chgm spl0
       sdens 10.00
       srad 1.40
       swin 0.30
       temp 300.0
       # calcenergy total
       calcenergy comps
       calcforce no
    end
    elec name vacuum
       mg-auto
       dime %i %i %i           """ % tuple(dime) + """
       cglen %5.1f %5.1f %5.1f """ % tuple(clen) + """
       fglen %5.1f %5.1f %5.1f """ % tuple(size) + """
       cgcent mol 1
       fgcent mol 1
       mol 1
       npbe
       bcfl sdh
       pdie 4.0
       sdie 4.0
       ion charge 1  conc 0.000 radius 2.0
       ion charge -1 conc 0.000 radius 2.0
       srfm mol
       chgm spl0
       sdens 10.00
       srad 1.40
       swin 0.30
       temp 300.0
       # calcenergy total
       calcenergy comps
       calcforce no
    end
    print elecEnergy water - vacuum end
    quit
    """
    # write pot dx "/scratch/scratch/tmeyer/md_pka/confE_tests/runs/pot.dx"
    # write pot dx "/scratch/scratch/tmeyer/projects/cavityfinder/pot"
    # write dielx dx "/scratch/scratch/tmeyer/projects/cavityfinder/dielx"
    f = open(run_folder + "/apbs.in", 'w')
    f.write(input_script)
    f.close()

    if cavity_parameter is not None:
        ### Prepare files for cavityfinder ###
        # PQR file in PDB style.
        if pqr_filename is None:
            cavityfinder_pqr = jobname + '_kb.pqr'
        else:
            cavityfinder_pqr = 'kb_' + pqr_filename
        structure.write_pqr(run_folder + cavityfinder_pqr, kb_style=True)

        if water_folder is not None:
            # File with waters.
            water_filename = water_folder + pqr_filename[:-4] + '.pdb'
            shutil.copy(water_filename, run_folder + 'water.pdb')

        #cavities 0.2 0.7 -1.0
        #focusedCavitySearch 0.15 0.7 -1.0

        ### generate input for cavityFinder ###
        input_script = """
pqr %s  """ % cavityfinder_pqr + """
output cavities

cavities %.2f 0.25 -1.0 """ % cavity_parameter + """
focusedCavitySearch %.2f 0.25 -1.0 """ % cavity_parameter + """

    """
        f = open(run_folder + "/cavityFinder.in", 'w')
        f.write(input_script)
        f.close()

        cavityfinder_command      = '/scratch/scratch/tmeyer/CHARMM_NAMD/cavityfinder < cavityFinder.in > cavityFinder.out'
        cavityfinder_copy_command = 'cp $JDIR/cavityFinder.out $SDIR/'
    else:
        cavityfinder_command      = ''
        cavityfinder_copy_command = ''


    ### Generate shell script ###
    tcsh_script = """#!/bin/tcsh

echo "Running on `hostname`"

set JOBNAME=%s """ % jobname    + """
set SDIR=%s    """ % run_folder + """
set JDIR=%s$JOBNAME  """ % dolly_run_folder + """


mkdir -p $JDIR
cp $SDIR/* $JDIR

cd $JDIR

%s                    """ %  cavityfinder_command   + """
%s apbs.in > apbs.out """ %  apbs_bin               + """
%s -e %s > coulomb.out   """ % (coulomb_bin, apbs_pqr) + """

cd ..

#mkdir $SDIR/output_cuda
cp $JDIR/apbs.out $SDIR/
cp $JDIR/coulomb.out $SDIR/
%s                    """ %  cavityfinder_copy_command   + """
#cp $JDIR/* $SDIR/

rm $JDIR/*
rmdir $JDIR
"""
    tcsh_script_name = run_folder + "run.sh"
    f = open(tcsh_script_name, 'w')
    f.write(tcsh_script)
    f.close()
    # Add right to execute file.
    st = os.stat(tcsh_script_name)
    os.chmod(tcsh_script_name, st.st_mode | stat.S_IEXEC)

def submit_job(run_folder, queues=None):
    if run_folder[-1] != '/':
        run_folder += '/'

    shell = subprocess.Popen('csh\n',\
            stdin=subprocess.PIPE,\
            stdout=subprocess.PIPE,\
            stderr=subprocess.PIPE,\
            shell=True\
            )

    #shell.stdin.write('qsub -o %s -e %s %srun.sh\n' % (run_folder, run_folder, run_folder) )

    if queues is None:
        # shell.stdin.write('qsub -q D47.q,D48.q,D49.q,D50.q,D51.q,D52.q,D53.q,D54.q,D55.q,D56.q,D57.q,D58.q,D59.q,D60.q,D64.q -o %s -e %s %srun.sh\n' % (run_folder, run_folder, run_folder) )
        # shell.stdin.write('qsub -q D61.q,D62.q,D63.q,D64.q -o %s -e %s %srun.sh\n' % (run_folder, run_folder, run_folder) )
        shell.stdin.write('qsub -q D61.q,D62.q,D64.q -o %s -e %s %srun.sh\n' % (run_folder, run_folder, run_folder) )
        # shell.stdin.write('qsub -q D47.q -o %s -e %s %srun.sh\n' % (run_folder, run_folder, run_folder) )
    else:
        shell.stdin.write('qsub -q %s -o %s -e %s %srun.sh\n' % (queues, run_folder, run_folder, run_folder) )
    #shell.stdin.write('qsub -q D64.q,D63.q,D62.q,D61.q -o %s -e %s %srun.sh\n' % (run_folder, run_folder, run_folder) )

    shell.stdin.write('exit\n')

def read_result(run_folder, epsilon=4.0):
    if run_folder[-1] != '/':
        run_folder += '/'

    is_valid = True
    is_done  = True
    ### Check if this is a finished or unfinished job ###
    validation_files = ['apbs.in', 'run.sh']
    for file in validation_files:
        if not os.path.exists(run_folder + file):
            is_valid = False
            break

    if not is_valid:
        print("Skipping folder '%s', since it does not contain a valid apbs job." % run_folder)
        return -1

    finished_files = ['apbs.out', 'coulomb.out']
    # finished_files = ['coulomb.out']
    for file in finished_files:
        if not os.path.exists(run_folder + file):
            is_done = False
            break

    if not is_done:
        return 0
    else:
        apbs_energy    = None
        coulomb_energy = None
        f = open(run_folder + 'apbs.out')
        reg = re.compile(r'^  Global net ELEC energy =\s+([-\d\.]+)E([-+\d]+)\s+kJ/mol$')
        for line in f:
            reg_m = reg.match(line)
            # print line
            if reg_m is not None:
                apbs_energy = float(reg_m.groups()[0])
                apbs_energy *= 10**int(reg_m.groups()[1])
                break
        f.close()
        # apbs_energy = 0.0

        f = open(run_folder + 'coulomb.out')
        reg = re.compile(r'^Total energy = ([-\d\.]+)e([-+\d]+) kJ/mol in vacuum.$')
        for line in f:
            reg_m = reg.match(line)
            # print line
            if reg_m is not None:
                coulomb_energy = float(reg_m.groups()[0])
                coulomb_energy *= 10**int(reg_m.groups()[1])
                # Epsilon in solvation calculations: 4
                coulomb_energy /= epsilon
                break
        f.close()

        return (apbs_energy, coulomb_energy)
        # if apbs_energy is not None and coulomb_energy is not None:
        #     finished_jobs[pdb] = [apbs_energy, coulomb_energy]
        # else:
        #     crashed_jobs.append(pdb)


def submit_jobs(base_folder, jobname, apbs_bin, coulomb_bin, restart=False, target_res=0.2, dolly_prefix='',
                use_water_template=False, subfolder_suffix='', get_pqr_from_md=False, par_files=[], conc=0.1,
                queues=None, dolly_run_folder=None, cavity_parameter=None):
    ### Check folder ###
    if base_folder[-1] != '/':
        base_folder += '/'

    if jobname != '':
        kbp_folder = base_folder + jobname + "kbp/"
        md_folder = base_folder + jobname + "md/"
    else:
        kbp_folder = base_folder
        md_folder = kbp_folder + '../md/'
        jobname = dolly_prefix


    frame_folder = kbp_folder + "frames/"
    if use_water_template:
        water_folder = kbp_folder + "water/"
    else:
        water_folder = None
    # frame_folder = base_folder + jobname + "/kbp/frames/"

    files = os.listdir(frame_folder)

    apbs_folder = kbp_folder + "confE%s/" % subfolder_suffix
    # apbs_folder = base_folder + jobname + "/kbp/confE/"
    run_folder  = apbs_folder + "run/"

    if restart:
        if os.path.isdir(apbs_folder):
            shutil.rmtree(apbs_folder)

    if not os.path.isdir(apbs_folder):
        os.mkdir(apbs_folder)

    pdbs = []
    reg = re.compile(r'^frame(\d+)\.pdb$')
    for file in files:
        reg_m = reg.match(file)
        # print file
        if reg_m is not None:
            number = reg_m.groups()[0]
            pdbs.append('frame' + number)

    shell = subprocess.Popen('csh\n',\
            stdin=subprocess.PIPE,\
            stdout=subprocess.PIPE,\
            stderr=subprocess.PIPE,\
            shell=True\
            )

    first = True
    for pdb in pdbs:
        ### Create a run folder. ###
        if not os.path.isdir(run_folder):
            os.mkdir(run_folder)

        run_f = run_folder + pdb + '/'
        if os.path.isdir(run_f):
            continue

        if get_pqr_from_md:
            #frame_folder
            pdb_filename = frame_folder + pdb + '.pdb'
            xplor_psf_filename = md_folder + 'protein_in_water.xplor.psf'
            s = file_parser.Simple_struct_parser()
            s.read_pdb(pdb_filename)
            s.read_xplor_psf(xplor_psf_filename, allow_missmatch=True)
            for par in par_files:
                s.read_par(par)

            atoms_to_delete = []
            for atom in s.atoms:
                if atom['resid'] == 1:
                    if atom['name'] not in ['C', 'O']:
                        atoms_to_delete.append(atom)
                if atom['resid'] == 3:
                    if atom['name'] not in ['N', 'HN', 'CA']:
                        atoms_to_delete.append(atom)
            s.del_atom(atoms_to_delete)

        else:
            ### Get Karlsberg+ PQR and convert in APBS PQR file. ###
            kbp_pqr  = "%sdone/%s/c_pH7_%s.reference.pqr" % (kbp_folder, pdb, pdb)
            #apbs_pqr = pdb + ".pqr"
            if not os.path.isfile(kbp_pqr):
                # print("Skipping '%s'. File not found in Karlsberg+ result folder:\n %s\n" % (pdb, kbp_pqr))
                continue

            os.mkdir(run_f)

            s = file_parser.Simple_struct_parser()
            s.read_pdb(kbp_pqr, is_pqr=True)


        pqr_filename = pdb + '.pqr'
        verbose = first
        first = False
        prepare_job(run_f, jobname+'_'+pdb, s, apbs_bin, coulomb_bin, pqr_filename=pqr_filename, verbose=verbose,
                    target_res=target_res, water_folder=water_folder, conc=conc, queues=queues,
                    dolly_run_folder=dolly_run_folder, cavity_parameter=cavity_parameter)


        if queues is None:
            #shell.stdin.write('qsub  -o %s -e %s %srun.sh\n' % (run_f, run_f, run_f) )
            shell.stdin.write('qsub -q D48.q,D49.q,D50.q,D51.q,D52.q,D53.q,D54.q,D55.q,D56.q,D57.q,D58.q,D59.q,D60.q,D64.q -o %s -e %s %srun.sh\n' % (run_f, run_f, run_f) )
            #shell.stdin.write('qsub -q D64.q,D63.q,D62.q,D61.q -o %s -e %s %srun.sh\n' % (run_f, run_f, run_f) )
        else:
            shell.stdin.write('qsub -q %s -o %s -e %s %srun.sh\n' % (queues, run_f, run_f, run_f) )

    shell.stdin.write('exit\n')


def read_results(base_folder, jobname='', subfolder_suffix='', epsilon=4.0, extendend=False):
    if extendend:
        return read_results2(base_folder, jobname, subfolder_suffix, epsilon, extendend)

    ### Check folder ###
    if base_folder[-1] != '/':
        base_folder += '/'

    if jobname != '':
        apbs_folder = base_folder + jobname + "/kbp/confE%s/" % subfolder_suffix
    else:
        apbs_folder = base_folder + "confE%s/" % subfolder_suffix
    run_folder  = apbs_folder + "run/"

    if os.path.exists(run_folder):
        folders = os.listdir(run_folder)
    else:
        folders = []

    pdbs = []
    reg = re.compile(r'^frame(\d+)$')
    for folder in folders:
        if not os.path.isdir(run_folder + folder):
            continue
        reg_m = reg.match(folder)
        if reg_m is not None:
            number = reg_m.groups()[0]
            pdbs.append('frame' + number)

    unfinished_jobs = []
    finished_jobs   = {}
    crashed_jobs    = []

    first = True
    for pdb in pdbs:
        run_f = run_folder + pdb + '/'

        is_valid = True
        is_done  = True
        ### Check if this is a finished or unfinished job ###
        validation_files = ['apbs.in', 'run.sh']
        for file in validation_files:
            if not os.path.exists(run_f + file):
                is_valid = False
                break

        if not is_valid:
            # print("Skipping folder '%s', since it does not contain a valid apbs job." % run_f)
            continue

        # finished_files = ['apbs.out', 'coulomb.out']
        finished_files = ['coulomb.out']
        for file in finished_files:
            if not os.path.exists(run_f + file):
                is_done = False
                break

        if not is_done:
            # unfinished_jobs.append(run_f)
            unfinished_jobs.append(pdb)
        else:
            apbs_energy    = None
            coulomb_energy = None
            if os.path.exists(run_f + 'apbs.out'):
                f = open(run_f + 'apbs.out')
                reg = re.compile(r'^  Global net ELEC energy =\s+([-\d\.]+)E([-+\d]+)\s+kJ/mol$')
                for line in f:
                    reg_m = reg.match(line)
                    if reg_m is not None:
                        apbs_energy = float(reg_m.groups()[0])
                        apbs_energy *= 10**int(reg_m.groups()[1])
                        break
                f.close()
            # else:
            #     apbs_energy = 0.0


            f = open(run_f + 'coulomb.out')
            reg = re.compile(r'^Total energy = ([-\d\.]+)e([-+\d]+) kJ/mol in vacuum.$')
            for line in f:
                reg_m = reg.match(line)
                # print line
                if reg_m is not None:
                    coulomb_energy = float(reg_m.groups()[0])
                    coulomb_energy *= 10**int(reg_m.groups()[1])
                    # Epsilon in solvation calculations: 4
                    coulomb_energy /= epsilon
                    break
            f.close()


            if apbs_energy is not None and coulomb_energy is not None:
                finished_jobs[pdb] = [apbs_energy, coulomb_energy]
            else:
                crashed_jobs.append(pdb)


    return (unfinished_jobs, crashed_jobs, finished_jobs)


def read_results_res(base_folder, residue_list, jobname='', subfolder_suffix='', epsilon=4.0, extendend=False):
    if extendend:
        return read_results2(base_folder, jobname, subfolder_suffix, epsilon, extendend)

    import MDAnalysis
    import numpy as np
    from MDAnalysis.analysis.distances import distance_array

    ### Check folder ###
    if base_folder[-1] != '/':
        base_folder += '/'

    if jobname != '':
        apbs_folder = base_folder + jobname + "/kbp/confE%s/" % subfolder_suffix
    else:
        apbs_folder = base_folder + "confE%s/" % subfolder_suffix
    run_folder  = apbs_folder + "run/"

    if os.path.exists(run_folder):
        folders = os.listdir(run_folder)
    else:
        folders = []

    pdbs = []
    reg = re.compile(r'^frame(\d+)$')
    for folder in folders:
        if not os.path.isdir(run_folder + folder):
            continue
        reg_m = reg.match(folder)
        if reg_m is not None:
            number = reg_m.groups()[0]
            pdbs.append('frame' + number)


    # Creat a selection containing all titratable residues.
    # residue_list_sel = ''
    # for i, residue in enumerate(residue_list):
    #     resname, resid, segname = re.split(r'[-_]', residue)
    #     sel_residue = '(resid %s and segid %s)' % (resid, segname)
    #     if residue_list_sel != '':
    #         sel_residue = ' or ' + sel_residue
    #     residue_list_sel += sel_residue
    # residue_list_sel = '(%s)' % residue_list_sel

    termini = {}
    residue_list_orderd = {}
    for i, residue in enumerate(residue_list):
        resname, resid, segname = re.split(r'[-_]', residue)
        if resname in ['CTE', 'NTE']:
            termini[(int(resid), segname)] = resname
        if not residue_list_orderd.has_key(segname):
            residue_list_orderd[segname] = []
        residue_list_orderd[segname].append(resid)
    residue_list_sel = ''
    for i, segname in enumerate(residue_list_orderd.keys()):
        sel_segname = 'segid %s' % segname

        sel_resid = ''
        for j, resid in enumerate(residue_list_orderd[segname]):
            if j > 0:
                sel_resid += ' or '
            sel_resid += 'resid %s' % resid

        if i > 0:
            residue_list_sel += ' or '
        residue_list_sel += '(%s and (%s))' % (sel_segname, sel_resid)

    # print residue_list_sel
    # sys.exit()


    # Get a PDB structure. (The .pqr file does not contain segment names)
    # Coordinates are not relevant here. The topology is needed.
    # original_pqr = base_folder + jobname + "/done/%s/c_pH7_%s.pqr" % (pdbs[0], pdbs[0])

    unfinished_jobs = []
    finished_jobs   = {}
    crashed_jobs    = []

    # first = True
    for pdb in pdbs:
        run_f = run_folder + pdb + '/'

        is_valid = True
        is_done  = True
        ### Check if this is a finished or unfinished job ###
        validation_files = ['apbs.in', 'run.sh']
        for file in validation_files:
            if not os.path.exists(run_f + file):
                is_valid = False
                break

        if not is_valid:
            # print("Skipping folder '%s', since it does not contain a valid apbs job." % run_f)
            continue

        # finished_files = ['apbs.out', 'coulomb.out']
        finished_files = ['coulomb.out']
        for file in finished_files:
            if not os.path.exists(run_f + file):
                is_done = False
                break

        if not is_done:
            # unfinished_jobs.append(run_f)
            unfinished_jobs.append(pdb)
        else:
            apbs_energy    = None
            coulomb_energy = None
            apbs_energy_per_atom = {}
            all_residues_in_range = {}
            if os.path.exists(run_f + 'apbs.out'):
                f = open(run_f + 'apbs.out')

                reg = re.compile(r'^  Global net ELEC energy =\s+([-\d\.]+)[eE]([-+\d]+)\s+kJ/mol$')
                reg_atom = re.compile(r'^\s+Atom \d+:\s+([-\d\.]+)[eE]([-+\d]+)\s+kJ/mol$')
                reg_calc = re.compile(r'^CALCULATION #(\d+) \((\w+)\): MULTIGRID$')

                # current_calc_nr = None
                current_calc = None
                for line in f:
                    reg_calc_m = reg_calc.match(line)
                    if reg_calc_m is not None:
                        # current_calc_nr = int(reg_calc_m.groups()[0]) - 1
                        current_calc = reg_calc_m.groups()[1]
                        if not apbs_energy_per_atom.has_key(current_calc):
                            apbs_energy_per_atom[current_calc] = []
                        apbs_energy_per_atom[current_calc].append([])
                        continue
                    reg_atom_m = reg_atom.match(line)
                    if reg_atom_m is not None:
                        energy = float(reg_atom_m.groups()[0])
                        energy *= 10**int(reg_atom_m.groups()[1])
                        apbs_energy_per_atom[current_calc][-1].append(energy)
                        continue
                    reg_m = reg.match(line)
                    if reg_m is not None:
                        apbs_energy = float(reg_m.groups()[0])
                        apbs_energy *= 10**int(reg_m.groups()[1])
                        break

                f.close()

                solv_energy_per_atom = []
                nr_of_atoms = len(apbs_energy_per_atom['water'][0])
                for atom_index in range(nr_of_atoms):
                    water_e = apbs_energy_per_atom['water'][-1][atom_index]
                    vacuum_e = apbs_energy_per_atom['vacuum'][-1][atom_index]
                    solv_e = water_e - vacuum_e
                    solv_energy_per_atom.append(solv_e)
            # else:
            #     apbs_energy = 0.0


            f = open(run_f + 'coulomb.out')
            reg = re.compile(r'^Total energy = ([-\d\.]+)[eE]([-+\d]+) kJ/mol in vacuum.$')
            reg_atom = re.compile(r'^\s+Atom \d+:  Energy  = ([-\d\.]+)[eE]([-+\d]+) kJ/mol$')
            coulomb_energy_per_atom = []
            for line in f:
                reg_atom_m = reg_atom.match(line)
                if reg_atom_m is not None:
                    energy = float(reg_atom_m.groups()[0])
                    energy *= 10**int(reg_atom_m.groups()[1])
                    # Epsilon in solvation calculations: 4
                    energy /= epsilon
                    coulomb_energy_per_atom.append(energy)
                    continue
                reg_m = reg.match(line)
                if reg_m is not None:
                    coulomb_energy = float(reg_m.groups()[0])
                    coulomb_energy *= 10**int(reg_m.groups()[1])
                    # Epsilon in solvation calculations: 4
                    coulomb_energy /= epsilon
                    break
            f.close()



            # Get a PDB structure. (The .pqr file does not contain segment names)
            # Coordinates are not relevant here. The topology is needed.
            original_pqr = base_folder + jobname + "/done/%s/c_pH7_%s.pqr" % (pdb, pdb)

            if len(solv_energy_per_atom) != len(coulomb_energy_per_atom):
                error = "Different number of atoms found in apbs.out (%i) and coulomb.out (%i)!" \
                        % (len(solv_energy_per_atom), len(coulomb_energy_per_atom))
                raise(AssertionError(error))
            # print np.sum(solv_energy_per_atom) - apbs_energy
            # print np.sum(coulomb_energy_per_atom) - coulomb_energy

            ### Assign energies to titratable residues. ###
            # Get a structure that has all atoms that the apbs input file had, but with segname entries.
            protein = MDAnalysis.Universe(original_pqr, format='PDB')

            # Assign solvation energy to mass entry and coulomb energy to charge entry of the atoms.
            protein.atoms.set_mass(solv_energy_per_atom)
            protein.atoms.set_charge(coulomb_energy_per_atom)

            cutoff = 0.0
            cutoff2 = 10.0
            # cutoff = 6.0
            # cutoff2 = 14.0
            # protein.atoms.write('/user/tmeyer/temp/test_prot.pdb')
            nr_of_residues = len(residue_list)
            residue_enery = np.zeros((nr_of_residues, 2), dtype=np.float)
            for res_nr, residue in enumerate(residue_list):

                resname, resid, segname = re.split(r'[-_]', residue)
                sel_residue = '(resid %s and segid %s)' % (resid, segname)
                # x = []
                # y = []
                # for cutoff in range(30, 400, 2):
                #     sel_around_res = '((around %.2f %s))' % (cutoff/10., sel_residue)
                #     # sphzone: uses the center of geometry as reference.
                #     # sel_around_res = '((sphzone %.2f %s))' % (cutoff/10., sel_residue)
                #     atoms_around_res = protein.select_atoms(sel_around_res)
                #     atoms_around_res.write('/user/tmeyer/temp/test_%s.pdb' % residue)
                #     de_solv = atoms_around_res.totalMass()
                #     e_coulomb = atoms_around_res.totalCharge()
                #     y.append(de_solv + e_coulomb)
                #     x.append(cutoff/10.)
                # from matplotlib import pyplot as plt
                # plt.figure()
                # plt.plot(x, y, 'rx-')
                # x = []
                # y = []
                # for cutoff in range(30, 400, 2):
                #     # sel_around_res = '((around %.2f %s))' % (cutoff/10., sel_residue)
                #     # sphzone: uses the center of geometry as reference.
                #     sel_around_res = '((sphzone %.2f %s))' % (cutoff/10., sel_residue)
                #     atoms_around_res = protein.select_atoms(sel_around_res)
                #     atoms_around_res.write('/user/tmeyer/temp/test_%s.pdb' % residue)
                #     de_solv = atoms_around_res.totalMass()
                #     e_coulomb = atoms_around_res.totalCharge()
                #     y.append(de_solv + e_coulomb)
                #     x.append(cutoff/10.)
                # from matplotlib import pyplot as plt
                # # plt.figure()
                # plt.plot(x, y, 'bx-')
                # plt.show()

                # sel_around_res = '((around %.2f %s))' % (cutoff/10., sel_residue)
                # sphzone: uses the center of geometry as reference.


                # sel_around_res = '((sphzone %.2f %s))' % (cutoff, sel_residue)
                #
                # atoms_around_res = protein.select_atoms(sel_around_res)
                # de_solv = atoms_around_res.totalMass()
                # e_coulomb = atoms_around_res.totalCharge()
                # residue_enery[res_nr, 0] = de_solv
                # residue_enery[res_nr, 1] = e_coulomb
                ref_coord = protein.select_atoms(sel_residue).centerOfGeometry()
                ref_coord = np.array([ref_coord], dtype=np.float32)


                # sel_around_res = '((sphlayer %.2f %.2f %s))' % (cutoff, cutoff2, sel_residue)
                # sel_around_res = '((sphzone %.2f %s))' % (cutoff2, sel_residue)
                sel_around_res = '((point %f %f %f %.4f))' \
                                 % (ref_coord[0][0], ref_coord[0][1], ref_coord[0][2], cutoff2)
                atoms_around_res = protein.select_atoms(sel_around_res)


                # comp_coord = atoms_around_res.coordinates()
                # distances = distance_array(ref_coord, comp_coord)[0]
                # weight = distances
                # weight -= cutoff
                # weight /= -(cutoff2 - cutoff)
                # weight += 1
                # weight = weight.clip(-np.inf, 1.0)
                # # assert weight.min() >= 0.0
                # # assert weight.max() <= 1.0

                # distances = distances.clip(cutoff, np.inf)
                # distances /= cutoff
                # distances -= 1
                # distances /= cutoff2
                # distances += 1
                # sys.exit()
                # from matplotlib import pyplot as plt
                # plt.figure()
                # plt.plot(distances, weight, 'xr')
                # plt.show()

                # de_solv = np.sum(atoms_around_res.masses() * weight)
                # e_coulomb = np.sum(atoms_around_res.charges() * weight)
                de_solv = np.sum(atoms_around_res.masses())
                e_coulomb = np.sum(atoms_around_res.charges())

                # residue_enery[res_nr, 0] = len(atoms_around_res.masses())
                # residue_enery[res_nr, 0] = de_solv
                # residue_enery[res_nr, 1] = e_coulomb
                nr_of_contributing_atoms = len(atoms_around_res.masses().nonzero()[0])
                # nr_of_contributing_atoms = len(atoms_around_res.masses())
                # nr_of_contributing_atoms = 1.0

                nr_of_residue_atoms = len(protein.select_atoms(sel_residue).atoms)
                # nr_of_residue_atoms = 1.0
                residue_enery[res_nr, 0] = de_solv / nr_of_contributing_atoms * nr_of_residue_atoms
                # / len(atoms_around_res.masses())
                residue_enery[res_nr, 1] = e_coulomb / nr_of_contributing_atoms * nr_of_residue_atoms
                # * nr_of_residue_atoms

                # Collect all other titratable residues within the range (with gc <-> gc distance)
                # residue_list_sel_in_sphere = '(%s and (point %f %f %f %.4f))' \
                #                              % (residue_list_sel, ref_coord[0][0], ref_coord[0][1], ref_coord[0][2],
                #                                 cutoff2)
                residue_list_sel_in_sphere1 = '(point %f %f %f %.4f)' \
                                             % (ref_coord[0][0], ref_coord[0][1], ref_coord[0][2], cutoff2)
                residue_list_sel_in_sphere2 = residue_list_sel
                # print residue_list_sel
                # residue_list_sel_in_sphere = "around %f resid 8" % cutoff2
                # residue_list_sel_in_sphere = "point -11.43438911 -17.85067177   3.11139154 %f" % cutoff2
                # residue_list_sel_in_sphere = '(around %f %s)' % (cutoff2, residue_list_sel)
                titr_residues_in_sphere = protein.select_atoms(residue_list_sel_in_sphere1)\
                                                 .select_atoms(residue_list_sel_in_sphere2).residues
                # print titr_residues_in_sphere
                # titr_residues_in_sphere = titr_residues_in_sphere.select_atoms(residue_list_sel_in_sphere2).residues
                # print titr_residues_in_sphere
                # titr_residues_in_sphere = protein.select_atoms(residue_list_sel_in_sphere).residues
                residues_in_range = []
                for res in titr_residues_in_sphere:
                    cog_coord = np.array(res.centerOfGeometry(), dtype=np.float32)
                    dist = np.linalg.norm(cog_coord - ref_coord[0])
                    if dist < cutoff:
                        weight = 1.0
                    else:
                        weight = (dist - cutoff) / -(cutoff2-cutoff) + 1.0
                        weight = weight.clip(0.0, 1.0)
                        if weight < 0.01:
                            continue

                    resid_comp = res.resnum
                    segname_comp = res.segments[0].name
                    if termini.has_key((resid_comp, segname_comp)):
                         resname_comp = termini[(resid_comp, segname_comp)]
                    else:
                        resname_comp = res.name
                    residue_comp = "%s-%i_%s" % (resname_comp, resid_comp, segname_comp)
                    residues_in_range.append((residue_comp, weight))
                    # if resid_comp == 141:
                    #     print "%.2f -> %.2f (%s)" % (dist, weight, residue_comp)

                all_residues_in_range[residue] = residues_in_range

                # if pdb == 'frame0':
                #     jobid = base_folder.split('/')[-3]
                #     atoms_around_res.write('/user/tmeyer/temp/test_%s_%s.pdb' % (jobid, residue))
                #     protein.atoms.write('/user/tmeyer/temp/test_prot_%s.pdb' % jobid)

            if apbs_energy is not None and coulomb_energy is not None:
                finished_jobs[pdb] = [apbs_energy, coulomb_energy, residue_enery, all_residues_in_range]
            else:
                crashed_jobs.append(pdb)

    return unfinished_jobs, crashed_jobs, finished_jobs


def parse_apbs_out(filename):
    apbs_energy = None
    f = open(filename)
    reg = re.compile(r'^  Global net ELEC energy =\s+([-\d\.]+)E([-+\d]+)\s+kJ/mol$')
    for line in f:
        reg_m = reg.match(line)
        if reg_m is not None:
            apbs_energy = float(reg_m.groups()[0])
            apbs_energy *= 10**int(reg_m.groups()[1])
            break
    f.close()
    return apbs_energy

def parse_coulomb_out(filename):
    coulomb_energy = None
    f = open(filename)
    reg = re.compile(r'^Total energy = ([-\d\.]+)e([-+\d]+) kJ/mol in vacuum.$')
    for line in f:
        reg_m = reg.match(line)
        # print line
        if reg_m is not None:
            coulomb_energy = float(reg_m.groups()[0])
            coulomb_energy *= 10**int(reg_m.groups()[1])
            # Epsilon in solvation calculations: 4
            #coulomb_energy /= epsilon
            break
    f.close()

    return coulomb_energy


def read_results2(base_folder, jobname='', subfolder_suffix='', epsilon=4.0, extended=False):

    raise AssertionError("This function outdated and should not be used. The development is continued in 'apbs_manager_exp2'.")

    ### Check folder ###
    if base_folder[-1] != '/':
        base_folder += '/'

    if jobname != '':
        apbs_folder = base_folder + jobname + "/kbp/confE%s/" % subfolder_suffix
    else:
        apbs_folder = base_folder + "confE%s/" % subfolder_suffix
    run_folder  = apbs_folder + "run/"

    if os.path.exists(run_folder):
        folders = os.listdir(run_folder)
    else:
        folders = []

    pdbs = []
    reg = re.compile(r'^frame(\d+)$')
    for folder in folders:
        if not os.path.isdir(run_folder + folder):
            continue
        reg_m = reg.match(folder)
        if reg_m is not None:
            number = reg_m.groups()[0]
            pdbs.append('frame' + number)

    unfinished_jobs = []
    finished_jobs   = {}
    crashed_jobs    = []

    first = True
    sites = None
    for pdb in pdbs:
        run_f = run_folder + pdb + '/'

        framenr = int(pdb[5:])

        is_valid = True
        is_done  = True
        ### Check if this is a finished or unfinished job ###
        validation_files = ['apbs.in', 'run.sh']
        for file in validation_files:
            if not os.path.exists(run_f + file):
                is_valid = False
                break

        if not is_valid:
            #print("Skipping folder '%s', since it does not contain a valid apbs job." % run_f)
            continue

        # finished_files = ['apbs.out', 'coulomb.out']
        finished_files = ['coulomb.out']
        for file in finished_files:
            if not os.path.exists(run_f + file):
                is_done = False
                break

        if not is_done:
            # unfinished_jobs.append(run_f)
            unfinished_jobs.append(pdb)
        else:
            crashed = False

            output_filename = run_f + 'apbs.out'
            apbs_energy = parse_apbs_out(output_filename)
            #if os.path.exists(run_f + 'apbs.out'):
            #    f = open(run_f + 'apbs.out')
            #    reg = re.compile(r'^  Global net ELEC energy =\s+([-\d\.]+)E([-+\d]+)\s+kJ/mol$')
            #    for line in f:
            #        reg_m = reg.match(line)
            #        if reg_m is not None:
            #            apbs_energy = float(reg_m.groups()[0])
            #            apbs_energy *= 10**int(reg_m.groups()[1])
            #            break
            #    f.close()
            # else:
            #     apbs_energy = 0.0

            output_filename = run_f + 'coulomb.out'
            coulomb_energy = parse_coulomb_out(output_filename) / epsilon
            #f = open(run_f + 'coulomb.out')
            #reg = re.compile(r'^Total energy = ([-\d\.]+)e([-+\d]+) kJ/mol in vacuum.$')
            #for line in f:
            #    reg_m = reg.match(line)
            #    # print line
            #    if reg_m is not None:
            #        coulomb_energy = float(reg_m.groups()[0])
            #        coulomb_energy *= 10**int(reg_m.groups()[1])
            #        # Epsilon in solvation calculations: 4
            #        coulomb_energy /= epsilon
            #        break
            #f.close()

            if apbs_energy is None or coulomb_energy is None:
                crashed = True

            if extended:
                output_filename = run_f + 'apbs_titr_zero.out'
                apbs_energy_titr_zero = parse_apbs_out(output_filename)

                output_filename = run_f + 'coulomb_titr_zero.out'
                coulomb_energy_titr_zero = parse_coulomb_out(output_filename) / epsilon

                if apbs_energy_titr_zero is None or coulomb_energy_titr_zero is None:
                    crashed = True

                ### Parse the sites (once per MD). ###
                if sites is None:
                    files = os.listdir(run_f + '/sites/')
                    sites = []
                    for file in files:
                        # site_1_A_frame1.pqr
                        reg = re.compile(r'^site_(.+)_frame(\d+).pqr$')
                        reg_m = reg.match(file)
                        if reg_m is not None:
                            site = reg_m.groups()[0]
                            #framenr = int(reg_m.groups()[1])
                            sites.append(site)

                abps_sites = {}
                coulomb_sites = {}
                for site in sites:
                    output_filename = run_f + 'sites_run/apbs_site_%s_frame%i.out' % (site, framenr)
                    abps_site = parse_apbs_out(output_filename) / epsilon
                    abps_sites[site] = abps_site

                    output_filename = run_f + 'sites_run/coulomb_site_%s_frame%i.out' % (site, framenr)
                    coulomb_site = parse_coulomb_out(output_filename) / epsilon
                    coulomb_sites[site] = coulomb_site

                    if abps_site is None or coulomb_site is None:
                        crashed = True

            if crashed:
                crashed_jobs.append(pdb)
            else:
                if not extended:
                    finished_jobs[pdb] = [apbs_energy, coulomb_energy]
                else:
                    finished_jobs[pdb] = [apbs_energy, coulomb_energy, apbs_energy_titr_zero, coulomb_energy_titr_zero,\
                        abps_sites, coulomb_sites]


    return (unfinished_jobs, crashed_jobs, finished_jobs)

def clean_crashed_jobs(base_folder, jobname, crashed_jobs_list, subfolder_suffix=''):
    if base_folder[-1] != '/':
        base_folder += '/'

    if jobname != '':
        apbs_folder = base_folder + jobname + "/kbp/confE%s/" % subfolder_suffix
    else:
        apbs_folder = base_folder + "confE%s/" % subfolder_suffix

    run_folder = apbs_folder + "run/"
    # run_folder = base_folder + jobname + "/kbp/confE/run/"

    for job in crashed_jobs_list:
        run_f = run_folder + job
        shutil.rmtree(run_f)


def clear_base_folder(base_folder):
    if base_folder[-1] != '/':
        base_folder += '/'

    jobnames = os.listdir(base_folder)
    for jobname in jobnames:
        apbs_folder = base_folder + jobname + "/kbp/confE/run/"
        if not  os.path.exists(apbs_folder):
            continue

        apbs_jobs = os.listdir(apbs_folder)
        for apbs_job in apbs_jobs:
            apbs_job_folder = apbs_folder + apbs_job + '/'
            files_in_apbs_folder = os.listdir(apbs_job_folder)
            if len(files_in_apbs_folder) < 3:
                print("Deleting empty folder: " + apbs_job_folder)
                shutil.rmtree(apbs_job_folder)


        (unfinished_jobs, crashed_jobs, finished_jobs) = read_results(base_folder, jobname)
        if len(unfinished_jobs) + len(crashed_jobs) == 0:
            continue

        print("Clearing folder: " + base_folder + jobname)
        clean_crashed_jobs(base_folder, jobname, crashed_jobs)
        clean_crashed_jobs(base_folder, jobname, unfinished_jobs)
        print("Folders deleted: %i/%i (unfinished/crashed)" % (len(crashed_jobs), len(unfinished_jobs)) )




if __name__ == '__main__':
    #apbs_bin    = "/scratch/scratch/tmeyer/CHARMM_NAMD/apbs"
    apbs_bin    = "/scratch/scratch/tmeyer/CHARMM_NAMD/apbs-1.3-amd64/bin/apbs"
    coulomb_bin = "/scratch/scratch/tmeyer/CHARMM_NAMD/coulomb"


    base_folder = "/scratch/scratch/tmeyer/md_pka/runs/general/"
    jobnames = []
    jobnames.append("c_2lzt_charged_his_large_longer")
    #jobnames.append("arg_switched")
    #jobnames.append("asp")
    #jobnames.append("asp_switched")
    #jobnames.append("cys")
    #jobnames.append("cys_switched")
    #jobnames.append("glu")
    #jobnames.append("glu_switched")
    #jobnames.append("hsd")
    #jobnames.append("hse")
    #jobnames.append("hsp")
    #jobnames.append("lys")
    #jobnames.append("lys_switched")
    #jobnames.append("tyr")
    #jobnames.append("tyr_switched")




    print "base folder: " + base_folder
    for jobname in jobnames:
        print "### " + jobname + " ###"

        kbp_folder = base_folder + jobname + '/kbp_std_noMin_new/'
        subfolder_suffix = ''
        #par_files = []
        #par_files.append('/scratch/scratch/tmeyer/karlsbergplus/par.inp')
        #par_files.append('/scratch/scratch/tmeyer/karlsbergplus/patches.prm')
        #queues = "D64.q"
        #dolly_run_folder = '/public/scratch/sakalli/abps/'

        submit_jobs(kbp_folder, '', apbs_bin, coulomb_bin, restart=False, target_res=0.3)
        #submit_jobs(kbp_folder, '', apbs_bin, coulomb_bin, restart=False, target_res=0.05, subfolder_suffix=subfolder_suffix,
        #            get_pqr_from_md=True, par_files=par_files, conc=0.0, queues=queues, dolly_run_folder=dolly_run_folder)
        (unfinished_jobs, crashed_jobs, finished_jobs) = read_results(kbp_folder, '', subfolder_suffix=subfolder_suffix, extendend=True)
        print "finished jobs: %i" % len(finished_jobs.keys())
        print "running jobs:  %i" % len(unfinished_jobs)
        print "crashed jobs:  %i" % len(crashed_jobs)

        #clean_crashed_jobs(kbp_folder, '', crashed_jobs, subfolder_suffix=subfolder_suffix)
        # -> Das script ausdrucken lassen auf welchem Dolly es lief?

        import numpy as np

        s_frames = []
        c_frames = []
        s_nt_frames = []
        c_nt_frames = []
        s_sites_sum_frames = []
        c_sites_sum_frames = []
        final_c_frames = []
        final_s_frames = []
        final_frames = []
        for (s, c, s_nt, c_nt, s_sites, c_sites) in finished_jobs.itervalues():
            s_frames.append(s)
            c_frames.append(c)
            s_nt_frames.append(s_nt)
            c_nt_frames.append(c_nt)

            s_sites_sum = np.sum(s_sites.values())
            c_sites_sum = np.sum(c_sites.values())
            s_sites_sum_frames.append(s_sites_sum)
            c_sites_sum_frames.append(c_sites_sum)

            final = s - s_sites_sum - s_nt + c - (c_nt + c_sites_sum)
            final_frames.append(final)
            final = s - s_sites_sum - s_nt
            final_s_frames.append(final)
            final = c - (c_nt + c_sites_sum)
            final_c_frames.append(final)

        print "Sum: %.1f +- %.1f (%.1f)" % (np.average(final_frames),      np.std(final_frames),     np.std(final_frames)/np.sqrt(len(final_frames)))
        #print "s: %.1f +- %.1f (%.1f)" % (np.average(s_frames),      np.std(s_frames),     np.std(s_frames)/np.sqrt(len(s_frames)))
        #print "s_nt: %.1f +- %.1f (%.1f)" % (np.average(s_nt_frames),      np.std(s_nt_frames),     np.std(s_nt_frames)/np.sqrt(len(s_nt_frames)))
        #print "s_sites: %.1f +- %.1f (%.1f)" % (np.average(s_sites_sum_frames),      np.std(s_sites_sum_frames),     np.std(s_sites_sum_frames)/np.sqrt(len(s_sites_sum_frames)))

        list = c_frames
        print "c:       %.1f +- %.1f (%.1f)" % (np.average(list),      np.std(list),     np.std(list)/np.sqrt(len(list)))
        list = c_nt_frames
        print "c_nt:    %.1f +- %.1f (%.1f)" % (np.average(list),      np.std(list),     np.std(list)/np.sqrt(len(list)))
        list = c_sites_sum_frames
        print "c_sites: %.1f +- %.1f (%.1f)" % (np.average(list),      np.std(list),     np.std(list)/np.sqrt(len(list)))
        print
        list = final_frames
        print "Sum:   %.1f +- %.1f (%.1f)" % (np.average(list),      np.std(list),     np.std(list)/np.sqrt(len(list)))
        list = final_s_frames
        print "Sum_s: %.1f +- %.1f (%.1f)" % (np.average(list),      np.std(list),     np.std(list)/np.sqrt(len(list)))
        list = final_c_frames
        print "Sum_c: %.1f +- %.1f (%.1f)" % (np.average(list),      np.std(list),     np.std(list)/np.sqrt(len(list)))


        import matplotlib.pyplot as plt

        plt.figure()
        #c_frames -= np.average(c_frames)
        #c_nt_frames -= np.average(c_nt_frames)
        #c_sites_sum_frames -= np.average(c_sites_sum_frames)
        #plt.plot(c_sites_sum_frames, linewidth=1)
        #plt.plot(c_frames, linewidth=1)
        #plt.plot(c_nt_frames, linewidth=1)
        #plt.legend(["sites", "c","c_nt"], loc='best', prop={'size': 7}).draw_frame(False)

        final_frames -= np.average(final_frames)
        final_s_frames -= np.average(final_s_frames)
        final_c_frames -= np.average(final_c_frames)
        plt.plot(final_frames, linewidth=3)
        plt.plot(final_s_frames, linewidth=2)
        plt.plot(final_c_frames, linewidth=2)
        plt.legend(["sum", "s","c"], loc='best', prop={'size': 7}).draw_frame(False)

        #plt.legend(["sites", "c","c_nt", "c_sum"], loc='best', prop={'size': 7}).draw_frame(False)
        #plt.hist(solv_data)
        plt.show()


        #c_data  = [c    for (s, c, s_nt, c_nt, s_sites, c_sites) in finished_jobs.itervalues()]
        #s_data  = [s for (s, c, s_nt, c_nt, s_sites, c_sites) in finished_jobs.itervalues()]
        #final_data = [c + s for (s, c, s_nt, c_nt, s_sites, c_sites) in finished_jobs.itervalues()]
        #
        #import numpy as np
        #
        #print "Coulomb: %.1f +- %.1f (%.1f)" % (np.average(c_data),      np.std(c_data),     np.std(c_data)/np.sqrt(len(c_data)))
        #print "Solvation: %.1f +- %.1f (%.1f)" % (np.average(s_data), np.std(s_data),  np.std(s_data)/np.sqrt(len(s_data)))
        #print "Sum: %.1f +- %.1f (%.1f)" % (np.average(final_data),      np.std(final_data), np.std(final_data)/np.sqrt(len(final_data)))
        #print


        import matplotlib.pyplot as plt

        # plt.figure()
        # plt.hist(solv_data)
# plt.show()

















    # # folder = "/scratch/scratch/tmeyer/md_pka/md_pka_manager/snase_runs/2oxp"
    # # folder = "/scratch/scratch/tmeyer/md_pka/md_pka_manager/snase_runs/2oxp_rois"
    # base_folder = "/scratch/scratch/tmeyer/md_pka/md_pka_manager/sb_runs"
    # # folder = "/scratch/scratch/mdkurs/nils/2lzt/2lzt_glu0"
    # jobnames = []
    # jobnames.append("2lzt")
    # jobnames.append("2lzt_charged_his")
    # jobnames.append("2lzt_neutral_acids")
    # jobnames.append("2lzt_neutral_bases_charged_tyr")
    # jobnames.append("2lzt_neutral_acids_charged_his")
    # jobnames.append("2lzt_neutral_bases_charged_tyr_arg")
    # jobnames.append("2lztno_pme")
    # jobnames.append("2lzt_no_pme_const_v")
    # # jobnames.append("2lzt_first_sim")
    # jobnames.append("2lzt_no_margin")
    # jobnames.append("2lzt_hse")
    # jobnames.append("2lzt_no_pme_const_v_1bar")


    # Barnase
    # jobnames = []
    # jobnames.append("1a2p")
    # # jobnames.append("1a2p_hse")
    # jobnames.append("1a2p_charged_his")
    # jobnames.append("1a2p_neutral_bases_charged_tyr")
    # jobnames.append("1a2p_neutral_acids_charged_his")
    # jobnames.append("1a2p_neutral_bases_charged_tyr_arg")

    # CD2
    # jobnames = []
    # jobnames.append("2hng")
    # #jobnames.append("2hng_hse")
    # jobnames.append("2hng_charged_his")
    # jobnames.append("2hng_neutral_bases_charged_tyr")
    # jobnames.append("2hng_neutral_acids_charged_his")
    # jobnames.append("2hng_neutral_bases_charged_tyr_arg")

    # jobnames = []
    # jobnames.append("2lzt")
    # jobnames.append("2lzt_hse")
    # jobnames.append("2lzt_charged_his")
    # jobnames.append("2lzt_neutral_bases_charged_tyr")
    # jobnames.append("2lzt_neutral_acids_charged_his")
    # jobnames.append("2lzt_neutral_bases_charged_tyr_arg")

    # jobnames = []
    # jobnames.append("2lzt_cooled")
    # jobnames.append("2lzt_cooled_hse")
    # jobnames.append("2lzt_cooled_charged_his")
    # jobnames.append("2lzt_cooled_n_bases_c_tyr")
    # jobnames.append("2lzt_cooled_n_acids_c_his")
    # jobnames.append("2lzt_cooled_n_bases_c_tyr_arg")

    #jobnames = []
    #jobnames.append("2lzt_constr")
    # jobnames.append("2lzt_constr_hse")
    # jobnames.append("2lzt_constr_charged_his")
    # # jobnames.append("2lzt_constr_n_bases_c_tyr")
    # jobnames.append("2lzt_constr_n_acids_c_his")
    # jobnames.append("2lzt_constr_n_bases_c_tyr_arg")

    # jobnames = []
    # jobnames.append("2lzt_cons001")
    # jobnames.append("2lzt_cons001_hse")
    # jobnames.append("2lzt_cons001_charged_his")
    # jobnames.append("2lzt_cons001_n_bases_c_tyr")
    # jobnames.append("2lzt_cons001_n_acids_c_his")
    # jobnames.append("2lzt_cons001_n_acids_c_his_e")
    # jobnames.append("2lzt_cons001_n_bases_c_tyr_arg")

    # jobnames = []
    # jobnames.append("1a2p_c")
    # # jobnames.append("1a2p_hse")
    # jobnames.append("1a2p_c_charged_his")
    # jobnames.append("1a2p_c_n_bases_c_tyr")
    # jobnames.append("1a2p_c_n_acids_c_his")
    # jobnames.append("1a2p_c_n_bases_c_tyr_arg")

    # jobnames = []
    # jobnames.append("2hng_c")
    # # jobnames.append("2hng_hse")
    # jobnames.append("2hng_c_charged_his")
    # jobnames.append("2hng_c_n_bases_c_tyr")
    # jobnames.append("2hng_c_n_acids_c_his")
    # jobnames.append("2hng_c_n_bases_c_tyr_arg")


    #base_folder = "/scratch/scratch/tmeyer/md_pka/md_pka_manager/sb_runs_2/"





    # pdbid = '1tqo'
    # exp_pkas = {'GLU-92_A' : 8.7}
    # pdbid = '2oeo'
    # exp_pkas = {'ASP-92_A' : 8.1}
    # pdbid = '2oxp'
    # exp_pkas = {'ASP-66_A' : 8.7}
    # pdbid = '2rbm'
    # exp_pkas = {'LYS-72_A' : 8.6}
    # pdbid = '2rks'
    # exp_pkas = {'LYS-38_A' : 10.4}
    # pdbid = '3c1e'
    # exp_pkas = {'LYS-125_A' : 6.2}



    # jobnames = []
    # jobnames.append(pdbid)
    # jobnames.append(pdbid + "_rois")
    #
    # base_folder = "/scratch/scratch/tmeyer/md_pka/md_pka_manager/snase_runs/"




    # jobnames = []
    # EPP-72_A:
    # jobnames.append("3ero")
    # jobnames.append("3ero_rois")
    # LYS-38_A:
    # jobnames.append("2rks")
    # jobnames.append("2rks_rois")

    # base_folder = "/scratch/scratch/tmeyer/md_pka/md_pka_manager/snase_runs/"




    # print "base folder: " + base_folder
    # for jobname in jobnames:
    #     print "### " + jobname + " ###"
    #
    #     submit_jobs(base_folder, jobname, apbs_bin, coulomb_bin, restart=False, target_res=0.3)
    #     (unfinished_jobs, crashed_jobs, finished_jobs) = read_results(base_folder, jobname)
    #     print "finished jobs: %i" % len(finished_jobs.keys())
    #     print "running jobs:  %i" % len(unfinished_jobs)
    #     print "crashed jobs:  %i" % len(crashed_jobs)
    #
    #     clean_crashed_jobs(base_folder, jobname, crashed_jobs)
    #     # -> Das script ausdrucken lassen auf welchem Dolly es lief?
    #
    #     c_data     = [c    for (solv, c) in finished_jobs.itervalues()]
    #     solv_data  = [solv for (solv, c) in finished_jobs.itervalues()]
    #     final_data = [c + solv for (solv, c) in finished_jobs.itervalues()]
    #
    #     import numpy as np
    #
    #     print "Coulomb: %.1f +- %.1f (%.1f)" % (np.average(c_data),      np.std(c_data),     np.std(c_data)/np.sqrt(len(c_data)))
    #     print "Solvation: %.1f +- %.1f (%.1f)" % (np.average(solv_data), np.std(solv_data),  np.std(solv_data)/np.sqrt(len(solv_data)))
    #     print "Sum: %.1f +- %.1f (%.1f)" % (np.average(final_data),      np.std(final_data), np.std(final_data)/np.sqrt(len(final_data)))
    #     print
    #
    #
    #
    #     # submit_apolar_jobs(base_folder, jobname, apbs_bin, coulomb_bin, restart=False)
    #     # (unfinished_jobs, crashed_jobs, finished_jobs) = read_apolar_results(base_folder, jobname)
    #
    #     import matplotlib.pyplot as plt
    #
    #     # plt.figure()
    #     # plt.hist(solv_data)
    #     # plt.show()



    #clear_base_folder("/scratch/scratch/tmeyer/md_pka/runs/general")

