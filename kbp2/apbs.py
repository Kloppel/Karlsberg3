# -*- coding: utf-8 -*-

import subprocess
import os
import stat
import re

import file_parser


def create_apbs_input(jobname, run_folder, structure, apbs_bin, coulomb_bin, target_res=0.3, verbose=True, ion_conc=0.1,
                tmp_run_folder=None, cavity_par=[]):
    """
    Preapares a apbs job. That includes a config file for APBS and a shell script to submit.

    Parameters:
    jobname:        The name used for the pqr file and for the workdir on the dollies.
    run_folder:     A folder to place the config script, the shell script and the structure.
    structure:      A file_parser.Simple_struct_parser() object for the structure.
    target_res:     The minimum resolution. The grid size is chosen to get a equal or larger resolution. Also the maximum
                    grid size is 289!
    verbose:        Print information about the chosen grid size in the shell.
    """

    if run_folder[-1] != '/':
        run_folder += '/'

    if tmp_run_folder is None:
        tmp_run_folder = run_folder

    ### Create pqr file ###
    pqr_filename = jobname + '.pqr'
    structure.write_pqr(run_folder + pqr_filename)

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
            print("WARNING (apbs): Resolution in dimension %i is larger than desired: %f" % (i, res))
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
        mol pqr %s """ % pqr_filename + """
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
       ion charge 1  conc %.3f radius 2.0 """ % ion_conc + """
       ion charge -1 conc %.3f radius 2.0 """ % ion_conc + """
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

    if cavity_par:
        ### Prepare files for cavityfinder ###
        # PQR file in PDB style.
        if pqr_filename is None:
            cavityfinder_pqr = jobname + '_kb.pqr'
        else:
            cavityfinder_pqr = 'kb_' + pqr_filename
        structure.write_pqr(run_folder + cavityfinder_pqr, kb_style=True)

        # if water_folder is not None:
        #     # File with waters.
        #     water_filename = water_folder + pqr_filename[:-4] + '.pdb'
        #     shutil.copy(water_filename, run_folder + 'water.pdb')

        # cavities 0.2 0.7 -1.0
        # focusedCavitySearch 0.15 0.7 -1.0

        ### generate input for cavityFinder ###
        input_script = """
pqr %s  """ % cavityfinder_pqr + """
output cavities

cavities %.2f %.2f %.2f """ % (cavity_par[0], cavity_par[1], cavity_par[2]) + """
focusedCavitySearch %.2f %.2f %.2f """ % (cavity_par[0], cavity_par[1], cavity_par[2]) + """

        """
        f = open(run_folder + "/cavityFinder.in", 'w')
        f.write(input_script)
        f.close()

        cavityfinder_command = '/scratch/scratch/tmeyer/CHARMM_NAMD/cavityfinder < cavityFinder.in > cavityFinder.out'
        cavityfinder_copy_command = 'cp $JDIR/cavityFinder.out $SDIR/'
    else:
        cavityfinder_command = ''
        cavityfinder_copy_command = ''



    ### Generate shell script ###
    tcsh_script = """#!/bin/tcsh
echo "Running on `hostname`"
"""

    if tmp_run_folder != run_folder:
        tcsh_script += """
set JOBNAME=%s """ % jobname    + """
set SDIR=%s    """ % run_folder + """
set JDIR=%s$JOBNAME  """ % tmp_run_folder + """

mkdir -p $JDIR
cp $SDIR/* $JDIR
cd $JDIR"""
    else:
        tcsh_script += """
cd %s """ % run_folder + """
        """

    tcsh_script += """
%s                    """ %  cavityfinder_command   + """
%s apbs.in > apbs.out """ %  apbs_bin               + """
%s -e %s > coulomb.out   """ % (coulomb_bin, pqr_filename)

    if tmp_run_folder != run_folder:
        tcsh_script += """
cd ..
cp $JDIR/apbs.out $SDIR/
cp $JDIR/coulomb.out $SDIR/
%s                    """ %  cavityfinder_copy_command   + """
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


def submitt_apbs_job(run_folder, qsub_parameter=''):
    if run_folder[-1] != '/':
        run_folder += '/'

    shell = subprocess.Popen('csh\n',\
            stdin=subprocess.PIPE,\
            stdout=subprocess.PIPE,\
            stderr=subprocess.PIPE,\
            shell=True\
            )

    print 'qsub %s -o %s -e %s %srun.sh\n' % (qsub_parameter, run_folder, run_folder, run_folder)

    shell.stdin.write('qsub %s -o %s -e %s %srun.sh\n' % (qsub_parameter, run_folder, run_folder, run_folder) )
    shell.stdin.write('exit\n')


#check again if it is ok
def run_local(run_folder, quiet=True):
    if run_folder[-1] != '/':
        run_folder += '/'

    # shell = subprocess.Popen('csh\n',\
    shell = subprocess.Popen('bash\n',\
            stdin=subprocess.PIPE,\
            stdout=subprocess.PIPE,\
            stderr=subprocess.PIPE,\
            shell=True\
            )

    shell.stdin.write('%srun.sh\n' % run_folder )
    shell.stdin.write('exit\n')

    if not quiet:
        print 'Waiting for APBS to finish..'
    while True:
        nextline = shell.stdout.readline()
        if not nextline:
            break


def read_result(run_folder, epsilon=4.0):

    if run_folder[-1] != '/':
        run_folder += '/'

    finished_files = ['apbs.out', 'coulomb.out']
    for file in finished_files:
        if not os.path.exists(run_folder + file):
            is_done = False
            break
    else:
        is_done = True

    if not is_done:
        # unfinished_jobs.append(run_f)
        return (None, )
    else:
        apbs_energy = None
        coulomb_energy = None

        f = open(run_folder + 'apbs.out')
        reg = re.compile(r'^  Global net ELEC energy =\s+([-\d\.]+)E([-+\d]+)\s+kJ/mol$')
        for line in f:
            reg_m = reg.match(line)
            if reg_m is not None:
                apbs_energy = float(reg_m.groups()[0])
                apbs_energy *= 10 ** int(reg_m.groups()[1])
                break
        f.close()

        f = open(run_folder + 'coulomb.out')
        reg = re.compile(r'^Total energy = ([-\d\.]+)e([-+\d]+) kJ/mol in vacuum.$')
        for line in f:
            reg_m = reg.match(line)
            if reg_m is not None:
                coulomb_energy = float(reg_m.groups()[0])
                coulomb_energy *= 10 ** int(reg_m.groups()[1])
                # Epsilon in solvation calculations: 4
                coulomb_energy /= epsilon
                break
        f.close()

    return apbs_energy, coulomb_energy

def start_apbs_job(jobname, run_folder, structure, apbs_bin, coulomb_bin, target_res=0.3, verbose=True, ion_conc=0.1,
                   qsub_parameter='', temp_run_folder=None, cavity_par=[]):

    if run_folder[-1] != '/':
        run_folder += '/'

    if temp_run_folder is not None and temp_run_folder[-1] != '/':
        temp_run_folder += '/'

    if os.path.exists(run_folder):
        # error = "Folder cannot be created since it already exists: %s" % run_folder
        # raise(AssertionError(error))
        return

        # create_apbs_input(jobname, run_folder, structure, target_res=target_res, verbose=verbose, conc=conc,
        #               qsub_parameter=qsub_parameter, tmp_run_folder=temp_run_folder)
        # run_local(run_folder)
    else:
        os.mkdir(run_folder)

    # Prepare pqr file
    create_apbs_input(jobname, run_folder, structure, apbs_bin, coulomb_bin, target_res=target_res, verbose=verbose, ion_conc=ion_conc,
                      tmp_run_folder=temp_run_folder, cavity_par=cavity_par)

    # Submit job
    # submitt_apbs_job(run_folder, qsub_parameter=qsub_parameter)

    run_local(run_folder)

def create_interaction_energy_apbs_input(jobname, run_folder, structure, apbs_bin, coulomb_bin, target_res=0.3, verbose=True,
                                         ion_conc=0.1, only_prepare_calcs = False):

    """
    Preapares a apbs interaction energy job. That includes a config file for APBS and a shell script to submit.

    Parameters:
    jobname:        The name used for the pqr file and for the workdir on the dollies.
    run_folder:     A folder to place the config script, the shell script and the structure.
    structure:      A file_parser.Simple_struct_parser() object for the structure.
    target_res:     The minimum resolution. The grid size is chosen to get a equal or larger resolution. Also the maximum
                    grid size is 289!
    verbose:        Print information about the chosen grid size in the shell.
    """

    if run_folder[-1] != '/':
        run_folder += '/'

    ### Create pqr file ###
    pqr_filename = jobname + '.pqr'
    structure.write_pqr(run_folder + pqr_filename)

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
        mol pqr %s """ % pqr_filename + """
    end

    elec name %s  """ % jobname + """
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
       ion charge 1  conc %.3f radius 2.0 """ % ion_conc + """
       ion charge -1 conc %.3f radius 2.0 """ % ion_conc + """
       srfm mol
       chgm spl0
       sdens 10.00
       srad 1.40
       swin 0.30
       temp 300.0
       # calcenergy total
       calcenergy comps
       calcforce no
       write pot dx %s """ % jobname + """
    end
    print elecEnergy %s end """ % jobname + """
    quit
    """
    f = open(run_folder + "/apbs.in", 'w')
    f.write(input_script)
    f.close()


    ### Generate shell script ###
    tcsh_script = """#!/bin/tcsh
echo "Running on `hostname`"
"""
    tcsh_script += """
cd %s """ % run_folder + """
    """
    tcsh_script += """
%s apbs.in > apbs.out """ %  apbs_bin               + """
%s -e %s > coulomb.out   """ % (coulomb_bin, pqr_filename)

    tcsh_script_name = run_folder + "run.sh"
    f = open(tcsh_script_name, 'w')
    f.write(tcsh_script)
    f.close()
    # Add right to execute file.
    st = os.stat(tcsh_script_name)
    os.chmod(tcsh_script_name, st.st_mode | stat.S_IEXEC)

    if not only_prepare_calcs:
        run_local(run_folder)

        # Submit job
        # submitt_apbs_job(run_folder, qsub_parameter=qsub_parameter)


def apply_potential(folder, jobname, multivalue_source, csv_file, dx_file):

    tcsh_script = """
%s %s %s %s & """ % (multivalue_source, csv_file, dx_file, folder + '%s_pot.phi' % jobname) + """
"""
    input_file =  open(folder + 'pot.sh', 'w')
    input_file.write('#!/bin/tcsh\n')
    input_file.write(tcsh_script)

    st = os.stat(folder + 'pot.sh')
    os.chmod(folder + 'pot.sh', st.st_mode | stat.S_IEXEC)

    # shell = subprocess.Popen('csh\n',\
    #         stdin=subprocess.PIPE,\
    #         stdout=subprocess.PIPE,\
    #         stderr=subprocess.PIPE,\
    #         shell=True\
    #         )
    # shell.stdin.write('%spot.sh & \n' % folder )
    # shell.stdin.write('exit\n')

def calculate_interaction_energy(folder, pqr_file, phi_file):

    potentials = []
    charges = []

    interaction_energy = 0

    f = open(phi_file, 'r')

    for line in f:
        line = re.split(r'[,]' ,line)
        potential = line[3]
        potential = float(potential)
        potentials.append(potential)

    q = open(pqr_file, 'r')
    for line in q:
        line = line.split()
        if line[0] == 'ATOM':
            charge = line[8]
            charge = float(charge)
            charges.append(charge)

    for charge, potential in zip(charges, potentials):
        multiplication = charge*potential
        interaction_energy+=multiplication

    KJ_ener =  interaction_energy*0.6*4.2
    kcal_ener = interaction_energy*0.6

    output_filepath = folder + 'interaction_energy.dat'
    out_file = open(output_filepath, 'w')
    out_file.write(str(kcal_ener))
    out_file.write(' Kcal/mol\n')
    out_file.write(str(KJ_ener))
    out_file.write(' KJ/mol\n')
    out_file.close()

    # print KJ_ener
    return KJ_ener

def pqr2csv(pqr_file):

    """    ATOM  12340  O1D PRD     3      -0.914   -0.770    9.784 -0.760 1.700"""

    f = open(pqr_file, 'r')
    q = open(pqr_file[:-4]+'.csv', 'a')

    for line in f:
        line_comps = line.split()
        if line_comps[0] == 'ATOM':
            string_to_write = '%s,%s,%s\n' % (line_comps[5], line_comps[6], line_comps[7])
            q.write(string_to_write)

    f.close()
    q.close()

    # st = os.stat(pqr_file[:-4]+'.csv')
    # os.chmod(pqr_file[:-4]+'.csv', st.st_mode | stat.S_IEXEC)

if __name__ == '__main__':

    apbs_bin    = "/scratch/scratch/tmeyer/CHARMM_NAMD/apbs"
    coulomb_bin = "/scratch/scratch/tmeyer/CHARMM_NAMD/coulomb"
    qsub_parameter = '-q D47.q,D48.q,D49.q,D50.q,D51.q,D52.q,D53.q,D54.q,D55.q,D56.q,D57.q,D58.q,D59.q,D60.q,D64.q'


    # file preparation for
    crd = '/user/jdragelj/python/CcO/ccobasic/consecutive_modelling_relax/basic_out.crd'
    psf = '/user/jdragelj/python/CcO/ccobasic/consecutive_modelling_relax/basic_out.xplor.psf'
    pdb_mod = file_parser.Simple_struct_parser()

    pdb_mod.read_par("/scratch/scratch/awoelke/md_cco/toppar/par_all22_prot_plus_heme_and_Cu.inp")
    pdb_mod.read_par("/scratch/scratch/awoelke/md_cco/toppar/par_all36_lipid.prm")

    pdb_mod.read_crd(crd)
    pdb_mod.read_xplor_psf(psf)

    new = pdb_mod.copy(segname='MEMB',exclude=True)

    # to run on dollys - jd
    # run_folder = '/scratch/scratch/jdragelj/tmp/apbs_test/run/testrun/'
    # temp_run_folder = '/public/scratch/jdragelj/kb/apbs_runs/'

    #to run on local - create whenever you want - jd
    run_folder = '/user/jdragelj/python/apbs/apbs_test/run/'
    temp_run_folder = '/user/jdragelj/python/apbs/apbs_test/tmp'

    # jobname = 'jdtest'
    jobname = 'cco_test'
    # jobname = 'mitja'


    start_apbs_job(jobname, run_folder, new, apbs_bin, coulomb_bin, target_res=0.5, verbose=True,
                   qsub_parameter=qsub_parameter, temp_run_folder=temp_run_folder)

    results = read_result(run_folder)

    print
    if len(results) == 1:
        print("No output files found. Job may not be finished yet, come back later.")
    else:
        if None in results:
            print("Job not finished yet or either the solvation or coulomb calculation failed!")
        else:
            print("Everything is fine! Here are the results:\n")
            apbs_energy, coulomb_energy = results
            print("Solvation energy: %5.3f" % apbs_energy)
            print("Coulomb energy:   %5.3f" % coulomb_energy)

