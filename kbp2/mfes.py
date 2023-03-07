# coding=utf-8

import re
import subprocess
import os
import stat
import shutil


#BRAINASTORMING
#
# 1. MODULE TO CREATE CHARMM INPUT WITH MFES
#     - chamrm structure is generated separately and mfes just addes commands
#     - copies files to the folder
#     - writes new files with t and meshsize modified into that specific folder
#
# 2. SUBMIT mfes claculation
#     - in progress
#
# 3. monitoring for jobs!
#      - since we are running separate external_scripts, best is to actually
#        just iterate over folders and check for conditions from mfes copy results

#@todo: WIP, incomplete module to be done in the future

def make_new_config_file(folder, t):

    config_filename = folder + "config.in"

    config_file = """
[general]
jobname = input
molecule = input.pqr
mode = energy

[experiment]
eps_in = 4
eps_out = 80
probe_radius = 1.4
cavity = no
ionc = 0.1
ionr = 2
exclusion_probe_radius = 2

[pka]
st_folder = ./ST/
sites_file = p2lzt_pH7_0.sites
calc_cte = no
calc_nte = no
explicit_models = yes

[model]
generator = standard
grid_resolution = 512
smoothing = t %i """ %t + """
boundary = boundary_very_coarse.vol
refine_file =
surface_stl = molecule_t30_512.stl
volume_vol  = protein.vol
debug = no

generator_residue = voxelizer
grid_residue_resolution = 64
smoothing_residue = t 10

[meshing]
molecule_surface = molecule_surface.opt
molecule_volume = molecule_volume.opt
boundary_volume = boundary_volume.opt
second_order_surface = no
residue_surface = residue_surface.opt
residue_volume = residue_volume.opt

[solver]
solver = mumps-inverse
solution_order = 2
maxsteps = 2

"""

    f = open(config_filename, 'w')
    f.write(config_file)
    f.close()

    return

def make_new_molecule_surface_file(folder, meshsize):

    surface_filename = folder + "molecule_surface.opt"

    surface_file = """
options.localh  1
options.meshsize  %.2f """ %meshsize+ """
options.minmeshsize  0
meshoptions.fineness  0.7
options.grading  1
options.optsteps2d  0
options.optsteps3d  0
"""


    f = open(surface_filename, 'w')
    f.write(surface_file)
    f.close()

    return

def start_mfes_with_charmm(charmm_structure, mfes_input_files, t=None, meshsize=None):

    '''mFES can be started from charmm and therefore this function just copies neeeded files and adds commdn to charmm structure! '''

    needed_folders = os.listdir(mfes_input_files)
    for filename in needed_folders:
        shutil.copy2(mfes_input_files+'/'+filename, charmm_structure.workdir + filename)
    charmm_structure.add_charmm_command('mfes sele all .and. .not. resname tip3 end', adj_task='minimize_modelled')

    if t is not None:
        # os.remove(charmm_structure.workdir + 'config.in')
        make_new_config_file(charmm_structure.workdir, t)

    if meshsize is not None:
        # os.remove(charmm_structure.workdir + 'config.in')
        make_new_molecule_surface_file(charmm_structure.workdir, meshsize)

    return charmm_structure

def check_mfes_job(folder):

    '''
    1. result.out exists and has a number -> calculation is all right!
    2. result.out exists and is empty -> calculation crashed
    3. result.out does not exist -> calculation crashed
    '''

    mfes_filename = folder + 'result.out'
    if os.path.exists(mfes_filename):
        if os.stat(mfes_filename).st_size <= 10:
            os.remove(mfes_filename)
            return False
        else:
            return True
    else:
        return False

def read_mfes_results(folder):

    ''' returns 1. real result, 2. None 3. intersected/unknown
        check for one folder where mfes runs'''

    result_status = check_mfes_job(folder)
    print folder
    if result_status:
        mfes_filename = folder + 'result.out'
        with open(mfes_filename, "r") as f:
            for line in f:
                line = line.split(" ")
                solvation_energy = float(line[0])
        return solvation_energy
    else:
        return None


def build_membrane(tickness, size_on_plane, shift, vdW, input_pqr, script_folder):

    """Usage: perl modify_pqr.pl <original pqr> <plane: xy, xz, yz>
    <membrane thickness [Angstroems]>
    <membrane size on plane [Angstroems]>
    <normal origin shift on membrane plane (absolute) [Angstroems]>
    <vdW radius of balls modelling membrane [Angstroems]> <output pqr>

(different plane choosing does not work yet)

An example of your computations is input.pqr in that same folder.
You run the script like:
perl modify_pqr.pl input.pqr xy 20 100 0 5 output.pqr

It takes input.pqr and generates output.pqr"""


    filename = script_folder + '/modify_pqr.pl'

    if os.path.exists(filename):
        os.remove(filename)
    if not os.path.exists(script_folder):
        os.mkdir(script_folder)

    scripta = open(filename, 'w')

    file_text = """
#!/usr/bin/perl
print "Hello World from modify_pqr!\n";

$num_args = $#ARGV + 1;
if ($num_args != 7) {
  print "\n Usage: perl modify_pqr.pl <original pqr> <plane: xy, xz, yz> <membrane thickness [Angstroems]> <membrane size on plane [Angstroems]> <normal origin shift on membrane plane (absolute) [Angstroems]> <vdW radius of balls modelling membrane [Angstroems]> <output pqr>\n";
  exit;
}

$origPQR   = $ARGV[0];
$plane     = $ARGV[1];
$mThick    = $ARGV[2];
$mSize     = $ARGV[3];
$mPos      = $ARGV[4];
$ballR     = $ARGV[5];
$outputPQR = $ARGV[6];

if ($plane ne 'xy' and $plane ne 'xz' and $plane ne 'yz' ){
    print "wrong plane chosen: $plane\n";
    print "options are: xy, xz or yz\n";
    exit;
}

$ballDist = $ballR;
$output = '';

for (my $x = -1*int($mSize*0.5); $x <= int($mSize*0.5+0.99); $x=$x+$ballDist){
    for (my $y = -1*int($mSize*0.5); $y <= int($mSize*0.5+0.99); $y=$y+$ballDist){
	for (my $z = -1*int($mThick*0.5); $z <= int($mThick*0.5+0.99); $z+=$ballDist){
	    $output .= sprintf("ATOM  12515 MEM  MEM A 999    %8.3f%8.3f%8.3f 0.000 %5.3f     A\n", $x, $y, $z+$mPos, $mThick);
	}
    }
}

$cmd =`cp $origPQR $outputPQR`;

my $OUTFILE;

open $OUTFILE, '>>', $outputPQR;

print { $OUTFILE } $output;

close $OUTFILE;

print "$outputPQR successfully written!";
exit;
    """
    scripta.write(file_text)
    scripta.close()

    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IEXEC)



    # TODO
    import subprocess
    shell = subprocess.Popen('bash\n',\
            stdin=subprocess.PIPE,\
            stdout=subprocess.PIPE,\
            stderr=subprocess.PIPE,\
            shell=True\
            )

    shell.stdin.write('cd ' + script_folder + '\n')
    shell.stdin.write("perl modify_pqr.pl %s xy %i %i %i %i %i_%i_%i_%i.pqr\n" % (input_pqr, tickness, size_on_plane, shift, vdW, tickness, size_on_plane, shift, vdW))
    shell.stdin.write('exit\n')
    # if not quiet:
    #     print 'Waiting for TAPBS to finish..'
    # while True:
    #     nextline = shell.stdout.readline()
    #     if not quiet:
    #         print nextline
    #     if not nextline:
    #         break



def insert_membrane_mfes(tickness, size_on_plane, shift, vdW, source_folder):

    """ not working properly"""

    filename = source_folder + '/mfes_script'

    if os.path.exists(filename):
        os.remove(filename)
    if not os.path.exists(source_folder):
        os.mkdir(source_folder)

    scripta = open(filename, 'w')

    file_text = """
#! /bin/bash
mfes-0.3c.x86_64 --ini config.in -c

cp input.pqr input_orig.pqr
perl /scratch/scratch/jdragelj/membrane_builder/modify_pqr.pl input_orig.pqr xy %i %i %i %i input.pqr """  % (tickness, size_on_plane, shift, vdW) + """"

mfes-0.3c.x86_64 --ini config.in
cat result.out | awk '{ split($0,a," "); print a[5]" "a[5]" "a[5]" "a[5] }' > r\
esult_temp
mv result_temp result.out
    """
    scripta.write(file_text)
    scripta.close()

    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IEXEC)

def read_conf_ener_result(run_folder, epsilon=4.0):

    if run_folder[-1] != '/':
        run_folder += '/'



    ### Parse mfes output ###
    mfes_energy = None


    ### Parse Coulomb output file ###\
    coulomb_energy = None
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

    return mfes_energy, coulomb_energy


def calc_prot_ener(run_folder, jobname, residue_list, structure, titratable_residues, mfes_settings, pka_cycle0, pdie=4, sdie=80,
                   ion_conc=0.1):
    """
    run mfes protonation calculation
    """

    #todo: idea: write config.ini out of mfes_settings
    raise AssertionError('Work in progress!')

    if run_folder[-1] != '/':
        run_folder += '/'

    if os.path.exists(run_folder):
        error = "Folder %s does exist." % run_folder
        raise AssertionError(error)
    else:
        os.mkdir(run_folder)

    # mfes_bin = mfes_settings['bin']
    # mfes_bin = 'mfes_pka_script'
    mfes_bin = 'mfes-0.3c.x86_64'

    ### Make sure residue_list entries are of the format ('ARG", 23, 'A') instead of ARG-23_A
    if type(residue_list[0]) is str:
        residue_list_descr = residue_list
        residue_list = []
        for residue_descr in residue_list_descr:
            resname, resid, segname = re.split(r'[-_]', residue_descr)
            resid = int(resid)
            residue_list.append((resname, resid, segname))


    ### Write .st and .sites files
    sites_file_str = ''
    resnames = []
    for resname, resid, segname in residue_list:
        st_filename = resname + '.st'
        if resname not in resnames:
            resnames.append(resname)
            # Create .st file for this residue
            st_file_str = ''
            res_def = titratable_residues[resname]

            if pka_cycle0:
                if resname not in pka_cycle0:
                    error = "Not all residues defined in 'pka_cycle0': No entry for residue %s" % resname
                    raise(AssertionError(error))
                if len(res_def) != len(pka_cycle0[resname]):
                    error = "Not all states are defined for residue %s in 'pka_cycle0'" % resname
                    raise(AssertionError(error))

            for state, state_def in enumerate(res_def):
                if pka_cycle0:
                    st_file_str += "%.2f pK %s %.f pK\n" % (state_def['pka'], state_def['name'], pka_cycle0[resname][state])
                else:
                    st_file_str += "%.2f pK %s\n" % (state_def['pka'], state_def['name'])
                for atom_nr, atom_name in enumerate(state_def['atoms']):
                    charge = state_def['atoms'][atom_name]
                    st_file_str += "ATOM   %4i %4s %3s A   1    9999.9999999.9999999.999%6.3f99.999      A\n" \
                        % (atom_nr, atom_name, resname, charge)
            f = open(run_folder + st_filename, 'w')
            f.write(st_file_str)
            f.close()
        sites_file_str += "%s %i %s %s\n" % (segname, resid, resname, st_filename)
    sites_filename = jobname + '.sites'
    f = open(run_folder + sites_filename, 'w')
    f.write(sites_file_str)
    f.close()

    ### Write pqr file
    pqr_filename = 'input.pqr'
    structure.write_pqr(run_folder + pqr_filename, kb_style=True)


    mfes_input_files = "/scratch/scratch/jdragelj/mfes_pka_files"
    needed_folders = os.listdir(mfes_input_files)
    for filename in needed_folders:
        shutil.copy2(mfes_input_files+'/'+filename, run_folder+filename)

    ### run mfes ###

    input_filename = 'config.in'
    output_filename  = jobname + '_mfes.out'

    shell = subprocess.Popen('bash\n',\
            stdin=subprocess.PIPE,\
            stdout=subprocess.PIPE,\
            stderr=subprocess.PIPE,\
            shell=True\
            )

    shell.stdin.write('cd ' + run_folder + '\n')
    shell.stdin.write("%s --ini %s > %s\n" % (mfes_bin, input_filename, output_filename))
    shell.stdin.write('exit\n')


    quiet = mfes_settings['quiet_mode']

    if not quiet:
        print 'Waiting for mFES to finish..'

    while True:
        nextline = shell.stdout.readline()
        if not quiet:
            print nextline
        if not nextline:
            break

    shutil.copy2(run_folder+'config.pkint', run_folder+jobname+'.pkint')
    shutil.copy2(run_folder+'config.g', run_folder+jobname+'.g')



if __name__ == '__main__':

    input_pqr = '/scratch/scratch/jdragelj/membrane_builder/input.pqr'
    # script_folder = '/scratch/scratch/jdragelj/membrane_builder/'
    # source_folder = '/scratch/scratch/jdragelj/mfes_vault/'

    # #test
    # tickness = 20
    # size_on_plane = 75
    # shift = -5
    # vdW = 4
    # build_membrane(tickness, size_on_plane, shift, vdW, input_pqr, script_folder)

    # build_membrane(16, 77, -5, 5, input_pqr, script_folder)
    # build_membrane(15, 75, -5, 5, input_pqr, script_folder)
    # insert_membrane_mfes(12, 75, -5, 5, source_folder)
    # insert_membrane_mfes(18, 75, -5, 5, source_folder)