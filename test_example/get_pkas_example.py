# coding=utf-8


import os
import kbp2

# note: "../" is not accepted in python, must be replced by a proper path

# full filepath of the PDB file that will be used for titration
# IMPORTANT NOTE: Karlsberg2+ replaces water molecules with a dielectric medium of eps = 80
# IMPORTANT NOTE: Please remove (crystal) water molecules before titration
pdb_filename = '../kb2plus_package/test_example/4pti.pdb1'
# working idrectory where all files related to titration will be made
workdir = '<user folder of choice unrelated to the package>'

# premade list that will append the toplogy files needed for CHARMM27 to build the protein structure
top = []
top.append('../kb2plus_package/test_example/toppar27/top.inp')
# files patches.rtf, sb_pat.rtf are needed for any Karlsberg2+ run with CHARMM27 force field
top.append('../kb2plus_package/test_example/toppar27/patches.rtf')
top.append('../kb2plus_package/test_example/toppar27/sb_pat.rtf')
par = []
par.append('../kb2plus_package/test_example/toppar27/par.inp')
# files patches.prm is needed for any Karlsberg2+ run with CHARMM27 force field
par.append('../kb2plus_package/test_example/toppar27/patches.prm')

# modelling decision given in a certain format
# in this case, water can be renamed to TIP3 with 'rename' or kept as HOH with 'keep'
# decision are used in the module kbp2.charmm.py
decisions = []
# in this decision water molecules will be renamed
decisions.append(('rename__HOH_TIP3', 'keep'))

################
### Settings ###
################

# the settings in this section are referring to the object (with a chosen name) kbp2_set of the class PkaCalcSettings
# the class PkaCalcSettings is located in the module kbp2.pka_calculation.py
# kbp2_set stores all needed info in order to titrate the structure, technical and modelling information

# object is created, initialized from the class, ready to save settings for the titration
kbp2_set = kbp2.pka_calculation.PkaCalcSettings()
# topology files stored from the list above
kbp2_set.top = top
# paramater files stored from the list above
kbp2_set.par = par
# modelling decision stored from the list above
kbp2_set.modelling_decisions = decisions
# working directory is set
kbp2_set.workdir = workdir
# a quick check if the workfolder exists and is created if not
if not os.path.exists(workdir):
    os.mkdir(workdir)
# IMPORTANT NOTE: when re-running the computation in the same folder, please manually delete all files previously made

# protocol is defining which pH adapted conformations (PACs) will be used
# more details on PACs in the publication (G. Kieseritzky, 2008)
# pH7 is usually enough for enzymes we usually investigate as it is their active pH value
# NOTE: comment out the PACs you want to exclude
protocol = []
protocol.append((-10, 'open_sb_acids'))
protocol.append((  7, 'h_min'))
protocol.append(( 20, 'open_sb_bases'))
kbp2_set.protocol = protocol

# structure (PDB file) is set for the titration
kbp2_set.set_structure(pdb_filename)
#the choice of folder keeping, 'keep', 'remove' - suggested for a big amount of computations
kbp2_set.remove_folders = 'keep'

# preoptimization of the structure - relaxation of carbonyl oxygens because of the crystal structure imprecisions,
# and/or minimization of the structure and hydrogens in a dielectric with eps = 4
# recommended for titrations of crystal structures
kbp2_set.preopt = {'carb_oxi_relax': True, \
                    'init_die_4': True}

######################################
### Settings the binary file paths ###
######################################
# Binaries filepaths can be set in the module pka_claulations.py within the PkaCalcSettings class definition
# in that case the lines below are not needed
# However, it is recommended to define filepaths in the input script for clarity
# Binary filepath definitions:
# 1. binary file for software TAPBS (computing protonation energy)
# replace with a full path
kbp2_set.tapbs_bin = 'LD_LIBRARY_PATH=../kb2plus_package/lib_fort_gnu ../kb2plus_package/bin/tapbs_1.3_cav_enere'
# 2. binary file for CHARMM
# replace with a full path
kbp2_set.charmm_bin = "<charmm binary filepath>"
# 3. binary file for software APBS (computing conformation energy)
# replace with a full path
kbp2_set.apbs_bin    = "../kb2plus_package/bin/apbs"
# 4. binary file for software TAPBS (computing electrostatic part of conformation energy)
# replace with a full path
kbp2_set.coulomb_bin = "../kb2plus_package/bin/coulomb"

# IMPORTANT NOTE: Binary for software Kalrsberg (Monte-Carlo algorithm) is hard coded in the module karlsberg.py
# IMPORTANT NOTE: Replace the line:
#         karlsberg_bin = "../kb2plus_package/bin/karlsberg2.x86_64_fixed"
#                 with a full filepath in karlsberg.py

###################################
##### TITRATABLE YAML CHOICE ######
###################################
# titratable yaml file containes definition of states of resdiues that are titrated,
# defined by patches and affected atomic partial charges
# in the folder aditional_files within kbp2 python package yaml files are kept
titratable_yaml = kbp2_set.titratable_yaml
kbp2_set.set_yaml(titratable_yaml)
titratable_definitions = kbp2.kbp_tools.parse_titratable_yaml(titratable_yaml)

#######################
### RUN KARLSBERG2+ ###
#######################
# final command, a function calc_pkas that will conduct the titration of a.a. residues for the given crystal structure
kbp2.pka_calculation.calc_pkas(kbp2_set, titratable_definitions)
