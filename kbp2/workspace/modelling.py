# -*- coding: utf-8 -*-

import kbp2
import os



# pdb_filename = '/scratch/scratch/tmeyer/md_pka/mod_bio/2ibx_m.pdb1'
# modelling_dir = '/scratch/scratch/tmeyer/md_pka/mod_bio/2ibx_m_charmm/'
# modelling_dir = '/scratch/scratch/tmeyer/md_pka/mod_bio/2ibx_m_charmm2/'
pdb_filename = '/scratch/scratch/pdb/pdb_bio_merged/zt/2zta.pdb1'
modelling_dir = '/scratch/scratch/tmeyer/md_pka/mod_bio/2zta_m_charmm3/'
# pdb_filename = '/scratch/scratch/tmeyer/md_pka/mod_bio/1hng.pdb1'
# modelling_dir = '/scratch/scratch/tmeyer/md_pka/mod_bio/1hng_charmm/'





if not os.path.exists(modelling_dir):
    os.mkdir(modelling_dir)

##########################################
### Read structure and model hydrogens ###
##########################################
s = kbp2.file_parser.Simple_struct_parser()
s.read_pdb(pdb_filename)

s = s.copy(residuelist=[('A', 0), ('B', 0)], exclude=True)
s.shift_segment('A', 1)
s.shift_segment('B', 1)


top = []
top.append('/scratch/scratch/tmeyer/karlsbergplus/top.inp')
top.append('/scratch/scratch/tmeyer/karlsbergplus/patches.rtf')
par = []
par.append('/scratch/scratch/tmeyer/karlsbergplus/par.inp')
par.append('/scratch/scratch/tmeyer/karlsbergplus/patches.prm')

c = kbp2.charmm.Charmm_manager(workdir=modelling_dir, top=top, par=par)
c.add_structure(s)

for segname in ['A', 'B']:
    first_resid = c.structure.struct[segname].resids()[0]
    first_residue = c.structure.struct[segname][first_resid]
    last_resid = c.structure.struct[segname].resids()[-1]
    last_residue = c.structure.struct[segname][last_resid]
    c.add_sequence(first_residue, 'GLY', is_last_residue=True)
    c.add_sequence(last_residue, 'GLU')

# 2zta
# The N terminus of the molecule is Gly0. This residue is an artifact of the engineered thrombin cleavage site used
# to purify the protein.26 The C terminus of the molecule is Glu32. The pKa value for Glu32 is for the side-chain.

# 2zta
# REMARK 465     GLU A    32
# REMARK 465     ARG A    33
# REMARK 465     GLU B    32
# REMARK 465     ARG B    33

# 2ibx
# REMARK 465     LEU A     1
# REMARK 465     VAL A     2
# REMARK 465     LYS A     3
# REMARK 465     SER A     4
# REMARK 465     GLN A   326
# REMARK 465     ARG A   327
# REMARK 465     GLU A   328

# c.add_decision('rename__HOH_TIP3', 'keep')
c.check_structures(quiet=False)
done = c.run_charmm()