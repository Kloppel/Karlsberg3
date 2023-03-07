# -*- coding: utf-8 -*-

import kbp2


project_name = '4pti'
protein = '/scratch/scratch/pdb/pdb_bio_merged/pt/4pti.pdb1'
exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/4pti.dat"

project_name = '1a2p_cm'
protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/1a2p_cm.pdb'
exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/barnase.dat"

project_name = '2lzt'
protein = '/scratch/scratch/pdb/pdb_bio_merged/lz/2lzt.pdb1'
exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/lysozyme.dat"

project_name = '3rn3'
protein = '/scratch/scratch/pdb/pdb_bio_merged/rn/3rn3.pdb1'
exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/rnasaa.dat"

project_name = '2rn2'
protein = '/scratch/scratch/pdb/pdb_bio_merged/rn/2rn2.pdb1'
exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/rnaseh.dat"

project_name = '1hng_99c2_a'
exp_pkas = {"ASP-2_A" : 3.50,
            "ASP-25_A" : 3.53,
            "ASP-26_A" : 3.58,
            "ASP-28_A" : 3.57,
            "ASP-62_A" : 4.18,
            "ASP-71_A" : 3.20,
            "ASP-72_A" : 4.14,
            "ASP-94_A" : 3.83,
            "CTE-99_A" : 3.11,
            "GLU-29_A" : 4.42,
            "GLU-33_A" : 4.20,
            "GLU-41_A" : 6.53,
            "GLU-56_A" : 3.95,
            "GLU-99_A" : 4.10}

project_name = '3icb'
protein = '/scratch/scratch/pdb/pdb_bio_merged/ic/3icb.pdb1'
exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/cabd.dat"

project_name = '1ppf_i'
# 1ppf is a complex. Only chain I is "Third domain of the turkey ovomucoid inhibitor (OMTKY3)"
protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/1ppf_i.pdb1'
exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/omtky3_mod.dat"

project_name = '1ert_a'
protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/1ert_a_min7.pdb'
exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/thioredoxin_tmod.dat"

# project_name = '1xnb'
# protein = '/scratch/scratch/pdb/pdb_bio_merged/xn/1xnb.pdb1'
# exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/extra_part_1xnb.dat"
# ASN O <-> N swapped
project_name = '1xnb_35sw'
protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/1xnb_35sw.pdb'
exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/extra_part_1xnb.dat"

project_name = '1pga'
protein = '/scratch/scratch/pdb/pdb_bio_merged/pg/1pga.pdb1'
exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/strepprotg.dat"

project_name = '2zta_cm3'
# resids shifted by one!
# Is a homo dimer, so all pKas also exist for chain B.
protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/2zta_cm3.pdb'
exp_pkas = {"NTE-1_A" : 7.88,
            "LYS-4_A" : 10.78,
            "GLU-7_A" :  4.6,
            "ASP-8_A" :  3.48,
            "LYS-9_A" : 11.26,
            "GLU-11_A" : 3.94,
            "GLU-12_A" :  4.05,
            "LYS-16_A" : 10.62,
            "TYR-18_A" :  9.82,
            "HIS-19_A" :  6.24,
            "GLU-21_A" :  4.38,
            "GLU-23_A" :  4.2,
            "LYS-28_A" : 11.1,
            "LYS-29_A" : 10.64,
            "GLU-33_A" :  4.62,
            "CTE-33_A" :  4.03}
exp_pkas.update( {"NTE-1_B" : 7.88,
            "LYS-4_B" : 10.78,
            "GLU-7_B" :  4.6,
            "ASP-8_B" :  3.48,
            "LYS-9_B" : 11.26,
            "GLU-11_B" : 3.94,
            "GLU-12_B" :  4.05,
            "LYS-16_B" : 10.62,
            "TYR-18_B" :  9.82,
            "HIS-19_B" :  6.24,
            "GLU-21_B" :  4.38,
            "GLU-23_B" :  4.2,
            "LYS-28_B" : 11.1,
            "LYS-29_B" : 10.64,
            "GLU-33_B" :  4.62,
            "CTE-33_B" :  4.03} )

project_name = '3bdc'
protein = '/scratch/scratch/pdb/pdb_bio_merged/bd/3bdc.pdb1'
exp_pkas = {"ASP-95_A" :  2.16,
            "ASP-19_A" :  2.21,
            "GLU-101_A" : 3.81,
            "GLU-129_A" : 3.75,
            "GLU-122_A" : 3.89,
            "GLU-67_A" :  3.76,
            "GLU-43_A" :  4.32,
            "GLU-135_A" : 3.76,
            "GLU-57_A" :  3.49,
            "ASP-21_A" :  6.54,
            "GLU-52_A" :  3.93,
            "ASP-40_A" :  3.87,
            "GLU-10_A" :  2.82,
            "GLU-75_A" :  3.26,
            "GLU-73_A" :  3.31}
            # ASP	77	<2.2
            # ASP	83	<2.2




##################################
### Get old Karlsberg+ results ###
##################################
kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_po1_ie80/' + project_name + '/'

kbp_results = kbp2.kbp_results.KbpResult()
kbp_results.read_kbp_results(kbp_job_folder)
old_pkas = kbp_results.pkas


#############################
### Get experimental pKas ###
#############################
if type(exp_pkas) == str:
    exp_pkas = kbp2.kbp_tools.parse_gernots_exp_pkas(exp_pkas)

