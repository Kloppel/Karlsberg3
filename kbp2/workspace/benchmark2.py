# -*- coding: utf-8 -*-

import cPickle as pickle
from collections import defaultdict
import os
import sys

from kbp2.workspace.tools import *
import kbp2.analyse_md_pkas
import kbp2.workspace.tools

from kbp2.cython.analyse_md_pkas_tools import analyse_tension

if __name__ == '__main__':

    ##################
    ### Parameters ###
    ##################
    pick_pdb = None
    quite = False
    write_suggested_prot_states = False
    skip_analyse_tension = True


# Specifically, the van der Waals radii used are as follows:
# tetrahedral carbon, tetrahedral nitrogen or sulfur with hydrogens attached, 2.0 A;
# trigonal carbon and trigonal NH, 1.7 A;
# trigonal CH, CH2, and sulfur, 1.85 A;
# hydroxyl and trigonal nitrogen, 1.5 A;
# trigonal NH2, 1.8 A;
# oxygen and water, 1.4 A (Rashin, 1984).


    max_iteration = 0

    # pick_pdb = '3bdc'

    selected_residue = ''
    # selected_residue = 'LYS-1_A'


    # ! std choices !
    # kbp_subfolders = ["kbp_noMin_refzero", "kbp_noMin_refzero_cav09", "kbp_noMin_refzero_cav08",
    #                   "kbp_noMin_refzero_cav07", "kbp_noMin"]




    # charmm_confe_subfolder = None
    # charmm_confe_subfolder = 'charmm_confE'
    charmm_confe_subfolder = 'charmm_confE_no14'

    # kbp_subfolders = ["kbp_noMin_refzero", "kbp_e4_refzero"]
    # kbp_subfolders = ["kbp_noMin", "kbp_noMin_all", "kbp_e4"]
    # kbp_subfolders = ["kbp_noMin", "kbp_e4"]
    # kbp_subfolders = ["kbp_noMin", "kbp_e4", "kbp_noMin_refzero_ter2", "kbp_e4_refzero_ter2"]
    # kbp_subfolders = ["kbp_noMin", "kbp_e1","kbp_e4", "kbp_noMin_refzero_ter2", "kbp_e1_refzero_ter2", "kbp_e4_refzero_ter2"]
    # kbp_subfolders = ["kbp_noMin_refzero_ter2", "kbp_e1_refzero_ter2", "kbp_e4_refzero_ter2"]
    # kbp_subfolders = ["kbp_noMin", "kbp_e1","kbp_e4", "kbp_noMin_refzero_ter2_np", "kbp_e4_refzero_ter2_np"]

    # kbp_subfolders = ["kbp_noMin_refzero_ter2_np", "kbp_noMin_refzero_rashin2_ter2_np", "kbp_e1_refzero_ter2_np", "kbp_e4_refzero_ter2_np"]
    # kbp_subfolders = ["kbp_noMin_refzero_ter2_np", "kbp_noMin_refzero_rashin2_ter2_np", "kbp_e1_refzero_ter2_np", "kbp_e4_refzero_ter2_np", "kbp_e4_refzero_pe8"]
    # kbp_subfolders = ["kbp_e4_refzero_ter2_np" , "kbp_e4_refzero_cav09"]
    kbp_subfolders = ["kbp_e4_refzero_ter2_np"]

    # kbp_subfolders = ["kbp_e4_refzero_ter2_np"]
    # kbp_subfolders = ["kbp_e4_refzero_rashin2_ter2_np"]
    # kbp_subfolders = ["kbp_e4_refzero_ter2_np", "kbp_e4_refzero_ter2_np_propka"]
    # kbp_subfolders = ["kbp_e4_refzero_ter2_np", "kbp_e4_refzero_ter2_np_fit"]

    # kbp_subfolders = ["kbp_e4_refzero_pe2"]


    # SNase prot_tasks:
    # kbp_subfolders = ["kbp_noMin", "kbp_noMin_cav09", "kbp_noMin_refzero"]
    # kbp_subfolders = ["kbp_noMin", "kbp_noMin_cav09", "kbp_noMin_refzero"]
    # kbp_subfolders = ["kbp_noMin", "kbp_noMin_cav09", "kbp_e4_refzero_ter2_np"]
    # SNase std_prot_tasks:
    # kbp_subfolders = ["kbp_noMin_refzero", "kbp_noMin_refzero_cav09"]
    # kbp_subfolders = ["kbp_e4_refzero_ter2_np", "kbp_e4_refzero_cav09", "kbp_e4_refzero_cav07"]
    # kbp_subfolders = ["kbp_e4_refzero_pe8"]

    # kbp_subfolders = ["kbp_e4_refzero_ter2_np", "kbp_e4_refzero_cav09", "kbp_e4_refzero_cav07", "kbp_e4_refzero_pe8"]
    # kbp_subfolders = ["kbp_e4_refzero_ter2_np", "kbp_e4_refzero_cav09", "kbp_e4_refzero_cav07", "kbp_e4_refzero_cav07_w2", "kbp_e4_refzero_pe8"]

    # kbp_subfolders = ["kbp_e4_refzero_ter2_np", "kbp_e4_refzero_cav09", "kbp_e4_refzero_cav07", "kbp_e4_refzero_cav07_w2"]
    # kbp_subfolders = ["kbp_e4_refzero_cav07_w2"]

    # kbp_subfolders = ["kbp_e4_refzero_cav09_w2"]
    # kbp_subfolders = ["kbp_e4_refzero_cav09"]
    # kbp_subfolders = ["kbp_e4_refzero_cav09_w2_pe8"]
    # kbp_subfolders = ["kbp_e4_refzero_cav09_w2", "kbp_e4_refzero_cav09_w2_pe8"]
    # kbp_subfolders = ["kbp_e4_refzero_ter2_np"]



    # kbp_subfolders = ["kbp_noMin_refzero_ter2_np"]

    # propka_folder = None
    propka_folder = "/scratch/scratch/tmeyer/md_pka/propka_runs/"

    benchmark_set = 'std'
    # benchmark_set = 'snase'

    protonation_scheme = 'prot_tasks'
    # protonation_scheme = 'std_prot_tasks'

    # import time
    # time.sleep(3600)

    sasa_factor = 0.0

    print("Maximum iteration: %i" % max_iteration)
    print("protonation scheme: %s" % protonation_scheme)
    print("CHARMM confE subfolder: %s" % charmm_confe_subfolder)
    # print("Karlsberg+ folder: %s\n" % kbp_subfolder)

    project_folders = []

    if benchmark_set == 'std':
        main_project_folder = '/scratch/scratch/tmeyer/md_pka/runs/general2/'

        project_name = '4pti'
        exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/4pti.dat"
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '1pga'
        exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/strepprotg.dat"
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '1a2p_cm'
        exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/barnase.dat"
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '2lzt'
        exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/lysozyme.dat"
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '3rn3'
        exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/rnasaa.dat"
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '2rn2'
        exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/rnaseh.dat"
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        # project_name = '1hng_cm'
        # project_name = '1hng_99c_a'
        project_name = '1hng_99c2_a'
        # exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/cd2.dat"
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
        # exp_pkas.update( {"ASP-2_B" : 3.50,
        #             "ASP-25_B" : 3.53,
        #             "ASP-26_B" : 3.58,
        #             "ASP-28_B" : 3.57,
        #             "ASP-62_B" : 4.18,
        #             "ASP-71_B" : 3.20,
        #             "ASP-72_B" : 4.14,
        #             "ASP-94_B" : 3.83,
        #             "CTE-99_B" : 3.11,
        #             "GLU-29_B" : 4.42,
        #             "GLU-33_B" : 4.20,
        #             "GLU-41_B" : 6.53,
        #             "GLU-56_B" : 3.95,
        #             "GLU-99_B" : 4.10})
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '3icb'
        exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/cabd.dat"
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '1ppf_i'
        exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/omtky3_mod.dat"
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '1ert_a'
        exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/thioredoxin_tmod.dat"
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        # project_name = '1xnb'
        # exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/extra_part_1xnb.dat"
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        project_name = '1xnb_35sw'
        exp_pkas = "/user/tmeyer/workspace/misc/benchmark_pkas/fromGernot/extra_part_1xnb.dat"
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '2zta_cm3'
        # shifted by one (+1)!
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
        # exp_pkas = {}
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
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        # project_name = '2zta_cm3_l'
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '3bdc'
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
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))




    if benchmark_set == 'snase':
        main_project_folder = '/scratch/scratch/tmeyer/md_pka/runs/snase_nocal2/'

        project_name = '3evq'
        exp_pkas = {'GLU-25_A' : 7.5}
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '3erq'
        exp_pkas = {'LYS-25_A' : 6.3}
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '3eji_m'
        exp_pkas = {'LYS-36_A' : 7.2}
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '3d6c'
        exp_pkas = {'GLU-38_A' : 7.2}
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '3dmu_m'
        exp_pkas = {'LYS-62_A' : 8.1}
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '2oxp_m'
        exp_pkas = {'ASP-66_A' : 8.7}
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '3ero'
        exp_pkas = {'GLU-72_A' : 7.3}
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '2rbm_m'
        exp_pkas = {'LYS-72_A' : 8.6}
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '3c1f'
        exp_pkas = {'LYS-104_A' : 7.7}
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        project_name = '3c1e'
        exp_pkas = {'LYS-125_A' : 6.2}
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))




        # # experiments say: global unfolding
        # project_name = '2oeo'
        # exp_pkas = {'ASP-92_A' : 8.1}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        # # experiments say: global unfolding
        # project_name = '1tt2'
        # exp_pkas = {'LYS-92_A' : 5.6}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        # # >= 10.4
        # project_name = '2rks_m'
        # exp_pkas = {'LYS-38_A' : 10.4}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        # project_name = '1u9r'
        # exp_pkas = {'GLU-66_A' : 8.5}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        # # project_name = '1tqo'
        # # exp_pkas = {'GLU-92_A' : 8.7}
        # # project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        # # # Dimer
        # project_name = '3d4d'
        # exp_pkas = {'GLU-91_A' : 7.1}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        # # Dimer
        # project_name = '3e5s'
        # exp_pkas = {'LYS-103_A' : 8.2 }
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))


    old_pkas_folders = []
    old_pkas_names = []
    old_pkas_folders.append('/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_po1_ie80/')
    old_pkas_names.append('KB+')
    # old_pkas_folders.append('/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_po1_refzero_ter2_ie80/')
    # old_pkas_names.append('KB+_ref0')


    nr_of_res = 0
    avg_def = 0

    buried = []
    project_results = []
    for project_folder, exp_pkas in project_folders:
        md_kbp2 = kbp2.jobs.MD_kbp2_project(None, None)
        pickle_filename = project_folder + 'md_kbp2.pickle'
        md_kbp2.restore_from_pickle(pickle_filename)


        ### To generate TCL code for RMSDs in VMD ###
        # structure_filename = md_kbp2.charmm_out_prefix
        # pdb_filename = md_kbp2.tasks[0].task_folder + 'modelling/' + structure_filename + '.pdb'
        # if not os.path.exists(pdb_filename):
        #     print '# File does not exist: ' + pdb_filename
        # # else:
        # #     print pdb_filename
        # print "mol new {%s} type {pdb} first 0 last -1 step 1 waitfor 1" % pdb_filename
        # for task in md_kbp2.tasks:
        #
        #     psf_filename = task.task_folder + 'md/' + 'protein_in_water.xplor.psf'
        #     dcd_filename = task.task_folder + 'md/' + structure_filename + '_' + task.taskname + '_md.dcd'
        #     if not os.path.exists(psf_filename):
        #         print '# File does not exist: ' + psf_filename
        #     # else:
        #     #     print psf_filename
        #     if not os.path.exists(dcd_filename):
        #         print '# File does not exist: ' + dcd_filename
        #     # else:
        #     #     print dcd_filename
        #
        #     print "set molid [mol new {%s} type {dcd} first 0 last -1 step 2 waitfor 1]" % dcd_filename
        #     print "mol addfile {%s} type {psf} first 0 last -1 step 1 waitfor 1 $molid" % psf_filename
        #
        #     # print "set molid [mol new {%s} type {psf} first 0 last -1 step 1 waitfor 1]" % psf_filename
        #     # print "mol addfile {%s} type {dcd} first 0 last -1 step 2 waitfor 1 $molid" % dcd_filename
        #     #
        # print
        #############################################




        # old_pkas = md_kbp2.kbp_result.pkas
        jobname = project_folder.split('/')[-2]
        old_pkas = []
        for old_pkas_folder in old_pkas_folders:
            # print "Karlsberg+ results taken from: " + kbp_job_folder
            kbp_job_folder = old_pkas_folder + jobname + '/'
            kbp_results = kbp2.kbp_results.KbpResult()
            kbp_results.read_kbp_results(kbp_job_folder)
            old_pkas.append(kbp_results.pkas)
        # old_pkas = None


        project_name = project_folder.strip('/').split('/')[-1]
        propka_pkas_file = propka_folder + project_name + '.pka'
        propka_pkas = kbp2.workspace.tools.parse_propka_pkas(propka_pkas_file)
        old_pkas.append(propka_pkas)

        # project_name = project_folder.strip('/').split('/')[-1]
        project_result_from_kbp_subfolders = []
        for kbp_subfolder in kbp_subfolders:
            print("\nKarlsberg+ folder: %s" % kbp_subfolder)

            project_result = analyse_md_pkas.ProjectResults(project_name)
            project_result.set_exp_pkas(exp_pkas)
            if old_pkas is not None:
                # project_result.set_old_pkas(zero_pkas)
                project_result.set_old_pkas(old_pkas)
            project_result.read_results(project_folder, kbp_subfolder, sasa_factor, quite, max_iteration,
                                        charmm_confe_subfolder=charmm_confe_subfolder,
                                        protonation_scheme=protonation_scheme)
            project_result_from_kbp_subfolders.append(project_result)


        # Show old pKas
        # project_result.detailed_project_results.combined_results = kbp_results
        # project_result.pkas = old_pkas
        # print "printing Karlsberg+ results"



        ### To get sasa of the crystal structure ###
        temp_workdir = '/tmp/charmm_manager1/'
        if os.path.exists(temp_workdir):
            import shutil
            shutil.rmtree(temp_workdir)
        os.mkdir(temp_workdir)

        cs = md_kbp2.charmm_struct
        cs.workdir = temp_workdir
        cs.add_decision('rename__HOH_TIP3', 'keep')
        cs.check_structures(quiet=True)
        cs.run_charmm()
        sasa_structure = cs.get_modelled_structure()
        for par_filename in cs.par:
            sasa_structure.read_par(par_filename)
        # s = cs.structure
        sasa_structure.sasa(estimate=False)

        ### Filter exp pKas according to sasa ###
        exp_pkas_dict = project_result.exp_pkas
        filtered_exp_pkas_dict = {}
        buried.append({})
        for residue, pKa in exp_pkas_dict.iteritems():
            resname, resid, segname = re.split(r'[-_]', residue)
            resid = int(resid)

            if resname not in ['CTE', 'NTE']:
                buried_sasa = sasa_structure[segname][resid].side_chain_sasa

                s = sasa_structure.copy(residuelist=[(segname, resid)])
                s.sasa(estimate=False)
                total_sasa = s[segname][resid].side_chain_sasa
            else:
                if resname == 'CTE':
                    ter_atoms = ['C', 'OT1', 'OT2', 'OXT']
                elif resname == 'NTE':
                    ter_atoms = ['N', 'HT1', 'HT2', 'HT3']

                buried_sasa = 0
                for atom in sasa_structure[segname][resid].iter_atoms():
                    if atom['name'] in ter_atoms:
                        buried_sasa += atom['sasa']

                s = sasa_structure.copy(residuelist=[(segname, resid)])
                s.sasa(estimate=False)
                total_sasa = 0
                for atom in s[segname][resid].iter_atoms():
                    if atom['name'] in ter_atoms:
                        total_sasa += atom['sasa']

            buried_fraction = (1.0 - (buried_sasa / total_sasa)) * 100.0
            buried[-1][residue] = buried_fraction


            # if buried_fraction > 80:
            #     print '%s : %.0f %%' % (residue, buried_fraction)

            # Comment in for filtering
        #     if buried_fraction <= 80:
        #         filtered_exp_pkas_dict[residue] = pKa
        # project_result.exp_pkas = filtered_exp_pkas_dict


        project_results.append(project_result_from_kbp_subfolders)




        ### Print avg shifts of experimental pKas
    #     mod_pkas = {'ASP' : 4.0,
    #                 'GLU' : 4.4,
    #                 'HIS' : 6.8,
    #                 'LYS' : 10.4,
    #                 'TYR' : 9.6,
    #                 'NTE' : 7.5,
    #                 'CTE' : 3.8,}
    #
    #     for residue, pka in filtered_exp_pkas_dict.iteritems():
    #         nr_of_res += 1
    #         mod_pka = mod_pkas[re.split(r'[-_]', residue)[0]]
    #
    #         avg_def += abs(pka - mod_pka)
    # print avg_def/nr_of_res
    # import sys
    # sys.exit()



    # lines = []
    # line = 'calc,'
    # for kbp_subfolder in kbp_subfolders:
    #     line += kbp_subfolder + ','
    # line.strip(',')
    # lines.append(line)
    # for project_result_from_kbp_subfolders in project_results:
    #     project_name = project_result_from_kbp_subfolders[0].name
    #     line = '%20s: ' % project_name
    #     for project_result, kbp_subfolder in zip(project_result_from_kbp_subfolders, kbp_subfolders):
    #         pkas = project_result.pkas
    #         exp_pkas = project_result.exp_pkas
    #         total_rmsd, rmsds, diffs = calc_rmsd(pkas, exp_pkas)
    #         rmsd = rmsds[0]
    #         diffs = diffs[0]
    #
    #         # sorted_residue_list = kbp2.tools.get_sorted_residue_list(diffs.keys())
    #         line += '%5.2f ' % rmsd
    #
    #         ouput_str=[]
    #         print_project_rmsds(project_results, verbose=False, plot=False, ouput_str=ouput_str)
    #
    #     lines.append(line)
    #
    # print
    # for line in lines:
    #     print line

    if project_folder is not None:
        old_pkas_names.append('PROPKA')


    lines = []
    line = ' calc,count'

    for old_pkas_name in old_pkas_names:
        line += ',' + old_pkas_name
    for kbp_subfolder in kbp_subfolders:
        line += ',' + kbp_subfolder
    line.strip(',')
    lines.append(line)
    project_results_by_kbp_subfolder = defaultdict(list)
    for project_result_from_kbp_subfolders in project_results:
        old_pkas_printed = False
        count_printed = False
        project_name = project_result_from_kbp_subfolders[0].name
        line = '%12s' % project_name
        for project_result, kbp_subfolder in zip(project_result_from_kbp_subfolders, kbp_subfolders):
            # pkas = project_result.pkas
            # exp_pkas = project_result.exp_pkas
            # total_rmsd, rmsds, diffs = calc_rmsd(pkas, exp_pkas)
            # rmsd = rmsds[0]
            # diffs = diffs[0]

            # sorted_residue_list = kbp2.tools.get_sorted_residue_list(diffs.keys())
            # line += '%5.2f ' % rmsd

            if not old_pkas_printed:
                old_pkas_printed = True
                old_pkas_list = project_result.old_pkas
                for old_pkas in old_pkas_list:
                    project_result.old_pkas = old_pkas
                    output_str = []
                    max_def_str_list = []
                    total_rmsd, rmsds, diffs, total_rmsd_old = print_project_rmsds([project_result], verbose=False,
                                                                                   plot=False, output_str=output_str)
                    project_name, nr_of_residues_in_project, total_rmsd_str, total_rmsd_old_str = output_str[0]

                    if not count_printed:
                        count_printed = True
                        line += ' (%2s): ' % nr_of_residues_in_project
                    line += '%11s ' % total_rmsd_old_str
                project_result.old_pkas = old_pkas_list


            output_str = []
            total_rmsd, rmsds, diffs, total_rmsd_old = print_project_rmsds([project_result], verbose=False, plot=False,
                                                           output_str=output_str)

            project_name, nr_of_residues_in_project, total_rmsd_str, total_rmsd_old_str = output_str[0]
            line += '%11s ' % total_rmsd_str

            project_results_by_kbp_subfolder[kbp_subfolder].append(project_result)

        lines.append(line)

    print
    for line in lines:
        print line

    line = "%12s" %  'all'
    max_dev_line = "%12s" %  'max_dev, '
    old_pkas_printed = False
    count_printed = False
    for i, kbp_subfolder in enumerate(kbp_subfolders):
        project_results_for_kbp_subfolder = project_results_by_kbp_subfolder[kbp_subfolder]

        if not old_pkas_printed:
            old_pkas_printed = True
            for i in range(len(old_pkas_names)):
                for project_result in project_results_for_kbp_subfolder:
                    project_result.old_pkas_list = project_result.old_pkas
                    project_result.old_pkas = project_result.old_pkas_list[i]

                max_def_str_list = []
                total_rmsd, rmsds, diffs, total_rmsd_old = print_project_rmsds(project_results_for_kbp_subfolder,
                                                                               verbose=False, plot=False,
                                                                               max_def_str_list=max_def_str_list)

                if not count_printed:
                    count_printed = True
                    nr_of_residues = 0
                    for diff in diffs:
                        nr_of_residues += len(diff)
                    line += ' (%2s): ' % nr_of_residues

                max_def_str, max_def_str_old = max_def_str_list
                max_dev_line += ',' + max_def_str_old

                line += '%11.2f ' % total_rmsd_old

                for project_result in project_results_for_kbp_subfolder:
                    project_result.old_pkas = project_result.old_pkas_list


        project_result.old_pkas = old_pkas_list

        max_def_str_list = []
        total_rmsd, rmsds, diffs, total_rmsd_old = print_project_rmsds(project_results_for_kbp_subfolder, verbose=False,
                                                                       plot=False, max_def_str_list=max_def_str_list)
        max_def_str, max_def_str_old = max_def_str_list
        max_dev_line += ',' + max_def_str

        line += '%11.2f ' % total_rmsd
    print line
    print max_dev_line

    ##################################################################################

    lines = []
    # line = ' calc,'
    #
    # for old_pkas_name in old_pkas_names:
    #     line += old_pkas_name + ','
    # for kbp_subfolder in kbp_subfolders:
    #     line += kbp_subfolder + ','
    # line.strip(',')
    # lines.append(line)

    for resname in ['ASP', 'GLU', 'LYS', 'HIS', 'TYR', 'NTE', 'CTE']:
        resname_filter = [resname]
        line = '      %s' % resname

        count_printed = False
        old_pkas_printed = False
        for i, kbp_subfolder in enumerate(kbp_subfolders):
            project_results_for_kbp_subfolder = project_results_by_kbp_subfolder[kbp_subfolder]

            if not old_pkas_printed:
                old_pkas_printed = True
                for i in range(len(old_pkas_names)):
                    for project_result in project_results_for_kbp_subfolder:
                        project_result.old_pkas_list = project_result.old_pkas
                        project_result.old_pkas = project_result.old_pkas_list[i]

                    total_rmsd, rmsds, diffs, total_rmsd_old = print_project_rmsds(project_results_for_kbp_subfolder,
                                                                                   verbose=False, plot=False,
                                                                                   resname_filter=resname_filter)
                    if not count_printed:
                        count_printed = True
                        nr_of_residues = 0
                        for diff in diffs:
                            nr_of_residues += len(diff)
                        line += ' (%2s): ' % nr_of_residues

                    line += '%11.2f ' % total_rmsd_old

                    for project_result in project_results_for_kbp_subfolder:
                        project_result.old_pkas = project_result.old_pkas_list


            project_result.old_pkas = old_pkas_list

            total_rmsd, rmsds, diffs, total_rmsd_old = print_project_rmsds(project_results_for_kbp_subfolder,
                                                                           verbose=False, plot=False,
                                                                           resname_filter=resname_filter)
            line += '%11.2f ' % total_rmsd
        lines.append(line)
    print
    for line in lines:
        print line









    ##################################################################################

    lines = []
    line = '%20s,' % ' calc'
    for kbp_subfolder in kbp_subfolders:
        line +=  '%20s,' % kbp_subfolder
    line.strip(',')
    lines.append(line)
    for project_no, project_result_from_kbp_subfolders in enumerate(project_results):
        project_name = project_result_from_kbp_subfolders[0].name

        sorted_residue_list = None
        residue_rmsds = []
        for project_result, kbp_subfolder in zip(project_result_from_kbp_subfolders, kbp_subfolders):
            pkas = project_result.pkas
            exp_pkas = project_result.exp_pkas
            total_rmsd, rmsds, diffs = calc_rmsd(pkas, exp_pkas)
            rmsd = rmsds[0]
            diffs = diffs[0]

            sorted_residue_list = kbp2.workspace.tools.get_sorted_residue_list(diffs.keys())
            residue_rmsds_in_kbp_subfolder = []
            for residue_descr in sorted_residue_list:
                residue_descr_kbp = kbp_tools.get_kbp_residue(residue_descr)
                residue_rmsds_in_kbp_subfolder.append(pkas[residue_descr_kbp])

            residue_rmsds.append(residue_rmsds_in_kbp_subfolder)


        for i, old_pkas_name in enumerate(old_pkas_names):
            old_pkas = project_result_from_kbp_subfolders[0].old_pkas[i]

            kbp_residue_rmsds_in_kbp_subfolder = []
            for residue_descr in sorted_residue_list:
                residue_descr_kbp = kbp_tools.get_kbp_residue(residue_descr)
                kbp_residue_rmsds_in_kbp_subfolder.append(old_pkas[residue_descr_kbp])
            residue_rmsds.insert(i, kbp_residue_rmsds_in_kbp_subfolder)


        residue_rmsds_in_kbp_subfolder = []
        for residue_descr in sorted_residue_list:
            # residue_descr_kbp = kbp_tools.get_kbp_residue(residue_descr)
            residue_rmsds_in_kbp_subfolder.append(exp_pkas[residue_descr])
        residue_rmsds.insert(0, residue_rmsds_in_kbp_subfolder)



        lines = []

        line = ' calc,'
        for kbp_subfolder in ["exp"] + old_pkas_names + kbp_subfolders + ["buried"]:
            line += kbp_subfolder + ','
        line.strip(',')

        lines.append(line)
        for i, residue_descr in enumerate(sorted_residue_list):
            resname, resid, segname = re.split(r'[-_]', residue_descr)
            resname = resname.replace('NTE', 'NTER').replace('CTE', 'CTER')
            resname = resname[0].upper() + resname[1:].lower()
            if project_name in ['2zta_cm3']:
                residue_str = resname + '%i' % int(resid) + ' (%s)' % segname
            else:
                residue_str = resname + '%i' % int(resid)
            line = "%10s: " % residue_str
            # line = "%10s: " % residue_descr
            for residue_rmsds_in_kbp_subfolder, kbp_subfolder in zip(residue_rmsds, ["exp"] + old_pkas_names + kbp_subfolders):
                line += "%6.2f," % residue_rmsds_in_kbp_subfolder[i]


            # print sorted_residue_list
            # print buried[project_no]
            line += "%3.0f" % buried[project_no][residue_descr]
            line.strip(',')
            lines.append(line)

        print "\n"
        print project_name + ":"
        for line in lines:
            print line
        print "\n"



    # All pkas
    # print_project_rmsds(project_results)


    # print("")
    # # Select one Project
    # selected_project = select_project(project_results, pick_pdb)
    #
    # if selected_project is not None:
    #     # Print pkas
    #     print_project_pkas(selected_project, all=True)
    #     # Plot MD weights
    #     plot_md_weight(selected_project)
    #
    #     if selected_residue:
    #         kbp2.workspace.tools.plot_titration_curve(selected_project, selected_residue)

