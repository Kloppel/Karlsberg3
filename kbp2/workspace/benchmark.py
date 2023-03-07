# -*- coding: utf-8 -*-

import cPickle as pickle
import os
import sys

from kbp2.workspace.tools import *
import kbp2.analyse_md_pkas

from kbp2.cython.analyse_md_pkas_tools import analyse_tension

if __name__ == '__main__':

    ##################
    ### Parameters ###
    ##################
    pick_pdb = None
    quite = False
    write_suggested_prot_states = False
    skip_analyse_tension = True

    # pick_pdb = '3bdc'
    # pick_pdb = '2lzt'
    # pick_pdb = '2lzt_l'
    # pick_pdb = '1a2p_cm'

    # pick_pdb = '3rn3'

    # pick_pdb = '3icb'
    # pick_pdb = '4pti'
    # pick_pdb = '1xnb'
    # pick_pdb = '1xnb_35sw'
    # pick_pdb = '1ppf_i'
    # pick_pdb = '1pga'

    # pick_pdb = '2rn2'

    # pick_pdb = '1hng_cm'
    # pick_pdb = '1hng_99c2_a'
    # pick_pdb = '1ert_a'

    # pick_pdb = '2zta_cm3'


# Specifically, the van der Waals radii used are as follows:
# tetrahedral carbon, tetrahedral nitrogen or sulfur with hydrogens attached, 2.0 A;
# trigonal carbon and trigonal NH, 1.7 A;
# trigonal CH, CH2, and sulfur, 1.85 A;
# hydroxyl and trigonal nitrogen, 1.5 A;
# trigonal NH2, 1.8 A;
# oxygen and water, 1.4 A (Rashin, 1984).

    # SNASE
    # i1
    pick_pdb = '3d6c'
    # bug
    # pick_pdb = '1u9r'
    # i1
    # pick_pdb = '3c1e'
    # ~ i1
    # pick_pdb = '1tqo'
    # # i1 ~
    # pick_pdb = '2oeo'
    # # i1 ~
    # pick_pdb = '1tt2'
    # ist ok
    # pick_pdb = '2rks'
    # evtl # i1
    # pick_pdb = '3erq'
    # i1
    # pick_pdb = '3evq'


    # Bug
    # pick_pdb = '2oxp'
    # Bug
    # pick_pdb = '3dmu'
    # pick_pdb = '2rks'

    # pick_pdb = '3c1f'
    # pick_pdb = '3eji'
    #pick_pdb = '3d4d'
    # pick_pdb = '3ero'




    max_iteration = 0

    selected_residue = ''
    # selected_residue = 'LYS-1_A'
    # selected_residue = 'DPP-21_A'
    # selected_residue = 'DPP-19_A'
    # selected_residue = 'NTE-1_A'
    # selected_residue = 'EPP-73_A'
    # selected_residue = 'NTE-1_A'
    # selected_residue = 'DPP-66_A'
    # selected_residue = 'EPP-68_A'
    # selected_residue = 'EPP-172_A'

    # charmm_confe_subfolder = None
    # charmm_confe_subfolder = 'charmm_confE'
    charmm_confe_subfolder = 'charmm_confE_no14'


    # kbp_subfolder = "kbp_e4"
    kbp_subfolder = "kbp_e4_refzero_ter2_np"
    # kbp_subfolder = "kbp_e1_refzero_ter2_np"
    # kbp_subfolder = "kbp_noMin_refzero_ter2_np"
    # kbp_subfolder = "kbp_noMin_refzero_rashin2_ter2_np"

    # kbp_subfolder = "kbp_noMin_refzero_cav08"

    # kbp_subfolder = "kbp_e4_refzero_rashin2_ter2_np"





    # benchmark_set = 'std'
    benchmark_set = 'snase'

    # protonation_scheme = 'prot_tasks'
    protonation_scheme = 'std_prot_tasks'
    # protonation_scheme = 'prot_tasks_172p'
    # protonation_scheme = 'prot_tasks_j172p'

    # 4.184kJ/mol= 1.0kcal/mol.
    sasa_factor = 0.0

    print("Maximum iteration: %i" % max_iteration)
    print("protonation scheme: %s" % protonation_scheme)
    print("Karlsberg+ folder: %s\n" % kbp_subfolder)

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
        # shifted by one!
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

        project_name = '3d6c'
        exp_pkas = {'GLU-38_A' : 7.2}
        project_folders.append((main_project_folder + project_name + '/', exp_pkas))

        # # project_name = '1u9r'
        # # exp_pkas = {'GLU-66_A' : 8.5}
        # # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # project_name = '2oxp_m'
        # exp_pkas = {'ASP-66_A' : 8.7}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # project_name = '3c1e'
        # exp_pkas = {'LYS-125_A' : 6.2}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # # project_name = '1tqo'
        # # exp_pkas = {'GLU-92_A' : 8.7}
        # # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # project_name = '2oeo'
        # exp_pkas = {'ASP-92_A' : 8.1}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # project_name = '1tt2'
        # exp_pkas = {'LYS-92_A' : 5.6}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # project_name = '2rks_m'
        # exp_pkas = {'LYS-38_A' : 10.4}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # project_name = '3erq'
        # exp_pkas = {'LYS-25_A' : 6.3}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # project_name = '3evq'
        # exp_pkas = {'GLU-25_A' : 7.5}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # project_name = '3dmu_m'
        # exp_pkas = {'LYS-62_A' : 8.1}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # # # # Dimer
        # # project_name = '3d4d'
        # # exp_pkas = {'GLU-91_A' : 7.1}
        # # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # project_name = '2rbm_m'
        # exp_pkas = {'LYS-72_A' : 8.6}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # project_name = '3ero'
        # exp_pkas = {'GLU-72_A' : 7.3}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # # # Dimer
        # # project_name = '3e5s'
        # # exp_pkas = {'LYS-103_A' : 8.2 }
        # # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # project_name = '3c1f'
        # exp_pkas = {'LYS-104_A' : 7.7}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))
        #
        # project_name = '3eji_m'
        # exp_pkas = {'LYS-36_A' : 7.2}
        # project_folders.append((main_project_folder + project_name + '/', exp_pkas))


    project_results = []
    for project_folder, exp_pkas in project_folders:
        md_kbp2 = kbp2.jobs.MD_kbp2_project(None, None)
        pickle_filename = project_folder + 'md_kbp2.pickle'
        md_kbp2.restore_from_pickle(pickle_filename)

        # old_pkas = md_kbp2.kbp_result.pkas
        jobname = project_folder.split('/')[-2]
        # kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache/' + jobname + '/'
        # kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi/' + jobname + '/'
        # kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_po1/' + jobname + '/'
        # kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_po1_refzero/' + jobname + '/'
        kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_po1_ie80/' + jobname + '/'
        # kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_po1_ie80/' + jobname + '/'
        # kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_po1_refzero_ie80/' + jobname + '/'
        # kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_po1_refzero_nter2_ie80/' + jobname + '/'
        # kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_po1_refzero_ter2_ie80/' + jobname + '/'

        # kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_po3_refzero_ter2/' + jobname + '/'
        # kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_alltitr/' + jobname + '/'
        # kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_refzero_alltitr/' + jobname + '/'


        # kbp_job_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_multi_po1_rashin_ie80/' + jobname + '/'
        print "Karlsberg+ results taken from: " + kbp_job_folder

        kbp_results = kbp2.kbp_results.KbpResult()
        kbp_results.read_kbp_results(kbp_job_folder)
        old_pkas = kbp_results.pkas
        # **************
        # old_pkas = None
        # **************


        #
        # zero_pkas = old_pkas
        # for residue_descr in old_pkas:
        #     resname, resid, segname = re.split(r'[-_]', residue_descr)
        #     titr_residues = kbp_results.descr.titratable_residues
        #     model_pka = titr_residues[resname][1]['pka']
        #     if model_pka < -50:
        #         model_pka = abs(model_pka + 100)
        #     else:
        #         model_pka = abs(model_pka)
        #     zero_pkas[residue_descr] = model_pka


        project_name = project_folder.strip('/').split('/')[-1]
        project_result = analyse_md_pkas.ProjectResults(project_name)
        project_result.set_exp_pkas(exp_pkas)
        if old_pkas is not None:
            # project_result.set_old_pkas(zero_pkas)
            project_result.set_old_pkas(old_pkas)
        project_result.read_results(project_folder, kbp_subfolder, sasa_factor, quite, max_iteration,
                                    charmm_confe_subfolder=charmm_confe_subfolder,
                                    protonation_scheme=protonation_scheme)

        # Show old pKas
        # project_result.detailed_project_results.combined_results = kbp_results
        # project_result.pkas = old_pkas
        # print "printing Karlsberg+ results"

        project_results.append(project_result)


        # # Recheck the protonation pattern
        # # Read MD_kbp2_project object for the project
        # md_kbp2 = kbp2.jobs.MD_kbp2_project(None, None)
        # pickle_filename = project_folder + 'md_kbp2.pickle'
        # md_kbp2.restore_from_pickle(pickle_filename)
        #
        # print "Checking folder " + project_folder
        # folders = os.listdir(project_folder)
        # for folder in folders:
        #     prot_pickel_filename = project_folder + folder + '/prot_task.pickle'
        #     job_pickel_filename = project_folder + folder + '/job_dump.pickle'
        #     if os.path.exists(job_pickel_filename):
        #         # This seems to be a task folder
        #
        #         md_folder = project_folder + folder + '/md/'
        #         residue_list = kbp_results.descr.sorted_residue_list
        #         titratable_yaml = \
        #             kbp_tools.parse_titratable_yaml('/scratch/scratch/tmeyer/kbplus2/titratable.yaml')
        #
        #         prot_state_yaml = kbp_tools.determine_md_protonation_pattern(md_folder, residue_list, titratable_yaml)
        #
        #         titr_residue_dict_md = md_kbp2.charmm_struct.get_titr_residue_dict()
        #         md_kbp2.charmm_struct.mod_titr_residue_dict_from_titrable_yaml(titr_residue_dict_md, prot_state_yaml,
        #                                                                        titratable_yaml)
        #
        #         # Write the prot_state
        #         f = open(prot_pickel_filename, 'w')
        #         pickle.dump((folder, titr_residue_dict_md), f)
        #         f.close()


    # All pkas
    print_project_rmsds(project_results)




    print("")
    # Select one Project
    selected_project = select_project(project_results, pick_pdb)

    if selected_project is not None:
        # Print pkas
        print_project_pkas(selected_project, all=True)
        # Plot MD weights
        plot_md_weight(selected_project)

        if selected_residue:
            kbp2.workspace.tools.plot_titration_curve(selected_project, selected_residue)



    if skip_analyse_tension:
        plt.show(True)
        sys.exit()


    # solv_ref = selected_project.detailed_project_results.kbp_results[0].conf_energies['solvation']
    # coul_ref = selected_project.detailed_project_results.kbp_results[0].conf_energies['coulomb']
    # for md in selected_project.detailed_project_results.kbp_results:
    #     solv = selected_project.detailed_project_results.kbp_results[md].conf_energies['solvation']
    #     coul = selected_project.detailed_project_results.kbp_results[md].conf_energies['coulomb']
    #     # print solv - solv_ref
    #     print coul - coul_ref


    user_selected_project = selected_project
    for selected_project in project_results:
        print selected_project.path
        #######################
        ### Analyse tension ###
        #######################
        ### Read md_kbp2 object for chosen project.
        project_folder = selected_project.path

        # Find latest iteration
        files_in_project_folder = os.listdir(project_folder)

        reg = re.compile(r'^prot_tasks_?i?(\d*).pickle$')
        last_iteration = 0
        for file in files_in_project_folder:
            reg_m = reg.match(file)
            if reg_m is not None:
                # if len(reg_m.groups()) > 0:
                iteration = reg_m.groups()[0]
                if iteration:
                    iteration = int(iteration)
                    if iteration > last_iteration:
                        last_iteration = iteration
        if max_iteration != -1:
            if last_iteration > max_iteration:
                last_iteration = max_iteration
        next_iteration = last_iteration + 1

        # Read protonation tasks
        if last_iteration != 0:
            pickle_prot_tasks = 'prot_tasks_i%i.pickle' % last_iteration
            if selected_project is user_selected_project:
                print("\nReading protonation from %s" % pickle_prot_tasks)
        else:
            pickle_prot_tasks = 'prot_tasks.pickle'
            if selected_project is user_selected_project:
                print("\nReading protonation from %s" % pickle_prot_tasks)
        prot_tasks_filename = project_folder + pickle_prot_tasks
        f = open(prot_tasks_filename, 'r')
        task_names, titr_residue_dicts = pickle.load(f)
        f.close()

        # Read MD_kbp2_project object for the project
        md_kbp2 = kbp2.jobs.MD_kbp2_project(None, None)
        pickle_filename = selected_project.path + 'md_kbp2.pickle'
        md_kbp2.restore_from_pickle(pickle_filename)




        # Analyse tension
        print("")
        # sasa_factor = 0.0
        min_delta_e = -20.0
        tension, best_states, best_mds = analyse_tension(selected_project.detailed_project_results, sasa_factor, min_delta_e)
        ph_values = selected_project.detailed_project_results.descr.ph_values
        residue_list = selected_project.detailed_project_results.descr.sorted_residue_list
        # residue_list_ext = selected_project.detailed_project_results.descr.residue_list_ext

        residue_list_filtered = []
        tension_filtered = []
        for residue_tension, residue_descr in zip(tension, residue_list):
            if np.sum(residue_tension) < -0.1:
                tension_filtered.append(residue_tension)
                residue_list_filtered.append(residue_descr)

        if selected_project is user_selected_project:
            interactive_plot(ph_values, tension_filtered, residue_list_filtered, show=False)




        # Remove identical states:
        # GLU: 1 and 3 (MD default: 3)
        # ASP: 1 and 3 (MD default: 3)
        # Set these residues do -1:
        filtered_resnames = ['ARG']
        replace_states = {'EPP' : {1 : 3},
                          'DPP' : {1 : 3}}
        for residue_best_state, residue_descr in zip(best_states, residue_list):
            resname = re.split(r'[-_]', residue_descr)[0]
            if resname in replace_states:
                res_replace_states = replace_states[resname]
                for ph_number, ph in enumerate(ph_values):
                    state = residue_best_state[ph_number]
                    if state in res_replace_states:
                        residue_best_state[ph_number] = res_replace_states[state]
            elif resname in filtered_resnames:
                # Set all entries to -1
                residue_best_state -= residue_best_state
                residue_best_state += -1



        # Make suggestions for MD simulations.
        # prot_states = []
        # current_prot_state = None
        # for ph_number, ph in enumerate(ph_values):
        #     if current_prot_state is None:
        #         # The first ph value
        #         current_prot_state = best_states[:, ph_number]
        #         prot_states.append(current_prot_state)
        #     else:
        #         new_prot_state = best_states[:, ph_number]
        #         # print type(new_prot_state)
        #
        #         if not all(new_prot_state == current_prot_state):
        #             prot_states.append(new_prot_state)
        #         else:
        #             prot_states.append(current_prot_state)

        # Make suggestions for MD simulations. Version II
        titratable_yaml = selected_project.detailed_project_results.descr.titratable_residues
        prot_states = []
        for ph_number, ph in enumerate(ph_values):
            md_weights = selected_project.detailed_project_results.weights
            chosen_md = np.argmax(md_weights[ph_number])

            current_prot = titr_residue_dicts[chosen_md]

            better_prot_kbp = {}
            for residue_nr, residue_descr in enumerate(residue_list):
                state = best_states[residue_nr, ph_number]
                if state != -1:
                    if residue_descr == "":
                        print residue_descr
                    better_prot_kbp[residue_descr] = state

            # Copy the protonation pattern, so it can be modified
            better_prot = kbp2.charmm.Charmm_manager.copy_titr_residue_dict(current_prot)
            md_kbp2.charmm_struct.mod_titr_residue_dict_from_titrable_yaml(better_prot, better_prot_kbp, titratable_yaml)

            # Convert better_prot back into a list of protonations
            prot_list = []
            for residue_nr, residue_descr in enumerate(residue_list):
                resname_kbp, resid, segname = re.split(r'[-_]', residue_descr)
                resname = kbp2.kbp_tools.get_real_resname(resname_kbp)
                resid = int(resid)
                residue_tuple = (resname, resid, segname)

                state = better_prot[residue_tuple]
                if state is None:
                    state = 0

                prot_list.append(state)

            prot_list = np.array(prot_list, dtype=np.int)
            prot_states.append(prot_list)



        unique_prot_states = []
        start_ph_values = []
        # Minimum length in terms of ph steps a protonation state has to live to be chosen as suggestion.
        # 2 -> 1 ph unit
        min_ph = 0.0
        max_ph = 12.0
        length_required = 1
        current_length = 0
        current_prot_state = None
        start_ph = None
        for new_prot_state, ph in zip(prot_states, ph_values):
            current_length += 1
            if current_prot_state is None:
                current_prot_state = new_prot_state
                start_ph = ph
                continue

            if not all(new_prot_state == current_prot_state):
                if current_length >= length_required:
                    if (ph > min_ph) and (start_ph < max_ph):
                        unique_prot_states.append(current_prot_state)
                        start_ph_values.append(start_ph)
                current_length = 0
                start_ph = ph
                current_prot_state = new_prot_state
        if current_length >= length_required:
            if (ph > min_ph) and (start_ph < max_ph):
                unique_prot_states.append(current_prot_state)
                start_ph_values.append(start_ph)


        if selected_project is user_selected_project:
            text = "%8s" % ''
            for start_ph in start_ph_values:
                text += "%7.1f" % start_ph
            print text
            for residue_descr in sorted(residue_list):
                residue_nr = residue_list.index(residue_descr)
                text = "%-10s" % residue_descr
                for prot_state in unique_prot_states:
                    text += "%7i" % prot_state[residue_nr]
                print text



        ### Compare with the currently chosen protonation states ###

        # Translate suggested states into titr_residue_dict dictionaries
        charmm_manager = kbp2.charmm.Charmm_manager()
        titratable_yaml = selected_project.detailed_project_results.descr.titratable_residues
        titr_residue_dicts_new = []
        for prot_state in unique_prot_states:
            titr_residue_dict_new = {}
            for residue_nr, residue_descr in enumerate(residue_list):
                resname_kbp, resid, segname = re.split(r'[-_]', residue_descr)
                resname = kbp2.kbp_tools.get_real_resname(resname_kbp)
                resid = int(resid)
                residue_tuple = (resname, resid, segname)

                state = prot_state[residue_nr]
                titr_residue_dict_new[residue_tuple] = state

            titr_residue_dicts_new.append(titr_residue_dict_new)

            # prot_state_yaml = {}
            # for residue_nr, residue_descr in enumerate(residue_list):
            #     resname, resid, segname = re.split(r'[-_]', residue_descr)
            #     # if prot_state[residue_nr] > -1:
            #     prot_state_yaml[residue_descr] = prot_state[residue_nr]
            #
            #     # # Replace entries with -1 (no tension or filtered) with the state used in the original MD.
            #     # resid = int(resid)
            #     # residue_tuple = (resname, resid, segname)
            #     # old_state = best_mds[residue_descr]
            #     # md_kbp2.charmm_struct.mod_titr_residue_dict_single(titr_residue_dict_new, residue_tuple)
            #
            # titr_residue_dict_new = md_kbp2.charmm_struct.get_titr_residue_dict()
            # charmm_manager.mod_titr_residue_dict_from_titrable_yaml(titr_residue_dict_new, prot_state_yaml, titratable_yaml)
            # titr_residue_dicts_new.append(titr_residue_dict_new)


        # Current: titr_residue_dicts
        # Suggested: titr_residue_dicts_new
        accepted_titr_residue_dicts = []
        accepted_start_ph_values = []
        for i, titr_residue_dict_new in enumerate(titr_residue_dicts_new):
            if selected_project is user_selected_project:
                print("\nSuggestion %i starting at ph %.1f:" % (i, start_ph_values[i]))
            any_match = False
            for j, titr_residue_dict in enumerate(titr_residue_dicts):
                match = True
                mismatches = []
                for residue_descr_kb in residue_list:
                    residue_descr = kbp_tools.get_real_residue(residue_descr_kb)
                    resname, resid, segname = re.split(r'[-_]', residue_descr)
                    residue_tuple = (resname, int(resid), segname)

                    if resname in filtered_resnames:
                        continue

                    suggested_state = titr_residue_dict_new[residue_tuple]
                    current_state = titr_residue_dict[residue_tuple]
                    if current_state is None:
                        current_state = 0

                    if current_state != suggested_state:
                        match = False
                        diff_str = "   %s: %i -> %i" % (residue_descr, current_state, suggested_state)
                        mismatches.append(diff_str)

                # assert(not (match and any_match))
                if match:
                    if selected_project is user_selected_project:
                        print("   Great!! -> Matches protonation state %i:" % j)
                    any_match = True
                elif len(mismatches) < 5:
                    if selected_project is user_selected_project:
                        print("   -> similar to protonation state %i:" % j)
                        for mismatch in mismatches:
                            print(mismatch)

            if not any_match:
                if selected_project is user_selected_project:
                    print("   ### NEW ###")
                accepted_titr_residue_dicts.append(titr_residue_dict_new)
                accepted_start_ph_values.append(start_ph_values[i])


        # generate name suggestions for the new tasks
        nr_of_new_states = len(accepted_titr_residue_dicts)
        next_tast_names = []
        for i in range(nr_of_new_states):
            ph_str = str(accepted_start_ph_values[i])
            ph_str = ph_str.replace('.', '_')
            suffix = ''
            c = 0
            while True:
                task_folder_name = 'i%i_ph%s' % (next_iteration, ph_str + suffix)
                if not os.path.exists(project_folder + task_folder_name):
                    break
                c += 1
                suffix = '__' + str(c)
            next_tast_names.append(task_folder_name)

        # Search the project folder for exiting tasks and change the name of the new tasks, if they have been calculated
        # already.
        new_tasks_already_calculated = 0
        folders = os.listdir(project_folder)
        for folder in folders:
            prot_pickel_filename = project_folder + folder + '/prot_task.pickle'
            if os.path.exists(prot_pickel_filename):
                f = open(prot_pickel_filename, 'r')
                name, titr_residue_dict = pickle.load(f)
                f.close()

                # Check if this task matches one of the new tasks
                residue_tuples = sorted(titr_residue_dict.keys())
                for i, accepted_titr_residue_dict in enumerate(accepted_titr_residue_dicts):
                    accepted_residue_tuples = sorted(accepted_titr_residue_dict.keys())
                    if residue_tuples != accepted_residue_tuples:
                        continue
                    for residue_tuple, state_ref in accepted_titr_residue_dict.iteritems():
                        state = titr_residue_dict[residue_tuple]
                        if state is None:
                            state = 0
                        if state != state_ref:
                            break
                    else:
                        # Match found! Use the existing task name
                        next_tast_names[i] = name
                        new_tasks_already_calculated += 1
        nr_of_new_states_not_done = nr_of_new_states - new_tasks_already_calculated
        print("\nResults of tension analysis for %s:" % selected_project.name)
        print("  %i new protonation states found, %i of them not yet calculated.") \
             % (nr_of_new_states, nr_of_new_states_not_done)

        if nr_of_new_states > 0 and write_suggested_prot_states:
            pickle_filename = project_folder + 'prot_tasks_i%i.pickle' % next_iteration
            f = open(pickle_filename, 'w')
            next_titr_residue_dicts = titr_residue_dicts + accepted_titr_residue_dicts
            all_task_names = task_names + next_tast_names
            pickle.dump((all_task_names, next_titr_residue_dicts), f)
            f.close()
            print("  Updated protonation states writen to: %s" % pickle_filename)
            print("  Number of protonations: %i" % len(all_task_names))
        else:
            print("  Number of protonations: %i" % len(task_names))



    plt.show(True)