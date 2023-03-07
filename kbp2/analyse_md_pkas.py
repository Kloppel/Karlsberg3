# -*- coding: utf-8 -*-

import cPickle as pickle
import os, re
from kbp2 import kbp_results, apbs_manager, kbp_tools


def read_project(project_folder, kbp_subfolder, sasa_factor=0.0, quiet=False, max_iteration=-1,
                 charmm_confe_subfolder=None, protonation_scheme='prot_tasks'):

    ################
    ### Settings ###
    ################

    pickle_kbp_results = 'kbp_results.pickle'
    pickle_kbp_results_avg_only = 'kbp_results_avg_only.pickle'
    pickle_prot_tasks_base = protonation_scheme
    apbs_suffix = ''
    epsilon = 4.0

    # Only working with read_all = True, since it will be overwritten by the content of the pickle file.
    frame_range = [10, -1]
    # frame_range = [30, -1]
    # frame_range = [10, 55]
    # frame_range = [55, -1]
    if frame_range != [10, -1]:
        print "frame range: " + str(frame_range)

    # The following flags are not active at the moment
    # Weighting of the pKas within an MD. Only for display and RMSDs of individual MDs. Not used for the total RMSD.
    # Choices are: 'avg', 'weighted'
    #weight_within_md = 'avg'

    # Weighting of pKas and MDs for the final pKas (and the total RMSD).
    # Choices are:
    # 'avg'      : pKas, protonation- and conformational energies are averaged within a MD. Total pKa is obtained with
    #              a Boltzmann sum over all MDs. The weight for a MD is given by the sum of protonation- and
    #              conformational energy.
    # 'weighted' : The pKas are obtained with a boltzmann sum over all frames in all MDs. The weight of each frame is
    #              given by the sum of protonation- and conformational energy.
    # 'mc'       : Use the Monte Carlo algorithm provided by 'Karlsberg' to obtain the pKas. All .g and .pkint files,
    #              as well as the conformational energies are provided to the Karlsberg binary. Takes a while.
    #weighting = 'avg'

    #############
    ### Flags ###
    #############
    read_all = False

    reread_results = False

    ################
    ### Analysis ###
    ################
    if project_folder[-1] != '/':
        project_folder += '/'

    # Find latest iteration
    files_in_project_folder = os.listdir(project_folder)
    reg = re.compile(r'^' + re.escape(pickle_prot_tasks_base) + r'_?i?(\d*).pickle$')
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
    # Read protonation tasks
    if last_iteration != 0:
        pickle_prot_tasks = '%s_i%i.pickle' % (pickle_prot_tasks_base, last_iteration)
        # print("\nReading protonation from %s" % pickle_prot_tasks)
    else:
        pickle_prot_tasks = pickle_prot_tasks_base + '.pickle'
        # print("\nReading protonation from %s" % pickle_prot_tasks)
    prot_tasks_filename = project_folder + pickle_prot_tasks
    f = open(prot_tasks_filename, 'r')
    task_names, titr_residue_dicts = pickle.load(f)
    f.close()
    # prot_tasks_filename = project_folder + pickle_prot_tasks
    # f = open(prot_tasks_filename, 'r')
    # task_names, titr_residue_dicts = pickle.load(f)
    # f.close()

    if not quiet:
        print("# Reading protonation %s from project: %s" % (pickle_prot_tasks, project_folder))

    all_task_results = []
    for task_name in task_names:
        if task_name == "ph7_1_26p":
            print "ph7_1_26p skipped"
            continue

        task_folder = project_folder + task_name + '/'

        kbp_folder = task_folder + kbp_subfolder + '/'
        if not os.path.exists(kbp_folder):
            print("Skipping task %s, it has not finished yet" % (project_folder + task_name,))
            continue

        pickel_filename_avg_only = task_folder + kbp_subfolder + '/' + pickle_kbp_results_avg_only
        if os.path.exists(pickel_filename_avg_only) and not read_all and not reread_results:
            if not quiet:
                print("Previous averaged results found for task %s -> reading pickle" % task_name)
            f = open(pickel_filename_avg_only, 'rb')
            task_results = pickle.load(f)
            f.close()
        else:
            pickel_filename = task_folder + kbp_subfolder + '/' + pickle_kbp_results
            if os.path.exists(pickel_filename) and not reread_results:
                if not quiet:
                    print("Previous results found for task %s -> reading pickle" % task_name)
                f = open(pickel_filename, 'rb')
                task_results = pickle.load(f)
                f.close()

                msg_printed = False
                for frame_nr, kbp_task in task_results:
                    if kbp_task.residue_sasa is None:
                        if not msg_printed:
                            print "calculating SASA"
                            msg_printed = True
                        kbp_task.calc_residue_sasa()
                f = open(pickel_filename, 'wb')
                pickle.dump(task_results, f, protocol=2)
                f.close()
            else:
                task_results = kbp_results.FrameKbpResults()
                if not quiet:
                    print("Reading Karlsberg+ results for: '%s'" % task_folder)

                ###################
                ### Karlsberg+  ###
                ###################
                main_frame_folder = task_folder + kbp_subfolder + '/done/'
                folders = os.listdir(main_frame_folder)
                reg = re.compile(r'^frame(\d+)$')
                for folder in folders:
                    reg_m = reg.match(folder)
                    if reg_m is None:
                        continue
                    frame_nr = int(reg_m.groups()[0])
                    frame_folder = main_frame_folder + folder + '/'

                    frame_results = kbp_results.KbpResult()
                    frame_results.read_kbp_results(frame_folder, detailed=True)
                    task_results.add_task(frame_results, frame_nr)

                    kb_energies = kbp_tools.parse_kb_energy(frame_folder)
                    task_results.set_misc_info_energies('kb_energy', frame_nr, kb_energies)

                if not quiet:
                    print(" -> %s frames found" % task_results.get_nr_of_frames())

                ############
                ### APBS ###
                ############
                residue_list = task_results.descr.sorted_residue_list
                kbp_folder = task_folder + kbp_subfolder + '/'
                if not quiet:
                    print("Reading APBS results for: %s" % kbp_folder)
                (unfinished_jobs, crashed_jobs, finished_jobs) = \
                    apbs_manager.read_results_res(kbp_folder, residue_list, subfolder_suffix=apbs_suffix, epsilon=epsilon)
                conf_energy_t = [(int(frame.replace('frame', '')), (s, c, sc_sites, near_sites))
                                 for (frame, (s, c, sc_sites, near_sites)) in finished_jobs.iteritems()]
                conf_energy_t.sort()

                nr_of_apbs_results = len(conf_energy_t)
                if nr_of_apbs_results != task_results.get_nr_of_frames():
                    error = "The number of APBS results (%i) does not match the number of Karlsberg+ results (%i)." % \
                            (nr_of_apbs_results, task_results.get_nr_of_frames())
                    raise AssertionError(error)

                for frame_nr, energies in conf_energy_t:
                    task_results.set_task_conf_energies('solvation', frame_nr, energies[0])
                    task_results.set_task_conf_energies('coulomb', frame_nr, energies[1])
                    task_results.set_task_conf_energies('residue_energy', frame_nr, energies[2])
                    # task_results.set_misc_info_energies('residue_energy', frame_nr, energies[2])
                    task_results.set_misc_info_energies('all_residues_in_range', frame_nr, energies[3])


                ############
                ### SASA ###
                ############
                # if not quiet:
                #     print("Calculating SASA for: %s" % kbp_folder)
                # # Calculate the SASA of all titratable residues
                for frame_nr, kbp_task in task_results:
                    kbp_task.calc_residue_sasa()

                ############################################
                ### Conformational energy for each frame ###
                ############################################
                # # Calculate and set the conformational energy for each frame.
                # for frame_nr, kbp_task in task_results:
                #     solv = kbp_task.conf_energies['solvation']
                #     coulomb = kbp_task.conf_energies['coulomb']
                #     sasa_correction = 0.0
                #     for i, reside_descr in enumerate(kbp_task.descr.sorted_residue_list):
                #         sasa_correction += -0.5 * kbp_task.residue_sasa[i]
                #     kbp_task.conf_energy = solv + coulomb + sasa_correction

                f = open(pickel_filename, 'wb')
                pickle.dump(task_results, f, protocol=2)
                f.close()

            #####################
            ### CHARMM energy ###
            #####################
            if charmm_confe_subfolder is not None:
                charmm_energy_folder = kbp_folder + charmm_confe_subfolder + '/'
                charmm_energies = read_charmm_results(charmm_energy_folder)
                for frame_nr in range(task_results.get_nr_of_frames()):
                    for energy_type, value in charmm_energies.iteritems():
                        task_results.set_task_conf_energies(energy_type, frame_nr, value[frame_nr])


            # Calculate averages and pickle the results
            task_results.average_over_frames(frame_range)

            task_results.avg_results.calc_protonation_energy()

            task_results.avg_results.find_pkas()

            task_results.avg_results.folder = task_folder

            # Store the averaged results in a separate file.
            tmp_kbp_results = task_results.kbp_results
            task_results.kbp_results = None
            f = open(pickel_filename_avg_only, 'wb')
            pickle.dump(task_results, f, protocol=2)
            f.close()
            task_results.kbp_results = tmp_kbp_results

        #####################################
        ### Average conformational energy ###
        #####################################
        solv = task_results.avg_results.conf_energies['solvation']
        coulomb = task_results.avg_results.conf_energies['coulomb']
        # sasa = task_results.avg_results.conf_energies['sasa']

        # sasa_correction = 0.0
        # sasa_factor = -0.4
        # for i, residue_descr in enumerate(task_results.avg_results.descr.sorted_residue_list):
        #     resname, resid, segname = re.split(r'[-_]', residue_descr)
        #
        #     # if resname in ['EPP', 'DPP']:
        #     # print residue_descr
        #     # if int(resid) == 38:
        #     #     print task_results.avg_results.residue_sasa[i]
        #     sasa_correction += sasa_factor * task_results.avg_results.residue_sasa[i]
        # # task_results.avg_results.conf_energy = solv + coulomb + sasa_correction


            # # import numpy as np
            # # R = 8.3144621                             # J/(K*mol)
            # RT_ln10 = R * 300.0 / 1000.0 * np.log(10) # kJ/mol
            # for state, residue_def  in enumerate(task_results.avg_results.descr.titratable_residues[resname]):
            #
            #     # if residue_def['name'] not in ['R', '0']:
            #     if residue_def['name'] not in ['R']:
            #         sasa_factor  = 0.5
            #         if resname in ['EPP', 'DPP', 'CTE']:
            #             sasa_factor *= 1.0
            #         elif resname in ['HSP', 'LYS', 'ARG', 'NTE']:
            #             sasa_factor *= -1.0
            #         else:
            #             sasa_factor = 0.0
            #         task_results.avg_results.pkaint[residue_descr][state] += \
            #             sasa_factor * task_results.avg_results.residue_sasa[i]
            #         # print residue_descr
            #         # print task_results.avg_results.pkaint[residue_descr][state] / RT_ln10 + 100


        charmm_elec = task_results.avg_results.conf_energies['ELEC']
        # charmm_vdw = task_results.avg_results.conf_energies['VDWaals']

        # task_results.avg_results.conf_energy = solv + coulomb
        task_results.avg_results.conf_energy = solv + charmm_elec
        # task_results.avg_results.conf_energy = solv + charmm_elec + sasa_correction
        # task_results.avg_results.conf_energy = charmm_elec
        # task_results.avg_results.conf_energy = solv + charmm_elec + charmm_vdw
        # task_results.avg_results.conf_energy = solv
        # task_results.avg_results.conf_energy = solv + charmm_vdw + charmm_elec

        # print task_results.avg_results.conf_energy
        # print sasa_correction
        # print "coulomb_elec"
        # print "charmm_elec_w1-4"
        # print "charmm_elec"
        # print "charmm_elec + vdW"

        all_task_results.append(task_results)


    combined_project_results = kbp_results.FrameKbpResults()
    for task_results in all_task_results:
        combined_project_results.add_task(task_results.avg_results)

    ### 1. Method ###
    # combined_project_results.combine_frames()

    ### 2. Method ###
    combined_project_results.combine_frames_karlsberg(cpus=3)

    ### 3. Method ###
    # for md_nr, avg_results in combined_project_results:
    #     avg_results.karlsberg_titration(cpus=3)
    # combined_project_results.combine_frames()


    pkas = combined_project_results.combined_results.pkas
    

    return pkas, combined_project_results


class ProjectResults(object):
    def __init__(self, project_name):
        self.name = project_name
        self.pkas = None
        self.exp_pkas = None
        self.old_pkas = None
        self.path = None
        self.sasa_factor = None

        # FrameKbpResults object
        self.detailed_project_results = None

    def set_exp_pkas(self, exp_pkas):
        if type(exp_pkas) == str:
            exp_pkas = kbp_tools.parse_gernots_exp_pkas(exp_pkas)



        # exp_pkas_all = exp_pkas
        # exp_pkas = {}
        # for residue_descr, exp_pka in exp_pkas_all.iteritems():
        #     resname, resid, segname = re.split(r'[-_]', residue_descr)
        #     mod_pkas = {'ASP' : 4.0,
        #                 'GLU' : 4.4,
        #                 'ARG' : 12.0,
        #                 'LYS' : 10.4,
        #                 'HIS' : 7.0,
        #                 'TYR' : 9.6,
        #                 'CYS' : 9.5,
        #                 'NTE' : 7.5,
        #                 'CTE' : 3.8}
        #     if abs(exp_pka - mod_pkas[resname]) > 1.0:
        #         exp_pkas[residue_descr] = exp_pka
            # exp_pkas[residue_descr] = mod_pkas[resname]



        self.exp_pkas = exp_pkas

    def set_old_pkas(self, old_pkas):
        self.old_pkas = old_pkas

    def read_results(self, path, kbp_subfolder, sasa_factor=0.0, quiet=False, max_iteration=-1,
                     charmm_confe_subfolder=None, protonation_scheme='prot_tasks'):
        self.path = path
        self.sasa_factor = sasa_factor
        pkas, detailed_project_results = read_project(path, kbp_subfolder, sasa_factor, quiet=quiet,
                                                      max_iteration=max_iteration, charmm_confe_subfolder=charmm_confe_subfolder,
                                                      protonation_scheme=protonation_scheme)
        self.pkas = pkas
        self.detailed_project_results = detailed_project_results

def read_charmm_results(workdir):

    if workdir[-1] != '/':
        workdir += '/'

    frame_workdirs = os.listdir(workdir)
    nr_and_folder = []
    for frame_workdir in frame_workdirs:
        if not 'frame' in frame_workdir:
            continue
        frame_nr = int(frame_workdir.strip('frame'))
        nr_and_folder.append((frame_nr, workdir + frame_workdir + '/'))

    charmm_energies = {}
    nr_and_folder.sort()
    for frame_nr, frame_workdir in nr_and_folder:
        charmm_out_filename = frame_workdir + "frame%i_charmm.out" % frame_nr
        # charmm_out_filename = frame_workdir + "c_pH7_frame%i.reference_charmm.out" % 59
        f = open(charmm_out_filename)
        for line in f:
            if 'ENER ENR:' in line:
                types_and_names = []
                types_and_values = []
                while not '----------' in line:
                    line = line.strip('\n')
                    types_and_names.append(line.split(':'))
                    line = f.next()
                line = f.next()
                while not '----------' in line:
                    types_and_values.append(line.split('>'))
                    line = f.next()
                break
        for (types_descr, names), (types, values) in zip(types_and_names, types_and_values):
            names_list = names.split()
            values_list = values.split()
            for name, value in zip(names_list, values_list):
                value = float(value)
                # Convert to kJ/mol
                value *= 4.184
                if not name in charmm_energies:
                    charmm_energies[name] = [value]
                else:
                    charmm_energies[name].append(value)
        f.close()
    return charmm_energies


    # energies_to_plot= ['ENERgy', 'ELEC', 'VDWaals', 'BONDs', 'ANGLes', 'DIHEdrals', 'IMPRopers']
    # from matplotlib import pyplot as plt
    # plt.figure()
    # for energy_to_plot in energies_to_plot:
    #     y = charmm_energies[energy_to_plot]
    #     plt.plot(y)
    # plt.legend(energies_to_plot)#
    # plt.show()

    # ENER ENR:  Eval#     ENERgy      Delta-E         GRMS
    # ENER INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
    # ENER CROSS:           CMAPs        PMF1D        PMF2D        PRIMO
    # ENER EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
    #  ----------       ---------    ---------    ---------    ---------    ---------
    # ENER>        0   1522.52937      0.00000     18.03137
    # ENER INTERN>      544.02813   1124.93499    169.23668    713.96413     66.85994
    # ENER CROSS>      -176.10529      0.00000      0.00000      0.00000
    # ENER EXTERN>     -468.39449   -451.99472      0.00000      0.00000      0.00000
    #  ----------       ---------    ---------    ---------    ---------    ---------




















