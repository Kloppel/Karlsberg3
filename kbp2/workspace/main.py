# -*- coding: utf-8 -*-
# coding=utf-8


# import sys
# sys.path.append('/user/tmeyer/workspace/script/protein_toolbox/packages')

import kbp2
import os, sys


# Todo: Proper check for CYS that are not in disu. bridges.

def start_md_pka_project(base_folder, jobname, protein, top, par, md_template_folder, roi=None):
    ### SETUP ###
    md_kbp2_template = kbp2.jobs.MD_kbp2_project(None, None)
    md_kbp2_template.top = top
    md_kbp2_template.par = par
    md_kbp2_template.roi = roi

    ####################################################
    ### 2) General modelling setup for the structure ###
    ####################################################
    project_folder = base_folder + jobname + '/'

    # DEBUG
    # import shutil, os
    # if os.path.exists(project_folder):
    #     shutil.rmtree(project_folder)

    ### SETUP ###
    # Create new md-pka project
    md_kbp2 = kbp2.jobs.MD_kbp2_project(project_folder, protein, template=md_kbp2_template)
    md_kbp2.md_template_folder = md_template_folder

    ### CHARMM ### <- This is the place to add any general modelling instructions.
    charmm_struct = md_kbp2.charmm_struct
    # Do not use calcium ions.
    charmm_struct.add_decision('rename__CA_CAL', 'keep')
    #############


    ################################
    ### Start the Karlsberg+ job ###
    ################################
    kbp_job_folder = kbp_cache_folder + jobname + '/'
    kbp_pdb = kbp_job_folder + md_kbp2.charmm_struct.charmm_out_prefix + '.pdb'
    print kbp_pdb
    status = ''
    if md_kbp2.is_written_prot_available():
        status = 'unnecessary'
    elif os.path.exists(kbp_pdb):
        # Check if it has finished
        kbp_result = kbp2.kbp_results.KbpResult()
        success = kbp_result.read_kbp_results(kbp_job_folder, allow_unfinished_job=True)

        # changed for PROPKA pKas
        # propka_folder = "/scratch/scratch/tmeyer/md_pka/propka_runs/"
        # propka_pkas_file = propka_folder + jobname + '.pka'
        # propka_pkas = kbp2.workspace.tools.parse_propka_pkas(propka_pkas_file)
        # kbp_result.pkas = propka_pkas

        if not success:
            # Karlsberg+ job is still running
            status = 'running'
            # Todo: Check for a crashed job
        else:
            # Karlsberg+ job done.
            status = 'done'
            md_kbp2.kbp_result = kbp_result
    else:
        # Model the structure
        kbp_modelling_folder = kbp_job_folder + 'modelling/'
        os.mkdir(kbp_job_folder)
        os.mkdir(kbp_modelling_folder)
        md_kbp2.charmm_struct.workdir = kbp_modelling_folder
        md_kbp2.charmm_struct.add_decision('rename__HOH_TIP3', 'keep')
        md_kbp2.charmm_struct.check_structures(quiet=True)

        md_kbp2.charmm_struct.run_charmm()
        kbp_struture = md_kbp2.charmm_struct.get_modelled_structure()

        kbp_struture.write_pdb(kbp_pdb)


    if status not in ['done', 'unnecessary']:
        # Prepare Karlsberg+ job
        kbp_job = kbp2.job_manager2.kbp_job(pdb=kbp_pdb, single_folder_run=True)
        kbp_job.kb_para['globals']['workDir'] = kbp_dolly_workdir
        kbp_job.parameter['run_dir'] = kbp_job_folder
        kbp_job.kb_para['sfc_globals']['preOpt'] = 1

        # Submit the job the the folder queue
        job = kbp2.queue.Folder_job('folder', kbp_queue_folder, jobname+'_pre_md', prefix='kbp_job_')
        job.set_job_object(kbp_job)
        job.submit()
        status = job.get_status()


    # print("Karlsberg+ pre MD run: %s" % status)
    if status not in ['done', 'unnecessary']:
        return status



    # if status == 'done':
    #     from matplotlib import pyplot as plt
    #     residue_nr = 2
    #     x = kbp_result.descr.ph_values
    #     y = kbp_result.deprot_curves[residue_nr]
    #     residue = kbp_result.descr.sorted_residue_list[residue_nr]
    #     print kbp_result.pkas[residue]
    #     plt.plot(x,y, 'x-r')
    #     # for i in range(len(kbp_result.occs[residue_nr])):
    #     #     plt.plot(x, kbp_result.occs[residue_nr][i], 'x-r')
    #     plt.title(residue)
    #     plt.ylim([-0.01, 1.01])
    #     plt.show()



    #####################################################################
    ### Define protonation states that should be used for the MD runs ###
    #####################################################################
    # If there are protonation instructions in the folder, then they are taken otherwise the default protocol is used.
    md_kbp2.setup_tasks()

    pickle_filename = project_folder + 'md_kbp2.pickle'
    md_kbp2.write_pickle(pickle_filename)

    pka_m = kbp2.md_pkas.pka_session()
    for task in md_kbp2.tasks:
        modelling_prefix = task.modelling_prefix
        task_folder = task.task_folder
        if task.md_prefix is None:
            task.md_prefix = modelling_prefix + '_' + task.taskname
        jobname = task.md_prefix
        structure_filename = md_kbp2.structure_filename

        task_md_template_folder = task.parent_project.md_template_folder
        msg = pka_m.add_job(modelling_prefix, task_folder, jobname=jobname, top=top, par=par,
                            md_template_folder=task_md_template_folder, structure_filename=structure_filename)

    pka_m.send_jobs_to_server()

    return pka_m



def print_job_status(job_status):
    for job in job_status:
        if type(job) == str:
            print job
            print('')
        elif job != None:
            print("PATH:        " + job.root_workdir)
            print("STATUS:      " + job.status)

            if job.md_status != 'done' and job.md_status != 'none':
                print(" MD:         %s -> %s" % (job.md_status, job.md_progress))
            else:
                print(" MD:         %s" % (job.md_status))

            if job.md_problem:
                print(" MD - problem : " + job.md_problem)

            if job.kbp_status != 'done' and job.kbp_status != 'none':
                print(" Karlsberg+: %s -> %s" % (job.kbp_status, job.kbp_progress))
            else:
                print(" Karlsberg+: %s" % (job.kbp_status))

            if job.kbp_problem:
                print(" Karlsberg+ - problem : " + job.kbp_problem)

            if job.apbs_status != 'done' and job.apbs_status != 'none':
                print(" APBS:        %s -> %s" % (job.apbs_status, job.apbs_progress))
            else:
                print(" APBS:       %s" % (job.apbs_status))

            print('')
        else:
            print "None returned."
            print ''

### Rewrite analysis script:
# Make it more efficient and storage model consistent.
# Create new guess from tension analysis.


if __name__ == '__main__':

    while True:

        #########################
        ### 1 ) General Setup ###
        #########################
        top = []
        top.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/top.inp")
        top.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.rtf")
        par = []
        par.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar_36/par_all22_prot.inp")
        par.append("/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.prm")

        md_template_folder_sb = '/user/tmeyer/workspace/projects/md_pkas/snase/templ_clarge_5_10ns_100ps/'
        md_template_folder_sb_long = '/user/tmeyer/workspace/projects/md_pkas/snase/templ_clarge_5_40ns_50ps/'
        md_template_folder_sb_3000 = '/user/tmeyer/workspace/projects/md_pkas/snase/templ_clarge_5_10ns_100ps_3000/'
        md_template_folder_sb_more = '/user/tmeyer/workspace/projects/md_pkas/snase/templ_clarge_5_10ns_50ps/'
        md_template_folder_snase = '/user/tmeyer/workspace/projects/md_pkas/snase/templ_clarge_20ns_200ps/'
        md_template_folder_hema = '/user/tmeyer/workspace/projects/md_pkas/snase/templ_20ns_200ps/'

        # folder_general = '/scratch/scratch/tmeyer/md_pka/runs/general_testing/'
        folder_general = '/scratch/scratch/tmeyer/md_pka/runs/general2/'
        folder_snase = '/scratch/scratch/tmeyer/md_pka/runs/snase_nocal2/'
        folder_hema = '/scratch/scratch/tmeyer/md_pka/runs/hema/'

        #Todo: Are set via global variables!!
        # kbp_cache_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache/'
        kbp_cache_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_po1_ie80/'
        kbp_queue_folder = '/scratch/scratch/tmeyer/md_pka/md_pka_manager/queue/kbp_jobs/'
        kbp_dolly_workdir = '/public/scratch/tmeyer/kb/pre_md/'



        std_projects = []

        project_name = '3bdc'
        protein = '/scratch/scratch/pdb/pdb_bio_merged/bd/3bdc.pdb1'
        std_projects.append((project_name, protein))

        project_name = '2lzt'
        protein = '/scratch/scratch/pdb/pdb_bio_merged/lz/2lzt.pdb1'
        std_projects.append((project_name, protein))

        project_name = '1a2p_cm'
        protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/1a2p_cm.pdb'
        std_projects.append((project_name, protein))

        project_name = '3rn3'
        protein = '/scratch/scratch/pdb/pdb_bio_merged/rn/3rn3.pdb1'
        std_projects.append((project_name, protein))

        project_name = '3icb'
        protein = '/scratch/scratch/pdb/pdb_bio_merged/ic/3icb.pdb1'
        std_projects.append((project_name, protein))

        project_name = '2rn2'
        protein = '/scratch/scratch/pdb/pdb_bio_merged/rn/2rn2.pdb1'
        std_projects.append((project_name, protein))

        project_name = '4pti'
        protein = '/scratch/scratch/pdb/pdb_bio_merged/pt/4pti.pdb1'
        std_projects.append((project_name, protein))

        # project_name = '1xnb'
        # protein = '/scratch/scratch/pdb/pdb_bio_merged/xn/1xnb.pdb1'
        # std_projects.append((project_name, protein))
        project_name = '1xnb_35sw'
        protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/1xnb_35sw.pdb'
        std_projects.append((project_name, protein))

        # 1ppf is a complex. Only chain I is "Third domain of the turkey ovomucoid inhibitor (OMTKY3)"
        project_name = '1ppf_i'
        protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/1ppf_i.pdb1'
        std_projects.append((project_name, protein))

        project_name = '1pga'
        protein = '/scratch/scratch/pdb/pdb_bio_merged/pg/1pga.pdb1'
        std_projects.append((project_name, protein))

        # # Wrong N-terminus and pKa calculation was done with preopt=0 what cased the asymmetry of the result.
        # # project_name = '2zta_cm2'
        # # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/2zta_cm2.pdb'
        # # std_projects.append((project_name, protein))
        project_name = '2zta_cm3'
        protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/2zta_cm3.pdb'
        std_projects.append((project_name, protein))

        project_name = '1ert_a'
        protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/1ert_a.pdb'
        std_projects.append((project_name, protein))

        # # project_name = '1hng_cm'
        # # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/1hng_cm.pdb'
        # # std_projects.append((project_name, protein))
        # # project_name = '1hng_99c_a'
        # # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/1hng_99c_a.pdb'
        # # std_projects.append((project_name, protein))
        project_name = '1hng_99c2_a'
        protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/1hng_99c2_a.pdb'
        std_projects.append((project_name, protein))

        for project_name, protein in std_projects:
            print("### Starting project %s: ###" % project_name)
            pka_m = start_md_pka_project(folder_general, project_name, protein, top, par, md_template_folder_sb)
            # pka_m = start_md_pka_project(folder_general, project_name, protein, top, par, md_template_folder_sb_3000)
            # pka_m = start_md_pka_project(folder_general, project_name, protein, top, par, md_template_folder_sb_more)
            if type(pka_m) != str:
                print_job_status(pka_m.get_status())
            else:
                print pka_m


        #############
        ### SNASE ###
        #############
        snase_projects = []

        # project_name = '3d6c'
        # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/3d6c.pdb1'
        # roi = [('GLU', 38, 'A')]
        # snase_projects.append((project_name, protein, roi))
        #
        # project_name = '2oxp_m'
        # #protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/2oxp.pdb1'
        # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/2oxp_m.pdb'
        # roi = [('ASP', 66, 'A')]
        # snase_projects.append((project_name, protein, roi))
        #
        # project_name = '3c1e'
        # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/3c1e.pdb1'
        # roi = [('LYS', 125, 'A')]
        # snase_projects.append((project_name, protein, roi))
        #
        # project_name = '3erq'
        # protein = '/scratch/scratch/pdb/pdb_bio_merged/er/3erq.pdb1'
        # roi = [('LYS', 25, 'A')]
        # snase_projects.append((project_name, protein, roi))
        #
        # project_name = '3evq'
        # protein = '/scratch/scratch/pdb/pdb_bio_merged/ev/3evq.pdb1'
        # roi = [('GLU', 25, 'A')]
        # snase_projects.append((project_name, protein, roi))
        #
        # project_name = '3dmu_m'
        # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/3dmu_m.pdb'
        # roi = [('LYS', 62, 'A')]
        # snase_projects.append((project_name, protein, roi))
        #
        # project_name = '2rbm_m'
        # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/2rbm_m.pdb'
        # roi = [('LYS', 72, 'A')]
        # snase_projects.append((project_name, protein, roi))
        #
        # project_name = '3ero'
        # protein = '/scratch/scratch/pdb/pdb_bio_merged/er/3ero.pdb1'
        # roi = [('GLU', 72, 'A')]
        # snase_projects.append((project_name, protein, roi))
        #
        # project_name = '3c1f'
        # protein = '/scratch/scratch/pdb/pdb_bio_merged/c1/3c1f.pdb1'
        # roi = [('LYS', 104, 'A')]
        # snase_projects.append((project_name, protein, roi))
        #
        # project_name = '3eji_m'
        # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/3eji_m.pdb'
        # roi = [('LYS', 36, 'A')]
        # snase_projects.append((project_name, protein, roi))



        # # experiments say: global unfolding
        # project_name = '2oeo'
        # protein = '/scratch/scratch/pdb/pdb_bio_merged/oe/2oeo.pdb1'
        # roi = [('ASP', 92, 'A')]
        # snase_projects.append((project_name, protein, roi))
        # # experiments say: global unfolding
        # project_name = '1tt2'
        # protein = '/scratch/scratch/pdb/pdb_bio_merged/tt/1tt2.pdb1'
        # roi = [('LYS', 92, 'A')]
        # snase_projects.append((project_name, protein, roi))

        # Dimer
        # # project_name = '3e5s'
        # # protein = '/scratch/scratch/pdb/pdb_bio_merged/e5/3e5s.pdb1'
        # # snase_projects.append((project_name, protein))

        # # Dimer
        # # # project_name = '3d4d'
        # # # protein = '/scratch/scratch/pdb/pdb_bio_merged/d4/3d4d.pdb1'
        # # # snase_projects.append((project_name, protein))
        # # #
        # # #

        # # project_name = '1u9r_m'
        # # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/1u9r_m.pdb'
        # # snase_projects.append((project_name, protein))

        # #
        # # project_name = '1tqo'
        # # protein = '/scratch/scratch/pdb/pdb_bio_merged/tq/1tqo.pdb1'
        # # snase_projects.append((project_name, protein))
        # #

        # # >= 10.4
        # # project_name = '2rks_m'
        # # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/2rks_m.pdb'
        # # roi = [('LYS', 38, 'A')]
        # # snase_projects.append((project_name, protein, roi))


        import time
        c = 0
        for project_name, protein, roi in snase_projects:
            print("### Starting project %s: ###" % project_name)
            pka_m = start_md_pka_project(folder_snase, project_name, protein, top, par, md_template_folder_snase,roi=roi)
            if type(pka_m) != str:
                print_job_status(pka_m.get_status())
            else:
                print pka_m
            print

            c += 1
            if c == 2:
                c = 0
                # time.sleep(180)
        print "All jobs submitted"


        std_long_projects = []
        # project_name = '2lzt_l'
        # protein = '/scratch/scratch/pdb/pdb_bio_merged/lz/2lzt.pdb1'
        # std_long_projects.append((project_name, protein))

        # project_name = '2zta_cm3_l'
        # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/2zta_cm3.pdb'
        # std_long_projects.append((project_name, protein))


        for project_name, protein in std_long_projects:
            print("### Starting project %s: ###" % project_name)
            pka_m = start_md_pka_project(folder_general, project_name, protein, top, par, md_template_folder_sb_long)
            if type(pka_m) != str:
                print_job_status(pka_m.get_status())
            else:
                print pka_m
            print



        ############
        ### HEMA ###
        ############
        # hama_projects = []
        #
        # project_name = '2ibx_cm2'
        # protein = '/scratch/scratch/tmeyer/md_pka/mod_bio/2ibx_cm2.pdb'
        # hama_projects.append((project_name, protein))
        #
        # for project_name, protein in hama_projects:
        #     print("### Starting project %s: ###" % project_name)
        #     pka_m = start_md_pka_project(folder_hema, project_name, protein, top, par, md_template_folder_hema)
        #     if type(pka_m) != str:
        #         print_job_status(pka_m.get_status())
        #     else:
        #         print pka_m
        #     print

        # import time
        # print "Sleeping for 10 min"
        # time.sleep(600)
        # print "Sleeping for 1 hour"
        # time.sleep(3600)
        break


