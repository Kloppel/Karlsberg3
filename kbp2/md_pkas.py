# -*- coding: utf-8 -*-
import re



import sys, os, shutil, time
sys.path.append('/user/tmeyer/workspace/script/protein_toolbox')

# from kbp2 import file_parser, charmm, job_manager2, folder_queue
from kbp2 import job_manager2, kbp_tools
from kbp2 import charmm
from kbp2 import file_parser
from kbp2 import apbs_manager
# import apbs_manager_exp2 as apbs_manager

import subprocess

import cPickle as pickle

from multiprocessing.connection import Client, Listener
from multiprocessing import Process, Queue, Manager

def tail(f, window=20):
    """
    Returns the last `window` lines of file `f` as a list.
    """
    BUFSIZ = 1024
    f.seek(0, 2)
    bytes = f.tell()
    size = window + 1
    block = -1
    data = []
    while size > 0 and bytes > 0:
        if bytes - BUFSIZ > 0:
            # Seek back one whole BUFSIZ
            f.seek(block * BUFSIZ, 2)
            # read BUFFER
            data.insert(0, f.read(BUFSIZ))
        else:
            # file too small, start from begining
            f.seek(0,0)
            # only read what was not read
            data.insert(0, f.read(bytes))
        linesFound = data[0].count('\n')
        size -= linesFound
        bytes -= BUFSIZ
        block -= 1
    return ''.join(data).splitlines()[-window:]


class MD_pKa_job(object):
    def __init__(self, modelling_prefix ,root_workdir, jobname, md_queue=None, kbp_queue=None, top=[], par=[]\
                 , md_template_folder=None, ph=7, structure_filename=''):
        """
        Initialize new job.
        """
        if root_workdir[-1] != '/':
            root_workdir += '/'

        # Not used, but is helpful to identify the job.
        self.pdb = structure_filename
        self.modelling_prefix = modelling_prefix
        self.root_workdir = root_workdir
        self.jobname      = jobname
        self.top          = top
        self.par          = par
        assert(type(ph) == int)
        self.ph           = ph

        # This are the Queue() endpoints to communicate with the main process. They are set externally if the
        # job is started.
        self.status_queue = None
        self.cmd_queue    = None


        ############################
        ### Modelling Parameters ###
        ############################
        self.modelling_folder = root_workdir + 'modelling/'
        # Can be: 'none' / 'crashed' / 'done'
        self.modelling_status = 'none'
        # Dict that stores information about the modelled structure.
        self.charmm_manager_data = {}

        #####################################
        ### Molecular Dynamics Parameters ###
        #####################################
        self.md_folder    = root_workdir + 'md/'
        # Path to the folder that contains the job queue for the MD jobs.
        if md_queue is None:
            self.md_queue = '/scratch/scratch/tmeyer/md_pka/md_pka_manager/queue/d58/'
        else:
            self.md_queue = md_queue
        # Can be: 'none' / 'preparing' / 'running' / 'crashed' / 'done'
        self.md_status    = 'none'
        # If something is wrong and the reason was found it is documented here.
        self.md_problem   = ''
        self.md_progress  = ''

        ############################
        ### Karlsberg+ Parameter ###
        ############################
        # old:
        # self.kbp_folder_suffix   = 'kbp_noMin_refzero/'
        # self.kbp_folder_suffix   = 'kbp_e1_refzero/'
        # self.kbp_folder_suffix   = 'kbp_e4_refzero/'

        # experimental:
        # self.kbp_folder_suffix   = 'kbp_noMin_all/'

        # std:
        # self.kbp_folder_suffix   = 'kbp_noMin/'
        # self.kbp_folder_suffix   = 'kbp_e1/'
        # self.kbp_folder_suffix   = 'kbp_e4/'


        ## self.kbp_folder_suffix   = 'kbp_noMin_refzero_ter2/'
        ## self.kbp_folder_suffix   = 'kbp_e1_refzero_ter2/'
        # self.kbp_folder_suffix   = 'kbp_e4_refzero_ter2/'
        # self.kbp_folder_suffix   = 'kbp_noMin_refzero_ter2_np/'
        # self.kbp_folder_suffix   = 'kbp_e1_refzero_ter2_np/'
        # self.kbp_folder_suffix   = 'kbp_e4_refzero_ter2_np/'

        # Partially calculated
        # self.kbp_folder_suffix   = 'kbp_e4_refzero_pe2/'


        # self.kbp_folder_suffix   = 'kbp_noMin_rashin_c/'
        ## self.kbp_folder_suffix   = 'kbp_noMin_refzero_rashin2_ter2/'
        ## self.kbp_folder_suffix   = 'kbp_e4_refzero_rashin2_ter2/'
        # self.kbp_folder_suffix   = 'kbp_noMin_refzero_rashin2_ter2_np/'
        # self.kbp_folder_suffix   = 'kbp_e4_refzero_rashin2_ter2_np/'

        # self.kbp_folder_suffix   = 'kbp_e4_refzero_ter2_np_propka/'
        # self.kbp_folder_suffix   = 'kbp_e4_refzero_ter2_np_fit/'

        # SNase:
        # self.kbp_folder_suffix   = 'kbp_noMin_refzero/'
        # self.kbp_folder_suffix   = 'kbp_noMin_refzero_cav09/'
        # self.kbp_folder_suffix   = 'kbp_noMin_refzero_cav08/'
        # self.kbp_folder_suffix   = 'kbp_noMin_refzero_cav07/'

        # self.kbp_folder_suffix   = 'kbp_e4_refzero_cav09/'
        # self.kbp_folder_suffix   = 'kbp_e4_refzero_cav07/'

        # self.kbp_folder_suffix   = 'kbp_e4_refzero_pe8/'

        # self.kbp_folder_suffix   = 'kbp_noMin/'
        # self.kbp_folder_suffix   = 'kbp_noMin_cav09/'
        # self.kbp_folder_suffix   = 'kbp_noMin_cav08/'
        # self.kbp_folder_suffix   = 'kbp_noMin_cav07/'

        # self.kbp_folder_suffix   = 'kbp_e4_refzero_cav09_w2/'
        self.kbp_folder_suffix   = 'kbp_e4_refzero_cav09/'
        # self.kbp_folder_suffix   = 'kbp_e4_refzero_cav07_w2/'

        # self.kbp_folder_suffix   = 'kbp_e4_refzero_cav09_w2_pe8/'



        # self.min = 0
        # self.min = 1
        self.min = 4

        # self.vdw = True
        self.vdw = False



        self.ref0 = True
        # self.ref0 = False


        # self.cavity_parameter = None
        self.cavity_parameter = 0.9
        # self.cavity_parameter = 0.8
        # self.cavity_parameter = 0.7

        self.water_template = False



        # self.kbp_folder_suffix   = 'kbp_e4_noMin_new2bb_wi/'
        # self.kbp_folder_suffix   = 'kbp_e4_noMin_new2bb_wi_3000/'
        # self.kbp_folder_suffix   = 'kbp_e1_noMin_new2bb_wi/'
        # self.kbp_folder_suffix   = 'kbp_noMin_new2bb_wi/'

        # self.kbp_folder_suffix   = 'kbp_e4_noMin_newi/'
        # self.kbp_folder_suffix   = 'kbp_e4_noMin_200/'

        # self.kbp_folder_suffix   = 'kbp_e4_cav_noMin_new2bb_wi/'
        # self.kbp_folder_suffix   = 'kbp_e4_noMin_new2bb_wiw/'
        # self.kbp_folder_suffix   = 'kbp_e4_noMin_pkinte/'
        # self.kbp_folder_suffix   = 'kbp_e4_noMin_new_cav/'
        # self.kbp_folder_suffix   = 'kbp_pe2_e4_noMin_new/'

        #self.kbp_folder_suffix   = 'kbp_std_new2/'
        #self.kbp_folder_suffix   = 'kbp/'
        #self.kbp_folder_suffix   = 'kbp_pe6/'
        #self.kbp_folder_suffix   = 'kbp_pe2/'
        #self.kbp_folder_suffix   = 'kbp_noMin/'
        #self.kbp_folder_suffix   = 'kbp_e4_cav/'
        #self.kbp_folder_suffix   = 'kbp_e4_new/'

        self.kbp_folder = root_workdir + self.kbp_folder_suffix

        # Submitted Karlsberg+ jobs: {<name of the submitted job> : <Karlsberg+ job object>}
        self.kbp_jobs_prepared = {}
        # Path to the folder that contains the job queue for the Karlsberg+ jobs.
        if kbp_queue is None:
            self.kbp_queue = '/scratch/scratch/tmeyer/md_pka/md_pka_manager/queue/kbp_jobs/'
        else:
            self.kbp_queue = kbp_queue
        # Can be: 'none' / 'preparing' / 'running' / 'crashed' / 'done'
        self.kbp_status    = 'none'
        self.kbp_problem   = ''
        self.kbp_progress  = ''
        # Contains information about failed jobs.
        # Syntax: [ ([str:error message], kbp-job object) ]
        self.kbp_errors = []


        ######################
        ### APBS Parameter ###
        ######################
        #self.apbs_subfolder_suffix = '_e2'
        self.apbs_subfolder_suffix = ''
        # Can be: 'none' / 'preparing' / 'running' / 'crashed' / 'done'
        self.apbs_status    = 'none'
        self.apbs_problem   = ''
        self.apbs_progress  = ''



        #########################
        ### Results Parameter ###
        #########################
        self.results_folder = root_workdir + 'results/'
        # Can be: 'none' / 'done'
        self.result_status  = 'none'


        # The jobs status. Can be set to: 'ready' / 'running' / 'stopping' / 'stopped' / 'crashed' / 'done'
        self.status = 'ready'

        # All actions are logged here.
        self.log = []

        # Is set to True if the job is new.
        # Is set to False if the job is dumped or has been restored from a file.
        self.new = True

        # The folder containing the external_scripts for the molecular dynamic simulation.
        if md_template_folder is None:
            self.md_template_folder = '/user/tmeyer/workspace/projects/md_pkas/snase/templ/'
        else:
            self.md_template_folder = md_template_folder
            if self.md_template_folder[-1] != '/':
                self.md_template_folder += '/'

    def replace_tags(self, filename, var_to_set):
        f = open(filename, 'r')
        lines = []
        for line in f:
            reg = re.compile(r'^.*?(\(#-[^(#-]+-#\)).*?$')
            reg_m = reg.match(line)
            if reg_m is not None:
                tag   = reg_m.groups()[0]
                value = var_to_set[tag]
                line = line.replace(tag, value)
            lines.append(line)
        f.close()

        f = open(filename, 'w')
        f.writelines(lines)
        f.close()


    @property
    def start(self):
        """
        Start the job.
        If the job has been started previously, it will be continued.
        If the job has crashed or did not finish, the last action will be restarted.
        """
        #-# print('Starting job ' + self.jobname)


        # Try to read the results from a previous run.
        # Some settings are allowed to be changed. (Karlsberg+ run folder)
        kbp_folder = self.kbp_folder
        apbs_subfolder_suffix = self.apbs_subfolder_suffix
        md_template_folder = self.md_template_folder
        pdb = self.pdb
        i_min = self.min
        ref0 = self.ref0
        vdw = self.vdw
        cavity_parameter = self.cavity_parameter
        water_template = self.water_template
        recovered = self.recover()
        self.kbp_folder = kbp_folder
        self.apbs_subfolder_suffix = apbs_subfolder_suffix
        self.md_template_folder = md_template_folder
        self.pdb = pdb
        self.min = i_min
        self.ref0 = ref0
        self.vdw = vdw
        self.cavity_parameter = cavity_parameter
        self.water_template = water_template




        if not recovered:
            self.status = 'running'
            self.log.append('Job %s started.' % self.jobname)
            #if os.path.exists(self.root_workdir):
            #    shutil.rmtree(self.root_workdir)
            #os.mkdir(self.root_workdir)
        else:
            self.log.append('Job %s continued.' % self.jobname)

        # This vairable was introduced later, so in some older jobs it does not exist.
        self.kbp_folder_suffix = self.kbp_folder.strip('/')[self.kbp_folder.strip('/').rfind('/')+1:]

        # Check the folder structure and reset a job if the corresponding folder is missing.
        if not os.path.exists(self.modelling_folder):
            print("Folder '" + self.modelling_folder + "' not found: resetting modelling status.")
            self.modelling_status = 'none'
            self.md_status        = 'none'
            self.kbp_status       = 'none'
        if not os.path.exists(self.md_folder):
            print("Folder '" + self.md_folder + "' not found: resetting MD status.")
            self.md_status  = 'none'
            self.md_problem  = ''
            self.md_progress = ''
            self.kbp_status = 'none'

        else:
            job_completed = False
            files_required = ['_md.dcd', '_md_min_e1.dcd', '_md_min_e4.dcd']
            for prefix in files_required:
                if not os.path.exists(self.md_folder + self.jobname + prefix):
                    break
            else:
                job_completed = True
            # if job_completed:
            #     self.md_status  = 'done'
            #     self.md_problem  = ''
            #     self.md_progress = ''
            # else:
            if not job_completed:
                self.md_status  = 'running'
                self.md_problem  = ''
                self.md_progress = ''
        if not os.path.exists(self.kbp_folder):
            print("Folder '" + self.kbp_folder + "' not found: resetting Karlsberg+ status.")
            self.kbp_status   = 'none'
            self.kbp_problem  = ''
            self.kbp_progress = ''
        if not os.path.exists(self.kbp_folder + 'done') or not os.path.exists(self.kbp_folder + 'kbp2_jobs'):
            print("Folder '" + self.kbp_folder + "' incomplete: resetting Karlsberg+ status.")
            self.kbp_status   = 'none'
            self.kbp_problem  = ''
            self.kbp_progress = ''


        if not os.path.exists(self.kbp_folder + 'confE%s/' % self.apbs_subfolder_suffix):
            print("Folder '" + self.kbp_folder + 'confE%s/' % self.apbs_subfolder_suffix + "' not found: resetting APBS status.")
            self.apbs_status = 'none'
            self.apbs_problem   = ''
            self.apbs_progress  = ''

        # Debug
        # self.kbp_queue = '/scratch/scratch/tmeyer/md_pka/md_pka_manager/queue/kbp_jobs/'
        # self.kbp_status = 'running'
        # self.md_status = 'running'
        # self.md_status = ''
        # self.kbp_status = ''
        # self.apbs_status = ''

        ####################
        ### 1) Modelling ###
        ####################
        if self.modelling_status != 'done':
            self.log.append('Checking the structure.')

            #################################
            ### Make decisions (optional) ###
            #################################
            # cm.add_decision('rename__CA_CAL', 'keep')
            # cm.add_decision('rename__HOH_TIP3', 'rename')
            # cm.add_decision('ligand_collisions', 'delete')

            # ToDo: Check if the modelling was successful
            if os.path.exists(self.modelling_folder):
                done = 1
            else:
                done = 0

            if done == 1:
                self.modelling_status = 'done'
            else:
                self.modelling_status = 'crashed'

            # Store what is important from Charmm_manager object and delete it, since it prevents any pickle dump.
            self.charmm_manager_data['prefix'] = self.modelling_prefix
            self.charmm_manager_data['top']    = self.top
            self.charmm_manager_data['par']    = self.par
            self.charmm_manager = None

        self.log.append('Modelling finished with status: %s' % self.modelling_status)

        self.dump()

        if self.modelling_status != 'done':
            reason = 'charmm_crash'
            self.status = 'crashed'
            self.dump()
            return reason

        if self.status == 'stopping':
            self.status = 'stopped'
            self.dump()
            return 'stop_requested'



        #############
        ### 2) MD ###
        #############
        # self.md_status = 'running'
        if self.md_status != 'done' and self.md_status != 'running':
            self.log.append('Preaparing MD external_scripts.')

            self.md_status  = 'preparing'
            self.md_problem = ''


            if not os.path.exists(self.md_template_folder):
                self.md_problem = "The folder containing the MD template external_scripts is not reachable: " \
                                  + self.md_template_folder
                self.md_status = 'crashed'
                self.status = 'crashed'
                self.dump()
                reason = 'md_crash'
                return reason

            # if os.path.exists(self.md_folder):
            #     # Restart the job.
            #     shutil.rmtree(self.md_folder)

            # Prepare folder
            if not os.path.exists(self.md_folder):
                os.mkdir(self.md_folder)

                # md_cache_folder = "/scratch/scratch/tmeyer/md_pka/runs/md_cache/"
                # files_to_copy = ['1_charmm_solvate.out', '2_namd_min.out', '3_namd_eq_heat_run.out', '4_namd_main.out',
                #                  'protein_in_water.xplor.psf', 'protein_in_water.pdb', 'protein_in_water.psf',
                #                  'protein_in_water_min_neutral.xplor.psf', 'protein_in_water_min_neutral.psf',
                #                  'protein_in_water_min_neutral.pdb',
                #                  'equ_constrains.pdb',
                #                  '%s_md.dcd' % self.jobname,
                #                  '%s_md.restart.coor' % self.jobname, '%s_md.restart.xsc' % self.jobname,
                #                  '%s_md.restart.vel' % self.jobname
                #                  ]
                #'%s_md_min_e1.dcd' % self.jobname, '%s_md_min_e4.dcd' % self.jobname

            # if os.path.exists(md_cache_folder + self.jobname):
            #     print("Cache found, copying files.")
            #     for file_to_copy in files_to_copy:
            #         src = md_cache_folder + self.jobname + '/md/' + file_to_copy
            #         dst = self.md_folder  + file_to_copy
            #         if os.path.exists(src):
            #             shutil.copy(src, dst)


            # Copy required files from modelling folder
            psf_src_file = self.modelling_folder + self.charmm_manager_data['prefix'] + '.psf'
            psf_input_file = self.md_folder      + self.charmm_manager_data['prefix'] + '.psf'
            shutil.copy2(psf_src_file, psf_input_file)
            crd_src_file = self.modelling_folder + self.charmm_manager_data['prefix'] + '.crd'
            crd_input_file = self.md_folder      + self.charmm_manager_data['prefix'] + '.crd'
            shutil.copy2(crd_src_file, crd_input_file)

            # Copy template external_scripts
            for fn in os.listdir(self.md_template_folder):
                if fn[0] != '.':
                    src = self.md_template_folder + fn
                    dst = self.md_folder          + fn
                    shutil.copy(src, dst)

            # Change variables in template external_scripts
            top_par_charmm_commands = ''
            par_namd_commands       = ''
            first = True
            for top in self.charmm_manager_data['top']:
                if first:
                    top_par_charmm_commands += "read rtf card name \""  + top + '\"\n'
                    first = False
                else:
                    top_par_charmm_commands += "read rtf card name \""  + top + '\" append\n'
            first = True
            for par in self.charmm_manager_data['par']:
                if first:
                    top_par_charmm_commands += "read para card name \"" + par + '\"\n'
                    par_namd_commands       += "parameters          "   + par + "\n"
                    first = False
                else:
                    top_par_charmm_commands += "read para card name \"" + par + '\" append\n'
                    par_namd_commands       += "parameters          "   + par + "\n"

            var_to_set = {}
            # 1_charmm_solvate.inp, Neutral.inp
            var_to_set['(#-toppar-#)']         = top_par_charmm_commands
            #                var_to_set['(#-inputname_psf-#)']  = psf_src_file
            #                var_to_set['(#-inputname_crd-#)']  = crd_src_file
            var_to_set['(#-inputname_psf-#)']  = self.charmm_manager_data['prefix'] + '.psf'
            var_to_set['(#-inputname_crd-#)']  = self.charmm_manager_data['prefix'] + '.crd'
            var_to_set['(#-waterboxoffset-#)'] = '10.0'
            # 2_namd_min.inp, 3_namd_eq_heat_run.inp get_frames.tcl
            var_to_set['(#-par_namd-#)']       = par_namd_commands
            md_jobname = self.jobname + '_md'
            var_to_set['(#-outputname-#)'] = md_jobname
            # prep.sh get_frames.tcl
            var_to_set['(#-path-#)']    = self.md_folder
            var_to_set['(#-jobname-#)'] = self.jobname
            # get_frames.tcl
            #var_to_set['(#-frame_folder-#)'] = self.kbp_folder + 'frames'
            #var_to_set['(#-water_folder-#)'] = self.kbp_folder + 'water'


            files_to_search = ['1_charmm_solvate.inp', '2_namd_min.inp', '3_namd_eq_heat_run.inp', \
                               'run_cuda.sh', 'run_parallel.sh', 'run_parallel_min.sh', 'Neutralize.inp', \
                               '4_namd_main.inp', '5_namd_minimize_e4.inp', '5_namd_minimize_e1.inp']
            # '5_namd_minimize_e4_3000.inp', '5_namd_minimize_e1_3000.inp']

            for fn in files_to_search:
                self.replace_tags(self.md_folder + fn, var_to_set)




            self.dump()

            # Submit jobs
            script =  "cd %s" % self.md_folder + "\n"
            script += "./prep.sh > prep.out \n"

            j = {'script' : script}
            #                queue_dir = '/scratch/scratch/tmeyer/md_pka/md_pka_manager/queue/d58/'
            number = ''
            done = False
            while not done:
                job_filename = "job_%s.pickle" % (self.jobname + number)
                if not os.path.exists(self.md_queue + job_filename):
                    f = open(self.md_queue + job_filename, 'w')
                    pickle.dump(j,f)
                    f.close()
                    done = True
                else:
                    if number == '':
                        number = '1'
                    else:
                        number = str(int(number) + 1)

            self.md_status = 'running'
            self.md_problem  = ''
            self.md_progress = 'Structure preparation started..'

            self.log.append('MD jobs submitted.')

            self.dump()


        if self.md_status != 'done':
            # Enter sleeping mode and wait for job to finish or crash.
            self.log.append('Waiting for MD jobs to finish.')
            if os.path.exists(self.md_folder + '4_namd_main.inp'):
                namd_output_filename = self.md_folder + '4_namd_main.out'
            else:
                namd_output_filename = self.md_folder + '3_namd_eq_heat_run.out'

            post_processing_output_filenames = [self.md_folder + '5_namd_minimize_e1.out',
                                                self.md_folder + '5_namd_minimize_e4.out']
            # self.md_folder + '5_namd_minimize_e4_3000.out',
            # self.md_folder + '5_namd_minimize_e4_3000.out']

            first_loop = True
            sleep = True
            main_run_finished = False
            while sleep:
                if first_loop:
                    first_loop = False
                else:
                    time.sleep(900)

                if not os.path.exists(self.md_folder):
                    self.md_status   = 'crashed'
                    self.md_problem  = "MD folder has not been created!"
                    self.md_progress = 'crashed'
                    break

                post_procesing_running = False
                for post_processing_output_filename in post_processing_output_filenames:
                    if os.path.exists(post_processing_output_filename):
                        post_procesing_running = True
                        break

                if not os.path.exists(namd_output_filename):
                    self.md_progress = 'queued'
                    if os.path.exists(self.md_folder + '1_charmm_solvate.out'):
                        self.md_progress = 'solvating protein'
                    if os.path.exists(self.md_folder + 'output_min/2_namd_min.out'):
                        self.md_progress = 'minimizing water'
                    if os.path.exists(self.md_folder + 'Neutralize.out'):
                        self.md_progress = 'adding ions'
                    if os.path.exists(self.md_folder + '3_namd_eq_heat_run.out'):
                        self.md_progress = 'heating up'
                    self.md_progress = 'Structure preparation pending.. (%s)' % self.md_progress
                if os.path.exists(namd_output_filename):
                    f = open(namd_output_filename, 'r')
                    lines = f.readlines()
                    f.close()

                    self.md_progress = "Main run started. (energy minimization)"

                    speed = 0.0
                    for line in lines[-1::-1]:
                        # Info: Initial time: 4 CPUs 0.0114801 s/step 0.0664358 days/ns 85.625 MB memory
                        if line.find("Info: Initial time:") > -1:
                            speed = float(line.split()[7])
                            speed = 1.0 / speed
                            break

                    for line in lines[-1::-1]:
                        if line.find("Program finished.") > -1:
                            #self.md_status   = 'done'
                            self.md_problem  = ''
                            #self.md_progress = 'done'
                            self.md_progress = 'Main run done.'
                            main_run_finished = True
                            #sleep = False
                            break
                        if line.find("FATAL ERROR") > -1:
                            self.md_status   = 'crashed'
                            self.md_problem  = line
                            self.md_progress = 'crashed'
                            sleep = False
                            break
                        if line[0:7].find('ENERGY:') > -1:
                            entries = line.split()
                            timestep = entries[1]
                            self.md_progress = "Main run at %.2f ns (%.1f ns/day)." \
                                               % (float(timestep) * 2.0 / 1e6, speed)
                            break
                #else:
                if post_procesing_running and main_run_finished:
                    #self.md_problem  = ''
                    self.md_progress = 'minimization running'

                    done = 0
                    for post_processing_output_filename in post_processing_output_filenames:
                        if not os.path.exists(post_processing_output_filename):
                            continue

                        f = open(post_processing_output_filename, 'r')
                        last_lines = tail(f)
                        f.close()

                        #f = open(post_processing_output_filename, 'w')
                        #last_lines = ['Manually repared:', 'Program finished.']
                        #for l in last_lines:
                        #    f.write(l + '\n')
                        #f.close()

                        for line in last_lines[-1::-1]:
                            if (line.find("Program finished.") > -1) or (line.find("WallClock:") > -1):
                                done += 1
                                # The log files are truncated here, since they are way too large.
                                f = open(post_processing_output_filename, 'w')
                                for l in last_lines:
                                    f.write(l + '\n')
                                f.close()
                                break
                            if line.find("FATAL ERROR") > -1:
                                self.md_status   = 'crashed'
                                self.md_problem  = line
                                self.md_progress = 'minimization crashed'
                                sleep = False
                                break

                    if done >= len(post_processing_output_filenames):
                        self.md_status   = 'done'
                        self.md_problem  = ''
                        self.md_progress = 'done'
                        sleep = False
                    else:
                        self.md_progress = 'Minimization step %i of %i.' % (done+1, len(post_processing_output_filenames))

                self.dump()

                stop = self.check_queue()
                if stop:
                    sleep = False
                    self.md_status   = 'stopped'
                    self.md_problem  = 'stop requested'
                    self.status = 'stopped'
                    self.dump()
                    reason = 'stop_requested'
                    return reason


        self.log.append('MD jobs done with status: ' + self.md_status)

        self.dump()

        if self.md_status != 'done':
            reason = 'md_crash'
            self.status = 'crashed'
            self.dump()
            return reason

        if self.status == 'stopping':
            reason = 'stop_requested'
            self.status = 'stopped'
            self.dump()
            return reason

        # self.kbp_status = ''
        # self.kbp_progress = ''
        # self.kbp_problem = ''

        #####################
        ### 3) Karlsberg+ ###
        #####################
        if self.kbp_status != 'done':
            self.log.append('Preparing Karlsberg+ jobs.')
            self.kbp_status  = 'preparing'
            self.kbp_problem = ''

            # Create folder for Karlsberg+ jobs if necessary.
            if not os.path.exists(self.kbp_folder):
                os.mkdir(self.kbp_folder)


            ### Extract frames if necessary. ###
            frame_folder = self.kbp_folder + 'frames/'
            water_folder = self.kbp_folder + 'water/'
            if not os.path.exists(frame_folder):
                os.mkdir(frame_folder)
                if not os.path.exists(water_folder):
                    os.mkdir(water_folder)

                if self.min == 0:
                    if not self.water_template:
                        shutil.copy2(self.md_template_folder + 'get_frames_now.tcl', self.kbp_folder + 'get_frames.tcl')
                    else:
                        shutil.copy2(self.md_template_folder + 'get_frames.tcl', self.kbp_folder + 'get_frames.tcl')
                elif self.min == 1:
                    if not water_template:
                        shutil.copy2(self.md_template_folder + 'get_frames_min_e1_now.tcl', self.kbp_folder + 'get_frames.tcl')
                    else:
                        shutil.copy2(self.md_template_folder + 'get_frames_min_e1.tcl', self.kbp_folder + 'get_frames.tcl')
                elif self.min == 4:
                    if not self.water_template:
                        shutil.copy2(self.md_template_folder + 'get_frames_min_e4_now.tcl', self.kbp_folder + 'get_frames.tcl')
                    else:
                        shutil.copy2(self.md_template_folder + 'get_frames_min_e4.tcl', self.kbp_folder + 'get_frames.tcl')



                # shutil.copy2(self.md_template_folder + 'get_frames_min_e4_now_3000min.tcl', self.kbp_folder + 'get_frames.tcl')
                # shutil.copy2(self.md_template_folder + 'get_frames_min_e1_now_3000min.tcl', self.kbp_folder + 'get_frames.tcl')

                # shutil.copy2(self.md_template_folder + 'get_frames_min_e4_nearw.tcl', self.kbp_folder + 'get_frames.tcl')
                # shutil.copy2(self.md_template_folder + 'get_frames.tcl', self.kbp_folder + 'get_frames.tcl')
                # shutil.copy2(self.md_template_folder + 'get_frames_min_e4.tcl', self.kbp_folder + 'get_frames.tcl')
                # shutil.copy2(self.md_template_folder + 'get_frames_min_e1.tcl', self.kbp_folder + 'get_frames.tcl')

                shutil.copy2(self.md_template_folder + 'get_frames.sh',  self.kbp_folder + 'get_frames.sh')

                # Change variables in template external_scripts
                var_to_set = {}
                md_jobname = self.jobname + '_md'
                var_to_set['(#-outputname-#)'] = md_jobname
                var_to_set['(#-path-#)']    = self.md_folder
                var_to_set['(#-frame_folder-#)'] = self.kbp_folder + 'frames'
                var_to_set['(#-water_folder-#)'] = self.kbp_folder + 'water'

                files_to_search = ['get_frames.tcl']

                for fn in files_to_search:
                    self.replace_tags(self.kbp_folder + fn, var_to_set)

                # Start get_frames.sh
                script = "cd %s\n./get_frames.sh\nexit\n" % self.kbp_folder
                process = subprocess.Popen("tcsh", \
                                           shell=True, \
                                           stdin=subprocess.PIPE, \
                                           stdout=subprocess.PIPE, \
                                           stderr=subprocess.PIPE, \
                                           )
                process.stdin.write(script)
                while True:
                    next_line = process.stdout.readline()
                    if not next_line:
                        break

                errs = process.stderr.readlines()
                #if errs:
                #    print self.jobname
                #    for err in errs:
                #        print err

                f = open(self.kbp_folder + 'get_frames.out', 'r')
                for line in f:
                    if line.find('ERROR') > -1:
                        if line.find('Stride'):
                            continue
                        self.status       = 'crashed'
                        self.kbp_status   = 'crashed'
                        self.kbp_problem  = "get_frames.tcl contains errors:\n%s" % (line)
                        self.kbp_progress = 'crashed'
                        self.dump()
                        f.close()
                        reason = 'kbp_crash'
                        return reason
                f.close()


            # Get list of pdb files.
            fd = os.listdir(frame_folder)
            pdb_list = []
            for entry in fd:
                if os.path.isfile( frame_folder + entry ):
                    if entry.find('~') > -1 or entry[0] == '.' or entry[-4:] != '.pdb':
                        continue
                    pdb_list.append(entry)


            # Prepare missing Karlsberg+ jobs and sleep until all runs have finished.
            self.log.append('Waiting for Karlsberg+ jobs to finish.')

            first_loop = True
            sleep = True
            while sleep:
                reset_crashed = False
                if first_loop:
                    first_loop = False
                    reset_crashed = True
                else:
                    # time.sleep(3600)
                    time.sleep(1800)

                reset_crashed_file = self.kbp_folder + 'restart_crashed'
                if os.path.exists(reset_crashed_file):
                    reset_crashed = True
                    # The crashed files are reset only once.
                    os.remove(reset_crashed_file)
                    self.log.append('Found a reset file for crashed Karlsberg+ jobs in:\n' + self.jobname)
                    self.kbp_status  = 'running'
                    self.kbp_problem = ''


                queued  = 0
                running = 0
                done    = 0
                crashed = 0

                # Check for existing jobs.
                pdb_list_submit = []
                for pdb in pdb_list:

                    run_dir  = self.kbp_folder + 'kbp2_jobs/'
                    done_dir = self.kbp_folder + 'done/'
                    pdb_job_folder = pdb[:-4] + '/'

                    # Did the job finish?
                    if os.path.exists(done_dir + pdb_job_folder):
                        # Was it successful?
                        success = False
                        result_file = done_dir + pdb_job_folder + "pka_results.dat"
                        if os.path.exists(result_file):
                            f = open(result_file)
                            first_line = f.readline()
                            f.close()
                            if first_line.find('Conformations:') == 0:
                                # Job finished successfully.
                                done += 1
                                success = True
                        if not success:
                            # Reset the job, if a file named "restart_crashed" is present in the kbp folder.
                            if reset_crashed:
                                shutil.rmtree(done_dir + pdb_job_folder)
                                # Submit the job.
                                pdb_list_submit.append(pdb)
                            else:
                                # Job crashed.
                                crashed += 1

                    # Is the job still running?
                    elif os.path.exists(run_dir + pdb_job_folder):
                        if os.path.isdir(run_dir + pdb_job_folder):

                            job_filename = self.kbp_queue + "kbp_job_%s_%s_%s.pickle" \
                                                            % (self.kbp_folder_suffix, self.jobname, pdb[:-4])
                            # Is the job still in queue folder?
                            if os.path.exists(job_filename):
                                running += 1
                            else:
                                # Reset and submit the job.
                                print("The file %s doeas not exist in the queue folder. Job will be reset." % job_filename)
                                shutil.rmtree(run_dir + pdb_job_folder)
                                pdb_list_submit.append(pdb)

                    else:
                        if os.path.isfile(run_dir + pdb_job_folder[:-1]):
                            # This was a bug and should not happen anymore.
                            print(run_dir + pdb_job_folder[:-1] + " is a file! Will be deleted.")
                            os.remove(run_dir + pdb_job_folder[:-1])

                        # Otherwise submit the job.
                        pdb_list_submit.append(pdb)


                if len(pdb_list_submit) > 0:

                    # Prepare Karlsberg+ template job.
                    kbp_template = job_manager2.kbp_job()

                    #            kbp_template.kb_para['sfc_globals']['rtfPatches'] = '/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.rtf'
                    #            kbp_template.kb_para['sfc_globals']['prmPatches'] = '/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/patches.prm'
                    #            kbp_template.kb_para['sfc_globals']['titratable'] = '/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/titratable.yaml'
                    #             kbp_template.kb_para['globals']['include'] = '/scratch/scratch/tmeyer/kbplus2_epsp1'
                    #             kbp_template.kb_para['globals']['include'] = '/scratch/scratch/tmeyer/kbplus2_epsw120'


                    kbp_template.kb_para['sfc'][0]['- name'] = "c_pH7"
                    kbp_template.kb_para['sfc'][0]['  pH']   = self.ph

                    # kbp_template.kb_para['sfc'][0]['unnamed_entries'][0] = \
                    #         "   - 'what=((RESNAME EPP .OR. RESNAME DPP) .AND. .NOT. BB) .OR. HYDROGENS,confs=0,kfc=0'"

                    # sfc = dict(kbp_template.kb_para['sfc'][0])
                    # sfc['- name'] = 'c_pH-10'
                    # sfc['  pH'] = -10
                    # kbp_template.kb_para['sfc'].append(sfc)
                    #
                    # sfc = dict(kbp_template.kb_para['sfc'][0])
                    # sfc['- name'] = 'c_pH20'
                    # sfc['  pH'] = 20
                    # kbp_template.kb_para['sfc'].append(sfc)


                    # Set up cavity calculation:
                    if self.cavity_parameter is not None:
                        kbp_template.parameter['cavity_parameter'] = 'cavity %.2f 0.25 0.0' % self.cavity_parameter
                        if self.water_template:
                            kbp_template.parameter['water_folder'] = self.kbp_folder + 'water/'

                    # Use an alternative version of kbp2:
                    # kbp_template.kb_para['sfc_globals']['tapbs'] = 'tapbs_1.3_enere'

                    # This one will create .ener files
                    kbp_template.kb_para['sfc_globals']['tapbs'] = 'tapbs_1.3_cav_enere'
                    # This one will use ion concentration of 0.0 for the model calculations
                    # kbp_template.kb_para['sfc_globals']['tapbs'] = 'tapbs_1.3_noionsref'

                    ###
                    ### Rashin 86 vdW radii
                    ###
                    # (ASP and GLU)
                    # kbp_template.kb_para['sfc_globals']['prm'] = '/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/par_rashin.inp'
                    # all
                    if self.vdw:
                        kbp_template.kb_para['sfc_globals']['prm'] = '/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/par_rashin2.inp'
                        # kbp_template.kb_para['sfc_globals']['prm'] = '/scratch/scratch/tmeyer/CHARMM_NAMD/toppar/md_mod/par_acids.inp'

                    if self.ref0:
                        # titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted.yaml'
                        # titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_nter2.yaml'
                        # titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter2.yaml'

                        titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter2_np.yaml'
                        # titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter2_np_propka.yaml'
                        # titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter2_np_fit.yaml'

                        kbp_template.kb_para['sfc_globals']['titratable'] = titratable_yaml
                    else:
                        # titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_ter.yaml'
                        pass

                        # titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_sidechains.yaml'
                        # titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter.yaml'

                        # kbp_template.kb_para['sfc_globals']['titratable'] = titratable_yaml


                    #kbp_template.kb_para['globals']['include'] = '/scratch/scratch/tmeyer/kbplus2_new/'
                    #kbp_template.parameter['exec'] = "perl -I/scratch/scratch/tmeyer/kbplus2_new /scratch/scratch/tmeyer/kbplus2_new/titrate.pl"

                    #kbp_template.kb_para['globals']['include'] = '/scratch/scratch/tmeyer/kbplus2_epsp6/'
                    #kbp_template.parameter['exec'] = "perl -I/scratch/scratch/tmeyer/kbplus2_epsp6 /scratch/scratch/tmeyer/kbplus2_epsp6/titrate.pl"
                    # kbp_template.kb_para['globals']['include'] = '/scratch/scratch/tmeyer/kbplus2_epsp2/'
                    # kbp_template.parameter['exec'] = "perl -I/scratch/scratch/tmeyer/kbplus2_epsp2 /scratch/scratch/tmeyer/kbplus2_epsp2/titrate.pl"
                    #kbp_template.kb_para['globals']['include'] = '/scratch/scratch/tmeyer/kbplus2_noMin/'
                    #kbp_template.parameter['exec'] = "perl -I/scratch/scratch/tmeyer/kbplus2_noMin /scratch/scratch/tmeyer/kbplus2_noMin/titrate.pl"



                    # -> Use one of these
                    # -> termini definition fixed in this one
                    kbp_template.kb_para['globals']['include'] = '/scratch/scratch/tmeyer/kbplus2_noMin/'
                    kbp_template.parameter['exec'] = "perl -I/scratch/scratch/tmeyer/kbplus2_noMin /scratch/scratch/tmeyer/kbplus2_noMin/titrate.pl"

                    # This one is for protein epsilon = 2
                    # kbp_template.kb_para['globals']['include'] = '/scratch/scratch/tmeyer/kbplus2_pe2_noMin_new/'
                    # kbp_template.parameter['exec'] = "perl -I/scratch/scratch/tmeyer/kbplus2_pe2_noMin_new /scratch/scratch/tmeyer/kbplus2_pe2_noMin_new/titrate.pl"

                    # This one is for protein epsilon = 8
                    # kbp_template.kb_para['globals']['include'] = '/scratch/scratch/tmeyer/kbplus2_pe8_noMin_new/'
                    # kbp_template.parameter['exec'] = "perl -I/scratch/scratch/tmeyer/kbplus2_pe8_noMin_new /scratch/scratch/tmeyer/kbplus2_pe8_noMin_new/titrate.pl"



                    kbp_template.kb_para['globals']['workDir'] = "/public/scratch/tmeyer/kb/%s/%s/" \
                                                                 % (self.kbp_folder_suffix, self.jobname)

                    #kbp_template.parameter['pdb_dir']  = frame_folder
                    kbp_template.parameter['pdb_dir']  = self.kbp_folder + 'frames/'
                    kbp_template.parameter['run_dir']  = self.kbp_folder + 'kbp2_jobs/'
                    kbp_template.parameter['done_dir'] = self.kbp_folder + 'done/'

                    if not os.path.exists(kbp_template.parameter['run_dir']):
                        os.mkdir(kbp_template.parameter['run_dir'])
                    if not os.path.exists(kbp_template.parameter['done_dir']):
                        os.mkdir(kbp_template.parameter['done_dir'])



                    # Replace all queue entries, if a file named "reset_queue" is present in the kbp folder.
                    reset_queue = False
                    reset_queue_file = self.kbp_folder + 'reset_queue'
                    if os.path.exists(reset_queue_file):
                        reset_queue = True

                    # Prepare the jobs.
                    kbp_jobs_prepared_list = []
                    for pdb in pdb_list_submit:
                        job_filename = self.kbp_queue +  "kbp_job_%s_%s_%s.pickle" \
                                                         % (self.kbp_folder_suffix, self.jobname, pdb[:-4])

                        job_filename_running = job_filename + ' - running'
                        job_filename_done    = job_filename + ' - done'
                        job_filename_result  = job_filename + ' - result'

                        submit = False

                        # Is the job running?
                        # if os.path.exists(job_filename_running):
                        #     pass
                        # This should not happen! If it is running there is a folder, that was detected above.
                        # if os.path.exists(job_filename_done):
                        #     # This is a problem. But since the Karlsberg+ folder do not exist, the job will
                        #     # crash anyway. Clean up, in case the listener had crashed.
                        #     os.remove(job_filename)
                        #     os.remove(job_filename_running)
                        #     crashed += 1
                        # else:
                        #     # Is counted above.
                        #     pass

                        # Is the job already queued?
                        if os.path.exists(job_filename):
                            if reset_queue:
                                # Clean up for a restart.
                                os.remove(job_filename)
                                if os.path.exists(job_filename_done):
                                    os.remove(job_filename_done)
                                if os.path.exists(job_filename_result):
                                    os.remove(job_filename_result)
                                # In case the listener crashed.
                                if os.path.exists(job_filename_running):
                                    os.remove(job_filename_running)
                                submit = True
                            else:
                                if os.path.exists(job_filename_done):
                                    done += 1
                                else:
                                    queued += 1
                        else:
                            submit = True
                        if submit:
                            # Submit the job.
                            job = job_manager2.kbp_job(pdb, kbp_template)
                            kbp_jobs_prepared_list.append(job)

                            f = open(job_filename, 'w')
                            pickle.dump(job, f)
                            f.close()
                            queued += 1

                    # Debug
                    sum = done + crashed + running + queued
                    if sum != len(pdb_list):
                        print("Oh, thats bad. %i jobs got lost!" % (sum - len(pdb_list),))

                if (done + crashed) == len(pdb_list):
                    sleep = False
                    self.kbp_status   = 'running'
                    self.kbp_problem  = ''
                    self.kbp_progress = 'checking Karlsberg+ jobs'
                else:
                    self.kbp_status   = 'running'
                    self.kbp_problem  = ''
                    self.kbp_progress = "Waiting for Karlsberg+ jobs to finish: %i/%i/%i/%i (done/crashed/running/queued)" \
                                        % (done, crashed, running, queued)

                stop = self.check_queue()
                if stop:
                    sleep = False
                    self.md_status   = 'stopped'
                    self.md_problem  = 'stop requested'
                    self.status = 'stopped'
                    self.dump()
                    reason = 'stop_requested'
                    return reason


            if crashed > 0:
                self.status       = 'crashed'
                self.kbp_problem  = "Some jobs crashed: %i of %i" % (crashed, len(self.kbp_jobs_prepared))
                self.kbp_progress = 'crashed'
            else:
                self.kbp_status   = 'done'
                self.kbp_problem  = ''
                self.kbp_progress = 'done'

        self.log.append('Karlsberg+ jobs done with status: ' + self.md_status)

        self.dump()

        if self.kbp_status != 'done':
            reason = 'kbp_crash'
            self.status = 'crashed'
            self.dump()
            return reason

        if self.status == 'stopping':
            reason = 'stop_requested'
            self.status = 'stopped'
            self.dump()
            return reason



        ###############
        ### 4) APBS ###
        ###############
        # self.apbs_status = 'none'
        if self.apbs_status != 'done':
            # self.apbs_status = 'done'
            # self.apbs_problem = ''
            # if False:
            self.log.append('Preparing APBS jobs.')
            self.apbs_status  = 'preparing'

            # Jobs are started in "confE<apbs_subfolder_suffix>/" subfolder of self.kbp_folder
            subfolder_suffix = self.apbs_subfolder_suffix

            # if self.apbs_problem == 'crashed':
            (unfinished_jobs, crashed_jobs, finished_jobs) = apbs_manager.read_results(self.kbp_folder, subfolder_suffix=subfolder_suffix)
            if len(crashed_jobs) > 0:
                apbs_manager.clean_crashed_jobs(self.kbp_folder, '', crashed_jobs, subfolder_suffix=subfolder_suffix)
            # if len(unfinished_jobs) > 0:
            #    apbs_manager.clean_crashed_jobs(self.kbp_folder, '', unfinished_jobs, subfolder_suffix=subfolder_suffix)

            self.apbs_problem = ''

            apbs_bin    = "/scratch/scratch/tmeyer/CHARMM_NAMD/apbs"
            coulomb_bin = "/scratch/scratch/tmeyer/CHARMM_NAMD/coulomb"

            dolly_run_folder = "/public/scratch/tmeyer/kb/apbs/" + self.kbp_folder_suffix + '/'

            # Use CAVITIES and water template?
            # apbs_manager.submit_jobs(self.kbp_folder, '', apbs_bin, coulomb_bin, restart=False, target_res=0.3,\
            #                          dolly_prefix=self.jobname, use_water_template=False,
            #                          subfolder_suffix=subfolder_suffix, dolly_run_folder=dolly_run_folder)

            apbs_manager.submit_jobs(self.kbp_folder, '', apbs_bin, coulomb_bin, restart=False, target_res=0.3, \
                                     dolly_prefix=self.jobname,
                                     use_water_template=water_template, cavity_parameter=cavity_parameter,
                                     subfolder_suffix=subfolder_suffix, dolly_run_folder=dolly_run_folder)

            # apbs_manager.submit_jobs(self.kbp_folder, '', apbs_bin, coulomb_bin, restart=False, target_res=0.3,\
            #                         dolly_prefix=self.jobname, use_water_template=True,
            #                         subfolder_suffix=subfolder_suffix, dolly_run_folder=dolly_run_folder,
            #                         cavity_parameter=0.8)

            self.apbs_status = 'running'

            first_loop = True
            sleep = True
            while sleep:
                if first_loop:
                    first_loop = False
                else:
                    time.sleep(3600)

                (unfinished_jobs, crashed_jobs, finished_jobs) = apbs_manager.read_results(self.kbp_folder, subfolder_suffix=subfolder_suffix)

                self.apbs_progress = "Waiting for APBS jobs to finish: %i/%i/%i (done/crashed/unfinished)" \
                                     % (len(finished_jobs), len(crashed_jobs), len(unfinished_jobs))

                if len(unfinished_jobs) == 0:
                    sleep = False

                stop = self.check_queue()
                if stop:
                    sleep = False
                    self.apbs_status   = 'stopped'
                    self.apbs_problem  = 'stop requested'
                    self.status = 'stopped'
                    self.dump()
                    reason = 'stop_requested'
                    return reason


            if len(crashed_jobs) > 0:
                self.apbs_status = 'crashed'
                self.apbs_problem = '%i jobs crashed.' % len(crashed_jobs)
            else:
                self.apbs_status = 'done'
                self.apbs_problem = ''


        if self.apbs_status != 'done':
            reason = 'apbs_crash'
            self.status = 'crashed'
            self.dump()
            return reason

        ##########################################
        ### Start CHARMM conformational energy ###
        ##########################################
        # Get the Karlsberg+ definition of states
        if not self.ref0:
            titratable_yaml = '/scratch/scratch/tmeyer/kbplus2/titratable.yaml'
        else:
            # titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter2.yaml'
            titratable_yaml = '/scratch/scratch/tmeyer/md_pka/runs/titratables/titratable_refsmall_shifted_ter2_np.yaml'
        titratable_residues = kbp_tools.parse_titratable_yaml(titratable_yaml)
        template_state = {'charge' : 0,
                          'patch' : None,
                          'external_patches' : None,
                          'rename' : None,
                          'special' : None}
        titratable_residues_kbp = {}
        for restype in titratable_residues:
            titratable_residues_kbp[restype] = []
            for state in titratable_residues[restype]:
                new_state = dict(template_state)
                if 'patch' in state:
                    new_state['patch'] = state['patch']
                titratable_residues_kbp[restype].append(new_state)

        # charmm_conf_folder = self.kbp_folder + 'charmm_confE/'
        charmm_conf_folder = self.kbp_folder + 'charmm_confE_no14/'

        # charmm_conf_folder = self.kbp_folder + 'charmm_confE_gbsw/'
        if not os.path.exists(charmm_conf_folder):
            os.mkdir(charmm_conf_folder)

        # Get the charmm object
        # f = open(self.modelling_folder + 'charmm_struct.pkl', 'rb')
        # charmm_struct = pickle.load(f)
        # f.close()
        # assert(isinstance(charmm_struct, charmm.Charmm_manager))
        top = ['/scratch/scratch/tmeyer/karlsbergplus/top.inp',
               '/scratch/scratch/tmeyer/karlsbergplus/patches.rtf']
        par = ['/scratch/scratch/tmeyer/karlsbergplus/par.inp',
               '/scratch/scratch/tmeyer/karlsbergplus/patches.prm']
        # charmm_bin = '/usr/local/dolly/bin/charmm36'
        # charmm_bin = 'charmm36b1_64'
        charmm_bin = '/public/local/scratch/bin/charmm36b1_64'
        charmm_struct = charmm.Charmm_manager(charmm_bin=charmm_bin, top=top, par=par)
        charmm_struct.set_titr_residues(titratable_residues_kbp)

        import tempfile
        tmp_folder_base = '/public/local/scratch/charmm_tmp/'
        tmp_charmm_conf_folder = tempfile.mkdtemp(prefix=tmp_folder_base) + '/'

        print "Starting CHARMM confE runs for " + self.jobname
        sleep = True
        while sleep:

            prepared_frame_folders = []
            kbp_jobfolders = os.listdir(self.kbp_folder + 'done/')
            for kbp_frame_jobfolder in kbp_jobfolders:
                frame_workdir_final = charmm_conf_folder + kbp_frame_jobfolder + '/'
                frame_workdir = tmp_charmm_conf_folder + kbp_frame_jobfolder + '/'

                if os.path.exists(frame_workdir_final):
                    charmm_output_filename = "%s%s_charmm.out" % (frame_workdir_final, kbp_frame_jobfolder)

                    if os.path.exists(charmm_output_filename):
                        f_out = open(charmm_output_filename, 'r')
                        charmm_normal_termination = False
                        for next_line in f_out:
                            if 'NORMAL TERMINATION BY NORMAL STOP' in next_line:
                                charmm_normal_termination = True
                        f_out.close()

                        if charmm_normal_termination:
                            continue
                        else:
                            print("resubmitting: " + frame_workdir_final)
                            shutil.rmtree(frame_workdir_final)
                            os.mkdir(frame_workdir)
                else:
                    os.mkdir(frame_workdir)

                pqr_filename = self.kbp_folder + 'done/' + kbp_frame_jobfolder \
                               + '/c_pH7_%s.reference.pqr' % kbp_frame_jobfolder
                pqr_in_ref_state = file_parser.Simple_struct_parser()
                pqr_in_ref_state.read_pdb(pqr_filename, is_pqr=True)

                if not charmm_struct.structures_checked:
                    charmm_struct.add_structure(pqr_in_ref_state)

                    charmm_struct.check_structures(quiet=True)

                    # Define the reference protonation vector
                    prot_vector = charmm_struct.get_titr_residue_dict(all=True)
                    ref_prot_vector = {}
                    for residue_tuple in prot_vector.keys():
                        ref_prot_vector[residue_tuple] = 0
                    charmm_struct.apply_titr_residue_dict(ref_prot_vector)

                    i = charmm_struct.charmm_config['tasks'].index('build_missing')
                    charmm_struct.charmm_config['tasks'].pop(i)
                    i = charmm_struct.charmm_config['tasks'].index('hbuild')
                    charmm_struct.charmm_config['tasks'].pop(i)
                    i = charmm_struct.charmm_config['tasks'].index('write_structure')
                    charmm_struct.charmm_config['tasks'].pop(i)
                    i = charmm_struct.charmm_config['tasks'].index('autogen')
                    charmm_struct.charmm_config['tasks'].pop(i)

                    charmm_struct.clean_up = True

                    # charmm_commands = 'ENERgy EPS 4.0\n'
                    charmm_commands = 'ENERgy EPS 4.0 E14Fac 0.0\n'
                    # charmm_commands = 'ENERgy EPS 8.0 E14Fac 0.0\n'
                    # charmm_commands = 'ENERgy EPS 2.0 E14Fac 0.0\n'
                    charmm_struct.add_charmm_command(charmm_commands, 'read_water_coordinates')
                    charmm_commands = 'BOMBlev -4\n'
                    charmm_struct.add_charmm_command(charmm_commands, 'read_structure_sequence')
                    charmm_commands = 'BOMBlev 0\n'
                    charmm_struct.add_charmm_command(charmm_commands, 'patches')

                    # Add a patch in case CT1 has been used to avoid problem with renamed atoms.
                    # The following is a copied light version of the corresponding section
                    # in Charmm_manager.__generate_charmm_input_script. Not a very elegant solution..
                    # The patch inserted here has to be the last one! This is happening, since all other patches are
                    # inserted after the section 'patches' in the script.
                    for segname in charmm_struct.charmm_instructions['ter'].keys():
                        nter, cter = charmm_struct.charmm_instructions['ter'][segname]
                        if cter is None:
                            pass
                        elif not ((('OXT' in cter) or ('OT1' in cter) and not ('HT1' in cter))):
                            if cter.resname == 'PRO':
                                pass
                            else:
                                charmm_commands = 'patch CT1N %s %i\n' % (cter.segname, cter.resid)
                                charmm_struct.add_charmm_command(charmm_commands, 'patches')

                else:
                    charmm_struct.transfer_coordinates(pqr_in_ref_state)
                charmm_struct.workdir = frame_workdir
                charmm_struct.title = kbp_frame_jobfolder

                # # Write charmm object into pickle
                # pickle_charmm_struct_filename = frame_workdir + 'charmm_struct.pkl'
                # f = open(pickle_charmm_struct_filename, 'wb')
                # pickle.dump(charmm_struct, f, protocol=2)
                # f.close()
                #
                # # Prepare and write python script
                # python_script = """#!/usr/bin/python
                # # -*- coding: utf-8 -*-
                # import cPickle as pickle
                # f = open('%scharmm_struct.pkl', 'rb') """ % frame_workdir + """
                # charmm_struct = pickle.load(f)
                # f.close()
                # charmm_struct.run_charmm()
                # """
                # f = open(frame_workdir + 'run.py', 'w')
                # for line in python_script.split('\n'):
                #     f.write(line.strip(' ') + '\n')
                # f.close()
                #
                # # Prepare and write shell script
                # shell_script = """#!/bin/tcsh
                # echo "Running on `hostname`"
                # cd %s""" % frame_workdir + """
                # python -u run.py
                # """
                # f = open(frame_workdir + 'run.sh', 'w')
                # for line in shell_script.split('\n'):
                #     f.write(line.strip(' ') + '\n')
                # f.close()

                prepared_frame_folders.append(frame_workdir)
                charmm_struct.run_charmm()
                if not charmm_struct.charmm_normal_termination:
                    print("Job crashed: " + frame_workdir)

                # os.remove(frame_workdir + 'charmm_struct.pkl')
                shutil.move(frame_workdir, frame_workdir_final)

            os.rmdir(tmp_charmm_conf_folder)
            # shell = subprocess.Popen('tcsh\n',\
            # stdin=subprocess.PIPE,\
            # stdout=subprocess.PIPE,\
            # stderr=subprocess.PIPE,\
            # shell=True\
            # )

            # script_folder = charmm_conf_folder
            # c = 0
            # n = 0
            # max_scripts_per_job = 999
            # shell_script = None
            # for i, prepared_frame_folder in enumerate(prepared_frame_folders):
            #     if shell_script is None:
            #         shell_script = "#!/bin/tcsh\n"
            #     shell_script += "chmod +x " + prepared_frame_folder + 'run.sh' + '\n'
            #     shell_script += prepared_frame_folder + 'run.sh' + '\n'
            #     c += 1
            #     if (c == max_scripts_per_job) or (i == len(prepared_frame_folders) - 1):
            #         script_name = script_folder + 'run%i.sh' % n
            #         f = open(script_name, 'w')
            #         for line in shell_script:
            #             f.write(line)
            #         f.write('\n')
            #         f.close()
            #         shell.stdin.write('qsub -o %s -e %s %s\n' % (script_folder, script_folder, script_name))
            #         shell_script = None
            #         n += 1
            #         c = 0

            # for prepared_frame_folder in prepared_frame_folders:
            #     shell.stdin.write('qsub -o %s -e %s %s\n' % (prepared_frame_folder, prepared_frame_folder,
            #                                                  prepared_frame_folder + 'run.sh'))

            # shell.stdin.write('exit\n')

            sleep = False

        print "Finished CHARMM confE runs done for " + self.jobname

        self.status = 'done'
        print("Job Done!")
        self.dump()

        return self.status

        # Todo: catch negative ion number in Neutralize.inp


    def stop(self):
        """
        Stop the job.
        """
        pass

    def dump(self):
        self.new = False
        job_backup_name = self.root_workdir + 'job_dump.pickle'
        f = open(job_backup_name, 'w')

        sq = self.status_queue
        cq = self.cmd_queue
        self.status_queue = None
        self.cmd_queue    = None
        pickle.dump(self.__dict__, f)
        self.status_queue = sq
        self.cmd_queue    = cq

        f.close()

        stop = self.check_queue()
        if stop:
            self.status = 'stopping'

        return 1

    def recover(self):
        job_backup_name = self.root_workdir + 'job_dump.pickle'
        if os.path.exists(job_backup_name):
            f = open(job_backup_name, 'r')
            recovered_job = pickle.load(f)
            f.close()

            sq = self.status_queue
            cq = self.cmd_queue
            self.__dict__ = recovered_job
            self.status_queue = sq
            self.cmd_queue    = cq

            self.check_queue()

            return 1
        else:
            return 0

    def check_queue(self):
        if self.cmd_queue is not None and self.status_queue is not None:
            stop = False
            # Look for commands.
            while not self.cmd_queue.empty():
                cmd = self.cmd_queue.get(block=False)
                if cmd == 'stop':
                    stop = True
            # Erase old status entry if available.
            while not self.status_queue.empty():
                try:
                    self.status_queue.get(block=False)
                except:
                    pass

            sq = self.status_queue
            cq = self.cmd_queue
            self.status_queue = None
            self.cmd_queue    = None
            sq.put(self)
            self.status_queue = sq
            self.cmd_queue    = cq



            return stop
        else:
            return False



def run_MD_pKa_job(job, status_queue, cmd_queue):
    job.status_queue = status_queue
    job.cmd_queue    = cmd_queue

    job.check_queue()

    job.start

    #    stop = False
#    while not stop:
#
#        while not cmd_queue.empty():
#            cmd = cmd_queue.get(block=False)
#            if cmd == 'stop':
#                stop = True
#        while not status_queue.empty():
#            try:
#                status_queue.get(block=False)
#            except:
#                pass
#
#        status_queue.put('I am still running')
#
#        if not stop:
#            time.sleep(10)

    # Clear command queue.
    while not cmd_queue.empty():
        cmd_queue.get(block=False)
#    while not status_queue.empty():
#        try:
#            status_queue.get(block=False)
#        except:
#            pass

    print(job.root_workdir + '   DONE')

# Wieder einsammeln (-> get falls fertig?)

class pka_manager(object):
    class pka_Manager_Error(Exception):
            def __init__(self, value):
                self.value = value
            def __str__(self):
                return repr(self.value)

#    def dump_job(self, job):
#        job.new = False
#        job_backup_name = job.root_workdir + 'job_dump.pickle'
#        f = open(job_backup_name, 'w')
#        pickle.dump(job, f)
#        f.close()
#
#    def recover_job(self, job):
#        job_backup_name = job.root_workdir + 'job_dump.pickle'
#        if os.path.exists(job_backup_name):
#            f = open(job_backup_name, 'r')
#            recovered_job = pickle.load(f)
#            f.close()
#            job.__dict__ = recovered_job
#            return 1
#        else:
#            return -1

    def __init__(self):
        """
        Initialize a new pka_manager object.
        """
#        self.status = {}
#        self.status['running'] = False

        self.queued_jobs  = []
        self.running_jobs = []
        self.done_jobs    = []

        self.processes = {}

    def start(self):
        """
        Starts the main loop to process the jobs in a separate thread.
        """

        for j in self.jobs:
            pass

#            j.start()
#            print("Job %s finished with status: %s" % (j.jobname, j.status))


    def status(self):
        """
        Reports the current status of the submitted jobs.
        """
        pass
    def stop(self):
        """
        Stops the main loop, as soon as possible.
        """
        pass
    def cancel(self):
        """
        Stops all jobs immediately. Eventually some external jobs must be stopped manually.
        """
        pass


    def start_server(self):
        address = ('localhost', 6000)
        server = Listener(address, authkey=b'MD_pKa_manager')

        print("MD_pKa_manager server started.")

        stop = False
        while not stop:
            # Wait for a new connection.
            conn = server.accept()

            try:
#                print('# Running Jobs:')
#                finished_processes = []
#                for root_folder in self.processes:
#                    print('  ' + root_folder)
#                    (p, s_q, c_q) = self.processes[root_folder]
#                    print('  ' + str(p.is_alive()))
#                    print('')
#                    if not p.is_alive():
#                        p.join()
#                        finished_processes.append(root_folder)
#                for root_folder in finished_processes:
#                    self.processes.pop(root_folder)


                # Receive data.
                msg = conn.recv()
                answer = 'received'

                # Process the received data.
                if type(msg) == str:
                    print("message received: " + msg)
                    if msg == 'stop':
                        stop = True
                elif type(msg) == tuple:
                    if len(msg) == 2 and type(msg[0]).__name__ == 'MD_pKa_job' and type(msg[1]) == str:
                        job = msg[0]
                        cmd = msg[1]
                        # The root work dir is used as unique ID
                        root_folder = job.root_workdir
                        start_job = False
                        if cmd == 'add':
                            if self.processes.has_key(root_folder):
                                print('Already running job received.')
                                answer = 'running'
                            else:
                                print('New job received.')
                                answer = 'added'
                                start_job = True

                        elif cmd == 'restart':
                            print('Restart request received.')
                            if self.processes.has_key(root_folder):
                                (p, s_q, c_q) = self.processes[root_folder]
                                if p.is_alive():
                                    print('Restart request rejected.')
                                    answer = 'Job is still running, please stop it first.'
                                else:
                                    print('Restart request accepted.')
                                    answer = 'Job restarted.'
                                    while not s_q.empty():
                                        # When the process is dead, there should be no risk in reading the queue.
                                        s_q.get()
                                    p.join()
                                    start_job = True
                            else:
                                print('New job received.')
                                answer = 'added'
                                start_job = True

                        elif cmd == 'get':
                            print('Get request received.')
                            if self.processes.has_key(root_folder):
                                (p, s_q, c_q) = self.processes[root_folder]
                                answer = s_q.get()
                                # Put pack the answer for the next time.
                                s_q.put(answer)
                            else:
                                answer = 'unknown job'

                        elif cmd == 'stop':
                            print('stop request received.')
                            (p, s_q, c_q) = self.processes[root_folder]
#                            answer = s_q.get()
                            c_q.put('stop')
                            answer = 'stop signal sent'

                        # Start a job in a new process.
                        if start_job:
                            manager = Manager()
                            s_q = manager.Queue()
                            c_q = manager.Queue()
                            p = Process(target=run_MD_pKa_job, args=(job,s_q, c_q))
                            p.start()
                            self.processes[root_folder] = (p, s_q, c_q)

                    elif len(msg) == 2 and type(msg[0]) == str and type(msg[1]) == str:
                        root_folder = msg[0]
                        cmd         = msg[1]
                        if cmd == 'kill':
                            print('Kill request received.')
                            if self.processes.has_key(root_folder):
                                (p, s_q, c_q) = self.processes[root_folder]
                                c_q.put('stop')
                                answer = 'job deleted'
                                try:
#                                    while not s_q.empty():
#                                        # When the process is dead, there should be no risk in reading the queue.
#                                        s_q.get()
                                    p.get(timeout=60)
                                    p.join(timeout=60)
                                except:
                                    answer = 'no answer'
                                self.processes.pop(root_folder)
                            else:
                                answer = 'Unknown Job.'
                            print('Kill request done with status: ' + answer)

                # Send an answer and close the connection.
                conn.send(answer)
                conn.close()

            except:
                print("An Error occured:")
                print(sys.exc_info())

        server.close()

class pka_session(object):
    def __init__(self):
        """
        Initialize a new pka_session object.
        """
        self.jobs = []

    def add_job(self, modelling_prefix ,root_workdir, jobname='', top=[], par=[], md_template_folder=None, ph=7,
                structure_filename=''):
        """
        Add a new job to the manager. This objects takes a PDB file as input, models the structure, starts
        MD simulations, start Karlsberg+ and finally it collects all data to calculate the pKa values.

        Parameters:
        modelling_prefix: prefix for the structure files in the modelling sub folder.
        root_workdir:     Path to the folder where all results will be stored. The folder must not exist.
        jobname:          A name used as prefix for output files.
        """
        if not jobname:
            jobname = modelling_prefix

        if not os.path.exists(root_workdir):
            os.mkdir(root_workdir)

        job = MD_pKa_job(modelling_prefix, root_workdir, jobname, top=top, par=par,
                         md_template_folder=md_template_folder, ph=ph, structure_filename=structure_filename)

        self.jobs.append(job)

        return 1

    def send(self, obj, cmd):
        address = ('localhost', 6000)
        conn = Client(address, authkey=b'MD_pKa_manager')
        conn.send( (obj, cmd) )
        msg = conn.recv()
        conn.close()
        return msg

    def send_jobs_to_server(self):
#        address = ('localhost', 6000)
        for job in self.jobs:
#            conn = Client(address, authkey=b'MD_pKa_manager')
#            conn.send( (job, 'add') )
#            msg = conn.recv()
            msg = self.send(job, 'add')
            # print(msg)
#            conn.close()

    def get_status(self):
#        address = ('localhost', 6000)
        jobs = []
        for job in self.jobs:
#            conn = Client(address, authkey=b'MD_pKa_manager')
#            conn.send( (job, 'get') )
#            msg = conn.recv()
            msg = self.send(job, 'get')
            jobs.append( msg )
#            conn.close()
        return jobs

    def stop(self):
#        address = ('localhost', 6000)
        for job in self.jobs:
#            conn = Client(address, authkey=b'MD_pKa_manager')
#            conn.send( (job, 'stop') )
#            msg = conn.recv()
            msg = self.send(job, 'stop')
            print(msg)
#            conn.close()

    def restart(self):
        for job in self.jobs:
            msg = self.send(job, 'restart')
            print(msg)

    def kill_job(self, root_folder):
        if root_folder[-1] != '/':
            root_folder += '/'
        for job in self.jobs:
            print job.root_workdir
            if job.root_workdir == root_folder:
                msg = self.send(job, 'kill')
                print(msg)

#    def stop_job(self, root_folder):
#        msg = self.send(root_folder, 'stop')
#        print(msg)


if __name__ == '__main__':
    pka_m = pka_manager()

    pka_m.start_server()

##    folder = '/user/tmeyer/workspace/projects/md_pkas/snase/prototypes/'
#    folder = '/scratch/scratch/tmeyer/md_pka/md_pka_manager/runs/'
#    pdb_folder    = '/user/tmeyer/workspace/projects/md_pkas/snase/pdbs/all/'
#
##    pdbs = ['1STN', '1TQO', '1TR5', '2OXP', '3P75', '3ITP']
##    pdbs = ['1TQO', '1STN']
#    pdbs = ['1TQO']
#
#
#    pka_m = pka_manager()
#
#    for pdb in pdbs:
#        if not pdb:
#            print("'" + pdb + "'")
#            continue
#        f = folder + pdb
##        if os.path.exists(f):
##            shutil.rmtree(f)
#
#        p = pdb_folder + pdb + '.pdb'
#        msg = pka_m.add_job(p, f, jobname=pdb)
#        if msg != 1:
#            print(msg)
#
#    pka_m.start()






















