# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:19:36 2011

@author: tmeyer
"""

import subprocess
import os
import time
#import pp
import shutil

class kbp_job(object):
    def __init__(self, pdb=None, template=None, single_folder_run=None):
        self.parameter = {}
        self.kb_para = {}

        # If set to a string, the string will be written into a file named cavities.conf in each job directory.
        # This will trigger cavity_tapbs_1.0 to include cavities in the calculations.
        # Syntax: cavity <cavity parameter> <grid resolution> <unused> 0.9 0.25 0.0
        # - cavity parameter: 0.7 - 0.9 does make sense
        # - grid resolution:  smaller that 0.25 will result in very long calculations.
        # e.g. 0.9 0.25 0.0
        self.parameter['cavity_parameter'] = -1

        # If set to a path, a water tmeplate file for the cavity calculations is copied from here.
        # Assumes that a file named <self.parameter['jobid']>.pdb exitst in this folder, containing water molecules.
        self.parameter['water_folder'] = ''
        
        ### set default general parameters ###
        #self.parameter['exec'] = 'titrate.pl'
#        self.parameter['exec'] = "perl -I/scratch/scratch/tmeyer/kbplus2 /scratch/scratch/tmeyer/kbplus2/titrate.pl"
        self.parameter['exec'] = "perl -I/scratch/scratch/tmeyer/kbplus2 /scratch/scratch/tmeyer/kbplus2/titrate.pl"
        self.parameter['jobid'] = None
        
        # For additional charmm commands
        self.parameter['charmm_commands'] = ''
        
        ### set default parameters for Karlsberg+2 ###
        self.kb_para = {}
        self.kb_para['key_order'] = ['globals', 'sfc_globals', 'sfc',\
                     'ce_globals', 'kb_globals']
        
        self.kb_para['globals'] = {\
            'homeDir':    None,\
            'include':   '/scratch/scratch/tmeyer/kbplus2/',\
            'workDir':   '/public/scratch/tmeyer/kb',\
            'cpus':      1,\
            'key_order': ['homeDir','include','workDir','cpus'],\
            'unnamed_entries' : []
            }
            
#            'charmm': 'charmm34b1_64',\            
        self.kb_para['sfc_globals'] = {\
            'pdb': None,\
            'rtf': '/scratch/scratch/tmeyer/karlsbergplus/top.inp',\
            'prm': '/scratch/scratch/tmeyer/karlsbergplus/par.inp',\
            'rtfPatches': '/scratch/scratch/tmeyer/karlsbergplus/patches.rtf',\
            'prmPatches': '/scratch/scratch/tmeyer/karlsbergplus/patches.prm',\
            'titratable': 'titratable.yaml',\
            'charmm': '/scratch/scratch/tmeyer/charmm35',\
            #'tapbs': 'tapbs_cavities',\
            'tapbs': 'cavity_tapbs_1.0',\
            'karlsberg': 'karlsberg2.x86_64',\
            'sasa': '/scratch/scratch/gernotf/sasa/sasa_static.bin',\
            'preOpt': 0,\
            'key_order': ['pdb', 'rtf','prm','rtfPatches','prmPatches',\
                          'titratable','charmm','tapbs','karlsberg','sasa',\
                          'preOpt'],\
            'unnamed_entries' : []\
            }
        
        # 'name' entry will be used as prefix für final 'name'
        self.kb_para['sfc'] = []
        self.kb_para['sfc'].append( {\
            '- name': 'c_pH7',\
            '  pH': 7,\
            '  redox': 0,\
            '  randomSeed': 1,\
            '  optimize': '',\
            'key_order': ['- name', '  pH', '  redox', '  randomSeed', '  optimize' ],\
            'unnamed_entries' : [\
                "   - 'what=HYDROGENS,confs=0,kfc=0'"]\
            } )

        # cexe: apbs
        # 'name' entry will be used as suffix für final 'name'
        self.kb_para['ce_globals'] = {\
            'aexe': 'apbs',\
            'name': 'ce',\
            'key_order': ['aexe', 'name'],\
            'unnamed_entries' : []\
            }
        
        self.kb_para['kb_globals'] = {\
            'name': 'kb',\
            'exe': 'karlsberg2.x86_64',\
            'pH_start': -10,\
            'pH_end': 20,\
            'pH_incr': 0.5,\
            'redox_start': 0,\
            'redox_end': 0,\
            'redox_incr': 25,\
            'temperatures': '',\
            'key_order': ['name', 'exe', 'pH_start', 'pH_end', 'pH_incr',\
                          'redox_start','redox_end','redox_incr','temperatures'],\
            'unnamed_entries' : [\
                ' - 300',\
                ' - 373',\
                ' - 465',\
                ' - 579',\
                ' - 722',\
                ' - 900'\
                ]\
            }
        
        
        ### overwrite parameters, if template is provided ###
        if template != None:
            self.parameter = dict( template.parameter )
            self.kb_para = dict( template.kb_para )
            self.single_folder_run = template.single_folder_run
        else:
            self.parameter['pdb_dir'] = None
            self.parameter['run_dir'] = None
            self.parameter['done_dir'] = None
            if single_folder_run is None:
                self.single_folder_run = False
            else:
                self.single_folder_run = single_folder_run

        if pdb != None:
            self.parameter['pdb'] = pdb
        else:
            self.parameter['pdb'] = None


        
        
        ### initialize misc. variables ###
        
        self.running = 0
        self.done = 0
        
        self.stdout = []
        self.stderr = []

        
    class jobManagerError(Exception):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)  
        
    def create_kb_script(self, job_dir, pdb=None, keep_job_names=False):
        self.kb_para['globals']['homeDir']  = job_dir
        if pdb is not None:
            self.kb_para['sfc_globals']['pdb']  = pdb
        else:
            self.kb_para['sfc_globals']['pdb']  = self.parameter['pdb']
        if not keep_job_names:
            for sfc in self.kb_para['sfc']:
                sfc['- name']                  += '_' + self.parameter['jobid']
            self.kb_para['ce_globals']['name'] += '_' + self.parameter['jobid']
            self.kb_para['kb_globals']['name'] += '_' + self.parameter['jobid']
        
        script = '---\n'
        for l0 in self.kb_para['key_order']:
            
            if l0 == 'sfc':
                
                script += l0 + ':\n'
                
                for l0_index in range(len(self.kb_para[l0])):
                    sfc = self.kb_para[l0][l0_index]
                    
                    for l1 in sfc['key_order']:
                        script += ' ' + str( l1 ) + ': ' \
                                + str( sfc[ l1 ] ) + '\n'
                
                    for l1 in sfc['unnamed_entries']:
                        script += ' ' + str( l1 ) + '\n'
                        
                    script += '\n'
                    
                script += '\n'
            else:
                script += l0 + ':\n'
                
                for l1 in self.kb_para[l0]['key_order']:
                    script += ' ' + str( l1 ) + ': ' \
                            + str( self.kb_para[l0][ l1 ] ) + '\n'
            
                for l1 in self.kb_para[l0]['unnamed_entries']:
                    script += ' ' + str( l1 ) + '\n'
                    
                script += '\n'                
                
        return script
        
    def start(self, keep_job_names=False, single_folder_run=None):
        # The single_folder_run parameter should not be used here. It can be set via the class constructor.
        # Kept for backward compatibility.

        if single_folder_run is None:
            single_folder_run = self.single_folder_run

        # Make sure slash convention is fulfilled.
        if single_folder_run:
            folder_list = ['run_dir']
            # The following two folder are not used with single_folder_run.
            self.parameter['pdb_dir']  = ''
            self.parameter['done_dir'] = ''
        else:
            folder_list = ['run_dir', 'pdb_dir', 'done_dir']
        for folder in folder_list:
            if self.parameter[folder][-1] != '/':
                self.parameter[folder] += '/'

        ### last chance to set default values ###
        if self.parameter['jobid'] == None:
            for suffix in ['.pdb', '.ent']:
                suffix_start = self.parameter['pdb'].find(suffix)
                if suffix_start != -1:
                    jobid = self.parameter['pdb'][:suffix_start]
                    break
            else:
                jobid = self.parameter['pdb']

            # Strip folder if necessary.
            folder_end = jobid.rfind('/')
            if folder_end != -1:
                jobid = jobid[folder_end+1:]

            # jobid = self.parameter['pdb'].strip('.pdb.ent')
            self.parameter['jobid'] = jobid
        
        for key in self.parameter.keys():
            if self.parameter[key] is None:
                print "ERROR in kbp_job: Parameter " + str(key) + " is not defined."
                return -1

        self.running = 1
        
        #print 'starting job: ' + str(self.parameter['jobid'])

        ### start Karlsberg+ ###
        self.run_karlsbergplus(keep_job_names, single_folder_run)
        
        self.done = 1
        
        return self
        
    def run_karlsbergplus(self, keep_job_names=False, single_folder_run=False):
        if single_folder_run:
            job_dir = self.parameter['run_dir']
            job_done_dir = self.parameter['done_dir']
            if job_done_dir == '':
                job_done_dir = job_dir
        else:
            job_dir = self.parameter['run_dir'] + self.parameter['jobid']
            job_done_dir = self.parameter['done_dir'] + self.parameter['jobid']

        if job_dir[-1] != '/':
            job_dir += '/'
        if job_done_dir != '' and job_done_dir[-1] != '/':
            job_done_dir += '/'

        # Get PDB filename without path.
        path_end = self.parameter['pdb'].rfind('/')
        if path_end != -1:
            pdb_wo_path = self.parameter['pdb'][path_end+1:]
        else:
            pdb_wo_path = self.parameter['pdb']

            
        ### create input script for Karlsberg+2 ###
        kb_script = self.create_kb_script(job_dir, pdb_wo_path, keep_job_names)
        
        
        # store current directory
        old_dir = os.getcwd()


        if not single_folder_run:
            ### create job dir ###
            if os.path.exists(job_dir):
#                s = "ERROR in job_manager:run_karlsbergplus: job directory exists"\
#                        + "already in run_dir."
#                raise self.jobManagerError(s)
                print "skipping " + self.parameter['jobid']\
                        + " since job directory exists already in run_dir"
                return -1
            
            if os.path.exists(job_done_dir):
#                s = "ERROR in job_manager:run_karlsbergplus: job directory exists"\
#                        + "already in done_dir."
#                raise self.jobManagerError(s)
                print "skipping " + self.parameter['jobid']\
                        + " since job directory exists already in done_dir"

                status = 1
                if os.path.exists(job_done_dir + '/pka_results.dat'):
                    f = open(job_done_dir + '/pka_results.dat', 'r')
                    self.stdout = f.readlines()
                    f.close()
                else:
                    self.stdout = ["The file 'pka_results.dat' does not exist."]
                    status = -1
                if os.path.exists(job_done_dir + '/pka_error.dat'):
                    f = open(job_done_dir + '/pka_error.dat', 'r')
                    self.stderr = f.readlines()
                    f.close()
                    status = 1
                else:
                    self.stdout = ["The file 'pka_error.dat' does not exist."]
                    status = -1


                return status

            os.mkdir( job_dir )

            timeouts = 0
            while True:
                time.sleep(5)
                if os.path.exists(job_dir):
                    break
                timeouts += 1
                print('Folder %s does not exist. Will wait 5s and try again.' % job_dir)
                if timeouts > 10:
                    raise AssertionError('Could not create folder %s' % job_dir)

            
            
        ### enter job dir ###
        os.chdir( job_dir )
        
        
        ### start external shell ###
        f = subprocess.Popen('csh\n',\
                stdin=subprocess.PIPE,\
                stdout=subprocess.PIPE,\
                stderr=subprocess.PIPE,\
                shell=True\
                )
        #        self.process = subprocess.Popen(self.parameter['exec'],\
#        f.stdin.write( 'pwd ' + '\n' )
        f.stdin.write( 'cd ' + job_dir + '\n' )
        
        ### copy PDB file to run_dir ###
        if single_folder_run:
            if self.parameter['pdb'].find(job_dir) == -1:
                command = 'cp ' + self.parameter['pdb'] + ' .'
            else:
                command = ''
        else:
            command = 'cp ' + self.parameter['pdb_dir']\
                                     + self.parameter['pdb'] + ' '\
                                     + pdb_wo_path
        if command:
            f.stdin.write( command + '\n' )
        
        
        filename = job_dir + '/titrate.yml'
        f_kb = open(filename, 'w')
        f_kb.write(kb_script)
        f_kb.close()



        if self.parameter['cavity_parameter'] > -1:
            #print("cavity calculation is activated.")
            # Copy pdb file with template waters if available.
            if self.parameter['water_folder'] != '':
                #print("Water folder is set up to: " + self.parameter['water_folder'])
                water_folder = self.parameter['water_folder']
                if water_folder[-1] != '/':
                    water_folder += '/'
                water_filename = water_folder + self.parameter['jobid'] + '.pdb'
                assert os.path.exists(water_filename)
                shutil.copy(water_filename, job_dir + 'water.pdb')

            # Set up cavity calculation.
            filename = job_dir + '/cavities.conf'
            f_cav = open(filename, 'w')
            f_cav.write(self.parameter['cavity_parameter'] + '\n')
            f_cav.close()

        #filename = job_dir + '/cavities.conf'
        #f_cav = open(filename, 'w')
        #f_cav.write('cavity 0.0 0.2 0.7\n')
        #f_cav.close()



        if self.parameter['charmm_commands'] != '':
            shutil.copy(self.parameter['charmm_commands'], job_dir + '/charmm_commands.inp')
        
        
        ### rename HSP and HSE in HIS
        f.stdin.write('sed -i s/HSE/HIS/g %s \n' % pdb_wo_path)
        f.stdin.write('sed -i s/HSP/HIS/g %s \n' % pdb_wo_path)
        f.stdin.write('sed -i s/HSD/HIS/g %s \n' % pdb_wo_path)
        
        
        ### run karlsberg ###
        f.stdin.write(self.parameter['exec'] + ' titrate.yml \n')


        if not single_folder_run:
            ### move results to done_dir ###
            command = 'mv ' + job_dir + ' ' + self.parameter['done_dir']
            f.stdin.write( command + '\n' )

        
        ### stop process ###
        f.stdin.write('exit \n')

        
        ### enter previous directory ###
        os.chdir( old_dir )
        
        
        ### get stdout ###
        while True:
            next_line = f.stdout.readline()
            
            # if karlsberg is terminates
            if not next_line:
                break
            
            next_line = next_line.strip('\n ')
            
            self.stdout.append( next_line)
            
        ### get stderr ###
        while True:
            next_line = f.stderr.readline()
            
            if not next_line:
                break
            
            next_line = next_line.strip('\n ')
            
            self.stderr.append( next_line)

        
        #f = open(self.parameter['done_dir'] + '/' + self.parameter['jobid']\
        #        + '/pka_results.dat', 'w')
        if single_folder_run:
            final_folder = job_dir
        else:
            final_folder = job_done_dir

        f = open(final_folder + '/pka_results.dat', 'w')
        for line in self.stdout:
            f.write(line + '\n')
        f.close()
#        f = open(self.parameter['done_dir'] + '/' + self.parameter['jobid']\
#                + '/pka_error.dat', 'w')
        f = open(final_folder + '/pka_error.dat', 'w')
        for line in self.stderr:
            f.write(line + '\n')
        f.close()
       
        
        return 1

 

    def read_occupancies(self):
        import storable
        import sys
        import pprint
    
        p = pprint.PrettyPrinter(indent=4)
        for f in sys.argv[1:]:
            p.pprint(storable.retrieve(f))         
   








     

#kbp_template = kbp_job()
#
##kbp_template.parameter['pdb'] = "2LZT.pdb"
#kbp_template.parameter['pdb_dir'] = "/user/tmeyer/workspace/script/python/karlsberg/pdbs"
#kbp_template.parameter['run_dir'] = "/scratch/scratch/tmeyer/kbp2jobs/tmp_testrun"
#kbp_template.parameter['done_dir'] = "/user/tmeyer/workspace/script/python/karlsberg/done"
#
#
##k = kbp_template
##k.start()
#
#job1 = kbp_job('2LZT.pdb' ,kbp_template)
#job2 = kbp_job('2LZT_1.pdb' ,kbp_template)
#
##job1.start()
##job2.start()
#
##print "# stdout:"
##for m in k.stdout:
##    print m
##print "# stderr:"
##for m in k.stderr:
##    print m
##print "---"
#
#
#### prepare job server ###
#
## number of processes
#njobs = 5
#
## tuple of all parallel python servers to connect with
##ppservers = ()
##ppservers = ("10.0.0.1",)
#
##job_server = pp.Server(njobs, ppservers=ppservers)
#job_server = pp.Server(njobs)
#
#print "Starting pp with", job_server.get_ncpus(), "workers"
#
#
##submitted_jobs = []
##submitted_jobs.append( job_server.submit(job1.start, (), modules=('os','subprocess')) )
##submitted_jobs.append( job_server.submit(job2.start, (), modules=('os','subprocess')) )
#
##submitted_jobs.append( job_server.submit(dosth, (), modules=("time",)) )
##submitted_jobs.append( job_server.submit(dosth, (), modules=("time",)) )
#
#
##### prepare jobs ###
#pdb_dir = '/user/tmeyer/workspace/script/python/karlsberg/pdbs'
#fd = os.listdir( pdb_dir )
#pdb_list = []
#for entry in fd:
#    if os.path.isfile( pdb_dir + '/' + entry ):
#        pdb_list.append( entry )
#
#kbpj_prepared = []
#for pdb in pdb_list:
#    job = kbp_job( pdb, kbp_template )
#    kbpj_prepared.append(job)
#
#    
#### submitt jobs ###
#kbj_submitted = []
#for job in kbpj_prepared:
#    kbj_submitted.append( job_server.submit(job.start, (), modules=('os','subprocess',)) )
#
#  
#### collect jobs (and wait for them) ###
#kbpj_finished = []
#for job in kbj_submitted:
#    kbpj_finished.append( job() )
#
#    
#for kbpj in kbpj_finished:
#    print '##################################'
#    for l in kbpj.stdout:
#        print l
#    print
#
#
#### print stats ###
#job_server.print_stats()



