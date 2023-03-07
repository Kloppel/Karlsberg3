# coding=utf-8
import subprocess


import os, re, time, sys
import cPickle as pickle
import subprocess
from multiprocessing import Pool
import datetime

# Needed if kbp2 package is not globally installed.
sys.path.append('/user/tmeyer/workspace/script/protein_toolbox/packages/')

# sys.path.append('/user/tmeyer/workspace/script/protein_toolbox')
# import job_manager2
from kbp2 import job_manager2

def job(path, job_file, job_data):
    script = job_data['script'] + '\n exit \n'
    print("starting job: " + script)
    process = subprocess.Popen("tcsh",\
                    shell=True,\
                    stdin=subprocess.PIPE,\
                    stdout=subprocess.PIPE,\
                    stderr=subprocess.PIPE,\
                    )
    process.stdin.write(script)

    output = ''
    while True:
        next_line = process.stdout.readline()
        if not next_line:
            break
        output += next_line
#            print next_line
#                    if next_line.find('NORMAL TERMINATION BY NORMAL STOP'):
#                        self.charmm_normal_termination = True

    if os.path.exists(path + job_file + ' - running'):
        os.remove(path + job_file + ' - running')
    f = open(path + job_file + ' - done', 'w')
    f.close()

def kbp_job(path, job_file, job_data):
    success = 1
    job_data.start()

    # try:
    #     finished_job = job_data.start
    # except:
    #     success = 0


    # if success == 0:
    #     finished_job = 'crashed'

    # f = open(path + job_file + ' - result', 'w')
    # pickle.dump(finished_job, f)
    # f.close()

    if os.path.exists(path + job_file):
        os.remove(path + job_file)

    c = 0
    while os.path.exists(path + job_file) and c < 10:
        print("File '%s' not deleted yet, waiting.." % path + job_file)
        c += 1
        time.sleep(30)
    if c == 10:
        print("Waring: File %s was not deleted!" % path + job_file)

    # c = 0
    # while os.path.exists(path + job_file) and c<7:
    #     time.sleep(10)
    #     c += 1
    if os.path.exists(path + job_file + ' - running'):
        os.remove(path + job_file + ' - running')



    # f = open(path + job_file + ' - done', 'w')
    # f.close()

    return success

#    f = open(path + job_file + ' - result', 'w')
#    pickle.dump(result, f)
#    f.close()

class listener(object):
    """
    Checks a folder for jobs and runs them.

    Jobs:
    A file in the format: "jobX.pickle" Where X is a number,

    Control statements
    "c - stop" : stopps the listener after finishing the last job.


    """
    def __init__(self, path):
        if path[-1] != '/':
            path += '/'
        self.path = path

    def run(self):
        run = True
        running_job_name = None

        if os.path.exists(self.path + 'c - stopped'):
            os.remove(self.path + 'c - stopped')
        f = open(self.path + 'c - running', 'w')
        f.close()

        while run:
            f = open(self.path + 'c - running', 'w')
            f.close()

            # Check existing files.
            jobs = []
            control = []
            files = os.listdir(self.path)
            for fn in files:
                # Look for jobs
                reg = re.compile(r'^job_([\w-]+).pickle$')
                reg_m = reg.match(fn)
                if reg_m is not None:
                    jobs.append( reg_m.groups()[0] )

                # Look for control statements
                reg = re.compile(r'^c - ([\w-]+)$')
                reg_m = reg.match(fn)
                if reg_m is not None:
                    control.append( reg_m.groups()[0] )

            for cmd in control:
                if cmd == 'stop':
                    run = False
                    os.remove(self.path + 'c - stop')
                    if os.path.exists(self.path + 'c - running'):
                        os.remove(self.path + 'c - running')


            if jobs and run:
                # Find the next unfinished job.
                for job in jobs:
                    running_job_name = 'job_%s.pickle' % job
                    if not os.path.exists(self.path + running_job_name + ' - done'):
                        break
                else:
                    running_job_name = None

                if running_job_name is not None:
                    print 'starting: ' + running_job_name
                    f = open(self.path + running_job_name + ' - running', 'w')
                    f.close()

                    # To make sure that the submitting process can finish writing.
                    time.sleep(5)

                    f = open(self.path + running_job_name)
                    running_job = pickle.load(f)
                    f.close()

                    script = running_job['script'] + '\n exit \n'
                    print script
                    process = subprocess.Popen("tcsh",\
                                    shell=True,\
                                    stdin=subprocess.PIPE,\
                                    stdout=subprocess.PIPE,\
                                    stderr=subprocess.PIPE,\
                                    )
                    process.stdin.write(script)

                    output = ''
                    while True:
                        next_line = process.stdout.readline()
                        if not next_line:
                            break
                        output += next_line
                        print next_line
    #                    if next_line.find('NORMAL TERMINATION BY NORMAL STOP'):
    #                        self.charmm_normal_termination = True


                    if os.path.exists(self.path + running_job_name):
                        os.remove(self.path + running_job_name)
                    # c = 0
                    # while os.path.exists(self.path + running_job_name) and c<7:
                    #     time.sleep(10)
                    #     c += 1
                    if os.path.exists(self.path + running_job_name + ' - running'):
                        os.remove(self.path + running_job_name + ' - running')



                    # f = open(self.path + running_job_name + ' - done', 'w')
                    # f.close()

                    running_job_name = None


            if not run:
                f = open(self.path + 'c - stopped', 'w')
                f.close()
                return 1
            else:
                print "waiting for somthing to do.."
                time.sleep(10)

    def run_parallel(self, max_processes):
        run = True
        # running_processes = 0
        running_jobs = []
        pool = Pool(processes=max_processes)

        if os.path.exists(self.path + 'c - stopped'):
            os.remove(self.path + 'c - stopped')
        if os.path.exists(self.path + 'c - stopping'):
            os.remove(self.path + 'c - stopping')
        f = open(self.path + 'c - running', 'w')
        f.close()

        while run:
            f = open(self.path + 'c - running', 'w')
            f.close()

            # Check existing files.
            jobs = []
            control = []
            files = os.listdir(self.path)
            files.sort()
            for fn in files:
                # Look for new jobs. (not running and not finished)
                reg = re.compile(r'^([^_]*)_?job_([\w-]+).pickle$')
                reg_m = reg.match(fn)
                if reg_m is not None:
                    # if os.path.exists(self.path + fn + ' - done'):
                    if (fn + ' - done') in files:
                        continue
                    # if os.path.exists(self.path + fn + ' - running'):
                    if (fn + ' - running') in files:
                        continue


                    if 'roifree' in fn:
                        continue


                    job_type = reg_m.groups()[0]
                    job_name = reg_m.groups()[1]
                    job_file = fn
                    jobs.append( (job_type, job_name, job_file) )

                # Look for control statements
                reg = re.compile(r'^c - ([\w-]+)$')
                reg_m = reg.match(fn)
                if reg_m is not None:
                    control.append( reg_m.groups()[0] )

            for cmd in control:
                if cmd == 'stop':
                    run = False
                    os.remove(self.path + 'c - stop')
                    if os.path.exists(self.path + 'c - running'):
                        os.remove(self.path + 'c - running')
                    f = open(self.path + 'c - stopping', 'w')
                    f.close()


            # Check the number of running processes.
            running_processes = 0
            running_jobs_new = []
            for rj in running_jobs:
                if not rj.ready():
                    running_jobs_new.append(rj)
                    running_processes += 1
                else:
                    rj.get()
            running_jobs = running_jobs_new


            if jobs and run:
                # To make sure that the submitting process can finish writing.
                time.sleep(5)

                while running_processes < max_processes and jobs:
                    running_processes += 1

                    (job_type, job_name, job_file) = jobs.pop(0)

                    print 'starting: ' + job_name

                    # Read the job data.
                    job_file_full = self.path + job_file
                    if not os.path.exists(job_file_full):
                        print("WARNING: Job file '%s' does not exist and is skipped!" % job_file)
                        continue

                    f = open(job_file_full)
                    job_data = pickle.load(f)
                    f.close()

                    f = open(self.path + job_file + ' - running', 'w')
                    f.close()

                    # Run the job depending on the job type.
                    if job_type == '':
                        running_jobs.append( \
                            pool.apply_async( \
                                job, (self.path, job_file, job_data) \
                            ) \
                        )
                    elif job_type == 'kbp':
                        running_jobs.append( \
                            pool.apply_async( \
                                kbp_job, (self.path, job_file, job_data) \
                            ) \
                        )


            if not run:
                pool.close()
                pool.join()
                f = open(self.path + 'c - stopped', 'w')
                f.close()
                return 1
            else:
                dt = datetime.datetime.now()
                if running_processes == max_processes:
                    print "Queue is full. (%i jobs running | %i jobs waiting) - %i/%i %i:%02i" \
                          % (running_processes, len(jobs), dt.month, dt.day, dt.hour, dt.minute)
                else:
                    print "waiting for new jobs.. (%i/%i slots used) - %i/%i %i:%02i" \
                          % (running_processes, max_processes, dt.month, dt.day, dt.hour, dt.minute)
                time.sleep(60)



if __name__ == '__main__':
    class job(object):
        script = 'sleep 5 \n echo "running..."\n'

#    o = job()
#    f = open('/scratch/scratch/tmeyer/md_pka/md_pka_manager/queue/d58/job1.pickle', 'w')
#    pickle.dump(o,f)
#    f.close()

#    l = listener('/scratch/scratch/tmeyer/md_pka/md_pka_manager/queue/d58')
#    l.run()

    l = listener('/scratch/scratch/tmeyer/md_pka/md_pka_manager/queue/kbp_jobs/')
    # l = listener('/scratch/scratch/tmeyer/md_pka/md_pka_manager/queue/kbp_jobs/queue2/')
    # l.run_parallel(100)
    l.run_parallel(140)






