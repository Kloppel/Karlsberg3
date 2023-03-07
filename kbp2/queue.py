# -*- coding: utf-8 -*-

from multiprocessing import Pool
import os
import cPickle as pickle


# For get_status in local jobs:
# if not rj.ready():
#     running_jobs_new.append(rj)
#     running_processes += 1
# else:
#     rj.get()


class Queue_manager(object):
    def __init__(self):
        self.jobs = []

        self.job_types = {}
        # The only jobs that are queued locally.
        self.job_types['local'] = {'max_jobs' : 3}
        # For the following jobs 'submit' and 'get_status' functions of the jobs are used.
        self.job_types['folder'] = {'max_jobs' : None}
        self.job_types['submit'] = {'max_jobs' : None}
        self.job_types['server'] = {'max_jobs' : None}

        # Prepare queues
        self.pools = {}
        for job_type in self.job_types.keys():
            max_processes = self.job_types[job_type]['max_jobs']
            self.pools[job_type] = Pool(processes=max_processes)

    def submit(self, job, wait=False):
        self.jobs.append(job)

        # Submit into queue.
        job_type = job.job_type
        pool = self.pools[job_type]
        pool.apply_async(job)

    def get_status(self):
        pass

    def submit_waiting(self):
        pass

    def terminate(self):
        for job_type in self.job_types.keys():
            pool = self.pools[job_type]
            pool.terminate()
            pool.join()

    def collect(self):
        for job_type in self.job_types.keys():
            pool = self.pools[job_type]
            pool.close()
            pool.join()


class Queue_manager_job(object):
    def __init__(self, job_type):
        self.job_type = job_type

        # 'queued', 'running' or 'done'
        self.status = ''
        self.msg = ''

    def submit(self):
        pass

    def get_status(self):
        pass

    def set_status(self):
        pass

    def __call__(self, *args, **kwargs):
        pass


class Folder_job(Queue_manager_job):
    def __init__(self, job_type, folder, name, prefix='', append=False):
        super(Folder_job, self).__init__(job_type)

        if folder[-1] != '/':
            folder += '/'

        self.folder = folder
        self.name = name
        self.prefix = prefix

        if not os.path.exists(folder):
            error = "The folder specified for the queue does not exist: %s" % folder
            raise AssertionError(error)

        if append:
           self.filename = self.__get_free_file_name()
        else:
            self.filename = self.__get_file_name()

        self.job_object = None

    def submit(self):
        filename = self.filename
        if not self.get_status():
            self.__dump_job(filename)

    def __get_file_name(self):
        prefix = self.prefix
        base_name = self.name
        suffix = '.pickle'
        return self.folder + prefix + base_name + suffix

    def __get_free_file_name(self):
        prefix = self.prefix
        base_name = self.name
        suffix = '.pickle'
        evading_suffix = self.__get_evading_suffix()
        return self.folder + prefix + base_name + evading_suffix + suffix

    def __get_evading_suffix(self):
        prefix = self.prefix
        base_name = self.name
        suffix = '.pickle'
        if not os.path.exists(self.folder + prefix + base_name + suffix):
            return ''
        else:
            c = 0
            while c < 9999:
                if c == 0:
                    evading_suffix = ''
                else:
                    evading_suffix = '_%i' % c
                filename = self.folder + prefix + base_name + evading_suffix + suffix
                if not os.path.exists(filename):
                    break
                c += 1
            return evading_suffix

    def get_status(self):
        if os.path.exists(self.filename):
            status = 'queued'
            if os.path.exists(self.filename + ' - running'):
                status = 'running'
            if os.path.exists(self.filename + ' - done'):
                status = 'done'
        else:
            status = ''
        return status

    def __dump_job(self, filename):
        f = open(filename, 'w')
        pickle.dump(self.job_object, f)
        f.close()

    def set_job_object(self, job_object):
        self.job_object = job_object


if __name__ == '__main__':
    # queue_folder = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_snase/test'
    queue_folder = '/scratch/scratch/tmeyer/md_pka/md_pka_manager/queue/kbp_jobs/'
    job = Folder_job('folder', queue_folder, 'protein', prefix='kbp_job_')



    # Start the Karlsberg+ job.
    kbp_pdb = '/scratch/scratch/pdb/pdb_bio_merged/bd/3bdc.pdb1'
    kbp_run_dir = '/scratch/scratch/tmeyer/md_pka/runs/kbp_cache_snase/run/single/'

    from kbp2 import job_manager2
    kbp_job = job_manager2.kbp_job(pdb=kbp_pdb)
    kbp_job.kb_para['globals']['workDir'] = '/public/scratch/tmeyer/kb/pre_md/'
    kbp_job.parameter['run_dir'] = kbp_run_dir

    job.set_job_object(kbp_job)



    job.submit()
    print job.get_status()



# job_filename = self.kbp_queue + "kbp_job_%s_%s_%s.pickle" \
#                % (self.kbp_folder_suffix, self.jobname, pdb[:-4])

# nur locale Jobs so starten
# bei allen anderen 'submit' aufrufen


# queue_manager = Queue_manager()
#
# for i in range(6):
#     job = Queue_manager_job('local')
#     queue_manager.submit(job)
#
# queue_manager.collect()

# job_status = queue_manager.submit(job)







