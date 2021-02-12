from subprocess import check_output
import time
import subprocess
import os


LSF_TEMPLATE = 'bsub -o {}.out -e {}-err.out -M {} -n {} '


def memory_supervised_lsf_run(job, prefix, workdir,
                              start_mem=1000, mem_delta=1000, mem_limit=80000, np=1):
    """
    Submits a job to the cluster LSF queue system. If the job is killed due
    to exceeded memory then the memory limit is increased.
    :param np: number of cores
    :param start_mem: default memory
    :param job: job command or bash file be run on cluster
    :param prefix: a string with the name of the program (used as prefix for stderr and stdout)
    :param workdir: directory where the job will be run
    :param mem_delta: memory increase steep
    :param mem_limit: maximum available memory
    """

    def run_monitored_job(command, work_dir):
        """
        Submits a job to the queue system and waits until it finishes
        :param command: Bsub command
        :param work_dir: directory where the job will be run
        """
        command = command.split()
        # Submit job to the queue ang get the job id
        out = str(check_output(command, cwd=work_dir))
        job_id = out.split('<')[1].split('>')[0]
        while True:
            # Check the job status every 20 seconds
            try:
                out1 = str(check_output(["bjobs", "-a", "-noheader", job_id]))
                status = out1.split()[2]
            except subprocess.CalledProcessError:
                continue
            if status == 'DONE' or status == 'EXIT':
                break
            elif 'not found' in out1:
                break
            else:
                time.sleep(20)

    def read_max_used_mem(bsub_out):
        """
        Reads the maximum and exceeded amount of memory during the last bsub run.
        :param bsub_out: bsub output file
        :return: maximum used memory, exceeded amount
        """
        memterm = False
        last = -1
        maxmem = 0
        try:
            with open(bsub_out) as f:
                lines = f.readlines()
                for nr, line in enumerate(lines):
                    if 'LSBATCH' in line:
                        last = nr
                for line in lines[last:]:
                    if 'TERM_MEMLIMIT' in line:
                        memterm = True
                if memterm:
                    for line in lines[last:]:
                        try:
                            if 'Max Memory' in line:
                                maxmem = int(float(line.split()[-2]))
                        except ValueError:
                            pass
        except IOError:
            pass
        return memterm, maxmem

    bsub_prefix = build_bsub_command(prefix, start_mem, np)
    jobstring = bsub_prefix + ' ' + job
    run_monitored_job(jobstring, workdir)
    current_mem = start_mem
    while True:
        # Check the job summary. If the memory limit was exceeded rerun with more memory
        mem_term, max_mem = read_max_used_mem(os.path.join(workdir, prefix + '.out'))
        if not mem_term:
            break
        elif max_mem > mem_limit:
            break
        elif max_mem == 0:
            current_mem += mem_delta
        elif max_mem < mem_limit:
            current_mem = max_mem + mem_delta

        bsub_prefix = build_bsub_command(prefix, current_mem, np)
        jobstring = bsub_prefix + ' ' + job
        run_monitored_job(jobstring, workdir)
    # Sometimes a file created by a node on cluster is not immediately visible to other cores
    # which results in os.in_dir.exists returning False on existing files. Waiting a bit helps
    time.sleep(10)
    return current_mem


def build_bsub_command(prefix, memory, np=1):
    """
    Creates the cluster job submission command
    Specific to LSF platform (needs to be modified for other cluster platforms
    :param np: number of cores
    :param prefix: prefix for the LSF output files
    :param memory: Memory limit
    :return: LSF command
    """
    command = LSF_TEMPLATE.format(prefix, prefix, memory, np)
    return command