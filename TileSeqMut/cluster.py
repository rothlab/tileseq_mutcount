#!/usr/bin/env python3.7

import pandas as pd
import os
import glob
import re
import subprocess
import time
import numpy as np

# other modules
from TileSeqMut import alignment

# The functions in this script are used to manage submitting jobs to different clusters
# This script is called by main.py

# GURU is depreciated
def alignment_sh_guru(fastq_map, ref_name, ref_seq, ref_path, sam_path, ds_sam_path, sh_output):
    """
    fastq_map: df contains paths to fastq files and downsamled fastq files
    ref_name: name for the reference sequence (same as project name)
    ref_seq: reference sequence to put in fasta file
    ref_path: build fasta in this path
    sam_path: path to sam files
    ds_sam_path: path to downsampled sam files
    sh_output: make sh files to submit to the cluster

    make reference for bowtie2
    make SGE sh files
    submit sh files to SGE
    """
    ## build reference
    ref = alignment.make_ref(ref_name, ref_seq, ref_path)

    for index, row in fastq_map.iterrows():
        sample_name = os.path.basename(row["R1"]).split("_")[0]

        shfile = os.path.join(sh_output, sample_name+"_aln.sh") # for each sample, the alignment is for both R1 and R2 (they are aligning separately)
        log_file = alignment.align_main(ref, row["R1"], row["R2"], sam_path, shfile)

        sub_cmd = f"qsub -cwd -N {'aln_'+sample_name} -e {log_file} {shfile}"
        os.system(sub_cmd)
        shfile_ds = os.path.join(sh_output, sample_name+"_aln_ds.sh")
        log_file_ds = alignment.align_main(ref, row["r1_ds"], row["r2_ds"], ds_sam_path, shfile_ds)
        sub_cmd = f"qsub -cwd -N {'aln_ds_'+sample_name} -e {log_file_ds} {shfile_ds}"
        os.system(sub_cmd)


def alignment_sh_galen(fastq_map, ref_name, ref_seq, ref_path, sam_path, sh_output, at, logging, rc, blacklist, queue):
    """
    fastq_map: df contains paths to fastq files and downsamled fastq files
    ref_name: name for the reference sequence (same as project name)
    ref_seq: reference sequence to put in fasta file
    ref_path: build fasta index in this path
    sam_path: path to sam files
    sh_output: make sh files to submit to the cluster
    at: alignment time, default to 8 hours
    logging: logging object
    rc: True if the user wants to use reverse complement as well

    make reference for bowtie2
    make SGE sh files
    submit sh files to SGE
    """
    ## build references
    ref = alignment.make_ref(ref_name, ref_seq, ref_path)
    phix = alignment.make_ref("phix", ref_seq, ref_path)

    if blacklist != "": 
        blacklistArg = f"#SBATCH --exclude={blacklist}\n"
    else:
        blacklistArg = ""

    if queue != "":
        queueArg = f"#SBATCH --partition={queue}\n"
    else:
        queueArg = ""

    # store sam paths
    fastq_map = pd.concat([fastq_map, pd.DataFrame(columns=["r1_sam", "r2_sam"])])
    all_job_id = []
    for index, row in fastq_map.iterrows():
        sample_name = os.path.basename(row["R1"]).split("_")[0]

        # for each sample, the alignment is for both R1 and R2 (they are aligned separately)
        shfile = os.path.join(sh_output, sample_name+"_aln.sh")

        # write header to sh file

        # create log file for alignment
        sam_log_f = os.path.join(sam_path, f"{sample_name}")

        if "Undetermined" in sample_name: # phix takes longer to align
            time_request = f"36:00:00"
            header = f"#!/bin/bash\n#SBATCH --time={time_request}\n#SBATCH --mem=2G\n#SBATCH --job-name=aln_{sample_name}\n#SBATCH " \
                 f"--error={sam_log_f}-%j.log\n{blacklistArg}{queueArg}#SBATCH --output={sam_log_f}-%j.log\n"
            # when align undetermined fastq files to phix, we consider reads in both direction, rc = True
            r1_sam, r2_sam, log_file = alignment.align_main(phix, row["R1"], row["R2"], sam_path, shfile, rc=True,
                                                            header=header)
        else:
            time_request = f"0{at}:00:00"
            header = f"#!/bin/bash\n#SBATCH --time={time_request}\n#SBATCH --mem=2G\n#SBATCH --job-name=aln_{sample_name}\n#SBATCH " \
                 f"--error={sam_log_f}-%j.log\n{blacklistArg}{queueArg}#SBATCH --output={sam_log_f}-%j.log\n"
            r1_sam, r2_sam, log_file = alignment.align_main(ref, row["R1"], row["R2"], sam_path, shfile, rc=rc, header=header)

        row["r1_sam"] = r1_sam
        row["r2_sam"] = r2_sam
        sub_cmd = ["sbatch", str(shfile)]
        jobs = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
        job_id = jobs.stdout.decode("utf-8").strip().split(".")[0]
        all_job_id.append(job_id.split()[-1])
        logging.info(f"{sample_name}: job id - {job_id}")
    return fastq_map, all_job_id


def alignment_sh_ccbr(fastq_map, ref_name, ref_seq, ref_path, sam_path, sh_output, at, logging, rc, cluster_name):
    """
    submit jobs to BC/BC2/DC
    return a df with columns: [R1, R2, r1_sam, r2_sam]
    return a list of job id that we just submited
    """

    # build reference
    ref = alignment.make_ref(ref_name, ref_seq, ref_path)
    phix = alignment.make_ref("phix", ref_seq, ref_path)
    # store sam paths
    fastq_map = pd.concat([fastq_map, pd.DataFrame(columns=["r1_sam", "r2_sam"])])
    time = at
    all_job_id = []
    for index, row in fastq_map.iterrows(): # go through all the fastq pairs
        sample_name = os.path.basename(row["R1"]).split("_")[0]

        shfile = os.path.join(sh_output, f"Aln_{sample_name}.sh")
        if "Undetermined" in sample_name:
            # when align undetermined fastq files to phix, we consider reads in both direction, rc = True
            r1_sam, r2_sam, log_file = alignment.align_main(phix, row["R1"], row["R2"], sam_path, shfile, rc=True)
        else:
            r1_sam, r2_sam, log_file = alignment.align_main(ref, row["R1"], row["R2"], sam_path, shfile, rc=rc)

        row["r1_sam"] = r1_sam
        row["r2_sam"] = r2_sam
        # create log file for alignment
        sam_log_f = os.path.join(sam_path, f"{sample_name}.log")
        if cluster_name == "BC2" or cluster_name == "BC":
            sub_cmd = ["submitjob2", "-w", str(time), "-c", "1", str(shfile), "2>", sam_log_f]
        else:
            sub_cmd = ["submitjob", "-w", str(time), "-c", "1", str(shfile), "2>", sam_log_f]
        jobs = subprocess.run(sub_cmd, stdout=subprocess.PIPE)

        job_id = jobs.stdout.decode("utf-8").strip().split(".")[0]
        all_job_id.append(job_id)
        logging.info(f"{sample_name}: job id - {job_id}")

    return fastq_map, all_job_id



def mut_count_sh_ccbr(sample_name, cmd, mt, mm, sh_output_dir, logger, cores, cluster_name):
    """
    Submit mutation count jobs to BC

    """
    # go through files df and submit jobs for each pair of sam files
    # counting mutations in raw sam output files
    shfile = os.path.join(sh_output_dir, f"Mut_count_{sample_name}.sh")
    log_f = os.path.join(sh_output_dir, f"Mut_count_{sample_name}.log")
    with open(shfile, "w") as sh:
        sh.write(cmd+"\n")
        os.system(f"chmod 755 {shfile}")
    # submit this to the cluster
    if cluster_name == "BC" or cluster_name == "BC2":
        sub_cmd = ["submitjob2","-w", str(mt), "-c", f"{cores}", "-m", f"{mm}", shfile, "&>>", log_f]
    else:  # DC
        sub_cmd = ["submitjob","-w", str(mt), "-c", f"{cores}", "-m", f"{mm}", shfile, "&>>", log_f]
    logger.debug(sub_cmd)
    job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
    job_id = job.stdout.decode("utf-8").strip().split(".")[0]
    # log sample name and job id
    logger.info(f"Sample {sample_name}: job id - {job_id}")
    return job_id


def mut_count_sh_galen(sample_name, cmd, mt, mm, sh_output_dir, logger, cores, blacklist, queue):
    """
    Submit mutation count jobs to DC
    """
    # go through files df and submit jobs for each pair of sam files
    # counting mutations in raw sam output files

    shfile = os.path.join(sh_output_dir, f"Mut_count_{sample_name}.sh")
    log_f = os.path.join(sh_output_dir, f"Mut_count_{sample_name}")

    time_request = f"{mt}:00:00"

    if blacklist != "": 
        blacklistArg = f"#SBATCH --exclude={blacklist}\n"
    else:
        blacklistArg = ""

    if queue != "":
        queueArg = f"#SBATCH --partition={queue}\n"
    else:
        queueArg = ""


    header = f"#!/bin/bash\n#SBATCH --time={time_request}\n#SBATCH --job-name=mtc_{sample_name}\n#SBATCH " \
             f"--cpus-per-task={cores}\n#SBATCH --error={log_f}-%j.log\n#SBATCH --mem={mm}G\n#SBATCH " \
             f"--output={log_f}-%j.log\n{blacklistArg}{queueArg}"
    
    # TODO: Delete existing .sh file for dropped samples and create new ones
    
    if os.path.isfile(shfile):
        os.remove(shfile)

    with open(shfile, "w") as sh:
        sh.write(header)
        sh.write(cmd+"\n")
        os.system(f"chmod 755 {shfile}")
    #sample_error_file = os.path.join(log_dir, f"sample_{sample_name}.log")
    # submit this to the cluster
    sub_cmd = ["sbatch", str(shfile)]
    logger.debug(sub_cmd)
    job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
    job_id = job.stdout.decode("utf-8").strip()
    # log sample name and job id
    logger.info(f"Sample {sample_name}: job id - {job_id}")
    return job_id.split()[-1]


#######################################################################################################################################################################################

def submit_mutcount_jobs(output_dir, param_json, sample, blacklist, sh_output, r1, r2, args, logger):
    """
    Make mutation counting jobs and submit for this sample
    """
    # submit job with main.py -r1 and -r2
    # run main.py with -r1 and -r2

    cmd = f"tileseq_mut -n {args.name} -r1 {r1} -r2 {r2} -o {output_dir} -p" \
                f" {param_json} --skip_alignment -log {args.log_level} -env {args.environment} -at {args.at} -mt {args.mt} -c {args.c} --logger_time {args.logger_time} " 
    if args.sr_Override:
        cmd = cmd + "--sr_Override "

    if args.wt_override:
        cmd = cmd + "--wt_override "
    
    if args.calibratePhredPhix:
        cmd = cmd + "--calibratePhredPhix "

    if args.calibratePhredWT:
        cmd = cmd + "--calibratePhredWT "

    if args.errorOverride:
        cmd = cmd + "--errorOverride "

    if args.environment == "BC2" or args.environment == "BC" or args.environment == "DC":
        logger.info("Submitting mutation counts jobs....")
        job_id = mut_count_sh_ccbr(sample, cmd, args.mt, args.mm, sh_output, logger,
                                            args.c, args.environment)
    # elif :
    #     logging.info("Submitting mutation counts jobs to DC...")
    #     # (sample_name, cmd, mt, sh_output_dir, logger)
    #     job_id = cluster.mut_count_sh_dc(sample, cmd, self._args.mt, self._args.mm, sh_output, self._log,
    #                                      self._args.c)
    elif args.environment == "GALEN":
        logger.info("Submitting mutation counts jobs to GALEN...")
        # (sample_name, cmd, mt, sh_output_dir, logger)
        job_id = mut_count_sh_galen(sample, cmd, args.mt, args.mm, sh_output, logger,
                                            args.c, blacklist, args.queue)  # this
    else:
        raise ValueError("Wrong environment")

    return job_id


def makejobs(output_dir, param_json, samples, blacklist, sh_output, sam_dir, args, logger):
        """
        For each pair of sam files in output/sam_files/
        submit mut count job to the cluster
        """
        if sam_dir == "":
            # get samples in param file
            sam_dir = os.path.join(output_dir, "sam_files/") # read sam file from sam_file
        # get sam files from parameter file
        if not os.path.isdir(sam_dir):
            logger(f"Directory: ./sam_files/ not found in {output_dir}")
            raise ValueError()

        logger.debug(f"Sam files are read from {sam_dir}")
        job_list = []
        for sample in samples:
            # assume all the sam files have the same name format (id_*.sam)
            sam_f_r1 = glob.glob(f"{sam_dir}/{sample}_*R1_*.sam")
            sam_f_r2 = glob.glob(f"{sam_dir}/{sample}_*R2_*.sam")
            if len(sam_f_r1) == 0 or len(sam_f_r2) == 0:
                logger.error(f"SAM file for sample {sample} not found. Please check your parameter file")
                raise ValueError()

            else:
                logger.info(f"Sample {sample}")
                logger.info(f"Read1: {sam_f_r1[0]}")
                logger.info(f"Read2: {sam_f_r2[0]}")
                sam_id_1 = os.path.basename(sam_f_r1[0]).split("_")[0]
                sam_id_2 = os.path.basename(sam_f_r2[0]).split("_")[0]
                if (sam_id_1 != sample) or (sam_id_1 != sam_id_2) or (sam_id_2 != sample):
                    logger.error("IDs in sam files don't match!")
                r1 = sam_f_r1[0]
                r2 = sam_f_r2[0]
            
            job_id = submit_mutcount_jobs(output_dir, param_json, sample, blacklist, sh_output, r1, r2, args, logger)
            job_list.append(job_id)
        logger.info("all samples resubmitted")


#######################################################################################################################################################################################
        # after resubmitting dropped jobs, track all resubmitted jobs on cluster
        jobs = ",".join(job_list)
        logger.debug(f"All jobs resubmitted: {jobs}")

#######################################################################################################################################################################################

def parse_jobs(job_list, env, logger):
    """
    return true if all the jobs in job list finished
    else wait for 10 mins and return how man jobs are running and queued
    job_list: list of job ids
    logger: logging object
    """
    qstat_cmd = ["qstat"] + job_list
    job = subprocess.run(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    qstat_out = job.stdout.decode("utf-8", errors="replace")
    qstat_err = job.stderr.decode("utf-8", errors="replace")

    f_id = []
    updated_list = []
    while True:
        running = []
        queued = []
        completed = []
        # jobs might be finished and no longer in queue
        # in this case the job ID are returned in sterr
        # extract job id from std err
        err = []
        if qstat_err != "":
            if env == "BC" or env == "BC2":
                err = qstat_err.split("\n")[:-1]
                id_regex = re.compile(r"(\d+).bc")
            elif env == "DC":
                err = qstat_err.split("\n")[:-1]
                id_regex = re.compile(r"(\d+).dc[0-9]+")

            f_id = [] # finished jobs
            for i in err:
                try:
                    match = id_regex.search(i)
                    job_id = match.group(1)
                    f_id.append(job_id)
                except:
                    logger.warning(i)
                    continue
            err_id = set(f_id)
            updated_list = [x for x in job_list if x not in err_id]

        if qstat_out != "":
            qstat_out = qstat_out.split("\n")[:-1]
            # R: running
            # Q: Queued
            # C: completed
            # E: error
            if env == "BC" or env == "BC2":
                id_regex = re.compile(r"(\d+).bc.+(R|Q|C|E)")
            elif env == "DC":
                id_regex = re.compile(r"(\d+).dc[0-9]+.+(R|Q|C|E)")

            for line in qstat_out:
                if ("---" in line) or ("Job ID" in line): continue
                match = id_regex.search(line)
                job_id = match.group(1)
                job_s = match.group(2)
                if job_s == "E" or job_s == "C":
                    completed.append(job_id)
                elif job_s == "R":
                    running.append(job_id)
                elif job_s == "Q":
                    queued.append(job_id)

        logger.info(f"{len(queued)} jobs queued")
        logger.info(f"{len(running)} jobs running")
        final_list = list(set(updated_list+running+queued))
        if final_list == []:
            return True
        else:
            # check in 10min
            time.sleep(600)
            job_list = final_list
            qstat_cmd = ["qstat"] + job_list
            job = subprocess.run(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            qstat_out = job.stdout.decode("utf-8", errors="replace")
            qstat_err = job.stderr.decode("utf-8", errors="replace")


def parse_jobs_galen(param_json, output_dir, sample_list, args, job_list, sleep_time, logger):
    """
    Galen uses slurm scheduler, different from BC and DC
    return true if all the jobs in job list finished
    else wait for 10 mins and return how man jobs are running and queued
    sample_list: list of sample ids
    job_list: list of job ids
    time: time (in seconds) to wait before checking jobs status
    logger: logging object
    """
    
    cmd = "squeue -u $USER -j {}".format(",".join(job_list))
    job = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # get status of jobs from job_list
    squeue_output = job.stdout.read().decode("utf-8", errors="replace") 
    squeue_error = job.stderr.read().decode("utf-8", errors="replace")
    #FIXME: Zombie process unless using job.communicate(): https://stackoverflow.com/questions/2760652/how-to-kill-or-avoid-zombie-processes-with-subprocess-module
    logger.debug(squeue_output)
    logger.debug(squeue_error)
    
    node_list = []

    while True:
        running_jobs, queued_jobs, completed_jobs, stuck_jobs, suspended_jobs = [], [], [], [], []

        if squeue_output != "":
            job_entries_df = pd.DataFrame([i.split() for i in squeue_output.split("\n")])
            job_entries_df = job_entries_df.rename(columns=job_entries_df.iloc[0]).dropna()
            logger.debug(job_entries_df)

            running_jobs = list(zip(job_entries_df[job_entries_df["ST"] == "R"]["JOBID"].tolist(), job_entries_df[job_entries_df["ST"] == "R"]["NAME"].tolist()))
            # running_jobs = job_entries_df[job_entries_df["ST"] == "R"]["JOBID"].tolist() # active job IDs
            queued_jobs = list(zip(job_entries_df[job_entries_df["ST"] == "PD"]["JOBID"].tolist(), job_entries_df[job_entries_df["ST"] == "PD"]["NAME"].tolist()))
            # queued_jobs = job_entries_df[job_entries_df["ST"] == "PD"]["JOBID"].tolist() # pending job IDs
            # completed_jobs = list(zip(job_entries_df[job_entries_df["ST"] == "CG"]["JOBID"].tolist(), job_entries_df[job_entries_df["ST"] == "CG"]["NAME"].tolist()))
            # completed_jobs = job_entries_df[job_entries_df["ST"] == "CG"]["JOBID"].tolist() # completing job IDs
            suspended_jobs = list(zip(job_entries_df[job_entries_df["ST"] == "CG"]["JOBID"].tolist(), job_entries_df[job_entries_df["ST"] == "CG"]["NAME"].tolist()))
            # suspended_jobs = job_entries_df[job_entries_df["ST"] == "S"]["JOBID"].tolist() # suspended job IDs
            stuck_jobs = list(zip(job_entries_df[job_entries_df["NODELIST(REASON)"].str.contains("launch failed requeued held")]["JOBID"].tolist(), job_entries_df[job_entries_df["NODELIST(REASON)"].str.contains("launch failed requeued held")]["JOBID"].tolist()))
            # stuck_jobs = job_entries_df[job_entries_df["NODELIST(REASON)"].str.contains("launch failed requeued held")]["JOBID"].tolist() # stuck job IDs

        logger.info(f"{len(queued_jobs)} jobs queued")
        logger.info(f"{len(running_jobs)} jobs running\n")
        # logger.info(f"{len(completed_jobs)} jobs completed\n")

        #attempt to requeue suspended jobs
        if len(suspended_jobs):
            logger.info(f"WARNING: Suspended jobs detected! {len(suspended_jobs)} jobs suspended_jobs")
            logger.info("Requeuing jobs {}...".format(" ".join(suspended_jobs)))
            requeue_cmd = "scontrol requeue {}".format(" ".join(suspended_jobs))

        #attempt to release stuck jobs
        if len(stuck_jobs):
            logger.info(f"WARNING: Failed/Held jobs detected! {len(stuck_jobs)} jobs stuck")
            logger.info("Attempting to release jobs {}...".format(" ".join(stuck_jobs)))
            release_cmd = "scontrol release {}".format(" ".join(stuck_jobs))

        job_list = running_jobs + queued_jobs + suspended_jobs + stuck_jobs
        if not len(job_list):
            return True
        
        job_samples_dict = {}
        for job, name in job_list: #populate a dictionary with job ids and sample names
            # logger.debug(f"job name: {name}")
            job_samples_dict[job] = re.search("(\d+)", name).group(1)
        # samples that have finished running     
        completed_jobs = list(set(sample_list)-set(list(job_samples_dict.values()))) 
        
        if completed_jobs:
            #if submitted jobs are alignment jobs
            if re.search("^(\w+)_", job_list[0][1]).group(1) == "aln": 
                # TODO: check the SAM directory for completed job's .SAM files
                logger.debug("Skipping alignment resubmissions")
            #if submitted jobs are mutcount jobs
            elif re.search("^(\w+)_", job_list[0][1]).group(1) == "mtc": 
                mutcount_list = glob.glob(os.path.join(output_dir, "counts_sample_*.csv"))
                # failed_samples = []
                all_samples_with_csv = []
                # go through all the mut count files
                for f in mutcount_list:
                    m = re.search('.*counts_sample_(.+?).csv', f)
                    # find sample ID in file name
                    if m:
                        sample_id = m.group(1)
                        all_samples_with_csv.append(sample_id)
                    else:
                        logger.error("Sample ID not found")
                        raise ValueError(f"Sample ID not found for {f}")
                    # check if the file is empty
                    # mut_n = self._checkoutput(f)
                    # if mut_n == 0:
                    #     failed_samples.append(sample_id)
                # if any sample does not have a csv file generated (jobs were completed but did not generate output files)
                dropped_samples = [i for i in completed_jobs if i not in all_samples_with_csv]
                # failed_samples += other_missing
                if dropped_samples:   
                    logger.info(f"Failed samples: {dropped_samples}")
                    logger.info("Resubmitting samples ...")
                    env_jobs_dir = os.path.join(output_dir, "GALEN_jobs")
                    parent_dir = os.path.abspath(os.path.join(output_dir, os.pardir))
                    sam_dir = os.path.join(parent_dir, "sam_files")

                    
                    dropped_job_list = [job_id for job_id, sample in job_samples_dict.items() if sample in dropped_samples]
                    

                    # determine failing nodes and resubmit sample
                    for job in dropped_job_list: 
                        cmd = ["sacct", "-j", job, "--format=NodeList"]
                        job = subprocess.run(cmd, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        sacct_output = job.stdout.read().decode("utf-8", errors="replace") 
                        node_df = pd.DataFrame([i.split() for i in sacct_output.split("\n")])
                        node_df = node_df.iloc[2:]
                        node = node_df.iloc[0, 0]
                        node_list.append(node)
                    blacklist = ",".join(node_list)

                    #resubmit failed jobs
                    makejobs(output_dir, param_json, dropped_samples, blacklist, env_jobs_dir, sam_dir, args, logger)

                    # finished = self._makejobs(env_jobs_dir, sam_dir, failed_samples)
                
        #check in 10 mins
        time.sleep(int(sleep_time))
        job_list = [job for job, name in job_list]
        cmd = "squeue -u $USER -j {}".format(",".join(job_list))
        job = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # get status of remaining jobs from job_list
        squeue_output = job.stdout.read().decode("utf-8", errors="replace")
        logger.debug(squeue_output)


def submit_given_jobs(shfile, logger, mt, mm, cores, env=""):
    """
    for a given sh file, submit to cluster, return job id
    """
    sample_name = shfile.split(".")[0].split("_")[-1]
    sh_dir = os.path.dirname(shfile)
    log_f = os.path.join(sh_dir, f"Mut_count_{sample_name}.log")
    if env == "GALEN":
        sub_cmd = ["sbatch", str(shfile)]
        job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
        job_id = job.stdout.decode("utf-8").strip().split()[-1]
    elif env == "BC2" or env == "BC":
        sub_cmd = ["submitjob2", "-w", str(mt), "-c", f"{cores}", "-m", f"{mm}", shfile, "&>>", log_f]
        job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
        job_id = job.stdout.decode("utf-8").strip()
    elif env == "DC":
        sub_cmd = ["submitjob", "-w", str(mt), "-c", f"{cores}", "-m", f"{mm}", shfile, "&>>", log_f]
        job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
        job_id = job.stdout.decode("utf-8").strip()
    else:
        raise ValueError("Wrong environment code")
    # log sample name and job id
    logger.info(f"Sample {sample_name}: job id - {job_id}")
    return job_id

#
# if __name__ == "__main__":
#     # test job list
#     job_list = ["352344", "348556"]
#     # parse_jobs(job_list, "DC", "")
#
#     parse_jobs_galen(job_list, "")
