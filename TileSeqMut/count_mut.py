#!/usr/bin/env python3.7

import pandas as pd
import re
import os
import math
import logging
import mmap
import time
import glob
import multiprocessing as mp
import argparse
# import datetime
from collections import deque

# modules in package
from TileSeqMut import help_functions
from TileSeqMut import locate_mut

## Read sam file (R1 and R2)
# count mutations in sam files
# output mutation counts


class readSam(object):

    def __init__(self, sam_r1, sam_r2, param, arguments, output_dir, cores, loggerobj):
        """
        sam_R1: read one of the sample
        sam_R2: read two of the sample
        param: parameter sheet
        arguments: user input arguments
        output_dir: output directory
        cores: number of processes to use
        loggerobj:  logger
        """
        self._r1 = sam_r1
        self._r2 = sam_r2
        self._param = param
        self._project, self._seq, self._cds_seq, self._tile_map, self._region_map, self._samples, self._var, self._relations = help_functions.parse_json(param)
        # settings
        self._qual = self._var["posteriorThreshold"]
        self._mutrate = self._var["mutRate"]
        self._output_counts_dir = output_dir
        self._cores = cores

        # sample information
        self._sample_id = os.path.basename(sam_r1).split("_")[0]
        self._sample_info = self._samples[self._samples["Sample ID"] == self._sample_id]

        # tile information
        self._sample_tile = self._sample_info["Tile ID"].values[0]
        # beginning position of this tile (cds position)
        self._tile_begins = (self._tile_map[self._tile_map["Tile Number"] == self._sample_tile]["Start AA"].values[0] *3)-2
        # ending position of this tile (cds position)
        self._tile_ends = self._tile_map[self._tile_map["Tile Number"] == self._sample_tile]["End AA"].values[0] *3
        self._tile_len = self._tile_ends - self._tile_begins
        # min cover is the minimum percentage of the tile that needs to be covered by a read
        # for example, mincover=0.9 means 90 percent of the tile needs to be covered
        min_cover = self._var["minCover"]
        self._min_map_len = math.ceil(self._tile_len * float(min_cover))
        self._cds_start = self._seq.cds_start.iloc[0]

        # get sample specific info
        self._sample_condition = self._sample_info["Condition"].values[0]
        self._sample_tp = self._sample_info["Time point"].values[0]
        self._sample_rep = int(self._sample_info["Replicate"].values[0])
        # seq_lookup is a dataframe with the following columns:
        # temp_pos, cds, dna_pos, protein, protein_pos
        self._seq_lookup = help_functions.build_lookup(self._cds_start, self._seq.cds_end.iloc[0], self._cds_seq)

        # path for the output files
        self._sample_counts_f = os.path.join(self._output_counts_dir,f"counts_sample_{self._sample_id}.csv")
        # r1 and r2 files are only used when posteriorQC is enabled
        self._sample_counts_r1_f = os.path.join(self._output_counts_dir, f"counts_sample_{self._sample_id}_r1.csv")
        self._sample_counts_r2_f = os.path.join(self._output_counts_dir, f"counts_sample_{self._sample_id}_r2.csv")
        # file used to track reads mapped outside of the tile
        self._sample_offreads_f = os.path.join(self._output_counts_dir,f"offreads_{self._sample_id}.csv")
        offreads = open(self._sample_offreads_f, "w")
        offreads.write("r1_start, r1_len, r2_start, r2_len\n")
        offreads.close()

        self._error_override = arguments.errorOverride
        self._base = arguments.base
        self._posteriorQC = arguments.posteriorQC
        self._mut_log = loggerobj
        self._mut_log.info(f"Counting mutations in sample-{self._sample_id}")
        self._mut_log.info(f"Sam file input R1:{sam_r1}")
        self._mut_log.info(f"Sam file input R2:{sam_r2}")

        output_csv = open(self._sample_counts_f, "w")
        # write log information to counts output
        output_csv.write(f"#Sample:{self._sample_id}\n#Tile:{self._sample_tile}\n#Tile Starts:{self._tile_begins}\n"
                         f"#Tile Ends:{self._tile_ends}\n#Condition:{self._sample_condition}\n"
                         f"#Replicate:{self._sample_rep}\n#Timepoint:{self._sample_tp}\n"
                         f"#Posterior cutoff:{self._qual}\n#min cover %:{min_cover}\n")
        output_csv.close()

        if self._posteriorQC:
            output_csv = open(self._sample_counts_r1_f, "w")
            # write log information to counts output
            output_csv.write(f"#Sample:{self._sample_id}\n#Tile:{self._sample_tile}\n#Tile Starts:{self._tile_begins}\n"
                             f"#Tile Ends:{self._tile_ends}\n#Condition:{self._sample_condition}\n"
                             f"#Replicate:{self._sample_rep}\n#Timepoint:{self._sample_tp}\n"
                             f"#Posterior cutoff:{self._qual}\n#min cover %:{min_cover}\n")
            output_csv.close()

            output_csv = open(self._sample_counts_r2_f, "w")
            # write log information to counts output
            output_csv.write(f"#Sample:{self._sample_id}\n#Tile:{self._sample_tile}\n#Tile Starts:{self._tile_begins}\n"
                             f"#Tile Ends:{self._tile_ends}\n#Condition:{self._sample_condition}\n"
                             f"#Replicate:{self._sample_rep}\n#Timepoint:{self._sample_tp}\n"
                             f"#Posterior cutoff:{self._qual}\n#min cover %:{min_cover}\n")
            output_csv.close()

        # build a df for tracking how many reads passed filter at certain position of the read 
        # df contains ["pos", "m_r1", "m_r2", "passed", "wt_passed"] columns
        # "pos" refers to all the nt positions in this tile 
        # "m_r1" refers to how many times the mutation was seen r1 
        # "m_r2" refers to how many times the mutation was seen r2 
        # "m_both" refers to how many times the mutation was seen in both r1 and r2 
        # "passed" refers to how many mut calls passed the filter
        # "wt_passed" refers to how many WT calls passed the filter
        self._track_reads = pd.DataFrame({}, columns=["pos", "m_r1", "m_r2", "m_both", "passed", "wt_passed"])
        
        self._track_reads["pos"] = range(self._tile_begins, self._tile_ends+1)
        self._track_reads = self._track_reads.set_index("pos")
        self._track_reads = self._track_reads.fillna(0)

        # single replicate override
        self._sroverride = arguments.sr_Override

    def adjust_er(self, wt_override=False):
        """
        Based on the wt samples given in the parameter sheet,
        calculate new error rate for phred scores
        this function is dependent on the tileseqMave package
        please make sure the correct script is installed before running the pipeline
        """
        r_df = pd.DataFrame(self._relations)
        phred_output = []
        
        if r_df.empty:
            if not wt_override:
                # no wt relationship is defined
                self._mut_log.warning("No WT relationship! No error rate correction can be applied")
                return phred_output

            else:  # treating EVERYTHING in the run as wt
                self._mut_log.warning("WT OVERRIDE ON, TREATING EVERYTHING AS WT")
                wt = pd.DataFrame({}, columns=["Condition 1", "Relationship", "Condition 2"])
                wt["Condition 1"] = self._samples["Condition"]
        else:
            # find out wt sample ID for this sample
            wt = r_df[r_df["Relationship"]=="is_wt_control_for"]
            if wt.empty:
                if not wt_override:
                    # no wt relationship is defined
                    self._mut_log.warning("No WT relationship! No error rate correction can be applied")
                    return phred_output

                else:  # treating EVERYTHING in the run as wt
                    self._mut_log.warning("WT OVERRIDE ON, TREATING EVERYTHING AS WT")
                    wt["Condition 1"] = self._samples["Condtion"]
        wt_id = ""

        if self._sample_condition in wt["Condition 1"].tolist():
            if self._sample_rep == 1:
                # this sample is wt 
                # run calibratePhred 
                # we only need 1 wt calibration for each tile (r1 and r2)
                # we are always using the first rep 
                wt_id = self._sample_id
            else:
                # this sample is rep 2 of the wt
                # find corresponding rep 1 of the wt 
                wt_id = self._samples[(self._samples["Condition"] == self._sample_condition) & (self._samples["Tile ID"] == self._sample_tile) & (self._samples["Replicate"] == 1)]["Sample ID"].values[0]

        elif self._sample_condition in wt["Condition 2"].tolist():
            # get corresponding wt 
            wt_name = wt[wt["Condition 2"] == self._sample_condition]["Condition 1"].values[0]
            # get wt sample for this tile and rep 1
            wt_id = self._samples[(self._samples["Condition"] == wt_name) & (self._samples["Tile ID"] == self._sample_tile) & (self._samples["Replicate"] == 1)]["Sample ID"].values[0]
        
        # check if calibrated
        self._mut_log.info(f"wt sample used to adjust phred scores: {wt_id}")
        phred_output_r1 = os.path.join(self._output_counts_dir, f"{wt_id}_T{self._sample_tile}_R1_calibrate_phred.csv")
        phred_output_r2 = os.path.join(self._output_counts_dir, f"{wt_id}_T{self._sample_tile}_R2_calibrate_phred.csv")
        # if phred output file is not generated
        if not os.path.isfile(phred_output_r1):
            # create an empty file as place holder
            # to prevent other jobs writing to the same file
            with open(phred_output_r1, 'w') as fp:
                pass
            log_f = os.path.join(self._output_counts_dir, f"{wt_id}_R1_phred.log")
            # find wt read 1
            # print(f"{os.path.dirname(self._r1)}/{wt_id}*_R1*.sam")
            wt_r1_sam  = glob.glob(f"{os.path.dirname(self._r1)}/{wt_id}*_R1*.sam")[0]
            self._mut_log.info("WT file used for calibration:\t" + wt_r1_sam)
            cmd_r1 = f"tsm calibratePhred {wt_r1_sam} -p {self._param} -o {phred_output_r1} -l {log_f} --silent --cores {self._cores}"
            if self._sroverride:
                cmd_r1 = cmd_r1 + " --srOverride"
            self._mut_log.info(cmd_r1)
            os.system(cmd_r1)

        if not os.path.isfile(phred_output_r2):
            # create an empty file as place holder 
            with open(phred_output_r2, 'w') as fp:
                pass
            log_f = os.path.join(self._output_counts_dir, f"{wt_id}_R2_phred.log")
            wt_r2_sam  = glob.glob(f"{os.path.dirname(self._r2)}/{wt_id}*_R2*.sam")[0]
            self._mut_log.info("WT file used for calibration:\t" + wt_r2_sam)
            cmd_r2 = f"tsm calibratePhred {wt_r2_sam} -p {self._param} -o {phred_output_r2} -l {log_f} --silent --cores {self._cores}"
            if self._sroverride:
                cmd_r2 = cmd_r2 + " --srOverride"
            self._mut_log.info(cmd_r2)
            os.system(cmd_r2)
        
        # check if both file has something in there
        # if they are empty, wait for them to finish (they might be running in other jobs
        t0 = time.time()  # check for timeout
        while os.stat(phred_output_r1).st_size == 0:
            self._mut_log.warning("phred file (R1) found, but empty.. Waiting for the job to finish...")
            os.system("sleep 300")
            t1 = time.time()
            total = float(t1-t0)
            if total > 10800:
                raise TimeoutError(f"Wating for Phred file {phred_output_r1} timeout")
        t0 = time.time()
        while os.stat(phred_output_r2).st_size == 0:

            self._mut_log.warning("phred file (R2) found, but empty.. Waiting for the job to finish...")
            os.system("sleep 300")
            t1 = time.time()
            total = float(t1-t0)
            if total > 10800:
                raise TimeoutError(f"Wating for Phred file {phred_output_r2} timeout")
        time.sleep(30)
        self._mut_log.info(f"Adjusted thred files generated: {phred_output_r1}, {phred_output_r2}")
        return [phred_output_r1, phred_output_r2]

    def adjust_er_phix(self):
        """
        If not phix thred adjusted, run tsm calibratePhred with phix reads
        Before running calibrate phred, shrink sam file size by selecting reads aligned to phix only
        """
        dir_path = os.path.dirname(os.path.realpath(__file__))
        phix_fasta = os.path.join(dir_path, "data/phix.fasta")
        tmp_header = os.path.join(self._output_counts_dir, "fakeheader")
        sam_dir = os.path.dirname(self._r1)
        phix_r1 = glob.glob(f"{sam_dir}/Undetermined*_R1_*.sam")[0]
        phix_r2 = glob.glob(f"{sam_dir}/Undetermined*_R2_*.sam")[0]
        with open(tmp_header, "w") as sam_header:
            sam_header.write("@SQ	SN:phiX	LN:5386")

        phred_output_r1 = os.path.join(self._output_counts_dir, f"phix_R1_calibrate_phred.csv")
        phred_output_r2 = os.path.join(self._output_counts_dir, f"phix_R2_calibrate_phred.csv")
        if not os.path.isfile(phred_output_r1):

            # create an empty file as place holder
            with open(phred_output_r1, 'w') as fp:
                pass
            # trim the sam file
            trimmed = phix_r1.replace(".sam", "_trimmed.sam")
            cmd = f"cat {tmp_header} {phix_r1} | samtools view -S -F 4 - -o {trimmed}"
            os.system(cmd)

            log_f = os.path.join(self._output_counts_dir, f"phix_phred.log")
            cmd = f"tsm calibratePhred {phix_r1} -p {self._param} -o {phred_output_r1} -l {log_f} --silent --cores {self._cores} --fastaref {phix_fasta}"
            if self._sroverride:
                cmd = cmd + " --srOverride"
            os.system(cmd)

        if not os.path.isfile(phred_output_r2):
            # create an empty file as place holder
            with open(phred_output_r2, 'w') as fp:
                pass
            # trim the sam file
            trimmed = phix_r2.replace(".sam", "_trimmed.sam")
            cmd = f"cat {tmp_header} {phix_r2} | samtools view -S -F 4 - -o {trimmed}"
            os.system(cmd)

            log_f = os.path.join(self._output_counts_dir, f"phix_phred.log")
            cmd = f"tsm calibratePhred {phix_r2} -p {self._param} -o {phred_output_r2} -l {log_f} --silent --cores" \
                  f" {self._cores} --fastaref {phix_fasta}"
            if self._sroverride:
                cmd = cmd + " --srOverride"
            os.system(cmd)

        # check if both file has something in there
        # if they are empty, wait for them to finish (they might be running in other jobs
        t0 = time.time()  # check for timeout
        while os.stat(phred_output_r1).st_size == 0:
            os.system("sleep 300")
            self._mut_log.warning("phred file (R1) found, but empty.. Waiting for the job to finish...")
            t1 = time.time()
            total = float(t1 - t0)
            if total > 10800:
                raise TimeoutError(f"Wating for Phred file {phred_output_r1} timeout")
        t0 = time.time()
        while os.stat(phred_output_r2).st_size == 0:
            os.system("sleep 300")
            self._mut_log.warning("phred file (R2) found, but empty.. Waiting for the job to finish...")
            t1 = time.time()
            total = float(t1 - t0)
            if total > 10800:
                raise TimeoutError(f"Wating for Phred file {phred_output_r2} timeout")
        time.sleep(30)
        self._mut_log.info(f"Adjusted thred files generated: {phred_output_r1}, {phred_output_r2}")
        return [phred_output_r1, phred_output_r2]


    def multi_core(self, adjusted_er):
        """
        Read two sam files at the same time, store mutations that passed filter
        """
        read_pair = 0  # total read pairs
        un_map = 0  # total number of unmapped reads
        read_nomut = 0  # read pairs that have no mutations

        hgvs_output = {}  # save all the hgvs in this pair of sam files
        r1_pop_hgvs = {}  # save all the hgvs that are 2/3nt changes based only on R1
        r2_pop_hgvs = {}  # save all the hgvs that are 2/3nt changes based only on R2
        off_mut = {}  # save all the positions that were mapped but outside of the tile
        row = {}  # create a dictionary to save all the reads information to pass on to locate_mut.py

        final_pairs = 0
        r1_popmut = 0
        r2_popmut = 0
        off_read = 0

        # used to track read length for read 1 and read 2
        # this information will be written in the counts file
        read_length_r1 = {}
        read_length_r2 = {}

        # used to track reads mapped outside of a tile
        offmap = open(self._sample_offreads_f, "a")

        # stream sam file 1 and sam file 2
        r1_f = open(self._r1, "r")
        r2_f = open(self._r2, "r")
        # init pool
        pool = mp.Pool(self._cores)
        jobs = []
        self._mut_log.info("Start reading SAM files...")
        for line_r1, line_r2 in zip(r1_f, r2_f):
            line_r1 = line_r1.split()
            line_r2 = line_r2.split()
            # if a line has less than 11 fields, skip this line
            if len(line_r1) < 11 or len(line_r2) < 11:
                # the read has no sequence
                self._mut_log.warning(line_r1)
                self._mut_log.warning(line_r2)
                self._mut_log.warning("Missing fields in read!")
                continue

            # count how many read pairs in this pair of sam files
            read_pair += 1
            mapped_name_r1 = line_r1[2]
            mapped_name_r2 = line_r2[2]
            # if one of the read didn't map to ref
            if mapped_name_r1 == "*" or mapped_name_r2 == "*":
                un_map += 1
                continue

            # check if read ID mapped
            read_name_r1 = line_r1[0]
            read_name_r2 = line_r2[0]
            if read_name_r1 != read_name_r2:
                self._mut_log.info(read_name_r1)
                self._mut_log.info(read_name_r2)
                self._mut_log.error("Read pair IDs did not map, please check fastq files")
                raise  ValueError("Read pair IDs did not map, please check fastq files")

            # get starting position for r1 and r2
            pos_start_r1 = line_r1[3]
            pos_start_r2 = line_r2[3]

            # record reads that are not mapped to this tile
            # get CIGAR string, sequence and quality
            CIGAR_r1 = line_r1[5]
            seq_r1 = line_r1[9]
            quality_r1 = line_r1[10]

            CIGAR_r2 = line_r2[5]
            seq_r2 = line_r2[9]
            quality_r2 = line_r2[10]

            r1_end = int(pos_start_r1) + len(seq_r1)
            r2_end = int(pos_start_r2) + len(seq_r2)
            # read 1 sequence must cover from start of the tile to 90% of the tile
            if ((r1_end - int(self._cds_start)) < (int(self._tile_begins) + int(self._min_map_len))) or \
                    ((r2_end - int(self._cds_start)) < (int(self._tile_begins) + int(self._min_map_len))):
                offmap.write(f"{int(pos_start_r1) - int(self._cds_start)}, {len(seq_r1)}, {int(pos_start_r2) - int(self._cds_start)}, {len(seq_r2)}\n")
                off_read += 1
                continue

            mdz_r1 = [i for i in line_r1 if "MD:Z:" in i]
            mdz_r2 = [i for i in line_r2 if "MD:Z:" in i]

            if len(mdz_r1) != 0:
                mdz_r1 = mdz_r1[0].split(":")[-1]
            else:
                mdz_r1 = ""

            if len(mdz_r2) != 0:
                mdz_r2 = mdz_r2[0].split(":")[-1]
            else:
                mdz_r2 = ""

            if ((not re.search('[a-zA-Z]', mdz_r1)) and ("I" not in CIGAR_r1) and ("D" not in CIGAR_r1)) and \
                    ((not re.search('[a-zA-Z]', mdz_r2)) and ("I" not in CIGAR_r2) and ("D" not in CIGAR_r2)):
                # if MDZ string only contains numbers
                # and no insertions shown in CIGAR string
                # means there is no mutation in this read
                # if both reads have no mutations in them, skip this pair
                read_nomut += 1
                continue

            # make the reads in the format of a dictionary
            # columns=["pos_start", "qual", "CIGAR", "mdz","seq"])
            row["pos_start_r1"] = pos_start_r1
            row["qual_r1"] = quality_r1
            row["cigar_r1"] = CIGAR_r1
            row["mdz_r1"] = mdz_r1
            row["seq_r1"] = seq_r1

            row["pos_start_r2"] = pos_start_r2
            row["qual_r2"] = quality_r2
            row["cigar_r2"] = CIGAR_r2
            row["mdz_r2"] = mdz_r2
            row["seq_r2"] = seq_r2

            # for all the reads that reached to this stage, check the read len for R1 and R2
            # read len = length of the actual reads
            read_len_r1 = len(seq_r1)
            read_len_r2 = len(seq_r2)
            # this info is for debug purpose
            if read_length_r1.get(read_len_r1, -1) == -1:
                read_length_r1[read_len_r1] = 1
            else:
                read_length_r1[read_len_r1] += 1

            if read_length_r2.get(read_len_r2, -1) == -1:
                read_length_r2[read_len_r2] = 1
            else:
                read_length_r2[read_len_r2] += 1

            jobs.append(pool.apply_async(process_wrapper, (row, self._seq, self._cds_seq, self._seq_lookup, self._tile_begins,
                                               self._tile_ends, self._qual, self._mut_log, self._mutrate,
                                                           self._base, self._posteriorQC, self._error_override, adjusted_er)))
            row = {} # flush
        offmap.close()
        r1_f.close()
        r2_f.close()
        self._mut_log.info("File streamed to subprocesses, waiting for jobs to finish...")
        # wait for all jobs to finish
        all_df = []
        for job in jobs:
            hgvs, outside_mut, all_posterior, hgvs_r1_clusters, hgvs_r2_clusters, track_df = job.get()
            if len(hgvs) != 0:
                final_pairs += 1
                if hgvs_output.get(hgvs, -1) != -1:
                    hgvs_output[hgvs] += 1
                else:
                    hgvs_output[hgvs] = 1

            if len(hgvs_r1_clusters) != 0:
                r1_popmut += 1
                if hgvs_r1_clusters in r1_pop_hgvs:
                    r1_pop_hgvs[hgvs_r1_clusters] += 1
                else:
                    r1_pop_hgvs[hgvs_r1_clusters] = 1

            if len(hgvs_r2_clusters) != 0:
                r2_popmut += 1
                if hgvs_r2_clusters in r2_pop_hgvs:
                    r2_pop_hgvs[hgvs_r2_clusters] += 1
                else:
                    r2_pop_hgvs[hgvs_r2_clusters] = 1

            if outside_mut != []:
                outside_mut = list(set(outside_mut))
                for i in outside_mut:
                    if not (i in off_mut):
                        off_mut[i] = 1
            if not all_posterior.empty:
                all_df.append(all_posterior)
            track_df = track_df.set_index("pos")
            track_df = track_df[track_df.index.isin(self._track_reads.index)]
            # add track df to track summary
            track_all = pd.concat([self._track_reads, track_df], axis=1).fillna(0)
            self._track_reads = track_all.groupby(by=track_all.columns, axis=1).sum()
        # clean up
        pool.close()
        self._mut_log.info("Pool closed")
        # track mutations that are not within the same tile
        # track sequencing depth for each sample
        # write this information to a tmp file
        # merge the tmp file (in main.py)

        # tmp_f = os.path.join(self._output_counts_dir, f"{self._sample_id}_tmp.csv")
        # off_mut = list(set(off_mut))
        # off_mut.sort()
        # f = open(tmp_f, "w")
        # f.write(
        #     f"sample,{self._sample_id},tile,{self._tile_begins}-{self._tile_ends},sequencing_depth,{read_pair},off_tile_reads,{off_read},off_tile_perc,{off_read / read_pair}\n")
        # f.close()

        output_csv = open(self._sample_counts_f, "a")
        output_csv.write(f"#Raw read depth:{read_pair}\n")
        output_csv.write(
            f"#Number of read pairs without mutations:{read_nomut}\n#Number of read pairs did not map to gene:{un_map}\n#Mutation positions outside of the tile:{off_mut}\n#Number of reads outside of the tile:{off_read}\n")
        output_csv.write(f"#Final read-depth:{read_pair - un_map - off_read}\n")
        output_csv.write(f"#Total read pairs with mutations:{final_pairs}\n")
        output_csv.write(
            f"#Comment: Total read pairs with mutations = Read pairs with mutations that passed the posterior threshold\n#Comment: Final read-depth = raw read depth - reads didn't map to gene - reads mapped outside of the tile\n")
        output_csv.write(f"#Read len for R1: {read_length_r1}\n")
        output_csv.write(f"#Read len for R2: {read_length_r2}\n")
        output_csv.close()

        self._mut_log.info(f"Raw sequencing depth: {read_pair}")
        self._mut_log.info(f"Number of reads without mutations:{read_nomut}")
        self._mut_log.info(f"Final read-depth: {read_pair - un_map - off_read}")
        # convert list to df with one col
        hgvs_df = pd.DataFrame.from_dict(hgvs_output, orient="index")
        hgvs_df = hgvs_df.reset_index()
        if not hgvs_df.empty:
            hgvs_df.columns = ["HGVS", "count"]
        else:
            hgvs_df = pd.DataFrame({"HGVS": [], "count": []})
        hgvs_df.to_csv(self._sample_counts_f, mode="a", index=False)
        del hgvs_df

        # save read coverage to csv file
        cov_file = os.path.join(self._output_counts_dir, f"coverage_{self._sample_id}.csv")
        self._track_reads.to_csv(cov_file)

        if self._posteriorQC:
            r1_df = pd.DataFrame.from_dict(r1_pop_hgvs, orient="index")
            r1_df = r1_df.reset_index()
            if not r1_df.empty:
                r1_df.columns = ["HGVS", "count"]
            r1_df.to_csv(self._sample_counts_r1_f, mode="a", index=False)
            del r1_f

            r2_df = pd.DataFrame.from_dict(r2_pop_hgvs, orient="index")
            r2_df = r2_df.reset_index()
            if not r2_df.empty:
                r2_df.columns = ["HGVS", "count"]
            r2_df.to_csv(self._sample_counts_r2_f, mode="a", index=False)
            del r2_f

            self._mut_log.info(f"Reading posterior unfiltered dfs ...")
            all_f = os.path.join(self._output_counts_dir, f"{self._sample_id}_posprob_all.csv")
            for df in all_df:
                df.to_csv(all_f, mode='a', header=False, index=False)
            self._mut_log.info(f"Posterior df saved to files ... {all_f}")
            del all_df


def process_wrapper(row, seq, cds_seq, seq_lookup, tile_begins, tile_ends, qual, locate_log, mutrate, base,
                    posteriorQC, errorOverride, adjusted_er):
    """
    Wrapper function to process each line (pair of reads)
    """
    mut_parser = locate_mut.MutParser(row, seq, cds_seq, seq_lookup, tile_begins, tile_ends, qual, locate_log,
                                      mutrate, base, posteriorQC, errorOverride, adjusted_er)
    hgvs, outside_mut, all_df, hgvs_r1_clusters, hgvs_r2_clusters, track_df = mut_parser._main()
    return hgvs, outside_mut, all_df, hgvs_r1_clusters, hgvs_r2_clusters, track_df
