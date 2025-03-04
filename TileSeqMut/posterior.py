#!/usr/bin/env python3.7

# This is used to calculate the posterior probability for variants
# the script is modified from varCallerSnippet.r

import math
import pandas as pd
import numpy as np
import time
from fractions import Fraction

def cluster(mut_cluster: dict, r1_qual: str, r2_qual: str, r1_mappos: dict,
                        r2_mappos: dict, mut_rate: float, cut_off: float, base: int, posteriorQC: bool, error_override: bool, adjustthred:list):
    """
    Parse cluster of mutations
    Mutations are clustered together if they are within
    @param mut_cluster dictionary where the values are dfs with mutations from both r1 and r2, the keys are
    group numbers, values are df of mutations
    @param r1_qual: quality of read 1
    @param r2_qual: quality of read 2
    @param r1_mappos: dictionary with position mapped between ref and read 1
    @param r2_mappos: dictionary with position mapped between ref and read 2
    @param mut_rate: mut rate defined by user
    @param cut_off: posterior cutoff defined by user
    @param base: phred base (default 33)
    @param posteriorQC: boolean, if posterior QC required
    @param error_override: boolean, if error probabilities should be obtained from observed error
    @param adjustthred: list with path to phred calibration files for R1 and R2
    """
    if posteriorQC:
        pos_df, wt_df, all_df, clustered_r1_mut, clustered_r2_mut = call_prob_withposterior(mut_cluster, r1_mappos,
                                                                                     r2_mappos, r1_qual, r2_qual,
                                                                                     mut_rate, base, cut_off, error_override, adjustthred)
    else:
        pos_df, wt_df, all_df, clustered_r1_mut, clustered_r2_mut = call_prob(mut_cluster, r1_mappos, r2_mappos, r1_qual,
                                                                       r2_qual, mut_rate, base, cut_off, error_override, adjustthred)

    return pos_df, wt_df, all_df, clustered_r1_mut, clustered_r2_mut


def call_prob(mut_cluster, r1_mappos, r2_mappos, r1_qual, r2_qual, mut_rate, base, cut_off, error_override, adjustthred):
    """
    @param mut_cluster dictionary where the values are dfs with mutations from both r1 and r2, the keys are
    group numbers, values are df of mutations
    @param r1_mappos: dictionary with position mapped between ref and read 1
    @param r2_mappos: dictionary with position mapped between ref and read 2
    @param r1_qual: quality of read 1
    @param r2_qual: quality of read 2
    @param mut_rate: mut rate defined by user
    @param base: phred base (default 33)
    @param cut_off: posterior cutoff defined by user
    @param error_override: boolean, if error probabilities should be obtained from observed error
    @param adjustthred: list with path to phred calibration files for R1 and R2
    """

    pos_prob = {"m": [], "prob": [], "read": []}
    wt_calls = {"m": [], "read": []} # to store confident WT calls

    for c in mut_cluster.keys():
        mutcall = mut_cluster[c]
        # each mutcall is a df with columns:
        # m_r1,read,pos,ref_r1,alt_r1,qual_r1,m_r2,ref_r2,alt_r2,qual_r2
        c_size = mutcall.shape[0]  # size of the cluster
        # iterate through the clusters
        for index, row in mutcall.iterrows():
            # if the mut is on read 2 only
            if (pd.isnull(row["m_r1"]) or row["alt_r1"] == "N") and not pd.isnull(row["ref_r2"]):
                # if the position is not on read 1
                if r1_mappos.get(int(row["pos"])) is None:
                    continue

                # get R1 wt qual with map pos and qual str
                r1_qual_base = r1_qual[r1_mappos[int(row["pos"])]]
                # calculate posterior prob
                pos = bayesian_variant_call([row["ref_r2"], row["alt_r2"]], [r1_qual_base, row["qual_r2"]],
                                            row["ref_r2"], mut_rate, base, c_size, error_override, adjustthred=adjustthred)
                # if posterior for alt base is greater than cutoff
                # AND greater than wt
                if pos[row["alt_r2"]] > cut_off and pos[row["ref_r2"]] <= pos[row["alt_r2"]]:
                    pos_prob["m"].append(row["m_r2"])
                    pos_prob["prob"].append(pos[row["alt_r2"]])
                    pos_prob["read"].append("r2")

                # if posterior for WT base is greater than cutoff AND greater than alt base
                elif pos[row["ref_r2"]] > cut_off and pos[row["alt_r2"]] < pos[row["ref_r2"]]:
                    wt_calls["m"].append(row["m_r2"])
                    wt_calls["read"].append("r2")

            # if r1 is not null
            elif (pd.isnull(row["m_r2"]) or row["alt_r2"] == "N") and not pd.isnull(row["ref_r1"]):

                if r2_mappos.get(int(row["pos"])) is None:
                    continue
                # get R1 wt qual with map pos and qual str
                r2_qual_base = r2_qual[r2_mappos[int(row["pos"])]]

                pos = bayesian_variant_call([row["alt_r1"], row["ref_r1"]], [row["qual_r1"], r2_qual_base],
                                            row["ref_r1"],
                                            mut_rate,
                                            base, c_size, error_override, adjustthred=adjustthred)
                if pos[row["alt_r1"]] > cut_off and pos[row["ref_r1"]] <= pos[row["alt_r1"]]:

                    pos_prob["m"].append(row["m_r1"])
                    pos_prob["prob"].append(pos[row["alt_r1"]])
                    pos_prob["read"].append("r1")

                # if posterior for WT base is greater than cutoff AND greater than alt base 
                elif pos[row["ref_r1"]] > cut_off and pos[row["alt_r1"]] < pos[row["ref_r1"]]:
                    wt_calls["m"].append(row["m_r1"])
                    wt_calls["read"].append("r1")

            elif (not pd.isnull(row["m_r2"])) and (not pd.isnull(row["ref_r1"])):

                basecall = [row["alt_r1"], row["alt_r2"]]
                qual = [row["qual_r1"], row["qual_r2"]]
                pos = bayesian_variant_call(basecall, qual, row["ref_r1"], mut_rate, base, c_size, error_override, adjustthred=adjustthred)

                # if (pos[row["ref_r1"]] > pos[row["alt_r1"]]) and (pos[row["ref_r2"]] > pos[row["alt_r2"]]): continue

                if pos[row["alt_r1"]] > pos[row["alt_r2"]] and pos[row["alt_r1"]] > cut_off:
                    pos_prob["m"].append(row["m_r1"])
                    pos_prob["prob"].append(pos[row["alt_r1"]])
                    pos_prob["read"].append("r1")

                elif pos[row["alt_r2"]] > pos[row["alt_r1"]] and pos[row["alt_r2"]] > cut_off:
                    pos_prob["m"].append(row["m_r2"])
                    pos_prob["prob"].append(pos[row["alt_r2"]])
                    pos_prob["read"].append("r2")

                elif pos[row["alt_r1"]] == pos[row["alt_r2"]] and pos[row["alt_r1"]] > cut_off:
                    pos_prob["m"].append(row["m_r1"])
                    pos_prob["prob"].append((pos[row["alt_r1"]], pos[row["alt_r2"]]))
                    pos_prob["read"].append(("r1", "r2"))

                # if posterior for WT base is greater than cutoff AND greater than alt base
                elif pos[row["ref_r1"]] > max(pos[row["alt_r1"]], pos[row["alt_r2"]]) and pos[row["ref_r1"]] > cut_off:
                    wt_calls["m"].append(row["m_r1"])
                    wt_calls["read"].append("r1")

    pos_df = pd.DataFrame(pos_prob)
    wt_df = pd.DataFrame(wt_calls)
    return pos_df, wt_df, pd.DataFrame({}), pd.DataFrame({}), pd.DataFrame({})


def call_prob_withposterior(mut_cluster, r1_mappos, r2_mappos, r1_qual, r2_qual, mut_rate, base, cut_off, error_override, adjustthred):
    """
    @param mut_cluster dictionary where the values are dfs with mutations from both r1 and r2, the keys are
    group numbers, values are df of mutations
    @param r1_mappos: dictionary with position mapped between ref and read 1
    @param r2_mappos: dictionary with position mapped between ref and read 2
    @param r1_qual: quality of read 1
    @param r2_qual: quality of read 2
    @param mut_rate: mut rate defined by user
    @param base: phred base (default 33)
    @param cut_off: posterior cutoff defined by user
    @param error_override: boolean, if error probabilities should be obtained from observed error
    @param adjustthred: list with path to phred calibration files for R1 and R2
    """

    pos_prob = {"m": [], "prob": [], "read": []}
    all_prob = {"m": [], "prob": [], "read": [], "pass": []}
    wt_prob = {"m": [], "read": []} # to store confident WT calls

    # this is used to record 2/3nt changes on r1 or r2
    clustered_r1_mut = [pd.DataFrame({"m": [], "prob": [], "read": []})]
    clustered_r2_mut = [pd.DataFrame({"m": [], "prob": [], "read": []})]
    for c in mut_cluster.keys():
        mutcall = mut_cluster[c]
        # each mutcall is a df with columns:
        # m_r1,read,pos,ref_r1,alt_r1,qual_r1,m_r2,ref_r2,alt_r2,qual_r2
        c_size = mutcall.shape[0]  # size of the cluster
        # iterate through the clusters
        pos = {}
        r = ""

        # use this df to save the temporary mutations in this cluster
        tmp_cluster_mut = {"m": [], "prob": [], "read":[]}
        for index, row in mutcall.iterrows():

            if (pd.isnull(row["m_r1"]) or row["alt_r1"] == "N") and not pd.isnull(row["ref_r2"]):
               
                if r1_mappos.get(int(row["pos"])) is None:
                    continue

                # get R1 wt qual with map pos and qual str
                r1_qual_base = r1_qual[r1_mappos[int(row["pos"])]]
                # calculate posterior prob
                pos = bayesian_variant_call([row["ref_r2"], row["alt_r2"]], [r1_qual_base, row["qual_r2"]],
                                            row["ref_r2"], mut_rate, base, c_size)

                # add this information to all prob df
                all_prob["m"].append(row["m_r2"])
                all_prob["prob"].append(pos[row["alt_r2"]])
                all_prob["read"].append("r2")

                # if pos[row["ref_r2"]] > pos[row["alt_r2"]]:
                #     # this means that wt has higher pos-prob, then we dump this mutation
                #     all_prob["pass"].append(-1)

                # if we only consider the probability of mutation that we got from r2
                # check if this is greater than the probability cutoff
                if pos[row["alt_r2"]] > cut_off and c_size > 1:
                    # add this information to tmp prob df
                    tmp_cluster_mut["m"].append(row["m_r2"])
                    tmp_cluster_mut["prob"].append(pos[row["alt_r2"]])
                    tmp_cluster_mut["read"].append("r2")

                if pos[row["alt_r2"]] > cut_off and pos[row["ref_r2"]] <= pos[row["alt_r2"]]:
                    all_prob["pass"].append(1)

                    pos_prob["m"].append(row["m_r2"])
                    pos_prob["prob"].append(pos[row["alt_r2"]])
                    pos_prob["read"].append("r2")
                else:
                    all_prob["pass"].append(-1)
                    # if posterior for WT base is greater than cutoff AND greater than alt base
                    if pos[row["ref_r2"]] > cut_off and pos[row["alt_r2"]] <= pos[row["ref_r2"]]:
                        wt_prob["m"].append(row["m_r2"])
                        wt_prob["read"].append("r2")

                r = "r2"

            elif (pd.isnull(row["m_r2"]) or row["alt_r2"] == "N") and not pd.isnull(row["ref_r1"]):

                if r2_mappos.get(int(row["pos"])) is None:
                    continue

                # get R2 wt qual with map pos and qual str
                r2_qual_base = r2_qual[r2_mappos[int(row["pos"])]]

                pos = bayesian_variant_call([row["alt_r1"], row["ref_r1"]], [row["qual_r1"], r2_qual_base],
                                            row["ref_r1"], mut_rate, base, c_size, error_override, adjustthred)

                # add this information to all prob df
                all_prob["m"].append(row["m_r1"])
                all_prob["prob"].append(pos[row["alt_r1"]])
                all_prob["read"].append("r1")

                # if we only consider the probability of mutation that we got from r2
                # check if this is greater than the probability cutoff
                if pos[row["alt_r1"]] > cut_off and c_size > 1:
                    # add this information to tmp prob df
                    tmp_cluster_mut["m"].append(row["m_r1"])
                    tmp_cluster_mut["prob"].append(pos[row["alt_r1"]])
                    tmp_cluster_mut["read"].append("r1")

                if pos[row["alt_r1"]] > cut_off and pos[row["ref_r1"]] <= pos[row["alt_r1"]]:
                    all_prob["pass"].append(1)

                    pos_prob["m"].append(row["m_r1"])
                    pos_prob["prob"].append(pos[row["alt_r1"]])
                    pos_prob["read"].append("r1")

                else:
                    all_prob["pass"].append(-1)
                    # if posterior for WT base is greater than cutoff AND greater than alt base
                    if pos[row["ref_r1"]] > cut_off and pos[row["alt_r1"]] <= pos[row["ref_r1"]]:
                        wt_prob["m"].append(row["m_r1"])
                        wt_prob["read"].append("r1")
                r = "r1"

            elif (not pd.isnull(row["m_r2"])) and (not pd.isnull(row["ref_r1"])):

                basecall = [row["alt_r1"], row["alt_r2"]]
                qual = [row["qual_r1"], row["qual_r2"]]
                pos = bayesian_variant_call(basecall, qual, row["ref_r1"], mut_rate, base, c_size, error_override, adjustthred)
                # add this information to all prob df
                all_prob["m"].append(row["m_r1"])
                all_prob["prob"].append((pos[row["alt_r1"]], pos[row["alt_r2"]]))
                all_prob["read"].append(("r1", "r2"))
                all_prob["pass"].append("-")

                # this step is used for QC
                # the idea is that to count number of mutations on each read that passed the cutoff (separately)
                # if we only consider the probability of mutation that we got from r1
                # check if this is greater than the probability cutoff
                if pos[row["alt_r1"]] > cut_off and c_size > 1:
                    # add this information to tmp prob df
                    tmp_cluster_mut["m"].append(row["m_r1"])
                    tmp_cluster_mut["prob"].append(pos[row["alt_r1"]])
                    tmp_cluster_mut["read"].append("r1")
                # if we only consider the probability of mutation that we got from r2
                # check if this is greater than the probability cutoff
                if pos[row["alt_r2"]] > cut_off and c_size > 1:
                    # add this information to tmp prob df
                    tmp_cluster_mut["m"].append(row["m_r2"])
                    tmp_cluster_mut["prob"].append(pos[row["alt_r2"]])
                    tmp_cluster_mut["read"].append("r2")

            # this means that there are mut on both reads
            if r == "":
                if pos[row["alt_r1"]] > pos[row["alt_r2"]] and pos[row["alt_r1"]] > cut_off:
                    pos_prob["m"].append(row["m_r1"])
                    pos_prob["prob"].append(pos[row["alt_r1"]])
                    pos_prob["read"].append("r1")

                elif pos[row["alt_r2"]] > pos[row["alt_r1"]] and pos[row["alt_r2"]] > cut_off:
                    pos_prob["m"].append(row["m_r2"])
                    pos_prob["prob"].append(pos[row["alt_r2"]])
                    pos_prob["read"].append("r2")

                elif pos[row["alt_r1"]] == pos[row["alt_r2"]] and pos[row["alt_r1"]] > cut_off:
                    pos_prob["m"].append(row["m_r1"])
                    pos_prob["prob"].append((pos[row["alt_r1"]], pos[row["alt_r2"]]))
                    pos_prob["read"].append(("r1", "r2"))

                # if posterior for WT base is greater than cutoff AND greater than alt base
                elif pos[row["ref_r1"]] > max(pos[row["alt_r1"]], pos[row["alt_r2"]]) and pos[row["ref_r1"]] > cut_off:
                    wt_prob["m"].append(row["m_r1"])
                    wt_prob["read"].append(("r1"))

            # check tmp prob df
            # if number of mutations passed on r1 or r2 greater than 1
            # output mutations, label them r1 and r2, also assign them to the same cluster
            tmp_cluster = pd.DataFrame(tmp_cluster_mut)
            r1_mut_cluster = tmp_cluster[tmp_cluster["read"] == "r1"]
            r2_mut_cluster = tmp_cluster[tmp_cluster["read"] == "r2"]
            if r1_mut_cluster.shape[0] > 1:
                clustered_r1_mut.append(r1_mut_cluster)
            if r2_mut_cluster.shape[0] > 1:
                clustered_r2_mut.append(r2_mut_cluster)

    # print(pos_prob)
    pos_df = pd.DataFrame(pos_prob)
    wt_df = pd.DataFrame(wt_prob)
    all_df = pd.DataFrame(all_prob)
    clustered_r1_mut = pd.concat(clustered_r1_mut)
    clustered_r2_mut = pd.concat(clustered_r2_mut)
    return pos_df, wt_df, all_df, clustered_r1_mut, clustered_r2_mut


def bayesian_variant_call(basecall, qual, wt, mut_rate, base, clusterSize, error_override, adjustthred):
    """
    Return a dict where the keys are bases and the values are posterior probabilities for that base.
    @param basecall: list of base calls (i.e R1 -> A R2 -> C :  ["A", "C"])
    @param phred: phred score for the base calls (in letters) ["!", "J"]
    @param wt: wild type base
    @param mut_rate: mutation rate
    @param error_override: boolean, if error probabilities should be obtained from observed error
    @param adjustthred: list with path to phred calibration files for R1 and R2
    return: dictinary with basecall as keys and post prob as values
    """
    # all possible hypo bases
    nt = list(set([wt]+basecall))
    # nt = ["A", "G", "C", "T"]
    # convert phred to int scores
    # convert phred to int scores
    # in the case of deletions or insertions, there are multiple phred scores for each mut call
    phred = []
    for i in qual:
        if len(i) == 1:
            phred.append(10 ** (-(ord(i) - int(base)) / 10))
        else:
            all_phred = [10 ** (-(ord(j) - int(base)) / 10) for j in i.split(",")]
            phred.append(np.prod(all_phred))
    
    if len(adjustthred) == 2:
        phred_r1_df = pd.read_csv(adjustthred[0], index_col=0)
        phred_r2_df = pd.read_csv(adjustthred[1], index_col=0)
        phred_r1_df["observed"] = phred_r1_df["observed"].fillna(phred_r1_df["specification"])
        phred_r2_df["observed"] = phred_r2_df["observed"].fillna(phred_r2_df["specification"])
        
        if error_override: # use the observed error rate
            phred_r1_df.observed[(phred_r1_df["observed"] == 0) | (phred_r1_df["observed"] == 1)] = phred_r1_df["specification"]
            phred_r2_df.observed[(phred_r2_df["observed"] == 0) | (phred_r2_df["observed"] == 1)] = phred_r2_df["specification"]
        else: # select the maximum error rate from observed error rate and specification error rate
            phred_r1_df.observed[(phred_r1_df["observed"] == 0) | (phred_r1_df["observed"] == 1) | (phred_r1_df["observed"] < phred_r1_df["specification"])] = phred_r1_df["specification"]
            phred_r2_df.observed[(phred_r2_df["observed"] == 0) | (phred_r2_df["observed"] == 1) | (phred_r2_df["observed"] < phred_r2_df["specification"])] = phred_r2_df["specification"]
        r1_qual = qual[0].split(",")
        phred_r1 = np.prod(phred_r1_df[phred_r1_df.index.isin(r1_qual)]["observed"])
        r2_qual = qual[1].split(",")
        phred_r2 = np.prod(phred_r2_df[phred_r2_df.index.isin(r2_qual)]["observed"])
        phred = [phred_r1, phred_r2]
    #print(f"adjusted phred: {phred}")
    # phred = [10**(-(ord(i) - 33) / 10) for i in phred]
    post_p = []
    for base in nt: # go through each nt
        log_odd = 0
        if base == wt:
            log_odd += math.log((1-mut_rate) ** clusterSize) - math.log(1-((1-mut_rate)**clusterSize)/3)
        elif base == "ins":
            log_odd += (math.log(mut_rate) - math.log(4)) - math.log(1-(mut_rate/4))
        else:
            log_odd += math.log(1-(1-mut_rate) ** clusterSize) - math.log(3) - math.log(1-(1-(1-mut_rate) ** clusterSize)/3)

        for j in range(len(basecall)):
            if basecall[j] == base:
                try:
                    log_odd += (math.log(1-phred[j]) - math.log(phred[j]) + math.log(3))
                except:
                    print(phred)
                    print(qual)
                    exit()
            else:
                log_odd += (math.log(phred[j]) - math.log(3) - math.log((1/3) -(phred[j]/9)))
        # print(log_odd)
        # print(basecall)
        if log_odd > 700: 
            logit_value = 1
        elif log_odd < -700:
            logit_value = 0
        else: 
            logit_value = math.exp(log_odd) / (1+math.exp(log_odd))
        post_p.append(logit_value)

    prob = dict(zip(nt, post_p))
    #output = dict(zip(basecall, [prob.get(base) for base in basecall]))
    #print(f"output: {output}")
    return prob


# if __name__ == "__main__":
#     basecall = ["T", "A"]
#     phred = ["I", "J"]
#     wt = "C"
#     mut_rate= 0.0025
#     prob = bayesian_variant_call(basecall, phred, wt, mut_rate)
#     print(prob)
