TITLE = "Kmer Cluster Finder"
DESC = "Given a fasta of nucleotide sequence(s), window size, and max mismatch threshold, find all sequence groups that meet the score threshold."
'''
Date: July 10, 2023
@author: sbrown
'''

## Import Libraries

import sys
import argparse
import os
import time
import math


## Declare global variables

DEBUG = False
VERB = False


## Classes and functions

class bcolors:
    CYAN = '\033[1;36;40m'
    BLUE = '\033[1;34;40m'
    GREEN = '\033[1;32;40m'
    YELLOW = '\033[1;33;40m'
    RED = '\033[1;31;40m'
    BOLDWHITE = '\033[1;37;40m'
    DARK = '\033[1;30;40m'
    PURPLE = '\033[1;35;40m'
    ENDC = '\033[0m'

def statprint(msg, msg_type = "STATUS"):
    typeColour = ""
    if msg_type == "ERROR":
        typeColour = bcolors.RED
    elif msg_type == "WARNING":
        typeColour = bcolors.YELLOW
    elif msg_type == "DEBUG":
        typeColour = bcolors.GREEN
    elif msg_type == "SUBPROCESS":
        typeColour = bcolors.GREEN
        msg_type = "     " + msg_type
    else:
        typeColour = bcolors.BOLDWHITE

    print("{message_color}{message_type}{end_color} {time_color}[{datetime}]{end_color}: {message}".format(message_color = typeColour, 
             message_type = msg_type,
             end_color = bcolors.ENDC, 
             time_color = bcolors.BLUE, 
             datetime = time.strftime("%Y/%m/%d %T"), 
             message = msg), flush=True)

COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
def revcomp(seq):
    return ''.join([COMPLEMENT[nt] for nt in seq[::-1]])

## Main

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = TITLE)
    parser.add_argument("--sequence", dest = "SEQFILE", help = "Input Sequence File", type = str)
    parser.add_argument("--window", dest = "WINDOW", help = "Window size", type = int)
    parser.add_argument("--mismatch", dest = "MISMATCH", help = "Number of mismatches allowed.", type = int)
    parser.add_argument("--output", dest = "OUTPUT", help = "Output tsv file to write", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    print("=======================================================")
    print("Python version: {}".format(sys.version))
    print("Python environment: {}".format(sys.prefix))
    print("Server: {}".format(os.uname()[1]))
    print("Current directory: {}".format(os.getcwd()))
    print("Command: {}".format(" ".join(sys.argv)))
    print("Time: {}".format(time.strftime("%Y/%m/%d %T")))
    print("=======================================================\n")


    ## Script content here

    statprint("Reading in sequence...")
    seq = ""
    seqs = []
    total_seq_len = 0
    for line in open(args.SEQFILE, "r"):
        if line.startswith(">"):
            ## write previous seq if it exists.
            if seq != "":
                seqs.append(seq)
                total_seq_len += len(seq)
                seq = ""
        else:
            seq += line.rstrip()
    ## write last sequence
    seqs.append(seq)
    total_seq_len += len(seq)

    total_estimated_comparisons = math.floor(((total_seq_len - args.WINDOW + 1) * (total_seq_len - args.WINDOW + 1)) / 2)
    statprint("There are {:,} estimated kmer comparisons to complete.".format(total_estimated_comparisons))

    progress_print_chunk_size = math.floor(total_estimated_comparisons / 1000)

    groups = []
    completed_comparisons = 0
    statprint("Calculating scores and grouping...")
    #for seq_a in seqs:
    for seq_a_i in range(len(seqs)):
        #for seq_b in seqs:
        for seq_b_i in range(seq_a_i, len(seqs)):   ## only do "upper triangle" of comparisons. ie if comparing A to B, don't also compare B to A.
            seq_a = seqs[seq_a_i]
            seq_b = seqs[seq_b_i]
            for r in range(len(seq_a) - args.WINDOW + 1):
                score_row = []
                if DEBUG and r%100==0: statprint("Checking all matches for position {}".format(r))
                for c in range(len(seq_b) - args.WINDOW + 1):
                    completed_comparisons += 1
                    if completed_comparisons % progress_print_chunk_size == 0:
                        percent_complete = 100 * completed_comparisons / total_estimated_comparisons
                        statprint("{} % complete ({} / {} comparions checked).".format(percent_complete, completed_comparisons, total_estimated_comparisons))

                    score_fwd = sum(seq_a[r:r+args.WINDOW][k] == seq_b[c:c+args.WINDOW][k] for k in range(args.WINDOW))
                    score_rev = sum(seq_a[r:r+args.WINDOW][k] == seq_b[c:c+args.WINDOW][::-1][k] for k in range(args.WINDOW))
                    score_comp = -1*sum(seq_a[r:r+args.WINDOW][k] == COMPLEMENT[seq_b[c:c+args.WINDOW][k]] for k in range(args.WINDOW))
                    score_revcomp = -1*sum(seq_a[r:r+args.WINDOW][k] == revcomp(seq_b[c:c+args.WINDOW])[k] for k in range(args.WINDOW))
                    ## keep the max magnitude value, retaining the original sign.
                    score = max([score_fwd, score_rev, score_comp, score_revcomp], key=lambda x: abs(x))

                    if abs(score) >= (args.WINDOW - args.MISMATCH):
                        ## is sufficiently similar.
                        ## track the max scoring pair of seqs.
                        #Also, to ensure only adding each sequence once per position it occurs at, appending "_[r/c]" for the index of the position and "_i" for OCR id (in case identical kmer happens to appear in different OCRs at the same position)
                        seq1 = "{}_{}_{}".format(seq_a[r:r+args.WINDOW], r, seq_a_i)
                        if score_fwd == score:
                            seq2 = "{}_{}_{}".format(seq_b[c:c+args.WINDOW], c, seq_b_i)
                        elif score_rev == score:
                            seq2 = "{}_{}_{}".format(seq_b[c:c+args.WINDOW][::-1], c, seq_b_i)
                        elif score_comp == score:
                            seq2 = "{}_{}_{}".format("".join([COMPLEMENT[k] for k in seq_b[c:c+args.WINDOW]]), c, seq_b_i)
                        elif score_revcomp == score:
                            seq2 = "{}_{}_{}".format(revcomp(seq_b[c:c+args.WINDOW]), c, seq_b_i)
                        else:
                            ## should not be here.
                            sys.exit("Something has gone terribly wrong...")

                        
                        ## now need to track these
                        ## ONLY IF NOT ON THE DIAGONAL (ie seq1==seq2, including the position ("seq1" includes the position and OCRind))
                        if not seq1 == seq2:
                            indexes_to_merge = []
                            seqs_added = False
                            for i in range(len(groups)):
                                if seq1 in groups[i]:
                                    ## add seq2 to this group
                                    groups[i].append(seq2)
                                    seqs_added = True
                                    ## potentially need to merge to this
                                    indexes_to_merge.append(i)
                                elif seq2 in groups[i]:
                                    ## add seq1 to this group
                                    groups[i].append(seq1)
                                    seqs_added = True
                                    ## potentially need to merge to this
                                    indexes_to_merge.append(i)
                            ## check if they were not grouped
                            if not seqs_added:
                                groups.append([seq1, seq2])
                            ## check if need to merge
                            if len(indexes_to_merge) > 1:
                                if len(indexes_to_merge) > 2:
                                    sys.exit("Oh dear, there should only ever be two sets that need merging...")
                                else:
                                    # groups[indexes_to_merge[0]].update(groups[indexes_to_merge[1]])
                                    groups[indexes_to_merge[0]] += groups[indexes_to_merge[1]]
                                    del groups[indexes_to_merge[1]]

        
    ## now have a list of sets, each set being a group of sequences.
    if len(groups) > 0:
        statprint("There are {} groups of sequences, with sizes {} (min), {} (max), and {} (mean)".format(len(groups),
                min([len(x) for x in groups]),
                max([len(x) for x in groups]),
                sum([len(x) for x in groups])/len(groups)))
    else:
        statprint("No groups")
    if DEBUG: statprint(groups, "DEBUG")



    ## write out groups of >=3 sequences.
    statprint("Writing output...")
    out = open(args.OUTPUT, "w")
    out.write("group\tkmer\n")
    i = 0
    for g in groups:
        if len(g) > 2:
            i += 1
            for k in g:
                out.write("{}\t{}\n".format(i, k.split("_")[0]))
    out.close()
    

    statprint("Done.")