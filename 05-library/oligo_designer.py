TITLE = "Oligo Designer"
DESC = "Given a list of candidates, generate all oligos to be synthesized"
'''
Date: October 30,2023
@author: sbrown
'''

## Import Libraries

import sys
import argparse
import os
import time
import Bio.SeqIO
from Bio.Seq import Seq
import math
import itertools


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


## Main

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = TITLE)
    parser.add_argument("--neutral_seq", dest = "NEUTRAL", help = "Neutral sequence, fasta.", type = str)
    parser.add_argument("--pad_re_5", dest = "PAD_RE_5", help = "5' padding and restriction enzyme site sequence", type = str)
    parser.add_argument("--re_pad_3", dest = "RE_PAD_3", help = "3' restriction enzyme site and padding sequence", type = str)
    parser.add_argument("--adaptor_5", dest = "ADAPTOR5", help = "5' Illumina adaptor sequence", type = str)
    parser.add_argument("--adaptor_3", dest = "ADAPTOR3", help = "3' Illumina adaptor sequence", type = str)
    parser.add_argument("--candidates", dest = "CANDIDATES", help = "Fasta of candidate sequences", type = str)
    parser.add_argument("--num_candidates", dest = "NUM_CANDIDATES", help = "Number of candidates per oligo", type = int)
    parser.add_argument("--candidate_gap", dest = "CANDIDATE_GAP", help = "Number of neutral bases between candidates", type = int)
    parser.add_argument("--forbidden_seqs", dest = "FORBIDDEN_SEQS", help = "Forbidden sequence(s) to not allow in oligos", type = str, nargs = "+")
    parser.add_argument("--spacer", action = "store_true", dest = "SPACER", help = "Flag for using a spacer element too.")
    parser.add_argument("--output", dest = "OUTPUT", help = "Output fasta file to write", type = str)
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

    ## load in 5' and 3' ends, adaptors.
    upstream = args.PAD_RE_5 + args.ADAPTOR5
    downstream = args.ADAPTOR3 + args.RE_PAD_3

    ## load in neutral background sequence (ONLY FIRST SEQUENCE PARSED)
    neutral = [str(record.seq.upper()) for record in Bio.SeqIO.parse(args.NEUTRAL, "fasta")][0]

    ## load in candidates
    candidates = {}
    longest_candidate_length = 0
    summed_candidate_length = 0
    for cand in Bio.SeqIO.parse(args.CANDIDATES, "fasta"):
        cname = cand.name
        cseq = str(cand.seq.upper())
        if len(cseq) > longest_candidate_length:
            longest_candidate_length = len(cseq)
        candidates[cname] = cseq
        summed_candidate_length += len(cseq)
        if str(Seq(cseq).reverse_complement()) != cseq:
            ## is not palindromic
            ## save the reverse compliment too.
            candidates["{}revc".format(cname)] = str(Seq(cseq).reverse_complement())
            summed_candidate_length += len(cseq)    ## length of reverse compliment is the same
    average_candidate_length = math.floor(summed_candidate_length / len(candidates))

    ## check if using a spacer element
    if args.SPACER:
        ## add spacer element too
        candidates["spacer"] = "X"*average_candidate_length

    ## determine oligo length based on longest candidate.
    ## adapter and RE bases + (num_candidates*max_candidate_length) + ((num_candidates+1)*candidate_gap_size)
    oligo_len = len(upstream) + len(downstream) + (args.NUM_CANDIDATES * longest_candidate_length) + ((args.NUM_CANDIDATES + 1) * args.CANDIDATE_GAP)
    statprint("Oligos will be made to be {} bases long to accomodate a candidate gap of {} and {} candidates (max candidate length is {})".format(oligo_len, args.CANDIDATE_GAP, args.NUM_CANDIDATES, longest_candidate_length))
    statprint("Average candidate length, and resulting length of spacer element(s), is {} bases.".format(average_candidate_length))
    statprint("There will be {} oligos created.".format(len(candidates)**args.NUM_CANDIDATES))

    # Generate all combinations 
    # for s in itertools.product(["a","b","c","d"], repeat=3):
    # ...     print(s)
    # ...
    # ('a', 'a', 'a')
    # ('a', 'a', 'b')
    # ('a', 'a', 'c')
    # ('a', 'a', 'd')
    # ('a', 'b', 'a')
    # ('a', 'b', 'b')
    # ('a', 'b', 'c')
    # ('a', 'b', 'd')
    # ...
    # ('d', 'c', 'c')
    # ('d', 'c', 'd')
    # ('d', 'd', 'a')
    # ('d', 'd', 'b')
    # ('d', 'd', 'c')
    # ('d', 'd', 'd')
    SUCCEED = False
    while not SUCCEED:
        SUCCEED=True
        oligos = {}
        num_failed = 0
        #candidates.keys() contains the names of the candidates in the oligo.
        for cand_set in itertools.product(candidates.keys(), repeat=args.NUM_CANDIDATES):
            ## while building these strings, will start with neutral sequence being spaces. Will fill in with actual sequence later.
            ## start with 5' adaptor and leading gap. creating spaces depending on the candidate gap length specified
            oligo_str = args.ADAPTOR5 + " "*args.CANDIDATE_GAP
            oligo_name = []
            ## combine the candidate sequences together and add the specified number of spaces based on the gap length. Do not add 3' adapter yet.
            for cand in cand_set:
                oligo_str += candidates[cand] + " "*args.CANDIDATE_GAP
                oligo_name.append(cand)

            ## create oligo name as string
            oligo_name = "-".join(oligo_name)
            
            #oligo currently includes the 5' adaptor, candidates and the gap sequences. there is no 3' adaptor and 3' RE added yet.
            num_padding_needed = oligo_len - len(args.PAD_RE_5) - len(oligo_str) - len(downstream)  ## desired length, minus 5'RE (not yet added) minus what is built so far, minus downstream (3' adapter and RE)
            oligo_str += " "*num_padding_needed

            ## add 3' adaptors
            oligo_str += args.ADAPTOR3
        

            ## go through the positions in the oligo and find the spaces/Xs. 
            ## if the position is a space or X, find the position index (i)
            ## fill in spaces and Xs (from spacer) with neutral sequence. Start neutral seq at index 0 immediately after adapter.
            oligo_array = [x for x in oligo_str]
            for i in range(len(oligo_array)):
                if oligo_array[i] in [" ", "X"]:
                    neutral_seq_indx=i-len(args.ADAPTOR5) ## this is to offset, so neutral starts immediately after adapter.
                    oligo_array[i] = neutral[neutral_seq_indx]  
            oligo_str = "".join(oligo_array)
            oligos[oligo_name] = args.PAD_RE_5 + oligo_str + args.RE_PAD_3            
            ## check that no forbidden sequence created
            PASS = True
            for fs in args.FORBIDDEN_SEQS:
                if fs in oligo_str:
                    PASS = False
                if str(Seq(fs).reverse_complement()) in oligo_str:
                    PASS = False
            
            if not PASS:
                # statprint("Generated oligo for {} has forbidden sequence {}: {}".format(oligo_name, fs, oligo_str), "ERROR")
                num_failed += 1
                SUCCEED=False
  
        if not SUCCEED:
            statprint(f"Forbidden seq found, shifting neutral by 1...there have been {num_failed} failed oligos so far... :(", "WARNING")
            #remove the first base of neutral sequence
            neutral=neutral[1:]
        


    
    ## write fasta output
    out = open(args.OUTPUT, "w")
    for oligo_name in oligos:
        out.write(">{}.{}\n{}\n".format(oligo_name, os.path.splitext(os.path.basename(args.NEUTRAL))[0].replace(".","_"), oligos[oligo_name]))
    out.close()    

    statprint("Done.")