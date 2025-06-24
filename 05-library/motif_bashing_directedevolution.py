TITLE = "Motif Bashing Directed Evolution"
DESC = "Repeatedly mutate random sequences until no motifs present."
'''
Date: November 6, 2023
@author: sbrown
'''

## Import Libraries

import sys
import argparse
import os
import time
import random
import pandas as pd
import subprocess
import warnings


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
    parser.add_argument("--motifs", dest = "MOTIFS", help = "Motif file to use", type = str,default='/projects/sbrown_prj/220318_DI/data/interim/motifs/JASPAR+novel+histone_motifs.meme')
    parser.add_argument("--tmpdir", dest = "TMP", help = "Temp dir to write to", type = str,default='./')
    parser.add_argument("--seqlen", dest = "SEQLEN", help = "Length of sequence to generate", type = int,default=125)
    parser.add_argument("--output", dest = "OUTPUT", help = "Output file to write", type = str,default='motif+novel-free_125_2.fasta')
    parser.add_argument("--seed", dest = "SEED", help = "Random Seed", type = int,default=59)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    random.seed(args.SEED)

    SCORE_THRESH = 11.16

    print("=======================================================")
    print("Python version: {}".format(sys.version))
    print("Python environment: {}".format(sys.prefix))
    print("Server: {}".format(os.uname()[1]))
    print("Current directory: {}".format(os.getcwd()))
    print("Command: {}".format(" ".join(sys.argv)))
    print("Time: {}".format(time.strftime("%Y/%m/%d %T")))
    print("=======================================================\n")


    fimo_outdir = os.path.join(args.TMP, "tmp") ## this will get repeatedly overwritten with each iteration

    motif_file = args.MOTIFS

    fasta_file = os.path.join(args.TMP, "sequence.fasta")

    num_sites = 900
    num_tested = 0

    ## randomly generate nucleotide string
    seq = "".join(random.choices(["A","C","T","G"], k=args.SEQLEN))


    while num_sites > 0:
        
        num_tested += 1
        out = open(fasta_file, "w")
        out.write(f">sequence\n{seq}\n")
        out.close()

        ## run fimo
        fimo_cmd='/projects/nknoetze_prj/miniconda3/envs/ocr_prj/bin/fimo --verbosity 1 --thresh 0.0001 --oc {} {} {}'.format(fimo_outdir,motif_file, fasta_file)
        #Run FIMO on the sequences
        VALID_COM = False
        while not VALID_COM:
            call = subprocess.Popen(fimo_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            (res, err) = call.communicate()
            if err.decode('ascii') == "":
                VALID_COM = True
            else:
                print("ERROR IN FIMO: {}\nError: {}".format(fimo_cmd, str(err.decode('ascii'))))
                raise Exception()
        
        try:
            fimo_results=pd.read_csv('{}/fimo.tsv'.format(fimo_outdir),sep='\t',index_col=False, comment="#")
        except:
            print("ERROR")
            print(fimo_outdir)
            print("fasta:")
            for line in open(fasta_file, "r"):
                print(line)
            print("fimo:")
            for line in open('{}/fimo.tsv'.format(fimo_outdir), "r"):
                print(line)
            sys.exit(1)
        #filtering for TFBS predictions with a score greater than 11.16
        fimo_results=(fimo_results[~fimo_results.motif_id.str.startswith('#')].rename(columns={'sequence_name':'gene_id'}).query('`score`>11.16').drop_duplicates())

        ## get number of result lines.
        num_sites = len(fimo_results.index)

        statprint(f"Sequence {num_tested} has {num_sites} motif site(s).")

        if num_sites > 0:
            ## mutate the sequence at these site(s)
            statprint("Mutating...")
            for index, row in fimo_results.iterrows():
                ## fimo result positions are 1-based
                start = int(row["start"])
                end = int(row["stop"])
                seq = seq[:start-1] + "".join(random.choices(["A","C","T","G"], k=end-start+1)) + seq[end:]
    
    ## hooray, if we are here that means we've found a sequence.
    statprint(f"Sequence: {seq}")
    out = open(args.OUTPUT, "w")
    out.write(f">seq_with_no_motif_sites\n{seq}\n")
    out.close()

    
    

    statprint("Done.")