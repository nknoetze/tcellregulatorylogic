TITLE = "Kmer Clusters to MEME"
DESC = "Given a file with kmers in groups, generate a MEME format file with one MEME motif for each cluster."
'''
Date: June 8, 2023
@author: sbrown
'''

## Import Libraries

import sys
import argparse
import os
import time


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
    parser.add_argument("--kmer_clusters", dest = "KMER", help = "Input kmer clusters file", type = str)
    parser.add_argument("--meme_output", dest = "MEME", help = "Output ", type = str)
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
    
    motifs = {}
    counts = {}

    HEADER = True
    for line in open(args.KMER, "r"):
        if HEADER:
            HEADER = False
        else:
            group, kmer = line.rstrip().split("\t")
            if group not in counts:
                counts[group] = [{"A":0, "C":0, "G":0, "T":0} for x in kmer]  ## make list of dictionaries, list length == kmer length.
            
            for i in range(len(kmer)):
                counts[group][i][kmer[i]] += 1

    ## write results, normalizing counts to fractions.

    out = open(args.MEME, "w")
    out.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
    for group in counts:
        ## get number of kmers that contributed
        num_kmer = sum([counts[group][0][x] for x in "ACGT"])
        kmer_len = len(counts[group])

        out.write("MOTIF CLUST{} CLUST{}.{}\n".format(group, group, num_kmer))
        out.write("letter-probability matrix: alength= 4 w= {} nsites= {} E= 0\n".format(kmer_len, num_kmer))
        for pos in range(len(counts[group])):
            out.write(" {} {} {} {}\n".format(counts[group][pos]["A"]/num_kmer,
                                              counts[group][pos]["C"]/num_kmer,
                                              counts[group][pos]["G"]/num_kmer,
                                              counts[group][pos]["T"]/num_kmer))
        out.write("\n")




    

    statprint("Done.")