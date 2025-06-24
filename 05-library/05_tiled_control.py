TITLE = "Tiled Control"
DESC = "Given control sequence, create all possible tiled sequences of a given length"
'''
Date: December 15,2023
@author: nknoetze
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
    parser.add_argument("--control_seq", dest = "CONTROL", help = "control sequence, fasta.", type = str)
    parser.add_argument("--pad_re_5", dest = "PAD_RE_5", help = "5' padding and restriction enzyme site sequence", type = str)
    parser.add_argument("--re_pad_3", dest = "RE_PAD_3", help = "3' restriction enzyme site and padding sequence", type = str)
    parser.add_argument("--adaptor_5", dest = "ADAPTOR5", help = "5' Illumina adaptor sequence", type = str)
    parser.add_argument("--adaptor_3", dest = "ADAPTOR3", help = "3' Illumina adaptor sequence", type = str)
    parser.add_argument("--tile_length", dest = "TILE_LENGTH", help = "length of tiled sequence to create", type = int)
    parser.add_argument("--forbidden_seqs", dest = "FORBIDDEN_SEQS", help = "Forbidden sequence(s) to not allow in oligos", type = str, nargs = "+")
    parser.add_argument("--output", dest = "OUTPUT", help = "Output fasta file to write", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB
    tile_length = args.TILE_LENGTH

    print("=======================================================")
    print("Python version: {}".format(sys.version))
    print("Python environment: {}".format(sys.prefix))
    print("Server: {}".format(os.uname()[1]))
    print("Current directory: {}".format(os.getcwd()))
    print("Command: {}".format(" ".join(sys.argv)))
    print("Time: {}".format(time.strftime("%Y/%m/%d %T")))
    print("=======================================================\n")


    ## Script content here
    ## load in control sequence (ONLY FIRST SEQUENCE PARSED)
    ## get the length of the control sequence, and the sequence name
    control = [str(record.seq.upper()) for record in Bio.SeqIO.parse(args.CONTROL, "fasta")][0]
    control_name=[str(record.description) for record in Bio.SeqIO.parse(args.CONTROL, "fasta")][0]   
    control_seq_length=len(control)
    
    tiled_seqs={}
    #obtain each subsequence for the control sequence based onthe specified tile length
    for i in range(0,control_seq_length-tile_length+1):
        succeed="TRUE"
        tiled_seq=control[i:i+tile_length]
        oligo_seq=args.ADAPTOR5+tiled_seq+args.ADAPTOR3
    # check if there are any forbidden sequences
        for fs in args.FORBIDDEN_SEQS:
            if fs in oligo_seq or str(Seq(fs).reverse_complement()) in oligo_seq:
                succeed="FALSE"
                print('forbidden sequence found.')
        if succeed=="TRUE":
            oligo_name=('{}.{}').format(i+1,control_name)
            tiled_seqs[oligo_name]=args.PAD_RE_5+oligo_seq+args.RE_PAD_3


# write fasta output
out = open(args.OUTPUT, "w")
for oligo_name in tiled_seqs:
    out.write(">{}\n{}\n".format(oligo_name, tiled_seqs[oligo_name]))
out.close()    

statprint("Done.")