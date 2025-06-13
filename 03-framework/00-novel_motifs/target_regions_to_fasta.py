TITLE = "Create FASTA file from Regions"
DESC = "Given a set of regions and target genes, generate fasta file."
'''
Date: July 4,2024
@author: nknoetze, sbrown
'''

## Import Libraries

import sys
import argparse
import os
import time
import pybedtools


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
    parser.add_argument("--regions", dest = "REGIONS", help = "Regions file containing OCRs/Peaks", type = str)
    parser.add_argument("--genes", dest = "GENES", help = "Set of genes to use for filtering", type = str)
    parser.add_argument("--output", dest = "OUTPUT", help = "Output fasta file to write", type = str)
    parser.add_argument("--ref_fasta", dest = "REF_FASTA", help = "Path to the genome reference file", type = str)
    parser.add_argument("--region_type", dest = "REGION_TYPE", help = "all, non-promoter, or promoter", type = str)
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

    ## load in gencode
    REF_GENOME_FASTA = pybedtools.BedTool(args.REF_FASTA)

    target_genes = []
    for line in open(args.GENES, "r"):
        target_genes.append(line.rstrip())

    regions = {}

    header = None
    for line in open(args.REGIONS, "r"):
        if not header:
            header = line.rstrip().split("\t")
        else:
            ## not header line
            # chrom   start   end     ocr_id  score   strand  ocr_type        query_gene      promoter_id     query_name      gene_group      data    dist_tss        promoter_name   promoter_type
            linedict = {x:y for x,y in zip(header, line.rstrip().split("\t"))}
            region_id = "{}\t{}".format(linedict["ocr_id"].replace(":","\t").replace("-","\t"), linedict["ocr_id"])
            ## region type - need to compare query_gene and promoter_id.
            promoter_id = linedict["promoter_id"]
            gene_id = linedict["query_gene"]
            region_type = linedict["ocr_type"]

            ## only use regions for genes of interest!
            if gene_id in target_genes:

                ## if promoter but not of gene in question
                if region_type == "promoter" and gene_id != promoter_id:
                    ## not promoter for this gene, so it is what we are calling a "linked promoter".
                    region_type = "linked_promoter"

                if args.REGION_TYPE == "all":
                    ## we will track all regions and their region_type for this set. Priority of region_type is promoter > linked_promoter > non-promoter
                    ## any region will either be non-promoter, or promoter/linked_promoter. 

                    ## NOTE April 29, 2024 - Updating to not include linked_promoters in the 'all' set. This should have been done originally.
                    # NOTE Sept 23, 2024 - Reverting back to including linked-promoters for comparing ATAC and DNAse datasets.
                    #if region_type != "linked_promoter":
                
                    if region_id not in regions:
                        regions[region_id] = region_type
                    else:
                        if regions[region_id] != "promoter":
                            ## update - could potentially be promoter now. Note this may result in no change.
                            regions[region_id] = region_type
                else:
                    ## only track regions of interest
                    if region_type == args.REGION_TYPE:
                        regions[region_id] = region_type

    ## only take unqiue regions, because if we are looking at the same OCR, we don't want to double count all the kmers in that OCR.
    regions = [x for x in set(regions)]

    ## convert to pybedtools object
    bed_regions = pybedtools.BedTool("\n".join([r for r in regions]), from_string=True)
    ## ONLY WRITE ONE SEQUENCE FOR DEBUGGING
    ##bed_regions = pybedtools.BedTool("\n".join([list(regions.keys())[0]]), from_string=True)

    ## get sequences
    a = bed_regions.sequence(fi=REF_GENOME_FASTA)     ## pybedtools sequence (NOT temp file)
    region_seqs = open(a.seqfn).read()
    ## write file.
    out = open(args.OUTPUT, "w")
    out.write(region_seqs)
    out.close()