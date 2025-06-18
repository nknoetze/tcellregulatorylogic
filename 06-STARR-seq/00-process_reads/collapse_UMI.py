TITLE = "Collapse UMIs"
DESC = "ollapse UMIs with a specified hamming distance into a single UMI. Count unique UMIs for each oligo"
'''
Date: January 2024
@author: nknoetze
'''

##########################################
###           IMPORT LIBRARIES         ###
##########################################
from umi_tools import UMIClusterer
import pandas as pd
import argparse

##########################################
###             ARGUMENTS              ###
##########################################
## Deal with command line arguments
parser = argparse.ArgumentParser(description = TITLE)
parser.add_argument("--UMI_bedfile", dest = "UMI_FILE", help = "Path to file containing UMI information and counts", type = str)
parser.add_argument("--mm_threshold", dest = "MM_THRESHOLD", help = "Mistmatch threshold to use in collapsing", type = int)
parser.add_argument("--outfile", dest = "OUTFILE", help = "Path and name of file to write output to", type = str)
args = parser.parse_args()

##########################################
###           Parse Input              ###
##########################################
UMI_dict={}
UMI_list=list()
with open(args.UMI_FILE,"r") as file:
  for line in file:
    #candidate, _, _, UMI, UMI_count,_ = line.rstrip().split('\t')
    candidate, UMI, UMI_count = line.rstrip().split('\t')
    #algorithm require bytes to calculate hamming distance
    UMI = bytes(UMI, "ascii")
    UMI_count = int(UMI_count)
    # Check if candidate is already in UMI_dict
    if candidate not in UMI_dict:
        # If not, create a new dictionary for the candidate
        UMI_dict[candidate] = {}
    # Store UMI and UMI count in the dictionary for the candidate
    UMI_dict[candidate][UMI] = UMI_count


##########################################
###             Collapse UMI           ###
##########################################
#set up cluster algorithm
print('Clustering UMIs with a distance of {}'.format(args.MM_THRESHOLD))
clusterer = UMIClusterer(cluster_method='cluster')
UMI_count_file=open(args.OUTFILE,'w')
UMI_count_file.write('candidate\tcount\n')
#perform collapsing for each candidate independently.
for candidate in UMI_dict.keys():
    candidate_dict=UMI_dict[candidate]
    clustered_umis=clusterer(candidate_dict,threshold=args.MM_THRESHOLD)
    n_unique_UMI=len(clustered_umis)
    UMI_count_file.write('{}\t{}\n'.format(candidate,n_unique_UMI))
UMI_count_file.close()

