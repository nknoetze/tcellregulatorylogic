TITLE = "OBTAIN UMI BED"
DESC = "Create a UMI BED file containing the UMI, candidate, and counts."
'''
Date: June 2024
@author: nknoetze
'''

##########################################
###           IMPORT LIBRARIES         ###
##########################################
import polars as pl
import argparse

##########################################
###             ARGUMENTS              ###
##########################################
## Deal with command line arguments
parser = argparse.ArgumentParser(description = TITLE)
parser.add_argument("--SAM_file", dest = "SAM_FILE", help = "Path to sam file containing alignments", type = str)
parser.add_argument("--outfile", dest = "OUTFILE", help = "Path and name of file to write output to", type = str)
args = parser.parse_args()

##########################################
###           Parse Input              ###
##########################################
print('reading in sam file')
sam_file = pl.read_csv(args.SAM_FILE,separator='\t',has_header=False,
    new_columns=["qname", "candidate"])

print('parsing sam file')
filtered_sam=(sam_file.with_columns([
            pl.col("qname").str.slice(-10).alias("UMI"),
            pl.col("qname").str.head(-11).alias("read")])
            .select(["candidate", "read", "UMI"]).unique())
            #pl.col("tlen").abs(),#.alias("stop"),
            #pl.lit("*").alias("strand"),
            #pl.lit("1").alias("start")])
            #.select(["candidate", "read", "UMI"]).unique())
            #.select(["candidate", "read", "start","tlen", "UMI", "strand"]).unique())
print('getting UMI counts')
umi_bedfile = (
    filtered_sam
    .group_by(["candidate", "UMI"])  # Group by UMI and candidate
    .len() # n_unique()????
    #.join(filtered_sam, on=["UMI", "candidate"]), how="inner")  # Merge with original DataFrame to retain original columns
    #.select(["candidate", "start", "tlen", "UMI", "len", "strand"])
    .unique()) 

#filtered_sam.write_csv("/projects/nknoetze_prj/ocr_prj/parsed_file.tsv",include_header=True,separator="\t")
umi_bedfile.write_csv(args.OUTFILE,include_header=False,separator="\t")
