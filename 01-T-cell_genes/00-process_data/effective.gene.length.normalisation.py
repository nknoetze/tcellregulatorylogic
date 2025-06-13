import argparse
import pandas as pd
import numpy as np

#----set args---

parser=argparse.ArgumentParser(description='normalise counts for size factor and/or gene length')

parser.add_argument('--count_matrix',metavar='-count_matrix',help='path to the matrix to normalise',type=str,required=True)
parser.add_argument('--metadata_file_path',metavar='-metadata_file_path',help='path to DICE and SRA metadata file',type=str,required=True)
parser.add_argument('--outdir',metavar='-outdir',help='output directory for the normalised counts',type=str,required=True)
parser.add_argument('--effective_lengths',metavar='-effective_lengths',help='file containing the effective gene lengths for each gene and sample"',type=str,required=True)
parser.add_argument('--prefix',metavar='-prefix',help='prefix to use for the output file. The prefix should describe the file being analysed.',type=str,required=True)
args=parser.parse_args()

#--- ASSIGN ARGS TO GLOBAL VARS 
count_matrix=args.count_matrix
metadata_file_path=args.metadata_file_path
outdir=args.outdir
effective_lengths=args.effective_lengths
prefix=args.prefix


#----------------
counts_matrix = pd.read_csv(count_matrix,sep='\t',compression='gzip')
effective_gene_lengths = pd.read_csv(effective_lengths, sep='\t', names=['gene_id','effective_gene_length','sample_id'])
metadata = pd.read_csv(metadata_file_path,sep='\t')

#for each sample and gene, normalise the expression value by the effective gene length
#Use a factor of 1000 to as a scaling factor. Helps with visualisation.
adjusted_counts = (pd.melt(counts_matrix,id_vars=['gene_id']).rename(columns={'variable':'sample_id'})
                  .merge(effective_gene_lengths,how='inner')
                  .assign(adjusted_count=lambda x:((x['value']/x['effective_gene_length'])*1000))
                  .drop(['effective_gene_length','value'],axis=1)
                  .pivot(index='gene_id',columns='sample_id',values='adjusted_count')
                  .reset_index().rename_axis(None,axis=1))

adjusted_counts_outfile = '{}/sra.dice.counts.{}.sf.genel.adj.tsv.gz'.format(outdir,prefix)
adjusted_counts.to_csv(adjusted_counts_outfile,sep='\t',doublequote=False,compression='gzip',index=False)
