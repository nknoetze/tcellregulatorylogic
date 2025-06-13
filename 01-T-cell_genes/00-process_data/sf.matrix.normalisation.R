suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(argparse))
suppressMessages(library(DESeq2))

#----set args---
parser <- ArgumentParser(description='normalise counts for size factor');

parser$add_argument(
  '--count_matrix', metavar='count_matrix',
  help='path to the matrix to normalise');

parser$add_argument(
  '--metadata_file_path', metavar='metadata_file_path',
  help='path to DICE and SRA metadata file');

parser$add_argument(
  '--outdir', metavar="outdir",
  help="output directory for the normalised counts");

parser$add_argument(
  '--prefix', metavar="prefix",
  help='prefix to use for the output file. The prefix should describe the file being analysed.')

parser$add_argument(
  '--tmp_outdir', metavar="tmp_outdir",
  help="scratch output directory for the vsd counts");

args=parser$parse_args();
count_matrix_file_path=args$count_matrix
metadata_file_path=args$metadata_file_path
outdir=args$outdir
prefix=args$prefix
tmp_outdir=args$tmp_outdir
#----------------
counts_file <- paste('zcat ',count_matrix_file_path,sep='')
counts_matrix <- fread(cmd=counts_file) %>%
  column_to_rownames(var='gene_id')
metadata <- fread(metadata_file_path)

dds <- DESeqDataSetFromMatrix(counts_matrix,colData = metadata,design=~1)
dds <- estimateSizeFactors(dds)

sf_norm_counts <- counts(dds, normalized = TRUE) %>% 
  as.data.frame() %>%
  rownames_to_column(var='gene_id')

sf_norm_counts_outfile <- paste(outdir,'sra.dice.counts.',prefix,'.sf.adj.tsv.gz',sep='')
fwrite(sf_norm_counts,sf_norm_counts_outfile,sep='\t',quote=FALSE)
