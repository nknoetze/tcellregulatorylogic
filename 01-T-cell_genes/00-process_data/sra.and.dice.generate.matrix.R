#----- Load packages
suppressMessages(library(tximport))
suppressMessages(library(tximportData))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(argparse))
suppressMessages(library(DESeq2))
#----set args---
parser <- ArgumentParser(description='clean up metadata and filter counts matrix');

parser$add_argument(
  '--gencode', metavar='gencode',
  help='path to the gencode file');

parser$add_argument(
  '--DICE_file_path', metavar='dice_file_path',
  help='path to DICE RSEM files');

parser$add_argument(
  '--SRA_file_path', metavar='sra_file_path',
  help='path to SRA RSEM files');

parser$add_argument(
  '--metadata_file', metavar='metadata_file',
  help='path to the dice and sra master metadata file');

parser$add_argument(
  '--outdir', metavar="outdir",
  help="output directory for the count files");

parser$add_argument(
  '--tmp_outdir', metavar="tmp_outdir",
  help="scratch output directory for the vsd counts");

args=parser$parse_args();
gencode_file_path=args$gencode
DICE_rsem_files=args$DICE_file_path
SRA_rsem_files=args$SRA_file_path
file_path=args$file_path
metadata_file=args$metadata_file
outdir=args$outdir
tmp_outdir=args$tmp_outdir

 # ------ READ IN THE DATA FILES
print('reading the data and preparing it with tximport')
metadata <- fread(metadata_file)
sra_metadata <- metadata %>% filter(authors!="DICE")
dice_metadata <- metadata %>% filter(authors=="DICE")

print('reading in the file paths')
file_path_sra <- c(paste0(SRA_rsem_files,"/",sra_metadata$sample_id,"/",sra_metadata$sample_id,'.genes.results',sep=''))

file_path_dice <- c(paste0(DICE_rsem_files,"/",dice_metadata$sample_id,"/",dice_metadata$sample_id,'.genes.results',sep=''))
all_paths <- c(file_path_dice,file_path_sra)
files <- file.path(all_paths)

gencode <- fread(gencode_file_path)
tx2gene <- gencode %>% dplyr::select(transcript_id, gene_id) %>% unique()


# prepare the txi import files. 
print('importing rsem files with tximport')
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

#edit the lengths that are 0 to 1 for it to run in deseq. this happens when the expression is 0 for all samples.
txi_matrix <- txi.rsem[["length"]]
txi_matrix[txi_matrix==0] <- 1 
txi.rsem[["length"]] <- txi_matrix

#------- RUNNING SAMPLES THROUGH DESEQ 
# preparing the sample metadata for DESeq
print('setting up the sample info for deseq')
sample_info<- data.frame(sample_id=metadata$sample_id) %>%
  mutate(sample=sample_id) %>%
  column_to_rownames(var="sample") 
# getting the size factors for the normalised counts 
print('importing data with DESeq')
DESeq.ds <- DESeqDataSetFromTximport(txi.rsem,colData = sample_info,design=~1)
# remove genes that have a count of 0 across all samples
keep <- rowSums(counts(DESeq.ds)) > 0
dds <- DESeq.ds[keep,]

# #get the raw counts
print('getting the un-normalised counts')
raw_counts <- counts(dds, normalized = FALSE) %>% 
  as.data.frame() %>%
  rownames_to_column(var='gene_id')

#----- ASSIGN OUTPUT VARIABLES
raw_counts_output <- paste(outdir,'sra.dice.counts.raw.tsv.gz',sep='')

# ----- WRITE THE OUTPUT FILES
# write the raw counts to a tsv file and update the metadata file
print('writing un-normalized counts')
fwrite(raw_counts, raw_counts_output,sep='\t',quote=FALSE)
