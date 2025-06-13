### Cutoff finder for maximum discrimination between T cell and non-T cell genes.
## June 1, 2023
# sbrown
suppressMessages(library(argparse))
### ----------------------------- ###
###        SET UP ARGUMENTS       ###
### ----------------------------- ###
# create parser object
parser <- ArgumentParser()
parser$add_argument("--gene_file", help="Path to the gene file to filter")
parser$add_argument("--additional_gene_file", help="Use only if using old method", default=FALSE)
parser$add_argument("--outdir", help="Path to the output directory")
parser$add_argument("--sigma_threshold", help="sigma threshold to use for target genes", type="integer")

args <- parser$parse_args()
additional_gene_file <- args$additional_gene_file
gene_file <- args$gene_file
outdir <- args$outdir
sigma_threshold <- args$sigma_threshold

##GENES THAT HAVE AT LEAST 2 OCR TYPES. LIST STILL INCLUDES GENES WITHOUT ANY OCRS :(
#genes <- read.table("/projects/nknoetze_prj/ocr_prj/data/processed/gene_lists/ranked_gene_list_missingocrgenes_removed.tsv", sep="\t", stringsAsFactors = F, header=T)
print('Reading in ranked gene list')
genes <- read.table(gene_file, sep="\t", stringsAsFactors = F, header=T)
genes <- genes[order(genes$t_o_rank),]

##rank is much more normally distributed than sum.
## we want genes with low trank and low orank
print('Getting thresholds')
mn <- mean(genes$t_o_rank)
sigma <- sd(genes$t_o_rank)
t_o_rank_thresh_targets <- mn - (sigma_threshold*sigma)
t_o_rank_thresh_background <- mn - (((sigma_threshold)-1)*sigma)
#### Write out gene lists ####
if(additional_gene_file!=FALSE){
    print('Requested old method')
    #genes <- read.table("/projects/nknoetze_prj/ocr_prj/data/processed/gene_lists/ranked_gene_list_missingocrgenes_removedv2.tsv", sep="\t", header=T, stringsAsFactors = F)
    additional_genes <- read.table(additional_gene_file, sep="\t", header=T, stringsAsFactors = F)
    target_genes <- subset(additional_genes, t_o_rank <= t_o_rank_thresh_targets)$gene_id
    background_genes <- subset(additional_genes, t_o_rank > t_o_rank_thresh_background)$gene_id
} else{
    target_genes <- subset(genes, t_o_rank <= t_o_rank_thresh_targets)$gene_id
    background_genes <- subset(genes, t_o_rank > t_o_rank_thresh_background)$gene_id
}
## 4 sigma targets
write.table(target_genes, paste(outdir,sigma_threshold,'sigma_target_genes.txt',sep=''), col.names = F, row.names=F, quote=F)
## 4 sigma ZoE 3 sigma background
write.table(background_genes, paste(outdir,(sigma_threshold)-1,'sigma_background_genes.txt',sep=''), col.names = F, row.names=F, quote=F)
