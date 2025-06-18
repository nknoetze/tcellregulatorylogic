### --------------------------------------- ###
###   Load source scripts&libraries         ###
### --------------------------------------- ###
source('/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/data/file.metrics.fx.R')
### ----------------------------- ###
###        SET UP ARGUMENTS       ###
### ----------------------------- ###
# create parser object
parser <- ArgumentParser()
parser$add_argument("-gen", "--gene_file", help="Path to the gene list")
parser$add_argument("-gc", "--gencode", help="Path to the parse gencode v19 file")
parser$add_argument("-ocr", "--linked_ocrs", help="Path to the linked OCR file")
parser$add_argument("-o", "--outdir", help="Path to the output directory")

args <- parser$parse_args()
gene_file <- args$gene_file
gencode_file <- args$gencode
linked_ocr_file <- args$linked_ocrs
outdir <- args$outdir

### ----------------------------- ###
###        READ IN FILES          ###
### ----------------------------- ###
print("Reading in the Data")

gene_list <- fread(gene_file) %>%
select(gene_id, gene_group,gene_name)
types <- c('protein_coding','TR_C_gene','IG_C_gene')

gencode <- fread(gencode_file) %>% 
  filter(type=='transcript',gene_type %in% types, transcript_type %in% types) %>% 
  filter(!(grepl("chrY|chrM",chrom))) 

#Get the gene name and id's for later :) 
gencode_genes <- gencode %>% select(gene_name,gene_id) %>% unique() %>% left_join(gene_list) %>%
  distinct(gene_id,gene_name)

#Get the TSS for each promoter based on the strand
#REMOVE TRANSCRIPT ID. For a single gene, different transcripts can have the same TSS
#So only keep the gene id, and the non-redundnant TSS
tss <- gencode %>% mutate(tss=ifelse(strand=='+',start,end)) %>%
  transmute('promoter_id'=gene_id,tss,'promoter_name'=gene_name) %>% unique()

#Get the linked OCRs and add the gencode gene information (gene name and gene group)
#based on the query gene (gene id)
linked_ocrs <- fread(linked_ocr_file) %>%
  rename('gene_id'=query_gene) %>% 
  inner_join(gencode_genes,by='gene_id')

##########################################
##             OCR-TSS DIST             ##
##########################################
print('Calculating distance of OCRs to nearest TSS')
ocr_tss_dist <- get_ocrs_dists(linked_ocrs,gencode_genes)

##########################################
##         WRITE OUTPUT FILES           ##
##########################################
fwrite(ocr_tss_dist,paste(outdir,'pulled_ocrs_tss_dist.tsv',sep=''),sep='\t',quote=FALSE)