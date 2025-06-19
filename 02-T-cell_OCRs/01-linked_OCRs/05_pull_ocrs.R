suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))

### ----------------------------- ###
###        SET UP ARGUMENTS       ###
### ----------------------------- ###
# create parser object
parser <- ArgumentParser()
parser$add_argument("-gc", "--gencode", help="Path to the parse gencode v19 file")
parser$add_argument("-ocr", "--filtered_ocr_file", help="Path to the filtered T-cell OCR dataset")
parser$add_argument("-int", "--interaction_file", help="Path to the annotated interaction file")
parser$add_argument("-gen", "--gene_file", help="Path to the gene list")
parser$add_argument("-o", "--outdir", help="Path to the output directory")


args <- parser$parse_args()
gencode <- args$gencode
filtered_ocr_file <- args$filtered_ocr_file
interaction_file <- args$interaction_file
gene_file <- args$gene_file
outdir <- args$outdir

### ----------------------------- ###
###        READ IN FILES          ###
### ----------------------------- ###
print("Reading in the Data")
#gencode <- fread('/projects/nknoetze_prj/references/annotations/gencode/gencode.v19.transc.type.1based.tsv') %>% select(gene_name,gene_id) %>% unique()
gencode <- fread(gencode) %>% select(gene_name,gene_id) %>% unique()
#gene_list <- fread('/projects/nknoetze_prj/ocr_prj/data/processed/gene_lists/ranked_gene_list.tsv')
gene_list <- fread(gene_file)

#ocr_file_name <- sub('\\..*$', '', basename('/projects/nknoetze_prj/ocr_prj/data/processed/meuleman_dnase/meuleman_tcell_2sample_1cd8_1cd4_ocrs_merged_annotated_tidy.tsv'))
ocr_file_name <- sub('\\..*$', '', basename(filtered_ocr_file))
ocr_file_name <- sub('_annotated', '', ocr_file_name)
#ocrs <- fread('/projects/nknoetze_prj/ocr_prj/data/processed/meuleman_dnase/meuleman_tcell_2sample_1cd8_1cd4_ocrs_merged_annotated_tidy.tsv') %>%
ocrs <- fread(filtered_ocr_file) %>%
  select(-width,-strand)
#interactions <- fread('/projects/nknoetze_prj/ocr_prj/data/processed/interactions/DICE/shared/Loop_2.5kb/connectivity8/merged_interactions_tss_within.tsv')
interactions <- fread(interaction_file)

### ----------------------------- ###
###     Prepare the Data          ###
### ----------------------------- ###
ocr_ranges <- makeGRangesFromDataFrame(ocrs,keep.extra.columns = TRUE,ignore.strand = TRUE,seqnames.field = c('chrom'))

### ------------------------------------------------- ###
###  Pull out the OCRs for genes with interactions    ###
### ------------------------------------------------- ###
print("Getting significant interactions for genes in the genome")
filtered_interactions <- interactions %>%
  #File is already filtered!
  #get the interactions where the gene's promoter is within one of the two bins
  filter(bin1_gene_id %in% gencode$gene_id|bin2_gene_id %in% gencode$gene_id) %>%
  separate_rows(bin1_gene_id, sep=";") %>%
  separate_rows(bin2_gene_id,sep=';') 

bin1_ranges <- makeGRangesFromDataFrame(filtered_interactions,seqnames.field = 'chrom1',start.field = 's1',end.field = 'e1',keep.extra.columns = TRUE,ignore.strand = TRUE)
bin2_ranges <- makeGRangesFromDataFrame(filtered_interactions,seqnames.field = 'chrom2',start.field = 's2',end.field = 'e2',keep.extra.columns = TRUE,ignore.strand = TRUE)

### ------------------------------------------------- ###
###           FIND OVERLAPPING HI-C BIN:OCRS          ###
### ------------------------------------------------- ###
###                INTERACTION BIN 1                  ###
### ------------------------------------------------- ###
#For each interaction bin, find it it overlaps with an OCR
#query:ocrs, subject: interactions
#keep all overlapping ocrs
#min overlap of 1
#any type of overlap ok
print('Finding overlaps between OCRs and Hi-C Bin 1')
ocr_overlaps_bin1 <- findOverlaps(query = ocr_ranges,subject = bin1_ranges,type = 'any',select='all',minoverlap = 1)
mapping_df1 <- DataFrame(interaction_indx=subjectHits(ocr_overlaps_bin1),ocr_indx=queryHits(ocr_overlaps_bin1)) %>% as.data.frame()
bin1_ocr_interactions <- mapping_df1 %>%
  mutate(ocr_id=ocr_ranges[ocr_indx]$ocr_id,interaction_id=bin1_ranges[interaction_indx]$interaction_id,ocr_type=ocr_ranges[ocr_indx]$ocr_type,
         ocr_gene_id=ocr_ranges[ocr_indx]$gene_id,interaction_bin1_gene_id=bin1_ranges[interaction_indx]$bin1_gene_id,ocr_bin1_promoter_id=ocr_ranges[ocr_indx]$gene_id,
         ocr_bin1_transcript_id=ocr_ranges[ocr_indx]$transcript_id) %>%
  transmute('bin1_ocr_id'=ocr_id,interaction_id,interaction_bin1_gene_id,ocr_gene_id,ocr_type,ocr_bin1_promoter_id,ocr_bin1_transcript_id) %>% unique()

### ------------------------------------------------- ###
###                INTERACTION BIN 2                  ###
### ------------------------------------------------- ###
print('Finding overlaps between OCRs and Hi-C Bin 2')
ocr_overlaps_bin2 <- findOverlaps(query = ocr_ranges,subject = bin2_ranges,type = 'any',select='all',minoverlap = 1)
mapping_df2 <- DataFrame(interaction_indx=subjectHits(ocr_overlaps_bin2),ocr_indx=queryHits(ocr_overlaps_bin2)) %>% as.data.frame()
bin2_ocr_interactions <- mapping_df2 %>%
  mutate(ocr_id=ocr_ranges[ocr_indx]$ocr_id,interaction_id=bin2_ranges[interaction_indx]$interaction_id,ocr_type=ocr_ranges[ocr_indx]$ocr_type,
         ocr_gene_id=ocr_ranges[ocr_indx]$gene_id,interaction_bin2_gene_id=bin2_ranges[interaction_indx]$bin2_gene_id,ocr_bin2_promoter_id=ocr_ranges[ocr_indx]$gene_id,
         ocr_bin2_transcript_id=ocr_ranges[ocr_indx]$transcript_id) %>%
  transmute('bin2_ocr_id'=ocr_id,interaction_id,interaction_bin2_gene_id,ocr_gene_id,ocr_type,ocr_bin2_promoter_id,ocr_bin2_transcript_id) %>% unique()

#Combine interactions & pulled ocrs for each interaction bin.
#NOTE: DOES NOT CONTAIN PROMOTER OCRS THAT ARE UNLINKED.
interacting_ocr_df <- inner_join(bin1_ocr_interactions,bin2_ocr_interactions,by='interaction_id') %>%
  rename('ocr_bin1_gene_id'=ocr_gene_id.x,'ocr_bin2_gene_id'=ocr_gene_id.y,'ocr_bin1_type'=ocr_type.x,'ocr_bin2_type'=ocr_type.y)
print('Done getting overlaps.')
### ------------------------------------------------- ###
###    OUTPUT FILE WITH LINKED OCRS FOR EACH GENE     ###
### ------------------------------------------------- ###
# For each gene in our list, get all the OCRs
#if there are no OCRs overlapping a hi-c bin , pull the gene's promoter OCR
#Create an output file that looks like a bedfile.
#output file contains extra information such as the ocr type, the query gene, and the name
#of the gene for the promoter if it is a promoter ocr
#file can be used as a bed file input for many downstream tools.
#Gets the genes that have ocrs annotated for it to reduce the amount of looping we have to do.
genes_w_ocrs <- gencode %>% filter(gene_id %in% ocrs$gene_id)
ocrs_by_gene <- list()

for(query_gene in genes_w_ocrs$gene_id){
  gene_interactions <- interacting_ocr_df %>% 
    filter(interaction_bin1_gene_id==query_gene | interaction_bin2_gene_id==query_gene) %>%
    mutate(query_gene=query_gene)
  if(nrow(gene_interactions)==0){
    print(paste("No interactions for ",query_gene,". Pulling promoter OCRs.",sep=''))
    ocr_df <- ocrs %>% filter(gene_id==query_gene) %>% 
      transmute(chrom,start,end,ocr_id,score='.',strand='*',ocr_type,query_gene,'promoter_id'=gene_id,transcript_id)
    ocrs_by_gene[[query_gene]] <- ocr_df
  } else{
    print(paste("Pulling OCRs for ",query_gene,sep=''))
    ocr_df <- gene_interactions %>% 
      mutate(bin1=paste(ocr_bin1_type,ocr_bin1_promoter_id,ocr_bin1_transcript_id,bin1_ocr_id,sep=";"),bin2=paste(ocr_bin2_type,ocr_bin2_promoter_id,ocr_bin2_transcript_id,bin2_ocr_id,sep=";")) %>% 
      select(bin1,bin2,query_gene) %>% unique() %>% pivot_longer(cols=contains('bin')) %>% 
      separate(value,into=c('ocr_type','promoter_id','transcript_id','ocr_id'),sep=";") %>% 
      mutate(ocr_coord=ocr_id) %>% select(-name) %>% unique() %>%
      separate(ocr_coord,into=c('chrom','start','end'),sep=":|-") %>%
      transmute(chrom,start,end,ocr_id,score='.',strand='*',ocr_type,query_gene,promoter_id,transcript_id)
    ocrs_by_gene[[query_gene]] <- ocr_df
  }
}

ocrs_by_gene_df_temp <- rbindlist(ocrs_by_gene) 
### ------------------------------------------------- ###
###       Subset Pulled OCRs to gene Lists            ###
### ------------------------------------------------- ###
print('Reducing master list to gene in our list of interest')
ocrs_by_gene_df_subset <- ocrs_by_gene_df %>% filter(query_gene %in% gene_list$gene_id)
#NOTE: DOES NOT CONTAIN PROMOTER OCRS THAT ARE UNLINKED
interacting_ocr_df_subset <- interacting_ocr_df %>% filter(interaction_bin1_gene_id %in% gene_list$gene_id | interaction_bin2_gene_id %in% gene_list$gene_id)

### ------------------------------------------------- ###
###           Write out the output files              ###
### ------------------------------------------------- ###
#NOTE: INTERACTIONS FILE DOES NOT CONTAIN PROMOTER OCRS THAT ARE UNLINKED.
fwrite(interacting_ocr_df,paste(outdir,'pulled_ocrs_master_interactions.tsv.gz',sep=''),sep='\t',quote=FALSE)
ocrs_by_gene_df %>% fwrite(paste(outdir,'pulled_ocrs_master_by_gene.tsv',sep=''),sep='\t',quote=FALSE)
#NOTE: INTERACTIONS FILE DOES NOT CONTAIN PROMOTER OCRS THAT ARE UNLINKED.
fwrite(interacting_ocr_df_subset,paste(outdir,'pulled_ocrs_interactions.tsv.gz',sep=''),sep='\t',quote=FALSE)
fwrite(ocrs_by_gene_df_subset,paste(outdir,'pulled_ocrs_by_gene.tsv.gz',sep=''),sep='\t',quote=FALSE)
print('done!')