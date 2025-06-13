suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

### ----------------------------- ###
###        SET UP ARGUMENTS       ###
### ----------------------------- ###
# create parser object
parser <- ArgumentParser()
parser$add_argument("-int", "--interaction_file", help="Path to the annotated unmerged interaction file")
parser$add_argument("-mergedint", "--merged_interaction_file", help="Path to the annotated merged interaction file")
parser$add_argument("-o", "--outdir", help="Path to the output directory")

args <- parser$parse_args()
interaction_file <- args$interaction_file
merged_interaction_file <- args$merged_interaction_file
outdir <- args$outdir

### ----------------------------- ###
###        READ IN FILES          ###
### ----------------------------- ###
#get interactions for cd4 or cd8 tcells
#merged_interactions <- fread('/projects/nknoetze_prj/ocr_prj/data/processed/interactions/DICE/shared/Loop_2.5kb/connectivity8/merged_interactions_tss_within.tsv')%>% 
merged_interactions <- fread(merged_interaction_file)%>% 
#Get interactions with a 'strong connection' in cd4 or cd8 t cells.
  filter(StrongConn_CD4Naive > 0.0 | StrongConn_CD8Naive > 0.0) %>% 
  select(-StrongConn_NaiveB,-StrongConn_Monocyte,-StrongConn_NaiveK)

#get unmerged interactions for cd4 or cd8 t cells
#unmerged_interactions <- fread('/projects/nknoetze_prj/ocr_prj/data/processed/interactions/DICE/shared/Loop_2.5kb/interactions_tss_within.tsv') %>% 
unmerged_interactions <- fread(interaction_file) %>% 
  filter(qval_CD4Naive < 0.01 | qval_CD8Naive < 0.01) %>%
  select(-qval_CD4Naive,-qval_CD8Naive,-qval_Monocyte,-qval_NaiveB,-qval_NaiveK)

### ----------------------------- ###
###        Find genes with        ###
###        merged interactions    ###
### ----------------------------- ###
merged_genes <- merged_interactions %>% 
  filter(StrongConn_CD4Naive > 0.5 | StrongConn_CD8Naive > 0.5) %>% 
  select(bin1_gene_id,bin2_gene_id) %>% 
  separate_rows(bin1_gene_id,sep=";") %>% 
  separate_rows(bin2_gene_id,sep=';') %>% 
  pivot_longer(cols=c('bin1_gene_id','bin2_gene_id')) %>% 
  transmute('gene_id'=value) %>% distinct() %>% filter(gene_id!='none')
  
### ----------------------------- ###
###        Find genes without     ###
###        merged interactions    ###
### ----------------------------- ###
#find the interactions that are not strong (<=0.5) for both cd8 and cd4 (these are genes without any interactions in merged dataset)
merged_missing_genes <- merged_interactions %>% 
  filter(StrongConn_CD4Naive <= 0.5 , StrongConn_CD8Naive <= 0.5) %>%
  select(bin1_gene_id,bin2_gene_id) %>% 
  separate_rows(bin1_gene_id,sep=";") %>% 
  separate_rows(bin2_gene_id,sep=';') %>% 
  pivot_longer(cols=c('bin1_gene_id','bin2_gene_id')) %>% 
  transmute('gene_id'=value) %>% distinct() %>% filter(gene_id!='none') %>%
  filter(!(gene_id %in% merged_genes$gene_id)) %>% unlist()

### ----------------------------- ###
###   Pull unmerged interactions  ###
###      for missing genes        ###
### ----------------------------- ###
raw_interactions <- list()
for (g in merged_missing_genes){
  df <- unmerged_interactions %>% filter(grepl(g,bin2_gene_id)|grepl(g,bin1_gene_id))
  raw_interactions[[g]] <- df
}

raw_interaction_df <- rbindlist(raw_interactions) %>% mutate(interaction_id=paste('unmerged',':',interaction_id,sep=''))

### ----------------------------- ###
###   Create final output file    ###
### ----------------------------- ###
all_interactions <- merged_interactions %>% 
  filter(StrongConn_CD4Naive > 0.5 | StrongConn_CD8Naive > 0.5) %>% 
  select(-StrongConn_CD4Naive,-StrongConn_CD8Naive) %>% 
  mutate(interaction_id=paste('merged',':',interaction_id,sep='')) %>% 
  rbind(raw_interaction_df) 

### ----------------------------- ###
###          Write file           ###
### ----------------------------- ###
fwrite(all_interactions,paste(outdir,'interactions_tss_within_tidy.tsv',sep=''),sep='\t',quote=FALSE)

