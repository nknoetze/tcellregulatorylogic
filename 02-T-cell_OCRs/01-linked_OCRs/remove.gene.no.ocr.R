suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))


### ----------------------------- ###
###        SET UP ARGUMENTS       ###
### ----------------------------- ###
# create parser object
parser <- ArgumentParser()
parser$add_argument("-gf", "--gene_file", help="Path gene file to filter")
parser$add_argument("-of", "--pulled_ocr_file", help="Path linked ocrs")
parser$add_argument("-o", "--outdir", help="Path to the output directory")

args <- parser$parse_args()
gene_file <- args$gene_file
pulled_ocr_file <- args$pulled_ocr_file
outdir <- args$outdir

### ----------------------------- ###
###        READ IN FILES          ###
### ----------------------------- ###

#gene_list <- fread('/projects/nknoetze_prj/ocr_prj/data/processed/gene_lists/ranked_gene_list.tsv')
gene_list <- fread(gene_file)
genes <- gene_list %>% distinct(gene_id,gene_name)

#pulled_ocrs <- fread('/projects/nknoetze_prj/ocr_prj/results/ocr_metric_files/meuleman_tcell_2sample_1cd8_1cd4_ocrs_merged/Loop_2.5kb/merged_interactions/connectivity8/tss_within/pulled_ocrs_tss_dist.tsv') 
pulled_ocrs <- fread(pulled_ocr_file) 

bin_size=strsplit(pulled_ocr_file,split='/')[[1]][8]
filtered_ocr_data_name=strsplit(pulled_ocr_file,split='/')[[1]][7]
merge_param=strsplit(pulled_ocr_file,split='/')[[1]][10]

exclude_genes <- pulled_ocrs %>% 
  mutate(promoter_type=case_when(ocr_type=='promoter' & promoter_id != query_gene ~ 'linked_promoter', ocr_type=='promoter' & promoter_id==query_gene ~ 'promoter',TRUE~'none')) %>% 
  filter(promoter_type !='linked_promoter') %>% 
  group_by(query_gene) %>% 
  summarise(n_ocrtype=n_distinct(ocr_type)) %>% ungroup() %>% filter(n_ocrtype==1) %>%
  left_join(genes,by=c('query_gene'='gene_id')) %>% filter(!(is.na(gene_name))) %>% transmute(gene_name,'gene_id'=query_gene) 

genes_no_ocrs <- gene_list %>% filter(!(gene_id %in% pulled_ocrs$query_gene)) %>% select(gene_name,gene_id)
genes_to_exclude <- rbind(exclude_genes,genes_no_ocrs)

filtered_genelist <- gene_list %>% filter(!(gene_id %in% genes_to_exclude$gene_id)) %>% arrange(t_o_rank)

outfile <- paste(outdir,sub('.tsv','',basename(gene_file)),'_filtered+',filtered_ocr_data_name,'+',bin_size,'+',merge_param,'.tsv',sep='')
print(outfile)
fwrite(filtered_genelist,outfile,sep='\t',quote=FALSE)
