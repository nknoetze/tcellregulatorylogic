library(data.table)
library(tidyverse)
options(scipen=999)

# Read in the median expression for each gene by cell subtype
median_expr <- fread('/projects/nknoetze_prj/ocr_prj/data/median_expr_nobatchcorr.tsv') %>% select(-cd4_n_stim,-cd8_n_stim)
gencode <- fread('/projects/nknoetze_prj/references/annotations/gencode/gencode.v19.transc.type.1based.tsv') %>% 
  filter(gene_type=="protein_coding" | gene_type=="TR_C_gene"| gene_type=="IG_C_gene") %>% 
  filter(transcript_type=="protein_coding" | transcript_type=="TR_C_gene"| transcript_type=="IG_C_gene") %>%
  select(gene_id,gene_name) %>% unique()
#list of tcell columns
tcells <- c('cd4_naive','cd8_naive','th17','th1-17','th1','tfh','th2','treg_memory','treg_naive','gene_id')

# select for t cell samples, calculate the sum of medians (across each cell type) for each gene
tcell_df <- median_expr %>% select(cd4_naive,cd8_naive,th17,`th1-17`,th1,tfh,th2,treg_memory,treg_naive,gene_id) %>% 
  mutate(tsum=rowSums(across(where(is.numeric)))) 

# select the non-t cell samples, calculate the sum of medians (across each cell type) for each gene
nontcell_df <- median_expr %>% select(-cd4_naive,-cd8_naive,-th17,-`th1-17`,-th1,-tfh,-th2,-treg_memory,-treg_naive,-gene_name)  %>% 
  mutate(osum=rowSums(across(where(is.numeric)))) 

#get the median of the sum of median for t cells. use as a filter to remove genes with insufficient expression in t cells
tcell_median <-median(tcell_df$tsum)

#combine the two dataframe, remove genes where their expression is below the median of tsum
#arrage in ascending order by the non-tcell sum (osum), rank.
#arrage in descending order by the tcell sum (tsum), rank.
#get the rank sums (combining the trank and orank)
#get the gene names
#order based on the rowsums in ascending order.
gene_df <- inner_join(tcell_df,nontcell_df) %>%
  filter(tsum > tcell_median) %>% 
  arrange(osum) %>% mutate(orank=1:nrow(.)) %>%
  arrange(desc(tsum)) %>% mutate(trank=1:nrow(.)) %>% 
  mutate(t_o_rank=trank+orank) %>% inner_join(gencode) %>%
  arrange(t_o_rank)  %>%
  mutate(rownumber=1:nrow(.)) %>%
  mutate(gene_group=ifelse(rownumber<41, 'tcell_specific','non_specific')) %>% 
  select(-rownumber)


fwrite(gene_df,'/projects/nknoetze_prj/ocr_prj/data/processed/gene_lists/ranked_gene_list.tsv',quote=FALSE,sep='\t')