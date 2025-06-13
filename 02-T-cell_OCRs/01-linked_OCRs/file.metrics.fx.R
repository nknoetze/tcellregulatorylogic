suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

##########################################
##             OCR DIST                 ##
##########################################
get_ocrs_dists <- function(linked_ocrs,gencode_gene_file){
  #remove gene group from gencode file
  gencode_tidy <- gencode_gene_file %>% distinct(gene_id,gene_name)
  ###################
  #  PROMOTER OCRs  #
  ###################
  #Get the promoter OCRs. If the OCRs overlaps the TSS of the query gene
  #Label as promoter, otherwise linked promoter
  promoter_ocrs <- linked_ocrs %>% rename('query_gene'=gene_id) %>% filter(ocr_type=="promoter") %>%
    mutate(promoter_type=ifelse(promoter_id !=query_gene,'Linked Promoter','Query Gene Promoter')) %>% select(-transcript_id) %>%
    unique()
  promoter_ocr_dist <- promoter_ocrs %>%
    #focus on the ocrs that overlap the tss of the query gene
    filter(promoter_type=='Query Gene Promoter') %>% unique() %>% 
    inner_join(tss,by=c('promoter_id')) %>% 
    #calculate the distance of the tss to the closest end of the ocr.
    #note if it is up or downstream
    mutate(dist_start=tss-start,dist_end=tss-end,
           dist_tss=ifelse(abs(dist_start)<abs(dist_end),dist_start,dist_end))%>%
    #for a given ocr and gene, only keep the nearest tss for the target gene. 
    group_by(ocr_id,query_gene) %>% 
    mutate(min_dist=min(abs(dist_tss))) %>%
    filter(abs(dist_tss)==min_dist) %>%
    select(-dist_start,-dist_end,-min_dist,-tss)%>% unique() %>%
    rename('query_name'=gene_name)
  
  distal_promoter_ocrs <- promoter_ocrs %>% 
    #focus on the ocrs that overlap a tss, but for an interacting gene
    filter(promoter_type=='Linked Promoter') 
  
  # we want to measure the distance of the distal promoter (promoter_id) to the query gene (query_gene)
  #so, get the TSS for the QUERY GENE
  target_promoter_tss <- distal_promoter_ocrs %>% 
    transmute('promoter_id'=query_gene) %>% 
    inner_join(tss) %>% transmute('query_gene'=promoter_id,tss) %>% unique()
  
  
  distal_promoter_ocr_dist <- distal_promoter_ocrs %>% rename('query_name'=gene_name) %>%
    inner_join(target_promoter_tss,by='query_gene') %>% 
    mutate(dist_start=tss-start,dist_end=tss-end,
           dist_tss=ifelse(abs(dist_start)<abs(dist_end),dist_start,dist_end)) %>%
    group_by(ocr_id,query_gene) %>% 
    mutate(min_dist=min(abs(dist_tss))) %>%
    filter(abs(dist_tss)==min_dist) %>% 
    select(-tss,-dist_start,-dist_end,-min_dist) %>% unique() %>% 
    #add the name of the distal promoter gene name
    left_join(gencode_tidy,by=c('promoter_id'='gene_id')) %>% 
    rename('promoter_name'=gene_name)
  
  promoter_ocr_dists <- full_join(distal_promoter_ocr_dist,promoter_ocr_dist) %>% unique()
  
  ######################
  #  NON-PROMOTER OCRs #
  ######################
  #get the non-promoter distal ocrs
  distal_ocrs <- linked_ocrs %>% filter(ocr_type=="non-promoter") %>% select(-transcript_id) %>% unique() %>%
    rename('query_gene'=gene_id)
  
  #remove transcript id, since different transcripts can have the same TSS for a gene
  distal_ocr_cognate_promoter_tss <- linked_ocrs %>% 
    filter(gene_id %in% distal_ocrs$query_gene,transcript_id !='none') %>% 
    transmute('query_gene'=gene_id,promoter_id) %>%
    left_join(tss) %>% rename('target_promoter_id'=promoter_id) %>% 
    select(-promoter_name) %>% unique()
  
  distal_ocr_dist <- distal_ocrs %>% inner_join(distal_ocr_cognate_promoter_tss,by='query_gene') %>% 
    #make sure the query gene and the target promoter id are the same
    #this will make sure we calculate the distance of the OCR to the target gene promoter
    filter(target_promoter_id==query_gene) %>% 
    mutate(dist_start=tss-start,dist_end=tss-end,
           dist_tss=ifelse(abs(dist_start)<abs(dist_end),dist_start,dist_end)) %>% 
    group_by(ocr_id,query_gene,target_promoter_id) %>%
    mutate(min_dist=min(abs(dist_tss))) %>% 
    filter(abs(dist_tss)==min_dist) %>% 
    ungroup() %>%
    select(-tss,-dist_start,-dist_end,-min_dist,-target_promoter_id) %>% unique() %>%
    mutate(promoter_name='none',promoter_type='none') %>%
    rename('query_name'=gene_name)
  ocr_dist <- full_join(distal_ocr_dist,promoter_ocr_dists) 
  return(ocr_dist)
}