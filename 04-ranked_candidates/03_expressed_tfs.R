suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

### ----------------------------- ###
###        SET UP ARGUMENTS       ###
### ----------------------------- ###
# create parser object
parser <- ArgumentParser()
parser$add_argument("-tf", "--tcell_tfs", help="Path list of TFs expressed in T cells")
parser$add_argument("-ff", "--feature_file", help="Path to the feature file obtained from framework")
parser$add_argument("-mt", "--motif_type", help="indicate type of motif (tf, novel, both)")
parser$add_argument("-ot", "--ocr_type", help="indicate type of ocr (promoter or non-promoter")
parser$add_argument("-ntf", "--n_tfs", help="Number of tfbs analysed at a time (2 or 3)",type='integer',required=TRUE)

args <- parser$parse_args()
tcell_tfs <- args$tcell_tfs
feature_file <- args$feature_file
motif_type <- args$motif_type
ocr_type <- args$ocr_type
n_tfs <- args$n_tfs

### ----------------------------- ###
###        READ IN THE DATA       ###
### ----------------------------- ###
print('reading in files')
tcell_expressed_tfs <- fread(tcell_tfs)
#tcell_expressed_tfs <- fread('/projects/nknoetze_prj/references/annotations/tf_fantom/fantom5_tf_hg19_tcellexpressed.tsv')


print('reading in feature file')
feature_results <- fread(feature_file) %>%
#feature_results <- fread('/projects/sbrown_prj/220318_DI/data/processed/tcell_ocr_metrics/OCR0040SB/results_230925.tsv',nThread=48) %>%
  mutate(fc=((target_value+0.000001)/(background_mean+0.000001)),
         feature=str_to_upper(feature)) %>%
  filter(region_set==ocr_type,datatype=='fimo_comotif_geneCount',p_val < 0.05,!(grepl("\\.",feature)),fc >2)


if(motif_type=='tfbs'){
  filtered_comotif_results <- feature_results %>% filter(str_count(feature,"_")==0) %>%  mutate(motif=feature) %>% separate_rows(motif,sep=";|::")
} else if(motif_type=='novel'){
  filtered_comotif_results <- feature_results %>% filter(str_count(feature,"PCA_")==(n_tfs)) 
  } else if(motif_type=='both'){
  filtered_comotif_results <- feature_results %>% filter(!(grepl("PCP_|PCNP_",feature))) %>% filter(str_count(feature, "PCA_")>0) %>% filter(str_count(feature, "PCA_")< n_tfs) %>%  mutate(motif=feature) %>% separate_rows(motif,sep=";|::")
}
print('done reading in and filtering initial file')


### ----------------------------- ###
###        FILTER RESULTS         ###
### ----------------------------- ###
if(motif_type=='both'){
  novel_sites <- filtered_comotif_results %>% 
    filter(grepl("_",motif)) %>% select(-motif) %>% unique()
  tf_sites <- filtered_comotif_results %>% 
    #only look at expression for tfs
    filter(!(grepl("_",motif))) %>% 
    group_by(feature,region_set) %>%
    #find the number of individual tfs for a given co-motif 
    # this can be more than the number of individual motifs for a given co-motif due to dimers
    mutate(n_motif=n_distinct(motif)) %>% 
    #remove motifs who's expression does not meet the threshold
    filter(motif %in% tcell_expressed_tfs$gene_name) %>% 
    #now check if the number of remaining motifs == the number of motifs for the co-motif
    #if so, keep it, since all the tfs need to be expressed!
    group_by(feature,region_set) %>% 
    filter(n_distinct(motif)==n_motif) %>% ungroup() %>% select(-n_motif,-motif) %>% unique() 
  expressed_comotifs <- rbind(novel_sites,tf_sites)
}

if(motif_type=='tfbs'){
  expressed_comotifs <- filtered_comotif_results %>% 
    group_by(feature,region_set) %>%
    #find the number of individual tfs for a given co-motif 
    # this can be more than the number of individual motifs for a given co-motif due to dimers
    mutate(n_motif=n_distinct(motif)) %>% 
    #remove motifs who's expression does not meet the threshold
    filter(motif %in% tcell_expressed_tfs$gene_name) %>% 
    #now check if the number of remaining motifs == the number of motifs for the co-motif
    #if so, keep it, since all the tfs need to be expressed!
    group_by(feature,region_set) %>% 
    filter(n_distinct(motif)==n_motif) %>% ungroup() %>% select(-n_motif,-motif) %>% unique() 
} 

if(motif_type=='novel'){
  expressed_comotifs <- filtered_comotif_results
}
print('done finding expressed TFs')

### ----------------------------- ###
###        WRITE FILES            ###
### ----------------------------- ###
outfile <- paste(gsub('.txt|.tsv','',feature_file),'.',ocr_type,'.',motif_type,'.filtered.tsv',sep='')
fwrite(expressed_comotifs,outfile,sep='\t',quote=FALSE)