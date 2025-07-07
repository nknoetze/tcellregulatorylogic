library(data.table)
library(pheatmap)
library(knitr)
library(tidytext)
#need to use this version of cmake in order to insall ggpubr/nloptr
#old_path <- Sys.getenv("PATH")
#Sys.setenv(PATH = paste("/gsc/software/linux-x86_64-centos7/cmake-3.22.2/bin", old_path, sep = ":")) 
#BiocManager::install('nloptr')

## NEED TO INSTALL OLD VERSION TO BE COMPATIBLE WITH CURRENT R VERSION..
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/ggpubr/ggpubr_0.2.tar.gz" 
#install.packages(packageurl, repos=NULL, type="source",dependencies = FALSE)
library(ggpubr)
## NEED TO INSTALL OLD VERSION TO BE COMPATIBLE WITH CURRENT R VERSION..
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.5.1.tar.gz" 
#install.packages(packageurl, repos=NULL, type="source",dependencies = FALSE)
library(rstatix)
library(tidyverse)
library(argpase)
### ----------------------------- ###
###        SET UP ARGUMENTS       ###
### ----------------------------- ###
# create parser object
parser <- ArgumentParser()
parser$add_argument("-tg", "--target_genes", help="Path to the list of target genes")
parser$add_argument("-od", "--outdir", help="Path to the framework output directory")
parser$add_argument("-res", "--framework_results", help="Path to the framework results file")

args <- parser$parse_args()
target_genes <- args$target_genes
outdir <- args$outdir
framework_results <- args$framework_results

### ----------------------------- ###
###     PROCESSING FUNCTION       ###
### ----------------------------- ###
process_data <- function(results_file,expression_file,tcell_expressed_tf_file,feature_location_dir){
  print('reading in file')
  results_files <- list.files(feature_location_dir,pattern='ranked',full.names = TRUE)
  
  #all results have been filtered already in the rank_tfs.smk pipeline. no need to filter for expressed tfs
  results_df <- map_dfr(results_files,fread,nThread=48)  %>% 
    distinct(feature,region_set,target_value,background_mean, p_val,fc,motif_type,rank) %>%
    mutate(feature_region=case_when(region_set=='promoter' ~ paste('p',feature,sep=''),
                                    region_set=='non-promoter' ~ paste('np',feature,sep=''),
                                    TRUE ~ paste('a',feature,sep=''))) 
  
  print('creating expression file')
  expressed_comotifs_df <- results_df %>% filter(motif_type !='novel') %>% 
    distinct(motif_type, feature, region_set) %>% 
    #separate the feature to get the individual motif/tf names
    separate_rows(feature,sep=";|::") %>%
    #focus just on tfs for expression
    filter(!(grepl("_",feature))) %>% unique() %>%
    #add the expression values
    inner_join(expression_file,by=c('feature'='gene_name')) %>% select(-gene_id) %>%
    pivot_longer(cols=c(-feature,-region_set,-motif_type)) %>% 
    #only keep TFs that are expressed in t cells
    filter(feature %in% tcell_expressed_tf_file$gene_name) %>% 
    #Only keep entries for cell types where expression meeets threshold
    filter(value>tcell_median) %>% 
    #cleanup names
    rename('cell_type'=name) %>% 
    mutate(cell_group=ifelse(cell_type %in% tcells,'T-cells','Non-T cells'),
           region_set=str_to_title(region_set))
  
  return(list(results=results_df, expressed_comotifs = expressed_comotifs_df)) 
}

### ----------------------------- ###
###     Read in files             ###
### ----------------------------- ###
gene_expression <- fread('/projects/nknoetze_prj/ocr_prj/data/median_expr_nobatchcorr.tsv')
tcells <- c('cd4_naive','cd8_naive','th17','th1-17','th1','tfh','th2','treg_memory','treg_naive','cd4_n_stim','cd8_n_stim')
tcell_median=237.04
tcell_expressed_tfs <- fread('/projects/nknoetze_prj/references/annotations/tf_fantom/fantom5_tf_hg19_tcellexpressed.tsv')

### ----------------------------- ###
###     PROCESS DATA              ###
### ----------------------------- ###
processed_zoe_files <- process_data(framework_results,expression_file = gene_expression,tcell_expressed_tf_file = tcell_expressed_tfs,feature_location_dir = outdir)
prefix <- list.files(outdir,pattern='ranked',full.names = TRUE) %>% basename() %>% str_extract("^[^.]+") %>% unique()
output_dir <- paste(outdir,'/',prefix,'.processedcomotif.zoe.RDS',sep='')

saveRDS(processed_zoe_files,output_dir)
