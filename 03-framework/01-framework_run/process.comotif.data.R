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
process_data <- function(results_file,gene_file, expression_file,tcell_expressed_tf_file,feature_location_dir){
  print('reading in file')
  results_files <- list.files(feature_location_dir,pattern='ranked',full.names = TRUE)
  #  names(results_files) <- basename(results_files) %>% sub(".*\\.(\\w+)\\.filtered\\.tsv", "\\1", .)
  
  #all results have been filtered already in the rank_tfs.smk pipeline. no need to filter for expressed tfs
  results_df <- map_dfr(results_files,fread,nThread=48)  %>% 
    distinct(feature,region_set,target_value,background_mean, p_val,fc,motif_type,rank) %>%
    #results_df <- map_dfr(results_files,fread,.id='motif_type',nThread=48)  %>% 
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
  
  #get the percent of t-cell subtypes and non-tcell cell types that express 
  #tfs for a given motif.
  percent_expressed_motifs_df <- expressed_comotifs_df %>% 
    group_by(feature,cell_group) %>% mutate(n_celltype=n_distinct(cell_type)) %>% ungroup() %>% 
    select(-cell_type,-value) %>% unique() %>% 
    pivot_wider(names_from=cell_group,values_from=n_celltype) %>% 
    mutate(`T Cell`=(`T-cells`/11)*100,`Non-T Cell`=(`Non-T cells`/13)*100) %>%
    replace(is.na(.), 0)
  
  #get stats. lots of results, so we can do the t-test for co-motifs
  stat.test_df <- percent_expressed_motifs_df %>% select(-`T-cells`,-`Non-T cells`) %>%
    pivot_longer(cols=c(`T Cell`,`Non-T Cell`), names_to="variable",values_to='value')%>% 
    group_by(region_set,motif_type) %>%
    t_test(value ~ variable, paired=FALSE) %>%
    adjust_pvalue(method='hochberg') %>%
    add_significance("p.adj") %>% ungroup()
  
  return(list(results=results_df, expressed_comotifs = expressed_comotifs_df,percent_expressed_motifs=percent_expressed_motifs_df,stat.test=stat.test_df)) 
}

### ----------------------------- ###
###     Read in files             ###
### ----------------------------- ###
gene_expression <- fread('/projects/nknoetze_prj/ocr_prj/data/median_expr_nobatchcorr.tsv')
types <- c('protein_coding','TR_C_gene','IG_C_gene')
gencode <- fread('/projects/nknoetze_prj/references/annotations/gencode/gencode.v19.transc.type.1based.tsv') %>% 
  filter(gene_type %in% types,transcript_type %in% types) %>% select(gene_id,gene_name) %>% unique()
gene_names <- gencode %>% distinct(gene_id,gene_name)
tf <- fread('/projects/nknoetze_prj/references/annotations/tf_fantom/fantom5_tf_hg19.tsv') %>% select(Symbol,EntrezGene)
tcells <- c('cd4_naive','cd8_naive','th17','th1-17','th1','tfh','th2','treg_memory','treg_naive','cd4_n_stim','cd8_n_stim')
tcell_median=237.04
tcell_expressed_tfs <- fread('/projects/nknoetze_prj/references/annotations/tf_fantom/fantom5_tf_hg19_tcellexpressed.tsv')

#gene lists
target_genes <- fread(target_genes,header = FALSE,col.names = c('gene_id')) %>%
  inner_join(gene_names) 

### ----------------------------- ###
###     PROCESS DATA              ###
### ----------------------------- ###
processed_zoe_files <- process_data(framework_results,gene_file = target_genes,expression_file = gene_expression,tcell_expressed_tf_file = tcell_expressed_tfs,feature_location_dir = outdir)
prefix <- list.files(outdir,pattern='ranked',full.names = TRUE) %>% basename() %>% str_extract("^[^.]+") %>% unique()
output_dir <- paste(outdir,'/',prefix,'.processedcomotif.zoe.RDS',sep='')

saveRDS(processed_zoe_files,output_dir)
