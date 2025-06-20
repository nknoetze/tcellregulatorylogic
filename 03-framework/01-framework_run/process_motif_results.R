library(data.table)
library(tidyverse)

##### FX##
process_data <- function(results_file, expression_file,tcell_expressed_tf_file){
  tcells <- tcells <- c('cd4_naive','cd8_naive','th17','th1-17','th1','tfh','th2','treg_memory','treg_naive','cd4_n_stim','cd8_n_stim')
  #filter results for fc >2 and p value < 0.05.
  results_df <- fread(results_file, nThread=48) %>% 
    filter(grepl("fimo_geneCount",datatype),p_val < 0.05,!(grepl("\\.",feature))) %>% 
    #not all results contain the background values,making dummy column so this column is removed for the next steps
    mutate(fc=((target_value+0.000001)/(background_mean+0.000001))) %>%
    filter(fc > 2) %>% 
    mutate(feature=str_to_upper(feature),
           motif_type=ifelse(grepl("PCA_AC", feature),"Novel Motif", "Known Motif"),
           feature_region=case_when(region_set=='promoter' ~ paste('p',str_to_upper(feature),sep=''),
                                    region_set=='non-promoter' ~ paste('np',str_to_upper(feature),sep=''),
                                    TRUE ~ paste('a',str_to_upper(feature),sep='')))# %>% select(-background_values)
  
  expressed_motifs_df <- results_df %>% 
    #separate the feature to get the individual motif/tf names
    mutate(tf=feature) %>% 
    separate_rows(tf,sep="::")%>% 
    #determine if feature is a dimer (contains two motifs)
    group_by(feature) %>% mutate(n_tf=n_distinct(tf)) %>% 
    select(tf,n_tf,feature, region_set, fc,feature_region) %>%
    #add the expression values
    inner_join(expression_file,by=c('tf'='gene_name')) %>% select(-gene_id) %>%
    pivot_longer(cols=c(-tf,-fc,-region_set,-n_tf,-feature,-feature_region)) %>% 
    #only keep TFs that are expressed in t cells
    filter(tf %in% tcell_expressed_tf_file$gene_name) %>% 
    #Only keep entries for cell types where expression meeets threshold
    filter(value>237.04) %>% 
    #Only keep entries for where both the tfs are expressed if part of a dimer
    group_by(feature) %>% filter(n_distinct(tf)==n_tf) %>% ungroup() %>% 
    #cleanup names
    rename('cell_type'=name) %>% 
    mutate(cell_group=ifelse(cell_type %in% tcells,'T-cells','Non-T cells'),
           cell_type=gsub("_n_"," ",cell_type),
           cell_type=gsub("_"," ",cell_type),
           cell_type=gsub("_n_"," ",cell_type),
           cell_type=str_to_title(cell_type),
           region_set=str_to_title(region_set))
  
  # list of motifs that are expressed and significantly enriched.
  #if novel motif, just get significant results.
  significant_motifs_df <- results_df %>% 
    filter(feature_region %in% expressed_motifs_df$feature_region | feature_region %in% results_df$feature_region & motif_type=="Novel Motif")
  
  return(significant_motifs_df) 
}

#####################################################################################################################
gene_expression <- fread('/projects/nknoetze_prj/ocr_prj/data/median_expr_nobatchcorr.tsv')
tcell_expressed_tfs <- fread('/projects/nknoetze_prj/references/annotations/tf_fantom/fantom5_tf_hg19_tcellexpressed.tsv')
zoe_results <- '/projects/sbrown_prj/220318_DI/data/processed/tcell_ocr_metrics/ZoE_single+allbutmotifs_10k/results_231025.tsv'
izoe_results <- '/projects/sbrown_prj/220318_DI/data/processed/tcell_ocr_metrics/iZoE_single+allbutmotifs_10k/results_231028.tsv'

processed_zoe_files <- process_data(zoe_results,expression_file=gene_expression,tcell_expressed_tf_file = tcell_expressed_tfs)
processed_izoe_files <- process_data(izoe_results,expression_file = gene_expression,tcell_expressed_tf_file = tcell_expressed_tfs)

fwrite(processed_zoe_files,'/projects/sbrown_prj/220318_DI/data/processed/tcell_ocr_metrics/ZoE_single+allbutmotifs_10k/results_231025.significant.motifs.tsv',sep='\t',quote=FALSE)
fwrite(processed_izoe_files,'/projects/sbrown_prj/220318_DI/data/processed/tcell_ocr_metrics/iZoE_single+allbutmotifs_10k/results_231028.significant.motifs.tsv',sep='\t',quote=FALSE)
