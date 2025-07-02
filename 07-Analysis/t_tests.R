library(data.table)
library(tidyverse)


process_counts <- function(data_file,type){
 if(type=='single'){
    oligos <- fread(data_file) %>% select(oligo,background,cell, m1,m2,m3) %>% filter(background!="strong") %>%
        pivot_longer(cols=c(-oligo,-background,-cell),names_to='motif_group',values_to='motif') %>% 
      group_by(oligo,background,cell) %>% filter(motif!="spacer") %>%
        filter(n_distinct(motif)==1) %>% ungroup()
  }
  else if(type=="pair"){
    oligos <- fread(data_file) %>% select(oligo,background,cell, m1,m2,m3) %>% filter(background!="strong") %>%
      pivot_longer(cols=c(-oligo,-background,-cell),names_to='motif_group',values_to='motif') %>% 
      group_by(oligo,background,cell) %>% filter(motif!="spacer") %>%
      filter(n_distinct(motif)==2) %>% ungroup()
  }
  final <- fread(data_file) %>% filter(background=="synthetic", oligo %in% oligos$oligo| oligo=="spacer-spacer-spacer") %>%
    rowwise() %>% mutate(num_na = sum(is.na(c(RNA1_log2fc, RNA2_log2fc, RNA3_log2fc)))) %>% filter(num_na < 2) %>% 
    select(oligo, background, cell,RNA1_log2fc_adj, RNA2_log2fc_adj, RNA3_log2fc_adj) %>% 
    pivot_longer(cols=c(-oligo,-background,-cell),names_to='sample_id',values_to='log2fc') %>%
    filter(!is.na(log2fc)) 
  return(final)
}

two_sided_t_test <- function(df,background_seqs_df) {
  # Split the data into two groups based on motif_type
  motif_values <- df$log2fc
  background_values <- background_seqs_df$log2fc
  t_test_result <- t.test(motif_values, background_values, alternative = "two.sided")
  results_df <- data.frame(
    statistic = c(t_test_result$statistic),
    p_value = c(t_test_result$p.value),
    conf_int = c(paste(t_test_result$conf.int, collapse = ", ")),
    method = c(t_test_result$method))
  return(results_df)
}

test_data <- function(dataframe){
  t_test_results_list <- list()
  orientations <- c('fwd','rev')
  cells <- c("Jurkat","K562")
  for(c in cells){
    for(o in orientations){
      for(m in 1:length(motifs)){
        motif <- motifs[m]
        motif_orientation <- paste(motif,o,sep='')
        motif_ignore <- paste0(motif, ifelse(grepl("rev", motif_orientation), "fwd", "rev"))
        
        with_mo <- dataframe %>% filter(str_detect(oligo,motif_orientation),!(str_detect(oligo,motif_ignore)),cell==c)
        with_mo_mean <- mean(with_mo$log2fc,na.rm = TRUE)
        with_mo_se <- sd(with_mo$log2fc, na.rm = TRUE) / sqrt(length(with_mo$log2fc))
        
        without_mo <- background_df %>% filter(!(str_detect(oligo,motif)),cell==c)
        without_mo_mean <- mean(without_mo$log2fc,na.rm = TRUE)
        without_mo_se <- sd(without_mo$log2fc, na.rm = TRUE) / sqrt(length(without_mo$log2fc))
        
        cohens <- as.data.frame(effectsize::cohens_d(with_mo$log2fc,without_mo$log2fc))
        t_test <- two_sided_t_test(with_mo,without_mo) %>% mutate(motif_orientation=motif_orientation,cell=c)
        df <- rbind(with_mo %>% mutate(group="Test Set"),without_mo %>% mutate(group="Background Set")) %>%
          mutate(cell=c,motif=motif,motif_orientation=motif_orientation,motif_cohens=cohens$Cohens_d,mean_log2fc_with_motif=with_mo_mean,mean_log2fc_without_motif=without_mo_mean,
                 se_without_motif=without_mo_se,se_with_motif=with_mo_se,orientation=o) %>% full_join(t_test,by = c("motif_orientation", "cell"))
        t_test_results_list[[paste(motif,c,o,sep='')]] <-  df
      }
    }
  }
  t_test_results <- rbindlist(t_test_results_list)
  return(t_test_results)
}

test_pairs <- function(dataframe){
  t_test_results_list <- list()
  orientations <- c('fwd-fwd','rev-rev','fwd-rev','rev-fwd')
  cells <- c("K562","Jurkat")
  for(c in cells){
    for(i in 1:(length(motifs) - 1)){  # Ensure not to exceed the last index
      motif1 <- motifs[i]
      for(j in (i+1):length(motifs)){
        motif2 <- motifs[j]
        
        for(o in orientations){
          orientation_split <- str_split_fixed(o,'-',2)
          motif1_orientation <- paste0(motif1, orientation_split[1])
          motif2_orientation <- paste0(motif2, orientation_split[2])
          motif1_ignore <- paste0(motif1, ifelse(grepl("rev", motif1_orientation), "fwd", "rev"))
          motif2_ignore <- paste0(motif2, ifelse(grepl("rev", motif2_orientation), "fwd", "rev"))
          
          with_mo <- dataframe %>% filter(str_detect(oligo,motif1),str_detect(oligo,motif2),!str_detect(oligo,motif1_ignore),!str_detect(oligo,motif2_ignore),cell==c) 
          with_mo_mean <- mean(with_mo$log2fc,na.rm = TRUE)
          with_mo_se <- sd(with_mo$log2fc, na.rm = TRUE) / sqrt(length(with_mo$log2fc))
          
          without_mo <- dataframe %>% filter(!str_detect(oligo,motif1),!str_detect(oligo,motif2),cell==c)
          without_mo_mean <- mean(without_mo$log2fc,na.rm = TRUE)
          without_mo_se <- sd(without_mo$log2fc, na.rm = TRUE) / sqrt(length(without_mo$log2fc))
          
          # Cohen's d for combined motifs (motif 1 and motif 2 vs no motif)
          pair_cohens <- as.data.frame(effectsize::cohens_d(with_mo$log2fc, without_mo$log2fc))
          t_test <- two_sided_t_test(with_mo,without_mo) %>% mutate(motif1=motif1_orientation,motif2=motif2_orientation,cell=c)
          df <-t_test %>% mutate(motif1=motif1_orientation,motif2=motif2_orientation,motif_pair=paste(motif1_orientation,motif2_orientation,sep='-'),
                                 pair_cohens=pair_cohens$Cohens_d,mean_log2fc_with_pair=with_mo_mean,mean_log2fc_without_pair=without_mo_mean,
                                 se_without_pair=without_mo_se,se_with_pair=with_mo_se,orientation=o)
          t_test_results_list[[paste(motif_pair=paste(motif1_orientation,motif2_orientation,sep='-'),c,sep='')]] <-  df
        }
      }
    }
  }
  t_test_results <- rbindlist(t_test_results_list)
  return(t_test_results)
}

annotate_synergistic_interactions <- function(motif_pair_list,cohen_single_df,cohen_pair_df){
  annotated_df <- list()
  for(c in c("K562","Jurkat")){
    for(p in 1:(nrow(motif_pair_list))){
      pair <- motif_pair_list$motif_pair[p]
      motif1 <- str_split(pair,'-')[[1]][1]
      motif2 <- str_split(pair,'-')[[1]][2]
      
      m1_mean <- cohen_single_df %>% filter(cell==c,motif_orientation==motif1) %>% transmute(motif_orientation,'motif1_log2fc'=mean_log2fc_with_motif,'cohens'=motif_cohens) %>% unique()
      m2_mean <- cohen_single_df %>% filter(cell==c,motif_orientation==motif2) %>% transmute(motif_orientation,'motif2_log2fc'=mean_log2fc_with_motif,'cohens'=motif_cohens) %>% unique()
      m1_m2_mean <- cohen_pair_df %>% filter(cell==c,motif_pair==pair) %>% transmute(motif_pair,'pair_mean_log2fc'=mean_log2fc_with_pair,'cohens'=pair_cohens)
      
      cohens <- rbind(m1_mean %>% select(motif_orientation,cohens),m2_mean %>% select(motif_orientation,cohens),m1_m2_mean %>% transmute("motif_orientation"=motif_pair,cohens)) 
      if(abs(m1_mean$motif1_log2fc)> abs(m2_mean$motif2_log2fc)){
        max_motif_log2fc <- m1_mean$motif1_log2fc
        nonmax_motif_log2fc <- m2_mean$motif2_log2fc
        max_motif <- motif1
        min_motif <- motif2
      } else{
        max_motif_log2fc <- m2_mean$motif2_log2fc
        nonmax_motif_log2fc <- m1_mean$motif1_log2fc
        max_motif <- motif2
        min_motif <- motif1
      }
      filtered <- m1_m2_mean  %>% select(-cohens) %>% mutate(max_motif_log2fc=max_motif_log2fc,nonmax_motif_log2fc=nonmax_motif_log2fc) %>% 
        group_by(motif_pair) %>% mutate(additive=sum(max_motif_log2fc,nonmax_motif_log2fc),cell=c) %>% ungroup()
      annotated_df[[paste(pair,c,sep='')]] <- filtered
      
    }
    annotated_interactions <- rbindlist(annotated_df)
  }
  return(annotated_interactions)
}


monotypic_counts <- process_counts('/projects/nknoetze_prj/ocr_prj/results/test_run/STARR-seq/candidate.master.count.file.tsv',type='single') 
ditypic_counts <- process_counts('/projects/nknoetze_prj/ocr_prj/results/test_run/STARR-seq/candidate.master.count.file.tsv',type='pair') 

################################################################################################################################################
motifs <- c('BATF::JUN','FOXO1::ELK3', 'PCA_AC1141','PCA_AC0518','PCA_AC2464' ,'PCA_AC0432','PCA_AC0437','PCA_AC2375','PCA_AC2374','PCA_AC1205', 'FOSB::JUN','ZNF274','HSF1','POU6F1','TCF7','GATA3','PCA_AC2802',"PBX2")
monotypic_results <- test_data(monotypic_counts) %>% mutate(signif=ifelse(p_value<(0.05/36)& abs(motif_cohens) >0.5,"Significant","Not Significant"))
ditypic_results <- test_pairs(ditypic_counts) %>% mutate(pair_signif=ifelse(p_value<(0.05/612)& abs(pair_cohens) >0.5,"Significant","Not Significant"))

######## INTERACTIONS
paired_interactions <- ditypic_results %>% 
  pivot_longer(cols=c(motif1,motif2)) %>%
  rename('motif_orientation'=value) %>% distinct(mean_log2fc_with_pair,cell,p_value,motif_pair,pair_cohens,motif_orientation,se_with_pair,pair_signif) %>% 
  left_join(monotypic_results %>% select(-p_value),by=c("motif_orientation","cell"),relationship='many-to-many') %>% 
  distinct(mean_log2fc_with_pair,pair_signif,p_value,cell,motif_pair,pair_cohens,motif_cohens,motif_orientation,se_with_pair,se_with_motif) %>%
  group_by(cell,motif_pair) %>% 
  filter(abs(motif_cohens)==max(abs(motif_cohens))) %>% ungroup() %>% 
  mutate(expected_cohens=motif_cohens,
         expected_se=ifelse(expected_cohens==motif_cohens,se_with_motif,"NA")) %>% 
  filter(expected_se!="NA") %>% 
  mutate(Z=(pair_cohens-expected_cohens)/sqrt(as.numeric(expected_se)^2+se_with_pair^2),
         interaction_p_value=2 * (1 - pnorm(abs(Z)))) %>% select(-expected_se) %>% ungroup() %>%
  mutate(interaction_signif=ifelse(pair_signif=="Significant" & interaction_p_value < 0.05,"Potential Interaction","No Interaction")) 

annotated_synergistic_interactions <- annotate_synergistic_interactions(paired_interactions %>% filter(interaction_signif=="Potential Interaction"),monotypic_results,ditypic_results) %>%
  mutate(synergistic=case_when(pair_mean_log2fc > 0 & additive > 0 & max_motif_log2fc >0 & nonmax_motif_log2fc >0 & pair_mean_log2fc > additive~ "Synergistic",
                               pair_mean_log2fc < 0 & additive < 0 & max_motif_log2fc <0 & nonmax_motif_log2fc <0 & pair_mean_log2fc < additive~ "Synergistic",
                               TRUE ~ "Not Synergistic"))

###WRITE RESULTS
fwrite(monotypic_results,'/projects/nknoetze_prj/ocr_prj/results/test_run/STARR-seq/monotypic.results.tsv',quote=FALSE,sep='\t')
paired_interactions %>% left_join(annotated_synergistic_interactions,by=c("motif_pair","cell")) %>% 
  fwrite('/projects/nknoetze_prj/ocr_prj/results/test_run/STARR-seq/ditypic.results.tsv',sep='\t',quote=FALSE)
