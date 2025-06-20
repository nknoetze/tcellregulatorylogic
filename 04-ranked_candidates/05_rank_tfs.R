suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

### ----------------------------- ###
###        SET UP ARGUMENTS       ###
### ----------------------------- ###
# create parser object
parser <- ArgumentParser()
parser$add_argument("-ff", "--feature_file", help="Path to the feature file to rank")

args <- parser$parse_args()
feature_file <- args$feature_file

### ----------------------------- ###
###          SCORE Fx             ###
### ----------------------------- ###
rank_score <- function(coenriched_motif_df){
  #For co-enriched promoter motifs, get the rank for the target value and the background value
  # we want a high target value and a low background value
  rank_weighted <- coenriched_motif_df %>% group_by(motif_type) %>% 
    mutate(target_value_rank=dense_rank(desc(as.numeric(target_value))),
           background_mean_rank=dense_rank(round(background_mean,digits=0)),
           target_value_rank_adj=(target_value_rank/max(target_value_rank)),
           background_mean_rank_adj=(background_mean_rank/max(background_mean_rank)),
           total_score_adj=(target_value_rank_adj+background_mean_rank_adj)) %>% 
    arrange(total_score_adj) %>% mutate(rank=dense_rank(total_score_adj)) %>% ungroup()
  return(rank_weighted)
}

### ----------------------------- ###
###        READ IN THE DATA       ###
### ----------------------------- ###
print('reading in file')
features <- fread(feature_file,nThread=48)
print('getting motif type')
features <- features %>% mutate(motif_type=case_when(str_count(feature,"_AC")==2 ~ "novel",str_count(feature,"_AC")==1 ~ "both",TRUE ~ "tfbs"))


### ----------------------------- ###
###        RANK THE FEATURES      ###
### ----------------------------- ###
print('ranking features')
ranked_features <-  rank_score(features)
print('done')
### ----------------------------- ###
###    WRITE THE OUTPUT FILE      ###
### ----------------------------- ###
outfile <- paste(gsub('.txt|.tsv','',feature_file),'.ranked.tsv',sep='')
fwrite(ranked_features,outfile,sep='\t',quote=FALSE)
 


