library(seqinr)
library(argparse)
library(data.table)
library(tidyverse)
library(tidytext)

###--------------------------------------------------------------------------------------------------------------------###
# SETUP ARGUMENTS
###--------------------------------------------------------------------------------------------------------------------###
# create parser object
parser <- ArgumentParser()
parser$add_argument("-ranked_features", "--ranked_features", help="Path to the file of ranked features")
parser$add_argument("-motif_archetypes", "--motif_archetypes", help="Path to the motif archetypes metadata file")
parser$add_argument("-n_novel", "--n_novel", help="Number of candidate novel motifs to select",type='integer')
parser$add_argument("-n_tfbs", "--n_tfbs", help="Number of candidate tfbs to select",type="integer")
parser$add_argument("-outdir", "--outdir", help="output directory for writing the fasta file and feature table")

args <- parser$parse_args()
ranked_features <- args$ranked_features
motif_archetypes <- args$motif_archetypes
n_novel <- args$n_novel
n_tfbs <- args$n_tfbs
outdir <- args$outdir

###--------------------------------------------------------------------------------------------------------------------###
# READ IN DATA
###--------------------------------------------------------------------------------------------------------------------###
print('reading in features file and archetypes')
#all_pairs_fulloverlap <- fread('/projects/sbrown_prj/220318_DI/data/processed/tcell_ocr_metrics/ZoE_motifs_10k_coenrich/results_231010.all.filtered.full.nonoverlapping.allsites.ranked.tsv',nThread=84) %>% 
all_pairs_fulloverlap <- fread(ranked_features,nThread=84) %>% 
  select(feature, motif, target_value,background_mean,rank, ocr_id,tfbs_pos,strand,matched_sequence,score,motif_type,motif_id,p_val) %>% unique()
filtered_pairs <- all_pairs_fulloverlap %>% 
  #Ensure motif is part of a tfbs:novel pair
  group_by(motif) %>% filter(n_distinct(motif_type)>1) %>% ungroup() %>% 
  #ensure feature is comprised of two distinct motifs
  group_by(feature) %>% filter(n_distinct(motif)==2) %>% ungroup()

#https://resources.altius.org/~jvierstra/projects/motif-clustering-v2.0beta/
#motif_archetypes <- fread('/projects/nknoetze_prj/ocr_prj/data/raw/motif_archetypes/metadata.tsv') 
motif_archetypes <- fread(motif_archetypes) 


###--------------------------------------------------------------------------------------------------------------------###
# FOR 18 CANDIDATE SEQUENCES, HALF TFBS AND HALF NOVEL.
###--------------------------------------------------------------------------------------------------------------------###
print('assigning motifs to archetypes')
filtered_pairs_archetypes <- filtered_pairs %>% rename('motif_class'=motif_type) %>% 
  # Sept 10 2024, remove novel:novel features since we dont use them. speed things up
  filter(motif_class!='novel') %>% 
  #get the cluster ID for each motif_id/binding site prediction
  #if the motif does not belong to a cluster, just use the motif name as the id.
  left_join(motif_archetypes) %>% select(-tf_name,-family_name,-motif_type,-PMID,-source_id) %>%
  mutate(cluster=ifelse(is.na(cluster),motif,cluster)) %>%
  #for each cluster, create a new motif name that combines ALL the motifs belonging to the cluster together.
  group_by(cluster) %>% mutate(cluster_motifs = paste(unique(motif), collapse = "|")) %>% ungroup() %>%
  #for each motif_class, sort by rank and target value. assign a unique id to each feature
  group_by(motif_class) %>% arrange(rank,desc(target_value),background_mean) %>%  mutate(id = match(feature, unique(feature))) %>% ungroup()

print('selecting candidate tfbs')
selected_tfbs_clusters <- filtered_pairs_archetypes %>% 
  #filter dataframe for the tfbs:tfbs pairs
  filter(motif_class=='tfbs') %>% 
  #obtain the first X ranked motifs, as designated by the NEW_MOTIF id
  distinct(cluster_motifs, .keep_all = TRUE) %>% head(n_tfbs) %>% 
  distinct(feature,motif,cluster,cluster_motifs,id) 

selected_tfbs_motif <- filtered_pairs_archetypes %>% 
  #filter dataframe for the tfbs:tfbs pairs
  filter(motif_class=='tfbs') %>% 
  #obtain the first X ranked motifs, as designated by the MOTIF name
  distinct(motif, .keep_all = TRUE) %>% head(n_tfbs) %>% 
  distinct(feature,motif,cluster,cluster_motifs,id) 

#Contains the ids for the features selected above. 
#If we only select X cluster_motifs, we can run into the case where  we have less than X MOTIFs since a MOTIF can belong to multiple cluster_motifs
#Selecting X cluster_motifs, and @ MOTIFs accounts for this! if we hit X cluster_motifs, we keep going until we hit another unique MOTIF
selected_tfbs_cluster_motifs <- rbind(selected_tfbs_motif,selected_tfbs_clusters) %>% 
  unique()

#get all features that are below the max id above. the max id is the feature that results in X unique motif clusters or Motifs.
selected_tfbs_clusters_features <- filtered_pairs_archetypes %>% 
  #filter(motif_class=='tfbs',id <= max(selected_tfbs_cluster_motifs$id))  %>%
  filter(motif_class=='tfbs') %>%
  filter(feature %in% selected_tfbs_cluster_motifs$feature | id <= max(selected_tfbs_cluster_motifs$id))  %>%
  distinct(feature,motif,cluster,cluster_motifs,id) 
  
print('getting novel motif candidates')
#Get all novel:tfbs co-enriched pairs 
# filter the co-enriched pairs. Only keep co-enriched pairs where the tfbs is one of the motifs in the 9 clusters generated above
tfbs_novel_pairs_features <- filtered_pairs_archetypes %>%  
  filter(motif_class=='both',motif %in% selected_tfbs_clusters_features$motif) %>%
  distinct(feature,cluster_motifs)

#Get the binding sites for the motifs/tfs for the features found above.
selected_novel_clusters<- filtered_pairs_archetypes %>% filter(motif_class=='both',feature %in% tfbs_novel_pairs_features$feature) %>%  
  #Since we want to select the top 9 unique NOVEL motifs, filter the results to only keep binding sites for NOVEL motifs
  filter(grepl("PCA_AC",cluster_motifs)) %>% 
  #obtain the first 9 ranked novel motifs, as designated by the NEW_MOTIF id
  distinct(cluster_motifs, .keep_all = TRUE) %>% head(n_novel) %>%
  distinct(feature,motif,cluster,cluster_motifs,id) 

#get all features that are below the max id above. the max id is the feature that results in X unique motif clusters.
#also make sure we are only looking at features where the motif is one of the motifs in the tfbs selected above
#THIS WILL NOT INCLUDE MOTIFS THAT ARE PART OF THE CLUSTERS CHOSEN ABOVE UNLESS THEY ARE FROM THE SAME MOTIF.
selected_novel_clusters_features <- filtered_pairs_archetypes %>% 
  filter(motif_class=='both',feature %in% tfbs_novel_pairs_features$feature,id <= max(selected_novel_clusters$id))
print('sequence selection')
# filter the results to obtain the features for the candidate novel motifs and tfbs 
final_candidates <- filtered_pairs_archetypes %>% 
  #get the selected features
  filter(feature %in% selected_tfbs_clusters_features$feature | feature %in% selected_novel_clusters_features$feature) %>%
  distinct(feature, motif,motif_class,cluster,cluster_motifs,matched_sequence,tfbs_pos,score) %>% 
  #create a motif name using only the motifs that are involved in selecting the candidate sequence for each cluster motif.
  #this will tell us what motifs the selected sequence could represent
  group_by(cluster_motifs) %>% mutate(motifs=paste(unique(motif),collapse='|')) %>% ungroup() %>%
  # for each motif, calculate the number of times a sequence is seen in the OCRs for our T-cell genes
  # Make sure you look at DISTINCT positions.
  #only use the sequences from the TFs that are part of the co-enriched features selected above.
  group_by(cluster_motifs,matched_sequence) %>% mutate(n=n_distinct(tfbs_pos)) %>% ungroup() %>% 
  #select all sequences that appear the most 
  #if there is a tie, select the sequence with the highest fimo score
  #do it based on cluster_motifs because it can represent multiple motifs selected at the beginning..
  group_by(cluster_motifs) %>% filter(n==max(n)) %>% filter(score==max(score)) %>% ungroup() %>% 
  distinct(motif,cluster_motifs, motifs, cluster, matched_sequence,score,n)
  
if(nrow(final_candidates)!=n_tfbs*2){
  final_candidates <- final_candidates %>% group_by(motifs) %>% filter(n==max(n)) %>% filter(score==max(score)) %>% ungroup()
  if(nrow(final_candidates)!=n_tfbs*2){
    final_candidates <- final_candidates %>% group_by(motif) %>% filter(n==max(n)) %>% filter(score==max(score)) %>% ungroup()
  }
}

#keeps all the information for qualified features (ranks, positions, target genes etc)
final_candidates_detailed <- filtered_pairs_archetypes %>% filter(feature %in% selected_tfbs_clusters_features$feature | feature %in% selected_novel_clusters_features$feature) %>%
  #create a motif name using only the motifs that are involved in selecting the candidate sequence for each cluster motif
  group_by(cluster_motifs) %>% mutate(motifs=paste(unique(motif),collapse='|')) %>% ungroup() %>%
  #Indicate if the feature was one of the 18 features initially selected for candidates.
  mutate(selected_feature=ifelse(feature %in% selected_tfbs_clusters$feature | feature %in% selected_novel_clusters$feature,'yes','no' ))

print('creating outputs')
#write candidate sequences to output file.
fasta_output <- paste(outdir,paste('candidate_seqs.',n_novel,'novel_',n_tfbs,'tfbs.fa', sep=""),sep='/')
candidate_table <- paste(outdir,paste('candidates.',n_novel,'novel_',n_tfbs,'tfbs.tsv', sep=""),sep='/')
candidates_detailed <- paste(outdir,paste('candidates.',n_novel,'novel_',n_tfbs,'tfbs.detailed.tsv', sep=""),sep='/')
fwrite(final_candidates,candidate_table,sep='\t',quote=FALSE)
fwrite(final_candidates_detailed,candidates_detailed,sep='\t',quote=FALSE)
write.fasta(sequences =as.list(final_candidates$matched_sequence), names = final_candidates$motifs, file.out = fasta_output)


