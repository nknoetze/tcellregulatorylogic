suppressMessages(library(argparse))
suppressMessages(library(doSNOW))
suppressMessages(library(foreach))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))

### ----------------------------- ###
###        SET UP ARGUMENTS       ###
### ----------------------------- ###
# create parser object
parser <- ArgumentParser()
parser$add_argument("-ff", "--feature_file", help="Path to the feature file obtained from framework",required=TRUE)
parser$add_argument("-tp", "--tfbs_pos", help="Path to the file containing tfbs/motif predictions",required=TRUE)
parser$add_argument("-ntf", "--n_tfs", help="Number of tfbs analysed at a time (2 or 3)",type='integer',required=TRUE)
parser$add_argument("-rs", "--ocr_type", help="ocr type analysed (promoter or non-promoter)",required=TRUE)
parser$add_argument("-nc", "--n_cores", help="number of cores to use",required=TRUE,type='integer')

args <- parser$parse_args()
feature_file <- args$feature_file
tfbs_pos <- args$tfbs_pos
n_tfs <- args$n_tfs
ocr_type <- args$ocr_type
n_cores <- args$n_cores

### ----------------------------- ###
###          READ IN FILES        ###
### ----------------------------- ###
print('reading feature file')
enriched_features <- fread(feature_file) %>% 
  select(feature,metric,region_set,target_value,background_mean,p_val,fc) %>%
  mutate(motif=feature) %>% 
  separate_rows(motif,sep=';') 

print('reading site predictions')
tfbs_locations <- fread(tfbs_pos) %>% 
  filter(score> 11.16) %>%
  rename('motif'=feature) %>%
  separate(ocr_id, into=c('chrom','ocr_start','ocr_end'),sep=':|-',remove=FALSE) %>%
  mutate(site_start=as.numeric(ocr_start)+start,site_end=as.numeric(ocr_start)+end,
         motif_pos_id=paste(motif,';',chrom,':',site_start,'-',site_end,sep='')) %>% select(-start,-end,-ocr_start,-ocr_end,-q_val) %>% unique()

print('filtering site predictions')
filtered_tfbs_locations <- tfbs_locations %>% filter(motif %in% enriched_features$motif)

enriched_features_locations <- as.data.table(enriched_features %>% 
  left_join(filtered_tfbs_locations))
print('finished reading files')

print('Splitting up dataframe')
start <- Sys.time()
# Split the data into overlapping and non-overlapping data frames
# data.table faster than tidyverse operations
enriched_features_locations[, n_motif := length(unique(motif)), by = .(feature, ocr_id)]
#Get df of features that have motifs in the same ocr to filter over.
overlapping_df <- enriched_features_locations[n_motif != 1, ]
#Keep the rest for later. These will be non-overlapping sites.
nonoverlapping_df <- enriched_features_locations[n_motif == 1, ]
print(Sys.time() - start )

### ----------------------------- ###
###         SET UP CLUSTER        ###
### ----------------------------- ###
#SET UP CLUSTER
cl <- makeCluster(n_cores)
registerDoSNOW(cl)
#Get number of features to iterate through
unique_motifs <- unique(overlapping_df$feature)

#set up progress bar
pb <- txtProgressBar(max = length(unique_motifs), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
start <- Sys.time()

### ----------------------------- ###
###       PERFORM FUNCTION        ###
### ----------------------------- ###
#finally, only keep features where there is a binding site for EACH motif in the feature.
#Use the dataframe that has features with motifs on the same ocr as the starting point
print('getting non overlapping sites')
nonoverlapping_sites <- foreach(i=1:length(unique_motifs),.options.snow=opts,.combine='rbind',.packages=c('data.table','GenomicRanges','tidyverse')) %dopar% {
  coenriched_motif=unique_motifs[i]
  filtered_df <- overlapping_df %>% filter(feature==coenriched_motif)
  #Split the motif pairs and get their locations in the ocrs
  motif1=stringr::str_split(coenriched_motif,';')[[1]][1]
  motif2=stringr::str_split(coenriched_motif,';')[[1]][2]
  
  #get motif locations
  motif1_loc <- filtered_df %>% filter(motif==motif1 & feature==coenriched_motif) %>% distinct(feature,motif,chrom,site_start,site_end) %>% makeGRangesFromDataFrame(ignore.strand = TRUE,start.field = 'site_start',end.field = 'site_end')
  motif2_loc <- filtered_df %>% filter(motif==motif2 & feature==coenriched_motif) %>% distinct(feature,motif,chrom,site_start,site_end) %>% makeGRangesFromDataFrame(ignore.strand = TRUE,start.field = 'site_start',end.field = 'site_end')
  
  #find overlapping sites
  overlapping_12 <- subsetByOverlaps(motif1_loc, motif2_loc, invert = FALSE,type='within') %>% as.data.frame()
  overlapping_21 <- subsetByOverlaps(motif2_loc, motif1_loc, invert = FALSE,type='within') %>% as.data.frame() 
  overlapping_sites <- rbind(overlapping_12,overlapping_21)
  
  overlapping_sites_df <- overlapping_sites %>%
    rename('chrom'=seqnames,'site_start'=start,'site_end'=end) %>%
    select(-width,-strand) %>% left_join(filtered_df,by=c("chrom", "site_start", "site_end"))
  
  #Remove any overlapping sites
  non_overlapping <- filtered_df %>% filter(!(motif_pos_id %in% overlapping_sites_df$motif_pos_id))
}

### ------------------------------- ###
### COMBINE NON-OVERLAPPING RESULTS ###
### ------------------------------- ###
nonoverlapping_sites <- nonoverlapping_sites %>%
  #add back the non-overlapping sites from the original df
  rbind(nonoverlapping_df) %>% 
  group_by(feature) %>% filter(n_distinct(motif)==n_tfs) %>% 
  mutate(tfbs_pos=paste(chrom,':',site_start,'-',site_end,sep='')) %>% 
 select(-motif_pos_id,-site_start,-site_end,-n_motif) %>% unique() 
stopCluster(cl)
print(Sys.time() - start )

outfile <- paste(gsub('.txt|.tsv','',feature_file),'.full.nonoverlapping.tsv',sep='')
print('writing output file')
print(outfile)
fwrite(nonoverlapping_sites,outfile,sep='\t',quote=FALSE)
