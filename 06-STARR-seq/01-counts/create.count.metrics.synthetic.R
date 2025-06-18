library(data.table)
library(tidyverse)
library(preprocessCore)

parse_file <- function(RNA_file_path,DNA_file_path,regex){
  RNA_UMI_count_files <- list.files(RNA_file_path,pattern=regex,full.names = TRUE)
  names(RNA_UMI_count_files) <- basename(RNA_UMI_count_files) %>% sub(regex, "", .)
  raw_RNA_count <- map_dfr(RNA_UMI_count_files,
                           ~fread(.x,sep='\t',col.names=c('candidate','count')), .id = 'sample_id',nThread=48)
  
  DNA_UMI_count_files <- list.files(DNA_file_path,pattern=regex,full.names = TRUE)
  names(DNA_UMI_count_files) <- basename(DNA_UMI_count_files) %>% sub(regex, "", .)
  raw_DNA_count <- map_dfr(DNA_UMI_count_files,
                           ~fread(.x,sep='\t',col.names=c('candidate','count')), .id = 'sample_id',nThread=48)
  
  df <- full_join(raw_DNA_count,raw_RNA_count) %>%
    mutate(candidate=gsub('neg_81_1-neg','neg81',candidate),
           candidate=gsub('synthetic_inert_1','synthetic',candidate),
           candidate=gsub('seq1305\\(\\+\\);strong','strong',candidate),
           candidate=gsub("ATF2\\|FOSB::JUN\\|FOSL1::JUND\\|JUN::JUNB\\|ATF3","FOSB::JUN",candidate),
           candidate=gsub("FOSL1::JUNB\\|FOS::JUNB\\|FOS::JUND\\|BATF\\|BATF::JUN\\|FOSL1::JUND\\|JUN::JUNB", "BATF::JUN",candidate))
  return(df)
}
# quantile_normalise <- function(df,c){
#   quantile_list <- list()
#   for(b in c("strong","random","strong|random")){
#     if(b=="strong|random"){
#       matrix <- df %>% filter(cell==c) %>% 
#         separate(candidate,into=c('candidate','background'),sep='\\.') %>% 
#         filter(!grepl(b,background)) %>% 
#         distinct(candidate,background,RNA1_log2fc,RNA2_log2fc,RNA3_log2fc) %>% 
#         pivot_longer(cols=c(-candidate,-background),names_to='sample_id',values_to='log2fc') %>% 
#         pivot_wider(names_from=c(background,sample_id),names_sep="@",values_from=log2fc) %>% 
#         column_to_rownames(var='candidate') %>% as.matrix()
#     } 
#     else{
#       matrix <- df %>% filter(cell==c) %>% 
#         separate(candidate,into=c('candidate','background'),sep='\\.') %>% 
#         filter(grepl(b,background)) %>% 
#         distinct(candidate,background,RNA1_log2fc,RNA2_log2fc,RNA3_log2fc) %>% 
#         pivot_longer(cols=c(-candidate,-background),names_to='sample_id',values_to='log2fc') %>% 
#         pivot_wider(names_from=c(background,sample_id),names_sep="@",values_from=log2fc) %>% 
#         column_to_rownames(var='candidate') %>% as.matrix()
#     }
#     quantile <- normalize.quantiles(matrix)
#     colnames(quantile) <- colnames(matrix)
#     rownames(quantile) <- rownames(matrix)
#     quantile_df <- quantile %>% as.data.frame() %>% rownames_to_column(var='candidate') %>% 
#       pivot_longer(cols=c(-candidate),values_to='qn_log2fc',names_to='sample') %>% 
#       separate(sample,into=c('background','sample_id'),sep="@") %>% 
#       mutate(candidate=paste(candidate,background,sep='.'),sample_id=gsub("log2fc","qn_log2fc",sample_id)) %>% select(-background) %>%
#       pivot_wider(names_from=sample_id,values_from=qn_log2fc) %>% 
#       group_by(candidate) %>% 
#       mutate(RNA123_qn_log2fc_mean=mean(c(RNA1_qn_log2fc,RNA2_qn_log2fc,RNA3_qn_log2fc),na.rm = TRUE),
#              RNA123_qn_log2fc_median=median(c(RNA1_qn_log2fc,RNA2_qn_log2fc,RNA3_qn_log2fc),na.rm = TRUE),
#              RNA123_qn_log2fc_sd=sd(c(RNA1_qn_log2fc,RNA2_qn_log2fc,RNA3_qn_log2fc),na.rm = TRUE),
#              RNA123_qn_log2fc_cv=(RNA123_qn_log2fc_sd/RNA123_qn_log2fc_mean)) %>% ungroup() %>%
#       mutate(cell=c)
#     quantile_list[[paste(c,b,sep="-")]] <- quantile_df
#   }
#   return(quantile_list)
# }

###### RAW COUNTS
raw_counts_df <- parse_file(RNA_file_path="/projects/holtlab_prj/scratch/23_miniP/data/processed/STARRseq/RNA/02-candidate_counts",
                            DNA_file_path = "/projects/holtlab_prj/scratch/23_miniP/data/processed/STARRseq/DNA/30-961947925/02-mapped_reads",
                            regex='.mismatch0.tsv.rawcount.tsv$')
raw_count <- raw_counts_df %>% filter(grepl("synthetic|strong",candidate)) %>%
  mutate(sample_id=paste(sample_id,"_raw_readcount",sep='')) %>% pivot_wider(names_from=sample_id,values_from=count) %>%
  mutate(Jurkat_RNA123_raw_readcount_mean = apply(across(starts_with("Jurkat_RNA")), 1,function(x) mean(x[x > 0], na.rm = TRUE)),
         Jurkat_RNA123_raw_readcount_median = apply(across(starts_with("Jurkat_RNA")), 1, function(x) median(x[x > 0], na.rm = TRUE)),
         Jurkat_RNA123_raw_readcount_sd = apply(across(starts_with("Jurkat_RNA")), 1, function(x) sd(x[x > 0], na.rm = TRUE)),
         K562_RNA123_raw_readcount_mean = apply(across(starts_with("K562_RNA")), 1,function(x) mean(x[x > 0],na.rm = TRUE)),
         K562_RNA123_raw_readcount_median = apply(across(starts_with("K562_RNA")), 1, function(x) median(x[x > 0], na.rm = TRUE)),
         K562_RNA123_raw_readcount_sd = apply(across(starts_with("K562_RNA")), 1, function(x) sd(x[x > 0], na.rm = TRUE)),
         DNA12345678_raw_readcount_mean = apply(across(starts_with("DNA")), 1, function(x) mean(x[x > 0], na.rm = TRUE)),
         DNA12345678_raw_readcount_median = apply(across(starts_with("DNA")), 1, function(x) median(x[x > 0], na.rm = TRUE)),
         DNA12345678_raw_readcount_sd = apply(across(starts_with("DNA")), 1, function(x) sd(x[x > 0], na.rm = TRUE))) %>%
  pivot_longer(cols = contains("RNA"),names_to = c("cell", ".value"), names_pattern = "^(Jurkat|K562)_(.*)")


####### DEDUP COUNTS
dedup_counts <- parse_file(RNA_file_path="/projects/holtlab_prj/23_miniP/data/processed/STARRseq/RNA/02-candidate_counts",
                           DNA_file_path = "/projects/holtlab_prj/23_miniP/data/processed/STARRseq/DNA/30-961947925/03-candidate_counts",
                           regex='.mismatch0.counts.0.tsv$')

dedup <- dedup_counts %>% filter(grepl("synthetic|strong",candidate)) %>% 
  mutate(sample_id=paste(sample_id,"_deduplicated",sep='')) %>% pivot_wider(names_from=sample_id,values_from=count) %>%
  mutate(Jurkat_RNA123_deduplicated_mean = apply(across(starts_with("Jurkat_RNA")), 1,function(x) mean(x[x > 0], na.rm = TRUE)),
         Jurkat_RNA123_deduplicated_median = apply(across(starts_with("Jurkat_RNA")), 1, function(x) median(x[x > 0], na.rm = TRUE)),
         Jurkat_RNA123_deduplicated_sd = apply(across(starts_with("Jurkat_RNA")), 1, function(x) sd(x[x > 0], na.rm = TRUE)),
         K562_RNA123_deduplicated_mean = apply(across(starts_with("K562_RNA")), 1,function(x) mean(x[x > 0],na.rm = TRUE)),
         K562_RNA123_deduplicated_median = apply(across(starts_with("K562_RNA")), 1, function(x) median(x[x > 0], na.rm = TRUE)),
         K562_RNA123_deduplicated_sd = apply(across(starts_with("K562_RNA")), 1, function(x) sd(x[x > 0], na.rm = TRUE)),
         DNA12345678_deduplicated_mean = apply(across(starts_with("DNA")), 1, function(x) mean(x[x > 0], na.rm = TRUE)),
         DNA12345678_deduplicated_median = apply(across(starts_with("DNA")), 1, function(x) median(x[x > 0], na.rm = TRUE)),
         DNA12345678_deduplicated_sd = apply(across(starts_with("DNA")), 1, function(x) sd(x[x > 0], na.rm = TRUE))) %>%
  pivot_longer(cols = contains("RNA"),names_to = c("cell", ".value"), names_pattern = "^(Jurkat|K562)_(.*)")

####### CPM COUNTS
cpm <- dedup_counts %>% filter(grepl("synthetic|strong",candidate)) %>% 
  pivot_wider(names_from = sample_id, values_from = count,values_fill=0) %>%
  #must have count >=10 in all dna replicates
  #remove candidates with 0 RNA in all replicates
  filter(if_all(starts_with("DNA"), ~ . >= 10), if_any(contains("RNA"), ~ . != 0)) %>% 
  pivot_longer(cols=c(-candidate),names_to='sample_id',values_to='count') %>% 
  mutate(sample_id=paste(sample_id,"_cpm",sep='')) %>% 
  group_by(sample_id) %>% mutate(library_size=sum(count)) %>% ungroup() %>%
  mutate(cpm=((count)/(library_size))*1000000)

cpm_library_sizes <- cpm %>% distinct(candidate,library_size,sample_id) %>%
  pivot_wider(names_from=sample_id,values_from=library_size,names_prefix="libsize_")

cpm_df <- cpm %>% select(-count,-library_size) %>% unique() %>% 
  pivot_wider(names_from=sample_id,values_from=cpm) %>% 
  #only use counts > 0 for the summary statistics
  mutate(Jurkat_RNA123_cpm_mean = apply(across(starts_with("Jurkat_RNA")), 1,function(x) mean(x[x > 0], na.rm = TRUE)),
         Jurkat_RNA123_cpm_median = apply(across(starts_with("Jurkat_RNA")), 1, function(x) median(x[x > 0], na.rm = TRUE)),
         Jurkat_RNA123_cpm_sd = apply(across(starts_with("Jurkat_RNA")), 1, function(x) sd(x[x > 0], na.rm = TRUE)),
         K562_RNA123_cpm_mean = apply(across(starts_with("K562_RNA")), 1,function(x) mean(x[x > 0],na.rm = TRUE)),
         K562_RNA123_cpm_median = apply(across(starts_with("K562_RNA")), 1, function(x) median(x[x > 0], na.rm = TRUE)),
         K562_RNA123_cpm_sd = apply(across(starts_with("K562_RNA")), 1, function(x) sd(x[x > 0], na.rm = TRUE)),
         DNA12345678_cpm_mean = apply(across(starts_with("DNA")), 1, function(x) mean(x[x > 0], na.rm = TRUE)),
         DNA12345678_cpm_median = apply(across(starts_with("DNA")), 1, function(x) median(x[x > 0], na.rm = TRUE)),
         DNA12345678_cpm_sd = apply(across(starts_with("DNA")), 1, function(x) sd(x[x > 0], na.rm = TRUE)),
  ) %>%
  #collapse DNA replicates into one sample
  #remove where the DNA sample ==0
  mutate(DNA_cpm = rowMeans(across(starts_with("DNA") & !matches("qn|sd|mean|median")), na.rm = TRUE)) %>% 
  filter(DNA_cpm > 0 ) %>% 
  pivot_longer(cols = contains("RNA"),names_to = c("cell", ".value"), names_pattern = "^(Jurkat|K562)_(.*)")

###LOG2FC CPM
log2fc <- cpm_df %>%
  mutate(RNA1_log2fc = ifelse(RNA1_cpm > 0 & DNA_cpm > 0, log2(RNA1_cpm / DNA_cpm), NA),
         RNA2_log2fc = ifelse(RNA2_cpm > 0 & DNA_cpm > 0, log2(RNA2_cpm / DNA_cpm), NA),
         RNA3_log2fc = ifelse(RNA3_cpm > 0 & DNA_cpm > 0, log2(RNA3_cpm / DNA_cpm), NA)) %>%
  #only calculate log2fc for oligos where there are at least two replicates for a given cell type.
  rowwise() %>%
  mutate(num_na = sum(is.na(c(RNA1_log2fc, RNA2_log2fc, RNA3_log2fc)))) %>%
  filter(num_na < 2) %>% select(-num_na) %>% 
  group_by(candidate,cell) %>% 
  #calculate summaries
  mutate(RNA123_log2fc_mean=mean(c(RNA1_log2fc,RNA2_log2fc,RNA3_log2fc),na.rm = TRUE),
         RNA123_log2fc_median=median(c(RNA1_log2fc,RNA2_log2fc,RNA3_log2fc),na.rm = TRUE),
         RNA123_log2fc_sd=sd(c(RNA1_log2fc,RNA2_log2fc,RNA3_log2fc),na.rm = TRUE),
         RNA123_log2fc_cv=(RNA123_log2fc_sd/RNA123_log2fc_mean)) %>% ungroup()

background_log2fc <- log2fc %>%
  transmute(candidate,'RNA1_background'=RNA1_log2fc,'RNA2_background'=RNA2_log2fc,'RNA3_background'=RNA3_log2fc, cell) %>%
  separate(candidate,into=c('oligo','background'),sep='\\.') %>% 
  filter(background=="synthetic",oligo=='spacer-spacer-spacer') %>% 
  select(-oligo)

log2fc_w_background <- log2fc %>% filter(!grepl('spacer-spacer-spacer',candidate)) %>% 
  separate(candidate,into=c('oligo','background'),sep='\\.',remove=FALSE) %>% 
  select(-oligo) %>% 
  left_join(background_log2fc,relationship='many-to-many',by=c("cell","background")) %>% 
  mutate(RNA1_log2fc_adj=RNA1_log2fc-RNA1_background,
         RNA2_log2fc_adj=RNA2_log2fc-RNA2_background,
         RNA3_log2fc_adj=RNA3_log2fc-RNA3_background) %>% 
  select(-RNA1_background,-RNA2_background,-RNA3_background) %>%
  full_join(log2fc %>% filter(grepl('spacer-spacer-spacer',candidate)))

# QN LOG2FC
#jurkat_norm <- quantile_normalise(log2fc,"Jurkat") %>% rbindlist()
#K562_norm <- quantile_normalise(log2fc,"K562") %>% rbindlist()
#quantile_normalised <- rbind(jurkat_norm,K562_norm)

full_df <- full_join(raw_count,dedup) %>% full_join(cpm_library_sizes) %>% full_join(log2fc_w_background) %>% #full_join(quantile_normalised) %>% 
  separate(candidate,into=c('oligo','background'),sep='\\.') %>%
  separate(oligo,into=c('m1','m2','m3'),sep='-',remove=FALSE) %>% 
  mutate(m1_orientation=ifelse(grepl("revc",m1),"rev","fwd"),
         m2_orientation=ifelse(grepl("revc",m2),"rev","fwd"),
         m3_orientation=ifelse(grepl("revc",m3),"rev","fwd"),
         m1=gsub("revc","",m1),
         m2=gsub("revc","",m2),
         m3=gsub("revc","",m3)) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))
df1 <- full_df %>% select(oligo,m1,m1_orientation,m2,m2_orientation,m3,m3_orientation,background,cell)
df2 <- full_df %>% select(-oligo,-m1,-m1_orientation,-m2,-m2_orientation,-m3,-m3_orientation,-background,-cell)
final <- cbind(df1,df2)



#MAKE SURE FILES ONLY CONTAINS OLIGOS WITH LOG2FC FOR 2+ REPLICATES FOR A GIVEN CELL TYPE
fwrite(final,'/projects/holtlab_prj/23_miniP/data/processed/STARRseq/candidate.master.count.file.reduced.tsv',sep='\t',quote=FALSE)  
monotypic <- final %>% select(oligo,background,cell, m1,m2,m3) %>% filter(background!="strong") %>%
  pivot_longer(cols=c(-oligo,-background,-cell)) %>% group_by(oligo,background,cell) %>% filter(value!="spacer") %>% 
  filter(n_distinct(value)==1) %>% ungroup()
final %>% filter(oligo %in% monotypic$oligo| oligo=="spacer-spacer-spacer",background!='random',background!='strong') %>%
 rowwise() %>% mutate(num_na = sum(is.na(c(RNA1_log2fc, RNA2_log2fc, RNA3_log2fc)))) %>% filter(num_na < 2) %>% select(-num_na) %>%
 fwrite('/projects/holtlab_prj/23_miniP/data/processed/STARRseq/candidate.master.count.file.monotypic.reduced.tsv',sep='\t',quote=FALSE)

# monotypic_random <- final %>% select(oligo,background,cell, m1,m2,m3) %>% filter(background=="random") %>%
#   pivot_longer(cols=c(-oligo,-background,-cell)) %>% group_by(oligo,background,cell) %>% filter(value!="spacer") %>% 
#   filter(n_distinct(value)==1) %>% ungroup()
# final %>% filter(oligo %in% monotypic_random$oligo| oligo=="spacer-spacer-spacer",background=='random') %>%
#  rowwise() %>% mutate(num_na = sum(is.na(c(RNA1_log2fc, RNA2_log2fc, RNA3_log2fc)))) %>% filter(num_na < 2) %>% select(-num_na) %>%
#  fwrite('/projects/holtlab_prj/23_miniP/data/processed/STARRseq/candidate.master.count.file.monotypic.random.tsv',sep='\t',quote=FALSE)



ditypic <- final %>% select(oligo,background,cell, m1,m2,m3) %>% filter(background!="strong") %>%
  pivot_longer(cols=c(-oligo,-background,-cell)) %>% group_by(oligo,background,cell) %>% filter(value!="spacer") %>% 
  filter(n_distinct(value)==2) %>% ungroup()
final %>% filter(oligo %in% ditypic$oligo| oligo=="spacer-spacer-spacer",background!='random',background!="strong") %>%
 rowwise() %>% mutate(num_na = sum(is.na(c(RNA1_log2fc, RNA2_log2fc, RNA3_log2fc)))) %>% filter(num_na < 2) %>% select(-num_na) %>%
 fwrite('/projects/holtlab_prj/23_miniP/data/processed/STARRseq/candidate.master.count.file.ditypic.reduced.tsv',sep='\t',quote=FALSE)

# ditypic_random <- final %>% select(oligo,background,cell, m1,m2,m3) %>% filter(background=="random") %>%
#   pivot_longer(cols=c(-oligo,-background,-cell)) %>% group_by(oligo,background,cell) %>% filter(value!="spacer") %>% 
#   filter(n_distinct(value)==2) %>% ungroup()
# final %>% filter(oligo %in% ditypic_random$oligo| oligo=="spacer-spacer-spacer",background=='random') %>%
#  rowwise() %>% mutate(num_na = sum(is.na(c(RNA1_log2fc, RNA2_log2fc, RNA3_log2fc)))) %>% filter(num_na < 2) %>% select(-num_na) %>%
#  fwrite('/projects/holtlab_prj/23_miniP/data/processed/STARRseq/candidate.master.count.file.ditypic.random.tsv',sep='\t',quote=FALSE)


tritypic <- final %>% select(oligo,background,cell, m1,m2,m3) %>% filter(background!="strong") %>%
  pivot_longer(cols=c(-oligo,-background,-cell)) %>% group_by(oligo,background,cell) %>% filter(value!="spacer") %>% 
  filter(n_distinct(value)==3) %>% ungroup()
final %>% filter(oligo %in% tritypic$oligo| oligo=="spacer-spacer-spacer",background!="strong") %>%
 rowwise() %>% mutate(num_na = sum(is.na(c(RNA1_log2fc, RNA2_log2fc, RNA3_log2fc)))) %>% filter(num_na < 2) %>% select(-num_na) %>%
fwrite('/projects/holtlab_prj/23_miniP/data/processed/STARRseq/candidate.master.count.file.tritypic.reduced.tsv',sep='\t',quote=FALSE)

# tritypic_random <- final %>% select(oligo,background,cell, m1,m2,m3) %>% filter(background=="random") %>%
#   pivot_longer(cols=c(-oligo,-background,-cell)) %>% group_by(oligo,background,cell) %>% filter(value!="spacer") %>% 
#   filter(n_distinct(value)==3) %>% ungroup()
# final %>% filter(oligo %in% tritypic_random$oligo| oligo=="spacer-spacer-spacer",background=='random') %>%
#  rowwise() %>% mutate(num_na = sum(is.na(c(RNA1_log2fc, RNA2_log2fc, RNA3_log2fc)))) %>% filter(num_na < 2) %>% select(-num_na) %>%
# fwrite('/projects/holtlab_prj/23_miniP/data/processed/STARRseq/candidate.master.count.file.tritypic.random.tsv',sep='\t',quote=FALSE)
