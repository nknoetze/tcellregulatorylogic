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
    mutate(candidate=gsub('synthetic_seq','synthetic',candidate),
           candidate=gsub('seq1305\\(\\+\\);strong','strong',candidate),
           candidate=gsub("ATF2\\|FOSB::JUN\\|FOSL1::JUND\\|JUN::JUNB\\|ATF3","FOSB::JUN",candidate),
           candidate=gsub("FOSL1::JUNB\\|FOS::JUNB\\|FOS::JUND\\|BATF\\|BATF::JUN\\|FOSL1::JUND\\|JUN::JUNB", "BATF::JUN",candidate))
  return(df)
}
####### DEDUP COUNTS
dedup_counts <- parse_file(RNA_file_path="/projects/nknoetze_prj/ocr_prj/results/test_run/RNA/02-candidate_counts",
                           DNA_file_path = "/projects/nknoetze_prj/ocr_prj/results/test_run/DNA/03-candidate_counts",
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

cpm_library_sizes <- cpm %>% filter(grepl("synthetic|strong",candidate)) %>% distinct(candidate,library_size,sample_id) %>%
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
  filter(num_na < 2) %>% select(-num_na) 

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

full_df <- full_join(dedup,cpm_library_sizes) %>% full_join(log2fc_w_background) %>%
  separate(candidate,into=c('oligo','background'),sep='\\.') %>%
  separate(oligo,into=c('m1','m2','m3'),sep='-',remove=FALSE) %>% 
  mutate(m1_orientation=ifelse(grepl("revc",m1),"rev","fwd"),
         m2_orientation=ifelse(grepl("revc",m2),"rev","fwd"),
         m3_orientation=ifelse(grepl("revc",m3),"rev","fwd"),
         m1=gsub("revc","",m1),
         m2=gsub("revc","",m2),
         m3=gsub("revc","",m3),
         m1_tmp=ifelse(m1=='spacer',m1,paste(m1,m1_orientation,sep='')),
         m2_tmp=ifelse(m2=='spacer',m2,paste(m2,m2_orientation,sep='')),
         m3_tmp=ifelse(m3=='spacer',m3,paste(m3,m3_orientation,sep='')),
         oligo=paste(m1_tmp,m2_tmp,m3_tmp,sep='-')) %>%
  select(-m1_tmp,-m2_tmp,-m3_tmp) %>% 
  mutate(across(where(is.numeric), ~ round(., 4)))
df1 <- full_df %>% select(oligo,m1,m1_orientation,m2,m2_orientation,m3,m3_orientation,background,cell)
df2 <- full_df %>% select(-oligo,-m1,-m1_orientation,-m2,-m2_orientation,-m3,-m3_orientation,-background,-cell)
final <- cbind(df1,df2)


#MAKE SURE FILES ONLY CONTAINS OLIGOS WITH LOG2FC FOR 2+ REPLICATES FOR A GIVEN CELL TYPE
fwrite(final,'candidate.master.count.file.tsv',sep='\t',quote=FALSE)  
# monotypic <- final %>% select(oligo,background,cell, m1,m2,m3) %>% filter(background!="strong") %>%
#   pivot_longer(cols=c(-oligo,-background,-cell)) %>% group_by(oligo,background,cell) %>% filter(value!="spacer") %>% 
#   filter(n_distinct(value)==1) %>% ungroup()
# final %>% filter(oligo %in% monotypic$oligo| oligo=="spacer-spacer-spacer",background!='random',background!='strong') %>%
#   rowwise() %>% mutate(num_na = sum(is.na(c(RNA1_log2fc, RNA2_log2fc, RNA3_log2fc)))) %>% filter(num_na < 2) %>% select(-num_na) %>%
#   fwrite('/projects/nknoetze_prj/ocr_prj/results/test_run/STARR-seq/candidate.master.count.monotypic.file.tsv',sep='\t',quote=FALSE)
# 
# 
# 
# ditypic <- final %>% select(oligo,background,cell, m1,m2,m3) %>% filter(background!="strong") %>%
#   pivot_longer(cols=c(-oligo,-background,-cell)) %>% group_by(oligo,background,cell) %>% filter(value!="spacer") %>% 
#   filter(n_distinct(value)==2) %>% ungroup()
# final %>% filter(oligo %in% ditypic$oligo| oligo=="spacer-spacer-spacer",background!='random',background!="strong") %>%
#   rowwise() %>% mutate(num_na = sum(is.na(c(RNA1_log2fc, RNA2_log2fc, RNA3_log2fc)))) %>% filter(num_na < 2) %>% select(-num_na) %>%
#   fwrite('/projects/nknoetze_prj/ocr_prj/results/test_run/STARR-seq/candidate.master.count.ditypic.file.tsv',sep='\t',quote=FALSE)