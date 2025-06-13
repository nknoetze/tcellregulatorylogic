library(data.table)
library(janitor)
library(stringr)
library(tidyverse)

sample_subject <- fread('/projects/dice_dart/metadata/sample_subject.txt') %>%
  filter(SAMPLE_USE=="Seq_RNA") %>% clean_names('all_caps') %>%
  rename("DBGAP_SAMPLE_ID"=DB_GA_P_SAMPLE_ID) %>% rename("DBGAP_SUBJECT_ID"=DB_GA_P_SUBJECT_ID) %>%
  rename("BIOSAMPLE_ACCESSION"=BIO_SAMPLE_ACCESSION) %>% rename("SUBJECT_SAMPLE_ID"=SAMPLE_ID)
subject_phentoypes <- fread('/projects/dice_dart/metadata/subject_phenotypes.txt',
                            col.names = c('DBGAP_SUBJECT_ID','SUBJECT_ID','AGE','RACE','SEX')) 
sample_attributes <- fread('/projects/dice_dart/metadata/sample_attributes.txt') %>%
  filter(ANALYTE_TYPE=="RNA") %>% clean_names('all_caps') %>% rename("DBGAP_SAMPLE_ID"=DB_GA_P_SAMPLE_ID) %>%
  rename("SUBJECT_SAMPLE_ID"=SAMPLE_ID)
sra <- fread('/projects/dice_dart/metadata/SraRunTable.txt') %>% 
  select(-LoadDate,-MBases,-MBytes,-AvgSpotLen,-DATASTORE_filetype,-DATASTORE_provider,-Consent,-is_tumor,-histological_type,-sex,
         -Platform,-ReleaseDate,-study_name,-study_design,-Sample_Name,-biospecimen_repository_sample_id,-InsertSize,-LibrarySource) %>%
  rename("SUBJECT_ID"=submitted_subject_id,"SUBJECT_SAMPLE_ID"=Library_Name) %>% clean_names('all_caps')

sample_metadata <- full_join(sample_subject,sample_attributes) %>%
  full_join(subject_phentoypes) %>% full_join(sra)
sample_metadata <- sample_metadata %>% mutate(SAMPLE_ID=1:nrow(sample_metadata)) %>%
  mutate(SAMPLE_ID=str_pad(SAMPLE_ID,4,pad=0)) %>% mutate(SAMPLE_ID=paste("DICE",SAMPLE_ID,sep='_')) %>%
  rename('sample_id'=SAMPLE_ID,'cell_type'=HISTOLOGICAL_TYPE)
fwrite(sample_metadata,'/projects/nknoetze_prj/promoter_prj/results/tidy/dice_master_metadata_rna.tsv',sep='\t',quote=FALSE)
