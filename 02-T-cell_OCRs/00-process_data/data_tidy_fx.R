library(GenomicRanges)
library(data.table)
library(tidyverse)

### --------------------------------------- ###
###           Clean up Metadata File        ###
### --------------------------------------- ###

tidy_meuleman_data_tcells <- function(meuleman_metadata){
  metadata_df <- fread(meuleman_metadata)
  meuleuman_metadata_tcell_tidy <- metadata_df %>% 
    filter(Cell_Group=="T-cell",!(grepl("induced",Description))) %>%
    mutate(Cell_Type_Detailed=gsub("_"," ",Cell_Type_Detailed),
           Cell_Type_Detailed=gsub("Tcell|T-cell","T cell",Cell_Type_Detailed),
           sample_id=paste(`Biosample name`,`Altius Biosample ID`,sep='.'),
           Cell_Group=case_when(grepl("CD8+",Cell_Type) ~ "CD8_Tcell",
                                Cell_Type=="T-cell" ~ "Tcell",
                                TRUE ~ "CD4_Tcell"))
  return(meuleuman_metadata_tcell_tidy)
}

### --------------------------------------- ###
###      Create Base T-cell OCR file        ###
### --------------------------------------- ###

tidy_meuleman_ocrs_tcell_samples <- function(hg19_ocrs,sample_ids,ocr_ids,tidy_metadata,min_sample){
  #read in files and make some adjustments
  ocrs <- fread(hg19_ocrs,nThread=18)
  sample_ids <- fread(ocr_sample_ids) %>% mutate(sample_id=paste(`Biosample name`,`Altius Biosample ID`,sep='.')) %>% filter(!is.na(`library order`)) %>% select(sample_id)
  ocr_idx <- fread(ocr_ids, nThread=18) %>% mutate(ocr_id=paste0(seqname,":",start,'-',end,sep=''))
  ocr_idx_id <- ocr_idx %>% select(identifier,ocr_id)
  # Assign the column and rownames (sample, and ocr identifier, respectively)
  #unforunately, hg19 has 5 duplicate ocr regions after being lifted over.. so we need to use the unique ocr identifier from meuleman as the rownames until later filtering.
  colnames(ocrs) <- sample_ids$sample_id
  rownames(ocrs) <- ocr_idx$identifier
  
  #Get the list of OCR ids/identifiers that are on a canonical chromosome (no non-canonical, or Y)
  filtered_ocr_ids <- ocr_idx %>% filter(!(grepl("_|chrY",ocr_id))) %>% select(ocr_id,identifier)

  # get the OCRs for T cells,and only keep those annotated in X+ T cell samples.
  # remove OCRs that are on non-canonical or Y chromosomes.
  tcell_ocrs_tidy <- ocrs %>% select(matches(tidy_metadata$sample_id,colnames(.))) %>% filter(rowSums(.)>=min_sample) %>% 
    filter(rownames(.) %in% filtered_ocr_ids$identifier) %>% rownames_to_column(var='identifier') %>% inner_join(ocr_idx_id)
  return(tcell_ocrs_tidy)
}

### --------------------------------------- ###
###  Create Base T-cell OCR file By Subtype ###
### --------------------------------------- ###

tidy_meuleman_ocrs_tcells_subtype <- function(tcell_ocrs_tidy, metadata_tidy,ocr_ids,min_sample){
  df_list <- list()
  filtered_metadata <- metadata_tidy %>% filter(Cell_Group != "Tcell")
  ocr_idx <- fread(ocr_ids, nThread=18) %>% mutate(ocr_id=paste0(seqname,":",start,'-',end,sep='')) %>% select(ocr_id,identifier)
  tcell_ocrs_subtype_tidy <- tcell_ocrs_tidy %>% select(-ocr_id) %>%
    reshape2::melt() %>%
    rename('sample_id'=variable) %>%
    inner_join(filtered_metadata) %>%
    # since the values are either 0/1 for each OCR:sample, we can't just count the rows unless we filter for values=1
    # instead, just summarise the value which will give you the number of samples that have the OCR annotated
    group_by(identifier,Cell_Group) %>% summarise(n_samples=sum(value)) %>%
    pivot_wider(names_from=Cell_Group, values_from=n_samples) %>% inner_join(ocr_idx)
  
  #OCRs need to be in at least two T cell samples (for CD8 and CD4 T cells collectively)
  #But the OCR only needs to be in at least one sample for each cell group
  tcell_ocrs_subtype <- tcell_ocrs_tidy %>% select(-ocr_id) %>% 
    reshape2::melt() %>% rename('sample_id'=variable) %>% inner_join(filtered_metadata) %>% 
    group_by(identifier,Cell_Group) %>% summarise(n_samples=sum(value)) %>% 
    pivot_wider(names_from=Cell_Group, values_from=n_samples) %>% inner_join(ocr_idx) %>% 
    mutate(sum=(CD4_Tcell+CD8_Tcell)) %>% filter(sum>=2) %>% 
    filter(CD4_Tcell>=1 & CD8_Tcell >=1) %>% select(-sum)
  #save files.
  df_list[[1]] <- tcell_ocrs_subtype
  return(df_list)
}

### --------------------------------------- ###
###       Tidy up Interaction Dataset       ###
### --------------------------------------- ###
# Function to create a single dataframe of the interactions from the shared DICE data
# Contains the X chromosome! :)
tidy_interactions <- function(interactions_data_dir){
  interactions_df <- data.frame()
  tidy_interactions_list <- list()
  bin_sizes <- c('Loop_1kb','Loop_2.5kb','Loop_5kb')
  #Makes sure we are looking at the Stringent files and significant interactions (Q < 0.01) for the 2.5 and 1kb bins
  interaction_files <- list.files(interactions_data_dir,pattern="FitHiChIP.interactions_FitHiC_Q0.01.bed",recursive=TRUE,full.names=TRUE)
  stringent_files_indx <- grepl(pattern = "Stringent",interaction_files)
  stringent_interaction_files <- interaction_files[stringent_files_indx]
  #Get the 5kb interaction files which use the unfiltered bed file!!!!!! (To be consistent with previous results??)
  interaction_files_5kb <- list.files(interactions_data_dir,pattern="FitHiChIP.interactions_FitHiC.bed",recursive=TRUE,full.names=TRUE)
  #combine to get full list.
  all_interaction_files <- append(interaction_files_5kb,stringent_interaction_files)
  
  
  for(file in all_interaction_files){
    # Read in the interaction files for each cell type. From the file name, extract the celltype and the bin size.
  # Select the useful columns, and rename some of them, append the cell type to the dataframe
  # Remove interactions on the Y chromosome.
    split_file=strsplit(file,"/")
    celltype=split_file[[1]][10]
    binsize=split_file[[1]][11]
    print(paste('Reading file for ',celltype,' Bin size:',binsize,sep=''))
    df <- fread(file) %>% transmute('chrom1'=chr1,s1,e1,'chrom2'=chr2,s2,e2,cc,'qval'=`Q-Value_Bias`) %>%
      filter(chrom1!="chrY"|chrom2!="chrY") %>% mutate(cell_type=celltype,bin_size=binsize)
    print(paste('There are ',nrow(df), ' interactions.',sep=''))
    interactions_df <- rbind(interactions_df,df)
  }
  for(bin in bin_sizes){
    print(paste('Parsing data for ',bin,sep=''))
    # After parsing all of the files, combine them together into a single dataframe
    # Do this by Binsize to create a list of dataframes. One for each binsize
    # keep the interactions that are significant (Qval < 0.01) in at least one of the cell types
    # For interactions that are missing in a cell-type (ie were never detected), their contact count should be changed to 0
    # Also, the corresponding q value, (which will be NA), should be labelled as 1 (insignificant)
    wider_interaction_df <- interactions_df %>% filter(bin_size==bin) %>%
      pivot_wider(names_from=cell_type,values_from=c(qval,cc)) %>%
      filter_at(vars(starts_with("qval_")),any_vars(.<0.01)) %>%
      select(chrom1, s1,e1,chrom2,s2,e2,cc_Monocyte,qval_Monocyte,cc_NaiveK,qval_NaiveK,cc_NaiveB,qval_NaiveB,cc_CD4Naive,qval_CD4Naive,cc_CD8Naive,qval_CD8Naive) %>%
      mutate_at(vars(starts_with("cc_")),~(ifelse(is.na(.),0,.))) %>% replace(is.na(.),1)
    tidy_interactions_list[[bin]] <- wider_interaction_df
  }
  return(tidy_interactions_list)
}

### --------------------------------------- ###
###  Tidy up Merged Interaction Datasets    ###
### --------------------------------------- ###
# Function to create a single dataframe of the interactions from the shared DICE data
# Contains the X chromosome! :)
tidy_merged_interactions <- function(interactions_data_dir,prefix){
  merged_interactions_df <- data.frame()
  tidy_merged_interactions_list <- list()
  bin_sizes <- c('Loop_2.5kb')
  interaction_files <- list.files(interactions_data_dir,pattern="FitHiChIP.interactions_FitHiC_Q0.01.merged.bed",recursive=TRUE,full.names=TRUE)
  file_indx <- grepl(pattern = prefix,interaction_files)
  interaction_files <- interaction_files[file_indx]

  for(file in interaction_files){
    # Read in the interaction files for each cell type. From the file name, extract the celltype and the bin size.
    # Select the useful columns, and rename some of them, append the cell type to the dataframe
    # Remove interactions on the Y chromosome.
    split_file=strsplit(file,"/")
    celltype=split_file[[1]][10]
    binsize=split_file[[1]][11]
    print(paste('Reading merged interaction file for ',celltype,' Bin size:',binsize,sep=''))
    df <- fread(file) %>% transmute('chrom1'=chr1,'s1'=bin1_low,'e1'=bin1_high,'chrom2'=chr2,'s2'=bin2_low,'e2'=bin2_high,sumCC,StrongConn,'qval'=`Q-Value_Bias`) %>%
      filter(chrom1!="chrY"|chrom2!="chrY",qval<0.01) %>%
      mutate(cell_type=celltype,bin_size=binsize) %>% select(-qval) %>% unique()
    print(paste('There are ',nrow(df), ' interactions.',sep=''))
    merged_interactions_df <- rbind(merged_interactions_df,df)
  }
  for(bin in bin_sizes){
    print(paste('Parsing data for ',bin,sep=''))
    # After parsing all of the files, combine them together into a single dataframe
    # Do this by Binsize to create a list of dataframes. One for each binsize
    # For merged interactions that are missing in a cell-type their contact count should be changed to 0
    # Also, the corresponding StrongConn(which will be NA), should be labelled as 0 as well.
    wider_interaction_df <- merged_interactions_df %>% filter(bin_size==bin) %>%
      pivot_wider(names_from=cell_type,values_from=c(StrongConn,sumCC)) %>%
      select(chrom1,s1,e1,chrom2,s2,e2,sumCC_Monocyte,StrongConn_Monocyte,sumCC_NaiveK,StrongConn_NaiveK,sumCC_NaiveB,StrongConn_NaiveB,sumCC_CD4Naive,StrongConn_CD4Naive,sumCC_CD8Naive,StrongConn_CD8Naive) %>%
      #If there merged interaction is not present in a cell type, assign CC as 0 and StrongConn as 0
      mutate_at(vars(starts_with("sumCC_")),~(ifelse(is.na(.),0,.))) %>% replace(is.na(.),0)
    tidy_merged_interactions_list[[bin]] <- wider_interaction_df
  }
  return(tidy_merged_interactions_list)
}

### --------------------------------------- ###
###     MERGE OVERLAPPING OCR PEAKS         ###
### --------------------------------------- ###
merge_ocrs <- function(unmerged_ocrs){
  #combine the ocrs together (reduce) to obtain non-overlapping intervals  
  merged_ocrs <- unmerged_ocrs %>% 
    separate(ocr_id,into=c('chr','start','end'),sep=':|-') %>% 
    makeGRangesFromDataFrame() %>% GenomicRanges::reduce() %>% as.data.frame() %>% 
    mutate(ocr_id=paste(seqnames,':',start,'-',end,sep='')) %>% filter(seqnames!="chrY", seqnames!="chrM")
  return(merged_ocrs)
}
