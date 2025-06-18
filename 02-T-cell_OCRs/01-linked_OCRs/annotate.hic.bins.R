suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(argparse))
### ----------------------------- ###
###        SET UP ARGUMENTS       ###
### ----------------------------- ###
# create parser object
parser <- ArgumentParser()
parser$add_argument("-interaction_file", "--interaction_file", help="Path to the interaction file to annotate")
parser$add_argument("-is_merged", "--is_merged", help="specify if the interactions are merged interactions (merged)",default='no')
parser$add_argument("-gencode", "--gencode", help="Path to the gencode file to annotate")
parser$add_argument("-o", "--outdir", help="Path to the output directory")

args <- parser$parse_args()
gencode <- args$gencode
interaction_file <- args$interaction_file
is_merged <- args$is_merged
outdir <- args$outdir


#read in gencode file. Exclude ChrY and MT genomes
#only include transcription start sites for transcripts that are protein coding and tr\ig.
types <- c('protein_coding','TR_C_gene','IG_C_gene')

#read in gencode file. Exclude ChrY and MT genomes
#only include transcription start sites for transcripts that are protein coding and tr\ig, or the ncRNA class
#gencode_v19 <- fread('/projects/nknoetze_prj/references/annotations/gencode/gencode.v19.transc.type.1based.tsv') %>% 
gencode_v19 <- fread(gencode) %>% 
  filter(type=='transcript',gene_type %in% types, transcript_type %in% types) %>% 
  mutate(tss=ifelse(strand=="+",start,end)) %>%
  filter(!(grepl("chrY|chrM",chrom)))

bin_size=str_split(interaction_file,pattern = '/')[[1]][10]
#interactions <- fread('/projects/nknoetze_prj/ocr_prj/data/processed/interactions/DICE/shared/Loop_2.5kb/DICE_interactions_merged_tidy.tsv')
interactions <- fread(interaction_file)

#assign bin ids 
bin1_ids <- interactions %>% select(chrom1,s1,e1) %>% unique() %>% mutate(bin1_id=1:nrow(.))
bin2_ids <- interactions %>% select(chrom2,s2,e2) %>% unique() %>% mutate(bin2_id=1:nrow(.))
#label the interaction id (combination of bin 1 and 2 ids)
interactions <- interactions %>% inner_join(bin1_ids) %>% inner_join(bin2_ids) %>% 
  mutate(interaction_id=paste0(bin1_id,'_',bin2_id,sep='')) 

## Get the bin coordinates for bin 1 and 2. Label the index as the row number
bin1 <- interactions %>%
  transmute('chr'=chrom1,'start'=s1,'end'=e1,'bin_id'=bin1_id) %>% unique() %>% 
  mutate(interaction_idx=1:nrow(.))
bin2 <- interactions %>%
  transmute('chr'=chrom2,'start'=s2,'end'=e2,'bin_id'=bin2_id) %>% unique() %>%
  mutate(interaction_idx=1:nrow(.))

#Create Grange Objects for bins and the TSS
bin1_ranges <- makeGRangesFromDataFrame(bin1,keep.extra.columns = TRUE,ignore.strand=TRUE)
bin2_ranges <- makeGRangesFromDataFrame(bin2,keep.extra.columns = TRUE,ignore.strand=TRUE)
tss_ranges <- gencode_v19 %>%
  transmute(chrom,'tss_start'=tss,'tss_end'=tss,gene_id, gene_name) %>% 
  makeGRangesFromDataFrame(.,start.field = c('tss_start'),end.field = c('tss_end'), seqnames.field = c('chrom'),ignore.strand = TRUE,keep.extra.columns = TRUE)

####
map_genes <- function(interaction_bin_df, bin_granges, tss_granges){
  #query-tss
  #subject-interactions
  #use within: ask: are there any TSS within the Interactions?
  tss_bin_overlap <- findOverlaps(query = tss_granges,subject = bin_granges,type='within') 
  
  # get the gene ids for each interaction
  # append it to the df containing the overlaps
  mapping_df <- DataFrame(interaction_idx=subjectHits(tss_bin_overlap),tss_idx=queryHits(tss_bin_overlap)) %>% as.data.frame()
  mapping_df <- mapping_df %>% mutate(gene_name=tss_granges[tss_idx]$gene_name,gene_id=tss_granges[tss_idx]$gene_id) %>% 
    select(gene_name,gene_id,interaction_idx) %>%
    mutate(interaction_idx=as.integer(interaction_idx))
  bin_genes <- interaction_bin_df %>% left_join(mapping_df) %>% replace(is.na(.), "none")
  #for each bin, create a column containing all of the gene ids whose TSS are within the bin
  bin_gene_overlaps <- bin_genes %>% select(-bin_id) %>% unique() %>% aggregate(gene_id~interaction_idx, data=., paste, collapse = ";")
  bin_annotations <- bin_genes %>% select(-gene_name,-gene_id) %>% unique() %>% merge(bin_gene_overlaps)
}

bin1_annotations <- map_genes(bin1,bin1_ranges,tss_ranges) %>% transmute('chrom1'=chr,'s1'=start,'e1'=end,'bin1_gene_id'=gene_id,'bin1_id'=bin_id)
bin2_annotations <- map_genes(bin2,bin2_ranges,tss_ranges) %>% transmute('chrom2'=chr,'s2'=start,'e2'=end,'bin2_gene_id'=gene_id,'bin2_id'=bin_id)

#get rid of extra columns. 
#merge the new annotations
if(is_merged=='merged'){
  interactions_tss_within <- interactions %>% select(bin1_id,bin2_id,chrom1,s1,chrom2,
                                                     s2,e2,interaction_id,StrongConn_Monocyte,StrongConn_NaiveK,StrongConn_NaiveB,StrongConn_CD4Naive,StrongConn_CD8Naive) %>% 
    merge(bin1_annotations) %>% 
    merge(bin2_annotations,by=c('bin2_id','chrom2','s2','e2')) %>%
    select(interaction_id,bin1_id,chrom1,s1,e1,bin1_gene_id, bin2_id,chrom2,s2,e2,bin2_gene_id,StrongConn_Monocyte,StrongConn_NaiveK,StrongConn_NaiveB,StrongConn_CD4Naive,StrongConn_CD8Naive) %>%
    unique()
  outfile <- paste(outdir,'merged_interactions_tss_within.tsv',sep='')
} else{
  interactions_tss_within <- interactions %>% select(bin1_id,bin2_id,chrom1,s1,chrom2,
                                                     s2,e2,interaction_id,qval_Monocyte,qval_NaiveK,qval_NaiveB,qval_CD4Naive,qval_CD8Naive) %>% 
    merge(bin1_annotations) %>% 
    merge(bin2_annotations,by=c('bin2_id','chrom2','s2','e2')) %>%
    select(interaction_id,bin1_id,chrom1,s1,e1,bin1_gene_id, bin2_id,chrom2,s2,e2,bin2_gene_id,qval_Monocyte,qval_NaiveK,qval_NaiveB,qval_CD4Naive,qval_CD8Naive) %>%
    unique()
  outfile <- paste(outdir,'interactions_tss_within.tsv',sep='')

}

fwrite(interactions_tss_within,paste(outdir,'interactions_tss_within.tsv',sep=''),sep='\t',quote=FALSE)

