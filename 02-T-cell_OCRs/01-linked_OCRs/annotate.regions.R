suppressMessages(library(argparse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
### AUTHOR: N. KNOETZE ###
### ----------------------------- ###
###        SET UP ARGUMENTS       ###
### ----------------------------- ###
# create parser object
parser <- ArgumentParser()
parser$add_argument("-gc", "--gencode", help="Path to the parse gencode v19 file")
parser$add_argument("-ocr", "--file_to_annotate", help="Path to the OCR or CHIP seq peak file to annotate")
parser$add_argument("-o", "--outdir", help="Path to the output directory")
parser$add_argument("-fn", "--file_name", help="output file name")
parser$add_argument("-gt", "--gene_type", help="specify if this run is for lincRNA or protein coding")

args <- parser$parse_args()
gencode <- args$gencode
file_to_annotate <- args$file_to_annotate
outdir <- args$outdir
file_name <-  args$file_name
GENE_TYPE <- args$gene_type

#read in gencode file. Exclude ChrY and MT genomes
# only use tss for transcripts that are protein coding
if(GENE_TYPE=="longncRNA"){
  types <- c('lincRNA','antisense', 'sense_intronic', 'sense_overlapping')
}else if(GENE_TYPE=='protein_coding'){
  types <- c('protein_coding','TR_C_gene','IG_C_gene')
}
outdir <- paste(outdir,'/',GENE_TYPE,sep='')
print(outdir)
if(!(dir.exists(outdir))){
  dir.create(outdir)
}


gencode_v19 <- fread(gencode) %>% 
  #fread('/projects/nknoetze_prj/references/annotations/gencode/gencode.v19.transc.type.1based.tsv') %>% 
  filter(type=='transcript',gene_type %in% types, transcript_type %in% types) %>% 
  mutate(tss_start=ifelse(strand=="+",start-500,end-500),
         tss_end=ifelse(strand=="+",start+500,end+500)) %>%
  filter(!(grepl("chrY|chrM",chrom)))

unannotated_file <- fread(file_to_annotate) %>%
    #fread('/projects/nknoetze_prj/ocr_prj/data/processed/meuleman_dnase/meuleman_tcell_2sample_1cd8_1cd4_ocrs_merged_tidy.tsv') %>%
    #fread('/projects/nknoetze_prj/ocr_prj/data/processed/interactions/DICE/shared/tcell_merged_chip_peaks.tsv') %>%
    rename('chrom'=seqnames)

#Create Grange Objects for OCRs/chip peaks and the TSS
region_ranges <- makeGRangesFromDataFrame(unannotated_file,keep.extra.columns = TRUE,seqnames.field = c('chrom'))

#Create Grange Objects for Gene TSS
tss_ranges <- gencode_v19 %>%
  select(chrom,tss_start,tss_end,gene_id, gene_name,transcript_id) %>% 
  makeGRangesFromDataFrame(.,start.field = c('tss_start'),end.field = c('tss_end'), seqnames.field = c('chrom'),ignore.strand = TRUE,keep.extra.columns = TRUE)

####
map_genes <- function(unannotated_file, region_granges, tss_granges){
  #query-tss
  #subject-interactions
  #use within: ask: are there any TSS within the OCRs/chipseq peaks?
  tss_ocr_overlap <- findOverlaps(query = tss_granges,subject = region_granges,minoverlap = 1,type='any',select='all') 
  # get the gene ids for each interaction
  # append it to the df containing the overlaps
  mapping_df <- DataFrame(region_idx=subjectHits(tss_ocr_overlap),tss_idx=queryHits(tss_ocr_overlap)) %>% as.data.frame()
  mapping_df <- mapping_df %>% mutate(gene_name=tss_granges[tss_idx]$gene_name,
                                      gene_id=tss_granges[tss_idx]$gene_id,
                                      ocr_id=region_granges[region_idx]$ocr_id,transcript_id=tss_granges[tss_idx]$transcript_id) %>% 
    select(gene_name,gene_id,ocr_id,transcript_id) %>% mutate(ocr_type='promoter')
  ocr_tss_genes <- unannotated_file %>% left_join(mapping_df) %>% replace(is.na(.), 'none') %>% 
    mutate(ocr_type=ifelse(ocr_type=='none','non-promoter','promoter')) %>% unique()
}

annotated_file <- map_genes(unannotated_file,region_ranges,tss_ranges)

#out_file='/projects/nknoetze_prj/ocr_prj/data/processed/meuleman_dnase/meuleman_tcell_2sample_1cd8_1cd4_ocrs_merged_annotated_tidy.tsv'
#out_file='/projects/nknoetze_prj/ocr_prj/data/processed/interactions/DICE/shared/tcell_merged_chip_peaks_merged_annotated_tidy.tsv'

out_file <- paste(outdir,'/',file_name,'.tsv',sep='')
print(out_file)
fwrite(annotated_file,out_file,sep='\t',quote=FALSE)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#		      __		        __		      __		      __	   	#
#	   (___()'`;	  (___()'`;	  (___()'`;	  (___()'`;	    #
#	   /,	  /`	    /,	  /`	  /,	  /`	  /,	  /`	    #
#	   \\"--\\	    \\"--\\	    \\"--\\	    \\"--\\	   	  #
#					  									                          #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #