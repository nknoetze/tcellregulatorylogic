options(width=60)
options(scipen = 999)
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(pheatmap))
suppressMessages(library(knitr))
suppressMessages(library(ggh4x))

#Goal: The goal of this analysis is to find a list of T-cell specific that have moderate - high expression in T-cells with little expression in non-T cells.
#
#The gene counts used in this analysis underwent were used with DESeq2 to normalise the counts using the median of ratios method. This method allows for the comparision of gene expression values between independent samples by correcting for RNA composition and library size. To account for gene length, the normalised reads obtained from DeSeq2 were divided by the effective gene length of the corresponding gene for each sample independently.
#-------------------------
#set up the metadata 
tcells <- c('cd4_naive','cd8_naive','cd4_n_stim','cd8_n_stim','th1','th2',
            'th17','th1-17','tfh','treg_memory', 'treg_naive')
cell_names <- c('Naive B-cells','Brain Endothelial Cells','CD34+ HSC - Cord Blood','CD34+ HSC - Peripheral',
                'CD4+ Naive T cells','Stimulated CD4+ Naive T cells','CD8+ Naive T cells',
                'Stimulated CD8+ Naive T cells','Classical Monocytes','Aortic Endothelial Cells','Brain Microvascular Endothelial Cells',
                'Hepatocyte Cell Line','Pulmonary Endothelial Cells','Non-Classical Monocytes','Oocyte Cell',
                'Sperm Cell', 'Differentiated Hepatocytes', 'T-Follicular Helper Cells', 'T-Helper 1 cells', 'T-Helper 1/17 cells',
                'T-Helper 17 cells', 'T-Helper 2 cells', 'Memory Regulatory T cells', 'Naive Regulatory T cells')
metadata_path='/projects/nknoetze_prj/promoter_prj/results/tidy/sra_dice_master_metadata.tsv'
metadata <- fread(metadata_path)%>%
  filter(cell_type!="nk_cd16pos") %>%
  mutate(cell_group=ifelse(cell_type %in% tcells,'Tcell',cell_type)) %>%
  select(cell_group,sample_id, cell_type)

#set up the names for the all the cell types
cell_types <- metadata %>% select(cell_type,cell_group) %>% 
  unique() %>% arrange(cell_type) %>%
  cbind(cell_names) 

offtarget_celltypes <- cell_types %>% filter(!(cell_type %in% tcells))
tcell_types <- cell_types %>% filter(cell_group=='Tcell') 

gene_types <- c('protein_coding','IG_C_gene','TR_C_gene')
#EXCLUDE the Y chromosome, remove mitochondrial genes
chromosomes <- paste('chr',c(1:22,'X'),sep='')
gencode_file <- fread('/projects/nknoetze_prj/references/annotations/gencode/gencode.v19.full.parsed.tsv',
                 col.names = c('chrom','type','start','end','strand','gene_id',
                               'transcript_id','gene_type','gene_name')) %>%
  filter(!grepl("^MT-",gene_name),chrom %in% chromosomes) %>% 
  select(gene_name,gene_id,gene_type) %>% unique()
#this list contains protein coding and the TR/IG genes
genes <- gencode_file %>%
  filter(gene_type %in% gene_types) %>%
  select(gene_name,gene_id) 

# converting read counts that were normalised by DESeq2 and adjusted for gene length
# to a dataframe and adding the metadata.
# Remove genes that have 0 expression across all samples
# includes protein coding AND TR/IG genes.********
counts_nobatchcorr <- fread('/projects/nknoetze_prj/promoter_prj/results/SRA-endothelial-hepatocyte-rsem/04-normalised-counts/sra.dice.counts.raw.sf.genel.adj.tsv.gz',nThread=12) %>%
  filter(rowSums(across(where(is.numeric)))>0) %>%
  reshape2::melt() %>%
  filter(gene_id %in% genes$gene_id) %>%
  rename('sample_id'=variable) %>%
  inner_join(metadata)

#-------- FUNCTION FOR CALCULATING MEDIAN EXPRESSION ------------------------
#--------------IN EACH CELL TYPE INDEPENDENTLY ------------------------------
get_median_expression <- function(counts_df){
# Get the median expression for each gene in each cell type, separately. Keep zeros
    median_expression <- counts_df %>%
      group_by(cell_type,gene_id) %>%
      summarise(median_expression=as.numeric(median(value))) %>%
      pivot_wider(names_from=cell_type,values_from = median_expression) 
    return(median_expression)
} 

## Step 1: Getting the Median Expression for Each Gene and cell type
#The median expression is calculated on a per-cell-type basis for each gene (protein coding, IGC, TRC genes).
median_expression_nobatchcorr <- get_median_expression(counts_nobatchcorr)
#median_expression_nobatchcorr %>% inner_join(genes) %>% fwrite('/projects/nknoetze_prj/ocr_prj/data/median_expr_nobatchcorr.tsv',sep='\t',quote=FALSE)