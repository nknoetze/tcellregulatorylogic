### --------------------------------------- ###
###          Load source scripts            ###
### --------------------------------------- ###
#source('/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/data/data_tidy_fx.R')

### --------------------------------------- ###
###        List of paths for tidying        ###
### --------------------------------------- ###

# METADATA FILES
meuleman_metadata_path="/projects/nknoetze_prj/ocr_prj/data/raw/meuleman_dnase/DHS_Index_and_Vocabulary_metadata_filtered.tsv"
# OCR FILES
meuleman_hg19ocr_path="/projects/nknoetze_prj/ocr_prj/data/raw/meuleman_dnase/dat_bin_FDR01_hg19.txt.gz"
ocr_sample_ids="/projects/nknoetze_prj/promoter_prj/data/metadata/DHS_Index_and_Vocabulary_metadata_modified.txt"
ocr_ids="/projects/nknoetze_prj/promoter_prj/data/raw/meuleman_dnase/DHS_Index_and_Vocabulary_hg19_WM20190703.txt.gz"
# GENCODE v19 FILE
gencodev19="/projects/nknoetze_prj/references/annotations/gencode/gencode.v19.transc.type.1based.tsv"
# INTERACTION DIRECTORY FOR THE SHARED, PROCESSED INTERACTIONS (CONTAINS CHRX)
interactions_dir ='/projects/nknoetze_prj/ocr_prj/data/processed/interactions/DICE/shared'

### --------------------------------------- ###
###             Tidy Metadata!              ###
### --------------------------------------- ###
print('creating OCR metadata for T cell samples')
tidy_metadata <- tidy_meuleman_data_tcells(meuleman_metadata_path)

### --------------------------------------- ###
###           Create OCR Datasets!          ###
### --------------------------------------- ###
print('creating OCR matrix for T cell samples. OCR in 1+ for any T cell sample')
tidy_tcell_1sample_ocrs <- tidy_meuleman_ocrs_tcell_samples(meuleman_hg19ocr_path,ocr_sample_ids,ocr_ids,tidy_metadata,min_sample=1)

tcell_ocr_subtypes_list <- tidy_meuleman_ocrs_tcells_subtype(tidy_tcell_1sample_ocrs,tidy_metadata,ocr_ids,min_sample=2)
print('creating OCR matrix. OCR in 2+ samples overall, and 1+ sample for CD8 AND CD4 T cells')
tcell_ocrs_subtype <- tcell_ocr_subtypes_list[1] %>% as.data.frame()

### --------------------------------------- ###
###       Tidy up Interaction Dataset       ###
### --------------------------------------- ###
#combine all the interaction files into one file. only keep interactions 
#that are significant in at least one cell type.
tidy_interaction_df_list <- tidy_interactions(interactions_dir)

### --------------------------------------- ###
###  Tidy up Merged Interaction Datasets    ###
### --------------------------------------- ###
#combine all the interaction files into one file. only keep interactions 
#that are significant in at least one cell type.
tidy_merged_interaction_connect8_df_list <- tidy_merged_interactions(interactions_dir,prefix='connectivity8')

### --------------------------------------- ###
###   Merge Overlapping OCRs                ###
### --------------------------------------- ###
merged_tcell_ocrs_subtype <- merge_ocrs(tcell_ocrs_subtype)

### --------------------------------------- ###
###         Write out the tidy files        ###
### --------------------------------------- ###
fwrite(merged_tcell_ocrs_subtype,'/projects/nknoetze_prj/ocr_prj/data/processed/meuleman_dnase/meuleman_tcell_2sample_1cd8_1cd4_ocrs_merged_tidy.tsv',sep="\t",quote=FALSE)

# Processed interactions
fwrite(tidy_interaction_df_list$Loop_2.5kb,'/projects/nknoetze_prj/ocr_prj/data/processed/interactions/DICE/shared/Loop_2.5kb/DICE_interactions_tidy.tsv',sep='\t',quote=FALSE)

# Processed Merged interactions
fwrite(tidy_merged_interaction_connect8_df_list$Loop_2.5kb,'/projects/nknoetze_prj/ocr_prj/data/processed/interactions/DICE/shared/Loop_2.5kb/connectivity8/DICE_interactions_merged_tidy.tsv',sep='\t',quote=FALSE)
