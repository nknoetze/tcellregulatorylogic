library(data.table)
library(tidyverse)

sra_metadata <- fread('/projects/nknoetze_prj/promoter_prj/results/tidy/sra_master_metadata.tsv') %>%
  rename(accession='geo_accession') %>%
  select(-library_selection)
dice_metadata <- fread('/projects/nknoetze_prj/promoter_prj/results/tidy/dice_master_metadata_rna_tidy.tsv') %>%
  janitor::clean_names(case='snake') %>%
  rename(accession='gap_accession',
         library_prep='library_selection') %>%
  mutate(age=as.character(age),
         subject_id=as.character(subject_id),
         organism='homo_sapiens',
         authors="DICE",
         avg_spot_len=50)

full_metadata <- full_join(dice_metadata,sra_metadata) %>%
  select(-bases,-bytes,-library_source,-body_site,-sample_use,-analyte_type,-is_tumor,-biospecimen_repository,-biosample_accession,-subject_sample_id) %>%
  mutate(cell_type=tolower(cell_type),
         cell_group=case_when(grepl("endothelial",cell_type) ~ "Endothelial Cells",
                              grepl("hepat",cell_type) ~ "Hepatocytes",
                              grepl("sperm",cell_type) ~ "Sperm Cell",
                              grepl("oocyte",cell_type) ~ "Oocyte Cell",
                              grepl("cd34+",cell_type) ~ "HSC",
                              TRUE ~ "Immune Cells"),
         broad_cell_type=case_when(grepl("monocytes",cell_type) ~ "Non/Classical Monocytes",
                                   grepl("nk_cd16pos", cell_type) ~ "NK Cells",
                                   grepl("b_naive", cell_type) ~ "Naive B Cells",
                                   grepl("^cd4|treg|tfh|^th", cell_type) ~ "Regulatory/Helper/CD4 T cells",
                                   grepl("cd8", cell_type) ~ "CD8 T cells",
                                   grepl("pulmonary", cell_type) ~ "Endothelial Cells - Pulmonary",
                                   grepl("brain", cell_type) ~ "Endothelial Cells - Brain",
                                   grepl("hepat", cell_type) ~ "Hepatocytes",
                                   grepl("aortic", cell_type) ~ "Endothelial Cells - Aortic",
                                   cell_type =="sperm_cell" ~ "Sperm Cells",
                                   cell_type =="oocyte_cell" ~ "Oocyte Cells",
                                   grepl("cordblood",cell_type) ~ "CD34+ HSC - Cord Blood",
                                   grepl("peripheral",cell_type) ~ "CD34+ HSC - Peripheral"),
         cell_label_name=case_when(grepl("monocytes",cell_type) ~ "Non/Classical Monocytes",
                                   grepl("nk_cd16pos", cell_type) ~ "NK Cells",
                                   grepl("b_naive", cell_type) ~ "Naive B Cells",
                                   grepl("cd4_naive", cell_type) ~ "Naive CD4 T Cells",
                                   grepl("cd8_naive", cell_type) ~ "Naive CD8 T Cells",
                                   grepl("cd4_n_stim", cell_type) ~ "Activated Naive CD4 T Cells",
                                   grepl("cd8_n_stim", cell_type) ~ "Activated Naive CD8 T Cells",
                                   grepl("^th|tfh", cell_type) ~ "Helper T Cells",
                                   grepl("^treg", cell_type) ~ "Regulatory T Cells",
                                   grepl("pulmonary", cell_type) ~ "Endothelial Cells - Pulmonary",
                                   grepl("brain", cell_type) ~ "Endothelial Cells - Brain",
                                   grepl("hepat", cell_type) ~ "Hepatocytes",
                                   grepl("aortic", cell_type) ~ "Endothelial Cells- Aortic",
                                   grepl("sperm", cell_type) ~ "Sperm Cell",
                                   grepl("oocyte", cell_type) ~ "Oocyte Cell",
                                   grepl("cordblood",cell_type) ~ "CD34+ HSC - Cord Blood",
                                   grepl("peripheral",cell_type) ~ "CD34+ HSC - Peripheral"))

fwrite(full_metadata, '/projects/nknoetze_prj/promoter_prj/results/tidy/sra_dice_master_metadata.tsv',sep='\t',quote=FALSE)
