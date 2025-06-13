### OUT DIR VARIABLES ###
out_dir=config['OUTDIR']

pulled_ocrs_out_dir=config['PULLED_OCRS']
pulled_ocrs_out_dir_1=f"{out_dir}/{pulled_ocrs_out_dir}"

ocr_metric_files_out_dir=config['OCR_METRIC_FILES']
ocr_metric_files_out_dir_1=f"{out_dir}/{ocr_metric_files_out_dir}"

### VARIABLES ###
conda_env=config['CONDA_ENV']
gencode=config['GENCODE']
gene_list=config['GENE_LIST']
gene_type=config['GENE_TYPE']
filtered_ocr_file_dir=config['TIDY_OCR_FILE_DIR']
filtered_gene_file_dir=config['FILTERED_GENE_DIR']
filtered_ocr_data_name=config['FILTERED_OCR_DATA_NAME']
bin_size=config['BIN_SIZE']
merge_param=config['MERGE_PARAM']
interaction_annotation=config['INTERACTION_ANNOTATION']
interaction_dir=config['INTERACTIONS_DIR']

# # # --------------------------------------------------------------
# # # SET RULE ALL
# # # --------------------------------------------------------------

rule all:
	input:
		expand(f"{filtered_ocr_file_dir}""/{gene_type}/{filtered_ocr_data_name}_annotated_tidy.tsv",filtered_ocr_file_dir=filtered_ocr_file_dir,gene_type=gene_type,filtered_ocr_data_name=filtered_ocr_data_name),
		expand(f"{interaction_dir}/""{bin_size}/{merge_param}/{gene_type}/merged_interactions_tss_within.tsv",interaction_dir=interaction_dir,merge_param=merge_param,bin_size=bin_size,gene_type=gene_type),
		expand(f"{interaction_dir}""/{bin_size}/{gene_type}/interactions_tss_within.tsv",interaction_dir=interaction_dir,bin_size=bin_size,gene_type=gene_type),
		expand(f"{interaction_dir}/""{bin_size}/{merge_param}/{gene_type}/interactions_tss_within_tidy.tsv",interaction_dir=interaction_dir,merge_param=merge_param,bin_size=bin_size,gene_type=gene_type),
		expand(f"{pulled_ocrs_out_dir_1}/""{filtered_ocr_data_name}/{bin_size}/merged_interactions/{merge_param}/tss_within/{gene_type}/pulled_ocrs_master_by_gene.tsv",pulled_ocrs_out_dir_1=pulled_ocrs_out_dir_1,filtered_ocr_data_name=filtered_ocr_data_name,bin_size=bin_size,gene_type=gene_type,merge_param=merge_param),
		expand(f"{filtered_gene_file_dir}/""ranked_{gene_type}_list_filtered+{filtered_ocr_data_name}+{bin_size}+{merge_param}.tsv",filtered_gene_file_dir=filtered_gene_file_dir,gene_type=gene_type,filtered_ocr_data_name=filtered_ocr_data_name,bin_size=bin_size,merge_param=merge_param),
		expand(f"{ocr_metric_files_out_dir_1}/""{filtered_ocr_data_name}/{bin_size}/merged_interactions/{merge_param}/tss_within/{gene_type}/pulled_ocrs_tss_dist.tsv",ocr_metric_files_out_dir_1=ocr_metric_files_out_dir_1,filtered_ocr_data_name=filtered_ocr_data_name,bin_size=bin_size,gene_type=gene_type,merge_param=merge_param)
		
# # # --------------------------------------------------------------
# # # ANNOTATE OCR REGIONS
# # # --------------------------------------------------------------
rule annotate_ocrs:
	input:
		script="/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/data/annotate.regions.R"
	params:
		file_to_annotate=f"{filtered_ocr_file_dir}""/{filtered_ocr_data_name}_tidy.tsv",
		gencode=gencode,
		outdir=filtered_ocr_file_dir,
		gene_type=gene_type,
		file_name="{filtered_ocr_data_name}_annotated_tidy"
	output:
		f"{filtered_ocr_file_dir}""/{gene_type}/{filtered_ocr_data_name}_annotated_tidy.tsv"
	shell:
		"Rscript {input.script} --file_to_annotate {params.file_to_annotate} --gencode {params.gencode} "
		"--outdir {params.outdir} --file_name {params.file_name} --gene_type {params.gene_type}"


# # # --------------------------------------------------------------
# # # ANNOTATE HIC BINS
# # # --------------------------------------------------------------
rule annotate_hic:
	input:
		script="/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/data/annotate.hic.bins.R"
	params:
		interaction_file=f"{interaction_dir}""/{bin_size}/{merge_param}/DICE_interactions_merged_tidy.tsv",
		gencode=gencode,
		outdir=interaction_dir,
		prefix=merge_param,
		gene_type=gene_type
	output:
		f"{interaction_dir}/""{bin_size}/{merge_param}/{gene_type}/merged_interactions_tss_within.tsv"
	shell:
		"Rscript {input.script} --gencode {params.gencode} --interaction_file {params.interaction_file} --outdir {params.outdir} --is_merged merged --prefix {params.prefix} --gene_type {params.gene_type}"

# # # --------------------------------------------------------------
# # # ANNOTATE UNMERGED HIC BINS
# # # --------------------------------------------------------------
rule annotate_unmerged_hic:
	input:
		script="/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/data/annotate.hic.bins.R"
	params:
		interaction_file=f"{interaction_dir}""/{bin_size}/DICE_interactions_tidy.tsv",
		gencode=gencode,
		outdir=interaction_dir,
		prefix="none",
		gene_type=gene_type
	output:
		f"{interaction_dir}""/{bin_size}/{gene_type}/interactions_tss_within.tsv"
	shell:
		"Rscript {input.script} --gencode {params.gencode} --interaction_file {params.interaction_file} --outdir {params.outdir} --is_merged no --prefix {params.prefix} --gene_type {params.gene_type}"

# # # --------------------------------------------------------------------------
# # # FILTER MERGED INTERACTIONS AND GET UNMERGED INTERACTIONS FOR MISSING GENES
# # # --------------------------------------------------------------------------
rule process_interaction_data:
	input:
		interaction_file=f"{interaction_dir}""/{bin_size}/{gene_type}/interactions_tss_within.tsv",
		merged_interaction_file=f"{interaction_dir}/""{bin_size}/{merge_param}/{gene_type}/merged_interactions_tss_within.tsv"
	params:
		script="/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/data/combine.interaction.sets.R",
		outdir=f"{interaction_dir}/""{bin_size}/{merge_param}/{gene_type}/"
	output:
		f"{interaction_dir}/""{bin_size}/{merge_param}/{gene_type}/interactions_tss_within_tidy.tsv"
	shell:
		"Rscript {params.script} --interaction_file {input.interaction_file} --merged_interaction_file {input.merged_interaction_file} --outdir {params.outdir}"


# # # --------------------------------------------------------------
# # # LINK OCRS TO GENES
# # # --------------------------------------------------------------
rule link_ocrs:
	input:
		filtered_ocr_file=f"{filtered_ocr_file_dir}""/{gene_type}/{filtered_ocr_data_name}_annotated_tidy.tsv",
		annotated_interaction_file=f"{interaction_dir}/""{bin_size}/{merge_param}/{gene_type}/interactions_tss_within_tidy.tsv",
		gene_file=f"{filtered_gene_file_dir}""/ranked_{gene_type}_list.tsv"
	params:
		script="/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/data/pull.ocrs.R",
		gencode=gencode,
		outdir=f"{pulled_ocrs_out_dir_1}/""{filtered_ocr_data_name}/{bin_size}/merged_interactions/{merge_param}/tss_within/{gene_type}/",
		is_merged="combined"
	output:
		pulled_ocrs=f"{pulled_ocrs_out_dir_1}/""{filtered_ocr_data_name}/{bin_size}/merged_interactions/{merge_param}/tss_within/{gene_type}/pulled_ocrs_master_by_gene.tsv"
	shell:
		"Rscript {params.script} --filtered_ocr_file {input.filtered_ocr_file} --interaction_file {input.annotated_interaction_file} "
		"--gene_file {input.gene_file} --gencode {params.gencode} "
		"--outdir {params.outdir} --is_merged {params.is_merged}"

# # # --------------------------------------------------------------
# # # REMOVE GENES WITH NO OCRS
# # # --------------------------------------------------------------
rule remove_genes:
	input:
		script="/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/data/remove.gene.no.ocr.R",
		pulled_ocr_file=f"{pulled_ocrs_out_dir_1}/""{filtered_ocr_data_name}/{bin_size}/merged_interactions/{merge_param}/tss_within/{gene_type}/pulled_ocrs_master_by_gene.tsv"
	params:
		gene_file=gene_list,
		outdir=f"{filtered_gene_file_dir}/",
		new_output=f"{filtered_gene_file_dir}/""ranked_{gene_type}_list_filtered.tsv"
	output:
		f"{filtered_gene_file_dir}/""ranked_{gene_type}_list_filtered+{filtered_ocr_data_name}+{bin_size}+{merge_param}.tsv"
	shell:
		"Rscript {input.script} --gene_file {params.gene_file} --pulled_ocr_file {input.pulled_ocr_file} "
		"--outdir {params.outdir}"

# # # --------------------------------------------------------------
# # # CREATE METRIC FILES
# # # --------------------------------------------------------------
rule create_metric_files:
	input:
		linked_ocrs=f"{pulled_ocrs_out_dir_1}/""{filtered_ocr_data_name}/{bin_size}/merged_interactions/{merge_param}/tss_within/{gene_type}/pulled_ocrs_master_by_gene.tsv"
	params:
		script="/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/data/create.metric.files.R",
		gencode=gencode,
		gene_file=gene_list,
		outdir=f"{pulled_ocrs_out_dir_1}/""{filtered_ocr_data_name}/{bin_size}/merged_interactions/{merge_param}/tss_within/{gene_type}/"
	output:
		pulled_ocrs=f"{ocr_metric_files_out_dir_1}/""{filtered_ocr_data_name}/{bin_size}/merged_interactions/{merge_param}/tss_within/{gene_type}/pulled_ocrs_tss_dist.tsv"
	shell:
		"Rscript {params.script} --linked_ocrs {input.linked_ocrs} --gene_file {params.gene_file} "
		"--gencode {params.gencode} --outdir {params.outdir}"