### OUT DIR VARIABLES ###
outdir=config['OUTDIR']

### VARIABLES ###
conda_env=config['CONDA_ENV']
gencode=config['GENCODE']
gene_expr_file=config['GENE_EXPR_FILE']
target_gene_list=config['TARGET_GENES_FILE']
expressed_tcell_tfs=config['EXPRESSED_TCELL_TFS']
n_tf=config['N_TF']
feature_file_name=config['FEATURE_FILE_NAME']
ocr_type=config['OCR_TYPE']
motif_type=config['MOTIF_TYPE']
motif_comparisons=config['MOTIF_COMPARISONS']
n_cores=config['N_CORES']
metric=config['METRIC']

# # # --------------------------------------------------------------
# # # SET RULE ALL
# # # --------------------------------------------------------------

rule all:
	input:
		expand(f"{outdir}""/{feature_file_name}.{ocr_type}.filtered.full.{metric}.nonoverlapping.{motif_comparisons}.ranked.tsv",
		outdir=outdir,feature_file_name=feature_file_name,ocr_type=ocr_type,motif_comparisons=motif_comparisons,metric=metric)

# # # --------------------------------------------------------------
# # # FILTER RESULTS FOR SIGNIFICANT RESULTS, RELEVANT RESULTS AND EXPRESSED TFS
# # # --------------------------------------------------------------
rule filter_results:
	input:
		filtered_features=f"{outdir}""/{feature_file_name}.tsv"
	params:
		script="03_expressed_tfs.R",
		tcell_tfs=expressed_tcell_tfs,
		n_tfs=n_tf,
		motif_type="{motif_type}",
		ocr_type="{ocr_type}"
	output:
		f"{outdir}""/{feature_file_name}.{ocr_type}.{motif_type}.filtered.tsv"
	shell:
		"Rscript {params.script} --feature_file {input.filtered_features} --ocr_type {params.ocr_type} "
		"--tcell_tfs {params.tcell_tfs} --motif_type {params.motif_type} --n_tfs {params.n_tfs}"

# --------------------------------------------------------------
# REMOVE OVERLAPPING TFBS FOR A GIVEN PAIR/TRIPLET
# --------------------------------------------------------------
rule get_nonoverlapping_sites:
	input:
		filtered_features=f"{outdir}""/{feature_file_name}.{ocr_type}.{motif_type}.filtered.tsv",
		tfbs_pos=f"{outdir}""/{feature_file_name}.tsv.fimo_comotif_geneCount+fimo_comotif_geneCount+{ocr_type}.RAW.featureLocations.tsv"
	params:
		script="04_nonoverlapping_sites.R",
		n_tfs=n_tf,
		ocr_type="{ocr_type}",
		n_cores=n_cores
	priority: 7
	log:
		f"{outdir}""/logs/{feature_file_name}.{ocr_type}.{motif_type}.filtered.full.nonoverlapping.log"
	output:
		f"{outdir}""/{feature_file_name}.{ocr_type}.{motif_type}.filtered.full.nonoverlapping.tsv"
	shell:
		"Rscript {params.script} --feature_file {input.filtered_features} --tfbs_pos {input.tfbs_pos} "
		"--ocr_type {params.ocr_type} --n_tfs {params.n_tfs} --n_cores {params.n_cores} 2>&1 | tee {log}"

# --------------------------------------------------------------
# COMBINE FILES
# --------------------------------------------------------------
rule combine_results:
	input:
		expand(f"{outdir}""/{feature_file_name}.{{ocr_type}}.{motif_type}.filtered.full.nonoverlapping.tsv",outdir=outdir,feature_file_name=feature_file_name,motif_type=motif_type)
	params:
		ocr_type="{ocr_type}"
	priority: 7
	output:
		temp=temp(f"{outdir}""/{feature_file_name}.{ocr_type}.full.tempfile.tsv"),
		outfile=f"{outdir}""/{feature_file_name}.{ocr_type}.filtered.full.nonoverlapping.allsites.tsv"
	#write header to file, take only the first line (since we need one header). Add remaining lines to file.
	shell:
		"""awk 'FNR==1' {input} > {output.temp} ; awk 'FNR==1' {output.temp} > {output.outfile}; awk 'FNR>1' {input} >> {output.outfile} """

# # # # --------------------------------------------------------------
# # # # RANK FEATURES
# # # # --------------------------------------------------------------

rule rank_features:
	input:
		nonoverlapping_features=f"{outdir}""/{feature_file_name}.{ocr_type}.filtered.full.nonoverlapping.allsites.tsv"
	params:
		script="05_rank_tfs.R",
		ocr_type="{ocr_type}"
	priority: 5
	output:
		f"{outdir}""/{feature_file_name}.{ocr_type}.filtered.full.nonoverlapping.{motif_comparisons}.ranked.tsv"
	shell:
		"Rscript {params.script} --feature_file {input.nonoverlapping_features}"
