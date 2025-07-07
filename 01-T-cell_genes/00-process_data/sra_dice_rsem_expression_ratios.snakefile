#PACKAGES
import pandas as pd

### OUT DIR VARIABLES ###
out_dir=config['OUTDIR']
tmp_dir=config['TEMP_OUTDIR']

star_out_dir=config['STAR']
rsem_out_dir=config['RSEM']
rsem_quants_out_dir=config['RSEM_FINAL']
matrices_out_dir=config['MATRICES']
corrected_matrices_out_dir=config['CORRECTED_MATRICES']
normalised_counts_out_dir=config['NORMALISED_COUNTS']

star_out_dir_1=f"{tmp_dir}/{star_out_dir}"
rsem_out_dir_2=f"{tmp_dir}/{rsem_out_dir}"

rsem_final_out_dir_1=f"{out_dir}/{rsem_quants_out_dir}"
raw_counts_out_dir_2=f"{out_dir}/{matrices_out_dir}"
corrected_matrices_out_dir_3=f"{out_dir}/{corrected_matrices_out_dir}"
normalised_counts_out_dir_4=f"{out_dir}/{normalised_counts_out_dir}"

### VARIABLES ###
trimmed_fastqs=config['TRIMMED_FASTQ']
star_index=config['STAR_INDEX']
rsem_index=config['RSEM_INDEX']
master_metadata_file=config['MASTER_METADATA']
gencode=config['GENCODE']
dice_rsem_files=config['DICE_RSEM_FILES']
conda_env_file=config['CONDA_CONFIG']
counts_file_prefix=config['COUNTS_FILE_PREFIX']
#--- GET SRR SAMPLE IDS FOR THE FASTQS
try:
    SAMPLE_DF = pd.read_table(master_metadata_file)
except FileNotFoundError:
    SAMPLE_DF = pd.DataFrame()
DICE_DF=SAMPLE_DF[SAMPLE_DF.authors=="DICE"]
SRA_DF=SAMPLE_DF[SAMPLE_DF.authors!="DICE"]
dice_sample_id_list=(DICE_DF['sample_id']).tolist()
sra_sample_id_list=(SRA_DF['sample_id']).tolist()

rule all:
	input:
		expand(f"{star_out_dir_1}""/{sra_sample}.Aligned.toTranscriptome.out.bam",sra_sample=sra_sample_id_list),
        expand(f"{rsem_final_out_dir_1}""/{sra_sample}/{sra_sample}.genes.results", sra_sample=sra_sample_id_list),
        expand(f"{rsem_final_out_dir_1}""/{sra_sample}/{sra_sample}.isoforms.results",sra_sample=sra_sample_id_list),
        f"{rsem_final_out_dir_1}""/dice_sra_effective_lengths.tsv",
        f"{raw_counts_out_dir_2}/sra.dice.counts.raw.tsv.gz",
        f"{corrected_matrices_out_dir_3}/sra.dice.counts.batchcorr.tsv.gz",
        expand(f"{normalised_counts_out_dir_4}/""sra.dice.counts.{counts_file_prefix}.sf.adj.tsv.gz",counts_file_prefix=counts_file_prefix),
        expand(f"{normalised_counts_out_dir_4}/""sra.dice.counts.{counts_file_prefix}.sf.genel.adj.tsv.gz",counts_file_prefix=counts_file_prefix)

# # # --------------------------------------------------------------
# # # ALIGN READS TO GENOME USING STAR
# # # --------------------------------------------------------------

rule run_star:
    input:
        script="gtex_runstar.py",
        fastq1=f"{trimmed_fastqs}""/{sra_sample}_1.trimmed.fastq.gz",
        index=star_index
    params:
        overhang="49",
        outdir=f"{star_out_dir_1}""/{sra_sample}"
    priority: 10
    output:
        f"{star_out_dir_1}""/{sra_sample}/{sra_sample}.Aligned.toTranscriptome.out.bam"
    run:
        shell("python {input.script} {input.index} {input.fastq1} "\
        "{wildcards.sra_sample} -o {params.outdir} --sjdbOverhang {params.overhang}")

rule move_star:
    input:
        bam=f"{star_out_dir_1}""/{sra_sample}/{sra_sample}.Aligned.toTranscriptome.out.bam"
    params:
        star_dir=f"{star_out_dir_1}""/{sra_sample}"
    priority: 9
    output:
        f"{star_out_dir_1}""/{sra_sample}.Aligned.toTranscriptome.out.bam"
    shell:
        "mv {input.bam} {output}; rm -rf {params.star_dir}"

# # # --------------------------------------------------------------
# # # GENERATE GENE COUNTS
# # # --------------------------------------------------------------

rule run_rsem:
    input:
        script="gtex_runrsem.py",
        bam=f"{star_out_dir_1}""/{sra_sample}.Aligned.toTranscriptome.out.bam",
        reference_index=rsem_index
    params:
        outdir=f"{rsem_out_dir_2}""/{sra_sample}/"
    priority: 8
    output:
        f"{rsem_out_dir_2}""/{sra_sample}/{sra_sample}.genes.results",
        f"{rsem_out_dir_2}""/{sra_sample}/{sra_sample}.isoforms.results"
    run:
        shell("python {input.script} {input.reference_index} {input.bam} {wildcards.sra_sample} -o {params.outdir} ")

rule move_rsem:
    input:
        gene_files=f"{rsem_out_dir_2}""/{sra_sample}/{sra_sample}.genes.results",
        isoform_files=f"{rsem_out_dir_2}""/{sra_sample}/{sra_sample}.isoforms.results"
    priority: 7
    output:
        gene_results=f"{rsem_final_out_dir_1}""/{sra_sample}/{sra_sample}.genes.results",
        isoform_results=f"{rsem_final_out_dir_1}""/{sra_sample}/{sra_sample}.isoforms.results"
    shell:
        "mv {input.gene_files} {output.gene_results}; mv {input.isoform_files} {output.isoform_results}"

# # # --------------------------------------------------------------
# # # GET EFFECTIVE GENE LENGTHS FOR EACH SAMPLE
# # # --------------------------------------------------------------

rule effective_lengths_sra:
    input:
        gene_results=f"{rsem_final_out_dir_1}""/{sra_sample}/{sra_sample}.genes.results"
    output:
        temp(f"{rsem_out_dir_2}""/{sra_sample}.effective.gene.lengths.tsv")
    priority: 6
    shell:
        """tail -n+2 {input.gene_results} | awk -v sample_id={wildcards.sra_sample} -v FS="\\t" -v OFS="\\t" "{{print \$1,\$4,sample_id}}" > {output}"""

rule effective_lengths_dice:
    input:
        dice_gene_results=f"{dice_rsem_files}""/{dice_sample}/{dice_sample}.genes.results"
    output:
        temp(f"{rsem_out_dir_2}""/{dice_sample}.effective.gene.lengths.tsv")
    priority: 6
    shell:
        """tail -n+2 {input.dice_gene_results} | awk -v sample_id={wildcards.dice_sample} -v FS="\\t" -v OFS="\\t" "{{print \$1,\$4,sample_id}}" > {output}"""

rule get_effective_lengths:
    input:
        dice_effective_lengths=expand(f"{rsem_out_dir_2}""/{dice_sample}.effective.gene.lengths.tsv",dice_sample=dice_sample_id_list),
        sra_effective_lengths=expand(f"{rsem_out_dir_2}""/{sra_sample}.effective.gene.lengths.tsv",sra_sample=sra_sample_id_list)
    output:
        merged_effective_lengths=f"{rsem_final_out_dir_1}""/dice_sra_effective_lengths.tsv"
    params:
        effective_lengths_dice_sra=f"{rsem_out_dir_2}""/*.effective.gene.lengths.tsv"
    priority: 5
    shell:
        "cat {params.effective_lengths_dice_sra} > {output.merged_effective_lengths}"

# # # --------------------------------------------------------------
# # # GENERATE THE RAW COUNTS MATRIX
# # # --------------------------------------------------------------

rule generate_raw_counts_matrix:
    input:
        script='sra.and.dice.generate.matrix.R',
        gencode=gencode,
        metadata_file=master_metadata_file,
        dice_rsem_files=dice_rsem_files,
        sra_rsem_files=rsem_final_out_dir_1
    params:
        outdir=f"{raw_counts_out_dir_2}/",
        log=f"{raw_counts_out_dir_2}/sra.dice.counts.log"
    priority: 4
    output:
        outfile=f"{raw_counts_out_dir_2}""/sra.dice.counts.raw.tsv.gz"
    conda: conda_env_file
    shell:
        "Rscript {input.script} --outdir {params.outdir} --gencode {input.gencode} --DICE_file_path {input.dice_rsem_files} "
        "--SRA_file_path {input.sra_rsem_files} --metadata_file {input.metadata_file} 2>&1 | tee {params.log}"

# # # --------------------------------------------------------------
# # # PERFORM SF AND GENE LENGTH NORMALISATION
# # # --------------------------------------------------------------

rule sf_normalisation:
    input:
        script="sf.matrix.normalisation.R",
        counts1=f"{raw_counts_out_dir_2}/sra.dice.counts.raw.tsv.gz",
        counts2=f"{corrected_matrices_out_dir_3}/""sra.dice.counts.batchcorr.tsv.gz",
        metadata=master_metadata_file
    params:
        counts_matrix=f"{out_dir}""/*/sra.dice.counts.{counts_file_prefix}.tsv.gz",
        outdir=f"{normalised_counts_out_dir_4}/",
        log=f"{normalised_counts_out_dir_4}""/logs/{counts_file_prefix}.sf.normalisation.log"
    priority: 2
    output:
        f"{normalised_counts_out_dir_4}/""sra.dice.counts.{counts_file_prefix}.sf.adj.tsv.gz"
    shell:
        "Rscript {input.script} --count_matrix {params.counts_matrix} --metadata_file_path {input.metadata} --outdir {params.outdir} "
        "--prefix {wildcards.counts_file_prefix} 2>&1 | tee {params.log}"

rule genelength_normalisation:
    input:
        script="effective.gene.length.normalisation.py",
        counts_matrix=f"{normalised_counts_out_dir_4}/""sra.dice.counts.{counts_file_prefix}.sf.adj.tsv.gz",
        metadata=master_metadata_file,
        effective_lengths=f"{rsem_final_out_dir_1}""/dice_sra_effective_lengths.tsv"
    params:
        outdir=f"{normalised_counts_out_dir_4}/",
        log=f"{normalised_counts_out_dir_4}""/logs/{counts_file_prefix}.genel.normalisation.log"
    priority: 2
    output:
        f"{normalised_counts_out_dir_4}/""sra.dice.counts.{counts_file_prefix}.sf.genel.adj.tsv.gz"
    shell:
        "python {input.script} --count_matrix {input.counts_matrix} --metadata_file_path {input.metadata} --outdir {params.outdir} "
        "--effective_lengths {input.effective_lengths} --prefix {wildcards.counts_file_prefix} 2>&1 | tee {params.log}"