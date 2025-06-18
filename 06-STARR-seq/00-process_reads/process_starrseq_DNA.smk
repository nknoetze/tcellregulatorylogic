# # Author
# # Nicole Knoetze
import pandas as pd
import os

### DIR VARIABLES ###
seq_reads_dir=config['SEQ_READS_DIR']
outdir=config['OUTPUT_DIR']
scratch_outdir=config['SCRATCH_OUTPUT_DIR']

### FILES ###
oligo_library_index=config['OLIGO_LIBRARY_INDEX']

### VARIABLES ###
cutadapt_env=config['CUTADAPT_ENV']
umitools_env=config['UMITOOLS_ENV']
mapping_env=config['MAPPING_ENV']
umi_collapsing_env=config['UMI_COLLAPSING_ENV']
datamash_env=config['DATAMASH_ENV']
library_names=config['LIBRARY_NAMES']
MM_thresholds=config['MM_THRESHOLDS']
umi_MMs=config['UMI_MM']
seed=config['SEED']
method=config['METHOD']
# --------------------------------------------------------------
# SET RULE ALL
# --------------------------------------------------------------
rule all:
    input:
        #expand(f"{scratch_outdir}""/03-candidate_counts/subsampled{seed}_{method}_DNA{library_name}.mismatch0.counts.0.tsv",outdir=outdir,library_name=library_names,seed=seed,method=method)
        expand(f"{outdir}""/03-candidate_counts/{method}_DNA{library_name}.counts.tsv",library_name=library_names,method=method)
    
# --------------------------------------------------------------
# EXTRACT UMI
# --------------------------------------------------------------
rule extract_UMI:
    input:
        R1=f"{seq_reads_dir}""/{library_name}_R1_001.fastq.gz",
        R2=f"{seq_reads_dir}""/{library_name}_R2_001.fastq.gz"
    conda: umitools_env
    log: f"{outdir}""/01-trimmed-fastq/logs/DNA{library_name}_umi_extract_log.txt"
    output:
        R1_umi=temp(f"{outdir}""/01-trimmed-fastq/DNA{library_name}_umi_R1.fastq.gz"),
        R2_umi=temp(f"{outdir}""/01-trimmed-fastq/DNA{library_name}_umi_R2.fastq.gz")
    shell:
        "umi_tools extract --extract-method=regex --bc-pattern='.*AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC(?P<umi_1>.{{10}}).*' -L {log} -I {input.R1} -S {output.R1_umi} --read2-in={input.R2} --read2-out={output.R2_umi}"

# --------------------------------------------------------------
# TRIM READS
# --------------------------------------------------------------
rule trim_reads:
    input:
        R1=f"{outdir}""/01-trimmed-fastq/DNA{library_name}_umi_R1.fastq.gz",
        R2=f"{outdir}""/01-trimmed-fastq/DNA{library_name}_umi_R2.fastq.gz"
    priority:1
    conda: cutadapt_env
    log: f"{outdir}""/01-trimmed-fastq/logs/DNA{library_name}.trimlog.txt"
    params:
        adapter_fwd="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        adapter_rev="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        quality_cutoff=10,
        min_length=80
    threads: 24
    output:
        R1_trimmed=f"{outdir}""/01-trimmed-fastq/DNA{library_name}_R1_trimmed.fastq.gz",
        R2_trimmed=f"{outdir}""/01-trimmed-fastq/DNA{library_name}_R2_trimmed.fastq.gz"
    shell:
        "cutadapt -a {params.adapter_fwd} -A {params.adapter_rev} -q {params.quality_cutoff} --minimum-length {params.min_length} -j {threads} -o {output.R1_trimmed} -p {output.R2_trimmed} {input.R1} {input.R2} 1> {log}"

# --------------------------------------------------------------
# RUN BOWTIE2
# --------------------------------------------------------------
#map the reads
# default output is sam file. 
# Filter sam file using -F 780 to exclude: read unmapped, mate unmapped, not primary alignment, read fails platform/vendor QC (keeps PCR duplicates)
#-f 2 give proper paired reads only
# final output is a bam file.
rule map_reads:
    input:
        R1=f"{outdir}""/01-trimmed-fastq/DNA{library_name}_R1_trimmed.fastq.gz",
        R2=f"{outdir}""/01-trimmed-fastq/DNA{library_name}_R2_trimmed.fastq.gz"
    params:
        oligo_library_index=oligo_library_index
    threads: 20
    conda: mapping_env
    log: f"{outdir}/""02-mapped_reads/logs/DNA{library_name}.mapping.log"
    output:
        f"{outdir}/""02-mapped_reads/DNA{library_name}.bam"
    shell:
        #"(bowtie2 -p {threads} -x {params.oligo_library_index} -1 {input.R1} -2 {input.R2} | samtools view -@ {threads} -Sb |samtools sort -@ {threads} -o {output} - ) 2> {log}"
        "(bowtie2 -p {threads} -x {params.oligo_library_index} -1 {input.R1} -2 {input.R2} | samtools view -@ {threads} -F 780 -f 2 -Sb > {output}) 2> {log}"

#--------------------------------------------------------------
# FILTER ALIGNMENTS (RNA AND DNA)
#--------------------------------------------------------------
#need to fix bam file after filtering to update pairing information
rule filter_alignments:
    input:
        f"{outdir}/""02-mapped_reads/DNA{library_name}.bam"
    conda: mapping_env
    params:
        MM_threshold=0
    output:
        temp(f"{scratch_outdir}""/02-mapped_reads/temp.DNA{library_name}.mismatch{MM_threshold}.bam")
    shell:
         """
        tmp_files=""
        for MM in {params.MM_threshold}; do 
            for M in $(seq 0 $MM); do
                bamtools filter -tag NM:$M -in {input} -out {output}.tmp_$M.bam 
                tmp_files="$tmp_files {output}.tmp_$M.bam"
            done
        done
        samtools cat -o {output} $tmp_files
        rm $tmp_files
        """  
        
rule fixmates:
    input:
        f"{scratch_outdir}""/02-mapped_reads/temp.DNA{library_name}.mismatch{MM_threshold}.bam"
    conda: mapping_env
    threads: 24
    output:
        temp(f"{scratch_outdir}""/02-mapped_reads/DNA{library_name}.mismatch{MM_threshold}.bam")
    shell:
        "samtools sort -@ {threads} -n {input} | samtools fixmate -@ {threads} - {output}"

#FINAL BAM FILE *****
rule filter_alignments2:
#keep read mapped in proper pair: this makes sure BOTH reads had a mismatch of 0 after bamtools filtering (-f 2)
#only keep the fwd read to reduce file size for processing (-f 64)
#Flag=66
    input:
        f"{scratch_outdir}""/02-mapped_reads/DNA{library_name}.mismatch{MM_threshold}.bam"
    conda: mapping_env
    threads: 24
    output:
        f"{outdir}""/02-mapped_reads/DNA{library_name}.mismatch{MM_threshold}.bam"
    shell:
        "samtools view -@ {threads} -f 66 -b {input} > {output}"

#--------------------------------------------------------------
# CREATE UMI BED FILE
#--------------------------------------------------------------
#only keep the fwd read to reduce file size for processing (-f 64)
#only keep read name_UMI, candidate, and candidate_seq_len for processing to reduce file size.
rule bam_to_tsv:
    input:
        f"{outdir}""/02-mapped_reads/DNA{library_name}.mismatch{MM_threshold}.bam"
    conda: mapping_env
    threads: 24
    output:
        f"{scratch_outdir}""/02-mapped_reads/DNA{library_name}.mismatch{MM_threshold}.tsv"
    shell:
        """samtools view -@ {threads} {input} | cut -f 1,3 | awk '{{if (substr($1, length($1)-9, 10) !~ /N/) print}}'> {output}"""

rule create_UMI_bed:
    input:
        f"{scratch_outdir}""/02-mapped_reads/DNA{library_name}.mismatch{MM_threshold}.tsv"
    params:
        script="/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/umi_starrseq/create_UMI_bed.py"
    output:
        f"{scratch_outdir}""/03-UMI_counts/DNA{library_name}.mismatch{MM_threshold}.UMI.bed"
    shell:
        "python {params.script} --SAM_file {input} --outfile {output}"

#--------------------------------------------------------------
# GET UMI COUNTS
#--------------------------------------------------------------
rule get_counts:
    input:
        f"{scratch_outdir}""/03-UMI_counts/DNA{library_name}.mismatch{MM_threshold}.UMI.bed"
    conda: umi_collapsing_env
    params: 
        script="/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/umi_starrseq/collapse_UMI.py",
        umi_MM="{umi_MM}"
    output:
        f"{outdir}""/03-candidate_counts/DNA{library_name}.mismatch{MM_threshold}.counts.{umi_MM}.tsv"
    shell:
        "python {params.script} --UMI_bedfile {input} --mm_threshold {params.umi_MM} --outfile {output}"

#--------------------------------------------------------------
# SUBSAMPLE COUNTS
#--------------------------------------------------------------
# rule subsample_counts:
#     input:
#         umi_counts_file=f"{outdir}""/03-candidate_counts/DNA{library_name}.mismatch0.counts.0.tsv",
#         wait=expand(f"{outdir}""/03-candidate_counts/DNA{library_name}.mismatch0.counts.0.tsv",library_name=library_names)
#     params: 
#         script="/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/umi_starrseq/subsample_library.R",
#         umi_counts_dir=f"{outdir}""/03-candidate_counts/",
#         seed="{seed}",
#         outdir=f"{scratch_outdir}""/03-candidate_counts/",
#         method="{method}"
#     output:
#         f"{scratch_outdir}""/03-candidate_counts/subsampled{seed}_{method}_DNA{library_name}.mismatch0.counts.0.tsv"
#     shell:
#         "Rscript {params.script} --umi_counts_dir {params.umi_counts_dir} --umi_counts_file {input.umi_counts_file} --seed_value {params.seed} --outdir {params.outdir} --method {params.method}"

#--------------------------------------------------------------
# Average subsampled COUNTS
#--------------------------------------------------------------
# rule average_counts:
#     input:
#         sumsampled_counts=expand(f"{scratch_outdir}""/03-candidate_counts/subsampled{seed}_{method}_DNA{{library_name}}.mismatch0.counts.0.tsv",seed=seed,method=method),
#         wait=expand(f"{scratch_outdir}""/03-candidate_counts/subsampled{seed}_{method}_DNA{library_name}.mismatch0.counts.0.tsv",library_name=library_names,seed=seed,method=method)
#     conda: datamash_env
#     output:
#         f"{outdir}""/03-candidate_counts/{method}_DNA{library_name}.counts.tsv"
#     shell:
#         """cat {input.sumsampled_counts} | grep -v "candidate" | datamash -sg 1 sum 2 | awk 'OFS=FS="\t" {{print $1,$2/10}}' > {output}"""
