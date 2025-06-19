# # Author
# # Nicole Knoetze
import pandas as pd
import os

### DIR VARIABLES ###
seq_reads_dir=config['SEQ_READS_DIR']
seq_reads_outdir=config['SEQ_READS_OUTDIR']
outdir=config['OUTPUT_DIR']
scratch_outdir=config['SCRATCH_OUTPUT_DIR']
### FILES ###
oligo_library_index=config['OLIGO_LIBRARY_INDEX']
dummy_seq_index=config['DUMMY_SEQ_INDEX']

### VARIABLES ###
cutadapt_env=config['CUTADAPT_ENV']
mapping_env=config['MAPPING_ENV']
datamash_env=config['DATAMASH_ENV']
umi_collapsing_env=config['UMI_COLLAPSING_ENV']
cell_line=config['CELL_LINE']
replicate=config['REPLICATE']
# --------------------------------------------------------------
# SET RULE ALL
# --------------------------------------------------------------
rule all:
    input:
        expand(f"{outdir}""/02-candidate_counts/{cell_line}_RNA{replicate}.mismatch0.counts.0.tsv",outdir=outdir,cell_line=cell_line,replicate=replicate)
# --------------------------------------------------------------
# EXTRACT UMI
# --------------------------------------------------------------
rule trim_reads:
    input:
        R1=f"{seq_reads_dir}""/{cell_line}_RNA{replicate}_1.fastq.gz",
        R2=f"{seq_reads_dir}""/{cell_line}_RNA{replicate}_2.fastq.gz"
    conda: cutadapt_env
    log: f"{seq_reads_outdir}""/01-trimmed-fastq/logs/{cell_line}_RNA{replicate}.trimlog.txt"
    params:
        adapter_fwd="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        adapter_rev="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        quality_cutoff=10,
        min_length=80
    threads: 24
    output:
        R1_trimmed=f"{seq_reads_outdir}""/01-trimmed-fastq/{cell_line}_RNA{replicate}_trimmed_1.fastq.gz",
        R2_trimmed=f"{seq_reads_outdir}""/01-trimmed-fastq/{cell_line}_RNA{replicate}_trimmed_2.fastq.gz"
    shell:
        "cutadapt -a {params.adapter_fwd} -A {params.adapter_rev} -q {params.quality_cutoff} --minimum-length {params.min_length} -j {threads} -o {output.R1_trimmed} -p {output.R2_trimmed} {input.R1} {input.R2} 1> {log}"

# --------------------------------------------------------------
# IDENTIFY DUMMY CONTAMINATION
# --------------------------------------------------------------
#map the reads to the dummy sequence
rule identify_contamination:
    input:
        R1_trimmed=f"{seq_reads_outdir}""/01-trimmed-fastq/{cell_line}_RNA{replicate}_trimmed_1.fastq.gz",
        R2_trimmed=f"{seq_reads_outdir}""/01-trimmed-fastq/{cell_line}_RNA{replicate}_trimmed_2.fastq.gz"
    params:
        dummy_seq_index=dummy_seq_index
    threads: 24
    conda: mapping_env
    log: f"{outdir}/""01-mapped_reads/logs/contamination_{cell_line}_RNA{replicate}.mapping.log"
    output:
        f"{outdir}/""01-mapped_reads/contamination/{cell_line}_RNA{replicate}.bam"
    shell:
        "(bowtie2 -p {threads} -x {params.dummy_seq_index} -1 {input.R1_trimmed} -2 {input.R2_trimmed} | samtools view -@ {threads} -F 780 -f 2 -Sb > {output}) 2> {log}"


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
        R1_trimmed=f"{seq_reads_outdir}""/01-trimmed-fastq/{cell_line}_RNA{replicate}_trimmed_1.fastq.gz",
        R2_trimmed=f"{seq_reads_outdir}""/01-trimmed-fastq/{cell_line}_RNA{replicate}_trimmed_2.fastq.gz"
    params:
        oligo_library_index=oligo_library_index
    threads: 24
    conda: mapping_env
    log: f"{outdir}/""01-mapped_reads/logs/{cell_line}_RNA{replicate}.mapping.log"
    output:
        f"{outdir}/""01-mapped_reads/{cell_line}_RNA{replicate}.bam"
    shell:
        "(bowtie2 -p {threads} -x {params.oligo_library_index} -1 {input.R1_trimmed} -2 {input.R2_trimmed} | samtools view -@ {threads} -F 780 -f 2 -Sb > {output}) 2> {log}"

#--------------------------------------------------------------
# FILTER ALIGNMENTS (RNA AND DNA)
#--------------------------------------------------------------
#need to do filtering for each mismatch independently, then merge the files. 
# Eg merge results for mismatch 0 and mismatch 1 is the MM_threshold is set to 1
rule filter_alignments:
    input:
        f"{outdir}/""01-mapped_reads/{cell_line}_RNA{replicate}.bam"
    conda: mapping_env
    params:
        MM_threshold=0
    output:
        temp(f"{scratch_outdir}""/01-mapped_reads/temp.{cell_line}_RNA{replicate}.mismatch0.bam")
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
#need fix bam file after filtering to update pairing information       
rule fixmates:
    input:
        f"{scratch_outdir}""/01-mapped_reads/temp.{cell_line}_RNA{replicate}.mismatch0.bam"
    conda: mapping_env
    threads: 24
    output:
        temp(f"{scratch_outdir}""/01-mapped_reads/{cell_line}_RNA{replicate}.mismatch0.bam")
    shell:
        "samtools sort -@ {threads} -n {input} | samtools fixmate -@ {threads} - {output}"

#FINAL BAM FILE *****
rule filter_alignments2:
#keep read mapped in proper pair: this makes sure BOTH reads had a mismatch of 0 after bamtools filtering (-f 2)
#only keep the fwd read to reduce file size for processing
    input:
        f"{scratch_outdir}""/01-mapped_reads/{cell_line}_RNA{replicate}.mismatch0.bam"
    conda: mapping_env
    threads: 24
    output:
        f"{outdir}""/01-mapped_reads/{cell_line}_RNA{replicate}.mismatch0.bam"
    shell:
        "samtools view -@ {threads} -f 66 -b {input} > {output}"

#--------------------------------------------------------------
# CREATE UMI BED FILE
#--------------------------------------------------------------
#only keep read name_UMI, candidate, for processing to reduce file size.
#Exclude reads where the UMI has at least one N.
rule bam_to_tsv:
    input:
        f"{outdir}""/01-mapped_reads/{cell_line}_RNA{replicate}.mismatch0.bam"
    conda: mapping_env
    threads: 24
    output:
        temp(f"{outdir}""/01-mapped_reads/{cell_line}_RNA{replicate}.mismatch0.tsv")
    shell:
        """samtools view -@ {threads} {input} | cut -f 1,3 | awk '{{if (substr($1, length($1)-9, 10) !~ /N/) print}}'> {output}"""

rule create_UMI_bed:
    input:
        f"{outdir}""/01-mapped_reads/{cell_line}_RNA{replicate}.mismatch0.tsv"
    params:
        script="/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/umi_starrseq/create_UMI_bed.py"
    output:
        f"{scratch_outdir}""/02-UMI_counts/{cell_line}_RNA{replicate}.mismatch0.UMI.bed"
    shell:
        "python {params.script} --SAM_file {input} --outfile {output}"

#--------------------------------------------------------------
# GET UMI COUNTS
#--------------------------------------------------------------
rule get_counts:
    input:
        f"{scratch_outdir}""/02-UMI_counts/{cell_line}_RNA{replicate}.mismatch0.UMI.bed"
    conda: umi_collapsing_env
    params: 
        script="/projects/nknoetze_prj/ocr_prj/src/tcell_ocr_prj/umi_starrseq/collapse_UMI.py",
        umi_MM=0
    output:
        f"{outdir}""/02-candidate_counts/{cell_line}_RNA{replicate}.mismatch0.counts.0.tsv"
    shell:
        "python {params.script} --UMI_bedfile {input} --mm_threshold {params.umi_MM} --outfile {output}"