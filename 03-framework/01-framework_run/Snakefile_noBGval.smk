## miniP - OCR feature testing
## Wrapper for random resampling framework

## March 10, 2023
## Scott Brown

VERSION = "1.0"
VERSION_DATE = "230310"

## run as:
## snakemake -p -j 32 --use-conda --configfile path/to/config.yaml


import os
import time
import glob
import logging


###############################
##  Define Global Variables  ##
###############################

#configfile: "config.yaml"

## Load Config
REF_GENOME_REFNAME = config["GENOME_REF"]
CHROM_SIZES = config["CHROM_SIZES"]
GENOME_BED_FILE = config["GENOME_BED_FILE"]
EXCL_REGIONS_BED = config["EXCL_REGIONS_BED_FILE"]
OCRMETRICS_DATA_DIR = config["OCRMETRICS_DATA_DIR"]

TARGET_GENES = config["TARGET_GENES"]
BACKGROUND_GENES = config["BACKGROUND_GENES"]
REGIONS = config["REGIONS"]
ALL_GENES = config["ALL_GENES"]

NUM_ITER = config["NUM_ITER"]
CORES = config["CORES"]
MAX_BG_WRITTEN = config["MAX_BACKGROUND_WRITTEN"]

FUNCTIONS = config["FUNCTIONS"]

FIMO_DATA = config["FIMO_DATA"]
EXPRESSED_TF_DATA = config["EXPRESSED_TFS"]

SCRAMBLE_METHOD = config["SCRAMBLE_METHOD"]

OUTPUT = config["OUTPUT"]

CONDA = config["ENV"]

#############################################
##  Define valid metrics for each datatype ##
#############################################

METRICS = {"fimo_geneCount": ["fimo_geneCount"],
           "fimo_geneCount_shuffling": ["fimo_geneCount_shuffling"],
           "fimo_geneCount_scrambling": ["fimo_geneCount_scrambling"],
           "fimo_comotif_geneCount": ["fimo_comotif_geneCount"],
           "fimo_comotif_geneCount_shuffling": ["fimo_comotif_geneCount_shuffling"],
           "fimo_comotif_geneCount_scrambling": ["fimo_comotif_geneCount_scrambling"]}
NUM_MOTIFS=config['NUM_MOTIFS']


################
##  Rule All  ##
################

list_files = [
    "{}.{}+{}+{}.RAW.tsv".format(OUTPUT, d, m, rt)
    for d in [f for f in FUNCTIONS.split(" ") if f in ["fimo_geneCount",
                                                        "fimo_geneCount_shuffling",
                                                        "fimo_geneCount_scrambling",
                                                        "fimo_comotif_geneCount",
                                                        "fimo_comotif_geneCount_shuffling",
                                                        "fimo_comotif_geneCount_scrambling"]]
    for m in METRICS[d]
    for rt in ["all", "promoter", "non-promoter"]
]

list_files_to_filter = [
    "{}.{}+{}+{}.FILTERED.tsv".format(OUTPUT, d, m, rt)
    for d in [x for x in FUNCTIONS.split(" ") if "fimo" in x]
    for m in METRICS[d]
    for rt in ["all", "promoter", "non-promoter"]
]

rule all: 
    input:
        #OUTPUT
        all_res = list_files,
        filtered_res = list_files_to_filter

########################
##  00 Assemble Data  ##
########################

rule cp_expressed_tfs:
    '''copy expressed tf file'''
    input: 
        expressed_tf = EXPRESSED_TF_DATA
    output:
        expressed_tf = os.path.join(OCRMETRICS_DATA_DIR, "fantom5_tf_hg19_expressed.tsv")
    shell:
        "cp {input.expressed_tf} {output.expressed_tf}"

rule cp_fimo:
    '''Copy Fimo Motif Data'''
    input:
        fimo = FIMO_DATA
    output:
        fimo = os.path.join(OCRMETRICS_DATA_DIR, "fimo_motifs_raw.txt")
    shell:
        "cp {input.fimo} {output.fimo};"

rule modify_fimo:
    '''Modify Meme file to separate TF name and motif id'''
    input:
        fimo = rules.cp_fimo.output.fimo
    output:
        fimo_modified = os.path.join(OCRMETRICS_DATA_DIR, "fimo_motifs.txt")
    shell:
        """awk '{{if($1=="MOTIF"){{sub(/.*\./, "", $3)}}; print}}' {input.fimo} > {output.fimo_modified}"""


rule expressed_fimo_id:
    '''Get file with motif Ids for those that have expressed TFs'''
    input:  
        #fimo = rules.cp_fimo.output.fimo,
        fimo = rules.modify_fimo.output.fimo_modified,
        expressed_tf = rules.cp_expressed_tfs.output.expressed_tf
    output:
        expressed_fimo_ids = os.path.join(OCRMETRICS_DATA_DIR, "fantom5_tf_hg19_expressed_ids.tsv")
    shell:
    #convert motif tfs to upper so we can filter by expressed TFs
    # modify the motif_id by replacing . with .\ to make it compatible with sed.
        """grep "MOTIF" {input.fimo} | cut -d " " -f 2-3 | tr [:lower:] [:upper:] | grep -f {input.expressed_tf} | cut -d " " -f 1 | sed 's/\./\\./g' > {output.expressed_fimo_ids}"""

rule filter_fimo:
    '''filter Fimo meme file for expressed tfs'''
    input:
        expressed_fimo_ids = rules.expressed_fimo_id.output.expressed_fimo_ids,
        fimo=os.path.join(OCRMETRICS_DATA_DIR, "fimo_motifs.txt")
    output:
        filtered_fimo=os.path.join(OCRMETRICS_DATA_DIR, "fimo_motifs_expressed.txt")
    shell:
    #print the header (first 9 lines) and filter the meme file based on the motif ids file. print up until the URL line for each motif id, echo to add new line to keep format
        """head -n 9 {input.fimo} > {output}; parallel -a {input.expressed_fimo_ids} 'sed -n "/{{}}/,/^URL/ p" {input.fimo} && echo' >> {output.filtered_fimo}"""

########################
##  01 Run Framework  ##
########################

rule run_framework:
    '''Run the random resampling for the selected functions and data'''
    input:
        target_genes = TARGET_GENES,
        background_genes = BACKGROUND_GENES,
        regions = REGIONS,
        fimo=rules.filter_fimo.output.filtered_fimo   ## EXPRESSION-FILTERED
    output:
        outfile = OUTPUT
    log: 
        log = "{}.log".format(OUTPUT)
    threads: CORES
    conda: CONDA
    params:
        framework_script = os.path.join(workflow.basedir,"test_features_noDistribution.py"),
        ocrmetrics_data_dir = OCRMETRICS_DATA_DIR,
        ref_fasta = REF_GENOME_REFNAME,
        ref_bed = GENOME_BED_FILE,
        chrom_sizes = CHROM_SIZES,
        all_genes = ALL_GENES,
        num_iter = NUM_ITER,
        max_bg = MAX_BG_WRITTEN,
        functions = FUNCTIONS,
        num_motifs = NUM_MOTIFS,
        excl_regions_bed = EXCL_REGIONS_BED,
        scramble_method = SCRAMBLE_METHOD
    shell:
        """python {params.framework_script} \
        --ocrmetrics_data_dir {params.ocrmetrics_data_dir} \
        --refgenome_fasta {params.ref_fasta} \
        --chrom_sizes_file {params.chrom_sizes} \
        --refgenome_bed {params.ref_bed} \
        --excl_regions_bed {params.excl_regions_bed} \
        --target_genes {input.target_genes} \
        --background_genes {input.background_genes} \
        --region_file {input.regions} \
        --all_genes {params.all_genes} \
        --num_iterations {params.num_iter} \
        --output {output.outfile} \
        --cores {threads} \
        --max_bg_written {params.max_bg} \
        --functions {params.functions} \
        --num_motifs {params.num_motifs} \
        --scramble_method {params.scramble_method} \
        2>&1 | tee {log};"""


################################
##  02 Extract Feature Scores ##
################################

rule extract_feature_scores:
    '''Extract positions of significant features'''
    input:
        resultfile = rules.run_framework.output.outfile,
        target_genes = TARGET_GENES,
        regions = REGIONS,
        fimo = rules.modify_fimo.output.fimo_modified,
    output:
        outfile = "{}.{{datatype}}+{{metric}}+{{region_type}}.RAW.tsv".format(OUTPUT),
        feat_outfile = "{}.{{datatype}}+{{metric}}+{{region_type}}.RAW.featureLocations.tsv".format(OUTPUT)
    conda: CONDA
    params:
        extraction_script = os.path.join(workflow.basedir, "generatePerBaseFeatureScores.py"),
        pthresh = 0.05,
        fcthresh = 2,
        ocrmetrics_data_dir = OCRMETRICS_DATA_DIR,
        ref_fasta = REF_GENOME_REFNAME,
    shell:
        """python {params.extraction_script} \
        --target_genes {input.target_genes} \
        --region_file {input.regions} \
        --result_file {input.resultfile} \
        --pthresh {params.pthresh} --fcthresh {params.fcthresh} \
        --ocrmetrics_data_dir {params.ocrmetrics_data_dir} \
        --datatype {wildcards.datatype} \
        --metric {wildcards.metric} \
        --region_set {wildcards.region_type} \
        --refgenome_fasta {params.ref_fasta} \
        --output {output.outfile} \
        --output_features {output.feat_outfile}"""
    

rule extract_feature_scores_tf_filtered:
    '''Extract positions of significant features'''
    input:
        resultfile = rules.run_framework.output.outfile,
        target_genes = TARGET_GENES,
        regions = REGIONS,
        fimo = rules.modify_fimo.output.fimo_modified,
        expressed_tf = rules.cp_expressed_tfs.output.expressed_tf
    output:
        outfile = "{}.{{datatype}}+{{metric}}+{{region_type}}.FILTERED.tsv".format(OUTPUT),
        feat_outfile = "{}.{{datatype}}+{{metric}}+{{region_type}}.FILTERED.featureLocations.tsv".format(OUTPUT)
    conda: CONDA
    params:
        extraction_script = os.path.join(workflow.basedir, "generatePerBaseFeatureScores.py"),
        pthresh = 0.05,
        fcthresh = 2,
        ocrmetrics_data_dir = OCRMETRICS_DATA_DIR,
        ref_fasta = REF_GENOME_REFNAME,
    shell:
        """python {params.extraction_script} \
        --target_genes {input.target_genes} \
        --region_file {input.regions} \
        --result_file {input.resultfile} \
        --pthresh {params.pthresh} --fcthresh {params.fcthresh} \
        --ocrmetrics_data_dir {params.ocrmetrics_data_dir} \
        --datatype {wildcards.datatype} \
        --metric {wildcards.metric} \
        --region_set {wildcards.region_type} \
        --refgenome_fasta {params.ref_fasta} \
        --expressed_features {input.expressed_tf} \
        --output {output.outfile} \
        --output_features {output.feat_outfile}"""