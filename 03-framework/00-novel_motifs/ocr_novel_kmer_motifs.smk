## Novel kmer motifs - OCR

## July 10, 2023
## Scott Brown

VERSION = "1.0"
VERSION_DATE = "230610"


import os
import time
import glob
import logging

###############################
##  Define Global Variables  ##
###############################
region_types = ["all"]

## Load Config
windows = config["WINDOWS"]
mismatches = config["MISMATCHES"]
RUN_ID = config["RUN_ID"]
WORKING_DIR=config["WORKING_DIR"]
REGIONS = config["REGIONS"]
TARGET_GENES = config["TARGET_GENES"]
REF_FASTA = config["REF_FASTA"]
JASPAR_DB=config["JASPAR_DB"]
FRAMEWORK_ENV=config["FRAMEWORK_ENV"]
MOTIF_ENV=config["MOTIF_ENV"]
ARCHETYPE_ENV=config["ARCHETYPE_ENV"]

################
##  Rule All  ##
################

rule all:
    input:
        expand(os.path.join(WORKING_DIR,RUN_ID, "05_clustered_motifs", "{region_type}+win{window}+thresh{mismatch}_archetypes.meme"), region_type=region_types, window=windows, mismatch=mismatches),
        expand("~/.local/src/genome-tools/setup_complete")



#############################
##  00 Generate OCR Fastas ##
#############################

rule make_fasta:
    input:
        regions = REGIONS,
        genes = TARGET_GENES,
        ref_fasta = REF_FASTA,
        script=os.path.join(workflow.basedir,"data/target_regions_to_fasta.py")
    conda: FRAMEWORK_ENV
    output:
        fasta = os.path.join(WORKING_DIR,RUN_ID,"00_ocr_seq","{region_type}+sequences.fasta")
    params:
        region_type = "all"
    shell:
        "python {input.script} --regions {input.regions} --genes {input.genes} --ref_fasta {input.ref_fasta} --region_type {params.region_type} --output {output}"
 
#############################
##  01 Run kmer clustering ##
#############################

rule run_kmercluster:
    input:
        fasta = rules.make_fasta.output.fasta
    conda: FRAMEWORK_ENV
    output:
        dot = os.path.join(WORKING_DIR,RUN_ID,"01_kmer_clusters","{region_type}+win{window}+thresh{mismatch}.tsv")
    params:
        script=os.path.join(workflow.basedir,"kmer_cluster_finder.py")
    shell:
        """python {params.script} --sequence {input.fasta} --window {wildcards.window} --mismatch {wildcards.mismatch} --output {output.dot} -v"""


############################
##  02 Generate MEME file ##
############################

rule generate_meme:
    input:
        dot = rules.run_kmercluster.output.dot
    conda: FRAMEWORK_ENV
    output:
        meme = os.path.join(WORKING_DIR,RUN_ID,"02_meme","{region_type}+win{window}+thresh{mismatch}.meme.txt")
    params:
        script=os.path.join(workflow.basedir,"kmer_clusters_to_meme.py")
    shell:
        """python {params.script} --kmer_clusters {input.dot} --meme_output {output.meme}"""


##########################
##  03 Search JASPAR DB ##
##########################

rule search_jaspar:
    input:
        meme = rules.generate_meme.output.meme,
        db=JASPAR_DB
    conda: MOTIF_ENV
    output:
        tomtomdir = directory(os.path.join(WORKING_DIR,RUN_ID,"03_tomtom","tomtom_{region_type}+win{window}+thresh{mismatch}_tomtom"))
    params:
        script=os.path.join(workflow.basedir,"kmer_clusters_to_meme.py")
    shell:
        """tomtom -o {output.tomtomdir} {input.meme} {input.db}"""


#########################
##  04 Filter clusters ##
#########################

rule filter_clusters:
    '''Write out clusters that are not found in the TOMTOM output file'''
    input:
        tomtom = rules.search_jaspar.output.tomtomdir,
        meme = rules.generate_meme.output.meme,
    #conda: FRAMEWORK_ENV
    output:
        novel_meme = os.path.join(WORKING_DIR,RUN_ID,"04_filtered_clusters","{region_type}+win{window}+thresh{mismatch}.novelmeme.txt")
    params:
        tomtom_file = os.path.join(rules.search_jaspar.output.tomtomdir, "tomtom.tsv")
    run:
        not_novel = set()
        for line in open(params.tomtom_file, "r"):
            clust = line.split("\t")[0]
            not_novel.add(clust)
        
        out = open(output.novel_meme, "w")

        WRITING = True
        for line in open(input.meme, "r"):
            if line.startswith("MOTIF"):
                if line.split(" ")[1] in not_novel:
                    WRITING = False
                else:
                    WRITING = True

            if WRITING:
                out.write(line)
        
        out.close()


rule get_novel_cluster_sizes:
    '''Write file with cluster sizes'''
    input:
        novel_meme = rules.filter_clusters.output.novel_meme
    output:
        cluster_counts = os.path.join(WORKING_DIR,RUN_ID,"04_filtered_clusters","{region_type}+win{window}+thresh{mismatch}.novelmeme.counts.txt")
    shell:
        "grep MOTIF {input.novel_meme} | awk -F '.' '{{print $2}}' > {output.cluster_counts}"


###################################
##  05 Cluster Filtered Clusters ##
###################################

rule create_single_motif_pfms:
    '''From meme file, create single PFM files for each motif'''
    input:
        novel_meme = rules.filter_clusters.output.novel_meme
    conda: FRAMEWORK_ENV
    output:
        pfm_dir = directory(os.path.join(WORKING_DIR,RUN_ID, "04_filtered_clusters", "{region_type}+win{window}+thresh{mismatch}_pfms"))
    params:
        script=os.path.join(workflow.basedir,"meme2jaspar.py")
    shell:
        "mkdir {output.pfm_dir} && python {params.script} {input.novel_meme} {output.pfm_dir}"

rule cluster_novel_motifs:
    '''Cluster the novel motifs'''
    input:
        novel_meme = rules.filter_clusters.output.novel_meme
    conda: MOTIF_ENV
    output:
        motif_clustering = os.path.join(WORKING_DIR,RUN_ID, "05_clustered_motifs", "{region_type}+win{window}+thresh{mismatch}_motifclust.txt")
    shell:
        "tomtom -dist kullback -motif-pseudo 0.1 -text -min-overlap 1 {input.novel_meme} {input.novel_meme} > {output.motif_clustering}"

rule install_genome_tools:
    output:
        touch("~/.local/src/genome-tools/setup_complete")
    conda: ARCHETYPE_ENV
    shell:
        """
        module load git
        mkdir -p ~/.local/src && cd ~/.local/src
        git clone https://github.com/jvierstra/genome-tools.git
        cd genome-tools
        python setup.py install --user
        touch ~/.local/src/genome-tools/setup_complete
        """

rule generate_archetypes:
    '''From clusters of novel motifs, make archetypic motif'''
    input:
        setup="~/.local/src/genome-tools/setup_complete",
        motif_clustering = rules.cluster_novel_motifs.output.motif_clustering,
        pfm_dir = rules.create_single_motif_pfms.output.pfm_dir
    conda: ARCHETYPE_ENV
    output:
        plot_dir = directory(os.path.join(WORKING_DIR,RUN_ID, "05_clustered_motifs", "{region_type}+win{window}+thresh{mismatch}_archetype-plots")),
        archetypes = os.path.join(WORKING_DIR,RUN_ID, "05_clustered_motifs", "{region_type}+win{window}+thresh{mismatch}_archetypes.uniprobe")
    params:
        script=os.path.join(workflow.basedir,"motif_clustering_Workflow_v2.1beta-human_240503.py")
    shell:
        "mkdir {output.plot_dir} && python {params.script} --clustered_motifs {input.motif_clustering} --pfm_dir {input.pfm_dir} --plot_dir {output.plot_dir} --output {output.archetypes}"

rule uniprobe2meme:
    '''Convert Uniprobe to meme format'''
    input:
        archetypes = rules.generate_archetypes.output.archetypes
    conda: MOTIF_ENV
    output:
        archetypes = os.path.join(WORKING_DIR,RUN_ID, "05_clustered_motifs", "{region_type}+win{window}+thresh{mismatch}_archetypes.meme")
    shell:
        "uniprobe2meme {input.archetypes} > {output.archetypes}"