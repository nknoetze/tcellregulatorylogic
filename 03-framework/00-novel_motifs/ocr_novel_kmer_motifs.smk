## Novel kmer motifs - OCR

## July 10, 2023
## Scott Brown

VERSION = "1.0"
VERSION_DATE = "230610"


import os
import time
import glob
import logging
#import pybedtools

###############################
##  Define Global Variables  ##
###############################
#windows = [8, 12, 16]
#windows = [12]
#mismatches = [1]
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
#WORKING_DIR = "/projects/sbrown_prj/220318_DI/data/processed/novel_kmer/"
#RUN_ID = "240430"

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
        # regions = "/projects/nknoetze_prj/ocr_prj/results/ocr_metric_files/meuleman_tcell_2sample_1cd8_1cd4_ocrs_merged/Loop_2.5kb/merged_interactions/connectivity8/tss_within/pulled_ocrs_tss_dist.tsv",
        regions = REGIONS,
        # genes = "/projects/sbrown_prj/sbrown_scratch/220318_DI/data/raw/ocrs/4sigma_target_genes.txt",
        genes = TARGET_GENES,
        # ref_fasta = "/projects/sbrown_prj/220318_DI/data/raw/ref/GRCh37.p13.genome.fa"
        ref_fasta = REF_FASTA,
        script=os.path.join(workflow.basedir,"data/target_regions_to_fasta.py")
    conda: FRAMEWORK_ENV
    output:
        fasta = os.path.join(WORKING_DIR,RUN_ID,"00_ocr_seq","{region_type}+sequences.fasta")
    params:
        region_type = "all"
        # region_type = lambda wildcards: wildcards.region_type
    shell:
        "python {input.script} --regions {input.regions} --genes {input.genes} --ref_fasta {input.ref_fasta} --region_type {params.region_type} --output {output}"
    # run:
    #     REF_GENOME_FASTA = pybedtools.BedTool(input.ref_fasta)

    #     target_genes = []
    #     for line in open(input.genes, "r"):
    #         target_genes.append(line.rstrip())

    #     regions = {}

    #     header = None
    #     for line in open(input.regions, "r"):
    #         if not header:
    #             header = line.rstrip().split("\t")
    #         else:
    #             ## not header line
    #             # chrom   start   end     ocr_id  score   strand  ocr_type        query_gene      promoter_id     query_name      gene_group      data    dist_tss        promoter_name   promoter_type
    #             linedict = {x:y for x,y in zip(header, line.rstrip().split("\t"))}
    #             region_id = "{}\t{}".format(linedict["ocr_id"].replace(":","\t").replace("-","\t"), linedict["ocr_id"])
    #             ## region type - need to compare query_gene and promoter_id.
    #             promoter_id = linedict["promoter_id"]
    #             gene_id = linedict["query_gene"]
    #             region_type = linedict["ocr_type"]

    #             ## only use regions for genes of interest!
    #             if gene_id in target_genes:

    #                 ## if promoter but not of gene in question
    #                 if region_type == "promoter" and gene_id != promoter_id:
    #                     ## not promoter for this gene, so it is what we are calling a "linked promoter".
    #                     region_type = "linked_promoter"

    #                 if params.region_type == "all":
    #                     ## we will track all regions and their region_type for this set. Priority of region_type is promoter > linked_promoter > non-promoter
    #                     ## any region will either be non-promoter, or promoter/linked_promoter. 

    #                     ## NOTE April 29, 2024 - Updating to not include linked_promoters in the 'all' set. This should have been done originally.
    #                     if region_type != "linked_promoter":
                    
    #                         if region_id not in regions:
    #                             regions[region_id] = region_type
    #                         else:
    #                             if regions[region_id] != "promoter":
    #                                 ## update - could potentially be promoter now. Note this may result in no change.
    #                                 regions[region_id] = region_type
    #                 else:
    #                     ## only track regions of interest
    #                     if region_type == params.region_type:
    #                         regions[region_id] = region_type

    #     ## only take unqiue regions, because if we are looking at the same OCR, we don't want to double count all the kmers in that OCR.
    #     regions = [x for x in set(regions)]

    #     ## convert to pybedtools object
    #     bed_regions = pybedtools.BedTool("\n".join([r for r in regions]), from_string=True)
    #     ## ONLY WRITE ONE SEQUENCE FOR DEBUGGING
    #     ##bed_regions = pybedtools.BedTool("\n".join([list(regions.keys())[0]]), from_string=True)

    #     ## get sequences
    #     a = bed_regions.sequence(fi=REF_GENOME_FASTA)     ## pybedtools sequence (NOT temp file)
    #     region_seqs = open(a.seqfn).read()
    #     ## write file.
    #     out = open(output.fasta, "w")
    #     out.write(region_seqs)
    #     out.close()

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
        #script = "/projects/sbrown_prj/220318_DI/src/codebase/tcell_ocr_prj/kmer_cluster_finder.py",
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
        #script = "/projects/sbrown_prj/220318_DI/src/analysis/kmer_clusters_to_meme.py",
        script=os.path.join(workflow.basedir,"kmer_clusters_to_meme.py")
    shell:
        """python {params.script} --kmer_clusters {input.dot} --meme_output {output.meme}"""


##########################
##  03 Search JASPAR DB ##
##########################

rule search_jaspar:
    input:
        meme = rules.generate_meme.output.meme,
        #db = "/projects/nknoetze_prj/ocr_prj/data/raw/jaspar/20230605205709_JASPAR2022_combined_matrices_2090_human_meme.txt"
        db=JASPAR_DB
    conda: MOTIF_ENV
    output:
        tomtomdir = directory(os.path.join(WORKING_DIR,RUN_ID,"03_tomtom","tomtom_{region_type}+win{window}+thresh{mismatch}_tomtom"))
    params:
        #script = "/projects/sbrown_prj/220318_DI/src/analysis/kmer_clusters_to_meme.py",
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
        # script = "/projects/sbrown_prj/220318_DI/src/analysis/meme2jaspar.py"
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
        #script = "/projects/sbrown_prj/220318_DI/src/analysis/motif_clustering_Workflow_v2.1beta-human_230725.py"
        # script=os.path.join(workflow.basedir,"motif_clustering_Workflow_v2.1beta-human_230725.py")
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