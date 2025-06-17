TITLE = "Test Features"
DESC = """Given a list of target genes of length n and a list of background genes length m, this script will: 
            - Obtain Target regions for the Target Genes
            - Test all defined Features within these regions.
            - Randomly select n genes from background genes, test features on these n genes.
            - Repeat the random resampling k times
            - Use the distribution of randomly resampled values for each feature to generate a p-value for each feature in the target gene list.
"""

'''
Date: January 05, 2023
@authors: sbrown, nknoetze
'''

## Import Stock Libraries

import sys
import argparse
import os
import time
import random
import multiprocessing as mp
import queue
import traceback
import math

import pybedtools
## Cite: 
''' 
Dale RK, Pedersen BS, and Quinlan AR. 2011. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. Bioinformatics 27(24):3423-3424.
Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26(6):841–842.
'''

import ocrmetrics

## Classes and functions

class bcolors:
    CYAN = '\033[1;36;40m'
    BLUE = '\033[1;34;40m'
    GREEN = '\033[1;32;40m'
    YELLOW = '\033[1;33;40m'
    RED = '\033[1;31;40m'
    BOLDWHITE = '\033[1;37;40m'
    DARK = '\033[1;30;40m'
    PURPLE = '\033[1;35;40m'
    ENDC = '\033[0m'

def statprint(msg, msg_type = "STATUS"):
    typeColour = ""
    if msg_type == "ERROR":
        typeColour = bcolors.RED
    elif msg_type == "WARNING":
        typeColour = bcolors.YELLOW
    elif msg_type == "DEBUG":
        typeColour = bcolors.GREEN
    elif msg_type == "SUBPROCESS":
        typeColour = bcolors.GREEN
        msg_type = "     " + msg_type
    else:
        typeColour = bcolors.BOLDWHITE

    print("{message_color}{message_type}{end_color} {time_color}[{datetime}]{end_color}: {message}".format(message_color = typeColour, 
             message_type = msg_type,
             end_color = bcolors.ENDC, 
             time_color = bcolors.BLUE, 
             datetime = time.strftime("%Y/%m/%d %T"), 
             message = msg), flush=True)


def processJob(in_q, out_q, i):
    '''Given an input queue, output queue, and index, process one unit'''
    # holder for results
    # [jobname, value]
    resHolder = []
    numInHolder = 0

    while True:
        dat = in_q.get()
        if dat is SENTINEL:
            break

        these_results = None

        jobname, func, feat, region_coords, region_sequences, region_gene_dists, set_genes, regions_genes, num_motifs, ref_fasta, chrom_sizes, ref_bed, excl_regions_bed, scramble_method, region_set  = dat

        these_results = FUNCTIONS_TO_RUN[func](region_coords = region_coords,
                                            region_sequences = region_sequences,
                                            feature = feat,
                                            region_gene_dists = region_gene_dists,
                                            set_genes = set_genes,
                                            regions_genes = regions_genes,
                                            jobname = jobname,
                                            num_motifs = num_motifs,
                                            ref_fasta = ref_fasta,
                                            chrom_sizes = chrom_sizes,
                                            ref_bed = ref_bed,
                                            excl_regions_bed = excl_regions_bed,
                                            scramble_method = scramble_method, 
                                            region_set = region_set)
        
        resHolder.append([jobname, these_results])

        numInHolder += 1

        if numInHolder == MAXHOLDERCOUNT:
            out_q.put(resHolder)
            resHolder = []
            numInHolder = 0
    # SENTINEL reached, write out remaining data
    out_q.put(resHolder)
    resHolder = []
    numInHolder = 0

    statprint("Processing in subprocess {} complete...".format(i), "SUBPROCESS")

def deleteTmpFile(in_q, i):
    '''given a pybedtools object, delete the associated tmp file'''
    while True:
        obj = in_q.get()
        if obj is SENTINEL:
            break
        os.remove(obj.fn)




## Declare global variables

DEBUG = False
VERB = False
SENTINEL = None
#MAXHOLDERCOUNT = 1000
MAXHOLDERCOUNT = 1



# This stores the functions in a dictionary so that then functions can be called in a loop. This also makes the code more amenable to multithreading if we need.
# Note that this must be in the form of a dictionary with key as the string and value as the ocrmetrics.function of the same name.
FUNCTIONS_TO_RUN = {"fimo_geneCount": ocrmetrics.fimo_geneCount,
                    "fimo_geneCount_shuffling": ocrmetrics.fimo_geneCount_shuffling,
                    "fimo_geneCount_scrambling": ocrmetrics.fimo_geneCount_scrambling,
                    "fimo_comotif_geneCount": ocrmetrics.fimo_comotif_geneCount,
                    "fimo_comotif_geneCount_shuffling": ocrmetrics.fimo_comotif_geneCount_shuffling,
                    "fimo_comotif_geneCount_scrambling": ocrmetrics.fimo_comotif_geneCount_scrambling
                    }

FEATURE_METRIC_FLIPPED_FUNCTIONS = ["fimo_geneCount","fimo_geneCount_shuffling","fimo_geneCount_scrambling","fimo_comotif_geneCount","fimo_comotif_geneCount_shuffling","fimo_comotif_geneCount_scrambling"]
SHUFFLING_FUNCTIONS = ["fimo_geneCount_shuffling","fimo_geneCount_scrambling","fimo_comotif_geneCount_shuffling","fimo_comotif_geneCount_scrambling"]
NON_SHUFFLING_FUNCTIONS = ["fimo_geneCount","fimo_comotif_geneCount"]

## Main

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = TITLE)
    parser.add_argument("--target_genes", dest = "TARGET_GENES", help = "Text file of target genes, ENSG, one per line.", type = str)
    parser.add_argument("--background_genes", dest = "BACKGROUND_GENES", help = "Text file of background genes, ENSG, one per line.", type = str)
    parser.add_argument("--all_genes", dest = "ALL_GENES", help = "Ranked gene list to map ENSG to gene name.", type = str)
    parser.add_argument("--refgenome_fasta", dest = "REF_FASTA", help = "Genome reference fasta file.", type = str)
    parser.add_argument("--chrom_sizes_file", dest = "CHROM_SIZES_FILE", help = "file of chromosome sizes for shuffling", type = str)
    parser.add_argument("--refgenome_bed", dest = "REF_GENOME_BED", help = "Bed file of chromosome sizes used for shuffling.", type = str)
    parser.add_argument("--excl_regions_bed", dest = "EXCL_REGIONS_BED", help = "Bed file of centromere and telomere regions.", type = str)
    parser.add_argument("--region_file", dest = "REGION_FILE", help = "File of regions linked to genes with additional info columns.", type = str)
    parser.add_argument("--ocrmetrics_data_dir", dest = "OCRMET_DAT_DIR", help = "Directory where ocrmetrics data files can be found.", type = str)
    parser.add_argument("--num_genes", dest = "NUM_GENES", help = "Number of genes to randomly select for each iteration. By default, this will be the size of 'target_genes'.", type = int)
    parser.add_argument("--num_iterations", dest = "NUM_ITER", help = "Number of iterations to repeatedly subsample background genes for background distribution.", type = int)
    parser.add_argument("--functions", dest = "SELECTED_FUNCTIONS", help = "List of function names to run, if only a subset is requested", type = str, nargs = "*", choices = FUNCTIONS_TO_RUN.keys())
    parser.add_argument("--num_motifs", dest = "NUM_MOTIFS", help = "List of the number of Motifs to explore for co-occuring motifs function", type = int, nargs = "*")
    parser.add_argument("--scramble_method", dest = "SCRAMBLE_METHOD", help = "dinucleotide, or mononucleotide scrambling", type = str)
    parser.add_argument("--distances", dest="DISTANCES", action = "store_true", help = "Flag to test gene-region distance.")
    parser.add_argument("--seed", dest = "SEED", help = "Random seed.", type = int, default = 56866345)
    parser.add_argument("--cores", dest = "maxNumberProcesses", help = "Number of cores to use.", type = int, default = 8)
    parser.add_argument("--max_bg_written", dest = "MAX_BACKGROUND_VALS_WRITTEN", help = "Max number of background values to write.", type = int, default = 1000)
    parser.add_argument("--sigonly", dest = "SIGONLY", action = "store_true", help = "Flag to only write significant (p<0.05) results", default = False)
    parser.add_argument("--output", dest = "OUTPUT", help = "File with results to write, .tsv format.", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    random.seed(args.SEED)

    ## set pybedtools tmp dir:
    pybedtools.set_tempdir('/projects/holtlab_prj/scratch/tmp')


    ## Welcome Message
    print("="*80)
    print("Python version: {}".format(sys.version))
    print("Python environment: {}".format(sys.prefix))
    print("Server: {}".format(os.uname()[1]))
    print("Current directory: {}".format(os.getcwd()))
    print("Command: {}".format(" ".join(sys.argv)))
    print("Time: {}".format(time.strftime("%Y/%m/%d %T")))
    print("{}\n".format("="*80))


    REF_GENOME_FASTA = pybedtools.BedTool(args.REF_FASTA)


    if args.SELECTED_FUNCTIONS and args.SELECTED_FUNCTIONS != "":
        ## subset FUNCTIONS_TO_RUN
        FUNCTIONS_TO_RUN = {x:FUNCTIONS_TO_RUN[x] for x in args.SELECTED_FUNCTIONS}
    

    ## initialize ocrmetrics
    statprint("Initializing ocrmetrics...")
    ocrmetrics_features = ocrmetrics.init(args.OCRMET_DAT_DIR, FUNCTIONS_TO_RUN.keys(),args.NUM_MOTIFS)

    statprint("The following feature types will be calculated:")
    [statprint("\t{}: {} feature(s)".format(func, len(ocrmetrics_features[func]))) for func in FUNCTIONS_TO_RUN]

    ## Reading in Gene Names for error reporting
    statprint("Reading in gene name annotations...")
    genes = {}
    header = None
    for line in open(args.ALL_GENES, "r"):
        if not header:
            header = line.rstrip().split("\t")
        else:
            linedict = {x:y for x,y in zip(header,line.rstrip().split("\t"))}
            genes[linedict["gene_id"]] = linedict["gene_name"]

    ## Read in Target Gene List
    statprint("Reading in target genes...")

    target_genes = []
    for line in open(args.TARGET_GENES, "r"):
        target_genes.append(line.rstrip())

    ## set NUM_GENES
    if not args.NUM_GENES:
        num_genes = len(target_genes)
        statprint("Setting num_genes to {}".format(num_genes))
    else:
        num_genes = args.NUM_GENES
        statprint("num_genes provided as {}".format(num_genes))


    ## Load in regions
    statprint("Loading all regions into memory...")
    ## initialize holder for region categories
    region_categories = ["all"]
    
    ## categorize regions by gene id
    regions = {}
    regions2genes = {}
    region_gene_distances = {}
    header = None
    for line in open(args.REGION_FILE, "r"):
        if not header:
            header = line.rstrip().split("\t")
        #if not line.startswith("chrom"):
        else:
            ## not header line
            # chrom   start   end     ocr_id  score   strand  ocr_type        query_gene      promoter_id     query_name      gene_group      data    dist_tss        promoter_name   promoter_type
            linedict = {x:y for x,y in zip(header, line.rstrip().split("\t"))}
            region_start = int(linedict["start"])
            region_end = int(linedict["end"])
            region_id = linedict["ocr_id"].replace(":","\t").replace("-","\t")
            ## region type - need to compare query_gene and promoter_id.
            promoter_id = linedict["promoter_id"]
            gene_id = linedict["query_gene"]
            region_type = linedict["ocr_type"]
            
            ## if promoter but not of gene in question
            if region_type == "promoter" and gene_id != promoter_id:
                ## not promoter for this gene, so it is what we are calling a "linked promoter".
                region_type = "linked_promoter"
                ## as of June 1, 2023, we no longer want to consider linked promoters, so dropping these.
            else:
                if region_type not in region_categories: region_categories.append(region_type)
                
                if gene_id not in regions:
                    regions[gene_id] = []
                ## store the tuple of region (broken into tab-separated chr start end) and region_type
                regions[gene_id].append((region_id, region_type))

                if region_id not in regions2genes:
                    regions2genes[region_id] = set()
                regions2genes[region_id].add(gene_id)

                dist = int(linedict["dist_tss"])
                region_gene_distances[(region_id, gene_id)] = dist

    try:
        ## set up multiprocessing
        statprint("Creating worker subprocesses...")
        metric_procs = []
        metric_iq = mp.Queue()
        metric_oq = mp.Queue()

        metric_num_submitted = 0
        metric_num_finished = 0

        for i in range(args.maxNumberProcesses):
            if VERB: statprint("\tProcess {} starting...".format(i+1))
            p = mp.Process(target=processJob, args=(metric_iq, metric_oq, i+1))
            metric_procs.append(p)
            p.start()
        ## now metric_iq is ready to receive jobs, like:
        ## [jobname, func, feat, region_coords, region_sequences, region_gene_dists]


        ## set up queue/job for deleting tmp files.
        statprint("Spinning up temp file deletion queue...")
        deletion_iq = mp.Queue()
        deletion_process = mp.Process(target=deleteTmpFile, args=(deletion_iq, 1))
        deletion_process.start()

        ## Fetch regions for Target Gene List (All regions, Distal regions, Promoter regions)
        statprint("Fetching Distal and Promoter regions for Target Genes...")
        
        target_regions = {region_cat: [] for region_cat in region_categories}
        ## each element is like "chr1\t1\t500\n..."
        target_region_gene_distances = {region_cat: [] for region_cat in region_categories}
        target_region_genes = {region_cat: [] for region_cat in region_categories}
        for tgene in target_genes:
            #if DEBUG: statprint("Target gene {} being processed...".format(genes[tgene]), "DEBUG")
            #this_gene_regions = {region_cat: [] for region_cat in region_categories}
            if tgene in regions:
                for tregion in regions[tgene]:
                    #if DEBUG: statprint("Processing {}--{}...".format(tregion[0], tregion[1]), "DEBUG")
                    target_regions["all"].append(tregion[0]+"\t"+tgene)
                    target_region_gene_distances["all"].append(region_gene_distances[(tregion[0], tgene)])
                    ## region-type specific lists:
                    target_regions[tregion[1]].append(tregion[0]+"\t"+tgene)
                    target_region_gene_distances[tregion[1]].append(region_gene_distances[(tregion[0], tgene)])
                    ## get all genes linked to this region
                    target_region_genes["all"].append(regions2genes[tregion[0]]) ## append the set of genes linked to this region.
                    target_region_genes[tregion[1]].append(regions2genes[tregion[0]])

            else:
                statprint("{} - {} did not have any regions associated with it".format(tgene, genes[tgene]), "WARNING")
        
        ## convert list of regions to pybedtools.BedTool object
        statprint("Convert region strings to one pybedtools.BedTool object")
        for region_type in region_categories:
            target_regions[region_type] = pybedtools.BedTool("\n".join(target_regions[region_type]), from_string=True)


        ## Get sequences of regions
        statprint("Fetching sequences of regions...")
        target_region_seqs = {region_cat: [] for region_cat in region_categories}

        for region_type in region_categories:
            a = target_regions[region_type].sequence(fi=REF_GENOME_FASTA, name=True)     ## pybedtools sequence (NOT temp file)
            target_region_seqs[region_type] = open(a.seqfn).read()            ## fasta entries

        ## Test Features
        statprint("Measuring features in target regions...")
        ## {"all": {feature1: value, feature2: value, ...},
        #   "promoter": {feature1: value, feature2: value, ...},
        #   "distal": {feature1: value, feature2: value, ...}}
        target_feature_results = {region_cat: {func: {feat: {} for feat in ocrmetrics_features[func]} for func in FUNCTIONS_TO_RUN} for region_cat in region_categories}

        # loop through region types
        for region_type in region_categories:
            # loop through functions/metrics
            for func in FUNCTIONS_TO_RUN:
                for feat in ocrmetrics_features[func]:
                    if VERB: statprint("\t...{}:{}".format(func, feat))
                    metric_iq.put(["target\t0\t{}\t{}\t{}".format(region_type, func, feat), func, feat, target_regions[region_type], target_region_seqs[region_type], target_region_gene_distances[region_type], target_genes, target_region_genes[region_type],args.NUM_MOTIFS,args.REF_FASTA,args.CHROM_SIZES_FILE,args.REF_GENOME_BED,args.EXCL_REGIONS_BED,args.SCRAMBLE_METHOD,"target"])
                    metric_num_submitted += 1

        ## wait for target results to finish so that we can calculate p-val as we go
        statprint("All 'target' jobs submitted to queue...")


        statprint("Monitoring queue to parse completed 'target' jobs...")
        while metric_num_finished < metric_num_submitted:
            ## slurp results
            res = metric_oq.get()
            if res is not SENTINEL:
                for item in res:    ## because I write to metric_oq in chunks to save IO.
                    ## completed job
                    metric_num_finished += 1
                    if metric_num_finished % (int(0.1*metric_num_submitted) + 1) == 0: statprint("{} jobs finished processing...".format(metric_num_finished))
                    jobname, valuedict = item
                    genegroup, repnum, region_type, func, feat = jobname.split("\t")

                    for metric in valuedict:
                        if genegroup == "target":
                            target_feature_results[region_type][func][feat][metric] = valuedict[metric]

        ###########
        ## dictionary to hold feature results for background sets
        background_feature_results = {region_cat: {func: {feat: {} for feat in ocrmetrics_features[func]} for func in FUNCTIONS_TO_RUN} for region_cat in region_categories}

        non_shuffling_flag=False
        for func in FUNCTIONS_TO_RUN:
            if func in SHUFFLING_FUNCTIONS:
                statprint("Measuring features in shuffled background regions...")
                for i in range(args.NUM_ITER):
                    ## write to console every 100 iterations
                    if i%100==0: statprint("\t{} iterations complete.".format(i))
                    ## Test Features
                    ## {"all": {feature1: value, feature2: value, ...},
                    #   "promoter": {feature1: value, feature2: value, ...},
                    #   "distal": {feature1: value, feature2: value, ...}}

                    # loop through region types
                    for region_type in region_categories:
                        # loop through functions/metrics
                        #no longer need to loop through functions to run.
                        #for func in FUNCTIONS_TO_RUN:
                        for feat in ocrmetrics_features[func]:
                            if VERB: statprint("\t...{}:{}".format(func, feat))
                            # pass all data to keep things standard. Function will receive all and only use what it needs.
                            #since we are shuffling the target gene regions, feed in target regions again!
                            metric_iq.put(["background\t0\t{}\t{}\t{}".format(region_type, func, feat), func, feat, target_regions[region_type], target_region_seqs[region_type], target_region_gene_distances[region_type], target_genes, target_region_genes[region_type],args.NUM_MOTIFS,args.REF_FASTA,args.CHROM_SIZES_FILE,args.REF_GENOME_BED,args.EXCL_REGIONS_BED,args.SCRAMBLE_METHOD,"background"])
                            metric_num_submitted += 1
            else:
                non_shuffling_flag=True
                
        #ONLY RUN THIS IF THERE ARE NON-SHUFFLING FUNCTIONS!
        #check if the flag was triggered. if atleast one of the requestion functions to run was a non-shuffling function. this will be True
        if non_shuffling_flag==True:
            ## Read in Background Gene List
            statprint("Reading in background gene list...")
            background_genes = []
            for line in open(args.BACKGROUND_GENES, "r"):
                background_genes.append(line.rstrip())

            statprint("Randomly subsampling background genes, calculating features, repeating {} times...".format(args.NUM_ITER))
            ## for i in 1:k:
            for i in range(args.NUM_ITER):
                ## write to console every 100 iterations
                if i%100==0: statprint("\t{} iterations complete.".format(i))

                ## Randomly select n genes from Background Gene List (without replacement)
                sampled_genes = random.sample(background_genes, num_genes)

                sampled_regions = {region_cat: [] for region_cat in region_categories}
                ## each element is like "chr1\t1\t500\n..."

                sampled_region_gene_distances = {region_cat: [] for region_cat in region_categories}

                sampled_region_genes = {region_cat: [] for region_cat in region_categories}

                ## make sure all sampled genes have regions.
                for sgene in sampled_genes.copy():
                    # check that this gene has regions. If not, repick (and make sure not double picking)
                    while sgene not in regions:
                        sampled_genes.remove(sgene)
                        sgene = random.sample(background_genes, 1)[0]
                        ## if gene already in subsampled set, pick again
                        while sgene in sampled_genes:
                            sgene = random.sample(background_genes, 1)[0]
                        sampled_genes.append(sgene)
                        ## now we have a new sgene that was not in the original set, and has regions.
                
                ## Fetch regions for Random Gene list
                for sgene in sampled_genes:
                    if sgene in regions:
                        for sregion in regions[sgene]:
                            sampled_regions["all"].append(sregion[0]+"\t"+sgene)
                            sampled_region_gene_distances["all"].append(region_gene_distances[(sregion[0], sgene)])
                            ## region-type specific lists:
                            sampled_regions[sregion[1]].append(sregion[0]+"\t"+sgene)
                            sampled_region_gene_distances[sregion[1]].append(region_gene_distances[(sregion[0], sgene)])
                            ## get all genes linked to this region
                            sampled_region_genes["all"].append(regions2genes[sregion[0]]) ## append the set of genes linked to this region.
                            sampled_region_genes[sregion[1]].append(regions2genes[sregion[0]])

                    else:
                        ## should never get here, as we ensure all genes have regions above.
                        statprint("{} - {} did not have any regions associated with it".format(sgene, genes[sgene]), "WARNING")

                ## convert list of regions to pybedtools.BedTool object
                #statprint("Convert region strings to one pybedtools.BedTool object")
                for region_type in region_categories:
                    sampled_regions[region_type] = pybedtools.BedTool("\n".join(sampled_regions[region_type]), from_string=True)


                ## Fetch region sequences
                sampled_region_seqs = {region_cat: [] for region_cat in region_categories}

                for region_type in region_categories:
                    a = sampled_regions[region_type].sequence(fi=REF_GENOME_FASTA, name=True)    ## pybedtools sequence (NOT temp file)
                    sampled_region_seqs[region_type] = open(a.seqfn).read()            ## fasta entries

                ## Test Features
                for region_type in region_categories:
                    # loop through functions/metrics
                    for func in FUNCTIONS_TO_RUN:
                        if func in NON_SHUFFLING_FUNCTIONS:
                            for feat in ocrmetrics_features[func]:

                                metric_iq.put(["background\t{}\t{}\t{}\t{}".format(i, region_type, func, feat), func, feat, sampled_regions[region_type], sampled_region_seqs[region_type], sampled_region_gene_distances[region_type], sampled_genes, sampled_region_genes[region_type],args.NUM_MOTIFS,args.REF_FASTA,args.CHROM_SIZES_FILE,args.REF_GENOME_BED,args.EXCL_REGIONS_BED,args.SCRAMBLE_METHOD,"background"])
                                metric_num_submitted += 1

                            if metric_num_submitted % (int(0.1*len(FUNCTIONS_TO_RUN)*args.NUM_ITER) + 1) == 0: statprint("...{} jobs submitted...".format(metric_num_submitted))
        
        
        ## Done submitting jobs to the multiprocessing
        statprint("All jobs submitted to queue...")
        ## Send END signal to queues (one for each subprocess)
        for i in range(args.maxNumberProcesses):
            metric_iq.put(SENTINEL)
        deletion_iq.put(SENTINEL)

        statprint("Monitoring queue to parse completed jobs...")
        while metric_num_finished < metric_num_submitted:
            ## slurp results
            res = metric_oq.get()
            if res is not SENTINEL:
                for item in res:    ## because I write to metric_oq in chunks to save IO.
                    ## completed job
                    metric_num_finished += 1
                    if metric_num_finished % (int(0.1*metric_num_submitted) + 1) == 0: statprint("{} jobs finished processing...".format(metric_num_finished))
                    jobname, valuedict = item
                    genegroup, repnum, region_type, func, feat = jobname.split("\t")

                    for metric in valuedict:
                        if genegroup == "target":
                            target_feature_results[region_type][func][feat][metric] = valuedict[metric]
                        elif genegroup == "background":
                            if metric not in background_feature_results[region_type][func][feat]:
                                background_feature_results[region_type][func][feat][metric] = [0,0] ## sum, p-value count
                            ## sum of values
                            background_feature_results[region_type][func][feat][metric][0] += valuedict[metric]
                            ## number values >= target value
                            if metric in target_feature_results[region_type][func][feat] and valuedict[metric] >= target_feature_results[region_type][func][feat][metric]:
                                background_feature_results[region_type][func][feat][metric][1] += 1

        


        ## make sure all metric_procs have finished.
        statprint("Shutting down all subprocesses...")
        for p in metric_procs:
            p.join()
        deletion_process.join()
    
    except:
        statprint("{}".format(sys.exc_info()[0]), "ERROR")
        statprint(traceback.format_exc(), "ERROR")
        statprint("Please wait while processes close...", "WARNING")
        while not metric_iq.empty():
            try:
                foo = metric_iq.get(block=False)
            except queue.Empty:
                pass
        for i in range(args.maxNumberProcesses):
            metric_iq.put(SENTINEL)
        for p in metric_procs:
            p.join()

        while not deletion_iq.empty():
            try:
                foo = deletion_iq.get(block=False)
            except queue.Empty:
                pass
        deletion_iq.put(SENTINEL)
        deletion_process.join()

        sys.exit("Exiting.")



    ## delete all pybedtools tmp files
    statprint("Deleting remaining pybedtools tmp files...")
    pybedtools.cleanup()
    
    ## dictionary to hold output
    ## {"all": {feature1: {"target_value": -1, "background_mean": -1, "p-value": -1, "background_values": []},
    #           feature2: {"target_value": -1, "background_mean": -1, "p-value": -1, "background_values": []},
    #           ...},
    #   "promoter": {feature1: {"target_value": -1, "background_mean": -1, "p-value": -1, "background_values": []},
    #           feature2: {"target_value": -1, "background_mean": -1, "p-value": -1, "background_values": []},
    #           ...},
    #   "distal": {feature1: {"target_value": -1, "background_mean": -1, "p-value": -1, "background_values": []},
    #           feature2: {"target_value": -1, "background_mean": -1, "p-value": -1, "background_values": []},
    #           ...}}
    statprint("Populating Metrics...")
    ## populate metrics based on set of feature keys in targets and background (as with kmers, target and/or background may not contain every single metric)
    out_data = {region_cat: {func: {feat: {metric: {"target_value": -1, "background_mean": -1, "p-value": -1} for metric in set(list(target_feature_results[region_cat][func][feat].keys()) + list(background_feature_results[region_cat][func][feat].keys()))} for feat in target_feature_results[region_cat][func]} for func in target_feature_results[region_cat]} for region_cat in region_categories}
    
    ## For each Feature and Target/Background region comparison
    statprint("Calculating stats...")
    for region_cat in region_categories:
        for func in out_data[region_cat]:
            for feat in out_data[region_cat][func]:
                for metric in out_data[region_cat][func][feat]:
                    ## get target value, 0 if not exists
                    ## NOTE that we are only populating missing metrics with 0. To date, only kmers and fimo comotifs may not report values for all metrics.
                    target_val = target_feature_results[region_cat][func][feat].get(metric,0)

                    
                    ## because we are not tracking an array and instead just a value, we need to check that the metric exists.
                    if metric in background_feature_results[region_cat][func][feat]:
                        background_average = background_feature_results[region_cat][func][feat][metric][0]/args.NUM_ITER
                        p_val = background_feature_results[region_cat][func][feat][metric][1] / args.NUM_ITER
                    else:
                        background_average = 0
                        p_val = 0

                    
                    #only keep significant results
                    if (args.SIGONLY and p_val < 0.05) or not args.SIGONLY:
                        out_data[region_cat][func][feat][metric]["target_value"] = target_val
                        out_data[region_cat][func][feat][metric]["background_mean"] = background_average
                        out_data[region_cat][func][feat][metric]["p-value"] = p_val

    output_data = {region_cat: {func: {feat: {metric: data for metric, data in out_data[region_cat][func][feat].items() if data["p-value"] != -1} for feat in out_data[region_cat][func]} for func in out_data[region_cat]} for region_cat in out_data}

    ## Output summary of Features, region comparison, Feature Score (Target, mean of Background?), p-values.
    statprint("Writing output file...")
    out = open(args.OUTPUT, "w")
    
    output_lines=[]
    output_lines.append("datatype\tfeature\tmetric\tregion_set\ttarget_value\tbackground_mean\tp_val\n")

    for region_cat in output_data:
        for func in output_data[region_cat]:
            for feat in output_data[region_cat][func]:
                for metric in output_data[region_cat][func][feat]:

                    ## to save output size, only write if target value > 0
                    if output_data[region_cat][func][feat][metric]["target_value"] > 0:

                        ## for some functions, "Feature" and "Metric" were flipped to allow a single analysis run outputting all features. Typically, a function returns a set of different metrics for one feature.
                        ## To "patch" this, we'll flip it if it is one of those functions.
                        if func in FEATURE_METRIC_FLIPPED_FUNCTIONS:
                            if "kmer" in func:
                                ## write metric as kmer length
                                lines="{datatype}\t{feature}\t{metric}\t{region_set}\t{target_val}\t{background_val}\t{p_val}\n".format(
                                    datatype = func,
                                    feature = metric,
                                    metric = "{}mer".format(len(metric)),
                                    region_set = region_cat,
                                    target_val = output_data[region_cat][func][feat][metric]["target_value"],
                                    background_val = output_data[region_cat][func][feat][metric]["background_mean"],
                                    p_val = output_data[region_cat][func][feat][metric]["p-value"],
                                )

                            else:
                                lines="{datatype}\t{feature}\t{metric}\t{region_set}\t{target_val}\t{background_val}\t{p_val}\n".format(
                                    datatype = func,
                                    feature = metric,
                                    metric = feat,
                                    region_set = region_cat,
                                    target_val = output_data[region_cat][func][feat][metric]["target_value"],
                                    background_val = output_data[region_cat][func][feat][metric]["background_mean"],
                                    p_val = output_data[region_cat][func][feat][metric]["p-value"],
                                )

                        else:
                            lines="{datatype}\t{feature}\t{metric}\t{region_set}\t{target_val}\t{background_val}\t{p_val}\n".format(
                                datatype = func,
                                feature = feat,
                                metric = metric,
                                region_set = region_cat,
                                target_val = output_data[region_cat][func][feat][metric]["target_value"],
                                background_val = output_data[region_cat][func][feat][metric]["background_mean"],
                                p_val = output_data[region_cat][func][feat][metric]["p-value"],
                            )    
                        output_lines.append(lines)
    out.write("".join(output_lines))
    out.close()
        
    

    statprint("Done.")
