TITLE = "Generate Per Base Feature Scores"
DESC = "Using results from feature enrichment test framework, generate per-base feature scores across OCRs"
'''
Date: April 19,2023
@author: sbrown, nknoetze
'''

## Import Libraries

import sys
import argparse
import os
import time
import pybedtools
import random
import string
import subprocess
import pandas as pd
import shutil
import warnings
warnings.simplefilter('error', pd.errors.DtypeWarning)


## Declare global variables

DEBUG = False
VERB = False


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


data = {}
feature_list = {}

## FUNCTION TO LOAD RELEVANT REFERENCE DATA
def data_init(source_data_dir, function_to_load):
    '''Function to initialize file paths, functions to run, data to load.'''
    ## define global variables to use.
    global data
    global feature_list

    ## get file paths
    repeatmasker_file = os.path.join(source_data_dir, "rmsk_filtered.txt")
    jaspar_file = os.path.join(source_data_dir, "jaspar_filtered.txt")
    fimo_file = os.path.join(source_data_dir, "fimo_motifs.txt")

    ## load data and populate features

    func = function_to_load
    ## data
    if "rmsk" in func and "rmsk" not in data:
        ## initial load of rmsk data
        data["rmsk"] = [{"chr": x.split("\t")[0],
            "start": x.split("\t")[1],
            "end": x.split("\t")[2],
            "repName": x.split("\t")[5],
            "repClass": x.split("\t")[6],
            "repFamily": x.split("\t")[7]} for x in open(repeatmasker_file)]

    elif func == "jaspar":
        data[func] = [{"chr": x.split("\t")[0],
            "start": x.split("\t")[1],
            "end": x.split("\t")[2],
            "tfbs": x.rstrip().split("\t")[6]} for x in open(jaspar_file)]

    elif "fimo" in func and "fimo" not in data:
        data["fimo"] = fimo_file


    ## features and feature-specific data
    if func == "rmsk_repName":
        feature_list[func] = list(set([x["repName"] for x in data["rmsk"]]))
        ## further process this data, getting regions for these specific features
        data[func] = {feat: [] for feat in feature_list["rmsk_repName"]}
        for element in data["rmsk"]:
            data[func][element["repName"]].append("{}\t{}\t{}".format(element["chr"], element["start"], element["end"]))
        ## convert these to BedTools objects.
        for feat in data[func]:
            data[func][feat] = pybedtools.BedTool("\n".join(data[func][feat]), from_string=True)

    elif func == "rmsk_repClass":
        feature_list[func] = list(set([x["repClass"] for x in data["rmsk"]]))
        ## further process this data, getting regions for these specific features
        data[func] = {feat: [] for feat in feature_list["rmsk_repClass"]}
        for element in data["rmsk"]:
            data[func][element["repClass"]].append("{}\t{}\t{}".format(element["chr"], element["start"], element["end"]))
        ## convert these to BedTools objects.
        for feat in data[func]:
            data[func][feat] = pybedtools.BedTool("\n".join(data[func][feat]), from_string=True)

    elif func == "rmsk_repFamily":
        feature_list[func] = list(set([x["repFamily"] for x in data["rmsk"]]))
        ## further process this data, getting regions for these specific features
        data[func] = {feat: [] for feat in feature_list["rmsk_repFamily"]}
        for element in data["rmsk"]:
            data[func][element["repFamily"]].append("{}\t{}\t{}".format(element["chr"], element["start"], element["end"]))
        ## convert these to BedTools objects.
        for feat in data[func]:
            data[func][feat] = pybedtools.BedTool("\n".join(data[func][feat]), from_string=True)

    elif func == "jaspar":
        feature_list[func] = list(set([x["tfbs"] for x in data["jaspar"]]))
        ## further process this data, getting regions for these specific features
        data["jaspar_tfbs"] = {feat: [] for feat in feature_list[func]}
        for element in data["jaspar"]:
            data["jaspar_tfbs"][element["tfbs"]].append("{}\t{}\t{}".format(element["chr"], element["start"], element["end"]))
        ## convert these to BedTools objects.
        for feat in data["jaspar_tfbs"]:
            data["jaspar_tfbs"][feat] = pybedtools.BedTool("\n".join(data["jaspar_tfbs"][feat]), from_string=True)
    
    # elif func == "kmer":
    #     feature_list[func] = list(set([x for x in data["kmer"]]))

    else:
        ## no special feature list, just set as single-length array with func name.
        feature_list[func] = [func]


## FUNCTIONS TO GET THE POSITIONS OF FEATURES ON OUR OCRS

def rmsk_repName(region_coords, feature, *args, **kwargs):
    '''Find locations of regions that overlap a RMSK repName feature'''
    ## return the overlapping regions of regions that have any overlapping bases with this featureset.
    a = region_coords.intersect(data["rmsk_repName"][feature])
    intersects = str(a).rstrip()
    
    os.remove(a.fn)

    ## add feature name to these positions and strand(*)
    intersects = "\n".join(["{}\t{}\t*".format(x, feature) for x in intersects.split("\n")])
    
    return intersects


def rmsk_repClass(region_coords, feature, *args, **kwargs):
    '''Find locations of regions that overlap a RMSK repClass feature'''
    ## return the overlapping regions of regions that have any overlapping bases with this featureset.
    a = region_coords.intersect(data["rmsk_repClass"][feature])
    intersects = str(a).rstrip()
    
    os.remove(a.fn)

    ## add feature name to these positionsand strand(*)
    intersects = "\n".join(["{}\t{}\t*".format(x, feature) for x in intersects.split("\n")])
    
    return intersects

def rmsk_repFamily(region_coords, feature, *args, **kwargs):
    '''Find locations of regions that overlap a RMSK repFamily feature'''
    ## return the overlapping regions of regions that have any overlapping bases with this featureset.
    a = region_coords.intersect(data["rmsk_repFamily"][feature])
    intersects = str(a).rstrip()
    
    os.remove(a.fn)

    ## add feature name to these positions and strand(*)
    intersects = "\n".join(["{}\t{}\t*".format(x, feature) for x in intersects.split("\n")])
    
    return intersects


def jaspar(region_coords, feature, *args, **kwargs):
    '''Find locations of regions that overlap a JASPAR TFBS'''
    ## return the fraction of regions that have any overlapping bases with this featureset.
    a = region_coords.intersect(data["jaspar_tfbs"][feature])
    intersects = str(a).rstrip()
    
    os.remove(a.fn)

    ## add feature name to these positions and strand(*)
    intersects = "\n".join(["{}\t{}\t*".format(x, feature) for x in intersects.split("\n")])
    
    return intersects


def fimo(region_sequences, features, *args, **kwargs):
    '''Run FIMO to find known motifs in sets of OCR sequences, and report their locations.'''
    #################################################################################
    ##                        Set up fasta file and directories                    ##
    #################################################################################  
    jobname = "fimorun_{}".format(''.join(random.choice(string.ascii_lowercase) for i in range(10)))
    #fimo_outdir='/projects/nknoetze_prj/nknoetze_scratch/ocr_prj/results/{}'.format(jobname)
    outdir='/tmp/{}'.format(jobname)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    fimo_outdir=os.path.join(outdir,"fimo")
    

    motif_file = data["fimo"]
    
    #Write out the sequences to the a fasta file. Keep redundant sequences
    fasta_file = os.path.join(outdir, "{}.fasta".format(jobname))
    fasta_num = 0
    out = open(fasta_file, "w")
    for seq in region_sequences.split("\n"):
        if seq.startswith(">"):
            out.write("{}\n".format(seq))
            fasta_num += 1
        else:
            out.write("{}\n".format(seq))
    out.close()
    #################################################################################
    ##                            Run FIMO to Scan for motifs                      ##
    #################################################################################
    fimo_cmd='/projects/nknoetze_prj/miniconda3/envs/ocr_prj/bin/fimo --verbosity 1 --thresh 0.0001 --o {} {} {}'.format(fimo_outdir,motif_file, fasta_file)
    #Run FIMO on the sequences
    VALID_COM = False
    while not VALID_COM:
        call = subprocess.Popen(fimo_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        (res, err) = call.communicate()
        if err.decode('ascii') == "":
            VALID_COM = True
        else:
            print("ERROR IN FIMO: {}\nError: {}".format(fimo_cmd, str(err.decode('ascii'))))
            print("FASTA FILE:")
            for line in open(fasta_file, "r"):
                print(line)
            raise Exception()
    #################################################################################
    ##               Read in List of all Motifs Scanned by FIMO                    ##
    #################################################################################   
    ### GET LIST OF ALL MOTIFS.
    with open(motif_file,'r') as fi:
        motif_ids=[]
        motif_alt_ids=[]
        for line in fi:
            if line.startswith("MOTIF"):
                line=line.split('\n')[0]
                split_line=line.split(' ')
                #get list of motif ids and alternate motif ids
                motif_ids.append(split_line[1])
                motif_alt_ids.append(split_line[2])
                #get non-redundant alternate motif ids
                motif_alt_ids=list(set(motif_alt_ids))
    #################################################################################
    ##                        Read and Parse FIMO Results                          ##
    #################################################################################
    try:
        fimo_results=pd.read_csv('{}/fimo.tsv'.format(fimo_outdir),sep='\t',index_col=False, comment="#")
    except:
        print("ERROR")
        print(fimo_outdir)
        print("fasta:")
        for line in open(fasta_file, "r"):
            print(line)
        print("fimo:")
        for line in open('{}/fimo.tsv'.format(fimo_outdir), "r"):
            print(line)
        sys.exit(1)
        #filter tfbs predictions for those with the fimo score threshold.
    fimo_results=(fimo_results[~fimo_results.motif_id.str.startswith('#')].query('`score`>11.16'))
    
    #remove the temp files
    shutil.rmtree(outdir)

    ## subset to significant features
    fimo_results_significant = fimo_results[fimo_results["motif_alt_id"].isin(list(features))][["sequence_name", "start", "stop","motif_id", "motif_alt_id", "p-value", "q-value", "strand","matched_sequence","score"]].drop_duplicates()

    fimo_results_all = fimo_results.filter(items=["sequence_name", "start", "stop","motif_id", "motif_alt_id", "p-value", "q-value", "strand","matched_sequence","score"]).drop_duplicates()

    return fimo_results_significant.to_csv(index=False), fimo_results_all.to_csv(index=False)



def kmer(region_sequences, features, *args, **kwargs):
    '''Find locations of regions that match a given kmer'''
    kmer_lengths = set([len(k) for k in features])
    kmer_res = ""
    for fasta_seq in region_sequences.split("\n"):
        if fasta_seq.startswith(">"):
            ## get ocr_id
            ocr_id = fasta_seq[1:]
            chrom = ocr_id.split(":")[0]
            start = int(ocr_id.split(":")[1].split("-")[0])
            end = int(ocr_id.split(":")[1].split("-")[1])
        else:
            fasta_seq = fasta_seq.upper()
            for kl in kmer_lengths:
                for i in range(len(fasta_seq)-kl):
                    if fasta_seq[i:(i+kl)] in features:
                        # add strand to result file (*)
                        kmer_res += "{}\t{}\t{}\t{}\t{}\t*\n".format(chrom, start+i, start+i+kl, ocr_id, fasta_seq[i:(i+kl)])

    return kmer_res.rstrip()














## Main

FUNCTIONS_TO_RUN = {"rmsk_repName": rmsk_repName,
                    "rmsk_repClass": rmsk_repClass,
                    "rmsk_repFamily": rmsk_repFamily,
                    "jaspar": jaspar,
                    "kmer": kmer,
                    "kmer_geneCount": kmer,
                    "kmer_revcomp": kmer,
                    "kmer_revcomp_geneCount": kmer,
                    "fimo": fimo,
                    "fimo_geneCount": fimo,
                    "fimo_comotif_geneCount" : fimo,
                    "fimo_geneCount_shuffling": fimo,
                    "fimo_comotif_geneCount_shuffling" : fimo,
                    "fimo_geneCount_scrambling": fimo,
                    "fimo_comotif_geneCount_scrambling" : fimo
                    }

FEATURE_METRIC_FLIPPED_FUNCTIONS = ["fimo", "fimo_comotif_geneCount","fimo_geneCount", "kmer", "kmer_geneCount", "kmer_revcomp", "kmer_revcomp_geneCount","fimo_geneCount_shuffling","fimo_comotif_geneCount_shuffling","fimo_geneCount_scrambling","fimo_comotif_geneCount_scrambling"]

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = TITLE)
    parser.add_argument("--target_genes", dest = "TARGET_GENES", help = "Text file of target genes, ENSG, one per line.", type = str)
    parser.add_argument("--region_file", dest = "REGION_FILE", help = "File of regions linked to genes with additional info columns.", type = str)
    parser.add_argument("--result_file", dest = "RESULT_FILE", help = "File of random resampling results.", type = str)
    parser.add_argument("--pthresh", dest = "P_THRESH", help = "p-value < threshold", type = float)
    parser.add_argument("--fcthresh", dest = "FC_THRESH", help = "Fold change > threshold", type = float)
    parser.add_argument("--ocrmetrics_data_dir", dest = "OCRMET_DAT_DIR", help = "Directory where ocrmetrics data files can be found.", type = str)
    parser.add_argument("--datatype", dest = "DATATYPE", help = "Datatype to use from result file - this was the 'function' in framework run.", type = str, choices = FUNCTIONS_TO_RUN.keys())
    parser.add_argument("--metric", dest = "METRIC", help = "METRIC to use from result file.", type = str)
    parser.add_argument("--region_set", dest = "REGION_SET", help = "Region set/type (all, promoter, linked_promoter, non-promoter).", type = str)
    parser.add_argument("--refgenome_fasta", dest = "REF_FASTA", help = "Genome reference fasta file.", type = str)
    parser.add_argument("--expressed_features", dest = "EXP_FEAT", help = "File with expressed features to subset to.", type = str)
    parser.add_argument("--output", dest = "OUTPUT", help = "Output tsv file to write", type = str)
    parser.add_argument("--output_features", dest = "OUTPUT_FEATURES", help = "Output tsv file to write with feature instances", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    print("="*80)
    print("Python version: {}".format(sys.version))
    print("Python environment: {}".format(sys.prefix))
    print("Server: {}".format(os.uname()[1]))
    print("Current directory: {}".format(os.getcwd()))
    print("Command: {}".format(" ".join(sys.argv)))
    print("Time: {}".format(time.strftime("%Y/%m/%d %T")))
    print("{}\n".format("="*80))


    ## Script content here

    ## Initialize reference data
    data_init(args.OCRMET_DAT_DIR, args.DATATYPE)

    ## Extract OCRs for each gene, combine into unique set of OCRs.
    REF_GENOME_FASTA = pybedtools.BedTool(args.REF_FASTA)

    ## read in expressed, if provided
    exp_feats = set()
    if args.EXP_FEAT:
        for line in open(args.EXP_FEAT, "r"):
            exp_feats.add(line.rstrip().upper())

    ## Read in Target Gene List
    statprint("Reading in target genes...")

    target_genes = []
    for line in open(args.TARGET_GENES, "r"):
        target_genes.append(line.rstrip())

    ## Load in regions
    statprint("Loading all relevant regions into memory...")

    
    ## categorize regions by gene id
    regions = {}
    regions2genes = {}

    header = None
    for line in open(args.REGION_FILE, "r"):
        if not header:
            header = line.rstrip().split("\t")
        #if not line.startswith("chrom"):
        else:
            ## not header line
            # chrom   start   end     ocr_id  score   strand  ocr_type        query_gene      promoter_id     query_name      gene_group      data    dist_tss        promoter_name   promoter_type
            linedict = {x:y for x,y in zip(header, line.rstrip().split("\t"))}
            region_id = "{}\t{}".format(linedict["ocr_id"].replace(":","\t").replace("-","\t"), linedict["ocr_id"])
            oid = linedict["ocr_id"]
            ## region type - need to compare query_gene and promoter_id.
            promoter_id = linedict["promoter_id"]
            gene_id = linedict["query_gene"]
            region_type = linedict["ocr_type"]

            ## only use regions for genes of interest!
            if gene_id in target_genes:
                ## if promoter not of gene in question, we are labelling it as a linked promoter
                if region_type == "promoter" and gene_id != promoter_id:
                    region_type = "linked_promoter"
                    ## as of June 1, 2023 no longer considering linked promoters
                else:
                    if args.REGION_SET == "all":
                        ## we will track all regions and their region_type for this set. Priority of region_type is promoter > linked_promoter > non-promoter
                        ## any region will either be non-promoter, or promoter/linked_promoter. 
                        if region_id not in regions:
                            regions[region_id] = region_type
                        else:

                            ### WHY IS THIS ONLY CAPTURING NON-PROMOTERS?????
                            if regions[region_id] != "promoter":
                                ## update - could potentially be promoter now. Note this may result in no change.
                                regions[region_id] = region_type

                        if oid not in regions2genes:
                            regions2genes[oid] = set()
                        regions2genes[oid].add(gene_id)
                    
                    else:
                        ## only track regions of interest
                        if region_type == args.REGION_SET:
                            regions[region_id] = region_type
                        
                            if oid not in regions2genes:
                                regions2genes[oid] = set()
                            regions2genes[oid].add(gene_id)

    
    statprint("Converting regions to pybedtools object...")
    ## convert to pybedtools object
    bed_regions = pybedtools.BedTool("\n".join([r for r in regions]), from_string=True)

    statprint("Getting actual sequences of regions...")
    ## get sequences
    a = bed_regions.sequence(fi=REF_GENOME_FASTA)     ## pybedtools sequence (NOT temp file)
    region_seqs = open(a.seqfn).read()
    ## convert fasta to dict
    region_seqs_dict = {}
    for reg in region_seqs.rstrip().split("\n"):
        if reg.startswith(">"):
            seqname = reg[1:]
        else:
            region_seqs_dict[seqname] = reg
    
    ## Initialize per-base array to track composite score.
    statprint("Initializing score array...")
    scores = {}
    for reg in region_seqs_dict:
        scores[reg] = [0 for x in range(len(region_seqs_dict[reg])+1)]


    ## Load in significant features.
    statprint("Loading significant features...")
    sig_features = set()
    header = None
    for line in open(args.RESULT_FILE, "r"):
        if not header:
            header = line.rstrip().split("\t")
        else:
            ## not header line
            # datatype        feature metric  region_set      target_value    background_mean p_val  background_vals
            linedict = {x:y for x,y in zip(header, line.rstrip().split("\t"))}
            
            datatype = linedict["datatype"]
            feature = linedict["feature"]#.upper()
            metric = linedict["metric"]
            region_set = linedict["region_set"]
            fc = (float(linedict["target_value"])+0.000001) / (float(linedict["background_mean"])+0.000001)
            pval = float(linedict["p_val"])

            #if feature == "MIRb" and region_set == args.REGION_SET and metric == args.METRIC and datatype == args.DATATYPE: print("\n".join([datatype, feature, metric, region_set, str(fc), str(pval)]))
            if datatype == args.DATATYPE and metric == args.METRIC and region_set == args.REGION_SET and pval < args.P_THRESH and fc > args.FC_THRESH:
                ## this is a relevant result.
                if not args.EXP_FEAT:
                    ## no expressed features provided, write all.
                    sig_features.add(feature)
                elif ":" not in feature and feature in exp_feats:
                    ## FIMO feature is straightforward, and is in expressed feature list
                    sig_features.add(feature)
                elif ":" in feature:
                    feats = feature.split("::")
                    if feats[0] in exp_feats and feats[1] in exp_feats:
                        ## both parts of dimer expressed
                        sig_features.add(feature)


    if DEBUG: statprint("sig_features are: {}...".format(sig_features))
    if len(sig_features) == 0:
        statprint("No significant features in result set", "WARNING")

        out = open(args.OUTPUT, "w")
        out_feat = open(args.OUTPUT_FEATURES, "w")
        if "fimo" in args.DATATYPE:
            out.write("datatype\tocr_id\tposition\tscore\tstrand\tgenes\tmatched_sequence\tmotif_id\tq_val\tscore\n")
            out_feat.write("datatype\tfeature\tocr_id\tstart\tend\tstrand\tgenes\tmatched_sequence\tmotif_id\tq_val\tscore\n")
            for reg in regions:
                oid = reg.split("\t")[3]
                out.write("{}\t{}\tNA\tNA\tNA\t{}\tNA\tNA\tNA\tNA\n".format(args.DATATYPE, oid, ",".join([x for x in regions2genes[oid]])))
                out_feat.write("{}\tNA\t{}\tNA\tNA\tNA\t{}\tNA\tNA\tNA\tNA\n".format(args.DATATYPE, oid, ",".join([x for x in regions2genes[oid]])))

        else:
            out.write("datatype\tocr_id\tposition\tscore\tstrand\tgenes\n")
            out_feat.write("datatype\tfeature\tocr_id\tstart\tend\tstrand\tgenes\n")
            for reg in regions:
                oid = reg.split("\t")[3]
                out.write("{}\t{}\tNA\tNA\tNA\t{}\n".format(args.DATATYPE, oid, ",".join([x for x in regions2genes[oid]])))
                out_feat.write("{}\tNA\t{}\tNA\tNA\tNA\t{}\n".format(args.DATATYPE, oid, ",".join([x for x in regions2genes[oid]])))
        out.close()
        out_feat.close()


    else:
        ## Identify positions that features map to on this OCR set.
        statprint("Finding where these features occur...")
        # run the function, for all, or the specific feature.
        if args.DATATYPE in FEATURE_METRIC_FLIPPED_FUNCTIONS:
            if "fimo" in args.DATATYPE:
            ## these functions do not receive a specific feature to run the analysis. They can receive 'features' to subset result before returning.
            # fimo will output the superset and filtered subset of results.
                feat_res, feat_res_all = FUNCTIONS_TO_RUN[args.DATATYPE](region_coords = bed_regions,
                                                region_sequences = region_seqs,
                                                features = sig_features)
                #if "fimo" in args.DATATYPE:
                #### I THINK CAN BE REMOVED? ALL OF THESE OPERATIONS COULD BE DONE IN THE SECTION BELOW, ESPECIALLY SINCE IT NEEDS A NEW
                ## RETURNS ALL RESULTS, NOT JUST SIGNIFICANT
                ## what is returned is chr, start, stop, ocr_id,motif_alt_id, STRAND and *matched_sequence*,motif_id for all sig features.
                # feat_res = "\n".join([
                #     "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                #         x.split(",")[0].split(":")[0], #chromosome
                
                #         #some functions return intersected bed files, so parse positions from from ocr id positions, and the
                #         #start and end columns returned from the intersection 
                #         int(x.split(",")[1]) + int(x.split(",")[0].split(":")[1].split("-")[0]) - 1,#start and end
                #         int(x.split(",")[2]) + int(x.split(",")[0].split(":")[1].split("-")[0]) - 1,
                #         x.split(",")[0], #ocr id
                #         x.split(",")[4],#motif_alt_id
                #         x.split(",")[7],#strand
                #         x.split(",")[8],#matched sequence
                #         x.split(",")[3]#motif id
                #     ) for x in feat_res_all.rstrip().split("\n") if not x.startswith("sequence_name")
                # ])
            #DO WE EVEN NEED THE FILTERED RESULTS????
                feat_res=feat_res_all
            else:
                feat_res = FUNCTIONS_TO_RUN[args.DATATYPE](region_coords = bed_regions,
                                                region_sequences = region_seqs,
                                                features = sig_features)


        else:
            ## these functions receive the feature, so return only the feature-specific results.
            feat_res_list = []
            for feat in sig_features:
                feat_res_list.append(FUNCTIONS_TO_RUN[args.DATATYPE](region_coords = bed_regions,
                                                        region_sequences = region_seqs,
                                                        feature = feat))
                # a bed file string with chrom start end, name and strand is returned with the overlapping regions
            ## merge
            feat_res = "\n".join(feat_res_list)




        ## walk through positions and update score array

        ###NEED TO ADD STRAND FOR SCORES
        statprint("Updating score array with these locations...and writing feature instances...")
        seen_ocrs = set()
        out = open(args.OUTPUT_FEATURES, "w")
        if "fimo" in args.DATATYPE:
            out.write("datatype\tfeature\tocr_id\tstart\tend\tstrand\tgenes\tmatched_sequence\tmotif_id\tq_val\tscore\n")
            for feat in feat_res.split("\n"):
                if feat.startswith("chr"):
                    feat = feat.split(",")
                    
                    ocr_id = feat[0]
                    seen_ocrs.add(ocr_id)

                    feature = feat[4]
                    start=int(feat[1])
                    end = int(feat[2])
                    strand=feat[7]
                    matched_sequence=feat[8]
                    motif_id=feat[3]
                    q_val=feat[6]
                    score=feat[9]

                    out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(args.DATATYPE, feature, ocr_id, start, end, strand,",".join([x for x in regions2genes[ocr_id]]),matched_sequence,motif_id,q_val,score))

                    ###SCOTT DOUBLE CHECK IF YOU NEED TO ADD 1 HERE?????
                    for idx in range(start, end+1):
                        scores[ocr_id][idx] += 1
            # write empty lines for OCRs that don't have any features.
            for ocr_id in regions:
                if ocr_id not in seen_ocrs:
                    out.write("{}\tNA\t{}\tNA\tNA\tNA\t{}\tNA\tNA\tNA\tNA\n".format(args.DATATYPE, ocr_id.split("\t")[3], ",".join([x for x in regions2genes[ocr_id.split("\t")[3]]])))
            out.close()
                
        else:
            out.write("datatype\tfeature\tocr_id\tstart\tend\tstrand\tgenes\n")
            for feat in feat_res.split("\n"):
                feat = feat.split("\t")
                #print(feat, flush=True)
                ocr_id = feat[3]
                seen_ocrs.add(ocr_id)

                feature = feat[4]

                ocr_start = int(ocr_id.split(":")[1].split("-")[0])
                #DO WE NEED TO ADJUST FOR -1????
                start = int(feat[1])
                end = int(feat[2])
                strand=feat[5]

                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(args.DATATYPE, feature, ocr_id, start-ocr_start, end-ocr_start, strand,",".join([x for x in regions2genes[ocr_id]])))
                ###DOUBLE CHECK RESULTS HERE????
                for idx in range(start-ocr_start, end-ocr_start+1):
                    scores[ocr_id][idx] += 1
            # write empty lines for OCRs that don't have any features.
            for ocr_id in regions:
                if ocr_id not in seen_ocrs:
                    out.write("{}\tNA\t{}\tNA\tNA\tNA\t{}\n".format(args.DATATYPE, ocr_id.split("\t")[3], ",".join([x for x in regions2genes[ocr_id.split("\t")[3]]])))
            out.close()

        ## write output
        ##NEEED TO ADD STRAND?????
        statprint("Writing output...")
        out = open(args.OUTPUT, "w")
        out.write("datatype\tocr_id\tposition\tscore\tgenes\n")
        for ocr_id in scores:
            for idx in range(len(scores[ocr_id])):
                out.write("{}\t{}\t{}\t{}\t{}\n".format(args.DATATYPE, ocr_id, idx, scores[ocr_id][idx], ",".join([x for x in regions2genes[ocr_id]])))
        out.close()




    

    statprint("Done.")