TITLE = "OCR Metrics"
DESC = """Functions to test OCR metrics."""

'''
Date: February 06, 2023
@authors: sbrown, nknoetze
'''

import sys
import os
import time
import pybedtools
from collections import Counter
from collections import defaultdict
import subprocess
import pandas as pd
import shutil
import warnings
import itertools
warnings.simplefilter('error', pd.errors.DtypeWarning)
import gc
import string
import random
import numpy as np
#from Bio import SeqIO 

#pybedtools.set_tempdir('/projects/sbrown_prj/sbrown_scratch/tmp')


data = {}
feature_list = {}

COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
def revcomp(seq):
    return ''.join([COMPLEMENT[nt] for nt in seq[::-1]])

def init(source_data_dir, functions_to_load,num_motifs):
    '''Function to initialize file paths, functions to run, data to load.'''
    ## define global variables to use.
    global data
    global feature_list

    ## get file paths
    #repeatmasker_file = os.path.join(source_data_dir, "rmsk_filtered.txt")
    #ocr_group_specific_file = os.path.join(source_data_dir, "group_specific_ocrs.tsv")
    #ocr_single_gene_specific_file = os.path.join(source_data_dir, "gene_specific_ocrs.tsv")
    jaspar_file = os.path.join(source_data_dir, "jaspar_filtered.txt")
    #kmer_file = os.path.join(source_data_dir, "kmers.txt")
    fimo_file = os.path.join(source_data_dir, "fimo_motifs.txt")
    ##not actually used...... may delete later lols.
    #filtered_fimo_file = os.path.join(source_data_dir, "fimo_motifs_expressed.txt")

    ## load data and populate features

    for func in functions_to_load:
        ## data
        if "rmsk" in func and "rmsk" not in data:
            ## initial load of rmsk data
            print('reading in repeat data')
            # initial load of rmsk data
            ## load data and populate features
            data["rmsk"] = pd.read_csv(repeatmasker_file,sep='\t',index_col=False,usecols=[0,1,2,5,6,7],names=['chr','start','end','.','strand','repName','repClass','repFamily','1','2','3','4'])

        elif func == "ocrs_group_specific":
            data[func] = pybedtools.BedTool("\n".join(
                list(
                    set(
                        ["\t".join(x.split("\t")[0:3]) for x in open(ocr_group_specific_file) if not x.startswith("chrom")]
                    )
                )
            ), from_string=True
            )
            
        elif func == "ocrs_single_gene_specific":
            data[func] = pybedtools.BedTool("\n".join(
                list(
                    set(
                        ["\t".join(x.split("\t")[0:3]) for x in open(ocr_single_gene_specific_file) if not x.startswith("chrom")]
                    )
                )
            ), from_string=True
            )

        elif func == "jaspar":
            data[func] = [{"chr": x.split("\t")[0],
              "start": x.split("\t")[1],
              "end": x.split("\t")[2],
              "tfbs": x.rstrip().split("\t")[6]} for x in open(jaspar_file)]
        
        elif "kmer" in func and "kmer" not in data:
            data["kmer"] = [x.rstrip() for x in open(kmer_file)]

        #if function is fimo co motifs, use filtered meme file. Otherwise regular.
        # elif func=="fimo_comotif_geneCount":
        #     if max(num_motifs)>2:
        #         data["fimo_comotif"] = filtered_fimo_file
        #     else:
        #         data["fimo_comotif"] = fimo_file
        # elif func=="fimo" or func =="fimo_geneCount" and "fimo" not in data:
        #     data["fimo"] = fimo_file
        elif "fimo" in func and "fimo" not in data:
            data["fimo"] = fimo_file

        ## features and feature-specific data
        if func == "rmsk_repName":
            feature_list[func]=data["rmsk"]['repName'].unique()
            #create a dictionary for each repeat class. each key (repeatName) will contain an empty list.
            data[func]={feat: [] for feat in feature_list[func]}
            #loop through the feature list, which contains the names of the unique repeat classes
            # returns a dataframe of coordinates for each repeat class
            for element in feature_list[func]:
                repeat_coords=data['rmsk'].loc[data['rmsk']['repName']==element].filter(items=['chr','start','end']).drop_duplicates()
                pybed=pybedtools.BedTool.from_dataframe(repeat_coords)
                data[func][element]=pybed

        elif func == "rmsk_repClass":
            feature_list[func]=data["rmsk"]['repClass'].unique()
            #create a dictionary for each repeat class. each key (repeatClass) will contain an empty list.
            data[func]={feat: [] for feat in feature_list[func]}
            #loop through the feature list, which contains the names of the unique repeat classes
            # returns a dataframe of coordinates for each repeat class
            for element in feature_list[func]:
                repeat_coords=data['rmsk'].loc[data['rmsk']['repClass']==element].filter(items=['chr','start','end']).drop_duplicates()
                pybed=pybedtools.BedTool.from_dataframe(repeat_coords)
                data[func][element]=pybed

        elif func == "rmsk_repFamily":
            feature_list[func]=data["rmsk"]['repFamily'].unique()
            #create a dictionary for each repeat class. each key (repeatFamily) will contain an empty list.
            data[func]={feat: [] for feat in feature_list[func]}
            #loop through the feature list, which contains the names of the unique repeat classes
            # returns a dataframe of coordinates for each repeat class
            for element in feature_list[func]:
                repeat_coords=data['rmsk'].loc[data['rmsk']['repFamily']==element].filter(items=['chr','start','end']).drop_duplicates()
                pybed=pybedtools.BedTool.from_dataframe(repeat_coords)
                data[func][element]=pybed

        elif func == "jaspar":
            feature_list[func] = list(set([x["tfbs"] for x in data["jaspar"]]))
            ## further process this data, getting regions for these specific features
            data["jaspar_tfbs"] = {feat: [] for feat in feature_list[func]}
            for element in data["jaspar"]:
                data["jaspar_tfbs"][element["tfbs"]].append("{}\t{}\t{}".format(element["chr"], element["start"], element["end"]))
            ## convert these to BedTools objects.
            for feat in data["jaspar_tfbs"]:
                data["jaspar_tfbs"][feat] = pybedtools.BedTool("\n".join(data["jaspar_tfbs"][feat]), from_string=True)
        
        else:
            ## no special feature list, just set as single-length array with func name.
            feature_list[func] = [func]


    ## return feature_list (instead of relying on class features)
    return feature_list




def gcContent(region_sequences, *args, **kwargs):
    '''Calculate fraction of bases that are G/C in input sequences (as a whole) and return this value'''
    total_bases = 0
    g_count = 0
    c_count = 0

    for seq in region_sequences.rstrip().split("\n"):
        if not seq.startswith(">"):
            total_bases += len(seq)
            g_count += seq.upper().count("G")
            c_count += seq.upper().count("C")
    
    return {"gc_frac": (g_count + c_count) / total_bases}


def rmsk_repName(region_coords, feature, *args, **kwargs):
    '''Calculate metrics for regions that overlap a RMSK repName feature'''
    ## return the fraction of regions that have any overlapping bases with this featureset.
    a = region_coords.window(data["rmsk_repName"][feature], w=0)

    ## get the unique set of regions that overlap
    overlapping_regions_set = set(["_".join(x.split("\t")[0:3]) for x in str(a).rstrip().split("\n")])
    ## Get count of each region in the input file (the number of times that region occurs in the region file)
    regions_counter = Counter("_".join(x.split("\t")[0:3]) for x in str(region_coords).rstrip().split("\n"))
    ## get corrected count of OCRs that overlap one or more feature regions (correcting count for number of occurrences of OCR in region file)
    overlapped_region_count = sum(regions_counter[r] for r in overlapping_regions_set)

    count_frac = overlapped_region_count / region_coords.count()

    ## return the fraction of bases in regions overlapping this featureset.
    ## note that .overlap() using start and end of regions, ignores chromosome. But .window() only pairs overlapping regions based on chromosome.
    # >>> a = pybedtools.BedTool("chr1\t100\t200\nchr1\t500\t600\nchr1\t800\t900", from_string=True)
    # >>> b = pybedtools.BedTool("chr1\t201\t301\nchr1\t550\t650\nchr1\t888\t9000", from_string=True)
    # >>> for x in a.window(b, w=0).overlap(cols=[2,3,5,6]):
    # ...     print("Interval is {} bases long, {} overlap with feature".format(len(x), x.fields[-1]))
    # ...
    # Interval is 100 bases long, 50 overlap with feature
    # Interval is 100 bases long, 12 overlap with feature

    overlap_data = a.overlap(cols=[2,3,6,7])    ## doing this stepwise to maintain the pointers to all temp files
    total_region_size = 0
    total_overlap_size = 0
    for overlap in overlap_data:
        total_region_size += len(overlap) # this is the length of the region inverval, not the overlap specifically
        total_overlap_size += int(overlap.fields[-1])
    
    # calculate fraction of total region size (bases) that contains feature being tested.
    if total_region_size > 0:
        overlap_frac = total_overlap_size / total_region_size
    else:
        overlap_frac = 0

    
    ## count of genes with feature
    gene_set = set()
    for hit in a:
        gene_set.add(hit.fields[3])
    
    gene_count = len(gene_set)

    os.remove(a.fn)
    os.remove(overlap_data.fn)

    return {"count_frac": count_frac, "overlap_frac": overlap_frac, "gene_count": gene_count}


def rmsk_repClass(region_coords, feature, *args, **kwargs):
    '''Calculate metrics for regions that overlap a RMSK repClass feature'''
    ## return the fraction of regions that have any overlapping bases with this featureset.
    a = region_coords.window(data["rmsk_repClass"][feature], w=0)
    ## get the unique set of regions that overlap
    overlapping_regions_set = set(["_".join(x.split("\t")[0:3]) for x in str(a).rstrip().split("\n")])
    ## Get count of each region in the input file (the number of times that region occurs in the region file)
    regions_counter = Counter("_".join(x.split("\t")[0:3]) for x in str(region_coords).rstrip().split("\n"))
    ## get corrected count of OCRs that overlap one or more feature regions (correcting count for number of occurrences of OCR in region file)
    overlapped_region_count = sum(regions_counter[r] for r in overlapping_regions_set)
    count_frac = overlapped_region_count / region_coords.count()

    ## return the fraction of bases in regions overlapping this featureset.
    overlap_data = a.overlap(cols=[2,3,6,7])    ## doing this stepwise to maintain the pointers to all temp files
    total_region_size = 0
    total_overlap_size = 0
    for overlap in overlap_data:
        total_region_size += len(overlap) # this is the length of the region inverval, not the overlap specifically
        total_overlap_size += int(overlap.fields[-1])

    # calculate fraction of total region size (bases) that contains feature being tested.
    if total_region_size > 0:
        overlap_frac = total_overlap_size / total_region_size
    else:
        overlap_frac = 0

    ## count of genes with feature
    gene_set = set()
    for hit in a:
        gene_set.add(hit.fields[3])
    
    gene_count = len(gene_set)

    os.remove(a.fn)
    os.remove(overlap_data.fn)

    return {"count_frac": count_frac, "overlap_frac": overlap_frac, "gene_count": gene_count}

def rmsk_repFamily(region_coords, feature, *args, **kwargs):
    '''Calculate metrics for regions that overlap a RMSK repFamily feature'''
    ## return the fraction of regions that have any overlapping bases with this featureset.
    a = region_coords.window(data["rmsk_repFamily"][feature], w=0)
    ## get the unique set of regions that overlap
    overlapping_regions_set = set(["_".join(x.split("\t")[0:3]) for x in str(a).rstrip().split("\n")])
    ## Get count of each region in the input file (the number of times that region occurs in the region file)
    regions_counter = Counter("_".join(x.split("\t")[0:3]) for x in str(region_coords).rstrip().split("\n"))
    ## get corrected count of OCRs that overlap one or more feature regions (correcting count for number of occurrences of OCR in region file)
    overlapped_region_count = sum(regions_counter[r] for r in overlapping_regions_set)
    count_frac = overlapped_region_count / region_coords.count()

    ## return the fraction of bases in regions overlapping this featureset.
    overlap_data = a.overlap(cols=[2,3,6,7])    ## doing this stepwise to maintain the pointers to all temp files
    total_region_size = 0
    total_overlap_size = 0
    for overlap in overlap_data:
        total_region_size += len(overlap) # this is the length of the region inverval, not the overlap specifically
        total_overlap_size += int(overlap.fields[-1])

    # calculate fraction of total region size (bases) that contains feature being tested.
    if total_region_size > 0:
        overlap_frac = total_overlap_size / total_region_size
    else:
        overlap_frac = 0

    ## count of genes with feature
    gene_set = set()
    for hit in a:
        gene_set.add(hit.fields[3])
    
    gene_count = len(gene_set)

    os.remove(a.fn)
    os.remove(overlap_data.fn)

    return {"count_frac": count_frac, "overlap_frac": overlap_frac, "gene_count": gene_count}

def ocrs_number(region_coords, *args, **kwargs):
    '''Number of regions linked to gene list'''
    return {"count": region_coords.count()}

def ocrs_group_specific(region_coords, *args, **kwargs):
    '''Number of non-shared regions linked to group'''
    a = region_coords.window(data["ocrs_group_specific"], w=0)
    val = a.count()
    os.remove(a.fn)
    return {"count": val}

def ocrs_single_gene_specific(region_coords, *args, **kwargs):
    '''Number of non-shared OCRs linked to gene list'''
    a = region_coords.window(data["ocrs_single_gene_specific"], w=0)
    val = a.count()
    os.remove(a.fn)
    return {"count": val}

def ocrs_set_specific(set_genes, regions_genes, *args, **kwargs):
    '''Number of OCRs that are only linked to genes in this set'''
    return {"count": sum([all(y in set_genes for y in x) for x in regions_genes])}


## OCR-Gene distances
def ocrs_distances(region_gene_dists, *args, **kwargs):
    '''summed and averaged distance and displacements'''
    res = {}
    res["summed-distance"] = sum([abs(x) for x in region_gene_dists])
    res["summed-displacement"] = sum(region_gene_dists)
    res["average-distance"] = sum([abs(x) for x in region_gene_dists]) / len(region_gene_dists)
    res["average-displacement"] = sum(region_gene_dists) / len(region_gene_dists)
    return res

## BLAT
def blat_similarity(region_sequences, jobname, *args, **kwargs):
    '''Run BLAT for the set of sequences and calculate sequence identity'''
    
    ###Need to get the region sequences, and write these to an output file that we can pass to BLAT.##
    
    jobname = "-".join(jobname.split("\t"))   ## switch tabs for dashes.

    outdir='/tmp/'
    blat_outfile=os.path.join(outdir,'blat.{}.pslx'.format(jobname))
    fasta_file = os.path.join(outdir, "{}.fasta".format(jobname))
    fasta_seqs = set()
    out = open(fasta_file, "w")
    ## reminder that region_sequences is a string of a fasta file. fasta headers are >gene::chr:start-end
    ## for checking for unique sequences, we want to check that the sequence has not been seen before, which is also encoded by the header after "::"
    ## new edit: if sequence is seen before, but genomic coords are unique, still write it.
    for seq in region_sequences.split("\n"):
        if seq.startswith(">"):
            if seq.split("::")[1] not in fasta_seqs:
                out.write("{}\n".format(seq))
                fasta_seqs.add(seq.split("::")[1])
        #elif seq not in fasta_seqs:
        else:
            out.write("{}\n".format(seq))
            #fasta_seqs.add(seq)
    out.close()

    #1. JOBNAME AS OUTPUT FILE
    #Use the same set of sequences for the database and the query.
    blat_cmd='/projects/nknoetze_prj/miniconda3/envs/ocr_prj/bin/pblat -noHead -threads=1 -minScore=0 -out=pslx -stepSize=1 -minIdentity=0 -tileSize=6 {} {} {}'.format(fasta_file,fasta_file,blat_outfile)
    #Run blat on the sequences
    VALID_COM = False
    while not VALID_COM:
        call = subprocess.Popen(blat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        (res, err) = call.communicate()
        if err.decode('ascii') == "":
            VALID_COM = True
        else:
            print("ERROR IN BLAT: {}\nError: {}".format(blat_cmd, str(err.decode('ascii'))))
            raise Exception()


    ## Read in the BLAT Results.
    ## Create column containing the total bp that match, and the total length of the query and target sequence
    #remove rows where the query and the target regions are the same
    #calculate the query and sequence identity for each comparison
    tidy_blat=(pd.read_csv(blat_outfile, sep='\t',header=None, names=['matches','mismatches','repmatches','ncount','qgapcount','qgapbases',
            'tgapcount','tgapbases','strand','qname','qsize','qstart','qend','tname','tsize',
            'tstart','tend','blockcount','blocksizes','qstarts','tstarts','seq1','seq2'],index_col=False)
              .query('qname!=tname')
              .assign(total_match=lambda x: (x['matches']+x['repmatches']),
                      total_length=lambda x: (x['qsize']+x['tsize']),
                      comparison=lambda x: (x['qname']+';'+x['tname']),
                      query_identity=lambda x: (x['total_match']/x['qsize']),
                      target_identity=lambda x: (x['total_match'])/x['tsize'])
          .filter(items=['comparison','total_match','total_length','target_identity','query_identity'])
          )

    # for each ocr comparison, get the BEST alignment (alignment with the highest sequence identities)
    # Also calculate the comparison identity, which captures the total number of matching bases for
    # the target AND database sequence, and the total length of both sequences.
    ### make sure you times total_matches by TWO because this number bases match in EACH sequence ###
    blat_metrics=(tidy_blat.assign(max_target_identity=(tidy_blat.groupby('comparison')['target_identity'].transform(lambda x: x.max())),
                               max_query_identity=(tidy_blat.groupby('comparison')['query_identity'].transform(lambda x: x.max())))
              .drop_duplicates(subset=['comparison','max_target_identity','max_query_identity']).assign(comparison_identity=lambda x: ((x['total_match']*2)/x['total_length']))
              .filter(items=['comparison','max_target_identity','max_query_identity','comparison_identity'])
             )
        
    #Calculate total number of possible OCR comparisons
    total_ocrs=len(fasta_seqs)
    total_ocr_comparisons = total_ocrs * (total_ocrs-1)
    #summarises the total number of bases matching across all comparisions / the total numberooked at across all comparisons / the total 
    avg_comparison_identity=sum(blat_metrics['comparison_identity'])/total_ocr_comparisons
    avg_overall_identity=sum(blat_metrics['max_query_identity']+blat_metrics['max_target_identity'])/(total_ocr_comparisons*2)
                             
    #remove the blat file
    os.remove(blat_outfile)
    os.remove(fasta_file)
    
    return {"average_comparison_identity": avg_comparison_identity, "average_overall_identity": avg_overall_identity}


def jaspar(region_coords, feature, *args, **kwargs):
    '''Calculate metrics for regions that overlap a JASPAR TFBS'''
    ## return the fraction of regions that have any overlapping bases with this featureset.
    a = region_coords.window(data["jaspar_tfbs"][feature], w=0)
    ## get the unique set of regions that overlap
    overlapping_regions_set = set(["_".join(x.split("\t")[0:3]) for x in str(a).rstrip().split("\n")])
    ## Get count of each region in the input file (the number of times that region occurs in the region file)
    regions_counter = Counter("_".join(x.split("\t")[0:3]) for x in str(region_coords).rstrip().split("\n"))
    ## get corrected count of OCRs that overlap one or more feature regions (correcting count for number of occurrences of OCR in region file)
    overlapped_region_count = sum(regions_counter[r] for r in overlapping_regions_set)
    count_frac = overlapped_region_count / region_coords.count()

    ## return the fraction of bases in regions overlapping this featureset.
    overlap_data = a.overlap(cols=[2,3,6,7])    ## doing this stepwise to maintain the pointers to all temp files
    total_region_size = 0
    total_overlap_size = 0
    for overlap in overlap_data:
        total_region_size += len(overlap) # this is the length of the region inverval, not the overlap specifically
        total_overlap_size += int(overlap.fields[-1])

    # calculate fraction of total region size (bases) that contains feature being tested.
    if total_region_size > 0:
        overlap_frac = total_overlap_size / total_region_size
    else:
        overlap_frac = 0

    ## count of genes with feature
    gene_set = set()
    for hit in a:
        gene_set.add(hit.fields[3])
    
    gene_count = len(gene_set)

    os.remove(a.fn)
    os.remove(overlap_data.fn)

    return {"count_frac": count_frac, "overlap_frac": overlap_frac, "gene_count": gene_count}


def kmer(region_sequences, *args, **kwargs):
    '''Return the metrics for kmer counting, counting all kmers at once'''
    ## this is different from the previous method because instead of looking for a specific kmer
    ## it scans through the input sequences and tallies each kmer it encounters,
    ## returning counts for all kmers at once (treating each kmer as a "metric")
    ## note that the kmer file can contain sequences of any length, but if a length is used,
    ## all possible kmers of that length that exist in the dataset must be included.
    
    #kmer_dict = {k: 0 for k in data["kmer"]}
    ## lets only track kmers in this set of region_sequences.
    kmer_dict = {}
    kmer_lengths = set([len(k) for k in data["kmer"]])
    num_ocrs = 0
    region_sequences = region_sequences.upper()
    for fasta_seq in region_sequences.split("\n"):
        if not fasta_seq.startswith(">"):
            num_ocrs += 1
            for kl in kmer_lengths:
                for i in range(len(fasta_seq)-kl):
                    if fasta_seq[i:(i+kl)] in kmer_dict:
                        kmer_dict[fasta_seq[i:(i+kl)]] += 1
                    else:
                        kmer_dict[fasta_seq[i:(i+kl)]] = 1

    return {k:kmer_dict[k]/num_ocrs for k in kmer_dict}

def kmer_geneCount(region_sequences, *args, **kwargs):
    '''Return the count of genes that each kmer is found in, counting all kmers that exist in this data at once'''
    kmer_dict = {}
    kmer_lengths = set([len(k) for k in data["kmer"]])
    region_sequences = region_sequences.upper()
    for fasta_seq in region_sequences.split("\n"):
        if fasta_seq.startswith(">"):
            ## get gene name
            gene = fasta_seq.split("::")[0][1:]

        else:
            for kl in kmer_lengths:
                for i in range(len(fasta_seq)-kl):
                    if fasta_seq[i:(i+kl)] in kmer_dict:
                        kmer_dict[fasta_seq[i:(i+kl)]].add(gene)
                    else:
                        kmer_dict[fasta_seq[i:(i+kl)]] = set([gene])

    return {k:len(kmer_dict[k]) for k in kmer_dict}

def kmer_revcomp(region_sequences, *args, **kwargs):
    '''Return the metrics for kmer counting, counting all kmers at once'''
    ## for every kmer encountered, the "min" of it and its reverse complement is used as the key,
    ## and both are used for enumeration.
    
    kmer_dict = {}
    kmer_lengths = set([len(k) for k in data["kmer"]])
    num_ocrs = 0
    region_sequences = region_sequences.upper()
    for fasta_seq in region_sequences.split("\n"):
        if not fasta_seq.startswith(">"):
            num_ocrs += 1
            for kl in kmer_lengths:
                for i in range(len(fasta_seq)-kl):
                    word = fasta_seq[i:(i+kl)]
                    revcomp_word = revcomp(word)
                    key = min(word, revcomp_word)
                    if key in kmer_dict:
                        kmer_dict[key] += 1
                    else:
                        kmer_dict[key] = 1

    return {k:kmer_dict[k]/num_ocrs for k in kmer_dict}

def kmer_revcomp_geneCount(region_sequences, *args, **kwargs):
    '''Return the count of genes that each kmer is found in, counting all kmers that exist in this data at once'''
    kmer_dict = {}
    kmer_lengths = set([len(k) for k in data["kmer"]])
    region_sequences = region_sequences.upper()
    for fasta_seq in region_sequences.split("\n"):
        if fasta_seq.startswith(">"):
            ## get gene name
            gene = fasta_seq.split("::")[0][1:]

        else:
            for kl in kmer_lengths:
                for i in range(len(fasta_seq)-kl):
                    word = fasta_seq[i:(i+kl)]
                    revcomp_word = revcomp(word)
                    key = min(word, revcomp_word)
                    if key in kmer_dict:
                        kmer_dict[key].add(gene)
                    else:
                        kmer_dict[key] = set([gene])

    return {k:len(kmer_dict[k]) for k in kmer_dict}


def fimo(region_sequences, jobname, *args, **kwargs):
    '''Run FIMO to find known motifs in sets of OCR sequences'''
    #################################################################################
    ##                        Set up fasta file and directories                    ##
    #################################################################################  
    jobname = "-".join(jobname.split("\t") + [''.join(random.choice(string.ascii_lowercase) for i in range(10))])   ## switch tabs for dashes, add random 10character string.
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
            out.write(">{}\n".format(fasta_num)) # doesn't it already have > in the header?
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
    #filtering for TFBS predictions with a score greater than 11.16
    fimo_results=(fimo_results[~fimo_results.motif_id.str.startswith('#')].query('`score`>11.16'))
    
    #for motif_ids, calculate the number of occurences
    fimo_motif_counts=(fimo_results.groupby(['motif_id'])['motif_id'].count().reset_index(name='n_occurrence'))
    #for motif_alt_ids, calculate the number of occurences
    fimo_motif_alt_counts=(fimo_results.groupby(['motif_alt_id'])['motif_alt_id'].count().reset_index(name='n_occurrence'))


    #Create a dictionary of each motif from the motif file
    #Add the n_occurence for each motif
    #If the motif is not found in the sequences, add the motif to the dictionary with a n_occurrence of 0
    motif_id_dict=dict()
    for mo_id in motif_ids:
        if mo_id in fimo_motif_counts['motif_id'].unique():
            fimo_df=fimo_motif_counts.loc[fimo_motif_counts['motif_id']==mo_id]
            motif_id_dict[mo_id]=fimo_df['n_occurrence'].values[0]
        else:
            motif_id_dict[mo_id]=0

    motif_alt_id_dict=dict()
    for mo_alt_id in motif_alt_ids:
        if mo_alt_id in fimo_motif_alt_counts['motif_alt_id'].unique():
            fimo_df=fimo_motif_alt_counts.loc[fimo_motif_alt_counts['motif_alt_id']==mo_alt_id]
            motif_alt_id_dict[mo_alt_id]=fimo_df['n_occurrence'].values[0]
        else:
            motif_alt_id_dict[mo_alt_id]=0

    
    #remove the temp files
    shutil.rmtree(outdir)
    #os.remove('{}/fimo.tsv'.format(fimo_outdir))
    #os.remove(fasta_file)
    
    #return {"alt_motif_occurrence": motif_alt_id_dict, "motif_occurrence": motif_id_dict}
    #return motif_id_dict.update(motif_alt_id_dict)
    return {**motif_id_dict, **motif_alt_id_dict}

def fimo_shuffling(region_coords, region_sequences, region_set, jobname, ref_fasta, chrom_sizes, ref_bed, excl_regions_bed, *args, **kwargs):
    '''Run FIMO to find known motifs in sets of OCR sequences'''
    #load required data
    motif_file = data["fimo"]
    
    jobname = "-".join(jobname.split("\t") + [''.join(random.choice(string.ascii_lowercase) for i in range(10))])   ## switch tabs for dashes, add random 10character string.
    
    outdir='/tmp/{}'.format(jobname)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    fimo_outdir=os.path.join(outdir,"fimo")
    fasta_file = os.path.join(outdir, "{}.fasta".format(jobname))

    if(region_set=='background'): #goes through scrambling
        #get non-redundant regions and convert to bedtools object
        regions_df=region_coords.to_dataframe(names=['chrom', 'start', 'end', 'query_gene'])
        #add the ocr id as the name column. this will correspond to the original ocr_id
        #only keep unique ocrs for shuffling
        regions_df['name']=regions_df.apply(lambda row: "{}:{}-{}".format(row['chrom'], row['start'], row['end']), axis=1)
        regions_to_shuffle=regions_df.filter(items=['chrom','start','end','name']).drop_duplicates()
        regions_bedtool=pybedtools.BedTool.from_dataframe(regions_to_shuffle)

        #randomly reposition regions
        shuffled=regions_bedtool.shuffle(g=chrom_sizes, noOverlapping=True,maxTries=5000, excl=excl_regions_bed).to_dataframe()
        #merge with the original regions file
        #since we kept the 'original' ocr/region_id, we can make sure that a repositioned region is
        #duplicated if it is linked to more than one gene
        shuffled_coords=regions_df.filter(items=['name', 'score', 'strand', 'query_gene']).merge(shuffled,on='name',how='left')
        #ocr id now corresponds to the SHUFFLED ocr.
        shuffled_coords['ocr_id'] = shuffled_coords.apply(lambda row: "{}:{}-{}".format(row['chrom'], row['start'], row['end']), axis=1)
        #get format for fimo
        #use gene_id as the name to keep it in the fasta file
        shuffled_coords=shuffled_coords.filter(items=['chrom', 'start', 'end', 'query_gene', 'score', 'strand', 'ocr_id'])

        ########Write fasta file for fimo
        #get sequences for shuffled coordinates
        shuffled_coords_bedtool=pybedtools.BedTool.from_dataframe(shuffled_coords)
        region_sequences=shuffled_coords_bedtool.sequence(fi=ref_fasta, name=True)
        region_sequences=open(region_sequences.seqfn).read().replace("::", " ")

    elif(region_set=="target"): # will not scramble sequences, just does fimo on the sequences as is.
        region_sequences = region_sequences.replace("::", " ")
        
    #Write out the sequences to the a fasta file. Keep redundant sequences
    fasta_file = os.path.join(outdir, "{}.fasta".format(jobname))
    fasta_num = 0
    out = open(fasta_file, "w")
    for seq in region_sequences.split("\n"):
        if seq.startswith(">"):
            out.write(">{}\n".format(fasta_num)) # doesn't it already have > in the header?
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
    #filtering for TFBS predictions with a score greater than 11.16
    fimo_results=(fimo_results[~fimo_results.motif_id.str.startswith('#')].query('`score`>11.16'))
    
    #for motif_ids, calculate the number of occurences
    fimo_motif_counts=(fimo_results.groupby(['motif_id'])['motif_id'].count().reset_index(name='n_occurrence'))
    #for motif_alt_ids, calculate the number of occurences
    fimo_motif_alt_counts=(fimo_results.groupby(['motif_alt_id'])['motif_alt_id'].count().reset_index(name='n_occurrence'))


    #Create a dictionary of each motif from the motif file
    #Add the n_occurence for each motif
    #If the motif is not found in the sequences, add the motif to the dictionary with a n_occurrence of 0
    motif_id_dict=dict()
    for mo_id in motif_ids:
        if mo_id in fimo_motif_counts['motif_id'].unique():
            fimo_df=fimo_motif_counts.loc[fimo_motif_counts['motif_id']==mo_id]
            motif_id_dict[mo_id]=fimo_df['n_occurrence'].values[0]
        else:
            motif_id_dict[mo_id]=0

    motif_alt_id_dict=dict()
    for mo_alt_id in motif_alt_ids:
        if mo_alt_id in fimo_motif_alt_counts['motif_alt_id'].unique():
            fimo_df=fimo_motif_alt_counts.loc[fimo_motif_alt_counts['motif_alt_id']==mo_alt_id]
            motif_alt_id_dict[mo_alt_id]=fimo_df['n_occurrence'].values[0]
        else:
            motif_alt_id_dict[mo_alt_id]=0

    
    #remove the temp files
    shutil.rmtree(outdir)
    return {**motif_id_dict, **motif_alt_id_dict}


def fimo_geneCount(region_sequences, jobname,*args, **kwargs):
    '''Run FIMO to find known motifs in sets of OCR sequences, and count the occurences by the number of genes'''

    #################################################################################
    ##                        Set up fasta file and directories                    ##
    #################################################################################  
    jobname = "-".join(jobname.split("\t") + [''.join(random.choice(string.ascii_lowercase) for i in range(10))])   ## switch tabs for dashes, add random 10character string.
    #fimo_outdir='/projects/nknoetze_prj/nknoetze_scratch/ocr_prj/results/{}'.format(jobname)
    outdir='/tmp/{}'.format(jobname)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    fimo_outdir=os.path.join(outdir,"fimo")
    

    motif_file = data["fimo"]

    ## replace :: with space in fasta
    ## currently, headers are >query_gene::chrX:start-end.
    region_sequences = region_sequences.replace("::", " ")
    
    #Write out the sequences to the a fasta file. Keep redundant sequences
    fasta_file = os.path.join(outdir, "{}.fasta".format(jobname))
    fasta_num = 0
    out = open(fasta_file, "w")
    for seq in region_sequences.split("\n"):
        if seq.startswith(">"):
            out.write(">{}.{}\n".format(seq,fasta_num)) #doesn't it already have > in the header?
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
    #filtering for TFBS predictions with a score greater than 11.16
    fimo_results=(fimo_results[~fimo_results.motif_id.str.startswith('#')].rename(columns={'sequence_name':'gene_id'}).query('`score`>11.16').filter(items=['motif_id','motif_alt_id','gene_id']).drop_duplicates())
    
    #for motif_ids, calculate the number of genes that have the motif
    fimo_motif_gene_counts=fimo_results.groupby('motif_id', as_index=False).agg({'gene_id': 'nunique'}).rename(columns={'gene_id':'n_gene'})

    #for motif_alt_ids, calculate the number of genes that have the motif
    fimo_motif_alt_gene_counts=fimo_results.groupby('motif_alt_id', as_index=False).agg({'gene_id': 'nunique'}).rename(columns={'gene_id':'n_gene'})

    #Create a dictionary of each motif from the motif file
    #Add the n_occurence for each motif
    #If the motif is not found in the sequences, add the motif to the dictionary with a n_occurrence of 0
    motif_id_gene_dict=dict()
    for mo_id in motif_ids:
        if mo_id in fimo_motif_gene_counts['motif_id'].unique():
            fimo_df=fimo_motif_gene_counts.loc[fimo_motif_gene_counts['motif_id']==mo_id]
            motif_id_gene_dict[mo_id]=fimo_df['n_gene'].values[0]
        else:
            motif_id_gene_dict[mo_id]=0

    motif_alt_id_gene_dict=dict()
    for mo_alt_id in motif_alt_ids:
        if mo_alt_id in fimo_motif_alt_gene_counts['motif_alt_id'].unique():
            fimo_df=fimo_motif_alt_gene_counts.loc[fimo_motif_alt_gene_counts['motif_alt_id']==mo_alt_id]
            motif_alt_id_gene_dict[mo_alt_id]=fimo_df['n_gene'].values[0]
        else:
            motif_alt_id_gene_dict[mo_alt_id]=0
    #remove the temp files
    shutil.rmtree(outdir)

    return {**motif_id_gene_dict, **motif_alt_id_gene_dict}

def fimo_geneCount_shuffling(region_coords,region_sequences,region_set,jobname,ref_fasta,chrom_sizes,ref_bed,excl_regions_bed,*args, **kwargs):
    '''Run FIMO to find known motifs in sets of OCR/region sequences that are randomly repositioned across the genome, and count the occurences by the number of genes'''

    #load required data
    motif_file = data["fimo"]

    jobname = "-".join(jobname.split("\t") + [''.join(random.choice(string.ascii_lowercase) for i in range(10))])   ## switch tabs for dashes, add random 10character string.

    #setup directories
    outdir='/tmp/{}'.format(jobname)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    fimo_outdir=os.path.join(outdir,"fimo")
    fasta_file = os.path.join(outdir, "{}.fasta".format(jobname))

    if(region_set=='background'): #goes through scrambling steps
        #################################################################################
        ##                                 Shuffle Sequences                           ##
        #################################################################################  
        #get non-redundant regions and convert to bedtools object
        regions_df=region_coords.to_dataframe(names=['chrom', 'start', 'end', 'query_gene'])
        #add the ocr id as the name column. this will correspond to the original ocr_id
        #only keep unique ocrs for shuffling
        regions_df['name']=regions_df.apply(lambda row: "{}:{}-{}".format(row['chrom'], row['start'], row['end']), axis=1)
        regions_to_shuffle=regions_df.filter(items=['chrom','start','end','name']).drop_duplicates()
        regions_bedtool=pybedtools.BedTool.from_dataframe(regions_to_shuffle)

        #randomly reposition regions
        #shuffled=regions_bedtool.shuffle(g=chrom_sizes, noOverlapping=True,incl=ref_bed,maxTries=5000).to_dataframe()
        shuffled=regions_bedtool.shuffle(g=chrom_sizes, noOverlapping=True,maxTries=5000, excl=excl_regions_bed).to_dataframe()
        #merge with the original regions file
        #since we kept the 'original' ocr/region_id, we can make sure that a repositioned region is
        #duplicated if it is linked to more than one gene
        shuffled_coords=regions_df.filter(items=['name', 'score', 'strand', 'query_gene']).merge(shuffled,on='name',how='left')
        #ocr id now corresponds to the SHUFFLED ocr.
        shuffled_coords['ocr_id'] = shuffled_coords.apply(lambda row: "{}:{}-{}".format(row['chrom'], row['start'], row['end']), axis=1)
        #get format for fimo
        #use gene_id as the name to keep it in the fasta file
        shuffled_coords=shuffled_coords.filter(items=['chrom', 'start', 'end', 'query_gene', 'score', 'strand', 'ocr_id'])
        #################################################################################
        ##                        Set up fasta file                                    ##
        #################################################################################  
        #get sequences for shuffled coordinates
        shuffled_coords_bedtool=pybedtools.BedTool.from_dataframe(shuffled_coords)
        region_sequences=shuffled_coords_bedtool.sequence(fi=ref_fasta, name=True)
        region_sequences=open(region_sequences.seqfn).read().replace("::", " ")

    elif(region_set=="target"): # will not scramble sequences, just does fimo on the sequences as is.
        region_sequences = region_sequences.replace("::", " ")

    #Write out the sequences to the a fasta file. Keep redundant sequences
    fasta_num = 0
    out = open(fasta_file, "w")
    for seq in region_sequences.split("\n"):
        if seq.startswith(">"):
            out.write(">{}.{}\n".format(seq,fasta_num)) # doesn't it already have > in the header?
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
    #filtering for TFBS predictions with a score greater than 11.16
    fimo_results=(fimo_results[~fimo_results.motif_id.str.startswith('#')].rename(columns={'sequence_name':'gene_id'}).query('`score`>11.16').filter(items=['motif_id','motif_alt_id','gene_id']).drop_duplicates())
        
    #for motif_ids, calculate the number of genes that have the motif
    fimo_motif_gene_counts=fimo_results.groupby('motif_id', as_index=False).agg({'gene_id': 'nunique'}).rename(columns={'gene_id':'n_gene'})

    #for motif_alt_ids, calculate the number of genes that have the motif
    fimo_motif_alt_gene_counts=fimo_results.groupby('motif_alt_id', as_index=False).agg({'gene_id': 'nunique'}).rename(columns={'gene_id':'n_gene'})

    #Create a dictionary of each motif from the motif file
    #Add the n_occurence for each motif
    #If the motif is not found in the sequences, add the motif to the dictionary with a n_occurrence of 0
    motif_id_gene_dict=dict()
    for mo_id in motif_ids:
        if mo_id in fimo_motif_gene_counts['motif_id'].unique():
            fimo_df=fimo_motif_gene_counts.loc[fimo_motif_gene_counts['motif_id']==mo_id]
            motif_id_gene_dict[mo_id]=fimo_df['n_gene'].values[0]
        else:
            motif_id_gene_dict[mo_id]=0

    motif_alt_id_gene_dict=dict()
    for mo_alt_id in motif_alt_ids:
        if mo_alt_id in fimo_motif_alt_gene_counts['motif_alt_id'].unique():
            fimo_df=fimo_motif_alt_gene_counts.loc[fimo_motif_alt_gene_counts['motif_alt_id']==mo_alt_id]
            motif_alt_id_gene_dict[mo_alt_id]=fimo_df['n_gene'].values[0]
        else:
            motif_alt_id_gene_dict[mo_alt_id]=0
    #remove the temp files
    shutil.rmtree(outdir)

    return {**motif_id_gene_dict, **motif_alt_id_gene_dict} 

def dinuc_shuffle(seq, num_shufs=None):
        # CODE ADAPTED FROM https://github.com/kundajelab/deeplift/blob/master/deeplift/dinuc_shuffle.py
        arr = np.frombuffer(bytearray(seq, "utf8"), dtype=np.int8)

        # Get the set of all characters, and a mapping of which positions have which
        # characters; use `tokens`, which are integer representations of the
        # original characters
        chars, tokens = np.unique(arr, return_inverse=True)

        # For each token, get a list of indices of all the tokens that come after it
        shuf_next_inds = []
        for t in range(len(chars)):
            mask = tokens[:-1] == t  # Excluding last char
            inds = np.where(mask)[0]
            shuf_next_inds.append(inds + 1)  # Add 1 for next token
    
        all_results = []

        for i in range(num_shufs if num_shufs else 1):
            # Shuffle the next indices
            for t in range(len(chars)):
                inds = np.arange(len(shuf_next_inds[t]))
                inds[:-1] = np.random.default_rng().permutation(len(inds) - 1)  # Keep last index same
                shuf_next_inds[t] = shuf_next_inds[t][inds]

            counters = [0] * len(chars)
        
            # Build the resulting array
            ind = 0
            result = np.empty_like(tokens)
            result[0] = tokens[ind]
            for j in range(1, len(tokens)):
                t = tokens[ind]
                ind = shuf_next_inds[t][counters[t]]
                counters[t] += 1
                result[j] = tokens[ind]
            all_results.append(chars[result].tostring().decode("ascii"))

        return all_results if num_shufs else all_results[0]

def fimo_geneCount_scrambling(region_sequences,region_set,jobname,ref_fasta,chrom_sizes,ref_bed,scramble_method,*args, **kwargs):
    '''Run FIMO to find known motifs in sets of OCR/region sequences that are randomly repositioned across the genome, and count the occurences by the number of genes'''
    #load required data
    motif_file = data["fimo"]
   
    jobname = "-".join(jobname.split("\t") + [''.join(random.choice(string.ascii_lowercase) for i in range(10))])   ## switch tabs for dashes, add random 10character string.

    #setup directories
    outdir='/tmp/{}'.format(jobname)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    fimo_outdir=os.path.join(outdir,"fimo")
    fasta_file = os.path.join(outdir, "{}.fasta".format(jobname))

    if(region_set=='background'): #goes through scrambling steps
        #################################################################################
        ##                            Scramble sequences                               ##
        #################################################################################
        gene_ocr_df=defaultdict(list)
        ocr_seq_df=defaultdict(str)
        seen=set()
        if scramble_method=="dinucleotide":    
            for line in region_sequences.split("\n"):
                if line.startswith(">"):## currently, headers are >query_gene::chrX:start-end.
                    gene_id=line.split("::")[0][1:]
                    ocr_id=line.split("::")[1]
                    gene_ocr_df[ocr_id].append(gene_id) 
                else:
                    seq=line
                    if ocr_id not in seen:#check if seq has been shuffled before
                        shuffled_seq=dinuc_shuffle(seq)
                        ocr_seq_df[ocr_id]=shuffled_seq
                        seen.add(ocr_id)
        elif scramble_method=="mononucleotide":
            for line in region_sequences.split("\n"):
                if line.startswith(">"):## currently, headers are >query_gene::chrX:start-end.
                    gene_id=line.split("::")[0][1:]
                    ocr_id=line.split("::")[1]
                    gene_ocr_df[ocr_id].append(gene_id) 
                else:
                    seq=line
                    if ocr_id not in seen:#check if seq has been shuffled before
                        seq_var=list(seq)
                        random.shuffle(seq_var)
                        shuffled_seq=''.join(seq_var)
                        ocr_seq_df[ocr_id]=shuffled_seq
                        seen.add(ocr_id)
        #Write out the sequences to the a fasta file. Keep redundant sequences
        
        fasta_num = 0
        out = open(fasta_file, "w")
        
        for ocr in ocr_seq_df: # get the scrambling sequence for each ocr id
            shuffled_seq=ocr_seq_df[ocr]
            for gene in gene_ocr_df[ocr]: #for each gene linked to the above OCR, write the fasta lines
                out.write(">{} {}.{}\n{}\n".format(gene, ocr,fasta_num,shuffled_seq)) # >gene_id ocr_id.num
                fasta_num+=1
        out.close()

    elif(region_set=="target"): # will not scramble sequences, just does fimo on the sequences as is.
        region_sequences = region_sequences.replace("::", " ") #currently, headers are >query_gene::chrX:start-end.
        
        #Write out the sequences to the a fasta file. Keep redundant sequences
        fasta_num = 0
        out = open(fasta_file, "w")
        for seq in region_sequences.split("\n"):
            if seq.startswith(">"):
                out.write(">{}.{}\n".format(seq,fasta_num))
                fasta_num += 1
            else:
                out.write("{}\n".format(seq))
        out.close()
    ################################################################################
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
    #filtering for TFBS predictions with a score greater than 11.16
    fimo_results=(fimo_results[~fimo_results.motif_id.str.startswith('#')].rename(columns={'sequence_name':'gene_id'}).query('`score`>11.16').filter(items=['motif_id','motif_alt_id','gene_id']).drop_duplicates())
        
    #for motif_ids, calculate the number of genes that have the motif
    fimo_motif_gene_counts=fimo_results.groupby('motif_id', as_index=False).agg({'gene_id': 'nunique'}).rename(columns={'gene_id':'n_gene'})

    #for motif_alt_ids, calculate the number of genes that have the motif
    fimo_motif_alt_gene_counts=fimo_results.groupby('motif_alt_id', as_index=False).agg({'gene_id': 'nunique'}).rename(columns={'gene_id':'n_gene'})

    #Create a dictionary of each motif from the motif file
    #Add the n_occurence for each motif
    #If the motif is not found in the sequences, add the motif to the dictionary with a n_occurrence of 0
    motif_id_gene_dict=dict()
    for mo_id in motif_ids:
        if mo_id in fimo_motif_gene_counts['motif_id'].unique():
            fimo_df=fimo_motif_gene_counts.loc[fimo_motif_gene_counts['motif_id']==mo_id]
            motif_id_gene_dict[mo_id]=fimo_df['n_gene'].values[0]
        else:
            motif_id_gene_dict[mo_id]=0

    motif_alt_id_gene_dict=dict()
    for mo_alt_id in motif_alt_ids:
        if mo_alt_id in fimo_motif_alt_gene_counts['motif_alt_id'].unique():
            fimo_df=fimo_motif_alt_gene_counts.loc[fimo_motif_alt_gene_counts['motif_alt_id']==mo_alt_id]
            motif_alt_id_gene_dict[mo_alt_id]=fimo_df['n_gene'].values[0]
        else:
            motif_alt_id_gene_dict[mo_alt_id]=0
    #remove the temp files
    shutil.rmtree(outdir)
    return {**motif_id_gene_dict, **motif_alt_id_gene_dict} 

def fimo_comotif_geneCount(region_sequences, jobname, num_motifs,*args, **kwargs):
    '''Run FIMO to find known motifs in sets of OCR sequences, and count the occurences by the number of genes'''

    #################################################################################
    ##                        Set up fasta file and directories                    ##
    #################################################################################  
    jobname = "-".join(jobname.split("\t") + [''.join(random.choice(string.ascii_lowercase) for i in range(10))])   ## switch tabs for dashes, add random 10character string.
    outdir='/tmp/{}'.format(jobname)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    fimo_outdir=os.path.join(outdir,"fimo")
    motif_file = data["fimo"]

    ## replace :: with space in fasta
    ## currently, headers are >query_gene::chrX:start-end.
    region_sequences = region_sequences.replace("::", " ")
    
    #Write out the sequences to the a fasta file. Keep redundant sequences
    fasta_file = os.path.join(outdir, "{}.fasta".format(jobname))
    fasta_num = 0
    out = open(fasta_file, "w")
    for seq in region_sequences.split("\n"):
        if seq.startswith(">"):
            out.write(">{}.{}\n".format(seq,fasta_num)) # doesn't it already have > in the header?
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
            raise Exception()

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
    fimo_results=(fimo_results[~fimo_results.motif_id.str.startswith('#')].rename(columns={'sequence_name':'gene_id'}).query('`score`>11.16').filter(items=['motif_id','motif_alt_id','gene_id']).drop_duplicates())
    #################################################################################
    ##                                 Co-Occuring Motifs                          ##
    #################################################################################
    comotif_alt_id_gene_dict=defaultdict(int)
    #Get all possible motif combinations
    for n in num_motifs:
        for gene in fimo_results['gene_id'].unique():
            gene_motifs=(motif[0] for motif in fimo_results.query('gene_id==@gene').filter(items=['motif_alt_id']).drop_duplicates().itertuples(index=False))
            gene_motif_combos=list(itertools.combinations(gene_motifs, n))

            for combo in gene_motif_combos:
                combo_sorted = tuple(sorted(combo))
                motif_str = ";".join(combo_sorted)
                comotif_alt_id_gene_dict[motif_str]+=1
    #remove the temp files
    shutil.rmtree(outdir)
    del fimo_results
    del gene_motif_combos
    del gene_motifs
    gc.collect()

    return comotif_alt_id_gene_dict

def fimo_comotif_geneCount_shuffling(region_coords, region_sequences,region_set, jobname,ref_fasta,chrom_sizes,ref_bed, num_motifs,excl_regions_bed,*args, **kwargs):
    '''Run FIMO to find known motifs in sets of OCR/region sequences, and count the occurences by the number of genes'''
    
    #load required data
    motif_file = data["fimo"]

    jobname = "-".join(jobname.split("\t") + [''.join(random.choice(string.ascii_lowercase) for i in range(10))])   ## switch tabs for dashes, add random 10character string.

    #setup directories
    outdir='/tmp/{}'.format(jobname)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    fimo_outdir=os.path.join(outdir,"fimo")
    fasta_file = os.path.join(outdir, "{}.fasta".format(jobname))


    if(region_set=='background'): # goes through shuffling steps
        #################################################################################
        ##                                 Shuffle Sequences                           ##
        #################################################################################  
        #get non-redundant regions and convert to bedtools object
        regions_df=region_coords.to_dataframe(names=['chrom', 'start', 'end', 'query_gene'])
        #add the ocr id as the name column. this will correspond to the original ocr_id
        #only keep unique ocrs for shuffling
        regions_df['name']=regions_df.apply(lambda row: "{}:{}-{}".format(row['chrom'], row['start'], row['end']), axis=1)
        regions_to_shuffle=regions_df.filter(items=['chrom','start','end','name']).drop_duplicates()
        regions_bedtool=pybedtools.BedTool.from_dataframe(regions_to_shuffle)

        #randomly reposition regions
        #shuffled=regions_bedtool.shuffle(g=chrom_sizes, noOverlapping=True,incl=ref_bed,maxTries=5000).to_dataframe()
        shuffled=regions_bedtool.shuffle(g=chrom_sizes, noOverlapping=True,maxTries=5000,excl=excl_regions_bed).to_dataframe()
        #merge with the original regions file
        #since we kept the 'original' ocr/region_id, we can make sure that a repositioned region is
        #duplicated if it is linked to more than one gene
        shuffled_coords=regions_df.filter(items=['name', 'score', 'strand', 'query_gene']).merge(shuffled,on='name',how='left')
        #ocr id now corresponds to the SHUFFLED ocr.
        shuffled_coords['ocr_id'] = shuffled_coords.apply(lambda row: "{}:{}-{}".format(row['chrom'], row['start'], row['end']), axis=1)
        #get format for fimo
        #use gene_id as the name to keep it in the fasta file
        shuffled_coords=shuffled_coords.filter(items=['chrom', 'start', 'end', 'query_gene', 'score', 'strand', 'ocr_id'])
        #################################################################################
        ##                        Set up fasta file                                    ##
        #################################################################################  
        #get sequences for shuffled coordinates
        shuffled_coords_bedtool=pybedtools.BedTool.from_dataframe(shuffled_coords)
        region_sequences=shuffled_coords_bedtool.sequence(fi=ref_fasta, name=True)
        region_sequences=open(region_sequences.seqfn).read().replace("::", " ")

    elif(region_set=='target'): # will not shuffle sequence. just does fimo on sequence as is.
        region_sequences = region_sequences.replace("::", " ")
    
    # #Write out the sequences to the a fasta file. Keep redundant sequences
    fasta_num = 0
    out = open(fasta_file, "w")
    for seq in region_sequences.split("\n"):
        if seq.startswith(">"):
            out.write(">{}.{}\n".format(seq,fasta_num)) # doesn't it already have > in the header?
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
            raise Exception()
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
    fimo_results=(fimo_results[~fimo_results.motif_id.str.startswith('#')].rename(columns={'sequence_name':'gene_id'}).query('`score`>11.16').filter(items=['motif_id','motif_alt_id','gene_id']).drop_duplicates())
    #################################################################################
    ##                                 Co-Occuring Motifs                          ##
    #################################################################################
    comotif_alt_id_gene_dict=defaultdict(int)
    #Get all possible motif combinations
    for n in num_motifs:
        for gene in fimo_results['gene_id'].unique():
            gene_motifs=(motif[0] for motif in fimo_results.query('gene_id==@gene').filter(items=['motif_alt_id']).drop_duplicates().itertuples(index=False))
            gene_motif_combos=list(itertools.combinations(gene_motifs, n))

            for combo in gene_motif_combos:
                combo_sorted = tuple(sorted(combo))
                motif_str = ";".join(combo_sorted)
                comotif_alt_id_gene_dict[motif_str]+=1
    #remove the temp files
    shutil.rmtree(outdir)
    del fimo_results
    del gene_motif_combos
    del gene_motifs
    gc.collect()

    return comotif_alt_id_gene_dict

def fimo_comotif_geneCount_scrambling(region_sequences,region_set, jobname,ref_fasta,chrom_sizes,ref_bed, scramble_method,num_motifs,*args, **kwargs):
    '''Run FIMO to find known motifs in sets of OCR/region sequences that are randomly repositioned across the genome, and count the occurences by the number of genes'''
    
    #load required data
    motif_file = data["fimo"]
    jobname = "-".join(jobname.split("\t") + [''.join(random.choice(string.ascii_lowercase) for i in range(10))])   ## switch tabs for dashes, add random 10character string.

    #setup directories
    outdir='/tmp/{}'.format(jobname)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    fimo_outdir=os.path.join(outdir,"fimo")
    fasta_file = os.path.join(outdir, "{}.fasta".format(jobname))

    if(region_set=='background'): # goes through scrambling steps
        #################################################################################
        ##                                 Scramble Sequences                           ##
        #################################################################################  
        gene_ocr_df=defaultdict(list)
        ocr_seq_df=defaultdict(str)
        seen=set()
        
        if scramble_method=="dinucleotide":    
            for line in region_sequences.split("\n"):
                if line.startswith(">"):## currently, headers are >query_gene::chrX:start-end.
                    gene_id=line.split("::")[0][1:]
                    ocr_id=line.split("::")[1]
                    gene_ocr_df[ocr_id].append(gene_id) 
                else:
                    seq=line
                    if ocr_id not in seen:#check if seq has been shuffled before
                        shuffled_seq=dinuc_shuffle(seq)
                        ocr_seq_df[ocr_id]=shuffled_seq
                        seen.add(ocr_id)
        elif scramble_method=="mononucleotide":
            for line in region_sequences.split("\n"):
                if line.startswith(">"):## currently, headers are >query_gene::chrX:start-end.
                    gene_id=line.split("::")[0][1:]
                    ocr_id=line.split("::")[1]
                    gene_ocr_df[ocr_id].append(gene_id) 
                else:
                    seq=line
                    if ocr_id not in seen:#check if seq has been shuffled before
                        seq_var=list(seq)
                        random.shuffle(seq_var)
                        shuffled_seq=''.join(seq_var)
                        ocr_seq_df[ocr_id]=shuffled_seq
                        seen.add(ocr_id)
        #Write out the sequences to the a fasta file. Keep redundant sequences
        
        fasta_num = 0
        out = open(fasta_file, "w")
        
        for ocr in ocr_seq_df:
            shuffled_seq=ocr_seq_df[ocr]
            for gene in gene_ocr_df[ocr]:
                out.write(">{} {}.{}\n{}\n".format(gene, ocr,fasta_num,shuffled_seq))# doesn't it already have > in the header?
                fasta_num+=1
        out.close()
    elif(region_set=="target"): # will not scramble sequences, just does fimo on the sequences as is.
        region_sequences = region_sequences.replace("::", " ")

        #Write out the sequences to the a fasta file. Keep redundant sequences
        fasta_num = 0
        out = open(fasta_file, "w")
        for seq in region_sequences.split("\n"):
            if seq.startswith(">"):
                out.write(">{}.{}\n".format(seq,fasta_num))
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
            raise Exception()
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
    fimo_results=(fimo_results[~fimo_results.motif_id.str.startswith('#')].rename(columns={'sequence_name':'gene_id'}).query('`score`>11.16').filter(items=['motif_id','motif_alt_id','gene_id']).drop_duplicates())
    #################################################################################
    ##                                 Co-Occuring Motifs                          ##
    #################################################################################
    comotif_alt_id_gene_dict=defaultdict(int)
    #Get all possible motif combinations
    for n in num_motifs:
        for gene in fimo_results['gene_id'].unique():
            gene_motifs=(motif[0] for motif in fimo_results.query('gene_id==@gene').filter(items=['motif_alt_id']).drop_duplicates().itertuples(index=False))
            gene_motif_combos=list(itertools.combinations(gene_motifs, n))

            for combo in gene_motif_combos:
                combo_sorted = tuple(sorted(combo))
                motif_str = ";".join(combo_sorted)
                comotif_alt_id_gene_dict[motif_str]+=1
    #remove the temp files
    shutil.rmtree(outdir)
    del fimo_results
    del gene_motif_combos
    del gene_motifs
    gc.collect()

    return comotif_alt_id_gene_dict