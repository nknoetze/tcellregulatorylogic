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
    repeatmasker_file = os.path.join(source_data_dir, "rmsk_filtered.txt")
    #ocr_group_specific_file = os.path.join(source_data_dir, "group_specific_ocrs.tsv")
    #ocr_single_gene_specific_file = os.path.join(source_data_dir, "gene_specific_ocrs.tsv")
    #jaspar_file = os.path.join(source_data_dir, "jaspar_filtered.txt")
    #kmer_file = os.path.join(source_data_dir, "kmers.txt")
    fimo_file = os.path.join(source_data_dir, "fimo_motifs.txt")
    ##not actually used...... may delete later lols.
    #filtered_fimo_file = os.path.join(source_data_dir, "fimo_motifs_expressed.txt")

    ## load data and populate features

    for func in functions_to_load:
        ## data
        if "fimo" in func and "fimo" not in data:
            data["fimo"] = fimo_file
        else:
            ## no special feature list, just set as single-length array with func name.
            feature_list[func] = [func]
    ## return feature_list (instead of relying on class features)
    return feature_list

def fimo(region_sequences, jobname, *args, **kwargs):
    '''Run FIMO to find known motifs in sets of OCR sequences'''
    #################################################################################
    ##                        Set up fasta file and directories                    ##
    #################################################################################  
    jobname = "-".join(jobname.split("\t") + [''.join(random.choice(string.ascii_lowercase) for i in range(10))])   ## switch tabs for dashes, add random 10character string.
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