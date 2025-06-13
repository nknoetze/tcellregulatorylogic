
## based on ipynb from https://github.com/vierstralab/motif-clustering/tree/master on July 25, 2023


import sys
import argparse
import os
import time
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt


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


DEBUG = False
VERB = False




## Main

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = "title")
    parser.add_argument("--clustered_motifs", dest = "CLUSTERED_MOTIFS", help = "TOMTOM Clustered Motif file.", type = str)
    parser.add_argument("--pfm_dir", dest = "PFM_DIR", help = "Directory of single PFMs", type = str)
    parser.add_argument("--plot_dir", dest = "PLOT_DIR", help = "Directory of archetypic plots", type = str)
    parser.add_argument("--output", dest = "OUTPUT", help = "Uniprobe Motif Output", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB


    tomtom = pd.read_table(args.CLUSTERED_MOTIFS).rename(columns={'#Query_ID': 'Query_ID'})
    # Remove the last 3 rows
    tomtom = tomtom.iloc[:-3]

    #convert to int
    tomtom[['Optimal_offset', 'Overlap']] = tomtom[['Optimal_offset', 'Overlap']].astype(int)
    
    sim = tomtom.pivot_table(index='Query_ID', columns='Target_ID', values='E-value', fill_value=np.nan)
    cols = sim.columns
    rows = sim.index

    sim = sim[cols].loc[rows]
    x = sim.values

    w = np.triu(x) +  np.triu(x, 1).T
    v = np.tril(x) + np.tril(x, -1).T

    sim.iloc[:,:] = np.nanmin(np.dstack([w, v]), axis=2)


    sim.fillna(100, inplace=True)
    sim = -np.log10(sim)
    sim[np.isinf(sim)] = 10


    # Cluster the square matrix



    Z = linkage(sim, method = 'complete', metric = 'correlation')

    cl = fcluster(Z, 0.7, criterion='distance')

    print(f'Number of motif clusters: {max(cl)}')

    motif_annot_df = pd.DataFrame({'motif_id':sim.index, 'cluster':cl}).set_index('motif_id')
    motif_annot_df['cluster'] = 'AC' + motif_annot_df['cluster'].astype(str).str.zfill(4)
    motif_annot_df.head()


    def relative_info_content(pwm):
        p = pwm/np.sum(pwm, axis = 1)[:,np.newaxis]
        ic = 2+np.sum(p*np.nan_to_num(np.log2(p)), axis = 1)
        ric = p*ic[:,np.newaxis]
        return ric


    def rev_compl(st):
        nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join(nn[n] for n in reversed(st))



    def abs_mean(x):
        return np.mean(np.abs(x))

    def process_cluster(df, tomtom_df):
        
        import os.path
        
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as mgridspec
        
        from genome_tools.plotting import sequence
        

        motifs = df.index

        rows = (tomtom_df['Query_ID'].isin(motifs)) & (tomtom_df['Target_ID'].isin(motifs))
        all_pairwise_df = tomtom_df[rows]        
        
        # Seed motif has the best optimal_offset in the group
        
        seed_motif = all_pairwise_df.groupby('Query_ID').agg({'Overlap': abs_mean}).sort_values('Overlap', ascending=False).index[0]
        seed_motif = all_pairwise_df.groupby('Query_ID').agg({'E-value': np.median}).sort_values('E-value', ascending=True).index[0]

        
        rows = (tomtom_df['Query_ID'] == seed_motif) & (tomtom_df['Target_ID'].isin(motifs))
        pairwise_df = tomtom_df[rows] 
        
        pivot_df = all_pairwise_df.pivot_table(index='Query_ID', columns='Target_ID', values='Optimal_offset')
        q = pivot_df.loc[seed_motif]
        qi = q[~np.isnan(q)].index
        
        
        query = []
        target = []
        offset = []
        orientation = []
        target_consensus = []
        query_consensus = []
        
        for m in motifs.difference(pairwise_df['Target_ID']):
            t = pivot_df[m]
            ti = t[~np.isnan(t)].index
            
            try:
                common_motif = (qi&ti)[0]
            except:
                print(f'ERROR: {m} no alignments available!')
                continue
                
            row_q = (all_pairwise_df['Query_ID'] == seed_motif) & (all_pairwise_df['Target_ID'] == common_motif)
            row_t = (all_pairwise_df['Query_ID'] == common_motif) & (all_pairwise_df['Target_ID'] == m)
            
            #print(all_pairwise_df[row_t|row_q])
            
            offset_q = all_pairwise_df[row_q]['Optimal_offset'].iloc[0]
            offset_t = all_pairwise_df[row_t]['Optimal_offset'].iloc[0]
            
            orientation_q = all_pairwise_df[row_q]['Orientation'].iloc[0]
            orientation_t = all_pairwise_df[row_t]['Orientation'].iloc[0]
            
            consensus_q = all_pairwise_df[row_q]['Query_consensus'].iloc[0]
            
            consensus_tq = all_pairwise_df[row_t]['Query_consensus'].iloc[0]
            consensus_tt = all_pairwise_df[row_t]['Target_consensus'].iloc[0]
            
            offset_p = len(consensus_tt) - offset_t - len(consensus_tq)
            
            target.append(m)
            query_consensus.append(consensus_q)

            if orientation_t == orientation_q:
                orientation.append('+')
                target_consensus.append(consensus_tt)
                
                if orientation_t == '+':
                    offset.append(offset_q+offset_t)
                else:
                    offset.append(offset_p+offset_q)
            else:            
                orientation.append('-')
                target_consensus.append(rev_compl(consensus_tt))
                
                if orientation_q == '-':        
                    offset.append(offset_p+offset_q)
                else:
                    offset.append(offset_q+offset_t)
            
        z = pd.DataFrame({
            'Query_ID': seed_motif,
            'Target_ID': target,
            'Optimal_offset': offset,
            'p-value': 0, 
            'E-value': 0,
            'q-value': 1,
            'Overlap': 0,
            'Query_consensus': query_consensus,
            'Target_consensus': target_consensus,
            'Orientation': orientation,
        })
        
        
        if len(z) > 0:
            pairwise_df = pd.concat([pairwise_df, z])
            #print(z)
        
        #print(len(motifs), z.shape[0], pairwise_df.shape[0])
        
        w = pairwise_df['Target_consensus'].str.len()
        left = min(-pairwise_df['Optimal_offset'])
        l_offset = -left - pairwise_df['Optimal_offset']
        right = max(l_offset + w)
        r_offset = right - w - l_offset
        
        alignment_df = pairwise_df.drop(['Query_ID', 'Optimal_offset', 'p-value', 'E-value', 'q-value', 'Overlap', 'Query_consensus'], axis=1)
        alignment_df.loc[:,'w'] = w
        alignment_df.loc[:,'l_offset'] = l_offset
        alignment_df.loc[:,'r_offset'] = r_offset
        alignment_df.columns = ['motif', 'consensus', 'strand', 'w', 'l_offset', 'r_offset']
        
        alignment_df.reset_index(drop=True, inplace=True)

        alignment_df = alignment_df.merge(df.reset_index(), left_on='motif', right_on='motif_id')
        
        alignment_df.reset_index(inplace=True)
        
        n = len(alignment_df)
        l = min(alignment_df['l_offset'])
        r = max(alignment_df['r_offset'] + alignment_df['w'])
        w = r - l
    
        summed_pwm = np.zeros((4, w, n))
        
        for i, row in alignment_df.iterrows():
            
            motif_id = row['motif']
            rc = row['strand'] == '-'
            left = row['l_offset']
            width = row['w']

            motif_pfm = os.path.join(args.PFM_DIR, motif_id + '.pfm')
            pwm = np.loadtxt(motif_pfm)
        
            if rc:
                pwm = pwm[::-1,::-1]
            

            extended_pwm = np.ones((4, int(w))) * 0.25
            extended_pwm[:,left:left+width] = pwm
                
            summed_pwm[:,:,i] += extended_pwm
            
        avg_pwm = np.nanmean(summed_pwm, axis=2).T
    
        ic = relative_info_content(avg_pwm)
        total_ic = ic.sum(axis=1)

        cdf = np.cumsum(total_ic)/np.sum(total_ic)
        s = np.where(cdf > 0.05)[0][0]
        e = np.where(cdf > 0.95)[0][0] + 1    

        avg_pwm = avg_pwm[s:e,:]
        
        ## plot
        
        fig = plt.figure()
        fig.set_size_inches((w+2)*.125+2, (n+1)*0.5+1)
        
        gs = mgridspec.GridSpec(n+1, 1)
        
        for i, row in alignment_df.iterrows():
            ax = fig.add_subplot(gs[i+1, :])
            
            motif_id = row['motif']
            rc = row['strand'] == '-'
            left = row['l_offset']
            width = row['w']

            motif_pfm = os.path.join(args.PFM_DIR, motif_id + '.pfm')
            pwm = np.loadtxt(motif_pfm)
        
            if rc:
                pwm = pwm[::-1,::-1]
            
            sequence.seq_plot(relative_info_content(pwm.T), ax=ax, offset=left)

            ax.axvspan(l-1, s, fc='lightgrey', alpha=0.5)
            ax.axvspan(e, r+1, fc='lightgrey', alpha=0.5)
            
            ax.set_xlim(left=l-1, right=r+1)
            ax.set_ylim(bottom=0, top=2.1)

            ax.xaxis.set_visible(False)
            ax.set_yticks([])
            
            source_id = motif_id
            tf_name = ""

            ax.set_ylabel(tf_name + '\n (' + source_id + ')', rotation=0, ha='right', va='center', fontname="Courier", fontsize='medium')

        # Archetype motif
        ax = fig.add_subplot(gs[0,:])

        sequence.seq_plot(relative_info_content(avg_pwm), ax=ax, offset=s)

        ax.set_xlim(left=l-1, right=r+1)
        ax.set_ylim(bottom=0, top=2.1)
        ax.xaxis.set_visible(False)
        ax.set_yticks([])

        ax.axvspan(s, e, fc='none', ec='r', lw=2, clip_on=False)
        [ax.spines[loc].set_visible(False) for loc in ['top', 'bottom', 'left', 'right']]

        ax.set_ylabel('Archetype\nconsensus', rotation=0, ha='right', va='center', fontname="Courier", fontsize='large', fontweight='bold', color='r')
        
        cluster_id = str(alignment_df['cluster'][0])
        gene_family = "gene_family"
        dbd = "dbd"
        cluster_name = cluster_id #+ ':' + gene_family + ':' + dbd

        figw, figh = fig.get_size_inches()
        height_frac = (figh-0.75)/figh
        
    
        
        gs.update(left=1-((figw-1.75)/figw), right=(figw-0.25)/figw, top=(figh-0.75)/figh, bottom=1-((figh-0.25)/figh))
        
        fig.suptitle(cluster_name.upper(), fontname="IBM Plex Mono", fontweight='bold', fontsize='large', y=1-(.5/figh), va='center')
        plt.savefig(os.path.join(args.PLOT_DIR, f'{cluster_id}.pdf'))
        
        header_line =  cluster_name + '\n'
        mat = pd.DataFrame(avg_pwm.T, index=['A:', 'C:', 'G:', 'T:']).to_string(header=False)
        return header_line + mat


    plt.ioff()

    with open(args.OUTPUT, 'w') as fh:
        
        for cluster, df in motif_annot_df.groupby('cluster'):
            print(cluster)
            consensus_pwm = process_cluster(df, tomtom)
            fh.write(consensus_pwm + '\n\n')
