>01-T-cell_genes
    >>00-process_data
        Includes code to process raw RNA seq data from the DICE consortium and the SRA, and the associated metadata.
    >>01-genes
        Code to obtain median expression values for protein-coding genes for each cell type. It then uses these values to rank genes based on their expression in T cells and non-T cells. Remove genes with no OCRs.
>02-T-cell_OCRS
    >>00_process_data
        Code to tidy up the OCR annotations. Only keeps OCRs annotated in both CD4 and CD8+ T cells. Merges overlapping OCRs. Cleans up the interaction datasets.
>03-framework
    >>00-novel_motifs
        Code obtains all OCRs for the top-ranked T cells genes and breaks the sequences down into 12-mers. 12-mers are clustered using a hamming distance of 1. PPMs are created and any resembling a TFBS PPM from JASPAR are remove. PPMs are once again clustered into a set of putative novel motifs
    >>01-framework_run
        Snakemake pipeline to run the framework. Framework uses FIMO to identify sites in OCRs for TFBS and novel motifs. Counts how many genes have a regulatory element with the given motif. Compares between gene sets for significance. Uses a monte carlo based framework.
>04-ranked_candidates
    Code filters the framework results to motifs that have a pval < 0.05 and a FC >2. If its a TFBS motif, the TF must be expressed in T cells. Gets all sites for the enriched motifs and removes any sites that fully overlap. Ranks motifs based on their target gene count and their non-target gene count.
>05-library
    Code uses the ranked feautres to select 18 candidate motifs. First starts with pairs comprised of two TFBS. Gets the TFs in the top-ranked pair and looks at which archetype they belong to. Iterate down the list until 9 TFBS are chosen representing 9 unique archetypes. To select novel motifs, filters pairs to those containing one of the 9 TFBS and selects the top 9 motifs. To select a representative sequence, choose the most common sequence, if there's a tie, select the one with the best score from FIMO. Code also generates a synthetic sequence in silico and ensures there are no motifs. Creates the oligo library and the tiled positive control.
>06-STARR-seq
    >>00-process_reads
        Process the DNA and RNA seq reads. When necessary, extracts the 10bp UMI and appends to the read name. Trims reads and removes adapters. maps reads to the oligo library. Only keep primary alignments, properly paired reads etc. Then filter to only keep read pairs where both reads have a perfect match to their respective oligo. Counts how many reads there are for each UMI. Final counts are de-duplicated
    >>01-counts
        Code to create the matrix file containing de-deuplicated, CPM normalised, and log2fc for the oligos across all samples. Only keep results for oligos where DNA_cpm > 0 and an RNA_cpm >0 for at least two of the three RNA samples. 


