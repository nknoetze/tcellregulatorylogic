### OUT DIR VARIABLES ###
outdir=config['OUTDIR']

### FILES ###
ranked_features=config['RANKED_FEATURES']
archetypes=config['MOTIF_ARCHETYPES']
### VARIABLES ###
conda_env=config['CONDA_ENV']
n_cores=config['N_CORES']
###
n_novel=config['N_NOVEL']
n_tfbs=config['N_TFBS']
###
num_candidates=config['NUM_CANDIDATES']
candidate_gap=config['CANDIDATE_GAP']
motifs=config['MOTIFS']
pad_re_5=config['PAD_RE_5']
re_pad_3=config['RE_PAD_3']
adaptor_5=config['ADAPTOR_5']
adaptor_3=config['ADAPTOR_3']
forbidden_seqs=config['FORIBIDDEN_SEQS']
output_file_name=config['OUTPUT_FILE_NAME']
control_seq=config['CONTROL_SEQ']
tile_length=int(config['TILE_LENGTH'])


# --------------------------------------------------------------
# SET RULE ALL
# --------------------------------------------------------------

rule all:
    input:
        expand([f"{{outdir}}/wetlab/bg-{{output_file_name}}.{nnovel}novel_{ntfbs}tfbs.oligolibrary.fa" for nnovel, ntfbs in zip(n_novel, n_tfbs)], outdir = outdir, output_file_name = output_file_name),
        expand(f"{outdir}/wetlab/tiledPosControl.seqs.fa",outdir=outdir),
        expand([f"{{outdir}}/wetlab/complete.{nnovel}novel_{ntfbs}tfbs.oligolibrary.fa" for nnovel, ntfbs in zip(n_novel, n_tfbs)], outdir = outdir),
        
# --------------------------------------------------------------
# SELECT 18 CANDIDATES. 9 TFBS AND 9 NOVEL MOTIFS
# --------------------------------------------------------------
rule select_candidates:
    input:
        ranked_features=ranked_features
    params:
        script="/projects/nknoetze_prj/ocr_prj/test_library/candidate.selection.R",
        archetypes=archetypes,
        outdir="{outdir}/wetlab"
    threads: 56
    output:
        table="{outdir}/wetlab/candidates.{nnovel}novel_{ntfbs}tfbs.tsv",
        fasta="{outdir}/wetlab/candidate_seqs.{nnovel}novel_{ntfbs}tfbs.fa"
    shell:
        "Rscript {params.script} --ranked_features {input.ranked_features} --motif_archetypes {params.archetypes} --n_novel {wildcards.nnovel} --n_tfbs {wildcards.ntfbs} --outdir {params.outdir}"

# --------------------------------------------------------------
# CREATE SYNTHETIC SEQUENCE
# --------------------------------------------------------------
rule generate_sequence:
    input:
        motifs=motifs
    params:
        script="/projects/nknoetze_prj/ocr_prj/test_library/motif_bashing_directedevolution.py",
        seqlen=125,
        seed=59,
        output='synthetic_seq.fasta'
    threads: 56
    output:
        synthetic_seq="{outdir}/synthetic_seq.fasta"
    shell:
        "python {params.script} --motifs {input.motifs} --seqlen {params.seqlen} --seed {params.seed} --output {params.output}"

# --------------------------------------------------------------
# GENERATE LIBRARY
# --------------------------------------------------------------
rule generate_library:
    input:
        fasta="{outdir}/wetlab/candidate_seqs.{nnovel}novel_{ntfbs}tfbs.fa",
        neutral_seq = "{outdir}/synthetic_seq.fasta"
    conda: conda_env
    params:
        script="/projects/nknoetze_prj/ocr_prj/test_library/oligo_designer.py",
       # neutral_seq = lambda wildcards: neutral_seqs[wildcards.output_file_name],
        pad_re_5=pad_re_5,
        re_pad_3=re_pad_3,
        adaptor_5=adaptor_5,
        adaptor_3=adaptor_3,
        num_candidates=num_candidates,
        candidate_gap=candidate_gap,
        forbidden_seqs= forbidden_seqs
    output:
        library = "{outdir}/wetlab/bg-{output_file_name}.{nnovel}novel_{ntfbs}tfbs.oligolibrary.fa"
    shell:
        "python {params.script} --candidates {input.fasta} --neutral_seq {input.neutral_seq} --pad_re_5 {params.pad_re_5} --re_pad_3 {params.re_pad_3} --adaptor_5 {params.adaptor_5} --adaptor_3 {params.adaptor_3} --num_candidates {params.num_candidates} --candidate_gap {params.candidate_gap} --forbidden_seqs {params.forbidden_seqs} --spacer --output {output.library}"
        # "python {params.script} --candidates {input.fasta} --neutral_seq {params.neutral_seq} --pad_re_5 {params.pad_re_5} --re_pad_3 {params.re_pad_3} --adaptor_5 {params.adaptor_5} --adaptor_3 {params.adaptor_3} --num_candidates {params.num_candidates} --candidate_gap {params.candidate_gap} --forbidden_seqs {params.forbidden_seqs} --spacer --output {output.library}"

# --------------------------------------------------------------
# CREATE TILED SEQUENCE CONTROL
# --------------------------------------------------------------
rule generate_tiled_control:
    input:
        control_seq=control_seq
    conda: conda_env
    params:
        script="/projects/nknoetze_prj/ocr_prj/test_library/tiled_control.py",
        pad_re_5=pad_re_5,
        re_pad_3=re_pad_3,
        adaptor_5=adaptor_5,
        adaptor_3=adaptor_3,
        forbidden_seqs=forbidden_seqs,
        tile_length=tile_length
    output:
        tiled_sequences = "{outdir}/wetlab/tiledPosControl.seqs.fa"
    shell:
        "python {params.script} --control_seq {input.control_seq} --pad_re_5 {params.pad_re_5} --re_pad_3 {params.re_pad_3} --adaptor_5 {params.adaptor_5} --adaptor_3 {params.adaptor_3} --forbidden_seqs {params.forbidden_seqs} --tile_length {params.tile_length} --output {output.tiled_sequences}"

# --------------------------------------------------------------
# Concatenate all files into one file for synthesis
# --------------------------------------------------------------
rule make_synthesis_file:
    input:
        libs = expand([f"{{outdir}}/wetlab/bg-{{output_file_name}}.{nnovel}novel_{ntfbs}tfbs.oligolibrary.fa" for nnovel, ntfbs in zip(n_novel, n_tfbs)], outdir = outdir, output_file_name = output_file_name),
        positive = expand(f"{outdir}/wetlab/tiledPosControl.seqs.fa",outdir=outdir)
    output:
        synthesize_file = "{outdir}/wetlab/complete.{nnovel}novel_{ntfbs}tfbs.oligolibrary.fa"
    shell:
        "cat {input.libs} {input.positive} > {output.synthesize_file}"
