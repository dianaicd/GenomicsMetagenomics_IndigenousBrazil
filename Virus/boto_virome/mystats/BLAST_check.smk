# BLAST check. Is the first hit against the virus of interest?
# non_viral_reads_files.txt & taxids ({sample}.{virus}.taxid) files required

configfile: 'config_BLAST.yaml'
include: "parse_resources.smk"

def strip_jump(line):
    line = str(line)
    line = line.rstrip('\n')
    return(line)

with open(config['non_viral_reads'], 'r') as nvr_file:
    nvr_names = [strip_jump(line) for line in nvr_file.readlines()]

rule all:
    input:
        expand("BLAST/non_viral_reads/{non_viral}", non_viral=nvr_names)

rule get_fasta:
    input:
        "../{sample}/{sample}.{virus}.bam"
    output:
        "BLAST/fasta/{sample}.{virus}.fasta"
    log:
        "BLAST/logs/fasta/{sample}.{virus}.log"
    shell:
        '''
        module add UHTS/Analysis/samtools/1.10;
        samtools fasta {input} > {output} 2>{log}
        '''

rule run_BLAST:
    input:
        "BLAST/fasta/{sample}.{virus}.fasta"
    output:
        nt="BLAST/blast_nt/{sample}.{virus}.nt_reads",
        refseq="BLAST/blast_refseq/{sample}.{virus}.refseq_reads"
    log:
        "BLAST/logs/blast/{sample}.{virus}.log"
    resources:
        runtime=lambda wildcards, attempt: get_runtime_alloc("BLAST_time", attempt, 24)
    shell:
        '''
        module add Blast/ncbi-blast/2.10.1+
        blastn -db nt/nt -query {input} -out {output.nt} -outfmt \
        "6 qaccver saccver staxid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        2> {log}

        blastn -db refseq -query {input} -out {output.refseq} -outfmt \
        "6 qaccver saccver staxid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore" \
        2>> {log}
        '''

rule run_count_reads:
    input:
        ids="BLAST/virus_taxid/{sample}.{virus}.taxid",
        nt_algnmts="BLAST/blast_nt/{sample}.{virus}.nt_reads"
        # refseq_algnmts="BLAST/blast_refseq/{sample}.{virus}.refseq_reads"
    output:
        nt_hits="BLAST/nt_hits/{sample}.{virus}_hits.csv"
        # refseq_hits="BLAST/refseq_hits/{sample}.{virus}_hits.csv"
    log:
        "BLAST/logs/hits/{sample}.{virus}.log"
    shell:
        '''
        python3 select_reads.py {input.ids} {input.nt_algnmts} {output.nt_hits} 2> {log}
        '''
        # python3 select_reads.py {input.ids} {input.refseq_algnmts} {output.refseq_hits} 2>> {log}

rule get_non_viral_reads:
    input:
        nt_algnmts="BLAST/blast_nt/{sample}.{virus}.nt_reads",
        # refseq_algnmts="BLAST/blast_refseq/{sample}.{virus}.refseq_reads",
        nt_hits="BLAST/nt_hits/{sample}.{virus}_hits.csv"
        # refseq_hits="BLAST/refseq_hits/{sample}.{virus}_hits.csv"
    output:
        non_viral_reads="BLAST/non_viral_reads/{sample}.{virus}.non_viral_reads.tsv",
        uniq_nt=temp("BLAST/non_viral_reads/{sample}.{virus}.uniq_reads_nt"),
        read_id_hits_nt=temp("BLAST/non_viral_reads/{sample}.{virus}.read_id_hits_nt")
        # uniq_refseq=temp("BLAST/non_viral_reads/{sample}.{virus}.uniq_reads_refseq"),
        # read_id_hits_refseq=temp("BLAST/non_viral_reads/{sample}.{virus}.read_id_hits_refseq")
    log:
        "BLAST/logs/non_viral_reads/{sample}.{virus}.log"
    shell:
        '''
        # NT
        set +e
        cut -f1 {input.nt_algnmts} | sort | uniq > {output.uniq_nt} 2> {log}
        cut -f1 -d',' {input.nt_hits} > {output.read_id_hits_nt} 2>> {log}

        if [[ -s {output.read_id_hits_nt} && -s {output.uniq_nt} ]]; then
            grep -vf {output.read_id_hits_nt} {output.uniq_nt} > {output.non_viral_reads} 2>> {log}
            exitcode=$?
            if [[ exitcode -eq 1 ]]; then
                echo '##### All the alignments are against the virus of interest #####' > {output.non_viral_reads} \
                2>> {log};
                exit 0
            fi
        else
            echo '##### NT: none of the alignments against virus of interest OR no alignments given by BLAST #####' > \
            {output.non_viral_reads} 2>> {log}
        fi
        '''
        # RefSeq
        # cut -f1 {input.refseq_algnmts} | sort | uniq > {output.uniq_refseq} 2>> {log}
        # cut -f1 -d',' {input.refseq_hits} > {output.read_id_hits_refseq} 2>> {log}
        # echo '########## RefSeq ##########' >> {output.non_viral_reads} 2>> {log}

        # if [[ -s {output.read_id_hits_refseq} && -s {output.uniq_refseq} ]]; then
        #     grep -vf {output.read_id_hits_refseq} {output.uniq_refseq} >> {output.non_viral_reads} 2>> log
        # else
        #     echo '##### Refseq: all alignments against virus OR non alignments given by BLAST #####' >> \
        #     {output.non_viral_reads} 2>> {log}
        # fi