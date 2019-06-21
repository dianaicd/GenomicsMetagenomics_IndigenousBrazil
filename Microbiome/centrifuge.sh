# database:
index=/scratch/cluster/monthly/sneuensc/sapfo/bin/centrifuge-1.0.4-beta/indices/nt
# bam
name=Bot15

# threads
nThreads=24

samtools fastq -@$nThreads $name/${name}_low_qual.bam |gzip >$name/$name.unmapped.fastq.gz

centrifuge -x $index -U $name/$name.unmapped.fastq.gz -S $name.centrifuge \
 --report-file ${name}_report.tsv -p $nThreads
