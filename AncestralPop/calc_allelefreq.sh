# Calculate allele frequencies (reference population)

# use outgroup to infer the ancestral state
anc=Yoruba
# minimum allele frequency
minmaf=0.05
# number of processors
p=32

 angsd -out alleles_maf${minmaf} -doMaf 3 \
       -doMajorMinor 1 -anc $anc.fa \
       -bam bam.filelist -GL 1 -minMaf $minmaf -p $p
