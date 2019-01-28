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

# Determine the variant sites in the reference population
gunzip -c alleles_maf${minmaf}.mafs.gz | cut -f 1,2 |sed 1d > sites_${minmaf}.txt
angsd sites index sites_${minmaf}.txt 
angsd sites print sites_${minmaf}.txt

# "Ancient" individual to test
ind=MN0008
angsd -dogeno 4 -doPost 2 -postCutoff 0.95 \
  -doMajorMinor 1 -doMaf 2 -i ${ind}.bam \
  -out genos_${ind} -GL 1