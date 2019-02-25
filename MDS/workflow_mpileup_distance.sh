bamfile=Botocudos_fromVictor
panel=Maanasa_americas_reheaded_filtered

ln -s ~/data/Git/Botocudos-scripts/AlleleCounts/count_and_sample.py ./
ln -s ~/data/Git/Botocudos-scripts/AlleleCounts/vcf_sample_random.py ./
ln -s ~/data/Git/Botocudos-scripts/MDS/counts2dist.py ./

CHR=($(seq 1 22))

# Get sites and reference, alternative alleles
cut -f 1,2,4,5 $panel.vcf |grep -v '#' |sed 's/\t/_/' >$panel.refalt

#-----------------------------------------------------------------------------#
# Run mpileup
for chr in ${chromosomes[@]}
do     
    cut -f 1,2 $panel.vcf |grep -v '#' |grep -P "^$chr\t"|\
        awk '{print($1"\t"$2-1"\t"$2)}' >${panel}_${chr}_sites.bed

    samtools mpileup -r $chr -Bl ${panel}_${chr}_sites.bed -b $bamfile \
     -a -o $bamfile_${chr}.mpileup -Q20 & 
done
  { sleep 5; echo waking up after 5 seconds; } &
  { sleep 1; echo waking up after 1 second; } &
  wait
  echo "Samtools mpileup done"

#-----------------------------------------------------------------------------#
# Count alleles and sample one base at random
for chr in ${CHR[@]}
do

    python count_and_sample.py $chr.mpileup \
        $chr.counts.gz $chr.sampled.gz $panel.refalt 

    gunzip -c ${chr}.sampled.gz >> $bamfile.counts.sampled.txt
done

#-----------------------------------------------------------------------------#
# Convert VCF to counts
python vcf_sample_random.py $panel.vcf $panel.counts $panel.counts.sampled.txt

paste $bamfile.counts.sampled.txt $panel.counts.sampled.txt > $bamfile.$panel.counts.sampled.txt

python counts2dist.py $bamfile.$panel.counts.sampled.txt $bamfile.$panel.dist