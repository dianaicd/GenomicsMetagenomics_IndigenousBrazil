### Written by Viridiana Villa Islas (villa.islas.vi@gmail.com)
#Your job name
#$-N Schmutzi
#Use current working directory
#$-cwd

#Join stdout and stderr
#$-j y
# Run job through bash shell
#$-S /bin/bash

#Send an email after the job has finished
#$-m e
#$-M villa.islas.vi@gmail.com

#If modules are needed, source modules environment (Do not delete the next line):
 . /etc/profile.d/modules.sh

### Arguments are 
### $1=bam file
### $2=path to the mt reference

#A report is sent to the email below

# change email to receive results there 

EMAIL='villa.islas.vi@gmail.com'

PETH="/cm/shared/apps/schmutzi/29sep2017"

## $1 is the file that is taken as the first parameter of the program (e.g.in:  ./program.sh file.txt ; file.txt is saved in the $1 variable)
file=$1 
base=`basename $file .bam`
 
####Preparing data: Extracting fastq reads from sort.rmdup.uniq.bam
module load samtools/1.2


###MD Tagging
###MD tags can be used by programs to perform e.g. SNP detection without a reference genome

samtools calmd -b $1  $2 > $base_A.MD.bam

####Index MD Tagged bam file

samtools index $base_A.MD.bam

###

module load  r/3.4.1 
module load schmutzi/29sep2017


#######
####ContDeam

$PETH/contDeam.pl --library double --out $file --uselength --ref $2 $base_A.MD.bam

#/cm/shared/apps/schmutzi/29sep2017/contDeam.pl --library double --out MA1GH1SS002_S2_L001-4_R1-2_001.collapsed_bamto.fastq.bwa_mt.sort.rmdup.uniq.bam --uselength --ref /mnt/Cromosoma/mavila/vvilla/rCRS/rCRS.fasta MA1GH1SS002_S2_L001-4_R1-2_001.collapsed_bamto.fastq.bwa_mt.sort.rmdup.uniq.bam_A.MD.bam


####Schmutzi 

$PETH/schmutzi.pl --ref $2 --t 8 $file $PETH/alleleFreqMT/197/freqs/ $base_A.MD.bam

#/cm/shared/apps/schmutzi/29sep2017/schmutzi.pl --ref /mnt/Cromosoma/mavila/vvilla/rCRS/rCRS.fasta --t 8 MA1GH1SS002_S2_L001-4_R1-2_001.collapsed_bamto.fastq.bwa_mt.sort.rmdup.uniq.bam /cm/shared/apps/schmutzi/29sep2017/alleleFreqMT/197/freqs/ MA1GH1SS002_S2_L001-4_R1-2_001.collapsed_bamto.fastq.bwa_mt.sort.rmdup.uniq.bam_A.MD.bam



exit;










