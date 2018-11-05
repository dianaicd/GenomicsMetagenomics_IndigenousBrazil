# MAKE PCA plot based on bam file and reference panel
# November 2012
# maricugh@gmail.com
# Update:
# This pipeline was optimized July 2013 by Maria to
# only filter out possible triallelic sites (most likely caused by sequence error)
# and ignore the possibility of damage
# ideally the bam files processed by this pipeline should have the base qualities
# recalibrated by mapDamage2 which reduces the quality of mismatches likely caused by damages
#
#  ****Modified Dec 2016 by Maria for Viridiana
#  added module loading, program names and few parsing commands
###
###
###Adapted to LIIGH cluster

### Arguments are
### $1=bam file
### $2=plink file basename
### $3=path to reference genome. IMPORTANT - this reference has to be in order, i.e. chr1 chr2 chr3 etc

file=$1
panel=$2
ref=$3

#samtools index $file
base=`basename $file .bam`

cut -f 1,4 $panel.map | sed 's/^/chr/' | sed 's/chr23/X/' | \
sed 's/chr24/Y/' > $panel.pos
##Convert the panel to homozygous
if [ ! -e $panel.rand.map ]
then
  echo "Converting panel to homozygous calls."
 awk '{printf $1" "$2" "$3" "$4" "$5" "$6; \
     for (i=7; i<=NF; i=i+2){ \
       x=int(.5 + rand()); \
       printf (" "$(i+x)" "$(i+x)); \
     } print "" \
   }' $panel.ped > $panel.rand.ped

 cp $panel.map $panel.rand.map
fi

#For bam file generates an mpileup file with the chr positions
#indicated in $2.pos and takes one allele randomly for each of the positions.
echo "Generating pileup of $base"
samtools mpileup -f $ref  -l $panel.pos $file | \
  get_random_allele_30.pl > $base.Q20.rand.mpileup


########################################Make tfam, tped and map files for sample############################################
############################################################################################################################

#tped file is generated with the chr and position coordinates followed
#by the homozygous allele.
#From the positions that overlap between the panel and the target.
echo "Generating tped of overlap between panel and sample."
overlap_mpileup_tped.pl  $panel.rand.map $base.Q20.rand.mpileup > $base.$panel.tped


###Prints a map file with only the coordinates of the overlapping
# positions in the mpileup file.
echo "Printing map file from pileup."
overlap_mpileup_map.pl $panel.rand.map $base.Q20.rand.mpileup  > $base.$panel.map

###Generates a fam file for the sample
echo "Generating fam file for the sample."
echo "$base $base 0 0 3 1" > $base.$panel.tfam


##Generates file tped.ids with the SNPs ids to be extracted from the reference
echo "Generating tped.ids with the SNPs ids to be extracted from the reference"
cut -f2 $base.$panel.tped > $base.$panel.tped.ids

###Extracts only the information for the SNPs indicated in tped.ids
echo "Extracting the information for the SNPs indicated in tped.ids"
plink --bfile $panel --extract $base.$panel.tped.ids --make-bed \
  --out $base.preFilt_sites --allow-no-sex


# SKIP TRIALLELIC SITES
###Generates a map file for the sites extracted using only
#the coordinates for homozygous
#paste $base.preFilt_sites.bim $base.$2.tped | awk '{if\
# ((($5 ~ /[A]/ || $6 ~ /[A]/) && $11 ~ /[A]/) || \
# (($5 ~ /[C]/ || $6 ~ /[C]/) && $11 ~ /C/) || \
# (($5 ~ /[G]/ || $6 ~ /[G]/) && $11 ~ /G/) || \
# (($5 ~ /[T]/ || $6 ~ /[T]/) && $11 ~ /[T]/)) \
# print $1,$2,$3,$4}' > $base.filt_sites.map

# SKIP TRIALLELIC SITES, C -> T and G -> A
echo "Skipping triallelic sites and C->T and G->A."
 paste $base.preFilt_sites.bim $base.$panel.tped | \
  awk '{if \
  (\
   (\
     (($5 ~ /[A]/ || $6 ~ /[A]/) && $11 ~ /[A]/) || \
     (($5 ~ /[C]/ || $6 ~ /[C]/) && $11 ~ /C/) || \
     (($5 ~ /[G]/ || $6 ~ /[G]/) && $11 ~ /G/) || \
     (($5 ~ /[T]/ || $6 ~ /[T]/) && $11 ~ /[T]/)\
   ) && \
   (\
    (($5 ~ /[C]/) && $11 !~ /T/) || \
    (($5 ~ /[G]/) && $11 !~ /A/)) \
   ) \
   print $1,$2,$3,$4}'  > $base.filt_sites.map

##Generates a file with the SNPs for the extraction of the overlapping
# positions in the reference panel and in the sample
echo "Generating $base.filt_sites.map.ids"
cut -f2 -d ' ' $base.filt_sites.map > $base.filt_sites.map.ids

###Extracts only the overlapping positions in the sample tfile
echo "Making $base.$panel.filt_sites"
plink --tfile $base.$panel --extract $base.filt_sites.map.ids \
  --recode --out $base.$panel.filt_sites


###Extracts only the overlapping positions in the reference file
echo "Extracting overlapping positions."
plink --file $panel.rand --extract $base.filt_sites.map.ids \
  --recode --out $base\_extract


##Merge the information of the sample and reference plink
#files that contain only the overlapping sites
echo "Making $base.filtSites\_MERGED"
plink --file $base.$panel.filt_sites --merge $base\_extract.ped $base\_extract.map \
  --recode --out $base.filtSites\_$panel\_MERGED --allow-no-sex


### BUILD PARFILES and run pca file

module add EcologyEvolution/eigensoft/6.1.4

echo "Creating files for smartpca."
mapfile=$base.filtSites\_$panel\_MERGED.map
mapbase=`basename $mapfile .map`
sed -i 's/ -9 / 1 /g' $mapbase.ped
###This depends on the name of your bam file.
samplename=`echo $base |cut -f1 -d '.'|cut -f1 -d'_'`
###This also depends on the name of your bam file and your sample.
sed -i 's/'$base'/'$samplename'/g' $mapbase.ped


###Generates a text file (.par) that indicates the file names where
#the information are contained:
echo "Creating -par file"
echo -e  genotypename: $mapbase.ped"\n"\
  snpname: $mapfile"\n"indivname: $mapbase.ped"\n"\
  evecoutname: $mapbase.evec"\n"\
  evaloutname: $mapbase.eval"\n"\
  numoutlieriter: 0 > $mapbase.par


##Prints the following line to keep track of the files processed
echo "smartpca -p $mapbase.par > $mapbase.out"

##Call smpartpca program, input file is .par generated in the command above.
#This command generates the files: .evec, .eval and .out.

smartpca -p $mapbase.par > $mapbase.out

###To make the plot of the PCA you will need the .evec and .eval files generated



exit;
