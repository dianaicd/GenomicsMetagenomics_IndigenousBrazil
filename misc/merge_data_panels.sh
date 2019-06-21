# Merge multiple datasets
source ~/data/Git/Botocudos-scripts/misc/source4merging.sh
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

#=============================================================================#
#-----------------------------------------------------------------------------#
 # Panels

#***************************#
# Maanasa
source ~/data/Git/Botocudos-scripts/misc/source4merging.sh

 # Prepare list with bam files
 panel=Maanasa_mask1_flip.bed

file=B.P.ASM.VMAnc
 if [ -e $file ] ; then rm $file ; fi
for BAMS in BotocudosAll Posth Malaspinas2014 fromVictorAncient
do
 make_bamlist --samples $BAMS --output $BAMS 
 cat $BAMS >> $file
done

mergeBAM2BED --panel $panel --bamlist $file 


# Compute genotype likelihoods once you are happy with the merged panel
geno_like --panel $panel --homozygous no --bamlist $bamlist --rmdamage no

 #**************************#
 # Human Origins
#***************************#
 # 1.2M
#***************************#
# Wollstein
source ~/data/Git/Botocudos-scripts/misc/source4merging.sh

panel=Jorde_Wollstein_hg19_final_noseconddegree_geno01.bed
file=B.P.ASM.VMAnc
 if [ -e $file ] ; then rm $file ; fi
for BAMS in BotocudosAll Posth Malaspinas2014 fromVictorAncient
do
 make_bamlist --samples $BAMS --output $BAMS 
 cat $BAMS >> $file
done

mergeBAM2BED --panel $panel --bamlist $file 

# Convert to eigenstratgeno once you have a panel you are happy with
# first, modify IDs so they are <39 characters
clus=/home/dcruzdva/archive/Panels/fromVictor/Maanasa_mask1_flip.clust 

if [ -e tmp.fam ] ; then rm tmp.fam ; fi  
while read line
do
   ind=$(echo $line |cut -f1 -d ' ' )
   match=$(cat $clus |sed 's/\t/ /g' |grep -P "^$ind " | cut -f2,3 -d ' ' \
      |sed 's/.variant/.v/g ; s/Great_Andamanese_//')
      
   if [ "$match" == "" ] ; then
      echo "$ind $ind 0 0 0 1" >>tmp.fam
   else
      echo "$match 0 0 0 1" >> tmp.fam
   fi 
done < $panel.$file.fam 

mv tmp.fam ${panel}.${file}.fam
# Do not forget to recode in order to change IDs for ind and pop!!!
plink --recode --bfile ${panel}.${file} --out ${panel}.${file}
make_par --panel $panel.$file --conversion ped2eigenstrat
convertf -p ${panel}.${file}_ped2eigenstrat.par

# Compute genotype likelihoods once you are happy with the merged panel

geno_like --panel $panel --homozygous no --bamlist $bamlist --rmdamage yes

#*****************************************************************************#
2# h.bed
source ~/data/Git/Botocudos-scripts/misc/source4merging.sh

panel=h.bed

file=B.P.ASM.VMAnc
 if [ -e $file ] ; then rm $file ; fi
for BAMS in BotocudosAll Posth Malaspinas2014 fromVictorAncient
do
 make_bamlist --samples $BAMS --output $BAMS 
 cat $BAMS >> $file
done

mergeBAM2BED --panel $panel --bamlist $file 

# Convert to eigenstratgeno once you have a panel you are happy with
# first, modify IDs so they are <39 characters
clus=/home/dcruzdva/archive/Panels/fromVictor/Genotypes/h_clust

if [ -e tmp.fam ] ; then rm tmp.fam ; fi  
while read line
do
   ind=$(echo $line |cut -f1 -d ' ' )
   match=$(grep -P "^$ind\t" $clus |sed 's/\t/ /g' | cut -f2,3 -d ' ' |sed 's/.variant/.v/g')
   if [ "$match" == "" ] ; then
      echo "$ind $ind 0 0 0 1" >>tmp.fam
   else
      echo "$match 0 0 0 1" >> tmp.fam
   fi 
done < $panel.$file.fam 

mv tmp.fam ${panel}.${file}.fam
# Do not forget to recode in order to change IDs for ind and pop!!!
plink --recode --bfile ${panel}.${file} --out ${panel}.${file}

make_par --panel $panel.$file --conversion ped2eigenstrat

convertf -p ${panel}.${file}_ped2eigenstrat.par

# Compute genotype likelihoods once you are happy with the merged panel
geno_like --panel $panel --homozygous no --bamlist $bamlist --rmdamage yes

#make_bamlist --samples Lindo --output Lindo 
# YanaLow
# MalTa
# TianYuan
# Scheib
bamlist=BPMVa
cat Botocudos Posth MaanasaAncient fromVictorAncient > $bamlist

mergeBAM2BED --panel $panel --bamlist Botocudos

panel=h.Botocudos.bed 
 
mergeBAM2BED --panel $panel --bamlist Posth

#***************************#
 # SGDP


 # BAM files
    # Ancient/low coverage
 # Botocudos
 # Posth
 # fromVictor
 # Lindo
 # Malaspinas2014
 
    # Modern
 # SGDP
 # Team A, Team B
 # Lindo2018, Aymara30x
 # fromVictor, ancient and modern
 # Maanasa, newworld, oldworld and ancient


 # Other
    # Ancient
 # Tianyuan
 # Mal'ta
    # Modern
 # Reich2012
  # Skoglund2015