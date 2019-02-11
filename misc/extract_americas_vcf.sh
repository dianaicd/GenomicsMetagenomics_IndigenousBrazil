# Filter samples from VCF file

# July 18, 2018
# Diana Cruz

# Input file from ~/archive/Panels/fromVictor
whole=Maanasa_mask1_flip.americas

# I filtered the names of the groups by hand and saved them in:
name=Maanasa_mask1_flip_whole

# Get .clust file
for i in $(cat $name)
  do
  grep -e " ${i}$"  ${whole}.clust
done > ${name}.clust

# Create .names file
# All this work is to solve the problem with the annoying header
cut -f 1,2 -d ' ' ${name}.clust |sed 's/ /_/' > ${name}.names
bcftools view -S ${name}.names ${whole}.vcf -o ${name}.vcf -O v
awk '{print($3," ",$1)}' ${name}.clust |sed 's/_/-/g'| sed 's/ \+/_/' > new_names.txt
bcftools reheader -s new_names.txt ${name}.vcf > ${name}_reheaded.vcf

# Make the panel file I use to plot later
awk '{print($1,$1,$3,$3,$3)}' ${name}.clust |sed 's/ \+/\t/g' >${name}.panel.txt

# Saved names from Moreno Mayar et al in:
#div=Maanasa.americas.Moreno_names.txt
div=Maanasa.whole_names.txt
# Should merge it with previous Panel
out=Maanasa_whole.panel
Rscript ~/data/Scripts/merge_panel_divisions.R ${name}.panel.txt $div $out

plink --make-bed --mind 0.1 --out ${name}_reheaded_filtered --vcf ${name}_reheaded.vcf
plink --recode vcf --bfile ${name}_reheaded_filtered --out ${name}_reheaded_filtered

bcftools view -h ${name}_reheaded_filtered.vcf |sed "s/\t/\n/g" |grep -v "#" |  \
 tail -n +9  > ${name}_reheaded_filtered.names

sed 's/_/ /' ${name}_reheaded_filtered.names | \
awk '{print $1 "\t" $1 "\t" $2}' > ${name}_reheaded_filtered.clust

awk '{print($1,$3,$1,$3,$3)}' ${name}_reheaded_filtered.clust | \
  sed 's/ \+/\t/g' >${name}_reheaded_filtered.panel.txt








# With that, you can run an MDS
# Select the individuals
bcftools view -h 90ind_nodamage.vcf >header.txt
sed 's/\t/\n/g' header.txt |grep -v "#" >names_included.txt
seq 0 89 |sed 's/\t/\n/g' -  >index.txt
paste index.txt names_included.txt >index_selected.txt
  # Then manually cut out the rows I won't need

index=($(cut -f 1 index_selected.txt))
index_1=($(echo "${index[@]}" |sed 's/ /\*3+4\n/g' | sed '$s/$/\*3+4/' |bc ))
index_2=($(echo "${index[@]}" |sed 's/ /\*3+6\n/g' | sed '$s/$/\*3+4/' |bc ))
l=$((${#index_1[@]}-1))
s=($(for i in $(seq 0 $l); do echo ${index_1[$i]}-${index_2[$i]}; done))
s=$(echo ${s[@]} | sed 's/ /,/g' - )

cut -f 1-3,$s,$((90*3+2))-$((99*3+6)) 90ind_nodamage_10boto.beagle > \
Han_Pima_Karitiana_Surui_10boto.beagle

nind=$(wc -l index_selected.txt| sed 's/ .*//' -)
nind=$(($nind - 1 + 10))
strind=($(for i in $(seq 0 $nind); do echo Ind$i Ind$i Ind$i; done))
strind=(marker allele1 allele2 ${strind[@]})
echo ${strind[@]} |sed 's/ \+/\t/g' | sed 's/\t\+/\t/g' >header_vcf.txt
tail -n+2 Han_Pima_Karitiana_Surui_10boto.beagle >tmp.beagle
cat header_vcf.txt  tmp.beagle >Han_Pima_Karitiana_Surui_10boto.beagle
rm tmp.beagle
