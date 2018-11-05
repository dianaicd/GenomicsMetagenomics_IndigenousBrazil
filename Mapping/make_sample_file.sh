ids=~/archive/Botocudo/database_names.csv
readFolder=~/archive/FASTQ/Botocudo
samples_name=samples_24Botocudos_Quack.txt

echo -e "ID\tData\tMAPQ\tLB\tPL\tSM" > tmp_samples.txt
    paste <(ls $readFolder/*/*gz | cut -d'/' -f8 | cut -d'.' -f1) \
          <(ls $readFolder/*/*gz) \
          <(ls $readFolder/*/*gz | rev | cut -d'/' -f1 | rev | cut -d'_' -f-3) \
          <(ls $readFolder/*/*gz | rev | cut -d'/' -f1 | rev | cut -d'_' -f3 ) | \
    awk '{print $1,$2,"30",$3,"ILLUMINA",$4}' OFS='\t' >> tmp_samples.txt

while read line
do
  ma=$(echo $line | cut -f 1 -d','  )
  mn=$(echo $line | cut -f4 -d',' )
  echo "I will change $ma for $mn"
  sed -i "s/$ma$/$mn/" tmp_samples.txt
done < $ids

mv tmp_samples.txt ${samples_name}
