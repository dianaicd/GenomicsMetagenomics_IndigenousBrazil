# Script to run D-stats with admixtools
# It includes personal communication with J. Victor Moreno-Mayar

panel=h_bamlist
## Plink .fam debe tener el nombre de la población en la primera columna
# Your plink file .fam must have the population name in the first column

 # 1. mientras todavía tienes tu plink en 
 # formato binario (bed bim fam), asegúrate de que 
 # la última columna de tu fam sean 1 (no 0 ni -9 ni otra cosa)
# 1. While you plink is in binary format, make sure the last column
# of the .fam family is 1 (not 9 nor any other value)
awk 'BEGIN{DOF="\t"} {print $1,$2,$3,$4,$5,1}' $panel.fam >tmp.fam 
mv tmp.fam $panel.fam 

 # 2. ya después lo conviertes a ped y aquí te va un ejemplo 
 # de mi archivo .par para el convertf
 # 2. Convert it to ped and specify .par file
# genotypename:    Panel.ped
# snpname:    Panel.map
# indivname:    Panel.ped
# outputformat:    EIGENSTRAT
# genotypeoutname:    Panel.eigenstratgeno
# snpoutname:    Panel.snp
# indivoutname:    Panel.ind
# familynames:    YES

if [ ! -e $panel.ped ]
then
    plink recode --bfile ${panel} --out $panel
fi

argument=(genotype snpname indivname outputformat genotypeoutname \
 snpoutname indivoutname familynames)
options=($panel.ped $panel.map $panel.ped EIGENSTRAT $panel.eigenstratgeno \
 $panel.snp $panel.ind YES)

for i in $(seq 0 7)
do 
  echo "${argument[$i]}: ${options[$i]}" >> ${panel}_convert.par
done

convertf -p ${panel}_convert.par 

# 3. después corres esto
R
a<-read.table("Panel.ind")
p<-as.character(a[,1])
p<-unlist(lapply(strsplit(p, ":"), "[[", 1))
a[,3]<-p
write.table(a, file="Panel.ind", sep=" ", quote=F, row.names=F, col.names=F)
q("no")

# cut -f 3 -d " " Panel.ind | sort | uniq > Panel.ind.list
pop=LagoaSanta
combFile=pop_combn_${pop}.txt
maxCPU=120
maxLines=10
total=$(wc -l $combFile |cut -f1 -d' ')
maxLoop=$(echo "$total / $maxLines" |bc)
parfile=maanasa_${pop}.par 

for i in $(seq 0 $maxLoop)
do 
    lo=$(echo "$i*$maxLines + 1" |bc)
    hi=$(echo "($i+1)*$maxLines" |bc)
    usedCPU=$(echo "$i % $maxCPU" |bc)

    qpDstat -p $parfile -l $lo -h $hi >Results/out_${lo}_${hi}_${pop}.txt 2>Results/err_${lo}_${hi}_${pop}.txt &
    #echo $lo $hi 
    if [ $usedCPU = 0 ]
    then
        echo "Waiting... $i"
        wait
    fi

done