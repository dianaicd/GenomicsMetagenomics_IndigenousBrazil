#!/bin/bash

date=2018_07_19
wd=~/scratch_monthly/Botocudos/BAM/Bonobo

ref=/home/dcruzdva/archive/References/Bonobo/panPan2.fa
f=($(cat fastq.txt))
cd $wd
ind=($(cat ma_ids.txt))
mn=($(cat mn_ids.txt))

#ind=(MA2384)
#mn=(MN00021)
# 0 3
#ind=(MA2385)
#mn=(MN00022)
# 3

#ind=(MA2391)
#mn=(MN00019)
#problems with MN00019/5

#ind=(MA2392)
#mn=(MN00056)
#problems with MN00056/1

#ind=(MA2394)
#mn=(MN00066)
#problems with MN00066/7

#ind=(MA2395)
#mn=(MN00067)
#problems with MN00067/3
#ind=(MA2396)
#mn=(MN00068)
#problems with MN00068/0
#ind=(MA2397)
#mn=(MN00069)
#problems with MN00069/0
ind=(MA2399)
mn=(MN00118)
#problems with MN00118/6
#problems with MN00118/8
#problems with MN00118/10
#problems with MN00118/11
#problems with MN00118/12

for i in 0 
do
  f=($(grep ${ind[$i]} fastq.txt | sed 's/,/ /g' ))
#  mkdir ${mn[$i]}
  cd ${mn[$i]}

  echo ${ind[$i]}
  echo ${f[$j]}
  
  libs=$(echo "${#f[@]} - 1" | bc)
  for j in 6 8 10 11 12

  do
    echo $j
   # mkdir $j 
    cd $j
    mapping_aDNA.sh --fastq1 ${f[$j]} --base ${j}_Bonobo \
      --ref $ref \
      --noClean \
      -p 9 -q 1 > out_Bonobo_${j}_Q1 2> err_Bonobo_${j}_Q1 \ &

    cd .. 
  done
  cd $wd

done
