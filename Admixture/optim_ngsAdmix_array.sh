#!/bin/bash

#BSUB -L /bin/bash
#BSUB -J Wollstein[1-900]
#BSUB -o log/out.Wollstein-%I.txt 
#BSUB -e log/err.Wollstein-%I.txt
#BSUB -R "rusage[mem=2000]"
#BSUB -M 2000000

#### Your shell commands below this line ####

panel=88ind_nodamage
name=24ind
kr=($(cat ks.txt))

echo $panel

k=$(echo ${kr[$LSB_JOBINDEX]} |sed 's/_.*//')
rep=$(echo ${kr[$LSB_JOBINDEX]} |sed 's/.*_//')
if [ ! -d $k ]
then
  mkdir $k
fi

echo "REP: $rep ; K: $k"

    ~/install/NGSadmix -likes ${panel}_${name}.beagle \
    -K ${k} -o ${k}/${panel}_${name}_k${k}_${rep}
    echo "Done"
  
