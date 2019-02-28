  # find the best replicate
  panel=$1 #88ind_nodamage
  ind=$2 #24ind
  
  n=$3 #10
  rep=$4 #10
  program=$5 #"ngsAdmix"
  pref=likelihoods

  if [ $program = "ADMIXTURE" ]
    then
      kr=($(cat ks.txt))
      for i in $(seq 1 $(echo $n*$rep - $rep | bc))
      do
        k=$(echo ${kr[$i]} |sed 's/_.*//')
        rep=$(echo ${kr[$i]} |sed 's/.*_//')
        grep ^Loglikelihood log/out.admixture-${i}.txt | \
          sed 's/Loglikelihood: /like=\t/' > \
          ${k}/${panel}_${ind}_k${k}_${rep}.log
      done
  fi

  for k in $(seq 2 $n)
  do
    for i in $(seq 1 $rep)
    do
      like=$(cat ${k}/${panel}_${ind}_k${k}_${i}.log | grep like | \
      sed 's/.*like=/\t/' |sed 's/ after.*//')
      echo "$i$like"
    done > ${pref}_k${k}_${panel}_${ind}.txt
#  mkdir k$k
#  mv ${panel}_${ind}_k${k}_* k$k
  done

Rscript ~/data/Scripts/find_likelihood.R $pref $panel ${ind}.txt $n $program
