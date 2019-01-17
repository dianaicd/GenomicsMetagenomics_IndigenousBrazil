path=~/Project/Simons 
rm missing_simons.txt
touch missing_simons.txt

while read line
do

  id=$(echo $line | cut -f1 -d ' ')
  sample=$(echo $line | cut -f3 -d ' ')
  pop=$(echo $line | cut -f 4 -d ' ')
  region=$(echo $line| cut -f5 -d ' ')
  echo "ID is $id"
  if [ -e ${path}/${id}.bam ]
  then
    echo "$id.bam will be renamed to ${region}_${pop}_${sample}.bam"
    ln -s ${path}/${id}.bam ./${region}_${pop}_${sample}.bam
    samtools index ${region}_${pop}_${sample}.bam &
  else
    echo "$id.bam not in your list"
    echo $line >> missing_simons.txt
  fi
done < Simons_sample_pop_region_country.txt

{ sleep 5; echo waking up after 5 seconds; } &
{ sleep 1; echo waking up after 1 second; } &
  wait
  echo all jobs are done!

