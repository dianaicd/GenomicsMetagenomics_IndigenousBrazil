while read ind
do
  echo "Working on $ind "
  name=$(basename $ind )
  samtools sort -o Sorted/${name} -@ 3 ${ind}.bam &
done < list_bam.txt
{ sleep 5; echo waking up after 5 seconds; } &
{ sleep 1; echo waking up after 1 second; } &
  wait
  echo all jobs are done!

