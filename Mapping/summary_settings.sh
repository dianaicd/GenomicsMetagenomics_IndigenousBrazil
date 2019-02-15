
# Script to work on .settings file
#sample=$1
library=$1
ids=( "$2" )
#name=$4

#pref=$(echo $sample$'\t'$library$'\t'$id)
seq_reads_se=0
seq_trash_se=0
seq_trash_se_frac=0
seq_retained_reads=0
seq_retained_nts=0
seq_retained_length=0

for id in ${ids[@]}
do
  seq_reads_se_id=$(grep 'Total number of reads:' $id.settings | sed 's/.*: //g')
  seq_reads_se=$(echo "$seq_reads_se + $seq_reads_se_id" |bc -l)

  seq_trash_se_id=$(grep 'Number of discarded mate 1 reads:' $id.settings | sed 's/.*: //g')
  seq_trash_se=$(echo "$seq_trash_se + $seq_trash_se_id" |bc -l)


  seq_retained_reads_id=$(grep 'Number of retained reads:' $id.settings | sed 's/.*: //g')
  seq_retained_reads=$(echo "$seq_retained_reads + $seq_retained_reads_id "|bc -l)

  seq_retained_nts_id=$(grep 'Number of retained nucleotides:' $id.settings | sed 's/.*: //g')
  seq_retained_nts=$(echo "$seq_retained_nts + $seq_retained_nts_id" |bc -l)

done

seq_trash_se_frac=$(echo "$seq_trash_se / $seq_reads_se" | bc -l)
seq_retained_length=$(echo "$seq_retained_nts / $seq_retained_reads" |bc -l)

# Fake the .settings file
echo "Fake $library.settings" >$library.settings
echo "Total number of reads: $seq_reads_se" >> $library.settings
echo "Number of discarded mate 1 reads: $seq_trash_se" >> $library.settings
echo "Number of retained reads: $seq_retained_reads" >> $library.settings
echo "Number of retained nucleotides: $seq_retained_nts" >> $library.settings

answer=$(echo $seq_reads_se$'\t'$seq_trash_se$'\t'$seq_trash_se_frac$'\t'$seq_retained_reads$'\t'$seq_retained_nts$'\t'$seq_retained_length)

echo ${answer}
