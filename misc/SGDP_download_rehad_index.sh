# Script to download files from the Simons Genome Diversity Project
#
dir_name=$1

is_incomplete=1
# 1. Try to download
while is_incomplete = 1
do
    while read line
    do
     wget -c $line &
    done < todownload.txt
    # Wait for downloads to finish
    { sleep 5; echo waking up after 5 seconds; } &
    { sleep 1; echo waking up after 1 second; } &
    wait
    echo "Files might or might not be fully downloaded"
    
    # 2. Compute MD5sum
    if [ -e validate.txt ]
    then
        rm validate.txt
    fi 
    touch validate.txt
    while read bam 
    do
        echo "validating $bam"
         md5sum $bam.bam >> validate.txt &
    done < ~/Project/Simons/${dir_name}_bam.txt
    # Wait for MD5sum to finish
    { sleep 5; echo waking up after 5 seconds; } &
    { sleep 1; echo waking up after 1 second; } &
    wait
    echo "Files might or might not be fully downloaded"
    rm todownload.txt
    is_incomplete=0

    # 3. Verify the MD5
    while read line
    do
        m=$(echo $line | cut -f1 -d ' ')
        file=$(echo $line | cut -f2 -d ' ')
        # Check if md5 value exists
        match=$(grep -c $m  ~/Project/Simons/prefix_md5.txt)
        if [ match = 0 ]
        then  # if 0 matches, then enter in the loop again
            is_incomplete=1
            grep $file ~/Project/Simons/bam_simons.txt | cut -f 16 >>todownload.txt
        fi 
    done < validate.txt
    rm validate.txt
done

# 4. Change the header and index
# To comper the SGDP bam files to those from the Americas,
# we need to look only at chromosomes 1 - 22.
# 
chroms=($(seq 1 22))

for bam in *bam
do
  ind=$(basename $bam .bam)
  name=$(grep $name ~/Project/Simons/Simons_sample_pop_region_country.txt | cut -f3)
  echo "$ind\t$name" >> bam_pop
  samtools view -b $ind.bam ${chroms[@]} > Reheaded/$name.bam &
done
{ sleep 5; echo waking up after 5 seconds; } &
{ sleep 1; echo waking up after 1 second; } &
  wait

cd Reheaded 
for bam in *bam 
do
    samtools index $bam &
done
{ sleep 5; echo waking up after 5 seconds; } &
{ sleep 1; echo waking up after 1 second; } &
  wait
  echo all jobs are done!