# Script to download files from the Simons Genome Diversity Project
#
dir_name=$1

is_incomplete=1

# Maybe we ran this script the other day but we do not remember it...
if [ -e all_complete.txt ]
then
    is_incomplete=0
fi 

# 1. Try to download
while [ $is_incomplete = 1 ]
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
    touch validated.txt
    while read bam 
    do
        echo "validating $bam"
         md5sum $bam.bam >> tmp_md5.txt &
    done < ~/Project/Simons/${dir_name}_bam.txt
    # Wait for MD5sum to finish
    { sleep 5; echo waking up after 5 seconds; } &
    { sleep 1; echo waking up after 1 second; } &
    wait
    echo "Files might or might not be fully downloaded"
    rm todownload.txt
    is_incomplete=0
    cat validated.txt tmp_md5.txt > tovalidate.txt

    # 3. Verify the MD5
    while read line
    do
        m=$(echo $line | cut -f1 -d ' ')
        file=$(echo $line | cut -f2 -d ' ')
        # Check if md5 value exists
        match=$(grep -c $m  ~/Project/Simons/prefix_md5.txt)
        if [ $match = 0 ]
        then  # if 0 matches, then enter in the loop again
            is_incomplete=1
            grep $file ~/Project/Simons/bam_simons.txt | cut -f 16 >>todownload.txt
        else
            if [ ! -e $file.bai ]
            then
                echo "indexing $file"
                samtools index $file & 
                echo $line >> validated.txt
            fi
        fi 
    done < tovalidate.txt
    { sleep 5; echo waking up after 5 seconds; } &
    { sleep 1; echo waking up after 1 second; } &
    wait

    if [ $is_incomplete = 0 ]
    then
        echo "well done" > all_complete.txt
    fi 

    rm validate.txt
done

# 4. Change the header and index
# To compare the SGDP bam files to those from the Americas,
# we need to look only at chromosomes 1 - 22.
# 
chroms=($(seq 1 22))
ls *bam > bams_to_rehead.txt
notready=1

if [ ! -d Reheaded ]
then 
    mkdir Reheaded
fi

while [ $notready = 1 ] 
do
    while read bam
    do
        ind=$(basename $bam .bam)
        name=$(grep $ind ~/Project/Simons/Simons_sample_pop_region_country.txt | cut -f3)
        echo "$ind\t$name" >> bam_pop
        echo " $ind.bam to Reheaded/$name.bam"
        #samtools view -b $ind.bam ${chroms[@]} > Reheaded/$name.bam &
        rehead_bammds.sh $ind $name &
    done < bams_to_rehead.txt

    samtools quickcheck -v Reheaded/*.bam > bad_bams.fofn   && echo 'all ok' || echo 'some files failed check, see bad_bams.fofn'
    n=$(wc -l bad_bams.fofn | cut -f1 -d ' ')
    if [ $n -gt 0 ]
    then 
        mv bad_bams.fofn bams_to_rehead.txt
    else 
        notready=0
    fi 
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