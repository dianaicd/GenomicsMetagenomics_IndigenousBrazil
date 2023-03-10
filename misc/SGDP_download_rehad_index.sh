# Script to download files from the Simons Genome Diversity Project
#
dir_name=$1

incomplete_downloads=1
continue_download=1
fully_retrieved=0
md5_found=0

# Maybe we ran this script the other day but we do not remember it...
if [ -e all_complete.txt ]
then
    incomplete_downloads=0
fi 

# 1. Try to download
while [ $incomplete_downloads = 1 ]
do
    while read line
    do
        if [ $continue_download = 1 ]
        then
            wget -c $line &
        else
            wget -r $line &
            $continue_download = 1
        fi
    done < todownload.txt
    # Wait for downloads to finish
    { sleep 5; echo waking up after 5 seconds; } &
    { sleep 1; echo waking up after 1 second; } &
    wait
    echo "Files might or might not be fully downloaded"
    
    # 2. Compute MD5sum
    touch validated.txt
    while read line 
    do
        bam=$(basename $line .bam)
        echo "validating $bam"
         md5sum $bam.bam >> tmp_md5.txt &
    done < todownload.txt
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
    done < tmp_md5.txt
    { sleep 5; echo waking up after 5 seconds; } &
    { sleep 1; echo waking up after 1 second; } &
    wait

    if [ $incomplete_downloads = 0 ]
    then
        echo "well done" > all_complete.txt
    fi 

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
    { sleep 5; echo waking up after 5 seconds; } &
    { sleep 1; echo waking up after 1 second; } &
    wait
    
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