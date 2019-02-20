# Script to download files from the Simons Genome Diversity Project
#
continent=$1
MD5_file=~/Project/Simons/prefix_md5.txt
project=~/Project/Simons
full_bam_list=~/Project/Simons/bam_simons.txt
# Must have the bamfile name and the region/pop/country in the
# third column. The field in the third column must also be unique
# E.g., Karitiana-2, Karitiana-3
names_pops=~/Project/Simons/Simons_sample_pop_region_country.txt

function download_file {
    # E.g., America
    continent=$1
    # E.g., LP6005441-DNA_A04.bam
    bamfile=$2
    is_incomplete=1

    cd $project/$continent
    # Try to download a file.

    url=$(grep $bam $full_bam_list | cut -f 16)
    while [ $is_incomplete = 1 ]
    do
        # First, try to continue the dowload
        fully_retrieved=0
        while [ $fully_retrieved = 0 ]
        do
            echo "Trying to download $bam."
            wget -vc $url 2>download_${bam}.txt 
            { sleep 5; echo waking up after 5 seconds; } &
            { sleep 1; echo waking up after 1 second; } &
            wait
            echo "File $bam might or might not be fully downloaded"
            fully_retrieved=$(grep -ci "fully retrieved" download_${bam}.txt)
        done
        
        #The file is apparently complete; 
        # try to do other more expensive (but useful) operations
        # that will tell you is the file is truncated/corrupted
        # We must try to download from zero if it is corrupted.

        # Quickcheck of bam
        echo "Quickcheck of $bam"
        samtools quickcheck -v $bam > bad_${bam}.fofn \
	        && error=0 \
	        || error=1
        if [ $error = 1 ] ; then rm $bam ; continue ; fi 

        # index
        echo "index $bam"
        samtools index $bam 2>index_${bam}.txt 
        error=$(wc -l index_${bam}.txt |cut -f1 -d' ')
        if [ $error > 0 ] ; then is_incomplete=1 ;rm $bam ; continue; fi
        
        # idxstats
        echo "idxstats $bam"
        samtools idxstats $bam 2>idxstats_${bam}.txt 
        error=$(wc -l idxstats_${bam}.txt |cut -f1 -d' ')
        if [ $error > 0 ] ; then is_incomplete=1 ; rm $bam ; continue; fi

        # md5
        echo "md5 $bam"
        md5sum $bam > md5sum_${bam}.txt  
        pattern=$(head -n1 md5sum_${bam}.txt | cut -f1 -d ' ')
        # Check if md5 value exists
        match=$(grep -c $pattern $MD5_FILE )
        if [ $match = 0 ]
        then 
            is_incomplete=1
            rm $bam 
            continue 
        else 
            is_incomplete = 0
            echo "Successfully downloaded" > complete_${bam}.txt 
        fi
    done
}


while read bam
do
    if [ ! -e complete_${bam}.txt ]
    then
        exists=$(grep -c $bam $full_bam_list)
        if [ $exists -gt 0 ]
        then
            download_file $continent "$bam.bam" &
        else
            echo "This file does not exist." > complete_${bam}.txt
        fi
    fi
done < $project/${continent}_bam.txt