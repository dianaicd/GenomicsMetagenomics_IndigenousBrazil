mkdir log
mkdir apparently_complete
maxCPU=20
todownload=download-links.txt
type=bam

awk -v type1=$type.* -v type2=$type '{sub(/.*\//, "", $1) ; sub(type1, type2, $1);print }' $todownload >names.txt
paste names.txt $todownload > name_link_${type}.txt


i=0
while read line
do
    usedCPU=$(echo "($i + 1) % $maxCPU" |bc)
    name=$(echo $line |cut -f1 -d' ')
    link=$(echo $line |cut -f2 -d' ')
    echo "Downloading $name"

    wget --no-use-server-timestamps -O $name $link >log/${name}.out 2> log/${name}.err & 

    if [ $usedCPU = 0 ]
    then
        echo "Waiting... $i"
        wait
        folder=$(echo $folder + 1 |bc)
        #mv *$type apparently_complete/
    fi
    i=$(echo $i + 1|bc )
done < name_link_${type}.txt