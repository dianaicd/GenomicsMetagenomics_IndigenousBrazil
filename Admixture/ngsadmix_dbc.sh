
enough_jobs(){
    local maxJobs=$1
    local i=$2
    launchedJobs=$(echo "($i + 1) % $maxJobs" |bc)
    if [ $launchedJobs = 0 ]
    then
        echo "Waiting... $i"
        wait
    fi
}

maxJobs=40
i=0
nThreads=4
if [ ! -e log ] ; then mkdir log ; fi
for k in $(seq 11 15)
do
    if [ ! -d $k ] ; then mkdir $k ; fi
    for rep in $(seq 1 100)
    do
        NGSadmix -likes ${panel}_${bamlist}.beagle \
            -K ${k} -P $nThreads \
            -o ${k}/${panel}_${bamlist}_k${k}_${rep} >\
            log/out.$panel.$bamlist.$rep.$k.txt 2> \
            log/err.$panel.$bamlist.$rep.$k.txt &

    enough_jobs $maxJobs $i
    i=$(echo $i +1 |bc)
    done
done
