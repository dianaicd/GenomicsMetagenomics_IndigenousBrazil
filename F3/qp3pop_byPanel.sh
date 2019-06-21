# Script to run F3-statistics, by panel
ln -s ~/data/Git/Botocudos-scripts/Dstat/switchPops_admixTools.r
#-----------------------------------------------------------------------------#
make_3poplist(){
    SHORTOPTS="p:b:"
    LONGOPTS="panel: bamlist: outgroup:"
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)     local panel=$2; shift;;
            -b|--bamlist)  local bamlist=$2; shift;;
            -o|--outgroup) local outgroup=$2; shift;;
        esac 
        shift
    done

    local pops=($(cut -f3 -d ' ' $panel.$bamlist.modified.ind |sort |uniq))

    if [ -e $panel.$bamlist.list_qp3test ]
    then rm $panel.$bamlist.list_qp3test
    fi

    nPops=$(echo ${#pops[@]} - 1 | bc)
    for i in $(seq 0 $nPops) ; do
        start=$(echo $i+1 |bc)
        for j in $(seq $start $nPops) ; do
            if [ ! "${pops[$i]}" = "$outgroup" ] ; then 
                if [ ! "${pops[$j]}" = "$outgroup" ] ; then 
                    echo "${pops[$i]} ${pops[$j]} $outgroup" >> $panel.$bamlist.list_qp3test
                fi
            fi 
        done
    done
}
#-----------------------------------------------------------------------------#
make_par(){
    SHORTOPTS="p:b:"
    LONGOPTS="panel: bamlist:"
    # This function uses my scripts
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)     local panel=$2; shift;;
             -b|--bamlist)  local bamlist=$2; shift;;
        esac 
        shift
    done

    local parops=(genotypename snpname indivname popfilename)
    local pargs=(eigenstratgeno snp modified.ind list_qp3test)

    if [ -e $panel.$bamlist.qp3test.par ]
    then rm $panel.$bamlist.qp3test.par 
    fi 

    for i in $(seq 0 3)
    do
        echo "${parops[$i]}: $panel.$bamlist.${pargs[$i]}" \
            >>$panel.$bamlist.qp3test.par
    done
}
#-----------------------------------------------------------------------------#
verify_length_IDs(){
    SHORTOPTS="p:b:"
    LONGOPTS="panel: bamlist:"
    # This function uses my scripts
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi

    eval set -- "$ARGS"
    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)     local panel=$2; shift;;
             -b|--bamlist)  local bamlist=$2; shift;;
        esac 
        shift
    done
    
    local i=0
    if [ -e $panel.$bamlist.longIDs ] ; then rm $panel.$bamlist.longIDs ; fi 
    while read line
    do
        local id=$(echo $line |cut -f1 -d ' ' )
        local lengthID=$(echo $id |wc -c)

        if [ $lengthID -gt 38 ] ; 
        then 
            echo "problem with $id of length $lengthID"
            echo "$i $id" >> $panel.$bamlist.longIDs
            i=$(echo $i+1 |bc)
            ind=$(echo $id |cut -f1 -d ':')
            sed -i "s/$ind:/$i:/" $panel.$bamlist.modified.ind 
        fi 
    done < $panel.$bamlist.modified.ind 

}
#-----------------------------------------------------------------------------#
# Maanasa
panel=Maanasa_mask1_flip
bamlist=B.P.ASM.VMAnc
outgroup=Yorubas
Rscript switchPops_admixtools.r ~/archive/Panels/Magic.panel $panel.$bamlist.ind
verify_length_IDs --panel $panel --bamlist $bamlist 

make_par $panel $bamlist 
make_3poplist --panel $panel --bamlist $bamlist --outgroup $outgroup & wait

qp3Pop -p $panel.$bamlist.qp3test.par > $panel.$bamlist.qp3test.results & wait

#-----------------------------------------------------------------------------#
# Wollstein
bamlist=B.P.ASM.VMAnc
panel=Jorde_Wollstein_hg19_final_noseconddegree_geno01
outgroup=YRI
sed -i 's/.*://' $panel.$bamlist.ind 
Rscript switchPops_admixtools.r ~/archive/Panels/Magic.panel $panel.$bamlist.ind
verify_length_IDs --panel $panel --bamlist $bamlist 

make_par $panel $bamlist 

make_3poplist --panel $panel --bamlist $bamlist --outgroup $outgroup
qp3Pop -p $panel.$bamlist.qp3test.par > $panel.$bamlist.qp3test.results & wait

#-----------------------------------------------------------------------------#
# h
bamlist=B.P.ASM.VMAnc
panel=h 
outgroup=Yoruba
sed -i 's/\.v:/\.variant:/' $panel.$bamlist.ind 
Rscript switchPops_admixTools.r ~/archive/Panels/Magic.panel $panel.$bamlist.ind
verify_length_IDs --panel $panel --bamlist $bamlist 

make_par $panel $bamlist 

make_3poplist --panel $panel --bamlist $bamlist --outgroup $outgroup 
qp3Pop -p $panel.$bamlist.qp3test.par > $panel.$bamlist.qp3test.results & wait 
