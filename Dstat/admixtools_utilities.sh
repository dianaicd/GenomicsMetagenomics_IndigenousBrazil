# Utilities for AdmixTools


make_DstatList(){
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
        start=$(echo $i+1 |bc)
        for j in $(seq $start $nPops) ; do
            if [ ! "${pops[$i]}" = "$outgroup" ] ; then 
                if [ ! "${pops[$j]}" = "$outgroup" ] ; then 
                    echo "Botocudos ${pops[$i]} ${pops[$j]} $outgroup" >> $panel.$bamlist.list_qpDstat
                fi
            fi 
        done
    done
}
#-----------------------------------------------------------------------------#
make_parDstat(){
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
    local pargs=(eigenstratgeno snp modified.ind list_qpDstat)

    if [ -e $panel.$bamlist.qpDstat.par ]
    then rm $panel.$bamlist.qpDstat.par 
    fi 

    for i in $(seq 0 3)
    do
        echo "${parops[$i]}: $panel.$bamlist.${pargs[$i]}" \
            >>$panel.$bamlist.qpDstat.par
    done
}

#-----------------------------------------------------------------------------#
manage_Dstat(){
    SHORTOPTS="p:b:c:l:"
    LONGOPTS="panel: bamlist: maxCPU: maxLines:"
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi
    eval set -- "$ARGS"

    local maxCPU=120
    local maxLines=10

    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)     local panel=$2; shift;;
            -b|--bamlist)   local bamlist=$2; shift;;
            -c|--maxCPU)    maxCPU=$2; shift;;
            -l|maxLines)    maxLines=$2; shift;;
        esac 
        shift
    done

    local combFile=$(grep popfilename $panel.$bamlist.qpDstat.par |sed 's/.*://')
    local total=$(wc -l $combFile |cut -f1 -d' ')
    local maxLoop=$(echo "$total / $maxLines" |bc)

    local lo=1
    local hi=10

    if [ ! -d Results ] ; then mkdir Results ; fi
    for i in $(seq 0 $maxLoop)
    do   
        lo=$(echo "$i*$maxLines + 1" |bc)
        hi=$(echo "($i+1)*$maxLines" |bc)
        usedCPU=$(echo "($i +1) % $maxCPU" |bc)

        qpDstat -p $panel.$bamlist.qpDstat.par -l $lo -h $hi > Results/$panel.$bamlist.$lo.$hi.qpDstat.results &
        if [ $usedCPU -eq 0 ]
        then
            echo "Waiting... $i"
            wait
        fi
    done
}

#-----------------------------------------------------------------------------#
wrap_Dstat(){
    SHORTOPTS="p:b:l:"
    LONGOPTS="panel: bamlist: maxLines:"
    ARGS=$(getopt -s bash -o "$SHORTOPTS" -l "$LONGOPTS" -n "ERROR" -- "$@")
    retVal=$?
    if [ $retVal -ne 0 ]; then echo "something went wrong with args"; exit; fi
    eval set -- "$ARGS"

    local maxCPU=120
    local maxLines=10

    while  [ $# -gt 0 ]; do
        case "$1" in
            -p|--panel)     local panel=$2; shift;;
            -b|--bamlist)   local bamlist=$2; shift;;
            -l|maxLines)    maxLines=$2; shift;;
        esac 
        shift
    done

    local combFile=$(grep popfilename $panel.$bamlist.qpDstat.par |sed 's/.*://')
    local total=$(wc -l $combFile |cut -f1 -d' ')
    local maxLoop=$(echo "$total / $maxLines" |bc)

    local lo=1
    local hi=10

    if [ -e $panel.$bamlist.qpDstat.results ]
    then rm $panel.$bamlist.qpDstat.results ; fi
    for i in $(seq 0 $maxLoop)
    do   
        lo=$(echo "$i*$maxLines + 1" |bc)
        hi=$(echo "($i+1)*$maxLines" |bc)
        cat Results/$panel.$bamlist.$lo.$hi.qpDstat.results >> $panel.$bamlist.qpDstat.results
    done
}
#-----------------------------------------------------------------------------#
make_3popList(){
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
make_par3pop(){
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