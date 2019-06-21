
# Script to sun Dstats
#-----------------------------------------------------------------------------#
panel=Maanasa_mask1_flip
file=B.P.ASM.VMAnc

ops=(genotypename snpname indivname popfilename)
values=(eigenstratgeno snp modified.ind list_qpDstat)

if [ -e ${panel}.$file.qp3test.par ] ; then rm ${panel}.$file.qp3test.par ; fi
for i in $(seq 0 3)
do
    echo "${ops[$i]}: ${panel}.${file}.${values[$i]}" >> ${panel}.$file.qp3test.par
done

qpDstat -p ${panel}.$file.qp3test.par > $panel.$bamlist.qpDstat.results 

#-----------------------------------------------------------------------------#
bamlist=B.P.ASM.VMAnc
panel=Jorde_Wollstein_hg19_final_noseconddegree_geno01
outgroup=YRI
sed -i 's/.*://' $panel.$bamlist.ind 
Rscript switchPops_admixtools.r ~/archive/Panels/Magic.panel $panel.$bamlist.ind
verify_length_IDs --panel $panel --bamlist $bamlist 

make_par $panel $bamlist 

make_3poplist --panel $panel --bamlist $bamlist --outgroup $outgroup
qpDstat -p $panel.$bamlist.qpDstat.par > $panel.$bamlist.qpDstat.results & wait

#-----------------------------------------------------------------------------#
# h
bamlist=B.P.ASM.VMAnc
panel=h 
outgroup=Yoruba
sed -i 's/\.v:/\.variant:/' $panel.$bamlist.ind 
Rscript switchPops_admixTools.r ~/archive/Panels/Magic.panel $panel.$bamlist.ind
verify_length_IDs --panel $panel --bamlist $bamlist 

make_par $panel $bamlist 

manage_Dstat --panel $panel --bamlist $bamlist