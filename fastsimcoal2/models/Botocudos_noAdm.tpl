//Number of population samples (demes)
4
//Population effective sizes (number of genes)
NAnc
NBotoc
NSurui
NKarit
//Sample sizes, ages and inbreeding
0 0 0
2 7 FIS_Botoc
4 0 FIS_Surui
4 0 FIS_Karit
//Growth rates : negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, growth rate, migr. matrix
7 historical event
TBotKarit 3 3 0 iBotKarit  0 0 instbot
TDivKarit 3 0 1 1          0 0
TBotSurui 2 2 0 iBotSurui  0 0 instbot
TDivSurui 2 0 1 1          0 0
TBotBotoc 1 1 0 iBotBotoc  0 0 instbot
TDivBotoc 1 0 1 1          0 0 
TBotAmr       0 0 0 iBotAnc    0 0 instbot
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 5.788249e-09 OUTEXP 