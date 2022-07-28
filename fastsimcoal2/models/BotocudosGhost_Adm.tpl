//Number of population samples (demes)
5
//Population effective sizes (number of genes)
NAnc
NBotoc
NSurui
NKarit
NGhost
//Sample sizes, ages and inbreeding
0 0 0
2 7 FIS_Botoc
4 0 FIS_Surui
4 0 FIS_Karit
0 0 0
//Growth rates : negative growth implies population expansion
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, growth rate, migr. matrix
14 historical event
TAdmKar   0 3 admAnc.Karit 1          0 0
TAdmSur   0 2 admAnc.Surui 1          0 0
7+TAdmBot 0 1 admAnc.Botoc 1          0 0
7+TAdmGh  4 1 admGhost     1          0 0
TBotKarit 3 3 0            iBotKarit  0 0 instbot
TDivKarit 3 0 1            1          0 0
TBotSurui 2 2 0            iBotSurui  0 0 instbot
TDivSurui 2 0 1            1          0 0
TBotBotoc 1 1 0            iBotBotoc  0 0 instbot
TDivBotoc 1 0 1            1          0 0 
TBotGhost 4 4 0            iBotGhost  0 0 instbot
TDivGhost 4 0 1            1          0 0
TBotAmr   0 0 0            iBotAnc    0 0 instbot
TBotBRA   0 0 0            iBotBRA    0 0 instbot
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 5.788249e-09 OUTEXP 