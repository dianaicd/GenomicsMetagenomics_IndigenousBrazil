//Number of population samples (demes)
4
//Population effective sizes (number of genes)
NeK
NeBot
NeMa
NeGh
//Samples sizes
4
2
4
0
//Growth rates : negative growth implies population expansion
0
//Number of migration matrices : 0 implies no migration between demes
3
// migration matrix
0.000 0.000 0.000 0.000
0.000 0.000 0.000 0.000
0.000 0.000 0.000 0.000
0.000 0.000 0.000 0.000
// migration matrix
0.000 0.000 0.000 0.000
0.000 0.000 0.000 0.005
0.000 0.000 0.000 0.000
0.000 0.005 0.000 0.000
// migration matrix
0.000 0.000 0.000 0.000
0.000 0.000 0.000 0.000
0.000 0.000 0.000 0.005
0.000 0.000 0.005 0.000
//historical event: time, source, sink, migrants, new size, growth rate, migr. matrix
5 historical event
migGBot         0 0 MigProp1 0           0 1
migGMay         0 0 MigProp2 0           0 2
Kar_Bot         1 0 1        rat_KB      0 0
Gh_Bra          3 0 1        rat_GBr     0 0
ancBr_ancMes    2 0 1        rat_anc_anc 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 2.5e-8
