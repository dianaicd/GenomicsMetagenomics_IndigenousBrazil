//Number of population samples (demes)
4
//Population effective sizes (number of genes)
NeKar
NeSur
NeBot
NeOut1
//Samples sizes
4
4
2
2
//Growth rates : negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, growth rate, migr. matrix
3 historical event
split_Kar_Sur            0 1 1 rat_KS       0 0
split_Bot_KarSur         2 1 1 rat_BotKS    0 0
split_Out_Bra            3 1 1 rat_OutBra   0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 2.5e-8
