//Number of population samples (demes)
7
//Population effective sizes (number of genes)
NeKar
NeSur
NeBot
NeAym
NeMay
NeUPopA
NeUPopM
//Samples sizes and age
4
4
2 8
2
4
0
0
//Growth rates : negative growth implies population expansion
0
0
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new size, growth rate, migr. matrix
10 historical event
split_Kar_Sur            0 1 1          rat_KS       0 0
split_Bot_KarSur         2 1 1          rat_BotKS    0 0
split_Aym                3 1 1          rat_Aym      0 0
split_May                4 1 1          rat_May      0 0
split_UPopM              6 4 1          rat_UPopM    0 0
split_UPopA              5 1 1          rat_UPopA    0 0
admix_MayUPopA           5 4 MigProp1   1            0 0
admix_BraUPopM1          6 1 MigProp2   1            0 0
admix_BraUPopM2          6 1 MigProp3   1            0 0
admix_BraAymara          3 1 MigProp4   1            0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 2.5e-8
