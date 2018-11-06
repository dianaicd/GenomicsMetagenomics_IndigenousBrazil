# Steps to follow in order to obtain the BAM files and statistics from mapping

The final outputs are stored under:

/home/dcruzdva/archive/BAM/Botocudos/2018_10_26

 For every SAMPLE, the script mapping_aDNA.sh was launched
  with the following parameters:
 
  mapping_aDNA.sh -i samples_.txt --ref /archive/unibe/eg/amalaspi/group/genomes/reference_human/hs.build37.1/hs.build37.1.fa -p 16 --depth --mapDamage


The 'sapfo' repository containin 'mapping_aDNA.sh' was cloned
from git@gitlab.isb-sib.ch:sneuensc/sapfo.git

commit a74a32ea6242c7be1a0cfbcbfd6851263e7338da
Author: Samuel Neuenschwander <samuel.neuenschwander@sib.swiss>
Date:   Thu Oct 25 09:01:45 2018 +0200