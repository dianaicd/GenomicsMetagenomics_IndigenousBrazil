#!/usr/bin/env python
#
# David Poznik
# 2015.7.2
# goQueryPhylogeny.py
#----------------------------------------------------------------------
import argparse, bz2, glob, os, subprocess, sys, re
import cPickle as pickle
from operator import attrgetter, itemgetter
from collections import Counter

DESCRIPTION = '''
Input:  observed Y genotypes
Output: intersection with 1000Y phylogeny
 
The script is designed with very low coverage aDNA in mind, 
but it would work fine with regular samples.

----------------------
Input details
----------------------
1. Query
   A directory of compressed (bz2) files, named and organized as in this example:
     stm/
        stm1.1000Y.genos.txt.bz2
        stm2.1000Y.genos.txt.bz2

   Specification:
   a. Indicate the ID of each individual in the file name, as above.
   b. Each file should be 2 columns.
      i.  hg19/b37 coordinate
      ii. genotype: a single character in { A, C, G, T }
   c. Include all confidently called genotypes (ref and alt) 
      at the 60,555 sites: 1000Y.b37.pos
      It's OK to include any superset of those 60,555 sites 
      (e.g. all those within Y.good2.b37.bed)
   d. For aDNA, use mapDamage2. 
      If depth > 1 at a given site and reads contradict one another, 
      it's better not to call the genotype than to randomly sample an allele.
          
2. Database
   a. (pos, geno) -> (ancSNPtupleList, derSNPtupleList)
          ancSNPtupleList is a list of snpTuples for which (pos, geno) is ancestral
          derSNPtupleList is a list of snpTuples for which (pos, geno) is derived
          snpTuple: (treeLabel, branch, anc, der)
      constructed by: 1000Y/ver5/tree/go9bSnpStats.py -b
   b. (treeLabel, branch) -> branchLabel
      constructed by: 1000Y/ver5/tree/go5aBuildBranchLabels.py

----------------------
Output details
----------------------
1. Summary file
   One line for each sample in each of three sections.
   Example line:
       stm1       A1-V168|1|211 BT-M42|0|79 CT-M168|0|45 F-M89|0|23 GHIJK-M3658|0|2 
                  b.11|39|2 K-M9|0|5 P-M45|0|29 R-M207|0|8 R1-M173|0|10 R1b-M343|0|9 
                  e2.434|8|7
   The data are a series of branchLabel|numAncestral|numDerived strings.
   branchLabel : either a canonical snp or a subtree.branchIndex, 
                 which refers to the 5-page 1000Y phylogeny figure in the supplement
   numDerived  : Number of snps derived genotypes observed from this branch

   Section 1: branches with 2+ observed derived
   Section 2: branches with 1+ observed derived
   Section 3: canonical branches with 1+ observed ancestral
   
2. Details files: 
   Two files for each individual: derived and ancestral hits.
   Columns:
       1. position
       2. ancestral allele
       3. derived allele
       4. observed genotype
       5. mutation type
          * (observed T and other allele is C) or (observed A and other allele is G)
          . otherwise
       6. tree label
       7. branch index
       8. branch label
'''

#----------------------------------------------------------------------
# arguments

parser = argparse.ArgumentParser(description=DESCRIPTION, 
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('inputDir', type=str, help='directory containing files to probe')
args = parser.parse_args()


#----------------------------------------------------------------------
# file names and constants

# input
inputFNstring = '../data/' + args.inputDir + '/*.1000Y.genos.*txt.bz2'
dbFN          = '../data/1000Y.snp.db.pkl'
index2labelFN = '../data/1000Y.index2label.pkl'

# output
outDir        = 'output/%s' % args.inputDir
summaryFN     = '%s/_%s.1000Y.txt' % (outDir, args.inputDir)
ancFNtemplate = '%s/{ID}.1000Y.ancestral.txt' % outDir
derFNtemplate = '%s/{ID}.1000Y.derived.txt' % outDir

# constants
damageTypeSet = set([('C', 'T'), ('G', 'A')])


#----------------------------------------------------------------------
# functions

def runCMD(cmd):
    'spawns an external process'
    if subprocess.Popen(cmd, shell=True).wait():
        sys.exit('\n'*3 + '!'*20 + '\nthe following command FAILED:\n' + cmd)

def mkdirP(dirName):
    'makes a directory'
    cmd = 'mkdir -p %s' % dirName
    runCMD(cmd)
    
def getBranchLabel(branchTuple):
    if branchTuple in branchTuple2labelDict:
        return branchTuple2labelDict[branchTuple]
    else:
        return '%s.%s' % branchTuple
    
def writeSNPlist(outFN, snpList):
    with open(outFN, 'w') as outFile:
        snpList = sorted(snpList, key=attrgetter('branch'), reverse=True)
        snpList = sorted(snpList, key=attrgetter('treeLabel'))
        for snp in snpList:
            outFile.write(str(snp) + '\n')
    print outFN

def compileHitString(ID, isDerived):
    'compile string of hits'
    
    someOutputList = ['%-15s :' % ID]
    allOutputList = list(someOutputList)

    if isDerived:
        branchTupleList = derBranchCounter.keys()
    else:
        branchTupleList = ancBranchCounter.keys()
    branchTupleList = sorted(branchTupleList, key=itemgetter(1), reverse=True) # branch
    branchTupleList = sorted(branchTupleList, key=itemgetter(0))               # tree
    for branchTuple in branchTupleList:
        numDer = derBranchCounter[branchTuple] \
            if branchTuple in derBranchCounter else 0
        numAnc = ancBranchCounter[branchTuple] \
            if branchTuple in ancBranchCounter else 0
        branchLabel = getBranchLabel(branchTuple)
        output = '%s|%d|%d' % (branchLabel, numAnc, numDer)
        
        allOutputList.append(output)
        if (isDerived and numDer > 1) or \
                                (not isDerived and re.search('-', branchLabel)):
            someOutputList.append(output)

    return someOutputList, allOutputList

#----------------------------------------------------------------------
# classes

class SNP:
    def __init__(self, pos, geno, snpTuple):
        self.pos = pos
        self.geno = geno
        self.treeLabel, self.branch, self.anc, self.der = snpTuple
        self.branchLabel = getBranchLabel((self.treeLabel, self.branch))
    def isDerived(self):
        return (self.geno == self.der)
    def isDamageType(self):
        if self.isDerived():
            return (self.anc, self.der) in damageTypeSet
        else:
            return (self.der, self.anc) in damageTypeSet
    def __str__(self):
        mutationType = '*' if self.isDamageType() else '.'
        return '%-8d %s %s %s %s %-3s %3d %s' % \
            (self.pos, self.anc, self.der, self.geno,
             mutationType,
             self.treeLabel, self.branch, self.branchLabel)
        
        
#----------------------------------------------------------------------
print 'Loading database...'

with open(dbFN, 'r') as inFile:
    snpDBdict = pickle.load(inFile)
with open(index2labelFN, 'r') as inFile:
    branchTuple2labelDict = pickle.load(inFile)


#----------------------------------------------------------------------
print '''Probing...\n
SNP files (see --help for format details)'''

mkdirP(outDir)
summaryFile = open(summaryFN, 'w')
derWSoutputList, globalAncOutputList = list(), list()

for inFN in glob.iglob(inputFNstring):
    ID = os.path.basename(inFN).split('.')[0]
    print '%-15s' % ID,

    # read data and build lists of hits    
    derSNPlist, ancSNPlist = list(), list()
    derBranchTupleList, ancBranchTupleList = list(), list()
    with bz2.BZ2File(inFN, 'rb') as inFile:
        for line in inFile:
            pos, geno = line.strip().split()
            pos = int(pos)
            geno = geno.upper()
            dbKeyTuple = (pos, geno)
            if dbKeyTuple in snpDBdict:
                ancSNPtupleList, derSNPtupleList = snpDBdict[dbKeyTuple]
                for snpTuple in ancSNPtupleList:
                    ancBranchTupleList.append(snpTuple[:2]) # (treeLabel, branch)
                    ancSNPlist.append(SNP(pos, geno, snpTuple))
                for snpTuple in derSNPtupleList:
                    derBranchTupleList.append(snpTuple[:2]) # (treeLabel, branch)
                    derSNPlist.append(SNP(pos, geno, snpTuple))

    # detail files
    derFN = derFNtemplate.replace('{ID}', ID)
    ancFN = ancFNtemplate.replace('{ID}', ID)
    writeSNPlist(derFN, derSNPlist)
    print ' ' * 15,
    writeSNPlist(ancFN, ancSNPlist)

    # count hits for summary file
    derBranchCounter = Counter(derBranchTupleList)
    ancBranchCounter = Counter(ancBranchTupleList)
    
    derOutputList, derAllOutputList = compileHitString(ID, isDerived = True)
    summaryFile.write(' '.join(derOutputList) + '\n')
    derWSoutputList.append(' '.join(derAllOutputList) + '\n')
    
    ancOutputList, _ = compileHitString(ID, isDerived = False)
    globalAncOutputList.append(' '.join(ancOutputList) + '\n')


#----------------------------------------------------------------------
summaryFile.write('\n')
for line in derWSoutputList:
    summaryFile.write(line)
summaryFile.write('\n')
for line in globalAncOutputList:
    summaryFile.write(line)
summaryFile.close()
print '\nSummary: %s' % summaryFN
