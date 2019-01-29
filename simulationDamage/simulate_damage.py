# !/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
Created on Thu Jan 24 14:17:43 2019

@author: fmichaud

"""

import pysam
import argparse
import numpy as np
import datetime


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type = str,
                    help="Path to the bam file from which we should include damage")
parser.add_argument("-b", "--bam_file", type = str, default="subsampling.bam",
                    help="Path to the bam file where will be stored the subsampling")
parser.add_argument("-f", "--fastq_file", type = str,default="damage.fastq",
                    help="Path to the fastq file where the damaged read will be stored ")

parser.add_argument("-ct", "--CTOT", type = str, default=None,
                    help="Path to the file giving the C to T frequency switch. No header allowed. First column should be the index and second the probability to have damage")
parser.add_argument("-ga", "--GTOA", type = str,default=None,
                    help="Path to the file giving the G to A frequency switch. No header allowed. First column should be the index and second the probability to have damage")
parser.add_argument("-m", "--mean_length", type = int,default=50,
                    help="Mean Length of the newly created read")
parser.add_argument("-s", "--sigma_length", type = int,default=30,
                    help="Standard deviation for the newly created read")
parser.add_argument("-min", "--minimum_length", type = int,default=30,
                    help="minimum length for the newly created read")
parser.add_argument("-max", "--maximum_length", type = int,default=90,
                    help="maximum length for the newly created read")
parser.add_argument("-n", "--number_read", type = int,default=0,
                    help="Number of read in the bamfile. If 0 (default) is given, will be recomputed, which is time consuming")
parser.add_argument("-ns", "--number_subsamble", type = int,default=None,
                    help="Number of read to subsamble. If None (default) is given, no subsambling will be performed")
parser.add_argument("-sk", "--skip_damage", default=False,
                    help="If True, skip the part about damaging DNA", action='store_true')

args = parser.parse_args()

def getTimeAndPrint():
    now = datetime.datetime.now()
    print(now.strftime("%d %b %H:%M:%S"))
    return(now)

def getTimeAndPrintDifference(previousTime: "datetime.datetime"):
    now = datetime.datetime.now()
    difference = now - previousTime
    print(f'Last operation took {str(difference)}')
    return(now)
def subSample(inputFile: str, outputFile: str, totalSize: int, subsampleSize: int):
    positionSelected = np.zeros(totalSize)
    position1 = np.random.choice(a=totalSize, size=subsampleSize, replace=False)
    positionSelected[position1] = 1
    samfileInput = pysam.AlignmentFile(inputFile, "rb")
    samfileOutput = pysam.AlignmentFile(outputFile, "wb", template = samfileInput)
    for i,read in enumerate(samfileInput):
        if (positionSelected[i]):
            samfileOutput.write(read)
    samfileInput.close()
    samfileOutput.close()


def damageDNA(sequence: 'np.array[str]', Frenquency: 'np.array[float]', inputBase: str, outputBase: str):
    positionT = np.array(np.where(sequence[0:len(Frenquency)] == inputBase)[0])
    subsetFrequency = Frenquency[positionT]
    randomFrequency = np.random.uniform(size=len(subsetFrequency))
    IndexpositionToBeChanged = np.array((subsetFrequency > randomFrequency))
    PositionToChange = positionT[IndexpositionToBeChanged]
    sequence[PositionToChange] = outputBase

def shortenDNA(sequence: 'np.array[str]', mean, sigma, hardcutMin, hardcutMax):
    nbBase = len(sequence)
    var = sigma**2
    theta = var / mean
    k = mean ** 2 / var
    newNbBase = np.random.gamma(theta,k)
    newNbBase = int(newNbBase)
    newNbBase = max([newNbBase,hardcutMin])
    newNbBase = min([newNbBase, hardcutMax,nbBase])
    return(sequence[0:newNbBase])


def shortenAndDamage(inputFile: str, outputFile: str, mean, sigma, hardcutMin, hardcutMax,
                     frequencyCtoT: 'np.array', frequencyGtoA: 'np.array'):
    samfileInput = pysam.AlignmentFile(inputFile, "rb")
    f = open(outputFile, "w")
    for read in samfileInput:
        seq = read.query_sequence
        inputArray = np.array(list(seq))
        inputArray = shortenDNA(inputArray, mean, sigma, hardcutMin, hardcutMax)
        damageDNA(inputArray, frequencyCtoT, "C", "T")
        inputArray = inputArray[::-1]
        damageDNA(inputArray, frequencyGtoA, "G", "A")
        inputArray = inputArray[::-1]
        seq = "".join(inputArray)
        readLength = len(inputArray)
        qual = read.qual[0:readLength]
        f.write(f"@{read.query_name}  length = {readLength}\n")
        f.write(f"{seq}\n")
        f.write(f"+\n")
        f.write(f"{qual}\n")
    f.close()
    samfileInput.close()

def countRead(intputFile: str):
    samfileInput = pysam.AlignmentFile(intputFile, "rb")
    counter = 0
    for _ in samfileInput:
        counter = counter + 1
    return(counter)

tinitial = tlastoperation = getTimeAndPrint()
if args.number_read == 0:
    print("counting the number of read")
    args.number_read = countRead(args.input)
    print(f"{args.number_read} read were detected")
    tlastoperation = getTimeAndPrintDifference(tinitial)

if args.number_subsamble is not None :
    print(f"Reducing the number of read to {args.number_subsamble} and printing to {args.bam_file}")
    subSample(args.input, args.bam_file, args.number_read, args.number_subsamble)
    print(f"Done reducing the number of read ")
    tlastoperation = getTimeAndPrintDifference(tlastoperation)
else:
    print("No subsampling will be performed")
    args.bam_file = args.input

if not args.skip_damage:
    if(args.CTOT == None or args.CTOT == None):
        print("You should provide input file for the frequency of mutation")
    else:
        print(f"Damaging the DNA. Printing result to {args.fastq_file}")
        allFrequency = np.loadtxt(args.CTOT, "float")
        FrenquencyCtoT = allFrequency[:, 1]
        allFrequency = np.loadtxt(args.GTOA, "float")
        FrenquencyGtoA = allFrequency[:, 1]
        shortenAndDamage(args.bam_file, args.fastq_file,
                         args.mean_length, args.sigma_length, args.minimum_length, args.maximum_length,
                         FrenquencyCtoT, FrenquencyGtoA)
        tlastoperation = getTimeAndPrintDifference(tlastoperation)
print("All done")