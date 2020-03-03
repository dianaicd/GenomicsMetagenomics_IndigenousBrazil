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
parser.add_argument("-f", "--fastq_file", type = str,default="duplicate.fastq",
                    help="Path to the fastq file where the duplicate read will be stored ")

parser.add_argument("-m", "--mean_length", type = int,default=0.5,
                    help="Mean number of base to cut")

parser.add_argument("-n", "--number_read", type = int,default=0,
                    help="Number of read in the bamfile. If 0 (default) is given, will be recomputed, which is time consuming")

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


def shortenDNA(sequence: 'np.array[str]', qualityScore: 'str', mean: float):
    positionCut = len(sequence) - np.random.poisson(mean)
    if(positionCut < 20):
        print("Some sequence were kept unchanged since too short")
        newSequence = sequence
        newQualityScore = qualityScore
    else:
        newSequence = sequence[0:positionCut]
        newQualityScore = qualityScore[0:positionCut]
    return((newSequence,newQualityScore))


def shortenAndDuplicateDNA(inputFile: str, outputFile: str, mean: float, nbOfRead: int):
    samfileInput = pysam.AlignmentFile(inputFile, "rb")
    f = open(outputFile, "w")
    nbOfDuplicate = np.random.poisson(1,nbOfRead)
    for index_read, read in enumerate(samfileInput):
        for _ in range(nbOfDuplicate[index_read]):
            seq = read.query_sequence
            inputArray = np.array(list(seq))
            (shrinkSequence,qual) = shortenDNA(inputArray,read.qual, mean)
            seq = "".join(shrinkSequence)
            readLength = len(shrinkSequence)
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


print(f"Creating duplicate. Printing result to {args.fastq_file}")

shortenAndDuplicateDNA(args.input, args.fastq_file,
                       args.mean_length,args.number_read)
tlastoperation = getTimeAndPrintDifference(tlastoperation)
print("All done")