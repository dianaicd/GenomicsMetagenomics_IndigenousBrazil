
import pysam
import argparse
import numpy as np
import datetime
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type = str,
                    help="Path to the bam file from which we should include damage")


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


samfileInput = pysam.AlignmentFile(args.input, "rb")
allPostion = []
beginm1 = 0
endm1 = 0
maxDuplicate = 30
allDuplicate = np.zeros(maxDuplicate)
count = 0
countDuplicate = 0
allEnd = []
allBeginning = []
allDuplicate = []
oldRef = 0
offset = 0
for indexRead, read in enumerate(samfileInput):
    if(read.reference_id > oldRef):
        offset = offset + read.positions[-1]
        oldRef = read.reference_id
    allEnd.append(read.positions[-1] + offset)
    allBeginning.append(read.positions[0] + offset)
    if read.is_duplicate:
        allDuplicate.append(indexRead)



for indexDuplicate in allDuplicate:
    readBegining = allBeginning[indexDuplicate]
    readEnd = allEnd[indexDuplicate]
    # windowSize = 1
    # foundParent = False
    # while not foundParent:
    #     if (readBegining ==allBeginning[indexDuplicate - windowSize]):
    #         foundParent = True
    #         allParent.append(indexDuplicate - windowSize)
    #
    #
    #                     (readBegining == allBeginning[indexDuplicate + windowSize]) or \
    #                     (readEnd == allEnd[indexDuplicate - windowSize]) or \
    #                     (readEnd == allEnd[indexDuplicate + windowSize])
    #     windowSize = windowSize + 1
    # print("found one parent")





_ , countRead = np.unique(allEnd,return_counts=True)
value, count = np.unique(countRead, return_counts=True)
print(value,count)
print(np.sum((value - 1)*count))
print(f"counter: {count}")
allPosition = np.array(allPostion)
print(f"I have {countDuplicate} duplicate")
print("printing stuff")
print(allDuplicate)
#plt.bar(range(maxDuplicate) ,allDuplicate)
#plt.yscale("log")
#plt.show()
plt.plot(allBeginning)
plt.show()
