import pysam, getopt, os, sys

bam_in = False
bam_out = False
trim = 0

options, remainder = getopt.getopt(sys.argv[1:],
                                   't:i:o:',
                                   ['trim=', 'bam_in=', 'bam_out='])

for opt, arg in options:
    if opt in ('-i', '--bam_in'):
        bam_in = arg
        if not os.path.isfile(bam_in):
            print("Input bam file does not exist.")
            sys.exist()
    elif opt in ('-o', '--bam_out'):
        bam_out = arg
    elif opt in ('-t', '--trim'):
        trim = int(arg)

if not bam_in:
    sys.exit("Must specify input bam with -b|--bam_in")

if not bam_out:
    bam_out = bam_in.replace(".bam", "_trim" + str(trim) + ".bam")


untrimmed_bam = pysam.AlignmentFile(bam_in, "rb")
trimmed_bam = pysam.AlignmentFile(bam_out, "wb", template = untrimmed_bam)

for read in untrimmed_bam:
    read.qual = '"'*trim + ''.join(list(read.qual)[trim:-trim]) + '"'*trim
    trimmed_bam.write(read)

untrimmed_bam.close()
trimmed_bam.close()