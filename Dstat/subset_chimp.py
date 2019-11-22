from Bio import SeqIO

myBim = "h.bim"
myChimp = "chimpHg19.fa"
myOutput = "chimpHg19.variants"

# %%
# Parse arguments
try:
    opts, args = getopt.getopt(sys.argv[1:],"b:a:o:",["bim=",
    "ancestral=", "out="])
except getopt.GetoptError:
    print('sys.argv[0] -i <inputfile> -o <outputfile>')
    sys.exit(2)
print(opts)
print(args)

for opt, arg in opts:
    # file with variants representing the ancestral alleles
    if opt in ("--ancestral", "-a"):
        myChimp = arg
    # tped with one individual representing H3
    elif opt in ("--bim", "-b"):
        myBim = arg
    elif opt in ("--out", "-o"):
        myOutput = arg


# Get all positions from panel's bim
myPos = {}
def add_key(line):
    chr,id,control,start,ref,alt = line.split("\t")
    myPos[chr+"_"+start] = 1


with open(myBim, "r") as bim:
    [ad_key(line) for line in bim.readlines()]


chimp_vars = {} 
chrs = set([coord.split("_")[0] for coord in list(myPos.keys())])


def write_line(pos):
    myLine = "\t".join([seq_record.id, str(pos),
             seq_record.seq[pos-1]]) + "\n"
    variants.write(myLine)
    #return(myLine)


with open(myOutput, 'w') as variants:
    # Go through fasta (ancestral states) and save variants
    # iterate on chromosomes
    for seq_record in SeqIO.parse(myChimp, "fasta"):
        if seq_record.id in chrs:
            print(seq_record.id)
            # get positions in a given chromosome
            positions = sorted([int(coord.split("_")[1]) for coord in list(myPos.keys()) if coord.split("_")[0] == seq_record.id])
            print(len(positions))
            
            [write_line(pos) for pos in positions]
            
    
    
    #print(len(positions))
    #print(positions[0])
    #print(positions[-1])

    #print(repr(seq_record.seq))

