#-----------------------------------------------------------------------------#
# Get chromosome size and break it in blocks

with open(ref_genome + ".fai", 'r') as index:
    
    chr_size = {}
    for line in index.readlines():
        chr,size = line.split()[0:2]
        chr_size[chr] = size
    chromosomes = list(chr_size.keys())
    chromosomes = [str(x) for x in range(1,23)] 
    lower_chr = {}
    upper_chr = {}
    for chr in chromosomes:
        chr = str(chr)
        lower_chr[chr] = [i for i in range(1, int(chr_size[chr]), int(float(blockSize)))]
        upper_chr[chr] = [i - 1 for i in lower_chr[chr][1:]]
        upper_chr[chr].append(chr_size[chr])

def get_range(chr = False):
    if(chr):
        ranges = ["_".join([str(chr), 
                    str(lower_chr[str(chr)][i]), 
                    str(upper_chr[str(chr)][i])]) for i in range(0, len(lower_chr[str(chr)]))]
    else:
        ranges = ["_".join([str(chr), 
                    str(lower_chr[str(chr)][i]), 
                    str(upper_chr[str(chr)][i])]) for chr in chromosomes for i in range(0, len(lower_chr[str(chr)]))]
    return(ranges)