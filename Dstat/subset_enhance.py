import getopt,sys
myChimp = "chimpHg19.variants"
myAfr = "africans_30_geno0.5.tped"
myH3 = "neanderthal.tped"
myOut = "enhanced_sites.txt"

# %%
# Parse arguments
try:
    opts, args = getopt.getopt(sys.argv[1:],"a:h:f:o:",["ancestral=",
    "H3=", "H3_free=", "out="])
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
    elif opt in ("--H3", "-h"):
        myH3 = arg
    # tped with populations that have not admixed with H3 since their split
    elif opt in ("--H3_free", '-f'):
        myAfr = arg
    elif opt in ("--out", "-o"):
        myOut = arg


chimp_vars = {}

def add_key(line):
    if not line.split()[2] == "N":
        chimp_vars["_".join(line.split()[0:2])] = line.split()[2]

print("Reading ancestral variants")
with open(myChimp, 'r') as file:
    [add_key(line) for line in file.readlines()]

# parse ped
def parse_ped(line):
    chr = line.split()[0]
    var_id = line.split()[1]
    pos = line.split()[3]
    allele1 = line.split()[4]
    allele2 = line.split()[5]
    return(chr,var_id,pos,allele1,allele2)

archaic_vars = {}
print("Extracting segregating sites in H3, in derived state")
# see if archaic carries alternative allele
def add_archaic(line):
    chr,var_id,pos,allele1,allele2 = parse_ped(line)
    # We want it to be homozygote
    if allele1 == allele2 and not allele1 == "N":
        id = "_".join([chr,pos])
        # but different from that of the chimp
        if id in chimp_vars.keys() and not allele1 == chimp_vars[id]:
            archaic_vars[id] = allele1

with open(myH3, 'r') as file:
    [add_archaic(line) for line in file.readlines()]

afr_vars = {}
print("Extracting ancestral states in population that has not admixed with H3 since their split")
# check if Africans carry ancestral allele
def is_ancestral(line):
    chr,var_id,pos,allele1,allele2 = parse_ped(line)
    id = "_".join([chr,pos])
    all_alleles = set(line.split()[4:])
    if 'N' in all_alleles:
        all_alleles.remove('N')

    # We want it to be homozygote
    if len(all_alleles) == 1:
        # but equal from that of the chimp
        if id in chimp_vars.keys() and chimp_vars[id] in all_alleles:
            #and different from the archaic
            if id in archaic_vars.keys() and not archaic_vars[id] in all_alleles:
                afr_vars[id] = var_id #all_alleles.pop()

with open(myAfr, 'r') as file:
    [is_ancestral(line) for line in file.readlines()]

print("Writing enhanced variants")
# Save sites
with open(myOut, 'w') as file:
    [file.write(variant + "\n") for variant in afr_vars.values()]
