import sys,getopt

options, remainder = getopt.getopt(sys.argv[1:], 'b:p:o:', ['sample_bim=',
                                                            'panel_bim=',
                                                            'output_alleles=',
                                                            'output_bim='])

for opt, arg in options:
    if opt in ('-b', '--sample_bim'):
        sample_bim = arg
    elif opt in ('-p', '--panel_bim'):
        panel_bim = arg
    elif opt in ('-o', '--output_alleles'):
        new_alleles = arg
    elif opt in ('-b', '--output_bim'):
        new_bim = arg

panel_alleles = {}

with open(panel_bim, 'r') as bimPanel:
    for line in bimPanel.readlines():
        chrom,ID,cM,pos,allele1,allele2 = line.replace("\n", "").split("\t")
        panel_alleles[chrom+":"+pos] = line #[allele1, allele2,ID]

def change_allele(allele, ref, alt):
    if allele in [ref,alt]:
        if allele == ref:
            return(ref)
        if allele == alt:
            return(alt)
    else:
        return('0')

with open(new_alleles, 'w') as outAlleles, open(sample_bim, 'r') as bimSample, open(new_bim, 'w') as newBimSample:

    for line in bimSample.readlines():
        chrom,old_ID,cM,pos,allele1,allele2 = line.replace("\n", "").split("\t")
        if chrom+":"+pos in panel_alleles:
            chrom, new_ID, cM, pos, ref_allele, alt_allele = panel_alleles[chrom+":"+pos].replace("\n", "").split("\t")
            new_allele1 = change_allele(allele1, ref_allele, alt_allele)
            new_allele2 = change_allele(allele2, ref_allele, alt_allele)

            myNewLine = "\t".join([new_ID, allele1, allele2, new_allele1, new_allele2]) + "\n"
            outAlleles.writelines(myNewLine)

            myNewLine = "\t".join([chrom, new_ID, cM, pos, allele1, allele2]) + "\n"
            newBimSample.write(myNewLine)
