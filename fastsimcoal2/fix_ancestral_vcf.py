from Bio import SeqIO
from pysam import VariantFile

vcf_path="test.vcf.gz"
vcf_path="A_Yoruba.chr22.vcf.gz"
out_vcf_path = 'test_anc.vcf'
ancestral_fasta = "chimpHg19.fa"
chr = "22"
chimp_chrs = SeqIO.to_dict(SeqIO.parse(ancestral_fasta, "fasta"))

# with open(ancestral_fasta) as handle:
#     for record in SeqIO.parse(handle, "fasta"):
#         print(f"This contig contains the ancestral alleles: {record.id}")


bcf_in = VariantFile(vcf_path)  # auto-detect input format
samples = [sample for sample in bcf_in.header.samples]

bcf_out = VariantFile(out_vcf_path, 'w', header=bcf_in.header)

for row in bcf_in.fetch():
    record = chimp_chrs[chr]
    anc_allele = record.seq[ row.pos -1 ].upper()
    new_alleles = row.alleles + (anc_allele,)
    # skip variants without ancestral allele information
    if anc_allele == ".":
        continue
    # skip variants from two or more mutations in homo lineage
    if len(set(new_alleles)) > 2:
        continue
    # 
    row.alleles = new_alleles
    if row.ref != anc_allele:
        # we have to swap alleles;
        # 01 -> 01
        # 00 -> 11
        # 11 -> 00
        for ind in samples:
            genotype = row.samples[ind]["GT"]
            if genotype == (0,0):
                genotype = (1,1)
            elif genotype == (1,1):    
                genotype = (0,0)
            row.samples[ind]["GT"] = genotype

    bcf_out.write(row)

        # print(f"not changing anc_allele: {anc_allele}  ref: {row.ref}  pos: {row.pos} alleles: {row.alleles}")    
    
    #print(f"anc_allele: {anc_allele}  ref: {row.ref}  pos: {row.pos}")


bcf_out.close()
#row.samples["MN0008"]["GT"]