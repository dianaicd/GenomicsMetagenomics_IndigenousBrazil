import pysam, os, sys, getopt
import numpy as np
import pandas as pd


# %%

def count_alleles_1D(in_vcf, region, inds_pop1, prefix_stats):
    nInd = len(inds_pop1)
    nChr = 2*nInd
    os.makedirs(os.path.dirname(prefix_stats), exist_ok=True)
    command = f'bcftools view {region} {in_vcf}  -i "AN={nChr} && AF=0" |bcftools stats > {prefix_stats}_all_sites.txt'
    os.system(command)
    command =  f'grep "number of records:" {prefix_stats}_all_sites.txt |rev|cut -f-1|rev'
    stream = os.popen(command)
    nHom = stream.read().replace("\n", "")
    der_allele_counts = {}
    der_allele_counts["0"] = nHom

    # subset sites with at least one derived allele among individuals (should be a small set)
    
    command = f' bcftools view --type snps -M2 {in_vcf} {region}  | bcftools query -i "AN={nChr} && AF>0" -f "%AC\n" | sort |uniq -c | sed "s/ \+/ /g"'
    stream = os.popen(command)
    x = stream.read()

    der_allele_counts.update(
        {y.split(" ")[2]:y.split(" ")[1] for y in x.split("\n")[0:nChr]}
    )

    return(der_allele_counts)


def count_alleles_2D(
    in_vcf, region,inds_pop1,inds_pop2,
    prefix_stats,prefix_tmp_files 
    ):
    nInd_pop1 = len(inds_pop1)
    nInd_pop2 = len(inds_pop2)
    nChr = 2 * (nInd_pop1 + nInd_pop2)
    # stats for the (0,0) entry
    os.makedirs(os.path.dirname(prefix_stats), exist_ok=True)
    command = f'bcftools view {region} {in_vcf}  -i "AN={nChr} && AF=0" |bcftools stats > {prefix_stats}_all_sites.txt'
    os.system(command)
    command =  f'grep "number of records:" {prefix_stats}_all_sites.txt |rev|cut -f-1|rev'
    stream = os.popen(command)
    nHom = stream.read().replace("\n", "")
    # initialize 2D sfs with (nrows, ncols)
    # nrows = nChrs in pop2; ncols = nChrs in pop1
    sfs_2D = np.zeros(((nInd_pop2 * 2) + 1, (nInd_pop1 * 2) + 1))
    sfs_2D[0,0] = nHom
    # subset sites with at least one derived allele among individuals (should be a small set)
    os.makedirs(os.path.dirname(prefix_tmp_files), exist_ok=True)
    command = f' bcftools view -Ob -o {prefix_tmp_files}.bcf --type snps -M2 {region}  -i "AN={nChr} && AF>0" {in_vcf}'
    os.system(command)
    os.system(f"tabix -f {prefix_tmp_files}.bcf")
    # simply count alleles
    bcf_in = pysam.VariantFile(f"{prefix_tmp_files}.bcf")  
    for row in bcf_in.fetch():
        n_der_pop1 = sum([allele for ind in inds_pop1 for allele in row.samples[ind]["GT"] ])
        n_der_pop2 = sum([allele for ind in inds_pop2 for allele in row.samples[ind]["GT"] ])
        sfs_2D[n_der_pop2, n_der_pop1] += 1
        #print(f"n_der_pop1: {n_der_pop1}  n_der_pop2: {n_der_pop2}")

    col_names = [f"d0_{i}" for i in range(0, 2*nInd_pop1 + 1)]
    row_names = [f"d1_{i}" for i in range(0, 2*nInd_pop2 + 1)]
    sfs_2D = pd.DataFrame(sfs_2D, columns = col_names, index = row_names)

    return(sfs_2D)

def sfs_to_realSFS_format(sfs, inds_pop1, inds_pop2):
    nInd_pop1 = len(inds_pop1)
    nInd_pop2 = len(inds_pop2)
    new_sfs = [sfs.iloc[i].iloc[j] for j in range(0, nInd_pop1 * 2 + 1) for i in range(0, nInd_pop2 *  2 + 1)]
    return(new_sfs)

def main():
    # %% --------------------------------------------------------------------------
    # Initialize some values
    chromosome = None
    # Parse args

    print('ARGV      :', sys.argv[1:])
    options, remainder = getopt.getopt(sys.argv[1:], 
                                        'i:s:t:r:p:o:',
                                        [
                                            'in_vcf=',
                                            'prefix_stats=',
                                            'tmp_vcf=',
                                            'region=',
                                            'in_popfile=',
                                            'out_file='
                                        ])
    print('OPTIONS   :', options)

    for opt, arg in options:
        if opt in ('-r', '--region'):
            chromosome = arg
        elif opt in ('-i', '--in_vcf'):
            in_vcf = arg
            prefix_stats = f"stats/{in_vcf}.chr{chromosome}.stats"
            prefix_tmp_files=f"tmp_vcf/{in_vcf}.chr{chromosome}"
            out_file=f'sfs/{in_vcf}.sfs'
        elif opt in ('-s', '--prefix_stats'):
            prefix_stats = arg
        elif opt in ('-t', '--tmp_vcf'):
            prefix_tmp_files = arg
        elif opt in ('-p', '--in_popfile'):
            in_popfile = arg
        elif opt in ('-o', '--out_file'):
            out_file = arg

    # %% --------------------------------------------------------------------------

    region = f"-r {chromosome}" if chromosome else ""
    popfile = pd.read_table(in_popfile, header=None, sep = "\s+", names = ["ind", "pop"])
    grouped = popfile.groupby("pop")
    populations =  list( {i:1 for i in popfile["pop"]} )
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    #print(f"popfile: [{popfile}]")
    if len(populations) > 1:
        inds_pop1,inds_pop2 = [list(grouped.get_group(population)["ind"]) for population in populations]
        sfs = count_alleles_2D(
                in_vcf, region, inds_pop1, inds_pop2,
                prefix_stats, prefix_tmp_files 
            )
        sfs.to_csv(out_file, index = True, sep = "\t")
        sfs_angsd_format = sfs_to_realSFS_format(sfs, inds_pop1, inds_pop2)
        with open(out_file+".realsfs", "w") as file:
            file.write(" ".join([str(x) for x  in  sfs_angsd_format]))
    else:
        inds_pop1 =  [list(grouped.get_group(population)["ind"]) for population in populations]
        nInd = len(inds_pop1)
        allele_counts = count_alleles_1D(in_vcf, region, inds_pop1, prefix_stats)
        with open(out_file, "w") as file:
            first_line = "\t".join([f"d0_{i}" for i in range(0, (nInd * 2) + 1)]) + "\n"
            second_line = "\t".join([allele_counts[str(i)] for i in range(0, (nInd * 2) + 1)]) + "\n"
            file.write(first_line)
            file.write(second_line)


# %%
if __name__ == "__main__":
    main()
