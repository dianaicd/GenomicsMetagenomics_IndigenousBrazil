#!/bin/bash


# Don't forget to pass the correct --account 
# for each job! (${PI}_${project})
#SBATCH --account=amalaspi_americas

# You must specify a valid email address!
#SBATCH --mail-user=diana-ivette.cruz-davalos@unil.ch

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=fail,end

# Job name
#SBATCH --job-name="shapeIIIITTTT please"

# Runtime and memory
# default: 12h (normal, ax-normal), 5d (long, ax-long)
# limit: 24h (normal, ax-normal), 10d (long, ax-long)
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=4G

# Partition
# Options are:
# cpu, gpu
#SBATCH --partition=cpu

# For parallel jobs
#SBATCH --cpus-per-task=1
##SBATCH --nodes=2
##SBATCH --ntasks=8
##SBATCH --ntasks-per-node=4


#SBATCH --output=logs/shapeit_%A-%a.out 
#SBATCH --error=logs/shapeit_%A-%a.err

# For array jobs
# Array job containing 100 tasks, run max 10 tasks at the same time
##SBATCH --array=1-100%10
#SBATCH --array=0-87
#### Your shell commands below this line ####

#inds=(MN0008 MN00056 KaritianaA KaritianaB MayaG MayaH SuruiA SuruiB)
inds=(A_Yoruba B_Yoruba A_French B_French)
individuals=($(for ind in ${inds[@]} ; do for chr in {1..22} ; do echo $ind; done; done))
chromosomes=($(for ind in ${inds[@]} ; do for chr in {1..22} ; do echo $chr; done; done))

module load gcc/9.3.0 shapeit4/4.2.2

#for CHR in {1..22}; do
ind=${individuals[$SLURM_ARRAY_TASK_ID]}
CHR=${chromosomes[$SLURM_ARRAY_TASK_ID]}



#-----------------------------------------------------------------------------#
# Call raw genotypes
echo "---- I AM CALLING RAW GENOTYPES -----"
minMapQ=30
minBaseQ=20
ref=/work/FAC/FBM/DBC/amalaspi/popgen/shared_ressources/reference_human/hs.build37.1/hs.build37.1.fa
threads=1
bam=bams/$ind.bam
raw_genos=Raw/$ind.chr$CHR.bcf
gq=30

bcftools mpileup -C 50 -q ${minMapQ} -Q ${minBaseQ} -a FMT/DP,SP \
                        -f ${ref} --threads ${threads} \
                        --regions ${CHR} \
                        -Ou  ${bam} | \
                    bcftools annotate -c RPB | \
                    bcftools call --threads ${threads} -c -V indels \
                        -Ob -o ${raw_genos} - 
#-----------------------------------------------------------------------------#
# Filter genotypes by depth and genotype quality
# 1/3 and twice
#less ../../Genotypes/2020_10_16/genomecov/MN0008_L3U_trim2.genomecov |grep -i genome | awk '{sum+=$2*$5;}END{print sum;}' |less
echo "---- I AM CALLING FILTERING GENOTYPES -----"
min_depth=$(grep -i genome genomecov/$ind.genomecov| awk '{sum+=$2*$5;}END{print sum/3;}') 
max_depth=$(grep -i genome genomecov/$ind.genomecov| awk '{sum+=$2*$5;}END{print sum*2;}') 

filtered_genos=Filtered/$ind.chr$CHR.depth.GQ$gq.vcf.gz



bcftools filter \
        --threads ${threads} \
        --exclude "QUAL <= ${gq} "\
        -Ob ${raw_genos} | \
        bcftools filter \
            --threads ${threads} \
            -Ob \
            --exclude "(SUM(DP4) < ${min_depth} | SUM(DP4) > ${max_depth} )" - \
            > ${filtered_genos}

#-----------------------------------------------------------------------------#
#Preparation
# UNPHASED_VCF=$ind.chr$CHR.vcf.gz
# UNPHASED_VCF_NOMAS=no_mas/$ind.chr$CHR.noMultiAllelicSites.vcf.gz
echo "---- I AM PHASING GENOTYPES -----"
UNPHASED_VCF=${filtered_genos}
UNPHASED_VCF_NOMAS=no_mas/$ind.chr$CHR.bcf

GEN_MAP=1000GP_Phase3/genetic_map_chr${CHR}_combined_b37.txt
PHASED_VCF=Phased/$ind.chr$CHR.onlyPhased.vcf

bcftools view -M 2 -O b $UNPHASED_VCF > $UNPHASED_VCF_NOMAS
bcftools index -f $UNPHASED_VCF_NOMAS
REF_VCF=1000Genomes/ALL.chr$CHR.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz

shapeit4.2 --input $UNPHASED_VCF_NOMAS \
    --map $GEN_MAP  --region $CHR \
    --output $PHASED_VCF \
    --sequencing --reference $REF_VCF

#-----------------------------------------------------------------------------#
#indexing
echo "---- I AM MERGING PHASED AND UNPHASED GENOTYPES -----"
bgzip $PHASED_VCF
bcftools index -f $UNPHASED_VCF
bcftools index -f $PHASED_VCF.gz

FINAL_VCF=Final_phased/$ind.chr$CHR.phased.vcf.gz

#Merging phased and unphased vcfs, keeping all unphased sites from the original vcf, but replacing the phased ones.
#bcftools merge --force-samples -Ov $UNPHASED_VCF $PHASED_VCF | awk 'BEGIN {ofs="\t"}    $0 ~ /^#CHROM/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}   $0 !~ /^#/ {  if(substr($11, 1, 3) != "./.") ; $10 = $11 ; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10 }' | bcftools view -O z > $FINAL_VCF

	#Merging phased and unphased vcfs, keeping all unphased sites from the original vcf, but replacing the phased ones.
	bcftools merge --force-samples $UNPHASED_VCF $PHASED_VCF.gz | awk 'BEGIN {OFS="\t"}
        $0 ~ /^#/ && $0 !~ /^#CHROM/ {print $0}
        $0 ~ /^#CHROM/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}
        $0 !~ /^#/ {
            if(substr($11, 1, 3) != "./.")
                $10 = $11
            print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
        }' | bcftools view -O z - > $FINAL_VCF
        
