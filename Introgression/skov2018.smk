localrules: cat_weights, concat_mut
chromosomes = [str(i) for i in range(1,23)]

rule all:
    input:
        #mut = "mutationrates.txt",
        #noAncestral = "SuruiA.chr22.observations_NOancestral.txt",
        #weights = "weights.txt",
        decoded = "SuruiA_decoded.Summary.txt"

rule get_callable:
    input:
        # repeatmaskfile = "masks/chr{chr}.fa.masked",
        # callabilityfile = "20140520.chr{chr}.strict_mask.fasta.gz"
        repeatmaskfile = "masks/chr{chr}.fa.masked",
        callabilityfile = "StrictMask/20140520.chr{chr}.strict_mask.fasta.gz"
    output:
        called_bases = "chr{chr}.txt",
        bed = "chr{chr}.bed"
    resources:
    params:
        windowsize = 1000,
        basename = "chr{chr}"
    shell:
        """
        module load gcc/9.3.0 python/2.7.18
        python MakeMaskfiles.py {input.repeatmaskfile} {input.callabilityfile} {params.windowsize} {wildcards.chr} {params.basename}

        # So for chromosome 17
        # python MakeMask.py chr17.fa.masked 20140520.chr17.strict_mask.fasta.gz 1000 17 chr17
        """

rule cat_weights:
    input:
        called_bases = expand("chr{chr}.txt", chr = chromosomes),
        bed = expand("chr{chr}.bed", chr = chromosomes)
    output:
        weightsTxt = "weights.txt",
        weightsBed = "weights.bed"
    shell:
        """
        cat {input.called_bases} > {output.weightsTxt}
        cat {input.bed} > {output.weightsBed} 
        """

rule variants_outgroup:
    input:
        #inds_ids = "outgroups.txt",
        #vcf = "vcfs/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz",
        vcf = "bcf/mixe_maya_huichol.chr{chr}.bcf",
        bed = "chr{chr}.bed"
    output:
        freqs = "chr{chr}.freq"
    shell:
        """
        #tabix -h ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -B chr17.bed | \
            #vcftools --vcf - --counts --stdout --keep outgroups.txt --remove-indels --min-alleles 2 --max-alleles 2 > chr17.freq

        #conda activate tabix 
        tabix  -h {input.vcf}  -R {input.bed} | \
            vcftools --vcf - --counts --stdout --remove-indels --min-alleles 2 --max-alleles 2 > {output.freqs}

        """

rule get_mu:
    input:
        called_bases = "chr{chr}.txt",
        freqs = "chr{chr}.freq"
    output:
        mut = "chr{chr}.mut"
    params:
        windowsize = 1000000,
        output_size = 1000
    shell:
        """
        module load gcc/9.3.0  python/2.7.18
        python Estimate_mutationrate.py {input.freqs} {params.windowsize} {params.output_size} {input.called_bases} {output.mut}

        # so for chromosome 17
        #python Estimate_mutationrate.py chr17.freq 1000000 1000 chr17.txt chr17.mut
        """

rule concat_mu:
    input:
        mut = expand("chr{chr}.mut", chr = chromosomes)
    output:
        "mutationrates.txt"
    shell:
        """
        cat {input} > {output}
        """

rule variants_training_ancestral:
    input:
        freqs = "chr{chr}.freq",
        train_vcf = "bcf/{ind}_maskStrict_repeatsUCSC_dp_GQ30_chr{chr}.bcf",
        bed = "chr{chr}.bed",
        weightsTxt = "chr{chr}.txt",
        ancestral = "homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_{chr}.fa" 
    output:
        "{ind}.chr{chr}.observations.txt"
    params:
        windowsize = 1000
    shell:
        """
        # tabix -fh 1000_genomes_phase3/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -B chr17.bed | \
        #    vcftools --vcf - --indv HG00096 --remove-indels --min-alleles 2 --max-alleles 2 --stdout --counts | \
        #    python FiltervariantsNOancestral.py chr17.freq 1000 chr17.txt HG00096.chr17.observations_NOancestral.txt

        tabix -fh {input.train_vcf} -R {input.bed} | \
            vcftools --vcf - --remove-indels --min-alleles 2 --max-alleles 2 --stdout --counts | \
                python Filtervariants.py {input.ancestral} {input.freqs} {params.windowsize} {input.weightsTxt} {output}
        """



rule variants_training_NOancestral:
    input:
        freqs = "chr{chr}.freq",
        train_vcf = "bcf/{ind}_maskStrict_repeatsUCSC_dp_GQ30_chr{chr}.bcf",
        bed = "chr{chr}.bed",
        weightsTxt = "chr{chr}.txt",
    output:
        "{ind}.chr{chr}.observations_NOancestral.txt"
    params:
        windowsize = 1000
    shell:
        """
        # tabix -fh 1000_genomes_phase3/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -B chr17.bed | \
        #    vcftools --vcf - --indv HG00096 --remove-indels --min-alleles 2 --max-alleles 2 --stdout --counts | \
        #    python FiltervariantsNOancestral.py chr17.freq 1000 chr17.txt HG00096.chr17.observations_NOancestral.txt

        tabix -fh {input.train_vcf} -R {input.bed} | \
            vcftools --vcf - --remove-indels --min-alleles 2 --max-alleles 2 --stdout --counts | \
                python FiltervariantsNOancestral.py {input.freqs} {params.windowsize} {input.weightsTxt} {output}
        """

rule cat_variants:
    input:
        #lambda wildcards: expand("{ind}.chr{chr}.observations_NOancestral.txt", ind = wildcards.ind, chr = chromosomes)
        lambda wildcards: expand("{ind}.chr{chr}.observations.txt", ind = wildcards.ind, chr = chromosomes)
    output:
        #"{ind}.obs_NOancestral.txt"
        "{ind}.obs.txt"
    shell:
        """
        cat {input} > {output}
        """

rule train_model:
    input:
        #train_obs = "{ind}.obs_NOancestral.txt",
        train_obs = "{ind}.obs.txt",
        hmm_params = "StartingParameters.hmm",
        weightsTxt = "weights.txt",
        mut = "mutationrates.txt"
    output:
        "{ind}_trained.hmm"
    conda:
        "python2.7"
    params:
        out_basename = "{ind}_trained"
    shell:
        """
        module load gcc/9.3.0  python/2.7.18
        python Train.py {input.train_obs} {params.out_basename}  {input.hmm_params} {input.weightsTxt} {input.mut}
        """

rule decode:
    input:
        train_obs = "{ind}.obs.txt",
        #train_obs = "{ind}.obs_NOancestral.txt",
        trained = "{ind}_trained.hmm",
        weightsTxt = "weights.txt",
        mut = "mutationrates.txt"
    output:
        "{ind}_decoded.Summary.txt"
    params:
        out_basename = "{ind}_decoded",
        windowsize = 1000
    shell:
        """
        module load gcc/9.3.0  python/2.7.18
        python Decode.py {input.train_obs} {params.out_basename} {input.trained} {input.weightsTxt} {input.mut} {params.windowsize}
        """
#  module load gcc/9.3.0  python/2.7.18; python Train.py SuruiA.chr22.observations_NOancestral.txt SuruiA_trained StartingParameters.hmm weights.txt mutationrates.txt
#python Decode.py SuruiA.chr22.observations_NOancestral.txt SuruiA_decoded SuruiA_trained.hmm weights.txt mutationrates.txt 1000
