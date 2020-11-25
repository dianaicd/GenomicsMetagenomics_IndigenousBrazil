#!/bin/bash
#SBATCH --account amalaspi_virome
#SBATCH --job-name DIAMOND_virome_complete
#SBATCH --partition ax-normal
#SBATCH --time 20:00:00
#SBATCH --mem-per-cpu 10G
#SBATCH --nodes 1
#SBATCH --cpus-per-task 4
#SBATCH --error DIAMOND_output/DIAMOND_job%j.txt
#SBATCH --output DIAMOND_output/DIAMOND_job%j.txt
#SBATC --constraint=AVX512

module add UHTS/Analysis/samtools/1.8
#diamond=/scratch/axiom/FAC/FBM/DBC/amalaspi/popgen/yarizmen/virome/DIAMOND/diamond-0.9.22/bin/diamond


function run_DIAMOND(){
	#Variables: mode (blastp/blastx), cores, name of the database, name of the query file, name of the output taxonomy file, and name of the output classification file
	mode=$1
        cpu=$2
	dbase=$3
        query=$4
        out_tax=$5
        out_alg=$6

	diamond $mode -p $cpu -d $dbase -q $query -o $out_tax -f 102
	diamond $mode -p $cpu -d $dbase -q $query -o $out_alg -f 6 qseqid sseqid staxids pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore
}

function classified_unclassified(){
        #Variables: names of the classification file and the report file.
        classif=$1

        #Classified reads.
        cut -f2 $classif | grep -wcv '0' > $(echo ${classif%_tax.txt}_ClassUnc.txt)
        #Unclassifed reads.
        cut -f2 $classif | grep -wc '0' >> $(echo ${classif%_tax.txt}_ClassUnc.txt)
        #Number of reads per taxid.
        cut -f2 $classif | grep -wv '0' | sort | uniq -c > $(echo ${classif%_tax.txt}_ReadsTaxon.txt)
}

mode=blastx
cpu=5
dbase=/scratch/axiom/FAC/FBM/DBC/amalaspi/virome/yarizmen/virome/Simulations/Viral_tax

#Run DIAMOND
for query in unmapped_fnas/*_unmapped.fna; do
	out_tax=$(echo ${query#unmapped_fnas/})
        out_alg=$(echo ${query#unmapped_fnas/})
        out_tax=$(echo DIAMOND_output/${out_tax%_unmapped.fna}_tax.txt)
        out_alg=$(echo DIAMOND_output/${out_alg%_unmapped.fna}_alg.txt)
        run_DIAMOND $mode $cpu $dbase $query $out_tax $out_alg
done

#wait

#Count
for classif in DIAMOND_output/*_tax.txt; do
        classified_unclassified $classif
done

#wait
