#!/usr/bin/env python
import pysam
from Bio import SeqIO

# %%
# poly = BED file with ascertained polymorphisms
poly = open('/Users/dcruz/Projects/Botocudos/Files/Quack/2018_10_03/mismatches.bed')
# samfile = target BAM
samfile = pysam.AlignmentFile("/Users/dcruz/Projects/Botocudos/Files/Quack/2018_10_03/Quack_MT_consensus.bam", "rb")
# consensus on the target individual
consensus = SeqIO.parse(open("/Users/dcruz/Projects/Botocudos/Files/Quack/2018_10_03/quack.fasta"),'fasta')
records = list(consensus)


# %%


def isDamaged(variant, ref):
    # Function to determine if a site might be damaged, as compared to the
    # reference base (i.e., the one found in the consensus)
    ref = ref.upper()
    variant = variant.upper()
    damage = {"C": "T", "G": "A"}
    if variant != ref:
        if ref in damage:
            if variant == damage[ref]:
                return ref + "_" + variant
            else:
                return "other"
        else:
            return "other"
    else:
        return False

# %%


def add_read(read, damaged_reads):
    if read in damaged_reads:
        damaged_reads[read] += 1
        return 0
    else:
        damaged_reads[read] = 1
        return 1

# %%


def add_positions(read, pos, damage_type):
    length_read = len(read)
    if pos > length_read/2:
        pos = -length_read + pos + 1

    if damage_type == "C_T":
        if pos in c_to_t:
            c_to_t[pos] += 1
        else:
            c_to_t[pos] = 1
    if damage_type == "G_A":
        if pos in g_to_a:
            g_to_a[pos] += 1
        else:
            g_to_a[pos] = 1
    if damage_type == "other":
        if pos in others:
            others[pos] += 1
        else:
            others[pos] = 1


# %%
counts = {'c_t': 0, 'g_a': 0, 'others':0}
COUNT_C_T = 0
COUNT_G_A = 0
COUNT_OTHERS = 0
gap = 0

v_start = []
v_end = []
v_chr = []
damaged_reads = {}
contam_reads = {}

c_to_t = {}
g_to_a = {}
others = {}
add_gap = {}

# %%
# Loop on the polymorphic sites
for r in poly:
    line = r.split()
    v_chr.append(line[0])
    v_start.append(int(line[1]))
    v_end.append(int(line[2]))
poly.close()

# %%
# Loop on the polymorphic sites
for i in range(0, len(records[0].seq)):

    add_gap[i-gap+1] = gap
    if records[0].seq[i] == "-":
        gap += 1

    add_gap[i] = gap

# %%


def count_mismatches(pileupread, ref, COUNT_C_T, COUNT_G_A, COUNT_OTHERS):
    # avoid weird sites
    if not pileupread.is_del and not pileupread.is_refskip:
        positions = pileupread.alignment.get_reference_positions()
        pos = [i for i, x in enumerate(positions) if x == BEGIN][0]
        variant = pileupread.alignment.query_sequence[
                pileupread.query_position]
        # avoid missing data
        if variant != "N":
            state = isDamaged(variant, ref)
            read = pileupread.alignment.query_name
            if(state == "C_T"):
                COUNT_C_T += add_read(read, damaged_reads)
                add_positions(read, pos, "C_T")
            elif(state == "G_A"):
                COUNT_G_A += add_read(read, damaged_reads)
                add_positions(read, pos, "G_A")
            elif(state == "other"):
                COUNT_OTHERS += add_read(read, contam_reads)
                add_positions(read, pos, "other")
            return [COUNT_C_T,
                    COUNT_G_A,
                    COUNT_OTHERS]

# %%


for i in range(0, len(v_chr)):
    iCHR = 0
    CHR = v_chr[i]
    BEGIN = v_start[i]
    BEGIN = BEGIN - add_gap[BEGIN]
    END = v_end[i]
    END = END - add_gap[END]

    # iterate (2 positions)
    # ATTENTION: pysam.pileup() is 0-based
    for pileupcolumn in samfile.pileup(CHR, BEGIN, END):
        if pileupcolumn.pos == BEGIN:
            # my records were 0-based
            ref = records[iCHR].seq[BEGIN+add_gap[BEGIN]]
            # loop on reads
            counts = [count_mismatches(pileupread, ref, COUNT_C_T, 
                                       COUNT_G_A, COUNT_OTHERS)
                for pileupread in pileupcolumn.pileups]

# %%

print("Number of reads with damage:")
print(len(damaged_reads))
print("Number of C to T:")
print(COUNT_C_T)
print("Number of G to A:")
print(COUNT_G_A)
print("Other reads with mismatches")
print(COUNT_OTHERS)

# %%
f = open("/Users/dcruz/Projects/Botocudos/Files/Quack/2018_10_03/c_to_t.txt",
         "w")
f.write("Position\tFrequency\n")
for key in c_to_t:
    chain = str(key) + "\t" + str(c_to_t[key]) + "\n"
    f.write(chain)
f.close()

# %%
f = open("/Users/dcruz/Projects/Botocudos/Files/Quack/2018_10_03/g_to_a.txt",
         "w")
f.write("Position\tFrequency\n")
for key in g_to_a:
    chain = str(key) + "\t" + str(g_to_a[key]) + "\n"
    f.write(chain)
f.close()

# %%
f = open("/Users/dcruz/Projects/Botocudos/Files/Quack/2018_10_03/others.txt",
         "w")
f.write("Position\tFrequency\n")
for key in others:
    chain = str(key) + "\t" + str(others[key]) + "\n"
    f.write(chain)
f.close()
# %%
for key in add_gap:
    print(add_gap[key])
# %%
samfile.close()
