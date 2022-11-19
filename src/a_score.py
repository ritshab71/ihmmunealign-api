import sys
from Bio import SeqIO, Align
from box import Box

COMMON_AREA_START_POS = 129
COMMON_AREA_END_POS = 243
COMMON_AREA_NUCL_LENGTH = COMMON_AREA_END_POS - COMMON_AREA_START_POS + 1
A_SLOPE = 0.0065
B_ADDITION = 0.0015
N_NUCLEOTIDE = 'N'
X_NUCLEOTIDE = 'X'

def count_mutations_v_gene(aln_ums, aln_v_gene):
    num_mutations_v_gene = 0
    num_nts_v_gene = 0

    print(f'aln_ums = {len(aln_ums)}')
    print(f'v_gene.aln_seq = {len(aln_v_gene)}')

    for i in range(0, len(aln_v_gene)):
        v_nucl = aln_v_gene[i]
        ums_nucl = aln_ums[i]

        if (v_nucl != ums_nucl):
            if (ums_nucl in [N_NUCLEOTIDE, X_NUCLEOTIDE]):
                num_nts_v_gene += 1
            else:
                num_mutations_v_gene += 1

    return Box({
        'num_mutations_v_gene': num_mutations_v_gene,
        'num_nts_v_gene': num_nts_v_gene
    })

def get_probability(aln_ums, v_gene):
    v_gene_info = count_mutations_v_gene(aln_ums, v_gene.aln_seq)
    cr_v_ratio = 115/len(v_gene.aln_seq) # todo: define 115 real length
    avg_mutations_cr = v_gene_info.num_mutations_v_gene * cr_v_ratio

    num_mutations_cr = 0
    num_nts_cr = 0

    start_index = COMMON_AREA_START_POS - 1 - v_gene.offset
    end_index = COMMON_AREA_END_POS - 1 - v_gene.offset

    for i in range(start_index, end_index):
        v_nucl = v_gene.aln_seq[i]
        ums_nucl = aln_ums[i]

        if (v_nucl != ums_nucl):
            if (ums_nucl in [N_NUCLEOTIDE, X_NUCLEOTIDE]):
                num_nts_cr += 1
            else:
                num_mutations_cr += 1

    prob_mutation_plus_b = (num_mutations_cr * A_SLOPE) * B_ADDITION
    nts_addition = (avg_mutations_cr/COMMON_AREA_NUCL_LENGTH) * A_SLOPE * num_nts_cr

    a_prob = prob_mutation_plus_b + nts_addition
    return a_prob