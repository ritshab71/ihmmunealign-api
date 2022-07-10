import sys
from Bio import SeqIO, Align
from box import Box
import src.alignment as Alignment

MIN_UMS_ALIGNMENT_END_OFFSET = 12

def get_best_alignment(sequence):
    j_genes = list(SeqIO.parse('src/files/IGHJ_repertoire.fa', 'fasta'))

    best_score = float('-inf')
    for gene in j_genes:
        j_sequence = gene.seq
        curr_alignment = Alignment.perform_local_alignment(sequence, j_sequence)

        if (curr_alignment.score > best_score):
            best_score = curr_alignment.score
            best_alignment = curr_alignment
            best_j_gene = j_sequence

    start1_index = best_alignment.start1
    start2_index = best_alignment.start2

    matching_j_gene_len = len(best_j_gene) - start2_index
    end_j_gene_match_position = matching_j_gene_len + start1_index

    if end_j_gene_match_position < len(sequence):
        ums_minus_trailing_j_gene = sequence[0:end_j_gene_match_position]
    else:
        ums_minus_trailing_j_gene = sequence

    return ums_minus_trailing_j_gene






