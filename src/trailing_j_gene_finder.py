import sys
from Bio import SeqIO, Align
from box import Box
from alignment import *
import json

MIN_UMS_ALIGNMENT_END_OFFSET = 12

def get_best_alignment(sequence):
    j_genes = list(SeqIO.parse('src/files/IGHJ_repertoire.fa', 'fasta'))

    best_score = float('-inf')
    for gene in j_genes:
        j_sequence = gene.seq
        # print(f'j= {gene.name}')
        curr_alignment_data = perform_local_alignment(sequence, j_sequence)
        # print(json.dumps(curr_alignment_data, sort_keys=True, indent=4))
        # print('\n\n')
        curr_alignment = curr_alignment_data.aln_object

        if (curr_alignment.score > best_score):
            best_score = curr_alignment.score
            best_alignment_data = curr_alignment_data
            best_alignment = curr_alignment
            best_j = gene
            best_j_gene = j_sequence

    # print(json.dumps(best_alignment_data, sort_keys=True, indent=4))

    start1_index = best_alignment.start1
    start2_index = best_alignment.start2

    matching_j_gene_len = len(best_j_gene) - start2_index
    end_j_gene_match_position = matching_j_gene_len + start1_index

    if end_j_gene_match_position < len(sequence):
        ums_minus_trailing_j_gene = sequence[0:end_j_gene_match_position]
    else:
        ums_minus_trailing_j_gene = sequence

    return Box({
        'ums_no_c': str(ums_minus_trailing_j_gene),
        'aln_ums': str(ums_minus_trailing_j_gene[start1_index:]),
        'j_sequence': f'{best_j.seq}',
        'j_seq_name': f'{best_j.name}',
        'aln_j_gene': f'{best_j.seq[start2_index:]}',
        'ums_offset': start1_index,
        'j_gene_offset': start2_index,
        'aln_visual': best_alignment_data.aln_format
    })






