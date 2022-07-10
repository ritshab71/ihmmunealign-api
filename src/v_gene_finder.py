import sys
from Bio import SeqIO, Align
import alignment as Alignment
from box import Box

MIN_UMS_ALIGNMENT_END_OFFSET = 12

def get_best_alignment(sequence):
    v_genes = list(SeqIO.parse('files/IGHV_repertoire.fa', 'fasta'))

    best_score = float('-inf')
    for gene in v_genes:
        v_sequence = gene.seq
        # returns the optimal alignment
        curr_alignment = Alignment.perform_local_alignment(sequence, v_sequence)

        if (curr_alignment.score > best_score):
            best_score = curr_alignment.score
            best_alignment = curr_alignment
            best_v_gene = gene

    start1_index = best_alignment.start1
    start2_index = best_alignment.start2

    return Box({
        'sequence': sequence,
        'aln_ums': sequence[start1_index:],
        'v_sequence': best_v_gene.seq,
        'v_seq_name': best_v_gene.name,
        'aln_v_gene': best_v_gene.seq[start2_index:],
        'ums_offset': start1_index,
        'v_gene_offset': start2_index,
        'alignment_info': best_alignment
    })





