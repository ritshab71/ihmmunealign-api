from numbers import Number
from operator import attrgetter
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.pairwise2 import format_alignment
from box import Box
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62

OPEN_GAP_COST = -100.0
EXTEND_GAP_COST = -5

# source: ---
def get_alignment_info(align1, align2, score, begin, end):
    gaps = 0
    matches = 0
    mismatches = 0

    start1 = (len(align1[:begin]) - align1[:begin].count("-") + 1)
    start2 = (len(align2[:begin]) - align2[:begin].count("-") + 1)

    for a, b in zip(align1[begin:end], align2[begin:end]):
        if (a == b):
            matches += 1
        elif a == "-" or b == "-":
            gaps += 1
        else:
            mismatches += 1

    return Box({
        'start1': start1 - 1,
        'start2': start2 - 1,
        'score': score,
        'matches': matches,
        'gaps': gaps,
        'mutations': mismatches
    })

def perform_local_alignment(sequence, v_sequence):
    # alignmentv1 = pairwise2.align.localds(sequence, v_sequence, matrix, open=OPEN_GAP_COST, extend=EXTEND_GAP_COST, one_alignment_only=True)
    # print(alignmentv1)
    alignment = pairwise2.align.localds(sequence, v_sequence, matrix, open=OPEN_GAP_COST, extend=EXTEND_GAP_COST)

    alignment_object = get_alignment_info(*alignment[0])
    alignment_format = format_alignment(*alignment[0])

    return Box({
        'aln_object': alignment_object,
        'aln_format': alignment_format
    })


# def main():
#     a = perform_local_alignment('acgctatggacgtctggggccaagggaccacggtcaccgtctcctca'.upper(), 'acggtatggacgtctggggccaagggaccacggtcaccgtctcctca'.upper())
#     # print(a.aln_format)
# if __name__ == "__main__":
#     main()