from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.pairwise2 import format_alignment
from box import Box

OPEN_GAP_COST = -10.0
EXTEND_GAP_COST = -0.5

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
    alignment = pairwise2.align.localxs(sequence, v_sequence, open=OPEN_GAP_COST, extend=EXTEND_GAP_COST, one_alignment_only=True)[0]
    alignment_object = get_alignment_info(*alignment)
    return alignment_object
