import sys
from Bio import SeqIO
from box import Box
import numpy as np
import v_gene_finder as VGeneFinder
import trailing_j_gene_finder as TrailingJGeneFinder
import a_score as AScore
import probability_holder as ProbabilityHolder
import hmm as HMM

D_GENE_ACCEPTANCE_TYPE = 1

def run_alignment(record):

    sequence = record.seq

    prob_holder = ProbabilityHolder.create()
    best_v_gene = VGeneFinder.get_best_alignment(sequence)

    aln_ums = best_v_gene.aln_ums
    aln_v_gene = best_v_gene.aln_v_gene
    v_seq_name = best_v_gene.v_seq_name
    v_gene_offset = best_v_gene.v_gene_offset

    v_gene = Box({
        'aln_seq': aln_v_gene,
        'name': v_seq_name,
        'offset': v_gene_offset,
        'length': len(aln_v_gene)
    })

    a_score = AScore.get_probability(aln_ums, v_gene)
    ums_no_c = TrailingJGeneFinder.get_best_alignment(sequence)

    hmm = HMM.create_model(v_gene, prob_holder, a_score)
    hmm.bake()
    viterbi_likelihood, viterbi_path = hmm.viterbi(ums_no_c)

    ihhmune_result = Box({
        'log_probability': viterbi_likelihood,
        'state_path': [s[1].name for s in viterbi_path[1:]],
        'v_aln_info': best_v_gene.alignment_info
    })

    return ihhmune_result

