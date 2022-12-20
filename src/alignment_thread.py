import sys
from Bio import SeqIO
from box import Box
import numpy as np
from a_score import get_probability
from hmm import create_model
import json
from probability_holder import create
from trailing_j_gene_finder import get_best_alignment as get_best_alignment_j
from v_gene_finder import get_best_alignment as get_best_alignment_v

D_GENE_ACCEPTANCE_TYPE = 1

def run_alignment_file(record):

    sequence = record.seq

    prob_holder = create()
    best_v_gene = get_best_alignment_v(sequence)


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

    a_score = get_probability(aln_ums, v_gene)
    ums_no_c = get_best_alignment_j(sequence)

    hmm = create_model(v_gene, prob_holder, a_score)
    hmm.bake()
    viterbi_likelihood, viterbi_path = hmm.viterbi(ums_no_c)

    ihhmune_result = Box({
        'log_probability': viterbi_likelihood,
        'state_path': [s[1].name for s in viterbi_path[1:]],
        'v_aln_info': best_v_gene.alignment_info
    })

    return ihhmune_result

def run_alignment_sequence(sequence):
    viterbi_likelihood=None
    a_score=None
    best_j_gene=None
    best_v_gene=None

    try:
        prob_holder = create()
        best_v_gene = get_best_alignment_v(sequence)

        # print(best_v_gene)

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
        a_score = get_probability(aln_ums, v_gene)

    # print(f'a_score={a_score}')
        best_j_gene = get_best_alignment_j(sequence)
        # print(json.dumps(best_j_gene, sort_keys=True, indent=4))
        ums_no_c = best_j_gene.ums_no_c

    # print(json.dumps(best_j_gene, sort_keys=True, indent=4))

        hmm = create_model(v_gene, prob_holder, a_score)
        hmm.bake()
        viterbi_likelihood, viterbi_path = hmm.viterbi(ums_no_c)

        ihhmune_result = Box({
            'log_probability': viterbi_likelihood,
            'a_score': a_score,
            # 'state_path': [s[1].name for s in viterbi_path[1:]],
            'v_aln_info': best_v_gene,
            'j_aln_info': best_j_gene
        })

        return ihhmune_result
    except RuntimeError as e:
        print(e)
        ihhmune_result = Box({
            'log_probability': viterbi_likelihood,
            'a_score': a_score,
            # 'state_path': [s[1].name for s in viterbi_path[1:]],
            'v_aln_info': best_v_gene,
            'j_aln_info': best_j_gene
        })

        return ihhmune_result

