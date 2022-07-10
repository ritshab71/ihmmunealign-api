import sys
from Bio import SeqIO
from box import Box
import v_gene_finder as VGeneFinder
import trailing_j_gene_finder as TrailingJGeneFinder
import a_score as AScore
import probability_holder as ProbabilityHolder
import exponential_decay as ExponentialDecay
from pomegranate import *

WAN_NUCLEOTIDE_COVERAGE = 0.0725
WAN_MUTATION_PROBABILITY = 0.138
WAN_MUTABILITY_SCORE = WAN_MUTATION_PROBABILITY/WAN_NUCLEOTIDE_COVERAGE

NO_HOTSPOT_NUCLEOTIDE_COVERAGE = 0.854
NO_HOTSPOT_MUTATION_PROBABILITY = 0.610
NO_HOTSPOT_MUTABILITY_SCORE = NO_HOTSPOT_MUTATION_PROBABILITY/NO_HOTSPOT_NUCLEOTIDE_COVERAGE
NAN_MUTABILITY_SCORE = (WAN_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2

RGYW_NUCLEOTIDE_COVERAGE = 0.0381
RGYW_MUTATION_PROBABILITY = 0.13
RGYW_MUTABILITY_SCORE = RGYW_MUTATION_PROBABILITY / RGYW_NUCLEOTIDE_COVERAGE

NGYW_MUTABILITY_SCORE = (RGYW_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2
RGNN_MUTABILITY_SCORE = (RGYW_MUTABILITY_SCORE * 0.25) + (NO_HOTSPOT_MUTABILITY_SCORE * 0.75)
RGYN_MUTABILITY_SCORE = (RGYW_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2

WRCY_NUCLEOTIDE_COVERAGE = 0.0353
WRCY_MUTATION_PROBABILITY = 0.122
WRCY_MUTABILITY_SCORE = WRCY_MUTATION_PROBABILITY/WRCY_NUCLEOTIDE_COVERAGE

WRCN_MUTABILITY_SCORE = (WRCY_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2
NNCY_MUTABILITY_SCORE = (WRCY_MUTABILITY_SCORE * 0.25) + (NO_HOTSPOT_MUTABILITY_SCORE * 0.75)
NRCY_MUTABILITY_SCORE = (WRCY_MUTABILITY_SCORE + NO_HOTSPOT_MUTABILITY_SCORE)/2

def initialise_mutability_dictionary():
    mutability_dict = {}
    with open('files/mutation_spectrum.txt') as mutation_spectrum:
        for line in mutation_spectrum:
            (key, val) = line.split()
            mutability_dict[key] = float(val)

    return mutability_dict

def test_wan(pent_nucl):
    second_nucl = pent_nucl[1]

    if (second_nucl == 'W'):
        return WAN_MUTABILITY_SCORE
    elif (second_nucl == 'N'):
        return NAN_MUTABILITY_SCORE
    else:
        return NO_HOTSPOT_MUTABILITY_SCORE

def test_rgyw(pent_nucl):
    second_nucl = pent_nucl[1]
    fourth_nucl = pent_nucl[3]
    fifth_nucl = pent_nucl[4]

    if (second_nucl == 'R' and fourth_nucl == 'Y' and fifth_nucl == 'W'):
        return RGYW_MUTABILITY_SCORE
    elif (second_nucl == 'N'):
        if (fourth_nucl == 'Y' and fifth_nucl == 'W'):
            return NGYW_MUTABILITY_SCORE
    elif (fourth_nucl == 'N'):
        if (second_nucl == 'R'):
            return RGNN_MUTABILITY_SCORE
    elif (second_nucl == 'R' and fourth_nucl == 'Y' and fifth_nucl == 'N'):
        return RGYN_MUTABILITY_SCORE
    else:
        return NO_HOTSPOT_MUTABILITY_SCORE

def test_wrcy(pent_nucl):
    first_nucl = pent_nucl[0]
    second_nucl = pent_nucl[1]
    fourth_nucl = pent_nucl[3]

    if (first_nucl == 'W' and second_nucl == 'R' and fourth_nucl == 'Y'):
        return WRCY_MUTABILITY_SCORE
    elif (fourth_nucl == 'N'):
        if (first_nucl == 'W' and second_nucl == 'R'):
            return WRCN_MUTABILITY_SCORE
    elif (second_nucl == 'N'):
        if (fourth_nucl == 'Y'):
            return NNCY_MUTABILITY_SCORE
    elif (first_nucl == 'N' and second_nucl == 'R' and fourth_nucl == 'Y'):
        return NRCY_MUTABILITY_SCORE
    else:
        return NO_HOTSPOT_MUTABILITY_SCORE

def get_mutability_score(nucl_pos, v_gene):
    nucl_index = nucl_pos - 1
    max_index = v_gene.length - 1
    v_seq = v_gene.aln_seq

    if (nucl_index == 0):
        penta_nucl = f'UU{v_seq[0:3]}'
    elif (nucl_index == 1):
        penta_nucl = f'U{v_seq[0:4]}'
    elif (nucl_index == max_index):
        penta_nucl = f'U{v_seq[max_index-2:max_index+1]}UU'
    elif (nucl_index == max_index - 1):
        penta_nucl = f'{v_seq[max_index-3:max_index+1]}U'
    else:
        penta_nucl = f'{v_seq[nucl_index-2:nucl_index+3]}'

    center_nucl = penta_nucl[2]

    if (center_nucl == 'A'):
        return test_wan(penta_nucl)
    if (center_nucl == 'G'):
        return test_rgyw(penta_nucl)
    if (center_nucl == 'C'):
        return test_wrcy(penta_nucl)
    else:
        return NO_HOTSPOT_MUTABILITY_SCORE

def get_match(tri_nucl, mutation_nucl):
    lookup_key = f'{tri_nucl[0]}({tri_nucl[1]}->{mutation_nucl}){tri_nucl[2]}'
    mutability_dict = initialise_mutability_dictionary()
    return mutability_dict[lookup_key]

def find_tri_nucleotide_probability(tri_nucl, mutation_nucl):
    curr_prob = get_match(tri_nucl, mutation_nucl)
    return curr_prob

def find_end_tri_nucleotide_probability(tri_nucl, mutation_nucl):
    total_prob = 0
    for nucl in ['A', 'C', 'G', 'T']:
        unique_tri_nucl = tri_nucl.replace('N', nucl)
        curr_prob = get_match(unique_tri_nucl, mutation_nucl)
        total_prob += curr_prob

    avg_prob = total_prob/4
    return avg_prob

def get_tri_nucleotide_probability(seq, nucl, nucl_pos, v_gene):
    nucl_index = nucl_pos - 1
    max_index = len(seq) - 1

    if (nucl_index == 0):
        tri_nucl = f'N{seq[0:2]}'
    elif (nucl_index == max_index):
        tri_nucl = f'{seq[nucl_index-1:nucl_index+1]}N'
    else:
        tri_nucl = f'{seq[nucl_index-1:nucl_index+2]}'

    if (tri_nucl[1] == nucl):
        print('Oh no its -1!')
        return -1

    result_prob = -1

    if (tri_nucl[0] == 'N' or tri_nucl[2] == 'N'):
        if (tri_nucl[1] == 'N'):
            return 0
        result_prob = find_end_tri_nucleotide_probability(tri_nucl, nucl)
    else:
        if (tri_nucl[1] == 'N'):
            return 0
        else:
            result_prob = find_tri_nucleotide_probability(tri_nucl, nucl)

    if (result_prob == -1):
        return -1
    else:
        return result_prob