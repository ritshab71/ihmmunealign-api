import math
import re
import sys
from Bio import SeqIO
from box import Box
import v_gene_finder as VGeneFinder
import trailing_j_gene_finder as TrailingJGeneFinder
import a_score as AScore

MIN_PROBABILITY_LIMIT = 9.9999999999999995E-7

def is_double(item):
    int_pattern = re.compile("^[0-9]*$")
    float_pattern = re.compile("^[0-9]*.[0-9]*$")
    if float_pattern.match(item) or int_pattern.match(item):
        return True
    else:
        return False

def calculate_normal_distribution_probability(mean, sd, x):
    var = float(sd) ** 2
    denom = (2 * math.pi * var) ** 0.5
    num = math.exp(-(float(x) - float(mean)) ** 2/(2 * var))
    return num/denom

def get_exo_probabilities(mean, std_dev, is_reverse_array):
    i = 0
    prob = 0
    exo_probs = []

    while(not(prob < MIN_PROBABILITY_LIMIT and i > mean)):
        prob = calculate_normal_distribution_probability(mean, std_dev, i)
        exo_probs.append(prob)
        i += 1

    sum_probs = sum(exo_probs)
    normalised_probs = [(prob/sum_probs) for prob in exo_probs]

    return normalised_probs

def create():
    with open('files/iHMMune_align_probabilities.PH') as probabilities_file:
        probabilities_input = probabilities_file.readlines()
        probabilities_input = [line.rstrip() for line in probabilities_input]
        probabilities_input = probabilities_input[1:]

    for entry in probabilities_input:
        entry_array = entry.split()
        identifier = entry_array[0]

        if identifier == 'VD_N':
            vd_n = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'DJ_N':
            dj_n = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'V_end_exo_mean':
            v_end_exo_mean = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'V_end_exo_stdDev':
            v_end_exo_std_dev = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'D_start_exo_mean':
            d_start_exo_mean = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'D_start_exo_stdDev':
            d_start_exo_std_dev = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'D_end_exo_mean':
            d_end_exo_mean = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'D_end_exo_stdDev':
            d_end_exo_std_dev = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'J_start_exo_mean':
            j_start_exo_mean = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'J_start_exo_stdDev':
            j_start_exo_std_dev = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'V_end_P':
            v_end_p = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'D_start_P':
            d_start_p = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'D_end_P':
            d_end_p = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'J_start_P':
            j_start_p = [float(d) for d in entry_array if is_double(d)]
        elif identifier == 'Gene_Mutation':
            gene_mutation = entry_array[1]

    return Box({
        'vd_n': vd_n,
        'dj_n': dj_n,
        'v_end_exo_mean': v_end_exo_mean,
        'v_end_exo_std_dev': v_end_exo_std_dev,
        'd_start_exo_mean': d_start_exo_mean,
        'd_start_exo_std_dev': d_start_exo_std_dev,
        'd_end_exo_mean': d_end_exo_mean,
        'd_end_exo_std_dev': d_end_exo_std_dev,
        'j_start_exo_mean': j_start_exo_mean,
        'j_start_exo_std_dev': j_start_exo_std_dev,
        'v_end_p': v_end_p,
        'd_start_p': d_start_p,
        'd_end_p': d_end_p,
        'j_start_p': j_start_p,
        'gene_mutation': gene_mutation
    })
