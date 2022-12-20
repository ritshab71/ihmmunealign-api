import re
from Bio import SeqIO
from box import Box
from exponential_decay import *
from mutability_score import *
from probability_holder import *
from pomegranate import *
from pomegranate import State

d_genes = list(SeqIO.parse('src/files/IGHD_repertoire.fa', 'fasta'))
j_genes = list(SeqIO.parse('src/files/IGHJ_repertoire.fa', 'fasta'))

def find_family_index(seq_name):
    regex_match = re.search(r"IGH[VDJ](\d{1})", seq_name)
    return int(regex_match.group(1))

def get_relative_mutation_probability(nucleotide, seq, nucl_position, v_gene, mutation_prob):
    six_percent_mutation = mutation_prob * 0.06
    two_percent_mutation = mutation_prob * 0.02

    reduced_mutation = mutation_prob - six_percent_mutation
    curr_mutation_fraction = get_tri_nucleotide_probability(seq, nucleotide, nucl_position, v_gene)
    relative_mutation_prob = curr_mutation_fraction * reduced_mutation * two_percent_mutation

    return relative_mutation_prob

def calculate_gene_emission_probability(nucl, nucl_position, seq, v_gene, a_score, gene_type):
    if (gene_type == 'V'):
        exp_mutation_prob = exponential_decay_v_gene(nucl_position, v_gene.offset)
    elif (gene_type == 'D'):
        exp_mutation_prob = 1.5 * exponential_decay_d_gene(nucl_position, v_gene.length)
    elif (gene_type == 'J'):
        exp_mutation_prob = exponential_decay_j_gene(nucl_position, v_gene.length)

    mutability_score = get_mutability_score(nucl_position, v_gene)
    mutation_prob = a_score * exp_mutation_prob * mutability_score
    no_mutation_prob = 1 - mutation_prob

    return Box({
      'A': no_mutation_prob if ('A' == nucl) else get_relative_mutation_probability('A', seq, nucl_position, v_gene, mutation_prob),
      'C': no_mutation_prob if ('C' == nucl) else get_relative_mutation_probability('C', seq, nucl_position, v_gene, mutation_prob),
      'G': no_mutation_prob if ('G' == nucl) else get_relative_mutation_probability('G', seq, nucl_position, v_gene, mutation_prob),
      'T': no_mutation_prob if ('T' == nucl) else get_relative_mutation_probability('T', seq, nucl_position, v_gene, mutation_prob),
        'N': no_mutation_prob if ('N' == nucl) else get_relative_mutation_probability('N', seq, nucl_position, v_gene, mutation_prob)

    })

def add_silent_states(hmm):
    silent_states = Box()
    for silent_state in ['X1a', 'X1b', 'X2', 'X3', 'X8', 'X9']:
        state = State(None, name=silent_state)
        hmm.add_state(state)
        state_key = f'{silent_state.lower()}'
        silent_states[state_key] = state

    x5a_states = []
    for i in range(len(d_genes)):
        x5a_state = State(None, name=f'X5a{i}')
        hmm.add_state(x5a_state)
        x5a_states.append(x5a_state)

    silent_states['x5a_states'] = x5a_states

    x5b_states = []
    for i in range(len(d_genes)):
        x5b_state = State(None, name=f'X5b{i}')
        hmm.add_state(x5b_state)
        x5b_states.append(x5b_state)

    silent_states['x5b_states'] = x5b_states

    x7a_states = []
    for i in range(len(d_genes)):
        x7a_state = State(None, name=f'X7a{i}')
        hmm.add_state(x7a_state)
        x7a_states.append(x7a_state)

    silent_states['x7a_states'] = x7a_states

    x7b_states = []
    for i in range(len(d_genes)):
        x7b_state = State(None, name=f'X7b{i}')
        hmm.add_state(x7b_state)
        x7b_states.append(x7b_state)

    silent_states['x7b_states'] = x7b_states

    x11a_states = []
    for i in range(len(j_genes)):
        x11a_state = State(None, name=f'X11a{i}')
        hmm.add_state(x11a_state)
        x11a_states.append(x11a_state)

    silent_states['x11a_states'] = x11a_states

    x11b_states = []
    for i in range(len(j_genes)):
        x11b_state = State(None, name=f'X11b{i}')
        hmm.add_state(x11b_state)
        x11b_states.append(x11b_state)

    silent_states['x11b_states'] = x11b_states

    return silent_states

def add_v_states(hmm, v_gene, a_score):
    v_states = []
    family_index = find_family_index(v_gene.name)

    for position, nucl in enumerate(v_gene.aln_seq):
        nucl_position = position + 1
        prob = calculate_gene_emission_probability(nucl, nucl_position, v_gene.aln_seq, v_gene, a_score, gene_type='V')
        dist = DiscreteDistribution({'A': prob.A, 'C': prob.C, 'G': prob.G, 'T': prob.T, 'N': prob.N})
        state = State(dist, name=f'V{position + 1}')

        hmm.add_state(state)
        v_states.append(Box({
            'state': state,
            'family_index': family_index
        }))

    return v_states

def add_d_states(hmm, v_gene, a_score):
    d_states_matrix = []

    for d_seq in d_genes:
        d_seq_name = d_seq.name
        family_index = find_family_index(d_seq_name)
        d_seq = d_seq.seq.upper()

        d_states = []
        for position, nucl in enumerate(d_seq):
            nucl_position = position + 1
            prob = calculate_gene_emission_probability(nucl, nucl_position, d_seq, v_gene, a_score, gene_type='D')
            dist = DiscreteDistribution({'A': prob.A, 'C': prob.C, 'G': prob.G, 'T': prob.T, 'N': prob.N})
            state = State(dist, name=f'D:{d_seq_name}{position + 1}')

            hmm.add_state(state)
            d_states.append(Box({
                'state': state,
                'family_index': family_index
            }))

        d_states_matrix.append(d_states)

    return d_states_matrix

def add_j_states(hmm, v_gene, a_score):
    j_states_matrix = []

    for j_seq in j_genes:
        j_seq_name = j_seq.name
        family_index = find_family_index(j_seq_name)
        j_seq = j_seq.seq.upper()

        j_states = []
        for position, nucl in enumerate(j_seq):
            nucl_position = position + 1
            prob = calculate_gene_emission_probability(nucl, nucl_position, j_seq, v_gene, a_score, gene_type='J')
            dist = DiscreteDistribution({'A': prob.A, 'C': prob.C, 'G': prob.G, 'T': prob.T, 'N': prob.N})
            state = State(dist, name=f'J:{j_seq_name}{position + 1}')

            hmm.add_state(state)
            j_states.append(Box({
                'state': state,
                'family_index': family_index
            }))

        j_states_matrix.append(j_states)

    return j_states_matrix

def add_n_states(hmm, prob_holder):
    vd_n_states = []
    dj_n_states = []

    vd_n_states_len = len(prob_holder.vd_n)
    for i in range(vd_n_states_len):
        dist = DiscreteDistribution({'A': 0.23599999999999999, 'C': 0.2387, 'G': 0.37130000000000002, 'T': 0.15379999999999999})
        vd_n_state = State(dist, name=f'VD_N{i}')
        hmm.add_state(vd_n_state)
        vd_n_states.append(vd_n_state)

    dj_n_states_len = len(prob_holder.dj_n)
    for i in range(dj_n_states_len):
        dist = DiscreteDistribution({'A': 0.2079, 'C': 0.3135, 'G': 0.33000000000000002, 'T': 0.14849999999999999})
        dj_n_state = State(dist, name=f'DJ_N{i}')
        hmm.add_state(dj_n_state)
        dj_n_states.append(dj_n_state)

    return Box({
        'vd_n': vd_n_states,
        'dj_n': dj_n_states
    })

def set_v_transitions(hmm, v_states, x1a, x1b, prob_holder):
    first_v_state = v_states[0].state
    family_index = v_states[0].family_index

    mean = prob_holder.v_end_exo_mean[family_index - 1]
    sd = prob_holder.v_end_exo_std_dev[family_index - 1]
    v_end_exo_probs = get_exo_probabilities(mean, sd, is_reverse_array=False)

    hmm.add_transition(hmm.start, first_v_state, 1)

    v_size = len(v_states) - 1
    curr, next = None, None
    for i in range(v_size):
        curr = v_states[i].state
        next = v_states[i + 1].state

        if (v_size - i <= len(v_end_exo_probs)):
            prob = v_end_exo_probs[v_size - i - 1]
            hmm.add_transition(curr, next, 1 - prob)
            hmm.add_transition(curr, x1b, prob)
        else:
            hmm.add_transition(curr, next, 1)

    hmm.add_transition(next, x1a, 1)

def set_vd_n_transitions_part_one(hmm, vd_n_states, prob_holder, x1a, x1b, x2, x3):
    # todo: toMagicalStateProb needs to be found out for now set to 0.01
    # magical_state_prob = 0.01
    x1a_to_x2_prob = 1
    x1b_to_x2_prob = 1

    hmm.add_transition(x1a, x2, x1a_to_x2_prob)
    hmm.add_transition(x1b, x2, x1b_to_x2_prob)

    x2_to_n1_prob = prob_holder.vd_n[0]
    x2_to_x3_prob = 1 - x2_to_n1_prob

    hmm.add_transition(x2, vd_n_states[0], x2_to_n1_prob)
    hmm.add_transition(x2, x3, x2_to_x3_prob)

    vd_n_size = len(vd_n_states) - 1
    curr, next = None, None
    for i in range(vd_n_size):
        curr = vd_n_states[i]
        next = vd_n_states[i + 1]

        n_to_n_prob = prob_holder.vd_n[i + 1]
        n_to_x3_prob = 1 - n_to_n_prob

        hmm.add_transition(curr, next, n_to_n_prob)
        hmm.add_transition(curr, x3, n_to_x3_prob)

    hmm.add_transition(next, x3, 1)

def set_vd_n_transitions_multi_part_two(hmm, x3, x5b_states):
    x3_to_x5b_prob = 1 / len(x5b_states)

    for x5b_state in x5b_states:
        hmm.add_transition(x3, x5b_state, x3_to_x5b_prob)

def set_dj_n_transitions_multi_part_two(hmm, x9, x11b_states):
    x9_to_x11b_prob = 1 / len(x11b_states)

    for x11b_state in x11b_states:
        hmm.add_transition(x9, x11b_state, x9_to_x11b_prob)

def set_d_transitions_multi(hmm, x5b_states, x5a_states, d_states_matrix, x7a_states, x7b_states, x8, prob_holder):
    for i in range(len(d_states_matrix)):
        d_seq = d_states_matrix[i]
        d_size = len(d_seq)

        d_state = d_states_matrix[i][0].state
        family_index = d_states_matrix[i][0].family_index

        mean_start = prob_holder.d_start_exo_mean[family_index - 1]
        sd_start = prob_holder.d_start_exo_std_dev[family_index - 1]
        d_start_exo_probs = get_exo_probabilities(mean_start, sd_start, is_reverse_array=False)

        mean_end = prob_holder.d_end_exo_mean[family_index - 1]
        sd_end = prob_holder.d_end_exo_std_dev[family_index - 1]
        d_end_exo_probs = get_exo_probabilities(mean_end, sd_end, is_reverse_array=False)

        x5bi_to_x5ai_prob = d_end_exo_probs[0]
        hmm.add_transition(x5b_states[i], x5a_states[i], x5bi_to_x5ai_prob)
        hmm.add_transition(x5b_states[i], d_state, 0)

        x5a_to_d1_prob = 1
        hmm.add_transition(x5a_states[i], d_state, x5a_to_d1_prob)

        d1_to_d2_prob = 1
        hmm.add_transition(d_state, d_states_matrix[i][1].state, d1_to_d2_prob)

        curr, next = None, None
        for j in range(len(d_states_matrix[i]) - 1):
            curr = d_states_matrix[i][j].state
            next = d_states_matrix[i][j + 1].state

            remainder = 1

            if (j < len(d_start_exo_probs)):
                exo_prob = d_start_exo_probs[j]
                hmm.add_transition(x5b_states[i], curr, exo_prob)

            if (j == len(d_start_exo_probs)):
                final_exo_to_d_prob = 1
                hmm.add_transition(x5b_states[i], curr, final_exo_to_d_prob)

            if (d_size - j <= len(d_end_exo_probs)):
                exo_prob = d_end_exo_probs[d_size - j - 1];
                remainder -= exo_prob
                hmm.add_transition(curr, x7b_states[i], exo_prob)

            hmm.add_transition(curr, next, remainder)

        hmm.add_transition(next, x7a_states[i], 1)
        hmm.add_transition(x7a_states[i], x8, 1)
        hmm.add_transition(x7b_states[i], x8, 1)

def set_dj_n_transitions_multi_part_one(hmm, dj_n_states, prob_holder, x7a_states, x7b_states, x8, x9):
    x8_to_n1_prob = prob_holder.dj_n[0]
    x8_to_x9_prob = 1 - x8_to_n1_prob
    hmm.add_transition(x8, dj_n_states[0], x8_to_n1_prob)
    hmm.add_transition(x8, dj_n_states[0], x8_to_x9_prob)

    curr, next = None, None
    for i in range(len(dj_n_states) - 1):
        curr = dj_n_states[i]
        next = dj_n_states[i + 1]
        n_to_n_prob = prob_holder.dj_n[i + 1]
        n_to_x9_prob = 1 - n_to_n_prob

        hmm.add_transition(curr, next, n_to_n_prob)
        hmm.add_transition(curr, x9, n_to_x9_prob)

    final_n_to_x9 = 1
    hmm.add_transition(next, x9, final_n_to_x9)

def set_j_transitions(hmm, x11a_states, x11b_states, j_states_matrix, prob_holder, is_c_region_removed):
    for i in range(len(j_states_matrix)):
        j_seq = j_states_matrix[i]

        j_state = j_states_matrix[i][0].state
        family_index = j_states_matrix[i][0].family_index

        mean_start = prob_holder.j_start_exo_mean[family_index - 1]
        sd_start = prob_holder.j_start_exo_std_dev[family_index - 1]
        j_start_exo_probs = get_exo_probabilities(mean_start, sd_start, is_reverse_array=False)

        x11bi_to_x11ai_prob = j_start_exo_probs[0]
        hmm.add_transition(x11b_states[i], x11a_states[i], x11bi_to_x11ai_prob)
        hmm.add_transition(x11b_states[i], j_state, 0)

        x11a_to_j1_prob = 1
        hmm.add_transition(x11a_states[i], j_state, x11a_to_j1_prob)

        j1_to_j2_prob = 1
        hmm.add_transition(j_state, j_states_matrix[i][1].state, j1_to_j2_prob)

        curr, next, last = None, None, None
        for j in range(len(j_states_matrix[i]) - 2):
            curr = j_states_matrix[i][j].state
            next = j_states_matrix[i][j + 1].state
            last = j_states_matrix[i][j + 2].state

            if (j < len(j_start_exo_probs)):
                exo_prob = j_start_exo_probs[j]
                hmm.add_transition(x11b_states[i], curr, exo_prob)

            if (j == len(j_start_exo_probs)):
                final_exo_to_j_prob = 1
                hmm.add_transition(x11b_states[i], curr, final_exo_to_j_prob)

            remainder = 1
            hmm.add_transition(curr, next, remainder)

        # todo: apply condition of is_c_removed
        hmm.add_transition(next, last, 1)
        hmm.add_transition(last, hmm.end, 1)

def create_model(v_gene, prob_holder, a_score):
    hmm = HiddenMarkovModel()

    silent_states = add_silent_states(hmm)
    x1a = silent_states.x1a
    x1b = silent_states.x1b
    x2 = silent_states.x2
    x3 = silent_states.x3
    x5a_states = silent_states.x5a_states
    x5b_states = silent_states.x5b_states
    x7a_states = silent_states.x7a_states
    x7b_states = silent_states.x7b_states
    x8 = silent_states.x8
    x9 = silent_states.x9
    x11a_states = silent_states.x11a_states
    x11b_states = silent_states.x11b_states

    n_states = add_n_states(hmm, prob_holder)
    vd_n_states = n_states.vd_n
    dj_n_states = n_states.dj_n

    v_states = add_v_states(hmm, v_gene, a_score)
    d_states_matrix = add_d_states(hmm, v_gene, a_score)
    j_states_matrix = add_j_states(hmm, v_gene, a_score)

    set_v_transitions(hmm, v_states, x1a, x1b, prob_holder)

    set_vd_n_transitions_part_one(hmm, vd_n_states, prob_holder, x1a, x1b, x2, x3)
    set_vd_n_transitions_multi_part_two(hmm, x3, x5b_states)

    set_d_transitions_multi(hmm, x5b_states, x5a_states, d_states_matrix, x7a_states, x7b_states, x8, prob_holder)

    set_dj_n_transitions_multi_part_one(hmm, dj_n_states, prob_holder, x7a_states, x7b_states, x8, x9)
    set_dj_n_transitions_multi_part_two(hmm, x9, x11b_states)

    set_j_transitions(hmm, x11a_states, x11b_states, j_states_matrix, prob_holder, is_c_region_removed=False)

    return hmm






