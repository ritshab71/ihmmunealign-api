import sys
from Bio import SeqIO
from box import Box
import math

EXPONENTIAL_DECAY_RATE = -0.0024;
AVERAGE_DGENE_LENGTH = 24;

def calculate_exponential_decay(position):
    return math.exp(EXPONENTIAL_DECAY_RATE * position)

def exponential_decay_v_gene(v_gene_nucl_position, v_gene_start_offset):
    seq_position = v_gene_nucl_position + v_gene_start_offset
    return calculate_exponential_decay(seq_position);

def exponential_decay_d_gene(d_gene_nucl_position, v_gene_length):
    seq_position = v_gene_length + d_gene_nucl_position
    return calculate_exponential_decay(seq_position);

def exponential_decay_j_gene(j_gene_nucl_position, v_gene_length):
    seq_position = v_gene_length + j_gene_nucl_position
    return calculate_exponential_decay(seq_position);

