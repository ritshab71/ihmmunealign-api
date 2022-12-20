import sys
from alignment_thread import *
from Bio import SeqIO

def multi_cell_align_file(input_file, is_cdrs):
    fasta_input = list(SeqIO.parse(input_file, 'fasta'))
    results = run_alignment_file(fasta_input[0])
    return results

def multi_cell_align_sequence(input_seq, is_cdrs):
    results = run_alignment_sequence(input_seq)
    return results
