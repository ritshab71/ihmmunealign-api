import sys
import alignment_thread as AlignmentThread
from Bio import SeqIO

def multi_cell_align(input_file, is_cdrs):
    fasta_input = list(SeqIO.parse(input_file, 'fasta'))
    results = AlignmentThread.run_alignment(fasta_input[0])
    return results

