from Bio import SeqIO
from Bio.Seq import Seq
import json
from alignment_thread import run_alignment_sequence
from box import Box

# def main():
#     input = list(SeqIO.parse('src/files/PW99.txt', 'fasta'))

#     a = []
#     for seq in input:

#         results = run_alignment_sequence(seq.seq)

#         res_json = Box({
#             'results': results,
#             'seq_name': seq.name
#         })
#         a.append(res_json)

#     print(json.dumps(a, sort_keys=True, indent=4))

def main2():
    results = run_alignment_sequence('tatgctatgcactgggtccgccaggctccaggcaaggggctagagtgggtggcagttatatcatatgatggaagtaataaatactacgcagactccgtgaagggccgattcaccatctccagagacaattccaagaacacgctgtatctgcaaatgaacagcctgagagctgaggacacggctgtgtattactgtgcgag'.upper())
    print(json.dumps(results, sort_keys=True, indent=4))

def main3():
    results = run_alignment_sequence('TCGGAGACCCTGTCCCTCACCTGCGCTGTCTATGGTGGGTCCTTCAGTGGTTACTACTGGAGCTGGATCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCAATCATAGTGGAAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCTGTGTATTACTGTGCGAGAGAAGTTATTAATGACTACGGTGACTTAATCGATGTCTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG')
    print(json.dumps(results, sort_keys=True, indent=4))


main2()
