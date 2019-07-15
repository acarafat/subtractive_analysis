"""
Subtractive Analysis Utilis
Requriement: 
    Python 2.7
Running from commandline: 
    python sa_utilis.py fasta_file sequence_length
Example:
    Running following code on sa_utilis.py will generate a fasta file named 
    CDHit_input.fasta containing sequences that has minimum length of 50 base-pairs.
    `python sa_utilis.py input.fasta 50`

Arafat Rahman
July 15 2019
arafat@nstu.edu.bd
"""
#!/usr/bin/env python

import sys
from Bio import SeqIO

def filter_by_length(file_in, thresehold=100):
    filtered_seqs = []
    for seq_record in SeqIO.parse(file_in, 'fasta'):
        if len(str(seq_record.seq)) >= 100:
            filtered_seqs.append(seq_record)
    SeqIO.write(filtered_seqs, open('../output/CDHit_input.fasta', 'w'), 'fasta')
    

if __name__ == '__main__':
    file_input = sys.argv[1]
    if len(sys.argv) == 3:
        thresh = sys.argv[2]
    filter_by_length(file_input, thresh)    
    
    