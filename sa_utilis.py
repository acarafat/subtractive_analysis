"""
Subtractive Analysis Utilis
"""

from Bio import SeqIO

def filter_by_length(file_in, thresehold=100):
    filtered_seqs = []
    for seq_record in SeqIO.parse(file_in, 'fasta'):
        if len(str(seq_record.seq)) >= 100:
            filtered_seqs.append(seq_record)
    SeqIO.write(filtered_seqs, open('CDHit_input.fasta', 'w'), 'fasta')
    

