#!/usr/bin/env python

# Paralog filtering for subtractive analysis

'''This script works on CD-HIT output.
CD-Hit is an web server for homolog sequence clustering. One output file of 
CD-Hit contains paralog sequence id in cluster. This script mine that output file,
find the unique sequence, and filter out those sequences from query in file.

Requriements: 
    Python 2.7
Running from commandline: 
    `python paralog_filtering.py input.fas.1.clstr.sorted CDHit_input.fasta`    
Example:
    CDHit analysis has several outputs. This file work on files with *.fas.1.clstr.sorted 
    extension. This output from CDHit report a cluster of fasta ids from fasta file 
    submited to CDHit. This code takes one sequence id from each paralog cluster, and then 
    extract the sequence of the selected id from fasta file submitted to CDHit.
    
    The following code mines on CDHit_input.fasta using the sequence ids from
    input.fas.1.clstr.sorted file and generate 1.no_paralogs.fasta file as output.
    `python subtractive_analysis.py input.fas.1.clstr.sorted CDHit_input.fasta` 
    
Arafat Rahman
July 15 2019
arafat@nstu.edu.bd
'''

from Bio import SeqIO
import sys

######################################
# Pick single sequence from a cluster#
######################################

def unique_in_cluster(filein):
    '''
    Input: Sorted cluster file of CD-HIT output
    Algo: Read raw data, by spliting, select unique sequence id in a cluster
    Output: List of unique sequence id
    '''
    cdhit_cluster = open(filein)
    raw_data = cdhit_cluster.read()
    
    # List to hold cluster representatvie sequences
    clust_rep = []
    for i in raw_data.split('>Cluster'):
        clust_rep.append(i.split('*')[0].split('>')[-1].split('...')[0])

    if '' in clust_rep:
        clust_rep.remove('')
    return clust_rep


# Filter sequences
def filter_seqs(file_in, seq_list):
    '''Given a fasta file and a list of fasta description.
    If a sequence description in fasta file have match in the list,
    this function will append the sequence record in a new filtered list.'''
    filtered_seqs = []
    for seq_record in SeqIO.parse(file_in, 'fasta'):
        if seq_record.id in seq_list:
            filtered_seqs.append(seq_record)
    return filtered_seqs
    
if __name__ == '__main__':
    clust_rep = unique_in_cluster(sys.argv[1])
    cdhit_input = sys.argv[2]
    unique_seqs = filter_seqs(cdhit_input, clust_rep)
    SeqIO.write(unique_seqs, open('../output/1.no_paralogs.fasta', 'w'), 'fasta')