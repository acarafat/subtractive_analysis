#!/usr/bin/env python

# Ortholog Check for Subtractive Analysis

"""
Subtractive Analysis Utilis
Requriements: 
    Python 2.7
    Biopython 1.68
Running from commandline: 
    `python subtractive_analysis.py input.fasta`
Example:
    Running following code on this script will run query of a set of amino acid
    sequences on NCBI database to find non-gut microbial and non-human orthologous
    sequences. The default parameters are env_nr database in blastp program with 
    0.0001 e-value. 
    `python subtractive_analysis.py input.fasta`

Arafat Rahman
July 15 2019
arafat@nstu.edu.bd
"""

import sys
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

# A generic function for doing BLAST
def do_qblast(query_file, prg = 'blastp', db = 'env_nr', org = '', e=0.0001):
    '''BLAST many sequences and save xml file'''
    query = open(query_file).read()
    blast_handle =  NCBIWWW.qblast(program=prg, database=db, 
                                     sequence=query, 
                                     entrez_query=org, expect=e)
    record_file = open('blast_tmp_record.xml','w')
    record_file.write(blast_handle.read())
    record_file.close()
    print 'BLAST successful for', query_file
    print 'Output saved in "blast_tmp_record.xml"'
    pass
                                             

# Check orthologous gut flora
def check_gut_orthologs(query_file):
    '''Do BLAST, find sequences with no orthologs in database,
    keep those non-orthologs in a list.'''
    do_qblast(query_file, 'blastp', 'env_nr', 'txid408170', 0.0001)
    blast_records = NCBIXML.parse(open('blast_tmp_record.xml'))
    not_gut_orthologs = []            
    for record in blast_records:
        print len(record.descriptions), 'hits for', record.query
        if len(record.descriptions) == 0:
            not_gut_orthologs.append(record.query)
    return not_gut_orthologs
    
    
# Check orthologus human proteins
def check_human_proteome(query_file):
    '''Do BLAST, find sequence with no orthologs in database, 
    keep those non-orthologs in a list.'''
    do_qblast(query_file, 'blastp', 'nr', 'txid9606', 0.0001)
    blast_records = NCBIXML.parse(open('blast_tmp_record.xml'))
    not_human_orthologs = []
    for record in blast_records:
        print len(record.descriptions), 'hits for', record.query
        if len(record.descriptions) == 0:
            not_human_orthologs.append(record.query)
    return not_human_orthologs            


# Filter sequences
def filter_seqs(file_in, seq_list):
    '''Given a fasta file and a list of fasta description.
    If a sequence description in fasta file have match in the list,
    this function will append the sequence record in a new filtered list.'''
    filtered_seqs = []
    for seq_record in SeqIO.parse(file_in, 'fasta'):
        if seq_record.description in seq_list:
            filtered_seqs.append(seq_record)
    return filtered_seqs



if __name__ == '__main__':
    query_in = sys.argv[1]
    
    # Gut ortholog screening
    print 'Checking for gut orthologs ...'
    non_orthologs_list = check_gut_orthologs(query_in)
    print 'Writing non-gut orthologs output ...'
    SeqIO.write(filter_seqs(query_in, non_orthologs_list), 
                open('../output/2.non_gut_orthologs.fasta','w'), 'fasta')
                
    # Human ortholog screening
    print 'Screening human orthologs ...'
    non_human_ortholog_list = check_human_proteome('../output/2.non_gut_orthologs.fasta')
    print 'Writing non-human orthologs output ...'
    SeqIO.write(filter_seqs(query_in, non_human_ortholog_list), 
                open('../output/3.non_human_orthologs.fasta','w'), 'fasta')
    print 'Done.'
    
    
