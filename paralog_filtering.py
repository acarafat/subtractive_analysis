# Paralog Filtering for Subtractive Analysis
# Coded by Arafat
# V0.1
# February 2017

from Bio import SeqIO

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
    clust_rep = unique_in_cluster('../Data/1488219032.result/1488219032.fas.1.clstr.sorted')
    unique_seqs = filter_seqs('../Data/GCF_000283715.1_ASM28371v1_protein.faa', clust_rep)
    SeqIO.write(unique_seqs, open('no_paralogs.fasta', 'w'), 'fasta')