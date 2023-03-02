#!/usr/bin/env python3
import re
import pandas as pd
import os
from Bio import SeqIO

''' ******************* BEGIN FUNCTIONS *************************************'''

def get_dict(contigs_txt):
    # get data from file
    df = pd.read_csv(contigs_txt, sep='\t', names=['header', 'species', 'score'])
    df = df.groupby(['species'])['header'].apply(','.join).reset_index()
    # prepare dictionary
    d = df.set_index('species').to_dict()['header']
    d2 = {key: [str(val) for val in value.split(',')] for key, value in d.items()}
    #print(d2['chicken'])
    print('Got dictionary with contigs.')
    return d2

def get_fasta(fasta_inp, dict):
    # make fasta for each species
    for key, value in dict.items():
        fasta_out = 'data_contig_based/tmp_' + key + '.fasta'
        #print(key, ' : ', value)
        with open(fasta_out, 'w+') as out:
            with open(fasta_inp, 'rt') as inp:
                for seq_record in SeqIO.parse(inp, "fasta"):
                    if (seq_record.id in value):
                        SeqIO.write(seq_record, out, 'fasta')
    print('Got fasta files.')

def get_gff(gtf_inp, dict, type):
    df = pd.read_csv(gtf_inp, sep='\t',
                     names=['header', 'tool', 'hit', 'score_aln', 'score_1', 'score_2', 'strand', 'num', 'info'])
    for key, value in dict.items():
        gtf_out = 'data_contig_based/tmp_' + type + '_' + key + '.gtf'
        #print(key, ' : ', value)
        df_tmp = df.loc[df['header'].isin(value)]
        #print(df_tmp.shape[0])  # no annotation for some contigs?!!???!?!
        df_tmp.to_csv(gtf_out, header=False, sep='\t', index=False)
    print('Got gtf files:', type)

''' ******************* END FUNCTIONS *************************************'''

fasta_inp = './genomes/fly_1500bp_genome_reduced.fasta'
contigs = './test/test_contigs.txt'
ann = './annotations/annot_mapped_1500.gtf'
pseudo = './pseudogenes/pseudo_mapped_15000.gtf'

directory = "data_contig_based"
if not os.path.exists(directory):
    os.makedirs(directory)

d = get_dict(contigs_txt=contigs)
get_fasta(fasta_inp=fasta_inp, dict=d)
get_gff(gtf_inp=ann, dict=d, type='annot')
get_gff(gtf_inp=pseudo, dict=d, type='pseudo')

# don't forget to remove directory if braker done

