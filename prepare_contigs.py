#!/usr/bin/env python3
import argparse
import pandas as pd
import os
from Bio import SeqIO

__author__ = "Natalia Nenasheva"
__email__ = "nenashen66@uni-greifswald.de"
__status__ = "development"

#####___USAGE_EXAMPLE___#####
# ./prepare_contigs.py -g ./genomes/fly_1500bp_genome_reduced.fasta -a ./annotations/annot_mapped_1500.gtf -p ./pseudogenes/pseudo_mapped_1500.gff3 -c ./test/test_contigs.txt -o ./data_contig_based


parser = argparse.ArgumentParser(description='After Augustify has picked/assigned parameters, it is possible to pull together all sequences with each parameter set & the annotations into separate files ' +
                                             '(3 files per parameter set that was ever assigned: genome, annotation, pseudogenes). ' +
                                             'Then these three sets of files might be used to run BRAKER. ' +
                                             'Thus, this script groups data from fasta and gtf files according to the species type assigned by the Augustify.')

parser.add_argument('-g', '--genome', required=True, type=str,
                    help='Genome fasta file (possibly softmasked).')
parser.add_argument('-a', '--annotation', required=True, type=str,
                    help='Annotation gtf file.')
parser.add_argument('-p', '--pseudogenes', required=True, type=str,
                    help='File with pseudogenes.')
parser.add_argument('-c', '--contigs', required=True, type=str,
                    help='parameters2contigs.lst from Augustify (parameters2contigs.lst may be got with flag -m).')
parser.add_argument('-o', '--out_dir', required=False, type=str,
                    help='Output directory to save all of the files.')
args = parser.parse_args()

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

if args.out_dir:
    directory = args.out_dir
else:
    directory = "data_contig_based"

if not os.path.exists(directory):
    os.makedirs(directory)

d = get_dict(contigs_txt=args.contigs)
lst = [key for key in d]
print('Species for braker:')
print(*lst, sep=' ')
get_fasta(fasta_inp=args.genome, dict=d)
get_gff(gtf_inp=args.annotation, dict=d, type='annot')
get_gff(gtf_inp=args.pseudogenes, dict=d, type='pseudo')

