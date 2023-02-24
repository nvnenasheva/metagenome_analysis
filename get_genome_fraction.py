#!/usr/bin/env python3

from Bio import SeqIO
import random
import argparse
import math
import subprocess


__author__ = "Natalia Nenasheva"
__email__ = "nenashen66@uni-greifswald.de"
__status__ = "development"


parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='This script reduce the genome of the model organism: take only 25% (--fraction/-f parameter) randomly selected from the genomes.' +
                                             'The program takes a genomic file with sequences as input. The output file contains a set (a fraction of the input set) of randomly' +
                                             ' selected sequences. The names of the selected sequences are saved to a intermediate file tmp_random_headers.txt')
parser.add_argument('-g', '--genome', required=True, type=str,
                    help='Genome fasta file (possibly softmasked)')
parser.add_argument('-o', '--output_fasta', required=True, type=str,
                    help='Genome fasta file that includes subsequences')
parser.add_argument('-f', '--fraction', required=True, type=float, default=0.25,
                    help='Proportion of genomic sequences to be selected.')
args = parser.parse_args()

#####___INFO___#####
# random_headers = 'tmp_random_headers.txt'
# output_fasta = 'random_sequences.fasta'
# genome = 'raw_genome.fasta'

#####___USAGE_EXAMPLE___#####
# ./get_genome_fraction.py -g /home/natalia/PycharmProjects/pythonProject/chop_up_genome/test_genome.fasta -o test_random_sequences.fasta -f=0.2

''' ******************* BEGIN FUNCTIONS *************************************'''
def save_random_headers(random_headers_list, random_headers_file):
    with open(random_headers_file, 'w') as fp:
        rand_headers = sum(random_headers_list, [])
        for item in rand_headers:
            fp.write(item+'\n')

def get_genome_fraction(fasta_inp, random_headers, fasta_out):
    with open(fasta_inp) as inp, open(fasta_out, 'w') as out:
        with open(random_headers, 'a+') as h:
            h.seek(0)
            lines = h.read().splitlines()
            for seq_record in SeqIO.parse(inp, "fasta"):
                if seq_record.id in lines:
                    SeqIO.write(seq_record, out, 'fasta')

def select_headers_randomly(fasta_inp, random_headers_file, fraction):
    source_names = []
    idxs_list = []

    last_index = int(subprocess.check_output(['grep', '-c', '>', fasta_inp]).decode('utf-8')) - 1     # save the index of the last sequence from input fasta file

    n = math.ceil(int(last_index + 1)*fraction)                                                       # take only 25% (fraction parameter) randomly selected from the genomes
    print('Input genome.fasta has {} sequences.'.format(last_index + 1))
    print('Select set of contigs, where fraction={}.'.format(fraction))
    print('-----------------------------------------------------------')

    with open(fasta_inp) as inp:
        idxs = []
        for index, seq_record in enumerate(SeqIO.parse(inp, "fasta")):
            source_name = seq_record.id.split("_")[0]
            idx = seq_record.id.split("_")[1]
            if source_name not in source_names:
                source_names.append(source_name)
                if idxs:
                    idxs_list.append(idxs)
                    idxs = []
                    idxs.append(idx)
                else:
                    idxs.append(idx)
            else:
                if int(index) == last_index:
                    idxs.append(idx)
                    idxs_list.append(idxs)
                else:
                    idxs.append(idx)

        rand_headers = []
        for i, name in zip(idxs_list, source_names):
            if n > len(i):
                n = len(i)
                rand_idxs = i
            else:
                rand_idxs = random.sample(i, n)

            rand_header = [name + '_' + s for s in rand_idxs]
            rand_headers.append(rand_header)

    save_random_headers(random_headers_list=rand_headers,
                        random_headers_file=random_headers_file)
    print("1. Save randomly selected sequences headers (n={}) to file: ".format(n), random_headers)

''' ******************* END FUNCTIONS *************************************'''

random_headers = 'tmp_random_headers.txt'
select_headers_randomly(args.genome, random_headers, args.fraction)
get_genome_fraction(args.genome, random_headers, args.output_fasta)
print("2. Got small fasta file with randomly selected sequences: ", args.output_fasta)
