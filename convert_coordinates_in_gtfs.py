#!/usr/bin/env python3

from Bio import SeqIO
import subprocess
import argparse
import json
import pandas as pd
import numpy as np
import os
import glob
import csv

#####___EXAMPLE___#####
# convert_coordinates_in_gtfs.py -headers /home/natalia/PycharmProjects/pythonProject/chop_up_genome/headers.txt -a ./annot.gtf -o ./annot_mapped.gtf -l=${N}


#####___INFO___#####
# headers = '/home/natalia/PycharmProjects/pythonProject/chop_up_genome/headers.txt'
# ann = './annot.gtf'
# out_bed = "tmp.bed"                                           # Output bed file (intermediate file)
# out_ann = 'annot_mapped.gtf'

__author__ = "Natalia Nenasheva"
__email__ = "nenashen66@uni-greifswald.de"
__status__ = "development"

parser = argparse.ArgumentParser(description='convert_coordinates_in_gtfs.py allows ' +
                                 'to transform coordinates of the initial annotation (-a/--annotation): ' +
                                 'the new file corresponds to the coordinate space of the fragmented genome. ' +
                                 'Fragment information is taken from the file -headers/--extended_headers.')

parser.add_argument('-headers', '--extended_headers', required=True, type=str,
                    help='It is an output file from chop_up_genome.py, that contains' +
                         ' extended information about subsequences into which the genome was chopped up.')
parser.add_argument('-a', '--annotation', required=True, type=str,
                    help='Annotation gtf file.')
parser.add_argument('-l', '--fragment_length', required=True, type=int,
                    help='the lengths of new sequences, that you specified in chop_up_genome.py script')
parser.add_argument('-o', '--out_annot', required=True, type=str,
                    help='Mapped annotation output file.')
args = parser.parse_args()


''' ******************* BEGIN FUNCTIONS *************************************'''


def find_pattern(s, start, end):
    #print(s)
    return (s.split(start))[1].split(end)[0]


def get_bed_file(inp_headers, out_bed):
    with open(inp_headers, "r") as inp, open(out_bed, 'w') as out:
        d = {}
        for line in inp:
            name_pat = 'source_name='
            coord_start_pat = 'sub_seq_start='
            key_pat = 'sub_seq_number='
            coord_end_pat = 'sub_seq_end='
            end = ';'

            source_name = find_pattern(line, name_pat, end)
            coord_start = find_pattern(line, coord_start_pat, end)
            coord_end = find_pattern(line, coord_end_pat, end)
            fragment = find_pattern(line, key_pat, end)

            out.write('%s\t%s\t%s' % (source_name, coord_start, coord_end))

            value = str(source_name) + '_' + fragment
            key = str(source_name) + '_' + str(coord_start) + '_' + str(coord_end)
            d[key] = value
            out.write("\n")

        print("Save dictionary with 'SEQNAME_STARTFRAGMENT_ENDFRAGMENT': 'SEQNAME_FRAGMENTNUMBER' information into the file tmp_dict.json.")
        json.dump(d, open("tmp_dict.json", 'w'))
        print("1. Make tmp_dict.json file for chopped up genome.")


def call_bedtools(out_bed, ann, ann_mapped):
    cmd = "bedtools intersect -wb -b "+ out_bed + " -a " + ann + " | sort -k17,17 -k18,18 > " + ann_mapped                                          # add -wb key,  sort -k17,17 -k18,18
    subprocess.call(cmd, shell=True)
    print("2. Intersect bed file (with fragments) and annotation for the genome and made annotation file for chopped up genome: " + str(ann_mapped))


def transform_coordinates(dict_child_names, ann_mapped, fragment_length):
    colnames = ['p_seq', 'p_tool', 'p_hit', 'p_start', 'p_end', 'p_score1', 'p_strand', 'p_score2', 'p_info', 'ch_seq', 'ch_start', 'ch_end']       # ch key is 'child', p is 'parent'
    df = pd.read_csv(ann_mapped, sep='\t', names=colnames)
    d = json.load(open(dict_child_names))

   # add some columns with supporting information
    source_col = df.columns.get_loc('ch_seq')
    df['additional_info'] = df.iloc[:,source_col:source_col+3].astype(str).agg('_'.join, axis=1)
    df['fragment_num'] = np.where(df['ch_end'] % fragment_length == 0, df['ch_end']/fragment_length, (df['ch_end']/fragment_length) + 1)
    df['ch_seq'] = df['additional_info'].map(d)

   # find coordinates for hints from annotation_mapped in fragments of sequences
    df['ch_hit_start'] = (df['p_start'] - ((df['fragment_num'].astype(int) - 1) * (df['ch_end'] - df['ch_start']))).astype(int)  # coordinates in annotations starts from 1, in fasta - from 0!!!
    df['ch_hit_end'] = (df['ch_hit_start'] + (df['p_end'] - df['p_start'] )).astype(int)

   # then save them into the dataframe
    df = df[['ch_seq', 'p_tool', 'p_hit', 'ch_hit_start', 'ch_hit_end', 'p_score1', 'p_strand', 'p_score2', 'p_info']]
    print(df.iloc[3][['p_info']])
    df.to_csv(ann_mapped, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
    print("3. Transform coordinates and rename headers for fragments sequences.")

'''
def rename_headers(input_fasta, output_fasta):
    with open(input_fasta) as inp, open(output_fasta, 'w') as out:
        for seq_record in SeqIO.parse(input_fasta, "fasta"):
            name_pat = 'source_name='
            fragment_pat = 'sub_seq_number='
            end = ';'
            source_name = find_pattern(seq_record.id, name_pat, end)
            fragment = find_pattern(seq_record.id, fragment_pat, end)

            name_to_add = source_name + '_' + fragment
            seq_record.id = name_to_add
            seq_record.description = name_to_add
            SeqIO.write(seq_record, out, 'fasta')
        print("5. Rename headers insight the fasta file with fragments of the initial genome sequences.")
'''

def remove_tmp():
    try:
        for tmp_files in glob.glob('./tmp*'):
            os.remove(tmp_files)
        print("4. Remove temporary files.")
    except:
        print("Error: %s file not found" % tmp_files)

''' ******************* END FUNCTIONS *************************************'''


get_bed_file(inp_headers=args.extended_headers,
             out_bed="tmp.bed")

call_bedtools(out_bed="tmp.bed",
              ann=args.annotation,
              ann_mapped=args.out_annot)


transform_coordinates(dict_child_names="tmp_dict.json",
                      ann_mapped=args.out_annot,
                      fragment_length=args.fragment_length)
remove_tmp()

# rename_headers(input_fasta=args.genome,
#                output_fasta="mapped_genome.fasta")

