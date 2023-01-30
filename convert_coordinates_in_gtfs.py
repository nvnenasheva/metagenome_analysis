#!/usr/bin/env python3

from Bio import SeqIO
import subprocess
import argparse

#####___EXAMPLE___#####
# convert_coordinates_in_gtfs.py -g /home/natalia/PycharmProjects/pythonProject/chop_up_genome/fly_15000bp_genome.fasta -a ./annot.gtf -o ./annot_mapped.gtf


#####___INFO___#####
# genome = '/home/natalia/PycharmProjects/pythonProject/chop_up_genome/output_example.fasta'
# ann = './annot.gtf'
# out_bed = "tmp.bed" # Output bed file (intermediate file)
# out_ann = 'annot_mapped.gtf'


parser = argparse.ArgumentParser()

parser.add_argument('-g', '--genome', required=True, type=str,
                    help='Genome fasta file (possibly softmasked)')
parser.add_argument('-a', '--annotation', required=True, type=str,
                    help='Annotation gtf file.')
parser.add_argument('-o', '--out_annot', required=True, type=str,
                    help='Mapped annotation output file.')

args = parser.parse_args()

def find_pattern(s, start, end):
    return (s.split(start))[1].split(end)[0]

def get_bed_file(fasta, bed):
    with open(fasta, "r") as inp, open(bed, 'w') as out:
        for record in SeqIO.parse(inp, "fasta"):
            name_pat = 'source_name='
            coord_start_pat = 'sub_seq_start='
            coord_end_pat = 'sub_seq_end='
            end = ';'
            out.write('%s\t%s\t%s' % (find_pattern(record.id, name_pat, end),
                                      find_pattern(record.id, coord_start_pat, end),
                                      find_pattern(record.id, coord_end_pat, end)))
            out.write("\n")
        print("1. Made bed file for chopped up genome.")


def call_bedtools(bed, ann):
    cmd = "bedtools intersect -b "+ bed + " -a " + ann + " > " + out_ann
    subprocess.call(cmd, shell=True)
    print("2. Intersect bed file (with fragments) and annotation for the genome.")

#=======================================================================================================================
genome = args.genome
ann = args.annotation
out_bed = "tmp.bed"
out_ann = args.out_annot

get_bed_file(genome, out_bed)
call_bedtools(out_bed, ann)
print("3. Made annotation file for chopped up genome: " + str(out_ann))
