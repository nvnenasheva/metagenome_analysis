#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import random

#####___TEST___#####
# chop_up_genome.py -g ./example.fasta -o ./output_example.fasta -l 100


#####___USAGE_EXAMPLE___#####
# use lengths = 1500, 3000, 6000, 15000, 50000, 100000, 500000
# export N=500000
# ./chop_up_genome.py -g ../augustify_modification/genome.fasta.masked -info ./fly_${N}bp_info.txt -headers ./headers_${N}.txt -o ./fly_${N}bp_genome.fasta -l=${N} -n=0
# and then map the genome annotation (annot.gtf) to this set of fragments:
# ./convert_coordinates_in_gtfs.py -headers ./headers_${N}.txt -a ./annot.gtf -o annot_mapped_${N}.gtf -l=${N}


#####___INFO___#####
# fragment lengths: 1500, 3000, 6000, 15000, 50000, 100000, 500000

# nucleotides per line = 60
# fragments = 10                      # how many fragments should one sequence be divided into
# n_overlapped = 30                   # the new sequences should overlap by N nucleotides
# fragment_length = 100               # the lengths of new sequences

__author__ = "Natalia Nenasheva"
__email__ = "nenashen66@uni-greifswald.de"
__status__ = "development"

parser = argparse.ArgumentParser(description='chop_up_sgenome.py allows to create an artificial genome from the input genomic sequence, consisting of short sequences. ' +
                                             'It is easy to specify the length of the fragments into which you want to split (-l/--fragment_length), ' +
                                             'randomly select several fragments for each sequence from the original file (-f/--n_fragments), ' +
                                             'and also select fragments taking into account overlappping (-n/--n_overlapped).'
                                             'In addition to the output_fasta file, a file with general information is also created (specified by -info/--output_info). ')

parser.add_argument('-g', '--genome', required=True, type=str,
                    help='Genome fasta file (possibly softmasked)')
parser.add_argument('-o', '--output_fasta', required=True, type=str,
                    help='Genome fasta file that includes subsequences')
parser.add_argument('-info', '--output_info', required=False, type=str, default='output_info.txt',
                    help='Output file contains additional information about input sequences.')
parser.add_argument('-headers', '--extended_headers', required=False, type=str, default='headers.txt',
                    help='Output file contains extended information about subsequences into which we split the genome.')
parser.add_argument('-l', '--fragment_length', required=True, type=int,
                    help='It is the lengths of new sequences')
parser.add_argument('-n', '--n_overlapped', required=False, type=int, default=0,
                    help='his key specifies that new sequences should overlap by N nucleotides.')
parser.add_argument('-f', '--n_fragments', required=False, type=int,
                    help='This key specifies how many fragments (sub sequences) per sequence have to be save to output.fasta; these fragments will be selected randomly.')
args = parser.parse_args()


''' ******************* BEGIN FUNCTIONS *************************************'''


def get_subsequences(sequence, sub_length, overlapping):
    res = []
    initial_name = sequence.id

    for idx in range(0, len(sequence) - sub_length + 1, sub_length - overlapping):
        sequence.id = initial_name
        name_to_add = 'sub_seq_start=' + str(idx) + ';' + 'sub_seq_end=' + str(idx+sub_length) + \
                      ';' + 'source_name=' + sequence.id + ';' + 'source_len=' + str(len(sequence)) + ';'
        sequence.id = name_to_add
        sequence.description = name_to_add
        res.append(sequence[idx: idx + sub_length])

        if overlapping == 0:                                # case without overalpping: n_overlapped = 0
            last_length = len(sequence) % sub_length        # how many nucleotides do we have in the last fragment?
            #print("Tail length: ", last_length)
            if last_length > 0 and len(sequence) - (idx + sub_length) < sub_length:
                name_to_add = 'sub_seq_start=' + str(len(sequence) - (len(sequence) % sub_length)) + ';' + \
                              'sub_seq_end=' + str(len(sequence)) + ';' + 'source_name=' + initial_name + \
                              ';' + 'source_len=' + str(len(sequence)) + ';'
                sequence.id = name_to_add
                sequence.description = name_to_add
                res.append(sequence[-last_length:])

        else:
            if idx + sub_length == len(sequence):
                break
            elif len(sequence) - (idx + sub_length) + overlapping < sub_length:
                new_ind = idx + sub_length - overlapping
                name_to_add = 'sub_seq_start=' + str(new_ind) + ';' + \
                              'sub_seq_end=' + str(len(sequence)) + ';' + 'source_name=' + initial_name + \
                              ';' + 'source_len=' + str(len(sequence)) + ';'
                sequence.id = name_to_add
                sequence.description = name_to_add
                res.append(sequence[new_ind:])

    print(str(len(res)) + " subsequences for " + initial_name)
    return res


''' ******************* END FUNCTIONS *************************************'''


n_overlapped = args.n_overlapped
fragment_length = args.fragment_length

with open(args.output_info, "w") as f:

    if args.n_fragments:
        n_fragments = args.n_fragments
    else:
        n_fragments = 'all'

    f.write('%s\t%s\t%s' % ('n_overlapped = ' + str(n_overlapped),
                        'fragment_length = ' + str(fragment_length),
                        'n_fragments = ' + str(n_fragments)))
    f.write("\n")
    f.write('========================================================' + "\n")
    f.write('seq_id' + '\t' + 'total_length' + "\t" + 'fragments' + "\n")

    print("Overlapping windows : " + str(n_overlapped))
    print("Fragment lengths : " + str(fragment_length))
    print(f"Try to write {n_fragments} for each sequence to output fasta file")
    print("========================================================")

    with open(args.genome) as inp, open(args.output_fasta, 'w') as out, open(args.extended_headers, 'w') as headers:

        for seq_record in SeqIO.parse(inp, "fasta"):
            if len(seq_record.seq) % fragment_length > 0:
                fragments_total_number = len(seq_record.seq) // fragment_length + 1                                     # added +1, because we have t take into account a tail_fragment too!
            else:
                fragments_total_number = len(seq_record.seq) // fragment_length                                         # length of the tail_fragment is equal to zero
            f.write('%s\t%i\t%i' % (seq_record.id, len(seq_record.seq), fragments_total_number))
            f.write("\n")

            if len(seq_record.seq) <= fragment_length:
                initial_name = seq_record.id
                name_to_add = 'sub_seq_number=1' + ';' + \
                              'sub_seq_start=0' + ';' + \
                              'sub_seq_end=' + str(len(seq_record.seq) - 1) + ';' + \
                              'source_name=' + seq_record.id + ';' + \
                              'source_len=' + str(len(seq_record.seq)) + ';'

                headers.write(name_to_add)
                headers.write('\n')
                name_to_add = initial_name + '_' + '1'

                seq_record.id = name_to_add
                seq_record.description = name_to_add
                SeqIO.write(seq_record, out, 'fasta')
                print("Only one subsequences for " + initial_name)

            else:
                sub_seqs = get_subsequences(sequence=seq_record,
                                            sub_length=fragment_length,
                                            overlapping=n_overlapped)
                for idx, sub_seq_record in enumerate(sub_seqs):

                    name_to_add = 'sub_seq_number=' + str(idx + 1) + ';' + sub_seq_record.id

                    sub_seq_record.id = sub_seq_record.id.split('source_name=')[1].split(';')[0] + '_' + str(idx + 1)
                    sub_seq_record.description = sub_seq_record.id
                    sub_seq_record.name = sub_seq_record.id

                    headers.write(name_to_add)
                    headers.write('\n')


                if n_fragments != "all" and int(n_fragments) > fragments_total_number:                                  # if the user tries to extract more fragments than there actually are, the program will process as many as there are available
                    n_fragments = "all"

                if n_fragments == 'all':
                    SeqIO.write(sub_seqs, out, 'fasta')
                else: # randomly select suq_sequences without repeats
                    sub_seqs = random.choices(list(sub_seqs), k=n_fragments)
                    SeqIO.write(sub_seqs, out, 'fasta')
