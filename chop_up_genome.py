#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import random


#chop_up_genome.py -g ./example.fasta -o ./output_example.fasta -l 100

# export N=500000
# ./chop_up_genome.py -g ../augustify_modification/genome.fasta.masked -info ./fly_${N}bp.txt -o ./fly_${N}bp_genome.fasta -l=${N} -n=0


# TEST
# 1500, 3000, 6000, 15000, 50000, 100000, 500000

# nucleotides per line = 60
# fragments = 10                      # how many fragments should one sequence be divided into
# n_overlapped = 30                   # the new sequences should overlap by N nucleotides
# fragment_length = 100               # the lengths of new sequences

parser = argparse.ArgumentParser()

parser.add_argument('-g', '--genome', required=True, type=str,
                    help='Genome fasta file (possibly softmasked)')
parser.add_argument('-o', '--output_fasta', required=True, type=str,
                    help='Genome fasta file that includes subsequences')
parser.add_argument('-info', '--output_info', required=False, type=str, default='output_info.txt',
                    help='Output file contains additional information about input sequences.')
parser.add_argument('-l', '--fragment_length', required=True, type=int,
                    help='the lengths of new sequences')
parser.add_argument('-n', '--n_overlapped', required=False, type=int, default=0,
                    help='new sequences should overlap by N nucleotides')
parser.add_argument('-f', '--n_fragments', required=False, type=int,
                    help='how many fragments (sub sequences) per sequence have to be save to output.fasta; these fragments will be selected randomly.')


args = parser.parse_args()

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
        #print("start_ind :", idx, "; end_idx :", idx + sub_length, " for ", initial_name)

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
    #print("__________________________________________________")
    return res


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

    with open(args.genome) as inp, open(args.output_fasta, 'w') as out:

        for seq_record in SeqIO.parse(args.genome, "fasta"):
            if len(seq_record.seq) % fragment_length > 0:
                fragments_total_number = len(seq_record.seq) // fragment_length + 1 # add +1, because we have t take into account a tail_fragment too!
            else:
                fragments_total_number = len(seq_record.seq) // fragment_length        # length of the tail_fragment is equal to zero
            f.write('%s\t%i\t%i' % (seq_record.id, len(seq_record.seq), fragments_total_number))
            f.write("\n")

            if len(seq_record.seq) <= fragment_length:
                initial_name = seq_record.id
                name_to_add = 'sub_seq_number=1' + ';' + \
                              'sub_seq_start=0' + ';' + \
                              'sub_seq_end=' + str(len(seq_record.seq) - 1) + ';' + \
                              'source_name=' + seq_record.id + ';' + \
                              'source_len=' + str(len(seq_record.seq)) + ';'
                seq_record.id = name_to_add
                seq_record.description = name_to_add
                SeqIO.write(seq_record, out, 'fasta')
                print("Only one subsequences for " + initial_name)

            else:
                sub_seqs = get_subsequences(sequence=seq_record,
                                            sub_length=fragment_length,
                                            overlapping=n_overlapped)
                for idx, sub_seq_record in enumerate(sub_seqs):
                    sub_seq_record.id = 'sub_seq_number=' + str(idx + 1) + ';' + sub_seq_record.id
                    sub_seq_record.description = sub_seq_record.id
                    sub_seq_record.name = sub_seq_record.id

                if n_fragments != "all" and int(n_fragments) > fragments_total_number:  # if the user tries to extract more fragments than there actually are, the program will process as many as there are available
                    n_fragments = "all"

                if n_fragments == 'all':
                    SeqIO.write(sub_seqs, out, 'fasta')
                else: # randomly select suq_sequences without repeats
                    sub_seqs = random.choices(list(sub_seqs), k=n_fragments)
                    SeqIO.write(sub_seqs, out, 'fasta')

########################################################################################################################
# 1. add function to transform/map global coordinates of annot.gtf & pseudo.gff3 to fragment coordinates (see drawing) - DONE (I decided to put it in a separate script)
# 2. catch when you try to generate more fragments than sanely possible (if we want to get 10 fragments, but in fact we just have 6, for example) - DONE
