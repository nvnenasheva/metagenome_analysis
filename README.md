# metagenome_analysis
The repository contains additional scripts for the analysis of metagenomic data.

## chop_up_genome.py

This program splits the input into a set of subsequences. 
There are several additional options:
1) set the length of the subsequences into which the genome should be divided (-l option)
2) splitting can be done taking into account the overlapping of subsequences on top of each other (-n option)
3) you can choose a randomly selected pool of subsequences (-f option is the number of subsequences for each sequence)

## get_genome_fraction.py
The goal of this script is to reduce the size of input genome => it helps to take only 25% (-f=0.25) of total amount of data randomly.

## convert_coordinates_in_gtfs.py
The main goal of this script is to map global coordinates from annot.gtf file to coordinates of genome sequence fragments.

## prepare_contigs.py
This script groups data from fasta and gtf files according to the species type assigned by the Augustify.

## slurm_braker.sh
It is the way to run BRAKER pipeline using all sequences that were pulled together (for each each parameter set) on the previous step.
Variable SPECIES has parameter sets from the _prepare_contigs.py_ output. 

* These three scripts could be used together:
```
export N=500000
# chop up the initial genome:
./chop_up_genome.py -g genome.fasta.masked -info fly_${N}bp_info.txt -o fly_${N}bp_genome.fasta -l=${N} -n=0 -headers headers_${N}bp.txt

# map the genome annotation (annot.gtf) to this set of fragments:
 ./convert_coordinates_in_gtfs.py -headers headers_${N}bp.txt -a annot.gtf -o annot_mapped_${N}.gtf -l=${N}
 
# extract 25% sequences randomly:
./get_genome_fraction.py -g fly_${N}bp_genome.fasta -o fly_${N}bp_genome_reduced.fasta -f=0.25 -headers=random_headers_${N}.txt

# intermediate result: annot_mapped_${N}.gtf and ffly_${N}bp_genome_reduced.fasta

# divide data according to the scpecies assigned by Augustify - useful to run braker pipeline
./prepare_contigs.py ...

# run BRAKER
sbatch ./slurm_braker.sh

```
