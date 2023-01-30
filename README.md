# metagenome_analysis
The repository contains additional scripts for the analysis of metagenomic data.

## chop_up_genome.py

This program splits the input into a set of subsequences. 
There are several additional options:
1) set the length of the subsequences into which the genome should be divided (-l option)
2) splitting can be done taking into account the overlapping of subsequences on top of each other (-n option)
3) you can choose a randomly selected pool of subsequences (-f option is the number of subsequences for each sequence)


## convert_coordinates_in_gtfs.py
The main goal of this script is to map global coordinates from annot.gtf file to coordinates of genome sequence fragments.

These two scripts coulb be used together:
```
export N=500000
chop_up_genome.py -g genome.fasta.masked -info fly_${N}bp.txt -o fly_${N}bp_genome.fasta -l=${N} -n=0

# and then mapp the genome annotation (annot.gtf) to this set of fragments:
convert_coordinates_in_gtfs.py -g fly_${N}bp_genome.fasta -a annot.gtf -o annot_mapped.gtf

# our results: annot_mapped.gtf and fly_${N}bp_genome.fasta
```
