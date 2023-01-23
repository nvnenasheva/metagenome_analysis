# metagenome_analysis
The repository contains additional scripts for the analysis of metagenomic data.

## Chop_up_genome.py script

This program splits the input into a set of subsequences. 
There are several additional options:
1) set the length of the subsequences into which the genome should be divided (-l option)
2) splitting can be done taking into account the overlapping of subsequences on top of each other (-n option)
3) you can choose a randomly selected pool of subsequences (-f option is the number of subsequences for each sequence)
