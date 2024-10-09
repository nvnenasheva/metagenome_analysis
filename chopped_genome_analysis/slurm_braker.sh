#!/bin/bash
#SBATCH -o br.%j.%N.out
#SBATCH -e br.%j.%N.err
#SBATCH -J braker
#SBATCH --get-user-env 
#SBATCH --time=72:00:00
#SBATCH -N 1 # number of nodes
#SBATCH -n 48
#SBATCH --mem=96000
#SBATCH -p pinky


export T=48
export SPECIES=$(echo 'E_coli_K12' 'aedes' 'amphimedon' 'arabidopsis' 'aspergillus_fumigatus' 'aspergillus_nidulans' 'aspergillus_oryzae' 'aspergillus_terreus' 'botrytis_cinerea' 'brugia' 'caenorhabditis' 'candida_albicans' 'candida_guilliermondii' 'candida_tropicalis' 'chaetomium_globosum' 'chicken' 'chlamydomonas' 'coccidioides_immitis' 'coprinus' 'coprinus_cinereus' 'coyote_tobacco' 'cryptococcus_neoformans_gattii' 'cryptococcus_neoformans_neoformans_B' 'cryptococcus_neoformans_neoformans_JEC21' 'debaryomyces_hansenii' 'encephalitozoon_cuniculi_GB' 'eremothecium_gossypii' 'fly' 'fusarium_graminearum' 'galdieria' 'histoplasma_capsulatum' 'kluyveromyces_lactis' 'laccaria_bicolor' 'leishmania_tarentolae' 'lodderomyces_elongisporus' 'magnaporthe_grisea' 'maize' 'nasonia' 'neurospora_crassa' 'phanerochaete_chrysosporium' 'pichia_stipitis' 'pneumocystis' 'rhizopus_oryzae' 's_aureus' 'saccharomyces_cerevisiae_S288C' 'saccharomyces_cerevisiae_rm11-1a_1' 'schistosoma' 'schizosaccharomyces_pombe' 'tetrahymena' 'thermoanaerobacter_tengcongensis' 'tomato' 'toxoplasma' 'trichinella' 'ustilago_maydis' 'volvox' 'wheat' 'yarrowia_lipolytica' 'zebrafish')
#PATH=~/BRAKER/scripts/braker.pl/:$PATH
#export PATH

cd ~/chopped_up_genomes/
if [ ! -d braker_fly_1500_output ]; then
  mkdir -p braker_fly_1500_output
fi
cd braker_fly_1500_output


export DIR=~/chopped_up_genomes/data_contig_based
for S in $SPECIES; do

  if [ -s ${DIR}/tmp_${S}.fasta ]; then
    echo Iteration for $S was skipped.
    continue
fi
  
  if [ ! -d $S ]; then
    echo Directory for $S was made.
  mkdir -p $S
fi
  
  cd $S
  braker.pl --species=${S} --skipAllTraining --esmode --softmasking --genome ${DIR}/tmp_${S}.fasta --annot ${DIR}/tmp_annot_${S}.gtf --eval_pseudo=${DIR}/tmp_pseudo_${S}.gtf --cores=${T} &> braker_${S}.out
  echo Done for ${S}.

done
