#!/bin/bash 
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=300G
#SBATCH -p normal
#SBATCH --qos=default
#SBATCH --job-name=clingo_batterie_tests
#SBATCH --mail-user=pgiraude@genoscope.cns.fr
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --error=log/%x-%J-%u.err
#SBATCH --output=log/%x-%J-%u.out

cd $SLURM_SUBMIT_DIR

module load conda
conda activate clingo

# Boucle sur tous les scripts Python commençant par main_lp
for script in main_lp*.py; do
    echo "Exécution de $script"
    python3 "$script"
done

# Déplacement des fichiers *1000_cluster*.txt
echo "Déplacement des fichiers *1000_cluster*.txt"
mv data/first/*1000_cluster*.txt retour_res/ret_pac/

# clingo <options> <predicate_files> > <output_file>
