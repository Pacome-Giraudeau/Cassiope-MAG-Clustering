#!/bin/bash 
#SBATCH -t 24:00:00 # temps demandé pour le travail (format HH:MM:SS ou J-HH:MM:SS)
#SBATCH -N 1 # nombre de noeuds demandés
#SBATCH -n 1 # nombre de tâches à exécuter
#SBATCH -c 32 # nombre de coeurs (cpus) 
#SBATCH --mem=300G # mémoire demandée (default: M)
#SBATCH -p normal # nom de la partition ("normal", "xlarge" ou "xxlarge")
#SBATCH --qos=default # nom de la queue de travail. pour normal, choisir default ou long, pour (x)xlarge, (x)xlarge ou (x)xlarge_month -sans parenthèse...)
#SBATCH --job-name=clingo_batterie_tests # nom du job
#SBATCH --mail-user=mlhermit@genoscope.cns.fr # email utilisateur: à remplir
#SBATCH --mail-type=BEGIN,FAIL,END #quand envoyer un mail (BEGIN, FAIL, END, ALL)
#SBATCH --error=log/%x-%J-%u.err #chemin du fichier où écrire l'erreur standard
#SBATCH --output=log/%x-%J-%u.out #chemin du fichier où écrire la sortie standard (facultatif, surtout si vous redirigez vous-même la sortie de clingo dans un fichier de sortie)


cd $SLURM_SUBMIT_DIR

module load conda #charge le module conda
conda activate clingo #active l'environnement clingo (modifiez pour le nom de votre environnement le cas échéant)

# commande clingo ou python: a vous de compléter !
python3 main_to_lp_copy.py

# clingo <options> <predicate_files> > <output_file>
