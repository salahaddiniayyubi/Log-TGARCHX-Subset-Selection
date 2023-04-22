#!/bin/bash -l
# L'argument '-l' est indispensable pour bénéficier des directives de votre .bashrc

# On peut éventuellement placer ici les commentaires SBATCH permettant de définir les paramètres par défaut de lancement :
#SBATCH --cpus-per-task=40
#SBATCH --partition=longrun
#SBATCH --mail-type=FAIL,END


# Activation de l'environnement anaconda python37
cd /share/home/orujov/Log-TGARCHX-Subset-Selection
export R_LIBS=""/share/castor/home/orujov/R/x86_64-pc-linux-gnu-library/4.2""
# Exécution du script habituellement utilisé, on utilise la variable CUDA_VISIBLE_DEVICES qui contient la liste des GPU logiques actuellement réservés (toujours à partir de 0)
# cd samir
Rscript main.R
