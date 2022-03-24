#!/bin/bash
#SBATCH -c 1
#SBATCH -t 7-00:00
#SBATCH -p holy-smokes
#SBATCH --mem=5000
#SBATCH -o traitR_%j.out
#SBATCH -e traitR_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=subirshakya@fas.harvard.edu

cd /n/holyscratch01/informatics/sshakya/avonet

module load R/4.1.0-fasrc01
export R_LIBS_USER=$HOME/apps/R_4-1-0:$R_LIBS_USER

Rscript /n/holyscratch01/informatics/sshakya/avonet/scripts/ou_plot_cmd.R $1
