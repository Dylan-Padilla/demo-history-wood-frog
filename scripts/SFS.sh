#!/bin/bash

#SBATCH --partition=day                    # which partition to run the job, options are in the Amarel guide
#SBATCH --job-name=aRanSyl_SFS                       # job name for listing in queue
#SBATCH --output=/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/results/logs/SFS/slurm-%j-%x.out
#SBATCH --error=/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/results/logs/SFS/slurm-%j-%x.err
#SBATCH --mem=100G                               # memory to allocate in Mb
#SBATCH -n 10                                   # number of cores to use
#SBATCH --time=1-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --mail-user=dylan.padilla@yale.edu           # email address to send status updates
#SBATCH --mail-type=FAIL,END                     # email for the following reasons

module --force purge
conda activate easySFS

## we will work with a subset of sites that are already filtered for low missing data proportions:

VCF="/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/data/aRanSyl_demographic_0.5pruned.vcf.gz"
POPS="/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/data/pops_demographic_file.txt"

## setting directory

cd /Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/SFS/aRanSyl_SFS/

## estimate projections

/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/SFS/easySFS/easySFS.py -i $VCF -p $POPS --preview -a

## running easySFS.py based on the projections recommended

/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/SFS/easySFS/easySFS.py -i $VCF -p $POPS -a --proj 50,32


