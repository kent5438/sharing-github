#!/usr/bin/env bash
#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#
# =============================================================================
# Setup Instructions
# =============================================================================
#
# 1. Add any other necessary Slurm arguments such as partition (-p) or account
#    (-A). If your system requires a walltime (-t), 24 hours (24:00:00) is
#    sufficient.  We recommend you do not remove any arguments below or Martian
#    may not run properly.
#
# 2. Change filename of slurm.template.example to slurm.template.
#
# =============================================================================
# Template
# =============================================================================
#
#SBATCH -J 10x_supernova
#SBATCH --partition=all
#SBATCH --export=ALL
#SBATCH --nodes=1 --ntasks=4 --cpus-per-task=8
#SBATCH --signal=2
#SBATCH --no-requeue
### Alternatively: --ntasks=1 --cpus-per-task=__MRO_THREADS__
###   Consult with your cluster administrators to find the combination that
###   works best for single-node, multi-threaded applications on your system.
#SBATCH --mem=300G
#SBATCH -o slurm.log
#SBATCH -e slurm.err

/export/EC1680U/kentchen/opt/supernova-2.1.1/supernova run \
        --id "supernova_run" \
        --fastqs="/mnt/NFS/EC2480U-P/kentchen/10x/1811KHX-0004/NTNU-SHLi/reads" \
        --localcores=32 \
        --localmem=300 \
        --maxreads=470000000
