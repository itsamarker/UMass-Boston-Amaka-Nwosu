#!/bin/bash
#SBATCH --job-name=trim # you can give your job a name
#SBATCH --nodes 1 # the number of processors or tasks
#SBATCH --cpus-per-task=2
#SBATCH --account=itcga # our account
#SBATCH --time=10:00:00 # the maximum time for the job
#SBATCH --mem=4gb # the amount of RAM
#SBATCH --partition=itcga # the specific server in chimera we are using
#SBATCH --reservation=2025_JUNE_ITCGA_WORKSHOP #This reservation lets us cut in line to use itcga cores
#SBATCH --error=%x-%A_%a.err   # a filename to save error messages into
#SBATCH --output=%x-%A_%a.out  # a filename to save any printed output into

# Module load
module load fastqc-0.11.9-gcc-10.2.0-osi6pqc

fastqc -o ~/project/results/untrimmed_fastqc/  ~/project/data/data/*.fastq
