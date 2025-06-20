#!/bin/bash
#SBATCH --job-name=sort # you can give your job a name
#SBATCH --nodes 1 # the number of processors or tasks
#SBATCH --cpus-per-task=16
#SBATCH --account=itcga # our account
#SBATCH --time=10:00:00 # the maximum time for the job
#SBATCH --mem=4gb # the amount of RAM
#SBATCH --partition=itcga # the specific server in chimera we are us$
#SBATCH --reservation=2025_JUNE_ITCGA_WORKSHOP #This reservation let$
#SBATCH --error=%x-%A_%a.err   # a filename to save error messages i$
#SBATCH --output=%x-%A_%a.out  # a filename to save any printed outp$

module load samtools-1.10-gcc-9.3.0-flukja5

input_dir=$1


for file in "$input_dir"/*.bam
do

name=$(basename "$file" .bam)

samtools sort "$file" -@ 4 -o "$input_dir/${name}_sorted.bam"
done
