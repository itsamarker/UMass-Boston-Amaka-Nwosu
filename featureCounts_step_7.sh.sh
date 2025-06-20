#!/bin/bash
#SBATCH --job-name=featureCounts # you can give your job a name
#SBATCH --nodes 1 # the number of processors or tasks
#SBATCH --cpus-per-task=16
#SBATCH --account=itcga # our account
#SBATCH --time=10:00:00 # the maximum time for the job
#SBATCH --mem=4gb # the amount of RAM
#SBATCH --partition=itcga # the specific server in chimera we are us$
#SBATCH --reservation=2025_JUNE_ITCGA_WORKSHOP #This reservation let$
#SBATCH --error=%x-%A_%a.err   # a filename to save error messages i$
#SBATCH --output=%x-%A_%a.out  # a filename to save any printed outp$

module load subread-2.0.2-gcc-10.2.0

annotation_dir=$1
input_dir=$2
output_dir=$3

for file in "$input_dir"/*_sorted.bam

do

name=$(basename "$file" _sorted.bam)

featureCounts -a "$annotation_dir/Homo_sapiens.GRCh38.111.gtf" \
	-o "$output_dir/${name}_counts.txt" \
	-T 24 \
	-p -B -C "$input_dir/${name}_sorted.bam"

done
