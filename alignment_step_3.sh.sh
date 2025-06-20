#!/bin/bash
#SBATCH --job-name=alignment # you can give your job a name
#SBATCH --nodes 1 # the number of processors or tasks
#SBATCH --cpus-per-task=16
#SBATCH --account=itcga # our account
#SBATCH --time=10:00:00 # the maximum time for the job
#SBATCH --mem=4gb # the amount of RAM
#SBATCH --partition=itcga # the specific server in chimera we are us$
#SBATCH --reservation=2025_JUNE_ITCGA_WORKSHOP #This reservation let$
#SBATCH --error=%x-%A_%a.err   # a filename to save error messages i$
#SBATCH --output=%x-%A_%a.out  # a filename to save any printed outp$

# Module load
module load gcc-10.2.0-gcc-9.3.0-f3oaqv7
module load python-3.8.12-gcc-10.2.0-oe4tgov
module load hisat2-2.1.0-gcc-9.3.0-u7zbyow

# Define variables
index_dir=$1
input_dir=$2 # takes this from the command line, first item after th$
output_dir=$3 # takes this from the command line, second item

for file in ${input_dir}/*_1_trim.fastq

do

# Pull basename \
name=$(basename ${file} _1_trim.fastq)

hisat2 -p 24 \
	-x "$index_dir" \
	-1 "$input_dir/${name}_1_trim.fastq" \
	-2 "$input_dir/${name}_2_trim.fastq" \
	-S "$output_dir/${name}.sam"

echo alignment is finished with $name

done
