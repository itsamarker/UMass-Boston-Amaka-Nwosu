#!/bin/bash
#SBATCH --job-name=trim2 # you can give your job a name
#SBATCH --nodes 1 # the number of processors or tasks
#SBATCH --cpus-per-task=16
#SBATCH --account=itcga # our account
#SBATCH --time=10:00:00 # the maximum time for the job
#SBATCH --mem=4gb # the amount of RAM
#SBATCH --partition=itcga # the specific server in chimera we are using
#SBATCH --reservation=2025_JUNE_ITCGA_WORKSHOP #This reservation lets us cut in line to use itcga cores
#SBATCH --error=%x-%A_%a.err   # a filename to save error messages into
#SBATCH --output=%x-%A_%a.out  # a filename to save any printed output into

# Module load
module load py-dnaio-0.4.2-gcc-10.2.0-gaqzhv4
module load py-xopen-1.1.0-gcc-10.2.0-5kpnvqq
module load py-cutadapt-2.10-gcc-10.2.0-2x2ytr5

# Define variables
input_dir=$1 # takes this from the command line, first item after the script
output_dir=$2 # takes this from the command line, second item

for file in ${input_dir}/*_1.fastq

do

# Pull basename \
name=$(basename ${file} _1.fastq)

# Run cutadapt
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim-n -m 25 \
-o ${output_dir}/${name}_1_trim.fastq \
-p ${output_dir}/${name}_2_trim.fastq \
${input_dir}/${name}_1.fastq \
${input_dir}/${name}_2.fastq

echo cutadapt is finished with $name

done
