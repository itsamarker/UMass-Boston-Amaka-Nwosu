## Magnetic Hyperthermia (MHT) in Ovarian Cancer: ITCGA Summer Workshop Project

Hello! This project is part of the **Summer ITCGA Workshop at UMass Boston**, in collaboration with **Dana-Farber: Integrated Training in Computational Genomics and Data Sciences**.

###  Project Overview
- Magnetic hyperthermia (MHT) is a thermal cancer treatment using **magnetic nanoparticles (MNPs)**.
- When exposed to an **alternating magnetic field (AMF)**, MNPs generate heat to:
  - Destroy cancer cells via apoptosis (programmed cell death)
  - Minimize damage to surrounding healthy tissue
- Effectiveness depends on:
  - High accumulation of MNPs within tumors
  - Tumors’ poor heat dissipation due to abnormal vasculature
- This study explores the genomic effects of MHT in **ovarian cancer cells**.
  - Focused on how the intracellular microenvironment impacts MHT outcomes
  - Cancer cells' disorganized cytoskeleton enhances MNP heat generation compared to normal cells

### RNA-seq Analysis Pipeline on Chimera Cluster

- **Download FastQ Files**
  - Used sratoolkit and a custom fastq-dump.sh script to download sequencing data from SRA

- **Quality Control**
  - Assessed read quality using FastQC

- **Adapter Trimming**
  - Removed Illumina adapters and low-quality bases using cutadapt
  - Employed SLURM scripts for batch trimming

- **Data Alignment & Expression Analysis**
  -Next steps include:
    - Aligning reads to a reference genome
    - Quantifying gene expression
    - Identifying differentially expressed genes between treatment groups
