#''''''''''''''''''''''''''''''''''''
#' Read Count Analysis with edgeR
#' @author Cooper Kimball-Rhines
#' @date 2025-06-10
#''''''''''''''''''''''''''''''''''''

# Install package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
install.packages("tidyverse")
install.packages("ggrepel")


# Load library
library(edgeR)
library(tidyverse)
library(ggrepel)

# Set working directory manually: Session -> Set Working Directory -> To Source File Location
# Load data from counts.txt file
counts <- read.table("../Downloads/all_counts.txt", header=TRUE, sep="\t", fill=TRUE) %>%
  select(,1:10)
  # rename_with(~str_remove(., '/itcgastorage/share_home/c.kimballrhines001/test_project/results/bam/')) %>%
  # rename_with(~str_remove(., '.bam')) |>
  
  #NOTE: We do not have any replicates in the data we've been working with, so I'm artificially creating
  #replicates based on the data we do have. DO NOT do this with real data or in your projects.
  #mutate(C2 = C1_S4_L001+3,
       #  C3 = C1_S4_L001+5,
       #  T2 = T1_S7_L001+3,
        # T3 = T1_S7_L001+5,
         #V2 = V1_S1_L001+3,
         #V3 = V1_S1_L001+5)

head(counts)

# Provide treatment metadata
condition <- c("MNPs", "MNPs", "MNPs", "MNPs-AMF", "MNPs-AMF", "MNPs-AMF", "PBS", "PBS", "PBS")
dim(counts)

# Filter out genes expressed below cutoff
# Convert column to numeric

counts[ , 2:10] <- lapply(counts[ , 2:10], as.numeric)

totalexp <- rowSums(counts[,2:10])
hist(totalexp)

counts <- filter(counts, totalexp > 10)

# Move annotation info to separate object
ann <- counts[,2:10]

# Combine gene counts and metadata into edgeR object
d <- DGEList(counts = counts[,2:10],
             group = factor(condition),
             genes = counts[,1:9])
str(d)
dim(d) # 12,558 genes across 9 samples

# Normalize the data
d <- calcNormFactors(d, method = "TMM")
d$samples

# Plot MDS
samples <- c("MNPs1", "MNPs-AMF1", "PBS1", "MNPs2", "MNPs_AMF2", "PBS2", "MNPs3", "MNPs_AMF3", "PBS3")

plotMDS(d, col = as.numeric(d$samples$group), label = d$samples$group)

# Generate design matrix
design <- model.matrix(~ 0 + d$samples$group)
design
colnames(design) <- levels(d$samples$group)

# Estimate dispersion
dg <- estimateGLMCommonDisp(d, design)
plotBCV(dg)

# Fit the model
fit <- glmFit(dg, design)

# This compares group 1 (Cs, 1) to group 2 (Ts, -1), ignoring group 3 (Vs, 0)
# and does a likelihood ratio test for each gene
fitCT <- glmLRT(fit, contrast=c(1, -1, 0))

# Sort out the differentially expressed genes
tabCT <- topTags(fitCT,n=Inf,adjust.method="BH", sort.by = "PValue")$table

# Make a significance column
tabCT <- tabCT %>% 
  mutate(significance = case_when((FDR < 0.05 & logFC > 1) ~ "Upregulated", (FDR < 0.05 & logFC < -1) ~ "Downregulated", .default = "Not significant"))

# Pull out the top differentially expressed genes
top_genes <- filter(tabCT, -log10(PValue) > 5)

# Visualize the expression data!
volcano_plot <- ggplot(tabCT, aes(x = logFC, y = -log10(PValue), color = significance)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Not significant" = "grey","Upregulated" = "red", "Downregulated" = "blue")) +
  geom_text_repel(data = top_genes, aes(label = Geneid), size = 3.5, fontface = 'bold') +
  geom_hline(yintercept = 5, linetype = "dashed") +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  theme_classic() +
  theme(legend.position = 'None',
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank()) +
  labs(x = "Log2 Fold Change", y = "-Log10 p-value")

volcano_plot

# Write out the most significantly differential genes
tabCT |>
  filter(significance != "Not significant") |>
  write_tsv(file = "topGenes.tsv")

# Let's investigate the top differentially expressed gene
topCount <- counts |>
  filter(Geneid == "ENSG00000275757") |>
  #select(-c(Geneid, Chr, Start, End, Strand,Length)) |>
  t()

topGene <- data.frame(cbind(as.character(d$samples$group), topCount))
colnames(topGene) <- c("Treatment","Reads")
topGene$Reads <- as.numeric(topGene$Reads)

genePlot <- ggplot(data = topGene, 
                   mapping = aes(x = Treatment, y = Reads, color = Treatment)) +
  geom_jitter(width = 0.25) +
  labs(y = "Read counts at HSPA6") +
  theme_bw()

genePlot # We can see that the counts of this gene are way higher in the TGF-b treatment
