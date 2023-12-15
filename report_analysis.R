# load-libraries ----------------------------------------------------------

library(tidyverse)
#install.packages("ggrepel")
library(ggrepel)
#install.packages("BiocManager")
library(BiocManager)

#BiocManager::install("tximport")
library(tximport)
#BiocManager::install("rhdf5")
library(rhdf5)

#install.packages("remotes")
library(remotes)
#remotes::install_github("pachterlab/sleuth")
#remotes::install_github("kevinblighe/EnhancedVolcano")
#BiocManager::install("fgsea")
library(ggrepel)
library(EnhancedVolcano)
library(fgsea)
#install.packages("devtools")
library(devtools)

# load-tables ----------------------------------------------------------
# Read outputs from Linux in R
allTPMs <- read.table("allTPMs.tsv", header=TRUE, sep="\t")
hypoxia_results <- read.table("hypoxiaDEA_results.tsv", header=TRUE, sep="\t")
res <- read.table("TPMs_and_DEA_results_file.tsv", header=TRUE, sep="\t")
genesets <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")

# Find the number of genes that are upregulated and downregulated over 150%
# with significant q-values < 0.001
sum(res$log2FC>=1.5 & res$qval<0.001)
#[1] 30
sum(res$log2FC<=-1.5 & res$qval<0.001)
#[1] 31

# Top 25 genes by product log2FC * -1log10(qval)
siggenes <- (res %>% arrange(desc(abs(pi))) %>% slice(1:50))$genes

# Genes in each geneset from siggenes
hypoxia_genes <- genesets[["HALLMARK_HYPOXIA"]]
sig_hypoxia_genes <- intersect(hypoxia_genes, siggenes)

cholesterol_genes <- genesets[["HALLMARK_CHOLESTEROL_HOMEOSTASIS"]]
sig_chol_genes <- intersect(siggenes, cholesterol_genes)

p53_genes <- genesets[["HALLMARK_P53_PATHWAY"]]
sig_p53_genes <- intersect(siggenes, p53_genes)

glycolysis_genes <- genesets[["HALLMARK_GLYCOLYSIS"]]
sig_glycolysis_genes <- intersect(siggenes, glycolysis_genes)


# HIF-1a target genes - 

hif_targets <- read_file("hif_targets.txt") %>% 
  str_split("\n") %>% 
  unlist() %>% 
  str_split(",") %>% 
  unlist() %>% 
  setdiff("") %>% 
  str_extract("[A-Z0-9]+") %>% 
  unique() %>%
  setdiff(c("V", "1", "CHEA", "HIF1A")) 

# plot tpms as bar chart for hif-1a target genes
hif1a_tpm <- res %>% filter(genes %in% hif_targets) %>%
  select(genes, normoxia_avg, hypoxia_avg, qval, log2FC) %>% 
  pivot_longer(cols=c(normoxia_avg, hypoxia_avg), names_to="condition", values_to="tpm") %>% 
  mutate(condition=ifelse(condition=="normoxia_avg", "Normoxia", "Hypoxia"))
  
hif1a_tpm %>% 
  filter(qval <= 0.001) %>% 
  ggplot() + 
  geom_bar(aes(x=genes, y=tpm, fill=condition), 
           stat="identity", 
           position="dodge") +
  xlab("Gene") +
  ylab("TPM") +
  labs(fill="Condition") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=c("red", "grey"))

# plot tpms as bar chart for urothelial differentiation markers
urothelial_markers <- c("GATA3", "UPK3A", "THBD", "TP63", "KRT5", "S100P", "UPK2")

uro_tpm <- res %>% filter(genes %in% urothelial_markers) %>% 
  select(genes, normoxia_avg, hypoxia_avg, qval) %>% 
  pivot_longer(cols=c(normoxia_avg, hypoxia_avg), names_to="condition", values_to="tpm") %>% 
  mutate(condition=ifelse(condition=="normoxia_avg", "Normoxia", "Hypoxia"))


uro_tpm %>% 
  filter(qval <= 0.01) %>% 
  ggplot() + 
  geom_bar(aes(x=genes, y=tpm, fill=condition), 
           stat="identity", 
           position="dodge") +
  xlab("Gene") +
  ylab("TPM") +
  labs(fill="Condition") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=c("red", "grey"))

# Volcanco plot
EnhancedVolcano(res, 
                x = "log2FC", 
                y = "qval", 
                title = "Urothelial Cells in Hypoxia",
                subtitle = "Differential expression analysis",
                FCcutoff = 1.5, 
                pCutoff = 0.001,
                lab = res$genes, 
                selectLab = siggenes,  
                labSize = 4.0, 
                max.overlaps = 1000, 
                pointSize = 2.0,
                legendPosition = 0,
                drawConnectors = TRUE,
                ylim = c(0, (max(-1*log10(res$qval)) + 0.5)), 
                xlim = c((max(abs(res$log2FC))*-1), max(abs(res$log2FC))+1),
                shape = c(1, 4, 23, 25),
                colAlpha = 1,
                hline = c(0.05),
                gridlines.major = FALSE, 
                gridlines.minor = FALSE)

# Gene set enrichment analysis (GSEA)
prerank <- res[c("genes", "pi")]
prerank <- setNames(prerank$pi, prerank$genes)

fgseaRes <- fgsea(pathways = genesets, stats = prerank, minSize=100, maxSize=500)

# Plot top enrichment pathways
topPathwaysUp <- fgseaRes[ES > 0.5][order(pval)][padj < 0.001][, pathway]
topPathwaysDown <- fgseaRes[ES < -0.5][order(pval)][padj < 0.001][, pathway]
  
topPathways <- c(topPathwaysUp, topPathwaysDown)

plotGseaTable(genesets[topPathways], prerank, fgseaRes, 
              gseaParam = 0.5) 

## Hypoxia
plotEnrichment(genesets[["HALLMARK_HYPOXIA"]], prerank) + 
  labs(title="Hypoxia pathway: Genes up-regulated in response to low oxygen levels (hypoxia)") +
  theme_classic() +
  xlab("Gene Rank") + 
  ylab("Enrichment score (ES)")

## Glycolysis
plotEnrichment(genesets[["HALLMARK_GLYCOLYSIS"]], prerank) + labs(title="Glycolysis pathway: Genes encoding proteins involved in glycolysis and gluconeogenesis") +
  theme_classic() +
  xlab("Gene Rank") + 
  ylab("Enrichment score (ES)")

## MTORC1
plotEnrichment(genesets[["HALLMARK_MTORC1_SIGNALING"]], prerank) + labs(title="MTORC1 Signalling Pathway: Genes up-regulated through activation of mTORC1 complex") +
  theme_classic() +
  xlab("Gene Rank") + 
  ylab("Enrichment score (ES)")

## Unfolded protein
plotEnrichment(genesets[["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]], prerank) + labs(title="Unfolded Protein Response: Genes up-regulated during unfolded protein response") +
  theme_classic() +
  xlab("Gene Rank") + 
  ylab("Enrichment score (ES)")

## Interferon gamma response
plotEnrichment(genesets[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]], prerank) + labs(title="Interferon Gamma Response: Genes up-regulated in response to interferon gamma") +
  theme_classic() +
  xlab("Gene Rank") + 
  ylab("Enrichment score (ES)")