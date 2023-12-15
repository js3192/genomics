### Load libraries ###

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
#BiocManager::install("scatterpie") 
#BiocManager::install("enrichplot")
library(enrichplot)

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


# HIF-1a target genes
#EPO, VEGF, HO-1, ADM, and Glut-1
hif1a_targets <- c("EPO", "VEGFA", "HMOX1", "ADM", "SLC2A1", "B2M", "HPRT")
hypoxia_pub_genes <- c("AKAP12", "ALDOB", "CASP6", "DTNA", "HS3ST1", "JUN", "KDELR3", "STC1")

# read hif_targets.txt
hif_targets <- read_file("hif_targets.txt") %>% 
  str_split("\n") %>% 
  unlist() %>% 
  str_split(",") %>% 
  unlist() %>% 
  setdiff("") %>% 
  # extract the gene names
  str_extract("[A-Z0-9]+") %>% 
  # remove "1" from list
  # remove duplicates
  unique() %>%
  setdiff("")

# plot tpms as bar chart for hif-1a target genes
hif1a_tpm <- res %>% filter(genes %in% sig_hypoxia_genes) %>%
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
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values=c("red", "grey"))

# plot tpms as bar chart for urothelial differentiation markers
urothelial_markers <- c("KRT5", "KRT14", "KRT20", "UPK1B", "UPK2", "UPK3A", "UPK3B", "UPK3BL", "UPK3BL1", "UPK3BL2")

uro_tpm <- res %>% filter(genes %in% urothelial_markers) %>% 
  select(genes, normoxia_avg, hypoxia_avg) %>% 
  pivot_longer(cols=c(normoxia_avg, hypoxia_avg), names_to="condition", values_to="tpm") %>% 
  mutate(condition=ifelse(condition=="normoxia_avg", "Normoxia", "Hypoxia"))


uro_tpm %>% 
  ggplot() + 
  geom_bar(aes(x=genes, y=tpm, fill=condition), 
           stat="identity", 
           position="dodge") +
  xlab("Gene") +
  ylab("TPM") +
  labs(fill="Condition") +
  theme_classic() 

# cell junction genes
cell_junction_genes <- c("")

# cell junction genes in urothelial cells
uro_cell_junction_genes <- res %>% filter(genes %in% cell_junction_genes) %>% 
  select(genes, normoxia_avg, hypoxia_avg) %>% 
  pivot_longer(cols=c(normoxia_avg, hypoxia_avg), names_to="condition", values_to="tpm") %>% 
  mutate(condition=ifelse(condition=="normoxia_avg", "Normoxia", "Hypoxia"))

uro_cell_junction_genes %>% 
  ggplot() + 
  geom_bar(aes(x=genes, y=tpm, fill=condition), 
           stat="identity", 
           position="dodge") +
  xlab("Gene") +
  ylab("TPM") +
  labs(fill="Condition") +
  theme_classic()

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

## Hypoxia
plotEnrichment(genesets[["HALLMARK_HYPOXIA"]], prerank) + 
  labs(title="Hypoxia pathway: Genes up-regulated in response to low oxygen levels (hypoxia)") +
  theme_classic() +
  xlab("Gene Rank") + 
  ylab("Enrichment score (ES)")
ggsave("hypoxiatargets.pdf")

## Glycolysis
plotEnrichment(genesets[["HALLMARK_GLYCOLYSIS"]], prerank) + labs(title="Glycolysis pathway: Genes encoding proteins involved in glycolysis and gluconeogenesis") +
  theme_classic() +
  xlab("Gene Rank") + 
  ylab("Enrichment score (ES)")
ggsave("glycolysistargets.pdf")

## MTORC1
plotEnrichment(genesets[["HALLMARK_MTORC1_SIGNALING"]], prerank) + labs(title="MTORC1 Signalling Pathway: Genes up-regulated through activation of mTORC1 complex") +
  theme_classic() +
  xlab("Gene Rank") + 
  ylab("Enrichment score (ES)")
ggsave("mtorc1targets.pdf")

## Unfolded protein
plotEnrichment(genesets[["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]], prerank) + labs(title="Unfolded Protein Response: Genes up-regulated during unfolded protein response, a cellular stress response related to the endoplasmic reticulum.") +
  theme_classic() +
  xlab("Gene Rank") + 
  ylab("Enrichment score (ES)")
ggsave("unfoldedptargets.pdf")

## P53
plotEnrichment(genesets[["HALLMARK_P53_PATHWAY"]], prerank) + labs(title="P53 Signalling Pathway: Genes involved in p53 pathways and networks") +
  theme_classic() +
  xlab("Gene Rank") + 
  ylab("Enrichment score (ES)")
ggsave("P53targets.pdf")

## DNA Repair 
plotEnrichment(genesets[["HALLMARK_DNA_REPAIR"]], prerank) + labs(title="DNA repair targets")
ggsave("DNArepair.pdf")

## Cholesterol homeostasis
plotEnrichment(genesets[["HALLMARK_CHOLESTEROL_HOMEOSTASIS"]], prerank) + labs(title="Cholesterol pathway: Genes involved in cholesterol homeostasis") +
  theme_classic() +
  xlab("Gene Rank") + 
  ylab("Enrichment score (ES)")
ggsave("cholesteroltargets.pdf")

## Increased UV response
plotEnrichment(genesets[["HALLMARK_UV_RESPONSE_UP"]], prerank) + labs(title="UV targets")
ggsave("uvtargets.pdf")

# gene set enrichment
gseaplot2(genesets[["HALLMARK_HYPOXIA"]], genesetID = 1, title="Hypoxia pathway: Genes up-regulated in response to low oxygen levels (hypoxia)")
gseaplot2(edo2, geneSetID = 1, title = edo2$Description[1])


# Plot top 10 pathways
topPathwaysUp <- fgseaRes[ES > 0.5][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < -0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(genesets[topPathways], prerank, fgseaRes, 
              gseaParam = 0.5)

allTPMs %>% 
  ggplot() +
  theme_minimal()

gseaplot2

gseaplot2(genesets, prerank, title = "Hypoxia pathway: Genes up-regulated in response to low oxygen levels (hypoxia)")

edo2 <- gseDO(geneList)
