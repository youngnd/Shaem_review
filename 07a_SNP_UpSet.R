#install.packages("vcfR")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("SNPRelate", quietly = TRUE))
  BiocManager::install("SNPRelate")
if (!requireNamespace("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2")
if (!requireNamespace("gridExtra", quietly = TRUE))
  BiocManager::install("gridExtra")
if (!requireNamespace("UpSetR", quietly = TRUE))
  BiocManager::install("UpSetR")
if (!requireNamespace("chromoMap", quietly = TRUE))
  BiocManager::install("chromoMap")
if (!requireNamespace("ggbio", quietly = TRUE))
  BiocManager::install("ggbio")
library(vcfR)
library(ggplot2)
library(gridExtra)
library(UpSetR)
library(chromoMap)
library(ggbio)
## Reading in positions of SNPs for each library. These are fixed (GN 1/1) SNPs in exonic regions
dataMi <- read.table(file="SMi.snp_positions.txt", sep="\t", header = FALSE)
dataMs <- read.table(file="SMs.snp_positions.txt", sep="\t", header = FALSE)
dataS4 <- read.table(file="SS4.snp_positions.txt", sep="\t", header = FALSE)
dataZ1 <- read.table(file="SZ1.snp_positions.txt", sep="\t", header = FALSE)
Mali <- dataMi$V1
Mauritius <- dataMs$V1
Senegal <- dataS4$V1
Zambia <- dataZ1$V1

read_sets = list(Mali = Mali,
                 Mauritius = Mauritius,
                 Senegal = Senegal,
                 Zambia = Zambia)

svglite::svglite(file="UpSet_SNPS.svg",width = 6,height=5)
upset(fromList(read_sets), 
      sets = c("Mali","Mauritius","Senegal","Zambia"), 
      keep.order = TRUE,order.by = "freq",
      sets.bar.color=c("#f99790ff","#7cad00c0","#00bdc3c0","#c57cffc0"),
#      sets.bar.color=c("#D73027","#FC8D59","#FEE090", "#91BFDB"),
      mb.ratio=c(0.7, 0.3), point.size=2, mainbar.y.label = "Number of SNPs",
      sets.x.label = "Total number of SNPs",show.numbers = "yes",
      text.scale = c(0.8, 0.8, 0.8, 0.8, 1, 1))
dev.off()