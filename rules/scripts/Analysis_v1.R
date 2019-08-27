# Author Paolo Uva
# 16 gennaio 2019 - Aggiornato script per usare dati normalizzati nei plot
#                   For details: http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-variance-stabilizing-transformation-and-the-rlog

# Species: human
# Protocol: mRNA (poly-A)
#os.getcwd()
# Per installare pacchetti:
# - CRAN: install.packages("NOME")
# - Bioconductor biocLite("NOME")


# RNA-Seq DE analysis with kallisto + DeSeq2
library("tximport")
library("genefilter")
library("pheatmap")
library("DESeq2")
library("gplots")
library("RColorBrewer")
library("ggplot2")
library("dplyr")

#setwd("/home/matteo/Scrivania/analysis")
setwd(snakemake@params[["current_dir"]])

#load("tx2gene_hsa.RData")
load(snakemake@params[["txgenes"]])
# File list
files <- file.path("kallisto", dir("kallisto"), "abundance.tsv")

names(files) <- gsub(".*/(.*)/.*", "\\1", files)

# Load and summarize at gene level

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene[,c(1,3)], ignoreTxVersion=TRUE)
save(txi, file = "txi.RData")

head(txi$counts)


# Run DeSeq2
# Create sampleTable from sample names
sampleTable <- data.frame(sample = factor(names(files)))
rownames(sampleTable) <- colnames(txi$counts)

# No design
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ 1)
dds <- DESeq(dds)



# Compare log2, vst and rlog transformations
# From http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-variance-stabilizing-transformation-and-the-rlog


dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = FALSE) # Suggested instead of rlog for > 30 samples because it's faster
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst")
)
colnames(df)[1:2] <- c("x", "y")  
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  



# Heatmap of most variable genes using rlog normalized data
# Select most variable genes by interquartile range
vsd.filt <- varFilter(assay(vsd),
				  var.func=IQR,
				  var.cutoff=0.9,
				  filterByQuantile=TRUE)
dim(vsd.filt)

#heatName_1 <- "Heatmap_Most_Var.png"
heatName_1 <- snakemake@output[["heatmap"]]
pheatmap(vsd.filt, # most variable
		 color = greenred(75),
		 #breaks = c(min(mat, na.rm = TRUE), seq(-3, 3, 6 / (pheatmap_n_colors - 2)), max(mat, na.rm = TRUE)),
		 cluster_rows = TRUE,
		 clustering_distance_rows = "correlation",
		 #clustering_distance_cols = "euclidean",
		 cluster_cols = TRUE,
		 #annotation = group_names,
		 #fontsize_row = 6,
		 #fontsize_col = 6,
		 scale = "row",
		 show_colnames = TRUE,
		 show_rownames = FALSE#,
		 #legend_breaks=c(-3,0,3),
		 #legend_labels=c("< -3","0","> 3"),
		 #filename=heatName_1
		 )

# PCA plot
plotPCA(vsd, intgroup=colnames(colData(dds)), ntop=1000)


# PCA plot with colors and shapes
library("ggplot2")
data <- plotPCA(vsd, intgroup=colnames(colData(dds)), ntop=1000, returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

### STOP HERE ###


#figName_3 <- "PCA.png"
figName_3 <- snakemake@output[["pca"]]
#png(filename=figName_3)
my_label=rep('',192)
my_label[names(files) %in% c("SAMPLE-IN-FC19-0319-R0001","UNDETERMINED")]="OK"

ggplot(data, aes(PC1, PC2, color="red")) +
	   geom_point(size=3) +
	   geom_text(data=data, mapping=aes(x=PC1+2, y=PC2+2, label=my_label), size=15) +
	   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
	   coord_fixed()
#dev.off()

