setwd("/Users/kalebvoight/Desktop/Kaleb\ Rotation")

packages = c("BiocManager") 
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

packages = c("DESeq2", "readxl", "dplyr", "utils", "ggplot2", "R.utils", "gplots", "ggrepel", "gprofiler2", "knitr", "kableExtra")

#if the packages do not exist, they will be installed:
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#Read in the table
data <- read.csv("Combined_DoD_Dataset.csv")
head(data)
data <- as_tibble(data)
data <- data  %>%
  dplyr::select(1:7)
rawcounts <- data
rawcounts$Gene_Name[is.na(rawcounts$Gene_Name)] <- paste0("Unknown_Gene_", which(is.na(rawcounts$Gene_Name)))
head(rawcounts)

###Edit treatment columns here
design1<-data.frame(experiment=colnames(rawcounts[,2:7]),  treatment = c("D14", "D14", "D21", "D21", "D28", "D28"))
head(design1)
# add the rownames to the design
rownames(design1)<-colnames(rawcounts[,2:7])
head(design1)

day <- DESeqDataSetFromMatrix(countData = (round(rawcounts[,2:7])), 
                                      colData = design1, design = ~  treatment)
day

dday <- DESeq(day)

###Edit treatment groups here
dday_res <- results(dday,contrast=c("treatment","D14", "D28"))
head(dday_res)
dim(dday_res)


# summary of the data
summary(dday_res)
res05 <- results(dday, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

dday_res@rownames <- rawcounts$Gene_Name
dday_res <- dday_res[order(dday_res$padj),]
head(dday_res)

# Volcano plot
rownames(dday_res) <- rawcounts$Gene_Name
dday_res$Gene_Name <- rawcounts$Gene_Name
length(dday_res$Gene_Name)
length(dday_res$diffexpressed)

p <- ggplot(data=dday_res, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()

p <- ggplot(data=dday_res, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal()

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

dday_res$diffexpressed <- "NO"

dday_res$diffexpressed[dday_res$log2FoldChange > 0.6 & dday_res$pvalue < 0.05] <- "UP"

dday_res$diffexpressed[dday_res$log2FoldChange < -0.6 & dday_res$pvalue < 0.05] <- "DOWN"

p <- ggplot(data=dday_res, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

dday_res$delabel <- NA
dday_res$delabel[dday_res$diffexpressed != "NO"] <- dday_res$Gene_Name[dday_res$diffexpressed != "NO"]

#Plot Volcano Plot
volcano_plot <- ggplot(data=dday_res, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  ggtitle(paste('Volcano Plot:', "D14", 'vs.', "D28"))

#Save Plot as a PNG file
ggsave("Volc_DoD.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

#GProfiler analysis
genes_to_query <- rownames(dday_res)[!is.na(dday_res$padj) & dday_res$padj < 0.05]

if (length(genes_to_query) == 0) {
  message("No genes to query.")
} else {
  #Perform GProfiler enrichment analysis
  gprofiler_results <- tryCatch(
    gprofiler2::gost(genes_to_query, organism = "hsapiens"),
    error = function(e) {
      message("Error during G:Profiler query: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(gprofiler_results) || nrow(gprofiler_results$result) == 0) {
    message("No results returned from G:Profiler.")
  } else {
    gprofiler_df <- as.data.frame(gprofiler_results$result)
    top_terms <- head(gprofiler_df[order(gprofiler_df$p_value), ], 10)
    
    if (nrow(top_terms) > 0) {
      enrich <- ggplot(top_terms, aes(x = reorder(term_name, -p_value), y = -log10(p_value))) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = "Enrichment Analysis of D14 vs D28",
             x = "Term", y = "-log10(P-value)") +
        theme_minimal()
      #Save Plot as a PNG file
      ggsave("Enrich_DoD.png", plot = enrich, width = 8, height = 6, dpi = 300)
    } else {
      message("No significant terms found to plot.")
    }
  }
}

#Select significant genes based on adjusted p-value
significant_genes <- dday_res[!is.na(dday_res$padj) & dday_res$padj < 0.05, ]

#Sort significant genes by absolute log2 fold change
significant_genes <- significant_genes[order(-abs(significant_genes$log2FoldChange)), ]

#Create a table of the top significant genes
top_significant_genes_table <- head(significant_genes, 10)

#Create a data frame
significant_genes_table <- data.frame(
  "Gene Name" = rownames(top_significant_genes_table),
  "Log2 Fold Change" = top_significant_genes_table$log2FoldChange,
  "Adjusted P-value" = format(top_significant_genes_table$padj, scientific = TRUE, digits = 3)
)

# Display the table
if (nrow(significant_genes_table) > 0) {
  kable(significant_genes_table, format = "html", caption = "Significant Genes between D14 and D28") %>%
    kable_styling(full_width = F)
} else {
  message("No significant genes found for display.")
}

#DoD Bar graph
vsd <- vst(dday, blind = FALSE)
treatment_info <- design1
normalized_counts <- assay(vsd)
normalized_counts_df <- as.data.frame(normalized_counts)
normalized_counts_df$Gene_Name <- rawcounts$Gene_Name

#Add Treatment
long_normalized_counts <- normalized_counts_df %>%
  tidyr::pivot_longer(-Gene_Name, names_to = "Sample", values_to = "Normalized_Expression") %>%
  left_join(treatment_info, by = c("Sample" = "experiment"))

#Filter for selected genes and the desired treatment groups
selected_genes <- c("NEUROG1", "NEUROG2", "MAP2", "TUBB3", "SOX2", "NANOG", "POU5F1")
selected_data <- long_normalized_counts %>%
  filter(Gene_Name %in% selected_genes & treatment %in% c("D14", "D21", "D28"))

###Edit colors for graph
custom_colors <- c("D14" = "#0072B2", "D21" = "#D95F02", "D28" = "#999999")

#Calculate mean and standard deviation
summary_data <- selected_data %>%
  group_by(Gene_Name, treatment) %>%
  summarize(Mean_Expression = mean(Normalized_Expression),
            SD = sd(Normalized_Expression))

#Create the bar plot with error bars
DoD <- ggplot(summary_data, aes(x = Gene_Name, y = Mean_Expression, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean_Expression - SD, ymax = Mean_Expression + SD),
                position = position_dodge(0.9), width = 0.25) +
  labs(title = "Relative Expression of Differentiation Genes",
       x = "Gene",
       y = "Normalized Expression Level") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Save Plot as a PNG file
ggsave("Bar_DoD.png", plot = DoD, width = 8, height = 6, dpi = 300)

#Heatmap
vsd <- vst(dday, blind=FALSE)
head(assay(vsd), 3)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$treatment
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library(pheatmap)

#Plot Heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Heatmap of DoD",
         filename = "HM_DoD.png")

#Plot PCA
pca <- plotPCA(vsd, intgroup=c("treatment")) + 
  ggtitle("PCA Modulation Day")

ggsave("PCA_DoD.png", plot = pca, width = 8, height = 6, dpi = 300)
