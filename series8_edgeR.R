getwd()
setwd("/home/behrad/Desktop/Series 8/")
library(edgeR)

#Reading data
Mock_1 <- read.table("expression_results/Series8_A549_Mock_1.genes.results", header = TRUE, sep ="\t", row.names = 1)
Mock_2 <- read.table("expression_results/Series8_A549_Mock_2.genes.results", header = TRUE, sep ="\t", row.names = 1)
Mock_3 <- read.table("expression_results/Series8_A549_Mock_3.genes.results", header = TRUE, sep ="\t", row.names = 1)

RSV_1 <- read.table("expression_results/Series8_A549_RSV_1.genes.results", header = TRUE, sep ="\t", row.names = 1)
RSV_2 <- read.table("expression_results/Series8_A549_RSV_2.genes.results", header = TRUE, sep ="\t", row.names = 1)
RSV_3 <- read.table("expression_results/Series8_A549_RSV_3.genes.results", header = TRUE, sep ="\t", row.names = 1)

#Creating expected counts table from RSEM expected counts of every sample  

Expected_Counts_Table <- cbind(Mock_1[,4], Mock_2[,4], Mock_3[,4], RSV_1[,4], RSV_2[,4], RSV_3[,4])

colnames(Expected_Counts_Table) <- c("Mock_1", "Mock_2", "Mock_3", "RSV_1", "RSV_2", "RSV_3")

rownames(Expected_Counts_Table) <- row.names(Mock_1)

head(Expected_Counts_Table) #let's see what it looks like


#Creating edgeR DGEList object from expected counts table
DG <- DGEList(counts = Expected_Counts_Table, group = c(rep("Mock", 3), rep("RSV", 3)),
              samples = colnames(Expected_Counts_Table), remove.zeros = TRUE)


#Having a look 
head(Expected_Counts_Table)
head(DG$counts)
head(cpm(DG))
dim(Expected_Counts_Table)
dim(cpm(DG))


#Removing genes with very low expression (keeping those with at least 10 cpm in at least 2 samples)
DG_keep <- rowSums(cpm(DG) > 10) >= 2
DG <- DG[DG_keep, ,keep.lib.sizes = FALSE]
dim(cpm(DG))


#Calculate scaling factors to normalize library sizes - TMM NORMALIZATION 
DG$samples$norm.factors  #Initial factors (all ones, we still need to normalize)


head(cpm(DG))            #Non normalized CPM 



DG <- calcNormFactors(DG)  #We compute the scaling factors for each library

DG$samples$norm.factors    #The scaling factors


head(cpm(DG)) 


#Multidimensional scaling plot of samples
plotMDS(DG, cex = 0.6, col = c(rep("blue", 3), rep("red", 3)))
DG <- estimateDisp(DG)
head(DG$tagwise.dispersion)
plotBCV(DG)

#Perform test for differential expression between the two conditions. 
Diff_expr_test <- exactTest(DG, c("Mock", "RSV"))

#Select genes with significant diff. expr. (FDR <= 0.01) 
Diff_expr_test_sign <- topTags(Diff_expr_test, n = nrow(Diff_expr_test), p.value = 0.01)

#Having a look at the DEGs

head(Diff_expr_test_sign$table)


tail(Diff_expr_test_sign$table)


#Number of upregulaated genes - that is genes that are more expressed in Infected w.r.t. Mock
sum(Diff_expr_test_sign$table$logFC > 0)


#Number of downregulated genes - that is genes that are less expressed in Infected w.r.t. Mock
sum(Diff_expr_test_sign$table$logFC < 0)

UP_regulated_genes <- row.names(Diff_expr_test_sign$table)[Diff_expr_test_sign$table$logFC > 0]
DOWN_regulated_genes <- row.names(Diff_expr_test_sign$table)[Diff_expr_test_sign$table$logFC < 0]

boxplot(cpm(DG)[UP_regulated_genes,], notch = TRUE, outline = FALSE, col = c(rep("blue", 3),
                                                                             rep("red", 3)), ylab = "CPM", xlab = "Samples", main = "CPM Distr. UP Regulated genes in RSV infected A549 cells")

boxplot(cpm(DG)[DOWN_regulated_genes,], notch = TRUE, outline = FALSE, col = c(rep("blue", 3), rep("red", 3)), ylab = "CPM", xlab = "Samples", main = "CPM Distr. DOWN Regulated genes in RSV infected A549 cells")


colors <- ifelse(Diff_expr_test$table$PValue <= 0.01, ifelse(Diff_expr_test$table$logFC > 0, "red", "blue"), "black")

plot(Diff_expr_test$table$logFC, -log10(Diff_expr_test$table$PValue), cex = 0.1, ylim = c(0,50), ylab = "-log10 Pvalue", xlab = "LogFC", main = "Vulcano Plot", col = colors)



library("vidger")

#Scatterplot
vsScatterPlot(x = "Mock", y = "RSV", data = DG, type = "edger")

#Volcano plot i.e. a scatterplot of logFC vs FDR (adh pv) 
vsVolcano(x = "Mock", y = "RSV", data = DG, type = "edger", padj = 0.01)

#MA plot i.e. a scatter plot of expression vs logFC
vsMAPlot(x = "Mock", y = "RSV", data = DG, type = "edger", padj = 0.01)


library(openxlsx)

wb <- createWorkbook()

addWorksheet(wb, "Expected Counts")
addWorksheet(wb, "CPM")
addWorksheet(wb, "Diff.Expr.All.Genes")
addWorksheet(wb, "Diff.Expr.SIGN")

writeData(wb = wb, sheet = "Expected Counts", x = Expected_Counts_Table, rowNames = TRUE)
writeData(wb = wb, sheet = "CPM", x = cpm(DG), rowNames = TRUE)
writeData(wb = wb, sheet = "Diff.Expr.All.Genes", x = Diff_expr_test$table, rowNames = TRUE)
writeData(wb = wb, sheet = "Diff.Expr.SIGN", x = Diff_expr_test_sign$table, rowNames = TRUE)

saveWorkbook(wb, "Diff_expr_analysis_series_8.xlsx", overwrite = TRUE)



