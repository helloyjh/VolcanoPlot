rank <- 'O'
save_path = paste0('Desktop/sub_count_ranked_v1.05/subject_count_rank_', rank, '.csv')
print(save_path)
cnt_data <- read.csv(file=save_path)
X_data <- subset(cnt_data, select = -c(0:3))
X_data <- t(X_data)
X_data <- X_data*1000000000
X_data <- round(X_data)

label <- rep(0L, length(cnt_data$diag))
label[which(cnt_data$diag=='ADD')] <- 1
Group <- as.factor(label)
colData <- as.data.frame(Group)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=X_data,
                              colData=colData, 
                              design=~Group)
                 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, pAdjustMethod="BH")
topT <- as.data.frame(res)

tsv_name <- paste0('Desktop/R_result/DESeq2/p_rank_', rank, '.tsv')
write.table(topT, file=tsv_name, quote=FALSE, sep='\t', col.names=NA)
library(EnhancedVolcano)
volc_name = paste0('Desktop/R_result/DESeq2/volcano_rank_', rank, '.pdf')
pdf(file=volc_name)
EnhancedVolcano(res,
                lab = rownames(res),
                title = paste0('Volcano Plot rank_', rank),
                pCutoff = 0.05,
                x = 'log2FoldChange',
                y = 'padj')
dev.off()

