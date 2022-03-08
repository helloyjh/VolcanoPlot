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
library(edgeR)

d <- DGEList(counts=X_data, group=Group)
d <- calcNormFactors(d)
d <- estimateDisp(d)

et <- exactTest(d)
res <- topTags(et, n = nrow(X_data), sort.by = "none")
res_tab <- res$table

tsv_name <- paste0('Desktop/R_result/edgeR/p_rank_', rank, '.tsv')
write.table(res_tab, file=tsv_name, quote=FALSE, sep='\t', col.names=NA)
library(EnhancedVolcano)
volc_name = paste0('Desktop/R_result/edgeR/volcano_rank_', rank, '.pdf')
pdf(file=volc_name)
EnhancedVolcano(res_tab, 
                lab = rownames(res_tab),
                title = paste0('Volcano Plot rank_', rank),
                x = 'logFC',
                y = 'FDR',
                pCutoff = 0.05,
                FCcutoff = 1)
dev.off()