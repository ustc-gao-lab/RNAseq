
library(edgeR)

express_data = read.table('C:/Users/gaoda/OneDrive - USTC/Gao lab/Scripts/example_heatmap_volcano_PZ/volcano/data_read_counts_HSV_vs_NI.txt', sep='\t', header=TRUE, row.names=1)
N = nrow(express_data)
sample = read.table('C:/Users/gaoda/OneDrive - USTC/Gao lab/Scripts/example_heatmap_volcano_PZ/volcano/data_samples_class.txt', sep='\t', header=FALSE, row.names=1)
colnames(sample) = c('class')


### CTL_HSV vs CTL_NI ###

group = factor(sample$class)
y = DGEList(counts=express_data, genes=rownames(express_data), group=group)
y$samples$group = relevel(y$samples$group, ref="CTL_NI")									
y = calcNormFactors(y, method="TMM")
design = model.matrix(~0+group, data=y$samples)
colnames(design) = levels(y$samples$group)
y = estimateDisp(y, design)
fit = glmFit(y, design)

comparison = makeContrasts(CTL_HSV-CTL_NI, levels=design)									
DEG_model = glmLRT(fit, contrast=comparison)
DEG_output = topTags(DEG_model, n=N)
write.table(DEG_output, "DEG_(CTL_HSV)_vs_(CTL_NI).txt", sep="\t", quote=FALSE)		
summary(de <- decideTestsDGE(DEG_model))

png("figure_volcano_(CTL_HSV)_vs_(CTL_NI).png", width=5, height=5, units='in', res=300)
plot(DEG_output$table$logFC, -log10(DEG_output$table$PValue), main="(CTL_HSV)_vs_(CTL_NI)", xlab="log2(FC)", ylab="-log10(pValue)", xlim=c(-10, 12), pch=16, cex=0.7)
with(subset(DEG_output$table, logFC >= 1 & PValue < 0.05 & FDR < 0.05), points(logFC, -log10(PValue), pch=16, cex=0.7, col="red"))
with(subset(DEG_output$table, logFC <= -1 & PValue < 0.05 & FDR < 0.05), points(logFC, -log10(PValue), pch=16, cex=0.7, col="blue"))
dev.off()


### TLR3_HSV vs TLR3_NI ###

group = factor(sample$class)
y = DGEList(counts=express_data, genes=rownames(express_data), group=group)
y$samples$group = relevel(y$samples$group, ref="TLR3_NI")									
y = calcNormFactors(y, method="TMM")
design = model.matrix(~0+group, data=y$samples)
colnames(design) = levels(y$samples$group)
y = estimateDisp(y, design)
fit = glmFit(y, design)

comparison = makeContrasts(TLR3_HSV-TLR3_NI, levels=design)									
DEG_model = glmLRT(fit, contrast=comparison)
DEG_output = topTags(DEG_model, n=N)
write.table(DEG_output, "DEG_(TLR3_HSV)_vs_(TLR3_NI).txt", sep="\t", quote=FALSE)	
summary(de <- decideTestsDGE(DEG_model))

png("figure_volcano_(TLR3_HSV)_vs_(TLR3_NI).png", width=5, height=5, units='in', res=300)
plot(DEG_output$table$logFC, -log10(DEG_output$table$PValue), main="(TLR3_HSV)_vs_(TLR3_NI)", xlab="log2(FC)", ylab="-log10(pValue)", xlim=c(-10, 12), pch=16, cex=0.7)
with(subset(DEG_output$table, logFC >= 1 & PValue < 0.05 & FDR < 0.05), points(logFC, -log10(PValue), pch=16, cex=0.7, col="red"))
with(subset(DEG_output$table, logFC <= -1 & PValue < 0.05 & FDR < 0.05), points(logFC, -log10(PValue), pch=16, cex=0.7, col="blue"))
dev.off()
