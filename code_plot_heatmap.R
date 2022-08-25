library(ComplexHeatmap)
library(circlize)

setwd("C:/Users/gaoda/OneDrive - USTC/Gao lab/Scripts/example_heatmap_volcano_PZ")
data = read.table('data_IFNa_ISG_zscore.txt', sep='\t', header=TRUE, row.names=1)

sample = read.table('data_IFNa_samples_annot.txt', sep='\t', header=FALSE, row.names=1)
colnames(sample) = c('genotype', 'condition')
sample_df = data.frame(genotype=sample$genotype, condition=sample$condition)
sample_annot = HeatmapAnnotation(df=sample_df, gap=unit(1,"mm"), width=unit(5,"mm"), 
								 col=list(genotype=c('CTL'='#00cc99','SNORA31_3H'='#66ccff','SNORA31_8H'='#3399ff','SNORA31_P'='#0066ff','STAT1'='#ff9933','TLR3'='#cc66ff'),
								 condition=c('IFNa'='#4d4d4d','NS'='#e6e6e6')))

DEG = read.table('data_IFNa_ISG_list_annot.txt', sep='\t', header=FALSE, row.names=1)
colnames(DEG) = c('CTL_DEG', 'SNORA31_DEG', 'STAT1_DEG')
DEG_df = data.frame(CTL_DEG=DEG$CTL_DEG, SNORA31_DEG=DEG$SNORA31_DEG, STAT1_DEG=DEG$STAT1_DEG)
DEG_annot = rowAnnotation(df=DEG_df, gap=unit(1,"mm"), width=unit(10,"mm"), 
								 col=list(CTL_DEG=c('Sig.'='#00cc99','Not Sig.'='#cccccc'),
								 		  SNORA31_DEG=c('Sig.'='#0066ff','Not Sig.'='#cccccc'),
								 		  STAT1_DEG=c('Sig.'='#ff9933','Not Sig.'='#cccccc')))

data_heatmap = Heatmap(data, name='z-score', col=colorRamp2(c(-2.1,0,2.1),c('#0000ff','#ffffff','#ff6600')),
				row_title='IFNa-Stimulated Genes', row_title_side='left', column_title_side='bottom',
				clustering_method_rows='complete', clustering_distance_rows='manhattan', row_dend_width=unit(10,'mm'),
				clustering_method_columns='complete', clustering_distance_columns='euclidean', column_dend_height=unit(20,'mm'),
				show_row_names=FALSE, show_column_names=TRUE, column_names_gp=gpar(fontsize=9), show_heatmap_legend=TRUE, top_annotation=sample_annot)
png(file = "C:/Users/gaoda/OneDrive - USTC/Gao lab/Scripts/example_heatmap_volcano_PZ/figure_IFNa_ISG_zscore_heatmap_w_DEG_annot2.png", width=9, height=9, units='in', res=300)
draw(data_heatmap + DEG_annot)
dev.off()

