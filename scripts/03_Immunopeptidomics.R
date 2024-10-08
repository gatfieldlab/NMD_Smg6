### Immunopeptidomics results and plots
### Virginie Ricci - October 2024

######################################################################
library(RColorBrewer)
library(randomcoloR)
library(DEXSeq)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(biomaRt)
library(stringr)
library(ggplot2)
library(tximport)
library(gridExtra)
library(ggrepel)
library(ggfortify)
library(ggpubr)
library(ggsignif)
library(cowplot)
library(ggh4x)
library(reshape2)
library(scales)
library(DESeq2)
library(Rsubread)
library(GenomicFeatures)
library(GenomicAlignments)
library(txdbmaker)
library(Rsamtools)
library(rtracklayer)
library(apeglm)
library(ashr)
library(ensembldb)
library(eisaR)
library(Repitools)
library(BiocParallel)
library(EnhancedVolcano)
library(corrplot)

######################################################################
### Immunopeptidomics results table
path_Immuno=''

df = read.csv(paste0(path_Immuno, 'immunopep_masterdf.csv'), sep=',', header=T)
head(df); dim(df)

df$sample=''
df[df$Genotype=='AlbCreERT2/CreERT2\nSmg6-wt/wt',]$sample = 'Ctrl_Smg6wt'
df[df$Genotype=='AlbCreERT2/CreERT2\nSmg6-mutant/mutant',]$sample = 'Ctrl_Smg6mut'
df[df$Genotype=='AlbCreERT2/CreERT2\nSmg6-wt/wt\nPten-ko/ko\nTsc1-ko/ko',]$sample = 'LdKO_Smg6wt'
df[df$Genotype=='AlbCreERT2/CreERT2\nSmg6-mutant/mutant\nPten-ko/ko\nTsc1-ko/ko',]$sample = 'LdKO_Smg6mut'
df$Background=''
df[grep('Ctrl', df$sample),]$Background = 'Ctrl'
df[grep('LdKO', df$sample),]$Background = 'LdKO'
df[df$pepsource == 'intronin',]$pepsource = 'intronic'

df$sampleID = paste0(df$sample, df$BR)


df = df[df$bindcat == 'binder', ]
df = df[df$sample_nr != 8, ]


Ctrl_Smg6wt1 = data.frame(psmed_peptide = df[df$sampleID == 'Ctrl_Smg6wt1',]$psmed_peptide, Ctrl_Smg6wt1 = df[df$sampleID == 'Ctrl_Smg6wt1',]$psmed_peptide_intensity)
Ctrl_Smg6wt2 = data.frame(psmed_peptide = df[df$sampleID == 'Ctrl_Smg6wt2',]$psmed_peptide, Ctrl_Smg6wt2 = df[df$sampleID == 'Ctrl_Smg6wt2',]$psmed_peptide_intensity)
Ctrl_Smg6wt3 = data.frame(psmed_peptide = df[df$sampleID == 'Ctrl_Smg6wt3',]$psmed_peptide, Ctrl_Smg6wt3 = df[df$sampleID == 'Ctrl_Smg6wt3',]$psmed_peptide_intensity)


Ctrl_Smg6mut1 = data.frame(psmed_peptide = df[df$sampleID == 'Ctrl_Smg6mut1',]$psmed_peptide, Ctrl_Smg6mut1 = df[df$sampleID == 'Ctrl_Smg6mut1',]$psmed_peptide_intensity)
Ctrl_Smg6mut2 = data.frame(psmed_peptide = df[df$sampleID == 'Ctrl_Smg6mut2',]$psmed_peptide, Ctrl_Smg6mut2 = df[df$sampleID == 'Ctrl_Smg6mut2',]$psmed_peptide_intensity)
Ctrl_Smg6mut3 = data.frame(psmed_peptide = df[df$sampleID == 'Ctrl_Smg6mut3',]$psmed_peptide, Ctrl_Smg6mut3 = df[df$sampleID == 'Ctrl_Smg6mut3',]$psmed_peptide_intensity)


LdKO_Smg6wt1 = data.frame(psmed_peptide = df[df$sampleID == 'LdKO_Smg6wt1',]$psmed_peptide, LdKO_Smg6wt1 = df[df$sampleID == 'LdKO_Smg6wt1',]$psmed_peptide_intensity)
LdKO_Smg6wt2 = data.frame(psmed_peptide = df[df$sampleID == 'LdKO_Smg6wt2',]$psmed_peptide, LdKO_Smg6wt2 = df[df$sampleID == 'LdKO_Smg6wt2',]$psmed_peptide_intensity)
LdKO_Smg6wt3 = data.frame(psmed_peptide = df[df$sampleID == 'LdKO_Smg6wt3',]$psmed_peptide, LdKO_Smg6wt3 = df[df$sampleID == 'LdKO_Smg6wt3',]$psmed_peptide_intensity)
LdKO_Smg6wt4 = data.frame(psmed_peptide = df[df$sampleID == 'LdKO_Smg6wt4',]$psmed_peptide, LdKO_Smg6wt4 = df[df$sampleID == 'LdKO_Smg6wt4',]$psmed_peptide_intensity)


LdKO_Smg6mut1 = data.frame(psmed_peptide = df[df$sampleID == 'LdKO_Smg6mut1',]$psmed_peptide, LdKO_Smg6mut1 = df[df$sampleID == 'LdKO_Smg6mut1',]$psmed_peptide_intensity)
LdKO_Smg6mut2 = data.frame(psmed_peptide = df[df$sampleID == 'LdKO_Smg6mut2',]$psmed_peptide, LdKO_Smg6mut2 = df[df$sampleID == 'LdKO_Smg6mut2',]$psmed_peptide_intensity)
LdKO_Smg6mut3 = data.frame(psmed_peptide = df[df$sampleID == 'LdKO_Smg6mut3',]$psmed_peptide, LdKO_Smg6mut3 = df[df$sampleID == 'LdKO_Smg6mut3',]$psmed_peptide_intensity)
LdKO_Smg6mut4 = data.frame(psmed_peptide = df[df$sampleID == 'LdKO_Smg6mut4',]$psmed_peptide, LdKO_Smg6mut4 = df[df$sampleID == 'LdKO_Smg6mut4',]$psmed_peptide_intensity)



full_df_matrix = data.frame(psmed_peptide=unique(df$psmed_peptide))
head(full_df_matrix); dim(full_df_matrix) # 2370

full_df_matrix = left_join(full_df_matrix, Ctrl_Smg6wt1, by='psmed_peptide')
full_df_matrix = left_join(full_df_matrix, Ctrl_Smg6wt2, by='psmed_peptide')
full_df_matrix = left_join(full_df_matrix, Ctrl_Smg6wt3, by='psmed_peptide')
full_df_matrix = left_join(full_df_matrix, Ctrl_Smg6mut1, by='psmed_peptide')
full_df_matrix = left_join(full_df_matrix, Ctrl_Smg6mut2, by='psmed_peptide')
full_df_matrix = left_join(full_df_matrix, Ctrl_Smg6mut3, by='psmed_peptide')

full_df_matrix = left_join(full_df_matrix, LdKO_Smg6wt1, by='psmed_peptide')
#full_df_matrix = left_join(full_df_matrix, LdKO_Smg6wt2, by='psmed_peptide') # excluded according to immunopeptidomics results
full_df_matrix = left_join(full_df_matrix, LdKO_Smg6wt3, by='psmed_peptide')
full_df_matrix = left_join(full_df_matrix, LdKO_Smg6wt4, by='psmed_peptide')

full_df_matrix = left_join(full_df_matrix, LdKO_Smg6mut1, by='psmed_peptide')
full_df_matrix = left_join(full_df_matrix, LdKO_Smg6mut2, by='psmed_peptide')
full_df_matrix = left_join(full_df_matrix, LdKO_Smg6mut3, by='psmed_peptide')
full_df_matrix = left_join(full_df_matrix, LdKO_Smg6mut4, by='psmed_peptide')

head(full_df_matrix); dim(full_df_matrix) # 2370


full_df_matrix[seq(2, ncol(full_df_matrix))] = apply(full_df_matrix[seq(2, ncol(full_df_matrix))], 2, as.numeric)

# matrix as log
full_df_matrix_log = full_df_matrix

# 0 for normal matrix
full_df_matrix[is.na(full_df_matrix)] = 0

# 1 for log matrix so that log(1) == 0
full_df_matrix_log[is.na(full_df_matrix_log)] = 1

head(full_df_matrix)
head(full_df_matrix_log)

# log transformation
full_df_matrix_log[seq(2, ncol(full_df_matrix_log))] = apply(full_df_matrix_log[seq(2, ncol(full_df_matrix_log))], 2, log10)


head(full_df_matrix); dim(full_df_matrix)
head(full_df_matrix_log); dim(full_df_matrix_log)

rownames(full_df_matrix) = full_df_matrix$psmed_peptide
rownames(full_df_matrix_log) = full_df_matrix_log$psmed_peptide

minimum = 1


### 
full_df_matrix_wt = full_df_matrix[grep('wt', names(full_df_matrix))]
head(full_df_matrix_wt); dim(full_df_matrix_wt)
full_df_matrix_wt$psmed_peptide = rownames(full_df_matrix_wt)


full_df_matrix_mut = full_df_matrix[grep('mut', names(full_df_matrix))]
head(full_df_matrix_mut); dim(full_df_matrix_mut)
full_df_matrix_mut$psmed_peptide = rownames(full_df_matrix_mut)

setdiff(rownames(full_df_matrix_wt), rownames(full_df_matrix_mut))
setdiff(rownames(full_df_matrix_mut), rownames(full_df_matrix_wt))


full_df_matrix = data.frame(psmed_peptide=sort(unique(c(rownames(full_df_matrix_wt), rownames(full_df_matrix_mut)))))
head(full_df_matrix); dim(full_df_matrix) # 2323
full_df_matrix = left_join(full_df_matrix, full_df_matrix_wt)
full_df_matrix = left_join(full_df_matrix, full_df_matrix_mut)

full_df_matrix = full_df_matrix[names(full_df_matrix)]
head(full_df_matrix); dim(full_df_matrix)
rownames(full_df_matrix) = full_df_matrix$psmed_peptide
full_df_matrix[is.na(full_df_matrix)] = 0


# Ctrl
Ctrl_Smg6wt_median = data.frame(Ctrl_Smg6wt_median = apply(full_df_matrix[grep('Ctrl_Smg6wt', colnames(full_df_matrix))], 1, median))
head(Ctrl_Smg6wt_median); dim(Ctrl_Smg6wt_median)

Ctrl_Smg6mut_median = data.frame(Ctrl_Smg6mut_median = apply(full_df_matrix[grep('Ctrl_Smg6mut', colnames(full_df_matrix))], 1, median))
head(Ctrl_Smg6mut_median); dim(Ctrl_Smg6mut_median)
# Ctrl


# LdKO
LdKO_Smg6wt_median = data.frame(LdKO_Smg6wt_median = apply(full_df_matrix[grep('LdKO_Smg6wt', colnames(full_df_matrix))], 1, median))
head(LdKO_Smg6wt_median); dim(LdKO_Smg6wt_median)

LdKO_Smg6mut_median = data.frame(LdKO_Smg6mut_median = apply(full_df_matrix[grep('LdKO_Smg6mut', colnames(full_df_matrix))], 1, median))
head(LdKO_Smg6mut_median); dim(LdKO_Smg6mut_median)
# LdKO


# Ctrl
Ctrls = cbind(Ctrl_Smg6wt_median, Ctrl_Smg6mut_median)
head(Ctrls); dim(Ctrls)
Ctrls$psmed_peptide = rownames(Ctrls)
Ctrls = left_join(Ctrls, unique(df[c('psmed_peptide', 'pepsource')]))


Ctrls_log = Ctrls
Ctrls_log[Ctrls_log == 0] = 1
Ctrls_log[c(1:2)] = apply(Ctrls_log[c(1:2)], 2, log10)
head(Ctrls_log); dim(Ctrls_log)
Ctrls_log = Ctrls_log[order(Ctrls_log$pepsource),]


ggplot(Ctrls_log) + geom_point(aes(Ctrl_Smg6wt_median, Ctrl_Smg6mut_median), alpha=0.3) + ggtitle('Ctrl - log transformed') + geom_abline(intercept = 0, slope = 1) + theme_bw() 
ggplot(Ctrls_log) + geom_point(aes(Ctrl_Smg6wt_median, Ctrl_Smg6mut_median, col=pepsource)) + ggtitle('Ctrl - log transformed') + geom_abline(intercept = 0, slope = 1) + theme_bw() 
ggplot(Ctrls_log) + geom_point(aes(Ctrl_Smg6wt_median, Ctrl_Smg6mut_median, col=pepsource)) + ggtitle('Ctrl - log transformed') + facet_grid(. ~ pepsource) + geom_abline(intercept = 0, slope = 1) + theme_bw() 
# Ctrl


# LdKO
LdKOs = cbind(LdKO_Smg6wt_median, LdKO_Smg6mut_median)
head(LdKOs); dim(LdKOs)
LdKOs$psmed_peptide = rownames(LdKOs)
LdKOs = left_join(LdKOs, unique(df[c('psmed_peptide', 'pepsource')]))


LdKOs_log = LdKOs
LdKOs_log[LdKOs_log == 0] = 1
LdKOs_log[c(1:2)] = apply(LdKOs_log[c(1:2)], 2, log10)
head(LdKOs_log); dim(LdKOs_log)
LdKOs_log = LdKOs_log[order(LdKOs_log$pepsource),]


ggplot(LdKOs_log) + geom_point(aes(LdKO_Smg6wt_median, LdKO_Smg6mut_median), alpha=0.3) + ggtitle('LdKO - log transformed') + geom_abline(intercept = 0, slope = 1) + theme_bw() 
ggplot(LdKOs_log) + geom_point(aes(LdKO_Smg6wt_median, LdKO_Smg6mut_median, col=pepsource)) + ggtitle('LdKO - log transformed') + geom_abline(intercept = 0, slope = 1) + theme_bw() 
ggplot(LdKOs_log) + geom_point(aes(LdKO_Smg6wt_median, LdKO_Smg6mut_median, col=pepsource)) + ggtitle('LdKO - log transformed') + facet_grid(. ~ pepsource) + geom_abline(intercept = 0, slope = 1) + theme_bw() 
# LdKO





### frequency table
head(full_df_matrix)

full_df_matrix_melt = melt(full_df_matrix)

head(full_df_matrix_melt); dim(full_df_matrix_melt)
full_df_matrix_melt$variable = as.character(full_df_matrix_melt$variable)

full_df_matrix_melt$sample = str_sub(full_df_matrix_melt$variable, end=-2)
head(full_df_matrix_melt)
full_df_matrix_melt = full_df_matrix_melt[full_df_matrix_melt$value != 0,]


# Ctrl
full_df_matrix_melt_Ctrl_wt = full_df_matrix_melt[grep('Ctrl_Smg6wt', full_df_matrix_melt$sample),]
head(full_df_matrix_melt_Ctrl_wt); dim(full_df_matrix_melt_Ctrl_wt)

Ctrl_wt_tbl = data.frame(psmed_peptide = names(table(full_df_matrix_melt_Ctrl_wt$psmed_peptide)), Freq_wt = cbind(table(full_df_matrix_melt_Ctrl_wt$psmed_peptide)))
head(Ctrl_wt_tbl); dim(Ctrl_wt_tbl)


full_df_matrix_melt_Ctrl_mut = full_df_matrix_melt[grep('Ctrl_Smg6mut', full_df_matrix_melt$sample),]
head(full_df_matrix_melt_Ctrl_mut); dim(full_df_matrix_melt_Ctrl_mut)

Ctrl_mut_tbl = data.frame(psmed_peptide = names(table(full_df_matrix_melt_Ctrl_mut$psmed_peptide)), Freq_mut = cbind(table(full_df_matrix_melt_Ctrl_mut$psmed_peptide)))
head(Ctrl_mut_tbl); dim(Ctrl_mut_tbl)


Ctrls_tbl = full_join(Ctrl_wt_tbl, Ctrl_mut_tbl)
head(Ctrls_tbl); dim(Ctrls_tbl) # 2213
Ctrls_tbl[is.na(Ctrls_tbl)] = 0
table(Ctrls_tbl$Freq_wt, Ctrls_tbl$Freq_mut)

Ctrls_df = data.frame(rbind(table(Ctrls_tbl$Freq_wt, Ctrls_tbl$Freq_mut)))
colnames(Ctrls_df) = seq(0, ncol(Ctrls_df)-1)
Ctrls_df
Ctrls_df$WT = rownames(Ctrls_df)
Ctrl_mdf <- gather(Ctrls_df, -WT, key="MUT", value="count")
Ctrl_mdf
Ctrl_mdf_copy = Ctrl_mdf
Ctrl_mdf$WT = factor(Ctrl_mdf$WT, levels=rev(seq(0, max(Ctrl_mdf$WT))))

pdf(file=paste0(path_Immuno, 'PeptideIntensity_Matrices_Ctrl.pdf'), width = 8, height=8)
ggplot(Ctrl_mdf, aes(x=MUT, y=WT, fill=log10(count), label=count))+
  geom_tile()+
  geom_text()+
  scale_fill_gradient(high="red2", low="white") + theme_classic() + ggtitle('Ctrl') + scale_x_discrete(position = "top") + ylab('WT') + coord_fixed()
dev.off()
# Ctrl


# LdKO
full_df_matrix_melt_LdKO_wt = full_df_matrix_melt[grep('LdKO_Smg6wt', full_df_matrix_melt$sample),]
head(full_df_matrix_melt_LdKO_wt); dim(full_df_matrix_melt_LdKO_wt)

LdKO_wt_tbl = data.frame(psmed_peptide = names(table(full_df_matrix_melt_LdKO_wt$psmed_peptide)), Freq_wt = cbind(table(full_df_matrix_melt_LdKO_wt$psmed_peptide)))
head(LdKO_wt_tbl); dim(LdKO_wt_tbl)


full_df_matrix_melt_LdKO_mut = full_df_matrix_melt[grep('LdKO_Smg6mut', full_df_matrix_melt$sample),]
head(full_df_matrix_melt_LdKO_mut); dim(full_df_matrix_melt_LdKO_mut)

LdKO_mut_tbl = data.frame(psmed_peptide = names(table(full_df_matrix_melt_LdKO_mut$psmed_peptide)), Freq_mut = cbind(table(full_df_matrix_melt_LdKO_mut$psmed_peptide)))
head(LdKO_mut_tbl); dim(LdKO_mut_tbl)


LdKOs_tbl = full_join(LdKO_wt_tbl, LdKO_mut_tbl)
head(LdKOs_tbl); dim(LdKOs_tbl) # 1712
LdKOs_tbl[is.na(LdKOs_tbl)] = 0
table(LdKOs_tbl$Freq_wt, LdKOs_tbl$Freq_mut)

LdKOs_df = data.frame(rbind(table(LdKOs_tbl$Freq_wt, LdKOs_tbl$Freq_mut)))
colnames(LdKOs_df) = seq(0, ncol(LdKOs_df)-1)
LdKOs_df
LdKOs_df$WT = rownames(LdKOs_df)
LdKO_mdf <- gather(LdKOs_df, -WT, key="MUT", value="count")
LdKO_mdf
LdKO_mdf_copy = LdKO_mdf
LdKO_mdf$WT = factor(LdKO_mdf$WT, levels=rev(seq(0, max(LdKO_mdf$WT))))

pdf(file=paste0(path_Immuno, 'PeptideIntensity_Matrices_LdKO.pdf'), width = 8, height=8)
ggplot(LdKO_mdf, aes(x=MUT, y=WT, fill=log10(count), label=count))+
  geom_tile()+
  geom_text()+
  scale_fill_gradient(high="red2", low="white") + theme_classic() + ggtitle('LdKO') + scale_x_discrete(position = "top") + ylab('WT') + coord_fixed()
dev.off()
# LdKO



Ctrls_OnlyInMut = Ctrls_tbl[Ctrls_tbl$Freq_wt == 0 & Ctrls_tbl$Freq_mut == 3,]
dim(Ctrls_OnlyInMut)

Ctrls_OnlyInWt = Ctrls_tbl[Ctrls_tbl$Freq_wt == 3 & Ctrls_tbl$Freq_mut == 0,]
dim(Ctrls_OnlyInWt)

Ctrls_log$top = ''
Ctrls_log[Ctrls_log$psmed_peptide %in% Ctrls_OnlyInMut$psmed_peptide,]$top = 'Ctrls_OnlyInMut'
Ctrls_log[Ctrls_log$psmed_peptide %in% Ctrls_OnlyInWt$psmed_peptide,]$top = 'Ctrls_OnlyInWt'

Ctrls_log = Ctrls_log[order(Ctrls_log$top),]


LdKOs_OnlyInMut = LdKOs_tbl[LdKOs_tbl$Freq_wt == 0 & LdKOs_tbl$Freq_mut == 4,]
dim(LdKOs_OnlyInMut)

LdKOs_OnlyInWt = LdKOs_tbl[LdKOs_tbl$Freq_wt == 3 & LdKOs_tbl$Freq_mut == 0,]
dim(LdKOs_OnlyInWt)


LdKOs_log$top = ''
LdKOs_log[LdKOs_log$psmed_peptide %in% LdKOs_OnlyInMut$psmed_peptide,]$top = 'LdKOs_OnlyInMut'
LdKOs_log[LdKOs_log$psmed_peptide %in% LdKOs_OnlyInWt$psmed_peptide,]$top = 'LdKOs_OnlyInWt'

LdKOs_log = LdKOs_log[order(LdKOs_log$top),]





pdf(file=paste0(path_Immuno, 'PeptideIntensity_Ctrl_median.pdf'), width = 8, height=8)
ggplot(Ctrls_log) + geom_point(aes(Ctrl_Smg6wt_median, Ctrl_Smg6mut_median, col=top), size=0.8) + ggtitle('Ctrl - log transformed', subtitle = paste0(nrow(Ctrls_log), ' peptides')) + facet_grid(. ~ pepsource) + geom_abline(intercept = 0, slope = 1) +
  theme_bw() + theme(title = element_text(face ="bold")) + coord_fixed() + scale_color_manual(values=c('black', 'lightblue', 'pink'))
dev.off()


pdf(file=paste0(path_Immuno, 'PeptideIntensity_LdKO_median.pdf'), width = 8, height=8)
ggplot(LdKOs_log) + geom_point(aes(LdKO_Smg6wt_median, LdKO_Smg6mut_median, col=top), size=0.8) + ggtitle('LdKO - log transformed', subtitle = paste0(nrow(LdKOs_log), ' peptides')) + facet_grid(. ~ pepsource) + geom_abline(intercept = 0, slope = 1) +
  theme_bw() + theme(title = element_text(face ="bold")) + coord_fixed() + scale_color_manual(values=c('black', 'lightblue', 'pink'))
dev.off()






#
Comparisons = data.frame(Group='Ctrls_OnlyInWT', psmed_peptide=Ctrls_OnlyInWT$psmed_peptide)
Comparisons = rbind(Comparisons, data.frame(Group='Ctrls_OnlyInMut', psmed_peptide=Ctrls_OnlyInMut$psmed_peptide))
Comparisons = rbind(Comparisons, data.frame(Group='LdKOs_OnlyInWT', psmed_peptide=LdKOs_OnlyInWT$psmed_peptide))
Comparisons = rbind(Comparisons, data.frame(Group='LdKOs_OnlyInMut', psmed_peptide=LdKOs_OnlyInMut$psmed_peptide))

head(Comparisons); dim(Comparisons) # 63

Comparions_table = table(Comparisons)


names(which(colSums(Comparions_table) != 1))


write.table(Comparisons, file = paste0(path_Immuno, 'PeptidesConditionComparison.txt'), quote = F, sep = '\t', row.names = F, col.names = T)




Comparisons_info = left_join(Comparisons, df[c('psmed_peptide', 'pepsource', 'mapped_proteins', 'coordinates_in_proteins')])
Comparisons_info = unique(Comparisons_info)
head(Comparisons_info); dim(Comparisons_info)

write.table(Comparisons_info, file = paste0(path_Immuno, 'PeptidesConditionComparisonInfo.txt'), quote = F, sep = '\t', row.names = F, col.names = T)

#








# Pearson's Chi-squared Test for Count Data
pdf(file=paste0(path_Immuno, 'PeptideIntensity_ChiSquaredTest.pdf'), width = 8, height=8)

# Ctrl
Ctrl_chisq = chisq.test(table(Ctrls_tbl$Freq_wt, Ctrls_tbl$Freq_mut))
Ctrl_chisq
Ctrl_chisq$observed
round(Ctrl_chisq$expected,2)
round(Ctrl_chisq$residuals, 3)

corrplot(Ctrl_chisq$residuals, is.cor = FALSE, title=paste0("\n\nPearson's r (", Ctrl_chisq$parameter, ") = ", round(Ctrl_chisq$statistic, 2) , ", p = ", Ctrl_chisq$p.value))


contrib <- 100*Ctrl_chisq$residuals^2/Ctrl_chisq$statistic
round(contrib, 3)

corrplot(contrib, is.cor = FALSE)


# LdKO
LdKO_chisq = chisq.test(table(LdKOs_tbl$Freq_wt, LdKOs_tbl$Freq_mut))
LdKO_chisq
LdKO_chisq$observed
round(LdKO_chisq$expected,2)
round(LdKO_chisq$residuals, 3)

corrplot(LdKO_chisq$residuals, is.cor = FALSE, title=paste0("\n\nPearson's r (", LdKO_chisq$parameter, ") = ", round(LdKO_chisq$statistic, 2) , ", p = ", LdKO_chisq$p.value))


contrib <- 100*LdKO_chisq$residuals^2/LdKO_chisq$statistic
round(contrib, 3)

corrplot(contrib, is.cor = FALSE)
dev.off()