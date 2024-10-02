#! R/4.4.0
### Read count
### Virginie Ricci - October 2024

######################################################################
print(version)

library(DEXSeq)
library(DESeq2)
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
library(BSgenome)
library(ribosomeProfilingQC)


######################################################################
# ADAPT THE CODE TO YOUR DATASET
# ADAPT THE CODE TO YOUR DATASET
# ADAPT THE CODE TO YOUR DATASET
######################################################################



######################################################################
### Variables to define

path_home=''
SCRIPT=paste0(path_home, 'scripts/')
DATA=paste0(path_home, 'data/')
RAW_DATA=paste0(path_home, 'raw_data/')
MAPPING=paste0(path_home, 'mapping_data/')
Results_path=paste0(path_home, 'results/')
DB=''
GENOME='Mmusculus.GRCm39.111'
GENOME_DB=paste0(DB, 'star/', GENOME, '/')
DB_GTF=paste0(DB, 'gtf/')
GTF_GENOME=paste0(DB_GTF, GENOME, '.sorted.gtf')

setwd(path_home)

Prefix='RNASeq'
Run_nb='2327'


system(paste0("mkdir -p ", Results_path))

setwd(Results_path)
######################################################################



######################################################################
### Functions
EnhancedVolcano_plot_fun = function(res, title_volcano, FC_cutoff=1, FDR_cutoff=0.05){
  keyvals <- ifelse(
    res$log2FoldChange < -log2(FC_cutoff), 'royalblue',
    ifelse(res$log2FoldChange > log2(FC_cutoff), 'gold',
           'black'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'gold'] <- 'upregulated'
  names(keyvals)[keyvals == 'black'] <- 'ns'
  names(keyvals)[keyvals == 'royalblue'] <- 'downregulated'
  
  volcano = EnhancedVolcano::EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='padj', pCutoff = FDR_cutoff, FCcutoff = FC_cutoff,
                                             legendPosition = 'bottom', legendLabSize = 12, legendIconSize = 4.0, drawConnectors = TRUE, widthConnectors = 0.75, labSize = 2,
                                             title = title_volcano, subtitle = paste0('FC cutoff: ', FC_cutoff, ' and FDR cutoff: ', FDR_cutoff), colCustom = keyvals) +
    theme(plot.tag.position = 'bottom') +
    labs(tag = paste0(length(names(keyvals)[keyvals == 'royalblue']), ' down, ', length(names(keyvals)[keyvals == 'black']), ' ns, ', length(names(keyvals)[keyvals == 'gold']), ' up'))
  # xlim = c(-xlims, xlims), 
  return(volcano)
}

EnhancedVolcano_plot_fun_rescale = function(res, title_volcano, FC_cutoff=1, FDR_cutoff=0.05, xlims=50){
  keyvals <- ifelse(
    res$log2FoldChange < -log2(FC_cutoff), 'royalblue',
    ifelse(res$log2FoldChange > log2(FC_cutoff), 'gold',
           'black'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'gold'] <- 'upregulated'
  names(keyvals)[keyvals == 'black'] <- 'ns'
  names(keyvals)[keyvals == 'royalblue'] <- 'downregulated'
  
  volcano = EnhancedVolcano::EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='padj', pCutoff = FDR_cutoff, FCcutoff = FC_cutoff,
                                             legendPosition = 'bottom', legendLabSize = 12, legendIconSize = 4.0, drawConnectors = TRUE, widthConnectors = 0.75, labSize = 2,
                                             title = title_volcano, subtitle = paste0('FC cutoff: ', FC_cutoff, ' and FDR cutoff: ', FDR_cutoff), colCustom = keyvals,
                                             xlim = c(-xlims, xlims)) +
    theme(plot.tag.position = 'bottom') +
    labs(tag = paste0(length(names(keyvals)[keyvals == 'royalblue']), ' down, ', length(names(keyvals)[keyvals == 'black']), ' ns, ', length(names(keyvals)[keyvals == 'gold']), ' up'))
  # xlim = c(-xlims, xlims), 
  return(volcano)
}

MA_plot = function(res, title_MA, tops=10, padj=0.05){
    res$ID = res$mgi_symbol

    res_signif = res[which(res$padj < 0.05),]

    MA = ggplot(data = res) +
    geom_point(
      aes(x = baseMean, y = log2FoldChange),
      size = 1.25,
      color = ifelse(is.na(res$padj) | res$padj >= 0.05, "grey50", "red"),
      alpha = 0.5) +
    geom_text_repel(data = res_signif[order(res_signif$pvalue),][c(1:tops),], aes(x = baseMean, y = log2FoldChange, label=ID), max.overlaps=Inf, force=2) +
    geom_hline(yintercept = 0) +
    scale_x_log10()  + theme_classic() + ggtitle(title_MA) +
    geom_smooth(data=res, aes(x = baseMean, y = log2FoldChange), method='loess', col='darkgreen', fill='darkgreen', se = F) +
    labs(subtitle=paste0('FDR cutoff: ', padj))

    return(MA)

}
######################################################################


######################################################################
# Input files

SampleBarcode = read.csv(paste0(RAW_DATA, 'Sample_Barcode.txt'), sep='\t', header = F)
names(SampleBarcode) = c('Sample', 'Barcode')
head(SampleBarcode); dim(SampleBarcode)
SampleBarcode$Library = ''
SampleBarcode$Background = ''
SampleBarcode$Genotype = '' 
SampleBarcode$name = ''

# .bam files
bams = list.files(path='', pattern='.bam$', full.names=T)
Library = ''
BC = ''
bams_info = data.frame(Library = Library, Barcode=BC, bam = bams)

sampleData = left_join(SampleBarcode, bams_info)
sampleData$RefSeq = 'mouse_genome'
######################################################################



######################################################################
# Reference genome annotation

# GTF files
GTF = rtracklayer::import(GTF_GENOME)
GTF_split = split(GTF, mcols(GTF)$gene_id)
GTF_exons = GTF[GTF$type == 'exon', ]

# Mmusculus.GRCm39.111.sorted.gtf
txdb = txdbmaker::makeTxDbFromGFF(file=GTF_GENOME, format='gtf') # 149076 transcripts
######################################################################


######################################################################
# Read count
print('Counts of reads (only filtered on rRNA and tRNA) mapped on genome ')
print('Each annotated element is a feature')
RNA_Counts_ALL <- GenomicAlignments::summarizeOverlaps(GTF_split, FILES$bam, mode="Union")
colnames(RNA_Counts_ALL) = FILES$Sample
colData(RNA_Counts_ALL)$Sample = FILES$Sample
colData(RNA_Counts_ALL)$Background = FILES$Background
colData(RNA_Counts_ALL)$Genotype = FILES$Genotype
colData(RNA_Counts_ALL)$name = FILES$name
colSums(assay(RNA_Counts_ALL))

RNA_Counts_ALL_df = as.data.frame(assay(RNA_Counts_ALL))
head(RNA_Counts_ALL_df); dim(RNA_Counts_ALL_df)

write.table(RNA_Counts_ALL_df, file = paste0(Results_path, 'NewStrategy_ReadCount.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
######################################################################



######################################################################
# Differential analysis

design_usage = formula( ~ Sample + exon + name:exon )
design = formula( ~ name)

RNA_Counts_ALL; dim(RNA_Counts_ALL)

RNA_Counts_ALL$condition <- factor(RNA_Counts_ALL$name, levels = unique(colData(RNA_Counts_ALL)$name))

listes=list()
listes$count = RNA_Counts_ALL

countings='RNA_Counts_ALL'

for (groups in c('Ctrl', 'LdKO', 'HCC')){
    print(groups)
    Selected_RNA_Counts_ALL = RNA_Counts_ALL[,grep(groups, colnames(RNA_Counts_ALL))]

    DESEQ = DESeqDataSet(se=Selected_RNA_Counts_ALL, design=design)

    smallestGroupSize <- 3
    keep <- rowSums(counts(DESEQ) >= 5) >= smallestGroupSize
    DESEQ <- DESEQ[keep,]

    DESEQ = DESeq(DESEQ, BPPARAM=BPPARAM)
    DESEQ_res <- results(DESEQ, contrast=c('name', levels(colData(DESEQ)$name)[1], levels(colData(DESEQ)$name)[2]), BPPARAM=BPPARAM)
    # resultsNames(DESEQ)[2] >> "name" "Ctrl_Smg6mut" "Ctrl_Smg6wt" 

    DESEQ_res_LFC <- lfcShrink(DESEQ, coef=2, type="ashr", BPPARAM=BPPARAM)

    DESEQ_counts = counts(DESEQ)
    DESEQ_vst = vst(DESEQ_counts)

    listes$DESEQ = DESEQ
    listes$DESEQ_res = DESEQ_res
    listes$DESEQ_res_LFC = DESEQ_res_LFC
    listes$DESEQ_counts = DESEQ_counts
    listes$DESEQ_vst = DESEQ_vst
    list_n = paste0('DESeq_', groups, '_', countings)
    assign(list_n, listes)
}





mart <- useEnsembl('ensembl', dataset = 'mmusculus_gene_ensembl', version = '111')
mart_all_info = getBM(filters = 'ensembl_gene_id', attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'gene_biotype', 'mgi_symbol', 'entrezgene_id'), values = rownames(assay(RNA_Counts_ALL)), mart = mart)
head(mart_all_info); dim(mart_all_info)
sub_mart_all_info = unique(subset(mart_all_info, select=c(ensembl_gene_id, gene_biotype, mgi_symbol)))
head(sub_mart_all_info); dim(sub_mart_all_info)

RNA_counts_all_info = left_join(RNA_counts_all, sub_mart_all_info)
head(RNA_counts_all_info); dim(RNA_counts_all_info)





DESeq_Ctrl_ALL = DESeq_Ctrl_RNA_Counts_ALL$DESEQ_res[!is.na(DESeq_Ctrl_RNA_Counts_ALL$DESEQ_res$log2FoldChange),]
head(DESeq_Ctrl_ALL); dim(DESeq_Ctrl_ALL)
str(DESeq_Ctrl_ALL)

DESeq_LdKO_ALL = DESeq_LdKO_RNA_Counts_ALL$DESEQ_res[!is.na(DESeq_LdKO_RNA_Counts_ALL$DESEQ_res$log2FoldChange),]
head(DESeq_LdKO_ALL); dim(DESeq_LdKO_ALL)
str(DESeq_LdKO_ALL)

DESeq_HCC_ALL = DESeq_HCC_RNA_Counts_ALL$DESEQ_res[!is.na(DESeq_HCC_RNA_Counts_ALL$DESEQ_res$log2FoldChange),]
head(DESeq_HCC_ALL); dim(DESeq_HCC_ALL)
str(DESeq_HCC_ALL)


head(sub_mart_all_info)
DS_Ctrl_ALL = data.frame(ensembl_gene_id=rownames(DESeq_Ctrl_ALL)); head(DS_Ctrl_ALL); dim(DS_Ctrl_ALL)
DS_Ctrl_ALL_final = left_join(DS_Ctrl_ALL, sub_mart_all_info); head(DS_Ctrl_ALL_final); dim(DS_Ctrl_ALL_final)
rownames(DESeq_Ctrl_ALL) = DS_Ctrl_ALL_final$mgi_symbol

DS_LdKO_ALL = data.frame(ensembl_gene_id=rownames(DESeq_LdKO_ALL)); head(DS_LdKO_ALL); dim(DS_LdKO_ALL)
DS_LdKO_ALL_final = left_join(DS_LdKO_ALL, sub_mart_all_info); head(DS_LdKO_ALL_final); dim(DS_LdKO_ALL_final)
rownames(DESeq_LdKO_ALL) = DS_LdKO_ALL_final$mgi_symbol

DS_HCC_ALL = data.frame(ensembl_gene_id=rownames(DESeq_HCC_ALL)); head(DS_HCC_ALL); dim(DS_HCC_ALL)
DS_HCC_ALL_final = left_join(DS_HCC_ALL, sub_mart_all_info); head(DS_HCC_ALL_final); dim(DS_HCC_ALL_final)
rownames(DESeq_HCC_ALL) = DS_HCC_ALL_final$mgi_symbol


# plots
pdf(file=paste0(Results_path, 'VolcanoPlots.pdf'), width = 8, height=8)
EnhancedVolcano_plot_fun(DESeq_Ctrl_ALL, 'DESeq: Control')
EnhancedVolcano_plot_fun(DESeq_LdKO_ALL, 'DESeq: LdKO')
EnhancedVolcano_plot_fun(DESeq_HCC_ALL, 'DESeq: HCC')
dev.off()


pdf(file=paste0(Results_path, 'MAPlots.pdf'), width = 8, height=8)
MA_plot(DS_Ctrl_ALL_final, 'DESeq: Control')
MA_plot(DS_LdKO_ALL_final, 'DESeq: LdKO_')
MA_plot(DS_HCC_ALL_final, 'DESeq: HCC')
dev.off()

######################################################################
