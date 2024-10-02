### QC plots - based on previous scripts available in the lab
### Virginie Ricci - October 2024

######################################################################
### Script that perform some QC on the new protocols used to generate
### ribo-seq libraries
######################################################################
library(ggplot2)
library(gridExtra)
library(cowplot)
library(DESeq2)

range01 <- function(x, minV, maxV){(x-minV)/(maxV-minV)}

######################################################################
### Variables to define
######################################################################
project='Enes'
path_home=''

Prefix='RiboSeq'
Run_nb=2371


path_plot=''
system(paste0('mkdir -p ', path_plot))

setwd(path_plot)

######################################################################
### Mapping stats
######################################################################

mappingStats_f= paste0(path_home, 'mapping_data/qc_', Prefix, '/', 'mappingStats.dat')
newMapStats <- read.table(mappingStats_f, header=T, sep="\t")
unique(newMapStats$Type)
newMapStats$Type <- factor(newMapStats$Type, levels=c('unmapped',
    'star-mouse-rrna', 
    'star-human-rrna', 
    'star-mouse-trna', , 
    'star-mouse-genome'))

newMapStats = newMapStats[!is.na(newMapStats$Type),]

maxs_1 = max(aggregate(newMapStats$Counts, by=list(Library=newMapStats$Library), FUN=sum)[,2])

p1 <- ggplot(newMapStats, aes(x=Library, y=Counts, fill=Type)) +
    geom_bar(stat="identity") +
    scale_fill_brewer(palette="Set1")


p2 <- ggplot(newMapStats, aes(x=Library, y=Counts, fill=Type)) +
    geom_bar(stat="identity", position="fill") +
    scale_fill_brewer(palette="Set1")

maxs = max(maxs_1, maxs_2)

p1 = p1 + ylim(0, maxs)

ggsave(paste0(path_plot, 'mappingStats.pdf'), p1, width=6)
ggsave(paste0(path_plot, 'mappingStatsPercentage.pdf'), p2, width=6)


######################################################################
### Check UMI filtering and de-duplication effects on samples total
### reads mapped to cDNA

## Global barcodes numbers:
barcodeCounts_f= paste0(path_home, 'mapping_data/qc_', Prefix, '/', 'globalBarcodeCounts.dat')
barcodeCounts <- read.table(barcodeCounts_f, header=FALSE, stringsAsFactors=TRUE)

colnames(barcodeCounts) <- c('Library', 'Barcode', 'Counts')
lnames <- levels(barcodeCounts[,1])
levels(barcodeCounts[,1]) <- unlist(strsplit(lnames, '_'))[seq(1, length(lnames)*6, by=6)]
gtot <- tapply(barcodeCounts$Counts, barcodeCounts$Library, sum)
barcodeCounts$Fraction <- barcodeCounts$Counts / gtot[barcodeCounts$Library]
head(barcodeCounts)

ssp <- ggplot(data=barcodeCounts, aes(x=Library, y=Counts, fill=Barcode)) +
    geom_bar(stat="identity") +
    ggtitle('Barcode splitting - all reads') +
    scale_fill_brewer(palette="Dark2") +
    theme_minimal()
ggsave(paste0(path_plot, 'globalSamplesSizes.pdf'), ssp, width=6)

ssp <- ggplot(data=barcodeCounts, aes(x=Library, y=Counts, fill=Barcode)) +
    geom_bar(stat="identity", position="fill") +
    ggtitle('Barcode splitting - percentage - all reads') +
    scale_fill_brewer(palette="Dark2") +
    theme_minimal()
ggsave(paste0(path_plot, 'globalSamplesSizesPercentage.pdf'), ssp, width=6)
