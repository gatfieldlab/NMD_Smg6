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

Prefix='RiboSeq'
Run_nb='2371'


system(paste0("mkdir -p ", Results_path))

setwd(Results_path)
######################################################################



######################################################################
### Functions
EnhancedVolcano_plot_fun_rescale = function(res, title_volcano, FC_cutoff=1, FDR_cutoff=0.05, xlims=10){
  keyvals <- ifelse(
    res$log2FoldChange < -FC_cutoff, 'royalblue',
    ifelse(res$log2FoldChange > FC_cutoff, 'gold',
           'black'))
  
  keyvals[which(res$padj > FDR_cutoff)] = 'black'

  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'gold'] <- 'upregulated'
  names(keyvals)[keyvals == 'black'] <- 'ns'
  names(keyvals)[keyvals == 'royalblue'] <- 'downregulated'

  res$labels = paste0(res$mgi_symbol , ':', res$ID)

  upreg = res[which(names(keyvals) == 'upregulated'),]
  upreg_labels1 = upreg[order(upreg$padj), ][c(1:20),]
  upreg_labels2 = upreg[order(upreg$log2FoldChange, decreasing=TRUE), ]
  upreg_labels2 = upreg_labels2[!is.na(upreg_labels2$padj),][c(1:20),]
  upreg_labels = unique(c(rownames(upreg_labels1), rownames(upreg_labels2)))

  downreg = res[which(names(keyvals) == 'downregulated'),]
  downreg_labels1 = downreg[order(downreg$padj), ][c(1:20),]
  downreg_labels2 = downreg[order(downreg$log2FoldChange, decreasing=FALSE), ]
  downreg_labels2 = downreg_labels2[!is.na(downreg_labels2$padj),][c(1:20),]
  downreg_labels = unique(c(rownames(downreg_labels1), rownames(downreg_labels2)))

  #labels = c(rownames(upreg_labels), rownames(downreg_labels))
  labels = c(upreg[upreg_labels,]$labels, downreg[downreg_labels,]$labels)


  
  volcano = EnhancedVolcano::EnhancedVolcano(res, lab=res$labels, selectLab = labels, x='log2FoldChange', y='padj', pCutoff = FDR_cutoff, FCcutoff = FC_cutoff,
                                             legendPosition = 'bottom', legendLabSize = 12, legendIconSize = 4.0, drawConnectors = TRUE, widthConnectors = 0.75, labSize = 2,
                                             title = title_volcano, subtitle = paste0('FC cutoff: ', FC_cutoff, ' and FDR cutoff: ', FDR_cutoff), colCustom = keyvals,
                                             xlim = c(-xlims, xlims)) +
    theme(plot.tag.position = 'bottom') +
    labs(tag = paste0(length(names(keyvals)[keyvals == 'royalblue']), ' down, ', length(names(keyvals)[keyvals == 'black']), ' ns, ', length(names(keyvals)[keyvals == 'gold']), ' up'))
  # xlim = c(-xlims, xlims), 
  return(volcano)
}



MA_plot = function(res, res_signif, title_MA, padj = 0.05){
  RES = data.frame(res)
  MA = ggplot(data = RES) +
    geom_point(
      aes(x = baseMean, y = log2FoldChange), #col=padj),
      size = 1.25,
      color = ifelse(is.na(RES$padj) | RES$padj >= padj, "grey50", "red"),
      alpha = 0.5,
      #shape = ifelse(abs(ds2data$log2FoldChange) >= 7.5, 17, 16)
    ) +
    #geom_text_repel(data = res_signif[order(res_signif$pvalue),][c(1:50),], aes(x = baseMean, y = log2FoldChange, label=ID), max.overlaps=Inf, force=2) +
    geom_hline(yintercept = 0) +
    scale_x_log10() + theme_classic() + ggtitle(title_MA) +
    geom_smooth(data=RES, aes(x = baseMean, y = log2FoldChange), method='loess', col='darkgreen', fill='darkgreen', se = F) +
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

#installed.genomes()
#available.genomes()[grep('Mmusculus', available.genomes())]

# if (!require("BiocManager"))
#         install.packages("BiocManager")
#     BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")

library('BSgenome.Mmusculus.UCSC.mm39')
genome = Mmusculus

CDS = prepareCDS(txdb)
CDS
######################################################################



######################################################################
# Get P site and then A site

yieldSize <- 10000000

sampleData$Asite_bam = str_replace(sampleData$bam, '.bam', '.Asite.bam')
for (bamFile in sampleData$bam){
  print(bamFile)

  files = sampleData[sampleData$bam == bamFile,]
  samples = files$Sample
  new_bamFile = files$Asite_bam
  print(samples)

  if (file.exists(new_bamFile) == FALSE){
    reads = GenomicAlignments::readGAlignments(bamFile, use.names = TRUE)

    ### part of getPsiteCoordinates() from ribosomeProfilingQC
    anchor='5end'

    BAMfile = BamFile(bamFile, yieldSize = yieldSize)
    shift = estimatePsite(BAMfile, CDS, genome)

    x = reads

    x <- x[qwidth(x)>shift & width(x)>shift & 
           cigarWidthAlongReferenceSpace(cigar(x), 
                                         N.regions.removed = TRUE)>shift]
    if(shift==0){
      return(x)
    }
    strds <- as.character(strand(x)) == "-"
    cigars <- cigar(x)
    cigars <- as.character(cigarNarrow(cigars))
    if(anchor=="5end"){
      cigars <- cigarQNarrow(cigars,
                            start=ifelse(strds, 1, shift+1),
                            end=ifelse(strds, -shift-1, -1))
    }else{
      l <- mcols(x)$qwidth
      shift <- l - shift
      cigars <- cigarQNarrow(cigars,
                            start=ifelse(strds, 1, shift),
                            end=ifelse(strds, -shift, -1))
    }

    x@cigar <- as.character(cigars)
    x@start <- x@start + attributes(cigars)$rshift
    x
    gr = GRanges(x)
    end(gr) = start(gr)
    gr
    ### part of getPsiteCoordinates()

    export(gr, BamFile(new_bamFile))
  }
}
Bam_Repeats_Asite = BamFileList(sampleData$Asite_bam)
names(Bam_Repeats_Asite) = sampleData$Sample
######################################################################



######################################################################
# Create counting bins >> bins can be assigned to several genes here (linked.to.single.gene.only=FALSE)

flattenedAnnot_Exons_all = GenomicFeatures::exonicParts(txdb, linked.to.single.gene.only=FALSE)
flattenedAnnot_Exons_all$tx_name = do.call(c, lapply(flattenedAnnot_Exons_all$tx_name, function(x){paste0(x, collapse='_')}))
flattenedAnnot_Exons_all$gene_id = do.call(c, lapply(flattenedAnnot_Exons_all$gene_id, function(x){paste0(x, collapse='_')}))
names(flattenedAnnot_Exons_all) = flattenedAnnot_Exons_all$gene_id
flattenedAnnot_Exons_all


# Get genes start and end coordinates
GTF_genes = unlist(transcriptsBy(txdb, 'gene'))
GTF_genes$gene_id = names(GTF_genes)
GTF_genes

GTF_genes_reduced = data.frame(GTF_genes) %>% 
  group_by(gene_id) %>% 
  summarise(  
    seqnames = unique(seqnames), 
    start = min(start),
    end   = max(end),
    strand = unique(strand), 
    gene_id = unique(gene_id))

GTF_genes_reduced = data.frame(GTF_genes_reduced)
GTF_genes_reduced = GRanges(GTF_genes_reduced)
names(GTF_genes_reduced) = GTF_genes_reduced$gene_id
GTF_genes_reduced

# Get transcripts start and end coordinates per genes
GTF_transcripts_reduced = data.frame(GTF_genes) %>% 
  group_by(tx_name) %>% 
  summarise(  
    seqnames = unique(seqnames), 
    start = min(start),
    end   = max(end),
    strand = unique(strand), 
    gene_id = unique(gene_id),
    tx_name = unique(tx_name))

GTF_transcripts_reduced = data.frame(GTF_transcripts_reduced)
GTF_transcripts_reduced = GRanges(GTF_transcripts_reduced)
GTF_transcripts_reduced


# Apply reduced() >> reduced ranges for each distinct (seqname, strand) pairing
# >> combined bins that are sequential/overlapping
flattenedAnnot_Exons_reduced = unlist(reduce(split(flattenedAnnot_Exons_all, elementMetadata(flattenedAnnot_Exons_all)$gene_id), with.revmap=TRUE))
flattenedAnnot_Exons_reduced
flattenedAnnot_Exons_reduced$type = 'exon'
# with.revmap >> Should the mapping from output to input ranges be stored in the returned object? If yes, then it is stored as metadata column revmap of type IntegerList.


# Get the 'reverse image' of each range
flattenedAnnot_Exons_reduced_reverse = gaps(flattenedAnnot_Exons_reduced)
flattenedAnnot_Exons_reduced_reverse$type = 'other'
flattenedAnnot_Exons_reduced_reverse

# Check which 'other' regions are within a gene
overlaps_reduced = findOverlaps(GTF_genes_reduced, flattenedAnnot_Exons_reduced_reverse)
overlaps_reduced # query, subject


flattenedAnnot_Exons_reduced_reverse_introns = flattenedAnnot_Exons_reduced_reverse[subjectHits(overlaps_reduced)]
flattenedAnnot_Exons_reduced_reverse_introns
flattenedAnnot_Exons_reduced_reverse_introns$type = 'intron'
names(flattenedAnnot_Exons_reduced_reverse_introns) = names(GTF_genes_reduced[queryHits(overlaps_reduced),])
flattenedAnnot_Exons_reduced_reverse_introns


flattenedAnnot_Exons_reduced_reverse_introns_unique = unique(flattenedAnnot_Exons_reduced_reverse_introns)
flattenedAnnot_Exons_reduced_reverse_introns_unique

# Combined exonic bins and intronic bins
flattenedAnnot_Exons_reduced_all = bindROWS(flattenedAnnot_Exons_reduced, objects=list(flattenedAnnot_Exons_reduced_reverse_introns_unique))

flattenedAnnot_Exons_reduced_all = sort(flattenedAnnot_Exons_reduced_all)
flattenedAnnot_Exons_reduced_all
flattenedAnnot_Exons_reduced_all$revmap = NULL
######################################################################



######################################################################
# A-site count in bins

print('Read count in reduced bins using Asite')
Counts_flat_reduced <- GenomicAlignments::summarizeOverlaps(flattenedAnnot_Exons_reduced_all, Bam_Repeats_Asite, mode="Union")
colData(Counts_flat_reduced)$Sample = sampleData_mini$Sample
colData(Counts_flat_reduced)$Background = sampleData_mini$Background
colData(Counts_flat_reduced)$Genotype = sampleData_mini$Genotype
colData(Counts_flat_reduced)$name = sampleData_mini$name
colSums(assay(Counts_flat_reduced))
assay(Counts_flat_reduced)[c(1:10),]

flattenedAnnot_Exons_reduced_all_df = data.frame(flattenedAnnot_Exons_reduced_all)
head(flattenedAnnot_Exons_reduced_all_df); dim(flattenedAnnot_Exons_reduced_all_df)
flattenedAnnot_Exons_reduced_all_df$id = paste0(flattenedAnnot_Exons_reduced_all_df$seqnames, ':', 
  flattenedAnnot_Exons_reduced_all_df$start, '-', flattenedAnnot_Exons_reduced_all_df$end, 
  '(', flattenedAnnot_Exons_reduced_all_df$strand, ')', ' ', flattenedAnnot_Exons_reduced_all_df$type)
head(flattenedAnnot_Exons_reduced_all_df)

names(Counts_flat_reduced) = flattenedAnnot_Exons_reduced_all_df$id
assay(Counts_flat_reduced)[c(1:10),]

Counts_flat_reduced_df = as.data.frame(assay(Counts_flat_reduced))
head(Counts_flat_reduced_df); dim(Counts_flat_reduced_df)


write.table(Counts_flat_reduced_df, file = paste0(Results_path, 'NewStrategy_ReadCount.txt'), quote = F, sep = '\t', row.names = T, col.names = T)
write.table(flattenedAnnot_Exons_reduced_all_df, file = paste0(Results_path, 'NewStrategy_Bins.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
######################################################################





######################################################################

DESeq_res_df = data.frame(assay(Counts_flat_reduced))
head(DESeq_res_df); dim(DESeq_res_df) # 535663

length(grep('exon', rownames(DESeq_res_df))) # 295535 exonic regions
length(grep('intron', rownames(DESeq_res_df))) # 240128 intronic regions


smallestGroupSize <- 3
keep <- rowSums(DESeq_res_df >= 5) >= smallestGroupSize
DESeq_res_df <- DESeq_res_df[keep,]
head(DESeq_res_df); dim(DESeq_res_df) # 119866

length(grep('exon', rownames(DESeq_res_df))) # 114600 exonic regions
length(grep('intron', rownames(DESeq_res_df))) # 5266 intronic regions


ctrl_wt = grep('Ctrl_Smg6wt', names(DESeq_res_df))
ctrl_mut = grep('Ctrl_Smg6mut', names(DESeq_res_df))
ldko_wt = grep('LdKO_Smg6wt', names(DESeq_res_df))
ldko_mut = grep('LdKO_Smg6mut', names(DESeq_res_df))



# MEDIAN # MEDIAN # MEDIAN
DESeq_res_df_median = data.frame(Ctrl_Smg6wt=rowMedians(as.matrix(DESeq_res_df[ctrl_wt])), Ctrl_Smg6mut=rowMedians(as.matrix(DESeq_res_df[ctrl_mut])),
                               LdKO_Smg6wt=rowMedians(as.matrix(DESeq_res_df[ldko_wt])), LdKO_Smg6mut=rowMedians(as.matrix(DESeq_res_df[ldko_mut])))

head(DESeq_res_df_median); dim(DESeq_res_df_median)
DESeq_res_df_median$ID=rownames(DESeq_res_df_median)

seqnames = str_split_fixed(rownames(DESeq_res_df_median), ':', 2)[,1]
coo = str_split_fixed(str_split_fixed(rownames(DESeq_res_df_median), ':', 2)[,2], '[(]', 2)[,1]
start = str_split_fixed(coo, '-', 2)[,1]
end = str_split_fixed(coo, '-', 2)[,2]
strand_types = str_split_fixed(str_split_fixed(rownames(DESeq_res_df_median), ':', 2)[,2], '[(]', 2)[,2]
strand = str_split_fixed(strand_types, '[)]', 2)[,1]
types = str_split_fixed(strand_types, '[)]', 2)[,2]
types = str_replace(types, ' ', '')


DESeq_res_df_median$seqnames = seqnames
DESeq_res_df_median$start = start
DESeq_res_df_median$end = end
DESeq_res_df_median$strand = strand
DESeq_res_df_median$type = types
head(DESeq_res_df_median); dim(DESeq_res_df_median)


DESeq_res_df_median_info_filtered_df = DESeq_res_df_median[DESeq_res_df_median$LdKO_Smg6mut > 0, ] # DESeq_res_df_median_info$Ctrl_Smg6mut_norm > 0 |
head(DESeq_res_df_median_info_filtered_df); dim(DESeq_res_df_median_info_filtered_df) # 119728

length(grep('exon', rownames(DESeq_res_df_median_info_filtered_df))) # 114523 exonic regions
length(grep('intron', rownames(DESeq_res_df_median_info_filtered_df))) # 5205 intronic regions


DESeq_res_df_median_info_filtered_df$width = as.integer(DESeq_res_df_median_info_filtered_df$end) - as.integer(DESeq_res_df_median_info_filtered_df$start) + 1

DESeq_res_df_median_info_filtered_df = DESeq_res_df_median_info_filtered_df[DESeq_res_df_median_info_filtered_df$width >=5, ] # to have sequences of the 3 reading frames

DESeq_res_df_median_info_filtered_df$score = '0'

DESeq_res_df_median_info_filtered_bed = DESeq_res_df_median_info_filtered_df[c('seqnames', 'start', 'end', 'ID', 'score', 'strand')]
head(DESeq_res_df_median_info_filtered_bed); dim(DESeq_res_df_median_info_filtered_bed)

DESeq_res_df_median_info_filtered_bed$ID = sprintf("%s%0.3d", DESeq_res_df_median_info_filtered_df$type, seq(1, nrow(DESeq_res_df_median_info_filtered_bed)))


head(DESeq_res_df_median_info_filtered_df); dim(DESeq_res_df_median_info_filtered_df) # 119724
DESeq_res_df_median_info_filtered_df$ID = DESeq_res_df_median_info_filtered_bed$ID

# bedtools getfasta is 0-based
# substract 1 to start position
DESeq_res_df_median_info_filtered_bed$start = as.numeric(DESeq_res_df_median_info_filtered_bed$start)
DESeq_res_df_median_info_filtered_bed$start = DESeq_res_df_median_info_filtered_bed$start-1
head(DESeq_res_df_median_info_filtered_bed); dim(DESeq_res_df_median_info_filtered_bed)


write.table(DESeq_res_df_median_info_filtered_df, file = paste0(Results_path, 'DEXseq_ExonsIntrons_NewStrategy.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
write.table(DESeq_res_df_median_info_filtered_bed, file = paste0(Results_path, 'DEXseq_ExonsIntrons_NewStrategy.bed'), quote = F, sep = '\t', row.names = F, col.names = F)


# Correct fasta headers (+1 to start position because bedtools getfasta is 0-based)
headers = DESeq_res_df_median_info_filtered_bed

headers$correct = paste0(headers$ID, '::', headers$seqnames, ':', headers$start + 1, '-', headers$end, '(', headers$strand, ')')
headers$wrong = paste0(headers$ID, '::', headers$seqnames, ':', headers$start, '-', headers$end, '(', headers$strand, ')')

head(headers)
headers$seqnames = NULL
headers$start = NULL
headers$end = NULL
headers$ID = NULL
headers$strand = NULL
headers$score = NULL
head(headers); dim(headers) # 119724

write.table(headers, file = paste0(Results_path, 'DEXseq_ExonsIntrons_NewStrategy_headers.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
######################################################################





######################################################################
# Get gene info

DESeq_res_df_median_info_filtered_gr = GRanges(DESeq_res_df_median_info_filtered_df)
DESeq_res_df_median_info_filtered_gr; length(DESeq_res_df_median_info_filtered_gr) # 119724

overlaps_final = findOverlaps(GTF_transcripts_reduced, DESeq_res_df_median_info_filtered_gr)
overlaps_final

DESeq_res_df_median_info_filtered_gr$ensembl_gene_id = ''
DESeq_res_df_median_info_filtered_gr$ensembl_transcript_id = ''

DESeq_res_df_median_info_filtered_gr[subjectHits(overlaps_final)]$ensembl_gene_id <- GTF_transcripts_reduced[queryHits(overlaps_final),]$gene_id
DESeq_res_df_median_info_filtered_gr[subjectHits(overlaps_final)]$ensembl_transcript_id <- GTF_transcripts_reduced[queryHits(overlaps_final),]$tx_name
DESeq_res_df_median_info_filtered_gr

DESeq_res_df_median_info_filtered = data.frame(DESeq_res_df_median_info_filtered_gr)


mart <- useEnsembl('ensembl', dataset = 'mmusculus_gene_ensembl', version = '111')
mart_deseq = getBM(filters = 'ensembl_gene_id', attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'gene_biotype', 'mgi_symbol', 'entrezgene_id', 'ensembl_peptide_id'), values = unique(DESeq_res_df_median_info_filtered$ensembl_gene_id), mart = mart)
head(mart_deseq); dim(mart_deseq)
mart_deseq = unique(subset(mart_deseq, select=c(ensembl_gene_id, gene_biotype, ensembl_transcript_id, ensembl_peptide_id, mgi_symbol)))
head(mart_deseq); dim(mart_deseq)

DESeq_res_df_median_info_filtered_final = left_join(DESeq_res_df_median_info_filtered, mart_deseq)
head(DESeq_res_df_median_info_filtered_final); dim(DESeq_res_df_median_info_filtered_final)
DESeq_res_df_median_info_filtered_final$score = NULL
DESeq_res_df_median_info_filtered_final$Full_ID = rownames(DESeq_res_df_median_info_filtered_df)
head(DESeq_res_df_median_info_filtered_final); dim(DESeq_res_df_median_info_filtered_final)


write.table(DESeq_res_df_median_info_filtered_final, file = paste0(Results_path, 'DEXseq_ExonsIntrons_NewStrategy_final.txt'), quote = F, sep = '\t', row.names = F, col.names = T)
######################################################################








######################################################################
# Differential analysis

Counts_flat_reduced; dim(Counts_flat_reduced)

dds = DESeqDataSetFromMatrix(countData = as.matrix(assay(Counts_flat_reduced)),
                              colData = colData(Counts_flat_reduced),
                              design= formula( ~ name))
dds

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 5) >= smallestGroupSize
dds <- dds[keep,]
head(counts(dds)); dim(dds)

dds = dds[DESeq_res_df_median_info_filtered_final$Full_ID, ] # keep bins with LdKO_Smg6mut > 0
dds; dim(dds)


# DESeq
design = formula( ~ name)


# Ctrl
DESEQ = dds[,grep('Ctrl', colnames(dds))]
DESEQ = DESeqDataSet(DESEQ, design=design)
DESEQ = DESeq(DESEQ)
DESEQ_res = results(DESEQ, contrast=c('name', levels(colData(DESEQ)$name)[grep('wt', levels(colData(DESEQ)$name))], levels(colData(DESEQ)$name)[grep('mut', levels(colData(DESEQ)$name))]))
DESEQ_res

DESEQ_Ctrl = DESEQ
DESEQ_Ctrl_res = DESEQ_res


DESEQ_Ctrl_res$Full_ID = rownames(DESEQ_Ctrl_res)
DESEQ_Ctrl_res_info = left_join(data.frame(DESEQ_Ctrl_res), DESeq_res_df_median_info_filtered_final)
head(DESEQ_Ctrl_res_info); dim(DESEQ_Ctrl_res_info)

first_row = DESEQ_Ctrl_res_info[1,]; first_row

if (nrow(first_row[(first_row$log2FoldChange > 0) & first_row$Ctrl_Smg6mut > first_row$Ctrl_Smg6wt, ]) == 0){
  if (nrow(first_row[(first_row$log2FoldChange < 0) & first_row$Ctrl_Smg6mut < first_row$Ctrl_Smg6wt, ]) == 0){
    print('Problem! log2FoldChange is positive and Ctrl_Smg6mut expression is lower than Ctrl_Smg6_wt expression or vice versa')
    print('Ctrl_Smg6mut expression should be bigger than Ctrl_Smg6_wt expression or vice versa')
  }
}


DESEQ_Ctrl_res_signif = DESEQ_Ctrl_res[which(DESEQ_Ctrl_res$padj < 0.05),]
DESEQ_Ctrl_res_signif; dim(DESEQ_Ctrl_res_signif)

res_Ctrl = DESeq_res_df_median_info_filtered_final[DESeq_res_df_median_info_filtered_final$Full_ID %in% rownames(DESEQ_Ctrl_res_signif), ]
head(res_Ctrl); dim(res_Ctrl)
length(unique(res_Ctrl$ensembl_gene_id))


# LdKO
DESEQ = dds[,grep('LdKO', colnames(dds))]
DESEQ = DESeqDataSet(DESEQ, design=design)
DESEQ = DESeq(DESEQ)
DESEQ_res = results(DESEQ, contrast=c('name', levels(colData(DESEQ)$name)[grep('mut', levels(colData(DESEQ)$name))], levels(colData(DESEQ)$name)[grep('wt', levels(colData(DESEQ)$name))]))
DESEQ_res

DESEQ_LdKO = DESEQ
DESEQ_LdKO_res = DESEQ_res


DESEQ_LdKO_res$Full_ID = rownames(DESEQ_LdKO_res)
DESEQ_LdKO_res_info = left_join(data.frame(DESEQ_LdKO_res), DESeq_res_df_median_info_filtered_final)
head(DESEQ_LdKO_res_info); dim(DESEQ_LdKO_res_info)

first_row = DESEQ_LdKO_res_info[1,]; first_row

if (nrow(first_row[(first_row$log2FoldChange > 0) & first_row$LdKO_Smg6mut > first_row$LdKO_Smg6wt, ]) == 0){
  if (nrow(first_row[(first_row$log2FoldChange < 0) & first_row$LdKO_Smg6mut < first_row$LdKO_Smg6wt, ]) == 0){
    print('Problem! log2FoldChange is positive and Ctrl_Smg6mut expression is lower than Ctrl_Smg6_wt expression or vice versa')
    print('Ctrl_Smg6mut expression should be bigger than Ctrl_Smg6_wt expression or vice versa')
  }
}


DESEQ_LdKO_res_signif = DESEQ_LdKO_res[which(DESEQ_LdKO_res$padj < 0.05),]
DESEQ_LdKO_res_signif; dim(DESEQ_LdKO_res_signif)

res_LdKO = DESeq_res_df_median_info_filtered_final[DESeq_res_df_median_info_filtered_final$Full_ID %in% rownames(DESEQ_LdKO_res_signif), ]
head(res_LdKO); dim(res_LdKO)
length(unique(res_LdKO$ensembl_gene_id))


# Plots
pdf(file=paste0(Results_path, 'VolcanoPlots.pdf'), width = 8, height=8)
EnhancedVolcano_plot_fun_rescale(DESEQ_Ctrl_res_info, 'DESeq: Control')
EnhancedVolcano_plot_fun_rescale(DESEQ_LdKO_res_info, 'DESeq: LdKO')
dev.off()


pdf(file=paste0(Results_path, 'MAPlots.pdf'), width = 8, height=8)
MA_plot(DESEQ_Ctrl_res, res_Ctrl, 'DESeq: Control')
MA_plot(DESEQ_LdKO_res, res_LdKO, 'DESeq: LdKO')
dev.off()
######################################################################
