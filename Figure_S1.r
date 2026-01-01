library(RColorBrewer)
library(tidyverse)
library(ggrepel)
library(ggrastr)
options(ggrastr.default.dpi=150)
library(patchwork)

library(SummarizedExperiment)                 # SummarizedExperiment
library(rtracklayer)                          # import/export BED files

library(csaw)                                 # windowing approach
library(edgeR)     

# annotation
library(org.Mm.eg.db)                         # ID conversion 
library(BSgenome.Mmusculus.UCSC.mm10)         # genome sequence
library(TxDb.Mmusculus.UCSC.mm10.knownGene)   # gene and transcript coordinates

# coverage tracks
library(Gviz) 
library(BiocParallel)



blacklist.gr <- import("~/Documents/TSC_ATAC/mm10-blacklist.v2.bed.gz",seqinfo=seqinfo(BSgenome.Mmusculus.UCSC.mm10))
# extension
readLength <- 150
blacklist.gr <- trim(resize(blacklist.gr,width=width(blacklist.gr) + 2*readLength, fix="center"))
# add mitochondrial genome to black list
chrM.gr <- GRanges("chrM", 
                IRanges(1,seqlengths(BSgenome.Mmusculus.UCSC.mm10)["chrM"]), 
                name="Mitochondrial genome")
blacklist.gr <- c(blacklist.gr, chrM.gr)


seqlevelsStyle(blacklist.gr) <- "NCBI"
seqlevelsStyle(BSgenome.Mmusculus.UCSC.mm10) <- "NCBI"
seqlevelsStyle(TxDb.Mmusculus.UCSC.mm10.knownGene) <- "NCBI"


# table with link to BAM files and sample annotation
alignments_ATAC <- read.table("~/Documents/TSC_ATAC/metadata.csv", header=TRUE, stringsAsFactors=FALSE,sep=',')
alignments_USF2 <- read.table("/Volumes/Extreme SSD/CUTRUN/USF2_metadata.csv", header=TRUE, stringsAsFactors=FALSE,sep=',')
alignments_H3K27Ac <- read.table("/Volumes/Extreme SSD/CUTRUN/H3K27Ac_metadata.csv", header=TRUE, stringsAsFactors=FALSE,sep=',')
alignments_HDAC1 <- read.table("/Volumes/Extreme SSD/CUTRUN/HDAC1_broad_metadata.csv", header=TRUE, stringsAsFactors=FALSE,sep=',')
alignments_HDAC2 <- read.table("/Volumes/Extreme SSD/CUTRUN/HDAC2_broad_metadata.csv", header=TRUE, stringsAsFactors=FALSE,sep=',')
alignments_H3K4me3 <- read.table("/Volumes/Extreme SSD/CUTRUN/HDAC2_broad_metadata.csv", header=TRUE, stringsAsFactors=FALSE,sep=',')
alignments_TFE3 <- read.table("/Volumes/Extreme SSD/CUTRUN/TFE3_metadata.csv", header=TRUE, stringsAsFactors=FALSE,sep=',')


# BAM files
bamFiles_ATAC <- setNames(alignments_ATAC$Filename, alignments_ATAC$ExternalSampleName)
bamFiles_USF2 <- setNames(alignments_USF2$Filename, alignments_USF2$ExternalSampleName)
bamFiles_H3K27Ac <- setNames(alignments_H3K27Ac$Filename, alignments_H3K27Ac$ExternalSampleName)
bamFiles_HDAC1 <- setNames(alignments_HDAC1$Filename, alignments_HDAC1$ExternalSampleName)
bamFiles_HDAC2 <- setNames(alignments_HDAC2$Filename, alignments_HDAC2$ExternalSampleName)
bamFiles_H3K4me3 <- setNames(alignments_H3K4me3$Filename, alignments_H3K4me3$ExternalSampleName)
bamFiles_TFE3 <- setNames(alignments_TFE3$Filename, alignments_TFE3$ExternalSampleName)


param <- readParam(pe="none", dedup=TRUE, discard=blacklist.gr,
                   restrict=standardChromosomes(BSgenome.Mmusculus.UCSC.mm10))


windows_ATAC <- windowCounts(bamFiles_ATAC, ext=0, width=150,
                     bin=TRUE, filter=0, param=param,
                     BPPARAM=MulticoreParam(12))
windows_USF2 <- windowCounts(bamFiles_USF2, ext=0, width=150,
                     bin=TRUE, filter=0, param=param,
                     BPPARAM=MulticoreParam(12))
windows_H3K27Ac <- windowCounts(bamFiles_H3K27Ac, ext=0, width=150,
                     bin=TRUE, filter=0, param=param,
                     BPPARAM=MulticoreParam(12))
windows_HDAC1 <- windowCounts(bamFiles_HDAC1, ext=0, width=150,
                     bin=TRUE, filter=0, param=param,
                     BPPARAM=MulticoreParam(12))
windows_HDAC2 <- windowCounts(bamFiles_HDAC2, ext=0, width=150,
                     bin=TRUE, filter=0, param=param,
                     BPPARAM=MulticoreParam(12))
windows_H3K4me3 <- windowCounts(bamFiles_H3K4me3, ext=0, width=150,
                     bin=TRUE, filter=0, param=param,
                     BPPARAM=MulticoreParam(12))
windows_TFE3 <- windowCounts(bamFiles_TFE3, ext=0, width=150,
                     bin=TRUE, filter=0, param=param,
                     BPPARAM=MulticoreParam(12))

# Add sample annotation
colData(windows_ATAC) <- cbind(colData(windows_ATAC), alignments_ATAC)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, rowRanges(windows_ATAC))
rowData(windows_ATAC) <- letterFrequency(seqs, letters=c('GC', 'N'), as.prob=TRUE) # G or C, N frequency
is.agap <- rowData(windows_ATAC)[,'N'] > 2/3 # i.e. more than 100 of 150 bases
table(is.agap)/nrow(windows_ATAC)
is.smaller <- width(rowRanges(windows_ATAC)) != 150
table(is.smaller)
windows_ATAC <- windows_ATAC[!is.agap & !is.smaller,]

colData(windows_USF2) <- cbind(colData(windows_USF2), alignments_USF2)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, rowRanges(windows_USF2))
rowData(windows_USF2) <- letterFrequency(seqs, letters=c('GC', 'N'), as.prob=TRUE) # G or C, N frequency
is.agap <- rowData(windows_USF2)[,'N'] > 2/3 # i.e. more than 100 of 150 bases
table(is.agap)/nrow(windows_USF2)
is.smaller <- width(rowRanges(windows_USF2)) != 150
table(is.smaller)
windows_USF2 <- windows_USF2[!is.agap & !is.smaller,]

colData(windows_H3K27Ac) <- cbind(colData(windows_H3K27Ac), alignments_H3K27Ac)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, rowRanges(windows_H3K27Ac))
rowData(windows_H3K27Ac) <- letterFrequency(seqs, letters=c('GC', 'N'), as.prob=TRUE) # G or C, N frequency
is.agap <- rowData(windows_H3K27Ac)[,'N'] > 2/3 # i.e. more than 100 of 150 bases
table(is.agap)/nrow(windows_H3K27Ac)
is.smaller <- width(rowRanges(windows_H3K27Ac)) != 150
table(is.smaller)
windows_H3K27Ac <- windows_H3K27Ac[!is.agap & !is.smaller,]

colData(windows_HDAC1) <- cbind(colData(windows_HDAC1), alignments_HDAC1)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, rowRanges(windows_HDAC1))
rowData(windows_HDAC1) <- letterFrequency(seqs, letters=c('GC', 'N'), as.prob=TRUE) # G or C, N frequency
is.agap <- rowData(windows_HDAC1)[,'N'] > 2/3 # i.e. more than 100 of 150 bases
table(is.agap)/nrow(windows_HDAC1)
is.smaller <- width(rowRanges(windows_HDAC1)) != 150
table(is.smaller)
windows_HDAC1 <- windows_HDAC1[!is.agap & !is.smaller,]

colData(windows_HDAC2) <- cbind(colData(windows_HDAC2), alignments_HDAC2)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, rowRanges(windows_HDAC2))
rowData(windows_HDAC2) <- letterFrequency(seqs, letters=c('GC', 'N'), as.prob=TRUE) # G or C, N frequency
is.agap <- rowData(windows_HDAC2)[,'N'] > 2/3 # i.e. more than 100 of 150 bases
table(is.agap)/nrow(windows_HDAC2)
is.smaller <- width(rowRanges(windows_HDAC2)) != 150
table(is.smaller)
windows_HDAC2 <- windows_HDAC2[!is.agap & !is.smaller,]

colData(windows_H3K4me3) <- cbind(colData(windows_H3K4me3), alignments_H3K4me3)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, rowRanges(windows_H3K4me3))
rowData(windows_H3K4me3) <- letterFrequency(seqs, letters=c('GC', 'N'), as.prob=TRUE) # G or C, N frequency
is.agap <- rowData(windows_H3K4me3)[,'N'] > 2/3 # i.e. more than 100 of 150 bases
table(is.agap)/nrow(windows_H3K4me3)
is.smaller <- width(rowRanges(windows_H3K4me3)) != 150
table(is.smaller)
windows_H3K4me3 <- windows_H3K4me3[!is.agap & !is.smaller,]

colData(windows_TFE3) <- cbind(colData(windows_TFE3), alignments_TFE3)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, rowRanges(windows_TFE3))
rowData(windows_TFE3) <- letterFrequency(seqs, letters=c('GC', 'N'), as.prob=TRUE) # G or C, N frequency
is.agap <- rowData(windows_TFE3)[,'N'] > 2/3 # i.e. more than 100 of 150 bases
table(is.agap)/nrow(windows_TFE3)
is.smaller <- width(rowRanges(windows_TFE3)) != 150
table(is.smaller)
windows_TFE3 <- windows_TFE3[!is.agap & !is.smaller,]


#TSS annotation
tss <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
tss <- resize(tss, width=1, fix="start")
overlapTSS <- rowRanges(windows_ATAC) %over% tss
rowData(windows_ATAC)$overlapTSS <- overlapTSS
overlapTSS <- rowRanges(windows_USF2) %over% tss
rowData(windows_USF2)$overlapTSS <- overlapTSS
overlapTSS <- rowRanges(windows_H3K27Ac) %over% tss
rowData(windows_H3K27Ac)$overlapTSS <- overlapTSS
overlapTSS <- rowRanges(windows_HDAC1) %over% tss
rowData(windows_HDAC1)$overlapTSS <- overlapTSS
overlapTSS <- rowRanges(windows_HDAC2) %over% tss
rowData(windows_HDAC2)$overlapTSS <- overlapTSS
overlapTSS <- rowRanges(windows_H3K4me3) %over% tss
rowData(windows_H3K4me3)$overlapTSS <- overlapTSS
overlapTSS <- rowRanges(windows_TFE3) %over% tss
rowData(windows_TFE3)$overlapTSS <- overlapTSS

combi <- list(c('WT','TSC'), c('WT','TFE3'), c('WT','dKO'),c('TFE3','dKO'))


plotMA <- function(df, s1, s2) {
  ggplot(df,
         mapping = aes(x = .data[[s1]] + .data[[s2]],
                       y = .data[[s1]] - .data[[s2]])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, overlapTSS), col='red', pch='.') ) +
    rasterize( geom_smooth(data = subset(df, overlapTSS), col='plum', se = FALSE) ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", s1, s2)) +
    theme_bw(base_size = 15)
}

# GC Plot
plotGC <- function(df, s1, s2) {
  ggplot(df,
         mapping = aes(x = GC,
                       y = .data[[s1]] - .data[[s2]])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, overlapTSS), col='red', pch='.') ) +
    rasterize( geom_smooth(data = subset(df, overlapTSS), col='plum', se = FALSE) ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "GC", y="M", title=sprintf("%s vs. %s", s1, s2)) +
    theme_bw(base_size = 15)
}


saveRDS(windows_ATAC,file='~/Documents/TSC_ATAC/windows_ATAC.rds')
saveRDS(windows_USF2,file='~/Documents/TSC_ATAC/windows_USF2.rds')
saveRDS(windows_H3K27Ac,file='~/Documents/TSC_ATAC/windows_H3K27Ac.rds')
saveRDS(windows_HDAC1,file='~/Documents/TSC_ATAC/windows_HDAC1.rds')
saveRDS(windows_HDAC2,file='~/Documents/TSC_ATAC/windows_HDAC2.rds')
saveRDS(windows_H3K4me3,file='~/Documents/TSC_ATAC/windows_H3K4me3.rds')
saveRDS(windows_TFE3,file='~/Documents/TSC_ATAC/windows_TFE3.rds')


#windows_ATAC <-readRDS(file='~/Documents/TSC_ATAC/windows_ATAC.rds')


# Filter windows from MACS2

peaks_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_summits.bed",format="bed")
peaks_USF2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/USF2_summits.bed",format="bed")
peaks_H3K27Ac=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/H3K27Ac_summits.bed",format="bed")
peaks_HDAC1=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC1_summits.bed",format="bed")
peaks_HDAC2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC2_summits.bed",format="bed")
peaks_H3K4me3=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K4me3_summits.bed",format="bed")
peaks_TFE3=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/TFE3_summits.bed",format="bed")



# count reads in windows specified by MACS2                                      
peak.counts <- regionCounts(bamFiles_ATAC, peaks_ATAC, param=param)
peak.abundances <- aveLogCPM(asDGEList(peak.counts))
peak.counts.filt <- peak.counts[peak.abundances > -3, ]
suppressWarnings(keep <- overlapsAny(rowRanges(windows_ATAC), peak.counts.filt))
filtered.windows_ATAC <- windows_ATAC[keep,]

peak.counts <- regionCounts(bamFiles_USF2, peaks_USF2, param=param)
peak.abundances <- aveLogCPM(asDGEList(peak.counts))
peak.counts.filt <- peak.counts[peak.abundances > -3, ]
suppressWarnings(keep <- overlapsAny(rowRanges(windows_USF2), peak.counts.filt))
filtered.windows_USF2 <- windows_USF2[keep,]

peak.counts <- regionCounts(bamFiles_H3K27Ac, peaks_H3K27Ac, param=param)
peak.abundances <- aveLogCPM(asDGEList(peak.counts))
peak.counts.filt <- peak.counts[peak.abundances > -3, ]
suppressWarnings(keep <- overlapsAny(rowRanges(windows_H3K27Ac), peak.counts.filt))
filtered.windows_H3K27Ac <- windows_H3K27Ac[keep,]

peak.counts <- regionCounts(bamFiles_HDAC1, peaks_HDAC1, param=param)
peak.abundances <- aveLogCPM(asDGEList(peak.counts))
peak.counts.filt <- peak.counts[peak.abundances > -3, ]
suppressWarnings(keep <- overlapsAny(rowRanges(windows_HDAC1), peak.counts.filt))
filtered.windows_HDAC1 <- windows_HDAC1[keep,]

peak.counts <- regionCounts(bamFiles_HDAC2, peaks_HDAC2, param=param)
peak.abundances <- aveLogCPM(asDGEList(peak.counts))
peak.counts.filt <- peak.counts[peak.abundances > -3, ]
suppressWarnings(keep <- overlapsAny(rowRanges(windows_HDAC2), peak.counts.filt))
filtered.windows_HDAC2 <- windows_HDAC2[keep,]

peak.counts <- regionCounts(bamFiles_H3K4me3, peaks_H3K4me3, param=param)
peak.abundances <- aveLogCPM(asDGEList(peak.counts))
peak.counts.filt <- peak.counts[peak.abundances > -3, ]
suppressWarnings(keep <- overlapsAny(rowRanges(windows_H3K4me3), peak.counts.filt))
filtered.windows_H3K4me3 <- windows_H3K4me3[keep,]

peak.counts <- regionCounts(bamFiles_TFE3, peaks_TFE3, param=param)
peak.abundances <- aveLogCPM(asDGEList(peak.counts))
peak.counts.filt <- peak.counts[peak.abundances > -3, ]
suppressWarnings(keep <- overlapsAny(rowRanges(windows_TFE3), peak.counts.filt))
filtered.windows_TFE3 <- windows_TFE3[keep,]

filtered.windows_ATAC <- normOffsets(filtered.windows_ATAC,span=.3,se.out=TRUE)
filtered.windows_USF2 <- normOffsets(filtered.windows_USF2,span=.3,se.out=TRUE)
filtered.windows_H3K27Ac <- normOffsets(filtered.windows_H3K27Ac,span=.3,se.out=TRUE)
filtered.windows_HDAC1 <- normOffsets(filtered.windows_HDAC1,span=.3,se.out=TRUE)
filtered.windows_HDAC2 <- normOffsets(filtered.windows_HDAC2,span=.3,se.out=TRUE)
filtered.windows_H3K4me3 <- normOffsets(filtered.windows_H3K4me3,span=.3,se.out=TRUE)
filtered.windows_TFE3 <- normOffsets(filtered.windows_TFE3,span=.3,se.out=TRUE)


saveRDS(filtered.windows_ATAC,file='~/Documents/TSC_ATAC/filtered.windows_ATAC.rds')
saveRDS(filtered.windows_USF2,file='~/Documents/TSC_ATAC/filtered.windows_USF2.rds')
saveRDS(filtered.windows_H3K27Ac,file='~/Documents/TSC_ATAC/filtered.windows_H3K27Ac.rds')
saveRDS(filtered.windows_HDAC1,file='~/Documents/TSC_ATAC/filtered.windows_HDAC1.rds')
saveRDS(filtered.windows_HDAC2,file='~/Documents/TSC_ATAC/filtered.windows_HDAC2.rds')
saveRDS(filtered.windows_H3K4me3,file='~/Documents/TSC_ATAC/filtered.windows_H3K4me3.rds')
saveRDS(filtered.windows_TFE3,file='~/Documents/TSC_ATAC/filtered.windows_TFE3.rds')



#filtered.windows_H3K27Ac <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_H3K27Ac.rds')
#filtered.windows_H3K4me3 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_H3K4me3.rds')

#############################################
######### Baseline Characterization #########
#############################################

### Behavior of chromatin generally by condition

peaks_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_peaks.narrowPeak")
summits_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_summits.bed")




peakAnno.edb <- annotatePeak(peaks_ATAC, tssRegion=c(-3000, 3000),
                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

anno <- as.data.frame(peakAnno.edb)$annotation
prom <- anno == 'Promoter (<=1kb)'


peaks_ATAC.filter <- peaks_ATAC[prom]

cutoff <- sort(peaks_ATAC.filter$score,decreasing=T)[5000]
idx <- peaks_ATAC.filter[peaks_ATAC.filter$score >= cutoff]$name


anchors <- summits_ATAC[summits_ATAC$name %in% idx]

ATAC_bw_WT1 <- import("~/Documents/TSC_ATAC/1706_WT.bigWig")
ATAC_bw_WT2 <- import("~/Documents/TSC_ATAC/1749_WT.bigWig")
ATAC_bw_WT3 <- import("~/Documents/TSC_ATAC/1750_WT.bigWig")
ATAC_bw_TSC1 <- import("~/Documents/TSC_ATAC/1747_TSC.bigWig")
ATAC_bw_TSC2 <- import("~/Documents/TSC_ATAC/1748_TSC.bigWig")
ATAC_bw_TSC3 <- import("~/Documents/TSC_ATAC/1755_TSC.bigWig")
ATAC_bw_TFE31 <- import("~/Documents/TSC_ATAC/1741_TFE3.bigWig")
ATAC_bw_TFE32 <- import("~/Documents/TSC_ATAC/1746_TFE3.bigWig")
ATAC_bw_TFE33 <- import("~/Documents/TSC_ATAC/1751_TFE3.bigWig")
ATAC_bw_dKO1 <- import("~/Documents/TSC_ATAC/1708_dKO.bigWig")
ATAC_bw_dKO2 <- import("~/Documents/TSC_ATAC/1744_dKO.bigWig")
ATAC_bw_dKO3 <- import("~/Documents/TSC_ATAC/1754_dKO.bigWig")

library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(ATAC_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(ATAC_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(ATAC_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(ATAC_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(ATAC_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(ATAC_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(ATAC_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(ATAC_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(ATAC_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(ATAC_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(ATAC_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(ATAC_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC2 + mat_TSC3) / 2
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3


col_fun<- circlize::colorRamp2(c(0, .5), c("white", "red"))

ht_list <- EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/ATAC_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/ATAC_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/ATAC_peaks_ALL_mean.pdf',height=2, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

maxl <- c(max(WT1_mean),max(WT2_mean),max(WT3_mean),max(TSC2_mean),max(TSC3_mean),
max(TFE31_mean),max(TFE32_mean),max(TFE33_mean),max(dKO1_mean),max(dKO2_mean),max(dKO3_mean))
maxl <- data.frame(maxl)
maxl$group <- c(rep('WT',3),rep('TSC',2),rep('TFE3',3),rep('dKO',3))

pdf('~/Documents/TSC_Paper/ATAC_peaks_box.pdf',height=3, width=4)
ggplot(maxl,aes(x=group,y=maxl)) + geom_boxplot() + theme_classic() + geom_point()
dev.off()


df <- data.frame(WT = rowMeans(mat_WT[partition == 'promoters',c(20:21)]),
  TSC1 = rowMeans(mat_TSC[partition == 'promoters',c(20:21)]),
  TFE3 = rowMeans(mat_TFE3[partition == 'promoters',c(20:21)]),
  dKO = rowMeans(mat_dKO[partition == 'promoters',c(20:21)]))
idx <- dim(df)[1]

df <- data.frame(unlist(df))
df$group <- c(rep('WT',idx),rep('TSC1',idx),rep('TFE3',idx),rep('dKO',idx))
colnames(df) <- c('val','group')

pdf('~/Documents/TSC_Paper/ATAC_promoter_peaks_dist.pdf',height=3, width=4)
ggplot(df, aes(x = val,color=group, fill = group)) +
  geom_density(alpha=0.1,bound=c(0,5),bw=.15) + xlim(0,2.5) + 
  scale_fill_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic()
dev.off()



### ATAC at TFE3 peaks


summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")

peaks_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_peaks.narrowPeak")
summits_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_summits.bed")

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"


peakAnno.edb <- annotatePeak(peaks_TFE3_CHIP_WT, tssRegion=c(-3000, 3000),
                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

anno <- as.data.frame(peakAnno.edb)$annotation
prom <- anno == 'Promoter (<=1kb)'

peaks_TFE3_CHIP_WT.filt <- peaks_TFE3_CHIP_WT[prom]

peakAnno.edb <- annotatePeak(peaks_ATAC, tssRegion=c(-3000, 3000),
                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

anno <- as.data.frame(peakAnno.edb)$annotation
prom <- anno == 'Promoter (<=1kb)'


peaks_ATAC.filter <- peaks_ATAC[prom]

cutoff <- sort(peaks_ATAC.filter$score,decreasing=T)[5000]
idx <- peaks_ATAC.filter[peaks_ATAC.filter$score >= cutoff]$name

#ATAC_overlap_TFE3_peaks<- subsetByOverlaps(peaks_TFE3_CHIP_FLCNKO, peaks_ATAC)
ATAC_overlap_TFE3_peaks<- subsetByOverlaps(peaks_ATAC.filter,peaks_TFE3_CHIP_WT.filt)


#anchors <- summits_TFE3_CHIP_FLCNKO[summits_TFE3_CHIP_FLCNKO$name %in% ATAC_overlap_TFE3_peaks$name]
#anchors <- anchors[anchors$name %in% idx]

anchors <- summits_ATAC[summits_ATAC$name %in% ATAC_overlap_TFE3_peaks$name]
anchors <- anchors[anchors$name %in% idx]

ATAC_bw_WT1 <- import("~/Documents/TSC_ATAC/1706_WT.bigWig")
ATAC_bw_WT2 <- import("~/Documents/TSC_ATAC/1749_WT.bigWig")
ATAC_bw_WT3 <- import("~/Documents/TSC_ATAC/1750_WT.bigWig")
ATAC_bw_TSC1 <- import("~/Documents/TSC_ATAC/1747_TSC.bigWig")
ATAC_bw_TSC2 <- import("~/Documents/TSC_ATAC/1748_TSC.bigWig")
ATAC_bw_TSC3 <- import("~/Documents/TSC_ATAC/1755_TSC.bigWig")
ATAC_bw_TFE31 <- import("~/Documents/TSC_ATAC/1741_TFE3.bigWig")
ATAC_bw_TFE32 <- import("~/Documents/TSC_ATAC/1746_TFE3.bigWig")
ATAC_bw_TFE33 <- import("~/Documents/TSC_ATAC/1751_TFE3.bigWig")
ATAC_bw_dKO1 <- import("~/Documents/TSC_ATAC/1708_dKO.bigWig")
ATAC_bw_dKO2 <- import("~/Documents/TSC_ATAC/1744_dKO.bigWig")
ATAC_bw_dKO3 <- import("~/Documents/TSC_ATAC/1754_dKO.bigWig")

library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(ATAC_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(ATAC_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(ATAC_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(ATAC_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(ATAC_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(ATAC_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(ATAC_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(ATAC_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(ATAC_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(ATAC_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(ATAC_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(ATAC_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC2 + mat_TSC3) / 2
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3


col_fun<- circlize::colorRamp2(c(0, 1.25), c("white", "red"))

ht_list <- EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/ATAC_TFE3_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/ATAC_TFE3_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/ATAC_TFE3_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

##############################
######### Figure SXA #########
##############################

# Show that chromatin at TFE3 targets is more accessible in TSC KO

keep <- overlapsAny(rowRanges(normed.windows_ATAC), summits_TFE3_CHIP_FLCNKO)
TFE3.filtered.windows.ATAC <- normed.windows_ATAC[keep,]
colData(TFE3.filtered.windows.ATAC)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.ATAC)$ExternalSampleName)
TFE3.filtered.windows.ATAC <- normOffsets(TFE3.filtered.windows.ATAC, se.out=TRUE)



#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.ATAC, log=TRUE, use.offsets= TRUE)

ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.ATAC))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.ATAC)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.ATAC)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.ATAC,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.ATAC)),
                        group = colData(TFE3.filtered.windows.ATAC)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.ATAC)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.ATAC))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)
#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.ATAC) <- cbind(rowData(TFE3.filtered.windows.ATAC), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.ATAC), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/peaks_WT_vs_TSC.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/peaks_WT_vs_TSC.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.05 & direction == 'down' & overlap != "" & name != 'Common')$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))


library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)

dbs <- c("ChEA_2022","Reactome_2024","GO_Biological_Process_2023")


enriched <- enrichr(genes_all, dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.05,]
p1<- ggplot(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),], (aes(x=-log10(Adjusted.P.value), y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + geom_point() + xlab('Adjusted P Value') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + ggtitle('ChEA') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"\\s+"), `[`, 1))) + theme(axis.text=element_text(colour="black"))


enriched <- enrichr(genes_proms, dbs)
enriched[[3]] <- enriched[[3]][enriched[[3]]$Adjusted.P.value < 0.05,]
p2<- ggplot(enriched[[3]][order(enriched[[3]]$Adjusted.P.value,decreasing=F),][rev(1:3),], (aes(x=-log10(Adjusted.P.value), y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + geom_point() + xlab('Adjusted P Value') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + ggtitle('GO_Biological_Process_2023') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[3]][order(enriched[[3]]$Adjusted.P.value,decreasing=F),][rev(1:3),]$Term,"\\(GO+"), `[`, 1))) + theme(axis.text=element_text(colour="black"))

pdf('~/Documents/TSC_Paper/ATAC_TFE3_Targets_WT_vs_TSC_UP.pdf',width=4.75,height=4)
p1
dev.off()

pdf('~/Documents/TSC_Paper/ATAC_TFE3_Targets_WT_vs_TSC_GO.pdf',width=4.75,height=4)
p2
dev.off()

######## Extract ATAC peaks which are DOPs

DOWs <- TFE3.filtered.windows.ATAC[rowData(TFE3.filtered.windows.ATAC)$FDR < 0.05 & rowData(TFE3.filtered.windows.ATAC)$logFC < 0]
#DOPs <- subsetByOverlaps(peaks_TFE3_CHIP_FLCNKO, DOWs)
#DOPs_WT <- subsetByOverlaps(peaks_TFE3_CHIP_WT, DOWs)
DOPs <- subsetByOverlaps(peaks_ATAC, DOWs)

#anchors <- summits_TFE3_CHIP_FLCNKO[summits_TFE3_CHIP_FLCNKO$name %in% DOPs$name]
anchors <- summits_ATAC[summits_ATAC$name %in% DOPs$name]


library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(ATAC_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT1
mat_WT2<- normalizeToMatrix(ATAC_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT2
mat_WT3<- normalizeToMatrix(ATAC_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT3

mat_TSC1<- normalizeToMatrix(ATAC_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC1
mat_TSC2<- normalizeToMatrix(ATAC_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC2
mat_TSC3<- normalizeToMatrix(ATAC_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC3

mat_TFE31<- normalizeToMatrix(ATAC_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE31
mat_TFE32<- normalizeToMatrix(ATAC_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE32
mat_TFE33<- normalizeToMatrix(ATAC_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE33

mat_dKO1<- normalizeToMatrix(ATAC_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO1
mat_dKO2<- normalizeToMatrix(ATAC_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO2
mat_dKO3<- normalizeToMatrix(ATAC_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO3


quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3


col_fun<- circlize::colorRamp2(c(0, .5), c("white", "red"))

ht_list <- EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/ATAC_TFE3_heatmap_DOPs_ALL.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/ATAC_TFE3_peaks_DOPs_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/ATAC_TFE3_peaks_DOPs_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()


### Show that specific chromatin closes in TFE3 KO

keep <- overlapsAny(rowRanges(normed.windows_ATAC), summits_TFE3_CHIP_FLCNKO)
TFE3.filtered.windows.ATAC <- normed.windows_ATAC[keep,]
colData(TFE3.filtered.windows.ATAC)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.ATAC)$ExternalSampleName)
TFE3.filtered.windows.ATAC <- normOffsets(TFE3.filtered.windows.ATAC, se.out=TRUE)



#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.ATAC, log=TRUE, use.offsets= TRUE)

ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.ATAC))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.ATAC)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.ATAC)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.ATAC,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.ATAC)),
                        group = colData(TFE3.filtered.windows.ATAC)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.ATAC)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.ATAC))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TFE3, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TFE3.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TFE3']]  + .data[['WT']],y = .data[['WT']] - .data[['TFE3']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TFE3.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TFE3')) +
    theme_bw(base_size = 15)
#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.ATAC) <- cbind(rowData(TFE3.filtered.windows.ATAC), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.ATAC), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TFE3', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/peaks_WT_vs_TFE3.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/peaks_WT_vs_TFE3.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "" & name != 'Common')$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

# Find no evidence of TFE3 chromatin closing



######## Extract ATAC peaks which are DOPs

DOWs <- TFE3.filtered.windows.ATAC[rowData(TFE3.filtered.windows.ATAC)$FDR < 0.1 & rowData(TFE3.filtered.windows.ATAC)$logFC > 0]
#DOPs <- subsetByOverlaps(peaks_TFE3_CHIP_FLCNKO, DOWs)
#DOPs_WT <- subsetByOverlaps(peaks_TFE3_CHIP_WT, DOWs)
DOPs <- subsetByOverlaps(peaks_ATAC, DOWs)

#anchors <- summits_TFE3_CHIP_FLCNKO[summits_TFE3_CHIP_FLCNKO$name %in% DOPs$name]
anchors <- summits_ATAC[summits_ATAC$name %in% DOPs$name]


library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(ATAC_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT1
mat_WT2<- normalizeToMatrix(ATAC_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT2
mat_WT3<- normalizeToMatrix(ATAC_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT3

mat_TSC1<- normalizeToMatrix(ATAC_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC1
mat_TSC2<- normalizeToMatrix(ATAC_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC2
mat_TSC3<- normalizeToMatrix(ATAC_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC3

mat_TFE31<- normalizeToMatrix(ATAC_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE31
mat_TFE32<- normalizeToMatrix(ATAC_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE32
mat_TFE33<- normalizeToMatrix(ATAC_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE33

mat_dKO1<- normalizeToMatrix(ATAC_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO1
mat_dKO2<- normalizeToMatrix(ATAC_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO2
mat_dKO3<- normalizeToMatrix(ATAC_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO3


quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3


col_fun<- circlize::colorRamp2(c(0, 1), c("white", "red"))

ht_list <- EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/ATAC_TFE3_heatmap_DOPs_TFE3KO_ALL.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/ATAC_TFE3_peaks_DOPs_TFE3KO_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/ATAC_TFE3_peaks_DOPs_TFE3KO_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

### Show WT vs dKO

keep <- overlapsAny(rowRanges(normed.windows_ATAC), summits_TFE3_CHIP_FLCNKO)
TFE3.filtered.windows.ATAC <- normed.windows_ATAC[keep,]
colData(TFE3.filtered.windows.ATAC)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.ATAC)$ExternalSampleName)
TFE3.filtered.windows.ATAC <- normOffsets(TFE3.filtered.windows.ATAC, se.out=TRUE)



#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.ATAC, log=TRUE, use.offsets= TRUE)

ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.ATAC))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.ATAC)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.ATAC)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.ATAC,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.ATAC)),
                        group = colData(TFE3.filtered.windows.ATAC)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.ATAC)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.ATAC))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - dKO, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.ATAC) <- cbind(rowData(TFE3.filtered.windows.ATAC), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.ATAC), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'dKO', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/peaks_WT_vs_dKO.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/peaks_WT_vs_dKO.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.5 & direction == 'up' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))






### Behavior of chromatin in relation to H3K27Ac

summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")
peaks_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_peaks.narrowPeak")
summits_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_summits.bed")
peaks_H3K27Ac=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K27Ac_broad_peaks.broadPeak")
summits_H3K27Ac=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K27Ac_broad_summits.bed")

summits_H3K27Ac$name = peaks_H3K27Ac$name


cutoff <- sort(peaks_ATAC$score,decreasing=T)[5000]
idx <- peaks_ATAC[peaks_ATAC$score >= cutoff]$name

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"


ATAC_overlap_TFE3_peaks <- subsetByOverlaps(peaks_ATAC,peaks_TFE3_CHIP_FLCNKO)
ATAC_overlap_TFE3_H3K27Ac_peaks <- subsetByOverlaps(ATAC_overlap_TFE3_peaks,peaks_H3K27Ac)

anchors <- summits_ATAC[summits_ATAC$name %in% ATAC_overlap_TFE3_H3K27Ac_peaks$name]
anchors <- anchors[anchors$name %in% idx]

ATAC_bw_WT1 <- import("~/Documents/TSC_ATAC/1706_WT.bigWig")
ATAC_bw_WT2 <- import("~/Documents/TSC_ATAC/1749_WT.bigWig")
ATAC_bw_WT3 <- import("~/Documents/TSC_ATAC/1750_WT.bigWig")
ATAC_bw_TSC1 <- import("~/Documents/TSC_ATAC/1747_TSC.bigWig")
ATAC_bw_TSC2 <- import("~/Documents/TSC_ATAC/1748_TSC.bigWig")
ATAC_bw_TSC3 <- import("~/Documents/TSC_ATAC/1755_TSC.bigWig")
ATAC_bw_TFE31 <- import("~/Documents/TSC_ATAC/1741_TFE3.bigWig")
ATAC_bw_TFE32 <- import("~/Documents/TSC_ATAC/1746_TFE3.bigWig")
ATAC_bw_TFE33 <- import("~/Documents/TSC_ATAC/1751_TFE3.bigWig")
ATAC_bw_dKO1 <- import("~/Documents/TSC_ATAC/1708_dKO.bigWig")
ATAC_bw_dKO2 <- import("~/Documents/TSC_ATAC/1744_dKO.bigWig")
ATAC_bw_dKO3 <- import("~/Documents/TSC_ATAC/1754_dKO.bigWig")

library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(ATAC_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT1
mat_WT2<- normalizeToMatrix(ATAC_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT2
mat_WT3<- normalizeToMatrix(ATAC_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT3

mat_TSC1<- normalizeToMatrix(ATAC_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC1
mat_TSC2<- normalizeToMatrix(ATAC_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC2
mat_TSC3<- normalizeToMatrix(ATAC_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC3

mat_TFE31<- normalizeToMatrix(ATAC_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE31
mat_TFE32<- normalizeToMatrix(ATAC_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE32
mat_TFE33<- normalizeToMatrix(ATAC_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE33

mat_dKO1<- normalizeToMatrix(ATAC_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO1
mat_dKO2<- normalizeToMatrix(ATAC_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO2
mat_dKO3<- normalizeToMatrix(ATAC_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO3

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3


col_fun<- circlize::colorRamp2(c(0, 1.25), c("white", "red"))

ht_list <- EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/ATAC_TFE3_H3K27Ac_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/ATAC_TFE3_H3K27Ac_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/ATAC_TFE3_H3K27Ac_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

### Behavior of chromatin in relation to USF2

summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")
peaks_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_peaks.narrowPeak")
summits_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_summits.bed")
peaks_USF2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/USF2_peaks.narrowPeak")
summits_USF2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/USF2_summits.bed")



cutoff <- sort(peaks_ATAC$score,decreasing=T)[5000]
idx <- peaks_ATAC[peaks_ATAC$score >= cutoff]$name

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"


ATAC_overlap_TFE3_peaks <- subsetByOverlaps(peaks_ATAC,peaks_TFE3_CHIP_FLCNKO)
ATAC_overlap_TFE3_USF2_peaks <- subsetByOverlaps(ATAC_overlap_TFE3_peaks,peaks_USF2)

anchors <- summits_ATAC[summits_ATAC$name %in% ATAC_overlap_TFE3_USF2_peaks$name]
anchors <- anchors[anchors$name %in% idx]


ATAC_bw_WT1 <- import("~/Documents/TSC_ATAC/1706_WT.bigWig")
ATAC_bw_WT2 <- import("~/Documents/TSC_ATAC/1749_WT.bigWig")
ATAC_bw_WT3 <- import("~/Documents/TSC_ATAC/1750_WT.bigWig")
ATAC_bw_TSC1 <- import("~/Documents/TSC_ATAC/1747_TSC.bigWig")
ATAC_bw_TSC2 <- import("~/Documents/TSC_ATAC/1748_TSC.bigWig")
ATAC_bw_TSC3 <- import("~/Documents/TSC_ATAC/1755_TSC.bigWig")
ATAC_bw_TFE31 <- import("~/Documents/TSC_ATAC/1741_TFE3.bigWig")
ATAC_bw_TFE32 <- import("~/Documents/TSC_ATAC/1746_TFE3.bigWig")
ATAC_bw_TFE33 <- import("~/Documents/TSC_ATAC/1751_TFE3.bigWig")
ATAC_bw_dKO1 <- import("~/Documents/TSC_ATAC/1708_dKO.bigWig")
ATAC_bw_dKO2 <- import("~/Documents/TSC_ATAC/1744_dKO.bigWig")
ATAC_bw_dKO3 <- import("~/Documents/TSC_ATAC/1754_dKO.bigWig")

library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(ATAC_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT1
mat_WT2<- normalizeToMatrix(ATAC_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT2
mat_WT3<- normalizeToMatrix(ATAC_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT3

mat_TSC1<- normalizeToMatrix(ATAC_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC1
mat_TSC2<- normalizeToMatrix(ATAC_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC2
mat_TSC3<- normalizeToMatrix(ATAC_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC3

mat_TFE31<- normalizeToMatrix(ATAC_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE31
mat_TFE32<- normalizeToMatrix(ATAC_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE32
mat_TFE33<- normalizeToMatrix(ATAC_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE33

mat_dKO1<- normalizeToMatrix(ATAC_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO1
mat_dKO2<- normalizeToMatrix(ATAC_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO2
mat_dKO3<- normalizeToMatrix(ATAC_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO3

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3


col_fun<- circlize::colorRamp2(c(0, 1.5), c("white", "red"))

ht_list <- EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/ATAC_TFE3_USF2_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/ATAC_TFE3_USF2_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/ATAC_TFE3_USF2_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

### IN TFE3 KO, at USF2 peaks, the chromatin is slightly more closed

### Behavior of chromatin in relation to HDAC2

summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")
peaks_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_peaks.narrowPeak")
summits_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_summits.bed")
peaks_HDAC2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC2_peaks.narrowPeak")
summits_HDAC2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC2_summits.bed")



cutoff <- sort(peaks_ATAC$score,decreasing=T)[5000]
idx <- peaks_ATAC[peaks_ATAC$score >= cutoff]$name

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"


ATAC_overlap_TFE3_peaks <- subsetByOverlaps(peaks_ATAC,peaks_TFE3_CHIP_FLCNKO)
ATAC_overlap_TFE3_HDAC2_peaks <- subsetByOverlaps(ATAC_overlap_TFE3_peaks,peaks_HDAC2)

anchors <- summits_ATAC[summits_ATAC$name %in% ATAC_overlap_TFE3_HDAC2_peaks$name]
anchors <- anchors[anchors$name %in% idx]


ATAC_bw_WT1 <- import("~/Documents/TSC_ATAC/1706_WT.bigWig")
ATAC_bw_WT2 <- import("~/Documents/TSC_ATAC/1749_WT.bigWig")
ATAC_bw_WT3 <- import("~/Documents/TSC_ATAC/1750_WT.bigWig")
ATAC_bw_TSC1 <- import("~/Documents/TSC_ATAC/1747_TSC.bigWig")
ATAC_bw_TSC2 <- import("~/Documents/TSC_ATAC/1748_TSC.bigWig")
ATAC_bw_TSC3 <- import("~/Documents/TSC_ATAC/1755_TSC.bigWig")
ATAC_bw_TFE31 <- import("~/Documents/TSC_ATAC/1741_TFE3.bigWig")
ATAC_bw_TFE32 <- import("~/Documents/TSC_ATAC/1746_TFE3.bigWig")
ATAC_bw_TFE33 <- import("~/Documents/TSC_ATAC/1751_TFE3.bigWig")
ATAC_bw_dKO1 <- import("~/Documents/TSC_ATAC/1708_dKO.bigWig")
ATAC_bw_dKO2 <- import("~/Documents/TSC_ATAC/1744_dKO.bigWig")
ATAC_bw_dKO3 <- import("~/Documents/TSC_ATAC/1754_dKO.bigWig")

library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(ATAC_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT1
mat_WT2<- normalizeToMatrix(ATAC_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT2
mat_WT3<- normalizeToMatrix(ATAC_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT3

mat_TSC1<- normalizeToMatrix(ATAC_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC1
mat_TSC2<- normalizeToMatrix(ATAC_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC2
mat_TSC3<- normalizeToMatrix(ATAC_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC3

mat_TFE31<- normalizeToMatrix(ATAC_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE31
mat_TFE32<- normalizeToMatrix(ATAC_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE32
mat_TFE33<- normalizeToMatrix(ATAC_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE33

mat_dKO1<- normalizeToMatrix(ATAC_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO1
mat_dKO2<- normalizeToMatrix(ATAC_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO2
mat_dKO3<- normalizeToMatrix(ATAC_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO3

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3


col_fun<- circlize::colorRamp2(c(0, 1.5), c("white", "red"))

ht_list <- EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/ATAC_TFE3_HDAC2_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/ATAC_TFE3_HDAC2_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/ATAC_TFE3_HDAC2_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

### IN TFE3 KO, at HDAC2 peaks, the chromatin is slightly more closed


### Behavior of chromatin in relation to HDAC1

summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")
peaks_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_peaks.narrowPeak")
summits_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_summits.bed")
peaks_HDAC1=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC1_peaks.narrowPeak")
summits_HDAC1=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC1_summits.bed")



cutoff <- sort(peaks_ATAC$score,decreasing=T)[5000]
idx <- peaks_ATAC[peaks_ATAC$score >= cutoff]$name

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"


ATAC_overlap_TFE3_peaks <- subsetByOverlaps(peaks_ATAC,peaks_TFE3_CHIP_FLCNKO)
ATAC_overlap_TFE3_HDAC1_peaks <- subsetByOverlaps(ATAC_overlap_TFE3_peaks,peaks_HDAC1)

anchors <- summits_ATAC[summits_ATAC$name %in% ATAC_overlap_TFE3_HDAC1_peaks$name]
anchors <- anchors[anchors$name %in% idx]


ATAC_bw_WT1 <- import("~/Documents/TSC_ATAC/1706_WT.bigWig")
ATAC_bw_WT2 <- import("~/Documents/TSC_ATAC/1749_WT.bigWig")
ATAC_bw_WT3 <- import("~/Documents/TSC_ATAC/1750_WT.bigWig")
ATAC_bw_TSC1 <- import("~/Documents/TSC_ATAC/1747_TSC.bigWig")
ATAC_bw_TSC2 <- import("~/Documents/TSC_ATAC/1748_TSC.bigWig")
ATAC_bw_TSC3 <- import("~/Documents/TSC_ATAC/1755_TSC.bigWig")
ATAC_bw_TFE31 <- import("~/Documents/TSC_ATAC/1741_TFE3.bigWig")
ATAC_bw_TFE32 <- import("~/Documents/TSC_ATAC/1746_TFE3.bigWig")
ATAC_bw_TFE33 <- import("~/Documents/TSC_ATAC/1751_TFE3.bigWig")
ATAC_bw_dKO1 <- import("~/Documents/TSC_ATAC/1708_dKO.bigWig")
ATAC_bw_dKO2 <- import("~/Documents/TSC_ATAC/1744_dKO.bigWig")
ATAC_bw_dKO3 <- import("~/Documents/TSC_ATAC/1754_dKO.bigWig")

library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(ATAC_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT1
mat_WT2<- normalizeToMatrix(ATAC_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT2
mat_WT3<- normalizeToMatrix(ATAC_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT3

mat_TSC1<- normalizeToMatrix(ATAC_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC1
mat_TSC2<- normalizeToMatrix(ATAC_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC2
mat_TSC3<- normalizeToMatrix(ATAC_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC3

mat_TFE31<- normalizeToMatrix(ATAC_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE31
mat_TFE32<- normalizeToMatrix(ATAC_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE32
mat_TFE33<- normalizeToMatrix(ATAC_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE33

mat_dKO1<- normalizeToMatrix(ATAC_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO1
mat_dKO2<- normalizeToMatrix(ATAC_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO2
mat_dKO3<- normalizeToMatrix(ATAC_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO3

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3


col_fun<- circlize::colorRamp2(c(0, 1.5), c("white", "red"))

ht_list <- EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/ATAC_TFE3_HDAC1_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/ATAC_TFE3_HDAC1_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/ATAC_TFE3_HDAC1_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

### IN TFE3 KO, at HDAC1 peaks, the chromatin is slightly more closed


### Behavior of chromatin in relation to HDAC1

summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")
peaks_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_peaks.narrowPeak")
summits_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_summits.bed")

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"


ATAC_overlap_TFE3_peaks <- subsetByOverlaps(peaks_ATAC,peaks_TFE3_CHIP_FLCNKO)
ATAC_overlap_TFE3_HDAC1_peaks <- subsetByOverlaps(ATAC_overlap_TFE3_peaks,peaks_HDAC1)
ATAC_overlap_TFE3_HDAC1_2_peaks <- subsetByOverlaps(ATAC_overlap_TFE3_HDAC1_peaks,peaks_HDAC2)
ATAC_overlap_TFE3_HDAC1_2_USF2_peaks <- subsetByOverlaps(ATAC_overlap_TFE3_HDAC1_2_peaks,peaks_USF2)


anchors <- summits_ATAC[summits_ATAC$name %in% ATAC_overlap_TFE3_HDAC1_2_USF2_peaks$name]
anchors <- anchors[anchors$name %in% idx]


ATAC_bw_WT1 <- import("~/Documents/TSC_ATAC/1706_WT.bigWig")
ATAC_bw_WT2 <- import("~/Documents/TSC_ATAC/1749_WT.bigWig")
ATAC_bw_WT3 <- import("~/Documents/TSC_ATAC/1750_WT.bigWig")
ATAC_bw_TSC1 <- import("~/Documents/TSC_ATAC/1747_TSC.bigWig")
ATAC_bw_TSC2 <- import("~/Documents/TSC_ATAC/1748_TSC.bigWig")
ATAC_bw_TSC3 <- import("~/Documents/TSC_ATAC/1755_TSC.bigWig")
ATAC_bw_TFE31 <- import("~/Documents/TSC_ATAC/1741_TFE3.bigWig")
ATAC_bw_TFE32 <- import("~/Documents/TSC_ATAC/1746_TFE3.bigWig")
ATAC_bw_TFE33 <- import("~/Documents/TSC_ATAC/1751_TFE3.bigWig")
ATAC_bw_dKO1 <- import("~/Documents/TSC_ATAC/1708_dKO.bigWig")
ATAC_bw_dKO2 <- import("~/Documents/TSC_ATAC/1744_dKO.bigWig")
ATAC_bw_dKO3 <- import("~/Documents/TSC_ATAC/1754_dKO.bigWig")

library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(ATAC_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT1
mat_WT2<- normalizeToMatrix(ATAC_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT2
mat_WT3<- normalizeToMatrix(ATAC_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_WT3

mat_TSC1<- normalizeToMatrix(ATAC_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC1
mat_TSC2<- normalizeToMatrix(ATAC_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC2
mat_TSC3<- normalizeToMatrix(ATAC_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TSC3

mat_TFE31<- normalizeToMatrix(ATAC_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE31
mat_TFE32<- normalizeToMatrix(ATAC_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE32
mat_TFE33<- normalizeToMatrix(ATAC_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_TFE33

mat_dKO1<- normalizeToMatrix(ATAC_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO1
mat_dKO2<- normalizeToMatrix(ATAC_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO2
mat_dKO3<- normalizeToMatrix(ATAC_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50) / max_dKO3

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3


col_fun<- circlize::colorRamp2(c(0, 1.5), c("white", "red"))

ht_list <- EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/ATAC_TFE3_HDAC_USF_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/ATAC_TFE3_HDAC_USF_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/ATAC_TFE3_HDAC_USF_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()


#######################################################################################
#So the paradox is that TSC KO does not affect chromatin opening at TFE3 promoter sites
#but increases H3K27Ac marks. Knocking out TFE3 dramatically closes chromatin at its own
#sites but does not affect H3K27Ac marks.
#######################################################################################



##############################
######### Figure SXB #########
##############################


#### Show TFE3/TFEB/TFEC/MITF chromatin accessibility

#TFE3 promoter peak is union_peak_118912
#TFEB promoter peak is union_peak_50854
#TFEC promoter is closed
#MITF promoter peak is union_peak_96985

peaks_ATAC_TFE = peaks_ATAC[peaks_ATAC$name %in% c('union_peak_118912','union_peak_50854',
  'union_peak_96985')]
summits_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_summits.bed")
summits_ATAC_TFE = summits_ATAC[summits_ATAC$name %in% c('union_peak_118912','union_peak_50854',
  'union_peak_96985')]


anchors <- summits_ATAC_TFE

mat_WT1<- normalizeToMatrix(ATAC_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(ATAC_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(ATAC_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(ATAC_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(ATAC_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(ATAC_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(ATAC_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(ATAC_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(ATAC_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(ATAC_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(ATAC_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(ATAC_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3


col_fun<- circlize::colorRamp2(c(0, 10), c("white", "red"))

ht_list <- EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/ATAC_MIT_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)



pdf('~/Documents/TSC_Paper/ATAC_TFEB_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(colMeans(mat_WT[1,]), colMeans(mat_TSC[1,]), colMeans(mat_TFE3[1,]), colMeans(mat_dKO[1,])) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()


pdf('~/Documents/TSC_Paper/ATAC_MITF_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(colMeans(mat_WT[2,]), colMeans(mat_TSC[2,]), colMeans(mat_TFE3[2,]), colMeans(mat_dKO[2,])) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

pdf('~/Documents/TSC_Paper/ATAC_TFE3_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(colMeans(mat_WT[3,]), colMeans(mat_TSC[3,]), colMeans(mat_TFE3[3,]), colMeans(mat_dKO[3,])) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()


pdf('~/Documents/TSC_Paper/ATAC_MIT_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/ATAC_MIT_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

#############



##############################
######### Figure SXC #########
##############################


#### Changes in H3K27Ac marks at TFE3 targets

#summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
#summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
#peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
#peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")
peaks_H3K27Ac=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K27Ac_broad_peaks.broadPeak")
summits_H3K27Ac=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K27Ac_broad_summits.bed")

summits_H3K27Ac$name = peaks_H3K27Ac$name

#peaks_USF2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/USF2_peaks.narrowPeak")


#seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
#seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
#seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
#seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"


#peakAnno.edb <- annotatePeak(peaks_TFE3_CHIP_WT, tssRegion=c(-3000, 3000),
#                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

#anno <- as.data.frame(peakAnno.edb)$annotation
#prom <- anno == 'Promoter (<=1kb)'

#peaks_TFE3_CHIP_WT.filt <- peaks_TFE3_CHIP_WT[prom]


#peakAnno.edb <- annotatePeak(summits_TFE3_CHIP_FLCNKO, tssRegion=c(-3000, 3000),
#                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")






#H3K27Ac_overlap_TFE3_peaks<- subsetByOverlaps(peaks_H3K27Ac,peaks_TFE3_CHIP_WT.filt)


extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

H3K4me1<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me1_ENCFF470XUZ.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K4me3<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me3_ENCFF692VTG.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K27ac<- import("~/Documents/TSC_ATAC/ExistingData/H3K27Ac_ENCFF279KZY.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
seqlevelsStyle(H3K4me1) <- "NCBI"
seqlevelsStyle(H3K4me3) <- "NCBI"
seqlevelsStyle(H3K27ac) <- "NCBI"


#enhancers<- subsetByOverlaps(H3K4me1, H3K4me3, invert=TRUE)
#promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

#H3K27Ac_enhancers<- subsetByOverlaps(H3K27Ac_overlap_TFE3_peaks, enhancers) 
#H3K27Ac_promoters<- subsetByOverlaps(H3K27Ac_overlap_TFE3_peaks, promoters) 

#Use ENCODE H3K27Ac peak calls
#H3K27Ac_overlap_TFE3_peaks<- subsetByOverlaps(peaks_H3K27Ac,peaks_TFE3_CHIP_WT.filt)
#H3K27Ac_overlap_TFE3_peaks<- subsetByOverlaps(H3K27Ac_overlap_TFE3_peaks,H3K27ac)
#H3K27Ac_overlap_TFE3_peaks<- subsetByOverlaps(H3K27Ac_overlap_TFE3_peaks,peaks_USF2)


H3K27Ac_keep<- subsetByOverlaps(peaks_H3K27Ac,H3K27ac)

cutoff <- sort(H3K27Ac_keep$score,decreasing=T)[5000]
idx <- H3K27Ac_keep[H3K27Ac_keep$score >= cutoff]$name

#anchors_enhancer <- summits_H3K27Ac[summits_H3K27Ac$name %in% H3K27Ac_enhancers$name]
#anchors_promoter <- summits_H3K27Ac[summits_H3K27Ac$name %in% H3K27Ac_promoters$name]
#anchors <- c(anchors_promoter,anchors_enhancer)
anchors <- summits_H3K27Ac[summits_H3K27Ac$name %in% idx]


H3K27Ac_bw_WT1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1706.spikein_normalized.bw")
H3K27Ac_bw_WT2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1749.spikein_normalized.bw")
H3K27Ac_bw_WT3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1750.spikein_normalized.bw")
H3K27Ac_bw_TSC1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1747.spikein_normalized.bw")
H3K27Ac_bw_TSC2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1748.spikein_normalized.bw")
H3K27Ac_bw_TSC3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1755.spikein_normalized.bw")
H3K27Ac_bw_TFE31 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1741.spikein_normalized.bw")
H3K27Ac_bw_TFE32 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1746.spikein_normalized.bw")
H3K27Ac_bw_TFE33 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1751.spikein_normalized.bw")
H3K27Ac_bw_dKO1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1708.spikein_normalized.bw")
H3K27Ac_bw_dKO2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1744.spikein_normalized.bw")
H3K27Ac_bw_dKO3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1754.spikein_normalized.bw")


library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(H3K27Ac_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(H3K27Ac_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(H3K27Ac_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(H3K27Ac_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(H3K27Ac_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(H3K27Ac_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(H3K27Ac_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(H3K27Ac_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(H3K27Ac_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(H3K27Ac_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(H3K27Ac_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(H3K27Ac_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3

col_fun<- circlize::colorRamp2(c(0, 1), c("white", "red"))


#partition<- c(rep("promoters", length(H3K27Ac_promoters)),
#              rep("enhancers", length(H3K27Ac_enhancers)))

# change the factor level so promoters come first
#partition<- factor(partition, levels=c("promoters", "enhancers"))

#partition_hp<- Heatmap(partition, col=structure(2:3, names = c("promoters", "enhancers")), 
#        name = "partition",
#        show_row_names = FALSE, width=unit(3,'mm'))




ht_list <-  EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/H3K27Ac_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/H3K27Ac_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/H3K27Ac_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



#### Changes in H3K27Ac marks across all TFE3 targets (even if no H3K27Ac peak there)

summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")

scores_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_peaks.bed",format="bed")

peaks_H3K27Ac=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K27Ac_broad_peaks.broadPeak")
summits_H3K27Ac=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K27Ac_broad_summits.bed")

summits_H3K27Ac$name = peaks_H3K27Ac$name



cutoff <- sort(scores_TFE3_CHIP_FLCNKO$score,decreasing=T)[5000]
idx <- scores_TFE3_CHIP_FLCNKO[scores_TFE3_CHIP_FLCNKO$score >= cutoff]$name

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"


extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

H3K4me1<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me1/ENCFF596ORE.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K4me3<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me3/ENCFF150KSB.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K27ac<- import("~/Documents/TSC_ATAC/ExistingData/H3K27Ac/ENCFF395EHX.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
seqlevelsStyle(H3K4me1) <- "NCBI"
seqlevelsStyle(H3K4me3) <- "NCBI"
seqlevelsStyle(H3K27ac) <- "NCBI"


enhancers<- subsetByOverlaps(H3K4me1, H3K4me3, invert=TRUE)
promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

H3K27Ac_enhancers<- subsetByOverlaps(peaks_TFE3_CHIP_FLCNKO, enhancers) 
H3K27Ac_promoters<- subsetByOverlaps(peaks_TFE3_CHIP_FLCNKO, promoters) 


anchors_enhancer <- summits_TFE3_CHIP_FLCNKO[summits_TFE3_CHIP_FLCNKO$name %in% H3K27Ac_enhancers$name]
anchors_promoter <- summits_TFE3_CHIP_FLCNKO[summits_TFE3_CHIP_FLCNKO$name %in% H3K27Ac_promoters$name]
anchors <- c(anchors_promoter,anchors_enhancer)
anchors <- anchors[anchors$name %in% idx]
anchors_enhancer <- anchors_enhancer[anchors_enhancer$name %in% idx]
anchors_promoter <- anchors_promoter[anchors_promoter$name %in% idx]

#anchors <- anchors[anchors$name %in% idx]


H3K27Ac_bw_WT1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1706.spikein_normalized.bw")
H3K27Ac_bw_WT2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1749.spikein_normalized.bw")
H3K27Ac_bw_WT3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1750.spikein_normalized.bw")
H3K27Ac_bw_TSC1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1747.spikein_normalized.bw")
H3K27Ac_bw_TSC2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1748.spikein_normalized.bw")
H3K27Ac_bw_TSC3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1755.spikein_normalized.bw")
H3K27Ac_bw_TFE31 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1741.spikein_normalized.bw")
H3K27Ac_bw_TFE32 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1746.spikein_normalized.bw")
H3K27Ac_bw_TFE33 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1751.spikein_normalized.bw")
H3K27Ac_bw_dKO1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1708.spikein_normalized.bw")
H3K27Ac_bw_dKO2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1744.spikein_normalized.bw")
H3K27Ac_bw_dKO3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1754.spikein_normalized.bw")


library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(H3K27Ac_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(H3K27Ac_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(H3K27Ac_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(H3K27Ac_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(H3K27Ac_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(H3K27Ac_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(H3K27Ac_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(H3K27Ac_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(H3K27Ac_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(H3K27Ac_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(H3K27Ac_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(H3K27Ac_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3

col_fun<- circlize::colorRamp2(c(0, 1), c("white", "red"))


partition<- c(rep("promoters", length(anchors_promoter)),
              rep("enhancers", length(anchors_enhancer)))

# change the factor level so promoters come first
partition<- factor(partition, levels=c("promoters", "enhancers"))

partition_hp<- Heatmap(partition, col=structure(2:3, names = c("promoters", "enhancers")), 
        name = "partition",
        show_row_names = FALSE, width=unit(3,'mm'))




ht_list <- partition_hp + EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/TFE3_anchored_H3K27Ac_TFE3_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, split= partition, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1[partition == 'promoters',])
WT2_mean<- colMeans(mat_WT2[partition == 'promoters',])
WT3_mean<- colMeans(mat_WT3[partition == 'promoters',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'promoters',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'promoters',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'promoters',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'promoters',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'promoters',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'promoters',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'promoters',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'promoters',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'promoters',])

WT_mean<- colMeans(mat_WT[partition == 'promoters',])
TSC_mean<- colMeans(mat_TSC[partition == 'promoters',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'promoters',])
dKO_mean<- colMeans(mat_dKO[partition == 'promoters',])


pdf('~/Documents/TSC_Paper/TFE3_anchored_H3K27Ac_TFE3_promoter_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/TFE3_anchored_H3K27Ac_TFE3_promoter_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



WT1_mean<- colMeans(mat_WT1[partition == 'enhancers',])
WT2_mean<- colMeans(mat_WT2[partition == 'enhancers',])
WT3_mean<- colMeans(mat_WT3[partition == 'enhancers',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'enhancers',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'enhancers',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'enhancers',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'enhancers',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'enhancers',])
#TFE33_mean<- colMeans(mat_TFE33[partition == 'enhancers',])#Bad
dKO1_mean<- colMeans(mat_dKO1[partition == 'enhancers',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'enhancers',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'enhancers',])

mat_TFE3 <- (mat_TFE31 + mat_TFE32) / 2


WT_mean<- colMeans(mat_WT[partition == 'enhancers',])
TSC_mean<- colMeans(mat_TSC[partition == 'enhancers',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'enhancers',])
dKO_mean<- colMeans(mat_dKO[partition == 'enhancers',])


pdf('~/Documents/TSC_Paper/TFE3_anchored_H3K27Ac_TFE3_enhancer_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/TFE3_anchored_H3K27Ac_TFE3_enhancer_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()


#### Changes in H3K27Ac marks

peaks_H3K27Ac=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K27Ac_broad_peaks.broadPeak")
summits_H3K27Ac=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K27Ac_broad_summits.bed")

summits_H3K27Ac$name = peaks_H3K27Ac$name

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"


peakAnno.edb <- annotatePeak(peaks_H3K27Ac, tssRegion=c(-3000, 3000),
                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

#anno <- as.data.frame(peakAnno.edb)$annotation

#SYMBOL
#CLEAR <- read.csv('~/Documents/TSC_Paper/Clear_genes.csv',header=FALSE)$V1


#peaks_H3K27Ac.keep <- subset(peakAnno.edb,distanceToTSS < 1000 & SYMBOL %in% CLEAR)
#anchors <- summits_H3K27Ac[summits_H3K27Ac$name %in% as.data.frame(peaks_H3K27Ac.keep)$name]



#idx <- tail(sort(peaks_H3K27Ac$score),n=5000)[1]
#peaks_H3K27Ac.keep <- peaks_H3K27Ac[peaks_H3K27Ac$score > idx]


extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

H3K4me1<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me1_ENCFF470XUZ.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K4me3<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me3_ENCFF692VTG.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K27ac<- import("~/Documents/TSC_ATAC/ExistingData/H3K27Ac_ENCFF279KZY.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
seqlevelsStyle(H3K4me1) <- "NCBI"
seqlevelsStyle(H3K4me3) <- "NCBI"
seqlevelsStyle(H3K27ac) <- "NCBI"


enhancers<- subsetByOverlaps(H3K4me1, H3K4me3, invert=TRUE)
promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

H3K27Ac_enhancers<- subsetByOverlaps(peaks_H3K27Ac.keep, enhancers) 
H3K27Ac_promoters<- subsetByOverlaps(peaks_H3K27Ac.keep, promoters) 

#Subset on ENCODE H3K27Ac peaks to clean up noise
H3K27Ac_enhancers <- subsetByOverlaps(H3K27Ac_enhancers, H3K27ac) 
H3K27Ac_promoters <- subsetByOverlaps(H3K27Ac_promoters, H3K27ac) 

anchors_enhancer <- summits_H3K27Ac[summits_H3K27Ac$name %in% H3K27Ac_enhancers$name]
anchors_promoter <- summits_H3K27Ac[summits_H3K27Ac$name %in% H3K27Ac_promoters$name]
anchors <- c(anchors_promoter,anchors_enhancer)

H3K27Ac_ENCODE <- subsetByOverlaps(peaks_H3K27Ac, H3K27ac) 
anchors <- summits_H3K27Ac[summits_H3K27Ac$name %in% H3K27Ac_ENCODE$name]


H3K27Ac_bw_WT1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1706.spikein_normalized.bw")
H3K27Ac_bw_WT2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1749.spikein_normalized.bw")
H3K27Ac_bw_WT3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1750.spikein_normalized.bw")
H3K27Ac_bw_TSC1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1747.spikein_normalized.bw")
H3K27Ac_bw_TSC2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1748.spikein_normalized.bw")
H3K27Ac_bw_TSC3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1755.spikein_normalized.bw")
H3K27Ac_bw_TFE31 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1741.spikein_normalized.bw")
H3K27Ac_bw_TFE32 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1746.spikein_normalized.bw")
H3K27Ac_bw_TFE33 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1751.spikein_normalized.bw")
H3K27Ac_bw_dKO1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1708.spikein_normalized.bw")
H3K27Ac_bw_dKO2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1744.spikein_normalized.bw")
H3K27Ac_bw_dKO3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1754.spikein_normalized.bw")


library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(H3K27Ac_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(H3K27Ac_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(H3K27Ac_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(H3K27Ac_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(H3K27Ac_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(H3K27Ac_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(H3K27Ac_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(H3K27Ac_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(H3K27Ac_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(H3K27Ac_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(H3K27Ac_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(H3K27Ac_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)


peak_mat <- cbind(apply(X=mat_TSC1, MARGIN=1, FUN=max),
  apply(X=mat_TSC2, MARGIN=1, FUN=max),
  apply(X=mat_TSC3, MARGIN=1, FUN=max),
  apply(X=mat_dKO1, MARGIN=1, FUN=max),
  apply(X=mat_dKO2, MARGIN=1, FUN=max),
  apply(X=mat_dKO3, MARGIN=1, FUN=max)
)

peak_mat <- cbind(apply(X=mat_WT1, MARGIN=1, FUN=max),
  apply(X=mat_WT2, MARGIN=1, FUN=max),
  apply(X=mat_WT3, MARGIN=1, FUN=max),
  apply(X=mat_TFE31, MARGIN=1, FUN=max),
  apply(X=mat_TFE32, MARGIN=1, FUN=max),
  apply(X=mat_TFE33, MARGIN=1, FUN=max)
)



idx <- peak_mat[rowSums(peak_mat[,1:3] > 1.5*peak_mat[,4:6]) >= 3,]

H3K27Ac_test <- H3K27Ac_ENCODE[rowSums(peak_mat[,1:3] < 2/3*peak_mat[,4:6]) >= 3]


peakAnno.edb <- annotatePeak(H3K27Ac_test, tssRegion=c(-3000, 3000),
                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")


anno <- as.data.frame(peakAnno.edb)$annotation

temp <- subset(peakAnno.edb,annotation == "Promoter (<=1kb)")

cat(unique(as.data.frame(peakAnno.edb)$SYMBOL),sep='\n')


quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3

col_fun<- circlize::colorRamp2(c(0, 1), c("white", "red"))


partition<- c(rep("promoters", length(H3K27Ac_promoters)),
              rep("enhancers", length(H3K27Ac_enhancers)))

# change the factor level so promoters come first
partition<- factor(partition, levels=c("promoters", "enhancers"))

partition_hp<- Heatmap(partition, col=structure(2:3, names = c("promoters", "enhancers")), 
        name = "partition",
        show_row_names = FALSE, width=unit(3,'mm'))




ht_list <- partition_hp + EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/H3K27Ac_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, split= partition, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1[partition == 'promoters',])
WT2_mean<- colMeans(mat_WT2[partition == 'promoters',])
WT3_mean<- colMeans(mat_WT3[partition == 'promoters',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'promoters',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'promoters',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'promoters',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'promoters',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'promoters',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'promoters',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'promoters',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'promoters',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'promoters',])

WT_mean<- colMeans(mat_WT[partition == 'promoters',])
TSC_mean<- colMeans(mat_TSC[partition == 'promoters',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'promoters',])
dKO_mean<- colMeans(mat_dKO[partition == 'promoters',])

#WT1_mean<- matrixStats::colMedians(as.matrix(mat_WT1[partition == 'promoters',]))
#WT2_mean<- matrixStats::colMedians(as.matrix(mat_WT2[partition == 'promoters',]))
#WT3_mean<- matrixStats::colMedians(as.matrix(mat_WT3[partition == 'promoters',]))
#TSC1_mean<- matrixStats::colMedians(as.matrix(mat_TSC1[partition == 'promoters',]))
#TSC2_mean<- matrixStats::colMedians(as.matrix(mat_TSC2[partition == 'promoters',]))
#TSC3_mean<- matrixStats::colMedians(as.matrix(mat_TSC3[partition == 'promoters',]))
#TFE31_mean<- matrixStats::colMedians(as.matrix(mat_TFE31[partition == 'promoters',]))
#TFE32_mean<- matrixStats::colMedians(as.matrix(mat_TFE32[partition == 'promoters',]))
#TFE33_mean<- matrixStats::colMedians(as.matrix(mat_TFE33[partition == 'promoters',]))
#dKO1_mean<- matrixStats::colMedians(as.matrix(mat_dKO1[partition == 'promoters',]))
#dKO2_mean<- matrixStats::colMedians(as.matrix(mat_dKO2[partition == 'promoters',]))
#dKO3_mean<- matrixStats::colMedians(as.matrix(mat_dKO3[partition == 'promoters',]))



pdf('~/Documents/TSC_Paper/H3K27Ac_promoter_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

pdf('~/Documents/TSC_Paper/H3K27Ac_promoter_peaks_ALL_mean.pdf',height=2, width=4)
bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean,TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(pooled = c(rep('WT',120),rep('TSC1',120),rep('WT',120),rep('TSC1',120)),origin = c(rep('WT',120),rep('TSC1',120),rep('TFE3',120),rep('dKO',120)), name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  group_by(name,origin) %>%
  summarize(mean_cl_boot(value,conf.int=0.9)) %>%
  ggplot(aes(x=name, y=y,color = origin, group = origin)) +
  geom_ribbon(aes(fill = origin, ymin = ymin, ymax = ymax), alpha = 0.1, color = NA) +
  geom_line(size=1) +  
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  scale_fill_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

pdf('~/Documents/TSC_Paper/H3K27Ac_promoter_peaks_ALL_mean_CI.pdf',height=3, width=4)
bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean,TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(origin = c(rep('WT',120),rep('TSC1',120),rep('WT',120),rep('TSC1',120)),pooled = c(rep('WT',120),rep('TSC1',120),rep('TFE3',120),rep('dKO',120)), name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  group_by(name,origin) %>%
  summarize(mean_cl_boot(value,conf.int=0.9)) %>%
  ggplot(aes(x=name, y=y,color = origin, group = origin)) +
  geom_ribbon(aes(fill = origin, ymin = ymin, ymax = ymax), alpha = 0.1, color = NA) +
  geom_line(size=1) +  
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  scale_fill_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()


df <- data.frame(WT = rowMeans(mat_WT[partition == 'promoters',c(20:21)]),
  TSC1 = rowMeans(mat_TSC1[partition == 'promoters',c(20:21)]),
  TFE3 = rowMeans(mat_TFE3[partition == 'promoters',c(20:21)]),
  dKO = rowMeans(mat_dKO[partition == 'promoters',c(20:21)]))
idx <- dim(df)[1]

df <- data.frame(unlist(df))
df$group <- c(rep('WT',idx),rep('TSC1',idx),rep('TFE3',idx),rep('dKO',idx))
colnames(df) <- c('val','group')

pdf('~/Documents/TSC_Paper/H3K27Ac_peaks_dist.pdf',height=2, width=4)
ggplot(df, aes(x = val,color=group, fill = group)) +
  geom_density(alpha=0.1,bound=c(0,5),bw=.15) + xlim(0,2.5) + 
  scale_fill_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic()
dev.off()


ggplot(df, aes(x=group,y = val,color=group, fill = group)) +
  geom_violin() + ylim(0,2.5)

  scale_fill_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic()


WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)


#WT_mean<- colMeans(mat_WT)
#TSC_mean<- colMeans(mat_TSC)
#TFE3_mean<- colMeans(mat_TFE3)
#dKO_mean<- colMeans(mat_dKO)


WT1_mean<- matrixStats::colMedians(as.matrix(mat_WT1))
WT2_mean<- matrixStats::colMedians(as.matrix(mat_WT2))
WT3_mean<- matrixStats::colMedians(as.matrix(mat_WT3))
TSC1_mean<- matrixStats::colMedians(as.matrix(mat_TSC1))
TSC2_mean<- matrixStats::colMedians(as.matrix(mat_TSC2))
TSC3_mean<- matrixStats::colMedians(as.matrix(mat_TSC3))
TFE31_mean<- matrixStats::colMedians(as.matrix(mat_TFE31))
TFE32_mean<- matrixStats::colMedians(as.matrix(mat_TFE32))
TFE33_mean<- matrixStats::colMedians(as.matrix(mat_TFE33))
dKO1_mean<- matrixStats::colMedians(as.matrix(mat_dKO1))
dKO2_mean<- matrixStats::colMedians(as.matrix(mat_dKO2))
dKO3_mean<- matrixStats::colMedians(as.matrix(mat_dKO3))


pdf('~/Documents/TSC_Paper/H3K27Ac_peaks_ALL_CLEAR.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

pdf('~/Documents/TSC_Paper/H3K27Ac_peaks_ALL_CLEAR.pdf',height=2, width=4)
bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean,TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(origin = c(rep('WT',120),rep('TSC1',120),rep('TFE3',120),rep('dKO',120)), name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  group_by(name,origin) %>%
  summarize(mean_cl_boot(value)) %>%
  ggplot(aes(x=name, y=y,color = origin, group = origin)) +
  geom_ribbon(aes(fill = origin, ymin = ymin, ymax = ymax), alpha = 0.1, color = NA) +
  geom_line(size=1) +  
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  scale_fill_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()







WT1_mean<- colMeans(mat_WT1[partition == 'promoters',])
WT2_mean<- colMeans(mat_WT2[partition == 'promoters',])
WT3_mean<- colMeans(mat_WT3[partition == 'promoters',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'promoters',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'promoters',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'promoters',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'promoters',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'promoters',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'promoters',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'promoters',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'promoters',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'promoters',])

WT_mean<- colMeans(mat_WT[partition == 'promoters',])
TSC_mean<- colMeans(mat_TSC[partition == 'promoters',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'promoters',])
dKO_mean<- colMeans(mat_dKO[partition == 'promoters',])


pdf('~/Documents/TSC_Paper/H3K27Ac_promoter_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

pdf('~/Documents/TSC_Paper/H3K27Ac_promoter_peaks_ALL_mean.pdf',height=3, width=4)
bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean,TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(origin = c(rep('WT',120),rep('TSC1',120),rep('TFE3',120),rep('dKO',120)), name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  group_by(name,origin) %>%
  summarize(mean_cl_boot(value,conf.int=0.9)) %>%
  ggplot(aes(x=name, y=y,color = origin, group = origin)) +
  geom_ribbon(aes(fill = origin, ymin = ymin, ymax = ymax), alpha = 0.1, color = NA) +
  geom_line(size=1) +  
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  scale_fill_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()





# See what is differentially marked by H3K27ac in TSC KO
filtered.windows_H3K27Ac <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_H3K27Ac.rds')
peaks_USF2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/USF2_peaks.narrowPeak")
summits_USF2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/USF2_summits.bed")



seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"

keep <- overlapsAny(rowRanges(filtered.windows_H3K27Ac), peaks_USF2)
TFE3.filtered.windows.H3K27Ac <- filtered.windows_H3K27Ac[keep,]
keep <- overlapsAny(rowRanges(TFE3.filtered.windows.H3K27Ac), peaks_TFE3_CHIP_WT)
TFE3.filtered.windows.H3K27Ac <- TFE3.filtered.windows.H3K27Ac[!keep,]




#TFE3.filtered.windows.H3K27Ac <- filtered.windows_H3K27Ac


#keep <- overlapsAny(rowRanges(TFE3.filtered.windows.H3K27Ac), anchors_enhancer)
#TFE3.filtered.windows.H3K27Ac.promoter <- TFE3.filtered.windows.H3K27Ac[keep,]

TFE3.filtered.windows.H3K27Ac.promoter <- TFE3.filtered.windows.H3K27Ac


colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
combi <- list(c('WT','TSC'), c('WT','TFE3'), c('WT','dKO'))

logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)
#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/H3K27Ac_enhancer_peaks_WT_vs_TSC.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/H3K27Ac_enhancer_WT_vs_TSC.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))


peaks_int <- subset(peaks,direction == 'down' & FDR < 0.1)

peakAnno.edb <- annotatePeak(peaks_int, tssRegion=c(-3000, 3000),
                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

anno <- as.data.frame(peakAnno.edb)$annotation

pdf('~/Documents/TSC_Paper/H3K27Ac_WT_vs_TSC_UP_anno_bar.pdf',width=3,height=1)
plotAnnoBar(peakAnno.edb)
dev.off()

TSC.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TSC.vs.WT.csv',sep=',')


idx <- match(peaks$gene,TSC.vs.WT$gene)

peaks$RNA_log2fc <- -TSC.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TSC.vs.WT$padj[idx]
df <- as.data.frame(peaks)

diff_open <- subset(df,direction == 'down' & FDR < 0.1 & RNA_p_adj < 0.1 & RNA_log2fc < -0.1)$gene




library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)

dbs <- c("ChEA_2022","Reactome_2024","GO_Biological_Process_2023")


enriched <- enrichr(genes_all, dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.05,]
p1<- ggplot(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),], (aes(x=-log10(Adjusted.P.value), y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + geom_point() + xlab('Adjusted P Value') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + ggtitle('ChEA') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"\\s+"), `[`, 1))) + theme(axis.text=element_text(colour="black"))


pdf('~/Documents/TSC_Paper/H3K27Ac_TFE3_Targets_WT_vs_TSC_UP.pdf',width=4.75,height=4)
p1
dev.off()


dbs <- c("ChEA_2022","Reactome_2024","GO_Biological_Process_2023")


enriched <- enrichr(diff_open, dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.05,]
p1<- ggplot(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),], (aes(x=-log10(Adjusted.P.value), y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + geom_point() + xlab('Adjusted P Value') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + ggtitle('ChEA') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"\\s+"), `[`, 1))) + theme(axis.text=element_text(colour="black"))


pdf('~/Documents/TSC_Paper/H3K27Ac_RNASeq_Targets_WT_vs_TSC_UP.pdf',width=4.75,height=4)
p1
dev.off()


# See what is differentially marked by H3K27ac in TSC KO vs dKO


#keep <- overlapsAny(rowRanges(filtered.windows_H3K27Ac), summits_TFE3_CHIP_FLCNKO)
TFE3.filtered.windows.H3K27Ac <- filtered.windows_H3K27Ac

#keep <- overlapsAny(rowRanges(TFE3.filtered.windows.H3K27Ac), anchors_enhancer)
#TFE3.filtered.windows.H3K27Ac.promoter <- TFE3.filtered.windows.H3K27Ac[keep,]

TFE3.filtered.windows.H3K27Ac.promoter <- TFE3.filtered.windows.H3K27Ac


colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(TSC - dKO, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.dKO.vs.TSC <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['dKO']],y = .data[['TSC']] - .data[['dKO']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.dKO.vs.TSC), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'TSC', 'dKO')) +
    theme_bw(base_size = 15)
#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .7, 'TSC',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .7, 'dKO', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/H3K27Ac_enhancer_peaks_TFE3_vs_dKO.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/H3K27Ac_enhancer_TFE3_vs_dKO.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.3 & direction == 'up' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))


library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)

dbs <- c("ChEA_2022","Reactome_2024","GO_Biological_Process_2023")


enriched <- enrichr(genes_all, dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.05,]
p1<- ggplot(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),], (aes(x=-log10(Adjusted.P.value), y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + geom_point() + xlab('Adjusted P Value') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + ggtitle('ChEA') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"\\s+"), `[`, 1))) + theme(axis.text=element_text(colour="black"))


pdf('~/Documents/TSC_Paper/H3K27Ac_TFE3_Targets_TFE3_vs_dKO_DOWN.pdf',width=4.75,height=4)
p1
dev.off()

# See what is differentially marked by H3K27ac in TSC KO vs dKO


keep <- overlapsAny(rowRanges(filtered.windows_H3K27Ac), summits_TFE3_CHIP_FLCNKO)
TFE3.filtered.windows.H3K27Ac <- filtered.windows_H3K27Ac[keep,]

#keep <- overlapsAny(rowRanges(TFE3.filtered.windows.H3K27Ac), anchors_enhancer)
#TFE3.filtered.windows.H3K27Ac.promoter <- TFE3.filtered.windows.H3K27Ac[keep,]

TFE3.filtered.windows.H3K27Ac.promoter <- TFE3.filtered.windows.H3K27Ac


colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(TSC - dKO, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.5)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .5, 'TSC',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .5, 'dKO', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/H3K27Ac_enhancer_peaks_TSC_vs_dKO.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/H3K27Ac_enhancer_TSC_vs_dKO.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.7,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.7 & direction == 'up' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))


library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)

dbs <- c("ChEA_2022","Reactome_2024","GO_Biological_Process_2023")


enriched <- enrichr(genes_all, dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.05,]
p1<- ggplot(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),], (aes(x=-log10(Adjusted.P.value), y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + geom_point() + xlab('Adjusted P Value') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + ggtitle('ChEA') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"\\s+"), `[`, 1))) + theme(axis.text=element_text(colour="black"))


pdf('~/Documents/TSC_Paper/H3K27Ac_TFE3_Targets_TSC_vs_dKO_UP.pdf',width=4.75,height=4)
p1
dev.off()


# See which peaks are globally differential H3K27ac in TSC KO


#keep <- overlapsAny(rowRanges(filtered.windows_H3K27Ac), c(enhancers,promoters))
#TFE3.filtered.windows.H3K27Ac.promoter <- filtered.windows_H3K27Ac[keep,]

TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_H3K27Ac

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)

pdf('~/Documents/TSC_Paper/H3K27Ac_DOPs_TSC_MA.pdf',height=5, width=3.5)
ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_point(col='gray',alpha=0.2) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red',alpha=0.9) ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15) + coord_flip()
dev.off()

#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/H3K27Ac_global_peaks_WT_vs_TSC.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/H3K27Ac_global_WT_vs_TSC.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

#MITF
MITF <- "MIDN;COMT;FADS2;ZMIZ1;PDK4;MACROD1;CCNL1;TBC1D10A;CARD10;ST6GAL1;ANKRD44;ACOT11;DYRK1A;PDIA4;PPP1R3B;MAD1L1;FKBP5;NRN1;EPB41;CTBP1;BRI3;UBC;EIF4EBP1;ATP6V1H;SIVA1;RBPMS;HPS4;HMGA2;GNG12;FBXO32;SSH1;GLUD1;NR4A1;NPLOC4;NR4A3;SP1;BHLHE40;UBE2N;PIR;OSBPL1A;BCAR3;COL18A1;CLIC4;CDKN1A;DOCK4;PLEKHF1;ITGB5;CHD9;BHLHE41;CHD6;TMEM51;SHB;AP2A2;PITPNC1;FHAD1;TUBA1C;HNF4A;MYO18A;EFHD2;TEAD1;MYBL1;CTNNBL1;TRIM62;ITPK1;ARHGEF19;NRG2;SLC39A14;KLF15;NCOR2;RNF144B;PCCA;ARHGEF3;RARA;SFT2D2;MGLL;COLEC12;HDAC5;CEBPB;HLF;SYNPO2;LRP5;SLC7A2;ATOH8;LRIG1;MGAT1;HES1;RALY;LRIG2;GABBR2;KCNJ10;PLEKHA5;KLF3;SULF2;TPCN1;CCS;WEE1;ETNK2;NFIC;RGS12;LTBR;GPD1L;BCAR1"
MITF <- strsplit(MITF,';')[[1]]



library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)

dbs <- c("ChEA_2022","Reactome_2024","GO_Biological_Process_2023")


enriched <- enrichr(genes_all, dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.05,]
p1<- ggplot(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),], (aes(x=-log10(Adjusted.P.value), y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + geom_point() + xlab('Adjusted P Value') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + ggtitle('ChEA') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"\\s+"), `[`, 1))) + theme(axis.text=element_text(colour="black"))


pdf('~/Documents/TSC_Paper/H3K27Ac_Global_TSC_vs_WT_UP.pdf',width=4.75,height=4)
p1
dev.off()

######## Extract H3K27Ac peaks which are DOPs

DOWs <- TFE3.filtered.windows.H3K27Ac.promoter[rowData(TFE3.filtered.windows.H3K27Ac.promoter)$FDR < 0.1 & rowData(TFE3.filtered.windows.H3K27Ac.promoter)$logFC < 0]
DOPs <- subsetByOverlaps(peaks_H3K27Ac, DOWs)

#anchors <- summits_TFE3_CHIP_FLCNKO[summits_TFE3_CHIP_FLCNKO$name %in% DOPs$name]
anchors <- summits_H3K27Ac[summits_H3K27Ac$name %in% DOPs$name]


library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(H3K27Ac_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(H3K27Ac_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(H3K27Ac_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(H3K27Ac_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(H3K27Ac_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(H3K27Ac_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(H3K27Ac_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(H3K27Ac_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(H3K27Ac_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(H3K27Ac_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(H3K27Ac_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(H3K27Ac_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3


col_fun<- circlize::colorRamp2(c(0, 1.5), c("white", "red"))

ht_list <- EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/H3K27Ac_DOPs_TSC.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/H3K27Ac_DOPs_TSC_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/H3K27Ac_DOPs_TSC_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()











##############################
######### Figure SXD #########
##############################


#### Changes in H3K4me3 marks at TFE3 targets

summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")
peaks_H3K4me3=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K4me3_peaks.narrowPeak")
summits_H3K4me3=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K4me3_summits.bed")


cutoff <- sort(peaks_H3K4me3$score,decreasing=T)[5000]
idx <- peaks_H3K4me3[peaks_H3K4me3$score >= cutoff]$name

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"


H3K4me3_overlap_TFE3_peaks<- subsetByOverlaps(peaks_H3K4me3,peaks_TFE3_CHIP_FLCNKO)


extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

H3K4me1<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me1/ENCFF596ORE.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K4me3<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me3/ENCFF150KSB.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K27ac<- import("~/Documents/TSC_ATAC/ExistingData/H3K27Ac/ENCFF395EHX.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
seqlevelsStyle(H3K4me1) <- "NCBI"
seqlevelsStyle(H3K4me3) <- "NCBI"
seqlevelsStyle(H3K27ac) <- "NCBI"


promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

H3K4me3_promoters<- subsetByOverlaps(H3K4me3_overlap_TFE3_peaks, promoters) 

anchors <- summits_H3K4me3[summits_H3K4me3$name %in% H3K4me3_promoters$name]



H3K4me3_bw_WT1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1706-H3K4me3.spikein_normalized.bw")
#H3K4me3_bw_WT2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1749-H3K4me3.spikein_normalized.bw")
H3K4me3_bw_WT3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1750-H3K4me3.spikein_normalized.bw")
H3K4me3_bw_TSC1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1747-H3K4me3.spikein_normalized.bw")
H3K4me3_bw_TSC2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1748-H3K4me3.spikein_normalized.bw")
H3K4me3_bw_TSC3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1755-H3K4me3.spikein_normalized.bw")
H3K4me3_bw_TFE31 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1741-H3K4me3.spikein_normalized.bw")
H3K4me3_bw_TFE32 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1746-H3K4me3.spikein_normalized.bw")
H3K4me3_bw_TFE33 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1751-H3K4me3.spikein_normalized.bw")
H3K4me3_bw_dKO1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1708-H3K4me3.spikein_normalized.bw")
H3K4me3_bw_dKO2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1744-H3K4me3.spikein_normalized.bw")
H3K4me3_bw_dKO3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1754-H3K4me3.spikein_normalized.bw")


library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(H3K4me3_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
#mat_WT2<- normalizeToMatrix(ATAC_bw_WT2, anchors, value_column = "score",
#  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(H3K4me3_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(H3K4me3_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(H3K4me3_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(H3K4me3_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(H3K4me3_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(H3K4me3_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(H3K4me3_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(H3K4me3_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(H3K4me3_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(H3K4me3_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
#quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT3) / 2
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3


col_fun<- circlize::colorRamp2(c(0, 6), c("white", "red"))

ht_list <- EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/H3K4me3_TFE3_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/H3K4me3_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5","#DA9195", "#DA9195")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/H3K4me3_ALL_mean.pdf',height=2, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()


####JUMp#####



# See which promoters are differentially marked by H3K4me3 in TSC KO


keep <- overlapsAny(rowRanges(filtered.windows_H3K4me3), summits_TFE3_CHIP_FLCNKO)
TFE3.filtered.windows.H3K4me3 <- filtered.windows_H3K4me3[keep,]

keep <- overlapsAny(rowRanges(TFE3.filtered.windows.H3K4me3), promoters)
TFE3.filtered.windows.H3K4me3.promoter <- TFE3.filtered.windows.H3K4me3[keep,]

colData(TFE3.filtered.windows.H3K4me3.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K4me3.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K4me3.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.H3K4me3.promoter)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K4me3.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K4me3.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K4me3.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K4me3.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K4me3.promoter)),
                        group = colData(TFE3.filtered.windows.H3K4me3.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K4me3.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K4me3.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)
#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K4me3.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K4me3.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K4me3.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/H3K4me3_promoter_peaks_WT_vs_TSC.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/H3K4me3_promoter_WT_vs_TSC.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.3 & direction == 'up' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

# See differential peaks by H3K4me3 in TSC KO


TFE3.filtered.windows.H3K4me3.promoter <- filtered.windows_H3K4me3

colData(TFE3.filtered.windows.H3K4me3.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K4me3.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K4me3.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.H3K4me3.promoter)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K4me3.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K4me3.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K4me3.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K4me3.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K4me3.promoter)),
                        group = colData(TFE3.filtered.windows.H3K4me3.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K4me3.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K4me3.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)
#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K4me3.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K4me3.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K4me3.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/H3K4me3_all_peaks_WT_vs_TSC.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/H3K4me3_all_WT_vs_TSC.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))


##############################
######### Figure SXE #########
##############################

#### Changes in HDAC2 binding at TFE3 targets

summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")
peaks_HDAC2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC2_peaks.narrowPeak")
summits_HDAC2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC2_summits.bed")

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"



#cutoff <- sort(summits_HDAC2$score,decreasing=T)[5000]
#idx <- summits_HDAC2[summits_HDAC2$score >= cutoff]$name

HDAC2_overlap_TFE3_peaks<- subsetByOverlaps(peaks_HDAC2,peaks_TFE3_CHIP_FLCNKO)


extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

H3K4me1<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me1/ENCFF596ORE.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K4me3<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me3/ENCFF150KSB.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K27ac<- import("~/Documents/TSC_ATAC/ExistingData/H3K27Ac/ENCFF395EHX.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
seqlevelsStyle(H3K4me1) <- "NCBI"
seqlevelsStyle(H3K4me3) <- "NCBI"
seqlevelsStyle(H3K27ac) <- "NCBI"


enhancers<- subsetByOverlaps(H3K4me1, H3K4me3, invert=TRUE)
promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

HDAC2_enhancers<- subsetByOverlaps(HDAC2_overlap_TFE3_peaks, enhancers) 
HDAC2_promoters<- subsetByOverlaps(HDAC2_overlap_TFE3_peaks, promoters) 

anchors_enhancer <- summits_HDAC2[summits_HDAC2$name %in% HDAC2_enhancers$name]
anchors_promoter <- summits_HDAC2[summits_HDAC2$name %in% HDAC2_promoters$name]
anchors <- c(anchors_promoter,anchors_enhancer)
#anchors <- anchors[anchors$name %in% idx]




HDAC2_bw_WT1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1706.spikein_normalized.bw")
HDAC2_bw_WT2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1749.spikein_normalized.bw")
HDAC2_bw_WT3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1750.spikein_normalized.bw")
HDAC2_bw_TSC1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1747.spikein_normalized.bw")
HDAC2_bw_TSC2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1748.spikein_normalized.bw")
HDAC2_bw_TSC3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1755.spikein_normalized.bw")
HDAC2_bw_TFE31 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC21741.spikein_normalized.bw")
HDAC2_bw_TFE32 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1746.spikein_normalized.bw")
HDAC2_bw_TFE33 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1751.spikein_normalized.bw")
HDAC2_bw_dKO1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1708.spikein_normalized.bw")
HDAC2_bw_dKO2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1744.spikein_normalized.bw")
HDAC2_bw_dKO3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1754.spikein_normalized.bw")


library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(HDAC2_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(HDAC2_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(HDAC2_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(HDAC2_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(HDAC2_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(HDAC2_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(HDAC2_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(HDAC2_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(HDAC2_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(HDAC2_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(HDAC2_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(HDAC2_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3

col_fun<- circlize::colorRamp2(c(0, 1), c("white", "red"))


partition<- c(rep("promoters", length(HDAC2_promoters)),
              rep("enhancers", length(HDAC2_enhancers)))

# change the factor level so promoters come first
partition<- factor(partition, levels=c("promoters", "enhancers"))

partition_hp<- Heatmap(partition, col=structure(2:3, names = c("promoters", "enhancers")), 
        name = "partition",
        show_row_names = FALSE, width=unit(3,'mm'))




ht_list <- partition_hp + EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/HDAC2_TFE3_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, split= partition, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1[partition == 'promoters',])
WT2_mean<- colMeans(mat_WT2[partition == 'promoters',])
WT3_mean<- colMeans(mat_WT3[partition == 'promoters',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'promoters',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'promoters',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'promoters',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'promoters',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'promoters',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'promoters',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'promoters',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'promoters',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'promoters',])

WT_mean<- colMeans(mat_WT[partition == 'promoters',])
TSC_mean<- colMeans(mat_TSC[partition == 'promoters',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'promoters',])
dKO_mean<- colMeans(mat_dKO[partition == 'promoters',])


pdf('~/Documents/TSC_Paper/HDAC2_TFE3_promoter_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/HDAC2_TFE3_promoter_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



WT1_mean<- colMeans(mat_WT1[partition == 'enhancers',])
WT2_mean<- colMeans(mat_WT2[partition == 'enhancers',])
WT3_mean<- colMeans(mat_WT3[partition == 'enhancers',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'enhancers',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'enhancers',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'enhancers',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'enhancers',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'enhancers',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'enhancers',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'enhancers',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'enhancers',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'enhancers',])



WT_mean<- colMeans(mat_WT[partition == 'enhancers',])
TSC_mean<- colMeans(mat_TSC[partition == 'enhancers',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'enhancers',])
dKO_mean<- colMeans(mat_dKO[partition == 'enhancers',])


pdf('~/Documents/TSC_Paper/HDAC2_TFE3_enhancer_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/HDAC2_TFE3_enhancer_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

##############################
######### Figure SXF #########
##############################

#### Changes in HDAC1 binding at TFE3 targets

summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")
peaks_HDAC1=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC1_peaks.narrowPeak")
summits_HDAC1=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC1_summits.bed")

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"



#cutoff <- sort(summits_HDAC2$score,decreasing=T)[5000]
#idx <- summits_HDAC2[summits_HDAC2$score >= cutoff]$name

HDAC1_overlap_TFE3_peaks<- subsetByOverlaps(peaks_HDAC1,peaks_TFE3_CHIP_FLCNKO)


extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

H3K4me1<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me1/ENCFF596ORE.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K4me3<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me3/ENCFF150KSB.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K27ac<- import("~/Documents/TSC_ATAC/ExistingData/H3K27Ac/ENCFF395EHX.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
seqlevelsStyle(H3K4me1) <- "NCBI"
seqlevelsStyle(H3K4me3) <- "NCBI"
seqlevelsStyle(H3K27ac) <- "NCBI"


enhancers<- subsetByOverlaps(H3K4me1, H3K4me3, invert=TRUE)
promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

HDAC1_enhancers<- subsetByOverlaps(HDAC1_overlap_TFE3_peaks, enhancers) 
HDAC1_promoters<- subsetByOverlaps(HDAC1_overlap_TFE3_peaks, promoters) 

anchors_enhancer <- summits_HDAC1[summits_HDAC1$name %in% HDAC1_enhancers$name]
anchors_promoter <- summits_HDAC1[summits_HDAC1$name %in% HDAC1_promoters$name]
anchors <- c(anchors_promoter,anchors_enhancer)
#anchors <- anchors[anchors$name %in% idx]




HDAC1_bw_WT1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1706-HDAC1.spikein_normalized.bw")
HDAC1_bw_WT2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1749-HDAC1.spikein_normalized.bw")
HDAC1_bw_WT3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1750-HDAC1.spikein_normalized.bw")
HDAC1_bw_TSC1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1747-HDAC1.spikein_normalized.bw")
HDAC1_bw_TSC2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1748-HDAC1.spikein_normalized.bw")
HDAC1_bw_TSC3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1755-HDAC1.spikein_normalized.bw")
HDAC1_bw_TFE31 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1741-HDAC1.spikein_normalized.bw")
HDAC1_bw_TFE32 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1746-HDAC1.spikein_normalized.bw")
HDAC1_bw_TFE33 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1751-HDAC1.spikein_normalized.bw")
HDAC1_bw_dKO1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1708-HDAC1.spikein_normalized.bw")
HDAC1_bw_dKO2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1744-HDAC1.spikein_normalized.bw")
HDAC1_bw_dKO3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1754-HDAC1.spikein_normalized.bw")


library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(HDAC1_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(HDAC1_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(HDAC1_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(HDAC1_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(HDAC1_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(HDAC1_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(HDAC1_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(HDAC1_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(HDAC1_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(HDAC1_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(HDAC1_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(HDAC1_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3

col_fun<- circlize::colorRamp2(c(0, 1.25), c("white", "red"))


partition<- c(rep("promoters", length(HDAC1_promoters)),
              rep("enhancers", length(HDAC1_enhancers)))

# change the factor level so promoters come first
partition<- factor(partition, levels=c("promoters", "enhancers"))

partition_hp<- Heatmap(partition, col=structure(2:3, names = c("promoters", "enhancers")), 
        name = "partition",
        show_row_names = FALSE, width=unit(3,'mm'))




ht_list <- partition_hp + EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/HDAC1_TFE3_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, split= partition, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1[partition == 'promoters',])
WT2_mean<- colMeans(mat_WT2[partition == 'promoters',])
WT3_mean<- colMeans(mat_WT3[partition == 'promoters',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'promoters',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'promoters',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'promoters',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'promoters',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'promoters',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'promoters',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'promoters',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'promoters',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'promoters',])

WT_mean<- colMeans(mat_WT[partition == 'promoters',])
TSC_mean<- colMeans(mat_TSC[partition == 'promoters',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'promoters',])
dKO_mean<- colMeans(mat_dKO[partition == 'promoters',])


pdf('~/Documents/TSC_Paper/HDAC1_TFE3_promoter_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/HDAC1_TFE3_promoter_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



WT1_mean<- colMeans(mat_WT1[partition == 'enhancers',])
WT2_mean<- colMeans(mat_WT2[partition == 'enhancers',])
WT3_mean<- colMeans(mat_WT3[partition == 'enhancers',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'enhancers',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'enhancers',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'enhancers',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'enhancers',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'enhancers',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'enhancers',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'enhancers',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'enhancers',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'enhancers',])



WT_mean<- colMeans(mat_WT[partition == 'enhancers',])
TSC_mean<- colMeans(mat_TSC[partition == 'enhancers',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'enhancers',])
dKO_mean<- colMeans(mat_dKO[partition == 'enhancers',])


pdf('~/Documents/TSC_Paper/HDAC1_TFE3_enhancer_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/HDAC1_TFE3_enhancer_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()


##############################
######### Figure SXF #########
##############################

#### Changes in TFE3 binding at TFE3 targets

summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")
peaks_TFE3=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/TFE3_peaks.narrowPeak")
summits_TFE3=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/TFE3_summits.bed")

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"



#cutoff <- sort(summits_HDAC2$score,decreasing=T)[5000]
#idx <- summits_HDAC2[summits_HDAC2$score >= cutoff]$name

TFE3_overlap_TFE3_peaks<- subsetByOverlaps(peaks_TFE3,peaks_TFE3_CHIP_FLCNKO)


extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

H3K4me1<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me1/ENCFF596ORE.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K4me3<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me3/ENCFF150KSB.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K27ac<- import("~/Documents/TSC_ATAC/ExistingData/H3K27Ac/ENCFF395EHX.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
seqlevelsStyle(H3K4me1) <- "NCBI"
seqlevelsStyle(H3K4me3) <- "NCBI"
seqlevelsStyle(H3K27ac) <- "NCBI"


enhancers<- subsetByOverlaps(H3K4me1, H3K4me3, invert=TRUE)
promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

TFE3_enhancers<- subsetByOverlaps(TFE3_overlap_TFE3_peaks, enhancers) 
TFE3_promoters<- subsetByOverlaps(TFE3_overlap_TFE3_peaks, promoters) 

anchors_enhancer <- summits_TFE3[summits_TFE3$name %in% TFE3_enhancers$name]
anchors_promoter <- summits_TFE3[summits_TFE3$name %in% TFE3_promoters$name]
anchors <- c(anchors_promoter,anchors_enhancer)
#anchors <- anchors[anchors$name %in% idx]




TFE3_bw_WT1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/TFE3s1706.spikein_normalized.bw")
TFE3_bw_WT2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/TFE3s1749.spikein_normalized.bw")
TFE3_bw_WT3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/TFE3s1750.spikein_normalized.bw")
TFE3_bw_TSC1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/TFE3s1747.spikein_normalized.bw")
TFE3_bw_TSC2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/TFE3s1748.spikein_normalized.bw")
TFE3_bw_TSC3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/TFE3s1755.spikein_normalized.bw")
TFE3_bw_TFE31 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/TFE3s1741.spikein_normalized.bw")
TFE3_bw_TFE32 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/TFE3s1746.spikein_normalized.bw")
TFE3_bw_TFE33 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/TFE3s1751.spikein_normalized.bw")
TFE3_bw_dKO1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/TFE3s1708.spikein_normalized.bw")
TFE3_bw_dKO2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/TFE3s1744.spikein_normalized.bw")
TFE3_bw_dKO3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/TFE3s1754.spikein_normalized.bw")



library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(TFE3_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(TFE3_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(TFE3_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(TFE3_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(TFE3_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(TFE3_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(TFE3_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(TFE3_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(TFE3_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(TFE3_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(TFE3_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(TFE3_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3

col_fun<- circlize::colorRamp2(c(0, 1.5), c("white", "red"))


partition<- c(rep("promoters", length(TFE3_promoters)),
              rep("enhancers", length(TFE3_enhancers)))

# change the factor level so promoters come first
partition<- factor(partition, levels=c("promoters", "enhancers"))

partition_hp<- Heatmap(partition, col=structure(2:3, names = c("promoters", "enhancers")), 
        name = "partition",
        show_row_names = FALSE, width=unit(3,'mm'))




ht_list <- partition_hp + EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/TFE3_TFE3_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, split= partition, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1[partition == 'promoters',])
WT2_mean<- colMeans(mat_WT2[partition == 'promoters',])
WT3_mean<- colMeans(mat_WT3[partition == 'promoters',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'promoters',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'promoters',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'promoters',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'promoters',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'promoters',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'promoters',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'promoters',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'promoters',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'promoters',])

WT_mean<- colMeans(mat_WT[partition == 'promoters',])
TSC_mean<- colMeans(mat_TSC[partition == 'promoters',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'promoters',])
dKO_mean<- colMeans(mat_dKO[partition == 'promoters',])


pdf('~/Documents/TSC_Paper/TFE3_TFE3_promoter_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/TFE3_TFE3_promoter_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



WT1_mean<- colMeans(mat_WT1[partition == 'enhancers',])
WT2_mean<- colMeans(mat_WT2[partition == 'enhancers',])
WT3_mean<- colMeans(mat_WT3[partition == 'enhancers',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'enhancers',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'enhancers',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'enhancers',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'enhancers',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'enhancers',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'enhancers',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'enhancers',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'enhancers',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'enhancers',])



WT_mean<- colMeans(mat_WT[partition == 'enhancers',])
TSC_mean<- colMeans(mat_TSC[partition == 'enhancers',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'enhancers',])
dKO_mean<- colMeans(mat_dKO[partition == 'enhancers',])


pdf('~/Documents/TSC_Paper/TFE3_TFE3_enhancer_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/TFE3_TFE3_enhancer_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()


##############################
######### Figure SXG #########
##############################

#### Changes in USF2 binding at TFE3 targets

summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")
peaks_USF2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/USF2_peaks.narrowPeak")
summits_USF2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/USF2_summits.bed")

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"



#cutoff <- sort(summits_HDAC2$score,decreasing=T)[5000]
#idx <- summits_HDAC2[summits_HDAC2$score >= cutoff]$name

USF2_overlap_TFE3_peaks<- subsetByOverlaps(peaks_USF2,peaks_TFE3_CHIP_FLCNKO)


extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

H3K4me1<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me1/ENCFF596ORE.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K4me3<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me3/ENCFF150KSB.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
H3K27ac<- import("~/Documents/TSC_ATAC/ExistingData/H3K27Ac/ENCFF395EHX.bed",
  format = "BED",extraCols = extraCols_narrowPeak)
seqlevelsStyle(H3K4me1) <- "NCBI"
seqlevelsStyle(H3K4me3) <- "NCBI"
seqlevelsStyle(H3K27ac) <- "NCBI"


enhancers<- subsetByOverlaps(H3K4me1, H3K4me3, invert=TRUE)
promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

USF2_enhancers<- subsetByOverlaps(USF2_overlap_TFE3_peaks, enhancers) 
USF2_promoters<- subsetByOverlaps(USF2_overlap_TFE3_peaks, promoters) 

anchors_enhancer <- summits_USF2[summits_USF2$name %in% USF2_enhancers$name]
anchors_promoter <- summits_USF2[summits_USF2$name %in% USF2_promoters$name]
anchors <- c(anchors_promoter,anchors_enhancer)
#anchors <- anchors[anchors$name %in% idx]




USF2_bw_WT1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1706.spikein_normalized.bw")
USF2_bw_WT2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1749.spikein_normalized.bw")
USF2_bw_WT3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1750.spikein_normalized.bw")
USF2_bw_TSC1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1747.spikein_normalized.bw")
USF2_bw_TSC2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1748.spikein_normalized.bw")
USF2_bw_TSC3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1755.spikein_normalized.bw")
USF2_bw_TFE31 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1741.spikein_normalized.bw")
USF2_bw_TFE32 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1746.spikein_normalized.bw")
USF2_bw_TFE33 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1751.spikein_normalized.bw")
USF2_bw_dKO1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1708.spikein_normalized.bw")
USF2_bw_dKO2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1744.spikein_normalized.bw")
USF2_bw_dKO3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1754.spikein_normalized.bw")



library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(USF2_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(USF2_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(USF2_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(USF2_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(USF2_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(USF2_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(USF2_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(USF2_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(USF2_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(USF2_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(USF2_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(USF2_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3

col_fun<- circlize::colorRamp2(c(0, 1), c("white", "red"))


partition<- c(rep("promoters", length(USF2_promoters)),
              rep("enhancers", length(USF2_enhancers)))

# change the factor level so promoters come first
partition<- factor(partition, levels=c("promoters", "enhancers"))

partition_hp<- Heatmap(partition, col=structure(2:3, names = c("promoters", "enhancers")), 
        name = "partition",
        show_row_names = FALSE, width=unit(3,'mm'))




ht_list <- partition_hp + EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/USF2_TFE3_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list, split= partition, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1[partition == 'promoters',])
WT2_mean<- colMeans(mat_WT2[partition == 'promoters',])
WT3_mean<- colMeans(mat_WT3[partition == 'promoters',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'promoters',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'promoters',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'promoters',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'promoters',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'promoters',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'promoters',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'promoters',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'promoters',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'promoters',])

WT_mean<- colMeans(mat_WT[partition == 'promoters',])
TSC_mean<- colMeans(mat_TSC[partition == 'promoters',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'promoters',])
dKO_mean<- colMeans(mat_dKO[partition == 'promoters',])


pdf('~/Documents/TSC_Paper/USF2_TFE3_promoter_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/USF2_TFE3_promoter_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



WT1_mean<- colMeans(mat_WT1[partition == 'enhancers',])
WT2_mean<- colMeans(mat_WT2[partition == 'enhancers',])
WT3_mean<- colMeans(mat_WT3[partition == 'enhancers',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'enhancers',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'enhancers',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'enhancers',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'enhancers',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'enhancers',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'enhancers',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'enhancers',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'enhancers',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'enhancers',])



WT_mean<- colMeans(mat_WT[partition == 'enhancers',])
TSC_mean<- colMeans(mat_TSC[partition == 'enhancers',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'enhancers',])
dKO_mean<- colMeans(mat_dKO[partition == 'enhancers',])


pdf('~/Documents/TSC_Paper/USF2_TFE3_enhancer_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/USF2_TFE3_enhancer_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()

# Look at USF2 binding more globally


anchors <- summits_USF2[summits_USF2$name %in% USF2_overlap_TFE3_peaks$name]

library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(USF2_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(USF2_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(USF2_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(USF2_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(USF2_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(USF2_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(USF2_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(USF2_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(USF2_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(USF2_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(USF2_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(USF2_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3

col_fun<- circlize::colorRamp2(c(0, 1), c("white", "red"))


ht_list <-  EnrichedHeatmap(mat_WT, pos_line = FALSE, column_title="WT", name = "WT", col=col_fun) + 
EnrichedHeatmap(mat_TSC, pos_line = FALSE, column_title="TSC", name = "TSC", col=col_fun) + 
EnrichedHeatmap(mat_TFE3, pos_line = FALSE, column_title="TFE3", name = "TFE3", col=col_fun) + 
EnrichedHeatmap(mat_dKO, pos_line = FALSE, column_title="dKO", name = "dKO", col=col_fun) 

pdf('~/Documents/TSC_Paper/USF2_nonpartition_TFE3_heatmap_ALL.pdf',height=4, width=5)
  draw(ht_list,, main_heatmap =2)
dev.off()

WT1_mean<- colMeans(mat_WT1)
WT2_mean<- colMeans(mat_WT2)
WT3_mean<- colMeans(mat_WT3)
TSC1_mean<- colMeans(mat_TSC1)
TSC2_mean<- colMeans(mat_TSC2)
TSC3_mean<- colMeans(mat_TSC3)
TFE31_mean<- colMeans(mat_TFE31)
TFE32_mean<- colMeans(mat_TFE32)
TFE33_mean<- colMeans(mat_TFE33)
dKO1_mean<- colMeans(mat_dKO1)
dKO2_mean<- colMeans(mat_dKO2)
dKO3_mean<- colMeans(mat_dKO3)

WT_mean<- colMeans(mat_WT)
TSC_mean<- colMeans(mat_TSC)
TFE3_mean<- colMeans(mat_TFE3)
dKO_mean<- colMeans(mat_dKO)


pdf('~/Documents/TSC_Paper/USF2_nonpartition_TFE3_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/USF2_nonpartition_TFE3_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()


# See which promoters are differentially marked by USF2 in TFE3 KO


#keep <- overlapsAny(rowRanges(filtered.windows_USF2), summits_TFE3_CHIP_FLCNKO)
#TFE3.filtered.windows.USF2 <- filtered.windows_USF2[keep,]
TFE3.filtered.windows.USF2 <- filtered.windows_USF2

#keep <- overlapsAny(rowRanges(TFE3.filtered.windows.USF2), promoters)
#TFE3.filtered.windows.USF2.promoter <- TFE3.filtered.windows.USF2[keep,]
TFE3.filtered.windows.USF2.promoter <- TFE3.filtered.windows.USF2
colData(TFE3.filtered.windows.USF2.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.USF2.promoter)$ExternalSampleName)


#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.USF2.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.USF2.promoter)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.USF2.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.USF2.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.USF2.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.USF2.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.USF2.promoter)),
                        group = colData(TFE3.filtered.windows.USF2.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.USF2.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.USF2.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)



contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.7)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.7))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)
#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.USF2.promoter) <- cbind(rowData(TFE3.filtered.windows.USF2.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.USF2.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TFE3', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/USF2_promoter_peaks_WT_vs_TFE3tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/USF2_promoter_WT_vs_TFE3.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,PValue < 0.1 & direction == 'down' & overlap != "")$overlap

areas <- subset(peaks,PValue < 0.1 & direction == 'up' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

# USF2 coloc with HDAC1 and HDAC2
peaks_USF2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/USF2_peaks.narrowPeak")
summits_USF2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/USF2_summits.bed")

peaks_HDAC1=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC1_peaks.narrowPeak")
summits_HDAC1=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC1_summits.bed")

peaks_HDAC2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC2_peaks.narrowPeak")
summits_HDAC2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/HDAC2_summits.bed")

peaks_H3K27Ac=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K27Ac_broad_peaks.broadPeak")
summits_H3K27Ac=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup/union/H3K27Ac_broad_summits.bed")

summits_H3K27Ac$name = peaks_H3K27Ac$name

anchors <- summits_USF2




USF2_bw_WT1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1706.spikein_normalized.bw")
USF2_bw_WT2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1749.spikein_normalized.bw")
USF2_bw_WT3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1750.spikein_normalized.bw")
USF2_bw_TSC1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1747.spikein_normalized.bw")
USF2_bw_TSC2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1748.spikein_normalized.bw")
USF2_bw_TSC3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1755.spikein_normalized.bw")
USF2_bw_TFE31 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1741.spikein_normalized.bw")
USF2_bw_TFE32 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1746.spikein_normalized.bw")
USF2_bw_TFE33 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1751.spikein_normalized.bw")
USF2_bw_dKO1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1708.spikein_normalized.bw")
USF2_bw_dKO2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1744.spikein_normalized.bw")
USF2_bw_dKO3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup_120/USF2s1754.spikein_normalized.bw")

HDAC1_bw_WT1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1706-HDAC1.spikein_normalized.bw")
HDAC1_bw_WT2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1749-HDAC1.spikein_normalized.bw")
HDAC1_bw_WT3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1750-HDAC1.spikein_normalized.bw")
HDAC1_bw_TSC1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1747-HDAC1.spikein_normalized.bw")
HDAC1_bw_TSC2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1748-HDAC1.spikein_normalized.bw")
HDAC1_bw_TSC3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1755-HDAC1.spikein_normalized.bw")
HDAC1_bw_TFE31 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1741-HDAC1.spikein_normalized.bw")
HDAC1_bw_TFE32 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1746-HDAC1.spikein_normalized.bw")
HDAC1_bw_TFE33 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1751-HDAC1.spikein_normalized.bw")
HDAC1_bw_dKO1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1708-HDAC1.spikein_normalized.bw")
HDAC1_bw_dKO2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1744-HDAC1.spikein_normalized.bw")
HDAC1_bw_dKO3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/1754-HDAC1.spikein_normalized.bw")

HDAC2_bw_WT1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1706.spikein_normalized.bw")
HDAC2_bw_WT2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1749.spikein_normalized.bw")
HDAC2_bw_WT3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1750.spikein_normalized.bw")
HDAC2_bw_TSC1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1747.spikein_normalized.bw")
HDAC2_bw_TSC2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1748.spikein_normalized.bw")
HDAC2_bw_TSC3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1755.spikein_normalized.bw")
HDAC2_bw_TFE31 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC21741.spikein_normalized.bw")
HDAC2_bw_TFE32 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1746.spikein_normalized.bw")
HDAC2_bw_TFE33 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1751.spikein_normalized.bw")
HDAC2_bw_dKO1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1708.spikein_normalized.bw")
HDAC2_bw_dKO2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1744.spikein_normalized.bw")
HDAC2_bw_dKO3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/HDAC2s1754.spikein_normalized.bw")

H3K27Ac_bw_WT1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1706.spikein_normalized.bw")
H3K27Ac_bw_WT2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1749.spikein_normalized.bw")
H3K27Ac_bw_WT3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1750.spikein_normalized.bw")
H3K27Ac_bw_TSC1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1747.spikein_normalized.bw")
H3K27Ac_bw_TSC2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1748.spikein_normalized.bw")
H3K27Ac_bw_TSC3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1755.spikein_normalized.bw")
H3K27Ac_bw_TFE31 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1741.spikein_normalized.bw")
H3K27Ac_bw_TFE32 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1746.spikein_normalized.bw")
H3K27Ac_bw_TFE33 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1751.spikein_normalized.bw")
H3K27Ac_bw_dKO1 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1708.spikein_normalized.bw")
H3K27Ac_bw_dKO2 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1744.spikein_normalized.bw")
H3K27Ac_bw_dKO3 <- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/BigWig/dup/H3K27Acs1754.spikein_normalized.bw")


library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(USF2_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(USF2_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(USF2_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(USF2_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(USF2_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(USF2_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(USF2_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(USF2_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(USF2_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(USF2_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(USF2_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(USF2_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)



quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT_USF2 <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC_USF2 <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3_USF2 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO_USF2 <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3

mat_WT1<- normalizeToMatrix(HDAC1_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(HDAC1_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(HDAC1_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(HDAC1_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(HDAC1_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(HDAC1_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(HDAC1_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(HDAC1_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(HDAC1_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(HDAC1_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(HDAC1_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(HDAC1_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT_HDAC1 <- (mat_WT1 + mat_WT2 + mat_WT3) / 3 / 1.5
mat_TSC_HDAC1 <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3 / 1.5
mat_TFE3_HDAC1 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3 / 1.5
mat_dKO_HDAC1 <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3 / 1.5


mat_WT1<- normalizeToMatrix(HDAC2_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(HDAC2_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(HDAC2_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(HDAC2_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(HDAC2_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(HDAC2_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(HDAC2_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(HDAC2_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(HDAC2_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(HDAC2_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(HDAC2_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(HDAC2_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT_HDAC2 <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC_HDAC2 <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3_HDAC2 <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO_HDAC2 <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3 


library(EnrichedHeatmap)
mat_WT1<- normalizeToMatrix(H3K27Ac_bw_WT1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT2<- normalizeToMatrix(H3K27Ac_bw_WT2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_WT3<- normalizeToMatrix(H3K27Ac_bw_WT3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TSC1<- normalizeToMatrix(H3K27Ac_bw_TSC1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC2<- normalizeToMatrix(H3K27Ac_bw_TSC2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TSC3<- normalizeToMatrix(H3K27Ac_bw_TSC3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_TFE31<- normalizeToMatrix(H3K27Ac_bw_TFE31, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE32<- normalizeToMatrix(H3K27Ac_bw_TFE32, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_TFE33<- normalizeToMatrix(H3K27Ac_bw_TFE33, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

mat_dKO1<- normalizeToMatrix(H3K27Ac_bw_dKO1, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO2<- normalizeToMatrix(H3K27Ac_bw_dKO2, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)
mat_dKO3<- normalizeToMatrix(H3K27Ac_bw_dKO3, anchors, value_column = "score",
  extend= 1000, mean_mode = "w0", w=50)

quantile(mat_WT1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_WT3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TSC3, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE31, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE32, c(0.1,0.25,0.5,0.9,1))
quantile(mat_TFE33, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO1, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO2, c(0.1,0.25,0.5,0.9,1))
quantile(mat_dKO3, c(0.1,0.25,0.5,0.9,1))

mat_WT_H3K27Ac <- (mat_WT1 + mat_WT2 + mat_WT3) / 3
mat_TSC_H3K27Ac <- (mat_TSC1 + mat_TSC2 + mat_TSC3) / 3
mat_TFE3_H3K27Ac <- (mat_TFE31 + mat_TFE32 + mat_TFE33) / 3
mat_dKO_H3K27Ac <- (mat_dKO1 + mat_dKO2 + mat_dKO3) / 3



col_fun<- circlize::colorRamp2(c(0, 1), c("white", "red"))



ht_list <-  EnrichedHeatmap(mat_WT_USF2, pos_line = FALSE, column_title="WT_USF2", name = "WT_USF2", col=col_fun) + 
EnrichedHeatmap(mat_TSC_USF2, pos_line = FALSE, column_title="TSC_USF2", name = "TSC_USF2", col=col_fun) + 
EnrichedHeatmap(mat_TFE3_USF2, pos_line = FALSE, column_title="TFE3_USF2", name = "TFE3_USF2", col=col_fun) + 
EnrichedHeatmap(mat_dKO_USF2, pos_line = FALSE, column_title="dKO_USF2", name = "dKO_USF2", col=col_fun) + 
EnrichedHeatmap(mat_WT_HDAC1, pos_line = FALSE, column_title="WT_HDAC1", name = "WT_HDAC1", col=col_fun) + 
EnrichedHeatmap(mat_TSC_HDAC1, pos_line = FALSE, column_title="TSC_HDAC1", name = "TSC_HDAC1", col=col_fun) +
EnrichedHeatmap(mat_TFE3_HDAC1, pos_line = FALSE, column_title="TFE3_HDAC1", name = "TFE3_HDAC1", col=col_fun) +
EnrichedHeatmap(mat_dKO_HDAC1, pos_line = FALSE, column_title="dKO_HDAC1", name = "dKO_HDAC1", col=col_fun) +
EnrichedHeatmap(mat_WT_HDAC2, pos_line = FALSE, column_title="WT_HDAC2", name = "WT_HDAC2", col=col_fun) + 
EnrichedHeatmap(mat_TSC_HDAC2, pos_line = FALSE, column_title="TSC_HDAC2", name = "TSC_HDAC2", col=col_fun) +
EnrichedHeatmap(mat_TFE3_HDAC2, pos_line = FALSE, column_title="TFE3_HDAC2", name = "TFE3_HDAC2", col=col_fun) +
EnrichedHeatmap(mat_dKO_HDAC2, pos_line = FALSE, column_title="dKO_HDAC2", name = "dKO_HDAC2", col=col_fun) +
EnrichedHeatmap(mat_WT_H3K27Ac, pos_line = FALSE, column_title="WT_H3K27Ac", name = "WT_H3K27Ac", col=col_fun) + 
EnrichedHeatmap(mat_TSC_H3K27Ac, pos_line = FALSE, column_title="TSC_H3K27Ac", name = "TSC_H3K27Ac", col=col_fun) +
EnrichedHeatmap(mat_TFE3_H3K27Ac, pos_line = FALSE, column_title="TFE3_H3K27Ac", name = "TFE3_H3K27Ac", col=col_fun) +
EnrichedHeatmap(mat_dKO_H3K27Ac, pos_line = FALSE, column_title="dKO_H3K27Ac", name = "dKO_H3K27Ac", col=col_fun)


pdf('~/Documents/TSC_Paper/USF2_HDAC12_heatmap_ALL.pdf',height=4, width=16)
  draw(ht_list, main_heatmap =2)
dev.off()





WT1_mean<- colMeans(mat_WT1[partition == 'promoters',])
WT2_mean<- colMeans(mat_WT2[partition == 'promoters',])
WT3_mean<- colMeans(mat_WT3[partition == 'promoters',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'promoters',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'promoters',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'promoters',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'promoters',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'promoters',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'promoters',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'promoters',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'promoters',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'promoters',])

WT_mean<- colMeans(mat_WT[partition == 'promoters',])
TSC_mean<- colMeans(mat_TSC[partition == 'promoters',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'promoters',])
dKO_mean<- colMeans(mat_dKO[partition == 'promoters',])


pdf('~/Documents/TSC_Paper/USF2_TFE3_promoter_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, TFE33_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","TFE33","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/USF2_TFE3_promoter_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



WT1_mean<- colMeans(mat_WT1[partition == 'enhancers',])
WT2_mean<- colMeans(mat_WT2[partition == 'enhancers',])
WT3_mean<- colMeans(mat_WT3[partition == 'enhancers',])
TSC1_mean<- colMeans(mat_TSC1[partition == 'enhancers',])
TSC2_mean<- colMeans(mat_TSC2[partition == 'enhancers',])
TSC3_mean<- colMeans(mat_TSC3[partition == 'enhancers',])
TFE31_mean<- colMeans(mat_TFE31[partition == 'enhancers',])
TFE32_mean<- colMeans(mat_TFE32[partition == 'enhancers',])
TFE33_mean<- colMeans(mat_TFE33[partition == 'enhancers',])
dKO1_mean<- colMeans(mat_dKO1[partition == 'enhancers',])
dKO2_mean<- colMeans(mat_dKO2[partition == 'enhancers',])
dKO3_mean<- colMeans(mat_dKO3[partition == 'enhancers',])



WT_mean<- colMeans(mat_WT[partition == 'enhancers',])
TSC_mean<- colMeans(mat_TSC[partition == 'enhancers',])
TFE3_mean<- colMeans(mat_TFE3[partition == 'enhancers',])
dKO_mean<- colMeans(mat_dKO[partition == 'enhancers',])


pdf('~/Documents/TSC_Paper/USF2_TFE3_enhancer_peaks_ALL.pdf',height=3, width=4)
  bind_rows(WT1_mean, WT2_mean, WT3_mean, TSC1_mean, TSC2_mean, TSC3_mean, TFE31_mean, TFE32_mean, dKO1_mean, dKO2_mean, dKO3_mean) %>%
  mutate(factor = c("WT1", "WT2", "WT3","TSC1","TSC2","TSC3","TFE31","TFE32","dKO1","dKO2","dKO3")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#DA9195", "#DA9195","#66BB6A","#66BB6A","#66BB6A","#9E9E9E","#9E9E9E","#9E9E9E","#039BE5","#039BE5","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()



pdf('~/Documents/TSC_Paper/USF2_TFE3_enhancer_peaks_ALL_mean.pdf',height=3, width=4)
  bind_rows(WT_mean, TSC_mean, TFE3_mean, dKO_mean) %>%
  mutate(factor = c("WT","TSC","TFE3","dKO")) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor)  %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#66BB6A","#9E9E9E","#039BE5")) +
  theme_classic(base_size = 14) +
  ylab("CPM") +
  xlab("")
dev.off()


########################################
############## RNA-SEQ #################
########################################



########### H3KA27c MARKS ##############

# Looks at correlation between RNASeq WT and TSC and H3K27Ac marks

TSC.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TSC.vs.WT.csv',sep=',')
filtered.windows_H3K27Ac <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_H3K27Ac.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_H3K27Ac

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)

#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TSC.vs.WT$gene)

peaks$RNA_log2fc <- -TSC.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TSC.vs.WT$padj[idx]
df <- as.data.frame(peaks)

#TSC.vs.WT <- subset(TSC.vs.WT,padj<0.05)
#TSC.vs.WT.down <- subset(TSC.vs.WT,log2FoldChange<0)
#TSC.vs.WT.up<- subset(TSC.vs.WT,log2FoldChange>0)
#rna_up <- intersect(dKO.vs.TFE3.down$gene,genes_all)
df.down <- as.data.frame(subset(peaks, direction == 'down' & FDR < 0.1 & overlap != ""))
dbs <- c("ChEA_2022","Reactome_Pathways_2024","GO_Biological_Process_2023")
library(enrichR)
enriched.down <- enrichr(df.down$gene, dbs)
human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)


MITF_genes <- str_split(enriched.down[[1]]$Genes[5],';')[[1]]
MITF_genes <- human2mouse$V2[match(MITF_genes,human2mouse$V1)]

LXR_genes <- str_split(enriched.down[[1]]$Genes[1],';')[[1]]
LXR_genes <- human2mouse$V2[match(LXR_genes,human2mouse$V1)]
RXR_genes <- str_split(enriched.down[[1]]$Genes[2],';')[[1]]
RXR_genes <- human2mouse$V2[match(RXR_genes,human2mouse$V1)]
PPARA_genes <- str_split(enriched.down[[1]]$Genes[3],';')[[1]]
PPARA_genes <- human2mouse$V2[match(PPARA_genes,human2mouse$V1)]
LRP_genes <- c(LXR_genes,RXR_genes,PPARA_genes)

Clear_genes <- read.csv('~/Documents/TSC_Paper/Clear_genes.csv',header=FALSE)$V1

df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_H3K27Ac_correlation.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_H3K27Ac_correlation_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()




# Looks at correlation between RNASeq WT and TSC and H3K27Ac marks

TFE3.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TFE3.vs.WT.csv',sep=',')
filtered.windows_H3K27Ac <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_H3K27Ac.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_H3K27Ac

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TFE3.vs.WT$gene)

peaks$RNA_log2fc <- -TFE3.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TFE3.vs.WT$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_H3K27Ac_correlation_WT_TFE3.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_H3K27Ac_correlation_WT_TFE3_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()




######### ATAC ############

# Looks at correlation between RNASeq WT and TSC and ATAC

TSC.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TSC.vs.WT.csv',sep=',')
filtered.windows_ATAC <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_ATAC.rds')
logCPM <- csaw::calculateCPM(filtered.windows_ATAC, log=TRUE, use.offsets= FALSE)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(filtered.windows_ATAC))) |>
  as_tibble()


### Calculate log2FC for ATAC
filtered.windows_ATAC.normed <- filtered.windows_ATAC[df$TSC > - 3 & df$WT > - 3]
filtered.windows_ATAC.normed <- normOffsets(filtered.windows_ATAC.normed)

colData(filtered.windows_ATAC.normed)$SampleType <- sub("^[^_]*_", "", colData(filtered.windows_ATAC.normed)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(filtered.windows_ATAC.normed, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(filtered.windows_ATAC.normed)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(filtered.windows_ATAC.normed))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(filtered.windows_ATAC.normed)$overlapTSS,#
         GC = rowData(filtered.windows_ATAC.normed)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(filtered.windows_ATAC.normed,
                        genes = as.data.frame(rowData(filtered.windows_ATAC.normed)),
                        group = colData(filtered.windows_ATAC.normed)$SampleType)
colnames(dgel) <- colnames(filtered.windows_ATAC.normed)
rownames(dgel) <- as.character(rowRanges(filtered.windows_ATAC.normed))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)

#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(filtered.windows_ATAC.normed) <- cbind(rowData(filtered.windows_ATAC.normed), tt)
merged <- mergeWindows(rowRanges(filtered.windows_ATAC.normed), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TSC.vs.WT$gene)

peaks$RNA_log2fc <- -TSC.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TSC.vs.WT$padj[idx]
df <- as.data.frame(peaks)

diff_open <- subset(df,direction == 'down' & FDR < 0.1 & RNA_p_adj < 0.1 & RNA_log2fc < -0.1)$gene


peakAnno.edb <- annotatePeak(peaks_TFE3_CHIP_WT, tssRegion=c(-3000, 3000),
                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

anno <- as.data.frame(peakAnno.edb)$annotation


peaks_int <- subset(peaks,direction == 'down' & FDR < 0.1 & RNA_p_adj < 0.1 & RNA_log2fc < -0.1)

#subsetByOverlaps(peaks_int,peaks_TFE3_CHIP_FLCNKO)

peakAnno.edb <- annotatePeak(peaks_int, tssRegion=c(-3000, 3000),
                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

anno <- as.data.frame(peakAnno.edb)$annotation


p1 <- plotAnnoBar(peakAnno.edb)
p2 <- plotDistToTSS(peakAnno.edb)

pdf('~/Documents/TSC_Paper/ATAC_RNA_Seq_int_TSC_vs_WT_dist.pdf',width=4.75,height=2.5)
p1 / p2
dev.off()

#1576 RNA-Seq DEGs
#379 ATAC DOPs
#Overlap 67

library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)
enriched.down <- enrichr(diff_open, dbs)

enriched.down[[1]] <- enriched[[1]][enriched.down[[1]]$Adjusted.P.value < 0.05,]
p1<- ggplot(enriched.down[[1]][order(enriched.down[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),], (aes(x=-log(Adjusted.P.value), y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + geom_point() + xlab('Adjusted P Value') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched.down[[1]][order(enriched.down[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"\\s+"), `[`, 1))) + theme(axis.text=element_text(colour="black"))

pdf('~/Documents/TSC_Paper/ATAC_RNA_Seq_int_TSC_vs_WT_UP.pdf',width=4.75,height=2.5)
p1
dev.off()

df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_ATAC_correlation.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_ATAC_correlation_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()




# Looks at correlation between RNASeq WT and TFE3 ATAC marks

TFE3.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TFE3.vs.WT.csv',sep=',')
filtered.windows_ATAC <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_ATAC.rds')
logCPM <- csaw::calculateCPM(filtered.windows_ATAC, log=TRUE, use.offsets= FALSE)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(filtered.windows_ATAC))) |>
  as_tibble()


### Calculate log2FC for ATAC
filtered.windows_ATAC.normed <- filtered.windows_ATAC[df$TSC > - 3 & df$WT > - 3]
filtered.windows_ATAC.normed <- normOffsets(filtered.windows_ATAC.normed)

colData(filtered.windows_ATAC.normed)$SampleType <- sub("^[^_]*_", "", colData(filtered.windows_ATAC.normed)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(filtered.windows_ATAC.normed, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(filtered.windows_ATAC.normed)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(filtered.windows_ATAC.normed))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(filtered.windows_ATAC.normed)$overlapTSS,#
         GC = rowData(filtered.windows_ATAC.normed)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(filtered.windows_ATAC.normed,
                        genes = as.data.frame(rowData(filtered.windows_ATAC.normed)),
                        group = colData(filtered.windows_ATAC.normed)$SampleType)
colnames(dgel) <- colnames(filtered.windows_ATAC.normed)
rownames(dgel) <- as.character(rowRanges(filtered.windows_ATAC.normed))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TFE3, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(filtered.windows_ATAC.normed) <- cbind(rowData(filtered.windows_ATAC.normed), tt)
merged <- mergeWindows(rowRanges(filtered.windows_ATAC.normed), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'up' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TFE3.vs.WT$gene)

peaks$RNA_log2fc <- -TFE3.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TFE3.vs.WT$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_ATAC_correlation_WT_TFE3.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_ATAC_correlation_WT_TFE3_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()




# Looks at correlation between RNASeq TSC and dKO ATAC marks

dKO.vs.TSC <- read.csv('~/Documents/TSC_RNASeq/dKO.vs.TSC.csv',sep=',')
filtered.windows_ATAC <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_ATAC.rds')
logCPM <- csaw::calculateCPM(filtered.windows_ATAC, log=TRUE, use.offsets= FALSE)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(filtered.windows_ATAC))) |>
  as_tibble()


### Calculate log2FC for ATAC
filtered.windows_ATAC.normed <- filtered.windows_ATAC[df$TSC > - 3 & df$WT > - 3]
filtered.windows_ATAC.normed <- normOffsets(filtered.windows_ATAC.normed)

colData(filtered.windows_ATAC.normed)$SampleType <- sub("^[^_]*_", "", colData(filtered.windows_ATAC.normed)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(filtered.windows_ATAC.normed, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(filtered.windows_ATAC.normed)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(filtered.windows_ATAC.normed))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(filtered.windows_ATAC.normed)$overlapTSS,#
         GC = rowData(filtered.windows_ATAC.normed)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(filtered.windows_ATAC.normed,
                        genes = as.data.frame(rowData(filtered.windows_ATAC.normed)),
                        group = colData(filtered.windows_ATAC.normed)$SampleType)
colnames(dgel) <- colnames(filtered.windows_ATAC.normed)
rownames(dgel) <- as.character(rowRanges(filtered.windows_ATAC.normed))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(TSC - dKO, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(filtered.windows_ATAC.normed) <- cbind(rowData(filtered.windows_ATAC.normed), tt)
merged <- mergeWindows(rowRanges(filtered.windows_ATAC.normed), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.01 & direction == 'up' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,dKO.vs.TSC$gene)

peaks$RNA_log2fc <- -dKO.vs.TSC$log2FoldChange[idx]
peaks$RNA_p_adj <- dKO.vs.TSC$padj[idx]
df <- as.data.frame(peaks)


df_signif <- subset(df, RNA_p_adj<0.05 & FDR < 0.2)
cat(subset(df_signif, direction == 'up' & RNA_log2fc > 0)$gene,sep='\n')
cat(unique(subset(df, direction == 'up' & FDR<0.05 & rep.logFC > 1 & gene!='NA')$gene),sep='\n')

write.table(df,'~/Documents/TSC_ATAC/ATAC_RNA_integrate_TSC_vs_dKO.csv')

df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_ATAC_correlation_TSC_dKO.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=rep.logFC,y=RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_ATAC_correlation_TSC_dKO_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()



########### H3K4me3 MARKS ##############

# Looks at correlation between RNASeq WT and TSC and H3K4me3 marks

TSC.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TSC.vs.WT.csv',sep=',')
filtered.windows_H3K4me3 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_H3K4me3.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_H3K4me3

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)

#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TSC.vs.WT$gene)

peaks$RNA_log2fc <- -TSC.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TSC.vs.WT$padj[idx]
df <- as.data.frame(peaks)

df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_H3K4me3_correlation.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_H3K4me3_correlation_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()




# Looks at correlation between RNASeq WT and TFE3 and H3K4me3 marks

TFE3.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TFE3.vs.WT.csv',sep=',')
filtered.windows_H3K4me3 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_H3K4me3.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_H3K4me3

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TFE3, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TFE3.vs.WT$gene)

peaks$RNA_log2fc <- -TFE3.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TFE3.vs.WT$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_H3K4me3_correlation_WT_TFE3.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_H3K4me3_correlation_WT_TFE3_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()



# Looks at correlation between RNASeq TSC and dKO and H3K4me3 marks

dKO.vs.TSC <- read.csv('~/Documents/TSC_RNASeq/dKO.vs.TSC.csv',sep=',')
filtered.windows_H3K4me3 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_H3K4me3.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_H3K4me3

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(TSC - dKO, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,dKO.vs.TSC$gene)

peaks$RNA_log2fc <- -dKO.vs.TSC$log2FoldChange[idx]
peaks$RNA_p_adj <- dKO.vs.TSC$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_H3K4me3_correlation_TSC_dKO.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_H3K4me3_correlation_TSC_dKO_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()


########### USF2 MARKS ##############

# Looks at correlation between RNASeq WT and TSC and USF2 marks

TSC.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TSC.vs.WT.csv',sep=',')
filtered.windows_USF2 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_USF2.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_USF2

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)

#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TSC.vs.WT$gene)

peaks$RNA_log2fc <- -TSC.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TSC.vs.WT$padj[idx]
df <- as.data.frame(peaks)

df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_USF2_correlation.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_USF2_correlation_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()




# Looks at correlation between RNASeq WT and TFE3 and USF2 marks

TFE3.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TFE3.vs.WT.csv',sep=',')
filtered.windows_USF2 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_USF2.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_USF2

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TFE3, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TFE3.vs.WT$gene)

peaks$RNA_log2fc <- -TFE3.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TFE3.vs.WT$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_USF2_correlation_WT_TFE3.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_USF2_correlation_WT_TFE3_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()



# Looks at correlation between RNASeq TSC and dKO and USF2 marks

dKO.vs.TSC <- read.csv('~/Documents/TSC_RNASeq/dKO.vs.TSC.csv',sep=',')
filtered.windows_USF2 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_USF2.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_H3K4me3

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(TSC - dKO, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,dKO.vs.TSC$gene)

peaks$RNA_log2fc <- -dKO.vs.TSC$log2FoldChange[idx]
peaks$RNA_p_adj <- dKO.vs.TSC$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_USF2_correlation_TSC_dKO.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_USF2_correlation_TSC_dKO_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()




########### HDAC1 MARKS ##############

# Looks at correlation between RNASeq WT and TSC and HDAC1 marks

TSC.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TSC.vs.WT.csv',sep=',')
filtered.windows_HDAC1 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_HDAC1.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_HDAC1

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)

#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TSC.vs.WT$gene)

peaks$RNA_log2fc <- -TSC.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TSC.vs.WT$padj[idx]
df <- as.data.frame(peaks)

df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_HDAC1_correlation.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_HDAC1_correlation_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()




# Looks at correlation between RNASeq WT and TFE3 and USF2 marks

TFE3.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TFE3.vs.WT.csv',sep=',')
filtered.windows_HDAC1 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_HDAC1.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_HDAC1

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TFE3, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TFE3.vs.WT$gene)

peaks$RNA_log2fc <- -TFE3.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TFE3.vs.WT$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_HDAC1_correlation_WT_TFE3.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_HDAC1_correlation_WT_TFE3_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()



# Looks at correlation between RNASeq TSC and dKO and USF2 marks

dKO.vs.TSC <- read.csv('~/Documents/TSC_RNASeq/dKO.vs.TSC.csv',sep=',')
filtered.windows_HDAC1 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_HDAC1.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_HDAC1

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(TSC - dKO, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,dKO.vs.TSC$gene)

peaks$RNA_log2fc <- -dKO.vs.TSC$log2FoldChange[idx]
peaks$RNA_p_adj <- dKO.vs.TSC$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_HDAC1_correlation_TSC_dKO.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=rep.logFC,y=RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_HDAC1_correlation_TSC_dKO_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()




########### HDAC2 MARKS ##############

# Looks at correlation between RNASeq WT and TSC and HDAC1 marks

TSC.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TSC.vs.WT.csv',sep=',')
filtered.windows_HDAC2 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_HDAC2.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_HDAC2

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)

#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TSC.vs.WT$gene)

peaks$RNA_log2fc <- -TSC.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TSC.vs.WT$padj[idx]
df <- as.data.frame(peaks)

df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_HDAC2_correlation.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_HDAC2_correlation_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()

# Looks at correlation between RNASeq WT and TFE3 and HDAC2 marks

TFE3.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TFE3.vs.WT.csv',sep=',')
filtered.windows_HDAC2 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_HDAC2.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_HDAC2

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TFE3, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TFE3.vs.WT$gene)

peaks$RNA_log2fc <- -TFE3.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TFE3.vs.WT$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_HDAC2_correlation_WT_TFE3.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_HDAC2_correlation_WT_TFE3_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()



# Looks at correlation between RNASeq TSC and dKO and USF2 marks

dKO.vs.TSC <- read.csv('~/Documents/TSC_RNASeq/dKO.vs.TSC.csv',sep=',')
filtered.windows_HDAC2 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_HDAC2.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_HDAC2

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(TSC - dKO, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,dKO.vs.TSC$gene)

peaks$RNA_log2fc <- -dKO.vs.TSC$log2FoldChange[idx]
peaks$RNA_p_adj <- dKO.vs.TSC$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_HDAC2_correlation_TSC_dKO.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=rep.logFC,y=RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_HDAC2_correlation_TSC_dKO_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()



########### TFE3 MARKS ##############

# Looks at correlation between RNASeq WT and TSC and HDAC1 marks

TSC.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TSC.vs.WT.csv',sep=',')
filtered.windows_TFE3 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_TFE3.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_TFE3

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TSC, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


### VIZ MA

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.1))

ggplot(df,mapping = aes(x = .data[['TSC']]  + .data[['WT']],y = .data[['WT']] - .data[['TSC']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'WT', 'TSC')) +
    theme_bw(base_size = 15)

#### END VIZ MA


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))


###


idx <- match(peaks$gene,TSC.vs.WT$gene)

peaks$RNA_log2fc <- -TSC.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TSC.vs.WT$padj[idx]
df <- as.data.frame(peaks)

df.down <- as.data.frame(subset(peaks, direction == 'down' & FDR < 0.1 & overlap != ""))
dbs <- c("ChEA_2022","Reactome_Pathways_2024","GO_Biological_Process_2023")
library(enrichR)
enriched.down <- enrichr(df.down$gene, dbs)
human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)


enriched.down[[1]] <- enriched[[1]][enriched.down[[1]]$Adjusted.P.value < 0.05,]
p1<- ggplot(enriched.down[[1]][order(enriched.down[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),], (aes(x=-log(Adjusted.P.value), y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + geom_point() + xlab('Adjusted P Value') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched.down[[1]][order(enriched.down[[1]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"\\s+"), `[`, 1))) + theme(axis.text=element_text(colour="black"))

pdf('~/Documents/TSC_Paper/TFE3_TSC_vs_WT_UP.pdf',width=4.75,height=4)
p1
dev.off()


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_TFE3_correlation.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_TFE3_correlation_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()


# Looks at correlation between RNASeq WT and TFE3 and TFE3 marks

TFE3.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TFE3.vs.WT.csv',sep=',')
filtered.windows_TFE3 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_TFE3.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_TFE3

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(WT - TFE3, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,TFE3.vs.WT$gene)

peaks$RNA_log2fc <- -TFE3.vs.WT$log2FoldChange[idx]
peaks$RNA_p_adj <- TFE3.vs.WT$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_TFE3_correlation_WT_TFE3.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=-rep.logFC,y=-RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(-rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(-RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_TFE3_correlation_WT_TFE3_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()



# Looks at correlation between RNASeq TSC and dKO and USF2 marks

dKO.vs.TSC <- read.csv('~/Documents/TSC_RNASeq/dKO.vs.TSC.csv',sep=',')
filtered.windows_TFE3 <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_TFE3.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_TFE3

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(TSC - dKO, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,dKO.vs.TSC$gene)

peaks$RNA_log2fc <- -dKO.vs.TSC$log2FoldChange[idx]
peaks$RNA_p_adj <- dKO.vs.TSC$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_TFE3_correlation_TSC_dKO.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=rep.logFC,y=RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_TFE3_correlation_TSC_dKO_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()




################################################
################ BIGWIG VIZ ATAC ###############
################################################



windows_ATAC_narrow <- windowCounts(bamFiles_ATAC, ext=0, width=10,
                     bin=TRUE, filter=0, param=param,
                     BPPARAM=MulticoreParam(12))

windows_ATAC_mod <- windowCounts(bamFiles_ATAC, ext=0, width=25,
                     bin=TRUE, filter=0, param=param,
                     BPPARAM=MulticoreParam(12))

# Add sample annotation
colData(windows_ATAC_narrow) <- cbind(colData(windows_ATAC_narrow), alignments_ATAC)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, rowRanges(windows_ATAC_narrow))
rowData(windows_ATAC_narrow) <- letterFrequency(seqs, letters=c('GC', 'N'), as.prob=TRUE) # G or C, N frequency
is.agap <- rowData(windows_ATAC_narrow)[,'N'] > 2/3 
table(is.agap)/nrow(windows_ATAC_narrow)
is.smaller <- width(rowRanges(windows_ATAC_narrow)) != 10
table(is.smaller)
windows_ATAC_narrow <- windows_ATAC_narrow[!is.agap & !is.smaller,]


colData(windows_ATAC_mod) <- cbind(colData(windows_ATAC_mod), alignments_ATAC)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, rowRanges(windows_ATAC_mod))
rowData(windows_ATAC_mod) <- letterFrequency(seqs, letters=c('GC', 'N'), as.prob=TRUE) # G or C, N frequency
is.agap <- rowData(windows_ATAC_mod)[,'N'] > 2/3 
table(is.agap)/nrow(windows_ATAC_mod)
is.smaller <- width(rowRanges(windows_ATAC_mod)) != 25
table(is.smaller)
windows_ATAC_mod <- windows_ATAC_mod[!is.agap & !is.smaller,]

#TSS annotation
tss <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
tss <- resize(tss, width=1, fix="start")
overlapTSS <- rowRanges(windows_ATAC_narrow) %over% tss
rowData(windows_ATAC_narrow)$overlapTSS <- overlapTSS
overlapTSS <- rowRanges(windows_ATAC_mod) %over% tss
rowData(windows_ATAC_mod)$overlapTSS <- overlapTSS

saveRDS(windows_ATAC_narrow,file='~/Documents/TSC_ATAC/windows_ATAC_narrow.rds')
saveRDS(windows_ATAC_mod,file='~/Documents/TSC_ATAC/windows_ATAC_mod.rds')

#windows_ATAC_narrow <- readRDS('~/Documents/TSC_ATAC/windows_ATAC_narrow.rds')


peaks_ATAC=import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/union/union_summits.bed")

# count reads in windows specified by MACS2                                      
peak.counts <- regionCounts(bamFiles_ATAC, peaks_ATAC, param=param)
peak.abundances <- aveLogCPM(asDGEList(peak.counts))
peak.counts.filt <- peak.counts[peak.abundances > -3, ]
#suppressWarnings(keep <- overlapsAny(rowRanges(windows_ATAC_narrow), peak.counts.filt))
#filtered.windows_ATAC <- windows_ATAC_narrow[keep,]
suppressWarnings(keep <- overlapsAny(rowRanges(windows_ATAC_mod), peak.counts.filt))
filtered.windows_ATAC <- windows_ATAC_mod[keep,]


colData(filtered.windows_ATAC)$SampleType <- sub("^[^_]*_", "", colData(filtered.windows_ATAC)$ExternalSampleName)



logCPM <- csaw::calculateCPM(filtered.windows_ATAC, log=TRUE, use.offsets= FALSE)
ct <- split(colnames(logCPM), colData(filtered.windows_ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(filtered.windows_ATAC))) |>
  as_tibble()

filtered.windows_ATAC.normed <- filtered.windows_ATAC[df$TSC > - 3 & df$WT > - 3]
windows_ATAC_narrow <- normOffsets(filtered.windows_ATAC.normed,se.out = windows_ATAC_narrow)

windows_ATAC_mod <- normOffsets(filtered.windows_ATAC.normed,se.out = windows_ATAC_mod)


colData(windows_ATAC_narrow)$SampleType <- sub("^[^_]*_", "", colData(windows_ATAC_narrow)$ExternalSampleName)
colData(windows_ATAC_mod)$SampleType <- sub("^[^_]*_", "", colData(windows_ATAC_mod)$ExternalSampleName)



#CPM <-  csaw::calculateCPM(windows_ATAC_narrow, log=FALSE, use.offsets = TRUE)
#gr <- rowRanges(windows_ATAC_narrow)
CPM <-  csaw::calculateCPM(windows_ATAC_mod, log=FALSE, use.offsets = TRUE)
gr <- rowRanges(windows_ATAC_mod)

for (cn in colnames(CPM)) {
  gr$score <- CPM[,cn]
  rtracklayer::export.bw(gr, sprintf('~/Documents/TSC_ATAC/%s.bigWig', cn),
                         format = "bigWig")
}



gr <- rowRanges(windows_ATAC_mod)
gr$score <- rowMeans(CPM[,c('1706_WT','1749_WT','1750_WT')])
rtracklayer::export.bw(gr, '~/Documents/TSC_ATAC/WT_mod.bigWig',
                         format = "bigWig")


gr <- rowRanges(windows_ATAC_mod)
gr$score <- rowMeans(CPM[,c('1748_TSC','1755_TSC')])
rtracklayer::export.bw(gr, '~/Documents/TSC_ATAC/TSC_mod.bigWig',
                         format = "bigWig")



gr <- rowRanges(windows_ATAC_mod)
gr$score <- rowMeans(CPM[,c('1741_TFE3','1746_TFE3','1751_TFE3')])
rtracklayer::export.bw(gr, '~/Documents/TSC_ATAC/TFE3_mod.bigWig',
                         format = "bigWig")


gr <- rowRanges(windows_ATAC_mod)
gr$score <- rowMeans(CPM[,c('1708_dKO','1744_dKO','1754_dKO')])
rtracklayer::export.bw(gr, '~/Documents/TSC_ATAC/dKO_mod.bigWig',
                         format = "bigWig")

source("https://github.com/PoisonAlien/trackplot/blob/master/R/trackplot.R?raw=true")

#Path to bigWig files
bigWigs = c('~/Documents/TSC_ATAC/WT_mod.bigWig', '~/Documents/TSC_ATAC/TSC_mod.bigWig', 
  '~/Documents/TSC_ATAC/TFE3_mod.bigWig', '~/Documents/TSC_ATAC/dKO_mod.bigWig',
  '~/Documents/TSC_ATAC/H3K27Ac_WT.bw',
  '~/Documents/TSC_ATAC/H3K27Ac_TSC1.bw',
'~/Documents/TSC_ATAC/H3K27Ac_TFE3.bw',
'~/Documents/TSC_ATAC/H3K27Ac_dKO.bw',
'~/Documents/TSC_ATAC/H3K4me3_WT.bw',
'~/Documents/TSC_ATAC/H3K4me3_TSC1.bw',
'~/Documents/TSC_ATAC/H3K4me3_TFE3.bw',
'~/Documents/TSC_ATAC/H3K4me3_dKO.bw')

#Make a table of bigWigs along with ref genome build
bigWigs = read_coldata(bws = bigWigs, build = "mm10")

rragd <- '4:32980998-33024180'
lamp1 <- '8:13157161-13177338'
ctsd <- '7:142373910-142389827'
hexa <- '9:59537540-59567109'
hexb <- '13:97174331-97200357'
gusb <- '5:129987011-130005049'
apoe <- '7:19694269-19703310'
apoc1 <- '7:19687480-19694659'
pck1 <- '2:173151073-173161274'
plin5 <- '17:56,109,597-56,119,606'
abca1 <- '4:53028796-53161988'
pcca <- '14:122482,365-122891944'
srebf1 <- '11:60197094-60226186'
rragc <- '4:123915433-123938997'
neu1 <- '17:34929253-34939297'
uap1l1 <- '2:25359489-25367682'

rragdt = track_extract(colData = bigWigs, loci = rragd)
lamp1t = track_extract(colData = bigWigs, loci = lamp1)
ctsdt = track_extract(colData = bigWigs, loci = ctsd)
hexat = track_extract(colData = bigWigs, loci = hexa)
hexabt = track_extract(colData = bigWigs, loci = hexb)
gusbt = track_extract(colData = bigWigs, loci = gusb)
apoet = track_extract(colData = bigWigs, loci = apoe)
apoc1 = track_extract(colData = bigWigs, loci = apoc1)
pck1t= track_extract(colData = bigWigs, loci = pck1)
plin5t= track_extract(colData = bigWigs, loci = plin5)
abca1t= track_extract(colData = bigWigs, loci = abca1)
pccat= track_extract(colData = bigWigs, loci = pcca)
srebf1t= track_extract(colData = bigWigs, loci = srebf1)
rragct = track_extract(colData = bigWigs, loci = rragc)
neu1t = track_extract(colData = bigWigs, loci = neu1)
uap1l1t = track_extract(colData = bigWigs, loci = uap1l1)


track_cols = c("#d35400","#2980b9", "#27ae60","#FA8072")

pdf('~/Documents/TSC_Paper/track_rragd.pdf',height=7, width=3)
track_plot(summary_list = rragdt,col = track_cols,y_max=c(5.2,5.2,5.2,5.2,0.12,0.12,0.12,0.12,0.6,0.6,0.6,0.6))
dev.off()

pdf('~/Documents/TSC_Paper/track_lamp1.pdf',height=7, width=3)
track_plot(summary_list = lamp1t,col = track_cols,y_max=c(2.4,2.4,2.4,2.4,0.45,0.45,0.45,0.45,2,2,2,2))
dev.off()

pdf('~/Documents/TSC_Paper/track_ctsd.pdf',height=7, width=3)
track_plot(summary_list = ctsdt,col = track_cols,y_max=c(2,2,2,2,0.9,0.9,0.9,0.9,3,3,3,3))
dev.off()

pdf('~/Documents/TSC_Paper/track_hexa.pdf',height=7, width=3)
track_plot(summary_list = hexat,col = track_cols,y_max=c(3.7,3.7,3.7,3.7,0.65,0.65,0.65,0.65,2,2,2,2))
dev.off()

pdf('~/Documents/TSC_Paper/track_hexb.pdf',height=7, width=3)
track_plot(summary_list = hexabt,col = track_cols)
dev.off()

pdf('~/Documents/TSC_Paper/track_gusb.pdf',height=7, width=3)
track_plot(summary_list = gusbt,col = track_cols)
dev.off()

pdf('~/Documents/TSC_Paper/track_apoe.pdf',height=7, width=3)
track_plot(summary_list = apoet,col = track_cols)
dev.off()

pdf('~/Documents/TSC_Paper/track_apoc1.pdf',height=7, width=3)
track_plot(summary_list = apoc1t,col = track_cols)
dev.off()

pdf('~/Documents/TSC_Paper/track_pck1.pdf',height=7, width=3)
track_plot(summary_list = pck1t,col = track_cols)
dev.off()

pdf('~/Documents/TSC_Paper/track_plin5.pdf',height=7, width=3)
track_plot(summary_list = plin5t,col = track_cols)
dev.off()

pdf('~/Documents/TSC_Paper/track_pcca.pdf',height=7, width=3)
track_plot(summary_list = pccat,col = track_cols)
dev.off()

pdf('~/Documents/TSC_Paper/track_abca1.pdf',height=7, width=3)
track_plot(summary_list = abca1t,col = track_cols)
dev.off()

pdf('~/Documents/TSC_Paper/track_srebf1.pdf',height=7, width=3)
track_plot(summary_list = srebf1t,col = track_cols)
dev.off()


pdf('~/Documents/TSC_Paper/track_rragc.pdf',height=7, width=3)
track_plot(summary_list = rragct,col = track_cols)
dev.off()


pdf('~/Documents/TSC_Paper/track_neu1.pdf',height=7, width=3)
track_plot(summary_list = neu1t,col = track_cols,y_max=c(1.7,1.7,1.7,1.7,0.9,0.9,0.9,0.9,3.5,3.5,3.5,3.5))
dev.off()

pdf('~/Documents/TSC_Paper/track_uap1l1.pdf',height=7, width=3)
track_plot(summary_list = uap1l1t,col = track_cols,y_max=c(1,1,1,1,0.4,0.4,0.4,0.4,1.5,1.5,1.5,1.5))
dev.off()


#bigwigAverage -b H3K27Acs1706.spikein_normalized.bw H3K27Acs1749.spikein_normalized.bw H3K27Acs1750.spikein_normalized.bw -o H3K27Ac_WT.bw -p max/2
#bigwigAverage -b H3K27Acs1747.spikein_normalized.bw H3K27Acs1748.spikein_normalized.bw H3K27Acs1755.spikein_normalized.bw -o H3K27Ac_TSC1.bw -p max/2
#bigwigAverage -b H3K27Acs1741.spikein_normalized.bw H3K27Acs1746.spikein_normalized.bw H3K27Acs1751.spikein_normalized.bw -o H3K27Ac_TFE3.bw -p max/2
#bigwigAverage -b H3K27Acs1708.spikein_normalized.bw H3K27Acs1744.spikein_normalized.bw H3K27Acs1754.spikein_normalized.bw -o H3K27Ac_dKO.bw -p max/2


#bigwigAverage -b 1706-H3K4me3.spikein_normalized.bw 1750-H3K4me3.spikein_normalized.bw -o H3K4me3_WT.bw -p max/2
#bigwigAverage -b 1747-H3K4me3.spikein_normalized.bw 1748-H3K4me3.spikein_normalized.bw 1755-H3K4me3.spikein_normalized.bw -o H3K4me3_TSC1.bw -p max/2
#bigwigAverage -b 1741-H3K4me3.spikein_normalized.bw 1746-H3K4me3.spikein_normalized.bw 1751-H3K4me3.spikein_normalized.bw -o H3K4me3_TFE3.bw -p max/2
#bigwigAverage -b 1708-H3K4me3.spikein_normalized.bw 1744-H3K4me3.spikein_normalized.bw 1751-H3K4me3.spikein_normalized.bw -o H3K4me3_dKO.bw -p max/2


################################################
################ ANNO TFE3 PEAKS ###############
################################################


library(ChIPseeker)

summits_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_summits.bed",format="bed")
summits_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_summits.bed",format="bed")
peaks_TFE3_CHIP_WT=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/WT/WT_mm10_peaks.bed",format="bed")
peaks_TFE3_CHIP_FLCNKO=import("/Volumes/Extreme SSD/CUTRUN/Normal_chow_TFE3_ChIP/Bam files/FLCNKO/FLCNKO_mm10_peaks.bed",format="bed")

seqlevelsStyle(summits_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(summits_TFE3_CHIP_FLCNKO) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_WT) <- "NCBI"
seqlevelsStyle(peaks_TFE3_CHIP_FLCNKO) <- "NCBI"

peakAnno.edb <- annotatePeak(peaks_TFE3_CHIP_WT, tssRegion=c(-3000, 3000),
                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

pdf('~/Documents/TSC_Paper/tfe_wt_peak_loc.pdf',height=3, width=6)
plotAnnoPie(peakAnno.edb)
dev.off()



peakAnno.edb <- annotatePeak(summits_TFE3_CHIP_FLCNKO, tssRegion=c(-3000, 3000),
                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb="org.Mm.eg.db")

pdf('~/Documents/TSC_Paper/tfe_flcnko_peak_loc.pdf',height=3, width=6)
plotAnnoPie(peakAnno.edb)
dev.off()



# Looks at correlation between RNASeq TSC and dKO and H3K27Ac marks

dKO.vs.TSC <- read.csv('~/Documents/TSC_RNASeq/dKO.vs.TSC.csv',sep=',')
filtered.windows_H3K27Ac <- readRDS(file='~/Documents/TSC_ATAC/filtered.windows_H3K27Ac.rds')

### Calculate log2FC for H3K27Ac
TFE3.filtered.windows.H3K27Ac.promoter<-filtered.windows_H3K27Ac

colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType <- sub("^[^_]*_", "", colData(TFE3.filtered.windows.H3K27Ac.promoter)$ExternalSampleName)

#### VIZ MA
logCPM <- csaw::calculateCPM(TFE3.filtered.windows.H3K27Ac.promoter, log=TRUE, use.offsets= TRUE)
ct <- split(colnames(logCPM), colData(TFE3.filtered.windows.ATAC)$SampleType)
df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(TFE3.filtered.windows.H3K27Ac.promoter))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(TFE3.filtered.windows.H3K27Ac.promoter)$overlapTSS,#
         GC = rowData(TFE3.filtered.windows.H3K27Ac.promoter)[,'G|C'])

gl <- lapply(combi, function(x) plotMA(df, x[1], x[2]))
patchwork::wrap_plots(gl, ncol=2)
#### END VIZ MA

dgel <- csaw::asDGEList(TFE3.filtered.windows.H3K27Ac.promoter,
                        genes = as.data.frame(rowData(TFE3.filtered.windows.H3K27Ac.promoter)),
                        group = colData(TFE3.filtered.windows.H3K27Ac.promoter)$SampleType)
colnames(dgel) <- colnames(TFE3.filtered.windows.H3K27Ac.promoter)
rownames(dgel) <- as.character(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter))

moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))

dgel <- estimateDisp(dgel, design=moma)
fit <- glmQLFit(dgel, moma,robust=TRUE)


contr <- makeContrasts(TSC - dKO, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.1)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(TFE3.filtered.windows.H3K27Ac.promoter) <- cbind(rowData(TFE3.filtered.windows.H3K27Ac.promoter), tt)
merged <- mergeWindows(rowRanges(TFE3.filtered.windows.H3K27Ac.promoter), tol=500L, max.width=5000L)

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .1, 'WT',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .1, 'TSC', 'Common'))

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.1,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()

peaks$gene <- str_extract(peaks$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.1 & direction == 'down' & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]
genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))

###


idx <- match(peaks$gene,dKO.vs.TSC$gene)

peaks$RNA_log2fc <- -dKO.vs.TSC$log2FoldChange[idx]
peaks$RNA_p_adj <- dKO.vs.TSC$padj[idx]
df <- as.data.frame(peaks)


df$color <- 'gray'
df$color[match(LRP_genes,df$gene)] <- 'blue'
df$color[match(setdiff(MITF_genes,LRP_genes),df$gene)] <- 'green'
df$color[match(Clear_genes,df$gene)] <- 'red'


pdf('~/Documents/TSC_Paper/RNASeq_H3K27Ac_correlation_TSC_dKO.pdf',height=3, width=4)
  ggplot(subset(df,color!='gray'),aes(x=rep.logFC,y=RNA_log2fc, label = gene, fill = color, color = color)) + 
  geom_point(alpha = 0.5) + geom_text_repel() + 
  scale_color_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) +
  theme_classic()
dev.off()


p1 <- ggplot(subset(df,color!='gray'),aes(rep.logFC,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-2,6)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()


p2 <- ggplot(subset(df,color!='gray'),aes(RNA_log2fc,fill=color)) + 
geom_density(alpha = 0.1,bounds=c(-1,2)) + 
scale_fill_manual(values = c("blue" = "blue", "green" = "green", "red" = "red")) + 
theme_classic()

pdf('~/Documents/TSC_Paper/RNASeq_H3K27Ac_correlation_TSC_dKO_gauss.pdf',height=3, width=4)
  p1 / p2
dev.off()



















# dKO open chromatin in RNA-Seq but transcripts aren't up

df.up <- as.data.frame(subset(peaks, direction == 'up' & FDR < 0.9 & overlap != ""))
df.down <- as.data.frame(subset(peaks, direction == 'down' & FDR < 0.9 & overlap != ""))
df <- as.data.frame(subset(peaks, FDR < 0.9 & overlap != ""))


library(EnhancedVolcano)

enriched.up <- enrichr(df.up$gene, dbs)
enriched.down <- enrichr(df.down$gene, dbs)



PPARA_genes <- human2mouse$V2[match(PPARA_genes,human2mouse$V1)]

df.up$color <- 'gray'
df.up$color[match(MITF_genes,df.up$gene)] <- 'red'

df.down$color <- 'gray'
df.down$color[match(LXR_genes,df.down$gene)] <- 'red'
df.down$color[match(RXR_genes,df.down$gene)] <- 'blue'
df.down$color[match(FOXA2_genes,df.down$gene)] <- 'orange'
df.down$color[match(PPARA_genes,df.down$gene)] <- 'green'

df$color <- 'gray'
df$color[match(MITF_genes,df$gene)] <- 'red'
df$color[match(union(LXR_genes,c(RXR_genes,FOXA2_genes,PPARA_genes)),df$gene)] <- 'blue'


ggplot(df.up,aes(x=rep.logFC,y=RNA_log2fc, label = gene, color=color)) + geom_point() + geom_text_repel()

ggplot(df.down,aes(x=rep.logFC,y=RNA_log2fc, label = gene, color=color)) + geom_point() + geom_text_repel()


ggplot(df,aes(x=rep.logFC,y=RNA_log2fc, label = gene, color=color)) + geom_point() + geom_text_repel()



df_RNA <- data.frame(peaks)
df_RNA <- df_RNA[!is.na(df_RNA$RNA_log2fc) & !is.na(df_RNA$gene),] 


df_RNA <- subset(df_RNA,FDR<0.9)
keyvals <- ifelse(df_RNA$FDR < 0.9 & df_RNA$direction == "up",'red',
	ifelse(df_RNA$FDR < 0.9 & df_RNA$direction == "down",'blue','gray'))

names(keyvals)[keyvals == 'red'] <- 'up'
names(keyvals)[keyvals == 'blue'] <- 'down'
names(keyvals)[keyvals == 'gray'] <- 'none'

df_RNA<-df_RNA[order(df_RNA$FDR,decreasing = TRUE),]

EnhancedVolcano(df_RNA,x = 'RNA_log2fc',y = 'RNA_p_adj',lab = df_RNA$gene,
	pCutoff = 0.05, FCcutoff=0.25,colCustom = keyvals,
	selectLab=df_RNA$gene[keyvals != 'gray'],colAlpha=.4)















### devtools::install_github("js229/Vennerable")


#### ATAC and RNA seq overlap
library(Vennerable)

n_ATAC_up <- length(unique(df.up$gene))  # Total peaks 
n_RNA_up <- length(unique(dKO.vs.TFE3.up$gene))  # Total peaks 

n_overlap <- length(intersect(unique(df.up$gene),unique(dKO.vs.TFE3.up$gene)))

venn_data <- Venn(SetNames = c("ATAC", "RNA"),
                  Weight = c(
                    "10" = n_ATAC_up, # Unique to A
                    "01" = n_RNA_up, # Unique to B
                    "11" = n_overlap         # Intersection
                  ))

# Plot the Venn diagram
plot(venn_data)



#####  Visualize overlap between USF2 peaks and H3K27Ac peaks in CUT&RUN

peaks_USF2=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/USF2_peaks.narrowPeak")
peaks_H3K27Ac=import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/H3K27Ac_peaks.narrowPeak")

USF2_overlap_H3K27Ac_peaks<- subsetByOverlaps(peaks_USF2, peaks_H3K27Ac)
length(peaks_USF2)
length(peaks_H3K27Ac)
length(USF2_overlap_H3K27Ac_peaks)

USF2_summit<- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/USF2_summits.bed")
H3K27Ac_summit<- import("/Volumes/Extreme SSD/CUTRUN/Trimmed/macs2_IgG/dup_120/union/H3K27Ac_summits.bed")


H3K27Ac_summit_500bp_window<- resize(H3K27Ac_summit, width = 500, fix="center")

hits<- findOverlaps(USF2_summit, H3K27Ac_summit_500bp_window)

summit_distance<- distance(USF2_summit[queryHits(hits)], H3K27Ac_summit[subjectHits(hits)])

table(summit_distance)

USF2_summit[queryHits(hits)][summit_distance ==0]

H3K27Ac_summit[subjectHits(hits)][summit_distance ==0]


# Compute signed distances
signed_distance <- function(A, B) {
  # Compute unsigned distance
  dist <- distance(A, B)
  
  # Determine signs based on whether A precedes or follows B
  sign <- ifelse(start(A) < start(B), -1, 1)
  
  # Apply sign to distance
  dist * sign
}


summit_distance<- signed_distance(USF2_summit[queryHits(hits)],
                                  H3K27Ac_summit[subjectHits(hits)])

distance_df<- table(summit_distance) %>%
  tibble::as_tibble() 

distance_df %>%
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  arrange(summit_distance) %>%
  ggplot(aes(x=summit_distance, y = n)) +
  geom_line()

df_binned <- distance_df %>%
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  arrange(summit_distance) %>%
  mutate(bin = floor(summit_distance / 5) * 5) %>%  # Create bins by grouping every 5 bp
  group_by(bin) %>%
  summarise(n = mean(n, na.rm = TRUE))  # Calculate average 'n' for each bin


df_binned %>%
  ggplot(aes(x=bin, y = n)) +
  geom_line() +
  scale_x_continuous(breaks = c(-250, 0, 250)) +
  xlab("distance to the summit \nof TAZ peaks (bp)") +
  ylab("peak density") +
  theme_classic(base_size = 14) 




#####  Visualize overlap between TFE3 peaks and dKO peaks in ATAC
TFE3_peaks<- import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/TFE3/TFE3_peaks.narrowPeak")
dKO_peaks<- import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/dKO/dKO_peaks.narrowPeak")

TFE3_overlap_dKO_peaks<- subsetByOverlaps(TFE3_peaks, dKO_peaks)
length(TFE3_peaks)
length(dKO_peaks)
length(TFE3_overlap_dKO_peaks)


TFE3_summit<- import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/TFE3/TFE3_summits.bed")
dKO_summit<- import("~/Documents/TSC_ATAC/Merged/Trimmed/macs2/dKO/dKO_summits.bed")


dKO_500bp_window<- resize(dKO_summit, width = 500, fix="center")

hits<- findOverlaps(TFE3_summit, dKO_500bp_window)

summit_distance<- distance(TFE3_summit[queryHits(hits)], dKO_summit[subjectHits(hits)])

table(summit_distance)

TFE3_summit[queryHits(hits)][summit_distance ==0]

dKO_summit[subjectHits(hits)][summit_distance ==0]


# Compute signed distances
signed_distance <- function(A, B) {
  # Compute unsigned distance
  dist <- distance(A, B)
  
  # Determine signs based on whether A precedes or follows B
  sign <- ifelse(start(A) < start(B), -1, 1)
  
  # Apply sign to distance
  dist * sign
}


summit_distance<- signed_distance(TFE3_summit[queryHits(hits)],
                                  dKO_summit[subjectHits(hits)])

distance_df<- table(summit_distance) %>%
  tibble::as_tibble() 

distance_df %>%
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  arrange(summit_distance) %>%
  ggplot(aes(x=summit_distance, y = n)) +
  geom_line()

df_binned <- distance_df %>%
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  arrange(summit_distance) %>%
  mutate(bin = floor(summit_distance / 5) * 5) %>%  # Create bins by grouping every 5 bp
  group_by(bin) %>%
  summarise(n = mean(n, na.rm = TRUE))  # Calculate average 'n' for each bin


df_binned %>%
  ggplot(aes(x=bin, y = n)) +
  geom_line() +
  scale_x_continuous(breaks = c(-250, 0, 250)) +
  xlab("distance to the summit \nof TAZ peaks (bp)") +
  ylab("peak density") +
  theme_classic(base_size = 14) 


### Plot TF peaks by activity
#promoters are TF peaks overlapping with H3K4me3 peaks
#active enhancers are TF peaks overlapping with H3K4me1 and H3K27ac peaks
#inactive enhancers are TF peaks overlapping with H3K4me1 but not H3K37ac peaks

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

H3K27ac<- import("~/Documents/TSC_ATAC/ExistingData/H3K27Ac/ENCFF395EHX.bed",
	format = "BED",extraCols = extraCols_narrowPeak)
H3K4me1<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me1/ENCFF596ORE.bed",
	format = "BED",extraCols = extraCols_narrowPeak)
H3K4me3<- import("~/Documents/TSC_ATAC/ExistingData/H3K4me3/ENCFF150KSB.bed",
	format = "BED",extraCols = extraCols_narrowPeak)
seqlevelsStyle(H3K27ac) <- "NCBI"
seqlevelsStyle(H3K4me1) <- "NCBI"
seqlevelsStyle(H3K4me3) <- "NCBI"


active_enhancers<- subsetByOverlaps(H3K4me1, H3K27ac)
inactive_enhancers<- subsetByOverlaps(H3K4me1, H3K27ac, invert=TRUE)
promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)



n_active_enhancers<- subsetByOverlaps(peaks_USF2,
                                      active_enhancers) %>% length()

n_inactive_enhancers<- subsetByOverlaps(peaks_USF2,
                                        inactive_enhancers) %>% length()

n_promoters<- subsetByOverlaps(peaks_USF2, 
                               promoters) %>% length()

n_unclassified<- length(peaks_USF2) - n_active_enhancers -
  n_inactive_enhancers - n_promoters


annotation_df<- data.frame(category = c("active_enhancers", "inactive_enhancers",
                        "promoters", "unclassified"),
peak_number = c(n_active_enhancers, n_inactive_enhancers, 
                n_promoters, n_unclassified))


annotation_df

annotation_df$category<- factor(annotation_df$category, 
                                levels = c("promoters", "active_enhancers",
                                           "inactive_enhancers", "unclassified"))

colors<- c("#8D1E0F", "#F57D2B", "#FADAC4", "#D4DADA")

ggplot(annotation_df, aes(x = "", y = peak_number, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Remove unnecessary axes
  labs(title = "YAP/TAZ/TEAD4 peaks") +
  scale_fill_manual(values = colors)


### HEATMAPS

enhancers<- subsetByOverlaps(H3K4me1, H3K4me3, invert=TRUE)

promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

TFE3_enhancers<- subsetByOverlaps(TFE3_peaks, enhancers) 

TFE3_promoters<- subsetByOverlaps(TFE3_peaks, promoters) 

TFE3_summit_enhancer<- TFE3_summit[TFE3_summit$name %in% TFE3_enhancers$name]
TFE3_summit_promoter<- TFE3_summit[TFE3_summit$name %in% TFE3_promoters$name]

# combine them
anchors<- c(TFE3_summit_promoter, TFE3_summit_enhancer) 

TFE3_1741_bw<- import("~/Documents/TSC_ATAC/1741_TFE3.bigWig")
TFE3_1746_bw<- import("~/Documents/TSC_ATAC/1746_TFE3.bigWig")
TFE3_1751_bw<- import("~/Documents/TSC_ATAC/1741_TFE3.bigWig")


# BiocManager::install("EnrichedHeatmap")
library(EnrichedHeatmap)

# extend 1000 bp on each side and use 50bp bin
mat1<- normalizeToMatrix(TFE3_1741_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

mat2<- normalizeToMatrix(TFE3_1746_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

mat3<- normalizeToMatrix(TFE3_1751_bw, anchors, value_column = "score",
                         extend= 1000, mean_mode = "w0", w=50)

quantile(mat1, c(0.1,0.25,0.5,0.9,1))
quantile(mat2, c(0.1,0.25,0.5,0.9,1))
quantile(mat3, c(0.1,0.25,0.5,0.9,1))

col_fun<- circlize::colorRamp2(c(0, 1.5), c("white", "red"))


partition<- c(rep("promoters", length(TFE3_promoters)),
              rep("enhancers", length(TFE3_enhancers)))

# change the factor level so promoters come first
partition<- factor(partition, levels=c("promoters", "enhancers"))

partition_hp<- Heatmap(partition, col=structure(2:3, names = c("promoters", "enhancers")), 
        name = "partition",
        show_row_names = FALSE, width=unit(3,'mm'))

partition_hp

ht_list<- partition_hp +
  EnrichedHeatmap(mat1, pos_line = FALSE, column_title="1", name = "1", col=col_fun) +
  EnrichedHeatmap(mat2, pos_line = FALSE, column_title="2", name = "2", col=col_fun) +
  EnrichedHeatmap(mat3, pos_line = FALSE, column_title="3", name = "3", col=col_fun)

draw(ht_list, split= partition, main_heatmap =2)




#H3K4me1 signal around peaks
YAP1_mean<- colMeans(mat1)
TAZ_mean<- colMeans(mat2)
TEAD4_mean<- colMeans(mat3)

bind_rows(YAP1_mean, TAZ_mean, TEAD4_mean) %>%
  mutate(factor = factor(c("YAP1", "TAZ", "TEAD4"), levels = c("YAP1", "TAZ", "TEAD4"))) %>%
  dplyr::select(factor, everything()) %>%
  tidyr::pivot_longer(-factor) %>%
  mutate(name = factor(name, levels = c(paste0("u",1:20), paste0("d", 1:20)))) %>%
  ggplot(aes(x=name, y=value)) +
  geom_line(aes(color = factor, group=factor)) +
  scale_x_discrete(breaks=c("u1", "d1", "d20"), labels = c("-1kb", "0", "1kb")) +
  scale_color_manual(values = c("#DA9195", "#E07B78", "#605D7D")) +
  theme_classic(base_size = 14) +
  ylab("RPKM") +
  xlab("")







































background <- windowCounts(bam_alt, bin=TRUE, width=10000L, param=param)
colData(background) <- cbind(colData(background), align_alt)


filterStat <- filterWindowsGlobal(windows_alt, background)
str(filterStat)

globalBG <- median( filterStat$back.abundances )
globalBG

head(filterStat$abundances - globalBG)


keep <- filterStat$filter > log2(3)
sum(keep)
mean(keep)
min(filterStat$abundances[keep])
globalBG + log2(3)

ggplot(data = tibble(x=filterStat$abundances),
       mapping = aes(x)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = globalBG, color='red') +
    geom_vline(xintercept = globalBG + log2(3), color='blue') +
    labs(x="Abundance", y="Frequency", title="Histogram of Window Abundances") +
    theme_bw(base_size=15)


rowData(windows_alt)$enrichment <- filterStat$filter

filtered.windows <- windows_alt[keep,]

rtracklayer::export(rowRanges(filtered.windows),
                    '~/Documents/TSC_ATAC/filtered_global.bed')












# Differential accessibility analysis
binned <- windowCounts(bam_alt, bin=TRUE, width=10000, param=param)
filtered.windows <- normFactors(binned, se.out=filtered.windows)


library(edgeR)
adj.counts <- cpm(asDGEList(binned), log=TRUE)
normfacs <- filtered.windows$norm.factors

par(mfrow=c(2, 3), mar=c(5, 4, 2, 1.5))
for (i in seq_len(length(bam_alt)-1)) {
    cur.x <- adj.counts[,2]
    cur.y <- adj.counts[,1+i]
    smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y,
        xlab="A", ylab="M", main=paste("1 vs", i+1))
    all.dist <- diff(log2(normfacs[c(i+1, 1)]))
    abline(h=all.dist, col="red")
}




logCPM <- csaw::calculateCPM(filtered.windows, log=TRUE)


dgel <- csaw::asDGEList(filtered.windows,
                        genes = as.data.frame(rowData(filtered.windows)),
                        group = colData(filtered.windows)$SampleType,
                        norm.factors=filtered.windows$norm.factors)
colnames(dgel) <- colnames(filtered.windows)
rownames(dgel) <- as.character(rowRanges(filtered.windows))
dgel

plotMDS(dgel, top=1000, gene.selection="common")





#filtered.windows <- normOffsetsNQ(filtered.windows)


filtered.windows <- normOffsets(filtered.windows,span=0.3)

filtered.windows <- normFactors(filtered.windows)


logCPM <- csaw::calculateCPM(filtered.windows, log=TRUE)

ct <- split(colnames(logCPM), colData(filtered.windows)$SampleType)

df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(filtered.windows))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(filtered.windows)$overlapTSS,
         GC = rowData(filtered.windows)[,'G|C'])

patchwork::wrap_plots(lapply(combi, function(x) plotMA(df, x[1], x[2])), ncol=2)
patchwork::wrap_plots(lapply(combi, function(x) plotGC(df, x[1], x[2])), ncol=2)

plotMA(df,'WT','TSC')
plotGC(df,'WT','TFE3')


dgel <- csaw::asDGEList(filtered.windows,
                        genes = as.data.frame(rowData(filtered.windows)),
                        group = colData(filtered.windows)$SampleType)
colnames(dgel) <- colnames(filtered.windows)
rownames(dgel) <- as.character(rowRanges(filtered.windows))
dgel

plotMDS(dgel, top=1000, gene.selection="common")


moma <- model.matrix(~ 0 + group, dgel$samples)
colnames(moma) <- sub('group','',colnames(moma))
moma


dgel <- estimateDisp(dgel, design=moma)



fit <- glmQLFit(dgel, moma)


contr <- makeContrasts(TSC - WT, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.05)) #Rec to set LFC ar default


### Plot MA plot
ct <- split(colnames(logCPM), colData(filtered.windows)$SampleType)

df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(filtered.windows))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(filtered.windows)$overlapTSS,
         GC = rowData(filtered.windows)[,'G|C'])

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.05))

ggplot(df,mapping = aes(x = .data[['TSC']] + .data[['WT']],y = .data[['TSC']] - .data[['WT']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'TSC', 'WT')) +
    theme_bw(base_size = 15)






###########
# TSC VS WT
###########
filtered.windows <- windows_alt[keep,]
filtered.windows <- normOffsets(filtered.windows)

logCPM <- csaw::calculateCPM(filtered.windows, log=TRUE, use.offsets = TRUE)


contr <- makeContrasts(TSC - WT, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.05)) #Rec to set LFC ar default


### Plot MA plot
ct <- split(colnames(logCPM), colData(filtered.windows)$SampleType)

df <- vapply(ct, function(x) {
  rowMeans(logCPM[,x,drop=FALSE])
}, rep(0, nrow(filtered.windows))) |>
  as_tibble() |>
  mutate(overlapTSS = rowData(filtered.windows)$overlapTSS,
         GC = rowData(filtered.windows)[,'G|C'])

df$signif.TSC.vs.WT <- as.logical(decideTests(res, p.value = 0.05))

ggplot(df,mapping = aes(x = .data[['WT']],y = .data[['TSC']] - .data[['WT']])) +
    geom_hex(bins = 100, aes(fill=after_stat(density)^(1/16)), show.legend = FALSE) +
    rasterize( geom_point(data = subset(df, signif.TSC.vs.WT), col='red', pch='.') ) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "A", y="M", title=sprintf("%s vs. %s", 'TSC', 'WT')) +
    theme_bw(base_size = 15)




tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(filtered.windows) <- cbind(rowData(filtered.windows), tt)


merged <- mergeWindows(rowRanges(filtered.windows), tol=150L)
str(merged,1)

merged$region

summary(width(merged$region))

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .05, 'TSC',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .05, 'WT', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/peaks_TSC_vs_WT.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/peaks_TSC_vs_WT.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.05,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()


#temp <- str_extract(subset(data,Significant == TRUE & str_sub(data$overlap,-1) == 'P')$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.05 & direction == 'up' & rep.logFC >0 & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]

genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))


tsc.vs.wt <- read.csv('~/Documents/TSC_RNASeq/TSC.vs.WT.csv',sep=',')
tsc.vs.wt <- subset(tsc.vs.wt,padj<0.05)
tsc.vs.wt.up <- subset(tsc.vs.wt,log2FoldChange>0)
rna_up <- intersect(tsc.vs.wt.up$gene,genes_all)






















###########
# TSC VS TFE3
###########


contr <- makeContrasts(TSC - TFE3, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.05)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(filtered.windows) <- cbind(rowData(filtered.windows), tt)


merged <- mergeWindows(rowRanges(filtered.windows), tol=150L)
str(merged,1)

merged$region

summary(width(merged$region))

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .05, 'TSC',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .05, 'TFE3', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/peaks_TSC_vs_TFE3.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/peaks_TSC_vs_TFE3.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.05,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()


#temp <- str_extract(subset(data,Significant == TRUE & str_sub(data$overlap,-1) == 'P')$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.05 & direction == 'up' & rep.logFC >0 & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]

genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))


tsc.vs.tfe3 <- read.csv('~/Documents/TSC_RNASeq/TSC.vs.TFE3.csv',sep=',')
tsc.vs.tfe3 <- subset(tsc.vs.tfe3,padj<0.05)
tsc.vs.tfe3.up <- subset(tsc.vs.tfe3,log2FoldChange>0)
rna_up <- intersect(tsc.vs.tfe3.up$gene,genes_all)




###########
# dKO VS TFE3
###########


contr <- makeContrasts(dKO - TFE3, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.05)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(filtered.windows) <- cbind(rowData(filtered.windows), tt)


merged <- mergeWindows(rowRanges(filtered.windows), tol=150L)
str(merged,1)

merged$region

summary(width(merged$region))

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .05, 'dKO',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .05, 'TFE3', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/peaks_dKO_vs_TFE3.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/peaks_dKO_vs_TFE3.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.05,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()


#temp <- str_extract(subset(data,Significant == TRUE & str_sub(data$overlap,-1) == 'P')$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.05 & direction == 'up' & rep.logFC >0 & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]

genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))


dKO.vs.TFE3 <- read.csv('~/Documents/TSC_RNASeq/dKO.vs.TFE3.csv',sep=',')
dKO.vs.TFE3 <- subset(tsc.vs.tfe3,padj<0.05)
dKO.vs.TFE3.up<- subset(dKO.vs.TFE3,log2FoldChange>0)
rna_up <- intersect(dKO.vs.TFE3.up$gene,genes_all)


###########
# dKO VS WT
###########


contr <- makeContrasts(dKO - WT, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.05)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(filtered.windows) <- cbind(rowData(filtered.windows), tt)


merged <- mergeWindows(rowRanges(filtered.windows), tol=150L)
str(merged,1)

merged$region

summary(width(merged$region))

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .05, 'dKO',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .05, 'WT', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/peaks_dKO_vs_WT.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/peaks_dKO_vs_WT.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.05,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()


#temp <- str_extract(subset(data,Significant == TRUE & str_sub(data$overlap,-1) == 'P')$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.05 & direction == 'up' & rep.logFC >0 & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]

genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))


dKO.vs.WT <- read.csv('~/Documents/TSC_RNASeq/dKO.vs.WT.csv',sep=',')
dKO.vs.WT <- subset(tsc.vs.tfe3,padj<0.05)
dKO.vs.WT.up<- subset(dKO.vs.WT,log2FoldChange>0)
rna_up <- intersect(dKO.vs.WT.up$gene,genes_all)


###########
# TFE3 VS WT
###########


contr <- makeContrasts(TFE3 - WT, levels = colnames(moma))

res <- glmQLFTest(fit, contrast = contr)

table(decideTests(res, p.value = 0.05)) #Rec to set LFC ar default


tt <- topTags(res, n=Inf, sort.by="none")$table[,c('logFC','logCPM','F','PValue','FDR')]

tt |>
    as_tibble() |>
    mutate(Windows = rownames(dgel), .before=1) |>
    arrange(FDR) |>
    head()

rowData(filtered.windows) <- cbind(rowData(filtered.windows), tt)


merged <- mergeWindows(rowRanges(filtered.windows), tol=150L)
str(merged,1)

merged$region

summary(width(merged$region))

anno <- detailRanges(merged$region,
                     txdb=TxDb.Mmusculus.UCSC.mm10.knownGene,
                     orgdb=org.Mm.eg.db, promoter=c(1500, 500), dist=10*1e3L)

merged$region$overlap <- anno$overlap
merged$region$left <- anno$left
merged$region$right <- anno$right
merged$region

region.stats <- combineTests(merged$id, tt) |>
    as_tibble() |>
    mutate(rep.logCPM = tt[rep.test,'logCPM'],
           Contrast = colnames(contr))

region.stats |>
    head()


peaks <- merged$region 

mcols(peaks) <- c(mcols(peaks), region.stats)

peaks$name <- ifelse(peaks$rep.logFC > 1 & peaks$FDR < .05, 'TFE3',
                     ifelse(peaks$rep.logFC < (-1) & peaks$FDR < .05, 'WT', 'Common'))

write.table(as.data.frame(peaks), file='~/Documents/TSC_ATAC/peaks_TFE3_vs_WT.tsv',
            row.names=FALSE, sep='\t')
rtracklayer::export(peaks, '~/Documents/TSC_ATAC/peaks_TFE3_vs_WT.bed')

table(peaks$name, peaks$direction)


data <- mcols(peaks) |>
    as_tibble() |>
    mutate(logFC = rep.logFC,
           Significant = FDR < 0.05,
           Annotated = overlap != "")

ggplot(data = data,
       mapping = aes(x = logFC,
                     y = -log10(PValue),
                     color = Significant,
                     shape = Annotated,
                     label = overlap)) +
    rasterize( geom_point() ) +
    geom_text_repel(data = data |> filter(Annotated & Significant & abs(logFC) > 1.5)) +
    geom_vline(xintercept = 0, lty = 2) +
    theme_bw()


#temp <- str_extract(subset(data,Significant == TRUE & str_sub(data$overlap,-1) == 'P')$overlap,"[^:]+")


areas <- subset(peaks,FDR < 0.05 & direction == 'up' & rep.logFC >0 & overlap != "")$overlap
proms <- str_sub(areas,-2) == 'PE' | str_sub(areas,-2) == 'PI' | str_sub(areas,-1) == 'P'
exons <- str_sub(areas,-1) == 'E'
introns <- str_sub(areas,-1) == 'I'
areas_proms <- areas[proms]
areas_introns <- areas[introns]
areas_exons <- areas[exons]

genes_proms <- str_extract(areas_proms,"[^:]+")
genes_exons <- str_extract(areas_exons,"[^:]+")
genes_introns <- str_extract(areas_introns,"[^:]+")
genes_all <- union(genes_proms,c(genes_exons,genes_introns))


dKO.vs.WT <- read.csv('~/Documents/TSC_RNASeq/TSC.vs.WT.csv',sep=',')
dKO.vs.WT <- subset(tsc.vs.tfe3,padj<0.05)
dKO.vs.WT.up<- subset(dKO.vs.WT,log2FoldChange>0)
rna_up <- intersect(dKO.vs.WT.up$gene,genes_all)







#LXR and RXR are low
cat(unique(genes_proms),sep='\n')

# Too few
cat(setdiff(genes_proms,c(genes_introns,genes_exons)),sep='\n')


# Contains LXR and RXR
cat(unique(genes_exons),sep='\n')

# Pretty clean enrichment of RXR, LXR, PPARA, CTCF
cat(setdiff(genes_exons,genes_proms),sep='\n')

# Pretty clean enrichment of RXR, LXR, PPARA, CTCF
cat(setdiff(genes_exons,c(genes_proms,genes_introns)),sep='\n')



# MITF stars popping up here, also LEF1, WT1, FOXA2; FOXA2 is a pioneer TF
# H3K36me3 and H3K27me3
cat(unique(genes_introns),sep='\n')
cat(setdiff(genes_introns,genes_proms),sep='\n')
cat(setdiff(genes_introns,c(genes_proms,genes_exons)),sep='\n')





