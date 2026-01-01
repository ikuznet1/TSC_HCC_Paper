library(reticulate)
library(ggfortify)
library(edgeR)
library(RColorBrewer)
library(EnhancedVolcano)
library(DESeq2)
library(tximport)
library(biomaRt)
library("sva")


##################
#### Fig 2B ######
##################

### Gene level

path = '~/Documents/TSC_RNASeq'
tx2gene <- read.table(paste0(path,'/t2g_clean.csv'),sep=',')
meta <- read.csv(paste0(path,'/metadata.csv'))
rownames(meta)<-meta$Number

meta <- meta[order(meta$Number),]


files <- paste(path,list.files(path,pattern = "\\.h5$",recursive=TRUE),sep='/')

txi.kallisto <- tximport(files, type = "kallisto", txOut = FALSE, tx2gene=tx2gene)
bulk <- txi.kallisto$counts

ddsSE <- DESeqDataSetFromTximport(txi.kallisto,meta,design=~Group + Batch)
ddsSE <- estimateSizeFactors(ddsSE)
idx <- rowSums( counts(ddsSE, normalized=TRUE) >= 5 ) >= 3
ddsSE <- ddsSE[idx,]

normalized_counts <- counts(ddsSE, normalized=TRUE)


ddsSE <- DESeq(ddsSE)
TSC.vs.WT<- lfcShrink(ddsSE,contrast=c('Group','TSC KO','WT'), type="ashr")
TFE3.vs.WT<- lfcShrink(ddsSE,contrast=c('Group','TFE3 KO','WT'), type="ashr")
TSC.vs.TFE3<- lfcShrink(ddsSE,contrast=c('Group','TSC KO','TFE3 KO'), type="ashr")
dKO.vs.WT<- lfcShrink(ddsSE,contrast=c('Group','dKO','WT'), type="ashr")
dKO.vs.TSC<- lfcShrink(ddsSE,contrast=c('Group','dKO','TSC KO'), type="ashr")
dKO.vs.TFE3<- lfcShrink(ddsSE,contrast=c('Group','dKO','TFE3 KO'), type="ashr")


vstSE <- vst(ddsSE,blind = FALSE)


mat <- assay(vstSE)
mm <- model.matrix(~Group, colData(vstSE))
mat <- limma::removeBatchEffect(mat, covariates=colData(vstSE)[,3], design=mm)
assay(vstSE) <- mat


pdf('~/Documents/TSC_Paper/RNASeq_pca.pdf',width=10,height=4)
plotPCA(vstSE,intgroup=c("Group"),ntop=5000) + theme_classic() + 
theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),legend.title=element_text(size=24),legend.text=element_text(size=24)) + 
labs(color='Disease') + geom_point(size=3.5)
dev.off()


##################
#### Fig 2C ######
##################


gene_an <- read.csv(paste0(path,'/t2g.txt'), sep='\t', na.strings="",header=FALSE)
gene_an <- gene_an[match(unique(gene_an$V2),gene_an$V2),]

rownames(gene_an) <- gene_an$V2

TSC.vs.WT$gene <- gene_an[rownames(TSC.vs.WT),]$V3
TFE3.vs.WT$gene <- gene_an[rownames(TFE3.vs.WT),]$V3
TSC.vs.TFE3$gene <- gene_an[rownames(TSC.vs.TFE3),]$V3
dKO.vs.WT$gene <- gene_an[rownames(dKO.vs.WT),]$V3
dKO.vs.TSC$gene <- gene_an[rownames(dKO.vs.TSC),]$V3
dKO.vs.TFE3$gene <- gene_an[rownames(dKO.vs.TFE3),]$V3


write.csv(TSC.vs.WT,'~/Documents/TSC_Paper/TSC.vs.WT.csv')
write.csv(TFE3.vs.WT,'~/Documents/TSC_Paper/TFE3.vs.WT.csv')
write.csv(TSC.vs.TFE3,'~/Documents/TSC_Paper/TSC.vs.TFE3.csv')
write.csv(dKO.vs.WT,'~/Documents/TSC_Paper/dKO.vs.WT.csv')
write.csv(dKO.vs.TSC,'~/Documents/TSC_Paper/dKO.vs.TSC.csv')
write.csv(dKO.vs.TFE3,'~/Documents/TSC_Paper/dKO.vs.TFE3.csv')

clear <- read.csv('~/Documents/TSC_Paper/CLEAR_genes.csv',header=F)
clear <- clear$V1

library(EnhancedVolcano)


# create custom key-value pairs for CLEAR components

TSC.vs.WT <- TSC.vs.WT[!is.na(TSC.vs.WT$padj),]
TSC.vs.WT[TSC.vs.WT$padj < 1e-20,]$padj = 1e-20
keyvals <- ifelse(TSC.vs.WT$gene %in% clear, 'red','blue')
names(keyvals)[keyvals == 'red'] <- 'CLEAR'
names(keyvals)[keyvals == 'gray'] <- 'Not CLEAR'
opacity <- ifelse(TSC.vs.WT$gene %in% clear, 1,0.1)

pdf('~/Documents/TSC_Paper/TSC.vs.WT_volcano.pdf',width=6,height=5)
EnhancedVolcano(TSC.vs.WT,lab=TSC.vs.WT$gene,FCcutoff=0.25,
	colCustom = keyvals,selectLab = clear,drawConnectors = TRUE,
	x='log2FoldChange',y='padj',pCutoff = 0.05,xlim=c(-2,2),colAlpha = opacity) + coord_flip()
dev.off()

dKO.vs.WT <- dKO.vs.WT[!is.na(dKO.vs.WT$padj),]
dKO.vs.WT[dKO.vs.WT$padj < 1e-20,]$padj = 1e-20
keyvals <- ifelse(dKO.vs.WT$gene %in% clear, 'red','blue')
names(keyvals)[keyvals == 'red'] <- 'CLEAR'
names(keyvals)[keyvals == 'gray'] <- 'Not CLEAR'
opacity <- ifelse(dKO.vs.WT$gene %in% clear, 1,0.1)

pdf('~/Documents/TSC_Paper/dKO.vs.WT_volcano.pdf',width=6,height=5)
EnhancedVolcano(dKO.vs.WT,lab=dKO.vs.WT$gene,FCcutoff=0.25,
	colCustom = keyvals,selectLab = clear,drawConnectors = TRUE,
	x='log2FoldChange',y='padj',pCutoff = 0.05,xlim=c(-4,4),colAlpha = opacity) + coord_flip()
dev.off()

dKO.vs.TSC <- dKO.vs.TSC[!is.na(dKO.vs.WT$padj),]
dKO.vs.TSC[dKO.vs.TSC$padj < 1e-20,]$padj = 1e-20
keyvals <- ifelse(dKO.vs.TSC$gene %in% clear, 'red','blue')
names(keyvals)[keyvals == 'red'] <- 'CLEAR'
names(keyvals)[keyvals == 'gray'] <- 'Not CLEAR'
opacity <- ifelse(dKO.vs.TSC$gene %in% clear, 1,0.1)

pdf('~/Documents/TSC_Paper/dKO.vs.TSC_volcano.pdf',width=6,height=5)
EnhancedVolcano(dKO.vs.TSC,lab=dKO.vs.TSC$gene,FCcutoff=0.25,
	colCustom = keyvals,selectLab = clear,drawConnectors = TRUE,
	x='log2FoldChange',y='padj',pCutoff = 0.05,xlim=c(-4,4),ylim=c(0,25),colAlpha = opacity) + coord_flip()
dev.off()


TFE3.vs.WT <- TFE3.vs.WT[!is.na(TFE3.vs.WT$padj),]
TFE3.vs.WT[TFE3.vs.WT$padj < 1e-20,]$padj = 1e-20
keyvals <- ifelse(TFE3.vs.WT$gene %in% clear, 'red','blue')
names(keyvals)[keyvals == 'red'] <- 'CLEAR'
names(keyvals)[keyvals == 'gray'] <- 'Not CLEAR'
opacity <- ifelse(TFE3.vs.WT$gene %in% clear, 1,0.1)

pdf('~/Documents/TSC_Paper/TFE3.vs.WT_volcano.pdf',width=6,height=5)
EnhancedVolcano(TFE3.vs.WT,lab=TFE3.vs.WT$gene,FCcutoff=0.25,
	colCustom = keyvals,selectLab = clear,drawConnectors = TRUE,
	x='log2FoldChange',y='padj',pCutoff = 0.05,xlim=c(-4,4),colAlpha = opacity) + coord_flip()
dev.off()

##################
#### Fig 2D ######
##################

library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)


TSC.vs.WT.names.up <- subset(TSC.vs.WT,padj < 0.05 & log2FoldChange > 0.25 & !is.na(gene))$gene
dKO.vs.TFE3.names.up <- subset(dKO.vs.TFE3,padj < 0.05 & log2FoldChange > 0.25 & !is.na(gene))$gene


websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("Enrichr") # Human genes   
}



dbs <- c("ChEA_2022","Reactome_2024","GO_Biological_Process_2023")


enriched <- enrichr(TSC.vs.WT.names.up, dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.05,]
p1<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:8),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:8),]$Term,"\\s+"), `[`, 1))) + theme(axis.text=element_text(colour="black"))

enriched <- enrichr(dKO.vs.TFE3.names.up, dbs)
enriched[[1]] <- enriched[[1]][enriched[[1]]$Adjusted.P.value < 0.05,]
p2<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:8),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:8),]$Term,"\\s+"), `[`, 1))) + theme(axis.text=element_text(colour="black"))


enriched <- enrichr(TSC.vs.WT.names.up, dbs)
enriched[[3]] <- enriched[[3]][enriched[[3]]$Adjusted.P.value < 0.05,]
p3<- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:8),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:8),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))

#enriched <- enrichr(dKO.vs.TFE3.names.up, dbs)
#enriched[[3]] <- enriched[[3]][enriched[[3]]$Adjusted.P.value < 0.05,]
#p4<- ggplot(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:5),], (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="P value",size="Overlap") + theme_classic()  + ggtitle('GO Biological Process Up') + scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[3]][order(enriched[[3]]$Combined.Score,decreasing=T),][rev(1:8),]$Term," \\(GO"), `[`, 1))) + theme(axis.text=element_text(colour="black"))

pdf('~/Documents/TSC_Paper/RNASeq_enrichments.pdf',width=6,height=8)

p1 / p2 / p3

dev.off()






