#######################################
#############  FIGURE 1B  #############
#######################################
mageck test -k MageckAnalysisStart_no_zero.txt -t 1,2,3 -c 0 -n chow_vs_lib --norm-method control --remove-zero control  --pdf-report --control-sgrna ./ctrl_sgrna_nonzero.txt --variance-estimation-samples 1,2,3 --gene-test-fdr-threshold  0.05 --additional-rra-parameters '--min-number-goodsgrna 2'
mageck test -k MageckAnalysisStart_no_zero.txt -t 4,5,6 -c 1,2,3 -n cdaa_vs_chow --norm-method control --remove-zero both --pdf-report --control-sgrna ./ctrl_sgrna_nonzero.txt --variance-estimation-samples 1,2,3,4,5,6 --gene-test-fdr-threshold  0.1 --additional-rra-parameters '--min-number-goodsgrna 2' --gene-lfc-method alphamedian
mageck test -k MageckAnalysisStart_no_zero_clean_cdaa.txt -t 4,5,6 -c 1,2,3 -n cdaa_vs_chow_clean --norm-method control --remove-zero both --pdf-report --control-sgrna ./ctrl_sgrna_nonzero.txt --variance-estimation-samples 1,2,3,4,5,6 --gene-test-fdr-threshold  0.1 --additional-rra-parameters '--min-number-goodsgrna 2' --gene-lfc-method alphamedian
mageck test -k MageckAnalysisStart_no_zero.txt -t 4,5,6 -c 0 -n cdaa_vs_lib --norm-method control --remove-zero control --pdf-report --control-sgrna ./ctrl_sgrna_nonzero.txt --variance-estimation-samples 4,5,6 --gene-test-fdr-threshold  0.05 --additional-rra-parameters '--min-number-goodsgrna 2'
mageck test -k MageckAnalysisStart_no_zero.txt -t 1,2,3,4,5,6 -c 0 -n combined_vs_lib --norm-method control --remove-zero control --pdf-report --control-sgrna ./ctrl_sgrna_nonzero.txt --variance-estimation-samples 1,2,3,4,5,6 --gene-test-fdr-threshold  0.05 --additional-rra-parameters '--min-number-goodsgrna 2'


mageck test -k MageckAnalysisStart_no_zero.txt -t 1 -c 0 -n chow1_vs_lib --norm-method control --remove-zero control  --pdf-report --control-sgrna ./ctrl_sgrna_nonzero.txt --variance-estimation-samples 1,2,3 --gene-lfc-method median
mageck test -k MageckAnalysisStart_no_zero.txt -t 2 -c 0 -n chow2_vs_lib --norm-method control --remove-zero control  --pdf-report --control-sgrna ./ctrl_sgrna_nonzero.txt --variance-estimation-samples 1,2,3 --gene-lfc-method median
mageck test -k MageckAnalysisStart_no_zero.txt -t 3 -c 0 -n chow3_vs_lib --norm-method control --remove-zero control  --pdf-report --control-sgrna ./ctrl_sgrna_nonzero.txt --variance-estimation-samples 1,2,3 --gene-lfc-method median

mageck test -k MageckAnalysisStart_no_zero.txt -t 4 -c 0 -n cdaa1_vs_lib --norm-method control --remove-zero control  --pdf-report --control-sgrna ./ctrl_sgrna_nonzero.txt --variance-estimation-samples 4,5,6 --gene-lfc-method median
mageck test -k MageckAnalysisStart_no_zero.txt -t 5 -c 0 -n cdaa2_vs_lib --norm-method control --remove-zero control  --pdf-report --control-sgrna ./ctrl_sgrna_nonzero.txt --variance-estimation-samples 4,5,6 --gene-lfc-method median
mageck test -k MageckAnalysisStart_no_zero.txt -t 6 -c 0 -n cdaa3_vs_lib --norm-method control --remove-zero control  --pdf-report --control-sgrna ./ctrl_sgrna_nonzero.txt --variance-estimation-samples 4,5,6 --gene-lfc-method median


# Pooled analysis


#cat(names(tail(a,n=100)),sep='\n') Top 100 genes give interesting results

#######################################
#############  FIGURE 1B  #############
#######################################

library(MAGeCKFlute)
library(clusterProfiler)
library(ggplot2)


file1 = file.path("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow_vs_lib.gene_summary.txt")
file2 = file.path("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow_vs_lib.sgrna_summary.txt")
genes<-read.csv(file1,sep='\t')

gdata = ReadRRA(file1,score='rra')
gdata$num <- genes$num
sdata = ReadsgRRA(file2)
gdata$HumanGene = TransGeneID(gdata$id, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")
sdata$HumanGene = TransGeneID(sdata$Gene, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")

## Remove missing or duplicate human genes
idx = duplicated(gdata$HumanGene)|is.na(gdata$HumanGene)
gdata = gdata[!idx, ]
#depmap_similarity = ResembleDepmap(gdata, symbol = "HumanGene", score = "Score")
#head(depmap_similarity)

gdata2<-OmitCommonEssential(gdata, symbol = "HumanGene")
alt_essential <- setdiff(gdata$id,gdata2$id)

essential <- read.table('~/Documents/PinAPL/mageck/CountFiles/EssentialGenes.csv',header=FALSE,sep='\t')
human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')
idx_match<- match(essential$V1,human2mouse$human_name)
essential_mouse <- human2mouse$mouse_name[idx_match]
essential_mouse <- essential_mouse[!is.na(essential_mouse)]


#TSG <- read.table('~/Documents/PinAPL/mageck/CountFiles/TSG_Genes.csv',header=TRUE,sep=',')
TSG <- read.table('~/Documents/TSC_Paper/cancerGeneList.csv',header=TRUE,sep=',')
TSG <- subset(TSG,Is.Tumor.Suppressor.Gene=='Yes')
idx_match<- match(TSG$Hugo.Symbol,human2mouse$human_name)
TSG_mouse <- human2mouse$mouse_name[idx_match]
TSG$mouse <- TSG_mouse
TSG <- TSG[!is.na(TSG$mouse),]

oncogene <- read.table('~/Documents/TSC_Paper/cancerGeneList.csv',header=TRUE,sep=',')
oncogene <- subset(oncogene,Is.Tumor.Suppressor.Gene=='No')
idx_match<- match(oncogene$Hugo.Symbol,human2mouse$human_name)
oncogene_mouse <- human2mouse$mouse_name[idx_match]
oncogene$mouse <- oncogene_mouse
oncogene <- oncogene[!is.na(oncogene$mouse),]

gdata$essential <- !is.na(match(gdata$id,essential_mouse))
gdata$essential_alt <- !is.na(match(gdata$id,alt_essential))
gdata$essential_combo <- gdata$essential | gdata$essential_alt

gdata$TSG <- !is.na(match(gdata$id,TSG$mouse))
gdata$oncogene <- !is.na(match(gdata$id,oncogene$mouse))

cutoff = -1
gdata <- gdata[order(gdata$Score,decreasing=F),]
gdata$rank <- 1: dim(gdata)[1]
gdata$color <- 'gray'
gdata$color[gdata$Score > 0 & gdata$FDR < 0.1] <- 'blue'
gdata$color[gdata$Score < 0 & gdata$FDR < 0.1] <- 'red'
gdata$alpha <- 0.1
gdata$alpha[gdata$color != 'gray'] = 1


p2 <- ggplot(gdata,aes(rank,Score)) + geom_point(aes(color=color,alpha=alpha)) + 
theme_classic() + scale_colour_manual(name="",values = c('gray'='gray','red'='red','blue'='blue'))

#cat(names(tail(sort(geneList),n=90)),sep='\n')

pdf(paste0('~/Documents/TSC_Paper/', 'CHOW_CRISPR_Screen_RankView.pdf'), width=5, height=2)
print(p2)
dev.off()

plot_frame <- data.frame(gene = subset(gdata,num>cutoff)$id, 
	stat = subset(gdata,num>cutoff)$Score,
	essential=subset(gdata,num>cutoff)$essential,
	essential_alt=subset(gdata,num>cutoff)$essential_alt,
	essential_combo=subset(gdata,num>cutoff)$essential | subset(gdata,num>cutoff)$essential_alt,
	TSG=subset(gdata,num>cutoff)$TSG,
  oncogene=subset(gdata,num>cutoff)$oncogene)

plot_frame <- plot_frame[order(plot_frame$stat),]

plot_frame$rank <- length(subset(gdata,num>cutoff)$Score):1

p1 <- ggplot(plot_frame, aes(x=rank,y=stat)) + 
geom_vline(xintercept= which(plot_frame$essential),size=0.1) + theme_void()

pdf(paste0('~/Documents/TSC_Paper/', 'CHOW_CRISPR_Screen_RankView_essential.pdf'), width=5, height=0.5)
print(p1)
dev.off()


#### CDAA
file1 = file.path("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa_vs_lib.gene_summary.txt")
file2 = file.path("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa_vs_lib.sgrna_summary.txt")
genes<-read.csv(file1,sep='\t')

gdata = ReadRRA(file1,score='rra')
gdata$num <- genes$num
sdata = ReadsgRRA(file2)
gdata$HumanGene = TransGeneID(gdata$id, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")
sdata$HumanGene = TransGeneID(sdata$Gene, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")

## Remove missing or duplicate human genes
idx = duplicated(gdata$HumanGene)|is.na(gdata$HumanGene)
gdata = gdata[!idx, ]
#depmap_similarity = ResembleDepmap(gdata, symbol = "HumanGene", score = "Score")
#head(depmap_similarity)

gdata2<-OmitCommonEssential(gdata, symbol = "HumanGene")
alt_essential <- setdiff(gdata$id,gdata2$id)

essential <- read.table('~/Documents/PinAPL/mageck/CountFiles/EssentialGenes.csv',header=FALSE,sep='\t')
human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')
idx_match<- match(essential$V1,human2mouse$human_name)
essential_mouse <- human2mouse$mouse_name[idx_match]
essential_mouse <- essential_mouse[!is.na(essential_mouse)]


#TSG <- read.table('~/Documents/PinAPL/mageck/CountFiles/TSG_Genes.csv',header=TRUE,sep=',')
TSG <- read.table('~/Documents/TSC_Paper/cancerGeneList.csv',header=TRUE,sep=',')
TSG <- subset(TSG,Is.Tumor.Suppressor.Gene=='Yes')
idx_match<- match(TSG$Hugo.Symbol,human2mouse$human_name)
TSG_mouse <- human2mouse$mouse_name[idx_match]
TSG$mouse <- TSG_mouse
TSG <- TSG[!is.na(TSG$mouse),]

oncogene <- read.table('~/Documents/TSC_Paper/cancerGeneList.csv',header=TRUE,sep=',')
oncogene <- subset(oncogene,Is.Tumor.Suppressor.Gene=='No')
idx_match<- match(oncogene$Hugo.Symbol,human2mouse$human_name)
oncogene_mouse <- human2mouse$mouse_name[idx_match]
oncogene$mouse <- oncogene_mouse
oncogene <- oncogene[!is.na(oncogene$mouse),]

gdata$essential <- !is.na(match(gdata$id,essential_mouse))
gdata$essential_alt <- !is.na(match(gdata$id,alt_essential))
gdata$essential_combo <- gdata$essential | gdata$essential_alt

gdata$TSG <- !is.na(match(gdata$id,TSG$mouse))
gdata$oncogene <- !is.na(match(gdata$id,oncogene$mouse))

cutoff = -1
gdata <- gdata[order(gdata$Score,decreasing=F),]
gdata$rank <- 1: dim(gdata)[1]
gdata$color <- 'gray'
gdata$color[gdata$Score > 0 & gdata$FDR < 0.1] <- 'blue'
gdata$color[gdata$Score < 0 & gdata$FDR < 0.1] <- 'red'
gdata$alpha <- 0.1
gdata$alpha[gdata$color != 'gray'] = 1


p2 <- ggplot(gdata,aes(rank,Score)) + geom_point(aes(color=color,alpha=alpha)) + 
theme_classic() + scale_colour_manual(name="",values = c('gray'='gray','red'='red','blue'='blue'))

#cat(names(tail(sort(geneList),n=90)),sep='\n')

pdf(paste0('~/Documents/TSC_Paper/', 'CDAA_CRISPR_Screen_RankView.pdf'), width=5, height=2)
print(p2)
dev.off()

plot_frame <- data.frame(gene = subset(gdata,num>cutoff)$id, 
  stat = subset(gdata,num>cutoff)$Score,
  essential=subset(gdata,num>cutoff)$essential,
  essential_alt=subset(gdata,num>cutoff)$essential_alt,
  essential_combo=subset(gdata,num>cutoff)$essential | subset(gdata,num>cutoff)$essential_alt,
  TSG=subset(gdata,num>cutoff)$TSG,
  oncogene=subset(gdata,num>cutoff)$oncogene)

plot_frame <- plot_frame[order(plot_frame$stat),]

plot_frame$rank <- length(subset(gdata,num>cutoff)$Score):1

p1 <- ggplot(plot_frame, aes(x=rank,y=stat)) + 
geom_vline(xintercept= which(plot_frame$essential),size=0.1) + theme_void()

pdf(paste0('~/Documents/TSC_Paper/', 'CDAACRISPR_Screen_RankView_essential.pdf'), width=5, height=0.5)
print(p1)
dev.off()

#######################################
#############  FIGURE 1C  #############
#######################################

chow1 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow1_vs_lib.gene_summary.txt",header=TRUE)
chow2 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow2_vs_lib.gene_summary.txt",header=TRUE)
chow3 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow3_vs_lib.gene_summary.txt",header=TRUE)
chow1 <- chow1[order(chow1$id),]
chow2 <- chow2[order(chow2$id),]
chow3 <- chow3[order(chow3$id),]
idx <- rowSums(data.frame(chow1$num,chow2$num,chow3$num) > 2) == 3
chow1 <- chow1[idx,]
chow2 <- chow2[idx,]
chow3 <- chow3[idx,]

a1 <- data.frame(chow1$neg.lfc,chow1$pos.lfc)
chow1$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(chow2$neg.lfc,chow2$pos.lfc)
chow2$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(chow3$neg.lfc,chow3$pos.lfc)
chow3$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])

chow1$goodsgrna <- apply(chow1,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
chow2$goodsgrna <- apply(chow2,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
chow3$goodsgrna <- apply(chow3,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))

scores <- data.frame(chow1$lfc,chow2$lfc,chow3$lfc)
rownames(scores) <- chow1$id
chow.log.fc <- robustbase::rowMedians(as.matrix(scores))
names(chow.log.fc) <- chow1$id
confidence <- data.frame(chow1$goodsgrna,chow2$goodsgrna,chow3$goodsgrna)
rownames(confidence) <- chow1$id

# Check for hogh logFC and discrepancy in direction
idy <- abs(rowSums(sign(scores))) == 3  | abs(chow.log.fc) < 1.5
scores <- scores[idy,]
chow.log.fc <- chow.log.fc[idy]
confidence <- confidence[idy,]

high_confidence <- rownames(confidence)[rowSums(confidence>=2) >= 2]


dis <- unlist(chow.log.fc)
p <- matrixTests::row_wilcoxon_twosample(scores,dis)
p.val <- p$pvalue
fdr <- p.adjust(p.val, method="BH")


library(EnhancedVolcano)


Knouse_TSGs <- read.csv('~/Documents/TSC_Paper/Knouse_Liver_TSGs.csv',header=FALSE)
Knouse_Broad_TSGs <- read.csv('~/Documents/TSC_Paper/Knouse_Broad_Liver_TSGs.csv',header=FALSE)
Knouse_depleted <- read.csv('~/Documents/TSC_Paper/Knouse_Liver_Depleted.csv',header=FALSE)
Tuson_TGs <- read.csv('~/Documents/TSC_Paper/Tuson_TGs.csv',header=TRUE)
Pan_Ca <- Tuson_TGs$PAN.Cancer_q.value
names(Pan_Ca) <- Tuson_TGs$Gene
TGs_TUSON <- names(Pan_Ca)[Pan_Ca<0.05]
h2m <- read.csv('~/Documents/TSC_Paper/human2mouse.csv',header=FALSE)
TGs_TUSON <- h2m$V2[match(TGs_TUSON, h2m$V1)]
TGs_TUSON <- TGs_TUSON[!is.na(TGs_TUSON)]
TOR_related <- c('Tsc1','Tsc2','Akt1','Akt2')
names(Pan_Ca) <- h2m$V2[match(names(Pan_Ca), h2m$V1)]
Pan_Ca <- Pan_Ca[!is.na(names(Pan_Ca))]

df <- data.frame(logFC = chow.log.fc, id = names(chow.log.fc), P = p.val, FDR = fdr)
df$Tuson <- Pan_Ca[match(rownames(df),names(Pan_Ca))]

keyvals <- ifelse(
    rownames(df) %in% unlist(Knouse_Broad_TSGs) , 'royalblue',
      ifelse(rownames(df) %in% unlist(Knouse_depleted), 'red',
        'gray'))
  keyvals[is.na(keyvals)] <- 'gray'
  names(keyvals)[keyvals == 'red'] <- 'high'
  names(keyvals)[keyvals == 'gray'] <- 'mid'
  names(keyvals)[keyvals == 'royalblue'] <- 'low'

alpha <- ifelse(
     rownames(df) %in% unlist(Knouse_TSGs) , 1,
      ifelse(rownames(df) %in% unlist(Knouse_depleted), 1,
        0.1))

labs <- rownames(df)
labs <- labs[alpha == 1]

EnhancedVolcano(df, lab = rownames(df),selectLab =labs ,x = 'logFC',
    y = 'FDR', pCutoff = 0.05,ylim = c(0,1),xlim=c(-8,3),
    colCustom = keyvals, colAlpha = alpha)

ggplot(subset(df,Tuson < 0.05),aes(logFC,-log10(Tuson))) + geom_point()

geneList <- chow.log.fc
p1 <- RankView(geneList,genelist = intersect(unlist(Knouse_Broad_TSGs),high_confidence),cutoff=c(quantile(dis,0.05),quantile(dis,0.95)),top = 0, bottom = 0,size=14,max.overlaps=Inf)



cdaa1 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa1_vs_lib.gene_summary.txt",header=TRUE)
cdaa2 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa2_vs_lib.gene_summary.txt",header=TRUE)
cdaa3 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa3_vs_lib.gene_summary.txt",header=TRUE)
cdaa1 <- cdaa1[order(cdaa1$id),]
cdaa2 <- cdaa2[order(cdaa2$id),]
cdaa3 <- cdaa3[order(cdaa3$id),]
idx <- rowSums(data.frame(cdaa1$num,cdaa2$num,cdaa3$num) > 2) == 3
cdaa1 <- cdaa1[idx,]
cdaa2 <- cdaa2[idx,]
cdaa3 <- cdaa3[idx,]

scores <- data.frame(cdaa1$neg.lfc,cdaa2$neg.lfc,cdaa3$neg.lfc)
rownames(scores) <- cdaa1$id
cdaa.log.fc <- robustbase::rowMedians(as.matrix(scores))
names(cdaa.log.fc) <- cdaa1$id

# Check for hogh logFC and discrepancy in direction
idy <- abs(rowSums(sign(scores))) == 3 | abs(cdaa.log.fc) < 2
scores <- scores[idy,]
cdaa.log.fc <- cdaa.log.fc[idy]

dis <- unlist(cdaa.log.fc)
p <- matrixTests::row_wilcoxon_twosample(scores,dis)
p.val <- p$pvalue
fdr <- p.adjust(p.val, method="BH")

geneList <- cdaa.log.fc
p1 <- RankView(geneList,cutoff=c(-1.3,1.3),top = 10, bottom = 10,size=14,max.overlaps=Inf)

df <- data.frame(logFC = cdaa.log.fc, id = names(cdaa.log.fc), P = p.val, FDR = fdr)
#df$FDR[combined.log.fc > 0] = pos.fdr[combined.log.fc > 0]

library(EnhancedVolcano)

Knouse_TSGs <- read.csv('~/Documents/TSC_Paper/Knouse_Liver_TSGs.csv',header=FALSE)
Knouse_depleted <- read.csv('~/Documents/TSC_Paper/Knouse_Liver_Depleted.csv',header=FALSE)
Tuson_TGs <- read.csv('~/Documents/TSC_Paper/Tuson_TGs.csv',header=TRUE)
Pan_Ca <- Tuson_TGs$PAN.Cancer_q.value
names(Pan_Ca) <- Tuson_TGs$Gene
TGs_TUSON <- names(Pan_Ca)[Pan_Ca<0.05]
h2m <- read.csv('~/Documents/TSC_Paper/human2mouse.csv',header=FALSE)
TGs_TUSON <- h2m$V2[match(TGs_TUSON, h2m$V1)]
TGs_TUSON <- TGs_TUSON[!is.na(TGs_TUSON)]
TOR_related <- c('Tsc1','Tsc2','Akt1','Akt2')

keyvals <- ifelse(
    rownames(df) %in% unlist(Knouse_TSGs) , 'royalblue',
      ifelse(rownames(df) %in% unlist(Knouse_depleted), 'red',
        'gray'))
  keyvals[is.na(keyvals)] <- 'gray'
  names(keyvals)[keyvals == 'red'] <- 'high'
  names(keyvals)[keyvals == 'gray'] <- 'mid'
  names(keyvals)[keyvals == 'royalblue'] <- 'low'

alpha <- ifelse(
     rownames(df) %in% unlist(Knouse_TSGs) , 1,
      ifelse(rownames(df) %in% unlist(Knouse_depleted), 1,
        0.1))

labs <- rownames(df)
labs[alpha != 1] = NA

EnhancedVolcano(df, lab = labs,x = 'logFC',
    y = 'FDR', pCutoff = 0.05,ylim = c(0,1),xlim=c(-8,3),
    colCustom = keyvals, colAlpha = alpha)



chow1 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow1_vs_lib.gene_summary.txt",header=TRUE)
chow2 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow2_vs_lib.gene_summary.txt",header=TRUE)
chow3 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow3_vs_lib.gene_summary.txt",header=TRUE)
chow1 <- chow1[order(chow1$id),]
chow2 <- chow2[order(chow2$id),]
chow3 <- chow3[order(chow3$id),]
cdaa1 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa1_vs_lib.gene_summary.txt",header=TRUE)
cdaa2 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa2_vs_lib.gene_summary.txt",header=TRUE)
cdaa3 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa3_vs_lib.gene_summary.txt",header=TRUE)
cdaa1 <- cdaa1[order(cdaa1$id),]
cdaa2 <- cdaa2[order(cdaa2$id),]
cdaa3 <- cdaa3[order(cdaa3$id),]
idx <- rowSums(data.frame(chow1$num,chow2$num,chow3$num,cdaa1$num,cdaa2$num,cdaa3$num) > 2) == 6
chow1 <- chow1[idx,]
chow2 <- chow2[idx,]
chow3 <- chow3[idx,]
cdaa1 <- cdaa1[idx,]
cdaa2 <- cdaa2[idx,]
cdaa3 <- cdaa3[idx,]


a1 <- data.frame(chow1$neg.lfc,chow1$pos.lfc)
chow1$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(chow2$neg.lfc,chow2$pos.lfc)
chow2$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(chow3$neg.lfc,chow3$pos.lfc)
chow3$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(cdaa1$neg.lfc,cdaa1$pos.lfc)
cdaa1$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(cdaa2$neg.lfc,cdaa2$pos.lfc)
cdaa2$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(cdaa3$neg.lfc,cdaa3$pos.lfc)
cdaa3$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])


chow1$goodsgrna <- apply(chow1,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
chow2$goodsgrna <- apply(chow2,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
chow3$goodsgrna <- apply(chow3,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
cdaa1$goodsgrna <- apply(cdaa1,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
cdaa2$goodsgrna <- apply(cdaa2,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
cdaa3$goodsgrna <- apply(cdaa3,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))


scores <- data.frame(chow1$lfc,chow2$lfc,chow3$lfc,cdaa1$lfc,cdaa2$lfc,cdaa3$lfc)
rownames(scores) <- cdaa1$id

confidence <- data.frame(chow1$goodsgrna,chow2$goodsgrna,chow3$goodsgrna,
  cdaa1$goodsgrna,cdaa2$goodsgrna,cdaa3$goodsgrna)
rownames(confidence) <- cdaa1$id

idx <- which(confidence == 0,arr.ind = TRUE)
mod_scores <- scores
mod_scores[cbind(idx[,1],idx[,2])] = NA

combined.log.fc <- robustbase::rowMedians(as.matrix(scores),na.rm = TRUE)
names(combined.log.fc) <- cdaa1$id

# Check for high logFC and discrepancy in direction
idy <- abs(rowSums(sign(scores))) >= 4 | abs(combined.log.fc) < 2
scores <- scores[idy,]
mod_scores <- mod_scores[idy,]
combined.log.fc <- combined.log.fc[idy]

dis <- unlist(combined.log.fc)
p <- matrixTests::row_wilcoxon_twosample(scores,dis)
p.val <- p$pvalue
fdr <- p.adjust(p.val, method="BH")

#pos.scores <- data.frame(chow1$pos.score,chow2$pos.score,chow3$pos.score,cdaa1$pos.score,cdaa2$pos.score,cdaa3$pos.score)
#dis <- unlist(pos.scores)
#pos.p <- matrixTests::row_wilcoxon_twosample(pos.scores,dis)
#pos.p.val <- pos.p$pvalue
#pos.fdr <- p.adjust(pos.p.val, method="BH")

geneList <- combined.log.fc
p1 <- RankView(geneList,cutoff=c(-1.3,1.3),top = 10, bottom = 10,size=14,max.overlaps=Inf)

df <- data.frame(logFC = combined.log.fc, id = names(combined.log.fc), P = p.val, FDR = fdr)
#df$FDR[combined.log.fc > 0] = pos.fdr[combined.log.fc > 0]

library(EnhancedVolcano)

Knouse_TSGs <- read.csv('~/Documents/TSC_Paper/Knouse_Liver_TSGs.csv',header=FALSE)
Knouse_depleted <- read.csv('~/Documents/TSC_Paper/Knouse_Liver_Depleted.csv',header=FALSE)
Tuson_TGs <- read.csv('~/Documents/TSC_Paper/Tuson_TGs.csv',header=TRUE)
Pan_Ca <- Tuson_TGs$PAN.Cancer_q.value
names(Pan_Ca) <- Tuson_TGs$Gene
TGs_TUSON <- names(Pan_Ca)[Pan_Ca<0.05]
h2m <- read.csv('~/Documents/TSC_Paper/human2mouse.csv',header=FALSE)
TGs_TUSON <- h2m$V2[match(TGs_TUSON, h2m$V1)]
TGs_TUSON <- TGs_TUSON[!is.na(TGs_TUSON)]
TOR_related <- c('Tsc1','Tsc2','Akt1','Akt2')
Knouse_hits <- read.csv('~/Documents/TSC_Paper/Knouse_Positive_Lib.csv',header=TRUE)
Knouse_hits_up <- subset(Knouse_hits,p.wilcox.bh < 0.05 & median.lfc.all>0)$id
Knouse_hits_down <- subset(Knouse_hits,p.wilcox.bh < 0.05 & median.lfc.all<0)$id


keyvals <- ifelse(
    rownames(df) %in% unlist(Knouse_hits_up) , 'royalblue',
      ifelse(rownames(df) %in% unlist(Knouse_hits_down), 'red',
        'gray'))
  keyvals[is.na(keyvals)] <- 'gray'
  names(keyvals)[keyvals == 'red'] <- 'low'
  names(keyvals)[keyvals == 'gray'] <- 'mid'
  names(keyvals)[keyvals == 'royalblue'] <- 'high'

alpha <- ifelse(
     rownames(df) %in% unlist(Knouse_hits_up) , 1,
      ifelse(rownames(df) %in% unlist(Knouse_hits_down), 1,
        0.1))

labs <- rownames(df)
labs[keyvals != 'royalblue'] = NA

p1 <- EnhancedVolcano(df, lab = labs,x = 'logFC',
    y = 'FDR', pCutoff = 0.1,FCcutoff = 0.25, ylim = c(0,2),xlim=c(-8,5),
    colCustom = keyvals, colAlpha = alpha,drawConnectors = TRUE)

pdf(paste0('~/Documents/TSC_Paper/', 'Merged_CRISPR_Screen_Volcano.pdf'), width=6, height=6)
p1
dev.off()


#######################################
#############  FIGURE 1D  #############
#######################################
chow1 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow1_vs_lib.gene_summary.txt",header=TRUE)
chow2 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow2_vs_lib.gene_summary.txt",header=TRUE)
chow3 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow3_vs_lib.gene_summary.txt",header=TRUE)
chow1 <- chow1[order(chow1$id),]
chow2 <- chow2[order(chow2$id),]
chow3 <- chow3[order(chow3$id),]
cdaa1 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa1_vs_lib.gene_summary.txt",header=TRUE)
cdaa2 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa2_vs_lib.gene_summary.txt",header=TRUE)
cdaa3 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa3_vs_lib.gene_summary.txt",header=TRUE)
cdaa1 <- cdaa1[order(cdaa1$id),]
cdaa2 <- cdaa2[order(cdaa2$id),]
cdaa3 <- cdaa3[order(cdaa3$id),]
idx <- rowSums(data.frame(chow1$num,chow2$num,chow3$num,cdaa1$num,cdaa2$num,cdaa3$num) > 2) == 6
chow1 <- chow1[idx,]
chow2 <- chow2[idx,]
chow3 <- chow3[idx,]
cdaa1 <- cdaa1[idx,]
cdaa2 <- cdaa2[idx,]
cdaa3 <- cdaa3[idx,]


a1 <- data.frame(chow1$neg.lfc,chow1$pos.lfc)
chow1$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(chow2$neg.lfc,chow2$pos.lfc)
chow2$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(chow3$neg.lfc,chow3$pos.lfc)
chow3$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(cdaa1$neg.lfc,cdaa1$pos.lfc)
cdaa1$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(cdaa2$neg.lfc,cdaa2$pos.lfc)
cdaa2$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(cdaa3$neg.lfc,cdaa3$pos.lfc)
cdaa3$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])


chow1$goodsgrna <- apply(chow1,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
chow2$goodsgrna <- apply(chow2,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
chow3$goodsgrna <- apply(chow3,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
cdaa1$goodsgrna <- apply(cdaa1,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
cdaa2$goodsgrna <- apply(cdaa2,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
cdaa3$goodsgrna <- apply(cdaa3,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))


scores <- data.frame(chow1$lfc,chow2$lfc,chow3$lfc,cdaa1$lfc,cdaa2$lfc,cdaa3$lfc)
rownames(scores) <- cdaa1$id

confidence <- data.frame(chow1$goodsgrna,chow2$goodsgrna,chow3$goodsgrna,
  cdaa1$goodsgrna,cdaa2$goodsgrna,cdaa3$goodsgrna)
rownames(confidence) <- cdaa1$id

idx <- which(confidence == 0,arr.ind = TRUE)
mod_scores <- scores
mod_scores[cbind(idx[,1],idx[,2])] = NA

combined.log.fc <- robustbase::rowMedians(as.matrix(scores),na.rm = TRUE)
names(combined.log.fc) <- cdaa1$id

# Check for high logFC and discrepancy in direction
idy <- (abs(rowSums(sign(scores))) >=4 | abs(combined.log.fc) < 2) & (rowSums(confidence >= 2) >= 3)
#idy <- abs(rowSums(sign(scores))) == 6 
scores <- scores[idy,]
mod_scores <- mod_scores[idy,]
combined.log.fc <- combined.log.fc[idy]

Tuson_TGs <- read.csv('~/Documents/TSC_Paper/Tuson_TGs.csv',header=TRUE)
Pan_Ca <- Tuson_TGs$PAN.Cancer_q.value
names(Pan_Ca) <- Tuson_TGs$Gene
TGs_TUSON <- names(Pan_Ca)[Pan_Ca<0.1]
h2m <- read.csv('~/Documents/TSC_Paper/human2mouse.csv',header=FALSE)
TGs_TUSON <- h2m$V2[match(TGs_TUSON, h2m$V1)]
TGs_TUSON <- TGs_TUSON[!is.na(TGs_TUSON)]

TSG <- read.table('~/Documents/TSC_Paper/cancerGeneList.csv',header=TRUE,sep=',')
TSG <- subset(TSG,Is.Tumor.Suppressor.Gene=='Yes')
idx_match<- match(TSG$Hugo.Symbol,human2mouse$human_name)
TSG_mouse <- human2mouse$mouse_name[idx_match]
TSG$mouse <- TSG_mouse
TSG <- TSG[!is.na(TSG$mouse),]

essential <- read.table('~/Documents/PinAPL/mageck/CountFiles/EssentialGenes.csv',header=FALSE,sep='\t')
human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')
idx_match<- match(essential$V1,human2mouse$human_name)
essential_mouse <- human2mouse$mouse_name[idx_match]
essential_mouse <- essential_mouse[!is.na(essential_mouse)]

oncogene <- read.table('~/Documents/TSC_Paper/cancerGeneList.csv',header=TRUE,sep=',')
oncogene <- subset(oncogene,Is.Tumor.Suppressor.Gene=='No')
idx_match<- match(oncogene$Hugo.Symbol,human2mouse$human_name)
oncogene_mouse <- human2mouse$mouse_name[idx_match]
oncogene$mouse <- oncogene_mouse
oncogene <- oncogene[!is.na(oncogene$mouse),]

df <- data.frame(gene = names(combined.log.fc), logFC = combined.log.fc)
df$class <- ifelse(
  df$gene %in% TSG$mouse, 'TSG',
  ifelse(
    df$gene %in% essential_mouse,'essential',
      'other'
    )
  )

pdf(paste0('~/Documents/TSC_Paper/', 'Merged_CRISPR_Screen_Volcano.pdf'), width=4, height=1.5)
ggplot(df,aes(logFC, fill = class)) + geom_density(alpha=.1,bw=1) + 
scale_fill_manual(values=c('TSG' = 'red','other' = 'gray','essential' = 'blue')) + theme_classic()
dev.off()

#######################################
#############  FIGURE 1E  #############
#######################################

df_sub <- subset(df,class == 'TSG')

bot <- quantile(df_sub$logFC,1/3)
top <- quantile(df_sub$logFC,2/3)

df_sub$sub_class[df_sub$logFC < bot] <- 'Down'
df_sub$sub_class[df_sub$logFC > top] = 'Up'

pdf(paste0('~/Documents/TSC_Paper/', 'Merged_CRISPR_Screen_Volcano_Up_vs_Down.pdf'), width=4, height=1.5)
ggplot(df_sub,aes(logFC, fill = sub_class)) + geom_density(alpha=.1,bw=0.5) + 
scale_fill_manual(values=c('red','blue')) + theme_classic()
dev.off()

library(enrichR)
library(forcats)
library(DOSE)


dbs <- c("ClinVar_2019")

#enriched <- enrichr(gset, dbs, background = background, include_overlap = TRUE)
enriched <- enrichr(df_sub$gene[df_sub$logFC < bot], dbs)
enriched_alt <- enrichr(df_sub$gene[df_sub$logFC > top], dbs)


pdf(paste0('~/Documents/TSC_Paper/', 'TSG_Up_vs_Down.pdf'), width=7, height=6)
p1<- ggplot(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + 
geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + 
scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Combined.Score,decreasing=T),][rev(1:5),]$Term,"_"), `[`, 1))) + 
theme(axis.text=element_text(colour="black"))

p2<- ggplot(enriched_alt[[1]][order(enriched_alt[[1]]$Combined.Score,decreasing=T),][rev(1:5),], 
  (aes(x=Combined.Score, y=fct_inorder(Term), color = as.numeric(Adjusted.P.value), size=parse_ratio(Overlap)))) + 
geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + 
scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched_alt[[1]][order(enriched_alt[[1]]$Combined.Score,decreasing=T),][rev(1:5),]$Term,"_"), `[`, 1))) + 
theme(axis.text=element_text(colour="black"))

p1/p2
dev.off()




cat(subset(df,class == 'TSG' & logFC > top)$gene,sep='\n')
#cat(subset(df,class == 'TSG' & logFC < -0.5)$gene,sep='\n')

#######################################
#############  FIGURE 1F  #############
#######################################

chow1 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow1_vs_lib.gene_summary.txt",header=TRUE)
chow2 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow2_vs_lib.gene_summary.txt",header=TRUE)
chow3 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow3_vs_lib.gene_summary.txt",header=TRUE)
chow1 <- chow1[order(chow1$id),]
chow2 <- chow2[order(chow2$id),]
chow3 <- chow3[order(chow3$id),]
cdaa1 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa1_vs_lib.gene_summary.txt",header=TRUE)
cdaa2 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa2_vs_lib.gene_summary.txt",header=TRUE)
cdaa3 <- read.table("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa3_vs_lib.gene_summary.txt",header=TRUE)
cdaa1 <- cdaa1[order(cdaa1$id),]
cdaa2 <- cdaa2[order(cdaa2$id),]
cdaa3 <- cdaa3[order(cdaa3$id),]
idx <- rowSums(data.frame(chow1$num,chow2$num,chow3$num,cdaa1$num,cdaa2$num,cdaa3$num) > 2) == 6
chow1 <- chow1[idx,]
chow2 <- chow2[idx,]
chow3 <- chow3[idx,]
cdaa1 <- cdaa1[idx,]
cdaa2 <- cdaa2[idx,]
cdaa3 <- cdaa3[idx,]


a1 <- data.frame(chow1$neg.lfc,chow1$pos.lfc)
chow1$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(chow2$neg.lfc,chow2$pos.lfc)
chow2$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(chow3$neg.lfc,chow3$pos.lfc)
chow3$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(cdaa1$neg.lfc,cdaa1$pos.lfc)
cdaa1$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(cdaa2$neg.lfc,cdaa2$pos.lfc)
cdaa2$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])
a1 <- data.frame(cdaa3$neg.lfc,cdaa3$pos.lfc)
cdaa3$lfc <- apply(a1, 1, function(x) x[1 + as.numeric(abs(x[2]) > abs(x[1]))])


chow1$goodsgrna <- apply(chow1,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
chow2$goodsgrna <- apply(chow2,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
chow3$goodsgrna <- apply(chow3,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
cdaa1$goodsgrna <- apply(cdaa1,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
cdaa2$goodsgrna <- apply(cdaa2,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))
cdaa3$goodsgrna <- apply(cdaa3,1,function(x) as.numeric(ifelse(as.numeric(x[15]) < 0,x[7],x[13])))


scores <- data.frame(chow1$lfc,chow2$lfc,chow3$lfc,cdaa1$lfc,cdaa2$lfc,cdaa3$lfc)
rownames(scores) <- cdaa1$id

confidence <- data.frame(chow1$goodsgrna,chow2$goodsgrna,chow3$goodsgrna,
  cdaa1$goodsgrna,cdaa2$goodsgrna,cdaa3$goodsgrna)
rownames(confidence) <- cdaa1$id

idx <- which(confidence == 0,arr.ind = TRUE)
mod_scores <- scores
mod_scores[cbind(idx[,1],idx[,2])] = NA

combined.log.fc <- robustbase::rowMedians(as.matrix(scores),na.rm = TRUE)
names(combined.log.fc) <- cdaa1$id

# Check for high logFC and discrepancy in direction
idy <- (abs(rowSums(sign(scores))) >=4 | abs(combined.log.fc) < 2) & (rowSums(confidence >= 2) >= 3)
#idy <- abs(rowSums(sign(scores))) == 6 
scores <- scores[idy,]
mod_scores <- mod_scores[idy,]
combined.log.fc <- combined.log.fc[idy]

Tuson_TGs <- read.csv('~/Documents/TSC_Paper/Tuson_TGs.csv',header=TRUE)
Pan_Ca <- Tuson_TGs$PAN.Cancer_q.value
names(Pan_Ca) <- Tuson_TGs$Gene
TGs_TUSON <- names(Pan_Ca)[Pan_Ca<0.1]
h2m <- read.csv('~/Documents/TSC_Paper/human2mouse.csv',header=FALSE)
TGs_TUSON <- h2m$V2[match(TGs_TUSON, h2m$V1)]
TGs_TUSON <- TGs_TUSON[!is.na(TGs_TUSON)]

TSG <- read.table('~/Documents/TSC_Paper/cancerGeneList.csv',header=TRUE,sep=',')
TSG <- subset(TSG,Is.Tumor.Suppressor.Gene=='Yes')
idx_match<- match(TSG$Hugo.Symbol,human2mouse$human_name)
TSG_mouse <- human2mouse$mouse_name[idx_match]
TSG$mouse <- TSG_mouse
TSG <- TSG[!is.na(TSG$mouse),]

essential <- read.table('~/Documents/PinAPL/mageck/CountFiles/EssentialGenes.csv',header=FALSE,sep='\t')
human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')
idx_match<- match(essential$V1,human2mouse$human_name)
essential_mouse <- human2mouse$mouse_name[idx_match]
essential_mouse <- essential_mouse[!is.na(essential_mouse)]

oncogene <- read.table('~/Documents/TSC_Paper/cancerGeneList.csv',header=TRUE,sep=',')
oncogene <- subset(oncogene,Is.Tumor.Suppressor.Gene=='No')
idx_match<- match(oncogene$Hugo.Symbol,human2mouse$human_name)
oncogene_mouse <- human2mouse$mouse_name[idx_match]
oncogene$mouse <- oncogene_mouse
oncogene <- oncogene[!is.na(oncogene$mouse),]

df <- data.frame(gene = names(combined.log.fc), logFC = combined.log.fc)
df$class <- ifelse(
  df$gene %in% TSG$mouse, 'TSG',
  ifelse(
    df$gene %in% essential_mouse,'essential',
      'other'
    )
  )


df$HumanGene = TransGeneID(df$gene, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")
## Remove missing or duplicate human genes
idx2 = duplicated(df$HumanGene)|is.na(df$HumanGene)
df = df[!idx2, ]



###
depmapDat = LoadDepmap()

Depmap = depmapDat$Depmap
sampleinfo = depmapDat$sampleinfo

#"bile duct cancer" , 'liver cancer'
sampleinfo$primary_disease = tolower(sampleinfo$primary_disease)
sampleinfo$subtype_disease = tolower(sampleinfo$subtype_disease)
sampleinfo$cell_line = tolower(sampleinfo$cell_line)
sampleinfo = sampleinfo[colnames(Depmap), ]
#idx = which(sampleinfo$primary_disease == 'liver cancer' & sampleinfo$subtype_disease == "hepatocellular carcinoma")

genes = intersect(df$HumanGene, rownames(Depmap))

dd = df[!duplicated(df[, 'HumanGene']), ]
rownames(dd) = dd[, 'HumanGene']


similarity = apply(Depmap[genes, ], 2, function(x) {
            tmp = cor.test(x, dd[genes, 'logFC'], method = 'pearson', 
                na.action = na.omit)
            c(tmp$estimate, tmp$p.value)
})

similarity = as.data.frame(t(similarity))
colnames(similarity) = c("estimate", "p.value")
#rownames(similarity) = sampleinfo[colnames(Depmap), 1]
similarity = similarity[order(-similarity$estimate), ]

###


library(matrixStats)
depmap_avg <- rowMedians(as.matrix(Depmap[genes,idx]))
#depmap_avg <- -log10(max(depmap_avg+1) - depmap_avg)*10
#depmap_avg <- depmap_avg - median(depmap_avg)
gdata_cc <- data.frame(Gene=genes,screen=dd[genes, 'logFC'],depmap=depmap_avg+1)


gdata_norm <- gdata_cc
gdata_norm$Gene <- c()
#gdata_norm <- as.data.frame(qn(gdata_norm))
#gdata_norm$Gene <- rownames(gdata_norm)


gdata_norm_new<-normalize.quantiles(as.matrix(gdata_norm),copy=TRUE)

rownames(gdata_norm_new) <- rownames(gdata_norm)
colnames(gdata_norm_new) <- colnames(gdata_norm)
gdata_norm_new <- as.data.frame(gdata_norm_new)
gdata_norm_new$Gene <- rownames(gdata_norm)
gdata_norm_new$label <- NA

cut_1 = -0.2998331
cut_2 = 0.8146948

gdata_norm_new$label[gdata_norm_new$screen < cut_1 & gdata_norm_new$depmap > cut_2] = gdata_norm_new$Gene[gdata_norm_new$screen < cut_1 & gdata_norm_new$depmap > cut_2]
gdata_norm_new$label[gdata_norm_new$screen > cut_2 & gdata_norm_new$depmap < cut_1] = gdata_norm_new$Gene[gdata_norm_new$screen > cut_2 & gdata_norm_new$depmap < cut_1]


p1 <- ggplot(gdata_norm_new,aes(screen,depmap)) + geom_point(aes(size = abs(screen - depmap)),alpha=0.2) + 
scale_size(range = c(0,3)) + 
geom_hline(yintercept=cut_1) + geom_hline(yintercept=cut_2) + 
geom_vline(xintercept=cut_1) + geom_vline(xintercept=cut_2) +
geom_text_repel(aes(label = label),max.overlaps = 50) + 
theme_classic()




pdf(paste0('~/Documents/TSC_Paper/', 'CRISPR_Screen_Depmap.pdf'), width=9, height=4)

print(p1)

dev.off()


#######################################
#############  FIGURE 1G  #############
#######################################

file1 = file.path("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow_vs_lib.gene_summary.txt")
file2 = file.path("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa_vs_lib.gene_summary.txt")
file3 = file.path("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa_vs_chow_clean.gene_summary.txt")
file4 = file.path("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa_vs_chow.gene_summary.txt")


gdata1 = ReadRRA(file1,score='lfc')
gdata2 = ReadRRA(file2,score='lfc')
gdata4 = ReadRRA(file3,score='lfc')
gdata3 = ReadRRA(file4,score='lfc')

cat(setdiff(gdata4$id,gdata3$id),sep='\n')
cat(intersect(subset(gdata2, Score < -5)$id,subset(gdata1, Score > -2)$id),sep='\n')


chow <- read.csv(file1,sep='\t')
cdaa <- read.csv(file2,sep='\t')
both <- read.csv(file4,sep='\t')

essential <- read.table('~/Documents/PinAPL/mageck/CountFiles/EssentialGenes.csv',header=FALSE,sep='\t')
human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')
idx_match<- match(essential$V1,human2mouse$human_name)
essential_mouse <- human2mouse$mouse_name[idx_match]
essential_mouse <- essential_mouse[!is.na(essential_mouse)]


clear <- read.table('~/Documents/TSC_Paper/CLEAR_genes.csv',header=FALSE,sep='\t')$V1


gdata1$Score[match(gdata3$id,gdata1$id)]

df <- data.frame(gene = gdata3$id, lfc = gdata3$Score, 
  lfc_chow = gdata1$Score[match(gdata3$id,gdata1$id)],
  lfc_cdaa = gdata2$Score[match(gdata3$id,gdata2$id)],
  FDR = gdata3$FDR,
  essential = !(gdata3$id %in% essential_mouse),
  clear = !(gdata3$id %in% clear),
  color = ifelse(gdata3$FDR < 0.1,'blue','gray'),
  color_alt = !(gdata3$id %in% intersect(subset(gdata1,FDR < 0.1 & Score>0)$id,subset(gdata2,FDR < 0.1 & Score>0)$id)),
  conf = both$num[match(gdata3$id,both$id)]
)
df$color[is.na(df$color)] = 'gray'
df$alpha = 0.2
df$alpha[!df$color_alt] = 1

df <- subset(df,conf >=4)

library(ggrepel)



p1 <- ggplot(df,aes(lfc_cdaa,lfc_chow)) + 
  geom_point(aes(size=10.385-abs(lfc_chow-lfc_cdaa),color=color_alt),alpha=df$alpha) +
  geom_abline(intercept = 0, slope = 1, size = 0.5) + 
  scale_size_continuous(range = c(0,2)) + 
  geom_text_repel(subset(df,!color_alt),max.overlaps = Inf, mapping=aes(lfc_cdaa,lfc_chow,label= gene)) +
  scale_colour_manual(values=c("blue", "gray")) + 
  theme_classic() + xlim(-5,5) + ylim(-5,5) + 
  geom_smooth(method='lm', se = FALSE)

pdf(paste0('~/Documents/TSC_Paper/', 'CDAA_vs_Chow_dot.pdf'), width=10, height=4)
p1
dev.off()





#######################################
#############  FIGURE 1H #############
#######################################

file1 = file.path("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/chow_vs_lib.gene_summary.txt")
file2 = file.path("/Users/ikuz/Documents/TSC_Paper/InVivoCRISPR/FinalAnalysis/cdaa_vs_lib.gene_summary.txt")

gdata1 = ReadRRA(file1,score='lfc')
gdata2 = ReadRRA(file2,score='lfc')

knouse <- read.csv('~/Documents/TSC_Paper/Knouse_Positive_Lib.csv')

a <- intersect(subset(gdata1,FDR < 0.1 & Score>0)$id,subset(gdata2,FDR < 0.1 & Score>0)$id)
b <- intersect(subset(gdata1,FDR < 0.1 & Score>0)$id,subset(knouse,p.wilcox.bh < 0.1 & median.lfc.all>0)$id)
c <- intersect(subset(gdata2,FDR < 0.1 & Score>0)$id,subset(knouse,p.wilcox.bh < 0.1 & median.lfc.all>0)$id)
unique(c(a,b,c))



#######################################
#############  FIGURE 1I #############
#######################################

human_muts <- read.csv('~/Documents/TSC_Paper/CLCA_muts.txt',sep='\t')
human_muts$Freq <- as.numeric(gsub("[^0-9.-]", "", human_muts$Freq))
human_muts <- subset(human_muts,Is.Cancer.Gene..source..OncoKB. == 'Yes' & Freq > 1)
#317 genes
human_muts$mouse <- human2mouse$mouse_name[match(human_muts$Gene,human2mouse$human_name)]

intersect(unique(c(a,b,c)),human_muts$mouse)

#Tuson_TGs$mouse <- human2mouse$mouse_name[match(Tuson_TGs$Gene,human2mouse$human_name)]
#Tuson_TGs[match(a,Tuson_TGs$mouse),]

a<-intersect(subset(gdata2,FDR < 0.1 & Score>0)$id,human_muts$mouse)
b<-intersect(subset(gdata1,FDR < 0.1 & Score>0)$id,human_muts$mouse)
c<-intersect(subset(knouse,p.wilcox.bh < 0.1 & median.lfc.all>0)$id,human_muts$mouse)
unique(c(a,b,c))
human_muts[human_muts$mouse %in% unique(c(a,b,c)),]



a<-intersect(subset(gdata2,Score>0)$id,human_muts$mouse)
b<-intersect(subset(gdata1,Score>0)$id,human_muts$mouse)
c<-intersect(subset(knouse,median.lfc.all>0)$id,human_muts$mouse)
unique(c(a,b,c))
human_muts[human_muts$mouse %in% unique(c(a,b,c)),]

#223 vs 317





file1 = file.path("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/chow_vs_lib_keepzero.gene_summary.txt")
file2 = file.path("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/chow_vs_lib_keepzero.sgrna_summary.txt")

all_reads = read.csv("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/FinalCounts.csv")

background <- unique(subset(all_reads,Base1 + Base2  + Base3 > 15)$Gene)

gdata = ReadRRA(file1,score='rra')
sdata = ReadsgRRA(file2)
gdata$HumanGene = TransGeneID(gdata$id, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")
sdata$HumanGene = TransGeneID(sdata$Gene, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")

gdata_s <- gdata %>% arrange(desc(Score))
sdata_s <- sdata %>% arrange(desc(LFC))

sgRNA_sig <- subset(sdata_s,FDR < 0.05 & LFC > 0)

counts_sg <- table(sgRNA_sig$Gene)
counts_sg <- counts_sg[names(counts_sg) %in% gdata_s$id]


gdata_s$sgCount <- 0


idx<-match(gdata_s$id,names(counts_sg))
idx <- idx[!is.na(idx)]
gdata_s[idx,]$sgCount <- counts_sg



gdata_filt <- subset(gdata_s,sgCount>1)


gset <- gdata_filt$id[1:125]

gset_tail <- tail(gdata_filt$id,n=125)

#Reactome enrichr

library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)

dbs <- c("ChEA_2022","Reactome_Pathways_2024","GO_Biological_Process_2023")


#enriched <- enrichr(gset, dbs, background = background, include_overlap = TRUE)
enriched <- enrichr(gset, dbs)

pdf(paste0('~/Documents/TSC_Paper/', 'CRISPR_Screen_Enrich.pdf'), width=7, height=4)
p1<- ggplot(enriched[[2]][order(enriched[[2]]$Adjusted.P.value,decreasing=F),][rev(1:5),], 
  (aes(x=Adjusted.P.value, y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + 
geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + 
ggtitle('Reactome 2024') + 
scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[2]][order(enriched[[2]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"_"), `[`, 1))) + 
theme(axis.text=element_text(colour="black"))
p1
dev.off()



#######################################
#############  FIGURE 1F  #############
#######################################





#######################################
#############  FIGURE 1G  #############
#######################################

file1 = file.path("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/cdaa_vs_chow_keepzero.gene_summary.txt")
file2 = file.path("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/cdaa_vs_chow_keepzero.sgrna_summary.txt")


genes<-read.csv(file1,sep='\t')

gdata = ReadRRA(file1,score='rra')
gdata$num <- genes$num
sdata = ReadsgRRA(file2)
gdata$HumanGene = TransGeneID(gdata$id, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")
sdata$HumanGene = TransGeneID(sdata$Gene, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")

## Remove missing or duplicate human genes
idx = duplicated(gdata$HumanGene)|is.na(gdata$HumanGene)
gdata = gdata[!idx, ]
#depmap_similarity = ResembleDepmap(gdata, symbol = "HumanGene", score = "Score")
#head(depmap_similarity)

gdata2<-OmitCommonEssential(gdata, symbol = "HumanGene")
alt_essential <- setdiff(gdata$id,gdata2$id)


TSG <- read.table('~/Documents/TSC_Paper/cancerGeneList.csv',header=TRUE,sep=',')
TSG <- subset(TSG,Is.Tumor.Suppressor.Gene=='Yes')
idx_match<- match(TSG$Hugo.Symbol,human2mouse$human_name)
TSG_mouse <- human2mouse$mouse_name[idx_match]
TSG$mouse <- TSG_mouse
TSG <- TSG[!is.na(TSG$mouse),]

oncogene <- read.table('~/Documents/TSC_Paper/cancerGeneList.csv',header=TRUE,sep=',')
oncogene <- subset(oncogene,Is.Tumor.Suppressor.Gene=='No')
idx_match<- match(oncogene$Hugo.Symbol,human2mouse$human_name)
oncogene_mouse <- human2mouse$mouse_name[idx_match]
oncogene$mouse <- oncogene_mouse
oncogene <- oncogene[!is.na(oncogene$mouse),]

gdata$essential <- !is.na(match(gdata$id,essential_mouse))
gdata$essential_alt <- !is.na(match(gdata$id,alt_essential))
gdata$essential_combo <- gdata$essential | gdata$essential_alt

gdata$TSG <- !is.na(match(gdata$id,TSG$mouse))
gdata$oncogene <- !is.na(match(gdata$id,oncogene$mouse))

geneList= gdata$Score
names(geneList) = gdata$id
p2 <- RankView(geneList,cutoff=c(-1.3,1.3),top = 10, bottom = 10,size=14,max.overlaps=Inf)
p2 <- p2 + theme(text = element_text(size = 14))

#cat(names(tail(sort(geneList),n=90)),sep='\n')

pdf(paste0('~/Documents/TSC_Paper/', 'CDAA_CRISPR_Screen_RankView.pdf'), width=5, height=5)
print(p2)
dev.off()

gene_of_int = c('Tsc1','Tsc2','Rraga','Rragb','Rragc','Rragd','Flcn',
  'Lamtor1','Lamtor2','Lamtor3','Lamtor4','Lamtor5','Tfe3','Tfeb',
  'Mitf','Mtor',"Rptor",'Deptor','Rictor',"Fnip1","Fnip2","Akt1","Akt2",'Akt3',
  'Rps6kb1','Eif4ebp1','Rheb','Akt1s1','Mlst8','Pik3ca','Pik3r1','Pten')
p2 <- RankView(geneList,cutoff=c(-1.3,1.3),top = 0, bottom = 0,size=14,max.overlaps=Inf,genelist=gene_of_int)
p2 <- p2 + theme(text = element_text(size = 14))

#######################################
#############  FIGURE 1G  #############
#######################################


file1 = file.path("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/cdaa_vs_chow_keepzero.gene_summary.txt")
file2 = file.path("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/cdaa_vs_chow_keepzero.sgrna_summary.txt")

all_reads = read.csv("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/FinalCounts.csv")

background <- unique(subset(all_reads,Base1 + Base2  + Base3 > 15)$Gene)

gdata = ReadRRA(file1,score='rra')
sdata = ReadsgRRA(file2)
gdata$HumanGene = TransGeneID(gdata$id, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")
sdata$HumanGene = TransGeneID(sdata$Gene, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")

gdata_s <- gdata %>% arrange(desc(Score))
sdata_s <- sdata %>% arrange(desc(LFC))

sgRNA_sig <- subset(sdata_s,FDR < 0.05 & LFC > 0)

counts_sg <- table(sgRNA_sig$Gene)
counts_sg <- counts_sg[names(counts_sg) %in% gdata_s$id]


gdata_s$sgCount <- 0


idx<-match(gdata_s$id,names(counts_sg))
idx <- idx[!is.na(idx)]
gdata_s[idx,]$sgCount <- counts_sg



gdata_filt <- subset(gdata_s,sgCount>0)


gset <- gdata_filt$id[1:125]


library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)

dbs <- c("ChEA_2022","Reactome_Pathways_2024","GO_Biological_Process_2023")


#enriched <- enrichr(gset, dbs, background = background, include_overlap = TRUE)
enriched <- enrichr(gset, dbs)

pdf(paste0('~/Documents/TSC_Paper/', 'CDAA_CRISPR_Screen_Enrich.pdf'), width=7, height=4)
p1<- ggplot(enriched[[2]][order(enriched[[2]]$Adjusted.P.value,decreasing=F),][rev(1:5),], 
  (aes(x=Adjusted.P.value, y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + 
geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + 
ggtitle('Reactome 2024') + 
scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[2]][order(enriched[[2]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"_"), `[`, 1))) + 
theme(axis.text=element_text(colour="black"))
p1
dev.off()



#######################################
#############  FIGURE S1C #############
#######################################

file1 = file.path("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/chow_vs_lib_keepzero.gene_summary.txt")
file2 = file.path("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/chow_vs_lib_keepzero.sgrna_summary.txt")

all_reads = read.csv("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/FinalCounts.csv")

background <- unique(subset(all_reads,Base1 + Base2  + Base3 > 15)$Gene)

gdata = ReadRRA(file1,score='rra')
sdata = ReadsgRRA(file2)
gdata$HumanGene = TransGeneID(gdata$id, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")
sdata$HumanGene = TransGeneID(sdata$Gene, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")

gdata_s <- gdata %>% arrange((Score))
sdata_s <- sdata %>% arrange((LFC))

sgRNA_sig <- subset(sdata_s,FDR < 0.05 & LFC > 0)

counts_sg <- table(sgRNA_sig$Gene)
counts_sg <- counts_sg[names(counts_sg) %in% gdata_s$id]


gdata_s$sgCount <- 0


idx<-match(gdata_s$id,names(counts_sg))
idx <- idx[!is.na(idx)]
gdata_s[idx,]$sgCount <- counts_sg



gdata_filt <- subset(gdata_s,sgCount>1)


gset <- gdata_filt$id[1:250]


library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)

dbs <- c("ChEA_2022","Reactome_Pathways_2024","GO_Biological_Process_2023")


#enriched <- enrichr(gset, dbs, background = background, include_overlap = TRUE)
enriched <- enrichr(gset, dbs)

pdf(paste0('~/Documents/TSC_Paper/', 'CRISPR_Screen_Deplete.pdf'), width=7, height=4)
p1<- ggplot(enriched[[3]][order(enriched[[2]]$Adjusted.P.value,decreasing=F),][rev(1:5),], 
  (aes(x=Adjusted.P.value, y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + 
geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + 
ggtitle('Reactome 2024') + 
scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[3]][order(enriched[[3]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"_"), `[`, 1))) + 
theme(axis.text=element_text(colour="black"))
p1
dev.off()


#######################################
#############  FIGURE S1D #############
#######################################

file1 = file.path("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/cdaa_vs_chow_keepzero.gene_summary.txt")
file2 = file.path("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/cdaa_vs_chow_keepzero.sgrna_summary.txt")

all_reads = read.csv("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/FinalCounts.csv")

background <- unique(subset(all_reads,Base1 + Base2  + Base3 > 15)$Gene)

gdata = ReadRRA(file1,score='rra')
sdata = ReadsgRRA(file2)
gdata$HumanGene = TransGeneID(gdata$id, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")
sdata$HumanGene = TransGeneID(sdata$Gene, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")

gdata_s <- gdata %>% arrange((Score))
sdata_s <- sdata %>% arrange((LFC))

sgRNA_sig <- subset(sdata_s,FDR < 0.05 & LFC > 0)

counts_sg <- table(sgRNA_sig$Gene)
counts_sg <- counts_sg[names(counts_sg) %in% gdata_s$id]


gdata_s$sgCount <- 0


idx<-match(gdata_s$id,names(counts_sg))
idx <- idx[!is.na(idx)]
gdata_s[idx,]$sgCount <- counts_sg



gdata_filt <- subset(gdata_s,sgCount>1)


gset <- gdata_filt$id[1:250]


library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)

dbs <- c("ChEA_2022","Reactome_Pathways_2024","GO_Biological_Process_2023")


#enriched <- enrichr(gset, dbs, background = background, include_overlap = TRUE)
enriched <- enrichr(gset, dbs)

pdf(paste0('~/Documents/TSC_Paper/', 'CDAA_CRISPR_Screen_Deplete.pdf'), width=7, height=4)
p2<- ggplot(enriched[[3]][order(enriched[[2]]$Adjusted.P.value,decreasing=F),][rev(1:5),], 
  (aes(x=Adjusted.P.value, y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + 
geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + 
ggtitle('Reactome 2024') + 
scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[3]][order(enriched[[3]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"_"), `[`, 1))) + 
theme(axis.text=element_text(colour="black"))
p2
dev.off()

#######################################
#############  FIGURE S1E #############
#######################################









gdata_s <- gdata %>% arrange(Score)
sdata_s <- sdata %>% arrange(LFC)

sgRNA_sig <- subset(sdata_s,FDR < 0.05 & LFC < 0)

counts_sg <- table(sgRNA_sig$Gene)
counts_sg <- counts_sg[names(counts_sg) %in% gdata_s$id]


gdata_s$sgCount <- 0


idx<-match(gdata_s$id,names(counts_sg))
idx <- idx[!is.na(idx)]
gdata_s[idx,]$sgCount <- counts_sg



gdata_filt <- subset(gdata_s,sgCount>0)


gset <- gdata_filt$id[1:125]


library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)

dbs <- c("ChEA_2022","Reactome_Pathways_2024","GO_Biological_Process_2023")


#enriched <- enrichr(gset, dbs, background = background, include_overlap = TRUE)
enriched <- enrichr(gset, dbs)

pdf(paste0('~/Documents/TSC_Paper/', 'CDAA_CRISPR_Screen_Deplete.pdf'), width=7, height=4)
p1<- ggplot(enriched[[2]][order(enriched[[2]]$Adjusted.P.value,decreasing=F),][rev(1:5),], 
  (aes(x=Adjusted.P.value, y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + 
geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + 
ggtitle('Reactome 2024') + 
scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[2]][order(enriched[[2]]$Adjusted.P.value,decreasing=F),][rev(1:5),]$Term,"_"), `[`, 1))) + 
theme(axis.text=element_text(colour="black"))
p1
dev.off()







plot_frame <- data.frame(gene = gdata$id, 
  stat = gdata$Score,
  rank=length(gdata$Score):1,
  essential=gdata$essential,
  essential_alt=gdata$essential_alt,
  essential_combo=gdata$essential | gdata$essential_alt,
  TSG=gdata$TSG)


p1 <- ggplot(plot_frame, aes(x=rank,y=stat)) + 
geom_vline(xintercept= which(plot_frame$essential),size=0.1) + theme_void()





























#mitf_gse61966 <- c('SHC4','COPB2','PDE3A','ATG9A','BRAT1',
#  'MAT2A','TMEM184B','ITPRIPL2','PPFIBP1','DDOST','RAB31','NCDN','ANKRD52')
#mitf_gse61966 <- human2mouse[match(mitf_gse61966,human2mouse$human_name),]$mouse_name

#gse61966 <- read.csv('~/Documents/TSC_Paper/GSE61966.top.table.tsv',sep='\t')
#gse61966 <- subset(gse61966,padj < 0.05)
#gse61966_Up <- subset(gse61966,log2FoldChange > 1)
#gse_up <- human2mouse[match(gse61966_Up$Symbol,human2mouse$human_name),]$mouse_name
#gse_up <- gse_up[!is.na(gse_up)]

#gse <- human2mouse[match(gse61966$Symbol,human2mouse$human_name),]$mouse_name
#gse61966$mouse <- gse


autophagy <- read.csv('~/Documents/TSC_Paper/Autophagy_Reactome.csv',sep=',')
autophagy <- autophagy$Gene
autophagy <- human2mouse[match(autophagy,human2mouse$human_name),]$mouse_name
autophagy <- autophagy[!is.na(autophagy)]

cat(names(geneList[autophagy][geneList[autophagy] > 1]),sep='\n')

gene_of_int = autophagy
p2 <- RankView(geneList,cutoff=c(-1.3,1.3),top = 0, bottom = 0,size=14,max.overlaps=Inf,genelist=gene_of_int)
p2 <- p2 + theme(text = element_text(size = 14))



gdata$RandomIndex = sample(1:nrow(gdata), nrow(gdata))
gdata = gdata[order(-gdata$Score), ]
gg = gdata[gdata$Score>0, ]
p1 = ScatterView(gg, x = "RandomIndex", y = "Score", label = "id",
                 y_cut = CutoffCalling(gdata$Score,1.3), 
                 groups = "top", top = 12, ylab = "-Log10RRA",max.overlaps=Inf)



gene_of_int = gse_up
p2 <- RankView(geneList,cutoff=c(-1.3,1.3),top = 0, bottom = 0,size=14,max.overlaps=Inf,genelist=gene_of_int)
p2 <- p2 + theme(text = element_text(size = 14))


gg = gdata[gdata$Score<0, ]

p2 = ScatterView(gg, x = "RandomIndex", y = "Score", label = "id",
                 y_cut = CutoffCalling(gdata$Score,-1.3), 
                 groups = "bottom", top = 12, ylab = "-Log10RRA",max.overlaps=Inf)

pdf(paste0('~/Documents/TSC_Paper/', 'CRISPR_Screen_RandomIndex.pdf'), width=3, height=7)
p1/p2
dev.off()


#######################################
#############  FIGURE 1D  #############
#######################################











geneList= gdata$Score
names(geneList) = gdata$id
p2 <- RankView(geneList,cutoff=c(-1.3,1.3),top = 10, bottom = 10,size=14,max.overlaps=Inf)
p2 <- p2 + theme(text = element_text(size = 14))


pdf(paste0('~/Documents/TSC_Paper/', 'CDAA_CRISPR_Screen_RankView.pdf'), width=5, height=5)
print(p2)
dev.off()
gene_of_int = c('Tsc1','Tsc2','Rraga','Rragb','Rragc','Rragd','Flcn',
  'Lamtor1','Lamtor2','Lamtor3','Lamtor4','Lamtor5','Tfe3','Tfeb',
  'Mitf','Mtor',"Rptor","Fnip1","Fnip2")
p2 <- RankView(geneList,cutoff=c(-1.3,1.3),top = 0, bottom = 0,size=14,max.overlaps=Inf,genelist=gene_of_int)
p2 <- p2 + theme(text = element_text(size = 14))

#cat(names(tail(sort(geneList),n=90)),sep='\n')

pdf(paste0('~/Documents/TSC_Paper/', 'CDAA_CRISPR_Screen_RankView_mtor.pdf'), width=5, height=5)
print(p2)
dev.off()



depmapDat = LoadDepmap()

Depmap = depmapDat$Depmap
sampleinfo = depmapDat$sampleinfo

#"bile duct cancer" , 'liver cancer'
sampleinfo$primary_disease = tolower(sampleinfo$primary_disease)
sampleinfo$subtype_disease = tolower(sampleinfo$subtype_disease)
sampleinfo$cell_line = tolower(sampleinfo$cell_line)
sampleinfo = sampleinfo[colnames(Depmap), ]
idx = which(sampleinfo$primary_disease == 'liver cancer' & sampleinfo$subtype_disease == "hepatocellular carcinoma")
Depmap = Depmap

genes = intersect(gdata$HumanGene, rownames(Depmap))

dd = gdata[!duplicated(gdata[, 'HumanGene']), ]
rownames(dd) = dd[, 'HumanGene']

similarity = apply(Depmap[genes, ], 2, function(x) {
            tmp = cor.test(x, dd[genes, 'Score'], method = 'pearson', 
                na.action = na.omit)
            c(tmp$estimate, tmp$p.value)
})

similarity = as.data.frame(t(similarity))
colnames(similarity) = c("estimate", "p.value")
#rownames(similarity) = sampleinfo[colnames(Depmap), 1]
similarity = similarity[order(-similarity$estimate), ]



gdata = ReadRRA(file1,score='lfc')
gdata$HumanGene = TransGeneID(gdata$id, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")
## Remove missing or duplicate human genes
idx2 = duplicated(gdata$HumanGene)|is.na(gdata$HumanGene)
gdata = gdata[!idx2, ]

dd = gdata[!duplicated(gdata[, 'HumanGene']), ]
rownames(dd) = dd[, 'HumanGene']

library(matrixStats)
depmap_avg <- rowMedians(as.matrix(Depmap[genes,idx]))
#depmap_avg <- -log10(max(depmap_avg+1) - depmap_avg)*10
#depmap_avg <- depmap_avg - median(depmap_avg)
gdata_cc <- data.frame(Gene=genes,screen=dd[genes, 'Score'],depmap=depmap_avg)

p1 = ScatterView(gdata_cc, x="depmap", y="screen", label = "Gene", 
                 model = "ninesquare", top = 10, slope=1,
                 display_cut = TRUE, y_cut = c(-1,1),
                 x_cut=c(-0.5,0.25),auto_cut_diag=1,max.overlaps=Inf)
pdf(paste0('~/Documents/TSC_Paper/', 'CRISPR_Screen_Depmap.pdf'), width=5, height=5)

print(p1)

dev.off()


#######################################
#############  FIGURE 1E  #############
#######################################
#mageck test -k Combined_GuideCounts.txt -t 3,4,5 -c 0,1,2 -n chow_vs_lib --norm-method total --remove-zero both --paired
#mageck test -k Combined_GuideCounts_with_CDAA.txt -t 6,7,8 -c 3,4,5 -n cdaa_vs_chow --norm-method control --remove-zero both --pdf-report --control-sgrna cntrl_sgrna.txt 
#cat(names(tail(a,n=100)),sep='\n') Top 100 genes give interesting results

library(MAGeCKFlute)
library(clusterProfiler)
library(ggplot2)


file1 = file.path("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/cdaa_vs_chow_keepzero.gene_summary.txt")
file2 = file.path("/Users/ikuz/Documents/PinAPL/mageck/CountFiles/cdaa_vs_chow_keepzero.sgrna_summary.txt")


gdata = ReadRRA(file1,score='rra')
sdata = ReadsgRRA(file2)
gdata$HumanGene = TransGeneID(gdata$id, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")
sdata$HumanGene = TransGeneID(sdata$Gene, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")

## Remove missing or duplicate human genes
idx = duplicated(gdata$HumanGene)|is.na(gdata$HumanGene)
gdata = gdata[!idx, ]
#depmap_similarity = ResembleDepmap(gdata, symbol = "HumanGene", score = "Score")
#head(depmap_similarity)

gdata2<-OmitCommonEssential(gdata, symbol = "HumanGene")
alt_essential <- setdiff(gdata$id,gdata2$id)

essential <- read.table('~/Documents/PinAPL/mageck/CountFiles/EssentialGenes.csv',header=FALSE,sep='\t')
human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
idx <- match(unique(human2mouse[,2]),human2mouse[,2])
human2mouse<-human2mouse[idx,]
colnames(human2mouse) <-c('human_name', 'mouse_name')
idx_match<- match(essential$V1,human2mouse$human_name)
essential_mouse <- human2mouse$mouse_name[idx_match]
essential_mouse <- essential_mouse[!is.na(essential_mouse)]


TSG <- read.table('~/Documents/PinAPL/mageck/CountFiles/TSG_Genes.csv',header=TRUE,sep=',')
idx_match<- match(TSG$Gene,human2mouse$human_name)
TSG_mouse <- human2mouse$mouse_name[idx_match]
TSG$mouse <- TSG_mouse
TSG <- TSG[!is.na(TSG$mouse),]


gdata$essential <- !is.na(match(gdata$id,essential_mouse))
gdata$essential_alt <- !is.na(match(gdata$id,alt_essential))
gdata$essential_combo <- gdata$essential | gdata$essential_alt

gdata$TSG <- !is.na(match(gdata$id,TSG$mouse))

geneList= gdata$Score
names(geneList) = gdata$id
p2 <- RankView(geneList,cutoff=c(-1.3,1.3),top = 10, bottom = 10,size=14,max.overlaps=Inf)
p2 <- p2 + theme(text = element_text(size = 14))


pdf(paste0('~/Documents/TSC_Paper/', 'CDAA_CRISPR_Screen_RankView.pdf'), width=5, height=5)
print(p2)
dev.off()
gene_of_int = c('Tsc1','Tsc2','Rraga','Rragb','Rragc','Rragd','Flcn',
  'Lamtor1','Lamtor2','Lamtor3','Lamtor4','Lamtor5','Tfe3','Tfeb',
  'Mitf','Mtor',"Rptor","Fnip1","Fnip2")
p2 <- RankView(geneList,cutoff=c(-1.3,1.3),top = 0, bottom = 0,size=14,max.overlaps=Inf,genelist=gene_of_int)
p2 <- p2 + theme(text = element_text(size = 14))

#cat(names(tail(sort(geneList),n=90)),sep='\n')

pdf(paste0('~/Documents/TSC_Paper/', 'CDAA_CRISPR_Screen_RankView_mtor.pdf'), width=5, height=5)
print(p2)
dev.off()

library("readxl")

Mouse.Mito <- read_excel("~/Downloads/hdWGCNA_TOM/Mouse.MitoCarta3.0.xls", sheet = "A Mouse MitoCarta3.0")

gdata$mito <- !is.na(match(gdata$id,Mouse.Mito$Symbol))

plot_frame <- data.frame(gene = gdata$id, 
  stat = gdata$Score,
  rank=length(gdata$Score):1,
  essential=gdata$essential,
  essential_alt=gdata$essential_alt,
  essential_combo=gdata$essential | gdata$essential_alt,
  TSG=gdata$TSG,
  Mito = gdata$mito)

p1 <- ggplot(plot_frame, aes(x=rank,y=stat)) + 
geom_vline(xintercept= which(plot_frame$essential),size=0.1) + theme_void()

#p1 <- ggplot(plot_frame, aes(x=rank,y=stat)) + 
#geom_vline(xintercept=which(plot_frame$essential_alt),size=0.1) + theme_void()

#p1 <- ggplot(plot_frame, aes(x=rank,y=stat)) + 
#geom_vline(xintercept=which(plot_frame$essential_combo) ,size=0.1) + theme_void()

pdf(paste0('~/Documents/TSC_Paper/', 'CDAA_CRISPR_Screen_RankView_essential.pdf'), width=5, height=0.5)
print(p1)
dev.off()



gene_of_int = c('Tsc1','Tsc2','Rraga','Rragb','Rragc','Rragd','Flcn',
  'Lamtor1','Lamtor2','Lamtor3','Lamtor4','Lamtor5','Tfe3','Tfeb',
  'Mitf','Mtor',"Rptor",'Deptor','Rictor',"Fnip1","Fnip2","Akt1","Akt2",'Akt3',
  'Rps6kb1','Eif4ebp1','Rheb','Akt1s1','Mlst8','Pik3ca','Pik3r1','Pten')
p2 <- RankView(geneList,cutoff=c(-1.3,1.3),top = 0, bottom = 0,size=14,max.overlaps=Inf,genelist=gene_of_int)
p2 <- p2 + theme(text = element_text(size = 14))

#cat(names(tail(sort(geneList),n=90)),sep='\n')

pdf(paste0('~/Documents/TSC_Paper/', 'CRISPR_Screen_RankView_mtor.pdf'), width=5, height=5)
print(p2)
dev.off()





