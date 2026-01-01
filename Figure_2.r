###################
### TCGA RNA-Seq###
###################


TCGA.RNA.drivers <- read.csv('~/Documents/TSC_Paper/HumanData/TCGA_RNASeq_TSC_Drivers.tsv',sep='\t')
RNA.up <- subset(TCGA.RNA.drivers,Higher.expression.in == '(B) TSC1_Driver' & q.Value < 0.05)$Gene
human2mouse <- read.csv('~/Downloads/hdWGCNA_TOM/human2mouse.csv',header=F)
TCGA.RNA.drivers$Log2.Ratio <- -TCGA.RNA.drivers$Log2.Ratio
MITF_all <- read.table('~/Documents/TSC_Paper/HumanData/MITF_geneset.txt')
MITF_all <- stringr::str_match_all(MITF_all, "symbol\":\"\\s*(.*?)\\s*\",\"href")
MITF_all<-MITF_all[[6]][,2]

set.seed(1)
MITF_all_subset <- MITF_all[sample(length(MITF_all), 1000)]


CLEAR <- read.csv('~/Documents/TSC_Paper/CLEAR_genes.csv',header=FALSE)$V1
CLEAR_human <- human2mouse$V1[match(CLEAR,human2mouse$V2)]
CLEAR_human <- CLEAR_human[!is.na(CLEAR_human)]

library(enrichR)
library(EnhancedVolcano)
library(DOSE)
library(forcats)
library(viridis)

dbs <- c("ChEA_2022","Reactome_Pathways_2024","GO_Biological_Process_2023")
enriched.up <- enrichr(RNA.up, dbs)

enriched.up[[1]] <- subset(enriched.up[[1]],Adjusted.P.value<0.05)

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

pdf('~/Documents/TSC_Paper/TCGA_TSC_RNA_enrichr.pdf',width=4,height=2.5)
p2<- ggplot(enriched.up[[1]][order(enriched.up[[1]]$Adjusted.P.value,decreasing=F),][1:8,], 
  (aes(x=-log(Adjusted.P.value), y=fct_inorder(Term), color = as.numeric(Combined.Score), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Adjusted P Value') + 
  ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched.up[[1]][order(enriched.up[[1]]$Adjusted.P.value,decreasing=F),][1:8,]$Term," \\s*"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256))) + coord_flip()
  p2
dev.off()

MITF <- enriched.up[[1]][1,]$Genes
MITF <- strsplit(MITF,';')[[1]]

keyvals <- ifelse(TCGA.RNA.drivers$Gene %in% CLEAR_human,'red','gray')
names(keyvals)[keyvals == 'red'] <- 'CLEAR'
names(keyvals)[keyvals == 'gray'] <- 'Other'

keyvals_alt <- ifelse(TCGA.RNA.drivers$Gene %in% MITF_all_subset,'red','gray')
names(keyvals_alt)[keyvals_alt == 'red'] <- 'MITF'
names(keyvals_alt)[keyvals_alt == 'gray'] <- 'Other'


colalp <- ifelse(TCGA.RNA.drivers$Gene %in% CLEAR_human,0.5,0.1)
TCGA.RNA.drivers$color <- keyvals
TCGA.RNA.drivers$colalp <- colalp

colalp <- ifelse(TCGA.RNA.drivers$Gene %in% MITF_all_subset,0.5,0.1)
TCGA.RNA.drivers$color_alt <- keyvals_alt
TCGA.RNA.drivers$colalp_alt <- colalp

TCGA.RNA.drivers.plot <- TCGA.RNA.drivers[order(keyvals),]

p1 <- EnhancedVolcano(TCGA.RNA.drivers.plot,x = 'Log2.Ratio',y = 'q.Value',lab = TCGA.RNA.drivers.plot$Gene,
  pCutoff = 0.05, FCcutoff=0.25,colCustom = sort(keyvals),
  selectLab=TCGA.RNA.drivers.plot$Gene[sort(keyvals) != 'gray'],
  colAlpha=TCGA.RNA.drivers.plot$colalp,ylim=c(0,5),
  xlim=c(-4,4))

TCGA.RNA.drivers.plot <- TCGA.RNA.drivers[order(keyvals_alt),]

p2 <- EnhancedVolcano(TCGA.RNA.drivers.plot,x = 'Log2.Ratio',y = 'q.Value',lab = TCGA.RNA.drivers.plot$Gene,
  pCutoff = 0.05, FCcutoff=0.25,colCustom = sort(keyvals_alt),
  selectLab=subset(TCGA.RNA.drivers.plot[sort(keyvals_alt) != 'gray',],q.Value<0.05)$Gene,
  colAlpha=TCGA.RNA.drivers.plot$colalp_alt,ylim=c(0,5),
  xlim=c(-4,4))

library(patchwork)
library(rasterpdf)

pdf('~/Documents/TSC_Paper/TCGA_TSC_RNA_volcano.pdf',width=6,height=6.5)
p1 + p2
dev.off()

#######################
### TCGA Proteomics ###
#######################
library(dplyr)
library(data.table)


TCGA.Prot.drivers <- read.csv('~/Documents/TSC_Paper/HumanData/TCGA_Proteomics_TSC_Drivers.tsv',sep='\t')
TOR_Genes <- c('RICTOR','MTOR','RPS6','TSC2','EIF4EBP1')
AKT_Genes <- c('AKT1','AKT2','AKT3')


idx <- TCGA.Prot.drivers$Gene %like% '_'
prot_split <- unlist(lapply(strsplit(TCGA.Prot.drivers$Gene,'_'),'[[',1))
TCGA.Prot.drivers$Name <- prot_split

keyvals <- ifelse(TCGA.Prot.drivers$Name %in% TOR_Genes,'red',ifelse(TCGA.Prot.drivers$Name %in% AKT_Genes,'green','gray'))
names(keyvals)[keyvals == 'red'] <- 'mTOR'
names(keyvals)[keyvals == 'green'] <- 'AKT'
names(keyvals)[keyvals == 'gray'] <- 'Other'

colalp <- ifelse(TCGA.Prot.drivers$Name %in% TOR_Genes,0.5,ifelse(TCGA.Prot.drivers$Name %in% AKT_Genes,0.5,0.1))

TCGA.Prot.drivers.Prot <- TCGA.Prot.drivers[!idx,] 
TCGA.Prot.drivers.Phos <- TCGA.Prot.drivers[idx,] 


pdf('~/Documents/TSC_Paper/TCGA_TSC_Prot_volcano.pdf',width=6,height=6.5)

p1 <- EnhancedVolcano(TCGA.Prot.drivers.Prot,x = 'Log2.Ratio',y = 'p.Value',lab = TCGA.Prot.drivers.Prot$Gene,
  pCutoff = 100, FCcutoff=100,ylim=c(0,2),xlim=c(-1.1,1.1),colCustom = keyvals[!idx],
  selectLab=TCGA.Prot.drivers.Prot$Name[TCGA.Prot.drivers.Prot$Name %in% c(TOR_Genes,AKT_Genes)],
  drawConnectors = TRUE,typeConnectors="closed")

p2 <- EnhancedVolcano(TCGA.Prot.drivers.Phos,x = 'Log2.Ratio',y = 'p.Value',lab = TCGA.Prot.drivers.Phos$Gene,
  pCutoff = 100, FCcutoff=100,ylim=c(0,2),xlim=c(-1.1,1.1),colCustom = keyvals[idx],
  selectLab=TCGA.Prot.drivers.Phos$Gene[TCGA.Prot.drivers.Phos$Name %in% c(TOR_Genes,AKT_Genes)],
  drawConnectors = TRUE,typeConnectors="closed")
p1 + p2

dev.off()

n_TSC = 11
n_other = 168

#######S6K
pRPS6 <- subset(TCGA.Prot.drivers.Phos,Name=='RPS6')
pRPS6.TSC.mean <- (pRPS6[1,3] + pRPS6[2,3])/2
pRPS6.TSC.ste <- (sqrt((pRPS6[1,5]^2 + pRPS6[2,5]^2)/2))/sqrt(n_TSC * 2)
pRPS6.Other.mean <- (pRPS6[1,4] + pRPS6[2,4])/2
pRPS6.Other.ste <- (sqrt((pRPS6[1,6]^2 + pRPS6[2,6]^2)/2))/sqrt(n_other * 2)
pRPS6.diff <- pRPS6.TSC.mean - pRPS6.Other.mean
pRPS6.diff.ste <- sqrt(pRPS6.TSC.ste^2 + pRPS6.Other.ste^2)
Z = pRPS6.diff/pRPS6.diff.ste
p = pnorm(q=Z, lower.tail=FALSE)*2

RPS6 <- subset(TCGA.Prot.drivers.Prot,Name=='RPS6')
RPS6.TSC.mean <- RPS6[1,3]
RPS6.TSC.ste <- pRPS6[1,5]/sqrt(n_TSC)
RPS6.Other.mean <- RPS6[1,4]
RPS6.Other.ste <- RPS6[1,6]/sqrt(n_other)
RPS6.diff <- RPS6.TSC.mean - RPS6.Other.mean
RPS6.diff.ste <- sqrt(RPS6.TSC.ste^2 + RPS6.Other.ste^2)
Z = RPS6.diff/RPS6.diff.ste
p = pnorm(q=Z, lower.tail=FALSE)*2

#RPS6.ratio <- (pRPS6.TSC.mean / RPS6.TSC.mean) - (pRPS6.Other.mean / RPS6.Other.mean)
boot_pRPS6.TSC <- rnorm(10000,pRPS6.TSC.mean,pRPS6.TSC.ste)
boot_RPS6.TSC <- rnorm(10000,RPS6.TSC.mean,RPS6.TSC.ste)
boot_pRPS6.Other <- rnorm(10000,pRPS6.Other.mean,pRPS6.Other.ste)
boot_RPS6.Other <- rnorm(10000,RPS6.Other.mean,RPS6.Other.ste)
boot_distro <- boot_pRPS6.TSC/boot_RPS6.TSC - boot_pRPS6.Other/boot_RPS6.Other
p = sum(boot_distro <= 0)/10000
#p = 0.0154


df1 <- data.frame(group = 'TSC',data='Phospho',mu_1 = pRPS6.TSC.mean + 1, se_1 = pRPS6.TSC.ste)
df2 <- data.frame(group = 'Other',data='Phospho',mu_1 = pRPS6.Other.mean + 1, se_1 = pRPS6.Other.ste)
df3 <- data.frame(group = 'TSC',data='Prot',mu_1 = RPS6.TSC.mean + 1, se_1 = RPS6.TSC.ste)
df4 <- data.frame(group = 'Other',data='Prot',mu_1 = RPS6.Other.mean + 1, se_1 = RPS6.Other.ste)
df <- rbind(df1,df2,df3,df4)

p1 <- ggplot(df,aes(x=group,y=mu_1,fill=group)) + 
geom_col() + facet_wrap(vars(data)) + ylim(c(0,2)) +
geom_errorbar(aes(x=group, ymin=mu_1-se_1, ymax=mu_1+se_1), width=0.4, position = position_dodge(.9)) + 
theme_classic()

TSC.distro <- boot_pRPS6.TSC/boot_RPS6.TSC 
Other.distro <- boot_pRPS6.Other/boot_RPS6.Other
plot.distro <- rbind(data.frame(group='TSC',val = TSC.distro),data.frame(group='Other',val = Other.distro))

p2 <- ggplot(plot.distro,aes(x=group ,y=val,fill=group)) + geom_violin() + 
ylim(-2,2) + theme_classic()

pdf('~/Documents/TSC_Paper/TCGA_TSC_Prot_RPS6.pdf',width=6,height=3.5)
p1 + p2
dev.off()

#######4EBP1
pRPS6 <- subset(TCGA.Prot.drivers.Phos,Name=='EIF4EBP1')
pRPS6.TSC.mean <- (pRPS6[1,3] + pRPS6[2,3] + pRPS6[3,3])/3
pRPS6.TSC.ste <- (sqrt((pRPS6[1,5]^2 + pRPS6[2,5]^2 + pRPS6[3,5]^2)/3))/sqrt(n_TSC * 3)
pRPS6.Other.mean <- (pRPS6[1,4] + pRPS6[2,4] + pRPS6[3,4])/3
pRPS6.Other.ste <- (sqrt((pRPS6[1,6]^2 + pRPS6[2,6]^2 + pRPS6[3,6]^2)/3))/sqrt(n_other * 3)
pRPS6.diff <- pRPS6.TSC.mean - pRPS6.Other.mean
pRPS6.diff.ste <- sqrt(pRPS6.TSC.ste^2 + pRPS6.Other.ste^2)
Z = pRPS6.diff/pRPS6.diff.ste
p = pnorm(q=Z, lower.tail=FALSE)*2

RPS6 <- subset(TCGA.Prot.drivers.Prot,Name=='EIF4EBP1')
RPS6.TSC.mean <- RPS6[1,3]
RPS6.TSC.ste <- pRPS6[1,5]/sqrt(n_TSC)
RPS6.Other.mean <- RPS6[1,4]
RPS6.Other.ste <- RPS6[1,6]/sqrt(n_other)
RPS6.diff <- RPS6.TSC.mean - RPS6.Other.mean
RPS6.diff.ste <- sqrt(RPS6.TSC.ste^2 + RPS6.Other.ste^2)
Z = RPS6.diff/RPS6.diff.ste
p = pnorm(q=Z, lower.tail=FALSE)*2

#RPS6.ratio <- (pRPS6.TSC.mean / RPS6.TSC.mean) - (pRPS6.Other.mean / RPS6.Other.mean)
boot_pRPS6.TSC <- rnorm(10000,pRPS6.TSC.mean,pRPS6.TSC.ste)
boot_RPS6.TSC <- rnorm(10000,RPS6.TSC.mean,RPS6.TSC.ste)
boot_pRPS6.Other <- rnorm(10000,pRPS6.Other.mean,pRPS6.Other.ste)
boot_RPS6.Other <- rnorm(10000,RPS6.Other.mean,RPS6.Other.ste)
boot_distro <- boot_pRPS6.TSC/boot_RPS6.TSC - boot_pRPS6.Other/boot_RPS6.Other
p = sum(boot_distro <= 0)/10000
#p = 0.7986


df1 <- data.frame(group = 'TSC',data='Phospho',mu_1 = pRPS6.TSC.mean + 1, se_1 = pRPS6.TSC.ste)
df2 <- data.frame(group = 'Other',data='Phospho',mu_1 = pRPS6.Other.mean + 1, se_1 = pRPS6.Other.ste)
df3 <- data.frame(group = 'TSC',data='Prot',mu_1 = RPS6.TSC.mean + 1, se_1 = RPS6.TSC.ste)
df4 <- data.frame(group = 'Other',data='Prot',mu_1 = RPS6.Other.mean + 1, se_1 = RPS6.Other.ste)
df <- rbind(df1,df2,df3,df4)

p1 <- ggplot(df,aes(x=group,y=mu_1,fill=group)) + 
geom_col() + facet_wrap(vars(data)) + ylim(c(0,2)) +
geom_errorbar(aes(x=group, ymin=mu_1-se_1, ymax=mu_1+se_1), width=0.4, position = position_dodge(.9)) + 
theme_classic()

TSC.distro <- boot_pRPS6.TSC/boot_RPS6.TSC 
Other.distro <- boot_pRPS6.Other/boot_RPS6.Other
plot.distro <- rbind(data.frame(group='TSC',val = TSC.distro),data.frame(group='Other',val = Other.distro))

p2 <- ggplot(plot.distro,aes(x=group ,y=val,fill=group)) + geom_violin() + 
ylim(0,4) + theme_classic()

pdf('~/Documents/TSC_Paper/TCGA_TSC_Prot_4EBP.pdf',width=6,height=3.5)
p1 + p2
dev.off()



#######AKT
pRPS6 <- subset(TCGA.Prot.drivers.Phos,Name %like% 'AKT')
pRPS6.TSC.mean <- (pRPS6[1,3] + pRPS6[2,3] + pRPS6[3,3] pRPS6[5,3] + pRPS6[6,3] + pRPS6[7,3])/6
pRPS6.TSC.ste <- (sqrt((pRPS6[1,5]^2 + pRPS6[2,5]^2 + pRPS6[3,5]^2 + pRPS6[5,5]^2 + pRPS6[6,5]^2 + pRPS6[7,5]^2)/6))/sqrt(n_TSC * 6)
pRPS6.Other.mean <- (pRPS6[1,4] + pRPS6[2,4] + pRPS6[3,4] + pRPS6[5,4] + pRPS6[6,4] + pRPS6[7,4])/6
pRPS6.Other.ste <- (sqrt((pRPS6[1,6]^2 + pRPS6[2,6]^2 + pRPS6[3,6]^2 + pRPS6[5,6]^2 + pRPS6[6,6]^2 + pRPS6[7,6]^2)/6))/sqrt(n_TSC * 6)
pRPS6.diff <- pRPS6.TSC.mean - pRPS6.Other.mean
pRPS6.diff.ste <- sqrt(pRPS6.TSC.ste^2 + pRPS6.Other.ste^2)
Z = pRPS6.diff/pRPS6.diff.ste
p = pnorm(q=Z, lower.tail=FALSE)*2

RPS6 <- subset(TCGA.Prot.drivers.Prot,Name %like% 'AKT')
RPS6.TSC.mean <- (RPS6[1,3] + RPS6[2,3] + RPS6[3,3])/3
RPS6.TSC.ste <- (sqrt((RPS6[1,5]^2 + RPS6[2,5]^2 + RPS6[3,5]^2)/3))/sqrt(n_TSC * 3)
RPS6.Other.mean <- (RPS6[1,4] + RPS6[2,4] + RPS6[3,4])/3
RPS6.Other.ste <- (sqrt((RPS6[1,6]^2 + RPS6[2,6]^2 + RPS6[3,6]^2)/3))/sqrt(n_other * 3)
RPS6.diff <- RPS6.TSC.mean - RPS6.Other.mean
RPS6.diff.ste <- sqrt(RPS6.TSC.ste^2 + RPS6.Other.ste^2)
Z = RPS6.diff/RPS6.diff.ste
p = pnorm(q=Z, lower.tail=TRUE)*2

#RPS6.ratio <- (pRPS6.TSC.mean / RPS6.TSC.mean) - (pRPS6.Other.mean / RPS6.Other.mean)
boot_pRPS6.TSC <- rnorm(10000,pRPS6.TSC.mean,pRPS6.TSC.ste)
boot_RPS6.TSC <- rnorm(10000,RPS6.TSC.mean,RPS6.TSC.ste)
boot_pRPS6.Other <- rnorm(10000,pRPS6.Other.mean,pRPS6.Other.ste)
boot_RPS6.Other <- rnorm(10000,RPS6.Other.mean,RPS6.Other.ste)
boot_distro <- boot_pRPS6.TSC/boot_RPS6.TSC - boot_pRPS6.Other/boot_RPS6.Other
p = sum(boot_distro > 0)/10000
#p = 0.0359


df1 <- data.frame(group = 'TSC',data='Phospho',mu_1 = pRPS6.TSC.mean + 1, se_1 = pRPS6.TSC.ste)
df2 <- data.frame(group = 'Other',data='Phospho',mu_1 = pRPS6.Other.mean + 1, se_1 = pRPS6.Other.ste)
df3 <- data.frame(group = 'TSC',data='Prot',mu_1 = RPS6.TSC.mean + 1, se_1 = RPS6.TSC.ste)
df4 <- data.frame(group = 'Other',data='Prot',mu_1 = RPS6.Other.mean + 1, se_1 = RPS6.Other.ste)
df <- rbind(df1,df2,df3,df4)

p1 <- ggplot(df,aes(x=group,y=mu_1,fill=group)) + 
geom_col() + facet_wrap(vars(data)) + ylim(c(-1,2)) +
geom_errorbar(aes(x=group, ymin=mu_1-se_1, ymax=mu_1+se_1), width=0.4, position = position_dodge(.9)) + 
theme_classic()

TSC.distro <- boot_pRPS6.TSC/boot_RPS6.TSC 
Other.distro <- boot_pRPS6.Other/boot_RPS6.Other
plot.distro <- rbind(data.frame(group='TSC',val = TSC.distro),data.frame(group='Other',val = Other.distro))

p2 <- ggplot(plot.distro,aes(x=group ,y=val,fill=group)) + geom_violin() + 
ylim(-6,0) + theme_classic()

pdf('~/Documents/TSC_Paper/TCGA_TSC_Prot_AKT.pdf',width=6,height=3.5)
p1 + p2
dev.off()


###################
### CLCA RNA-Seq###
###################

CLCA.AKT.LOW <- read.csv('~/Documents/TSC_Paper/HumanData/CLCA_RNASeq_AKT2_low_TSC_Mut_Better.tsv',sep='\t')
RNA.up <- subset(CLCA.AKT.LOW,Higher.expression.in == '(A) <=-0.50' & q.Value < 0.1)$Gene

library(enrichR)
library(EnhancedVolcano)
library(DOSE)
library(forcats)
library(viridis)

dbs <- c("ChEA_2022","Reactome_Pathways_2024","GO_Biological_Process_2023")
enriched.up <- enrichr(RNA.up, dbs)

enriched.up[[1]] <- subset(enriched.up[[1]],Adjusted.P.value<0.05)

wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

pdf('~/Documents/TSC_Paper/CLCA_TSC_RNA_enrichr.pdf',width=3,height=2)
p2<- ggplot(enriched.up[[1]][order(enriched.up[[1]]$Adjusted.P.value,decreasing=F),][1:3,], 
  (aes(x=-log(Adjusted.P.value), y=fct_inorder(Term), color = as.numeric(Combined.Score), 
  size=parse_ratio(Overlap)))) + geom_point() + xlab('Adjusted P Value') + 
  ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + 
  scale_y_discrete(labels= fct_inorder(
    wrapText(sapply(
      strsplit(enriched.up[[1]][order(enriched.up[[1]]$Adjusted.P.value,decreasing=F),][1:3,]$Term," \\s*"),
         `[`, 1),35))) + 
  theme(axis.text=element_text(colour="black"))+
  scale_color_stepsn(colors=rev(magma(256))) + coord_flip()
  p2
dev.off()




