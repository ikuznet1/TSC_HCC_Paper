library(PhosR)
library(stringr)
library(org.Mm.eg.db)
library(UniprotR)

library(seqinr)

translations <- read.fasta(file ="~/Documents/TSC_Paper/gencode.vM36.pc_translations.fa",
  seqtype = c("AA"))
names(translations) <- unlist(lapply(str_split(names(translations),'\\|'),'[[',7))


#######################
###### LOAD DATA ######
#######################

data <- read.csv('~/Documents/TSC_Paper/Phosphoproteomics_peptides.csv')
data <- data[data$Modifications.in.Master.Proteins != "",]
data <- data[!is.na(data$Abundances..Normalized...F1..131..Sample..KO),]


# Solit multi-phos peptides
mods <- data$Modifications.in.Master.Proteins
mods <- unlist(lapply(str_split(mods,'];'),'[[',1))
phos_sites_peptide <- unlist(lapply(str_split(mods,'xPhospho'),'[[',1))
phos_sites_peptide <- as.numeric(unlist(lapply(str_split(phos_sites_peptide,'\\ '),'[[',2)))
mods <- unlist(lapply(str_split(mods,'xPhospho '),'[[',2))
mods <- gsub("\\[","", mods)
mods <- gsub("\\]","", mods)
mods <- str_split(mods,';\\ ')
temp<-unlist(lapply(mods,length))
to_pad <- phos_sites_peptide - temp

for (i in which(to_pad>0)){
	mods[[i]] <- c(mods[[i]],rep('S',to_pad[i]))
}

mods <- unlist(mods)
site <- as.numeric(gsub("[A-Z]","",mods))
residue <- gsub("[0-9]","",mods)

abundances <- as.matrix(data[,20:29])
peptide <- data$Annotated.Sequence
master.protein <- data$Positions.in.Master.Proteins

abundances <- abundances[rep(seq_along(phos_sites_peptide), phos_sites_peptide), ]
peptide <- peptide[rep(seq_along(phos_sites_peptide), phos_sites_peptide)]
master.protein <- master.protein[rep(seq_along(phos_sites_peptide), phos_sites_peptide)]


abundances <- abundances[!is.na(site),]
residue <- residue[!is.na(site)]
peptide <- peptide[!is.na(site)]
master.protein <- master.protein[!is.na(site)]
site <- site[!is.na(site)]



sample_names <- unlist(lapply(str_split(colnames(abundances),'F1'),'[[',2))
sample_names <- gsub("Sample..", "_", sample_names)
sample_names <- gsub("\\..", "", sample_names)
sample_group <- unlist(lapply(str_split(sample_names,'_'),'[[',2))
colnames(abundances) <- sample_names


peptide <- gsub("\\[", "", peptide)
peptide <- gsub("\\]", "", peptide)
peptide <- gsub("\\.", "", peptide)
peptide <- gsub("-", "", peptide)


master.protein <- unlist(lapply(str_split(master.protein,'\\ '),'[[',1))
names <- read.csv('~/Documents/TSC_Paper/Phospho_omics_mapping.csv')
gene_symbols<- names[match(master.protein,names$Entry),'Entry.Name']
gene_symbols <- unlist(lapply(str_split(gene_symbols,'_'),'[[',1))

alt_symbols <- names[match(master.protein,names$Entry),'Gene.Names']
alt_symbols <- unlist(lapply(str_split(alt_symbols,'\\ '),'[[',1))

labels <- paste0(gene_symbols,';',residue,site,';',peptide)

n_occur <- data.frame(table(labels))
#idx <- n_occur[n_occur$Freq > 1,]

keep <- list()
abundances_final <- matrix(, nrow =  dim(n_occur)[1], ncol = dim(abundances)[2])
j=0
for (i in n_occur$labels) {
	j = j+1
	idx<-which(labels %in% i)
	keep <- c(keep,idx[1])
	if (length(idx) > 1){abundances_final[j,] = colSums(abundances[idx,])}
	else {abundances_final[j,] = abundances[idx,]}
}
keep <- unlist(keep)
colnames(abundances_final) <- colnames(abundances)
rownames(abundances_final) <- labels[keep]
residue <- residue[keep]
peptide <- peptide[keep]
site <- site[keep]

symbols <- gene_symbols[keep]
alt_symbols <- alt_symbols[keep]

flanking <- c()
failed <- c()
#Flanking sequences
for (i in 1:length(site)){
	gene_to_pull <- str_to_title(alt_symbols[i])
	idx <- which(names(translations) == gene_to_pull)
	if (sum(idx) == 0){
		print(paste0('Not present in dataset: ',gene_to_pull))
		flanking <- c(flanking,paste(rep('-',31),collapse=''))
		failed <- c(failed,gene_to_pull)
		next
	}
	candidate_translations <- translations[idx]
	candidate_translations <- lapply(candidate_translations,'paste',collapse='')
	peptide_contained_translations <- unlist(candidate_translations[grepl(peptide[i],candidate_translations, fixed = TRUE)])
	if (length(peptide_contained_translations) == 0){
		partial_pep <- substr(peptide[i],2,nchar(peptide[i])-2)
		peptide_contained_translations <- unlist(candidate_translations[grepl(partial_pep,candidate_translations, fixed = TRUE)])
	}

	if (length(peptide_contained_translations) == 0){
		print(paste0('Unresolvable error (peptide absence in sequence) for ',gene_to_pull))
		flanking <- c(flanking,paste(rep('-',31),collapse=''))
		failed <- c(failed,gene_to_pull)
		next
	}

	idx <- substr(peptide_contained_translations,site[i],site[i]) == residue[i]
	if (sum(idx) != 0){
		translation_pull <- peptide_contained_translations[idx][[1]]
		translation_pull<- paste0(paste(rep('-',15),collapse=''),translation_pull,paste(rep('-',15),collapse=''),collapse='')
		flanking <- c(flanking,substr(translation_pull,site[i],site[i] + 30))
	}
	else{
		print(paste0('Checking for off by 1 error in translation for ',gene_to_pull))
		site_t = site[i] + 1
		idx <- substr(peptide_contained_translations,site_t,site_t) == residue[i]
		if (sum(idx) == 0 ){
			idy <- gregexpr(peptide[i],peptide_contained_translations,fixed=TRUE)
			offset <- gregexpr('S',peptide[i])[[1]][1] - 1
			idy <- unlist(lapply(idy,'[[',1)) + offset
			idz <- which.min(abs(idy - site[i]))
			site_t <- idy[idz]
			translation_pull <- peptide_contained_translations[idz][[1]]
			print(paste0('Manually realigned peptides for ',gene_to_pull))
			translation_pull<- paste0(paste(rep('-',15),collapse=''),translation_pull,paste(rep('-',15),collapse=''),collapse='')
			flanking <- c(flanking,substr(translation_pull,site_t,site_t + 30))
		}
		else if (!grepl(peptide[i],peptide_contained_translations[idx][[1]], fixed = TRUE)){
			print(paste0('Unresolvable error (peptide absence in sequence) for ',gene_to_pull))
			print(substr(peptide_contained_translations,site[i] - 10,site[i] + 10))
			flanking <- c(flanking,paste(rep('-',31),collapse=''))
			failed <- c(failed,gene_to_pull)
			next
		}
		else {
			translation_pull <- peptide_contained_translations[idx][[1]]
			print('Resolved off by 1 error')
			translation_pull<- paste0(paste(rep('-',15),collapse=''),translation_pull,paste(rep('-',15),collapse=''),collapse='')
			flanking <- c(flanking,substr(translation_pull,site_t,site_t + 30))
		}
	}
}



#######################
####### SET UP ########
#######################


ppe <- PhosphoExperiment(assays = list(Quantification = as.matrix(log(abundances_final + 1))),
	Sequence=flanking,
	GeneSymbol=symbols,
	Residue=residue,
	Site=site)

ppe_filtered <- selectGrps(ppe, sample_group, 0.5, n=1) 

#### Don't need to impute
#ppe_imputed_tmp <- scImpute(ppe_filtered, 0.5, sample_group)[,colnames(ppe_filtered)]


p1 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
    labels=colnames(ppe_filtered), 
    panel = "quantify", grps = sample_group)

p1 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
       labels=colnames(ppe_filtered), panel = "dendrogram", 
       grps = sample_group)


plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), grps=sample_group, 
    labels = colnames(ppe_filtered), panel = "pca") +
    ggplot2::ggtitle("Before batch correction")


###############################
####### BATCH CORRECTION ######
###############################


sites = paste(sapply(GeneSymbol(ppe_filtered), function(x)x),";",
                 sapply(Residue(ppe_filtered), function(x)x),
                 sapply(Site(ppe_filtered), function(x)x),
                 ";", sep = "")


design = model.matrix(~ sample_group - 1)
design


data("SPSs")

ctl = which(sites %in% SPSs)
ppe_norm = RUVphospho(ppe_filtered, M = design, k = 2, ctl = ctl)


plotQC(SummarizedExperiment::assay(ppe_norm, "normalised"), grps=sample_group, 
    labels = colnames(ppe), panel="dendrogram")

g1 <- plotQC(SummarizedExperiment::assay(ppe_norm,"normalised"), grps=sample_group, 
    labels = colnames(ppe_filtered), panel = "pca") +
    ggplot2::ggtitle("After batch correction")

library(ggplot2)

pdf(paste0('~/Documents/TSC_Paper/', 'Phospho_omics_PCA.pdf'), width=5, height=4)
g1 + theme_classic()
dev.off()

###################################################
####### Differential Phosphorylation Testing ######
###################################################

# Toss out TSC1 since differential phosphorylation is artifactual

retained_symbols <- symbols
not_TSC1 <- which(retained_symbols != 'TSC1')


ppe <- PhosphoExperiment(assays = list(Quantification = as.matrix(log(abundances_final[not_TSC1,] + 1))),
	Sequence=flanking[not_TSC1],
	GeneSymbol=retained_symbols[not_TSC1],
	Residue=residue[not_TSC1],
	Site=site[not_TSC1])

sites = paste(sapply(GeneSymbol(ppe), function(x)x),";",
                 sapply(Residue(ppe), function(x)x),
                 sapply(Site(ppe), function(x)x),
                 ";", sep = "")


design = model.matrix(~ sample_group - 1)



data("SPSs")

ctl = which(sites %in% SPSs)
ppe_norm = RUVphospho(ppe, M = design, k = 2, ctl = ctl)


plotQC(SummarizedExperiment::assay(ppe_norm, "normalised"), grps=sample_group, 
    labels = colnames(ppe), panel="dendrogram")
plotQC(SummarizedExperiment::assay(ppe_norm,"normalised"), grps=sample_group, 
    labels = colnames(ppe_filtered), panel = "pca") +
    ggplot2::ggtitle("After batch correction")


suppressPackageStartupMessages({
  library(calibrate)
  library(limma)
  library(directPA)
  library(org.Mm.eg.db)
  library(reactome.db)
  library(annotate)
  library(PhosR)
})


data("PhosphoSitePlus")


sites = paste(sapply(GeneSymbol(ppe_norm), function(x)x),";",
                 sapply(Residue(ppe_norm), function(x)x),
                 sapply(Site(ppe_norm), function(x)x),
                 ";", sep = "")

X <- model.matrix(~ sample_group - 1)
fit <- lmFit(SummarizedExperiment::assay(ppe_norm, "normalised"), X)

table.Control <- topTable(eBayes(fit), number=Inf,coef=1)
table.KO <- topTable(eBayes(fit), number=Inf,coef=2)


DE1.RUV <- c(sum(table.Control[,"adj.P.Val"] < 0.05),sum(table.KO[,"adj.P.Val"] < 0.05))

# extract top-ranked phosphosites for each group comparison
contrast.matrix1 <- makeContrasts(sample_groupKO-sample_groupControl, levels=X)  # defining group comparisons

fit1 <- contrasts.fit(fit, contrast.matrix1)
table.KOvsCNTRL <- topTable(eBayes(fit1), number=Inf)

DE2.RUV <- c(sum(table.KOvsCNTRL[,"adj.P.Val"] < 0.05))



library(EnhancedVolcano)


label_plot_a <- paste0(lapply(str_split(rownames(table.KOvsCNTRL),';'),'[[',1),"_",lapply(str_split(rownames(table.KOvsCNTRL),';'),'[[',2))

pdf(paste0('~/Documents/TSC_Paper/', 'Phospho_omics_vocano_without_TSC1.pdf'), width=4, height=8)
EnhancedVolcano(table.KOvsCNTRL, lab = label_plot_a, 
	x = 'logFC', y = 'adj.P.Val',pCutoff = 0.05,
	FCcutoff=0.1, ylim = c(0, 7), xlim=c(-2.5,2.5))
dev.off()

signif <- subset(table.KOvsCNTRL,adj.P.Val<0.1)


#FNIP1;S907;RNSLSILVPHGDKESSDKKN 
#FNIP1;S907;RNSLSILVPHGDKESSDKK 
#FNIP1;S907;RNSLSILVPHGDKE 
#FNIP1;S220;RAFSEQGPLRL
#FLCN;S62;RAHSPAEGASSESSSPGPKK 
#TFE3;S545;RAASDPLLSSVSPAVSKA 
###################################################
############ SIGNALLOME CONSTRUCTION ##############
###################################################

phospho = SummarizedExperiment::assay(ppe_norm, "normalised")
phospho.mean <- meanAbundance(phospho, grps = sample_group)
aov <- matANOVA(mat=phospho, grps=sample_group)
idx <- (aov < 0.05) & (rowSums(phospho.mean > 0.5) > 0)
phospho.reg <- phospho[idx, ,drop = FALSE]
phos.std <- standardise(phospho.reg)
rownames(phos.std) <- paste0(GeneSymbol(ppe_norm), ";", Residue(ppe_norm), Site(ppe_norm), ";")[idx]

phos.seq <- Sequence(ppe_norm)[idx]

library(viridis)
library(pheatmap)
library(RColorBrewer)

phos.matrices <- kinaseSubstrateScore(substrate.list = PhosphoSite.mouse, 
                                    mat = phos.std, seqs = phos.seq, 
                                    numMotif = 5, numSub = 1, verbose = FALSE)

phosScoringMatrices <- phos.matrices
top <- 3
sites <- c()
for (i in seq_len(ncol(phosScoringMatrices$combinedScoreMatrix))) {
        sites <- union(sites, names(sort(phosScoringMatrices$combinedScoreMatrix[, 
            i], decreasing = TRUE)[seq_len(top)]))
}
o <- intersect(colnames(phosScoringMatrices$combinedScoreMatrix), 
        rownames(KinaseFamily))
annotation_col = data.frame(group = KinaseFamily[o, "kinase_group"], 
        family = KinaseFamily[o, "kinase_family"])
rownames(annotation_col) <- o

color_map = colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdYlBu")))(1000)
#color_map = viridis(100,option="magma")
pheatmap(phosScoringMatrices$combinedScoreMatrix[sites, 
            ], color=color_map, annotation_col = annotation_col, cluster_rows = TRUE, 
            cluster_cols = TRUE, fontsize = 7, main = paste("Top", 
                top, "phosphosite(s) for each kinase"))

pdf(paste0('~/Documents/TSC_Paper/', 'Kinase_substrate_relationship_without_TSC1.pdf'), width=8, height=8)
phos.matrices <- kinaseSubstrateScore(substrate.list = PhosphoSite.mouse, 
                                    mat = phos.std, seqs = phos.seq, 
                                    numMotif = 5, numSub = 1, verbose = FALSE)


dev.off()


pdf(paste0('~/Documents/TSC_Paper/', 'Kinase_substrate_relationship_without_TSC1_stringent.pdf'), width=8, height=8)
phos.matrices.stringent <- kinaseSubstrateScore(substrate.list = PhosphoSite.mouse, 
                                    mat = phos.std, seqs = phos.seq, 
                                    numMotif = 10, numSub = 2, verbose = FALSE)
dev.off()

pdf(paste0('~/Documents/TSC_Paper/', 'Kinase_substrate_relationship_without_TSC1_loose.pdf'), width=8, height=8)
phos.matrices.loose <- kinaseSubstrateScore(substrate.list = PhosphoSite.mouse, 
                                    mat = phos.std, seqs = phos.seq, 
                                    numMotif = 3, numSub = 1, verbose = FALSE)
dev.off()

pdf(paste0('~/Documents/TSC_Paper/', 'Kinase_activity_by_sample.pdf'), width=8, height=8)
heatmap(phos.matrices$ksActivityMatrix)
dev.off()




set.seed(1)
phos.predMat <- kinaseSubstratePred(phos.matrices, top=30, inclusion=10, verbose = TRUE) 

#Delete non-unique rows
phos.predMat<-phos.predMat[unique(rownames(phos.predMat)),]


a<-sort(phos.matrices$combinedScoreMatrix[,'AKT1'])
akt_sites <- unlist(c(PhosphoSite.mouse["AKT1"],PhosphoSite.mouse["AKT2"],
	PhosphoSite.mouse["AKT3"],PhosphoSite.mouse["Humphrey.Akt"],
	PhosphoSite.mouse["Yang.Akt"]))
idx <- names(a) %in% akt_sites
a[idx]

#FOXO1 S253





pdf(paste0('~/Documents/TSC_Paper/', 'Kinase_module_circle.pdf'), width=8, height=8)

Signalomes_results <- Signalomes(KSR=phos.matrices, 
                                predMatrix=phos.predMat, 
                                exprsMat=phos.std, module_res = 8,
                                KOI=c('MTOR','AKT1'),verbose = TRUE)
dev.off()



#phos.matrices$motifScoreMatrix["FLCN;S62;",]


### generate palette
my_color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))
kinase_all_color <- my_color_palette(ncol(phos.matrices$combinedScoreMatrix))
names(kinase_all_color) <- colnames(phos.matrices$combinedScoreMatrix)
kinase_signalome_color <- kinase_all_color[colnames(phos.predMat)]

g1 <- plotSignalomeMap(signalomes = Signalomes_results, color = kinase_signalome_color)

pdf(paste0('~/Documents/TSC_Paper/', 'Kinase_bubble.pdf'), width=8, height=4)
g1 + scale_size_continuous(range = c(0, 10))
dev.off()


pdf(paste0('~/Documents/TSC_Paper/', 'Kinase_network.pdf'), width=8, height=8)
plotKinaseNetwork(KSR = phos.matrices, predMatrix = phos.predMat, 
	threshold = 0.9, color = kinase_all_color, type = 'chord')
dev.off()


pdf(paste0('~/Documents/TSC_Paper/', 'Kinase_graph.pdf'), width=7.3, height=6)
plotKinaseNetwork(KSR = phos.matrices, predMatrix = phos.predMat, 
	threshold = 0.9, color = kinase_all_color, type = 'graph')
dev.off()



Signalomes_results$Signalomes$AKT1$exprs
Signalomes_results$Signalomes$AKT1$annotation
Signalomes_results$kinaseSubstrates$AKT1


subset(Signalomes_results$Signalomes$AKT1$annotation,kinase=='AKT1')


Signalomes_results$proteinModules

comparison_table <- table.KOvsCNTRL
new_names <- paste0(unlist(lapply(str_split(rownames(table.KOvsCNTRL),';'),'[[',1)),';',unlist(lapply(str_split(rownames(table.KOvsCNTRL),';'),'[[',2)),';')
AKT_changes <- c()
for (i in rownames(Signalomes_results$Signalomes$AKT1$annotation)){
	idx <- which(i == new_names)
	AKT_changes <- c(AKT_changes,mean(comparison_table$logFC[idx]))
}
names(AKT_changes) <- rownames(Signalomes_results$Signalomes$AKT1$annotation)

MTOR_changes <- c()
for (i in rownames(Signalomes_results$Signalomes$MTOR$annotation)){
	idx <- which(i == new_names)
	MTOR_changes <- c(MTOR_changes,mean(comparison_table$logFC[idx]))
}
names(MTOR_changes) <- rownames(Signalomes_results$Signalomes$MTOR$annotation)

df1 <- data.frame(group = 'AKT',value = AKT_changes)
df2 <- data.frame(group = 'MTOR',value = MTOR_changes)
df <- rbind(df1,df2)

pdf(paste0('~/Documents/TSC_Paper/', 'Kinase_density.pdf'), width=7.5, height=3)
	ggplot(df,aes(value, colour = group, fill = group)) + 
	geom_density(alpha = 0.1, bw = .1) + theme_classic() + xlim(-1.25,1.25)
dev.off()

#### HELPER FUNCTIONS

.kinaseNetwork <- function(predMatrix, KSR, threskinaseNetwork,
    kinase_signalome_color) {

    kinase_cor <- stats::cor(KSR$combinedScoreMatrix)

    cor_kinase_mat <- kinase_cor
    diag(cor_kinase_mat) <- 0
    kinase_network <- lapply(seq_len(ncol(cor_kinase_mat)), function(x)
        names(which(cor_kinase_mat[, x] > threskinaseNetwork)))
    names(kinase_network) <- colnames(cor_kinase_mat)

    cor_kinase_mat <- apply(cor_kinase_mat, 2,
                            function(x) x > threskinaseNetwork)
    cor_kinase_mat[cor_kinase_mat == FALSE] <- 0
    cor_kinase_mat[cor_kinase_mat == TRUE] <- 1

    network <- igraph::graph_from_adjacency_matrix(cor_kinase_mat,
        mode = "undirected", diag = FALSE)

    kinaseNetwork.res <- list(kinaseNetwork = kinase_network,
        kinaseCor = cor_kinase_mat)

    return(kinaseNetwork.res)

}

.phosphositeClusters <- function(KSR, verbose = TRUE) {
    substrate_cor <- stats::cor(t(KSR$combinedScoreMatrix))
    substrate_hclust <- stats::hclust(stats::dist(KSR$combinedScoreMatrix),
        method = "ward.D")
    if (verbose)
        message("calculating optimal number of clusters...")
    res <- lapply(seq(2,10,1), function(x) {
        substrate_clusters <- stats::cutree(substrate_hclust, k = x)
        cor.res <- lapply(seq_len(x), function(y) {
            substrate_cor = substrate_cor[substrate_clusters == y,
                substrate_clusters == y]
            diag(substrate_cor) <- 0
            return(substrate_cor)
        })
        cor.res <- lapply(cor.res, function(x) median(x))
        cor.res <- unlist(cor.res)
        cor.logic <- sum(cor.res >= 0.5) == x
        if (cor.logic) {
            cluster = mean(cor.res)
            names(cluster) = x
            return(cluster)
        }
    })
    if (isTRUE(is.null(unlist(res)))) {
        res <- lapply(seq(2,10,1), function(x) {
            substrate_clusters <- cutree(substrate_hclust, k = x)
            cor.res <- lapply(seq_len(x), function(y) {
                substrate_cor = substrate_cor[substrate_clusters == y,
                    substrate_clusters == y]
                diag(substrate_cor) <- 0
                return(substrate_cor)
            })
            cor.res <- unlist(lapply(cor.res, function(x) median(x)))
            cor.logic <- sum(cor.res >= 0.1) == x
            if (cor.logic) {
                cluster = median(cor.res)
                names(cluster) = x
                return(cluster)
            }
        })
        res <- as.numeric(names(which(unlist(res) == max(unlist(res)))[1]))
    } else {
        res <- as.numeric(names(which(unlist(res) == max(unlist(res)))[1]))
    }
    if (verbose)
        message(paste0("optimal number of clusters = ", res))
    substrate_clusters <- stats::cutree(substrate_hclust, k = res)
    return(substrate_clusters)
}

#' @import circlize
#' @importFrom utils stack
#' @importFrom rlang .data
.phosRsignalome <- function(predMatrix, signalomeCutoff, kinase_signalome_color,
    modules) {

    signalomeKinase <- colnames(predMatrix)

    signalomeSubstrates <- list()
    for (i in seq_len(length(signalomeKinase))) {
        signalomeSubstrates[[i]] =
            mapply(function(x) x[1],
                strsplit(names(which(
                    predMatrix[,signalomeKinase[[i]]] > signalomeCutoff)),
                    ";"))
    }
    names(signalomeSubstrates) <- signalomeKinase

    ############## generate circlize plot
    dftoPlot_signalome <- stack(signalomeSubstrates)
    dftoPlot_signalome$modules <- modules[dftoPlot_signalome$values]
    #adjacencyData <- with(dftoPlot_signalome, table(ind, modules))
    d = table(dftoPlot_signalome$ind, dftoPlot_signalome$modules)
    adjacencyData <- matrix(table(dftoPlot_signalome$ind, 
        dftoPlot_signalome$modules), nrow = nrow(d), ncol = ncol(d))
    rownames(adjacencyData) = rownames(d)
    colnames(adjacencyData) = colnames(d)

    m = sort(as.integer(unique(dftoPlot_signalome$modules)))
    grid.col <- c(kinase_signalome_color, rep("grey", length(unique(m))))
    names(grid.col) <- c(rownames(adjacencyData), as.character(unique(m)))
    adjacencyData = adjacencyData[,!grepl("noModule",colnames(adjacencyData))]
    
    n = length(grid.col)
    circos.clear()
    circos.par(start.degree = 180)
    circos.initialize(factors = "a", xlim = c(0, n))
    chordDiagram(adjacencyData, transparency = 0.2,
                order = c(rownames(adjacencyData),
                            rev(unique(m))),
                grid.col = grid.col, #big.gap = 15,
                annotationTrack = c("name", "grid"), scale = TRUE)
    title("Signalomes")

    return(signalomeSubstrates)
}

#' @importFrom dplyr count %>%
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
.getSignalomes <- function(predMatrix, exprsMat, KOI, resKinaseNetwork,
    signalomeSubstrates, modules, kinaseGroup, verbose = TRUE) {
    
    protein = mapply("[[", strsplit(rownames(exprsMat), ";"),
        MoreArgs = list(1))
    KinaseSubstrateList <- resKinaseNetwork$kinaseNetwork

    ############## annotate kinase-substrate relationship
    kinaseAnnot = annoKinaseSubstrateRelation(predMatrix)

    ############## proportion of regulation by kinase
    signalomeMatrix <- stack(signalomeSubstrates)
    signalomeMatrix$cluster <- modules[signalomeMatrix$values]

    balloon_bycluster <- signalomeMatrix
    balloon_bycluster <- na.omit(balloon_bycluster) %>%
        dplyr::count(.data$cluster, .data$ind)
    balloon_bycluster$ind <- as.factor(balloon_bycluster$ind)
    balloon_bycluster$cluster <- as.factor(balloon_bycluster$cluster)
    balloon_bycluster <- as.data.frame(
        tidyr::pivot_wider(balloon_bycluster, 
            names_from = .data$ind, values_from = .data$n))

    rownames(balloon_bycluster) = balloon_bycluster$cluster
    balloon_bycluster = balloon_bycluster[,-1]

    balloon_bycluster[is.na(balloon_bycluster)] <- 0
    balloon <- do.call(rbind, lapply(seq_len(nrow(balloon_bycluster)),
        function(x) {
            res = mapply(function(y, balloon_bycluster, x) {
                y/sum(balloon_bycluster[x, ]) * 100
            }, balloon_bycluster[x, ],
            MoreArgs = list(balloon_bycluster = balloon_bycluster, x = x))
        }))
    rownames(balloon) = rownames(balloon_bycluster)
    colnames(balloon) = colnames(balloon_bycluster)
    kinaseProportions <- round(balloon, 3)

    ############## generate kinase specific signalomes
    res = generateSignalome(kinaseAnnot, kinaseGroup, predMatrix, KOI,
                            KinaseSubstrateList, kinaseProportions,
                            signalomeSubstrates, exprsMat, protein, modules, 
                            verbose = verbose)
    return(res)
}

annoKinaseSubstrateRelation = function(predMatrix) {
    kinaseAnnot <- lapply(seq_len(nrow(predMatrix)), function(x) {
        scores <- predMatrix[x, ]
        siteName <- rownames(predMatrix)[[x]]
        kinaseTop <- names(which(scores == max(scores)))
        score <- as.numeric(max(scores))
        res <- c(siteName, kinaseTop, score)
        return(res)
    })
    kinaseAnnot <- do.call(rbind, kinaseAnnot)
    rownames(kinaseAnnot) <- kinaseAnnot[, 1]
    kinaseAnnot <- data.frame(kinaseAnnot[, -c(1)], stringsAsFactors = FALSE)
    colnames(kinaseAnnot) <- c("kinase", "score")
    kinaseAnnot$score <- as.numeric(kinaseAnnot$score)
    kinaseAnnot
}

generateSignalome = function(kinaseAnnot, kinaseGroup, predMatrix, KOI,
                            KinaseSubstrateList, kinaseProportions,
                            signalomeSubstrates, exprsMat, protein, modules,
                            verbose = TRUE) {
    annotation <- data.frame(kinase = as.factor(kinaseAnnot$kinase),
        kinaseFamily = as.factor(kinaseGroup[kinaseAnnot$kinase]),
        score = as.numeric(kinaseAnnot$score))
    rownames(annotation) <- rownames(predMatrix)
    
    kinaseProp = kinaseProportions[!grepl("noModule",
                                        rownames(kinaseProportions)),]
    m = modules[!grepl("noModule", modules)]
    
    res <- lapply(KOI, function(x) {
        if (x %in% names(KinaseSubstrateList)) {
            regModule <- which(kinaseProp[, x] > 1)
            kinases <- unique(c(x, KinaseSubstrateList[[x]]))
            substrates <- unique(unlist(
                lapply(kinases, function(x) signalomeSubstrates[[x]])))
            exprs_dat <- lapply(regModule, function(y) {
                exprsMat[protein %in%
                        names(m[m == y]) &
                        protein %in% substrates, ]
            })
            exprs_dat <- do.call(rbind, exprs_dat)
            annotation_dat <- annotation[rownames(annotation) %in%
                    rownames(exprs_dat), ]
            kinaseSignalome <- list(exprs = exprs_dat,
                                    annotation = annotation_dat)
            return(kinaseSignalome)
        } else {
            if (verbose)
                message(paste0(x, " is not found"))
        }
    })
    names(res) <- KOI
    res
}


Signalomes_Custom <- function(KSR, predMatrix, exprsMat, KOI, threskinaseNetwork = 0.9, 
    signalomeCutoff = 0.5, module_res = NULL, filter = FALSE, 
    verbose = TRUE) 
{
    if (!is.null(module_res)) {
        if (module_res < 20) {
            module_res = as.integer(module_res)
        }
        else {
            stop("module resolution should be an integer lower than 20")
        }
    }
    protein_assignment = mapply("[[", strsplit(rownames(KSR$combinedScoreMatrix), 
        ";"), MoreArgs = list(1))
    utils::data("KinaseFamily", envir = environment())
    kinaseGroup <- KinaseFamily[, "kinase_group"]
    names(kinaseGroup) <- KinaseFamily[, "gene_symbol"]
    my_color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, 
        "Accent"))
    kinase_all_color <- my_color_palette(ncol(KSR$combinedScoreMatrix))
    names(kinase_all_color) <- colnames(KSR$combinedScoreMatrix)
    kinase_signalome_color <- kinase_all_color[colnames(predMatrix)]
    my_color_palette_kinaseGroup <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, 
        "Set2"))
    kinaseGroup_color <- my_color_palette_kinaseGroup(length(unique(kinaseGroup)))
    names(kinaseGroup_color) <- unique(kinaseGroup)
    resKinaseNetwork <- .kinaseNetwork(predMatrix, KSR, threskinaseNetwork, 
        kinase_signalome_color)
    substrate_clusters <- .phosphositeClusters(KSR, verbose)
    cluster_assignment <- as.factor(substrate_clusters)
    dat.long <- data.frame(table(cluster_assignment, protein_assignment))
    dftoHeatmap <- tidyr::pivot_wider(dat.long, names_from = .data$protein_assignment, 
        values_from = .data$Freq)[, -1]
    dftoHeatmap[is.na(dftoHeatmap)] <- 0
    dftoHeatmap[dftoHeatmap > 0] <- 1
    hclust_res <- stats::hclust(stats::dist(t(dftoHeatmap)), 
        method = "ward.D")
    tree_height <- as.numeric(names(table(hclust_res$height)))
    branching <- as.numeric(table(hclust_res$height))
    tree_height_calc = unlist(lapply(seq(2, length(tree_height), 
        1), function(x) {
        h <- tree_height[[x]]
        m <- stats::cutree(hclust_res, h = h)
        return(length(table(m)))
    }))
    if (!is.null(module_res)) {
        hcutree = which(tree_height_calc <= module_res) + 1
    }
    else {
        hcutree <- min(tree_height[tree_height > 0])
    }
    modules <- stats::cutree(hclust_res, h = tree_height[[hcutree[[1]]]])
    if (filter) {
        filter_modules = modules %in% which(table(modules) < 
            10)
        modules[filter_modules] = "noModule"
    }
    signalomeSubstrates <- .phosRsignalome(predMatrix, signalomeCutoff, 
        kinase_signalome_color, modules)
    signalomes_of_KOI <- .getSignalomes(predMatrix, exprsMat, 
        KOI, resKinaseNetwork, signalomeSubstrates, modules, 
        kinaseGroup, verbose = verbose)
    signalome_res <- list(Signalomes = signalomes_of_KOI, proteinModules = modules, 
        kinaseSubstrates = signalomeSubstrates)
    return(signalome_res)
}

















o <- rownames(signif)

#Tc <- cbind(table.Control[o,"logFC"],table.KO[o,"logFC"])
Tc <- signif[o,c("logFC")]
#rownames(Tc) <- sites[match(o, rownames(ppe_norm))]
#rownames(Tc) <- gsub("(.*)(;[A-Z])([0-9]+)(;)", "\\1;\\3;", rownames(Tc))
#colnames(Tc) <- c('Cntrl','KO')
names(Tc) <- sites[match(o, rownames(ppe_norm))]
names(Tc) <- gsub("(.*)(;[A-Z])([0-9]+)(;)", "\\1;\\3;", names(Tc))

Tc.gene <- phosCollapse(Tc, id=gsub(";.+", "", names(Tc)), stat=abs(Tc), by = "max")
colnames(Tc.gene) <- 'logFC'
up <- rownames(Tc.gene)[Tc.gene[,'logFC'] >0]
down <- rownames(Tc.gene)[Tc.gene[,'logFC'] <0]


#Tc.gene <- phosCollapse(Tc, id=gsub(";.+", "", rownames(Tc)), 
                        #stat=apply(abs(Tc), 1, max), by = "max")


#geneSet <- names(sort(Tc.gene[,1], 
#                        decreasing = TRUE))[seq(round(nrow(Tc.gene) * 0.1))]

geneset <- names(sort(Tc.gene[,1], decreasing = TRUE))[seq(round(nrow(Tc.gene) * 0.1))]


pathways = as.list(reactomePATHID2EXTID)

path_names = as.list(reactomePATHID2NAME)
name_id = match(names(pathways), names(path_names))
names(pathways) = unlist(path_names)[name_id]

pathways = pathways[which(grepl("Mus musculus", names(pathways), ignore.case = TRUE))]

pathways = lapply(pathways, function(path) {
    gene_name = unname(getSYMBOL(path, data = "org.Mm.eg.db"))
    toupper(unique(gene_name))
})


path1 <- pathwayOverrepresent(geneSet, annotation=pathways, 
                                universe = names(Tc.gene), alter = "greater")
path2 <- pathwayRankBasedEnrichment(Tc.gene[,1], 
                                    annotation=pathways, 
                                    alter = "greater")

lp1 <- -log10(as.numeric(path2[names(pathways),1]))
lp2 <- -log10(as.numeric(path1[names(pathways),1]))
plot(lp1, lp2, ylab="Overrepresentation (-log10 pvalue)", xlab="Rank-based enrichment (-log10 pvalue)", main="Comparison of 1D pathway analyses", xlim = c(0, 10))


sel <- which(lp1 > 1.5 & lp2 > 0.9)
textxy(lp1[sel], lp2[sel], gsub("_", " ", gsub("REACTOME_", "", names(pathways)))[sel])

Tc.temp <- Tc
rownames(Tc.temp) <- toupper(rownames(Tc.temp))



sort(Tc.gene[,1],decreasing = TRUE)

path1 <- pathwayOverrepresent(rownames(Tc), annotation=lapply(PhosphoSite.mouse, function(x){gsub(";[STY]", ";", x)}), 
                                universe = unique(gene_symbols), alter = "greater")

path1$pathways[1:5,]

directPA(Tc.gene, direction=0, annotation=lapply(PhosphoSite.mouse, function(x){gsub(";[STY]", ";", x)}))






