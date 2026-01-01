python3 estimateTissueBoundary.py /Users/ikuz/Downloads/SeqScope/script/spatialcoordinates.txt /Users/ikuz/Downloads/SeqScope/script/HDMI_SeqScope_2nd.txt 750 ../tissueboundary


Rscript getSimpleGrid.R "MiSeq" "./STAR_outSolo.out/GeneFull/raw/" "/Users/ikuz/Downloads/SeqScope/script/spatialcoordinates.txt" "2107,2108,2109,2110,2111,2112,2113,2114,2115,2116,2117,2118,2119" "1300" "1300" "300" "/Users/ikuz/Downloads/SeqScope/script/binning" "./collapse.cpp"

Rscript getSimpleGrid.R "HiSeq" "./STAR_outSolo.out/GeneFull/raw/" "/Users/ikuz/Downloads/SeqScope/script/spatialcoordinates.txt" "2117" "100" "100" "300" "/Users/ikuz/Downloads/SeqScope/script/binning" "./collapse.cpp"


Rscript getSlidingGrid.R "MiSeq" "./STAR_outSolo.out/GeneFull/raw/" "/Users/ikuz/Downloads/SeqScope/script/spatialcoordinates.txt" "2117" "100" "100" "300" "/Users/ikuz/Downloads/SeqScope/script/binning" 60 0 30000 0 30000 "./collapse.cpp"



FC = /opt/homebrew/Cellar/gcc/14.2.0_1/bin/gfortran
F77 = /opt/homebrew/Cellar/gcc/14.2.0_1/bin/gfortran
FLIBS = -L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/14


## Function
seurat_to_spe <- function(seu, sample_id, img_id) {
    ## Convert to SCE
    sce <- Seurat::as.SingleCellExperiment(seu)
    
    ## Extract spatial coordinates
    spatialCoords <- as.matrix(
        seu@images[[img_id]]@coordinates[, c("Y_expand", "X_expand")])
    
    ## Extract and process image data
    img <- SpatialExperiment::SpatialImage(
        x = as.raster(seu@images[[img_id]]@image))
    
    imgData <- DataFrame(
        sample_id = sample_id,
        image_id = img_id,
        data = I(list(img)),
        scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
    
    # Convert to SpatialExperiment
    spe <- SpatialExperiment(
        assays = assays(sce),
        rowData = rowData(sce),
        colData = colData(sce),
        metadata = metadata(sce),
        reducedDims = reducedDims(sce),
        altExps = altExps(sce),
        sample_id = sample_id,
        spatialCoords = spatialCoords,
        imgData = imgData
    )
    # indicate all spots are on the tissue
    spe$in_tissue <- 1
    spe$sample_id <- sample_id
    # Return Spatial Experiment object
    spe
}



library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(celda)


zone <- read.csv('~/Downloads/SeqScope/Zonation.csv')

prior <- readRDS("/Users/ikuz/Downloads/SeqScope/script/prior/Liver_normal_10um_annotated.RDS")
prior$names <- prior@active.ident


prior_td <- readRDS("/Users/ikuz/Downloads/SeqScope/script/prior/Liver_TD_10um_annotated.RDS")
prior_td$names <- prior_td@active.ident


data <- readRDS("/Users/ikuz/Downloads/SeqScope/script/binning/10um_SimpleSqureGrids.RDS")
data[["Spatial"]] <- JoinLayers(data[["Spatial"]])
counts <- GetAssayData(data, layer="counts", assay="Spatial") 
spot_counts <- colSums(counts)

genes.percent.expression <- rowMeans(counts[,spot_counts > 100]>0 )*100   
genes.filter <- names(genes.percent.expression[genes.percent.expression>0])  #select genes expressed in at least 5% of cells
counts.sub <- counts[genes.filter,]

genes <- read.csv("/Users/ikuz/Downloads/SeqScope/script/STAR_outSolo.out/GeneFull/raw/features.tsv",sep='\t',header = F)
rownames(genes) <- genes$V1
rename_genes <- genes$V2[match(rownames(counts.sub),genes$V1)]
rownames(counts.sub) <- rename_genes

idx <- intersect(rownames(counts.sub),rownames(prior))


counts.sub <- counts.sub[rownames(counts.sub) %in% idx,]

dup <- rownames(counts.sub)[duplicated(rownames(counts.sub))]
for (x in dup){
	idx <- which(rownames(counts.sub) == x)
	counts.sub[idx[1],] <- colSums(counts.sub[idx,])
	counts.sub <- counts.sub[-idx[-1],]
}

data_new <- CreateSeuratObject(counts=counts.sub,assay="Spatial")
data_new@meta.data <- data@meta.data
data_new$image <- data$image
data_new$nCount_Spatial <- colSums(counts.sub)
data_new$nFeature_Spatial <- colSums(counts.sub>0)

VlnPlot(data_new, features = "nFeature_Spatial", pt.size = 0.1)
SpatialFeaturePlot(data_new, features = "nCount_Spatial",shape=22,pt.size.factor=3.5)

data_new <- subset(data_new,nFeature_Spatial > 50)


anchors <- FindTransferAnchors(reference = prior, query = data_new, dims = 1:30, reference.reduction = "pca", normalization.method = 'SCT')

predictions <- TransferData(anchorset = anchors, refdata = prior$names, dims = 1:30)

data_new <- AddMetaData(data_new,metadata = predictions)


table(data_new$predicted.id)
table(prior$names)

prior <- RunUMAP(prior, dims = 1:30, reduction = "pca", return.model = TRUE)

data_new <- MapQuery(anchorset = anchors, reference = prior, query = data_new,
    refdata = list(celltype = "names"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(prior, reduction = "umap", group.by = "names", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(data_new, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2


data_new <- SetIdent(data_new, value = "predicted.celltype")


SpatialDimPlot(prior, label = TRUE, label.size = 3,pt.size.factor=3.5)
SpatialDimPlot(data_new, label = TRUE, label.size = 3,pt.size.factor=3.5)

data_new <- SCTransform(data_new, assay = "Spatial", verbose = FALSE, variable.features.n = 3000)

cat(rownames(subset(FindMarkers(data_new,only.pos=T,ident.1='Macrophage'),p_val_adj<0.05)),sep='\n')
cat(rownames(subset(FindMarkers(data_new,only.pos=T,ident.1='ENDO'),p_val_adj<0.05)),sep='\n')
cat(rownames(subset(FindMarkers(data_new,only.pos=T,ident.1='RBC'),p_val_adj<0.05)),sep='\n')
cat(rownames(subset(FindMarkers(data_new,only.pos=T,ident.1='Hep_Z1'),p_val_adj<0.05)),sep='\n')


pseudo_data <- AggregateExpression(data_new, assays = "Spatial", return.seurat = T, group.by='predicted.celltype')
pseudo_prior <- AggregateExpression(prior, assays = "Spatial", return.seurat = T, group.by='names')

### Clean up


data_new$orig.ident <- 'new'
data_new$names <- data_new$predicted.id
prior$orig.ident <- 'old'

DefaultAssay(data_new) <- 'Spatial'

sce.raw <- as.SingleCellExperiment(CreateSeuratObject(counts=counts.sub,assay="Spatial"))
sce <- as.SingleCellExperiment(data_new)
z <- data_new$names
sce <- decontX(sce, background = sce.raw, z=z)

markers <- list(Endo = c("Aqp1","Kdr"), Hep = c("Cyp2f2", "B2m", "Apoe"), Macro = c("Cd74","Cd5l","H2-Ab1"))

cellTypeMappings <- list(Hep_Z1 = 1, Endo = 7, Macro = 8)

umap <- reducedDim(sce, "REF.UMAP")
plotDimReduceCluster(x = sce$decontX_clusters,
    dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXMarkerExpression(sce,
    markers = markers[['Hep']],
    groupClusters = cellTypeMappings,
    ncol = 2)


data_clean <- CreateSeuratObject(decontXcounts(sce),assay="Spatial")
data_clean@meta.data <- data_new@meta.data
data_clean$image <- data_new$image
data_clean$nCount_Spatial <- colSums(decontXcounts(sce))
data_clean$nFeature_Spatial <- colSums(decontXcounts(sce)>0)


DefaultAssay(prior) <- 'Spatial'

sce <- as.SingleCellExperiment(prior)
z <- prior$names
sce <- decontX(sce, z=z)

markers <- list(Endo = c("Aqp1","Kdr"), Hep = c("Cyp2f2", "B2m", "Apoe"), Macro = c("Cd74","Cd5l","H2-Ab1"))

cellTypeMappings <- list(Hep_Z1 = 1, Endo = 5, Macro = 6)

umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(x = sce$decontX_clusters,
    dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXMarkerExpression(sce,
    markers = markers[['Macro']],
    groupClusters = cellTypeMappings,
    ncol = 2)


prior_clean <- CreateSeuratObject(decontXcounts(sce),assay="Spatial")
prior_clean@meta.data <- prior@meta.data
prior_clean$image <- prior$image
prior_clean$nCount_Spatial <- colSums(decontXcounts(sce))
prior_clean$nFeature_Spatial <- colSums(decontXcounts(sce)>0)



DefaultAssay(prior_td) <- 'Spatial'

sce <- as.SingleCellExperiment(prior_td)
z <- prior_td$names
sce <- decontX(sce, z=z)

markers <- list(Endo = c("Aqp1","Kdr"), Hep = c("Cyp2f2", "B2m", "Apoe"), Macro = c("Cd74","Cd5l","H2-Ab1"))

cellTypeMappings <- list(Hep_Z1 = 2, Endo = 11, Macro = 6)

umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(x = sce$decontX_clusters,
    dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXMarkerExpression(sce,
    markers = markers[['Macro']],
    groupClusters = cellTypeMappings,
    ncol = 2)


prior_td_clean <- CreateSeuratObject(decontXcounts(sce),assay="Spatial")
prior_td_clean@meta.data <- prior_td@meta.data
prior_td_clean$image <- prior_td$image
prior_td_clean$nCount_Spatial <- colSums(decontXcounts(sce))
prior_td_clean$nFeature_Spatial <- colSums(decontXcounts(sce)>0)


saveRDS(data_clean,'~/Downloads/SeqScope/data_clean.rds')
saveRDS(prior_clean,'~/Downloads/SeqScope/prior_clean.rds')
saveRDS(prior_td_clean,'~/Downloads/SeqScope/prior_td_clean.rds')



### Integrate


combo <- merge(data_clean,prior_clean)
combo <- SCTransform(combo, assay = "Spatial", verbose = FALSE, variable.features.n = 3000)
combo <- RunPCA(combo)
combo <- RunUMAP(combo, dims = 1:30)

combo <- IntegrateLayers(object = combo, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE, normalization.method = "SCT")


# re-join layers after integration
combo[["Spatial"]] <- JoinLayers(combo[["Spatial"]])

combo <- RunUMAP(combo, dims = 1:30, reduction = "integrated.cca")
DimPlot(combo,reduction='umap',group.by='names')


aggregate<- AggregateExpression(combo, group.by = c("orig.ident", "names"), return.seurat = TRUE)

idx <- intersect(intersect(rownames(data_clean),rownames(prior_clean)),rownames(combo))

combo$celltype.ident <- paste(combo$names, combo$orig.ident, sep = "_")
Idents(combo) <- "celltype.ident"

slot(object = combo@assays$SCT@SCTModel.list[[2]], name="umi.assay")<-"Spatial"
combo <- PrepSCTFindMarkers(combo)

hep_marks <- FindMarkers(combo, ident.1 = c("Hep_Z1_new","Hep_Z3_new","Hep_Nuc_new","Hep_Injured_new","Hep_MT_new","Hep_Eef1a1_new","Hep_Mup10_new","Hep_Z1_old","Hep_Z3_old","Hep_Nuc_old","Hep_Injured_old","Hep_MT_old","Hep_Eef1a1_old","Hep_Mup10_old"),only.pos=T)
hep_marks <- rownames(subset(hep_marks,p_val_adj<0.05))

hep_marks_new <- FindMarkers(combo, ident.1 = c("Hep_Z1_new","Hep_Z3_new","Hep_Nuc_new","Hep_Injured_new","Hep_MT_new","Hep_Eef1a1_new","Hep_Mup10_new"),only.pos=T)
hep_marks_new <- rownames(subset(hep_marks_new,p_val_adj<0.05))

hep_marks_old <- FindMarkers(combo, ident.1 = c("Hep_Z1_old","Hep_Z3_old","Hep_Nuc_old","Hep_Injured_old","Hep_MT_old","Hep_Eef1a1_old","Hep_Mup10_old"),only.pos=T)
hep_marks_old <- rownames(subset(hep_marks_old,p_val_adj<0.05))

hep_marks_all <- union(union(hep_marks,hep_marks_new),hep_marks_old)

idx_non_hep <- setdiff(idx,hep_marks_all)

p1 <- CellScatter(aggregate, "old_Macrophage", "new_Macrophage",features=idx_non_hep) + geom_text(label=idx_non_hep) #cd68 Fgb is up
p1 <- CellScatter(aggregate, "old_ENDO", "new_ENDO",features=idx_non_hep) + geom_text(label=idx_non_hep) #FLT4, Gpihbp1, Esam



saveRDS(combo,'~/Downloads/SeqScope/clean_normal_HCC_integrated.rds')


macro.resp <- FindMarkers(combo, ident.1 = "Macrophage_new", ident.2 = "Macrophage_old", verbose = FALSE,recorrect_umi=T,features=idx_non_hep,only.pos=F)
macro.resp.up <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
macro.resp.down <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
cat(rownames(subset(macro.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(macro.resp.down,p_val_adj<0.05)),sep='\n')

ec.resp <- FindMarkers(combo, ident.1 = "ENDO_new", ident.2 = "ENDO_old", verbose = FALSE,recorrect_umi=T,features=idx_non_hep,only.pos=F)
ec.resp.up <- subset(ec.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
ec.resp.down <- subset(ec.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
cat(rownames(subset(ec.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(ec.resp.down,p_val_adj<0.05)),sep='\n')

hsc.resp <- FindMarkers(combo, ident.1 = "HSC_new", ident.2 = "HSC_old", verbose = FALSE,recorrect_umi=T,features=idx_non_hep,only.pos=F)
hsc.resp.up <- subset(hsc.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
hsc.resp.down <- subset(hsc.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
cat(rownames(subset(hsc.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(hsc.resp.down,p_val_adj<0.05)),sep='\n')

hepz1.resp <- FindMarkers(combo, ident.1 = "Hep_Z1_new", ident.2 = "Hep_Z1_old", verbose = FALSE,recorrect_umi=T,features=idx,only.pos=F)
hepz1.resp.up <- subset(hepz1.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
hepz1.resp.down <- subset(hepz1.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
cat(rownames(subset(hepz1.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(hepz1.resp.down,p_val_adj<0.05)),sep='\n')

hepz1.resp <- FindMarkers(combo, ident.1 = c("Hep_Z1_new","Hep_Z3_new","Hep_Nuc_new","Hep_Injured_new","Hep_MT_new","Hep_Eef1a1_new","Hep_Mup10_new"), ident.2 = c("Hep_Z1_old","Hep_Z3_old","Hep_Nuc_old","Hep_Injured_old","Hep_MT_old","Hep_Eef1a1_old","Hep_Mup10_old"), verbose = FALSE,recorrect_umi=T,features=idx,only.pos=F)
hepz1.resp.up <- subset(hepz1.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
hepz1.resp.down <- subset(hepz1.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
cat(rownames(subset(hepz1.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(hepz1.resp.down,p_val_adj<0.05)),sep='\n')



# Integrate with TD data

combo_all <- merge(combo,prior_td_clean)
combo_all <- SCTransform(combo_all, assay = "Spatial", verbose = FALSE, variable.features.n = 3000)
combo_all <- RunPCA(combo_all)
combo_all <- RunUMAP(combo_all, dims = 1:30)

combo_all <- IntegrateLayers(object = combo_all, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE, normalization.method = "SCT")

# re-join layers after integration
combo_all[["Spatial"]] <- JoinLayers(combo_all[["Spatial"]])

combo_all$orig.ident[combo_all$orig.ident == 'new'] = 'HCC'
combo_all$orig.ident[combo_all$orig.ident == 'old'] = 'normal'
combo_all$orig.ident[combo_all$orig.ident == '2116'] = 'td'
combo_all$orig.ident[combo_all$orig.ident == '2117'] = 'td'
combo_all$orig.ident[combo_all$orig.ident == '2118'] = 'td'
combo_all$orig.ident[combo_all$orig.ident == '2119'] = 'td'

combo_all$new_names <- combo_all$names
combo_all$new_names[combo_all$new_names == 'Hep_Z1'] = 'Hep_PP'
combo_all$new_names[combo_all$new_names == 'Hep_Z3'] = 'Hep_PC'
combo_all$new_names[combo_all$new_names == 'HSC-A'] = 'HSC'
combo_all$new_names[combo_all$new_names == 'HSC-N'] = 'HSC'
combo_all$new_names[combo_all$new_names == 'Macrophage_Kupffer'] = 'Macrophage'
combo_all$new_names[combo_all$new_names == 'Macrophage_Inflamed'] = 'Macrophage'

combo_all$celltype.ident <- paste(combo_all$new_names, combo_all$orig.ident, sep = "_")
Idents(combo_all) <- "celltype.ident"

slot(object = combo_all@assays$SCT@SCTModel.list[[2]], name="umi.assay")<-"Spatial"
combo_all <- PrepSCTFindMarkers(combo_all)


macro.resp <- FindMarkers(combo_all, ident.1 = "Macrophage_td", ident.2 = "Macrophage_normal", verbose = FALSE,recorrect_umi=T,features=idx_non_hep,only.pos=F)
macro.resp.up <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
macro.resp.down <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
cat(rownames(subset(macro.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(macro.resp.down,p_val_adj<0.05)),sep='\n')

hepz1.resp <- FindMarkers(combo_all, ident.1 = "Hep_PP_td", ident.2 = "Hep_PP_HCC", verbose = FALSE,recorrect_umi=T,features=idx,only.pos=F)
hepz1.resp.up <- subset(hepz1.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
hepz1.resp.down <- subset(hepz1.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
cat(rownames(subset(hepz1.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(hepz1.resp.down,p_val_adj<0.05)),sep='\n')

hepz1.resp <- FindMarkers(combo_all, ident.1 = c("Hep_PP_td","Hep_PC_td","Hep_Nuc_td","Hep_Injured_td","Hep_MT_td","Hep_Eef1a1_td"), ident.2 = c("Hep_PC_normal","Hep_PP_normal","Hep_Nuc_normal","Hep_Injured_normal","Hep_MT_normal","Hep_Eef1a1_normal","Hep_Mup10_normal"), verbose = FALSE,recorrect_umi=T,features=idx,only.pos=F)
hepz1.resp.up <- subset(hepz1.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
hepz1.resp.down <- subset(hepz1.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
cat(rownames(subset(hepz1.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(hepz1.resp.down,p_val_adj<0.05)),sep='\n')


DefaultAssay(prior_td) <- 'SCT'
prior_td <- RunUMAP(prior_td, dims = 1:30, reduction = "pca", return.model = TRUE)

anchors <- FindTransferAnchors(reference = prior_td, query = combo_all, dims = 1:30, reference.reduction = "pca", normalization.method = 'SCT')

combo_all <- MapQuery(anchorset = anchors, reference = prior_td, query = combo_all,
    refdata = list(celltype = "names"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(prior_td, reduction = "umap", group.by = "names", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(combo_all, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

#saveRDS(combo_all,'~/Downloads/SeqScope/clean_normal_td_HCC_integrated.rds')


table(combo_all$predicted.celltype,combo_all$orig.ident)

clust_proj_norm <- subset(combo_all,orig.ident == 'normal')$predicted.celltype
table(prior$names,clust_proj_norm)
#Hep_Z1 is mapped to Hep_MT and Hep_PP


clust_proj_HCC <- subset(combo_all,orig.ident == 'HCC')$predicted.celltype

data_clean$names <- clust_proj_HCC
data_clean$predicted.celltype.score <- subset(combo_all,orig.ident == 'HCC')$predicted.celltype.score


SpatialDimPlot(subset(data_clean,names=='Hep_PP'), label = TRUE, label.size = 3,pt.size.factor=2)
VlnPlot(data_clean,'predicted.celltype.score',group.by='names',pt.size=0)
hist(subset(data_clean,names=='Hep_PP')$nCount_Spatial)
SpatialDimPlot(subset(data_clean,names %in% c('Hep_PP','Hep_PC','Hep_MT')), label = TRUE, group.by='names',label.size = 3,pt.size.factor=1.5)
SpatialFeaturePlot(data_clean,'predicted.celltype.score',pt.size.factor=1.5)


combo_all$celltype.ident <- paste(combo_all$predicted.celltype, combo_all$orig.ident, sep = "_")
Idents(combo_all) <- "celltype.ident"


macro.resp <- FindMarkers(combo_all, ident.1 = c("Macrophage_Kupffer_HCC"), ident.2 = c("Macrophage_Kupffer_normal"), verbose = FALSE,recorrect_umi=T,features=idx_non_hep,only.pos=F)
macro.resp.up <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
macro.resp.down <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
cat(rownames(subset(macro.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(macro.resp.down,p_val_adj<0.05)),sep='\n')
#Kuppfer upregulate CD68

macro.resp <- FindMarkers(combo_all, ident.1 = c("HSC-N_HCC"), ident.2 = c("HSC-N_normal"), verbose = FALSE,recorrect_umi=T,features=idx_non_hep,only.pos=F)
macro.resp.up <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
macro.resp.down <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
cat(rownames(subset(macro.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(macro.resp.down,p_val_adj<0.05)),sep='\n')
#Stellate cells express SPARC <- https://www.nature.com/articles/s41598-017-18981-9


macro.resp <- FindMarkers(combo_all, ident.1 = c("HSC-A_HCC"), ident.2 = c("HSC-A_normal"), verbose = FALSE,recorrect_umi=T,features=idx_non_hep,only.pos=F)
#macro.resp.up <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
#macro.resp.down <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
macro.resp.up <- subset(macro.resp,avg_log2FC > 0)
macro.resp.down <- subset(macro.resp,avg_log2FC < 0)
cat(rownames(subset(macro.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(macro.resp.down,p_val_adj<0.05)),sep='\n')


macro.resp <- FindMarkers(combo_all, ident.1 = "ENDO_HCC", ident.2 = "ENDO_normal", verbose = FALSE,recorrect_umi=T,features=idx_non_hep,only.pos=F)
macro.resp.up <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
macro.resp.down <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
cat(rownames(subset(macro.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(macro.resp.down,p_val_adj<0.05)),sep='\n')

macro.resp <- FindMarkers(combo_all, ident.1 = "Hep_PP_HCC", ident.2 = "Hep_PP_normal", verbose = FALSE,recorrect_umi=T,only.pos=F)
macro.resp.up <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC > 0)
macro.resp.down <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01 & avg_log2FC < 0)
cat(rownames(subset(macro.resp.up,p_val_adj<0.05)),sep='\n')
cat(rownames(subset(macro.resp.down,p_val_adj<0.05)),sep='\n')


x <- data_clean@images$image@coordinates$x
y <- data_clean@images$image@coordinates$y
cell.names <- rownames(data_clean@images$image@coordinates)

#Niches

coords <- as.matrix(cbind(x,y))
rownames(coords) <- cell.names
neighbors    <- FindNeighbors(coords, k.param = 30)
neighbors$nn <- neighbors$nn[cell.names,cell.names]

# Continuouing on the BuildNicheAssay function
# build cell x cell type matrix
ct.mtx <- matrix(
    data = 0,
    nrow = length(cell.names),
    ncol = length(unlist(unique(data_clean[['names']])))
)
rownames(ct.mtx) <- cell.names
colnames(ct.mtx) <- unique(unlist(data_clean[['names']]))

cts <- data_clean[['names']]
for (i in 1:length(cell.names)) {
    ct <- as.character(cts[cell.names[[i]], ])
    ct.mtx[cell.names[[i]], ct] <- 1
}
  
# create niche assay
sum.mtx <- as.matrix(neighbors$nn %*% ct.mtx)
niche.assay <- CreateAssayObject(counts = t(sum.mtx))
data_clean[['niche']] <- niche.assay
DefaultAssay(data_clean) <- 'niche'
  
  
# cluster niches assay
data_clean <- ScaleData(data_clean)
results <- kmeans(
    x = t(data_clean[['niche']]@scale.data),
    centers = 5,
    nstart = 30
)
data_clean$niches <- results[["cluster"]]
  
SpatialDimPlot(data_clean,group.by='niches', label = TRUE, label.size = 3,pt.size.factor=2)
SpatialDimPlot(subset(data_clean,niches=='3'),group.by='niches', label = TRUE, label.size = 3,pt.size.factor=2)#injured heps



table(data_clean$niches,data_clean$names) / rowSums(table(data_clean$niches,data_clean$names))*100


#Preprocessing Complete

pdf(paste0('~/Downloads/SeqScope/', 'nCount_Spatial.pdf'), width=2, height=13)
SpatialFeaturePlot(data_clean, features = "nCount_Spatial",shape=22,pt.size.factor=3.5) + theme(aspect.ratio = 13)
dev.off()

pdf(paste0('~/Downloads/SeqScope/', 'CellType.pdf'), width=5, height=13)
SpatialDimPlot(data_clean, group.by = "names",shape=22,pt.size.factor=3.5,label=F) + theme(aspect.ratio = 13)
dev.off()

pdf(paste0('~/Downloads/SeqScope/', 'Niches.pdf'), width=5, height=13)
SpatialDimPlot(data_clean, group.by = "niches",shape=22,pt.size.factor=3.5,label=F) + theme(aspect.ratio = 13)
dev.off()

pdf(paste0('~/Downloads/SeqScope/', 'Injured_niche.pdf'), width=5, height=13)
SpatialDimPlot(subset(data_clean,niches=='3' | names == 'Macrophage_Inflamed'), group.by = "names",shape=22,pt.size.factor=3.5,label=F) + theme(aspect.ratio = 13)
dev.off()



data <- readRDS("/Users/ikuz/Downloads/SeqScope/script/binning/SlidingSquareGrids.RDS")
counts <- GetAssayData(data, layer="counts", assay="RNA") 
spot_counts <- colSums(counts)

genes.percent.expression <- rowMeans(counts[,spot_counts > 100]>0 )*100   
genes.filter <- names(genes.percent.expression[genes.percent.expression>0])  #select genes expressed in at least 5% of cells
counts.sub <- counts[genes.filter,]

dup <- rownames(counts.sub)[duplicated(rownames(counts.sub))]
for (x in dup){
	idx <- which(rownames(counts.sub) == x)
	counts.sub[idx[1],] <- colSums(counts.sub[idx,])
	counts.sub <- counts.sub[-idx[-1],]
}



# High resolution mapping
data_sliding <- CreateSeuratObject(counts=counts.sub,assay="Spatial")
data_sliding@meta.data <- data@meta.data
data_sliding$image <- data$image
data_sliding$nCount_Spatial <- colSums(counts.sub)
data_sliding$nFeature_Spatial <- colSums(counts.sub>0)

data_sliding <- subset(data_sliding,nFeature_Spatial>50)


SpatialFeaturePlot(data_sliding, features = "nCount_Spatial",shape=22,pt.size.factor=1.75) 


anchors <- FindTransferAnchors(reference = prior_td, query = data_sliding, dims = 1:30, reference.reduction = "pca", normalization.method = 'SCT')

data_sliding <- MapQuery(anchorset = anchors, reference = prior_td, query = data_sliding,
    refdata = list(celltype = "names"), reference.reduction = "pca", reduction.model = "umap")

data_sliding$names <- data_sliding$predicted.celltype


p1 <- DimPlot(prior_td, reduction = "umap", group.by = "names", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(data_sliding, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

SpatialDimPlot(subset(data_sliding,names %in% c('Hep_Injured','Macrophage_Inflamed')), group.by='names',shape=22,pt.size.factor=1.75) 

x <- data_sliding@images$image@coordinates$x
y <- data_sliding@images$image@coordinates$y
cell.names <- rownames(data_sliding@images$image@coordinates)

#Niches

coords <- as.matrix(cbind(x,y))
rownames(coords) <- cell.names
neighbors    <- FindNeighbors(coords, k.param = 100)
neighbors$nn <- neighbors$nn[cell.names,cell.names]

# Continuouing on the BuildNicheAssay function
# build cell x cell type matrix
ct.mtx <- matrix(
    data = 0,
    nrow = length(cell.names),
    ncol = length(unlist(unique(data_sliding[['names']])))
)
rownames(ct.mtx) <- cell.names
colnames(ct.mtx) <- unique(unlist(data_sliding[['names']]))

cts <- data_sliding[['names']]
for (i in 1:length(cell.names)) {
    ct <- as.character(cts[cell.names[[i]], ])
    ct.mtx[cell.names[[i]], ct] <- 1
}
  
# create niche assay
sum.mtx <- as.matrix(neighbors$nn %*% ct.mtx)
niche.assay <- CreateAssayObject(counts = t(sum.mtx))
data_sliding[['niche']] <- niche.assay
DefaultAssay(data_sliding) <- 'niche'
  
  
# cluster niches assay
data_sliding <- ScaleData(data_sliding)
results <- kmeans(
    x = t(data_sliding[['niche']]@scale.data),
    centers = 10,
    nstart = 30
)
data_sliding$niches <- results[["cluster"]]
  
table(data_sliding$niches,data_sliding$names) / rowSums(table(data_sliding$niches,data_sliding$names))*100
#2 Inflamed macs
#3 HEP MT
#4 Injured Heps
#5 Kuppfer macs
#6 Hep PP
#9 Endo

SpatialDimPlot(data_sliding,group.by='niches', label = TRUE, label.size = 3,pt.size.factor=1.75)
SpatialDimPlot(subset(data_sliding,niches=='5'),group.by='niches', label = TRUE, label.size = 3,pt.size.factor=1.75)

SpatialDimPlot(subset(data_sliding,niches %in% c('5','4')),group.by='niches', label = TRUE, label.size = 3,pt.size.factor=1.75)


results <- kmeans(
    x = t(data_sliding[['niche']]@scale.data),
    centers = 3,
    nstart = 30
)
data_sliding$niches <- results[["cluster"]]
table(data_sliding$niches,data_sliding$names) / rowSums(table(data_sliding$niches,data_sliding$names))*100

#1 Hep_Injured, Hep_PP, HPC, HSC-N, Macrophage_Inflamed
#2
#3
#5 HSC-A, ENDO
SpatialDimPlot(subset(data_sliding,niches %in% c('5')),group.by='names', label = F, label.size = 3,pt.size.factor=1.75)

SpatialDimPlot(subset(data_sliding,niches %in% c('1') & names %in% c('Macrophage_Inflamed','HPC','Hep_Injured')),group.by='names', label = F, label.size = 3,pt.size.factor=1.75)

SpatialDimPlot(subset(data_sliding,niches %in% c('5') & names %in% c('HSC-A','ENDO')),group.by='names', label = F, label.size = 3,pt.size.factor=1.75)

SpatialDimPlot(subset(data_sliding,niches %in% c('4') & names %in% c('Macrophage_Kupffer')),group.by='names', label = F, label.size = 3,pt.size.factor=1.75)




####################
####### S2B ########
####################
pdf(paste0('~/Documents/TSC_Paper/', 'Spatial_Feats.pdf'), width=5, height=5)
p1 <- SpatialFeaturePlot(data_sliding, features = "nFeature_Spatial",shape=22,pt.size.factor=1.75) 
ggrastr::rasterise(p1,layers='Spatial',dpi = 300)
dev.off()

####################
####### S2C ########
####################

pdf(paste0('~/Documents/TSC_Paper/', 'Spatial_UMAP.pdf'), width=6, height=4)
DimPlot(combo_all,group.by='predicted.celltype',reduction='ref.umap',label=T,raster=TRUE,pt.size = 2)
dev.off()

####################
####### S2D ########
####################

tbl.all <- t(table(combo_all$predicted.celltype,combo_all$orig.ident)) / rowSums(t(table(combo_all$predicted.celltype,combo_all$orig.ident)))
tbl.all <- tbl.all['HCC',]
df <- data.frame(grp = as.factor(names(tbl.all)),val = tbl.all,x=1)
df <- df[order(df$val),]
df$grp <- factor(df$grp,levels=df$grp)

pdf(paste0('~/Documents/TSC_Paper/', 'Cell_abundance_pie.pdf'), width=3, height=4)
ggplot(df,aes(y=val,fill=grp,x=x)) + geom_bar(position='stack',stat='identity') + coord_polar("y", start=0) + theme_classic()
dev.off()

idx <- rownames(df) %in% c("HSC-A","HSC-N","ENDO","Macrophage_Inflamed","Macrophage_Kupffer")
df <- df[idx,]
df$val <- df$val / sum(df$val)

pdf(paste0('~/Documents/TSC_Paper/', 'Cell_abundance_pie_of_int.pdf'), width=3, height=4)
ggplot(df,aes(y=val,fill=grp,x=x)) + geom_bar(position='stack',stat='identity') + coord_polar("y", start=0) + theme_classic()
dev.off()


####################
####### S2E ########
####################

#HSC_N <- c('Ecm1','Dcn','Sod3','Prelp')
#HSC_A <- c('Col3a1','Col1a1','Col1a2','Acta2')
#M_Kupffer <- c('Clec4f','Cd5l','Marco','C1qc')
#M_inflammed <- c('Cd74','H2-Aa','H2-Ab1','H2-Eb1')
#H_injured <- c('Saa1','Saa2','Hp','Lrg1')
#HPC <- c('Spp1','Mmp7','Clu','Epcam')
#Hep_MT <- c('mt-Co1','mt-Co2','mt-Atp6','mt-Cytb')
#Hep_PP <- c('Gpx3','Cyp2f2','Apoe','Alb')
#Hep_PC <- c('Cyp2e1','Cyp2c69','Cyp2a5','Cyp1a2')
#Hep_Nuc <- c('Neat1','Malat1','Mlxipl','Errfi1')



clust_markers <- read.csv('~/Documents/TSC_Paper/SeqScope_Cluster_marks.csv')
ir <- clust_markers %>% group_by(cluster) %>% select(gene)
ir <- group_split(ir)
gene_scores <- lapply(ir,function(x){ x$gene})
ir <- unlist(lapply(ir,function(x){unique(x$cluster)}))

combo_all <- AddModuleScore(combo_all,gene_scores,name='marks')
combo_all$predicted.celltype <- factor(combo_all$predicted.celltype,levels=ir)

pdf(paste0('~/Documents/TSC_Paper/', 'Spatial_Dot.pdf'), width=8, height=4)
DotPlot(subset(combo_all,orig.ident=='HCC'),c('marks1','marks2','marks3','marks4','marks5','marks6','marks7','marks8','marks9','marks10','marks11','marks12','marks13'),group.by = 'predicted.celltype',col.min = 0, col.max=2)
dev.off()

####################
####### S2F ########
####################

TFE3_genes <- read.csv('~/Documents/TSC_Paper/TFE3.vs.WT.csv')
TFE3_genes <- subset(TFE3_genes,log2FoldChange < 0 & padj < 0.05)$gene
TFE3_dgenes <- read.csv('~/Documents/TSC_Paper/dKO.vs.TSC.csv')
TFE3_dgenes <- subset(TFE3_dgenes,log2FoldChange < 0 & padj < 0.05)$gene
clear_genes <- read.csv('~/Documents/TSC_Paper/CLEAR_genes.csv',header=FALSE)$V1

combo_all <- AddModuleScore(combo_all,list(TFE3_dgenes,clear_genes))
combo_all$orig.ident <- factor(combo_all$orig.ident,levels=c('normal','td','HCC'))

p1 <- VlnPlot(combo_all,'Cluster1',pt.size=0,group.by='orig.ident',assay='SCT')

pdf(paste0('~/Documents/TSC_Paper/', 'TFE3_targets_all.pdf'), width=3, height=4)
p1
dev.off()


hep <- c("Hep_PP","Hep_MT","Hep_Injured","Hep_Eef1a1","Hep_Nuc","Hep_PC")
prog <- c("HPC")
stellate <- c("HSC-A","HSC-N")
mac <- c("Macrophage_Inflamed","Macrophage_Kupffer")
endo <- c("ENDO")
rbc <- c("RBC")
combo_all$sup_names <- as.character(combo_all$predicted.celltype)
combo_all$sup_names[combo_all$sup_names %in% hep] = 'hep'
combo_all$sup_names[combo_all$sup_names %in% prog] = 'prog'
combo_all$sup_names[combo_all$sup_names %in% stellate] = 'stellate'
combo_all$sup_names[combo_all$sup_names %in% mac] = 'mac'
combo_all$sup_names[combo_all$sup_names %in% endo] = 'endo'
combo_all$sup_names[combo_all$sup_names %in% rbc] = 'rbc'
combo_all$sup_names <- factor(combo_all$sup_names,levels=c('hep','prog','stellate','mac','endo','rbc'))

combo_all$group.sup_names <- paste0(combo_all$orig.ident,'_',combo_all$sup_names)

lvs <- c(paste0('normal','_',c('hep','prog','stellate','mac','endo','rbc')),
paste0('td','_',c('hep','prog','stellate','mac','endo','rbc')),
paste0('HCC','_',c('hep','prog','stellate','mac','endo','rbc')))


combo_all$group.sup_names <- factor(combo_all$group.sup_names,levels=lvs)

pdf(paste0('~/Documents/TSC_Paper/', 'TFE3_targets_celltype.pdf'), width=4, height=5)
DotPlot(combo_all,'Cluster1',group.by='group.sup_names',col.min=0,col.max=1,dot.min=0.5)
dev.off()


####################
####### S2G ########
####################

off_int <- subset(combo_all,predicted.celltype %in% c('HSC-A','HSC-N'))
off_int$predicted.celltype <- factor(off_int$predicted.celltype,levels=c('HSC-N','HSC-A'))
tbl.i <- t(table(off_int$predicted.celltype,off_int$orig.ident)) / rowSums(t(table(off_int$predicted.celltype,off_int$orig.ident)))
df <- data.frame(class=factor(c('Normal','TD','HCC'),levels=c('Normal','TD','HCC')),active=tbl.i[,'HSC-A'])
pdf(paste0('~/Documents/TSC_Paper/', 'Stellate.pdf'), width=3, height=4)
ggplot(df,aes(x=class,y=active)) + geom_bar(stat="identity") + theme_classic()
dev.off()


off_int <- subset(combo_all,predicted.celltype %in% c('Macrophage_Inflamed','Macrophage_Kupffer'))
off_int$predicted.celltype <- factor(off_int$predicted.celltype,levels=c('Macrophage_Inflamed','Macrophage_Kupffer'))
tbl.i <- t(table(off_int$predicted.celltype,off_int$orig.ident)) / rowSums(t(table(off_int$predicted.celltype,off_int$orig.ident)))
df <- data.frame(class=factor(c('Normal','TD','HCC'),levels=c('Normal','TD','HCC')),active=tbl.i[,'Macrophage_Inflamed'])
pdf(paste0('~/Documents/TSC_Paper/', 'Macro.pdf'), width=3, height=4)
ggplot(df,aes(x=class,y=active)) + geom_bar(stat="identity") + theme_classic()
dev.off()

off_int <- subset(combo_all,predicted.celltype %in% c('HPC',"Hep_PP","Hep_MT","Hep_Injured","Hep_Eef1a1","Hep_Nuc","Hep_PC"))
off_int$predicted.celltype <- factor(off_int$predicted.celltype,levels=c('HPC',"Hep_PP","Hep_MT","Hep_Injured","Hep_Eef1a1","Hep_Nuc","Hep_PC"))
tbl.i <- t(table(off_int$predicted.celltype,off_int$orig.ident)) / rowSums(t(table(off_int$predicted.celltype,off_int$orig.ident)))
df <- data.frame(class=factor(c('Normal','TD','HCC'),levels=c('Normal','TD','HCC')),active=tbl.i[,'HPC'])
pdf(paste0('~/Documents/TSC_Paper/', 'HPC.pdf'), width=3, height=4)
ggplot(df,aes(x=class,y=active)) + geom_bar(stat="identity") + theme_classic()
dev.off()

####################
####### S2H ########
####################

hep <- c("Hep_PP","Hep_MT","Hep_Injured","Hep_Eef1a1","Hep_Nuc","Hep_PC")
prog <- c("HPC")
stellate <- c("HSC-A","HSC-N")
mac <- c("Macrophage_Inflamed","Macrophage_Kupffer")
endo <- c("ENDO")
rbc <- c("RBC")
data_sliding$sup_names <- as.character(data_sliding$names)
data_sliding$sup_names[data_sliding$sup_names %in% hep] = 'hep'
data_sliding$sup_names[data_sliding$sup_names %in% prog] = 'prog'
data_sliding$sup_names[data_sliding$sup_names %in% stellate] = 'stellate'
data_sliding$sup_names[data_sliding$sup_names %in% mac] = 'mac'
data_sliding$sup_names[data_sliding$sup_names %in% endo] = 'endo'
data_sliding$sup_names[data_sliding$sup_names %in% rbc] = 'rbc'
data_sliding$sup_names <- factor(data_sliding$sup_names,levels=c('hep','prog','stellate','mac','endo','rbc'))

pdf(paste0('~/Documents/TSC_Paper/', 'Spatial_Cell_Dist.pdf'), width=5, height=5)
p1 <- SpatialDimPlot(subset(data_sliding,sup_names != 'hep'),group.by='sup_names', label = F, label.size = 3,shape=22,pt.size.factor=3) +  theme(aspect.ratio = 1)
ggrastr::rasterise(p1,layers='Spatial',dpi = 300)
dev.off()

####################
####### S2I ########
####################

pdf(paste0('~/Documents/TSC_Paper/', 'Spatial_Niche_Dist.pdf'), width=5, height=5)
p1 <- SpatialDimPlot(data_sliding,group.by='niches', label = F, label.size = 3,shape=22,pt.size.factor=1.35,stroke=0) +  theme(aspect.ratio = 1)
ggrastr::rasterise(p1,layers='Spatial',dpi = 300)
dev.off()


df <- table(data_sliding$niches,data_sliding$names) / rowSums(table(data_sliding$niches,data_sliding$names))*100

df <- df[,c('ENDO','Hep_Injured','HPC','HSC-A','HSC-N','Macrophage_Inflamed','Macrophage_Kupffer')]
df <- as.data.frame(t(t(df) / rowSums(t(df))))

pdf(paste0('~/Documents/TSC_Paper/', 'Spatial_Niche_Num.pdf'), width=7, height=4)
ggplot(df,aes(Var1,Freq,fill=Var1)) + geom_bar(stat='identity') + facet_wrap(~Var2,ncol=7) + theme_classic()
dev.off()

####################
####### S2J ########
####################

idx <- intersect(intersect(rownames(subset(combo_all,orig.ident=='HCC')),rownames(subset(combo_all,orig.ident=='normal'))),rownames(subset(combo_all,orig.ident=='td')))

macro.resp <- FindMarkers(combo_all, ident.1 = c("Hep_PP_HCC","Hep_MT_HCC","Hep_Injured_HCC","Hep_Eef1a1_HCC","Hep_Nuc_HCC","Hep_PC_HCC"), ident.2 = c("Hep_PP_normal","Hep_MT_normal","Hep_Injured_normal","Hep_Eef1a1_normal","Hep_Nuc_normal","Hep_PC_normal"), verbose = FALSE,recorrect_umi=T,only.pos=F,features=idx,logfc.threshold = 0)
macro.resp <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01)
cat(rownames(macro.resp),sep='\n')
cat(rownames(subset(macro.resp,avg_log2FC>0 & p_val_adj < 0.05)),sep='\n')

# Cleanly MITF
macro.resp <- FindMarkers(combo_all, ident.1 = c("Hep_PP_HCC"), ident.2 = c("Hep_PP_normal"), verbose = FALSE,recorrect_umi=T,only.pos=F,features=idx,logfc.threshold = 0)
macro.resp <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01)
cat(rownames(macro.resp),sep='\n')
cat(rownames(subset(macro.resp,avg_log2FC>0 & p_val_adj < 0.05)),sep='\n')


library(enrichR)
library(forcats)
library(DOSE)
library(patchwork)

dbs <- c("ChEA_2022","Reactome_Pathways_2024","GO_Biological_Process_2023")
enriched <- enrichr(rownames(subset(macro.resp,avg_log2FC>0 & p_val_adj < 0.05)), dbs, background = rownames(macro.resp), include_overlap = TRUE)

enriched[[1]] <- subset(enriched[[1]],Adjusted.P.value < 0.05)
pdf(paste0('~/Documents/TSC_Paper/', 'Hep_PP_Enrich_CHEA.pdf'), width=7, height=4)
p1<- ggplot(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:3),], 
  (aes(x=Adjusted.P.value, y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + 
geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + 
ggtitle('ChEA') + 
scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:3),]$Term,"_"), `[`, 1))) + 
theme(axis.text=element_text(colour="black"))
p1
dev.off()


enriched <- enrichr(rownames(macro.resp), dbs, include_overlap = TRUE)

enriched[[1]] <- subset(enriched[[1]],Adjusted.P.value < 0.05)
pdf(paste0('~/Documents/TSC_Paper/', 'Hep_PP_back_Enrich_CHEA.pdf'), width=7, height=4)
p2<- ggplot(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:3),], 
  (aes(x=Adjusted.P.value, y=fct_inorder(Term), color = as.numeric(Combined.Score), size=parse_ratio(Overlap)))) + 
geom_point() + xlab('Combined Score') + ylab('Term') + labs(color="Combined Score",size="Overlap") + theme_classic()  + 
ggtitle('ChEA') + 
scale_y_discrete(labels= fct_inorder(sapply(strsplit(enriched[[1]][order(enriched[[1]]$Adjusted.P.value,decreasing=F),][rev(1:3),]$Term,"_"), `[`, 1))) + 
theme(axis.text=element_text(colour="black"))
p2
dev.off()

pdf(paste0('~/Documents/TSC_Paper/', 'Hep_PP_Enrich_CHEA_overlay.pdf'), width=7, height=4)
p1 / p2
dev.off()

####################
####### S2K ########
####################

pdf(paste0('~/Documents/TSC_Paper/', 'TFE3_targets_in_HCC.pdf'), width=7, height=4)
DotPlot(subset(combo_all,orig.ident=='HCC'),'Cluster1',group.by='predicted.celltype',col.min=0,col.max=2,dot.min=0.5)
dev.off()

####################
####### S2L ########
####################


macro.resp <- FindMarkers(combo_all, ident.1 = c("Hep_PP_HCC"), ident.2 = c("Hep_PP_normal"), verbose = FALSE,recorrect_umi=T,only.pos=F,features=idx,logfc.threshold = 0)
HCC.vs.normal <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01)


macro.resp <- FindMarkers(combo_all, ident.1 = c("Hep_PP_HCC"), ident.2 = c("Hep_PP_td"), verbose = FALSE,recorrect_umi=T,only.pos=F,features=idx,logfc.threshold = 0)
HCC.vs.td <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01)

overlap <- intersect(rownames(HCC.vs.normal),rownames(HCC.vs.td))

df <- data.frame(gene = overlap, x = HCC.vs.normal[overlap,]$avg_log2FC, y = HCC.vs.td[overlap,]$avg_log2FC)

df$color <- 'black'
df$color[df$gene %in% clear_genes] <- 'red'

df$label <- df$gene
df$label[!(df$gene %in% clear_genes)] <- NA

pdf(paste0('~/Documents/TSC_Paper/', 'HCC_vs_norm_and_td.pdf'), width=7, height=5.5)
ggplot(df,aes(x,y,label=label,color = color)) + geom_point() + geom_text_repel(max.overlaps = Inf) + scale_color_manual(values = c('black','red')) + theme_classic()
dev.off()



a<-FindMarkers(subset(combo_all,orig.ident=='HCC'),ident.1='Macrophage_Inflamed_HCC',only.pos=T)
mac_marks <- rownames(subset(a,p_val_adj<0.05))


data_sliding <- SetIdent(data_sliding,value = 'niches')
DefaultAssay(data_sliding) <- 'Spatial'
data_sliding_new <- SCTransform(data_sliding,assay='Spatial',conserve.memory=T)
data_sliding_new<-PrepSCTFindMarkers(data_sliding_new)

a<-FindMarkers(subset(data_sliding_new,predicted.celltype=='Macrophage_Inflamed'),ident.1='2',only.pos=F,features=intersect(idx_non_hep,rownames(data_sliding_new)))
cat(rownames(subset(a,p_val_adj<0.05 & avg_log2FC>0)),sep='\n')
cat(rownames(subset(a,p_val_adj<0.05 & avg_log2FC<0)),sep='\n')

cat(rownames(subset(a)),sep='\n')

SpatialFeaturePlot(data_sliding_new,'Cd74', shape=22,pt.size.factor=3) +  theme(aspect.ratio = 1)








library(EnhancedVolcano)

keyvals <- ifelse(rownames(macro.resp) %in% TFE3_dgenes, 'red','gray')
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'gray'] <- 'OTHER'
names(keyvals)[keyvals == 'red'] <- 'CLEAR'

EnhancedVolcano(macro.resp, x = 'avg_log2FC', y = 'p_val_adj', 
lab = rownames(macro.resp),pCutoff = 0.05, 
selectLab = rownames(macro.resp)[which(names(keyvals) %in% c('red'))],
colCustom = keyvals,
FCcutoff = 0.1,xlim = c(-4, 4))





macro.resp <- FindMarkers(combo_all, ident.1 = "Macrophage_Kupffer_HCC", ident.2 = "Macrophage_Kupffer_normal", verbose = FALSE,recorrect_umi=T,only.pos=F,features=idx_non_hep)
macro.resp <- subset(macro.resp,pct.1 > 0.01 & pct.2 > 0.01)

EnhancedVolcano(macro.resp, x = 'avg_log2FC', y = 'p_val_adj', lab = rownames(macro.resp),pCutoff = 0.05, xlim = c(-4, 4))

______________



combo_all <- RunUMAP(combo_all, dims = 1:30, reduction = "integrated.cca")
DimPlot(combo_all,reduction='umap',group.by='names')

combo_all$celltype.ident <- paste(combo_all$names, combo_all$orig.ident, sep = "_")
Idents(combo_all) <- "celltype.ident"

slot(object = combo_all@assays$SCT@SCTModel.list[[2]], name="umi.assay")<-"Spatial"
combo_all <- PrepSCTFindMarkers(combo_all)

saveRDS(combo_all,'~/Downloads/SeqScope/clean_normal_td_HCC_integrated.rds')


_______________________




sce.raw <- as.SingleCellExperiment(data_new)
data_new <- subset(data_new,nFeature_Spatial > 50)
z <- 1 * (data_new$nFeature_Spatial > 100)
sce <- as.SingleCellExperiment(data_new)
sce <- decontX(sce,background = sce.raw, z=z)


umap <- reducedDim(sce, "decontX_UMAP")
plotDimReduceCluster(x = sce$decontX_clusters,
    dim1 = umap[, 1], dim2 = umap[, 2])


data_new2 <- CreateSeuratObject(counts=decontXcounts(sce),assay="Spatial")
data_new2@meta.data <- data_new@meta.data
data_new2$image <- data_new$image
data_new2$nCount_Spatial <- colSums(decontXcounts(sce))
data_new2$nFeature_Spatial <- colSums(decontXcounts(sce)>0)
data_new2 <- subset(data_new2,nFeature_Spatial > 100)


VlnPlot(data_new2, features = "nFeature_Spatial", pt.size = 0.1)
SpatialFeaturePlot(data_new2, features = "nCount_Spatial",shape=22,pt.size.factor=3.5)


data_new2 <- SCTransform(data_new2, assay = "Spatial", verbose = FALSE, variable.features.n = 3000)
#data_new2 <- AddModuleScore(data_new2,list(unique(unlist(zone[,1:3])),unique(unlist(zone[,4:6]))),ctrl = 10)
#SpatialFeaturePlot(data_new2, features = c("Cluster1"),shape=22,pt.size.factor=3.5,min.cutoff=0)
#VariableFeatures(data_new2) <- unique(unlist(zone))[unique(unlist(zone)) %in% rownames(data_new2)]

data_new2 <- RunPCA(data_new2, assay = "SCT", verbose = FALSE)
data_new2 <- FindNeighbors(data_new2, reduction = "pca", dims = 1:30)
data_new2 <- FindClusters(data_new2, verbose = FALSE)
data_new2 <- RunUMAP(data_new2, reduction = "pca", dims = 1:30)

DimPlot(data_new2, reduction = "umap", label = TRUE)
FeaturePlot(data_new2, 'Cluster1', reduction = "umap", label = TRUE)


SpatialDimPlot(data_new2, label = TRUE, label.size = 3,pt.size.factor=3.5)
VlnPlot(data_new,'nFeature_Spatial')

#0 is Hgfac, Sardh, Por
#1 is Hbb-bs, Coq8a, Sephs2
#2 is Cdh23, Acbd6, Nccrp1
#3 is Cdrt4os2, Znhit1, Ppm1h <- ECs
#4 is Nxpe4, Pard3, Egfem1
cat(rownames(subset(FindMarkers(data_new2,only.pos=T,ident.1=4),p_val_adj<0.05)),sep='\n')


head(FindMarkers(data_new,ident.1=0))

data0 <- data_new[,data_new@active.ident == 0]

SpatialFeaturePlot(data_new2, features = "Nxpe4",shape=22,pt.size.factor=3.5)
data0 <- SCTransform(data0, assay = "Spatial", verbose = FALSE)
data0 <- RunPCA(data0, assay = "SCT", verbose = FALSE)
data0 <- FindNeighbors(data0, reduction = "pca", dims = 1:30)
data0 <- FindClusters(data0, verbose = FALSE)
data0 <- RunUMAP(data0, reduction = "pca", dims = 1:30)
DimPlot(data0, reduction = "umap", label = TRUE
SpatialDimPlot(data0, label = TRUE, label.size = 3,pt.size.factor=3.5)










