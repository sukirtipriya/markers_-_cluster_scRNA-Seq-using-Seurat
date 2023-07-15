
# script to identify cluster identity........................
# Finding markers in every cluster
# finding conserved markes
# finding marked DE between conditions

library(Seurat)
library(tidyverse)

# load data

ifnb_harmony <- readRDS('../ifnb_harmony.rds')
str(ifnb_harmony)
view(ifnb_harmony@meta.data)

# visualize data

DimPlot(ifnb_harmony, reduction = 'umap', group.by = 'seurat_clusters', label =T)
condition <- DimPlot(ifnb_harmony, reduction = 'umap', group.by ='stim')

conditionalclusters

# findAll markers-----------------------------------------------

FindAllMarkers(ifnb_harmony,
                   logfc.threshold = 0.25,
                   min.pct = 0.1,
                   only.pos = TRUE,
                   test.use = 'DESeq2',
                   slot = 'counts')

# findConserved markers------------------------------------------
 
# Notes:
# slot depends on the type of the test used,
# default is data slot that stores normalized data
# DefaultAssay(ifnb_harmony)

DefaultAssay(ifnb_harmony)

markers_cluster3 <- FindConservedMarkers(ifnb_harmony,
                                     ident.1 = 3,
                                     grouping.var = 'stim')

head(markers_cluster3)

# let's visualize top features                                                        

FeaturePlot(ifnb_harmony, features = c('FCGR3A'), min.cutoff = 'q10')

# min-cut off explanation:
seq(1,5)
SetQuantile('q50', seq(1,5))
SetQuantile('q10', seq(1,5))

# rename cluster 3 ident

Idents(ifnb_harmony)
ifnb_harmony <- RenameTdents(ifnb_harmony, '3' = 'CD16 Mono')

DimPlot(ifnb_harmony, reduction = 'umap', label = T)

# cells already have annotations provided in the metadata

view(ifnb_harmony@meta.data)

# setting cluster identities is an iterative step
# multiple approaches could be taken - automatic/manual annotation (sometimes both)
# need to make sure each cell type cell type forms a separate cluster

# setting Idents as seurat annotations provided (also a sanity check!)

Idents(ifnb_harmony) <- ifnb_harmony@meta.data$seurat_annotations
Idents(ifnb_harmony)

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)

# findMarkers between conditions------------------------------------------

ifnb_harmony$celltype.cnd <- paste0(ifnb_harmony$seurat_annotations, '_', ifnb_harmony$stim)
view(ifnb_harmony@meta.data)
Idents(ifnb_harmony) <- ifnb_harmony$celltype.cnd

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)

# find markers

b.interferon.response <- FindMarkers(ifnb_harmony, ident.1 = 'CD16 Mono_STIM',
                                     ident.2 = 'CD16 Mono_CTRL')

# plotting conserved features vs Features between conditions 

head(markers_cluster3)

FeaturePlot(ifnb_harmony, features = c('FCGR3A','AIF1','IFIT1'), split.by = 'stim',
                       min.cutoff = 'q10')






