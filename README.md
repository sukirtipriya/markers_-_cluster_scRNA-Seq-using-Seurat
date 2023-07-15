# markers_-_cluster_scRNA-Seq-using-Seurat

Single-cell RNA-seq marker identification

Identifying gene markers for each cluster

Seurat has the functionality to perform a variety of analyses for marker identification; for instance, we can identify markers of each cluster relative to all other clusters by using the FindAllMarkers() function. This function essentially performs a differential expression test of the expression level in a single cluster versus the average expression in all other clusters.

To be identified as a cluster or cell type marker, within the FindAllMarkers() function, we can specify thresholds for the minimum percentage of cells expressing the gene in either of the two groups of cells (min.pct) and minimum difference in expression between the two groups (min.dff.pct).

![image](https://github.com/sukirtipriya/markers_-_cluster_scRNA-Seq-using-Seurat/assets/88479900/18c44d68-eb43-41ea-8f28-ca3015f3e87d)

Single-cell RNA-seq clustering analysis


![image](https://github.com/sukirtipriya/markers_-_cluster_scRNA-Seq-using-Seurat/assets/88479900/10120fc6-cb76-417f-a0e5-d706234c413e)

To identify clusters, the following steps will be performed:

    Normalization and transformation of the raw gene counts per cell to account for differences in sequencing depth per cell.
    Identification of high variance genes.
    Regression of sources of unwanted variation (e.g. number of UMIs per cell, mitochondrial transcript abundance, cell cycle phase).
    Identification of the primary sources of heterogeneity using principal component (PC) analysis and heatmaps.
    Clustering cells based on significant PCs (metagenes).


