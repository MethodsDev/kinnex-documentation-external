3.2 Vignette for Tertiary processing for SC-Kinnex
===================================================

(The subsamples dataset used in this vignette has since been archived, and is currently in the process of being replaced with larger complete PBMC dataset.)


This vigentte leverages various parts of the `Seurat package <https://satijalab.org/seurat/>`_ and follows along in parts the `"Seurat - Guided Clustering Tutorial" <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>`_
The standalone utility scIsoseqUtil.py, developed at MDL, creates sparce matrices from Isoquant Outs namely, transcript_model_reads and transcript_models_gtfs. 
The script is provided in source repo here. 


Creating sparse matrices for use with Seurat
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setting the environment:

.. code:: bash

    sudo apt install r-base
    sudo R -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
    sudo R -e 'BiocManager::install("argparse")'
    conda create -n scIsoseqUtil
    conda activate scIsoseqUtil
    conda install bioconda::r-argparse
    pip install pysam



.. code:: bash

    # from:
    # https://github.com/MethodsDev/scIsoquantMatrixBuilder

    wget https://github.com/MethodsDev/kinnex-documentation-external/archive/refs/heads/main.zip

    scIsoseqUtil.py --sample_id ${sample_id} \
                    --bam ${sample_id}.aligned.sorted.bam \
                    --transcript_model_reads ${sample_id}.transcript_model_reads.tsv.gz \
                    --transcript_models_gtf ${sample_id}.transcript_models.gtf.gz

Analysing sparse matrices created above
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The code below is an ``R`` code, blocks can be copied to ``Rmd`` to excute locally:

.. code:: bash

    install_if_missing <- function(packages) {
    if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
        install.packages(setdiff(packages, rownames(installed.packages())))
        }
    }

    install_if_missing(c('tidyverse','stringr','dplyr','edgeR','ggrepel','DESeq2','Seurat','clustermole'))



.. code:: bash

    # {r setup, include=FALSE}
        knitr::opts_chunk$set(echo = TRUE)
        library(tidyverse)
        library(Seurat)


Input counts matrix created above from step1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    #{r}
        data_dir = "scKinnex.genes-sc_matrix_from_isoquant/"
        output_prefix = "scKinnex.genes"


Reading data in using Read10x()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    #{r}
        data = Read10X(data.dir=data_dir,
               gene.column = 1,
               cell.column = 2,
               unique.features = TRUE,
               strip.suffix = FALSE)


UMI counts per cell:
~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    #{r}
        umi_counts_per_cell = colSums(data)


sorting:
~~~~~~~~

.. code:: bash

    #{r}
        umi_counts_per_cell = sort(umi_counts_per_cell, decreasing = T)

plotting :
~~~~~~~~~~~

.. code:: bash

    #{r}
        plot(umi_counts_per_cell, log='xy')
        ggsave(filename='PBMC_complete_umi_counts_per_cell.png',path=data_dir, plot = last_plot())


.. image:: ../_images/umi_counts_per_cell.png
   :align: center


Creating seurat object from counts matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    #{r}
    seurat_obj <- CreateSeuratObject(counts = data, project = "project", min.cells = 3, min.features = 200)
    seurat_obj


Terminal Out:

11390 features across 500 samples within 1 assay 
Active assay: RNA (11390 features, 0 variable features)
1 layer present: counts

.. code:: bash

    #{r}
    # before filtering
    seurat_obj@meta.data %>% summarize(median(nCount_RNA), median(nFeature_RNA))    


Terminal Out:

median(nCount_RNA)        median(nFeature_RNA)
<dbl>                     <dbl>
2794.17                   799

PercentageFeatureSet - Calculate the percentage of all counts that belong to a given set of features

.. code:: bash

    #{r}
        seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")


Exploring seurat object:


.. code:: bash

    #{r}
        seurat_obj
        seurat_obj@meta.data %>% head()



UMI counts per cell
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    #{r}
    seurat_obj@meta.data %>% dplyr::select(nCount_RNA) %>% 
    arrange(desc(nCount_RNA)) %>% 
    dplyr::mutate(i=row_number()) %>%
    ggplot(aes(x=i, y=nCount_RNA)) + geom_point() + theme_bw() + 
    scale_y_continuous(trans='log10') +
    scale_x_continuous(trans='log10') +
    ggtitle("nCount_RNA: UMI counts per cell")


.. image:: ../_images/nCount_RNA-umi_counts_per_cell.png
   :align: center


Feature counts per cell:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    #{r}
        seurat_obj@meta.data %>% dplyr::select(nFeature_RNA) %>% arrange(desc(nFeature_RNA)) %>% dplyr::mutate(i=row_number()) %>%
        ggplot(aes(x=i, y=nFeature_RNA)) + geom_point() + theme_bw() + 
        scale_y_continuous(trans='log10') +
        scale_x_continuous(trans='log10') +
        ggtitle("nFeature_RNA: gene count per cell") 


.. image:: ../_images/nFeature_RNA-gene_count_per_cell.png
   :align: center


Visualize QC metrics as a violin plot

VlnPlot: Draws a violin plot of single cell data (gene expression, metrics, PC scores, etc.)

.. code:: bash

    #{r}
    # Visualize QC metrics as a violin plot
    VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

.. image:: ../_images/VlnPlot.png
   :align: center


.. code:: bash

    #{r}
    plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    plot1 + plot2

.. image:: ../_images/FeatureScatter.png
   :align: center



Saving object.RDS

.. code:: bash

    #{r}
    #saveRDS(seurat_obj, file = paste0(data_dir, output_prefix, "-seurat_obj-preCellFiltering.rds"))


filtering cells on percent Mitochondria:

.. code:: bash

    #{r}
    #### filtering cells on 
    length(seurat_obj$percent.mt < 15)
    seurat_obj <- subset(seurat_obj, 
                     percent.mt < 15)
    seurat_obj


An object of class Seurat 
11390 features across 447 samples within 1 assay 
Active assay: RNA (11390 features, 0 variable features)
1 layer present: counts

Summarize:

.. code:: bash

    #{r}
        seurat_obj@meta.data %>% summarize(median(nCount_RNA), median(nFeature_RNA))


Terminal Output:
median(nCount_RNA)      median(nFeature_RNA)
<dbl>                   <int>
2870.08	                817


NormalizeData : Normalize the count data present in a given assay.
Normalization methods =
“LogNormalize”: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. 
This is then natural-log transformed using log1p.

.. code:: bash

    #{r}
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)



FindVariableFeatures: Identifies features that are outliers on a 'mean variability plot'.

selection.method =
“vst”: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). 
Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). 
Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).

.. code:: bash

    #{r}
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(seurat_obj), 10)

    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(seurat_obj)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    plot1 + plot2

.. image:: ../_images/VariableFeatures.png
   :align: center

ScaleData: 
Scales and centers features in the dataset. 
If variables are provided in vars.to.regress, they are individually regressed against each feature, and the resulting residuals are then scaled and centered.

.. code:: bash

    #{r}
    all.features <- rownames(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj, features = all.features)


Performing  PCA :
~~~~~~~~~~~~~~~~~

RunPCA: Run Principal Component Analysis on gene expression using IRLBA. For details about stored PCA calculation parameters, see `PrintPCAParams`.
VizDimLoadings: Visualize top genes associated with reduction components
DimPlot:
Graphs the output of a dimensional reduction technique (PCA by default). Cells are colored by their identity class.

.. code:: bash

    #{r}
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
    VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
    DimPlot(seurat_obj, reduction = "pca") + NoLegend()
    DimHeatmap(seurat_obj, dims = 1:3, cells = 500, balanced = TRUE)
    ElbowPlot(seurat_obj)


.. figure:: ../_images/DimPlot.png
    :height: 500px
    :width: 1000px
    :align: center


.. figure:: ../_images/ElbowPlot.png
   :height: 500px
   :width: 1000px
   :align: center


Generating UMAP : 
~~~~~~~~~~~~~~~~~

.. code:: bash

    #{r}
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
    DimPlot(seurat_obj, reduction = "umap")

    FeaturePlot(seurat_obj, features = c("nFeature_RNA"))

    FeaturePlot(seurat_obj, features = c("nCount_RNA"))

    FeaturePlot(seurat_obj, features = c("percent.mt"))

Feature Count plots from terminal out:

.. list-table:: 
    :widths: 50 50 

    * - .. figure:: ../_images/UMAP_DimPlot.png
           :alt: UMAP_DimPlot.png

           UMAP_DimPlot

      - .. figure:: ../_images/nFeature_RNA_FeaturePlot.png
           :alt: nFeature_RNA_FeaturePlot.png

           nFeature_RNA_FeaturePlot

Feature Count plots from terminal out:

.. list-table:: 
    :widths: 50 50

    * - .. figure:: ../_images/nCount_RNA_FeaturePlot.png
           :alt: nFeature_RNA_FeaturePlot

           nFeature_RNA_FeaturePlot

      - .. figure:: ../_images/percent_mt_FeaturePlot.png
           :alt: percent_mt_FeaturePlot

           percent_mt_FeaturePlot


.. code:: bash

    #{r}
    # counts and fractions of cells

    cluster_counts_n_fracs = seurat_obj@meta.data %>% group_by(seurat_clusters) %>% tally() %>%  mutate(frac=prop.table(n))

    cluster_counts_n_fracs

    saveRDS(seurat_obj, file = paste0(output_prefix, "-seurat_obj.rds"))


Terminal Out:

seurat_clusters n frac
<fctr> <int> <dbl>
0	219	0.52771084	
1	128	0.30843373		
2	45	0.10843373	
3	23	0.05542169	


DE, find markers:
~~~~~~~~~~~~~~~~~

find markers for every cluster compared to all remaining cells, report only the positive ones

.. code:: bash

    #{r}
    # find markers for every cluster compared to all remaining cells, report only the positive
    # ones
    seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
    seurat_obj.markers %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 1)   

.. code:: bash

    #{r}
    top_20_markers = seurat_obj.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>% slice_head(n=20) %>% ungroup()


    top_20_markers

.. code:: bash

    #{r}
    max_cluster <- max(as.numeric(top_20_markers$cluster)) - 1

    for (clnum in 0:max_cluster) {
        cluster = top_20_markers %>% filter(cluster == clnum)
  
  
            gene.symbols = sapply(cluster$gene, function(x) { str_split(x, "\\^")[[1]][1] })
  
            gene.symbols = grep("ENSG|ENST|novel", gene.symbols, value=T, invert=T)
  
        cat(paste0(clnum,":"))
        cat(gene.symbols, sep=",")
        cat("\n")
    }


.. code:: bash

    Terminal Out:

    0:IL7R,LTB,PRKCQ-AS1,RPL34,RCAN3,GAS5,TCF7,LEF1,MAL,CD27,CCR7,ANKRD44-AS1,RGCC,RGS10,NOSIP,TMEM123,CAMK4
    1:NKG7,GZMH,CST7,GZMA,GNLY,FGFBP2,CCL5,CCL4,PRF1,EFHD2,PLEK,HOPX,PFN1,GZMM,CALM1,GZMB,SH3BGRL3,CTSW,XCL2,TRGC2
    2:CD79A,IGHM,CD79B,BANK1,HLA-DQA1,BCL11A,HLA-DRA,TCL1A,TNFRSF13C,HLA-DMB,HLA-DRB1,SWAP70,VPREB3,RALGPS2
    3:CSTA,SERPINA1,CFD,VCAN,RGS2,MNDA,CD68,CYP27A1,RETN,CPVL,CLEC12A,LMO2,GRN,LST1,CYBB,NCF2,LILRA5,FCN1


Run above list through: http://xteam.xbio.top/ACT to get cell type predictions.
To read more on the ACT tool, the publication can be found here, 
`"Annotation of cell types (ACT): a convenient web server for cell type annotation". <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01249-5>`_
The detailed report can be navigated using `"Help" <http://xteam.xbio.top/ACT/help.jsp>`_ page for ACT. 


.. figure:: ../_images/cluster0_genelist.png
   :scale: 90%
   :align: center

.. figure:: ../_images/Treeofcellontology_cluster0.png
   :scale: 100%
   :align: center

.. figure:: ../_images/cluster1_genelist.png
   :scale: 90%
   :align: center

.. figure:: ../_images/Treeofcellontology_cluster1.png
   :scale: 100%
   :align: center


.. figure:: ../_images/cluster2_genelist.png
   :scale: 90%
   :align: center

.. figure:: ../_images/Treeofcellontology_cluster2.png
   :scale: 100%
   :align: center

.. figure:: ../_images/cluster3_genelist.png
   :scale: 90%
   :align: center

.. figure:: ../_images/Treeofcellontology_cluster3.png
   :scale: 100%
   :align: center

.. code:: bash

    #{r}
    # save files for later read/cell tracking

    write.table( Idents(seurat_obj), paste0(output_prefix, "-cell_cluster_assignments.tsv"), quote=F, row.names=T, sep="\t")

.. code:: bash

    #{r}
    saveRDS(seurat_obj, file = paste0(output_prefix, "-seurat_obj.rds"))


Examining specific gene sets example
Note, this helps to have the gene-symbol annotated gene features.

.. code:: bash

    #{r}
    # example definition of marker genes for certain cell types

    marker_genes = list()

    marker_genes[["CD8"]] = c("CD8A","CD8B","GZMB","TRB","TRA","PRF1", "GZMB")


.. code:: bash

    #{r}

    # function to extract gene ids with the relevant gene symbols


    feature_names = rownames(seurat_obj@assays$RNA$counts)

    get_feature_names_with_gene_symbols = function(gene_symbols) {
  
    gene_ids = c()
  
    for (gene_symbol in gene_symbols) {
        found_genes = grep(paste0(gene_symbol,"\\^"), feature_names, value=T) 
     if (length(found_genes) > 0) {
            gene_ids = c(gene_ids, found_genes)
        }
    }
    return(gene_ids)
    }

.. code:: bash

    #{r}
    # paint umaps according to the features of interest

    feature_ids = get_feature_names_with_gene_symbols(marker_genes[["CD8"]])

    VlnPlot(seurat_obj, features = feature_ids)
    FeaturePlot( seurat_obj, features = feature_ids)


.. figure:: ../_images/foi_VlnPlot.png
    :height: 500px
    :width: 1000px
    :align: center


.. figure:: ../_images/foi_FeaturePlot.png
   :height: 500px
   :width: 1000px
   :align: center