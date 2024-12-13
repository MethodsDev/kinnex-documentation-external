3.2 Vignette for creating pseudo-bulk counts matrix 
====================================================

In order to perform isoform switching analysis downstream, we would create psuedo samples by breaking/sub-sampling from each cluster into 3 to create replicates.
By doing so, each cluster becomes 'condition' and each of the sub-clusters becomes a 'replicate'

For easy sub-sampling, we will annotate sample names as 'Sample_R+cluster.id+replicate{1-3}', for ex: 'Sample_R01' corresponds to a sub-cluster within cluster index 0.

We follow the pseudo-bulk vignette here closely:

Pseudo-bulk ref : https://hbctraining.github.io/scRNA-seq_online/lessons/pseudobulk_DESeq2_scrnaseq.html

Other Refs:

- Seurat Cheatsheet: https://hbctraining.github.io/scRNA-seq_online/lessons/seurat_cheatsheet.html
- Seurat commands: https://satijalab.org/seurat/articles/essential_commands.html#fetchdata



Installing BiocStyle:

.. code:: bash

    #{r}
    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
    BiocManager::install("BiocStyle")


.. code:: bash

    #{r}
    remotes::install_github("cvarrichio/Matrix.utils") 

loading required libraries:

.. code:: bash

    #{r}
    library(tidyverse)
    library(DESeq2)
    library(Seurat)
    library(SingleCellExperiment)
    library(BiocStyle)
    library(data.table)
    library(Matrix.utils)
    knitr::opts_chunk$set(warning=FALSE, error=FALSE, message=FALSE)
    library(dplyr)


Part1: Reading and Exploring Seurat object: 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`kinnex_sc/complete_dataset/PBMC_BioIVT_10x3p_complete_refGuided.genes-sc_matrix_from_isoquant/PBMC_complete_scKinnex.genes-seurat_obj.rds`

.. code:: bash

    #{r}
    seurat_readObj <- readRDS('kinnex_sc/complete_dataset/PBMC_BioIVT_10x3p_complete_refGuided.genes-sc_matrix_from_isoquant/PBMC_complete_scKinnex.genes-seurat_obj.rds')
    seurat_readObj
    length(unique(rownames(seurat_readObj@meta.data)))

Terminal Output:

An object of class Seurat:
30015 features across 12850 samples within 1 assay 
Active assay: RNA (30015 features, 2000 variable features)
3 layers present: counts, data, scale.data
2 dimensional reductions calculated: pca, umap


12850

.. code:: bash

    #{r}
    #library(dplyr)
    Features(seurat_readObj) %>% head()
    rownames(seurat_readObj) %>% head()
    all(rownames(seurat_readObj@meta.data) == Cells(seurat_readObj))
    Idents(seurat_readObj) %>% head()

Terminal Out: 

.. code:: bash
    
    [1] "novel-gene-chr10-13364^novel-gene-chr10-13364" "novel-gene-chr10-2287^novel-gene-chr10-2287"  
    [3] "novel-gene-chr10-26782^novel-gene-chr10-26782" "novel-gene-chr10-29791^novel-gene-chr10-29791"
    [5] "novel-gene-chr10-30176^novel-gene-chr10-30176" "novel-gene-chr10-31726^novel-gene-chr10-31726"
    [1] "novel-gene-chr10-13364^novel-gene-chr10-13364" "novel-gene-chr10-2287^novel-gene-chr10-2287"  
    [3] "novel-gene-chr10-26782^novel-gene-chr10-26782" "novel-gene-chr10-29791^novel-gene-chr10-29791"
    [5] "novel-gene-chr10-30176^novel-gene-chr10-30176" "novel-gene-chr10-31726^novel-gene-chr10-31726"
    [1] TRUE
    AAACAACGAAAGAATC AAACAACGACAGTCTA AAACAACGAGTTAGAA AAACACCTGCTTCCAC AAACACCTGGAGGAGG 
               0                3                0                1                1 
    AAACACCTGGTACCTA 
               2 
    Levels: 0 1 2 3 4 5 6 7 8 9 10 11



Part2: Dividing clusters into sub-clusters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

seperating cluster info

.. code:: bash

    #{r}
    seurat_readObj@meta.data
    seurat_clusters_info <- FetchData(object = seurat_readObj, vars = c("seurat_clusters"), layer = "meta.data")


randomly sampling indexes from each cluster:

.. code:: bash

    #{r}
    seurat_clusters_info$cell_bc <- rownames(seurat_clusters_info)


Dividing cluster into 3 subsets random sampling:

.. code:: bash

    #{r}
    for (id in unique(seurat_readObj$seurat_clusters)) {
        cluster_name <- paste("cluster_",id, sep="")
        print(cluster_name)
        df_temp <- seurat_clusters_info[seurat_clusters_info$seurat_clusters==id,]
        df_temp$sample_id <- sample(factor(rep(1:3, length.out=nrow(seurat_clusters_info[seurat_clusters_info$seurat_clusters==id,])), 
                          labels=paste0("sample_R",id,1:3)))
        assign(cluster_name,df_temp)
    }

adding samples ids to metadata objects:

.. code:: bash

    #{r}
    #summary(cluster_0$sample_id)
    #summary(cluster_2$sample_id)
    #summary(cluster_10$sample_id)
    #summary(cluster_11$sample_id)
    ss <- seurat_readObj
    w_sample_ids <- rbind(cluster_0, cluster_1,cluster_2,cluster_3, cluster_4, cluster_5, cluster_6,    
                      cluster_7,cluster_8,cluster_9,cluster_10,cluster_11)

    #ss@meta.data
    #w_sample_ids

Adding sample names in Sample:

.. code:: bash

    #{r}
    ss@meta.data <- merge(ss@meta.data, dplyr::select(w_sample_ids, sample_id), by=0, all=TRUE)
    ss@meta.data

Assigning back to Seurat Object:

.. code:: bash

    #{r}
    seurat_readObj <- ss

  
Extracting metadata:

.. code:: bash

    #{r}
    metadata <- seurat_readObj@meta.data


Extract raw counts and metadata to create SingleCellExperiment object

.. code:: bash

    #{r}
    counts <- seurat_readObj@assays$RNA$counts


Set up metadata as desired for aggregation and DE analysis

.. code:: bash

    #{r}
    metadata$cluster_id <- factor(seurat_readObj@active.ident)

Part3 - Create single cell experiment object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. code:: bash

    #{r}
    sce <- SingleCellExperiment(assays = list(counts = counts), 
                           colData = metadata)


Exploring the raw counts for the dataset
Checking the assays present

.. code:: bash

    #{r}
    assays(sce)

Terminal Out:

List of length 1
names(1): counts

Check the counts matrix

.. code:: bash

    #{r}
    dim(counts(sce))
    counts(sce)[1:6, 1:6]


Part4: Preparing the single-cell dataset for pseudobulk analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extracting necessary metrics for aggregation by cell type in a sample:

.. code:: bash

    #{r}
    # Extract unique names of clusters (= levels of cluster_id factor variable)
    cluster_names <- levels(colData(sce)$cluster_id)
    cluster_names

    # Total number of clusters
    length(cluster_names)


Number of cells in each cluster:

.. code:: bash

    #{r}
    for (i in cluster_names) {
     print(paste(i, length(colData(sce)$cluster_id[colData(sce)$cluster_id==i]), sep = ":"))
    }

Terminal Out:

"0:2314"
"1:2123"
"2:1940"
"3:1822"
"4:1604"
"5:1062"
"6:1035"
"7:281"
"8:280"
"9:235"
"10:101"
"11:53"


.. code:: bash

    #{r}
    # Extract unique names of samples (= levels of sample_id factor variable)
    sample_names <- levels(colData(sce)$sample_id)
    sample_names

    # Total number of samples
    length(sample_names)   


Part5: Subset metadata
~~~~~~~~~~~~~~~~~~~~~~~

Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)

.. code:: bash

    #{r}
    #colData(sce)
    groups <- colData(sce)[, c("cluster_id", "sample_id")]
    head(groups)

Aggregate across cluster-sample groups
- transposing row/columns to have cell_ids as row names matching those of groups

.. code:: bash

    #{r}
    aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups, fun = "sum") 

Exploring aggregated output matrix

.. code:: bash

    #{r}
    class(aggr_counts)
    dim(aggr_counts)
    aggr_counts[1:6, 1:6]

Transpose aggregated matrix to have genes as rows and samples as columns

.. code:: bash

    #{r}
    aggr_counts <- t(aggr_counts)
    aggr_counts[1:6, 1:6]

Understanding tstrsplit()

.. code:: bash

    #{r}
    ## Exploring structure of function output (list)
    tstrsplit(colnames(aggr_counts), "_") %>% str()

    ## Comparing the first 10 elements of our input and output strings
    head(colnames(aggr_counts), n = 10)
    head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)

    aggr_counts


.. code:: bash

    #{r}
    # As a reminder, we stored our cell types in a vector called cluster_names
    cluster_names


    # Loop over all cell types to extract corresponding counts, and store information in a list

    ## Initiate empty list
    counts_ls <- list()

    for (i in 1:length(cluster_names)) {

        ## Extract indexes of columns in the global matrix that match a given cluster
        column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
        ## Store corresponding sub-matrix as one element of a list
        counts_ls[[i]] <- aggr_counts[, column_idx]
        names(counts_ls)[i] <- cluster_names[i]

    }

    # Explore the different components of the list
    str(counts_ls)  


Part6: Generating matching metadata at the sample-level
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. code:: bash

    #{r}
    # Reminder: explore structure of metadata
    head(colData(sce))

    # Extract sample-level variables
    metadata <- colData(sce) %>% 
    as.data.frame() %>% 
    dplyr::select(seurat_clusters,sample_id)

    dim(metadata)
    head(metadata)

    # Exclude duplicated rows
    metadata <- metadata[!duplicated(metadata), ]

    dim(metadata)
    head(metadata)


Rename rows:

.. code:: bash

    #{r}
    rownames(metadata) <- metadata$sample_id
    head(metadata)

Number of cells per sample and cluster

.. code:: bash

    #{r}
    t <- table(colData(sce)$sample_id,
           colData(sce)$cluster_id)
    t 


.. code:: bash

    #{r}
    temp <- '11_sample_R113'
    tstrsplit(temp, "_")[[1]]
    paste(tstrsplit(temp, "_")[[2]],tstrsplit(temp, "_")[[3]],sep='_')

Creating metadata list

.. code:: bash

    #{r}
    ## Initiate empty list
    metadata_ls <- list()

    for (i in 1:length(counts_ls)) {
  
        ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
        df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
        head(df)
        ## Use tstrsplit() to separate cluster (cell type) and sample IDs
        df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
        df$sample_id  <- paste(tstrsplit(temp, "_")[[2]],tstrsplit(temp, "_")[[3]],sep='_')

    
        ## Retrieve cell count information for this cluster from global cell count table
        idx <- which(colnames(t) == unique(df$cluster_id))
        cell_counts <- t[, idx]
        ## Remove samples with zero cell contributing to the cluster
        cell_counts <- cell_counts[cell_counts > 0]
    
        ## Match order of cell_counts and sample_ids
        sample_order <- match(df$sample_id, names(cell_counts))
        cell_counts <- cell_counts[sample_order]
    
        ## Append cell_counts to data frame
        df$cell_count <- cell_counts
    
    
        ## Join data frame (capturing metadata specific to cluster) to generic metadata
        df <- plyr::join(df, metadata, 
                     by = intersect(names(df), names(metadata)))
    
        ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
        rownames(df) <- df$cluster_sample_id
    
        ## Store complete metadata for cluster i in list
        metadata_ls[[i]] <- df
        names(metadata_ls)[i] <- unique(df$cluster_id)
        }

    # Explore the different components of the list
    str(metadata_ls)


we have matching lists of counts matrices and sample-level metadata for each cell type, and we are ready to proceed with pseudobulk differential expression analysis.

.. code:: bash

    #{r}
    # Double-check that both lists have same names
    all(names(counts_ls) == names(metadata_ls))
    #counts_ls$`0`
    #counts_ls[[idx]]

In absence of 'group_id', we can assign cluster names as groups

.. code:: bash

    #{r}
    colnames(counts_ls[[1]])


merging the matrices - one for each cluster - corresponding to counts for 3 replicates - to get gene counts:

.. code:: bash

    #{r}
    merged_sm <- RowMergeSparseMatrices(counts_ls[[1]],counts_ls[[2]])

    for (i in 3:length(counts_ls)) {
        print(i)
        print(colnames(counts_ls[[i]]))
        merged_sm <- RowMergeSparseMatrices(merged_sm, counts_ls[[i]])
    }

    colnames(merged_sm)


Writing combined counts to tsv:
Note - the tedious code below, which can use an R proficient R, worksaround the structure of the dgCsparse matrix object to assign rownames as 'isoquant_id' to thee final counts table.


.. code:: bash

    #{r}
    #head(merged_sm, 3)
    write.table(as.matrix(merged_sm), 
            file ="kinnex_sc/complete_dataset/pseudo_bulk_counts.tsv",
            row.names=TRUE,
            sep="\t")

    # temp <-
    # read.table(file ="kinnex_sc/complete_dataset/pseudo_bulk_counts.tsv",
    #            sep="\t")
    # 
    # 
    # temp$isoform_id <- rownames(temp)
    # head(temp)
    # 
    # write.table(temp,
    #             file ="kinnex_sc/complete_dataset/pseudo_bulk_counts.tsv",col.names = TRUE,
    #             sep="\t")

