---
title: "S2. TGFbeta Timecourse C1 CAGE Analysis: Bulk Data Analysis"
author: "Andrew T. Kwon (taejun.kwon@riken.jp)"
date: "May 30, 2018"
output: 
  html_document: 
    keep_md: yes
    number_sections: yes
---



# Overview

Complementing the C1 CAGE data, we also have the corresponding bulk nAnT-iCAGE data. This will be used in later stages to aid the C1 CAGE data analysis. We process the bulk data using the edgeR package and perform differential expression.

---------

# Data Preparation

## Library loading and parameter set up


```r
# Load the required libraries

library(tidyverse)
library(edgeR)

# load the s1 data for annotation, and the prepared raw counts data
load("~/Projects/Molecular_Network/TGFbeta_Timecourse/results/R_sessions/s1_scran.RData")
bulk.exp <- readRDS(file.path(dirs$results, "R_sessions/GENCODEv25.cage_cluster.coord.mask.dpi_sum.Rds"))
```

## Helper functions for visualizing bulk data expression.

```r
dotplot_expression <- function(ids, dat=dge, log=TRUE, id.type='clusterID')
{
  if (id.type == 'clusterName') {
    ids <- dplyr::filter(annot, clusterName %in% ids)$clusterID
  }
  
  gtab <- dplyr::filter(annot, clusterID %in% ids) %>% dplyr::group_by(geneNameStr) %>% dplyr::arrange(clusterName, .by_group=TRUE)
  ids <- gtab$clusterID
  ids <- ids[ids %in% rownames(dat)]
  
  if (length(ids) > 0)
  {
    ylim <- c(min(edgeR::cpm(dge, log=log)[ids,]), 
              max(edgeR::cpm(dge, log=log)[ids,]))
    plots <- map(ids, function(id) {
      tab <- data.frame(Exp=edgeR::cpm(dge, log=log)[id,], Timepoint=rep(c('t00','t06','t24'), each=3))
      tab$Time <- factor(tab$Timepoint)
      title <- dplyr::filter(annot, clusterID == id)$clusterName
      ggplot(tab, aes(x=Timepoint, y=Exp, color=Timepoint, fill=Timepoint)) + geom_point(size=3)  + ylim(ylim) + fontsize + ggtitle(title, subtitle=id) + scale_color_manual(values=sample.colors$Timepoint[tab$Timepoint])
    })
    cols <- ifelse(length(plots) >= 6, 3, ifelse(length(plots) < 6 & length(plots) > 1, 2, 1))
    scater::multiplot(plotlist=plots, cols=cols)
  } else {
    print('No CAGE clusters found')
  }
}

dotplot_expression_by_genes <- function(geneSymbols, dat=dge, log=TRUE)
{    
    ids <- dplyr::filter(annot, geneNameStr %in% geneSymbols)$clusterID
    dotplot_expression(ids, dat=dat, log=log)
}
```

---------

# Main Analysis

We filter lowly expressed promoters and perform differential expression analysis.


```r
# from the object, get the raw counts
n <- bulk.exp$dpi
bulk.exp <- as.matrix(bulk.exp[-1,2:10])
rownames(bulk.exp) <- n[2:length(n)]
groups <- factor(rep(c('t00','t06','t24'), each=3))
dge <- DGEList(counts=bulk.exp, group=groups)

# Filtering for lowly expressed tag clusters
keepRows <- rowSums(cpm(dge) > 1) >= 3 & rowSums(bulk.exp > 10) >= 1
dge <- dge[keepRows,] # 28031 rows left
dge$samples$lib.size <- colSums(dge$counts)

# normalization and common dispersion calculation
# biological coefficient of variation: should be around 0.4
dge <- calcNormFactors(dge)

design <- model.matrix(~groups)
dge <- estimateDisp(dge, design)

rm(bulk.exp, n, groups, keepRows, design)
```

