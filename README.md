# Ragas <img src="man/figs/Ragas.logo.png" align="right" width="300"/>

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Introduction to RAGAS
===
Jinghua Gu, Uthra Balaji
11/10/2023

Ragas (**R** **A**dvanced **G**allery for **A**nalysis of
**S**ingle-cell Data) is an R package that provides enhanced analysis
and visualization for single-cell RNA-Seq. Developed under a unique
consideration for [subcluster
analysis](https://github.com/jig4003/Ragas/blob/main/vignettes/subcluster.md),
Ragas offers functions that seamlessly integrates essential components
of the subcluster analysis. Ragas has the following unique advantages:
-   **Unified data structure**  
    All data from a variety of analysis are assembled under a unified
    structure called the “post-integration” (Pi) object
    
-   **Subcluster re-projection**  
    Objects from subclustering can be integrated and re-projected to the
    main object so that cell sub-populations from different clusters can
    be systematically analyzed side-by-side.
    
-   **Enhanced Visualization**  
    We have developed and re-developed many popular tools for single
    cell RNA-Seq, including matrix plot, annotated dot plot, stacked
    violin plot, as well as more basic box plot and bar plot to
    visualize expression and cell proportion changes acorss cluster or
    groups

Documentation
===
Read the following tutorials for a jump-start:
-   [**Tutorial: A quick start guide**](https://github.com/jig4003/Ragas/blob/main/vignettes/QuickStart.md) 
-   [**Tutorial: Ragas for subcluster analysis**](https://github.com/jig4003/Ragas/blob/main/vignettes/subcluster.md)

Installation
===
1. Install R
2. Install most of the dependencies using BiocManager::install
```r
BiocManager::install(c("Seurat", "ggplot2", "randomcoloR", "SingleCellExperiment", "Matrix", "future", "dplyr", "reshape2", "scales", "ComplexHeatmap", "grid", "circlize", "gplots", "muscat", "limma", "ggtree", "patchwork", "ggprism", "rstatix", "cowplot", "aplot", "crayon"), update = FALSE)
```
3. Install the ggsankey package from github
```r
devtools::install_github("davidsjoberg/ggsankey", upgrade = "never")
```
4. Install the Ragas package from github
```r
devtools::install_github("jig4003/Ragas", upgrade = "never")
```
5. In case users receive a version inconsistency warning about TMB (e.g., "Warning in checkDepPackageVersion(dep_pkg = "TMB"): Package version inconsistency detected...), reinstall glmmTMB
```r
install.packages("glmmTMB", type="source")
```
Docker
===
We have created a Ragas docker image based on the Jupyter docker (https://jupyter-docker-stacks.readthedocs.io/en/latest/). To access the Ragas docker, simply run the following:
```bash
docker pull uthrabalaji/ragas-docker:latest
```
