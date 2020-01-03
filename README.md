# BIOINF590_Code

Authors:  
* Stephanie The _University of Michigan_
* Elisa Warner _University of Michigan_

Final Project for BIOINF 590  

Goal is creation of a tool-kit for analyzing sparse single-cell data. The project is designed for people with single-cell data of several patients with many cells for a given cell type. The project gives various Euclidean and non-Euclidean distance matrices that can be used to compare patients based on their single-cell data. For best results, reduce the number of comparative genes so that n > p (_n_ number of cells is greater than _p_ number of gene expression values).

* `Stephanie_DE_correlation_gene_selection.R` : This R file uses Seurat DE to narrow down a set of genes to those with only 0.9 correlation. Results in $C$ number of datasets for $C$ different cell types. Each dataset contains many patients with many cells. One row of the dataset is characterized by one patient cell and its gene expression values over $k$ number of genes.
* `jblogdet.ipynb` : This function calculates the Jensen-Bregman Log Determinant Divergence for sparse matrices. Additional regularization of your sparse matrix may be necessary
* `Elisa_Similarity_Score.ipynb` : Analysis of the matrices (cells v gene expression). First, we conduct a 2-D PCA, then a 3-D PCA. Finally, we calculate distance metrics of patients against each other. The first is characterized by an $\ell_2$-norm, where we average the column expresson values of the gene expressions for each patient to get one mean expression vector for each patient. Then we calculate the $\ell_2$ distance between each patient vector to see how distanced they are from each other. Next, we try a covariance method, where we create a $k \times k$ matrix for each patient of each cell type. Then patients within a cell type are compared by the distances/similarities between their covariance matrices. We try three divergence measures: 1) Jensen-Bregman Log Determinant Divergence, 2) Affine-Invariant Riemannian Metric, 3) Log-Euclidean Riemannian Metric.
