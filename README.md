# Agglomerative-Hierarchical-Clustering

Applying your AHC program to the NCI microarray dataset
This dataset has 64 columns and 6830 rows, where each column is an observation (a cell line), and each row represents a feature (a gene). Therefore, the dataset is represented via its transposed **data matrix**.
Implemented AHC with the following linkage functions: **single linkage, complete linkage, average linkage and centroid linkage**. Output is a data structure that represents a **dendrogram**.
Implemented a function getClusters that takes a dendrogram and a positive integer K as arguments, and its output is the **K clusters** obtained by cutting the dendrogram at an **appropriate height**
