## Escape Data Correlations and Principal Component Analysis

Using our .csv escape data from [Individual_Antibody_Escape_Profiles](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Individual_Antibody_Escape_Profiles), we generated Pearson coefficients and a correlogram for all pairs of antibodies using the R script, [correlogram_and_PCA.R](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/blob/main/Correlogram_and_PCA/correlogram_and_PCA.R).

[correlogram_and_PCA.R](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/blob/main/Correlogram_and_PCA/correlogram_and_PCA.R) also runs a principal component analysis (PCA). To ensure appropriate clustering, we ran elbow, silhouette and gap analysis, and picked the optimal *k* based on the results from these. We chose *k* = 12 because it was the result given by the gap method, and it fell between the numbers suggested by the elbow (7) and silhouette (19) methods.

To run [correlogram_and_PCA.R](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/blob/main/Correlogram_and_PCA/correlogram_and_PCA.R) with our .csv escape data files, we also required [reference.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/blob/main/Correlogram_and_PCA/reference.csv) file, which is also included here. I also wished to order my correlogram in a particular way, so included a predefined order, [corr_order.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/blob/main/Correlogram_and_PCA/corr_order.csv).

Output is given in the [Output](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Correlogram_and_PCA/Output) subfolder.
