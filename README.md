# Network-Based Clustering Unveils Interconnected Landscapes of Genomic and Clinical Features Across Myeloid Malignancies

This repository provides supplementary information, data, and code associated with the manuscript:

Bayer et al. 2023, "Network-based clustering unveils interconnected landscapes of genomic and clinical features across myeloid malignancies"

https://doi.org/10.1101/2023.10.25.563992

### Online Tool

We have also developed an online tool for the classification of new patients across myeloid malignancies. Access it here:

[Myeloid Malignancies Prediction Tool](https://cbg.bsse.ethz.ch/myeloid-prediction/)

### Software

The software is available as an R package, which is hosted on CRAN. For more information on how to use this software, including installation, system requirements and application, please refer to our detailed documentation.
- CRAN Package: [Software Package](https://CRAN.R-project.org/package=clustNet)
- Installation and Example: [Code Demo](https://github.com/cbg-ethz/clustNet/blob/main/README.md)
- Documentation:  [Software Documentation (PDF)](https://cran.r-project.org/web/packages/clustNet/clustNet.pdf/)

### Reproducibility

The numerical simulations can be reproduced via the [numerical_simulations](numerical_simulations) folder. The results on the public TCGA datasets can be reproduced by running the files in the [tcga_analysis](tcga_analysis) folder. Given the full set of MDA data, all results and figures of the pan-myeloid analysis can be reproduced by running the files provied in the [analysis](analysis) folder. The cluster results from the pan-myeloid clustering analysis can be reproduced by running the code provided in the [euler_cluster](euler_cluster) folder.

![Screenshot 2023-10-31 at 08 46 43](https://github.com/cbg-ethz/myeloid-clustering/assets/38718986/b3a6c5a6-0bf1-491a-a496-e86c5dd58b76)

### License
GNU GPL (see [LICENSE](LICENSE) file for more details)
