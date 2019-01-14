# MouseCortex
Welcome to the GitHub repository associated with "**Single-cell transcriptomic analysis of mouse neocortical development**"  
**Loo, Simon _et al._**, _Nature Communications_, 2019:  
https://www.nature.com/articles/s41467-018-08079-9

## Abstract
The development of the mammalian cerebral cortex depends on careful orchestration of proliferation, maturation, and migration events, ultimately giving rise to a wide variety of neuronal and non-neuronal cell types. To better understand cellular and molecular processes that unfold during late corticogenesis, we perform single-cell RNA-seq on the mouse cerebral cortex at a progenitor driven phase (embryonic day 14.5) and at birth—after neurons from all six cortical layers are born. We identify numerous classes of neurons, progenitors, and glia, their proliferative, migratory, and activation states, and their relatedness within and across age. Using the cell-type-specific expression patterns of genes mutated in neurological and psychiatric diseases, we identify putative disease subtypes that associate with clinical phenotypes. Our study reveals the cellular template of a complex neurodevelopmental process, and provides a window into the cellular origins of brain diseases.

## What you'll find here
* `class.R`: Additions and modifications to **Shekhar _et al._'s** `class.R` file inititally published [here](https://github.com/broadinstitute/BipolarCell2016) 
* `E14_processing.R`: A detailed, fully documented summary of the processing steps performed for the E14.5 cortex.
	* Requires `class.R` to access the functions necessary for data analysis
	* All processing steps and parameters for analysis of the P0 cortex were identical to those illustrated in this template
* `E14_combined_matrix.txt.gz`: The main gene expression matrix for E14.5 cortex.
	* Rows are genes, columns are cells. 
	* Each biological replicate is represented here.
* `E14_tSNE.py`: Code for visualizing the data using `t-SNE`, in Python
* `P0_combined_matrix.txt.gz`: The main gene expression matrix for P0 cortex.
	* Rows are genes, columns are cells.
	* Each biological replicate is represented here.

We also created web-based tools for visualizing our data, these can be found at:
http://zylkalab.org/data

## Required R libraries
* genefilter
* sva
* igraph
* ggplot2
* Matrix
* gmodels
* RANN
* reshape
* cluster
* SLICER
* gplots

## Required Python libraries (for t-SNE visualization)
* scikit-learn
* matplotlib
* pandas

## Change log compared to Shekhar version:
* Added several new slots to S4 object, including `pca.eigenvalues`, `sils`, and `numclust`
* Added new function called `perform_refined_clustering`, which iteratively runs `doGraph_clustering` to optimize the number of nearest neighbors, then iteratively runs `doGraph_clustering` again to refine the cluster assignments using computed Silhouette widths
* Computes eigenvalues as part of `doPCA` and stores these values for the permutation test
* Adds a legend to t-SNE plots
* Modifies the `binomcount.test` function to compute log-fold-change as `(x+1)/(y+1)` to avoid NA/Inf values
* Adds checks to the `binomcount.test` to ensure that `effect.size` isn't NA
* Modifies cluster naming on `dot.plot`
* Prevents the usage of "Dingbats" font class in `ggplot` calls

## Questions/Issues
Please direct any questions or issues with the code published here to `jeremy [underscore] simon (at) med [dot] unc [dot] edu`
