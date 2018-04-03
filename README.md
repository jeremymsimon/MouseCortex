# MouseCortex
Welcome to the GitHub repository associated with Loo et al., "Single-cell transcriptomic catalog of mouse cortical development", bioRxiv 2017: https://www.biorxiv.org/content/early/2017/10/25/208744

## Abstract
We generated a single-cell transcriptomic catalog of the developing mouse cerebral cortex that includes numerous classes of neurons, progenitors, and glia, their proliferation, migration, and activation states, and their relatedness within and across timepoints. Cell expression profiles stratified neurological disease-associated genes into distinct subtypes. Complex neurodevelopmental processes can be reconstructed with single-cell transcriptomics data, permitting a deeper understanding of cortical development and the cellular origins of brain diseases.

## What you'll find here
The code maintained here includes additions and modifications to Shekhar et al.'s class.R file published here: https://github.com/broadinstitute/BipolarCell2016 that we used to identify cell types of the developing mouse cortex

The key steps are implemented in R, and are summarized and documented in the main R markdown file `E14_processing.R`. This requires the user to load `class.R` to access a wide variety of functions for data analysis. 

We also include code for visualizing the data using t-SNE, in Python: `E14_tSNE.py`

We also created web-based tools for visualizing our data, these can be found at:
http://zylkalab.org/data


## Change log compared to Shekhar version:
* Added several new slots to S4 object, including `pca.eigenvalues`, `sils`, and `numclust`
* Added new function called `perform_refined_clustering`, which iteratively runs `doGraph_clustering` to optimize the number of nearest neighbors, then iteratively runs `doGraph_clustering` again to refine the cluster assignments using computed Silhouette widths
* Computes eigenvalues as part of `doPCA` and stores these values for the permutation test
* Adds a legend to tSNE plots
* Modifies the `binomcount.test` function to compute log-fold-change as `(x+1)/(y+1)` to avoid NA/Inf values
* Adds checks to the `binomcount.test` to ensure that `effect.size` isn't NA
* Modifies cluster naming on `dot.plot`
* Prevents the usage of "Dingbats" font class in `ggplot` calls

