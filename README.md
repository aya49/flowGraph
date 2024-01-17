# flowGraph

[![DOI](https://zenodo.org/badge/DOI/10.1101/837765.svg)](https://doi.org/10.1101/837765)

flowGraph is an R package used to identify candidate biomarkers for disease diagnosis in flow cytometry data. It does so by identifying driver cell populations whose abundance changes significantly and independently given a disease.

flowGraph takes cell counts as input and outputs SpecEnr values for each cell population within a flow cytometry sample, based on their expected proportion. SpecEnr accounts for dependencies between cell populations such that we can use it to flag only cell populations whose abundance change is incurred wholly or in part because of its association with a sample class (e.g. healthy vs sick).

## Citation

The theory, proof, and algorithm behind the SpecEnr statistic used in the flowGraph package can be found in the following [paper](https://doi.org/10.1002/cyto.a.24503). Please consider citing if you found it helpful.

bibtex:
```
@article{yue2022automated,
  title={Automated identification of maximal differential cell populations in flow cytometry data},
  author={Yue, Alice and Chauve, Cedric and Libbrecht, Maxwell W and Brinkman, Ryan R},
  journal={Cytometry Part A},
  volume={101},
  number={2},
  pages={177--184},
  year={2022},
  publisher={Wiley Online Library}
}
```

The scripts and data from the paper can be downloaded on [Zenodo](https://zenodo.org/record/3991166).


## Installation

**flowGraph** can be installed via Bioconductor.

You can also intall the development version directly from Github using BiocManager:

```{r}
if (!require("BiocManager")) install.packages('BiocManager') 
BiocManager::install("aya49/flowGraph")
```

This package has been tested on R version 3.5+. If you have R version 3.5 - 4.0, download this package as a zip file, go to the `DESCRIPTION` file, change the R version requirement, and install the package locally.

## Usage

See our [vignette](vignettes/flowGraph.Rmd) for different use cases of the package: generating features, calculating summary statistics, visualizing results, and data specific use cases.

flowGraph takes as input a (sample x cell population) raw cell count matrix. See and example of how to obtain this matrix for one file (one row) [here](https://github.com/aya49/flowGraph/fcs_processing_guide.R)

flowGraph gives as output a flowGraph object that contains the SpecEnr abundance statistic and optionally its q-values. The flowGraph object can then be used to calculate more statistics, plot and visualize q-values on a cell hierarchy, etc.

```{r}
# load the package
library("flowGraph)
```
