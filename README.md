# mldr.resampling

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/mldr.resampling)](https://CRAN.R-project.org/package=mldr.resampling)
[![Downloads](http://cranlogs.r-pkg.org/badges/mldr.resampling)](https://cran.r-project.org/package=mldr.resampling)
[![TotalDownloads](http://cranlogs.r-pkg.org/badges/grand-total/mldr.resampling?color=yellow)](https://cran.r-project.org/package=mldr.resampling)

Collection of the state of the art multilabel resampling algorithms. The objective of these algorithms is to achieve balance in multilabel datasets.

## Installation

Use `install.packages` to install *mldr.resampling* and its dependencies:

```R
install.packages("mldr.resampling")
```

Alternatively, you can install it via `install_github` from the
[devtools](https://github.com/r-lib/devtools) package.

```R
devtools::install_github("madr0008/mldr.resampling")
```

## Building from source

Use `devtools::build` from [devtools](https://github.com/r-lib/devtools)
to build the package:

```R
devtools::build(args = "--compact-vignettes=gs+qpdf")
```

## Usage and examples

This package has an interface function that can be called in order to execute the desired algorithms, on the desired *mldr* datasets. This function can be called as follows:

```R
library(mldr.resampling)

resample(birds, c("MLSOL", "MLeNN"), P=30, k=5, TH=0.4)
```

For more examples and detailed explanation on available functions, please refer to the documentation.

## Citation

Please, cite *mldr.resampling* as follows:

```
@article{RIVERA2023126806,
  title = {mldr.resampling: Efficient reference implementations of multilabel resampling algorithms},
  journal = {Neurocomputing},
  pages = {126806},
  year = {2023},
  issn = {0925-2312},
  doi = {https://doi.org/10.1016/j.neucom.2023.126806},
  url = {https://www.sciencedirect.com/science/article/pii/S0925231223009293},
  author = {Antonio J. Rivera and Miguel A. Dávila and D. Elizondo and María J. {del Jesus} and Francisco Charte},
  keywords = {Multilabel learning, Imbalanced learning, Resampling algorithms, R software package},
  abstract = {Resampling algorithms are a useful approach to deal with imbalanced learning in multilabel scenarios. These methods have to deal with singularities in the multilabel data, such as the occurrence of frequent and infrequent labels in the same instance. Implementations of these methods are sometimes limited to the pseudocode provided by their authors in a paper. This Original Software Publication presents mldr.resampling, a software package that provides reference implementations for eleven multilabel resampling methods, with an emphasis on efficiency since these algorithms are usually time-consuming.}
}
```
