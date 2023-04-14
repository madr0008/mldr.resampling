# mldr.resampling

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

# Data sets birds, emotions and genbase are
# provided within the package mldr
resample(bibtex, c("MLSOL", "MLeNN"), P=30, k=5, TH=0.4)
```

For more examples and detailed explanation on available functions, please refer to the documentation.
