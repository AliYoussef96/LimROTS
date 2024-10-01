# LimROTS
[![issues](https://img.shields.io/github/issues/AliYoussef96/LimROTS)](https://github.com/AliYoussef96/LimROTS/issues)
[![pulls](https://img.shields.io/github/issues-pr/AliYoussef96/LimROTS)](https://github.com/AliYoussef96/LimROTS/pulls)
[![R-CMD-check](https://github.com/AliYoussef96/LimROTS/workflows/rworkflows/badge.svg)](https://github.com/AliYoussef96/LimROTS/actions)
<!--[![codecov](https://codecov.io/gh/AliYoussef96/LimROTS/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/AliYoussef96/LimROTS?branch=devel)-->
<!--[![codefactor](https://www.codefactor.io/repository/github/AliYoussef96/LimROTS/badge)](https://www.codefactor.io/repository/github/AliYoussef96/LimROTS)-->

LimROTS: A Hybrid Method Integrating Empirical Bayes and Reproducibility-Optimized Statistics for Robust Analysis of Protein and Metabolite Data

Differential expression analysis is a prevalent method utilised in the examination of diverse biological data.
The reproducibility-optimized test statistic (ROTS) ([Tomi Suomi et al.,](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005562)) has been developed with a modified t-statistic based on the data's intrinsic characteristics and ranks features according to their statistical significance for differential expression between two or more groups, as shown by the f-statistic. Focusing on proteomics and metabolomics, the current ROTS implementation cannot account for technical or biological covariates such as MS batches or gender differences among the samples. Consequently, we developed LimROTS, which employs a reproducibility-optimized test statistic utilizing the limma ([Ritchie ME et al.,](https://academic.oup.com/nar/article/43/7/e47/2414268)) methodology to simulate more complex experimental designs.


## Installation instructions

The devel version of LimROTS can be installed from GitHub as follows:

``` r
# Install LimROTS if not already installed
if (!requireNamespace("LimROTS", quietly = TRUE)) {
  remotes::install_github("AliYoussef96/LimROTS")
}
```

``` r
remotes::install_github("AliYoussef96/LimROTS")
```

## Code of Conduct

Please note that the LimROTS project is released with a
[Contributor Code of Conduct](https://bioconductor.org/about/code-of-conduct/).
By contributing to this project, you agree to abide by its terms. Contributions
are welcome in the form of feedback, issues and pull requests. You can find the
contributor guidelines of the LimROTS
[here](https://github.com/AliYoussef96/LimROTS/blob/main/CONTRIBUTING.md).

## Acknowledgements
Please note that LimROTS was only made possible thanks to many other R and
rOpenGov software authors, which are cited in the vignettes describing
this package.

This package was developed using the following resources:

- [_usethis_](https://cran.r-project.org/web/packages/usethis/) to generate an
  initial template.
- Continuous code testing is performed on
  [GitHub actions](https://github.com/features/actions) and include R CMD check,
- Code coverage assessment is possible thanks to
  [codecov](https://app.codecov.io/gh/).
- The documentation website is automatically updated thanks to
  [_pkgdown_](https://cran.r-project.org/web/packages/pkgdown/).
- The documentation is formatted thanks to
  [_devtools_](https://cran.r-project.org/web/packages/devtools/) and
  [_roxygen2_](https://cran.r-project.org/web/packages/roxygen2/).

