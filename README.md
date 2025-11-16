# LimROTS: A Hybrid Method Integrating Empirical Bayes and Reproducibility-Optimized Statistics for Robust Differential Expression Analysis

[![issues](https://img.shields.io/github/issues/AliYoussef96/LimROTS)](https://github.com/AliYoussef96/LimROTS/issues)
[![pulls](https://img.shields.io/github/issues-pr/AliYoussef96/LimROTS)](https://github.com/AliYoussef96/LimROTS/pulls)
[![R-CMD-check](https://github.com/AliYoussef96/LimROTS/workflows/rworkflows/badge.svg)](https://github.com/AliYoussef96/LimROTS/actions)
<!--[![codecov](https://codecov.io/gh/AliYoussef96/LimROTS/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/AliYoussef96/LimROTS?branch=devel)-->
<!--[![codefactor](https://www.codefactor.io/repository/github/AliYoussef96/LimROTS/badge)](https://www.codefactor.io/repository/github/AliYoussef96/LimROTS)-->

Differential expression analysis is a prevalent method utilised in the
examination of diverse biological data. The reproducibility-optimized test
statistic (ROTS) ([Tomi Suomi et
al.,](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005562))
has been developed with a modified t-statistic based on the data's intrinsic
characteristics and ranks features according to their statistical significance
for differential expression between two or more groups, as shown by the
f-statistic. Focusing on proteomics and metabolomics, the current ROTS
implementation cannot account for technical or biological covariates such as MS
batches or gender differences among the samples. Consequently, we developed
LimROTS, which employs a reproducibility-optimized test statistic utilizing the
limma empirical bayes ([Ritchie ME et
al.,](https://academic.oup.com/nar/article/43/7/e47/2414268)) methodology to
simulate more complex experimental designs.

## Installation instructions

### Option 1: Install from Bioconductor (recommended)

The package is available on Bioconductor release version. To install it, follow these steps,

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("LimROTS")
```

The package is also available on Bioconductor as a development (devel) version. To install it, follow these steps,

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version='devel')

BiocManager::install("LimROTS")
```

### Option 2: Install from GitHub

You can install the package directly from GitHub,

``` r
if (!requireNamespace("LimROTS", quietly = TRUE)) {
  remotes::install_github("AliYoussef96/LimROTS")
}
```
or

``` r
remotes::install_github("AliYoussef96/LimROTS" , ref  = "devel")
```

## Code of Conduct

Please note that the LimROTS project is released with a [Contributor Code of
Conduct](https://bioconductor.org/about/code-of-conduct/). By contributing to
this project, you agree to abide by its terms. Contributions are welcome in the
form of feedback, issues and pull requests. You can find the contributor
guidelines of the LimROTS
[here](https://github.com/AliYoussef96/LimROTS/blob/main/CONTRIBUTING.md).

## Acknowledgements

Please note that LimROTS was only made possible thanks to many other R and
rOpenGov software authors, which are cited in the vignettes describing this
package.

This package was developed using the following resources:

-   [*usethis*](https://cran.r-project.org/web/packages/usethis/) to generate an
    initial template.
-   Continuous code testing is performed on [GitHub
    actions](https://github.com/features/actions) and include R CMD check,
-   Code coverage assessment is possible thanks to
    [codecov](https://app.codecov.io/gh/).
-   The documentation website is automatically updated thanks to
    [*pkgdown*](https://cran.r-project.org/web/packages/pkgdown/).
-   The documentation is formatted thanks to
    [*devtools*](https://cran.r-project.org/web/packages/devtools/) and
    [*roxygen2*](https://cran.r-project.org/web/packages/roxygen2/).

## Disclaimer

LimROTS is an independent implementation inspired by the ROTS method. It is not affiliated with, 
endorsed by, or maintained by the original ROTS or limma development teams. 
This implementation was developed based on the methodologies described in:

ROTS: Elo LL, Filén S, Lahesmaa R, Aittokallio T. Reproducibility-optimized test statistic 
for ranking genes in microarray studies. IEEE/ACM Transactions on Computational Biology and Bioinformatics. 2008;5(3):423–31.

limma: Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. 
limma powers differential expression analyses for RNA-sequencing and microarray studies. 
Nucleic Acids Research. 2015;43(7):e47. https://doi.org/10.1093/nar/gkv007

Users are advised to cite the original publications when referencing the ROTS or limma methods.

## Citation

If you use `LimROTS` in your research, please cite our publication:

> Anwar, A. M., Jeba, A., Lahti, L., & Coffey, E. (2025). LimROTS: A Hybrid Method Integrating Empirical Bayes and 
Reproducibility-Optimized Statistics for Robust Differential Expression Analysis. *Bioinformatics*, btaf570. 
[https://doi.org/10.1093/bioinformatics/btaf570](https://doi.org/10.1093/bioinformatics/btaf570)