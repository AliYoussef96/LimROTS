# LimROTS: A Hybrid Method Integrating Empirical Bayes and Reproducibility-Optimized Statistics for Robust Differential Expression Analysis

[![issues](https://img.shields.io/github/issues/AliYoussef96/LimROTS)](https://github.com/AliYoussef96/LimROTS/issues)
[![pulls](https://img.shields.io/github/issues-pr/AliYoussef96/LimROTS)](https://github.com/AliYoussef96/LimROTS/pulls)
[![R-CMD-check](https://github.com/AliYoussef96/LimROTS/workflows/rworkflows/badge.svg)](https://github.com/AliYoussef96/LimROTS/actions)
<!--[![codecov](https://codecov.io/gh/AliYoussef96/LimROTS/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/AliYoussef96/LimROTS?branch=devel)-->
<!--[![codefactor](https://www.codefactor.io/repository/github/AliYoussef96/LimROTS/badge)](https://www.codefactor.io/repository/github/AliYoussef96/LimROTS)-->


Differential expression analysis is commonly used to study diverse biological datasets. The reproducibility-optimized test statistic (ROTS) ([Elo et al., 2008](https://ieeexplore.ieee.org/document/4359873/)) uses a modified t-statistic adapted to the intrinsic characteristics of the data and ranks features by their statistical significance between two or more groups, as measured by the F-statistic. However, The ROTS publication ([Suomi et al., 2017](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005562)) does not accommodate technical or biological covariates. LimROTS ([Anwar et al., 2025](https://doi.org/10.1093/bioinformatics/btaf570)) addresses this limitation by combining a reproducibility-optimized test statistic with the limma empirical Bayes approach ([Ritchie et al., 2015](https://academic.oup.com/nar/article/43/7/e47/2414268)), enabling the analysis of more complex experimental designs. These validated solutions have been available since December 16, 2024, in the Bioconductor development version, and since April 16, 2025, in Bioconductor release 3.21. Although similar linear modeling features were later incorporated into ROTS in Bioconductor 3.21, the implementation differs, and, to our knowledge, no formal publication describing it is available. Survival analysis with covariates represents a natural extension of the linear framework in LimROTS, with a development version available in Bioconductor devel 2.3.8 for ROTS and version 1.3.17 for LimROTS.


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

LimROTS is a standalone package that extends ROTS ([Elo et al., 2008](https://ieeexplore.ieee.org/document/4359873); [Tomi Suomi et
al., 2017](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005562)) and the limma ([Ritchie ME et
al., 2015](https://academic.oup.com/nar/article/43/7/e47/2414268)) framework. 
It is a newly developed, independently implemented method that incorporates covariates in reproducibility-optimized testing. 
LimROTS is not affiliated with, endorsed by, or maintained by the original ROTS or limma developers. Users are advised to cite the original publications when referencing these methods.

ROTS: Elo LL, Filén S, Lahesmaa R, Aittokallio T. 
Reproducibility-optimized test statistic for ranking genes in microarray studies. 
IEEE/ACM Transactions on Computational Biology and Bioinformatics. 2008;5(3):423–31.

Suomi T, Seyednasrollah F, Jaakkola MK, Faux T, Elo LL (2017) ROTS: 
An R package for reproducibility-optimized statistical testing. 
PLOS Computational Biology 13(5): e1005562. https://doi.org/10.1371/journal.pcbi.1005562

limma: Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. 
limma powers differential expression analyses for RNA-sequencing and microarray studies. 
Nucleic Acids Research. 2015;43(7):e47. https://doi.org/10.1093/nar/gkv007 

## Citation

If you use `LimROTS` in your research, please cite our publication:

> Anwar, A. M., Jeba, A., Lahti, L., & Coffey, E. (2025). LimROTS: A Hybrid Method Integrating Empirical Bayes and 
Reproducibility-Optimized Statistics for Robust Differential Expression Analysis. *Bioinformatics*, btaf570. 
[https://doi.org/10.1093/bioinformatics/btaf570](https://doi.org/10.1093/bioinformatics/btaf570)