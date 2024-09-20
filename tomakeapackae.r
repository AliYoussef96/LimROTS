
gitcreds::gitcreds_set()

usethis::create_package(".")

############################

roxygen2::roxygenize()

Rcpp::compileAttributes()


library(devtools)
#test()

roxygen2::roxygenize()

devtools::check(vignettes = FALSE)

devtools::check(vignettes = TRUE)

devtools::build(vignettes = FALSE)
devtools::build_vignettes()

library(BiocCheck)
BiocCheck()
