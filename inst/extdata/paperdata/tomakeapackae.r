
gitcreds::gitcreds_set()

usethis::create_package(".")

############################

roxygen2::roxygenize()

#Rcpp::compileAttributes()


usethis::use_testthat()

library(devtools)
#test()

roxygen2::roxygenize()

devtools::check(vignettes = FALSE)

devtools::test()


devtools::check(vignettes = TRUE)


devtools::build(vignettes = FALSE)

devtools::build_vignettes(quiet = FALSE)

library(BiocCheck)
BiocCheck(checkDir = "../LimROTS paper/")
