
<!-- README.md is generated from README.Rmd. Please edit that file -->

The REPAC R/Bioconductor package implements the method described in
Imada et al. to test for differential Polyadenylation Usage (DPU).

## Instalation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then run the following code:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

devtools::install_github('eddieimada/REPAC')
```

## Citation

Below is the citation output from using `citation('REPAC')` in R. Please
run this yourself to check for any updates on how to cite **REPAC**.

``` r
print(citation('REPAC'), bibtex = TRUE)
#> 
#> Imada EL, Wilks C, Langmead B, Marchionni L (2022). "Unleashing
#> alternative polyadenylation analyses with REPAC." _bioRxiv_. <URL:
#> https://www.biorxiv.org/content/10.1101/TODO>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Unleashing alternative polyadenylation analyses with REPAC},
#>     author = {Eddie L. Imada and Christopher Wilks and Ben Langmead and Luigi Marchionni},
#>     year = {2022},
#>     journal = {bioRxiv},
#>     url = {https://www.biorxiv.org/content/10.1101/TODO},
#>   }
```
