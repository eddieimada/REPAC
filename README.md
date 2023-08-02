
<!-- README.md is generated from README.Rmd. Please edit that file -->

The REPAC R package implements the method described in
[Imada](https://www.biorxiv.org/content/10.1101/2022.03.14.484280v1) et
al. to test for differential Polyadenylation Usage (DPU).

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
#> alternative polyadenylation analyses with REPAC." _bioRxiv_. doi:
#> 10.1101/2022.03.14.484280 (URL:
#> https://doi.org/10.1101/2022.03.14.484280), <URL:
#> https://www.biorxiv.org/content/10.1101/2022.03.14.484280>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {REPAC: analysis of alternative polyadenylation from RNA-sequencing data},
#>     author = {Eddie L. Imada and Christopher Wilks and Ben Langmead and Luigi Marchionni},
#>     year = {2023},
#>     journal = {Genome Biology},
#>     doi = {10.1186/s13059-023-02865-5},
#>     url = {https://link.springer.com/article/10.1186/s13059-023-02865-5},
#>   }
```
