#' Obtain an isometric log ratio transformed matrix from REPAC results
#'
#' `ilr_counts` generates an sometric log ratio transformed matrix based on the
#'  References and Sites from a REPAC result.
#'
#'
#' @name ilr_counts
#'
#' @param se The (Ranged)SummarizedExperiment object used on the analysis.
#' @param results A data.frame or tibble with results from REPAC

#' @return An isometric log ratio transformed matrix.
#'
#'
#' @author Eddie Imada
#'
#' @examples
#' ## fit model and test for DPU
#' groups <- rep(c("A", "B"), each=5)
#' dMat <- model.matrix(~ 0+groups)
#' colnames(dMat) <-  gsub("groups", "",colnames(dMat))
#' cMat <- makeContrasts(
#'     levels=colnames(dMat),
#'     BvsA= B - A
#' )
#' results <- fit_repac(se, dMat, cMat)
#' results <- resolve_comparisons(results)
#' mat <- ilr_counts(se, results$BvsA)
#' @export
ilr_counts <- function(se, results) {
    pa <- SummarizedExperiment::assays(se)[[1]]
    mat <- matrix(nrow=nrow(results), ncol = ncol(pa))
    for (i in 1:nrow(results)){
        ref <-results[,"Ref"][i, ,drop=T]
        site <-results[,"Site"][i, ,drop=T]
        mat[i,] <- compositions::ilr(t(pa[c(ref,site),])+1)[,1]
    }
    rownames(mat) <- results[,"Site", drop=T]
    return(mat)
}
