#' Fit model and test for differential polyA usage (DPU)
#'
#' `fit_repac` transforms the counts obtained from `create_pa_rse()` into compositions
#' and performs an isometric log ratio transformation prior to fitting the model and
#' testing via ANOVA.
#'
#' @name fit_repac
#'
#' @param se A (Ranged)SummarizedExperiment object.
#' @param gene_name A character vector containing a gene ID to be assigned to each
#' PA site.
#' @param group A character(1) with the name of the variable indicating the groups
#' to be compared. The variable must be a `colData()` column of the rse.
#' @param covariates A character vector with the name of the variables indicating
#' the co-variates to be included in the model. The variables must be a `colData()`
#' column of the rse.
#' @param method A character(1) with the method to be used to adjust the p-values for
#' multiple hypothesis testing. See [p.adjust][stats::p.adjust]
#'
#' @return A tibble with the following columns:
#' gene_name - Name of the gene associated with the PAS (e.g., SYMBOL)
#' Ref - ID of the reference PAS
#' ID - ID of the PAS tested in comparison to `Ref`
#' cFC - Compositional fold-change in the ilr-space
#' mDiff - Mean difference between the coverages in percentage across groups
#' t - T-statistic value
#' p.val - p-value
#' adj.p.val - Adjusted p-value
#'
#' @author Eddie Imada
#'
#' @examples
#' ## fit model and test for DPU
#' results <- fit_repac(se, group="groups", covariates = NULL)
#'
#' #' ## fit model and test for DPU, correcting for batch effects
#' results <- fit_repac(se, group="groups", covariates = "batches")
#'
#' @export
fit_repac <- function(se, group, covariates=NULL, method="BH"){
    # check inputs
    stopifnot(
        "'group' must be a variable in colData()" =
            group %in% names(colData(se)),
        "'covariates' must be a variable in colData()" =
            all(group %in% names(colData(se)))
    )
    n <- table(rowData(se)$gene_name) > 1
    keep <- names(n[n == TRUE])
    se <- se[rowData(se)$gene_name %in% keep,]
    genes <- unique(rowData(se)$gene_name)
    results <- furrr::future_map_dfr(genes, function(id){
        # Get gene count + pa sites
        mat <- SummarizedExperiment::assays(se)[[1]][which(rowData(se)$gene_name == id),]
        mat <- t(mat+1)
        ### Select reference pa site
        ref <- 1
        cd <- SummarizedExperiment::colData(se)
        res <- purrr::map_dfr(seq_along(colnames(mat)[-1]), function(idx){
            nm <- colnames(mat)[idx+1]
            ref.nms <- colnames(mat)[ref]
            #weight <- gene.weight[colnames(mat)[ref],]
            counts <- cbind(mat[,ref],mat[,idx+1])
            comp <- compositions::acomp(counts)
            icomp <- compositions::ilr(comp)[,1]
            if (length(covariates) != 0) {
                fmla <- reformulate(c(covariates, group), response = "icomp")
                fit <- lm(fmla, data=model.frame(fmla, cd))
            } else {
                fmla <- reformulate(group, response = "icomp")
                fit <- lm(fmla, data=model.frame(fmla, cd))
            }
            mod <- car::Anova(fit)
            p.val <- mod[[group,"Pr(>F)"]]
            if (p.val > 0.1) {
                ref <<- idx+1
            }
            ci <- grep(group, names(fit$coefficients))
            du <- compositions::ilrInv(fit$coefficients[1] + fit$coefficients[ci])[2] - ilrInv(fit$coefficients[1])[2]
            cFC <- fit$coefficients[ci]
            Ts <- round(summary(fit)$coefficients[[ci, "t value"]], 3)
            tibble::tibble(gene_name=id, Ref=ref.nms, ID=nm, cFC=cFC, mDiff=du, "t"=Ts, p.val=p.val)
        })
        res$adj.p.val <- p.adjust(res$p.val, method = method)
        res
    },.progress = TRUE)

    return(results)
}
