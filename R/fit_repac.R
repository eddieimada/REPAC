#' Fit model and test for differential polyA usage (DPU)
#'
#' `fit_repac` transforms the counts obtained from `create_pa_rse()` into compositions
#' and performs an isometric log ratio transformation prior to fitting the model
#'
#' @name fit_repac
#'
#' @param se A (Ranged)SummarizedExperiment object.
#' @param design A design matrix indicating the groups
#' to be compared. See [model.matrix][stats::model.matrix].
#' @param contrasts A numeric matrix with rows corresponding to coefficients in the
#' design matrix and columns containing contrasts. See [makeContrasts][limma::makeConstrasts].
#'
#' @return A list of tibbles with length equal the number of contrasts.
#' Each tibble contains the following columns:
#' gene_name - Name of the gene associated with the PAS (e.g., SYMBOL)
#' Ref - ID of the reference PAS
#' ID - ID of the PAS tested in comparison to `Ref`
#' cFC - Compositional fold-change in the ilr-space
#' mcDiff - Mean difference between the coverages in percentage across groups
#' t - T-statistic value
#' p.val - p-value
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
#'     BvsS= B - A
#' )
#' results <- fit_repac(se, dMat, cMat)
#'
#' @export
fit_repac <- function(se, design, contrasts=NULL){
    # check inputs
    dMat <- design
    cMat <- contrasts
    n <- table(SummarizedExperiment::rowData(se)$gene_name) > 1
    keep <- names(n[n == TRUE])
    se <- se[SummarizedExperiment::rowData(se)$gene_name %in% keep,]
    genes <- unique(SummarizedExperiment::rowData(se)$gene_name)
    results <- furrr::future_map_dfr(genes, function(id){
        # Get gene count + pa sites
        mat <- SummarizedExperiment::assays(se)[[1]][
            which(SummarizedExperiment::rowData(se)$gene_name == id),
            ]
        mat <- t(mat+1)
        res <- purrr::map(data.frame(combn(seq_along(colnames(mat)), 2)), function(idx) {
            nm <- colnames(mat)[idx[2]]
            ref.nms <- colnames(mat)[idx[1]]
            counts <- cbind(mat[, idx[1]], mat[, idx[2]])
            comp <- compositions::acomp(counts)
            icomp <- compositions::ilr(comp)[, 1]
            ### Extract results
            if (length(contrasts) == 0) {
                fit <- limma::lmFit(icomp, dMat)
                ref.mu <- compositions::ilrInv(mean(icomp[which(rowSums(dMat)==1)]))[2]
                mDiff <- apply(dMat[,as.logical(attributes(dMat)$assign)],2, function(l){
                    ref.mu-compositions::ilrInv(mean(icomp[as.logical(l)]))[2]
                })
            } else {
                fit <- limma::lmFit(icomp, dMat)
                fit <- limma::contrasts.fit(fit, cMat)
                mu <- apply(dMat,2, function(l){
                    compositions::ilrInv(mean(icomp[as.logical(l)]))[2]
                })
                mDiff <- apply(cMat,2,function(c){
                    sum(mu*c)
                })
            }
            fit <- limma::eBayes(fit)
            coefs <- colnames(fit$coefficients)
            tG <- map(coefs, function(x, y) {
                results <- limma::topTable(y, coef=x, n=Inf)
                names(results)[1] <- "cFC"
                results <- tibble::add_column(results, Ref=ref.nms, Site=nm, .before = 1 )
                results <- tibble::add_column(results, mcDiff=mDiff[x], .after = "cFC" )
                results <- tibble::add_column(results, Contrast=x)

            }, y=fit)
            names(tG) <- coefs
            return(tG)
        })
        res <- dplyr::bind_rows(flatten(res))
        res <- tibble::add_column(res, gene_name=gsub("_.+", "", res$Ref), .before = 1)
        res
    },.progress = TRUE)
    results <- tibble::tibble(results[-c(5,8,9)])
    results <- split(results[,-7], results$Contrast)
    return(results)
}
