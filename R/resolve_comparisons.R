#' Resolve which comparisons to keep based on the selected strategy
#'
#' The `fit_repac()` function will perform all pairwise comparisons among the PAS of a
#' gene. This is done to mantain flexibility to decide how the comparisons among PAS
#' are made.
#' `resolve_comparisons()` will resolve the comparisons according to the strategy selected
#' and output a reduced set of results.
#'
#' The "significance"is the default strategy to resolve the comparisons. This strategy will
#' simply keep the most significant comparisons for each reference site. The "closest"
#' strategy updates the reference site to the current site when these sites have
#' significantly differences in compositions. Finally, the "longest" strategy is the
#' strategy described in the REPAC manuscript, in which the reference is kept until a
#' it no longer has significant differences in compositions with the site being tested,
#' at which point it is updated to the current site.
#'
#'
#' @name resolve_comparisons
#'
#' @param data A list of tibbles outputed by `fit_repac()`.
#' @param threshold Significance cut-off for the "closest" and "longest" strategies.
#' @param resolve.by The strategy to be applied.
#'
#'
#' @return A tibble with adjusted p-values.
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
#' results <- resolve_comparisons(results)
#'
#' @export
resolve_comparisons <- function(data, threshold=0.05, resolve.by = c("significance", "closest", "longest")) {
    method <- match.arg(resolve.by)
    purrr::map(data, function(results) {
    results$Ref.idx <- gsub(".+_", "", results$Ref)
    results$Site.idx <- gsub(".+_", "", results$Site)
    l <- split(results, results$gene_name)
    resolved <- purrr::map_dfr(l, function(gene){
        idxs <- unique(gene$Site.idx)
        ref <- gene$Ref.idx[1]
        if (method == "significance") {
            tmp <- purrr::map_dfr(unique(gene$Ref), function(idx){
                sel <- gene[gene$Ref == idx,]
                sel <- sel[which.min(sel$P.Value),]
                sel[,-c(8,9)]
            })
            tmp$adjusted.P.Value <- p.adjust(tmp$P.Value)
        }
        else if (method == "closest") {
            tmp <- purrr::map_dfr(idxs, function(idx){
                site <- dplyr::filter(gene, Ref.idx == ref & Site.idx == idx)
                if (site$P.Value <= threshold) {
                    ref <<- idx
                    site[,-c(8,9)]
                } else {
                    site[,-c(8,9)]
                }
            })
            tmp$adjusted.P.Value <- p.adjust(tmp$P.Value)
        }
        else {
            tmp <- purrr::map_dfr(idxs, function(idx){
                site <- dplyr::filter(gene, Ref.idx == ref & Site.idx == idx)
                if (site$P.Value >= threshold) {
                    ref <<- idx
                    site[,-c(8,9)]
                } else {
                    site[,-c(8,9)]
                }
            })
            tmp$adjusted.P.Value <- p.adjust(tmp$P.Value)
        }
        return(tmp)
    })
    return(resolved)
    })
}
