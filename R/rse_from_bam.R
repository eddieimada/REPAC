#' Create a RangedSummarizedExperiment from bam files
#'
#' `rse_from_bam` creates a RangedSummarizedExperiment object by extracting
#'  counts from bam files of annotated polyA sites provided through an annotation.
#'
#'
#' @name rse_from_bam
#'
#' @param bams A character list with the path for the bam files.
#' @param annotation A GenomicRanges object with the coordinates of annotated
#' PA-sites. This object must have an elementMetadata column named "gene_name" with
#' gene-level IDs.
#' @param colData A data frame containing sample metadaata
#' @param ... Additional parameters passed to `Rsubread::featureCounts`
#'
#' @return A
#'  [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'  object.
#'
#' @author Eddie Imada
#'
#' @examples
#' ## Get path to bam files
#' bams <- list.file("~", full.names=T, pattern = "*.bam$")
#'
#' ## Load and quantify human hg38 PAS sites
#' data(hg38_pa)
#' pheno <- data.frame(groups=rep(c("A", "B"), each=5))
#' rse_from_bam(bams, hg38_pa, pheno, strandSpecific = 0,  allowMultiOverlap = F,
#' fraction = F, countMultiMappingReads=F, isPairedEnd = T, nthreads = 8,
#' countReadPairs = T,  requireBothEndsMapped = F)
#'
#' @export
rse_from_bam <- function(bams, annotation, colData,...){
    apa_gr <- annotation
    apa_gr <- sort(apa_gr)
    apa_gr <- c(apa_gr[strand(apa_gr) == "+"], rev(apa_gr[strand(apa_gr) == "-"]))
    apa_gr$geneID <- unlist(tapply(apa_gr$gene_name, apa_gr$gene_name,
                                   function(x){paste(x, sprintf("%02d", 1:length(x)),
                                                     sep="_")})[unique(apa_gr$gene_name)])
    apa_gr <- apa_gr[order(apa_gr$geneID)]
    names(apa_gr) <- apa_gr$geneID
    saf <- data.frame(apa_gr)[,c("geneID","seqnames", "start", "end", "strand")]
    names(saf) <- c("geneID","Chr", "Start", "End", "Strand")
    pa_quant <- Rsubread::featureCounts(files=bams, annot.ext = saf,
                                        isGTFAnnotationFile = F, ...)

    pa.counts <- pa_quant$counts
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=pa.counts),
                                                     rowRanges = apa_gr,
                                                     colData = DataFrame(colData)
    )
    return(se)
}
