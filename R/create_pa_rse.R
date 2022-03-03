#' Create a RangedSummarizedExperiment from recount3
#'
#' `create_pa_rse creates a RangedSummarizedExperiment object by extracting
#'  expression information from recount3 bigWig files from annotated polyA
#'  sites provided through a BED file. It will also pull relevant metadata
#'  for the project through recount3.
#'
#' @name create_pa_rse
#'
#' @param organism A character(1) specifying which organism you want to download
#' data from. Supported options are "human" or "mouse".
#' @param project A character(1) with the ID for a given study. See
#' [available_projects][recount3::available_projects] from the recount3 package.
#' @param annotation A character(1) with the path to a BED file containing
#' the coordinates of annotated PA-sites.
#' @param sample_id A character vector with the ID(s) for a given sample inside
#' the study. See
#' [available_samples][recount3::available_samples] from the recount3 package.
#' @param bed_cols A character vector with the column names for the BED file.
#' By default it uses the standard fields for a BED file ("chr", "start",
#' "end", "ID", "score", "strand").
#'
#' @return A
#'  [RangedSummarizedExperiment-class][SummarizedExperiment::RangedSummarizedExperiment-class]
#'  object.
#'
#' @author Eddie Imada
#'
#' @examples
#' ## Explore available projects in recount3
#' human_projects <- recount3::available_projects()
#'
#' ## Find the project you are interested in and obtain PA expression
#' create_pa_rse(organism = "mouse", project="SRP048707", annotation="mysites.bed")
#'
#' @export
create_pa_rse <- function(organism=c("human", "mouse"), project=NULL, annotation=NULL, sample_id=NULL, bed_cols= c("chr", "start", "end", "ID", "score", "strand")) {
    print("Retrieving metadata...")
    human_samples <- recount3::available_samples(organism = organism)
    metadata <- purrr::map_dfr(project, function(id){
    samples <- human_samples[human_samples$project == id,]$external_id
    home <- human_samples[human_samples$project == id, "project_home"][1]
    metadata <- recount3::locate_url(id,
                           project_home = home,
                           organism = organism,
                           type="metadata",
                           sample = sample)
    metadata <- recount3::file_retrieve(metadata)
    metadata <- recount3::read_metadata(metadata)
    metadata$BigWigURL <- recount3::locate_url(id, project_home = home, type="bw", sample = metadata$external_id, organism = organism)
    return(metadata)
    })
    keep <- purrr::map_lgl(metadata, ~ all(!is.na(.x)))
    metadata <- metadata[,keep]
    if (length(sample_id) == 0) {
        urls <- metadata$BigWigURL
    } else {
        metadata <- metadata[metadata$external_id %in% sample_id,]
        urls <- metadata$BigWigURL
    }
    print("Loading PA sites...")
    apa_gr <- GenomicRanges::makeGRangesFromDataFrame(readr::read_tsv(annotation, col_names = bed_cols), keep.extra.columns = T)
    names(apa_gr) <- paste0(seqnames(apa_gr), ":", start(apa_gr), "-", end(apa_gr))
    ord <- apa_gr
    strand(ord) <- "*"
    apa_gr <- apa_gr[names(sort(ord))]
    print("Retrieving counts...")
    mat <- furrr::future_map_dfc(seq_along(urls), function(x) {
        y <- megadepth::get_coverage(urls[[x]], op = "sum",
                          annotation = annotation,
                          prefix=file.path(tempdir(), basename(urls[[x]]))
        )
        y <- sort(y)
        y$score

    },.progress = T)
    mat <- as.matrix(mat)
    mat <- round(sweep(mat, 2,
                       STATS=metadata$recount_qc.star.average_mapped_length,
                       FUN = "/"))
    se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=mat),
                                colData = DataFrame(metadata),
                                rowRanges = apa_gr
    )
    rownames(se) <- se$ids
    se <- se[sort(rownames(se))]
}
