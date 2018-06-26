#' read the results from runnning FusionInspector
#'
#' This function is called when \code{\link{import_star_fusion}} is called
#' with \code{use_fusion_inspector=TRUE}. It will override coordinates
#' (chromosome, strand, and breakpoint) with those corresponding to
#' FusionInspector's "mini genome". The resulting list of fusion objects
#' can currently only be plotted with \code{\link{plotFusionReadsSimple}}.
#'
#' Note: there may not be perfect match between the fusions found by
#' STAR-Fusion and FusionInspector (since the latter reruns STAR on the
#' mini contig; also, selfie fusions are ignored).  If no match can be
#' found, the whole fusion is skipped. 
#'
#' @param filename Which file to import. Note that 
#'
#' @return The report is returned as \code{data.frame}
#'
#' @export
import_fusion_inspector <- function (filename='FusionInspector-inspect/finspector.igv.FusionJuncSpan',
                                     limit) {
    ## Try to read the FusionInspector report.
    ## Within one 'scaffold', only keep the leftmost and rightmost coordinate
    report <- withCallingHandlers({
        col_types_fusioninspector = readr::cols_only(
          "#scaffold" = col_character(),
          "fusion_break_name" = col_skip(),
          "break_left" = col_integer(),
          "break_right" = col_integer(),
          "num_junction_reads" = col_skip(),
          "num_spanning_frags" = col_skip(),
          "spanning_frag_coords" = col_skip()
          )
        if (missing(limit)) {
            readr::read_tsv(
              file = filename,
              col_types = col_types_fusioninspector)
        } else {
            readr::read_tsv(
              file = filename,
              col_types = col_types_fusioninspector,
              n_max = limit)
        }
    },
      error = function(cond) {
          message(paste0("Reading ", filename, " caused an error: ", cond[[1]]))
          stop(cond)
      },
      warning = function(cond) {
          message(paste0("Reading ", filename, " caused a warning: ", cond[[1]]))
          warning(cond)
      })
    colnames(report)[1] <- 'id'
    d <-
      data.frame(id=with(report, id[ !duplicated(id) ]),
                 start=with(report, tapply(break_left, id, min)),
                 end=with(report, tapply(break_right, id, max)))
    rownames(d) <- d$id
    d
}                                        #import_fusion_inspector
