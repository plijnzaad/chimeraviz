#' Simple version of plotFusionReads
#'
#' @export
plotFusionReadsSimple <- function(fusion, leftFlank=10000, rightFlank=leftFlank) {
  Gviz::displayPars(fusion@fusion_reads_alignment) <- list(
    showTitle = FALSE,
    showMismatches = FALSE # Show mismatched reads?
  )

  nf <- grid::grid.layout(
    nrow = 1,
    ncol = 1,
    heights = grid::unit(c(6), "null"),
    widths = grid::unit(c(1), "null"))
  
  grid::grid.newpage()
  Gviz::plotTracks(
    fusion@fusion_reads_alignment,
    from = fusion@gene_upstream@breakpoint-leftFlank,
    to = fusion@gene_downstream@breakpoint+rightFlank,
    chromosome = fusion@fusion_reads_alignment@chromosome,
    type = "pileup",
    add = TRUE)
}                                       #plotFusionReadsSimple
