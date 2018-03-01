#' Simple version of plotFusionReads
#'
#' @export
plotFusionReadsSimple <- function(fusion, leftFlank=10000, rightFlank=leftFlank) {
  Gviz::displayPars(fusion@fusionReadsAlignment) <- list(
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
    fusion@fusionReadsAlignment,
    from = fusion@geneA@breakpoint-leftFlank,
    to = fusion@geneB@breakpoint+rightFlank,
    chromosome = fusion@fusionReadsAlignment@chromosome,
    type = "pileup",
    add = TRUE)
}                                       #plotFusionReadsSimple


#' Plot fusion with separate junction and spanning reads
#'
#' @export
plotContigReads <- function(fusion, bam.junction, bam.spanning, 
                            leftFlank=1000, rightFlank=leftFlank) {
    ## as plotFusionReads, but with separate tracks for spanning and junction
    need <- "starfusion+fusioninspector"
    if (fusion@fusionTool != need)
      stop("Wrong object type, need: ", need, " got: ",  fusion@fusionTool)

    options(ucscChromosomeNames=FALSE)
    spanning <- addFusionReadsAlignment(fusion, bam.spanning,
                                        chromosome=fusion@geneA@chromosome)
    spanning@fusionReadsAlignment@name <- 'spanning'
    junction <- addFusionReadsAlignment(fusion, bam.junction,
                                        chromosome=fusion@geneA@chromosome)
    junction@fusionReadsAlignment@name <- 'junction'
    
    Gviz::displayPars(spanning@fusionReadsAlignment) <- list(
      showTitle = TRUE,
      showMismatches = FALSE,
      col='red',
      fill='red'
      )

    Gviz::displayPars(junction@fusionReadsAlignment) <- list(
      showTitle = TRUE,
      showMismatches = FALSE,
      col='blue',
      fill='blue'
      )
    
    nf <- grid::grid.layout(
      nrow = 2,
      ncol = 1,
      heights = grid::unit(c(2, 2), "null"),
      widths = grid::unit(c(1, 1), "null"))
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = nf))

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    Gviz::plotTracks(
      junction@fusionReadsAlignment,
      from = junction@geneA@breakpoint - leftFlank,
      to = junction@geneB@breakpoint + rightFlank,
      chromosome = junction@fusionReadsAlignment@chromosome,
      type = "pileup",
      add = TRUE)
    grid::popViewport(1)

    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    Gviz::plotTracks(
      spanning@fusionReadsAlignment,
      from = spanning@geneA@breakpoint - leftFlank,
      to = spanning@geneB@breakpoint + rightFlank,
      chromosome = spanning@fusionReadsAlignment@chromosome,
      type = "pileup",
      add = TRUE)
    grid::popViewport(1)
}                                       #plotContigReads

