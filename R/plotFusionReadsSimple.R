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
#' @param fusion Fusion object
#' @param bam.split name of bam file with junction reads in mini contig coords
#' @param bam.spanning name of bam file with spanning reads in mini contig coords
#' @param leftFlank,rightFlank how much to show left(right) of the split
#' @export
plotContigReads <- function(fusion,
                            bam.split,
                            bam.spanning, 
                            leftFlank=1000,
                            rightFlank=leftFlank) {
    ## as plotFusionReads, but with separate tracks for spanning and split
    ## 
    need <- "starfusion+fusioninspector"
    if (fusion@fusionTool != need)
      stop("Wrong object type, need: ", need, " got: ",  fusion@fusionTool)

    options(ucscChromosomeNames=FALSE)
    spanning <- addFusionReadsAlignment(fusion, bam.spanning,
                                        chromosome=fusion@geneA@chromosome)
    spanning@fusionReadsAlignment@name <- 'spanning'

    split <- addFusionReadsAlignment(fusion, bam.split,
                                        chromosome=fusion@geneA@chromosome)
    split@fusionReadsAlignment@name <- 'split'

    nspanning <- fusion@spanningReadsCount
    nsplit <- fusion@splitReadsCount

    rm(fusion)
    
    Gviz::displayPars(spanning@fusionReadsAlignment) <- list(
      showTitle = TRUE,
      showMismatches = TRUE,
      col='red',
      fill='red'
      )

    Gviz::displayPars(split@fusionReadsAlignment) <- list(
      showTitle = TRUE,
      showMismatches = TRUE,
      col='blue',
      fill='blue'
      )

    nf <- grid::grid.layout(
      nrow = 2,
      ncol = 1,
      heights = grid::unit(c(nsplit, nspanning), "null"),
      widths = grid::unit(c(1, 1), "null"))
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = nf))

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    Gviz::plotTracks(
      split@fusionReadsAlignment,
      from = split@geneA@breakpoint - leftFlank,
      to = split@geneB@breakpoint + rightFlank,
      chromosome = split@fusionReadsAlignment@chromosome,
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

