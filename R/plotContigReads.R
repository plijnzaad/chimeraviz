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
    if (fusion@fusion_tool != need)
      stop("Wrong object type, need: ", need, " got: ",  fusion@fusion_tool)

    options(ucscChromosomeNames=FALSE)
    spanning <- add_fusion_reads_alignment(fusion, bam.spanning,
                                        chromosome=fusion@gene_upstream@chromosome)
    spanning@fusion_reads_alignment@name <- 'spanning'

    split <- add_fusion_reads_alignment(fusion, bam.split,
                                        chromosome=fusion@gene_upstream@chromosome)
    split@fusion_reads_alignment@name <- 'split'

    nspanning <- fusion@spanning_reads_count
    nsplit <- fusion@split_reads_count

    rm(fusion)
    
    Gviz::displayPars(spanning@fusion_reads_alignment) <- list(
      showTitle = TRUE,
      showMismatches = TRUE,
      col='red',
      fill='red'
      )

    Gviz::displayPars(split@fusion_reads_alignment) <- list(
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
      split@fusion_reads_alignment,
      from = split@gene_upstream@breakpoint - leftFlank,
      to = split@gene_downstream@breakpoint + rightFlank,
      chromosome = split@fusion_reads_alignment@chromosome,
      type = "pileup",
      add = TRUE)
    grid::popViewport(1)

    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    Gviz::plotTracks(
      spanning@fusion_reads_alignment,
      from = spanning@gene_upstream@breakpoint - leftFlank,
      to = spanning@gene_downstream@breakpoint + rightFlank,
      chromosome = spanning@fusion_reads_alignment@chromosome,
      type = "pileup",
      add = TRUE)
    grid::popViewport(1)
}                                       #plotContigReads

