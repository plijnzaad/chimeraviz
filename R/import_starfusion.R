#' Import results from a STAR-Fusion run into a list of Fusion objects.
#'
#' A function that imports the results from a STAR-Fusion run, typically from
#' a star-fusion.fusion_candidates.final.abridged file, into a list of Fusion
#' objects.
#'
#' @param filename Filename for the STAR-Fusion
#' star-fusion.fusion_candidates.final.abridged results file.
#' @param genome_version Which genome was used in mapping (hg19, hg38, etc.).
#' @param limit A limit on how many lines to read.
#' @param use_fusioninspector If STAR-Fusion was run with the FusionInpspector
#'    option, specifying \code{FALSE} will read the fusion objects
#'    from the 'mini genome' it creates. See \code{\link{import_fusioninspector}} and
#'    \code{\link{plotFusionReadsSimple}}
#' @return A list of Fusion objects.
#'
#' @examples
#' starfusionData <- system.file(
#'   "extdata",
#'   "star-fusion.fusion_candidates.final.abridged.txt",
#'   package = "chimeraviz")
#' fusions <- import_starfusion(starfusionData, "hg19", 3)
#' # This should import a list of 3 fusions described in Fusion objects.
#'
#' @export
import_starfusion <- function (filename, genome_version, limit,
                               use_fusioninspector=FALSE) {

  # Is the genome version valid?
  valid_genomes <- c("hg19", "hg38", "mm10")
  if (is.na(match(tolower(genome_version), tolower(valid_genomes)))) {
    stop("Invalid genome version given")
  }

  # If the limit is set, is the value valid?
  if (missing(limit) == FALSE) {
    if (is.numeric(limit) == FALSE || limit <= 0) {
      stop("limit must be a numeric value bigger than 0")
    }
  }

  # Try to read the fusion report
  report <- withCallingHandlers({
      col_types_starfusion <- readr::cols_only(
##        "#FusionName" = col_skip(),
        "#FusionName" = col_character(), 
        "JunctionReadCount" = col_integer(),
        "SpanningFragCount" = col_integer(),
        "SpliceType" = col_skip(),
        "LeftGene" = col_character(),
        "LeftBreakpoint" = col_character(),
        "RightGene" = col_character(),
        "RightBreakpoint" = col_character(),
        "LargeAnchorSupport" = col_character(),
        "LeftBreakDinuc" = col_character(),
        "LeftBreakEntropy" = col_number(),
        "RightBreakDinuc" = col_character(),
        "RightBreakEntropy" = col_number(),
        "FFPM" = col_number(),
        "PROT_FUSION_TYPE" = col_character(),
        "annots" = col_character(), # e.g. '[SELFIE]'
        "FUSION_CDS" = col_character() # fusion sequence?
      )
      if (missing(limit)) {
        # Read all lines
        readr::read_tsv(
          file = filename,
          col_types = col_types_starfusion
        )
      } else {
        # Only read up to the limit
        readr::read_tsv(
          file = filename,
          col_types = col_types_starfusion,
          n_max = limit
        )
      }
    },
    error = function(cond) {
      message(paste0("Reading ", filename, " caused an error: ", cond[[1]]))
      stop(cond)
    },
    warning = function(cond) {
      # Begin Exclude Linting
      #message(paste0("Reading ", filename, " caused a warning: ", cond[[1]]))
      # End Exclude Linting
    }
  )

  if(use_fusioninspector)
      fi.table <- import_fusioninspector(limit=limit+10)

  # Set variables
  id                   <- NA
  inframe              <- NA
  fusion_tool          <- ifelse(use_fusioninspector, "starfusion+fusioninspector", "starfusion")
  spanning_reads_count <- NA
  split_reads_count    <- NA

  # List to hold all Fusion objects
  fusion_list <- vector("list", dim(report)[1])

  # Iterate through each line in the .tsv file
  for (i in 1:dim(report)[1]) {

    # Import starfusion-specific fields
    fusion_tool_specific_data <- list()
    fusion_tool_specific_data[["FusionName"]] = report[[i, "#FusionName"]]
    fusion_tool_specific_data[["LargeAnchorSupport"]] <-
      report[[i, "LargeAnchorSupport"]]
    fusion_tool_specific_data[["LeftBreakDinuc"]] <-
      report[[i, "LeftBreakDinuc"]]
    fusion_tool_specific_data[["LeftBreakEntropy"]] <-
      report[[i, "LeftBreakEntropy"]]
    fusion_tool_specific_data[["RightBreakDinuc"]] <-
      report[[i, "RightBreakDinuc"]]
    fusion_tool_specific_data[["RightBreakEntropy"]] <-
      report[[i, "RightBreakEntropy"]]
    fusion_tool_specific_data[["FFPM"]] <-
      report[[i, "FFPM"]]
    fusion_tool_specific_data[["annots"]] = report[[i, "annots"]]
    seq <- report[[i, "FUSION_CDS"]] # looks like '...aacgatcgtgaACTGATCGATG...'

    if(!is.null(report$PROT_FUSION_TYPE)) {
        ## only available when loading from coding_effect file (i.e.
        ## star-fusion.fusion_predictions.abridged.annotated.coding_effect.tsv)
        ft <- report[[i, "PROT_FUSION_TYPE"]]
        fusion_tool_specific_data [["inframe"]] <- ft
        if (ft == 'INFRAME')
          inframe <- TRUE
        if (ft == 'FRAMESHIFT')
          inframe <- FALSE
    }
    
    # id for this fusion
    id <- as.character(i)

    left_breakpoint <-
      unlist(strsplit(report[[i, "LeftBreakpoint"]], split = ":"))
    righ_breakpoint <-
      unlist(strsplit(report[[i, "RightBreakpoint"]], split = ":"))

    # Chromosome names
    chromosome_upstream <- left_breakpoint[1]
    chromosome_downstream <- righ_breakpoint[1]

    # Breakpoints
    breakpoint_upstream <- as.numeric(left_breakpoint[2])
    breakpoint_downstream <- as.numeric(righ_breakpoint[2])

    # Strand
    strand_upstream <- left_breakpoint[3]
    strand_downstream <- righ_breakpoint[3]

    # Number of supporting reads
    split_reads_count <- report[[i, "JunctionReadCount"]]
    spanning_reads_count <- report[[i, "SpanningFragCount"]]

    # STAR-Fusion doesn't provide the fusion sequence
    junction_sequence_upstream <- Biostrings::DNAString()
    junction_sequence_downstream <- Biostrings::DNAString()

    # Gene names
    gene_names_1 <- unlist(strsplit(report[[i, "LeftGene"]], split = "\\^"))
    gene_names_2 <- unlist(strsplit(report[[i, "RightGene"]], split = "\\^"))
    name_upstream <- gene_names_1[1]
    name_downstream <- gene_names_2[1]

    # Ensembl ids
    ensembl_id_upstream <- gene_names_1[2]
    ensembl_id_downstream <- gene_names_2[2]

    if (use_fusioninspector) {
        ## override things to go with the FI's virtual 'mini genome', but
        ## first save them:
        fusion_tool_specific_data[['orig_chromosome_upstream']] <- chromosome_upstream
        fusion_tool_specific_data[['orig_breakpoint_upstream']] <- breakpoint_upstream
        fusion_tool_specific_data[['orig_strand_upstream']] <- strand_upstream
        fusion_tool_specific_data[['orig_chromosome_downstream']] <- chromosome_downstream
        fusion_tool_specific_data[['orig_breakpoint_downstream']] <- breakpoint_downstream
        fusion_tool_specific_data[['orig_strand_downstream']] <- strand_downstream

        fusion.name <- report[[i, "#FusionName"]]
        if(fusion.name %in% fi.table$id) {
            chromosome_upstream <- fusion.name
            chromosome_downstream <- fusion.name
            strand_upstream <- '+'
            strand_downstream <- '+'
            breakpoint_upstream <- fi.table[fusion.name, 'start']
            breakpoint_downstream <- fi.table[fusion.name, 'end']
        } else {
            warning(sprintf("Could not find location for fusion %s
on virtual contig, ignoring it (line %d)", fusion.name, i+1))
            next
        }
    }
    
    # PartnerGene objects
    gene_upstream <- new(
      Class = "PartnerGene",
      name = name_upstream,
      ensembl_id = ensembl_id_upstream,
      chromosome = chromosome_upstream,
      breakpoint = breakpoint_upstream,
      strand = strand_upstream,
      junction_sequence = junction_sequence_upstream,
      transcripts = GenomicRanges::GRangesList()
    )
    gene_downstream <- new(
      Class = "PartnerGene",
      name = name_downstream,
      ensembl_id = ensembl_id_downstream,
      chromosome = chromosome_downstream,
      breakpoint = breakpoint_downstream,
      strand = strand_downstream,
      junction_sequence = junction_sequence_downstream,
      transcripts = GenomicRanges::GRangesList()
    )

    fusion_list[[i]] <- new(
      Class = "Fusion",
      id = id,
      fusion_tool = fusion_tool,
      genome_version = genome_version,
      spanning_reads_count = spanning_reads_count,
      split_reads_count = split_reads_count,
      fusion_reads_alignment = Gviz::AlignmentsTrack(),
      gene_upstream = gene_upstream,
      gene_downstream = gene_downstream,
      inframe = inframe,
      fusion_tool_specific_data = fusion_tool_specific_data
    )
  }

  # Return the list of Fusion objects (without the NULL elements)
  Filter(function(elt)!is.null(elt), fusion_list)
}
