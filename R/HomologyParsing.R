#' A function for determining nucleotide overlap that links genes.
#' @param ResultMatrix A matrix where with the upper triangle containing matrices of genes linked by nucleotide overlap
#' @keywords Homology, Orthology, Core Genome
#' @export
#' @examples
#' Catalog()

Catalog <- function(ResultMatrix,
                    Verbose = FALSE) {

  if (Verbose == TRUE) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  Cycler <- function(StartPosition,
                     TotalPositions) {
    if (TotalPositions != 1L &
        StartPosition != 1L) {
      return(c(StartPosition:TotalPositions, 1L:(StartPosition - 1L)))
    } else if (TotalPositions == 1L) {
      return(1L)
    } else if (StartPosition == 1L) {
      return(StartPosition:TotalPositions)
    }
  }
  ######
  # Collect size of matrix
  # Initialize structures
  ######
  Size <- dim(ResultMatrix)[1]
  Pairing <- matrix(data = vector("list",
                                  length = 1L),
                    nrow = dim(ResultMatrix)[1L],
                    ncol = dim(ResultMatrix)[2L])
  ######
  # Pull pairings from initial data structure
  ######
  for (m1 in 1:(Size - 1L)) {
    for (m2 in (m1 + 1L):Size) {
      Pairing[m1, m2][[1]] <- ResultMatrix[m1, m2][[1]]
    }
  }
  ######
  # Get dimensions
  # Initialize intermediate output structure
  ######
  FilledPositions <- Pairing[upper.tri(Pairing)]
  DimMatrix <- sapply(FilledPositions,
                      function(x) dim(x))
  PairsMatrix <- matrix(data = NA_integer_,
                        nrow = sum(DimMatrix[1, ]),
                        ncol = Size)
  ######
  # Fill in all pairs into intermediate output
  ######
  Count <- 1L
  for (m1 in seq_len(Size - 1L)) {
    for (m2 in (m1 + 1L):Size) {
      CurrentRows <- dim(Pairing[m1, m2][[1]])[1]
      PairsMatrix[(Count:(Count + CurrentRows - 1L)), m1] <- as.integer(Pairing[m1, m2][[1]][, 1])
      PairsMatrix[(Count:(Count + CurrentRows - 1L)), m2] <- as.integer(Pairing[m1, m2][[1]][, 2])
      Count <- Count + CurrentRows
    }
  }
  ######
  # Return all sets of genes that are correlated
  # These will be both transitive and intransitive
  ######
  GeneSets <- vector("list",
                     length = nrow(PairsMatrix))
  CoreSize <- ((ncol(PairsMatrix) - 1L) * ncol(PairsMatrix)) / 2L
  for (k1 in seq_along(GeneSets)) {
    MatchVector <- PairsMatrix[k1, ]
    Cycle <- Cycler(StartPosition = k1,
                    TotalPositions = length(GeneSets))
    GeneSetIDs <- vector("integer",
                         length = CoreSize)
    AddCount <- 1L
    GeneSetIDs[AddCount] <- k1
    for (k2 in seq_along(Cycle)) {
      Objective <- PairsMatrix[Cycle[k2], ]
      if (all(MatchVector == Objective, na.rm = TRUE) &
          !is.na(all(MatchVector == Objective, na.rm = TRUE)) &
          !all(is.na(MatchVector == Objective))) {
        MatchVector[which(is.na(MatchVector) &
                            !is.na(Objective))] <- Objective[which(is.na(MatchVector) &
                                                                     !is.na(Objective))]
        AddCount <- AddCount + 1L
        GeneSetIDs[AddCount] <- Cycle[k2]
      }
    }
    GeneSetIDs <- GeneSetIDs[which(GeneSetIDs > 0L)]
    GeneSetIDs <- sort(unique(GeneSetIDs))
    GeneSets[[k1]] <- PairsMatrix[GeneSetIDs,
                                  ,
                                  drop = FALSE]
    if (Verbose == TRUE) {
      setTxtProgressBar(pb = pBar,
                        value = k1/length(GeneSets))
    }
  }
  GeneSets <- unique(GeneSets)
  if (Verbose == TRUE) {
    TimeStop <- Sys.time()
    cat("\n")
    print(TimeStop - TimeStart)
  }
  return(GeneSets)
}
