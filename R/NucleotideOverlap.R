#' A function for determining nucleotide overlap that links genes.
#' @param SyntenyObject An object of class "Synteny"
#' @param PATH A character string identifying a sqlite database
#' @param GeneCalls A list of dataframes, one for each genome, with named columns "Start", "Stop", and "Strand"
#' @keywords Homology
#' @export
#' @examples
#' NucleotideOverlap()

NucleotideOverlap <- function(SyntenyObject,
                              GeneCalls,
                              PATH,
                              Verbose = FALSE) {

  L <- nrow(SyntenyObject)
  stopifnot("DECIPHER" %in% .packages(),
            is(SyntenyObject, "Synteny"),
            L > 1L,
            L == length(GeneCalls))

  if (Verbose == TRUE) {
    TotalTimeStart <- Sys.time()
  }
  ResultMatrix <- matrix(data = vector("list",
                                       length = 1L),
                         nrow = L,
                         ncol = L)
  TotalLength <- L^2 - L
  TotalCounter <- 0L
  ######
  # Function to extend matrices
  ######
  Ext.Check <- function(CurrentMatrix,
                        PositionCounter,
                        AdditionalRows) {
    if (PositionCounter > nrow(CurrentMatrix) * 0.95) {
      CurrentMatrix <- rbind(CurrentMatrix,
                             AdditionalRows)
    }
    return(CurrentMatrix)
  }
  ######
  # scroll through every hit table in the synteny object
  ######
  pBar <- txtProgressBar(style = 1L)
  for (m1 in seq_len(L - 1L)) {
    for (m2 in (m1 +1L):L) {
      ######
      # Collect Start Stop and Strand from current subject and query
      ######
      Q.Start <- as.integer(GeneCalls[[m1]]$Start)
      Q.Stop <- as.integer(GeneCalls[[m1]]$Stop)
      S.Start <- as.integer(GeneCalls[[m2]]$Start)
      S.Stop <- as.integer(GeneCalls[[m2]]$Stop)
      Q.Length <- Q.Stop - Q.Start + 1L
      S.Length <- S.Stop - S.Start + 1L
      ######
      # Select current hit table, subset to only single index
      ######
      CurrentHitTable <- SyntenyObject[m1, m2, drop = FALSE][[1]]
      CurrentHitTable <- CurrentHitTable[which(CurrentHitTable[, "index1"] == 1L), ]
      CurrentHitTable <- CurrentHitTable[which(CurrentHitTable[, "index2"] == 1L), ]
      CurrentHitTable <- CurrentHitTable[order(CurrentHitTable[,
                                                               "start1",
                                                               drop = FALSE]),
                                         ,
                                         drop = FALSE]
      ######
      # Collect the starts and stops for all hits correct for strandedness
      ######
      HitWidths <- CurrentHitTable[, "width"]
      Q.HitStarts <- CurrentHitTable[, "start1"]
      S.HitStarts <- CurrentHitTable[, "start2"]
      Q.HitEnds <- Q.HitStarts + HitWidths - 1L
      S.HitEnds <- S.HitStarts + HitWidths - 1L
      Strand <- CurrentHitTable[, "strand"]
      for (i in seq_along(HitWidths)) {
        if (Strand[i] == 0L) {
          S.HitStarts[i] <- CurrentHitTable[i, "start2"]
          S.HitEnds[i] <- CurrentHitTable[i, "start2"] + CurrentHitTable[i, "width"] - 1L
        } else if (Strand[i] == 1L) {
          S.HitStarts[i] <- CurrentHitTable[i, "start2"] - CurrentHitTable[i, "width"] + 1L
          S.HitEnds[i] <- CurrentHitTable[i, "start2"]
        }
      }
      ######
      # Record into Q.Match every hit index that lands in a gene in the query
      ######
      QueryMatrix <- matrix(NA_integer_,
                            ncol = 5L,
                            nrow = nrow(CurrentHitTable))
      ExtraRows <- QueryMatrix
      HitCounter <- 1L
      AddCounter <- 1L
      # QueryMatrix <- matrix(rep(NA_integer_, times = 5L), ncol = 5L)
      for (z1 in seq_along(Q.Start)) {
        CurrentGene <- HitCounter
        Q.NucOverLapL <- NA_integer_
        Q.NucOverLapR <- NA_integer_
        S.NucPositionL <- NA_integer_
        S.NucPositionR <- NA_integer_
        while (HitCounter <= length(HitWidths)) {
          if (Q.HitEnds[HitCounter] < Q.Start[z1]) {
            # Hit ends before current query gene begins
            HitCounter <- HitCounter + 1L
          } else if (Q.HitStarts[HitCounter] < Q.Start[z1] &
                     Q.HitEnds[HitCounter] >= Q.Start[z1] &
                     Q.HitEnds[HitCounter] <= Q.Stop[z1]) {
            # Hit overlaps left bound of current query gene
            CurrentGene <- z1
            Q.NucOverLapL <- Q.Start[z1]
            Q.NucOverLapR <- Q.HitEnds[HitCounter]
            NewWidth <- Q.NucOverLapR - Q.NucOverLapL + 1L
            S.NucPositionR <- S.HitEnds[HitCounter]
            S.NucPositionL <- S.HitEnds[HitCounter] - NewWidth + 1L
            # Add To Vector!
            QueryMatrix[AddCounter, ] <- c(CurrentGene,
                                           Q.NucOverLapL,
                                           Q.NucOverLapR,
                                           S.NucPositionL,
                                           S.NucPositionR)
            QueryMatrix <- Ext.Check(CurrentMatrix = QueryMatrix,
                                     PositionCounter = AddCounter,
                                     AdditionalRows = ExtraRows)
            AddCounter <- AddCounter + 1L
            HitCounter <- HitCounter + 1L
          } else if (Q.HitStarts[HitCounter] >= Q.Start[z1] &
                     Q.HitEnds[HitCounter] <= Q.Stop[z1]) {
            # Hit occurs entirely within current query gene
            CurrentGene <- z1
            Q.NucOverLapL <- Q.HitStarts[HitCounter]
            Q.NucOverLapR <- Q.HitEnds[HitCounter]
            S.NucPositionL <- S.HitStarts[HitCounter]
            S.NucPositionR <- S.HitEnds[HitCounter]
            # Add To Vector!
            QueryMatrix[AddCounter, ] <- c(CurrentGene,
                                           Q.NucOverLapL,
                                           Q.NucOverLapR,
                                           S.NucPositionL,
                                           S.NucPositionR)
            QueryMatrix <- Ext.Check(CurrentMatrix = QueryMatrix,
                                     PositionCounter = AddCounter,
                                     AdditionalRows = ExtraRows)
            AddCounter <- AddCounter + 1L
            HitCounter <- HitCounter + 1L
          } else if (Q.HitStarts[HitCounter] >= Q.Start[z1] &
                     Q.HitStarts[HitCounter] <= Q.Stop[z1] &
                     Q.HitEnds[HitCounter] > Q.Stop[z1]) {
            # Hit overlaps right bound of current query gene
            CurrentGene <- z1
            Q.NucOverLapL <- Q.HitStarts[HitCounter]
            Q.NucOverLapR <- Q.Stop[z1]
            NewWidth <- Q.NucOverLapR - Q.NucOverLapL + 1L
            S.NucPositionL <- S.HitStarts[HitCounter]
            S.NucPositionR <- S.HitStarts[HitCounter] + NewWidth - 1L
            # Add To Vector!
            QueryMatrix[AddCounter, ] <- c(CurrentGene,
                                           Q.NucOverLapL,
                                           Q.NucOverLapR,
                                           S.NucPositionL,
                                           S.NucPositionR)
            QueryMatrix <- Ext.Check(CurrentMatrix = QueryMatrix,
                                     PositionCounter = AddCounter,
                                     AdditionalRows = ExtraRows)
            AddCounter <- AddCounter + 1L
            HitCounter <- HitCounter + 1L
            break
          } else if (Q.HitStarts[HitCounter] <= Q.Start[z1] &
                     Q.HitEnds[HitCounter] >= Q.Stop[z1]) {
            # Hit eclipses current query gene
            CurrentGene <- z1
            Q.NucOverLapL <- Q.Start[z1]
            Q.NucOverLapR <- Q.Stop[z1]
            S.NucPositionL <- S.HitStarts[HitCounter] + (Q.HitStarts[HitCounter] - Q.Start[z1]) + 1L
            S.NucPositionR <- S.HitEnds[HitCounter] - (Q.HitEnds[HitCounter] - Q.Stop[z1]) - 1L
            # Add To Vector!
            QueryMatrix[AddCounter, ] <- c(CurrentGene,
                                           Q.NucOverLapL,
                                           Q.NucOverLapR,
                                           S.NucPositionL,
                                           S.NucPositionR)
            QueryMatrix <- Ext.Check(CurrentMatrix = QueryMatrix,
                                     PositionCounter = AddCounter,
                                     AdditionalRows = ExtraRows)
            AddCounter <- AddCounter + 1L
            break
          } else if (Q.HitStarts[HitCounter] > Q.Stop[z1]) {
            # Hit occurs after current query gene
            break
          }
        }
      }
      if (Verbose == TRUE) {
        TotalCounter <- TotalCounter + 1L
        setTxtProgressBar(pb = pBar,
                          value = TotalCounter/TotalLength)
      }
      QueryMatrix <- QueryMatrix[apply(QueryMatrix,
                                       1L,
                                       function(x) !all(is.na(x))),
                                 ,
                                 drop = FALSE]
      if (dim(QueryMatrix)[1] == 0L) {
        OutPutMatrix <- matrix(NA_integer_,
                               nrow = 1L,
                               ncol = 3L)
        OverLapMatrix <- matrix(NA_integer_,
                                nrow = 1L,
                                ncol = 3L)
      } else {
        HitCounter <- 1L
        AddCounter <- 1L
        OverLapMatrix <- matrix(NA_integer_,
                                nrow = nrow(QueryMatrix),
                                ncol = 3L)
        QueryMatrix <- QueryMatrix[order(QueryMatrix[,
                                                     4L,
                                                     drop = FALSE]),
                                   ,
                                   drop = FALSE]
        ExtraRows <- OverLapMatrix
        QueryMap <- QueryMatrix[, 1L]
        S.HitStarts <- QueryMatrix[, 4L]
        S.HitEnds <- QueryMatrix[, 5L]
        for (z2 in seq_along(S.Start)) {
          while (HitCounter <= length(QueryMap)) {
            if (S.HitEnds[HitCounter] < S.Start[z2]) {
              # Hit ends before current subject gene begins
              HitCounter <- HitCounter + 1L
            } else if (S.HitStarts[HitCounter] < S.Start[z2] &
                       S.HitEnds[HitCounter] >= S.Start[z2] &
                       S.HitEnds[HitCounter] <= S.Stop[z2]) {
              # Hit overlaps left bound of current subject gene
              CurrentGene <- z2
              QueryGenePosition <- QueryMap[HitCounter]
              ExactOverLap <- S.HitEnds[HitCounter] - S.Start[z2] + 1L
              # Add to vector !
              OverLapMatrix[AddCounter, ] <- c(QueryGenePosition,
                                               CurrentGene,
                                               ExactOverLap)
              OverLapMatrix <- Ext.Check(CurrentMatrix = OverLapMatrix,
                                         PositionCounter = AddCounter,
                                         AdditionalRows = ExtraRows)
              AddCounter <- AddCounter + 1L
              HitCounter <- HitCounter + 1L
            } else if (S.HitStarts[HitCounter] >= S.Start[z2] &
                       S.HitEnds[HitCounter] <= S.Stop[z2]) {
              # Hit occurs entirely within current subject gene
              CurrentGene <- z2
              QueryGenePosition <- QueryMap[HitCounter]
              ExactOverLap <- S.HitEnds[HitCounter] - S.HitStarts[HitCounter] + 1L
              # Add to vector !
              OverLapMatrix[AddCounter, ] <- c(QueryGenePosition,
                                               CurrentGene,
                                               ExactOverLap)
              OverLapMatrix <- Ext.Check(CurrentMatrix = OverLapMatrix,
                                         PositionCounter = AddCounter,
                                         AdditionalRows = ExtraRows)
              AddCounter <- AddCounter + 1L
              HitCounter <- HitCounter + 1L
            } else if (S.HitStarts[HitCounter] >= S.Start[z2] &
                       S.HitStarts[HitCounter] <= S.Stop[z2] &
                       S.HitEnds[HitCounter] > S.Stop[z2]) {
              # Hit overlaps right bound of current subject gene
              CurrentGene <- z2
              QueryGenePosition <- QueryMap[HitCounter]
              ExactOverLap <- S.Stop[z2] - S.HitStarts[HitCounter] + 1L
              # Add to vector !
              OverLapMatrix[AddCounter, ] <- c(QueryGenePosition,
                                               CurrentGene,
                                               ExactOverLap)
              OverLapMatrix <- Ext.Check(CurrentMatrix = OverLapMatrix,
                                         PositionCounter = AddCounter,
                                         AdditionalRows = ExtraRows)
              AddCounter <- AddCounter + 1L
              break
            } else if (S.HitStarts[HitCounter] <= S.Start[z2] &
                       S.HitEnds[HitCounter] >= S.Stop[z2]) {
              # Hit eclipses current subject gene
              CurrentGene <- z2
              QueryGenePosition <- QueryMap[HitCounter]
              ExactOverLap <- S.Stop[z2] - S.Start[z2] + 1L
              # Add to vector !
              OverLapMatrix[AddCounter, ] <- c(QueryGenePosition,
                                               CurrentGene,
                                               ExactOverLap)
              OverLapMatrix <- Ext.Check(CurrentMatrix = OverLapMatrix,
                                         PositionCounter = AddCounter,
                                         AdditionalRows = ExtraRows)
              AddCounter <- AddCounter + 1L
              break
            } else if (S.HitStarts[HitCounter] > S.Stop[z2]) {
              # Hit occurs after current subject gene
              break
            }
          }
        }
      }
      OverLapMatrix <- OverLapMatrix[apply(OverLapMatrix,
                                           1L,
                                           function(x) !all(is.na(x))),
                                     ,
                                     drop = FALSE]
      if (dim(OverLapMatrix)[1] == 0L) {
        OverLapMatrix <- matrix(NA_integer_,
                                nrow = 1L,
                                ncol = 3L)
        OutPutMatrix <- matrix(NA_integer_,
                               nrow = 1L,
                               ncol = 3L)
      } else if (dim(OverLapMatrix)[1] > 1L) {
        OutPutMatrix <- OverLapMatrix
        for (z4 in seq_len(nrow(OutPutMatrix))) {
          z5 <- which(OverLapMatrix[, 1L] == OverLapMatrix[z4, 1L] &
                        OverLapMatrix[, 2L] == OverLapMatrix[z4, 2L])
          OutPutMatrix[z4, 3L] <- sum(OverLapMatrix[z5, 3L])
        }
      }
      OutPutMatrix <- unique(OutPutMatrix)
      OutPutMatrix <- OutPutMatrix[order(OutPutMatrix[,
                                                      1L,
                                                      drop = FALSE]),
                                   ,
                                   drop = FALSE]
      if (Verbose == TRUE) {
        TotalCounter <- TotalCounter + 1L
        setTxtProgressBar(pb = pBar,
                          value = TotalCounter/TotalLength)
      }
      colnames(OutPutMatrix) <- c("QueryIndex", "SubjectIndex", "NucleotideOverlap")
      ResultMatrix[m1, m2] <- list(OutPutMatrix)
    }
  }
  if (Verbose == TRUE) {
    TotalTimeStop <- Sys.time()

    cat("\n")
    print(TotalTimeStop - TotalTimeStart)
  }
  return(ResultMatrix)
}










