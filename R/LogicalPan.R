#' A function for parsing genomic General Format Files
#' @param HomologList a list of matrices that represents sets of homologs
#' @param GeneCalls a list of matrices that represent gene calls
#' @keywords Homologs
#' @export
#' @examples
#' LogicalPan()

LogicalPan <- function(HomologList, 
                       GeneCalls,
                       Verbose = FALSE,
                       Plot = FALSE) {
  InitialMatrix <- vector("list",
                          length = length(HomologList))
  CatchMatrix <- vector("list",
                        length = length(HomologList))
  if (Verbose == TRUE) {
    pBar <- txtProgressBar(style = 1L)
    TimeStart <- Sys.time()
  }
  for (i in seq_along(InitialMatrix)) {
    InitialMatrix[[i]] <- vector("logical",
                                 length = ncol(HomologList[[i]]))
    CatchMatrix[[i]] <- vector("integer",
                               length = ncol(HomologList[[i]]))
    for (j in seq_along(InitialMatrix[[i]])) {
      if (any(!is.na(HomologList[[i]][, j]))) {
        InitialMatrix[[i]][j] <- TRUE
        CatchMatrix[[i]][j] <- unique(HomologList[[i]][which(!is.na(HomologList[[i]][, j])), j])
      }
    }
    if (Verbose == TRUE) {
      setTxtProgressBar(pb = pBar,
                        value = i/length(InitialMatrix))
    }
  }
  InitialMatrix <- do.call(cbind,
                           InitialMatrix)
  CatchMatrix <- do.call(cbind,
                         CatchMatrix)
  GeneRows <- sapply(GeneCalls,
                     function(x) nrow(x))
  MatRows <- vector("integer",
                    length = length(GeneCalls))
  for (i in seq_along(GeneCalls)) {
    MatRows[i] <- length(unique(CatchMatrix[i, which(!is.na(CatchMatrix[i, ]))]))
    if (is.na(MatRows[i])) {
      MatRows[i] <- 0L
    }
  }
  AddRows <- GeneRows - MatRows
  AdditionalMatrix <- matrix(FALSE,
                             ncol = sum(AddRows),
                             nrow = length(GeneCalls))
  Count <- 1L
  for (i in seq_along(AddRows)) {
    CurrentRows <- AddRows[i]
    AdditionalMatrix[i, (Count:(Count + CurrentRows - 1L))] <- TRUE
    Count <- Count + CurrentRows
  }
  PanMatrix <- cbind(InitialMatrix,
                     AdditionalMatrix)
  PanMatrix <- PanMatrix[, order(colSums(PanMatrix),
                                 decreasing = TRUE)]
  if (Verbose == TRUE) {
    cat("\n")
    TimeStop <- Sys.time()
    TotalTime <- (TimeStop - TimeStart)
    print(TotalTime)
  }
  if (Plot == TRUE) {
    image(t(PanMatrix),
          col = c("white", "blue"),
          main = "Presence / Absence")
  }
  return(PanMatrix)
}









