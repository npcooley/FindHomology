#' A function for parsing genomic General Format Files
#' @param GFFAddress a matrix of lists, the upper triange of which is filled with matrices of putative homologs
#' @keywords GFF
#' @export
#' @examples
#' GFFParser()


GFFParser <- function(GFFAddress,
                      Verbose = FALSE) {
  GeneCalls <- vector("list",
                      length = length(GFFAddress))
  if (Verbose == TRUE) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  for (i in seq_along(GFFAddress)) {
    z1 <- gzcon(url(GFFAddress[i]))
    z2 <- textConnection(readLines(z1))
    z3 <- readLines(z2)
    z4 <- strsplit(z3,
                   split = "\t")
    Start <- sapply(z4,
                    function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                       yes = x[4],
                                       no = NA),
                    USE.NAMES = FALSE)
    Stop <- sapply(z4,
                   function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                      yes = x[5],
                                      no = NA),
                   USE.NAMES = FALSE)
    Strand <- sapply(z4,
                     function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                        yes = x[7],
                                        no = NA),
                     USE.NAMES = FALSE)
    Product <- sapply(z4,
                      function(x) ifelse(test = length(x) == 9L & x[3] == "CDS",
                                         yes = x[9],
                                         no = NA),
                      USE.NAMES = FALSE)
    Product <- str_extract(Product,
                           "(?<=product=)(.*)(?=;protein_id)")
    z5 <- as.integer(Start[which(!is.na(Start) &
                                   !is.na(Stop) &
                                   !is.na(Strand) &
                                   !is.na(Product))])
    z6 <- as.integer(Stop[which(!is.na(Start) &
                                  !is.na(Stop) &
                                  !is.na(Strand) &
                                  !is.na(Product))])
    z7 <- Strand[which(!is.na(Start) &
                         !is.na(Stop) &
                         !is.na(Strand) &
                         !is.na(Product))]
    z8 <- Product[which(!is.na(Start) &
                          !is.na(Stop) &
                          !is.na(Strand) &
                          !is.na(Product))]
    GeneCalls[[i]] <- data.frame("Start" = z5,
                                 "Stop" = z6,
                                 "Strand" = z7,
                                 "Annotation" = z8,
                                 stringsAsFactors = FALSE)
    if (Verbose == TRUE) {
      setTxtProgressBar(pb = pBar,
                        value = i/length(GFFAddress))
    }
  }
  if (Verbose == TRUE) {
    TimeStop <- Sys.time()
    TotalTime <- (TimeStop - TimeStart)
    cat("\n")
    print(TotalTime)
  }
  return(GeneCalls)
}

