# DataClean.R
# This function cleans data sets.  For now, this means removing bad columns.

#' @export
DataClean <- function(dat, clippy = FALSE, rm.uniform.cols = TRUE) {

  num.unq.val <- rep(0, ncol(dat))

  for(i in 1:ncol(dat)){
    num.unq.val [i] <- length(unique(dat[, i]))
  }

  bad.cols <- which(num.unq.val == 1 | num.unq.val == nrow(dat))

  if(length(bad.cols) > 0) {
    if (clippy == TRUE) {
      cat("Data has completely uniform and/or uniquely-valued columns.")
      cat("Would you like to remove these columns?\ny/n: ")
      x <- readLines(con = stdin(), 1)
      
      while ( x != "y" & x != "n") {
        cat("Type only 'y' or 'n': ")
        x <- readLines(con = stdin(), 1)
      }
      
      if (x == "y") {
        dat <- dat[, -bad.cols]
      }
    } else {
      if (rm.uniform.cols == TRUE) {
        dat <- dat[, -bad.cols]
      }
    }
  }

  dat
}
