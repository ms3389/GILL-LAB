is_X_between_Y_n_Z <- function(x, y, z) {
  #dependencies
  require(tidyverse)
  require(magrittr)
  require(IRanges)
  #create single-digit IRanges object
  x = IRanges(x)
  #create IRanges object from start and end numbers
  r = IRanges(start = y, end = z)
  #return only y and z pairs that correspond to positive matches with x
  Ans = IRanges::as.data.frame(IRanges::findOverlapPairs(x, r, type = "within"))
  return(Ans)
}