#' @title Normalize miRNA According to GMPR
#'
#' @description
#' 1st step: calculate pairwise r_jk for two samples j and k
#' only look at those genes in the samples for which
#' count_ij * count_ik != zero
#' for those were it is not zero, look at the ratio
#' count_ij / count_ik. Then take the median of those.
#' 
#' @details
#' 
#' @param counts A data.frame with counts
#' 
#' @return MISSING
#' @author Nikos Bosse \email{nikosbosse@gmail.com}
#' @examples
#'  
#' counts <- mirmeth::example_data_wide
#' s_gmpr(counts)
#'
#' @export


gmpr_normalization <- function(counts) {
  counts <- counts[, !(colnames(counts) == "id")]
  n = ncol(counts)
  s_gmpr <- vector(mode = 'numeric', length = n)
  for (j in 1:n) {
    r <- vector(mode = 'numeric', length = n)
    for (k in 1:n) {
      v <- counts[,j]
      w <- counts[,k]
      nonzero <- v != 0 & w
      r[k] <- median(v[nonzero] / w[nonzero])
    }
    s_gmpr[j] <- exp((1/n) * sum(log(r)))
  }
  return(s_gmpr)
}




# function not yet working


deseq_normalization <- function(countsmatrix) {
  
  ##############################################
  # alternative: DESeq2 normalization approach
  rowprodsfun <- function(x) {
    prod(x[x !=0]) ^ (1/n) # over sum(ip)?
  }
  Ki_R <- apply(countsmatrix[, ip == 0], MARGIN = 1, FUN = rowprodsfun)
  # Ki_R <- matrixStats::rowProds(counts[, ip == 0])
  s_deseq <- counts[] / Ki_R # counts[, ip == 0]?
  s_deseq <- matrixStats::colMedians(as.matrix(s_deseq)) 
  ##############################################
}


