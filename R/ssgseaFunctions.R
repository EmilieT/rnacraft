
#
# SSP-GSEA functions
#

#
#' SSP-GSEA Empirical Cumulative Distribution Function
#' @description
#' Compute the difference between a weighted Empirical Cumulative Distribution Function
#' (ECDF) of the genes in the signature and the ECDF of the remaining genes.
#' @param signature character vector of gene names
#' @param data numeric vector of expression data, named by genes
#' @param alpha weighting factor
#' @return vector (same length as data), containing the ECDF difference for each gene
#' @note this function is internally used by \link{ssgsea.score}. usually, you should not
#' call it directly.
#' @seealso \link{ssgsea.score}
#'
#
ssgsea.ecdf <- function(signature, data, alpha=0.25) {
  signature <- unique(signature)
  signok <- intersect(signature, names(data))
  if (length(signok) != length(signature)) {
    dif = setdiff(signature, signok)
    warning(length(dif), " genes not found : ", paste(dif, collapse=","))
  }
  rnk <- rank(-data)
  pin <- rnk[names(sort(-data))]^alpha
  pin[! (names(pin) %in% signok)] <- 0
  pin <- cumsum(pin) / sum(pin[signok])

  pout <- rep(1, length(data)) / (length(data) - length(signok))
  names(pout) <- names(pin)
  pout[names(pout) %in% signok] <- 0
  pout <- cumsum(pout)

  pin - pout
}

#
#' SSP-GSEA score of single sample on signature
#' @param signature character vector of gene names
#' @param data numeric vector of expression data, named by genes
#' @param alpha weighting factor
#' @return numeric score of sample on signature
#'
#
ssgsea.score <- function(signature, data, alpha=0.25) {
  sum(ssgsea.ecdf(signature, data, alpha))
}


#
#' SSP-GSEA calculation for all samples/ all signatures
#

ssGSEAcalculation <- function(sig.lst,rna.tab){
  score.part <- lapply(sig.lst, function(sig) {
    cat(length(sig), "genes\n")
    apply(rna.tab, 2, function(col) ssgsea.score(sig, col))
  })

  sig.nam <- names(sig.lst)
  sig.nam <- sig.nam[substr(sig.nam,nchar(sig.nam)-2,nchar(sig.nam)) == "_UP"]
  sig.nam <- substr(sig.nam,1,nchar(sig.nam)-3)

  score.comb <- lapply(sig.nam, function(nam) {
    score.part[[paste0(nam, "_UP")]] - score.part[[paste0(nam, "_DN")]]
  })
  names(score.comb) <- sig.nam
  score.es <- c(score.part, score.comb)
  return(do.call("rbind",score.es))
}


