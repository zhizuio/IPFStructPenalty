#' IPFStructPenalty
#' @title Function producing stratified/ balanced folds for cross validation
#' @description
#' Get balanced folds for cross validation, which are used for tuning penalization parameters. This function is mainly used within the function \code{epsgo}. See the \code{R} package \pkg{c060} for details.
#' @export
balancedFolds <- function(class.column.factor, cross.outer)
{
  #stolen from MCREstimate package
  # used for stratified(balanced) classification or regression
  # get balanced folds from pamr
  sampleOfFolds  <- get("balanced.folds",envir=asNamespace("pamr"))(class.column.factor, nfolds=cross.outer)
  permutated.cut <- rep(0,length(class.column.factor))
  for (sample in 1:cross.outer)
  {
    cat(sample,"\n")
    permutated.cut[sampleOfFolds[[sample]]] <- sample
  }
  return(permutated.cut)
}
