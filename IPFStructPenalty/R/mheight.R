#' IPFStructPenalty
#' @title Structured penalized regression
#' @description
#' Subfunction for the tree-lasso and IPF-tree-lass methods to compute the minimun height of the dendrogram.
#' @param y response matrix
#' @param foldid an optional vector identifying what fold each observation is in. The \code{mheigth} returns the maximum of the leave-one-fold data's minimun height. Default is \code{NULL}.
#' @export
mheight <- function(y, foldid=NULL){
  if(is.null(foldid)){
    myDist0 <- 1 - abs(fastCorr(y))
    myDist <- myDist0[lower.tri(myDist0)]
    a0 <- dist(t(y))
    a0[1:length(a0)] <- myDist
    # hierarchical clustering for multivariate responses
    myCluster <- hclust(a0, method = "complete")
    min.height <- min(myCluster$height/max(myCluster$height))
  }else{
    min.height <- 0
    for(i in 1:max(foldid)){
      y0 <- y[!foldid==i,]
      myDist0 <- 1 - abs(fastCorr(y0))
      myDist <- myDist0[lower.tri(myDist0)]
      a0 <- dist(t(y0))
      a0[1:length(a0)] <- myDist
      # hierarchical clustering for multivariate responses
      myCluster <- hclust(a0, method = "complete")
      min.height <- max(min(myCluster$height/max(myCluster$height)), min.height)
    }
  }
  return(min.height)
}
