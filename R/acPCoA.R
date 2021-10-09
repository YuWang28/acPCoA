#' Perform AC-PCoA for confounding factor adjustement based on Principle Coordinate Analysis
#'
#' @param DistanceMatrix the n by n data distance matrix, where n is the number of samples. The (i,j)-th entry stands for the distance (dissimilarity) between the i-th sample and the j-th sample.
#' @param ConfounderMatrix the n by q confounder matrix, where n is the number of samples, q is the number of confounding factors.
#' @param nPC number of principal components to compute
#' @param lambdas the tuning parameter, non-negative.
#' @param kernel the kernel to use: "linear", "gaussian".
#' @param bandwidth bandwidth h for Gaussian kernel. Optional.
#' @param anov True or False. Whether the penalty term has the between groups sum of squares interpretation. Default is True.
#' @param perc the best lambda is defined to be the smallest lambda with R(lambda)<=perc (if anov=T), or R(lambda)<=perc*R(lambda=0) (if anov=F) in the nPC principal components.
#'
#' @return The principal components and the projected data
#' \item{v}{the principal components, p by nPC matrix}
#' \item{Xv}{the projected data, i.e. X times v}
#' \item{eigenvalueX}{eigenvalues for the PCs}
#' \item{varianceX}{variance explained by the PCs}
#' \item{varianceX_perc}{percentage of total variance in X explained by the PCs. If eval=F, NA is returned.}
#' \item{eigenvalueXperm}{eigenvalues for the PCs, permutation. If eval=F, NA is returned.}
#' \item{varianceXperm}{variance explained by the PCs, permutation. If eval=F, NA is returned.}
#' \item{sigPC}{the significant PCs. If eval=F, NA is returned.}
#' \item{...}{Input parameters for the function}
#' @export
#'
#' @import acPCA
#' @examples
#' \dontrun{
#' X <- data_mbqc_groupA$DistMat.BC;
#' Y <- data_mbqc_groupA$ConfounderMat;
#' result_acPCoA <- acPCoA(DistanceMatrix=X, ConfounderMatrix=Y, nPC=2, lambda=seq(0, 20, 0.05), kernel="linear")
#' ggplot(as.data.frame(result_acPCoA$Xv),aes(x=V1,y=V2,color=data_mbqc_groupA$Specimen))+geom_point()
#' }
acPCoA <- function (DistanceMatrix, ConfounderMatrix, nPC=2, lambdas=seq(0, 20, 0.05),kernel="linear",bandwidth=NULL,anov=T,perc=0.05)
{
  D=as.matrix(DistanceMatrix)
  centre <- function(D, n) {
    One <- matrix(1, n, n)
    mat <- diag(n) - One/n
    mat.cen <- mat %*% D %*% mat
  }
  D <- as.matrix(D)
  n <- nrow(D)
  epsilon <- sqrt(.Machine$double.eps)*10000
  names <- rownames(D)

  delta1 <- centre((-0.5 * D^2), n)
  trace <- sum(diag(delta1))
  D.eig <- eigen(delta1,symmetric=TRUE)
  #D.eig$values=Re(D.eig$values)
  min.eig <- min(D.eig$values)
  zero.eig <- which(abs(D.eig$values) < epsilon)
  D.eig$values[zero.eig] <- 0
    eig <- D.eig$values
    k <- length(which(eig > epsilon))
    vectors <- sweep(D.eig$vectors[, 1:k], 2, sqrt(eig[1:k]),
                     FUN = "*")
    values <- data.frame(eig[1:k])
    colnames(values) <- c("Eigenvalues")
    rownames(values) <- 1:nrow(values)
    rownames(vectors) <- names
    colnames(vectors) <- colnames(vectors, do.NULL = FALSE,
                                  prefix = "Axis.")
    result_tune <- acPCAtuneLambda(X=vectors, Y=ConfounderMatrix, nPC, lambdas, centerX=T, centerY=T, scaleX=F, scaleY=F, kernel, bandwidth, anov, perc, quiet=T)
    result <- acPCA(X=vectors, Y=ConfounderMatrix, nPC, lambda=result_tune$best_lambda, eval=F, numPerm=100, alpha=0.05, plot=T, centerX=T, centerY=T, scaleX=F, scaleY=F, kernel, bandwidth)
    return(result)

}
