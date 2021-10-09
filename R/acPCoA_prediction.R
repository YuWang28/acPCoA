#' Perform AC-PCoA on new data points for classification and prediction
#'
#' @param training.dist the n by n data distance matrix, where n is the number of training samples.
#' @param testing.dist the n by m data distance matrix, where n is the number of training samples and m is the number of testing samples.
#' @param training.confounder the n by q confounder matrix of training samples, where n is the number of samples, q is the number of confounding factors.
#' @param training.label Label of traning samples for classification. Optional.
#' @param testing.label Label of testing samples for classification. Optional.
#' @param nPC number of principal components to compute
#' @param kernel the kernel to use: "linear", "gaussian".
#' @param bandwidth bandwidth h for Gaussian kernel. Optional.
#' @param lambdas the tuning parameter, non-negative.
#' @param anov True or False. Whether the penalty term has the between groups sum of squares interpretation. Default is True.
#' @param perc the best lambda is defined to be the smallest lambda with R(lambda)<=perc (if anov=T), or R(lambda)<=perc*R(lambda=0) (if anov=F) in the nPC principal components.
#'
#' @return principal components and the projected data of training and testing samples
#' \item{v}{the principal components, p by nPC matrix}
#' \item{Xv_training}{the projected data of training samples, i.e. X times v}
#' \item{Xv_testing}{the projected data of testing samples}
#' \item{eigenvalueX}{eigenvalues for the PCs}
#' \item{varianceX}{variance explained by the PCs}
#' \item{PredictingAccuracy}{the accuracy of predicting labels of testing samples using Random Forest, if training and testing labels are provided.}
#' @export
#'
#' @examples
#' \dontrun{
#' training.dist=data_mbqc_groupA$DistMat.BC[1:679,1:679]
#' testing.dist=data_mbqc_groupA$DistMat.BC[1:679,680:848]
#' training.confounder=data_mbqc_groupA$ConfounderMat[1:679,]
#' label=as.factor(data_mbqc_groupA$Specimen)
#' training.label=label[1:679]
#' testing.label=label[680:848]
#' result_prediction=acPCoA_prediction(training.dist,testing.dist,training.confounder,training.label,testing.label,nPC=2,kernel="linear",bandwidth=NULL,lambdas=seq(0, 5, 0.05),anov=T,perc=0.05)
#' }

acPCoA_prediction <- function(training.dist,testing.dist,training.confounder,training.label,testing.label,nPC=2,kernel="linear",bandwidth=NULL,lambdas=seq(0, 5, 0.05),anov=T,perc=0.05){
  D_star <- as.matrix(testing.dist)
  D <- as.matrix(training.dist)
  n <- nrow(D)
  n_star <- ncol(D_star)
  centre <- function(D, n) {
    One <- matrix(1, n, n)
    mat <- diag(n) - One/n
    mat.cen <- mat %*% D %*% mat
  }
  centre_testing <- function(K,K_star,n,n_star){
    One1 <- matrix(1/n, n, n)
    One2 <- matrix(1/n, n, n_star)
    mat.cen_testing <- K_star - One1%*%K_star - K%*%One2 + One1%*%K%*%One2;
  }
  names <- rownames(D)
  D_colMean0=D-t(matrix(rep(colMeans(D),n),n,n))
  D_bothMean0=D_colMean0-matrix(rep(rowMeans(D_colMean0),n),n,n)
  epsilon <- sqrt(.Machine$double.eps)*1000000
  ## begin
  delta1 <- centre((-0.5 * D^2), n)
  Gram_testing <- centre_testing((-0.5 * D^2),(-0.5 * D_star^2), n,n_star)
  trace <- sum(diag(delta1))
  D.eig <- eigen(delta1,symmetric=TRUE)
  min.eig <- min(D.eig$values)
  zero.eig <- which(abs(D.eig$values) < epsilon)
  #zero.eig <- which(D.eig$values < epsilon)
  D.eig$values[zero.eig] <- 0
  eig <- D.eig$values
  k <- length(which(eig > epsilon))
  pcoa_vectors <- sweep(D.eig$vectors[, 1:k], 2, sqrt(eig[1:k]),
                   FUN = "*")
  vectors_nor <- sweep(D.eig$vectors[, 1:k], 2, 1/sqrt(eig[1:k]),
                       FUN = "*")
  pcoa_vector_testing = t(Gram_testing)%*%vectors_nor
  values <- data.frame(eig[1:k])
  colnames(values) <- c("Eigenvalues")
  rownames(values) <- 1:nrow(values)
  rownames(pcoa_vectors) <- names
  colnames(pcoa_vectors) <- colnames(pcoa_vectors, do.NULL = FALSE, prefix = "Axis.")
  #out <- (list(vector_testing=vector_testing,values = values, vectors = vectors, vectors_nor = vectors_nor))

  if(dim(pcoa_vectors)[2]<nPC){
    nPC=dim(pcoa_vectors)[2]
    }
  result_tune <- acPCAtuneLambda(X=pcoa_vectors, Y=training.confounder, nPC, lambdas, kernel, bandwidth,perc, centerX=T, centerY=T, scaleX=F, scaleY=F,  anov,  quiet=T)
  result <- acPCA(X=pcoa_vectors, Y=training.confounder, nPC, lambda=result_tune$best_lambda, kernel, bandwidth=NULL)
  training.input <- result$Xv
  validation.input <- pcoa_vector_testing%*%result$v
  result_prediction=list("Xv_training"=training.input,"v"=result$v,"Xv_testing"=validation.input,"eigenvalueX"=result$eigenX,"varianceX"=result$varX)
  if(is.null(training.label) || is.null(testing.label)){
    return(result_prediction)
  }else{
    rf_classifier <- randomForest(x=training.input,y=training.label)
    prediction <- predict(rf_classifier,validation.input)
    accuracy <- mean(testing.label==prediction)
    result_prediction$PredictingAccuracy=accuracy
    return(result_prediction)
  }
}
