#' Tune the lambda parameter in function acPCA
#'
#' @param X the n by p data matrix, where n is the number of samples, p is the number of variables. Missing values in X should be labeled as NA. If a whole sample in X is missing, it should be removed.
#' @param Y the n by q confounder matrix, where n is the number of samples, q is the number of confounding factors. Missing values in Y should be labeled as NA.
#' @param nPC number of principal components to compute.
#' @param lambdas a vector with the tuning parameters, non-negative values. If 0 is not in lambdas, it will be added to lambdas.
#' @param centerX center the columns in X. Default is True.
#' @param centerY center the columns in Y. Default is True.
#' @param scaleX scale the columns in X to unit standard deviation. Default is False.
#' @param scaleY scale the columns in Y to unit standard deviation. Default is False.
#' @param kernel the kernel to use: "linear", "gaussian".
#' @param bandwidth bandwidth h for Gaussian kernel. Optional.
#' @param anov True or False. Whether the penalty term has the between groups sum of squares interpretation. Default is True.
#' @param perc the best lambda is defined to be the smallest lambda with R(lambda)<=perc (if anov=T), or R(lambda)<=perc*R(lambda=0) (if anov=F) in the nPC principal components.
#' @param quiet True or False. Output the progress of the program. Default is False.
#' @return Results for tuning lambda
#' \item{ratio}{R(lambda): a vector with the ratios. Same length as lambdas}
#' \item{best_lambda}{the best lambda after cross-validation}
#' \item{...}{Input parameters for the function}
#' @export
#'
acPCAtuneLambda <- function(X, Y, nPC, lambdas, centerX=T, centerY=T, scaleX=F, scaleY=F, kernel=c("linear", "gaussian"), bandwidth=NULL, anov=T, perc=0.05, quiet=F){
  ####check whether a whole row in X is mis
  Xmis <- apply(X, 1, function(row){sum(!is.na(row))})
  if (sum(Xmis==0)){
    stop(paste("The following samples in X is missing, please remove them in X and Y: rows ", paste(which(Xmis==0), collapse =  " "), sep=""))
  }
  ####check whether the number of samples in X and Y match
  if (dim(X)[1]!=dim(Y)[1]){
    stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
  }
  ####check whether the elements in lambdas are non-negative
  if (sum(lambdas<0)){
    stop("All elements in lambdas should be non-negative")
  }
  lambdas <- sort(lambdas)
  ####add 0 to lambdas
  if (sum(lambdas==0)==0){
    lambdas <- c(0, lambdas)
  }

  nsam <- dim(X)[1]
  p <- dim(X)[2]

  ####if Y is a vector, change it to a matrix
  if (is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
  }

  ####missing data
  X[is.na(X)] <- mean(X, na.rm=T)
  Y[is.na(Y)] <- mean(Y, na.rm=T)
  X <- scale(X, center = centerX, scale = scaleX)
  Y <- scale(Y, center = centerY, scale = scaleY)

  K <- calkernel(Y, kernel, bandwidth)

  ratio <- c()
  for (i in 1:length(lambdas)){
    if (!quiet){
      if (i==length(lambdas)){
        cat("\r", round(i/length(lambdas)*100), "%", "completed\n")
      } else {
        cat("\r", round(i/length(lambdas)*100), "%", "completed")
      }
    }
    lambda <- lambdas[i]
    #result_acpca1 <- eigs_sym(calAv, k = nPC, which = "LA", n = p, args = list(X = X, K = K, lambda = lambda))
    #result_acpca1 <- eigs_sym(crossprod(X, (diag(dim(K)[1])-lambda*K)%*%X), k = nPC, which = "LA", n = p)
    A = calAv(args = list(X = X, K = K, lambda = lambda))
    result_acpca = eigen(A,symmetric=TRUE)
    result_acpca$values = result_acpca$values[1:nPC]
    result_acpca$vectors = result_acpca$vectors[,1:nPC]
    v <- matrix(result_acpca$vectors, ncol = nPC)
    Xv <- X%*%v
    ratio <- rbind(ratio, diag(crossprod(Xv, K%*%Xv))/diag(crossprod(Xv)) )
  }

  if (anov){
    ### the analysis of variance interpretation
    thres <- perc
  } else {
    ### otherwise
    thres <- max(ratio[1,])*perc
  }
  tmp <- which(apply(ratio <= thres, 1, sum)==ncol(ratio))
  if (length(tmp)==0){
    best_lambda <- lambdas[which(rowMeans(ratio)==min(rowMeans(ratio)))]
    warning("lambda is not large enough, increase max(lambdas) or increase perc")
  } else {
    best_lambda <- lambdas[min(tmp)]
  }
  return(list(ratio = ratio, best_lambda = best_lambda, lambdas = lambdas, kernel=kernel, bandwidth=bandwidth))
}
