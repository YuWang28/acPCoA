#' @title Tree data.
#' @description NGS  whole  genome  shotgun (WGS) sequencing data of white oak trees.
#'
#' @format A list with multiple elements, which are
#' \describe{
#'   \item{DistMat.d2star}{the n by n distance matrix. n is the number of samples. d2star distance is applied.}
#'   \item{DistMat.d2shepp}{the n by n distance matrix. n is the number of samples. d2shepp distance is applied.}
#'   \item{DistMat.d2}{the n by n distance matrix. n is the number of samples. d2 distance is applied.}
#'   \item{DistMat.cvtree}{the n by n distance matrix. n is the number of samples. Cvtree distance is applied.}
#'   \item{DistMat.euclidean}{the n by n distance matrix. n is the number of samples. Euclidean distance is applied.}
#'   \item{DistMat.manhattan}{the n by n distance matrix. n is the number of samples. Manhattan distance is applied.}
#'   \item{ConfounderMat}{the n by q confounder matrix}
#'   \item{ContinentalOri}{Samples were divided into three geographic categories according to their continental origins, which are NorthAmerica (NA), West Europe (WE), and East Europe and Asia (EEA).}
#'   \item{batch}{Samples were divided into four batches according to the NCBI BioProject from which they were downloaded and the analysis platforms they used}
#' }
"data_tree"
