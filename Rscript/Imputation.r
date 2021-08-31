#!/usr/bin/env Rscript

# K-nearest neighbor smoothing for high-throughput single-cell RNA-Seq data
# (R implementation.)

# Authors:
#   Yun Yan <yun.yan@nyumc.org>
#   Florian Wagner <florian.wagner@nyu.edu>
# Copyright (c) 2017, 2018 New York University

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(rsvd))

randomized_pca <- function(tmat, d, seed){
  # @param tmat A non-negative matrix with samples by features
  # @return A matrix with features by samples projected on PCA space
  set.seed(seed)
  #rsvd_obj <- rsvd(scale(tmat, center = TRUE, scale = FALSE), k=d)
  #rsvd_obj$u %*% diag(rsvd_obj$d)
  rpca_obj <- rpca(tmat, k=d, center=T, scale=F, retx=T, p=10, q=7)
  rpca_obj$x
}

normalization_median <- function(mat){
  # Median normalization
  # @param mat A non-negative matrix with genes by samples
  num_transcripts <- Matrix::colSums(mat)
  size_factor <- median(num_transcripts, na.rm = T) / num_transcripts
  t(t(mat) * size_factor)
}

freeman_tukey_transform <- function(mat){
  sqrt(mat) + sqrt(mat + 1)
}

pdist <- function(tmat){
  # @param tmat A non-negative matrix with samples by features
  # @reference http://r.789695.n4.nabble.com/dist-function-in-R-is-very-slow-td4738317.html
  mtm <- Matrix::tcrossprod(tmat)
  sq <- rowSums(tmat^2)
  out0 <- outer(sq, sq, "+") - 2 * mtm
  out0[out0 < 0] <- 0

  sqrt(out0)
}

smoother_aggregate_nearest_nb <- function(mat, D, k){
  # @param mat A matrix in a shape of #genes x #samples.
  # @param D A predefined distance matrix in a shape of #samples x #samples.
  # @param k An integer to choose \code{k} nearest samples (self-inclusive) to
  #  aggregate based on the distance matrix \code{D}. If \code{k} is greater than
  #  #samples, \code{k} is forced to be #samples to continue aggregation.
  sapply(seq_len(ncol(mat)), function(cid){
    nb_cid <- head(order(D[cid, ]), k)
    closest_mat <- mat[, nb_cid, drop=FALSE]
    return(Matrix::rowSums(closest_mat))
  })
}

knn_smoothing <- function(mat, k, d=10, seed=42){
  #' KNN-smoothing on UMI-filtered single-cell RNA-seq data
  #'
  #' @param mat A numeric matrix with gene names on rows and cell names on columns.
  #' @param k Number of nearest neighbours to aggregate.
  #' @param d Number of Principal components.
  #' @param seed Seed number. (default=42)
  #' @return A smoothed numeric matrix.
  #' @examples
  #' X <- matrix(abs(sin(seq(from=1, to=1000, length.out = 1000))),
  #' nrow = 25, byrow = T)
  #' y <- rep(1:4, each=10)
  #' dim(X)
  #' colnames(X) <- as.character(paste0("s", seq_len(ncol(X))))
  #' rownames(X) <- as.character(paste0("g", seq_len(nrow(X))))
  #' S <- knn_smoother(X, k=5)
  #' plot(X[1, ], X[3, ], col=factor(y), main="original")
  #' plot(S[1, ], S[3, ], col=factor(y), main="smoothed")
  #' @export

  cname <- colnames(mat)
  gname <- rownames(mat)

  num_steps <- ceiling(log2(k + 1))
  S <- mat
  for (p in seq(1, num_steps)){
    k_step <- min(2^p - 1, k)
    message(paste0('Step ', p, '/', num_steps, ': ',
                   'Smoothing using k=', k_step))
    Y <- freeman_tukey_transform(normalization_median(S))
    if (! is.null(d)) {
      Y <- t(randomized_pca(t(Y), d=d, seed=seed))
    }
    D <- pdist(t(Y))
    S <- smoother_aggregate_nearest_nb(mat, D, k_step + 1)
  }
  if (! is.null(cname)) colnames(S) <- cname
  if (! is.null(gname)) rownames(S) <- gname

  S
}

#############################################
#Perform MAGIC and Knnimpute

library(Rmagic)
library(Seurat)
library(reticulate)
library(SAVER)
library(RaceID)
# data files
data_dir <- "."
expr_files <- list.files(".", pattern = "ExpressionData.csv", recursive = TRUE, include.dirs = TRUE)

for(file in expr_files){
  dat <- read.csv(paste0(data_dir, "/", file), row.names = 1)
  dat_name <- dirname(file)
  
  print(file)
  
  # remove erccs
  #ercc.idx <- grep("^ERCC-", rownames(dat))
  #if( length(ercc.idx) ) { dat <- dat[-ercc.idx,] }
  
  # Create RaceID object
  sc <- SCseq(dat)
  sc <- filterdata(sc,mintotal=1000,CGenes=rownames(intestinalData)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(intestinalData))])
  Y <- sc@ndata * 10^4

  # run magic
  cat('\nRunning MAGIC')
  MAGIC_data <- magic(as.data.frame(as.matrix(t(Y))))
  MAGIC_data <- log2(t(MAGIC_data$result)+1)
  new_file <- gsub(".csv", "_magic_imputed.csv", file)
  new_file <- paste0(data_dir, "/", new_file)
  write.csv(MAGIC_data,
            file = new_file, quote=F, row.names = T)

  rm(MAGIC_data);gc()
  # run SAVER
  cat('\nRunning SAVER')
  SAVER_data <- saver(dat[,colnames(sc@ndata)], ncores = 12, estimates.only = TRUE)
  SAVER_data <- as.data.frame(log2(SAVER_data+1))
  new_file <- gsub(".csv", "_saver_imputed.csv", file)
  new_file <- paste0(data_dir, "/", new_file)
  write.csv(SAVER_data,
            file = new_file, quote=F, row.names = T)
  
  rm(SAVER_data);gc()
  # run knn-smoothing
  cat('\nRunning knn-smoothing')
  KNN_data <- knn_smoothing(dat[,colnames(sc@ndata)],k=30)
  KNN_data <- as.matrix(KNN_data) %*% diag(10^4/colSums(KNN_data))
  colnames(KNN_data) <- colnames(sc@ndata)
  KNN_data <- as.data.frame(log2(KNN_data+1))
  new_file <- gsub(".csv", "_knn_imputed.csv", file)
  new_file <- paste0(data_dir, "/", new_file)
  write.csv(KNN_data,
            file = new_file, quote=F, row.names = T) 
  
  rm(KNN_data);gc()
}