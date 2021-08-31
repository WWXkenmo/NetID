#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
# loading Parameters
if(length(args)!=10) {
    stop("Rscript GRN_build_pip.R ExpressionProfile Number_of_K(KNN), Number_of_Seed_Cell, If_Identity, Threshold_num, id, TI_method, varID_file, out_dir")
} else if (length(args == 10)) {
    cat(paste("Input: ",args[1]," \n", "The Number of NN(K): ", args[2], " \n", "The Number of Sampled Cell: ", args[3], " \n", "Whether Running Identity Mode: ", args[4], " \n", "Threshold Number of Cells: ",args[5]," \n","id: ",args[6]," \n", "TI method: ", args[7], " \n","Species: ", args[8],"\n","VarID_file: ",args[9],"\n","Save: ",args[10],sep=""))
}
input <- args[1];
k_n <- as.integer(args[2]);
sample_n <- as.integer(args[3]);
ident <- as.logical(args[4]);
Threshold_Num <- as.integer(args[5]);
id <- args[6];
TI <- args[7];
species <- as.integer(args[8]);
VarID_file <- args[9];
save <- args[10];



cat("Loading Require Packages\n")
####################################################
###Try to figure out while using Seurat produced KNN generate better KNN compare to the RaceID
require(Rcpp)
suppressPackageStartupMessages(require(igraph))
suppressPackageStartupMessages(require(SingleCellExperiment))
require(reticulate)
require(rARPACK)
sourceCpp("/mnt/data1/weixu/NetID/RaceID_VarID_code/Rcpp_code/VarID.cpp")
source("/mnt/data1/weixu/NetID/RaceID_VarID_code/R_code/RaceID_utils.R")
require(RaceID)
suppressPackageStartupMessages(require(Matrix))
suppressPackageStartupMessages(require(scater))
require(parallel)
require(compiler)
require(quadprog)
require(irlba)
suppressPackageStartupMessages(require(FNN))
require(MASS)
source("/mnt/data1/weixu/NetID/RaceID_VarID_code/R_code/VarID_functions.R")
require(runner)
require(matrixStats)

FilterOverLap <- function(nn,threshold=0.7,id=id){
    ####Filter the seed node based on shared NN
    JacDist <- function(Dist.m){
            R1 <- Dist.m %*% t(Dist.m)
            B <- rowSums(Dist.m)
            A <- matrix(rep(B,length(B)),nrow=length(B))
            A <- A + t(A)
            Dist <- R1/(A-R1)
            Dist
        }
        ####Transform nn into binary network
    binN <- matrix(0, nrow=ncol(nn),ncol=ncol(nn))
    for(i in 1:nrow(t(nn))){
        binN[i,t(nn)[i,]] <- 1 
    }
    colnames(binN) <- rownames(binN) <- colnames(nn)
    binN <- as(binN,"dgCMatrix")
    d = JacDist(binN)

    d_seed = d[id,id]
    
    ###take the upper tri matrix
    d_seed[!upper.tri(d_seed, diag = FALSE)] <- 0
    d_seed[d_seed<threshold] <- 0
    d_seed[d_seed!=0] <- 1

    ###select out the index
    d_seed <- as(as.matrix(d_seed), "dgCMatrix")
    d_seed <- summary(d_seed)
    d_seed <- d_seed[,-3]
    
    ###filter out seed cell with multiple_overlap
    #ind1 <- as.integer(names(table(d_seed[,1][d_seed[,1]>=2])))
    #ind2 <- as.integer(names(table(d_seed[,2][d_seed[,2]>=2])))
    seed_out <- c(d_seed[,1],d_seed[,2])
    seed_out <- as.integer(names(table(seed_out)))
    d_seed_out <- d[id,id][seed_out,seed_out]
    seed_out <- colnames(d_seed_out)
    
    ###Cluster the cells
    d_g <- igraph::graph_from_adjacency_matrix(d_seed_out,weighted = TRUE,mode="undirected")
    mc <- igraph::walktrap.community(d_g,steps=5)
    mc <- igraph::communities(mc)
    
    ###order the connection of seed cell with other cells
    con_d <- colMeans(d[id,id])
    
    #Filtering
    seed_in <- NULL
    for(i in 1:length(mc)){
        con_d_sub <- con_d[mc[[i]]]
        seed_in <- c(seed_in,names(con_d_sub)[which.min(con_d_sub)])
    }
    seed_out <- seed_out[seed_out %in% seed_in == FALSE]

    ####Return Seed Cell
    id <- id[id %in% seed_out == FALSE]
    id
}

pruneKnn2 <- function(expData,distM=NULL,large=TRUE,regNB=TRUE,bmethod=NULL,batch=NULL,regVar=NULL,offsetModel=TRUE,thetaML=FALSE,theta=10,ngenes=2000,span=.75,pcaComp=NULL,algorithm="kd_tree",metric="pearson",genes=NULL,knn=25,alpha=1,nb=3,SNN=TRUE,no_cores=NULL,FSelect=FALSE,seed=12345,...){

    expData <- as.matrix(expData)
    
    rs <- rowSums(expData > 0)
    cs <- colSums(expData)
    expData <- expData[rs>0,cs>0]

    
    if (!is.null(batch) ) batch <- batch[colnames(expData)]    
    if ( is.null(genes) ) genes <- rownames(expData)
    hflag <- FALSE
    if ( !is.null(batch) & !is.null(bmethod) ){
        if ( bmethod == "harmony" ){
            hflag  <- TRUE
            hbatch <- batch
            batch  <- NULL
        }
    }
    
    expData <- expData[genes,]
    colS    <- cs[ cs > 0 ]
    Xpca    <- NULL
    bg      <- NULL

    if ( large ){
        distM <- NULL
        
        if ( regNB ){
            if ( offsetModel ){
                regData <- compResiduals0(expData[genes,],batch=batch,regVar=regVar,span=span,no_cores=no_cores,ngenes=ngenes,seed=seed,thetaML=thetaML,theta=theta)
            }else{
                regData <- compResiduals(expData[genes,],batch=batch,regVar=regVar,span=span,no_cores=no_cores,ngenes=ngenes,seed=seed)
            }
            z <- regData$pearsonRes
        }else{
            regData <- NULL
            z <- t(t(expData)/colS*min(colS))
        }
        
        if ( FSelect ){
            bg <- fitBackVar(expData[genes,])
            backModel <- bg$fit
            genes     <- bg$genes
            expData   <- expData[genes,]
        }
        
        z <- z[genes,]
        f <- rowSums(is.na(z)) == 0 
        z <- z[f,]
        z <- t(apply(z,1,scale))
        set.seed(seed)

        if ( !is.null(pcaComp) ){
            pcaComp <- min( pcaComp, ncol(expData) - 1)
            pcaComp <- min( pcaComp, nrow(expData) - 1)
            Xpca <- irlba(A = t(z), nv = pcaComp)
        }else{
            pcaComp <- 100
            pcaComp <- min( pcaComp, ncol(expData) - 1)
            pcaComp <- min( pcaComp, nrow(expData) - 1)
            Xpca <- irlba(A = t(z), nv = pcaComp)
            
            
            g <- Xpca$d/sum(Xpca$d)
            g <- mean_run(g,3)
            y <- g[ -length(g) ] - g[-1]
            mm <- numeric(length(y))
            nn <- numeric(length(y))
            for ( i in 1:length(y)){
                mm[i] <- mean(y[i:length(y)]) 
                nn[i] <- sqrt(var(y[i:length(y)]))
            }
            ind <- which( y - (mm + nn) < 0 )
            for ( i in ind ) { if ( sum( i:(i + 3) %in% ind ) >= 2 ){ pcaComp <- i; break} }
            pcaComp <- max( 15, pcaComp )
        }
        dimRed <- Xpca$u[,1:pcaComp]%*%diag(Xpca$d[1:pcaComp])
        if ( hflag ){
            dimRed <- t(dimRed)
            colnames(dimRed) <- colnames(expData)
            dimRed <- HarmonyMatrix( dimRed, hbatch ,do_pca=FALSE,...)
            nn     <- get.knn(t(dimRed), k=knn, algorithm=algorithm)
            nn     <- t( cbind( 1:ncol(expData),nn$nn.index) )
            colnames(nn) <- colnames(expData)
        }else{
            nn     <- get.knn(dimRed, k=knn, algorithm=algorithm)
            nn     <- t( cbind( 1:ncol(expData),nn$nn.index) )
            dimRed <- t(dimRed)
            colnames(nn) <- colnames(dimRed) <- colnames(expData)
        }
    }else{
        if ( FSelect ){
            genes <- bg$genes
            expData <- expData[genes,]
        }
        dimRed <- NULL
        if ( is.null(distM) ) distM <- dist.gen(t(as.matrix(expData[genes,])), method = metric)
        maxD <- 2
        if ( metric == "euclidean" ) maxD <- max(distM)
        nn <- apply(distM,1,function(x){ j <- order(x,decreasing=FALSE); head(j,knn + 1); } )
        regData <- NULL
    } 
    
    cQP      <- cmpfun(QP)
    if(SNN == TRUE){
        ###calculating SNN and pruning on the SNN
        JacDist <- function(Dist.m){
            R1 <- Dist.m %*% t(Dist.m)
            B <- rowSums(Dist.m)
            A <- matrix(rep(B,length(B)),nrow=length(B))
            A <- A + t(A)
            Dist <- R1/(A-R1)
            Dist
        }
        ####Transform nn into binary network
        binN <- matrix(0, nrow=ncol(nn),ncol=ncol(nn))
        for(i in 1:nrow(t(nn))){
            binN[i,t(nn)[i,]] <- 1 
        }
        colnames(binN) <- rownames(binN) <- colnames(nn)
        binN <- as(binN,"dgCMatrix")
        d = JacDist(binN)

        ##Reconstruct NN object
        snn <- NULL
        for(i in 1:ncol(d)){
            ind <- order(d[,i],decreasing=TRUE)[1:(k+1)]
            snn <- cbind(snn, ind)
        }
        colnames(snn) <- colnames(nn)
        nn <- snn   
    }
    
    localFUN <- function(x,expData,colS,alpha,nb,cQP){
        y    <- as.matrix( expData[,x] )
        fit  <- fitLogVarLogMean(y)
        #fit <- backModel
        FNData <- t(t(y)/colS[x]*min(colS))
        k  <- FNData[,1]
        m  <- FNData[,-1]
        k0 <- as.vector(y[,1])
        m0 <- y[,-1]
        weights <- tryCatch(round(cQP(k,m,TRUE)$w,5), error = function(err){ rep(1/ncol(m),ncol(m)) } )
        
        if ( is.null(alpha) ){
            u    <- apply(y[,-1],1,function(x,w){sum(x * w)},w = weights)
            v    <- sqrt( lvar(y[,1],fit) )
            W    <- sum(weights)
            b1   <- max( (  u - ( y[,1] + v ) * W )/v, na.rm=TRUE)
            b2   <- max( ( -u + ( y[,1] - v ) * W )/v, na.rm=TRUE)
            lb   <- 0  
            Dmat <- matrix(1)
            dvec <- 0
            Amat <- matrix( c(1,1,1), nrow=1)
            bvec <- c( b1, b2, lb)
            
            suppressWarnings( opt <- tryCatch( {
                rs <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 0, factorized=FALSE )
                TRUE
            }, error = function(err){ FALSE } ))
            
            if ( opt ) alpha <- rs$solution else alpha <- 1
            if ( alpha == Inf ) alpha <- 1 
        }
    
        weights <- c(alpha,weights)
        weights <- weights/sum(weights)
        z <- applyProb(y + 1,coefficients(fit),weights)
        rownames(z) <- rownames(y)
        colnames(z) <- colnames(y)
        
        p <- apply(z,2,function(x){ exp( mean(log( p.adjust(x[order(x,decreasing=FALSE)],method="bonferroni")[1:nb] + 1e-16 ) ) ) })[-1]
        names(p) <- colnames(m)
        c(alpha,p)
    }

    localFUN2 <- function(x,expData,colS,alpha,nb,cQP){
        y    <- as.matrix( expData[,x] )
        fit  <- fitLogVarLogMean(y)
        #fit <- backModel
        FNData <- t(t(y)/colS[x]*min(colS))
        k  <- FNData[,1]
        m  <- FNData[,-1]
        k0 <- as.vector(y[,1])
        m0 <- y[,-1]
        weights <- tryCatch(round(cQP(k,m,TRUE)$w,5), error = function(err){ rep(1/ncol(m),ncol(m)) } )
        
        if ( is.null(alpha) ){
            u    <- apply(y[,-1],1,function(x,w){sum(x * w)},w = weights)
            v    <- sqrt( lvar(y[,1],fit) )
            W    <- sum(weights)
            b1   <- max( (  u - ( y[,1] + v ) * W )/v, na.rm=TRUE)
            b2   <- max( ( -u + ( y[,1] - v ) * W )/v, na.rm=TRUE)
            lb   <- 0  
            Dmat <- matrix(1)
            dvec <- 0
            Amat <- matrix( c(1,1,1), nrow=1)
            bvec <- c( b1, b2, lb)
            
            suppressWarnings( opt <- tryCatch( {
                rs <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 0, factorized=FALSE )
                TRUE
            }, error = function(err){ FALSE } ))
            
            if ( opt ) alpha <- rs$solution else alpha <- 1
            if ( alpha == Inf ) alpha <- 1 
        }
    
        weights <- c(alpha,weights)
        weights <- weights/sum(weights)
        z <- applyProb(y + 1,coefficients(fit),weights)
        rownames(z) <- rownames(y)
        colnames(z) <- colnames(y)
        
        p <- apply(z,2,function(x){ exp( mean(log( x[order(x,decreasing=FALSE)][1:nb] + 1e-16 ) ) ) })[-1]
        names(p) <- colnames(m)
        c(alpha,p)
    }

    out <- apply(t(nn),1,localFUN,expData=expData,colS=colS,alpha=alpha,nb=nb,cQP=cQP)
    out2 <- apply(t(nn),1,localFUN2,expData=expData,colS=colS,alpha=alpha,nb=nb,cQP=cQP)

    colnames(out) <- colnames(nn)
    colnames(out2) <- colnames(nn)
    pars <- list(large=large,regNB=regNB,offsetModel=offsetModel,thetaML=thetaML,theta=theta,ngenes=2000,span=.75,pcaComp=pcaComp,algorithm=algorithm,metric=metric,genes=genes,knn=knn,alpha=alpha,nb=nb,no_cores=no_cores,FSelect=FSelect,seed=seed)
    return(list(distM=distM,dimRed=dimRed,pvM=out[-1,],pvM.raw=out2[-1,],NN=nn,B=bg,regData=regData,pars=pars,alpha=out[1,],pca=Xpca))
}
getPPI_String <- function (object = NULL, species = 9606, score_threshold = 600, 
    save = FALSE) 
{
	suppressPackageStartupMessages(require("data.table"))
	suppressPackageStartupMessages(require("igraph"))
	
    linkFiles <- paste("https://stringdb-static.org/download/protein.links.v11.0/", 
        species, ".protein.links.v11.0.txt.gz", sep = "")
    if (!file.exists(sub(pattern = ".gz", replacement = "", x = basename(linkFiles)))) {
        if (!file.exists(basename(linkFiles))) 
            download.file(linkFiles, destfile = basename(linkFiles))
        gf <- gzfile(basename(linkFiles), "rt")
    }
    PPI <- read.table(gf, header = T, sep = "")
	PPI[,1] <- as.factor(PPI[,1])
	PPI[,2] <- as.factor(PPI[,2])
    close(gf)
    infoFiles <- paste("https://stringdb-static.org/download/protein.info.v11.0/", 
        species, ".protein.info.v11.0.txt.gz", sep = "")
    if (!file.exists(sub(pattern = ".gz", replacement = "", x = basename(infoFiles)))) {
        if (!file.exists(basename(infoFiles))) 
            download.file(infoFiles, destfile = basename(infoFiles))
        gf <- gzfile(basename(infoFiles), "rt")
    }
    Pinfo <- read.table(gf, header = T, sep = "\t", colClasses = c("character", 
        "character", "NULL", "NULL"), quote = "", row.names = 1)
    close(gf)
    PPI <- subset(PPI, combined_score > score_threshold)
    ENSP1 <- levels(PPI[, 1])
    levels(PPI[, 1]) <- toupper(Pinfo[ENSP1, ])
    ENSP2 <- levels(PPI[, 2])
    levels(PPI[, 2]) <- toupper(Pinfo[ENSP2, ])
    if (!is.null(object)) {
        gene_data <- rownames(object)
        gene_data_upper <- toupper(gene_data)
        gene_data <- as.data.frame(unique(as.data.table(data.frame(gene_data, 
            gene_data_upper)), by = "gene_data_upper"))
        rownames(gene_data) <- gene_data[, 2]
        PPI <- PPI[which(is.element(PPI[, 1], gene_data[, 2])), 
            ]
        PPI <- PPI[which(is.element(PPI[, 2], gene_data[, 2])), 
            ]
        levels(PPI[, 1]) <- gene_data[levels(PPI[, 1]), 1]
        levels(PPI[, 2]) <- gene_data[levels(PPI[, 2]), 1]
    }
    nodes <- union(PPI[, 1], PPI[, 2])
    links <- PPI[, 1:2]
    net <- graph_from_data_frame(d = links, vertices = nodes, 
        directed = FALSE)
    net <- igraph::simplify(net)
    if (save) {
        saveRDS(as_adj(net), paste(species, "_ppi_matrix_STRING-11.0.Rda", 
            sep = ""))
    }
    file.remove(paste(species, ".protein.links.v11.0.txt.gz", 
        sep = ""))
    file.remove(paste(species, ".protein.info.v11.0.txt.gz", 
        sep = ""))
    return(as_adj(net))
}

scoreMeans <- function(score, indx){
    s <- NULL
    for(i in names(table(indx))){
      s <- c(s, mean(score[indx == i]))  
    }
    names(s) <- names(table(indx))
    s
}
####load seurat object for comparsion latter
cat('Running RaceID and VarID to process dataset and construct cell KNN graph...')
require(RaceID)
## load input data from https://github.com/dgrun/VarID_analysis/blob/master/inputData_hematopoiesis.rds
d  <- count <- readRDS(input)
sc <- SCseq(d)
sc <- filterdata(sc,mintotal=1000,CGenes=rownames(intestinalData)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(intestinalData))])
expData  <- getExpData(sc)

if("VarID_res.rds" %in% list.files(VarID_file)){
    res <- readRDS(paste(VarID_file,"VarID_res.rds",sep=""))}else{
    cat("\nNo VarID file, running VarID...")
    res   <- pruneKnn2(expData,knn=(k_n-1),no_cores=1,SNN=FALSE)
}
cl    <- graphCluster(res,pvalue=0.01)
probs <- transitionProbs(res,cl)
sc <- updateSC(sc,res=res,cl=cl)
sc <- compumap(sc)

###calculate pesudotime
#sce <- SingleCellExperiment(assays = list(logcounts = res$regData$pearsonRes))
#diff <- DiffusionMap(t(as.matrix(logcounts(sce))), verbose = TRUE,n_pcs=20,n_local=20)
#dpt <- DPT(diff,tips = 4686) # the root cell is determined by differentiation potency
#plotmap(sc)
#plotmap(sc,um=TRUE)

#reducedDims(sce) <- list(umap=sc@umap)
#sce$label <- as.factor(sc@cluster$kpart)
#sce$dpt <- dpt$dpt
#gg1 <- plotReducedDim(sce, dimred = "umap", colour_by = "label",text_by = "label")
#gg2 <- plotReducedDim(sce, dimred = "umap", colour_by = "dpt")
#gg3 <- plotReducedDim(sce, dimred = "umap", colour_by = "scent_score")

#pt <- data.frame(Pesudotime=dpt$dpt)
#rownames(pt) <- colnames(sce)

####Create folder
dir.create(paste(save,id,sep=""))
setwd(paste(save,id,sep=""))
dir.create("raw")
dir.create("netID")

##save VarID object
saveRDS(res, paste("VarID_res.rds",sep=""))
###calculate pesudotime
if(TI != "NO_TI"){
if(TI == "DiffusionMap"){
    cl <- sc@cluster$kpart
    cat('\nUsing SCENT to determine the root cell...\n')
    ###Integrate SCENT to determine the start cell.
    suppressPackageStartupMessages(library(SCENT))
    PPI <- getPPI_String(sc@ndata * min(sc@counts),species=species)
    score <- CompCCAT(log2(as.matrix(sc@ndata[rownames(PPI),] * min(sc@counts))+1),PPI)
    root <- scoreMeans(score, cl)
    cat(paste("\nSCENT predict cluster ",names(root)[which.max(root)]," as the root: ",max(root),"\n",sep="" ))
    
    cat('\nUsing DiffushionMap to model the trajectory...')
    suppressPackageStartupMessages(library(destiny))
    sce <- SingleCellExperiment(assays = list(logcounts = res$regData$pearsonRes))
    diff <- DiffusionMap(t(as.matrix(logcounts(sce))), verbose = TRUE,n_pcs=20,n_local=20)
    dpt <- DPT(diff,tips = which.max(score)) # the root cell is determined by differentiation potency
    #plotmap(sc)
    #plotmap(sc,um=TRUE)
    UMAP <- as.matrix(sc@umap)
    colnames(UMAP) <- c("UMAP1","UMAP2")
    reducedDims(sce) <- list(umap=UMAP)
    sce$label <- as.factor(sc@cluster$kpart)
    sce$dpt <- dpt$dpt
    sce$scent_score <- score
    gg1 <- plotReducedDim(sce, dimred = "umap", colour_by = "label",text_by = "label")
    gg2 <- plotReducedDim(sce, dimred = "umap", colour_by = "dpt")
    gg3 <- plotReducedDim(sce, dimred = "umap", colour_by = "scent_score")
    
    Path <- paste("DiffusionMap.tiff",sep="")
    tiff(Path,res=300,height=3000,width=3000)
    gg2
    dev.off()
    
    pt <- data.frame(Pesudotime=dpt$dpt)
    rownames(pt) <- colnames(sce)
    
    t <- dpt$dpt
    Y <- sc@ndata * 10^4
    Y <- Y[rowMeans(Y)>0.01,]
    #Y <- Y[VariableFeatures(hp),]
    #gam.pval <- apply(Y,1,function(z){
    #    d <- data.frame(z=z, t=t)
    #    suppressWarnings({
    #        tmp <- mgcv::gam(z ~ s(t), data=d)
    #        })
    #        r <- summary(tmp)
    #        p <- r$s.table[4]
    #        p
    #    })
    gamFit <- function(z,pt){
    d <- data.frame(x=z, t=pt)
    suppressWarnings({
        tmp <- mgcv::gam(z ~ s(t), data=d)
    })
    r <- summary(tmp)
    p <- r$s.table[4]
    p
    }
    suppressPackageStartupMessages(library(parallel))
    cat("\nPerform GAM Fitting...")
	cat(paste("\nUsing ",6," cores...",sep=""))
    cl <- makeCluster(6)
    clusterExport(cl, 'Y')
    clusterExport(cl, 't')
    clusterExport(cl,'gamFit')
    gam.pval <- parApply(cl,X=Y,1,gamFit,pt=t)
	stopCluster(cl)

    gam.res <- cbind(gam.pval,apply(Y,1,var))
    gam.res <- gam.res[is.nan(gam.pval)==FALSE,]
    colnames(gam.res) <- c("GAMpvalue","Var")
    gam.res <- as.data.frame(gam.res)
    gam.res <- gam.res[order(gam.res[,2],decreasing=TRUE),]
}
#############################
#running slingshot algorithm to predict trajectory
#rd <- sc@umap
if(TI == "Slingshot"){
    cl <- sc@cluster$kpart
    cat('\nUsing SCENT to determine the root cell...\n')
    ###Integrate SCENT to determine the start cell.
    suppressPackageStartupMessages(library(SCENT))
    PPI <- getPPI_String(sc@ndata * min(sc@counts),species=species)
    score <- CompCCAT(log2(as.matrix(sc@ndata[rownames(PPI),] * min(sc@counts))+1),PPI)
    root <- scoreMeans(score, cl)
    cat(paste("\nSCENT predict cluster ",names(root)[which.max(root)]," as the root: ",max(root),"\n",sep="" ))
    
    cat('\nUsing slingshot to model the trajectory...')
    suppressPackageStartupMessages(library(slingshot))

    sce$label <- as.factor(cl)
    sce <- slingshot(sce, clusterLabels = 'label', reducedDim = 'umap',start.clus = names(root)[which.max(root)])
    
    Path <- paste("slingshot.tiff",sep="")
    tiff(Path,res=300,height=3000,width=3000)
    plot(reducedDims(sce)$umap, col = as.factor(cl), pch=16, asp = 1)
    lines(SlingshotDataSet(sce), lwd=2, col = 'black')
    dev.off()

    cat('\nFitting GAM...')
    suppressPackageStartupMessages(library(vgam))
    t <- pt <- slingPseudotime(sce)
    # for time, only look at the 100 most variable genes 
    Y <- sc@ndata * 10^4
    Y <- Y[rowMeans(Y)>0.01,]
    # fit a GAM with a loess term for pseudotime
    vgam.res <- NULL
    for(ind in 1:nrow(Y)){
    z <- Y[ind,]
    d <- data.frame(cbind(t,z))
    if(length(table(na.omit(d)$z))>7){
    name <- NULL
    for(i in 1:(ncol(d)-1)){
        name <- c(name,paste("pt_",i,sep=""))
    };name_full <- c(name,"gene")
    colnames(d) <- name_full
    formula <- "cbind(pt_1"
    if(length(name)>1){
        for(i in 2:length(name)){
            formula <- paste(formula,",",name[i],sep="")
        }
    };formula <- as.formula(paste(formula,") ~ s(gene)",sep=""))
    suppressWarnings({
      tmp <- vgam( formula, gaussianff,data=d)
    })
    r <- summaryvgam(tmp)
    p <- min(na.omit(r@anova$`P(Chi)`))
    var <- var(z)
    line <- c(p,var)
    line <- matrix(line,nrow=1,ncol=2)
    rownames(line) <- rownames(Y)[ind]
    vgam.res <- rbind(vgam.res,line)
    }
}
colnames(vgam.res) <- c("VGAMpvalue","Variance")
gam.res <- vgam.res
gam.res <- gam.res[order(gam.res[,2],decreasing=TRUE),]
}
if(TI != "DiffusionMap" & TI != "Slingshot"){
    cat("\nError: TI method need to be diffusionmap or slingshot")
    break
}
}
if(TI != "NO_TI") write.csv(pt,"./raw/Pseudotime.csv")
if(TI != "NO_TI") write.csv(gam.res,"./raw/GeneOrdering.csv")
if(TI != "NO_TI") write.csv(pt,"./netID/Pseudotime.csv")
if(TI != "NO_TI") write.csv(gam.res,"./netID/GeneOrdering.csv")

## prepare normalized expression data for GRN inference
cat('\nPerform geosketching...')
X <- apply(res$regData$pearsonRes,1,scale)
geosketch <- import('geosketch')
s <- svds(X, k=50)
X.pcs <- s$u

sketch.size <- as.integer(sample_n)
sketch.indices <- geosketch$gs(X.pcs, sketch.size, one_indexed = TRUE)

## prepare normalized expression data for GRN inference
exp <- sc@ndata * 10^4

require(Matrix)
pvalue <- 0.01
id <- unlist(sketch.indices)
x  <- t(res$NN)[id,]
y  <- Matrix(rep(0,ncol(res$NN)*length(id)), ncol=ncol(res$NN))
rownames(y) <- rownames(x)
colnames(y) <- colnames(res$NN)

## prune sampled neighbourhoods
for ( i in rownames(y) ){
    p <- res$pvM[,i]
    p[p < pvalue] <- 0
    y[i,res$NN[,i]] <- c(1,p)
}
## pruned p-value matrix. p-values of pruned neighbours are set to 0.
y <- as.matrix( t(y) )
y.prune <- y
y.prune[y.prune!=0] <- 1

pvalue <- 0
x  <- t(res$NN)[id,]
y  <- Matrix(rep(0,ncol(res$NN)*length(id)), ncol=ncol(res$NN))
rownames(y) <- rownames(x)
colnames(y) <- colnames(res$NN)

## prune sampled neighbourhoods
for ( i in rownames(y) ){
    p <- res$pvM.raw[,i]
    p[p < pvalue] <- 0
    y[i,res$NN[,i]] <- c(1,p)
}

## pruned p-value matrix. p-values of pruned neighbours are set to 0.
y <- as.matrix( t(y) )
y.unprune <- y

y.weighted <- y.unprune * y.prune

cs <- rowSums(y.weighted)
f <- cs > 0
y.weighted <- y.weighted[f,]
count <- apply(y.weighted,2,function(x){sum(x!=0)})
y.final <- matrix(0,nrow=nrow(y.weighted),ncol=ncol(y.weighted))
colnames(y.final) <- colnames(y.weighted)
rownames(y.final) <- rownames(y.weighted)
id <- rownames(y.weighted)[rownames(y.weighted) %in% colnames(y.weighted) == FALSE]
for(i in id){
    vec <- which(y.weighted[i,] == max(y.weighted[i,]))
    vec_count <- count[vec]
    y.final[i,vec[which(vec_count == min(vec_count))]] <- 1
}

if(ident){y.final[colnames(y.final),colnames(y.final)] <- diag(ncol(y.final))}else{
    y.final[colnames(y.final),colnames(y.final)] <- y.weighted[colnames(y.final),colnames(y.final)]
}
count <- function(x) { sum(x!=0)}
y.final <- y.final[,apply(y.final,2,count)>Threshold_Num]
aggregate <- function(ind, exp){
    exp.ind <- exp[,which(ind!=0)]
    ind <- ind[which(ind!=0)]/sum(ind)
    agre <- exp.ind %*% as.matrix(ind)
    agre
}
ks.final <- apply(y.final,2,function(x){  f <- colnames(exp) %in% rownames(y.final)[ x > 0]; rowMeans(as.matrix(exp[,f]))} )
colnames(ks.final) <- colnames(y.final)

####Save The File
write.csv(log2(ks.final+1),"./netID/ExpressionData.csv")
if(TI != "NO_TI") write.csv(log2(ks.final[rownames(gam.res),]+1),"./netID/ExpressionData.csv")
write.csv(log2(as.matrix(exp)+1),"./raw/ExpressionData.csv")
if(TI != "NO_TI") write.csv(log2(as.matrix(exp)[rownames(gam.res),]+1),"./raw/ExpressionData.csv")
cat("\nDone!\n")

#tf <- geneS[geneS %in% TF$TF]
## highlight central cell of each neighbourhood in the map
#types <- sc@cpart * NA
#types[colnames(zr)] <- 1
#plotsymbolsmap(sc,types,um=TRUE)

## ...or Pearson residuals
## exp <- res$regData$pearsonRes

## avergae expression across all neighbourhoods
## keep only transcription factors
#ks <- ks[intersect(rownames(ks),d$tf$Symbol),]
#tf <- geneS[geneS %in% TF$TF]
#system.time({
#g3 <- GENIE3(ks.final,nCores=12,verbose=TRUE,nTrees=500,regulators = rownames(ks.final), targets = rownames(ks.final) )
#}
#)

#system.time({
#g4 <- GENIE3(exp,nCores=12,verbose=TRUE,nTrees=500,regulators = rownames(ks.final), targets = rownames(ks.final))
#}
#)

#ChIP_net_filter <- subset(ChIP_net, Gene1 %in% tf & Gene2 %in% tf == FALSE)
#ChIP_net_filter <- ChIP_net
#S <- NULL
#E <- NULL
#for(i in 1:nrow(ChIP_net_filter)){
#    if(sum(tf == ChIP_net_filter[i,1]) == 1){
#        S <- c(S, ChIP_net_filter[i,1])
#        E <- c(E, ChIP_net_filter[i,2])
#    }
#   if(sum(tf == ChIP_net_filter[i,2]) == 1){
#        S <- c(S, ChIP_net_filter[i,2])
#        E <- c(E, ChIP_net_filter[i,1])
#    }
#}
#ChIP_tf_net <- data.frame(source=S,end=E)

#g4[g4== 0] <- min(g4[g4!= 0])
#g3[g3== 0] <- min(g3[g3!= 0])

#g3_pair <- as(as.matrix(g3), "dgCMatrix")
#g3_pair <- summary(g3_pair)
#g3_pair <- as.data.frame(g3_pair)
#g3_pair$i_name <- rownames(g3)[g3_pair$i]
#g3_pair$j_name <- colnames(g3)[g3_pair$j]

#g3_pair$link <- rep(0, nrow(g3_pair))
#for(i in 1:nrow(ChIP_tf_net)){
#    i_cor <- which(rownames(g3) == ChIP_tf_net[i,1])
#    j_cor <- which(colnames(g3) == ChIP_tf_net[i,2])
#    g3_pair[g3_pair$i == i_cor & g3_pair$j == j_cor,"link"] <- 1
#}

#g4_pair <- as(as.matrix(g4), "dgCMatrix")
#g4_pair <- summary(g4_pair)
#g4_pair <- as.data.frame(g4_pair)
#g4_pair$i_name <- rownames(g4)[g4_pair$i]
#g4_pair$j_name <- colnames(g4)[g4_pair$j]

#g4_pair$link <- g3_pair$link
##############################
#g3_pair <- g3_pair[order(g3_pair$x,decreasing=TRUE),]
#g4_pair <- g4_pair[order(g4_pair$x,decreasing=TRUE),]
