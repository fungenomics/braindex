# Selin Jessa
# April 2020
#
# This set of functions is adapted from Marie Forest and Claudia Kleinman
# and is based off the file creationFDA_GeneLevel.R.
# I have customized the functions and started to document them.

# data_dir <- "/Volumes/My Passport/braindex/brainspan/"
data_dir <- "."
load(file.path(data_dir, "data/CleanedData.genes.RData"))

library(fda)

prepGene <- function(gene) {

    num_age_label <- c(8, 12, 16, 21, 37, 92, 612, 1132, 1600, 2120)

    # Get the index of the gene in the BrainSpan data
    tmp <- which(geneName == gene)

    geneExons <- gene
    GE_matrix <- list(geneOfInterest=gene,
                      GE_matrix=lnGE_BS[tmp,],
                      reducedBrain=reducedBrain,
                      brainRegionsNamesReduced=brainRegionsNamesReduced,
                      num_age_reduced=num_age_reduced,
                      num_age_label=num_age_label,
                      ageNames=ageNames,
                      donorID_ls=donorID_ls,
                      geneExons=geneExons)

    return(GE_matrix)

}


#' Reform gene expression (gene x sample) into a matrix (brain regions x donor),
#' with NAs for missing combinations
#'
#' @param lnGE_matGene One-row matrix of gene expression (gene x sample) from Brainspan
#' @param reducedBrain Named list, where each element is a consolidated brain region,
#' and contains a vector of indices for samples belonging to that region
#' @param donors List where each element is a donor, and the corresponding
#' element contains a vector of indices for samples from that donor
#'
#' @return A matrix of gene expression (brain region x donor)
createSubMatrix <- function(lnGE_matGene, reducedBrain, donors){

    ret <- matrix(NA, ncol=length(donors), nrow=length(reducedBrain))

    for(i in 1:length(reducedBrain)){
        for(j in 1:length(donors)){
            tmpset <- intersect(reducedBrain[[i]], donors[[j]])
            if( length(tmpset) == 1 ){
                ret[i,j] <- lnGE_matGene[, tmpset ]
            }else if( length(tmpset) > 1 ){
                print(tmpset)
            }
        }
    }

    return(ret)
}



#' Impute NA values in gene expression matrix
#'
#' @param ge_mat Numeric matrix of gene expression as returned by createSubMatrix
#' @param num_age Numeric, ages in post-conception weeks (PCW)
#' @param use.log Logical, whether to use the natural log of the age in PCW
#' @param k Numeric, value of k for k-nearest-neighbours imputation. Default: 6.
#'
#' @return A matrix of gene expression (brain region x donor),
#' with missing values imputed
imputeNA <- function(ge_mat, num_age, use.log = TRUE, k = 6){
    # Not Idiot proof... will work only with matrices wihtout too many NA...
    # It uses the K-nearest neighbor algorithm, with K = 6
    # rows: Brain regions ; cols: age
    ret <- ge_mat

    for(i in 1:nrow(ret)){
        for(j in 1:ncol(ret)){
            if(is.na(ret[i,j])){
                if(use.log){
                    tmp <- idx_KNN(j, log(num_age), ge_mat[i,], k)
                }else{
                    tmp <- idx_KNN(j, num_age, ge_mat[i,], k)
                }
                ret[i,j] <- mean(ge_mat[i,tmp])
            }
        }
    }

    return(ret)
}



#' Perform imputation using k-nearest-neighbours
idx_KNN <- function(idx, x, yNA, K){
    idx_x <- x[idx]
    x_NA <- x
    x_NA[is.na(yNA)]<- NA
    ret <- 1:length(x)
    ret <- ret[order(abs(idx_x - x_NA ))]
    return(ret[1:K])
}



colorsNA <- function(ge_mat, na_col = 'red'){
    ret <- matrix( 'black', nrow = nrow(ge_mat), ncol=ncol(ge_mat))
    ret[is.na(ge_mat)] <- na_col
    return(ret)
}



#' Create functional data objects and performs smoothing
#'
#' @return The element fd from the object of class fdSmooth which contains a
#' smooth of the data (see ?smooth.basis for details)
fdaBrainSpan <- function( logAge, y,  ordr, ordr_rough ,lambda1, fd_names = NULL )  {

    GEdata <- t(y)
    norder <- ordr
    nbasis <- length(logAge) + norder - 2

    # create FD objects using a B-spline basis
    bsp <- create.bspline.basis( rangeval = c(min(logAge),max(logAge)),breaks = logAge, nbasis = nbasis, norder=norder )

    # Define parameters for FD analysis
    sGEfdPar <- fdPar(bsp, ordr_rough, lambda=lambda1)

    # Smooth data using FD analysis
    sGEfd <- smooth.basis(logAge, GEdata, sGEfdPar, fdnames = fd_names)

    # Add column names
    colnames(sGEfd$fd$coefs) <- fd_names$region

    # Return only the smooth
    return(sGEfd$fd)

}



#' Not used here?
gcv_fdaBS <- function(lAge, y, ord, ordr_rough, range_logGCV  ){
    ## Code for the GCV to determine lambda
    loglam <- seq(range_logGCV[1],range_logGCV[2],0.25)
    nlam <- length(loglam)
    dfsave <- rep(NA, nlam)
    gcvsave <- rep(NA, nlam)
    for(ilam in 1:nlam){
        cat(paste('log10 lambda = ', loglam[ilam],'\n'))
        lambda <- 10^loglam[ilam]
        fd_test <- fdaBrainSpan(lAge, y, 5, 3, lambda )
        dfsave[ilam ] <- fd_test$df
        gcvsave[ilam] <- sum(fd_test$gcv)
    }
    ret <- list( loglambda= loglam, gcvVal=gcvsave, dfval=dfsave)
    return(ret)
}



# Plotting ----

# plotFit.wNA <- function( xs, ge_NA, ge_Imputed, ge_fd, titles, xlabel='Log(Weeks after conception)', ylabel="Gene Expression"){
#   # Plot fitted data with imputed NA in red
#   col_points <- colorsNA(ge_NA)
#   for(j in 1:nrow(ge_Imputed)){
#     plot(  xs, ge_Imputed[j,], col=col_points[j,], main = titles[j], lwd=2,  xlab=xlabel, ylab=ylabel )
#     lines(  ge_fd[j] , lwd=2 )
#   }
# }

# plotFDA <- function(lst_GE_wNA, lst_GE_IMP, lst_GE_FD,lst_GE_FD_Deriv, exonOfInt, brNamesReduced, age_label,numAgeReduced, plotsPath){
#
#   if (!dir.exists(file.path(plotsPath, "Gene_FDA_curves/"))) {
#     dir.create(file.path(plotsPath, "Gene_FDA_curves/"))
#   }
#   if (!dir.exists(paste0(plotsPath, "Gene_fitted_curves/"))) {
#     dir.create(file.path(plotsPath, "Gene_fitted_curves/"))
#   }
#   if (!dir.exists(paste0(plotsPath, "Gene_derivatives/"))) {
#     dir.create(file.path(plotsPath, "Gene_derivatives/"))
#   }
#
#   for (i in 1:length(exonOfInt)) {
#
#     # For the plot of smooth curves, we can see what the method for plot.fd()
#     # is doing here: https://github.com/cran/fda/blob/master/R/plot.fd.R#L18
#     pdf( file.path(plotsPath, paste0('Gene_FDA_curves/FDA_plot_gene_', exonOfInt[i],".pdf")), width = 15, height = 10)
#     plot( lst_GE_FD[[i]] , main=paste("Gene:", exonOfInt[i] ,sep=" "), lwd=2, xaxt="n", xlab='Weeks after conception', ylab="Gene Expression")
#     legend("topright", brNamesReduced, lwd=2,col = 1:6, lty = 1:5)
#     axis(1, at=log(numAgeReduced), labels=rep('',length(numAgeReduced)))
#     num_age_label <- c(8, 12, 16, 21, 37, 92, 612, 1132, 1600, 2120)
#     axis(1, at=log(num_age_label), labels=num_age_label)
#     dev.off()
#
#     pdf( file.path(plotsPath, paste0('Gene_fitted_curves/Fitted_curve_gene_',exonOfInt[i], ".pdf")), width = 14, height = 14)
#     par(mfrow=c(4,4))
#     titres <- c()
#     for(j in 1:length(brNamesReduced)){
#       titres <- c(titres,paste("Gene:",exonOfInt[i],"--- B.R.:",brNamesReduced[j] ,sep=" ") )
#     }
#     plotFit.wNA(log(numAgeReduced), lst_GE_wNA[[i]], lst_GE_IMP[[i]], lst_GE_FD[[i]],  titles = titres )
#     dev.off()
#
#     pdf( file.path(plotsPath, paste0('Gene_derivatives/FDA_derivative_gene_',exonOfInt[i], ".pdf")), width = 15, height = 10)
#     plot(lst_GE_FD_Deriv[[i]], lwd=2, xaxt="n", xlab='Weeks after conception', ylab="Gene Expression", main=paste("FDA derivative Gene:",exonOfInt[i] ,sep=" "))
#     legend("topright", brNamesReduced, lwd=2,col = 1:6, lty = 1:5)
#     axis(1, at=log(numAgeReduced), labels=rep('',length(numAgeReduced)))
#     axis(1, at=log(age_label), labels=age_label)
#     dev.off()
#   }
# }



# Wrapper ----

createFDAcurves <- function(info, lambda_opt = 10^-2){

    # Gene name
    geneExonsNames <- info$geneExons

    # Initialize empty vectors for gene expression
    lst_GE_wNA <- vector("list", length = length(geneExonsNames))
    lst_GE_IMP <- vector("list", length = length(geneExonsNames))
    lst_GE_FD <- vector("list", length = length(geneExonsNames))
    lst_GE_FD_Deriv <- vector("list", length = length(geneExonsNames))

    for(i in 1:length(geneExonsNames)){

        # Get age labels, brain region names, and gene name
        fdNames <- list( age = info$ageNames,
                         region = info$brainRegionsNamesReduced,
                         gene = geneExonsNames[i] )

        # Transform gene expression row into matrix with NAs for missing brain regions / timepoints
        lst_GE_wNA[[i]] <- createSubMatrix(info$GE_matrix[i,], info$reducedBrain,  info$donorID_ls)

        # Impute NA values using k-nearest-neighbours algorithm
        lst_GE_IMP[[i]] <- imputeNA(lst_GE_wNA[[i]], info$num_age_reduced)

        # Perform functional data analysis and get a smooth of the data
        lst_GE_FD[[i]] <- fdaBrainSpan(  log(info$num_age_reduced), lst_GE_IMP[[i]], 5, 3, lambda_opt, fd_names = fdNames )

        # Compute the derivative of the functional data object
        lst_GE_FD_Deriv[[i]] <- deriv.fd(lst_GE_FD[[i]] )

    }

    fdNames$age_pcw <- info$num_age_reduced

    # plotFDA(lst_GE_wNA,
    #         lst_GE_IMP,
    #         lst_GE_FD,
    #         lst_GE_FD_Deriv,
    #         exonOfInt = geneExonsNames,
    #         info$brainRegionsNamesReduced,
    #         info$num_age_label,
    #         info$num_age_reduced,
    #         plotsPath)

    return(  list( lst_GE_FD_Deriv = lst_GE_FD_Deriv,
                   lst_GE_FD       = lst_GE_FD,
                   lst_GE_IMP      = lst_GE_IMP,
                   lst_GE_wNA      = lst_GE_wNA,
                   fdNames         = fdNames) )

}
