# Selin Jessa
# April 2020
#
# This set of functions is adapted from Marie Forest and Claudia Kleinman
# and is based off the file creationFDA_GeneLevel.R.
# I have customized the functions and started to document them.

# data_dir <- "/Volumes/My Passport/braindex/brainspan/"
data_dir <- "."
load(file.path(data_dir, "data/CleanedData.genes.RData"))
load(file.path(data_dir, "data/palette_brainspan.Rda"))

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

    return(  list( lst_GE_FD_Deriv = lst_GE_FD_Deriv,
                   lst_GE_FD       = lst_GE_FD,
                   lst_GE_IMP      = lst_GE_IMP,
                   lst_GE_wNA      = lst_GE_wNA,
                   fdNames         = fdNames) )

}
