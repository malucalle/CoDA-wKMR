################################################################################
# It includes the posibility of using cmultRepl() instead of adding a
# pseudocount.
# It also includes the possibility of working with continuous responses.
################################################################################


#------------------------------------------------------------------------------#
#                           AUXILIARY FUNCTIONS
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# NAME: rename_OTU
# FUNCTION: it renames the rows of an phyloseq object not to have repeated
#           names. The idea is to write at genus level for example,
#                       f_Lactobacillales_g_unclassified

# INPUT:    a "phyloseq" object with @tax_table information and the rank
#           given as a number of the column of @tax_table.
# OUTPUT:   the "taxonomyTable" object with non-repeated modified names.
#------------------------------------------------------------------------------#


#' Rename each taxon
#'
#' \code{rename_OTU} assigns a non - repeated name to each bacteria for a given
#'  taxonomic level.
#'
#'
#' @param phy the phyloseq object where the information is contained.
#' @param rank a \code{numeric} value indicating the taxonomic level where the
#' names should be taken. It corresponds with the column number of
#' \code{phy@tax_table} associated to the desired taxonomic rank.
#' @param db the bacterial database used to name each OTU of the phyloseq
#' object: SILVA (\emph{default}) or GreenGenes.
#'
#'
#'
#' @return A vector with the names for each bacteria. It has not repeated names
#' and none of them starts with \emph{Incertae_Sedis} or \emph{unclassified}.
#'
#' @export rename_OTU


  rename_OTU <- function(phy,rank, db ="SILVA"){

    # Test if the objects are from the expected class
      stopifnot(class(phy) == "phyloseq")
      stopifnot(class(rank) == "numeric")
    # Rank abbreviation (Kingdom, Phylum, Class, Order, Family, Genus, Specie)
      Ranking<-c("k","p","c","o","f","g","s")


  ##############################################################################
  # AUXILIAR FUNCTION
  ##############################################################################

    replace_rare <- function(j, rk, tax_table, Nam, Ranking, db){

    # The column of tax_table we are working on
      u <- rk[j]
    # The first previous column which is not "unclassified"
    # If the name is unclassified . . .
      if (tax_table[j,u]=="unclassified"){
        while(tax_table[j,u] %in% c("unclassified", "Incertae_Sedis")){
          u <- u - 1}
      # Modify Nam
        if (db=="SILVA"){
          V<-paste(Ranking[u], tax_table[j,u],
                   paste(unlist(strsplit(Nam[j],split="_"))[-c(1:2)],
                         collapse = "_"),
                   sep="_")
        } else {
          v<-paste(tax_table[j,u],
                   paste(unlist(strsplit(Nam[j],split="_"))[-c(1:2)],
                         collapse = "_"),
                   sep="_")}

    # . . . else if the name is Incertae_Sedis
      } else if (tax_table[j,u]=="Incertae_Sedis"){
        while(tax_table[j,u] %in% c("unclassified", "Incertae_Sedis")){
          u <- u - 1}
        if (db=="SILVA"){
        # Modify Nam
          V<-paste(Ranking[u], tax_table[j,u],
                   paste(unlist(strsplit(Nam[j],split="_")), collapse = "_"),
                   sep="_")
        } else {
          V<-paste(tax_table[j,u],
                   paste(unlist(strsplit(Nam[j],split="_")), collapse = "_"),
                   sep="_")
        }
      }

    # Return V
      return(V)
  }

  ##############################################################################

  # Vector with initial names
    Nam<-phy@tax_table[,rank]
    if (db=="SILVA"){Nam <- paste(Ranking[rank],Nam,sep="_")}
  # Make a copy to work with (Nam2)
    Nam2<-Nam
  # Initial r value (counting the maximum repeated values for a certain name)
    r<-2
  # Initial value for the number of categories used
    i<-1
  # A vector with the rank of the name (initially rank value)
    rk <- rep(rank, length(Nam))

  # While r>1 (while there are repeated names)
    while (r>1){
    # Repeated names
      Rep.Nam<-names(table(Nam)[(table(Nam)>1)])
    # Indices with repeated names
      Rep.Idx<-which(Nam %in% Rep.Nam)
    # Modify rk
      rk[Rep.Idx]<- rk[Rep.Idx] - 1
    # Modify Nam for Rep.Idx
      if(db=="SILVA"){
        Nam[Rep.Idx] <- paste(Ranking[rank-i],phy@tax_table[Rep.Idx,rank-i],
                              Nam2[Rep.Idx],sep="_")
      }else{
        Nam[Rep.Idx] <- paste(phy@tax_table[Rep.Idx,rank-i],
                              Nam2[Rep.Idx],sep="_")
      }

    # Modify r as the number of the maximum repeated name
      r<-max(table(Nam))
      i<-i+1

    }

  # Load library
    library(qdapRegex)
  # Extract the first value for each name
    First.NAM <- unlist(lapply(Nam, function(x)
      unlist(rm_between(x, "_","_", extract = T))[1]))
  # Indices for the names to replace (with the first name as unclassified
  # or Incertae_Sedis)
    IDX.Rep <- which(First.NAM %in% c("unclassified", "Incertae"))
  # If there are unclassified, . . .
    if (length(IDX.Rep) !=0){
      for (i in 1:length(IDX.Rep)){
        Nam[IDX.Rep[i]] <-  replace_rare(IDX.Rep[i],
                                         rk,
                                         phy@tax_table,
                                         Nam,
                                         Ranking,
                                         db)
    }

  }

  return(Nam)
}

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# FUNCTION: w.geoM


# INPUT:
#       w: vector of weights
#       x: compositional vector
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# NAME: w.geoM
# FUNCTION: it computes the weighted - Geometric mean as defined in
#           "Changing the Reference Measure in the Simplex and Its Weighting
#            Effects" - Egozcue and Pawlowsky (2016).
# INPUT:
#       w: vector of weights
#       x: compositional vector

# OUTPUT: the weighted geometric mean of x
#
#------------------------------------------------------------------------------#


#' @title Weighted geometric mean
#' @description  Weighted geometric mean
#'
#' \code{w.geoM} computes the weighted geometric mean as defined in Egozcue and
#' Pawlowsky (2016)
#'
#'
#' @param w vector of weights of the same length as x.
#' @param x vector for which the weighted geometric mean is computed.
#'
#'
#' @return the weighted geometric mean value.
#'
#' @export w.geoM
#'
#' @examples
#'
#' # Build compositional vectors
#'   x <- prop.table(runif(10))

#' # Vector of weights
#'   w <- c(2,rep(1,9))
#'
#' # Weighted geometric mean
#'   w.geoM(w,x)


  w.geoM <- function(w, x){
    # Redefine y
      y <- x/w
    # Sum of weights
      sw <- sum(w)
    # Weighted geometric mean
      w.g <- exp(sum(w*log(y))/sw)

    # Return the value
      return(w.g)

  }

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# NAME: w.dist
# FUNCTION: computes the weighted distance between two vectors
# INPUT:
#       w: vector of weights
#      x1: composition for sample 1
#      x2: composition for sample 2
#
#------------------------------------------------------------------------------#

#' Computes the weighted distance between two compositions
#'
#'
#' @param w vector of weights
#' @param x1 first composition
#' @param x2 second composition
#'
#' @return a value with the weighted distance between x1 and x2
#'
#' @examples
#'
#' # Build compositional vectors
#'   x1 <- prop.table(runif(10))
#'   x2 <- prop.table(runif(10))
#'
#' # Vector of weights
#'   w <- c(rep(1,4),3,3,rep(1,4))
#'
#' # The distance
#'    w.dist(w,x1,x2)
#' @export plot.tab


  w.dist <- function(w, x1, x2){
    # Define y1 and y2
      y1 <- x1/w
      y2 <- x2/w
    # Define the distance directly
      w.d <- sqrt(sum(w*(log(y1/w.geoM(w,x1)) - log(y2/w.geoM(w,x2)))^2))
    # Return the value
      return(w.d)
  }

#------------------------------------------------------------------------------#

################################################################################
#                                 MAIN FUNCTIONS
################################################################################

#------------------------------------------------------------------------------#
# NAME: wcomp.MiRKAT
# FUNCTION: it computes the p-value for the incidence of the microbiome over the
#           response variable, given a vector of weights.
# INPUT:
#           x: matrix with the composition for each individual
#           y: response variable: dichotomous or continuous
#           w: vector of weights for each taxa
#         cov: covariates to have into account as a possible confounding
#              factors



# NOTE: MiRKAT library should be downloaded from:

#               http://research.fhcrc.org/wu/en/software.html
#------------------------------------------------------------------------------#

#' @title Weighted MiRKAT
#' @description  Runs MiRKAT function using the specified weights for each taxon.
#'
#'
#' @param x matrix with the composition for each individual
#' @param y response variable: dichotomous or continuous
#' @param w vector of weights for each taxa
#' @param cov covariates to have into account as a possible confounding factors
#'
#' @return the p-value associated to the compositional data with the corresponding
#' vector of weights.
#'
#' @examples
#'
#'# Matrix of samples
#'  x <- matrix(runif(100), nrow = 10)
#'# Response vector
#'  y <- rep(c(0,1),5)
#'# Vector of weights
#'  w <- c(2,rep(1,9))
#'
#'# Run wcomp.MiRKAT function
#'  wcomp.MiRKAT(x,y,w)
#' @export wcomp.MiRKAT


  wcomp.MiRKAT <- function(x,
                           y,
                           w = rep(1, ncol(x)),
                           cov = NULL){


#------------------------------------------------------------------------------#
# Build distance matrix
#------------------------------------------------------------------------------#

# Load library
  suppressMessages(library(MiRKAT))
  suppressMessages(library(proxy))

# Definte dist as dist{proxy}
#  dist <- proxy::dist

# wDM: weigthed Distance Matrix
  # Build the matrix
    wDM <- as.matrix(dist(x,method = w.dist, w=w ))

  #----------------------------------------#
  # P-value for weighted data
  #----------------------------------------#

  # Distance to Kernel
    K1<-D2K(wDM)
  # Extract the weighted p-value
    wp <- MiRKAT(y = as.numeric(y),
                 X = cov,
                 Ks = K1)
    return(wp)

  }

#------------------------------------------------------------------------------#





#------------------------------------------------------------------------------#
# FUNCTION: multiwcomp.MiRKAT
# INPUT:
#           x: matrix with compositional information
#           y: response variable: dichotomous or continuous
#       w.seq: sequence of weights that are going to be applied for each taxa
#         cov: data.fraem with covariates to have into account as a possible
#              confounding factors
# pseudocount: number of counts added to the information not to have zeros
#------------------------------------------------------------------------------#


# NOTE: MiRKAT library should be downloaded from:

#               http://research.fhcrc.org/wu/en/software.html
#------------------------------------------------------------------------------#



#' @title Taxa importance evaluation through weighted MiRKAT
#' @description  For different weight values, computes for each variable the p-value using
#' `wcomp.MiRKAT` function and implements a linear regression model giving back
#' its slope and the p-value for the slope.
#'
#'
#' @param x matrix with the composition for each individual
#' @param y response variable: dichotomous or continuous
#' @param w.seq vector of weights to consider in the analysis
#' @param cov covariates to hace into account as a possible confounding factors
#' @param zero.rep logical value indicating if cmultRepl() frunction from
#' zCompositions should be used in order to replace null values.
#' @param pseudocount a value to add to the whole matrix in order to avoid zeros.
#' Only if zero.rep == FALSE.
#' @return a data.frame with the slope and its p-value for the linear regression
#' model of wcomp.MiRKAT p-values considering the given weights.
#'
#' @examples
#'
#'# Matrix of samples
#'  x <- matrix(runif(100), nrow = 10)
#'# Response vector
#'  y <- rep(c(0,1),5)
#'
#'# Run multiwcomp.MiRKAT function
#'  multiwcomp.MiRKAT(x,y)
#' @export wcomp.MiRKAT


  multiwcomp.MiRKAT <- function(x, y, w.seq = c(.1, .4, .7, 1, 2, 3, 4, 5),
                                cov=NULL, zero.rep = T, pseudocount = 1){

#------------------------------------------------------------------------------#
# Auxiliar function
#------------------------------------------------------------------------------#

    weight.pval <- function(w){

      # Define a vector where p-values are going to be added
      pw <- vector()
      # Complete the vector
      for (i in 1:ncol(OTU)){
        # Build the weigths vector
        w.vec <- rep(1,ncol(OTU)); w.vec[i] <- w
        # Fill the table with the corresponding p.value (pseudocount made, so 0)
        pw <-  c(pw, wcomp.MiRKAT(OTU, y = res, w = w.vec, cov = cov))

      }
      return(pw)
    }


#------------------------------------------------------------------------------#
# Required packages
#------------------------------------------------------------------------------#

# List of required packages
  list.of.packages <- c("MiRKAT", "zCompositions")
# Not installed packages from the list.of.packages
  new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
# Install not installed ones
  if(length(new.packages)) install.packages(new.packages)

#------------------------------------------------------------------------------#
# Extract information
#------------------------------------------------------------------------------#

# Compositional information called OTU
  if (zero.rep){
    if (sum(x==0)>0) {
      library(zCompositions)
      OTU <- cmultRepl(x)
    } else { OTU <- x}
  } else {OTU <- x + pseudocount}
# Response variable (if it is a factor, the levels 0 and 1)
  if (class(y)=="factor"){
  res <- factor(y,labels=c(0,1))
  } else {res <- y}


#------------------------------------------------------------------------------#
# Tables to save the information
#------------------------------------------------------------------------------#

# Table to save the slopes and its p-values
  slope.tab <- matrix(0, nrow = ncol(OTU), ncol = 2)
# Add names to columns and rows
  row.names(slope.tab) <- colnames(OTU)
  colnames(slope.tab) <- c("Slope", "p-val Slope")


#------------------------------------------------------------------------------#
# Complete the information
#------------------------------------------------------------------------------#

# Build a parallelization scenario
  suppressMessages(library(foreach))
  suppressMessages(library(doParallel))
# Number of cores of the computer but one
  no_cores <- detectCores() - 2
# Register the number of cores
  registerDoParallel(no_cores)

# Variables of interest computed in parallel
  pw.tab <- foreach(h=w.seq,
                    .export=c("w.dist","w.geoM","wcomp.MiRKAT", "weight.pval",
                              "res","OTU","cov"),
                    .combine='cbind',
                    .multicombine=TRUE,
                    .init=vector()) %dopar% {
                      weight.pval(h)
                    }
# Stop the parallelization
  stopImplicitCluster()


# Table with the slopes for each variable and the p-values
  for (i in 1:nrow(slope.tab)){
    sum.FIT <- summary(lm(-log10(pw.tab[i,] + 1e-10) ~ w.seq))$coefficients
    slope.tab[i,] <-c(sum.FIT[2,1], sum.FIT[2,4])
  }



# Return slope.tab table
  return(list(slope.tab, pw.tab))

  }

#------------------------------------------------------------------------------#




#' @title Plot the results in wcomp.MiRKAT
#' @description Graphical representation of the slopes presented in the output
#' of wcomp.MiRKAT function.
#'
#'
#' @param tab the output of wcomp.MiRKAT function
#' @param col two colors for the positive and negatives slopes
#'
#' @return a barplot representing the slopes for each variable
#'
#' @examples
#'
#'# Matrix of samples
#'  x <- matrix(runif(100), nrow = 10)
#'# Response vector
#'  y <- rep(c(0,1),5)
#'
#'# Run multiwcomp.MiRKAT function
#'  A <- multiwcomp.MiRKAT(x,y)
#'# Represent the results
#'  wM.barplot(A[[1]])
#' @export wM.barplot


  wM.barplot <- function(tab, col= c("red", "blue")){

  # Load library ggplot2
    suppressMessages(library(ggplot2))
  # Tab as data.frame
    tab<-as.data.frame(tab)
  # Add a column with the name of each row
    tab$Name <- row.names(tab)
  # Add a column with the color for the bar
    tab$color <- ifelse(tab$Slope>0, "red", "blue")
  # Order tab
    tab <- tab[order(tab$Slope),]
  # Order the levels of Name according to the order presented in the table
    tab$Name <- factor(tab$Name, levels = tab$Name)


  # IMP.plot
    Slope.plot <- ggplot(tab, aes(x=factor(Name), y=Slope)) +
      geom_bar(stat="identity", aes(fill = color),
               size=1) +
      guides(size = FALSE) + # Not to show the legend of the size
      scale_fill_manual(name = "Group of . . .",
                        values = c(col[1], col[2]),
                        breaks = c(col[1], col[2]),
                        labels=c("NUM","DEN")) +

      ylab("Slope size") +
      xlab("") + theme_bw() +
      coord_flip() +
      ggtitle("Taxa contribution") +
      theme(strip.text.x = element_text(size=12, angle=0,
                                        face="bold",colour="white"),
            strip.text.y = element_text(size=12, face="bold"),
            strip.background = element_rect(colour="black",
                                            fill="black"),
            plot.title = element_text(size=20, vjust=2.25, hjust=0.5,
                                      face = "bold"),
            legend.position="none")

    return(Slope.plot)

  }

#------------------------------------------------------------------------------#
