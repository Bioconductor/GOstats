##########################################
##
##  Copyright 2005 James W. MacDonald 
##
##  A function to output a table from running
##  GOHyperG() on a set of Affymetrix probe IDs
##
##
##  Modified 6-9-05 to output the table
##  Added hyperG2Affy 10-10-2005
##   - This function takes GOIDs form hyperGtable and outputs
##     a list giving the Affy IDs associated with each GOID
##
##
##  Moved hyperGtable and hyperG2Affy to GOstats 2-23-06
##
##  Added universe argument 3-22-06
###########################################

hyperGtable <- function(probids, lib, type="MF", pvalue=0.05,
                        min.count=10, save = FALSE, output = TRUE,
                        filename = NULL, universe = NULL){
  require(lib, quietly = TRUE, character.only = TRUE) || stop(paste("The ", lib, " package is required"))
  lls <- getLL(probids, lib)
  lls <- unique(lls)
  lls <- lls[!is.na(lls)]
  tmp <- GOHyperG(lls, lib, type, universe)
  index <- tmp$pvalues < pvalue & tmp$goCounts > min.count
  if(sum(index) == 0)
    stop(paste("There are no significant GO terms using a p-value of ",pvalue,
               " and a minimum count of ",  min.count,".\n", sep = ""), call. = FALSE)
  wh <- mget(names(tmp$pvalues[index]), GOTERM)
  tmp.terms <- sapply(wh, Term)
  out <- data.frame(names(tmp$pvalues[index]), tmp.terms, round(tmp$pvalues[index], 3), 
                       tmp$intCounts[index], tmp$goCounts[index])
  names(out) <- c("GO Term","Description","p-values",
                     "Number significant","Number on chip")
  if(output){
    if(!is.null(filename)){
      write.table(out, paste(filename, ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
    }else{
      stop("You need to give a filename in order to output this table!\n")
    }
  }
  if(save)
    out
}

hyperG2Affy <- function(probids, lib, type="MF", pvalue=0.05,
                        min.count=10, universe = NULL){
  require(lib, quietly = TRUE, character.only = TRUE) || stop(paste("The ", lib, " package is required"))
  lls <- getLL(probids, lib)
  lls <- unique(lls)
  lls <- lls[!is.na(lls)]
  tmp <- GOHyperG(lls, lib, type, universe)
  index <- tmp$pvalues < pvalue & tmp$goCounts > min.count
  if(sum(index) == 0)
    stop(paste("There are no significant GO terms using a p-value of ",pvalue,
               " and a minimum count of ",  min.count,".\n", sep = ""), call. = FALSE)
  wh <- mget(names(tmp$pvalues[index]), GOTERM)
  tmp.terms <- sapply(wh, Term)
  index2 <- match(names(tmp.terms), names(tmp$go2Affy))
  out <- lapply(tmp$go2Affy[index2], function(x) x[x %in% probids])
  names(out) <- names(tmp.terms)
  out
}
 
