##Copyright R. Gentleman 2004, all rights reserved
##some functions to do the shortest path analysis

shortestPath <- function(g, GOnode, mapfun=NULL, chip=NULL)
{
    if( !is(g, "graph") )
        stop("first argument must be a graph")

    if( edgemode(g) != "undirected" )
        stop("only undirected graphs for now, sorry")

    if( !is.character(GOnode) || length(GOnode) > 1 )
        stop("bad GO node description")

    ##obtain the LLIDs at the GO term
    go2egfun <- .get_eg_to_go_fun(mapfun, chip, reverse=TRUE)
    LLs <- unique(go2egfun(GOnode)[[1]])

    m1 <- match(LLs, nodes(g))
    notthere <- LLs[is.na(m1)]
    there <- LLs[!is.na(m1)]
    if( length(there) <= 1 )
        stop("no nodes correspond to the specified term")

    nthere <- length(there)
    rval <- vector("list", (nthere*(nthere-1))/2)
    nms <- rep("", (nthere*(nthere-1))/2)
    k <- 1
    there <- as.character(there)
    for(i in 1:(nthere-1))
        for(j in (i+1):nthere) {
            nms[k] <- paste(there[i], there[j], sep=":")
            rval[[k]] <- sp.between(g,there[i], there[j])
            k <- k+1
        }
    names(rval) <- nms
    return(list(shortestpaths=rval, nodesUsed=there,
                nodesNotUsed=notthere))
}

compGdist <- function(g, whNodes, verbose=FALSE) {
    nNodes = length(whNodes)
    if(nNodes <= 1) 
        stop("can only compute distances between two or more nodes")
    if( any( !(whNodes %in% nodes(g))) )
        stop("specified nodes are not in the graph")
    tfmatrix <- matrix(0, nrow= nNodes, ncol=nNodes)
    for(i in 1:nNodes) {
        vvX <- dijkstra.sp(g, whNodes[i])
        tfmatrix[i,] <- vvX$distance[whNodes]
        if( verbose)
            print(paste("Processing node", i, "of", nNodes))
    }
    dimnames(tfmatrix) <- list(whNodes, whNodes)
    tfmatrix
}

compCorrGraph <- function(eSet, k=1, tau=0.6) {
    if( tau < 0 || tau > 1 ) stop("bad tau value")
    cB <- abs(cor(t(exprs(eSet))))
    cB[cB < tau] <- 0
    whE = cB>0
    ##FIXME: any two genes whose exprssion values are perfectly
    ##correlated will get dropped - seems pretty unlikely
    cB[whE] <- ((1-cB)[whE])^k
    ##set diag to zero, since we do not want self-edges
    diag(cB) = 0
    requireNamespace("SparseM", quietly=TRUE) || stop("need SparseM package")
    v1<-SparseM::as.matrix.csr(cB, nrow=dim(cB)[1], ncol=dim(cB)[2])
    sparseM2Graph(v1, featureNames(eSet))
}

 ##given a matrix and a logical vector return two lists
 ##the row names and the column names for each selected entry

 idx2dimnames = function(x, idx) {
     if(length(idx) != length(x) )
         stop("x and idx different lengths")
     rowV = row(x)[idx]
     colV = col(x)[idx]
     return(list(rowNames = dimnames(x)[[1]][rowV],
                 colNames = dimnames(x)[[2]][colV]))
 }

##take as input a distance matrix
##if the number of infs is 1 less than the number of rows
##then that nod is not connected
 notConn = function(dists) {
   allInf = (dists==Inf)
   dns = idx2dimnames(dists, allInf)
   tab1 = table(dns[[1]])
   whAll = tab1[tab1 == (nrow(dists)-1)]
   if( any(tab1 != (nrow(dists)-1) & tab1 != length(whAll)) )
      warning("inf dists seem to be wrong")
   return(names(whAll))
}

