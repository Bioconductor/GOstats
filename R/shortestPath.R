##Copyright R. Gentleman 2004, all rights reserved
##some functions to do the shortest path analysis

shortestPath <- function(g, GOnode)
{
    if( !is(g, "graph") )
        stop("first argument must be a graph")

    if( edgemode(g) != "undirected" )
        stop("only undirected graphs for now, sorry")

    if( !is.character(GOnode) || length(GOnode) > 1 )
        stop("bad GO node description")

    require("GO", character.only=TRUE) ||
                      stop("GO library not available")

    ##obtain the LLIDs at the GO term

    LLs <- get(GOnode, GOGO2LL)

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
    rval <- NULL
    nNodes = length(whNodes)
    for(i in 1:nNodes) {
        vvX <- dijkstra.sp(g, whNodes[i])
        rval[[i]] <- vvX$distance[whNodes[-i]]
        if( verbose)
            print(paste("Processing node", i, "of", nNodes))
    }
    tfmatrix <- matrix(0, nr= nNodes, nc=nNodes)
    dimnames(tfmatrix) <- list(whNodes, whNodes)
    for(i in 1:nNodes){
        tfmatrix[i,names(rval[[i]])] <- rval[[i]]
    }
    return(tfmatrix)
}

compCorrGraph <- function(eSet, k=1, tau=0.6) {
    if( tau < 0 || tau > 1 ) stop("bad tau value")
    cB <- abs(cor(t(exprs(eSet))))
    cB[cB < tau] <- 0
    whE = cB[cB>0]
    cB[whE] <- ((1-cB)[whE])^k
    diag(cB) = 0
    ##FIXME: cB=1?
    require("SparseM", quietly=TRUE) || stop("need SparseM package")
    v1<-as.matrix.csr(cB, nr=dim(cB)[1], nc=dim(cB)[2])
    rv <- sparseM2Graph(v1, geneNames(eSet))
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

