##given a set of gene/probe identifiers obtain the GO graph that has all
##GO ids that those genes are annotated at, and
makeGOGraph <- function (x, what = "MF", lib = "hgu95av2",
                         removeRoot=TRUE)
{
    require(GO) || stop("no GO library")
    require(lib, character.only=TRUE, quietly=TRUE) ||
                     stop("data library not available")
    match.arg(what, c("MF", "BP", "CC"))
    wh <- paste("GO", what, "PARENTS", sep = "")
    dataenv <- get(wh, mode = "environment")
    GOpkgname<- paste(lib, "GO", sep="")
    GOpkgenv <- get(GOpkgname, mode="environment")
    newNodes <- mget(x, env = GOpkgenv, ifnotfound=NA)
    if( length(newNodes) == 1)
       bd = is.na(newNodes[[1]])
    else
       bd = is.na(newNodes)
    newNodes <- newNodes[!bd]

    newNodes <- lapply(newNodes, function(x) x[sapply(x, hasGOannote,
        what)])
    oldEdges <- vector("list", length = 0)
    oldNodes <- vector("character", length = 0)
    for (i in 1:length(newNodes)) {
        newN <- sapply(newNodes[[i]], function(x) x$GOID)
        done <- FALSE
        while (!done) {
            newN <- newN[!(newN %in% oldNodes)]
            if (length(newN) == 0)
                done <- TRUE
            else {
                oldNodes <- c(oldNodes, newN)
                numE <- length(newN)
                nedges <- mget(newN, env = dataenv, ifnotfound=NA)
                nedges <- nedges[!is.na(nedges)]
                oldEdges <- c(oldEdges, nedges)
                if (length(nedges) > 0)
                  newN <- sort(unique(unlist(nedges)))
                else newN <- NULL
            }
        }
    }
    rE <- vector("list", length = length(oldNodes))
    names(rE) <- oldNodes
    rE[names(oldEdges)] <- oldEdges
    rE <- lapply(rE, function(x) x <- which(oldNodes %in% x))
    names(rE) <- oldNodes

    rval = new("graphNEL", nodes = oldNodes, edgeL = lapply(rE,
        function(x) list(edges = x)), edgemode = "directed")

    if( removeRoot )
        rval = removeNode("GO:0003673", rval)
    rval
}


 ## helper function, determines if there is a go annotation for the
 ## desired mode
 hasGOannote <- function(x, which="MF") {
     if( !is.null(x$Ontology) ) {
         if( !is.na(x$Ontology) && x$Ontology == which )
            return(TRUE) else return(FALSE)
     }
     if( !is.character(x) )
         stop("wrong argument")
     tm <- get(x, env=GOTERM)
     if( names(tm)[1] == which)
         return(TRUE)
     return(FALSE)
 }


##start with one specific term and find it's more general terms
## if A has edges to B, C and D
##then at step 2, the nodes are A,B,C,D; new nodes B,C,D
##we need to find edges for each of them
##we're going to build a graphNEL
oneGOGraph <- function(x, dataenv) {
    require(GO) || stop("no GO library")
    if( length(x) > 1 )
        stop("wrong number of GO terms")
    if( !is.environment(dataenv) )
        stop("second argument must be an environment")

    oldEdges <- vector("list", length=0)
    oldNodes <- vector("character", length=0)
    newN <- x
    done <- FALSE
    while( !done ) {
        newN <- newN[!(newN %in% oldNodes)]
        if( length(newN) == 0 )
            done <- TRUE
        else {
            oldNodes <- c(oldNodes, newN)
            numE <- length(newN)
            nedges <- mget(newN, env=dataenv, ifnotfound=NA)
            nedges <- nedges[!is.na(nedges)]
            oldEdges <- c(oldEdges, nedges)
            if( length(nedges) > 0 )
                newN <- sort(unique(unlist(nedges)))
            else
                newN <- NULL
        }
    }
    rE <- vector("list", length=length(oldNodes))
    names(rE) <- oldNodes
    ##DON'T CHANGE: we must have the indices in the edgeL slot be
    ##   integer offsets  to the node list
    rE[names(oldEdges)] <- oldEdges
    rE <- lapply(rE, function(x) match(x, oldNodes))
    names(oldNodes) = oldNodes
    return(new("graphNEL", nodes=oldNodes, edgeL = lapply(rE,
               function(x) list(edges=x)),edgemode="directed" ))
}

##given two GO graphs, g1 and g2 join them together into a single new
##GO graph
combGOGraph <- function(g1, g2)
{
    if( !is(g1, "graphNEL") || !is(g2, "graphNEL") )
        stop("both arguments must be GO graphs")

    newn <- union(nodes(g1), nodes(g2))
    e1 <- edges(g1)
    e2 <- edges(g2)
    nnodes <- length(newn)
    rval <- vector("list", length=nnodes)
    names(rval) <- newn
    for( node in newn ) {
        el1 <- e1[[node]]
        el2 <- e2[[node]]
        rval[[node]] <- list(edges=match(union(el1, el2), newn))
    }
    return(new("graphNEL", nodes = newn, edgeL=rval, edgemode="directed"))
}

 ##GOleaves: the leaves of the GO graph are those nodes that have no
 ##inedges
  GOLeaves <- function(inG) {
      nG <- nodes(inG)
      iE <- inEdges(nG, inG)
      nG[sapply(iE, length) == 0]
  }

 ##distance between GO terms as described in B. Ding/R. Gentleman
 distDGGO <- function(term1, term2, dataenv) {
     g1 <- oneGOGraph(term1, dataenv)
     g2 <- oneGOGraph(term2, dataenv)
     n1 <- nodes(g1)
     n2 <- nodes(g2)
     ##drop the GO root -- if there
     m1 <- match("GO:0003673", n1)
     if( !is.na(m1) ) n1 <- n1[-m1]
     m2 <- match("GO:0003673", n2)
     if( !is.na(m2) ) n2 <- n2[-m2]
     length(intersect(n1,n2))/length(union(n1,n2))
 }

 ##the distance suggested by Cheng et al
 distCGO <- function(term1, term2, dataenv) {
     if( is(term1, "graph") )
         g1 <- term1
     else
         g1 <- oneGOGraph(term1, dataenv)
     if( is(term2, "graph") )
         g2 <- term2
     else
         g2 <- oneGOGraph(term2, dataenv)
     ig <- intersection(g1, g2)
     lfi <- GOleaves(ig)
     ##need a pathlength
 }

##three functions to get all the GO information for a set of GO terms
 getGOCategory <- function(x) {
     if( !is.character(x) )
         stop("need a character argument")
     if(length(x) == 0 )
         return( character(0))
     wh <- mget(x, env=GOTERM, ifnotfound=NA)
     return( sapply(wh, function(x) names(x)) )
 }

 getGOParents <- function(x) {
     if( !is.character(x) )
         stop("need a character argument")
     if(length(x) == 0 )
         return(list())
     hasMF <- mget(x, env=GOMFPARENTS, ifnotfound=NA)
     hasBP <- mget(x, env=GOBPPARENTS, ifnotfound=NA)
     hasCC <- mget(x, env=GOCCPARENTS, ifnotfound=NA)
     lenx <- length(x)
     rval <- vector("list", length=lenx)
     names(rval) <- x
     rval <- vector("list", length=lenx)
     names(rval) <- x
     for(i in 1:lenx) {
         if( (length(hasMF[[i]]) > 1 ) || !is.na(hasMF[[i]]) )
             rval[[i]] <- list(Ontology="MF", Parents=hasMF[[i]])
         else if( (length(hasMF[[i]]) > 1 ) || !is.na(hasBP[[i]]) )
             rval[[i]] <- list(Ontology="BP", Parents=hasBP[[i]])
         else if( (length(hasMF[[i]]) > 1 ) || !is.na(hasCC[[i]]) )
             rval[[i]] <- list(Ontology="CC", Parents=hasCC[[i]])
         else
             stop(paste(x[i], "is not a member of any ontology"))
     }
     return(rval)
 }

 getGOChildren <- function(x) {
     if( !is.character(x) )
         stop("need a character argument")
     if(length(x) == 0 )
         return(list())
     hasMF <- mget(x, env=GOMFCHILDREN, ifnotfound=NA)
     hasBP <- mget(x, env=GOBPCHILDREN, ifnotfound=NA)
     hasCC <- mget(x, env=GOCCCHILDREN, ifnotfound=NA)
     lenx <- length(x)
     rval <- vector("list", length=lenx)
     names(rval) <- x
     rval <- vector("list", length=lenx)
     names(rval) <- x
     for(i in 1:lenx) {
         if( (length(hasMF[[i]]) > 1 ) || !is.na(hasMF[[i]]) )
             rval[[i]] <- list(Ontology="MF", Children=hasMF[[i]])
         else if( (length(hasMF[[i]]) > 1 ) || !is.na(hasBP[[i]]) )
             rval[[i]] <- list(Ontology="BP", Children=hasBP[[i]])
         else if( (length(hasMF[[i]]) > 1 ) || !is.na(hasCC[[i]]) )
             rval[[i]] <- list(Ontology="CC", Children=hasCC[[i]])
         else
             rval[[i]] <- list()
     }
     return(rval)
 }

 getGOTerm <- function(x) {
     if( !is.character(x) )
         stop("need a character argument")
     if(length(x) == 0 )
         return(list())
     terms <- mget(x, env=GOTERM, ifnotfound=NA)
     ontology <- sapply(terms, function(x) names(x))
     return(split(terms, ontology))
 }
