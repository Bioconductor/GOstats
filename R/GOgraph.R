##given a set of LOCUSLINK IDs obtain the GO graph that has all
##GO ids that those genes are annotated at, and
makeGOGraph <- function (x, Ontology = "MF", removeRoot=TRUE)
{
    require(GO) || stop("no GO library")
    match.arg(Ontology, c("MF", "BP", "CC"))
    wh <- paste("GO", Ontology, "PARENTS", sep = "")
    dataenv = get(wh, mode="environment")
    newNodes <- mget(x, env = GOLOCUSID2GO, ifnotfound=NA)
    if( length(newNodes) == 1)
       bd = is.na(newNodes[[1]])
    else
       bd = is.na(newNodes)
    newNodes <- newNodes[!bd]

    newNodes <- lapply(newNodes, function(x) x[sapply(x, function(x)
                                                      x$Ontology == Ontology)])
    oldEdges <- vector("list", length = 0)
    oldNodes <- vector("character", length = 0)
    for (i in 1:length(newNodes)) {
        newN <- unique(sapply(newNodes[[i]], function(x) x$GOID))
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


 ## helper function, determines if there is a GO annotation for the
 ## desired mode
 hasGOannote <- function(x, which="MF") {
     if( is(x, "GOTerms") ) {
         cat = Category(x)
         if( !is.na(cat) && cat == which )
            return(TRUE) else return(FALSE)
     }
     if( is.list(x) ) {
         gT = sapply(x, function(y) is(y, "GOTerms"))
         if( any(gT) ) {
             if( all(gT) ) {
                 cats = sapply(x, Category)
                 return(cats == which)
             }
             else
                 stop("mixed arguments not allowed")
         }
     }
     if( !is.character(x) )
         stop("wrong argument")
     tm <- getGOCategory(x)
     return(tm == which)
 }


##start with one specific term and find it's more general terms
## if A has edges to B, C and D
##then at step 2, the nodes are A,B,C,D; new nodes B,C,D
##we need to find edges for each of them
##we're going to build a graphNEL
GOGraph = function(x, dataenv) {
    require(GO) || stop("no GO library")
    if (!is.environment(dataenv))
        stop("second argument must be an environment")
    ##this is the old oneGOGraph code - but it just works for
    ##multiple inputs
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

##to be deprecated
oneGOGraph <- function(x, dataenv) {
    require(GO) || stop("no GO library")
    if( length(x) > 1 )
        stop("wrong number of GO terms")
    if (!is.environment(dataenv))
        stop("second argument must be an environment")
    GOGraph(x, dataenv)
}
##given two GO graphs, g1 and g2 join them together into a single new
##GO graph
combGOGraph <- function(g1, g2)
{
    .Deprecated("join", "graph")
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

 ##similarity functions - based on graphs
 simUI = function(g1, g2) {
     if(!is(g1, "graph") || !is(g2, "graph") )
         stop("only works for graphs")
     n1 = nodes(g1); n2 = nodes(g2)
    length(intersect(n1, n2))/length(union(n1, n2))
 }

 simLP = function(g1, g2) {
    if(!is(g1, "graph") || !is(g2, "graph") )
         stop("only works for graphs")
    require("RBGL") || stop("need RBGL for this similarity")
    ig <- intersection(g1, g2)
    lfi <- GOLeaves(ig)
    degs = degree(ig)
    root = names(degs$outDegree)[degs$outDegree == 0]
    paths = sp.between(ig, lfi, root)
    plens = sapply(paths, function(x) x$length)
    max(plens)
}

 ##a helper function to get the right GOIDs
 .getWHEC = function(llid, wh, eCodes) {
     x = GOLOCUSID2GO[[llid]]
     if( !is.null(eCodes) )
         x = dropECode(x, eCodes)
     unique(unlist(getOntology(x, wh)))
  }

 simLL = function(ll1, ll2, Ontology="MF", measure = "LP",
                      dropCodes=NULL) {
    wh = match.arg(Ontology, c("MF", "BP", "CC"))
    ll1GO = .getWHEC(ll1, wh, dropCodes)
    ll2GO = .getWHEC(ll2, wh, dropCodes)
    dataenv = get(paste("GO", wh, "PARENTS", sep=""),
                    mode="environment")
    g1 = GOGraph(ll1GO, dataenv)
    g2 = GOGraph(ll2GO, dataenv)
    sm = match.arg(measure, c("LP", "UI"))
    sim = switch(sm,
           LP = simLP(g1, g2),
           UI = simUI(g1, g2))
    return(list(sim=sim, measure=measure, g1 = g1, g2 =g2))
  }


##three functions to get all the GO information for a set of GO terms
##FIXME: these need to be renovated - probably removed even..
 getGOCategory <- function(x) {
     if( !is.character(x) )
         stop("need a character argument")
     if(length(x) == 0 )
         return( character(0))
     wh <- mget(x, env=GOTERM, ifnotfound=NA)
     return( sapply(wh, Category) )
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
     ontology <- sapply(terms, Category)
     terms = sapply(terms, Ontology)
     return(split(terms, ontology))
 }
