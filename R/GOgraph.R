 combGOGraph = function(g1, g2) .Defunct("join", package = "GOstats")

##given a set of LOCUSLINK IDs obtain the GO graph that has all
##GO ids that those genes are annotated at, and
makeGOGraph <- function (x, Ontology = "MF", removeRoot=TRUE)
{
    match.arg(Ontology, c("MF", "BP", "CC"))
    wh <- paste("GO", Ontology, "PARENTS", sep = "")
    dataenv = get(wh, mode="environment")
    newNodes <- mget(x, env = GOENTREZID2GO, ifnotfound=NA)
    if( length(newNodes) == 1)
       bd = is.na(newNodes[[1]])
    else
       bd = is.na(newNodes)
    newNodes <- newNodes[!bd]

    newNodes <- lapply(newNodes, function(x) x[sapply(x,
                         function(x) {if (is.na(x$Ontology) )
                                          return(FALSE)
                                      else
                                          x$Ontology == Ontology})])
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

    if (removeRoot) {
        ## For compatibility with old GO packages, check for "GO:0003673" if
        ## "all" is not found.
        if (is.element("all", nodes(rval)))
          rval = removeNode("all", rval)
        else if (is.element("GO:0003673", nodes(rval)))
          rval = removeNode("GO:0003673", rval)
    }
    rval
}


##start with one specific term and find it's more general terms
## if A has edges to B, C and D
##then at step 2, the nodes are A,B,C,D; new nodes B,C,D
##we need to find edges for each of them
##we're going to build a graphNEL
GOGraph = function(x, dataenv) {
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
    if( length(x) > 1 )
        stop("wrong number of GO terms")
    if (!is.environment(dataenv))
        stop("second argument must be an environment")
    GOGraph(x, dataenv)
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
    ig <- intersection(g1, g2)
    if( numNodes(ig) < 2 ) return(0)
    lfi <- GOLeaves(ig)
    degs = degree(ig)
    root = names(degs$outDegree)[degs$outDegree == 0]
    paths = sp.between(ig, lfi, root)
    plens = sapply(paths, function(x) x$length)
    max(plens)
}

 ##a helper function to get the right GOIDs
 .getWHEC = function(llid, wh, eCodes) {
     x = GOENTREZID2GO[[llid]]
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
    if( length(nodes(g1)) == 0 || length(nodes(g2)) == 0 )
      return(NA)
    sm = match.arg(measure, c("LP", "UI"))
    sim = switch(sm,
           LP = simLP(g1, g2),
           UI = simUI(g1, g2))
    return(list(sim=sim, measure=measure, g1 = g1, g2 =g2))
  }


