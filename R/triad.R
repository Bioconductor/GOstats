 reduce2Degreek <- function(graph, k) {
    if(edgemode(graph) != "undirected")
       stop("only undirected graphs for now")
    done <- FALSE
    while(!done) {
        degs <- degree(graph)
        dgsub <- degs[degs> (k-1)]
        if( length(dgsub) == 0 )
           stop("no nodes left :-(")
        if( length(degs) == length(dgsub))
           done <- TRUE
        else
           graph <- subGraph(names(dgsub), graph)
   }
   return(graph)
  }

 isTriad <- function(x, y, z, elz, ely) {
   if( all(c(x,y) %in% elz ) && all(c(x,z) %in%  ely) )
     return(TRUE)
   return(FALSE)
 }

 enumPairs <- function(iVec) {
   leni <- length(iVec)
   if(leni < 2) return(vector(mode(iVec), length=0))
   eP <- vector("list", length=choose(leni, 2)/2)
   k<-1
   for( i in 1:(leni-1) ) {
     for(j in (i+1):leni) {
       eP[[k]] <- c(iVec[i], iVec[j])
       k <- k+1
     }
   }
   return(eP)
  }


 triadCensus <- function(graph) {
   g1 <- reduce2Degreek(graph, 2) ##all members have to have degree 2
   triads <- NULL
   k <- 1
   el <- edges(g1)
   for( n1 in nodes(g1) ) {
      for(j in enumPairs(el[[n1]])) {
        if( isTriad(n1, j[1], j[2], el[[j[2]]], el[[j[1]]]) ) {
            cand <- sort(c(n1,j))
            dupd <- sapply(triads, function(x) all(x==cand))
            if( length(triads)>1 && any(dupd) )
                next
            triads[[k]] <- cand
            k<- k+1
        }
      }
   }
   return(triads)
 }



