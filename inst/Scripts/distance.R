##Copyright B. Ding, 2004, all rights reserved

##FIXME: I think some of Jianhua's new environments remove all need
##for these functions
## function to get all the parent nodes for "go" (including "go" itself)
##
##NAME of Beiying's function:
##go.all.parents <- function(go)

## function to get all the children nodes for "go" (including itself)
##go.all.children <- function(go)

## function to find all genes annotated at GO including its children nodes
## env1: mapping from LL to GO
## env2: mapping from GO to LL

##go2allgenes <- function(gos,env1=NULL,env2=NULL)
## calculate the distance between two GO nodes by the cardinality of their GO graph intersection divided by the cardinality of their graph union.
##go.dist <- function(go1,go2)


## calculate distance between GO graphs for geneslist1 and geneslist2
## using all GO nodes the genes can be mapped to
## env3: mapping from GO to ALL LL

geneslist.dist.allnode.ll <- function(geneslist1,geneslist2,ontology="BP",
env1=NULL, env2=NULL, env3=NULL)
  {
    geneslist1.graph <- NULL

    for (j in geneslist1)
      {
        geneslist1.gos <- unlist(mget(j,get(env1)))
        geneslist1.ontology <- geneslist1.gos[!is.na(mget(geneslist1.gos,get(paste("GO",ontology,"ID2TERM",sep=""))))]

    #  print(geneslist.ontology)
        if(length(geneslist1.ontology)==0)
          geneslist1.ontology <- NA

        geneslist1.graph <- union(geneslist1.graph,names(go.all.parents(geneslist1.ontology)))
      }

    geneslist2.graph <- NULL

    for (j in geneslist2)
      {
        geneslist2.gos <- unlist(mget(j,get(env1)))
        geneslist2.ontology <- geneslist2.gos[!is.na(mget(geneslist2.gos,get(paste("GO",ontology,"ID2TERM",sep=""))))]

    #  print(geneslist.ontology)
        if(length(geneslist2.ontology)==0)
          geneslist2.ontology <- NA

        geneslist2.graph <- union(geneslist2.graph,names(go.all.parents(geneslist2.ontology)))
      }

    geneslist1.geneslist2.intersect.graph <- intersect(geneslist1.graph,geneslist2.graph)
    geneslist1.geneslist2.intersect.graph <- geneslist1.geneslist2.intersect.graph[!(geneslist1.geneslist2.intersect.graph %in% c("GO:0008150","GO:0003674","GO:0005575","GO:0003673"))]
    #    print(gene.geneslist.intersect.graph)

    if(length(geneslist1.geneslist2.intersect.graph)>0)
      {

        geneslist1.geneslist2.intersect.graph.lls <- list()

        for (j in geneslist1.geneslist2.intersect.graph)
          {
              geneslist1.geneslist2.intersect.graph.lls[[j]] <-
              intersect(unlist(mget(j, get(env3), ifnotfound=NA)), geneslist2)
          }

        geneslist1.geneslist2.intersect.graph.lls <-
        geneslist1.geneslist2.intersect.graph.lls[sapply(
        geneslist1.geneslist2.intersect.graph.lls, function(x) length(x)>0)]

        score <- sum(sapply(geneslist1.geneslist2.intersect.graph.lls,
        function(x)
        length(x))/sapply(mget(names(geneslist1.geneslist2.intersect.graph.lls),
        get(env3), ifnotfound=NA),function(x)
                          length(x)))
    }else{
        score <- 0
    }
    return(score)
}

## significance test for gene.dist.allnode.ll where geneslist1 is the gene set of interest and geneslist2 is the candidate gene list
geneslist.dist.allnode.ll.sig <- function(geneslist1, geneslist2,
                                          geneslist2.pool=NULL,
                                          ontology="BP", n.samp=100,
                                          env1=NULL, env2=NULL,
                                          env3=NULL)
{
    score <- geneslist.dist.allnode.ll(geneslist1, geneslist2,
                                       ontology=ontology, env1, env2, env3)
                                        # print(score)
    score.null <- NULL

    for ( i in 1:n.samp) {
        print(i)
        geneslist2.null <- sample(geneslist2.pool, length(geneslist2),
                                  FALSE)
                                        # print(geneslist2.null)
        score.null[i] <- geneslist.dist.allnode.ll(geneslist1,
                      geneslist2.null, ontology=ontology, env1, env2, env3)
    }

    pval <- mean(score<=score.null)
    return(list(score=score, pval=pval))
  }


