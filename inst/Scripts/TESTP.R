##FIRST LET'S SEE WHAT THE P-VALUE PROBLEM LOOKS LIKE

##FIXME: just ran the first simulation set1 - need to do set2 again
##and start looking at the appropriate set of adjustments
 library("GOstats")
 library("hgu95av2.db")

 w1 <- as.list(hgu95av2LOCUSID)

 w2 <- unique(unlist(w1))

 set.seed(123)

 myLL <- sample(w2, 100)

 xx <- hyperGTest(myLL)

 xx$numLL

 xx$numInt
[1] 60

 sum(xx$pvalues < 0.01)
[1] 9

 rval = NULL
 for(i in 1:10) {
     mx = sample(w2, 100)
     xx = hyperGTest(mx)
     rval[[i]] = list(nI = xx$numInt,gC= xx$goCounts,pV= xx$pvalues)
 }

 nI = sapply(rval, function(x) x$nI)
 gC = sapply(rval, function(x) length(x$gC))

 countpvs = function(il, pv=0.05)
    sapply(il, function(x) sum(x$pV < pv))

 nPV = countpvs(rval)

 sapply(rval, function(x) max(x$pV))


##analysis

 analyzeRun = function(x, pv=0.05) {
     nI = sapply(x, function(y) y$nI)
     gC = sapply(x, function(y) length(y$gC))
     nPV = countpvs(x, pv)
     return(list(nI=nI, gC=gC, pvs = nPV))
 }

##for set1
set1 = genPVs(100,100)

s1a = analyzeRun(set1)
mean(s1a$pvs/s1a$gC) #should be about 0.05, but it isn't

s2a = analyzeRun(set1, 0.01)
mean(s2a$pvs/s2a$gC)  ##is pretty close to 0.01

s3a = analyzeRun(set1, 0.1)
mean(s3a$pvs/s3a$gC)  ##about 3 times what it should be

#for set2
set2 = genPVs(100,100)

s2.05 = analyzeRun(set2)
mean(s2.05$pvs/s2.05$gC) #should be about 0.05, but it isn't

s2.01 = analyzeRun(set2, 0.01)
mean(s2.01$pvs/s2.01$gC)  ##is pretty close to 0.01

s2.1 = analyzeRun(set2, 0.1)
mean(s2.1$pvs/s2.1$gC)  ##about 3 times what it should be

##look at int g size

set3 = genPVs(500, 75)


##run a simulation

 genPVs = function(nsim, nsamp, verbose=TRUE) {
      rval = NULL
      for(i in 1:nsim) {
          mx = sample(w2, nsamp)
          rval[[i]] = hyperGTest(mx)
          if( verbose )
              print(paste("finished", i, "of", nsim))
      }
      return(rval)
  }

 set.seed(455)

 set1 = genPVs(100, 100)


#################
## compare dists
#################

 data(Ndists)
 data(Bdists)

 bothBad = (Ndists==Inf & Bdists==Inf)

 badN = (Ndists==Inf & Bdists !=Inf)

 badB = (Bdists==Inf & Ndists !=Inf)

