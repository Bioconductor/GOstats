.First.lib <- function(libname, pkgname, where) {
    ## require's version specificity doesn't allow for patterns.
    ## These will simulate it, in a hackish manner
    installed <- installed.packages()

    ## Make sure the metadata is at least 1.6.x
    ## Metadata is only 'Suggests' tho, so ok if they
    ## don't have it.
    cur <- match("hgu95av2", installed[,"Package"])
    if (!is.na(cur)) {
        if (compareVersion(installed[cur,"Version"], "1.6.0") < 0)
            stop("Your version of hgu95av2 must be at least 1.6.x")
        else
            require("hgu95av2")
    }

    ## The other packages are dependencies, must have them.  Make
    ## sure the BioC packages are the devel versions (as of 30/06/04)

    cur <- match("GO", installed[,"Package"])
    if ((!is.na(cur))&&(compareVersion(installed[cur,"Version"], "1.6.0")
                      < 0))
        stop("You must have package GO version >= 1.6.0")
    else
        require("GO")

    cur <- match("graph", installed[,"Package"])
    if ((!is.na(cur))&&(compareVersion(installed[cur,"Version"], "1.4.1")
                      < 0))
        stop("You must have package graph version >= 1.4.1")
    else
        require("graph")

    cur <- match("annotate", installed[,"Package"])
    if ((!is.na(cur))&&(compareVersion(installed[cur,"Version"], "1.4.0")
                      < 0))
        stop("You must have package annotate version >= 1.4.0")
    else
        require("annotate")

    cur <- match("Biobase", installed[,"Package"])
    if ((!is.na(cur))&&(compareVersion(installed[cur,"Version"], "1.4.15")
                      < 0))
        stop("You must have package Biobase version >= 1.4.15")
    else
        require("Biobase")

    cur <- match("RBGL", installed[,"Package"])
    if ((!is.na(cur))&&(compareVersion(installed[cur,"Version"], "1.2.2")
                      < 0))
        stop("You must have package RBGL version >= 1.2.2")
    else
        require("RBGL")

    cur <- match("genefilter", installed[,"Package"])
    if ((!is.na(cur))&&(compareVersion(installed[cur,"Version"], "1.4.0")
                      < 0))
        stop("You must have package genefilter version >= 1.4.0")
    else
        require("genefilter")

    cur <- match("multtest", installed[,"Package"])
    if ((!is.na(cur))&&(compareVersion(installed[cur,"Version"], "1.4.1")
                      < 0))
        stop("You must have package multtest version >= 1.4.1")
    else
        require("multtest")


    if(.Platform$OS.type == "windows" && require(Biobase) && interactive()
        && .Platform$GUI ==  "Rgui"){
        addVigs2WinMenu("GOstats")
    }


}
