.First.lib <- function(libname, pkgname, where) {
    if( !require("Biobase")) stop("can only load GOstats with Biobase")
    if( !require("graph")) stop("can only load GOstats with graph")
    if( !require("GO")) stop("can only load GOstats with GO")

    if(.Platform$OS.type == "windows" && require(Biobase) && interactive()
        && .Platform$GUI ==  "Rgui"){
        addVigs2WinMenu("GOstats")
    }


}
