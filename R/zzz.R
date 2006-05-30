.onAttach <- function(libname, pkgname) {
    if (.Platform$OS.type == "windows" && require("Biobase")
        && interactive() && .Platform$GUI == "Rgui") {
        addVigs2WinMenu("GOstats")
    }
}
