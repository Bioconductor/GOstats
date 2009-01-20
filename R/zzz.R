.onAttach <- function(libname, pkgname) {
    if (.Platform$OS.type == "windows" && require("Biobase")
        && interactive() && .Platform$GUI == "Rgui") {
        addVigs2WinMenu("GOstats")
    }
}

GOenv <- function(what) {
    getAnnMap(what, "GO", load=TRUE,
              type=c("db", "env"))
}

.get_eg_to_go_fun <- function(mapfun, chip, reverse=FALSE) {
    if (!is.null(chip) && is.character(chip)) {
        if (!is.null(mapfun) && is.function(mapfun))
          warning("ignoring 'chip' argument in favor of 'mapfun'")
        else {
            ## create mapfun from chip
            mapfun <- tryCatch({
                eg2go <- make_eg_to_go_map(chip)
                if (reverse)
                  eg2go <- revmap(eg2go)
                function(x) mget(x, eg2go, ifnotfound=NA)
            }, error=function(e) {
                ## if only we had classed exceptions!  and we can't
                ## reliable grep the condition message since it may be
                ## localized.
                msg <- paste(conditionMessage(e),
                             "\nDB-based version of ", chip, " not found.",
                             "\nReverting to use of environment-based GO")
                warning(msg, call.=FALSE)
                NULL
            })
        }
    }

	if (!is.function(mapfun)) {
		stop("either mapfun or chip must be specified")
	}

	mapfun
}

