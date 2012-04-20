require("GOstats") || stop("unable to load GOstats")
BiocGenerics:::testPackage("GOstats", "UnitTests", ".*_test\\.R$")
