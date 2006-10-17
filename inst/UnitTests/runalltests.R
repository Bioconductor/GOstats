require("RUnit", quietly=TRUE) || stop("RUnit package not found")
require("GOstats")


TEST_DATA_DIR <- "data"
runitPat <- ".*_test\\.R$"
runitDirs <- c(".")
suite <- defineTestSuite(name="GOstats Test Suite",
                         dirs=runitDirs,
                         testFileRegexp=runitPat,
                         rngKind="default",
                         rngNormalKind="default")
result <- runTestSuite(suite)
printTextProtocol(result, showDetails=FALSE)


