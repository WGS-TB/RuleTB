source("Interpretable.R")

applyRules = function(filename, broad = TRUE, nameCol = 1, numLoci = 24) {
  if (broad) {
    load("BroadRules.RData")
  } else {
    load("NarrowRules.RData")
  }
  Database = prepareDatabase(filename, relCol = nameCol, numLoci = numLoci, unique = FALSE)
  testData = Database[[1]]
  correctLoci = match(Loci, substr(colnames(testData), 2, nchar(colnames(testData))))
  if (any(is.na(correctLoci))) {
    print(paste("Warning:", sum(is.na(correctLoci)), "loci are missing! This may reduce prediction accuracy!"))
  }
  testData = testData[,correctLoci]
  predictions = applyAllRules(testData, Rules, NULL, hard = FALSE)
  rownames(predictions) = Database[[2]]
  predictions
}