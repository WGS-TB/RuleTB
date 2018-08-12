source("Utilities.R")
source("Alternative.R")

### This function sequentially extracts the most successful rules for predicting the lineages of
### a given set of strains. If exact =  TRUE, the optimization uses integer linear programming.
sequentialDriver = function(strains, lineages, exact, baseFile = NULL) {
  redLines = gsub(" ", "", lineages)
  redLines = gsub("/", "", redLines)
  uLineages = unique(redLines)
  L = length(uLineages)
  curStrains = strains
  curLineages = redLines
  allRules = vector("list", L)
  for (ind in 1:(L - 1)) {
    print(paste("Processing iteration", ind))
    curBaseFile = gsub(".csv", paste0(ind, ".csv"), baseFile)
    curRules = interpretDriver(curStrains, curLineages, exact = exact, baseFile = curBaseFile)
    curBest = which.min(sapply(curRules, function(x) {x[[2]]}))
    curBestLine = names(curBest)
    if (ind < (L - 1)) {
      allRules[[ind]] = curRules[[curBestLine]]
      names(allRules)[ind] = curBestLine
      curStrains = curStrains[curLineages != curBestLine, , drop = FALSE]
      curLineages = curLineages[curLineages != curBestLine]
    }
    else {
      if (curBest == 2) {
        curRules = curRules[c(2,1)]
      }
      allRules[ind:L] = curRules
      names(allRules)[ind:L] = names(curRules)
    }
  }
  allRules
}

### This function computes the confusion matrix of predictions and correct labels (names must match)
getConfusion = function (labels, predictions) {
  cNames = sort(unique(c(labels, predictions)))
  L = length(cNames)
  confusion = matrix(0, L, L, dimnames = list(cNames, cNames))
  Tab = table(labels, predictions)
  confusion[rownames(Tab), colnames(Tab)] = Tab
  # commonNames = intersect(rownames(confusion), colnames(confusion))
  # stopifnot(all(length(commonNames) == dim(confusion)))
  Diag = confusion[cbind(cNames, cNames)]
  accuracy = sum(Diag)/sum(confusion)
  cSum = colSums(confusion[,cNames])
  rSum = rowSums(confusion[cNames,])
  confusion = cbind(rbind(confusion, Specificity = Diag/cSum), Sensitivity = c(Diag/rSum, accuracy))
  confusion
}

### This function finds the rule best explaining the features defining each lineage in given strains  
interpretDriver = function(strains, lineages, exact = TRUE, baseFile = NULL) {
  redLines = gsub(" ", "", lineages)
  redLines = gsub("/", "", redLines)
	uLineages = unique(redLines)
	allRules = vector("list", length(uLineages))
	names(allRules) = uLineages
	for (line in uLineages) {
	  curClasses = (redLines == line)
	  curFile = gsub(".csv", paste0(line, ".csv"), baseFile)
	  print(paste("Processing lineage", line))
	  allRules[[line]] = findRules(strains, curClasses, exact = exact, baseFile = curFile)
	}
	allRules
}

### This auxiliary function creates concatenated ranges specified by given start and stop positions 
makeRanges = function(starts, stops) {
  L = length(starts)
  stopifnot(length(stops) == L)
	allLens = stops - starts
	result = rep(NA, sum(allLens))
	curPos = 1
	for (ind in 1:L) {
		curLen = allLens[ind]
		curRange = 0:(curLen - 1)
		result[curPos + curRange] = starts[ind] + curRange
		curPos = curPos + curLen
	}
	result
}

### This function is a driver for finding the rules predicting the classes of a given input matrix
### If exact = TRUE, integer linear programming is used to solve the optimization problem
### lambda specifies the misclassification (slack) penalty, and baseFile is used for the LP files
findRules = function(inputMatrix, classes, exact, baseFile, lambda = 100) {
	Input = createRuleILP(inputMatrix, classes = classes, exact = exact, lambda = lambda)
	Res = runAndParse(baseFile, Input, baseID = 0, prefix1 = "Y", prefix2 = "X", mat = NULL)	
	Res
}

### This function creates the (integer) linear program formulation for predicting the input classes
### based on Malioutov and Varshney's paper, "Exact Rule Learning via Boolean Compressed Sensing"
createRuleILP = function(inputMatrix, classes, exact, lambda) {
	numSamples = length(classes)
	numPos = sum(classes == TRUE)
	numNeg = sum(classes == FALSE)
	ranges = apply(inputMatrix, 2, range)
	n = ncol(inputMatrix)
	allVars = ranges[2,] - ranges[1,]
	Ind = rbind(rep(1:n, allVars), makeRanges(ranges[1,], ranges[2,]))
	numVars = sum(allVars)
	numConstraints = numSamples + 2 * numVars
	if (!exact) {
		numConstraints = numConstraints + numNeg
	}
	Lines = vector("character", numConstraints + 5)
	Lines[1] = "Minimize"
	obj1 = paste0("X", Ind[1,], "U", Ind[2,], collapse = " + ")
	obj2 = paste0("X", Ind[1,], "L", Ind[2,], collapse = " + ")
	obj  = paste0(obj1, " + ", obj2)
	if (!exact) {
		obj = paste0(obj, " + ", paste0(lambda, " Y", 1:numSamples, collapse = " + "))
	}
	Lines[2] = obj
	Lines[3] = "Subject to"
	for (ind in 1:numSamples) {
		curInput = inputMatrix[ind,]
		curDelta = curInput - ranges[1,]
		curDir = lapply(1:n, function(x) {c(rep("U", curDelta[x]), rep("L", allVars[x] - curDelta[x]))})
		curLine = paste0("X", Ind[1,], unlist(curDir), Ind[2,], collapse = " + ")
		extra = rep(paste0("Y", ind), !exact)
		curStatus = classes[ind]
		if (curStatus) {
			Lines[ind + 3] = paste0(curLine, rep(" - ", !exact), extra, " = 0")
		}
		else {
			Lines[ind + 3] = paste0(curLine, rep(" + ", !exact), extra, " >= 1")
		}
	}
	Lines[4 + numSamples] = "Binary"
	Lines[4 + numSamples + (1:numVars)] = paste0("X", Ind[1,], "U", Ind[2,])
	Lines[4 + numSamples + numVars + (1:numVars)] = paste0("X", Ind[1,], "L", Ind[2,])
	if (!exact) {
		Lines[4 + numSamples + 2 * numVars + (1:numNeg)] = paste0("Y", which(!classes))
	}
	Lines[length(Lines)] = "End"
	Lines
}

### This function applies a collection of rules to a given set of strains, in the specified order.
### If correct lineages are specified, a confusion matrix for the predictions is also computed.
### If hard = TRUE, the fractions are converted to a TRUE or FALSE format for each of the rules.
applyAllRules = function(strains, rules, lineages = NULL, hard = FALSE) {
  ruleStrains = names(rules)
  redRules = lapply(rules, function(x) {x[[3]]})
  M = length(redRules)
  N = nrow(strains)
  predLineages = matrix(NA, N, M)
  colnames(predLineages) = ruleStrains
  for (ind in 1:N) {
    for (index in 1:M) {
      curVal = checkRule(strains[ind,], redRules[[index]], hard = hard)
      if (is.na(curVal)) {
        next
      }
      else {
        predLineages[ind, index] = curVal
      }
    }
  }
  output = predLineages
  if (!is.null(lineages)) {
    bestLineages = ruleStrains[apply(predLineages, 1, function(x) {ifelse(all(is.na(x)), NA, which.max(x))})]
    confusion = getConfusion(lineages, bestLineages)
    output = list(predLineages, confusion)
  }
  output
}

### This function computes the fraction of a set of rules that holds for a particular strain
### If hard = TRUE, converts the result into a Boolean; otherwise, returns the fraction itself.
checkRule = function(strain, rule, hard) {
  curLoci = as.numeric(rule[1,])
  curBounds = as.numeric(rule[2,])
  curSenses = rule[3,]
  upper = which(curSenses == "U")
  lower = which(curSenses == "L")
  hold = mean(c(strain[curLoci[upper]] <= curBounds[upper], strain[curLoci[lower]] > curBounds[lower]))
  if (hard) {
    hold = (hold == 1)
  }
  hold
}

### This function prepares a set of rules in a tabular format, with respect to given loci names 
tableRules = function(rules, lociNames) {
  Lens = sapply(rules, function(x) {length(unique(x[1,]))})
  M = max(Lens)
  N = length(rules)
  Table = matrix(NA, M + 1, N + 1)
  Table[,1] = paste0("Rule", c("", 1:M))
  Table[1,-1] = names(rules)
  for (ind in 1:N) {
    curRule = rules[[ind]]
    curUnique = sort(as.integer(unique(curRule[1,])))
    for (index in 1:length(curUnique)) {
      curName = lociNames[curUnique[index]]
      cur = curRule[-1, (curRule[1,] == as.character(curUnique[index]))]
      if (class(cur) == "matrix") { # an upper bound and a lower bound
        curBounds = sort(as.integer(cur[1,]))
        curEntry = paste0(curBounds[1] + 1, " <= ", curName, " <= ", curBounds[2])
      }
      else { # a single bound
        curBound = as.integer(cur[1])
        if (cur[2] == "L") {
          curEntry = paste0(curName, " >= ", curBound + 1)
        }
        else {
          curEntry = paste0(curName, " <= ", curBound)
        }
      }
      Table[index + 1, ind + 1] = curEntry
    } 
  }
  Table
}

makeTextRules = function(redRules) {
  textRules = lapply(redRules, function(x) {
    y=x[1,]
    curLoci=altLoci[as.numeric(y)]
    curSenses=rep("<=",length(y))
    curSenses[x[3,]=="L"]=">="
    curCutoff=as.numeric(x[2,])
    curCutoff[curSenses==">="]=curCutoff[curSenses==">="]+1
    curRules=paste(curLoci, curSenses, curCutoff)
  })
  textRules
}

### This function computes the similarity of two sets of rules based on their average Jaccard index
compareRules = function(rules1, rules2) {
  Names = sort(names(rules1))
  stopifnot(sort(names(rules2)) == Names)
  indices = rep(NA, length(Names))
  names(indices) = Names
  for (name in Names) {
    curRule1 = rules1[[name]]
    curRule2 = rules2[[name]]
    jointRule = cbind(curRule1, curRule2)
    overlap = sum(duplicated(t(jointRule)))
    indices[name] = overlap/(ncol(jointRule) - overlap)
  }
  output = list(indices, mean(indices))
  output
}

### This function prepares a database from an input csv file; relCol is the column with the labels
### numLoci is the number of loci to consider; sep is the separator used for complex values
### If unique = TRUE, only one of the strains containing each unique pattern will be extracted 
### if cutoff > 0, only lineages with a number of examples exceeding the cutoff will be returned
prepareDatabase = function(dataFile = "Database.csv", relCol = 3, numLoci = 24, sep = "\\+", unique = TRUE, cutoff = 0) {
  Database = as.matrix(read.csv(dataFile))
  DatabaseR = Database[,relCol + (1:numLoci)]
  strainsDB = Database[,relCol]
  for (ind in 1:length(transTable)) {
    symbol = transTable[ind]
    DatabaseR[DatabaseR %in% c(symbol, toupper(symbol))] = ind - 1
  }
  for (ind in 1:nrow(DatabaseR)) {
    curStrain = DatabaseR[ind, ]
    complexPos = getComplexPositions(curStrain, sep = sep)
    if (length(complexPos) > 1) {
      stop("Error: the split is not unique!")
    }
    else if (length(complexPos) == 1) {
      strain1 = pickPositions(curStrain, complexPos, 1, sep = sep)
      strain2 = pickPositions(curStrain, complexPos, 2, sep = sep)
      DatabaseR[ind,] = strain1
      DatabaseR = rbind(DatabaseR, strain2)
      strainsDB[nrow(DatabaseR)] = strainsDB[ind]
    }
  }
  rownames(DatabaseR) = NULL
  DatabaseR = matrix(as.numeric(DatabaseR), nrow = nrow(DatabaseR), dimnames = dimnames(DatabaseR))
  if (unique) {
    output = extractUniqueRows(DatabaseR, repeats = FALSE, extra = strainsDB)
  }
  else {
    output = list(DatabaseR, strainsDB)
  }
  strains = output[[1]]
  lineages = output[[2]]
  Tab = table(lineages)
  goodLineages = names(Tab)[Tab >= cutoff]
  L = length(goodLineages)
  goodInds = lineages %in% goodLineages
  strains = strains[goodInds,]
  lineages = lineages[goodInds]
  output = list(strains, lineages)
  output
}