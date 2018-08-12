### Change this to the directory where your CPLEX software executable, ./cplex, is located
cplexDir = '/Users/admin/Applications/IBM/ILOG/CPLEX_Studio1261/cplex/bin/x86-64_osx/'

testDataset = rbind(c("1","2","4","5"), c("1","2","4","5"), c("2","3","5","6"), c("1","2+3","4","5+6"), c("1+2","3","4+5","6"), c("1","2","3+4","5"), c("1+2","2+3","4","5+6"), c("1+2","2+3","4+5","5+6"))

testClassification = function(Res, Dir='.', clonal=333:374, mixed=375:415, sep="\\+", metric = "constant", W = FALSE, useAll = TRUE, duplicate = FALSE, K = 1, reciprocal = TRUE, ind = "") {
	initDir = getwd()
	setwd(Dir)
	L = length(Res)
	N = length(clonal) + length(mixed)
	offset = length(clonal)
	allFiles = list.files(pattern = paste0('Simulation', ind, 'Results.*.RData'))
	redFiles = substr(allFiles, nchar('SimulationResults'), nchar(allFiles))
	matches = regexpr("[0-9]+", redFiles)
	fileNums = as.integer(substr(redFiles,as.integer(matches),as.integer(matches)+attr(matches,"match.length")-1))
	stopifnot(length(allFiles) == L)
	classification = matrix(NA, L, N)
	for (index in 1:L) {
		print(index)
		curFile = allFiles[fileNums == index] 
		load(curFile)
		finalPool = curResult[[3]]
		curStrains = Res[[index]][[1]]
		simpleStrains = matrix(as.numeric(curStrains[setdiff(1:nrow(curStrains), c(clonal, mixed)) ,]), ncol = 24)
		complexStrains = curStrains[c(clonal, mixed), ]
		weights = NULL
		if (W) {
			weights = computeWeights(simpleStrains)
		}
		uniqueS = extractUniqueRows(simpleStrains, repeats = FALSE)[[1]]
		initialNum = nrow(uniqueS)
		redPool = finalPool[1:initialNum,]
		if(any(redPool[do.call(order,as.data.frame(redPool)),]!=uniqueS[do.call(order,as.data.frame(uniqueS)),])){
			stop("Error: this should never happen!")
		}
		curClassification = computeClassification(finalPool, initialNum, complexStrains, metric = metric, weights = weights, sep = sep, useAll = useAll, duplicate = FALSE, K = 1, reciprocal = reciprocal)
		classification[index,] = curClassification[[1]]
	}
	setwd(initDir)
	classification
}

testResolutions = function(Res, clonalRes, mixedRes, Dir = '.', clonal = 333:374, mixed = 375:415, sep = "\\+", ind = "") {
	initDir = getwd()
	setwd(Dir)
	L = length(Res)
	N = length(clonal) + length(mixed)
	offset = length(clonal)
	correctness = matrix(NA, L, N)
	fractionalD = matrix(NA, L, N)
	predictions = matrix(NA, L, N)
	allFiles = list.files(pattern = paste0('Simulation', ind, 'Results.*.RData'))
	redFiles = substr(allFiles, nchar('SimulationResults'), nchar(allFiles))
	stopifnot(length(allFiles) == L)
	matches = regexpr("[0-9]+", redFiles)
	fileNums = as.integer(substr(redFiles,as.integer(matches),as.integer(matches)+attr(matches,"match.length")-1))
	for (index in 1:L) {
		curFile = allFiles[fileNums == index]
		load(curFile)
		predictions[index,] = curResult[[1]]
		allRes = curResult[[2]]
		finalPool = curResult[[3]]
		curStrains = Res[[index]][[1]]
		curLoci = apply(curStrains[c(clonal, mixed),], 1, function(x) {getComplexPositions(x, sep = sep)})
		for (ind1 in 1:length(clonal)) {
			curClonalRes = finalPool[allRes[[ind1]], ]
			correctness[index, ind1] = computeCorrectness(curClonalRes, clonalRes[[index]][[ind1]])
			fractionalD[index, ind1] = computeDistances(curClonalRes, clonalRes[[index]][[ind1]], curLoci[[ind1]])
		}
		for (ind2 in 1:length(mixed)) {
			ind2o = offset + ind2
			curMixedRes = finalPool[allRes[[ind2o]], ]
			correctness[index, ind2o] = computeCorrectness(curMixedRes, mixedRes[[index]][[ind2]])
			fractionalD[index, ind2o] = computeDistances(curMixedRes, mixedRes[[index]][[ind2]], curLoci[[ind2o]])
		}
	}
	setwd(initDir)
	output = list(correctness, fractionalD, predictions)
	output
}

computeCorrectness = function(resolvedStrains, trueResolutions) {
	n1 = nrow(resolvedStrains)/2
	n2 = nrow(trueResolutions)/2
	resolvedRed = resolvedStrains[2 * (1:n1), , drop = FALSE]
	commonStrains = sum(duplicated(rbind(trueResolutions, resolvedRed)))
	correctness = commonStrains/(n1 * n2)
	if (correctness > 1) {
		print("This should never happen")
		print(resolvedStrains)
		print(trueResolutions)
	}
	correctness
}

computeDistances = function(resolvedStrains, trueResolutions, relLoci) {
	L = length(relLoci)
	if (L == 1) {
		return(0)
	}
	n1 = nrow(resolvedStrains)/2
	n2 = nrow(trueResolutions)/2
	resolvedRed = resolvedStrains[2 * (1:n1), , drop = FALSE]
	trueRed = trueResolutions[2 * (1:n2), , drop = FALSE]
	Mat = matrix(NA, n1, n2)
	for (ind1 in 1:n1) {
		for (ind2 in 1:n2) {
			Mat[ind1, ind2] = computeSwapDistance(resolvedRed[ind1,], trueRed[ind2,], relLoci)
		}
	}
	result = mean(Mat) / (L - 1)
	result
}

computeSwapDistance = function(strain1, strain2, relLoci) {
	stopifnot(all(strain1[-relLoci] == strain2[-relLoci]))
	Len = length(relLoci)
	subStrain1 = strain1[relLoci]
	subStrain2 = strain2[relLoci]
	numDiffs = sum(subStrain1 != subStrain2)
	dist = min(numDiffs, Len - numDiffs)
	dist
}

findComplement = function(complexStrain, simpleStrain, sep = "\\+") {
	compLoci = which(complexStrain != as.character(simpleStrain))
	fullLoci = matrix(unlist(strsplit(complexStrain[compLoci], sep)), nrow = 2)
	otherLoci = fullLoci[fullLoci != matrix(simpleStrain[compLoci],ncol = length(compLoci),nrow = 2,byrow = TRUE)]
	otherStrain = simpleStrain
	otherStrain[compLoci] = otherLoci
	result = as.numeric(otherStrain)
	result
}

exactDriver = function(Directory) {
	Results = list()
	initDir = getwd()
	index = substr(Directory, nchar(Directory), nchar(Directory))
	setwd(Directory)
	for (metric in c("hamming")) {
		for (W in c(FALSE)) {
			for (byDist in c(FALSE)) {
				curResult = testHeuristic(metric = metric, sep = "\\+", W = W, byDist = byDist, index = index)
				Results[[paste0(capitalize(metric), index, rep("W", W), rep("U", !W), 1 + byDist)]] = curResult
			}
		}
	}
	for (metric in c("constant", "linear")) {
		for (W in c(TRUE, FALSE)) {
			for (byDist in c(TRUE)) {
				curResult = testHeuristic(metric = metric, sep = "\\+", W = W, byDist = byDist, index = index)
				Results[[paste0(capitalize(metric), index, rep("W", W), rep("U", !W), 1 + byDist)]] = curResult
			}
		}
	}
	setwd(initDir)
	Results
}

testHeuristic = function(metric = "linear", inds = 1:100, clonal = 333:374, mixed = 375:415, sep = "\\+", W = TRUE, exact = TRUE, integral = TRUE, byDist = FALSE, index = "0") {
	myBaseFile = paste0(substr(baseFile, 1, 10), index, substr(baseFile, 11, nchar(baseFile)))
	L = length(inds)
	complex = sort(union(clonal, mixed))
	N = length(complex)
	distMatrix = matrix(NA, L, N)
	rownames(distMatrix) = inds
	colnames(distMatrix) = complex
	for (ind in inds) {
		print(ind)
		curFile = gsub('\\.', paste0(ind, '\\.'), myBaseFile)
		if (exact) {
			curID = as.numeric(index) * length(inds) + ind - 1
			curResult = splitExact(curFile, metric = metric, weight = W, sep, integral, byDist, baseID = curID)
			curFileX = paste0(substr(curFile,1,nchar(curFile)-4),capitalize(metric),rep("W",W),1+byDist,".RData")
			save(curResult, file = curFileX)
			distMatrix[ind,] = curResult[[1]]
		}
		else {
			fullTab = processFile(curFile, ncol = 24, skip = 0, sep = sep, split = FALSE)
			curResult = splitHeuristic(fullTab, metric = metric, weight = W, sep = sep)
			curFileH = paste0(substr(curFile,1,nchar(curFile)-4),capitalize(metric),rep("W",W),"Heuristic.RData")
			save(curResult, file = curFileH)
			distMatrix[ind,] = curResult[[1]]
		}
	}
	clonalNames = as.character(clonal)
	mixedNames = as.character(mixed)
	spec = apply(distMatrix, 1, function(class) {sum(class[clonalNames])/length(clonalNames)})
	sens = apply(distMatrix, 1, function(class) {sum(1 - class[mixedNames])/length(mixedNames)})
	output = list(distMatrix, spec, sens)
	output
}

getPairwiseIntersections = function(complexStrains, sep = "\\+") {
	m = nrow(complexStrains)
	n = ncol(complexStrains)
	if (m == 1) {
		output = list(c(), c(), matrix(0, 1, 1))
		return(output)
	}
	numIntersections = matrix(0, m, m)
	allIntersections = vector("list", m - 1)
	names(allIntersections) = 1:(m - 1)
	strainsAsLists = vector("list", m)
	for (ind in 1:m) {
		if (ind %% 100 == 0) {
			print(paste("Processed", ind, "complex strains so far"))
		}
		strainsAsLists[[ind]] = lapply(strsplit(complexStrains[ind,], sep), as.numeric)
	}
	for (ind1 in 1:(m - 1)) {
		if (ind1 %% 100 == 0) {
			print(paste("Processed", ind1, "intersections so far"))
		}
		List1 = strainsAsLists[[ind1]]
		allIntersections[[ind1]] = vector("list", m - ind)
		for (ind2 in (ind1 + 1):m) {
			List2 = strainsAsLists[[ind2]]
			curIntersection = getIntersectionLists(List1, List2, specialNA = FALSE)
			curLens = sapply(curIntersection, length)
			curNumber = prod(curLens)
			numIntersections[ind1, ind2] = curNumber
			numIntersections[ind2, ind1] = curNumber
			if (curNumber > 0) {
				allIntersections[[ind1]][[as.character(ind2)]] = curIntersection
			}
		}
	}
	totalNumber = sum(numIntersections)/2
	Mat = matrix(NA, totalNumber, n + 2)
	curPos = 0
	for (ind1 in 1:(m - 1)) {
		for (ind2 in (ind1 + 1):m) {
			curLength = numIntersections[ind1, ind2]
			if (curLength > 0) {
				curRange = curPos + (1:curLength)
				curIntersection = allIntersections[[as.character(ind1)]][[as.character(ind2)]]
				Mat[curRange, 1] = ind1
				Mat[curRange, 2] = ind2
				Mat[curRange, -(1:2)] = expandOptions(curIntersection) 
				curPos = tail(curRange, 1)	
			}
		}
	}
	uniqueInts = extractUniqueRows(Mat[, -(1:2)], repeats = TRUE)
	uniqueIntersections = uniqueInts[[1]]
	uniqueSets = uniqueInts[[2]]
	uniqueIndices = sapply(uniqueSets, function(x) {sort(unique(as.vector(Mat[x, 1:2])))})
	output = list(uniqueIntersections, uniqueIndices, numIntersections)
	output
}

JaccardIndex = function(Set1, Set2, wildcard = 0) {
	if ((wildcard %in% Set1) && !(wildcard %in% Set2)) {
		Set1[Set1 == wildcard] = setdiff(Set2, Set1)[1]
	}
	if ((wildcard %in% Set2) && !(wildcard %in% Set1)) {
		Set2[Set2 == wildcard] = setdiff(Set1, Set2)[1]
	}
	num = length(intersect(Set1, Set2))
	den = length(union(Set1, Set2))
	output = num/den
	output
}

splitExact = function(filename, metric = "linear", weight = FALSE, sep = "\\+", integral = TRUE, byDist = TRUE, baseID = 0, attrTable = NULL, useFreqs = FALSE, minNew = FALSE, numLoci = 24) {
	dataset = processFile(filename, ncol = numLoci, skip = 0, sep = sep, split = FALSE)
	rownames(dataset) = 1:nrow(dataset) # added on June 7, 2015
	Distr = getDistributions(dataset, sep = sep)								# stage 1: create the splits
	n = length(Distr[[3]])
	simpleStrains = matrix(as.integer(dataset[-Distr[[3]],]), ncol = numLoci)
	W = NULL
	if (weight) {
		W = computeWeights(simpleStrains)
	}
	allComplexStrains = dataset[Distr[[3]], , drop = FALSE]
	freqs = NULL
	if (useFreqs) {
		simpleStrainsU = extractUniqueRows(simpleStrains, TRUE)
		simpleStrains = simpleStrainsU[[1]]
		freqs = sapply(simpleStrainsU[[2]], length)
	}
	numSimple = nrow(simpleStrains)
	print("Constructing the strain intersections")
	pairwiseInts = getPairwiseIntersections(allComplexStrains, sep = sep)
	adjMatrix = (pairwiseInts[[3]] > 0)
	pairwiseGraph = graph.adjacency(adjMatrix, mode = "undirected")
	comps  = clusters(pairwiseGraph)
	numClusters = comps$no
	print(paste("There are", numClusters, "connected components in this problem"))
	Results = vector("list", numClusters)
	attrTableFull = attrTable
	prevPos = 0
	for (clusterInd in 1:numClusters) {
		splitFurther = FALSE
		curSubset = (comps$membership == clusterInd)
		curSubgraph = induced.subgraph(pairwiseGraph, V(pairwiseGraph)[curSubset])
		stopifnot(is.connected(curSubgraph))
		complexStrains = allComplexStrains[curSubset, , drop = FALSE]
		attrTable = attrTableFull[curSubset]
		n = nrow(complexStrains)
		print(paste("Processing cluster number", clusterInd, "containing", n, "complex strains"))
		if (n == 1) { # SPECIAL TREATMENT			
			output = processSingleton(complexStrains, simpleStrains, metric, W, sep, filename, integral, byDist, baseID, attrTable, freqs, minNew)
		}
		else {
			splitFurther = (n > 5)
			output = processMultiples(complexStrains, simpleStrains, metric, W, sep, filename, integral, byDist, baseID, attrTable, freqs, minNew, splitFurther = splitFurther, curSubgraph = curSubgraph)
		}
		numClusters = 1
		if (splitFurther) {
			numClusters = length(output) 
		}
		if (numClusters == 1) {
			output = list(output)
		}
		for (ind in 1:numClusters) {
			curOut = output[[ind]]
			numOptima = length(curOut)
			allClassifications = sapply(curOut, function(x) {x[[1]]})
			meanClassification = rowMeans(matrix(allClassifications, ncol = numOptima))
			names(meanClassification) = names(curOut[[1]][[1]])
			baseID = baseID + numOptima + 1
			curComplexStrains = complexStrains[names(meanClassification), , drop = FALSE]
			Results[[prevPos + ind]] = list(meanClassification, curOut, curComplexStrains)
		}
		prevPos = prevPos + numClusters
	}
	fullMeanClassification = unlist(sapply(Results, function(x) {x[[1]]}, simplify = FALSE, USE.NAMES = TRUE))
	output = list(fullMeanClassification, Results)
	output
}

makeConstraint = function(indices, prefix = "U", Sum = NULL) {
	if (is.null(Sum)) {
		Sum = length(indices) - 1
	}
	Line = paste(paste0(prefix, indices), collapse = " + ")
	Line = paste(Line, "<=", Sum)
	Line
}

runFile = function(filename, Lines, baseID) {
	CPLEXfile = gsub(".csv", paste0("I", baseID, ".lp"), filename)
	fw = file(CPLEXfile, "w")
	writeLines(Lines, fw)
	close(fw)
	createScript(scriptName = "auxScript", baseID = baseID, filename = CPLEXfile)
	shortFilename = unlist(strsplit(as.character(CPLEXfile), "\\."))[1]
	scriptFile = paste0(shortFilename, "auxScript", baseID)
	system(paste0("./", scriptFile))
	# file.remove(scriptFile)
	file.remove(c(CPLEXfile, scriptFile))
}

readFile = function(filename) {
	f = file(filename, "r")
	resultLines = readLines(f, warn = FALSE)
	close(f)
	resultLines
}

parseLines = function(resultLines, prefix1 = "U", prefix2 = "C", mat = TRUE, integer = TRUE) {
	# NOTE: use "X", "Y", FALSE for DBs and "Y", "X", NULL for rule extracts!
	# NOTE: use "X", "Y", "color" for coloring programs like GridColoring.R!
	failLines = grep("^MIP - Integer infeasible", resultLines)
	if (length(failLines) > 0) {
		return(NULL)
	}
	successStart = ifelse(integer, "^MIP - Integer optimal solution", "^Dual simplex - Optimal")
	optLines = grep(successStart, resultLines)
	stopifnot(length(optLines) == 1)
	optLine = resultLines[optLines[1]]
	optPos = regexpr("=", optLine)
	optimum = as.numeric(substr(optLine, optPos + 1, nchar(optLine)))
 	usedLines = grep(paste0("^", prefix1, ".*0$"), resultLines)
 	wrongLines = grep("-", resultLines)
 	usedLines = setdiff(usedLines, wrongLines)
	usedVars = resultLines[usedLines]
	match = regexpr("[0-9]+", usedVars)
	stopifnot(all(as.integer(match) == 2))
	if (!is.null(mat) && (mat == "color" || mat == "binary")) {
		usedRows = as.integer(substr(usedVars, 2, 2 + attr(match,"match.length") - 1)) + 1
		usedVars = substr(usedVars, 2 + attr(match,"match.length"), nchar(usedVars))
		match1 = regexpr("[0-9]+", usedVars)
		stopifnot(all(as.integer(match1) == 2))
		usedCols = as.integer(substr(usedVars, 2, 2 + attr(match1,"match.length") - 1)) + 1
		usedVars = substr(usedVars, 2 + attr(match1,"match.length"), nchar(usedVars))
		match2 = regexpr("[0-9]+", usedVars)
		stopifnot(all(as.integer(match2) == 2))
		usedStrains = as.integer(substr(usedVars, 2, 2 + attr(match2,"match.length") - 1))
	}
	else {	
		usedStrains = as.integer(substr(usedVars, 2, 2 + attr(match,"match.length") - 1))
		coverLines = grep(paste0("^", prefix2, ".*0$"), resultLines)
		coverLines = setdiff(coverLines, wrongLines)
		coverVars = resultLines[coverLines]
		match1 = regexpr("[0-9]+", coverVars)
		stopifnot(all(as.integer(match1) == 2))
		compStrains = as.integer(substr(coverVars, 2, 2 + attr(match1,"match.length") - 1))
		coverVars = substr(coverVars, 2 + attr(match1,"match.length"), nchar(coverVars))
		match2 = regexpr("[0-9]+", coverVars)
		stopifnot(all(as.integer(match2) == 2))	
		coverStrains = as.integer(substr(coverVars, 2, 2 + attr(match2,"match.length") - 1))
	}
	if (is.null(mat)) {
		coverMatrix = rbind(locus = compStrains, cutoff = coverStrains, sense = substr(coverVars, 1, 1))
	}
	else if (mat == "binary") { # this will make a mistake if the last row and column are fully planted!
		allIndices = cbind(usedRows, usedCols)
		Dim = max(allIndices)
		coverMatrix = matrix(1, Dim, Dim)
		Max = max(usedStrains)
		for (ind in 1:Max) {
			curPos = (usedStrains == ind)
			coverMatrix[allIndices[curPos,]] = coverMatrix[allIndices[curPos,]] + 2 ^ (ind - 1)
		}
	}
	else if (mat == "color") { # this will make a mistake if the last row or column are fully planted!
		coverMatrix = matrix(NA, max(usedRows), max(usedCols))
		coverMatrix[cbind(usedRows, usedCols)] = usedStrains
	}
	else if (mat) {
		coverMatrix = matrix(coverStrains[order(compStrains)], nrow = 2)
	}
	else {
		coverMatrix = rbind(unused = coverStrains, closest = compStrains)
	}
	output = list(usedStrains, optimum, coverMatrix)
	output
}

runAndParse = function(filename, Lines, baseID = 0, prefix1 = "U", prefix2 = "C", mat = TRUE, integer = TRUE) {
	initDir = getwd()
	setwd(cplexDir)
	runFile(filename, Lines, baseID = baseID)
	CPLEXfile = gsub(".csv", paste0("I", baseID, ".lp"), filename)
	LOGfile = paste0(CPLEXfile, ".log")
	resultLines = readFile(LOGfile)
	setwd(initDir)
	output = parseLines(resultLines, prefix1, prefix2, mat = mat, integer = integer)
	output
}

makeCPLEXfile = function(Q, q, n, bounds, repList, closestDists, byDist, integral, freqs = NULL, free = vector()) {
	# print("Preparing a CPLEX file")
	# print(freqs)
	numConstraints = 2.5 * Q + 2 * q + n
	Lines = vector("character", numConstraints + 5)
	mapBack = convertToMap(repList, Q)
	allLengths = diff(bounds)
	repBounds = rep(1:n, allLengths)
	CInds = paste0(repBounds, "S", mapBack)
	halfBounds = rep(1:n, allLengths/2)
	halfQLower = rep(bounds[-length(bounds)],allLengths/2) + unlist(sapply(allLengths, function(x){1:(x/2)}))
	halfQUpper = rep(bounds[-length(bounds)],allLengths/2) + unlist(sapply(allLengths, function(x){x:((x/2)+1)}))
	CInds1 = paste0(halfBounds, "S", mapBack[halfQLower])
	CInds2 = paste0(halfBounds, "S", mapBack[halfQUpper])
	Lines[1] = "Minimize"
	if (byDist) {
		redRange = setdiff(1:q, free)
		Lines[2] = paste(closestDists[redRange], paste0("U", redRange), collapse = " + ")		
	}
	else {
		Lines[2] = paste(paste0("U", setdiff(which(closestDists != 0), free)), collapse = " + ")
	}
	if (!is.null(freqs)) {
		present = which(closestDists < eps)
		stopifnot(length(present) == length(freqs))
		goodFreqs = (freqs != 1)
		if (any(goodFreqs)) {
			bonus = present[goodFreqs]
			matches = match(mapBack, bonus, 0)
			redCInds = paste0(repBounds[matches != 0], "S", mapBack[matches != 0])
			Lines[2] = paste(Lines[2], "-", paste0(log(freqs[goodFreqs][matches], 2), " C", redCInds, collapse = " - "))
			stopifnot(length(grep("- 0", Lines[2])) == 0)
		}
	}
	# stage 3: create all the constraints
	Lines[3] = "Subject to"
	Range2 = 3 + 1:Q
	Lines[Range2] = paste0("C", CInds, " - U", mapBack, " <= 0")
	Range3 = 3 + Q + 1:(Q/2)
	Lines[Range3] = paste0("C", CInds1, " - C", CInds2, " = 0")
	start4 = 3 + 1.5 * Q
	for (ind in 1:q) {
		curR = repList[[ind]]
		Lines[start4 + ind] = paste0(paste0("C", repBounds[curR], "S", ind, collapse=" + ")," - U", ind, " >= 0")
	}
	start5 = 3 + 1.5 * Q + q
	for (index in 1:n) {
		curBounds = mapBack[bounds[index] + 1:allLengths[index]] # was: [bounds[index] + 1:(allLengths[index]/2)]
		Lines[start5 + index] = paste0(paste0("C", index, "S", curBounds, collapse = " + "), " = 2")
	}
	# integrality constraints or upper bounds at the end
	URange = 4 + 1.5 * Q + q + n + 1:q
	CRange = 4 + 1.5 * Q + 2 * q + n + 1:Q
	if (integral) {
		Lines[4 + 1.5 * Q + q + n] = "Binary"
		Lines[URange] = paste0("U", 1:q)
		Lines[CRange] = paste0("C", CInds)
	}
	else {
		Lines[4 + 1.5 * Q + q + n] = "Bounds"
		Lines[URange] = paste0("U", 1:q, " <= 1")
		Lines[CRange] = paste0("C", CInds, " <= 1")
	}
	Lines[length(Lines)] = "End"
	Lines
}

createScript = function(scriptName = "auxScript", baseID = 0, filename = NULL) {
	Lines = rep("", 9)
	Lines[1] = "#!/bin/bash"
	if (is.null(filename)) {
		Lines[2] = "for file in *.lp; do"
	}
	else {
		Lines[2] = paste0("for file in ", filename, "; do")
	}
	Lines[3] = "  cat <<EOF | ./cplex | tee $file.log"
	Lines[4] = "set log *"
	Lines[5] = "read $file"
	Lines[6] = "opt"
	Lines[7] = "disp sol var *"
	Lines[8] = "EOF"
	Lines[9] = "done"
	shortFilename = unlist(strsplit(as.character(filename), "\\."))[1]
	fullScriptName = paste0(shortFilename, scriptName, baseID)
	fw = file(fullScriptName, "w")
	writeLines(Lines, fw)
	close(fw)
	system(paste("chmod", "+x", fullScriptName, collapse = " "))
}

extractUniqueRows = function(Table, repeats = TRUE, extra = NULL) { 
	# Partly based on the solution at r.789695.n4.nabble.com/Count-unique-rows-columns-in-a-matrix-td844731.html
	stopifnot(!any(is.na(Table)))
	Q = nrow(Table)
	if (Q == 1) {
		output = list(Table)
		if (repeats) {
			repList = list(1)
			output = c(output, list(repList))
		}
		if (!is.null(extra)) {
			output = c(output, list(extra))
		}
		return(output)
	}
	myOrder = do.call(order, as.data.frame(Table))
	Table = Table[myOrder, , drop = FALSE]
	equalToPrevious = (rowSums(Table[-1, , drop = FALSE] != Table[-Q, , drop = FALSE]) == 0)
	lowerBounds = c(0, which(!equalToPrevious)) + 1
	output = list(Table[lowerBounds, , drop = FALSE])
	if (repeats) {
		N = nrow(output[[1]])
		repList = vector("list", N)
		upperBounds = c(lowerBounds[-1] - 1, Q)
		repList = lapply(1:N, function(x) {myOrder[lowerBounds[x] : upperBounds[x]]})
		output = c(output, list(repList))
	}
	if (!is.null(extra)) {
		stopifnot(length(extra) == Q)
		extra = extra[myOrder]
		uniqueExtra = extra[lowerBounds]
		output = c(output, list(uniqueExtra))
	}
	output
}

convertToMap = function(repList, Max = NULL) {
	if (is.null(Max)) {
		Max = max(unlist(repList))
	}
	Map = rep(NA, Max)
	for (ind in 1:length(repList)) {
		Map[repList[[ind]]] = ind
	}
	Map
}

findAllSolutions = function(complexStrains, closestDists, closestLists, Lines, filename, baseID, result, redTable) {
	n = nrow(complexStrains)
	usedStrains = result[[1]]					
	usedStrains = usedStrains[closestDists[usedStrains] >= eps]
	optimum = result[[2]]
	curCover = result[[3]]
	curCost = optimum
	curUsed = usedStrains
	# print(paste("The optimal cost is", curCost))
	iter = 1
	output = list()
	while (curCost - optimum < eps) {
		curProbs = apply(curCover, 2, function(x) {JaccardIndex(closestLists[[x[1]]], closestLists[[x[2]]])})
		curResolutions = lapply(1:n, function(x) {redTable[curCover[,x], ]})
		names(curResolutions) = rownames(complexStrains) 
		names(curProbs) = rownames(complexStrains)
		output[[iter]] = list(curProbs, curResolutions)
		# print(paste("Processing iteration", iter, "to determine additional optima"))
		Lines = c(Lines[1:3], makeConstraint(curUsed), Lines[4:length(Lines)])
		curResult = runAndParse(filename, Lines, baseID = baseID + iter)
		curStrains = curResult[[1]]
		iter = iter + 1
		if (is.null(curStrains)) {
			break
		}
		curUsed = curStrains[closestDists[curStrains] >= eps]
		curCost = curResult[[2]]
		curCover = curResult[[3]]
		# print(paste("The optimal cost is", curCost))
	}
	output
}

findForcedSplits = function(filename, Lines, initialResults, baseID) {
	# step 1: check flexibility of assignment
	optimum  = initialResults[[2]]
	curCover = initialResults[[3]]
	n = ncol(curCover)
	forced = rep(FALSE, n)
	for (ind in 1:n) {
		curLines = c(Lines[1:3], paste0("C", ind, "S", curCover[	1, ind], " = 0"), Lines[4:length(Lines)])
		curResult = runAndParse(filename, curLines, baseID = baseID + ind)
		if (is.null(curResult) || curResult[[2]] - optimum > eps) {
			forced[ind] = TRUE
		}
	}
	forced
}

processSingleton = function(complexStrains, simpleStrains, metric, W, sep, filename, integral, byDist, baseID, attrTable, freqs, minNew, freebies = NULL) {
	extraLength = 0
	if (!is.null(freebies)) {
		simpleStrains = rbind(simpleStrains, freebies)
		repeats = which(duplicated(simpleStrains))
		if (length(repeats) > 0) {
			simpleStrains = simpleStrains[-repeats, , drop = FALSE]
		}
		extraLength = nrow(freebies) - length(repeats)
		if (!is.null(freqs)) {
			freqs = c(freqs, rep(1, extraLength))
			stopifnot(length(freqs) == nrow(simpleStrains))
		}
	}
	curRow = complexStrains[1,]
	metricD = metric
	if (!byDist) {
		metricD = "discrete"
	}
	curPos = getComplexPositions(curRow, sep)
	curRes = attrTable[[1]]
	if (is.null(curRes)) {
		curRes = c(1, rep(0, length(curPos) - 1))
	}
	Mat = getCover(curRow, "min", simpleStrains, metric = metricD, weights = W, restrict = curRes, matrix = TRUE, attr = TRUE, extraLength = extraLength)
	Bests = attr(Mat, "bests")
	if (minNew && any(Mat < eps)) {
		if (any(Mat[,3] < eps)) {
			bestRes = which(Mat[,3] < eps)
			if (!is.null(freqs)) {
				if (extraLength > 0) {
					fullBests  = attr(Mat, "fullBests")
					extraBests = fullBests[bestRes]
					allSavings = sapply(extraBests, function(x) {freqs[x[[1]]] * freqs[x[[2]]]})	
				}
				else {
					allBests = Bests[bestRes]
					allSavings = sapply(allBests,   function(x) {freqs[x[[1]]] * freqs[x[[2]]]})
				}
				goodRes = (max(allSavings) - allSavings < eps)
			}
		}	
		else {
			bestRes = which(Mat[,1] < eps | Mat[,2] < eps)
			if (!is.null(freqs)) {
				Left0  = which(Mat[bestRes,1] < eps)
				Right0 = which(Mat[bestRes,2] < eps)
				if (extraLength > 0) {
					fullBests  = attr(Mat, "fullBests")
					extraBests = fullBests[bestRes]
					allSavingsL = sapply(extraBests[Left0 ], function(x) {log(freqs[x[[1]]], 2)})
					allSavingsR = sapply(extraBests[Right0], function(x) {log(freqs[x[[2]]], 2)})
				}
				else {
					allBests = Bests[bestRes]
					allSavingsL = sapply(allBests[Left0 ],   function(x) {log(freqs[x[[1]]], 2)})
					allSavingsR = sapply(allBests[Right0],   function(x) {log(freqs[x[[2]]], 2)})
				}
				bestSums = Mat[bestRes,3]
				if (any(Left0)) {
					bestSums[Left0 ] = bestSums[Left0 ] - allSavingsL
				}
				if (any(Right0)) {
					bestSums[Right0] = bestSums[Right0] - allSavingsR
				}
				goodRes = (bestSums - min(bestSums) < eps)
			}
		}
		bestRes = bestRes[goodRes]
	}
	else {
		bestRes = which(Mat[,3] < min(Mat[,3]) + eps)
	}
	Opts = rownames(Mat)[bestRes]
	numOptima = length(Opts)
	fullOpts = matrix(as.numeric(unlist(lapply(Opts, function(x) {strsplit(x,",")}))), nrow = numOptima, byrow = TRUE)
	if (!is.null(attrTable[[1]])) {
		expandOpts = matrix(as.numeric(curRes), nrow = numOptima, ncol = length(curRes), byrow = TRUE)
		expandOpts[, (curRes == 0)] = fullOpts
		curCovers1 = t(apply(expandOpts, 1, function(x) {pickPositions(curRow, curPos,     x, sep)}))
		curCovers2 = t(apply(expandOpts, 1, function(x) {pickPositions(curRow, curPos, 3 - x, sep)}))
	}
	else {
		if (nrow(Mat) == 1 && length(fullOpts) == 1) {
			curCovers1 = matrix(pickPositions(curRow, curPos, 1, sep), nrow = 1)
			curCovers2 = matrix(pickPositions(curRow, curPos, 2, sep), nrow = 1)
		}
		else {
			curCovers1 = t(apply(fullOpts, 1, function(x) {pickPositions(curRow, curPos, c(1,     x), sep)}))
			curCovers2 = t(apply(fullOpts, 1, function(x) {pickPositions(curRow, curPos, c(2, 3 - x), sep)}))
		}
	}
	redBests = Bests[bestRes]
	output = vector("list", numOptima)
	for (ind in 1:numOptima) {
		allRes = list(rbind(curCovers1[ind,], curCovers2[ind,]))
		names(allRes) = rownames(complexStrains)
		allProbs = JaccardIndex(redBests[[ind]][[1]], redBests[[ind]][[2]])
		names(allProbs) = rownames(complexStrains)
		output[[ind]] = list(allProbs, allRes)
	}
	output
}

processMultiples = function(complexStrains, simpleStrains, metric, W, sep, filename, integral, byDist, baseID, attrTable, freqs, minNew, splitFurther = FALSE, freebies = NULL, curSubgraph = NULL) {
	n = nrow(complexStrains)
	if (!is.null(attrTable)) {
		compLoci = sapply(attrTable, function(x) {sum(x == 0) + !(all(x == 0))})
	}
	else {
		compLoci = apply(complexStrains, 1, function(x) {length(grep(sep, x))})
	}
	Q = sum(2 ^ compLoci)
	numLoci = ncol(complexStrains)
	fullTable = matrix(NA, Q, numLoci)
	print(paste("There are", n, "complex strains to process"))
	bounds = 0
	for (ind in 1:n) {
		curRow = complexStrains[ind,]
		curPos = getComplexPositions(curRow, sep = sep)
		if (!is.null(attrTable)) {
			curRes = attrTable[[ind]]
			curLen = compLoci[ind]
		}
		else {
			curLen = length(curPos)
		}
		curRange = tail(bounds,1) + 1:(2^curLen)
		if (!is.null(attrTable) && !all(curRes == 0)) {
			curIndices = matrix(curRes, 2^(curLen - 1), length(curPos), byrow = TRUE)
			curIndices = rbind(curIndices, matrix(3 - curRes, 2^(curLen - 1), length(curPos),byrow=TRUE))
			if (any(curRes == 0)) {
				curFree = which(curRes == 0)
				curCube = hcube(rep(2, curLen - 1))
				curIndices[, curFree] = rbind(curCube, curCube)
			}
		}
		else {
			curIndices = hcube(rep(2, curLen))
		}
		fullTable[curRange,] = t(apply(curIndices,1,function(x){pickPositions(curRow,curPos,x,sep=sep)}))
		bounds = c(bounds, tail(curRange, 1))
	}
	resultU = extractUniqueRows(fullTable, repeats = TRUE)					
	redTable = resultU[[1]]
	repList = resultU[[2]]
	q = nrow(redTable) # indexes the U variables
	# print(paste("There are", q, "unique simple strains to process"))
	if (nrow(simpleStrains) > 0) {
		closestDistsRed = !duplicated(rbind(simpleStrains, redTable))
		closestDistsRed = as.numeric(closestDistsRed[-(1:nrow(simpleStrains))])
		stopifnot(length(closestDistsRed) == q)
	}
	else {
		closestDistsRed = rep(1, q)
	}
	if (!byDist) {
		closestDists = closestDistsRed
	}
	else {
		closestDists = rep(NA, q)
		closestLists = vector("list", q)
		if (!is.null(freqs)) {
			matchPos = rep(0, q)
		}
		print(paste("There are", q, "unique strains to process"))
		for (ind in 1:q) {
			if (ind %% 1000 == 0) {
				print(paste("Processing unique strain number", ind))
			}
			curClosest = getClosestDistance(redTable[ind,], simpleStrains, metric = metric, weights=W)
			closestDists[ind]  = curClosest$distance
			closestLists[[ind]] = curClosest$bests
			if (!is.null(freqs) && (curClosest$distance < eps)) {
				matchPos[ind] = curClosest$index
			}
		}
		if (!is.null(freqs)) {
			curFreqs = freqs[matchPos]
		}
	}
	free = vector()
	if (!is.null(freebies)) {
		numFreebies = nrow(freebies)
		free = which(duplicated(rbind(freebies, redTable))) - numFreebies
	}
	if (minNew) {
		Lines0  = makeCPLEXfile(Q, q, n, bounds, repList, closestDistsRed, FALSE, integral, NULL, free = free)
		result0 = runAndParse(filename, Lines0, baseID = baseID)
		optimum = result0[[2]]
		newStrains = setdiff(which(closestDistsRed > eps), free)
		Lines = makeCPLEXfile(Q, q, n, bounds, repList, closestDists, byDist, integral, curFreqs)
		Lines = c(Lines[1:3], makeConstraint(newStrains, Sum = optimum), Lines[4:length(Lines)])
	}
	else {
		Lines = makeCPLEXfile(Q, q, n, bounds, repList, closestDists, byDist, integral, curFreqs, free = free)
	}
	result = runAndParse(filename, Lines, baseID = baseID)
	if (splitFurther) {
		forced = findForcedSplits(filename, Lines, result, baseID)
		baseID = baseID + n + 1
		if (all(forced)) {
			output = findAllSolutions(complexStrains,closestDists,closestLists,Lines,filename,baseID,result,redTable)
		}
		else {
			freebies = redTable[unique(as.vector(result[[3]][,forced])), , drop = FALSE]	
			redComplexStrains  = complexStrains[!forced, , drop = FALSE]
			redAttrTable = attrTable[!forced]
			redSubgraph = induced.subgraph(curSubgraph, V(curSubgraph)[!forced])
			redComps = clusters(redSubgraph)
			numClusters = redComps$no
			output = vector("list", numClusters)
			print(paste("Further subdividing the problem into", numClusters, "components"))
			for (clusterInd in 1:numClusters) {
				miniSubset = (redComps$membership == clusterInd)
				miniSubgraph = induced.subgraph(redSubgraph, V(redSubgraph)[miniSubset])
				stopifnot(is.connected(miniSubgraph))
				miniComplexStrains = redComplexStrains[miniSubset, , drop = FALSE]
				miniAttrTable = redAttrTable[miniSubset]
				if (nrow(miniComplexStrains) == 1) {
					output[[clusterInd]] = processSingleton(miniComplexStrains, simpleStrains, metric, W, sep, filename, integral, byDist, baseID = baseID, attrTable = miniAttrTable, freqs = freqs, minNew = minNew, freebies = freebies)
				}
				else {
					output[[clusterInd]] = processMultiples(miniComplexStrains, simpleStrains, metric, W, sep, filename, integral, byDist, baseID = baseID, attrTable = miniAttrTable, freqs = freqs, minNew = minNew, freebies = freebies)
				}
			}
			if (any(forced)) {
				forcedRes = lapply(which(forced), function(x) {redTable[result[[3]][,x], ]})
				names(forcedRes) = rownames(complexStrains)[forced]
				forcedCovers = result[[3]][ ,forced, drop = FALSE]
				forcedP = apply(forcedCovers, 2, function(x) {JaccardIndex(closestLists[[x[1]]], closestLists[[x[2]]])})
				names(forcedP) = rownames(complexStrains)[forced]
				output[[numClusters + 1]] = list(list(forcedP, forcedRes))
			}
		}
	}
	else {
		output = findAllSolutions(complexStrains, closestDists, closestLists, Lines, filename, baseID, result, redTable)
	}
	output
}