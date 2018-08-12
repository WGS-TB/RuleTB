require(combinat)
require(Hmisc)
#source("NeighborJoin.R")
#source("CNVImpute.R")
#source("BeastInterface.R")
#source("Alternative.R")

MISSING = 0
SIMPLE  = 1
MIXED   = 2
CLONAL  = 3
MALFORMED  = -1
TOOMISSING = -2

badSymbols  = c("", "_", "-", "?", "c", "m", "x", "EMPTY", "%") # missing symbols
maxCNV 		= 15 # largest number of copy number variants - change if needed
Reciprocals = 1/(1:maxCNV)
Harmonics   = c(cumsum(Reciprocals), 0)
eps = 1e-8

Loci = as.character(c(154, 580, 960, 1644, 2059, 2531, 2687, 2996, 3007, 3192, 4348, 802, 424, 577, 1955, 2163))
Loci[length(Loci)] = paste0(Loci[length(Loci)], "b")
Loci = c(Loci, as.character(c(2165, 2347, 2401, 2461, 3171, 3690, 4156, 4052)))
altLoci = c("MIRU02","ETRD","MIRU10","MIRU16","MIRU20","MIRU23","MIRU24","MIRU26","MIRU27","ETRE","MIRU39","MIRU40", "Mtub04","ETRC","Mtub21","QUB11b","ETRA","Mtub29","Mtub30","ETRB","Mtub34","Mtub39","QUB4156","QUB26")
transTable = c(as.character(0:9), letters[1:6])

computeConsensus = function(directory, clonalNames=as.character(333:374), mixedNames=as.character(375:415), size=8){
	complex = sort(union(clonal, mixed))
	N = length(complex)
	allFiles = list.files('.', pattern = 'BeastSimulation.*.RData')
	match = regexpr("[0-9]+", allFiles)
	inds = as.integer(substr(allFiles, as.integer(match), as.integer(match) + attr(match, "match.length") - 1))
	Tab = table(inds)
	goodInds = sort(as.integer(names(Tab)[Tab==size]))
	L = length(goodInds)
	allCalls = matrix(0, L, N)
	rownames(allCalls) = goodInds
	colnames(allCalls) = c(clonalNames, mixedNames)
	scorer = function(class) {Res = sum(class[clonalNames]) + sum(1 - class[mixedNames]); Res}
	for (index in 1:L) {
		print(index)
		curFiles = allFiles[inds == goodInds[index]]
		for (file in curFiles) {
			load(file)
			curRow = output[[2]]
			curNames = names(curRow)
			allCalls[index, curNames] = allCalls[index, curNames] + curRow
		}
	}
	allClasses = (allCalls >= size/2)
	allScores = apply(allClasses, 1, scorer)
	output = list(allCalls, allScores)
	output
}

computeStats = function(clonalNames=as.character(333:374), mixedNames=as.character(375:415), N=100, extra=0) {
	List1 = vector("list", N)
	List2 = vector("list", N)
	allFiles = list.files('.', pattern = paste0('^Simulation.*.RData'))
	scorer1 = function(class) {Res = sum(class[clonalNames]); Res}
	scorer2 = function(class) {Res = sum(1 - class[mixedNames]); Res}
	for (file in allFiles) {
		print(file)
		load(file)
		match = regexpr("[0-9]+", file)
		ind = as.integer(substr(file, as.integer(match), as.integer(match) + attr(match, "match.length") - 1))
		curScores1 = matrix(NA, length(allResults[[1]]), length(allResults))
		rownames(curScores1) = c("Standard", "Baseline", "Beast")
		colnames(curScores1) = names(allResults)
		curScores2 = curScores1
		for (index in 1:length(allResults)) {
			curSub = allResults[[index]]
			curScores1[,index] = c(scorer1(curSub[[1]]), scorer1(curSub[[2]]), scorer1(curSub[[3]][[2]]))
			curScores2[,index] = c(scorer2(curSub[[1]]), scorer2(curSub[[2]]), scorer2(curSub[[3]][[2]]))
		}
		List1[[ind - extra]] = cbind(List1[[ind - extra]], curScores1)
		List2[[ind - extra]] = cbind(List2[[ind - extra]], curScores2)
	}
	bigMatrix1 = matrix(NA, N, length(List1[[1]]))
	colnames(bigMatrix1) = as.vector(outer(rownames(List1[[1]]), colnames(List1[[1]]), paste0))
	bigMatrix2 = bigMatrix1
	for (ind in 1:N) {
		curItem1 = List1[[ind]]
		curRow1 = as.vector(curItem1)
		curNames1 = as.vector(outer(rownames(curItem1), colnames(curItem1), paste0))
		bigMatrix1[ind, curNames1] = curRow1
		curItem2 = List2[[ind]]
		curRow2 = as.vector(curItem2)
		curNames2 = as.vector(outer(rownames(curItem2), colnames(curItem2), paste0))
		bigMatrix2[ind, curNames2] = curRow2
	}
	output = list(bigMatrix1, bigMatrix2)
	output
}

extraDriver = function(clonalNames=as.character(333:374), mixedNames=as.character(375:415), List=NULL, once=FALSE) {
	numbers = List[1,]
	lengths = List[2,]
	L = ncol(List)
	List = vector("list", L)
	scorer = function(class) {Res = sum(class[clonalNames]) + sum(1 - class[mixedNames]); Res}
	for (ind in 1:L) {
		curNumber = numbers[ind]
		curLength = lengths[ind]
		curFile = paste0('SimulationResults', curNumber, '.csv')
		curIndex = 9 - curLength # there are 8 settings total!
		print(curFile)
		if (once) {
			print(paste("There is", 1, "setting to process"))
			Result = metaDriver(curFile, curIndex, curIndex)
		}
		else {
			print(paste("There are", curLength, "settings to process"))
			Result = metaDriver(curFile, curIndex)
		}
		curScores = matrix(NA, length(Result[[1]]), length(Result))
		rownames(curScores) = c("Standard", "Baseline", "Beast")
		colnames(curScores) = names(Result)
		for (index in 1:length(Result)) {
			curSub = Result[[index]]
			curScores[,index] = c(scorer(curSub[[1]] - 2), scorer(curSub[[2]]), scorer(curSub[[3]][[2]]))
		}
		List[[ind]] = curScores
	}
	save(List, file = "Scores.RData")
	List	
}

superDriver = function(directory, clonalNames=as.character(333:374), mixedNames=as.character(375:415), done=0, upto = 100) {
	setwd(directory)
	scorer = function(class) {Res = sum(class[clonalNames]) + sum(1 - class[mixedNames]); Res}
	allFiles = list.files('.', pattern = 'Simulation.*.csv')
	allFiles = allFiles[(done + 1):min(upto,length(allFiles))]
	L = length(allFiles)
	List = vector("list", L)
	for (ind in 1:L) {
		print(paste("Processing file number", ind))
		file = allFiles[ind]
		print(file)
		Result = metaDriver(file)
		curScores = matrix(NA, length(Result[[1]]), length(Result))
		rownames(curScores) = c("Standard", "Baseline", "Beast")
		colnames(curScores) = names(Result)
		for (index in 1:length(Result)) {
			curSub = Result[[index]]
			curScores[,index] = c(scorer(curSub[[1]] - 2), scorer(curSub[[2]]), scorer(curSub[[3]][[2]]))
		}
		List[[ind]] = curScores
	}
	save(List, file = "Scores.RData")
	List
}

metaDriver = function(filename, minIndex = 1, maxIndex = 8, hybrid = TRUE) {
	modes = c("min", "max")
	metrics = c("constant", "linear")
	boths = c(FALSE)
	weighs = c(TRUE, FALSE)
	analyses = c(FALSE)
	conservatives = c(TRUE)
	trees = c(FALSE)
	repeateds = c(FALSE)
	skipBoth = c(0)
	skipOne = c(0)
	allResults = list()
	Groupings = list()
	index = 1
	for (mode in modes) {
		for (metric in metrics) {
			for (both in boths) {
				for (weigh in weighs) {
					for (analysis in analyses) {
						for (conservative in conservatives) {
							for (tree in trees) {
								for (repeated in repeateds) {
									if (both) {
										for (skip in skipBoth) {
											curString = paste0(capitalize(mode), capitalize(metric), "Both", weigh, skip)
											redString = curString # paste0("Both", weigh, skip)
											print(curString)
											if (hybrid) {
												allResults[[curString]] = splitHybrid(filename, metric = metric, weight = weigh, integral = TRUE, byDist = (mode == "max"))
											}
											else {
												allResults[[curString]] = driver(filename, skip = skip, metric = metric, mode = mode, both = both, weigh = weigh, analysis = analysis, conservative = conservative, tree = tree, repeated = repeated)
											}
											Groupings[[redString]] = c(Groupings[[redString]], index)
											index = index + 1
										}
									}
									else {
										for (skip in skipOne) {
											curString = paste0(capitalize(mode), capitalize(metric), "Solo", weigh, skip)
											redString = curString # paste0("Solo", weigh, skip)
											if (index >= minIndex && index <= maxIndex) {
												print(curString)
												if (hybrid) {
													allResults[[curString]] = splitHybrid(filename, metric = metric, weight = weigh, integral = TRUE, byDist = (mode == "max"))
												}
												else {
													allResults[[curString]] = driver(filename, skip = skip, metric = metric, mode = mode, both = both, weigh = weigh, analysis = analysis, conservative = conservative, tree = tree, repeated = repeated)
												}
												Groupings[[redString]] = c(Groupings[[redString]], index)
											}
											else {
												print(paste("Skipping", curString))
											}
											index = index + 1
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	baseFilename = substr(filename, 1, nchar(filename) - 4)
	if (minIndex > 1 || maxIndex < 8) {
		baseFilename = paste0(baseFilename, 'L', minIndex - 1)
	}
	save(allResults, file = paste0(baseFilename, ".RData"))
	# extFilename = paste0(baseFilename, "Results.txt")
	# allCalls = processCalls(allResults, Groupings, filename = extFilename)
	allResults
}

processAnalysis = function(Results, Calls) {
	L = length(Results)
	Output = vector("list", L)
	names(Output) = names(Results)
	for (ind in 1:L) {
		print(ind)
		curResult = Results[[ind]]
		curCalls  = Calls[, ind]
		N = length(curCalls)
		curMat = matrix(NA, N, 3, dimnames = list(rownames(Calls), c("Call", "Imputed Simple", "Imputed Complex")))
		curMat[,1] = curCalls
		impSimple = curResult[[1]] > 0
		curMat[names(impSimple), 2] = impSimple
		impComplex = curResult[[2]] > 0
		curMat[names(impComplex), 3] = impComplex
		goodRows = which(curMat[,1] == "mixed" | curMat[,1] == "clonal")
		goodNames = rownames(curMat)[goodRows]
		redMat = curMat[goodRows, , drop = FALSE]
		Tab1 = table(redMat[,1], redMat[,2])
		Tab2 = table(redMat[,1], redMat[,3])
		if (length(curResult) == 5) { # includes both cultures
			G = length(goodRows)
			compareMat = matrix(NA, 4 * G, 24)
			compareMat[4 * (1:G) - 3, ]  = curResult[[4]][goodNames,]
			compareMat[4 * (1:G) - 2, ] = curResult[[5]][goodNames,]
			curCollection = curResult[[3]]
			allNames = rownames(curCollection)
			for (index in 1:G) {
				curName = goodNames[index]
				goodColls = which(allNames == curName)
				compareMat[4 * index - 1, ] = curCollection[goodColls[1], ]
				if (length(goodColls) > 1) {
					compareMat[4 * index, ] = curCollection[goodColls[2], ]
				}
			}
			Output[[ind]] = list(redMat, Tab1, Tab2, compareMat)
		}
		else {
			Output[[ind]] = list(redMat, Tab1, Tab2)
		}
	}
	Output
}

processCalls = function(Results, Groupings, filename) {
	L = length(Groupings)
	for (ind in 1:L) {
		curGroup = Groupings[[ind]]
		curResults = Results[curGroup]
		curStandard = lapply(curResults, function(x) {x[[1]]})
		if (ind == 1) {
			numPatients = length(curStandard[[1]])
			Tab = matrix(MISSING, numPatients, L * 3)
			colnames(Tab) = outer(c("Standard", "Baseline", "Beast"), names(Groupings), paste0)
			rownames(Tab) = names(curStandard[[1]])
		}
		standardColumn = 3 * (ind - 1) + 1
		Tab[, standardColumn] = curStandard[[1]]
		if (sum(!duplicated(t(matrix(unlist(curStandard), ncol = length(curStandard))))) > 1) {
			stop(paste("Error: not all standard calls agree for", names(Groupings)[ind]))
		}
		baselineColumn = 3 * (ind - 1) + 2
		Tab[Tab[, standardColumn] == SIMPLE, baselineColumn] = SIMPLE
		curBaseline = lapply(curResults, function(x) {x[[2]]})
		curBaseline = lapply(curBaseline, function(x) {y = x; y[x < 0.5] = MIXED; y[x >= 0.5] = CLONAL; y})
		for (index in 1:length(curBaseline)) { # Use the maximum value, so that CLONAL overrides MIXED
			curCalls = curBaseline[[index]]
			Tab[names(curCalls), baselineColumn] = pmax(Tab[names(curCalls), baselineColumn], curCalls)
		}
		beastColumn = 3 * (ind - 1) + 3
		Tab[Tab[, standardColumn] == SIMPLE, beastColumn] = SIMPLE
		curBeast = lapply(curResults, function(x) {x[[3]][[2]]})
		curBeast = lapply(curBeast, function(x) {y = x; y[x < 0.5] = MIXED; y[x >= 0.5] = CLONAL; y})
		for (index in 1:length(curBeast)) { # Use the maximum value, so that CLONAL overrides MIXED
			curCalls = curBeast[[index]]
			Tab[names(curCalls), beastColumn] = pmax(Tab[names(curCalls), beastColumn], curCalls)
		}
	}
	TabR = Tab
	TabR[TabR == MISSING] = "missing"
	TabR[TabR == SIMPLE]  = "simple"
	TabR[TabR == MIXED]   = "mixed"
	TabR[TabR == CLONAL]  = "clonal"
	write.table(TabR, file = filename, quote = FALSE, sep = "\t")
	Tab
}

computeDiversity = function(Vector) { # computes the inverse of the Gaston-Hunter diversity index
	Tab = table(Vector)
	N = sum(Tab)
	if (N <= 1) {
		stop("Error: the vector has too few non-NA entries for a meaningful diversity measure to be computed!")
	}
	Num = sum(Tab * (Tab - 1))
	Den = N * (N - 1)
	index = Num / Den
	index
}

computeMeanStep = function(Vector) { # computes the mean number of steps separating two entries of a vector
	Tab = table(Vector)
	N = sum(Tab)
	if (N <= 1) {
		stop("Error: the vector has too few non-NA entries for a meaningful diversity measure to be computed!")
	}
	Vals = as.numeric(names(Tab))
	L = length(Vals)
	Sum = 0
	if (L == 1) {
		return(Sum)
	}
	for (ind1 in 1:(L - 1)) {
		for (ind2 in (ind1 + 1):L) {
			Sum = Sum + Tab[ind1] * Tab[ind2] * (Vals[ind2] - Vals[ind1])
		}
	}
	meanStep = 2 * Sum / (N * (N - 1))
	meanStep
}

computeWeights = function(Matrix, generalized = FALSE) { # computes the column weights of a given matrix
	if (generalized) {
		weights = apply(Matrix, 2, function(x) {1/computeMeanStep(x)})
	}
	else {
		weights = apply(Matrix, 2, computeDiversity)
	}
	weights
}

driver = function(filename = "MIRU_Latest.csv", ncol = 24, skip = 0, sep = "\\+", metric = "linear", mode = "min", both = TRUE, weigh = FALSE, analysis = FALSE, conservative = FALSE, tree = FALSE, repeated = FALSE) {
	print("Parsing the input file")
	if (both) {
		firstCulture  = processFile(filename = filename, ncol = ncol, skip = skip, sep = sep, split = FALSE)
		skip2 = skip + ncol # shift right!
		secondCulture = processFile(filename = filename, ncol = ncol, skip = skip2, sep = sep, split = FALSE)
		N = nrow(firstCulture)
		if (nrow(secondCulture) != N) {
			stop(paste("Error: the row numbers", nrow(firstCulture), "and", nrow(secondCulture), "do not match!"))
		}
		standardCalls = rep(NA, N)
		names(standardCalls) = rownames(firstCulture)
		List = vector("list", N)
		names(List) = rownames(firstCulture)
		default = NA
		if (conservative) {
			default = findMostCommonLoci(rbind(firstCulture, secondCulture))
		}
		for (ind in 1:N) {
			strain1 =  firstCulture[ind, ]
			strain2 = secondCulture[ind, ]
			standardCalls[ind] = standardCall(rbind(strain1, strain2), sep = sep)
			List[[ind]] = findPairwiseCover(strain1, strain2, sep = sep, conservative = conservative, default = default)
		}
		redList = lapply(List, function(x) {x[[1]]})
		redLengths = sapply(redList, length)
		redList = redList[redLengths > 0]
		fullTable = lapply(redList, function(x) {matrix(unlist(x), ncol = ncol, byrow = TRUE)})
		attrTable = lapply(redList, function(x) {lapply(x, function(y) {attr(y, "restrict")})})
		attrTable = do.call(c, attrTable)
		attrTable = attrTable[sapply(attrTable, length) > 0]
		fullTable = do.call(rbind, fullTable)
		dimnames(fullTable) = list(rep(names(redList), redLengths[redLengths > 0]), c())
	}
	else {
		fullTable = processFile(filename = filename, ncol = ncol, skip = skip, sep = sep, split = FALSE)
		standardCalls = sapply(split(fullTable, 1:nrow(fullTable)), function(x) {standardCall(matrix(x, nrow = 1))})
		names(standardCalls) = rownames(fullTable)
		fullTable = fullTable[rowSums(is.na(fullTable)) < ncol(fullTable), ] # removes completely missing rows
	}
	splitTable = separateTable(fullTable, makeAttr = !both)
	if (!both) {
		attrTable = splitTable[[3]]
	}
	complexStrains = splitTable[[2]]
	Collection = matrix(as.integer(splitTable[[1]]), ncol = ncol, dimnames = dimnames(splitTable[[1]]))
	if (analysis) {
		simpleMissing  = rowSums(is.na(Collection))
		complexMissing = rowSums(is.na(complexStrains))
	}
	# print("Performing imputation")
	weights = NULL
	if (weigh) {
		weights = computeWeights(Collection)
	}
	if (!conservative || !both) {
		imputedIntervals = imputeTree(Collection, metric = metric, weights = weights)
		imputedPos = imputedIntervals[[1]]
		imputedVals = round(rowMeans(imputedIntervals[[2]]))
		Collection[imputedPos] = imputedVals
	}
	n = nrow(Collection)
	# print(paste("There are", n, "simple strains in the dataset"))
	N = nrow(complexStrains)
	coverStrains = matrix(NA, nrow = 2 * N, ncol = ncol)
	# print("Identifying the optimal covers for complex strains")
	print(paste("There are", N, "complex strains to be processed"))
	# print(paste("They have a total of", sum(is.na(complexStrains)), "missing entries"))
	for (ind in 1:N) {
		if (ind %% 10 == 0) {
			print(paste("Processed", ind, "strains so far"))
		}
		currentVector = complexStrains[ind,]
		currentAttr   = attrTable[[ind]]
		currentCover = getCover(currentVector, mode = mode, Collection = Collection, sep = sep, metric = metric, weights = weights, restrict = currentAttr)
		coverStrains[2 * ind - 1,]  = currentCover[[1]]
		coverStrains[2 * ind,] 		= currentCover[[2]]
	}
	rownames(coverStrains) = rep(rownames(complexStrains), each = 2)
	fullCollection = rbind(Collection, coverStrains)
	if (analysis) {
		result = list(simpleMissing, complexMissing, fullCollection)
		if (both) {
			result = c(result, list(firstCulture, secondCulture))
		}
	}
	else {
		Rownames = rownames(fullCollection)
		# print(paste("There are", length(unique(Rownames)), "usable patients out of a total of", nrow(firstCulture)))
		complexTable = makeComplexTable(Rownames)
		rownames(fullCollection) = 1:nrow(fullCollection)
		Input = createPhylipInput(fullCollection)
		redFilename = substr(filename, 1, nchar(filename) - 4)
		phylipFilename = paste0("Phylip", redFilename, capitalize(mode), capitalize(metric), skip, rep("Both", both),  rep("Weighted", weigh), ".txt")
		writePhylipFile(Input[[1]], paste(phylipDir, phylipFilename, sep = "/"))
		write.table(as.matrix(Input[[2]]), file = phylipFilename, sep = "\t", row.names = FALSE)
		# print("Computing the distance matrix")
		distanceMatrix = getPairwiseDistances(fullCollection, metric = metric, weights = weights)
		baseFilename = paste0("Beast", redFilename, capitalize(mode), capitalize(metric), skip, rep("Both", both),  rep("Weighted", weigh)) # the three items added at the end of the filename is the new way!
		distFile = paste0(baseFilename, '.distances')
		write.table(distanceMatrix, file = distFile, row.names = Rownames, col.names = Rownames, sep = "\t", quote = FALSE)
		print("Constructing the neighbor-joining tree")
		NJtree = neighborJoin(distanceMatrix)
		# print("Reconstructing optimal assignments")
		reconstruction = intervalReconstruct(fullCollection, NJtree)
		# print("Computing probabilities of mixed infection")
		numComplex = length(complexTable)
		# print(paste("There are", numComplex, "complex strains to be processed"))
		ProbVector = rep(NA, numComplex)
		for (ind in 1:numComplex) {
			indices = complexTable[[ind]]
			trueIndices = which(names(NJtree) %in% as.character(indices))
			ProbVector[ind] = sum(getMultiProbabilityEvolved(reconstruction, NJtree, trueIndices))
		}
		outputT = ProbVector
		names(outputT) = names(complexTable)
		outputB = driverBeast(fullCollection, complexTable, baseFilename, NJtree, analysis = tree, repeated = repeated)
		# outputB = c(NULL, NULL) # temporary!!
		if (tree) {
			result = outputB
		}
		else {
			result = list(standardCalls, outputT, outputB)
		}
	}
	result
}

standardCall = function(strains, sep = "\\+") {
	if (all(is.na(strains))) {
		return(MISSING)
	}
	diverse = unlist(apply(strains, 1, function(x) {getComplexPositions(x, sep = sep)}))
	variable = which(apply(strains, 2, function(x) {length(unique(na.omit(x))) > 1}))
	allLoci = union(diverse, variable)
	L = length(allLoci)
	if (L == 0) {
		return(SIMPLE)
	}
	if (L == 1) {
		return(CLONAL)
	}
	return(MIXED)
}

makeComplexTable = function(Rownames) {
	uniqueNames = unique(Rownames)
	N = length(uniqueNames)
	Match = match(Rownames, uniqueNames)
	Table = split(1:length(Rownames), Match)
	names(Table) = uniqueNames
	Lengths = sapply(Table, length)
	redTable = Table[Lengths > 1]
	redTable
}

testSimilarity = function(Table1, Table2) {
	similar = c()
	for (ind in 1:nrow(Table1)) {
		strain1 = Table1[ind,]
		strain2 = Table2[ind,]
		curSim = TRUE
		for (p in 2:(ncol(Table1)-1)) {
			if (!(is.na(strain1[p])) && !(is.na(strain2[p])) && !(strain1[p] == strain2[p])) {
				if (!(is.na(strain1[p-1])) && !(is.na(strain2[p-1])) && !(strain1[p] == strain2[p-1] && strain1[p-1] == strain2[p])) {
					if (!(is.na(strain1[p+1])) && !(is.na(strain2[p+1])) && !(strain1[p] == strain2[p+1] && strain1[p+1] == strain2[p])) {
						curSim = FALSE
						break
					}
				}
			}
		}
		if (curSim && !(all(strain1 == strain2, na.rm = TRUE))) {
			similar = c(similar, ind)
		}
	}
	similar
}

processFile = function(filename, ncol, skip, sep, split = TRUE, filterNA = FALSE, hex = FALSE) {
	Table = read.csv(filename, stringsAsFactors = FALSE, check.names = FALSE)
	if (hex) {
		for (ind in 1:6) {
			Table[Table == LETTERS[ind]] = 9 + ind
			Table[Table == letters[ind]] = 9 + ind 
		}
	}
	for (symbol in badSymbols) {
		Table[Table == symbol] = NA
	}
	for (i in 1:9) {
		Table[Table == paste0(i,"s")] = as.character(i)
	}
	Table[Table == ">1400bp"] = as.character(maxCNV)
	allNA = apply(Table, 1, function(x) {all(is.na(x))}) # filter out all-NA rows
	Table = Table[!allNA, ]
	rownames(Table) = Table[,1]
	allSkip = 1:(skip + 1)
	if (skip + ncol + 2 <= ncol(Table)) {
		allSkip = c(allSkip, (skip + ncol + 2):ncol(Table))
	}
	Table = Table[,-allSkip]
	if (filterNA) {
		allNA = apply(Table, 1, function(x) {all(is.na(x))}) # filter out all-NA rows again!
		Table = Table[!allNA, ]
	}
	colOrder = order(colnames(Table))
	Table = Table[, colOrder]
	Table = as.matrix(Table)
	coll = substr(sep, 2, 2)
	for (ind1 in 1:nrow(Table)) {
		for (ind2 in 1:ncol(Table)) {
			if (!is.na(Table[ind1,ind2])) {
				# NOTE: the separated vales are now forced to be in increasing order (June 7, 2015)
				cur = try(paste(sort(as.integer(unlist(strsplit(Table[ind1,ind2],sep)))),collapse = coll))
				if (class(cur) == "try-error") {
					return(paste('{"code":"-1","message":"unrecognized symbol"',Table[ind1,ind2],'}'))
				}
				else {
					Table[ind1, ind2] = cur
				}
			}
		}
	}
	if (split) {
		output = separateTable(Table)
	}
	else {
		output = Table
	}
	output
}

separateTable = function(Table, makeAttr = FALSE) {
	# This function separates a table into simple and complex strains
	mixedInds = which(nchar(Table) >= 3, arr.ind = TRUE) # NOTE: Assumes there are never more than 99 CNVs per locus
	mixedTab = table(mixedInds[,1])
	mixedRows = as.numeric(names(mixedTab))
	if (makeAttr) {
		L = length(mixedRows)
		attrTable = vector("list", L)
		for (ind in 1:L) {
			attrTable[[ind]] = rep(0, mixedTab[ind])
		}
	}
	complexTable = Table[ mixedRows,]
	simpleTable  = Table[-mixedRows,]
	output = list(simpleTable, complexTable)
	if (makeAttr) {
		output = c(output, list(attrTable))
	}
	output
}

getDistance = function(Vector, Matrix, metric, weights = NULL, coeffs = NULL) {
	# This function computes the constant or linear distance metric between 2 vectors
	# Note that the second argument can also be a matrix where each row is a vector;
	# in this case, the output is a vector with entry i containing distance to row i.
	# If a vector of weights (of the same length as Vector) is given, it will be used.
	m = nrow(Matrix)
	n = ncol(Matrix)
	distance = rep(NA, m)
	repVector = matrix(Vector, nrow = m, ncol = n, byrow = TRUE)
	repWeights = 1
	if (!is.null(weights)) {
		repWeights = matrix(weights, nrow = m, ncol = n, byrow = TRUE)
	}
	if (metric == "discrete") {
		distance = apply(repVector != Matrix, 1, any)
	}
	else if (metric == "hamming" || metric == "categorical") {
		distance = rowSums((repVector != Matrix) * repWeights, na.rm = TRUE)
	}
	else if (metric == "goldstein") {
		distance = rowSums((repVector - Matrix)^2 * repWeights / n, na.rm = TRUE)
	}
	else if (metric == "constant" || metric == "hybrid") {
		if (is.null(coeffs)) {
			distance = rowSums(abs(repVector - Matrix) * repWeights, na.rm = TRUE)
		}
		else {
			distance = rowSums(matrix(coeffs[abs(repVector - Matrix) + 1], nrow = m) * repWeights, na.rm = TRUE)
		}
		if (metric == "hybrid") {
			MULT = 1 + log(ncol(Matrix), 2)
			distance = distance + MULT * rowSums((repVector != Matrix) * repWeights, na.rm = TRUE)
		}
	}
	else if (metric == "linear") {
		minMatrix = pmin(repVector, Matrix, na.rm = FALSE)
		minMatrix[minMatrix == 0] = maxCNV + 1 # Technical trick
		maxMatrix = pmax(repVector, Matrix, na.rm = FALSE)
		maxMatrix[maxMatrix == 0] = maxCNV + 1 # Technical trick
		preDistance = Harmonics[maxMatrix] - Harmonics[minMatrix]
		distance = rowSums(matrix(preDistance, nrow  = m, ncol = n) * repWeights, na.rm = TRUE)
	}
	else if (metric == "maximum") {
		distance = apply(abs(repVector - Matrix) * repWeights, 1, function(x) {max(x, na.rm = TRUE)})
	}
	else {
		print(paste("Error: invalid metric", metric))
	}
	distance
}

getPairwiseDistances = function(Matrix, metric, weights, coeffs = NULL, sep = "\\+") {
	# This function computes the pairwise distances between every row of the matrix.
	# NOTE: It is now extended to work with complex strains using a splitting trick (June 6, 2015).
	allComplex = apply(Matrix, 1, function(x) {getComplexPositions(x, sep = sep)})
	allLens = sapply(allComplex, length)
	K = nrow(Matrix)
	N = K + sum(allLens > 0)
	someComplex = (N > K)
	Distances = matrix(0, N, N)
	if (someComplex) {
		complexInds = which(allLens > 0)
		numComplex = N - K
		subMat = Matrix[complexInds,]
		simpleInds = setdiff(1:K, complexInds)
		numSimple = K - numComplex
		Matrix = Matrix[simpleInds,]
		Map = vector("list", K)
		Map[simpleInds]  = split(1:numSimple, 1:numSimple)
		Map[complexInds] = split(numSimple + 1:(2 * numComplex), rep(1:numComplex, each = 2))
		redComplex = allComplex[complexInds]
		redLens = allLens[complexInds]
		allComp = 1:numComplex
		res1 = lapply(allComp,function(x){pickPositions(subMat[x,],redComplex[[x]],rep(1,redLens[x]),sep)})
		res2 = lapply(allComp,function(x){pickPositions(subMat[x,],redComplex[[x]],rep(2,redLens[x]),sep)})
		res  = lapply(allComp,function(x){rbind(res1[[x]], res2[[x]])})
		Matrix = matrix(as.numeric(rbind(Matrix, do.call(rbind, res))), nrow = numSimple + 2 * numComplex)
	}
	# print(paste("There are", N, "strains to be processed"))
	for (ind in 1:(N-1)) {
		follow = (ind + 1):N
		curDists = getDistance(Matrix[ind, ], Matrix[follow, , drop = FALSE], metric = metric, weights = weights, coeffs = coeffs)
		Distances[ind, follow] = curDists
		Distances[follow, ind] = curDists
	}
	if (someComplex) {
		redDistances = matrix(0, K, K)
		for (ind1 in 1:(K-1)) {
			cur1 = Map[[ind1]]
			for (ind2 in (ind1 + 1):K) {
				cur2 = Map[[ind2]]
				meanDists = mean(Distances[cbind(cur1, cur2)])
				redDistances[ind1, ind2] = meanDists
				redDistances[ind2, ind1] = meanDists
			}
		}
		Distances = redDistances
	}
	Distances
}

getCover = function(complexVector, mode, Collection = NULL, sep = "\\+", metric = "linear", weights = NULL, restrict = NULL, matrix = FALSE, attr = FALSE, extraLength = 0) {
	if (extraLength > 0) {
		avoid = nrow(Collection) + 1 - (1:extraLength)
		redCollection = Collection[-avoid, ]
	}
	complexPositions = getComplexPositions(complexVector, sep = sep)
	NAPositions = which(is.na(complexVector))
	q = length(NAPositions)
	curNA0 = rep(NA, q)
	curNA1 = rep(NA, q)
	if (!is.null(restrict) && length(restrict) != length(complexPositions)) {
		stop("Error: misspecified restriction vector!")
	}
	openPositions = which(restrict == 0)
	L = length(openPositions)
	if (L == 0 || (L == 1 && length(complexPositions) == 1)) {
		if (L == 0) {
			Vector0 = pickPositions(complexVector, complexPositions, restrict, sep = sep)
			Vector1 = pickPositions(complexVector, complexPositions, 3 - restrict, sep = sep)
		}
		else { # there is only one possibility
			Vector0 = pickPositions(complexVector, complexPositions, 1, sep = sep)
			Vector1 = pickPositions(complexVector, complexPositions, 2, sep = sep)
		}
		# infer missing positions from the closest vectors in the collection!
		closest0 = getClosestDistance(Vector0, Collection, metric = metric, weights = weights)
		closest1 = getClosestDistance(Vector1, Collection, metric = metric, weights = weights)
		Vector0[NAPositions] = Collection[closest0$index, NAPositions]
		Vector1[NAPositions] = Collection[closest1$index, NAPositions]
		if (matrix) {
			Matrix = matrix(c(closest0$distance, closest1$distance, closest0$distance + closest1$distance), nrow=1)
			rownames(Matrix) = "1"
			if (attr) {
				bests0 = closest0$bests
				bests1 = closest1$bests
				if (extraLength > 0) {
					fullBests = list(list(bests0, bests1))
					bests0 = setdiff(bests0, avoid)
					if (length(bests0) == 0) {
						closest0 = getClosestDistance(Vector0, redCollection, metric = metric, weights = weights)
						bests0 = closest0$bests
					}
					bests1 = setdiff(bests1, avoid)
					if (length(bests1) == 0) {
						closest1 = getClosestDistance(Vector1, redCollection, metric = metric, weights = weights)
						bests1 = closest1$bests
					}
				}
				Bests = list(list(bests0, bests1))
			}
		}
	}
	else {
		curArg = c()
		if (mode == "random") {
			curArg = c(1 + floor(runif(L) * 2))
		}
		else if (mode == "min" || mode == "max"){
			if (mode == "min") {
				curOpt = Inf
			}
			else {
				curOpt = -Inf
			}
			allIndices = hcube(rep(2, L))
			N = nrow(allIndices)
			if (matrix) {
				Matrix = matrix(NA, N, 3)
				rownames(Matrix) = apply(allIndices, 1, function(x) {paste0(x, collapse = ",")})
				if (attr) {
					Bests = vector("list", N)
					if (extraLength > 0) {
						fullBests = vector("list", N)
					}
				}
			}
			for (index in 1:N) {
				firstChoice = restrict
				firstChoice[openPositions] = allIndices[index,]
				Vector0  = pickPositions(complexVector, complexPositions, firstChoice, sep = sep)
				closest0 = getClosestDistance(Vector0, Collection, metric = metric, weights = weights)
				distance0 = closest0$distance
				secondChoice = 3 - firstChoice
				Vector1  = pickPositions(complexVector, complexPositions, secondChoice, sep = sep)
				closest1 = getClosestDistance(Vector1, Collection, metric = metric, weights = weights)
				distance1 = closest1$distance
				sumDist = distance0 + distance1
				if (matrix) {
					Matrix[index,] = c(distance0, distance1, sumDist)
					if (attr) {
						bests0 = closest0$bests
						bests1 = closest1$bests
						if (extraLength > 0) {
							fullBests[[index]] = list(bests0, bests1)
							bests0 = setdiff(bests0, avoid)
							if (length(bests0) == 0) {
								closest0 = getClosestDistance(Vector0, redCollection, metric = metric, weights = weights)
								bests0 = closest0$bests
							}
							bests1 = setdiff(bests1, avoid)
							if (length(bests1) == 0) {
								closest1 = getClosestDistance(Vector1, redCollection, metric = metric, weights = weights)
								bests1 = closest1$bests
							}
						}
						Bests[[index]] = list(bests0, bests1)
					}
				}
				if ((mode == "min" && sumDist < curOpt) || (mode == "max" && sumDist > curOpt)) {
					curArg = firstChoice
					curNA0 = Collection[closest0$index, NAPositions]
					curNA1 = Collection[closest1$index, NAPositions]
					curOpt = sumDist
					if (closest0$index == closest1$index) {
						WARN = TRUE
					}
					else {
						WARN = FALSE
					}
				}
			}
		}
		else {
			stop(paste("Error: invalid mode", mode))
		}
		Vector0 = pickPositions(complexVector, complexPositions,     curArg, sep = sep)
		Vector0[NAPositions] = curNA0
		Vector1 = pickPositions(complexVector, complexPositions, 3 - curArg, sep = sep)
		Vector1[NAPositions] = curNA1
	}
	if (matrix) {
		if (attr) {
			attr(Matrix, "bests") = Bests
			if (extraLength > 0) {
				attr(Matrix, "fullBests") = fullBests
			}
		}
		return(Matrix)
	}
	return(list(Vector0, Vector1))
}

getComplexPositions = function(complexVector, sep) {
	# This function returns the positions in a vector that are complex, indicated by sep
	complexPositions = grep(sep, complexVector)
	complexPositions
}

pickPositions = function(complexVector, complexPositions, choice, sep) {
	pickedVector = complexVector
	preOptions = strsplit(complexVector[complexPositions], sep)
	allLens = sapply(preOptions, length)
	if (any(allLens != 2)) {
		print("Warning: some loci have more than 2 copy numbers; only the first two are kept!")
		preOptions = lapply(preOptions, function(x) {x[1:2]})
	}
	Options = matrix(unlist(preOptions), nrow = 2)
	pickedVector[complexPositions] = Options[cbind(choice, 1:ncol(Options))]
	pickedVector = as.integer(pickedVector)
	pickedVector
}

getClosestDistance = function(Vector, Collection, metric, weights, coeffs = NULL) {
	# This function computes the closest row vector in a given collection to a given vector using the given metric.
	distance = getDistance(Vector, Collection, metric = metric, weights = weights, coeffs = coeffs)
	optimum = min(distance)
	best = which.min(distance)
	bests = which(distance - optimum < eps)
	output = list(distance = optimum, index = best, bests = bests)
	output
}

expandOptions = function(List) { # transforms from CNF to DNF form
	L = length(List)
	allOptions = hcube(sapply(List, length))
	Matrix = apply(allOptions, 1, function(x) {Vec = rep(NA, L); for (ind in 1:L) {Vec[ind] = List[[ind]][x[ind]]}; Vec})
	Matrix = t(Matrix)
	Matrix
}

getIntersectionLists = function(List1, List2, specialNA = FALSE) {
	if (specialNA) {
		Res = mapply(function(x,y) 
		{if (any(is.na(x)) || any(is.na(y))) {as.vector(na.omit(union(x,y)))} else {intersect(x,y)}}, 
		List1, List2)
	}
	else {
		Res = mapply(function(x,y) {intersect(x,y)}, List1, List2)
	}
	Res
}

getIntersection = function(strain1, strain2, sep = "\\+", specialNA = FALSE) {
	List1 = lapply(strsplit(strain1, sep), as.numeric)
	List2 = lapply(strsplit(strain2, sep), as.numeric)
	Res12 = getIntersectionLists(List1, List2, specialNA = specialNA)
	Res12
}

imputeTwoPairs = function(strain11, strain12, strain21, strain22, sep, conservative = FALSE, default = NA){
	dataMatrix = rbind(strain11, strain12, strain21, strain22)
	allNACols = which(colSums(is.na(dataMatrix)) == 4)
	if (conservative) {
		dataMatrix[ , allNACols] = rep(default[allNACols], each = 4)
	}
	complexPositions1 = union(getComplexPositions(strain11, sep), getComplexPositions(strain12, sep))
	complexPositions2 = union(getComplexPositions(strain21, sep), getComplexPositions(strain22, sep))
	complexPositions  = union(complexPositions1, complexPositions2)
	goodPositions = setdiff(1:length(strain11), complexPositions)
	dataMatrixR = matrix(as.numeric(dataMatrix[, goodPositions]), nrow = 4)
	Tree = createEmptyTree(numNodes = 7)
	Tree = createNewNode(Tree, children = c(1,2), branchLengths = c(1,1), 5)
	Tree = createNewNode(Tree, children = c(3,4), branchLengths = c(1,1), 6)
	Tree = createNewNode(Tree, children = c(5,6), branchLengths = c(1,1), 7)
	reconstruction = intervalReconstruct(dataMatrixR, Tree)
	lowerBounds = reconstruction[1:4,   1:ncol(dataMatrixR)]
	upperBounds = reconstruction[1:4, -(1:ncol(dataMatrixR))]
	if (any(lowerBounds != upperBounds)) {
		print("Warning: could not impute unambiguously!")
		print(reconstruction)
	}
	dataMatrix[ , goodPositions] = lowerBounds
	# Also impute the positions whose complexity is limited to only one subtree
	complexPositions2Not1 = setdiff(complexPositions2, complexPositions1)
	submatrix1 = dataMatrix[1:2, complexPositions2Not1, drop = FALSE]
	if (!(any(submatrix1[1,] != submatrix1[2,], na.rm = TRUE))) {
		submatrix1[1, is.na(submatrix1[1,])] = submatrix1[2, is.na(submatrix1[1,])]
		submatrix1[2, is.na(submatrix1[2,])] = submatrix1[1, is.na(submatrix1[2,])]
		dataMatrix[1:2, complexPositions2Not1] = submatrix1
	}
	complexPositions1Not2 = setdiff(complexPositions1, complexPositions2)
	submatrix2 = dataMatrix[3:4, complexPositions1Not2, drop = FALSE]
	if (!(any(submatrix2[1,] != submatrix2[2,], na.rm = TRUE))) {
		submatrix2[1, is.na(submatrix2[1,])] = submatrix2[2, is.na(submatrix2[1,])]
		submatrix2[2, is.na(submatrix2[2,])] = submatrix2[1, is.na(submatrix2[2,])]
		dataMatrix[3:4, complexPositions1Not2] = submatrix2
	}
	output = list(dataMatrix, complexPositions)
	output
}

testCompatible = function(restrict1, restrict2) {
	# Compatibility means not having a match at one position and a mismatch at another
	if (length(restrict1) != length(restrict2)) {
		print("The lengths of the restriction vectors differ")
		return(FALSE)
	}
	pos12 = which((restrict1 == 1 & restrict2 == 2) | (restrict1 == 2 & restrict2 == 1))
	pos21 = which((restrict1 == 1 & restrict2 == 1) | (restrict1 == 2 & restrict2 == 2))
	return(length(pos12) * length(pos21) == 0)
}

findPairwiseCover = function(strain1, strain2, sep, conservative = FALSE, default = NA) {
	### NOTE: This will fail if two strains have complex loci with partial as well as complete overlaps!
	L = length(strain1)
	stopifnot(length(strain2) == L)
	missing1  = which(is.na(strain1))
	missing2  = which(is.na(strain2))
	missing12 = intersect(missing1, missing2)
	complex1  = getComplexPositions(strain1, sep = sep)
	complex2  = getComplexPositions(strain2, sep = sep)
	complex12 = intersect(complex1, complex2)
	complexStatus1 = (length(complex1) > 0)
	if (complexStatus1) {
		Options1  = matrix(unlist(strsplit(strain1[complex1], sep)), nrow = 2)
		restrict1 = rep(0, ncol(Options1))
	}
	else {
		restrict1 = NULL
	}
	complexStatus2 = (length(complex2) > 0)
	if (complexStatus2) {
		Options2  = matrix(unlist(strsplit(strain2[complex2], sep)), nrow = 2)
		restrict2 = rep(0, ncol(Options2))
	}
	else {
		restrict2 = NULL
	}
	if (length(missing12) == L) {
		output = list(list(), complexStatus1, complexStatus2)
		return(output)
	}
	else if (length(missing1) == L) {
		strain2[missing2] = default[missing2]
		attr(strain2, "restrict") = restrict2
		output = list(list(strain2), complexStatus1, complexStatus2)
		return(output)
	}
	else if (length(missing2) == L) {
		strain1[missing1] = default[missing1]
		attr(strain1, "restrict") = restrict1
		output = list(list(strain1), complexStatus1, complexStatus2)
		return(output)
	}
	allComplex = union(complex1, complex2)
	restrictC  = rep(0, length(allComplex))
	if (!(complexStatus1 || complexStatus2)) {
		numDiff = sum(strain1 != strain2, na.rm = TRUE)
		if ((!is.na(numDiff)) && numDiff > 0) { # new special case!
			# print("Special case: two unequal simple strains")
			restrictS = rep(0, numDiff)
			indS = 1
		}
	}
	indC = 1
	compatible = TRUE
	ind1 = 1
	ind2 = 1
	allUnions = vector("list", L)
	for (ind in 1:L) {
		if (ind %in% missing12) {
			allUnions[[ind]] = NA
		}
		else if (ind %in% missing1) {
			if (ind %in% complex2) {
				allUnions[[ind]] = Options2[,ind2]
				ind2 = ind2 + 1
				indC = indC + 1
			}
			else {
				allUnions[[ind]] = strain2[ind]
			}
		}
		else if (ind %in% missing2) {
			if (ind %in% complex1) {
				allUnions[[ind]] = Options1[,ind1]
				ind1 = ind1 + 1
				indC = indC + 1
			}
			else {
				allUnions[[ind]] = strain1[ind]
			}
		}
		else {
			if (ind %in% complex12) {
				opts1 = Options1[,ind1]
				opts2 = Options2[,ind2]
				intersection = intersect(opts1, opts2)
				if (length(intersection) == 1) {
					warning("Found two complex loci intersecting in one value! Processing strains separately.")
					compatible = FALSE
					break
				}
				else if (length(intersection) == 0) {
					compatible = FALSE
					break
				}
				allUnions[[ind]] = union(opts1, opts2)
				ind1 = ind1 + 1
				ind2 = ind2 + 1
				indC = indC + 1
			}
			else if (ind %in% complex1) {
				opts1 = Options1[,ind1]
				if (strain2[ind] %in% opts1) {
					allUnions[[ind]] = opts1
					restrictC[indC] = which(opts1 != strain2[ind])
				}
				else {
					compatible = FALSE
					break
				}
				ind1 = ind1 + 1
				indC = indC + 1
			}
			else if (ind %in% complex2) {
				opts2 = Options2[,ind2]
				if (strain1[ind] %in% opts2) {
					allUnions[[ind]] = opts2
					restrictC[indC] = which(opts2 == strain1[ind])
				}
				else {
					compatible = FALSE
					break
				}
				ind2 = ind2 + 1
				indC = indC + 1
			}
			else {
				if (strain1[ind] == strain2[ind]) {
					allUnions[[ind]] = strain1[ind]
				}
				else {
					if (!(complexStatus1 || complexStatus2)) {
						allUnions[[ind]] = sort(c(strain1[ind], strain2[ind]))
						if (strain1[ind] == allUnions[[ind]][1]) {
							restrictS[indS] = 1
						}
						else {
							restrictS[indS] = 2
						}
						indS = indS + 1
					}
					else {
						compatible = FALSE
						break
					}
				}
			}
		}
	}
	if (compatible) { # use allUnions for the combined strain
		mutualCover = unlist(lapply(allUnions, pasteFunction))
		if (conservative) {
			missingC = is.na(mutualCover)
			mutualCover[missingC] = default[missingC]
		}
		if (!(complexStatus1 || complexStatus2)) {
			if ((!is.na(numDiff)) && numDiff > 0) {
				attr(mutualCover, "restrict") = restrictS
			}
		}	
		else {
			attr(mutualCover, "restrict") = restrictC
		}
		outCovers = list(mutualCover)
	}
	else { # cover each strain separately
		if (conservative) {
			strain1[missing1] = default[missing1]
			strain2[missing2] = default[missing2]
		}
		attr(strain1, "restrict") = restrict1
		attr(strain2, "restrict") = restrict2
		outCovers = list(strain1, strain2)
	}
	output = list(outCovers, complexStatus1, complexStatus2)
	output
}

pasteFunction = function(x, sep = "\\+") {
	if (all(is.na(x))) {
		x
	} 
	else {
		paste(x, collapse = substr(sep, 2, 2))
	}
}

findMode = function(x) {
	Tab = table(x)
	mode = as.numeric(names(Tab)[which.max(Tab)])
	mode
}

findMostCommonLoci = function(Table) {
	loci = apply(Table, 2, findMode)
	loci
}

rankByDistance = function(distanceMatrix) {
	rankMatrix = apply(distanceMatrix, 1, function(x) {rank(x, ties.method = "min") - 1})
	rankMatrix
}

processDistance = function(file) {
	Table = read.table(file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, header = TRUE, row.names = NULL)
	Matrix = as.matrix(Table[,-1])
	rankMatrix = rankByDistance(Matrix)
	complexTable = makeComplexTable(colnames(Table)[-1])
	L = length(complexTable)
	Tab = matrix(NA, L, 2)
	rownames(Tab) = names(complexTable)
	colnames(Tab) = c("meanRank", "size")
	for (ind in 1:L) {
		curInds = complexTable[[ind]]
		curMat = rankMatrix[curInds, curInds]
		Q = length(curInds)
		Tab[ind,] = c(sum(curMat)/(Q*(Q-1)), Q)
	}
	Tab
}

processDistances = function(optionsList, baseline = "BeastMIRU_Latest", extension = ".distances", special = 3) {
	FUN = function(x, y) {outer(x, y, paste0)}
	L = length(optionsList)
	fullOptions = optionsList[[L]]
	for (ind in (L-1):1) {
		fullOptions = as.vector(do.call(FUN, list(optionsList[[ind]], fullOptions)))
	}
	Q = length(fullOptions)
	for (opt in 1:Q) {
		curOpt = fullOptions[opt]
		print(curOpt)
		curFile = paste0(baseline, curOpt, extension)
		curTab = processDistance(curFile)
		if (opt == 1) {
			Mat = matrix(NA, nrow(curTab), Q + 1)
			rownames(Mat) = rownames(curTab)
			colnames(Mat) = c("size", fullOptions)
			Mat[, 1:2] = curTab[, 2:1]
		}
		else {
			Mat[, opt + 1] = curTab[, 1]
		}
	}
	outFilename = paste0(baseline, optionsList[[special]][1], extension)
	write.table(Mat, file = outFilename, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
}

testLocalIndependence = function(Table, sep = "\\+") {
	L = ncol(Table)
	pVals = rep(NA, L)
	names(pVals) = colnames(Table)
	for (ind in 1:L) {
		curCol = Table[,ind]
		redCol = curCol[grep(sep, curCol)]
		if (length(redCol) > 0) {
			splitCol = matrix(unlist(strsplit(redCol, sep)), ncol = 2, byrow = TRUE)
			splitCol = matrix(as.numeric(splitCol), ncol = 2)
			splitCol = rbind(splitCol, splitCol[,2:1])
			curTab = table(splitCol[,1], splitCol[,2])
			if (length(curTab) >= 2) {
				print(curTab)
				pVals[ind] = chisq.test(curTab, simulate.p.value = TRUE, B = 20000)$p.value
			}
		}
	}
	pVals
}