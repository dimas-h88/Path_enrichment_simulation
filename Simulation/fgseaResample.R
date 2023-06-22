library(BiocParallel)
calcGseaStatShuffles <- function(stats,
                                 selectedStats,
                                 reranks, 
                                 gseaParam = 1,
                                 scoreType = c("std", "pos", "neg")) {
    scoreType <- match.arg(scoreType)
    
    statsAdj <- abs(stats)^gseaParam
    
    N <- length(stats)
    K <- length(selectedStats)
    SH <- length(reranks)
    
    X <- matrix(-1/(N-K), nrow=SH, ncol = N)
    NS <- sum(statsAdj[selectedStats])
    for (i in seq_len(SH)) {
        X[i, reranks[[i]][selectedStats]] <- statsAdj[selectedStats]/NS
    }
    X1 <- colMeans(X)
    
    ES <- c(0, cumsum(X1))
    minP <- min(ES)
    maxP <- max(ES)
    
    res <- switch(scoreType,
                  std = ifelse(maxP == -minP, 0, ifelse(maxP > -minP, maxP, minP)),
                  pos = maxP,
                  neg = minP)
    list(res=res)
}


fgseaMultilevelSE <- function(pathways,
                              stats,
                              se,
                              sampleSize  = 101,
                              minSize     = 1,
                              maxSize     = length(stats)-1,
                              eps         = 1e-50,
                              scoreType   = c("std", "pos", "neg"),
                              nproc       = 0,
                              gseaParam   = 1,
                              BPPARAM     = NULL,
                              nPermSimple = 1000,
                              nResample   = 11,
                              absEps      = NULL)
{
    scoreType <- match.arg(scoreType)
    pp <- fgsea:::preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, scoreType)
    pathwaysFiltered <- pp$filtered
    pathwaysSizes <- pp$sizes
    
    reranks <- replicate(nResample, {
        stats1 <- stats + rnorm(length(stats), sd=se)
        stats1 <- stats1[names(pp$stats)]
        stats1Ord <- rank(-stats1, ties.method = "random")
        stats1Ord
    }, simplify=FALSE)
    
    
    stats <- pp$stats
    m <- length(pathwaysFiltered)
    
    if (m == 0) {
        return(data.table(pathway=character(),
                          pval=numeric(),
                          padj=numeric(),
                          log2err=numeric(),
                          ES=numeric(),
                          NES=numeric(),
                          size=integer(),
                          leadingEdge=list()))
    }
    
    
    
    # Warning message for deprecated absEps parameter
    if (!is.null(absEps)){
        warning("You are using deprecated argument `absEps`. ",
                "Use `eps` argument instead. ",
                "`absEps` was assigned to `eps`.")
        eps <-  absEps
    }
    
    # Warning message for to small value for sampleSize
    if (sampleSize < 3){
        warning("sampleSize is too small, so sampleSize = 3 is set.")
        sampleSize <- max(3, sampleSize)
    }
    
    #To avoid warnings during the check
    log2err=nMoreExtreme=pathway=pval=padj=NULL
    nLeZero=nGeZero=leZeroMean=geZeroMean=nLeEs=nGeEs=isCpGeHalf=NULL
    ES=NES=size=leadingEdge=NULL
    .="damn notes"
    
    minSize <- max(minSize, 1)
    eps <- max(0, min(1, eps))
    
    
    if (sampleSize %% 2 == 0){
        sampleSize <-  sampleSize + 1
    }
    
    
    gseaStatRes <- do.call(rbind,
                           lapply(pathwaysFiltered,
                                  calcGseaStat,
                                  stats             = stats,
                                  returnLeadingEdge = TRUE,
                                  scoreType         = scoreType))
    
    gseaStatResShuffles <- do.call(rbind,
                                   lapply(pathwaysFiltered,
                                          calcGseaStatShuffles,
                                          stats             = stats,
                                          reranks           = reranks,
                                          scoreType         = scoreType))
    
    leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
    pathwayScores <- unlist(gseaStatResShuffles[, "res"])
    
    
    seeds <- sample.int(10^9, 1)
    BPPARAM <- fgsea:::setUpBPPARAM(nproc=nproc, BPPARAM=BPPARAM)
    
    simpleFgseaRes <- fgsea:::fgseaSimpleImpl(pathwayScores=pathwayScores, pathwaysSizes=pathwaysSizes,
                                              pathwaysFiltered=pathwaysFiltered, leadingEdges=leadingEdges,
                                              permPerProc=nPermSimple, seeds=seeds, toKeepLength=m,
                                              stats=stats, BPPARAM=SerialParam(), scoreType=scoreType)
    
    switch(scoreType,
           std = simpleFgseaRes[, modeFraction := ifelse(ES >= 0, nGeZero, nLeZero)],
           pos = simpleFgseaRes[, modeFraction := nGeZero],
           neg = simpleFgseaRes[, modeFraction := nLeZero])
    
    
    simpleFgseaRes[, leZeroMean := NULL]
    simpleFgseaRes[, geZeroMean := NULL]
    simpleFgseaRes[, nLeEs := NULL]
    simpleFgseaRes[, nGeEs := NULL]
    simpleFgseaRes[, nLeZero := NULL]
    simpleFgseaRes[, nGeZero := NULL]
    
    simpleFgseaRes[modeFraction < 10, pval := as.numeric(NA)]
    simpleFgseaRes[modeFraction < 10, padj := as.numeric(NA)]
    simpleFgseaRes[modeFraction < 10, NES := as.numeric(NA)]
    
    if (any(simpleFgseaRes$modeFraction < 10)){
        warning("There were ",
                paste(sum(simpleFgseaRes$modeFraction < 10)),
                " pathways for which P-values were not calculated properly due to ",
                "unbalanced (positive and negative) gene-level statistic values. ",
                "For such pathways pval, padj, NES, log2err are set to NA. ",
                "You can try to increase the value of the argument nPermSimple (for example set it nPermSimple = ",
                paste0(format(nPermSimple * 10, scientific = FALSE), ")"))
    }
    
    # Storing NA fgseaSimple results in a separate data.table
    naSimpleRes <- simpleFgseaRes[is.na(pval)]
    naSimpleRes[, padj := as.numeric(NA)]
    naSimpleRes[, log2err := as.numeric(NA)]
    naSimpleRes[, modeFraction := NULL]
    
    simpleFgseaRes <- simpleFgseaRes[!is.na(pval)]
    
    leftBorder <- log2(qbeta(0.025,
                             shape1 = simpleFgseaRes$nMoreExtreme,
                             shape2=nPermSimple - simpleFgseaRes$nMoreExtreme + 1))
    
    rightBorder <- log2(qbeta(1 - 0.025,
                              shape1 = simpleFgseaRes$nMoreExtreme + 1,
                              shape2=nPermSimple - simpleFgseaRes$nMoreExtreme))
    
    crudeEstimator <- log2((simpleFgseaRes$nMoreExtreme + 1) / (nPermSimple + 1))
    
    simpleError <- 0.5 * pmax(crudeEstimator - leftBorder, rightBorder - crudeEstimator)
    
    
    multError <- sapply((simpleFgseaRes$nMoreExtreme + 1) / (nPermSimple + 1), multilevelError, sampleSize)
    
    
    if (all(multError >= simpleError)){
        simpleFgseaRes[, log2err := 1/log(2)*sqrt(trigamma(nMoreExtreme + 1) - trigamma((nPermSimple + 1)))]
        simpleFgseaRes[, modeFraction := NULL]
        simpleFgseaRes <- rbindlist(list(simpleFgseaRes, naSimpleRes), use.names = TRUE)
        
        setorder(simpleFgseaRes, pathway)
        simpleFgseaRes[, "nMoreExtreme" := NULL]
        setcolorder(simpleFgseaRes, c("pathway", "pval", "padj", "log2err",
                                      "ES", "NES", "size", "leadingEdge"))
        
        simpleFgseaRes <- simpleFgseaRes[]
        return(simpleFgseaRes)
    }
    
    
    dtSimpleFgsea <- simpleFgseaRes[multError >= simpleError]
    dtSimpleFgsea[, log2err := 1/log(2)*sqrt(trigamma(nMoreExtreme + 1) - trigamma(nPermSimple + 1))]
    dtSimpleFgsea[, modeFraction := NULL]
    
    
    
    dtMultilevel <- simpleFgseaRes[multError < simpleError]
    # Probability estimation in the denominator of the multilevel algortihm
    # (s_r(q) >= 0 in the notation of the article):
    dtMultilevel[, "denomProb" := (modeFraction + 1) / (nPermSimple + 1)]
    
    multilevelPathwaysList <- split(dtMultilevel, by="size")
    # In most cases, this gives a speed increase with parallel launches.
    indxs <- sample(1:length(multilevelPathwaysList))
    multilevelPathwaysList <- multilevelPathwaysList[indxs]
    
    seed=sample.int(1e9, size=1)
    
    sign <- if (scoreType %in% c("pos", "neg")) TRUE else FALSE
    cpp.res <- fgsea:::multilevelImpl(multilevelPathwaysList, stats,
                                      sampleSize, seed, eps, sign=sign,
                                      BPPARAM=BPPARAM)
    cpp.res <- rbindlist(cpp.res)
    
    
    result <- rbindlist(multilevelPathwaysList)
    
    # `cppMpval` - P-values that are computed in cpp code
    # `isCpGeHalf` is a flag that mathces: whether the conditional probability
    # is greater than or equal to 0.5 (see article for details)
    result[, pval := pmin(1, cpp.res$cppMPval / denomProb)]
    result[, isCpGeHalf := cpp.res$cppIsCpGeHalf]
    result[, log2err := multilevelError(pval, sampleSize = sampleSize)]
    result[isCpGeHalf == FALSE, log2err:= NA]
    
    if (!all(result$isCpGeHalf)){
        warning("For some of the pathways the P-values were likely overestimated. ",
                "For such pathways log2err is set to NA.")
    }
    
    result[, isCpGeHalf := NULL]
    result[, modeFraction := NULL]
    result[, denomProb := NULL]
    
    
    result <- rbindlist(list(result, dtSimpleFgsea, naSimpleRes), use.names = TRUE)
    result[, nMoreExtreme := NULL]
    
    result[pval < eps, c("pval", "log2err") := list(eps, NA)]
    result[, padj := p.adjust(pval, method = "BH")]
    
    if (nrow(result[pval==eps & is.na(log2err)])){
        warning("For some pathways, in reality P-values are less than ",
                paste(eps),
                ". You can set the `eps` argument to zero for better estimation.")
    }
    
    setcolorder(result, c("pathway", "pval", "padj", "log2err",
                          "ES", "NES", "size", "leadingEdge"))
    
    setorder(result, pathway)
    
    result <- result[]
    result
}
