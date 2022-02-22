# JCI-FCI

###############################################################################
## JCI-FCI
###############################################################################

JCI_fci <- function(suffStat, indepTest, alpha, labels, p,
                  skel.method = c("stable", "original", "stable.fast"),
                  type = c("normal", "anytime", "adaptive"),
                  fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                  m.max = Inf, pdsep.max = Inf, rules = rep(TRUE, 10),
                  doPdsep = TRUE, biCC = FALSE, conservative = FALSE,
                  maj.rule = FALSE, numCores = 1, verbose = FALSE, F_tester=NULL, F_node_names = list(),orientFnode=TRUE)
{
  ## Purpose: Perform FCI-Algorithm, i.e., estimate PAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat, indepTest: info for the tests
  ## - p: number of nodes in the graph
  ## - alpha: Significance level of individual partial correlation tests
  ## - verbose: 0 - no output, 1 - detailed output
  ## - fixedGaps: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical); gaps fixed here are not changed
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximum size of conditioning set
  ## - pdsep.max: maximaum size of conditioning set for Possible-D-SEP
  ## - rules: array of length 10 wich contains TRUE or FALSE corresponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - doPdsep: compute possible dsep
  ## - biCC: TRUE or FALSE variable containing if biconnected components are
  ##         used to compute pdsep
  ## - conservative: TRUE or FALSE defining if
  ##          the v-structures after the pdsep
  ##          have to be oriented conservatively or nor
  ## - maj.rule: TRUE or FALSE variable containing if the majority rule is
  ##             used instead of the normal conservative
  ## - labels: names of the variables or nodes
  ## - type: it specifies the version of the FCI that has to be used.
  ##         Per default it is normal, the normal FCI algorithm. It can also be
  ##         anytime for the Anytime FCI and in this cas m.max must be specified;
  ##         or it can be adaptive for Adaptive Anytime FCI and in this case
  ##         m.max must not be specified.
  ## - numCores: handed to skeleton(), used for parallelization
  
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Dec 2009; update: Diego Colombo, 2012; Martin Maechler, 2013
  
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  
  ## Check that the type is a valid one
  type <- match.arg(type)
  if (type == "anytime" && m.max == Inf)
    stop("To use the Anytime FCI you must specify a finite 'm.max'.")
  if (type == "adaptive" && m.max != Inf)
    stop("To use the Adaptive Anytime FCI you must not specify 'm.max'.")
  
  if (conservative && maj.rule)
    stop("Choose either conservative FCI or majority rule FCI")
  
  cl <- match.call()
  if (verbose) cat("Compute Skeleton\n================\n")
  
  skel <- JCI_skeleton(suffStat, indepTest, alpha, labels = labels, method = skel.method,
                     fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                     NAdelete=NAdelete, m.max=m.max, numCores=numCores, verbose=verbose, F_node_names=F_node_names,F_tester=F_tester)
  skel@call <- cl # so that makes it into result
  G <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  pMax <- skel@pMax
  n.edgetestsSKEL <- skel@n.edgetests
  max.ordSKEL <- skel@max.ord
  allPdsep <- NA
  tripleList <- NULL
  
  if (doPdsep) {
    if (verbose) cat("\nCompute PDSEP\n=============\n")
    pc.ci <- pc.cons.intern(skel, suffStat, indepTest,
                            alpha = alpha, version.unf = c(1,1),
                            maj.rule = FALSE, verbose = verbose)
    ## Recompute (sepsets, G, ...):
    pdsepRes <- pdsep(skel@graph, suffStat, indepTest = indepTest, p = p,
                      sepset = pc.ci$sk@sepset, alpha = alpha, pMax = pMax,
                      m.max = if (type == "adaptive") max.ordSKEL else m.max,
                      pdsep.max = pdsep.max, NAdelete = NAdelete,
                      unfVect = pc.ci$unfTripl, # "tripleList.pdsep"
                      biCC = biCC, verbose = verbose)
    
    ## update the graph & sepset :
    G <- pdsepRes$G
    sepset <- pdsepRes$sepset
    pMax <- pdsepRes$pMax
    allPdsep <- pdsepRes$allPdsep
    n.edgetestsPD <- pdsepRes$n.edgetests
    max.ordPD <- pdsepRes$max.ord
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      tmp.pdsep <- new("pcAlgo", graph = as(G, "graphNEL"), call = cl,
                       n = integer(0), max.ord = as.integer(max.ordSKEL),
                       n.edgetests = n.edgetestsSKEL, sepset = sepset,
                       pMax = pMax, zMin = matrix(NA, 1, 1))
      sk. <- pc.cons.intern(tmp.pdsep, suffStat, indepTest, alpha,
                            verbose = verbose, version.unf = c(1, 1),
                            maj.rule = maj.rule)
      tripleList <- sk.$unfTripl
      ## update the sepsets
      sepset <- sk.$sk@sepset
    }
  }
  else {## !doPdsep : "do not Pdsep"
    n.edgetestsPD <- 0
    max.ordPD <- 0
    allPdsep <- vector("list", p)
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      nopdsep <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                                verbose = verbose, version.unf = c(2, 1),
                                maj.rule = maj.rule)
      tripleList <- nopdsep$unfTripl
      ## update the sepsets
      sepset <- nopdsep$sk@sepset
    }
  }
  if (verbose)
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules),
        "\nCompute collider:\n")
  
  #res <- udag2pag(pag = G, sepset, rules = rules, unfVect = tripleList,
  #               verbose = verbose)
  res <- udag2pag_JCI(pag = G, sepset, rules = rules, unfVect = tripleList,verbose = verbose,F_node_names= F_node_names, orientFnode = orientFnode)
  
  colnames(res) <- rownames(res) <- labels
  new("fciAlgo", amat = res, call = cl, n = integer(0),
      max.ord = as.integer(max.ordSKEL),
      max.ordPDSEP = as.integer(max.ordPD),
      n.edgetests = n.edgetestsSKEL, n.edgetestsPDSEP = n.edgetestsPD,
      sepset = sepset, pMax = pMax, allPdsep = allPdsep)
} ## {JCI-fci}


udag2pag_JCI <- function(pag, sepset, rules = rep(TRUE,10), unfVect = NULL, verbose = FALSE, orientCollider = TRUE,F_node_names=list(),orientFnode=TRUE)
{
  ## Purpose: Transform the Skeleton of a pcAlgo-object to a PAG using
  ## the rules of Zhang. The output is an adjacency matrix.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - pag: adjacency matrix of size pxp
  ## - sepset: list of all separation sets
  ## - rules: array of length 10 wich contains TRUE or FALSE corrsponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - unfVect: Vector with ambiguous triples (coded as number using triple2numb)
  ## - verbose: 0 - no output, 1 - detailed output
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 6 Mar 2009; cleanup: Martin Maechler, 2010
  ## update: Diego Colombo, 2012
  
  ## Notation:
  ## ----------------------------------------------------------------------
  ## 0: no edge
  ## 1: -o
  ## 2: -> (arrowhead)
  ## 3: - (tail)
  ## a=alpha
  ## b=beta
  ## c=gamma
  ## d=theta
  
  stopifnot(is.logical(rules), length(rules) == 10)
  if(!is.numeric(pag))  storage.mode(pag) <- "numeric"
  
  if (any(pag != 0)) {
    p <- as.numeric(dim(pag)[1])
    
    ## orient collider
    if (orientCollider) {
      ind <- which(pag == 1, arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        x <- ind[i, 1]
        y <- ind[i, 2]
        allZ <- setdiff(which(pag[y, ] != 0), x)
        for (z in allZ) {
          if (pag[x, z] == 0 && !((y %in% sepset[[x]][[z]]) ||
                                  (y %in% sepset[[z]][[x]]))) {
            if (length(unfVect) == 0) {
              if (verbose) {
                cat("\n", x, "*->", y, "<-*", z, "\n")
                cat("Sxz=", sepset[[z]][[x]], "and",
                    "Szx=", sepset[[x]][[z]], "\n")
              }
              pag[x, y] <- pag[z, y] <- 2
            }
            else {
              if (!any(unfVect == triple2numb(p, x, y, z), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, z, y, x), na.rm = TRUE)) {
                if (verbose) {
                  cat("\n", x, "*->", y, "<-*", z, "\n")
                  cat("Sxz=", sepset[[z]][[x]], "and",
                      "Szx=", sepset[[x]][[z]], "\n")
                }
                pag[x, y] <- pag[z, y] <- 2
              }
            }
          }
        }
      }
    } ## end: Orient collider
    ids_of_F_nodes= (p-length(F_node_names)+1):p # the last nodes are F nodes by assumption
    ids_non_F_nodes= 1:(p-length(F_node_names))
    #Orient all arrows out of F nodes
    if (orientFnode) {
      ind <- which(pag != 0, arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        x <- ind[i, 1]
        y <- ind[i, 2]
        if ((x %in% ids_of_F_nodes) & (y %in% ids_non_F_nodes)) {
          pag[x, y] <- 2
          pag[y, x] <- 3
        } 
        if ((x %in% ids_non_F_nodes) & (y %in% ids_of_F_nodes)) {
          pag[x, y] <- 3
          pag[y, x] <- 2
        }
        if ((x %in% ids_of_F_nodes) & (y %in% ids_of_F_nodes)) {
          pag[x, y] <- 2
          pag[y, x] <- 2
        }        
      }
    }
    
    old_pag1 <- matrix(0, p, p)
    while (any(old_pag1 != pag)) {
      old_pag1 <- pag
      ##-- R1 ------------------------------------------------------------------
      if (rules[1]) {
        ind <- which((pag == 2 & t(pag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pag[b, ] != 0 & pag[, b] == 1) & (pag[a, ] == 0 & pag[, a] == 0))
          indC <- setdiff(indC, a)
          if (length(indC) > 0) {
            if (length(unfVect) == 0) {
              pag[b, indC] <- 2
              pag[indC, b] <- 3
              if (verbose)
                cat("\nRule 1",
                    "\nOrient:", a, "*->", b, "o-*", indC,
                    "as:", b, "->", indC, "\n")
            }
            else {
              for (c in indC) {
                if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) &&
                    !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                  pag[b, c] <- 2
                  pag[c, b] <- 3
                  if (verbose)
                    cat("\nRule 1",
                        "\nConservatively orient:", a, "*->", b, "o-*",
                        c, "as:", b, "->", c, "\n")
                }
              } ## for( c )
            }
          }
        } ## for( i )
      }
      ##-- R2 ------------------------------------------------------------------
      if (rules[2]) {
        ind <- which((pag == 1 & t(pag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          c <- ind[i, 2]
          indB <- which((pag[a, ] == 2 & pag[, a] == 3 & pag[c, ] != 0 & pag[, c] == 2) | (pag[a, ] == 2 & pag[, a] != 0 & pag[c, ] == 3 & pag[, c] == 2))
          if (length(indB) > 0) {
            pag[a, c] <- 2
            if (verbose) {
              cat("\nRule 2","\n")
              cat("Orient:", a, "->", indB, "*->", c, "or", a, "*->", indB, "->", c, "with", a, "*-o", c, "as:", a, "*->", c, "\n")
            }
          }
        }
      }
      ##-- R3 ------------------------------------------------------------------
      if (rules[3]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          d <- ind[i, 2]
          indAC <- which((pag[b, ] != 0 & pag[, b] == 2) & (pag[, d] == 1 & pag[d, ] != 0))
          if (length(indAC) >= 2) {
            if (length(unfVect) == 0) {
              counter <- 0
              while ((counter < (length(indAC) - 1)) && (pag[d, b] != 2)) {
                counter <- counter + 1
                ii <- counter
                while (ii < length(indAC) && pag[d, b] != 2) {
                  ii <- ii + 1
                  if (pag[indAC[counter], indAC[ii]] == 0 && pag[indAC[ii], indAC[counter]] == 0) {
                    if (verbose) {
                      cat("\nRule 3","\n")
                      cat("Orient:", d, "*->", b, "\n")
                    }
                    pag[d, b] <- 2
                  }
                }
              }
            }
            else {
              comb.indAC <- combn(indAC, 2)
              for (j in 1:dim(comb.indAC)[2]) {
                a <- comb.indAC[1, j]
                c <- comb.indAC[2, j]
                if (pag[a, c] == 0 && pag[c, a] == 0 && c != a) {
                  if (!any(unfVect == triple2numb(p, a, d, c), na.rm = TRUE) &&
                      !any(unfVect == triple2numb(p, c, d, a), na.rm = TRUE)) {
                    pag[d, b] <- 2
                    if (verbose) {
                      cat("\nRule 3","\n")
                      cat("Conservatively orient:", d, "*->", b, "\n")
                    }
                  }
                }
              }
            }
          }
        }
      }
      ##-- R4 ------------------------------------------------------------------
      if (rules[4]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)## b o-* c
        while (length(ind) > 0) {
          b <- ind[1, 1]
          c <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all a s.t. a -> c and a <-* b
          indA <- which((pag[b, ] == 2 & pag[, b] != 0) &
                          (pag[c, ] == 3 & pag[, c] == 2))
          ## chose one a s.t. the initial triangle structure exists and the edge hasn't been oriented yet
          while (length(indA) > 0 && pag[c,b] == 1) {
            a <- indA[1]
            indA <- indA[-1]
            ## path is the initial triangle
            ## abc <- c(a, b, c)
            ## Done is TRUE if either we found a minimal path or no path exists for this triangle
            Done <- FALSE
            ### MM: FIXME?? Isn't  Done  set to TRUE in *any* case inside the following
            ### while(.), the very first time already ??????????
            while (!Done && pag[a,b] != 0 && pag[a,c] != 0 && pag[b,c] != 0) {
              ## find a minimal discriminating path for a,b,c
              md.path <- minDiscrPath(pag, a,b,c, verbose = verbose)
              ## if the path doesn't exists, we are done with this triangle
              if ((N.md <- length(md.path)) == 1) {
                Done <- TRUE
              }
              else {
                ## a path exists
                ## if b is in sepset
                if ((b %in% sepset[[md.path[1]]][[md.path[N.md]]]) ||
                    (b %in% sepset[[md.path[N.md]]][[md.path[1]]])) {
                  if (verbose)
                    cat("\nRule 4",
                        "\nThere is a discriminating path between",
                        md.path[1], "and", c, "for", b, ",and", b, "is in Sepset of",
                        c, "and", md.path[1], ". Orient:", b, "->", c, "\n")
                  pag[b, c] <- 2
                  pag[c, b] <- 3
                }
                else {
                  ## if b is not in sepset
                  if (verbose)
                    cat("\nRule 4",
                        "\nThere is a discriminating path between",
                        md.path[1], "and", c, "for", b, ",and", b, "is not in Sepset of",
                        c, "and", md.path[1], ". Orient:", a, "<->", b, "<->",
                        c, "\n")
                  pag[a, b] <- pag[b, c] <- pag[c, b] <- 2
                }
                Done <- TRUE
              }
            }
          }
        }
      }
      ##-- R5 ------------------------------------------------------------------
      if (rules[5]) {
        ind <- which((pag == 1 & t(pag) == 1), arr.ind = TRUE) ## a o-o b
        while (length(ind) > 0) {
          a <- ind[1, 1]
          b <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all c s.t. a o-o c and c is not connected to b
          indC <- which((pag[a, ] == 1 & pag[, a] == 1) & (pag[b, ] == 0 & pag[, b] == 0))
          ## delete b since it is surely in indC
          indC <- setdiff(indC, b)
          ## find all d s.t. b o-o d and d is not connected to a
          indD <- which((pag[b, ] == 1 & pag[, b] == 1) & (pag[a, ] == 0 & pag[, a] == 0))
          ## delete a since it is surely in indD
          indD <- setdiff(indD, a)
          if (length(indC) > 0 && length(indD) > 0) {
            counterC <- 0
            while ((counterC < length(indC)) && pag[a, b] == 1) {
              counterC <- counterC + 1
              c <- indC[counterC]
              counterD <- 0
              while ((counterD < length(indD)) && pag[a, b] == 1) {
                counterD <- counterD + 1
                d <- indD[counterD]
                ## this is the easiest one
                if (pag[c, d] == 1 && pag[d, c] == 1) {
                  if (length(unfVect) == 0) { ## normal version
                    pag[a, b] <- pag[b, a] <- 3
                    pag[a, c] <- pag[c, a] <- 3
                    pag[c, d] <- pag[d, c] <- 3
                    pag[d, b] <- pag[b, d] <- 3
                    if (verbose)
                      cat("\nRule 5",
                          "\nThere exists an uncovered circle path between", a, "and", b,
                          ". Orient:", a, "-", b, "and", a, "-", c, "-", d, "-", b, "\n")
                  }
                  else { ## conservative: check that every triple on the circle is faithful
                    path2check <- c(a,c,d,b)
                    if (faith.check(path2check, unfVect, p)) {
                      pag[a, b] <- pag[b, a] <- 3
                      pag[a, c] <- pag[c, a] <- 3
                      pag[c, d] <- pag[d, c] <- 3
                      pag[d, b] <- pag[b, d] <- 3
                      if (verbose)
                        cat("\nRule 5",
                            "\nThere exists a faithful uncovered circle path between",
                            a, "and", b, ". Conservatively orient:",
                            a, "-", b, "and", a, "-", c, "-", d, "-", b, "\n")
                    }
                  }
                }
                ## search with a breitensuche a minimal uncovered circle path
                else {
                  ## Find a minimal uncovered circle path for these a,b,c, and d.
                  ## This path has already been checked to be uncovered and
                  ## to be faithful for the conservative case
                  ucp <- minUncovCircPath(p, pag = pag, path = c(a,c,d,b),
                                          unfVect = unfVect, verbose = verbose)
                  ## there is a path ---> orient
                  if (length(ucp) > 1) {
                    ## orient every edge on the path as --
                    n <- length(ucp)
                    pag[ucp[1], ucp[n]] <- pag[ucp[n], ucp[1]] <- 3 ## a--b
                    for (j in 1:(length(ucp)-1)) ## each edge on the path --
                      pag[ucp[j], ucp[j + 1]] <- pag[ucp[j + 1], ucp[j]] <- 3
                  }
                }
              }
            }
          }
        }
      }
      ##-- R6 ------------------------------------------------------------------
      if (rules[6]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)## b o-* c
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          c <- ind[i, 2]
          if (any(pag[b, ] == 3 & pag[, b] == 3)) {
            pag[c, b] <- 3
            if (verbose)
              cat("\nRule 6",
                  "\nOrient:", b, "o-*", c, "as", b, "-*", c, "\n")
          }
        }
      }
      ##-- R7 ------------------------------------------------------------------
      if (rules[7]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          c <- ind[i, 2]
          indA <- which((pag[b, ] == 3 & pag[, b] == 1) & (pag[c, ] == 0 & pag[, c] == 0))
          indA <- setdiff(indA, c)
          if (length(indA) > 0) {
            if (length(unfVect) == 0) {
              pag[c, b] <- 3
              if (verbose)
                cat("\nRule 7",
                    "\nOrient:", indA, "-o", b, "o-*",
                    c, "as", b, "-*", c, "\n")
            }
            else for (a in indA)
              if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                pag[c, b] <- 3
                if (verbose)
                  cat("\nRule 7",
                      "\nConservatively orient:", a, "-o", b, "o-*",
                      c, "as", b, "-*", c, "\n")
              }
          }
        }
      }
      ##-- R8 ------------------------------------------------------------------
      if (rules[8]) {
        ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          c <- ind[i, 2]
          indB <- which(pag[, a] == 3 & (pag[a, ] == 2 | pag[a, ] == 1) &
                          pag[c, ] == 3 & pag[, c] == 2)
          if (length(indB) > 0) {
            pag[c, a] <- 3
            if (verbose)
              cat("\nRule 8",
                  "\nOrient:", a, "->", indB, "->", c,
                  "or", a, "-o", indB, "->", c, "with", a,
                  "o->", c, "as", a, "->", c, "\n")
          }
        }
      }
      ##-- R9 ------------------------------------------------------------------
      if (rules[9]) {
        ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE) ## a o-> c
        while (length(ind) > 0) {
          a <- ind[1, 1]
          c <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all b s.t. a (o-)--(o>) b and b and c are not connected
          indB <- which((pag[a, ] == 2 | pag[a, ] == 1) &
                          (pag[, a] == 1 | pag[, a] == 3) &
                          (pag[c, ] == 0 & pag[, c] == 0))
          ## delete c from indB since it is surely inside
          indB <- setdiff(indB, c)
          ## chose one b s.t. the initial structure exists and the edge hasn't been oriented yet
          while ((length(indB) > 0) && (pag[c,a] == 1)) {
            b <- indB[1]
            indB <- indB[-1]
            ## find a minimal uncovered pd path from initial (a,b,c) :
            upd <- minUncovPdPath(p, pag, a,b,c,
                                  unfVect = unfVect, verbose = verbose)
            ## there is a path ---> orient it
            if (length(upd) > 1) {
              pag[c, a] <- 3
              if (verbose)
                cat("\nRule 9",
                    "\nThere exists an uncovered potentially directed path between", a, "and", c,
                    ". Orient:", a, " ->",c, "\n")
            }
          }
        }
      }
      ##-- R10 ------------------------------------------------------------------
      if (rules[10]) {
        ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE) ## a o-> c
        while (length(ind) > 0) {
          a <- ind[1, 1]
          c <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all b s.t. b --> c
          indB <- which((pag[c, ] == 3 & pag[, c] == 2))
          if (length(indB) >= 2) {
            counterB <- 0
            while (counterB < length(indB) && (pag[c, a] == 1)) {
              counterB <- counterB + 1
              b <- indB[counterB]
              indD <- setdiff(indB, b)
              counterD <- 0
              while ((counterD < length(indD)) && (pag[c, a] == 1)) {
                counterD <- counterD + 1
                d <- indD[counterD]
                ## this is the easiest one
                if ((pag[a, b] == 1 || pag[a, b] == 2) &&
                    (pag[b, a] == 1 || pag[b, a] == 3) &&
                    (pag[a, d] == 1 || pag[a, d] == 2) &&
                    (pag[d, a] == 1 || pag[d, a] == 3) && pag[d, b] == 0 && pag[b, d] == 0) {
                  if (length(unfVect) == 0) { ## normal version
                    pag[c, a] <- 3
                    if (verbose)
                      cat("\nRule 10 [easy]",
                          "\nOrient:", a, "->", c, "\n")
                  }
                  else ## conservative version: check faithfulness of b-a-d
                    if (!any(unfVect == triple2numb(p,b,a,d), na.rm = TRUE) &&
                        !any(unfVect == triple2numb(p,d,a,b), na.rm = TRUE)) {
                      pag[c, a] <- 3
                      if (verbose)
                        cat("\nRule 10 [easy]",
                            "\nConservatively orient:", a, "->", c, "\n")
                    }
                }
                ## search with a breitensuche two minimal uncovered circle paths
                else {
                  ## find all x s.t. a (o-)--(o>) x
                  indX <- which((pag[a, ] == 1 | pag[a, ] == 2) &
                                  (pag[, a] == 1 | pag[, a] == 3), arr.ind = TRUE)
                  indX <- setdiff(indX, c)
                  if (length(indX >= 2)) {
                    counterX1 <- 0
                    while (counterX1 < length(indX) && pag[c, a] == 1) {
                      counterX1 <- counterX1 + 1
                      first.pos <- indA[counterX1]
                      indX2 <- setdiff(indX, first.pos)
                      counterX2 <- 0
                      while (counterX2 < length(indX2) && pag[c, a] == 1) {
                        counterX2 <- counterX2 + 1
                        sec.pos <- indX2[counterX2]
                        t1 <- minUncovPdPath(p, pag, a, first.pos, b,
                                             unfVect = unfVect, verbose = verbose)
                        if (length(t1) > 1) { # otherwise, can skip next minUnc..()
                          t2 <- minUncovPdPath(p, pag, a, sec.pos, d,
                                               unfVect = unfVect, verbose = verbose)
                          if (length(t2) > 1 &&
                              first.pos != sec.pos && pag[first.pos, sec.pos] == 0) {
                            ## we found 2 uncovered pd paths
                            if (length(unfVect) == 0) { ## normal version
                              pag[c, a] <- 3
                              if (verbose)
                                cat("\nRule 10", "\nOrient:", a, "->", c, "\n")
                            }
                            else if(!any(unfVect == triple2numb(p,first.pos, a, sec.pos), na.rm = TRUE) &&
                                    !any(unfVect == triple2numb(p,sec.pos, a, first.pos), na.rm = TRUE)) {
                              ## conservative version
                              pag[c, a] <- 3
                              if (verbose)
                                cat("\nRule 10",
                                    "\nConservatively orient:", a, "->", c, "\n")
                            }
                          }
                        }
                      } #  # while ( counterX2 .. )
                    }
                  }
                } # else
              } # while ( counterD .. )
            } # while ( counterB .. )
          } # if (length(indB) .)
        }
      } ## if (rules[10] ..)
    }
  }
  pag
}



JCI_skeleton <- function(suffStat, indepTest, alpha, labels, p,
                       method = c("stable", "original", "stable.fast"), m.max = Inf,
                       fixedGaps = NULL, fixedEdges = NULL,
                       NAdelete = TRUE, numCores = 1, verbose = FALSE,F_node_names=list(),F_tester=NULL)
{
  ## Purpose: Perform undirected part of PC-Algorithm, i.e.,
  ## estimate skeleton of DAG given data
  ## Order-independent version! NEU
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat: List containing all necessary elements for the conditional
  ##             independence decisions in the function "indepTest".
  ## - indepTest: predefined function for testing conditional independence
  ## - alpha: Significance level of individual partial correlation tests
  ## - fixedGaps: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical); gaps fixed here are not changed
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - numCores: number of cores to be used for calculation if
  ##   method = "stable.fast"
  ## ----------------------------------------------------------------------
  ## Value:
  ## - G, sepset, pMax, ord, n.edgetests
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 09.12.2009
  ## Modification: Diego Colombo; Martin Maechler; Alain Hauser
  
  ## x,y,S konstruieren
  ##-   tst <- try(indepTest(x,y,S, obj))
  ##-   if(inherits(tst, "try-error"))
  ##-     stop("the 'indepTest' function does not work correctly with 'obj'")
  
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p))
      p <- length(labels)
    else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    ## Don't want message, in case this is called e.g. from fciPlus():
    ## else
    ##   message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  ids_of_F_nodes= (p-length(F_node_names)+1):p # the last nodes are F nodes by assumption
  method <- match.arg(method)
  ## C++ version still has problems under Windows; will have to check why
  ##  if (method == "stable.fast" && .Platform$OS.type == "windows") {
  ##    method <- "stable"
  ##    warning("Method 'stable.fast' is not available under Windows; using 'stable' instead.")
  ##  }
  
  ## G := !fixedGaps, i.e. G[i,j] is true  iff  i--j  will be investigated
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps)) )
    stop("fixedGaps must be symmetric")
  else
    G <- !fixedGaps
  
  diag(G) <- FALSE
  
  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges)) )
    stop("fixedEdges must be symmetric")
  
  ## Check number of cores
  stopifnot((is.integer(numCores) || is.numeric(numCores)) && numCores > 0)
  if (numCores > 1 && method != "stable.fast") {
    warning("Argument numCores ignored: parallelization only available for method = 'stable.fast'")
  }
  if (method == "stable.fast") {
    ## Do calculation in C++...
    if (identical(indepTest, gaussCItest))
      indepTestName <- "gauss"
    else
      indepTestName <- "rfun"
    options <- list(
      verbose = as.integer(verbose),
      m.max = as.integer(ifelse(is.infinite(m.max), p, m.max)),
      NAdelete = NAdelete,
      numCores = numCores)
    res <- .Call("estimateSkeleton", G, suffStat, indepTestName, indepTest, alpha, fixedEdges, options);
    G <- res$amat
    ## sepset <- res$sepset
    sepset <- lapply(seq_p, function(i) c(
      lapply(res$sepset[[i]], function(v) if(identical(v, as.integer(-1))) NULL else v),
      vector("list", p - length(res$sepset[[i]])))) # TODO change convention: make sepset triangular
    pMax <- res$pMax
    n.edgetests <- res$n.edgetests
    ord <- length(n.edgetests) - 1L
  }
  else {
    ## Original R version
    pval <- NULL
    sepset <- lapply(seq_p, function(.) vector("list",p))# a list of lists [p x p]
    ## save maximal p value
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0L
    n.edgetests <- numeric(1)# final length = max { ord}
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord+1L] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      if (verbose)
        cat("Order=", ord, "; remaining edges:", remEdges,"\n",sep = "")
      if(method == "stable") {
        ## Order-independent version: Compute the adjacency sets for any vertex
        ## Then don't update when edges are deleted
        G.l <- split(G, gl(p,p))
      }
      for (i in 1:remEdges) {
        if(verbose && (verbose >= 2 || i%%100 == 0)) cat("|i=", i, "|iMax=", remEdges, "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (G[y, x] && !fixedEdges[y, x]) {
          if (x %in% ids_of_F_nodes && y %in% ids_of_F_nodes){
            # two F-nodes are NOT separable. Bidirected in-between
            pMax[x, y]<-1L
          }
          else{
            # at most one of x,y is an F node
            # sepset is augmented to include all the remaining F-nodes.
            nbrsBool <- if(method == "stable") G.l[[x]] else G[,x]
            nbrsBool[y] <- FALSE
            nbrs <- seq_p[nbrsBool]
            length_nbrs <- length(nbrs)
            if (length_nbrs >= ord) {
              if (length_nbrs > ord)
                done <- FALSE
              S <- seq_len(ord)
              #print(suffStat$dm)
              repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
                n.edgetests[ord1] <- n.edgetests[ord1] + 1
                # Check if either x, y are in ids_of_F_nodes
                if (x %in% ids_of_F_nodes || y %in% ids_of_F_nodes){
                  # IF GAUSSIAN, THERE IS AN ISSUE. I WAS ASSUMING DISC TEST CHECK INV WRT F NODE
                  if (x %in% ids_of_F_nodes){
                    # create a new suffStat using only pair of environments.
                    # || returns a scalar. | returns the vector which is what we need.
                    #newdf <-  suffStat$dm[suffStat$dm[labels[x]]==0 | suffStat$dm[labels[x]]==1,]
                    newdf<-suffStat$dm
                    # remove F-node columns
                    #other_F_node_names=F_node_names[F_node_names!=labels[y]] # remains a list
                    # update columns so that there is no constant. ow tester complains
                    for (local_name in colnames(newdf)){
                      if (length(unique(newdf[[local_name]]))==1){
                        newdf[[1,local_name]]<-newdf[[1,local_name]]+1 # this will do
                      }
                    }
                    #newdf <- newdf[,!(colnames(newdf)%in%drops)]                    
                    #print(nrow(newdf))
                    newsuffStat <- list(dm = newdf, adaptDF=FALSE)
                    #newsuffStat <- suffStat
                    other_F_nodes <- ids_of_F_nodes[ids_of_F_nodes!=x]
                  }
                  else{# y is an F-node
                    #print('we are here')
                    #newdf <-  suffStat$dm[suffStat$dm[labels[y]]==0 | suffStat$dm[labels[y]]==1,]
                    newdf<-suffStat$dm
                    # remove the other F-node columns. Not necessary but tester complains they are constant.
                    # REMARK: Dropping some of the F ndoes creates issues with labeling order changing. Instead just add 1 to an element in the F ndoes outside the one used in testing to make the test happy.
                    #other_F_node_names=F_node_names[F_node_names!=labels[y]] # remains a list
                    for (local_name in colnames(newdf)){
                      if (length(unique(newdf[[local_name]]))==1){
                        newdf[[1,local_name]]<-newdf[[1,local_name]]+1 # this will do
                      }
                      #print('After:')
                      #print(unique(newdf[[local_name]]))
                      #cat('\n')
                    }
                    newsuffStat <- list(dm = newdf, adaptDF=FALSE)
                    #newsuffStat <- suffStat
                    other_F_nodes <- ids_of_F_nodes[ids_of_F_nodes!=y]
                    
                  }
                  W <- nbrs[S]
                  #W <- W[!(W%in%other_F_nodes)] # drop F-nodes from conditioning
                }
                else{
                  # both are observational, we need to look at one of the environments - assuming soft intervention.
                  # Pick environment 1 for the first F node. This can be chosen in a smarter way if the data is not symmetric.
                  #newdf <-  suffStat$dm[suffStat$dm[,F_node_names[[1]]]==0,]  # THIS COMMA WASNT NEEDED IN TEST.R!!?
                  # in F_09, both env. don't have constant variables
                  #print(newdf)
                  # also need to drop this F-node column because otherwise indep test complains it has a single value
                  # also need to drop other F-nodes since they also may take single value
                  #drops=F_node_names # remains a list
                  #newdf <- newdf[,!(colnames(newdf)%in%drops)]
                  #print(names(newdf))
                  #print(newdf)
                  #print(nrow(newdf))
                  #newsuffStat <- list(dm = newdf, adaptDF=FALSE)
                  newsuffStat <- suffStat
                  # drop F-nodes from conditioning set, to be added later if indep.
                  W <- nbrs[S]
                  #W <- W[!(W%in%ids_of_F_nodes)] 
                }
                #pval <- indepTest(x, y, nbrs[S], newsuffStat)
                #print(names(newsuffStat$dm))
                #print(unique(newsuffStat$dm))
                #cat('We are at ')
                #cat(x,y,'\n')
                #print(nrow(newsuffStat$dm))
                #print(newsuffStat$dm)
                # print(x)
                # print(y)
                # print(W)
                # print(corr)
                #pval<-CCIT( np$array(df['X']),np$array(df['Y']),np$array(df['Z']))
                # GAUSSIAN NOT CODED                
                # if not Gaussian and discrete:
                #print(unique(newsuffStat$dm[['Jnk']]))
                cat(x,y,'\n')
                cat(W,'\n')
                pval <- indepTest(x, y, W, newsuffStat)
                # if gaussian, need to provide the test corr matrix
                #pval <- indepTest(x, y, W, newsuffStat)
                if (verbose)# && pval>=alpha)
                {cat("x=", names(df)[x], " y=", names(df)[y], " S=", names(df)[W], ": pval =", pval, "\n")#nbrs[S]
                }
                if(is.na(pval))
                  pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                if (pMax[x, y] < pval)
                  pMax[x, y] <- pval
                if(pval >= alpha) { 
                  G[x, y] <- G[y, x] <- FALSE
                  
                  if (x %in% ids_of_F_nodes || y %in% ids_of_F_nodes){
                    if (x %in% ids_of_F_nodes){
                      # add all the remaining F-nodes to the separating set
                      sepset[[x]][[y]] <- union(nbrs[S],ids_of_F_nodes[ids_of_F_nodes!=x])
                    }
                    else{
                      sepset[[x]][[y]] <- union(nbrs[S],ids_of_F_nodes[ids_of_F_nodes!=y])
                    }
                  } 
                  else{
                    sepset[[x]][[y]] <- union(nbrs[S],ids_of_F_nodes)
                  }
                  #print(names(df)[sepset[[x]][[y]]])
                  #cat("\n")
                  break
                }
                else {
                  nextSet <- getNextSet(length_nbrs, ord, S)
                  if (nextSet$wasLast)
                    break
                  S <- nextSet$nextSet
                }
              } ## repeat
            }
          }
        }
      }# for( i )
      ord <- ord + 1L
    } ## while()
    for (i in 1:(p - 1)) {
      for (j in 2:p)
        pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j,i])
    }
  }
  
  ## transform matrix to graph object :
  Gobject <-
    if (sum(G) == 0) {
      new("graphNEL", nodes = labels)
    } else {
      colnames(G) <- rownames(G) <- labels
      as(G,"graphNEL")
    }
  
  ## final object
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}## end{ JCI_skeleton }

