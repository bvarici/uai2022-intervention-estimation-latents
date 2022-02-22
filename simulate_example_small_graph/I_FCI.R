###############################################################################
## I-FCI
###############################################################################

I_fci <- function(suffStat, indepTest, alpha, labels, p,
                skel.method = c("stable", "original", "stable.fast"),
                type = c("normal", "anytime", "adaptive"),
                fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                m.max = Inf, pdsep.max = Inf, rules = rep(TRUE, 10),
                doPdsep = TRUE, biCC = FALSE, conservative = FALSE,
                maj.rule = FALSE, numCores = 1, verbose = FALSE, F_tester=NULL, F_node_names = list())
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
  
  skel <- I_skeleton(suffStat, indepTest, alpha, labels = labels, method = skel.method,
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
  res <- udag2pag(pag = G, sepset, rules = rules, unfVect = tripleList,
                  verbose = verbose)
  colnames(res) <- rownames(res) <- labels
  new("fciAlgo", amat = res, call = cl, n = integer(0),
      max.ord = as.integer(max.ordSKEL),
      max.ordPDSEP = as.integer(max.ordPD),
      n.edgetests = n.edgetestsSKEL, n.edgetestsPDSEP = n.edgetestsPD,
      sepset = sepset, pMax = pMax, allPdsep = allPdsep)
} ## {I-fci}

I_skeleton <- function(suffStat, indepTest, alpha, labels, p,
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
            # two F-nodes are separable given emptyset
            G[x, y] <- G[y, x] <- FALSE
            sepset[[x]][[y]] <- integer(0) # represents emptyset
            pMax[x, y]<-1L
          }
          else{
            # at most one of x,y is an F node
            # procedure for both cases are almost identical except the data used for one is F node is different than the whole dataset
            # furthermore, sepset is augmented to include all the remaining F-nodes.
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
                    newdf <-  suffStat$dm[suffStat$dm[labels[x]]==0 | suffStat$dm[labels[x]]==1,]
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
                    other_F_nodes <- ids_of_F_nodes[ids_of_F_nodes!=x]
                  }
                  else{# y is an F-node
                    #print('we are here')
                    newdf <-  suffStat$dm[suffStat$dm[labels[y]]==0 | suffStat$dm[labels[y]]==1,]
                    # remove the other F-node columns. Not necessary but tester complains they are constant.
                    # REMARK: Dropping some of the F ndoes creates issues with labeling order changing. Instead just add 1 to an element in the F ndoes outside the one used in testing to make the test happy.
                    #other_F_node_names=F_node_names[F_node_names!=labels[y]] # remains a list
                    for (local_name in colnames(newdf)){
                      #print(unique(newdf[local_name]))
                      #print(length(unique(newdf[local_name])))
                      #print('Before:')
                      #print(unique(newdf[[local_name]]))
                      if (length(unique(newdf[[local_name]]))==1){
                        newdf[[1,local_name]]<-newdf[[1,local_name]]+1 # this will do
                      }
                      #print('After:')
                      #print(unique(newdf[[local_name]]))
                      #cat('\n')
                    }
                    # print('FINALLY MAIN F NODE:')
                    # print(labels)
                    # print(y)
                    # print(labels[y])
                    # print(names(newdf))
                    # print(newdf)
                    # print(unique(newdf[labels(y)]))
                    #newdf <- newdf[,!(colnames(newdf)%in%drops)]
                    #print(newdf)
                    #print(y)
                    #print(labels[y])
                    #print(suffStat$dm[labels[y]]==1)
                    #print(newdf)
                    newsuffStat <- list(dm = newdf, adaptDF=FALSE)
                    other_F_nodes <- ids_of_F_nodes[ids_of_F_nodes!=y]
                    #print(other_F_nodes)
                    # local_data<-newdf
                    # for (j in names(local_data)){
                    #   print(j)
                    #   print(unique(local_data[j]))
                    #   cat('\n')
                    # }
                  
                    
                  }
                W <- nbrs[S]
                W <- W[!(W%in%other_F_nodes)] # drop F-nodes from conditioning
                }
                else{
                  # both are observational, we need to look at one of the environments - assuming soft intervention.
                  # Pick environment 1 for the first F node. This can be chosen in a smarter way if the data is not symmetric.
                  newdf <-  suffStat$dm[suffStat$dm[,F_node_names[[1]]]==0,]  # THIS COMMA WASNT NEEDED IN TEST.R!!?
                  # in F_09, both env. don't have constant variables
                  #print(newdf)
                  # also need to drop this F-node column because otherwise indep test complains it has a single value
                  # also need to drop other F-nodes since they also may take single value
                  drops=F_node_names # remains a list
                  newdf <- newdf[,!(colnames(newdf)%in%drops)]
                  #print(names(newdf))
                  #print(newdf)
                  #print(nrow(newdf))
                  newsuffStat <- list(dm = newdf, adaptDF=FALSE)
                  # drop F-nodes from conditioning set, to be added later if indep.
                  W <- nbrs[S]
                  W <- W[!(W%in%ids_of_F_nodes)] 
                }
                #pval <- indepTest(x, y, nbrs[S], newsuffStat)
                #print(names(newsuffStat$dm))
                #print(unique(newsuffStat$dm))
                #cat('We are at ')
                #cat(x,y,'\n')
                #print(nrow(newsuffStat$dm))
                #print(newsuffStat$dm)
                if (identical(indepTest, gaussCItest)){
                  if (x %in% ids_of_F_nodes || y %in% ids_of_F_nodes){
                    # always mapping F node to x for internal workings of CCIT
                    if (x %in% ids_of_F_nodes){
                      #print('hereX')
                      if(length(W)==0){
                        print('here!')
                        # need to feed NONE
                        #pval<-F_tester( np$array(newsuffStat$dm[y]),np$array(newsuffStat$dm[x]),NULL,num_iter = as.integer(30), bootstrap = TRUE, nthread = as.integer(20))
                        #pval<-F_tester( np$array(newsuffStat$dm[y]),np$array(newsuffStat$dm[x]),NULL,max_depths=c(as.integer(2),as.integer(3)),n_estimators=c(as.integer(10),as.integer(20)))
                        all_p_values_local=c()
                        # subsample and take the median p-value
                        for (a in 1:20){
                          ranges<-sample(nrow(newsuffStat$dm),1000)
                          #print(ranges)
                          #b=np$array(newsuffStat$dm[ranges,y])
                          pval<-F_tester( np$array(data.frame(newsuffStat$dm[ranges,x])),np$array(data.frame(newsuffStat$dm[ranges,y])),NULL,max_depths=c(as.integer(2),as.integer(3)),n_estimators=c(as.integer(10),as.integer(20)))
                          #pval<-F_tester( np$array(data.frame(newsuffStat$dm[ranges,y])),np$array(data.frame(newsuffStat$dm[ranges,x])),NULL)
                          #pval<-F_tester( np$array(newsuffStat$dm[y]),np$array(newsuffStat$dm[x]),NULL)
                          all_p_values_local<-c(all_p_values_local,pval)
                        }
                        pval<-median(all_p_values_local)
                        print(pval)
                      }
                      else{
                      pval<-F_tester( np$array(newsuffStat$dm[x]),np$array(newsuffStat$dm[y]),np$array(newsuffStat$dm[W]))
                    }}
                    else{
                      #print(ncol(newsuffStat$dm[y]))
                      #print(ncol(newsuffStat$dm[x]))
                      #print((newsuffStat$dm[W]))
                      
                      if(length(W)==0){
                        print('here!')
                        # need to feed NONE
                        #pval<-F_tester( np$array(newsuffStat$dm[y]),np$array(newsuffStat$dm[x]),NULL,num_iter = as.integer(30), bootstrap = TRUE, nthread = as.integer(20))
                        #pval<-F_tester( np$array(newsuffStat$dm[y]),np$array(newsuffStat$dm[x]),NULL,max_depths=c(as.integer(2),as.integer(3)),n_estimators=c(as.integer(10),as.integer(20)))
                        all_p_values_local=c()
                        # subsample and take the median p-value
                        for (a in 1:20){
                          ranges<-sample(nrow(newsuffStat$dm),1000)
                          #print(ranges)
                          #b=np$array(newsuffStat$dm[ranges,y])
                          pval<-F_tester( np$array(data.frame(newsuffStat$dm[ranges,y])),np$array(data.frame(newsuffStat$dm[ranges,x])),NULL,max_depths=c(as.integer(2),as.integer(3)),n_estimators=c(as.integer(10),as.integer(20)))
                          #pval<-F_tester( np$array(data.frame(newsuffStat$dm[ranges,y])),np$array(data.frame(newsuffStat$dm[ranges,x])),NULL)
                          #pval<-F_tester( np$array(newsuffStat$dm[y]),np$array(newsuffStat$dm[x]),NULL)
                          all_p_values_local<-c(all_p_values_local,pval)
                        }
                        pval<-median(all_p_values_local)
                        print(pval)
                      }
                      else{
                        #pval<-F_tester( np$array(newsuffStat$dm[y]),np$array(newsuffStat$dm[x]),np$array(newsuffStat$dm[W]),num_iter = as.integer(30), bootstrap = TRUE, nthread = as.integer(20))
                        pval<-F_tester( np$array(newsuffStat$dm[y]),np$array(newsuffStat$dm[x]),np$array(newsuffStat$dm[W]))
                      }
                    }}
                  else{# then we can use Gaussian test
                      n_local<-nrow(newsuffStat$dm)
                      corr<-cor(newsuffStat$dm)
                      newsuffStat <- list(C = corr, n = n_local)
                      pval <- indepTest(x, y, W, newsuffStat)
                  }
                }
                  # print(x)
                  # print(y)
                  # print(W)
                  # print(corr)
                  #pval<-CCIT( np$array(df['X']),np$array(df['Y']),np$array(df['Z']))
                
                else{
                  # if not Gaussian and discrete:
                  #print(unique(newsuffStat$dm[['Jnk']]))
                  cat(x,y,'\n')
                  cat(W,'\n')
                  pval <- indepTest(x, y, W, newsuffStat)
                }
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
}## end{ I_skeleton }

