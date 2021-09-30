  GBM <- function(input.data, ref.data, weight.matrix = NULL, method = "BT", n=NULL)
{
  input.data$Z = qchisq(input.data$P,df = 1,lower = F)
  Z = input.data$Z
  p = input.data$P
  M = length(Z)
  if (is.null(weight.matrix)) {W = rep(1, M)} else {W = weight.matrix}
  ldmat = cor(as.matrix(ref.data), use = "p")
  ldmat = 0.99 * ldmat + diag(0.01, nrow=dim(ldmat)[2], ncol=dim(ldmat)[2])
  
  if (length(Z) <= 1) stop("less than 1 SNP.")
  if (!(length(Z) == dim(ldmat)[1] & dim(ldmat)[2] == dim(ldmat)[1])) {
    stop("dimension do not match. Check dimension of Z and R.")
  }
  if(method %in% c("TPM","FCP","Simes","GM")){
    sumstat.GBM <- as.function(get(paste("sumstat",test, sep = ".")))
    return(sumstat.GBM(p))
  }else if(method %in% c("ExtSimes","RTP","SimpleM","ARTP",
                         "DOT","TQ","GATES","VEGAS",
                         "aSPUs","HYST","COMBAT","ART","ART.A",
                         "ACATO","BJ","GBJ","HC","GHC")){
    sumstat.GBM <- as.function(get(paste("sumstat",test, sep = ".")))
    return(sumstat.GBM(p,ldmat))
  }else if(method %in% c("MLR","FLM","PCA")){
    sumstat.GBM <- as.function(get(paste("sumstat",test, sep = ".")))
    if(is.null(n)) {
      stop("require the sample number n")
    }else{
      return(sumstat.GBM(p,ldmat,n))
    }
  }else if(method %in% c("BT","SKAT","SKATO","ACAT")){
    sumstat.GBM <- as.function(get(paste("sumstat",test, sep = ".")))
    return(sumstat.GBM(p,ldmat,W))
  }else {
    stop("Please use the right method!")
  }
}

  
  
  KATpval <- function (Q.all, lambda, acc = 0.01, lim = 1e+07) 
  {
    if (all(abs(lambda - lambda[1])/max(abs(lambda)) < 1e-10)) {
      return(pchisq(Q.all/mean(lambda), length(lambda), lower = FALSE))
    }
    pval = Liu.pval(Q.all, lambda)
    i1 = which(is.finite(Q.all))
    for (i in i1) {
      tmp = Wdavies(Q.all[i], lambda, acc = acc * pval[i], 
                    lim = lim)
      pval[i] = tmp$Qq
      if ((tmp$ifault > 0) | (pval[i] <= 0) | (pval[i] >= 1)) {
        sad = rep(1, length(Q.all))
        if (sum(Q.all > 0) > 0) {
          sad[Q.all > 0] = sapply(Q.all[Q.all > 0], saddle, lambda = lambda)
        }
        id = which(is.na(sad))
        if (length(id) > 0) {
          sad[id] = Liu.pval(Q.all[id], lambda)
        }
        pval[i] = sad
      }
    }
    return(pval)
  }
  
  Wdavies <- function (qx, lambda, delta = NULL, lim = 1e+09, acc = 1e-06) 
  {
    m = length(lambda)
    if (is.null(delta))  delta = rep(0, m)
    h = rep(1, m)
    sigma = 0
    out = .C("qfc", as.double(lambda), as.double(delta), as.integer(h), 
             as.integer(m), as.double(sigma), as.double(qx), as.integer(lim), 
             as.double(acc), trace = as.double(rep(0, 7)), ifault = as.integer(0), 
             res = as.double(0), PACKAGE = "mkatr")
    out$res = 1 - out$res
    return(list(Qq = out$res, trace = out$trace, ifault = out$ifault))
  }
  
  Liu.pval <- function (Q.all, lambda) 
  {
    c1 = rep(0, 4)
    for (i in 1:4) {
      c1[i] = sum(lambda^i)
    }
    muQ = c1[1]
    sigmaQ = sqrt(2 * c1[2])
    s1 = c1[3]/c1[2]^(3/2)
    s2 = c1[4]/c1[2]^2
    if (s1^2 > s2) {
      a = 1/(s1 - sqrt(s1^2 - s2))
      d = s1 * a^3 - a^2
      l = a^2 - 2 * d
    }else {
      l = 1/s2
      a = sqrt(l)
      d = 0
    }
    muX = l + d
    sigmaX = sqrt(2) * a
    param = list(l = l, d = d, muQ = muQ, muX = muX, sigmaQ = sigmaQ, 
                 sigmaX = sigmaX)
    Qx = (Q.all - param$muQ)/param$sigmaQ * param$sigmaX + param$muX
    return(pchisq(Qx, df = param$l, ncp = param$d, lower.tail = FALSE))
  }
  
  LD_R <- function (cor_G) 
  {
    if (is.positive.definite(cor_G) == FALSE) {
      cor_G <- make.positive.definite(cor_G)
    }
    if (is.positive.definite(cor_G) == FALSE) {
      cor_G <- cor(as.matrix(x), use = "p")
      diag(x) <- 1.0001
    }
    if (is.positive.definite(cor_G) == FALSE) {
      diag(cor_G) <- 1.001
    }
    if (is.positive.definite(cor_G) == FALSE) {
      diag(cor_G) <- 1.01
    }
    return(cor_G)
  }
  
  nsp <- function (ldmat, eps = NULL, ...) 
  {
    if (is.null(eps)) 
      eps <- sqrt(.Machine$double.eps)
    CC <- eigen(ldmat, TRUE)
    u <- CC$vectors
    d <- CC$values
    CC <- d > d[1] * eps
    if (!all(CC)) {
      d <- d[CC]
      u <- u[,CC]
    }
    L <- length(d)
    d <- sqrt(1/d)
    H <- u %*% (d * t(u))
    H <- 0.5 * (H + t(H))
    list(H = H, L = L)
  }
  
  vegas.call <- function (x, cor_G, vegas.pct, n_simul) 
  {
    stopifnot(length(x) == ncol(cor_G))
    vegas_vec <- ceiling(vegas.pct * ncol(cor_G))
    vegas_vec <- sort(vegas_vec)
    if (vegas_vec[1] > 1) {
      vegas.pct <- c(0, vegas.pct)
      vegas_vec <- c(1, vegas_vec)
    }
    chisq_vec <- qchisq(x, 1, lower.tail = FALSE)
    chisq_vec[chisq_vec == Inf] <- 60
    n_snps <- length(x)
    n_tests <- length(vegas_vec)
    TS_obs <- rep(NA, n_tests)
    TS_obs[1] <- max(chisq_vec, na.rm = TRUE)
    chisq_vec <- sort(chisq_vec, decreasing = TRUE)
    for (j in 2:n_tests) TS_obs[j] <- sum(chisq_vec[1:vegas_vec[j]])
    rd <- rmvnorm(n_simul, mean = rep(0, n_snps), sigma = cor_G)
    rd2 <- rd^2
    rd2 <- apply(rd2, 1, sort, decreasing = TRUE)
    pPerm0 <- rep(NA, n_tests)
    T0s <- apply(rd2, 2, max)
    pPerm0[1] <- (sum(T0s >= TS_obs[1]) + 1)/(length(T0s) + 1)
    for (j in 2:n_tests) {
      for (i in 1:n_simul) T0s[i] <- sum(rd2[1:vegas_vec[j], 
                                             i])
      pPerm0[j] <- (sum(T0s >= TS_obs[j]) + 1)/(length(T0s) + 
                                                  1)
    }
    v1 <- paste0("VEGAS.p", vegas.pct)
    v1[vegas_vec == ncol(cor_G)] <- "VEGAS.all"
    v1[vegas_vec == 1] <- "VEGAS.max"
    names(pPerm0) <- v1
    pPerm0
  }
  
  `%^%` <-  function (U, k) 
  {
    UUU <- eigen(U, symmetric = TRUE)
    Uval <- UUU$val
    Uvec <- UUU$vec
    Uvec <- Uvec[, Uval > 1e-07]
    Uval <- Uval[Uval > 1e-07]
    Uvec %*% (t(Uvec) * (Uval^k))
  }
  
  getbasismatrix <- function(evalarg, basisobj, nderiv=0, returnMatrix=FALSE) 
  {
    if (is.numeric(basisobj) && inherits(evalarg, "basisfd")) {
      temp     <- basisobj
      basisobj <- evalarg
      evalarg  <- temp
    }
    ##
    ##  check EVALARG
    ##
    #  if (!(is.numeric(evalarg)))  stop("Argument EVALARG is not numeric.")
    if(is.null(evalarg)) stop('evalarg required;  is NULL.')
    Evalarg <- evalarg
    # turn off warnings in checking if argvals can be converted to numeric
    op <- options(warn=-1)
    evalarg <- as.numeric(Evalarg)
    options(op)
    nNA <- sum(is.na(evalarg))
    if(nNA>0)
      stop('as.numeric(evalarg) contains ', nNA,
           ' NA', c('', 's')[1+(nNA>1)],
           ';  class(evalarg) = ', class(Evalarg))
    
    #  check basisobj
    
    if (!(inherits(basisobj, "basisfd")))
      stop("Second argument is not a basis object.")
    
    #  search for stored basis matrix
    
    if (!(length(basisobj$basisvalues) == 0 || is.null(basisobj$basisvalues))) {
      #  one or more stored basis matrices found,
      #  check that requested derivative is available
      if (!is.vector(basisvalues)) stop("BASISVALUES is not a vector.")
      basisvalues <- basisobj$basisvalues
      nvalues     <- length(basisvalues)
      #  search for argvals match
      N  <- length(evalarg)
      OK <- FALSE
      for (ivalues in 1:nvalues) {
        basisvaluesi <- basisvalues[ivalues]
        if (!is.list(basisvaluesi)) stop("BASISVALUES does not contain lists.")
        argvals <- basisvaluesi[[1]]
        if (!length(basisvaluesi) < nderiv+2) {
          if (N == length(argvals)) {
            if (all(argvals == evalarg)) {
              basismat <- basisvaluesi[[nderiv+2]]
              OK <- TRUE
            }
          }
        }
      }
      #   dimnames
      dimnames(basismat) <- list(NULL, basisobj$names)
      
      if (OK){
        if((!returnMatrix) && (length(dim(basismat)) == 2)){
          return(as.matrix(basismat))
        }
        return(basismat)
      }
    }
    
    #  Extract information about the basis
    
    type     <- basisobj$type
    nbasis   <- basisobj$nbasis
    params   <- basisobj$params
    rangeval <- basisobj$rangeval
    dropind  <- basisobj$dropind
    
    #  -----------------------------  B-spline basis  -------------------
    
    if (type == "bspline") {
      if (length(params) == 0) {
        breaks   <- c(rangeval[1], rangeval[2])
      } else {
        breaks   <- c(rangeval[1], params, rangeval[2])
      }
      norder   <- nbasis - length(breaks) + 2
      basismat <- bsplineS(evalarg, breaks, norder, nderiv, returnMatrix)
      
      #  -----------------------------  Constant basis  --------------------
      
    } else if (type == "const") {
      basismat  <- matrix(1,length(evalarg),1)
      
      #  -------------------------------  Fourier basis  -------------------
      
    } else if (type == "fourier") {
      period   <- params[1]
      basismat <- fourier(evalarg, nbasis, period, nderiv)
      
      #  -----------------------  Unrecognizable basis  --------------------
      
    } else {
      stop("Basis type not recognizable")
    }
    #  dimnames
    dimnames(basismat) <- list(NULL, basisobj$names)
    
    #  remove columns for bases to be dropped
    
    if (length(dropind) > 0) basismat <- basismat[,-dropind]
    if (length(evalarg) == 1) {
      basismat = matrix(basismat,1,length(basismat))
    }
    
    
    if((!returnMatrix) && (length(dim(basismat)) == 2)){
      #  coerce basismat to be nonsparse
      return(as.matrix(basismat))
    } else {
      #  allow basismat to be sparse if it already is
      return(as.matrix(basismat))
    }
    
  }
  
  int2Lfd <- function(m=0)
  {
    
    if (inherits(m, "Lfd")) {
      Lfdobj <- m
      return(Lfdobj)
    }
    
    if (!is.numeric(m)) 
      stop("Argument not numeric and not a linear differential operator.")
    
    
    if (length(m) != 1) stop("Argument is not a scalar.")
    
    if (round(m) != m)  stop("Argument is not an integer.")
    
    if (m < 0)   stop("Argument is negative.")
    
    #  all the checks passed, set up a functional data object
    #  The range will not be used in this case, and can be set
    #  to [0, 1]
    
    #  set up the list object for the homogeneous part
    
    if (m==0) {
      #  if derivative is zero, BWTLIST is empty
      bwtlist <- NULL
    } else {
      basisobj <- create.constant.basis(c(0,1))
      bwtlist  <- vector("list", m)
      for (j in 1:m) bwtlist[[j]] <- fd(0, basisobj)
    }
    
    #  define the Lfd object
    
    Lfdobj <- Lfd(m, bwtlist)
    
    return(Lfdobj)
    
  }
  
  Lfd <- function(nderiv=0, bwtlist=vector("list",0))
  {
    
    if (!is.numeric(nderiv))
      stop("Order of operator is not numeric.")
    if (nderiv != round(nderiv))
      stop("Order of operator is not an integer.")
    if (nderiv < 0)
      stop("Order of operator is negative.")
    
    #  check that bwtlist is either a list or a fd object
    
    if (!inherits(bwtlist, "list") && !inherits(bwtlist, "fd") &&
        !is.null(bwtlist) && !missing(bwtlist))
      stop("BWTLIST is neither a LIST or a FD object")
    
    #  if bwtlist is missing or NULL, convert it to a constant basis FD object
    
    if (is.null(bwtlist)) {
      bwtlist <- vector("list", nderiv)
      if (nderiv > 0) {
        conbasis <- create.constant.basis()
        for (j in 1:nderiv) bwtlist[[j]]  <- fd(0, conbasis)
      }
    }
    
    #  if BWTLIST is a fd object, convert to a list object.
    
    if (inherits(bwtlist, "fd")) bwtlist <- fd2list(bwtlist)
    
    #  check size of bwtlist
    
    nbwt <- length(bwtlist)
    
    if (nbwt != nderiv & nbwt != nderiv + 1)
      stop("The size of bwtlist inconsistent with NDERIV.")
    
    #  check individual list entries for class
    #  and find a default range
    
    if (nderiv > 0) {
      rangevec <- c(0,1)
      for (j in 1:nbwt) {
        bfdj <- bwtlist[[j]]
        if (inherits(bfdj, "fdPar")) {
          bfdj <- bfdj$fd
          bwtlist[[j]] <- bfdj
        }
        if (!inherits(bfdj, "fd") && !inherits(bfdj, "integer"))
          stop(paste("An element of BWTLIST contains something other ",
                     " than an fd object or an integer"))
        if (inherits(bfdj, "fd")) {
          bbasis   <- bfdj$basis
          rangevec <- bbasis$rangeval
        } else {
          if (length(bfdj) == 1) {
            bwtfd <- fd(bfdj, conbasis)
            bwtlist[[j]] <- bwtfd
          }
          else stop("An element of BWTLIST contains a more than one integer.")
        }
      }
      
      #  check that the ranges are compatible
      
      for (j in 1:nbwt) {
        bfdj    <- bwtlist[[j]]
        if (inherits(bfdj, "fdPar")) bfdj <- bfdj$fd
        bbasis <- bfdj$basis
        btype  <- bbasis$type
        #  constant basis can have any range
        if (!btype == "const") {
          brange = bbasis$rangeval
          if (any(rangevec != brange)) stop(
            "Ranges are not compatible.")
        }
      }
    }
    
    #  Save call
    
    Lfd.call <- match.call()
    
    #  set up the Lfd object
    
    #  S4 definition
    # Lfdobj <- new("Lfd", call=Lfd.call, nderiv=nderiv, bwtlist=bwtlist)
    
    #  S3 definition
    
    Lfdobj <- list(call=Lfd.call, nderiv=nderiv, bwtlist=bwtlist)
    oldClass(Lfdobj) <- "Lfd"
    
    Lfdobj
    
  }
  
  fourier <- function(x, nbasis = n, period = span, nderiv = 0)
  {
    xNames <- names(x)
    #
    x      <- as.vector(x)
    n      <- length(x)
    onen   <- rep(1,n)
    xrange <- range(x)
    span   <- xrange[2] - xrange[1]
    
    #  check period and set up omega
    
    if (period <= 0) stop("PERIOD not positive.")
    omega  <- 2*pi/period
    omegax <- omega*x
    
    #  check nbasis
    
    if (nbasis <= 0) stop("NBASIS not positive")
    
    #  check nderiv
    
    if (nderiv <  0) stop("NDERIV is negative.")
    
    #  if nbasis is even, add one
    
    if ((nbasis %/% 2)*2 == nbasis) nbasis <- nbasis + 1
    
    #  compute basis matrix
    
    basismat <- matrix(0,n,nbasis)
    if (nderiv == 0) {
      #  The fourier series itself is required.
      basismat[,1] <- 1/sqrt(2)
      if(nbasis>1){
        j    <- seq(2,nbasis-1,2)
        k    <- j/2
        args <- outer(omegax,k)
        basismat[,j]   <- sin(args)
        basismat[,j+1] <- cos(args)
      }
    } else {
      #  The derivative of order nderiv is required.
      basismat[,1] <- 0.0
      if(nbasis>1){
        if (nderiv == (nderiv %/% 2)*2) {
          mval  <- nderiv/2
          ncase <- 1
        } else {
          mval <- (nderiv-1)/2
          ncase <- 2
        }
        j    <- seq(2,nbasis-1,2)
        k    <- j/2
        fac  <- outer(onen,((-1)^mval)*(k*omega)^nderiv)
        args <- outer(omegax,k)
        if (ncase == 1) {
          basismat[,j]   <-  fac * sin(args)
          basismat[,j+1] <-  fac * cos(args)
        } else {
          basismat[,j]   <-  fac * cos(args)
          basismat[,j+1] <- -fac * sin(args)
        }
      }
    }
    basismat <- basismat/sqrt(period/2)
    #
    fNames <- "const"
    n2 <- (nbasis %/%2)
    if(n2>0){
      SC <- outer(c("sin", "cos"), 1:n2, paste, sep="")
      fNames <- c(fNames, as.vector(SC))
    }
    #
    dimnames(basismat) <- list(xNames, fNames)
    #
    return(basismat)
  }
  
  eval.basis <- function(evalarg, basisobj, Lfdobj=0, returnMatrix=FALSE) 
  {
    
    if (is.numeric(basisobj) && inherits(evalarg, "basisfd")) {
      temp     <- basisobj
      basisobj <- evalarg
      evalarg  <- temp
    }
    
    #  check EVALARG
    
    #  if (!(is.numeric(evalarg))){# stop("Argument EVALARG is not numeric.")
    # turn off warnings in checking if argvals can be converted to numeric.
    if(is.numeric(evalarg)){
      if(!is.vector(evalarg))stop("Argument 'evalarg' is not a vector.")
      Evalarg <- evalarg
    } else {
      op <- options(warn=-1)
      Evalarg <- as.numeric(evalarg)
      options(op)
      nNA <- sum(is.na(Evalarg))
      if(nNA>0)
        stop('as.numeric(evalarg) contains ', nNA,
             ' NA', c('', 's')[1+(nNA>1)],
             ';  class(evalarg) = ', class(evalarg))
      #  if(!is.vector(Evalarg))
      #      stop("Argument EVALARG is not a vector.")
    }
    
    #  check basisobj
    
    if (!(inherits(basisobj, "basisfd"))) stop(
      "Second argument is not a basis object.")
    
    #  check LFDOBJ
    
    Lfdobj <- int2Lfd(Lfdobj)
    ##
    ## 2.  set up
    ##
    #  determine the highest order of derivative NDERIV required
    
    nderiv <- Lfdobj$nderiv
    
    #  get weight coefficient functions
    
    bwtlist <- Lfdobj$bwtlist
    ##
    ## 3.  Do
    ##
    #  get highest order of basis matrix
    
    basismat <- getbasismatrix(evalarg, basisobj, nderiv, returnMatrix)
    
    #  Compute the weighted combination of derivatives is
    #  evaluated here if the operator is not defined by an
    #  integer and the order of derivative is positive.
    
    if (nderiv > 0) {
      nbasis    <- dim(basismat)[2]
      oneb      <- matrix(1,1,nbasis)
      nonintwrd <- FALSE
      for (j in 1:nderiv) {
        bfd    <- bwtlist[[j]]
        bbasis <- bfd$basis
        if (bbasis$type != "constant" || bfd$coefs != 0) nonintwrd <- TRUE
      }
      if (nonintwrd) {
        for (j in 1:nderiv) {
          bfd   <- bwtlist[[j]]
          if (!all(c(bfd$coefs) == 0.0)) {
            wjarray   <- eval.fd(evalarg, bfd, 0, returnMatrix)
            Dbasismat <- getbasismatrix(evalarg, basisobj, j-1,
                                        returnMatrix)
            basismat  <- basismat + (wjarray %*% oneb)*Dbasismat
          }
        }
      }
    }
    
    if((!returnMatrix) && (length(dim(basismat)) == 2)){
      return(as.matrix(basismat))
    } else {
      return(basismat)
    }
    
  }
  
  basisfd <- function (type, rangeval, nbasis, params, dropind = vector("list", 0), 
                       quadvals = vector("list", 0), values = vector("list",0), basisvalues = vector("list", 0)) 
  {
    if (nargs() == 0) {
      type <- "bspline"
      rangeval <- c(0, 1)
      nbasis <- 2
      params <- vector("list", 0)
      dropind <- vector("list", 0)
      quadvals <- vector("list", 0)
      values <- vector("list", 0)
      basisvalues <- vector("list", 0)
      basisobj <- list(type = type, rangeval = rangeval, nbasis = nbasis, 
                       params = params, dropind = dropind, quadvals = quadvals, 
                       values = values, basisvalues = basisvalues)
      oldClass(basisobj) <- "basisfd"
      return(basisobj)
    }
    type = "fourier"
    
    if (missing(quadvals)) 
      quadvals <- vector("list", 0)
    else if (!(length(quadvals) == 0 || is.null(quadvals))) {
      nquad <- dim(quadvals)[1]
      ncol <- dim(quadvals)[2]
      if ((nquad == 2) && (ncol > 2)) {
        quadvals <- t(quadvals)
        nquad <- dim(quadvals)[1]
        ncol <- dim(quadvals)[2]
      }
      if (nquad < 2) 
        stop("Less than two quadrature points are supplied.")
      if (ncol != 2) 
        stop("'quadvals' does not have two columns.")
    }
    if (!(length(values) == 0 || missing(values) || is.null(values))) {
      n <- dim(values)[1]
      k <- dim(values)[2]
      if (n != nquad) 
        stop(paste("Number of rows in 'values' not equal to number of", 
                   "quadrature points."))
      if (k != nbasis) 
        stop(paste("Number of columns in 'values' not equal to number of", 
                   "basis functions."))
    }
    else values <- vector("list", 0)
    if (!(length(basisvalues) == 0 || missing(basisvalues) || 
          !is.null(basisvalues))) {
      if (!is.list(basisvalues)) 
        stop("BASISVALUES is not a list object.")
      sizevec <- dim(basisvalues)
      if (length(sizevec) != 2) 
        stop("BASISVALUES is not 2-dimensional.")
      for (i in 1:sizevec[1]) {
        if (length(basisvalues[[i, 1]]) != dim(basisvalues[[i, 
                                                            2]])[1]) 
          stop(paste("Number of argument values not equal number", 
                     "of values."))
      }
    }
    else basisvalues <- vector("list", 0)
    if (missing(dropind)) 
      dropind <- vector("list", 0)
    if (length(dropind) > 0) {
      ndrop = length(dropind)
      if (ndrop >= nbasis) 
        stop("Too many index values in DROPIND.")
      dropind = sort(dropind)
      if (ndrop > 1 && any(diff(dropind)) == 0) 
        stop("Multiple index values in DROPIND.")
      for (i in 1:ndrop) {
        if (dropind[i] < 1 || dropind[i] > nbasis) 
          stop("A DROPIND index value is out of range.")
      }
      nvalues = length(values)
      if (nvalues > 0 && length(values[[1]] > 0)) {
        for (ivalue in 1:nvalues) {
          derivvals = values[[ivalue]]
          derivvals = derivvals[, -dropind]
          values[[ivalue]] = derivvals
        }
      }
    }
    if (type == "fourier") {
      paramvec <- rangeval[2] - rangeval[1]
      period <- params[1]
      if (period <= 0) 
        stop("Period must be positive for (a Fourier basis")
      params <- period
      if ((2 * floor(nbasis/2)) == nbasis) 
        nbasis <- nbasis + 1
    }
    else if (type == "bspline") {
      if (!missing(params)) {
        nparams <- length(params)
        if (nparams > 0) {
          if (params[1] <= rangeval[1]) 
            stop("Smallest value in BREAKS not within RANGEVAL")
          if (params[nparams] >= rangeval[2]) 
            stop("Largest value in BREAKS not within RANGEVAL")
        }
      }
    }
    else if (type == "expon") {
      if (length(params) != nbasis) 
        stop("No. of parameters not equal to no. of basis fns for (exponential basisobj$")
    }
    else if (type == "polyg") {
      if (length(params) != nbasis) 
        stop("No. of parameters not equal to no. of basis fns for (polygonal basisobj$")
    }
    else if (type == "power") {
      if (length(params) != nbasis) 
        stop("No. of parameters not equal to no. of basis fns for (power basisobj$")
    }
    else if (type == "const") {
      params <- 0
    }
    else if (type == "monom") {
      if (length(params) != nbasis) 
        stop("No. of parameters not equal to no. of basis fns for (monomial basisobj$")
    }
    else stop("Unrecognizable basis")
    obj.call <- match.call()
    basisobj <- list(call = obj.call, type = type, rangeval = rangeval, 
                     nbasis = nbasis, params = params, dropind = dropind, 
                     quadvals = quadvals, values = values, basisvalues = basisvalues)
    oldClass(basisobj) <- "basisfd"
    basisobj
  }
  
  GBJ_pvalue <- function(observed_gbj, d, pairwise_cors, times_to_try=5)
  {
    if (observed_gbj<=0)
    {
      return ( list(GBJ_corp=1, eFlag=0) )
    } else
    {
      # Try to rerun multiple times if fails at first.
      times_tried <- 0
      eFlag <- 1
      while (times_tried<times_to_try & eFlag==1)
      {
        # Error flag, if we don't manage to throw this, then break out of while loop
        eFlag <- 0
        
        # Make gBJ_BB bounds
        GBJ_z_bounds <- rep(NA,d)
        
        # Sometimes starting too high gives an error
        if (times_tried == 2) {
          prev_bound <- 5
        } else if (times_tried == 3) {
          prev_bound <- 3
        } else if (times_tried == 4) {
          prev_bound <- 1
        } else {
          prev_bound <- 8.2
        }
        
        for ( kkk in 1:(ceiling(d/2)) )
        {
          temp_gbj <- tryCatch(uniroot(GBJ_objective, interval=c(0, prev_bound), d=d, k_vec=kkk, pairwise_cors=pairwise_cors, offset=observed_gbj), error=function(e) e, warning=function(w) w)
          
          # Sometimes, we can't go high enough in t, because pnorm,etc will just round to 0, and thus
          # the signs on both sides of the interval will be the same.  In this case, we will try again
          # a few times and then give up.
          if(length(class(temp_gbj))>1) {
            eFlag <- 1
            times_tried <- times_tried + 1
            break
          } else {
            GBJ_z_bounds[kkk] <- temp_gbj$root
          }
          prev_bound <- GBJ_z_bounds[kkk]
        }
      }
      
      # If eFlag still 1, then we tried multiple times times with no success.
      if(eFlag==1)
      {
        return ( list(GBJ_corp=NA, eFlag=eFlag) )
      }
      
      # Only need the first half
      # Make sure to sort in increasing order for crossprob_cor
      GBJ_z_bounds[(ceiling(d/2)+1):d] <- GBJ_z_bounds[ceiling(d/2)]
      GBJ_z_bounds <- sort(GBJ_z_bounds, decreasing=FALSE)
      
      # qnorm can't handle more precision than 10^-16
      # Also crossprob_cor can only handle Z up to 8.2
      GBJ_z_bounds[which(GBJ_z_bounds > 8.2)]= 8.2
      
      # Send it to the C++.
      if (sum(abs(pairwise_cors)) == 0) {
        # For the independence flag in the c++, just have to send a number < -1.
        GBJ_corp <- ebb_crossprob_cor_R(d=d, bounds=GBJ_z_bounds, correlations=rep(-999,2))
      } else {
        GBJ_corp <- ebb_crossprob_cor_R(d=d, bounds=GBJ_z_bounds, correlations=pairwise_cors)
      }
      
      return ( list(GBJ_corp=GBJ_corp, eFlag=eFlag) )
    }
  }
  
  GBJ_objective <- function (t_vec, d, k_vec = NULL, pairwise_cors, offset = 0) 
  {
    k_vec = NULL
    offset = 0
    t_vec <- sort(abs(t_vec), decreasing = TRUE)
    if (is.null(k_vec)) {
      k_vec <- 1:d
    }
    p_values <- 1 - pchisq(t_vec^2, df = 1)
    GBJ_indicator <- which(p_values < k_vec/d)
    first_half <- which(k_vec <= ceiling(d/2))
    non_zero <- intersect(GBJ_indicator, first_half)
    if (length(non_zero) == 0) {
      return(rep(0, length(t_vec)) - rep(offset, length(t_vec)))
    }
    mean_null <- rep(NA, length(t_vec))
    sigsq_null <- rep(NA, length(t_vec))
    mean_null[non_zero] <- 2 * d * surv(t_vec[non_zero])
    sigsq_null <- calc_var_nonzero_mu(d = d, t = t_vec, mu = rep(0, 
                                                                 length(t_vec)), pairwise_cors = pairwise_cors)
    mu_null <- rep(NA, length(t_vec))
    rho_null <- rep(NA, length(t_vec))
    gamma_null <- rep(NA, length(t_vec))
    mu_null[non_zero] <- mean_null[non_zero]/d
    rho_null[non_zero] <- (sigsq_null[non_zero] - d * mu_null[non_zero] * 
                             (1 - mu_null[non_zero]))/(d * (d - 1) * mu_null[non_zero] * 
                                                         (1 - mu_null[non_zero]))
    gamma_null[non_zero] = rho_null[non_zero]/(1 - rho_null[non_zero])
    Z_common_means <- rep(NA, length(t_vec))
    for (iii in 1:length(t_vec)) {
      if (iii %in% non_zero) {
        Z_common_means[iii] <- uniroot(qnorm_mu, lower = 0, 
                                       upper = 100, t = t_vec[iii], kkk = k_vec[iii], 
                                       d = d)$root
      }
    }
    alt_mean_vec <- k_vec
    alt_var_vec <- rep(NA, length(t_vec))
    alt_var_vec[non_zero] <- calc_var_nonzero_mu(d = d, t = t_vec[non_zero], 
                                                 mu = Z_common_means[non_zero], pairwise_cors = pairwise_cors)
    mu_alt <- rep(NA, length(t_vec))
    rho_alt <- rep(NA, length(t_vec))
    gamma_alt <- rep(NA, length(t_vec))
    mu_alt[non_zero] <- alt_mean_vec[non_zero]/d
    rho_alt[non_zero] <- (alt_var_vec[non_zero] - d * mu_alt[non_zero] * 
                            (1 - mu_alt[non_zero]))/(d * (d - 1) * mu_alt[non_zero] * 
                                                       (1 - mu_alt[non_zero]))
    gamma_alt[non_zero] = rho_alt[non_zero]/(1 - rho_alt[non_zero])
    pq_mat_null <- matrix(data = NA, nrow = length(t_vec), ncol = 2)
    pq_mat_alt <- matrix(data = NA, nrow = length(t_vec), ncol = 2)
    pq_mat_null[non_zero, 1] <- -mu_null[non_zero]/(d - 1)
    pq_mat_null[non_zero, 2] <- -(1 - mu_null[non_zero])/(d - 
                                                            1)
    gamma_check_null <- apply(pq_mat_null, 1, max)
    pq_mat_alt[non_zero, 1] <- -mu_alt[non_zero]/(d - 1)
    pq_mat_alt[non_zero, 2] <- -(1 - mu_alt[non_zero])/(d - 
                                                          1)
    gamma_check_alt <- apply(pq_mat_alt, 1, max)
    non_zero <- which(gamma_null >= gamma_check_null & gamma_alt >= 
                        gamma_check_alt)
    if (length(non_zero) == 0) {
      return(rep(0, length(t_vec)) - rep(offset, length(t_vec)))
    }
    null_loglik <- rep(0, length(t_vec))
    null_param_mat <- cbind(k_vec[non_zero], mu_null[non_zero], 
                            gamma_null[non_zero])
    null_loglik[non_zero] <- apply(null_param_mat, 1, ebb_loglik, 
                                   d = d)
    alt_loglik <- rep(0, length(t_vec))
    alt_param_mat <- cbind(k_vec[non_zero], mu_alt[non_zero], 
                           gamma_alt[non_zero])
    alt_loglik[non_zero] <- apply(alt_param_mat, 1, ebb_loglik, 
                                  d = d)
    BB_GBJ_stats <- alt_loglik - null_loglik
    return(BB_GBJ_stats - rep(offset, length(t_vec)))
  }
  
  TPM <- function (p, trunc = 0.2) 
{
  stopifnot((trunc > 0) & (trunc <= 1))
  stopifnot(is.vector(p))
  stopifnot((min(p) >= 0) & (max(p) <= 1))
  w <- prod(p^(p <= trunc))
  M <- length(p)
  if (w > trunc) {
    return(1)
  }
  else {
    pr <- 0
    for (k in 1:M) {
      s <- 0:(k - 1)
      term1 <- sum(w * (w <= (trunc^k)) * (((k * log(trunc)) - 
                                              log(w))^s)/factorial(s))
      term2 <- (trunc^k) * (w > (trunc^k))
      pr <- pr + (choose(M, k) * ((1 - trunc)^(M - k)) * 
                    (term1 + term2))
    }
    return(pr)
  }
}

FCP <- function (p) 
{
  y = log(p)
  y = -2 * sum(y)
  pr = pchisq(y, df = 2 * length(p), lower.tail = FALSE)
  return(pr)
}

Simes <- function(p)
{
  p = p[order(p,decreasing = F)]
  M = length(p)
  P_Simes = NULL
  for(i in 1:M){
    P_Simes = c(P_Simes,(M*p[i]/i))
  }
  return(min(P_Simes))
}

ExtSimes <- function (p, ldmat) 
{
  p <- p[order(p,decreasing = F)]
  eff.snpcount.fun <- function(ldmat) {
    ldmat <- as.matrix(ldmat)
    snpcount.local <- dim(ldmat)[1]
    if (snpcount.local <= 1) 
      return(1)
    ev <- eigen(ldmat, only.values = TRUE)$values
    if (sum(ev < 0) != 0) {
      ev <- ev[ev > 0]
      ev <- ev/sum(ev) * snpcount.local
    }
    ev <- ev[ev > 1]
    snpcount.local - sum(ev - 1)
  }
  eff.snpcount.global <- eff.snpcount.fun(ldmat)
  n_values <- length(p)
  candid <- sapply(1:n_values, function(i) {
    (eff.snpcount.global * p[i])/eff.snpcount.fun(ldmat[1:i, 
                                                        1:i])
  })
  p_ext_simes <- min(candid)
  return(p_ext_simes)
}

RTP <- function(p,ldmat, k=NULL) 
{
  Z <- qchisq(p,df = 1,lower = F)
  M <- length(Z)
  if(is.null(k)) k <- round(0.5*M)
  ## transform Z-scores to chi-squares
  y <- Z^2
  ## calculate TQ-statistic
  yy <- y %*% rep(1, M)
  ## eigen decomposition of the LD matrix
  ee <- eigen(ldmat); eivec <- ee$vectors; eigva <- ee$values
  pc <- eivec %*% diag(sqrt(1/eigva)) %*% t(eivec)
  ## calculate decorreated statistics
  x <- (Z %*% pc)^2
  k <- round(0.5*M)
  #k <- M
  px <- 1-pchisq(x, df=1)
  P <- sort(px)
  Z <- sum(-log(P[1:k]))
  p.RTP <- integrate(function(x,y,m,n) 1-pgamma(log(qbeta(x,m+1,n-m))*m+y,m),0,1,Z,k,M)$va
  return(p.RTP)
}

SimpleM <- function(p,ldmat, pca_cut_perc = 0.995)
{
  if (is.positive.definite(ldmat) == FALSE) {
    ldmat <- make.positive.definite(ldmat)
  }
  if (is.positive.definite(ldmat) == FALSE) {
    ldmat <- cor(as.matrix(ref.data), use = "p")
    diag(ref.data) <- 1.0001
  }
  if (is.positive.definite(ldmat) == FALSE) {
    diag(ldmat) <- 1.001
  }
  if (is.positive.definite(ldmat) == FALSE) {
    diag(ldmat) <- 1.01
  }
  if (is.positive.definite(ldmat) == FALSE) 
    stop("ldmat is not positive definite. Please re-calculate with ld.Rsquare function.\n")
  min_p_obs <- min(p)
  num_of_snps <- length(p)
  cor_r <- ldmat
  eigen_values <- eigen(cor_r, only.values = TRUE)$values
  eigen_values_sorted <- sort(eigen_values, decreasing = TRUE)
  sum_eigen_values <- sum(eigen_values_sorted)
  M_eff_G <- 1
  for (k in 1:num_of_snps) {
    temp <- sum(eigen_values_sorted[1:k])/sum_eigen_values
    if (temp >= pca_cut_perc) {
      M_eff_G <- k
      break
    }
  }
  return(1 - (1 - min_p_obs)^M_eff_G)
}

GM <- function(p,a=0.0383)
{
  M <- length(p)
  Y <- sum(qgamma(p=1-p,shape = a,scale = 1))
  p <- 1-pgamma(Y,shape = M*a,scale = 1)
}

ARTP <- function (p, ldmat, k = NULL, w = NULL, ...) 
{
  tol.cor <- sqrt(.Machine$double.eps)
  ldmat <- abs(ldmat)
  ldmat[upper.tri(ldmat, TRUE)] <- 0
  m <- apply(ldmat, 2, function(.) all(. < 1 - tol.cor))
  M <- sum(m)
  P <- DOT(p,ldmat)
  P <- sort(P)
  if (is.null(k)) 
    k <- round(M * 0.5)
  k <- min(k, L)
  if (is.null(w)) 
    w <- rep(1, k)
  w <- w[1:k]
  z <- ((1 - P)/c(1, 1 - P[-M]))^(M:1)
  p <- (1 - z)[1:k]
  q <- qnorm(p) * w
  if (max(q) == Inf) 
    return(1)
  sumQ <- cumsum(qnorm(p) * w)
  pSg <- matrix(cumsum(w^2), k, k)
  pSg[lower.tri(pSg)] <- t(pSg)[lower.tri(pSg)]
  pCr <- cov2cor(pSg)
  sQ <- sumQ/sqrt(diag(pSg))
  Y <- max(sQ)
  P <- try(mvtnorm::pmvnorm(lower = rep(-Inf, k), upper = rep(max(sQ), 
                                                              k), sigma = pCr))
  if (inherits(P, "try-error")) {
    print("bad")
    P <- 1
  }
  else P <- P[1]
  return(P)
}

DOT <- function (p, ldmat, tol.cor = NULL, tol.egv = NULL, ...) 
{
  Z <- qchisq(p,df = 1,lower = F)
  d <- nsp(ldmat, eps = tol.egv, ...)
  H <- d$H
  L <- d$L
  X <- H %*% Z
  ssq <- sum(X^2)
  p.DOT <- 1 - pchisq(ssq, df=L)
  return(p.DOT)
}

TQ <- function(p, ldmat, approx.method="chi-square")
{
  ## the number of statistics
  Z <- qchisq(p,df = 1,lower = F)
  n <- length(Z)
  S <- ldmat
  ## transform Z-scores to chi-squares
  y <- Z^2
  ## calculate TQ-statistic
  yy <- y %*% rep(1, n)
  ## eigen decomposition of the LD matrix
  ee <- eigen(S); eivec <- ee$vectors; eigva <- ee$values
  pc <- eivec %*% diag(sqrt(1/eigva)) %*% t(eivec)
  ## calculate decorreated statistics
  x <- (Z %*% pc)^2
  ## calculate DOT-statistic
  xx <- x %*% rep(1, n)
  ## calculate DOT P-value
  DOT <- 1 - pchisq(as.numeric(xx), df = n, ncp=0)
  ## calculate TQ O-value
  TQ <- CompQuadForm::imhof(yy, lambda = eigva, delta = rep(0, n), epsrel = 1e-11, limit = 2e5)$Qq
  ## calculate average pair-wise LD
  rho.av2 <- sqrt(mean(S[lower.tri(S)]^2))
  ## two way to approaximate TQ P-value
  if(approx.method=="chi-square"){
    P.TQ <- 1-pchisq( (yy - (1 - rho.av2)*(n - 1) ) / (rho.av2*(n - 1) + 1), df=1)
  }else if(approx.method=="equicorrelation"){
    Sgm.e <- (1-rho.av2)*diag(1,n,n) + rho.av2 * matrix(1,n,n)
    ev <- eigen(Sgm.e)$values
    P.TQ <- imhof(yy, lambda = eigen(Sgm.e)$values, delta = rep(0, n), epsrel = 1e-11, limit = 2e5)$Qq
  }
  return(P.TQ)
}

GATES <- function(p, ldmat, keyGeneLoc=FALSE)
{
  pval_sort <- sort(p)
  pval_order <- order(p)
  n_snps <- length(p)
  cor_G <- LD_R(ldmat)
  cor_P <- cor_G[pval_order, pval_order]
  cor_P <- 0.2982 * cor_P^6 - 0.0127 * cor_P^5 + 0.0588 *  cor_P^4 +
    0.0099 * cor_P^3 + 0.6281 * cor_P^2 - 9e-04 * cor_P
  eff.snpcount.fun <- function(ldmat) {
    ldmat <- as.matrix(ldmat)
    snpcount.local <- dim(ldmat)[1]
    if (snpcount.local <= 1) 
      return(1)
    ev <- eigen(ldmat, only.values = TRUE)$values
    if (sum(ev < 0) != 0) {
      ev <- ev[ev > 0]
      ev <- ev/sum(ev) * snpcount.local
    }
    ev <- ev[ev > 1]
    snpcount.local - sum(ev - 1)
  }
  eff.snpcount.global <- eff.snpcount.fun(cor_P)
  n_values <- length(pval_sort)
  candid <- sapply(1:n_values, function(i)(eff.snpcount.global * pval_sort[i])/eff.snpcount.fun(cor_P[1:i,1:i]))
  if(keyGeneLoc==FALSE){
    return(min(candid))
  }else{
    Pg <- min(candid)
    keyGloc <- which(min(candid) == candid)[1]
    out <- c(Pg, keyGloc)
    names(out) <- c("Pg", "keyGeneLoc")
    return(out)
  }
}

VEGAS <- function(p, ldmat, vegas.pct = c(0.1, 0.2, 0.3, 0.4, 1),seed = 12345)
{
  set.seed(seed)
  x = p
  cor_G <- LD_R(ldmat)
  vegas.call <- function (x, cor_G, vegas.pct, n_simul) 
  {
    stopifnot(length(x) == ncol(cor_G))
    vegas_vec <- ceiling(vegas.pct * ncol(cor_G))
    vegas_vec <- sort(vegas_vec)
    if (vegas_vec[1] > 1) {
      vegas.pct <- c(0, vegas.pct)
      vegas_vec <- c(1, vegas_vec)
    }
    chisq_vec <- qchisq(x, 1, lower.tail = FALSE)
    chisq_vec[chisq_vec == Inf] <- 60
    n_snps <- length(x)
    n_tests <- length(vegas_vec)
    TS_obs <- rep(NA, n_tests)
    TS_obs[1] <- max(chisq_vec, na.rm = TRUE)
    chisq_vec <- sort(chisq_vec, decreasing = TRUE)
    for (j in 2:n_tests) TS_obs[j] <- sum(chisq_vec[1:vegas_vec[j]])
    rd <- rmvnorm(n_simul, mean = rep(0, n_snps), sigma = cor_G)
    rd2 <- rd^2
    rd2 <- apply(rd2, 1, sort, decreasing = TRUE)
    pPerm0 <- rep(NA, n_tests)
    T0s <- apply(rd2, 2, max)
    pPerm0[1] <- (sum(T0s >= TS_obs[1]) + 1)/(length(T0s) + 
                                                1)
    for (j in 2:n_tests) {
      for (i in 1:n_simul) T0s[i] <- sum(rd2[1:vegas_vec[j], 
                                             i])
      pPerm0[j] <- (sum(T0s >= TS_obs[j]) + 1)/(length(T0s) + 
                                                  1)
    }
    v1 <- paste0("VEGAS.p", vegas.pct)
    v1[vegas_vec == ncol(cor_G)] <- "VEGAS.all"
    v1[vegas_vec == 1] <- "VEGAS.max"
    names(pPerm0) <- v1
    pPerm0
  }
  if (is.positive.definite(cor_G) == FALSE) 
    stop("cor_G is not positive definite. Please re-calculate with ld.Rsqure function.\n")
  pval_vegas <- vegas.call(x = x, cor_G = cor_G, vegas.pct = vegas.pct, 
                           n_simul = 1000)
  if (any(pval_vegas <= 0.005)) {
    pval_vegas <- vegas.call(x = x, cor_G = cor_G, vegas.pct = vegas.pct, 
                             n_simul = 10000)
  }
  if (any(pval_vegas <= 5e-04)) {
    pval_vegas <- vegas.call(x = x, cor_G = cor_G, vegas.pct = vegas.pct, 
                             n_simul = 1e+05)
  }
  if (any(pval_vegas <= 5e-05)) {
    pval_vegas <- vegas.call(x = x, cor_G = cor_G, vegas.pct = vegas.pct, 
                             n_simul = 1e+06)
  }
  if (any(pval_vegas <= 5e-06) && max.simulation > 1e+06) {
    pval_vegas <- vegas.call(x = x, cor_G = cor_G, vegas.pct = vegas.pct, 
                             n_simul = max.simulation)
  }
  return(pval_vegas)
}

aSPUs <- function(p, ldmat, pow = c(1:8, Inf), n.perm = 1000, Ps = FALSE)
{
  Z <- qchisq(p,df = 1,lower = F)
  n <- length(Z)
  k <- n
  if (Ps == TRUE) 
    Z <- qnorm(1 - Z/2)
  U <- Z
  CovS <- ldmat
  eS <- eigen(CovS, symmetric = TRUE)
  ev <- eS$values
  CovSsqrt <- eS$vectors %*% diag(sqrt(pmax(ev, 0)), k)
  Ts = rep(NA, length(pow))
  for (j in 1:length(pow)) {
    if (pow[j] < Inf) 
      Ts[j] = sum(U^pow[j])
    else Ts[j] = max(abs(U))
  }
  pPerm0 = rep(NA, length(pow))
  T0s = numeric(n.perm)
  s <- sample(1:10^5, 1)
  for (j in 1:length(pow)) {
    set.seed(s)
    for (b in 1:n.perm) {
      U00 <- rnorm(k, 0, 1)
      U0 <- CovSsqrt %*% U00
      if (Ps == TRUE) 
        U0 <- abs(U0)
      if (pow[j] < Inf) {
        T0s[b] = round(sum(U0^pow[j]), digits = 8)
      }
      if (pow[j] == Inf) {
        T0s[b] = round(max(abs(U0)), digits = 8)
      }
    }
    pPerm0[j] = sum(abs(Ts[j]) <= abs(T0s))/n.perm
    P0s = ((n.perm - rank(abs(T0s))) + 1)/(n.perm)
    if (j == 1) 
      minp0 = P0s
    else minp0[which(minp0 > P0s)] = P0s[which(minp0 > P0s)]
  }
  Paspu <- (sum(minp0 <= min(pPerm0)) + 1)/(n.perm + 1)
  pvs <- c(pPerm0, Paspu)
  Ts <- c(Ts, min(pPerm0))
  return(pvs[10])
}

MLR <- function(p,ldmat,n)
{
  Z <- qchisq(p,df = 1,lower = F)
  m = length(Z)
  U05 = ldmat %^% (-.5)
  U05Z = U05 %*% Z        # U^-.5 %*% Z
  R2 <- sum(U05Z^2) / n       # assimp m/n
  Fstat = ((n - m) / m) * R2 / (1-R2)  # F-statistic
  p.MLR = as.double(pf(Fstat, m, n - m, lower.tail = FALSE)) # F-test
  return(p.MLR)
}

FLM <- function(p,ldmat,POS,w=NULL,n)
{
  basis <- basisfd(type = "fourier", rangeval = c(0, 1), nbasis = 25, 
                   params = 1, dropind = NULL, quadvals = NULL, 
                   values = NULL, basisvalues = NULL)
  basis$names <- c("const",as.vector(outer(c("sin", "cos"),1:12, paste, sep = "")))
  
  Z <- qchisq(p,df = 1,lower = F)
  M <- length(Z)
  if (is.null(w))  w <- rep(1, M)
  B <- eval.basis(POS, basis)  # as matrix (m x Kb)
  WB <- B * w
  Mat <- t(WB) %*% ldmat %*% WB
  BUB_1 <- Mat %^% (-1)
  
  BZstat <- as.vector(t(WB) %*% Z)
  RSS <- n - sum(BZstat * (BUB_1 %*% BZstat))
  
  Fstat <- ((n - 25) / 25) * (n - RSS) / RSS   # F-statistic
  p <- pf(Fstat, 25, n - 25, lower.tail = FALSE)
  
  return(p)
}

PCA <- function(p,ldmat,w=NULL,n) 
{
  Z <- qchisq(p,df = 1,lower = F)
  if (is.null(w))  w <- rep(1, M)
  WZ  <- as.vector(w * Z) * sqrt(n)
  WUW <- as.matrix(t(ldmat * w) * w) * n
  eX1 <- eigen(WUW, symmetric = TRUE)
  pCA <- with (eX1, {
    values[values < 0] <- 0
    prop.var <- values / sum(values)
    cum.var <- cumsum(prop.var)
    list(scores = vectors, importance = cum.var)
  })
  CPV <- pCA$importance   # Cumulative Proportion of Variance (CPV)
  M <- min(which(CPV >= 0.85))  # components for which Explained variance fraction is about 85%
  BBB <- as.matrix(pCA$scores[,1:M])
  GY <- as.vector(t(BBB) %*% WZ)
  CC <- as.matrix(t(BBB) %*% WUW %*% BBB)
  m <- qr(BBB)$rank
  if (m > 1) {
    RSS <- (n - sum(GY * as.vector((CC %^% (-1)) %*% GY)))
  } else { RSS <- (n - GY * GY / CC) }
  Fstat <- ((n - m) / m) * (n - RSS) / RSS    # F-statistic
  p <- pf(Fstat, m, n - m, lower.tail = FALSE)
  minP <- 100;
  minM <- 0
  if (p < minP) minM <- M
  minP <- min(p, minP)
  c(minP)
}

HYST <- function (p, ldmat) 
{
  M <- length(p)
  o.pv <- order(p)
  ldmat2 <- ldmat[o.pv, o.pv]
  out <- GATES(ldmat = ldmat2, p = sort(p),keyGeneLoc=TRUE)
  PGs <- out[1]
  keyGs <- out[2]
  Hyst <- -2 * sum(log(PGs))
  p.hyst = 1 - pchisq(Hyst, df = 2)
  return(p.hyst)
}

COMBAT <- function (p, ldmat, vegas.pct = c(0.1, 0.2, 0.3, 0.4, 1), 
                    pca_cut_perc = 0.995, nperm = 100, seed = 12345)
{
  pvalues <- as.numeric(p)
  n_snps <- length(pvalues)
  cor_G <- LD_R(ldmat)
  set.seed(seed)
  pval_gates <- GATES(p, ldmat)
  pval_vegas <- VEGAS(p, ldmat)
  pval_simpleM <- SimpleM(p, ldmat)
  gene_pvals <- c(GATES = pval_gates, pval_vegas, simpleM = pval_simpleM)
  rd <- rmvnorm(nperm, mean = rep(0, n_snps), sigma = cor_G)
  rd2 <- rd^2
  simul_pval_mat <- pchisq(rd2, 1, lower.tail = FALSE)
  func1 = function(p,ldmat) {
    p_gates <- GATES(p, ldmat)
    p_vegas <- VEGAS(p, ldmat)
    p_simpleM <- SimpleM(p, ldmat)
    c(p_gates, p_vegas, p_simpleM)
  }
  gene_pval_mat = apply(simul_pval_mat, 1, func1, ldmat = ldmat)
  gene_pval_mat = t(gene_pval_mat)
  method_cor <- cor(gene_pval_mat)
  order_pvals <- order(gene_pvals)
  sort_pvals <- gene_pvals[order_pvals]
  method_cor <- method_cor[order_pvals, order_pvals]
  p_combat_simes <- ext_simes(sort_pvals, method_cor)
  return(p_combat_simes)
}

BT <- function (p, ldmat, W=NULL) 
{
  M = length(p)
  Z <- qchisq(p,df = 1,lower = F)
  if (is.null(W)) W = rep(1, M)
  Zw = Z * W
  Rw = t(ldmat * W) * W
  eR = eigen(Rw, sym = TRUE)
  lamR = abs(eR$val)
  eta = colSums(eR$vec) * sqrt(lamR)
  R1 = sum(eta^2)
  Qb = sum(Zw)^2
  p.BT = pchisq(Qb/R1, 1, lower = FALSE)
  return(p.BT)
}

SKAT <- function(p, ldmat, W = NULL, rho = c((0:5/10)^2, 0.5, 1)) 
{
  M = length(p)
  Z <- qchisq(p,df = 1,lower = F)
  if (is.null(W))  W = rep(1, M)
  Zw = Z * W
  Rw = t(R * W) * W
  eR = eigen(Rw, sym = TRUE)
  lamR = abs(eR$val)
  Qv = sum(Zw^2)
  p.SKAT = KATpval(Qv, lamR)
  return(p.SKAT)
}

SKATO <- function(p, ldmat, W = NULL, rho = c((0:5/10)^2, 0.5, 1))
{
  M = length(p)
  Z <- qchisq(p,df = 1,lower = F)
  L = length(rho)
  if (L <= 2) return(NA)
  if (is.null(W))
    W = rep(1, M)
  Zw = Z * W
  Rw = t(ldmat * W) * W
  eR = eigen(Rw, sym = TRUE)
  lamR = abs(eR$val)
  eta = colSums(eR$vec) * sqrt(lamR)
  R1 = sum(eta^2)
  R2 = sum(eta^2 * lamR)
  c2 = outer(eta, eta)
  Lamq = eigen(diag(lamR) - R2/R1^2 * c2, symmetric = TRUE, 
               only.values = TRUE)$val
  Qb = sum(Zw)^2
  pvalb = pchisq(Qb/R1, 1, lower = FALSE)
  Qv = sum(Zw^2)
  pvalv = KATpval(Qv, lamR)
  
  L1 = L - 1
  rho1 = rho[-L]
  Qw = (1 - rho) * Qv + rho * Qb
  pval = rep(1, L)
  pval[1] = pvalv
  pval[L] = pvalb
  Lamk = vector("list", L)
  Lamk[[L]] = R1
  Lamk[[1]] = lamR
  for (k in 2:L1) {
    mk = rho[k] * c2
    diag(mk) = diag(mk) + (1 - rho[k]) * lamR
    aak = zapsmall(abs(eigen(mk, sym = TRUE, only.val = TRUE)$val))
    Lamk[[k]] = aak[aak > 0]
    pval[k] = KATpval(Qw[k], Lamk[[k]])
  }
  minP = min(pval)
  L = length(rho)
  qval = rep(0, L1)
  for (k in 1:L1) qval[k] = Liu0.qval(minP, Lamk[[k]])
  q1 = qchisq(minP, 1, lower = FALSE)
  tauk = (1 - rho1) * R2/R1 + rho1 * R1
  katint = function(xpar) {
    eta1 = sapply(xpar, function(eta0) min((qval - tauk * 
                                              eta0)/(1 - rho1)))
    KATpval(eta1, Lamq) * dchisq(xpar, 1)
  }
  prec = 1e-04
  p.value = try({
    minP + integrate(katint, 0, q1, subdivisions = 1000, 
                     abs.tol = minP * prec)$val
  }, silent = TRUE)
  while (class(p.value) == "try-error") {
    prec = prec * 2
    p.value = try({
      minP + integrate(katint, 0, q1, abs.tol = minP * 
                         prec)$val
    }, silent = TRUE)
  }
  p.value = min(p.value, minP * L)
  return(list(p.value = c(A = p.value, S2 = pvalv, S = pvalb), 
              pval = pval, rho.est = rho[which.min(pval)]))
}

ART <- function(p,ldmat) 
{
  Z <- qchisq(p,df = 1,lower = F)
  M <- length(Z)
  ## transform Z-scores to chi-squares
  y <- Z^2
  ## calculate TQ-statistic
  yy <- y %*% rep(1, M)
  ## eigen decomposition of the LD matrix
  ee <- eigen(ldmat); eivec <- ee$vectors; eigva <- ee$values
  pc <- eivec %*% diag(sqrt(1/eigva)) %*% t(eivec)
  ## calculate decorreated statistics
  x <- (Z %*% pc)^2
  k <- round(0.5*M)
  #k <- M
  px <- 1-pchisq(x, df=1)
  P <- sort(px)
  Lw <- sum(log(P[1:(k-1)]))
  Pk <- P[k]
  d = (k-1)*(digamma(L+1) - digamma(k))
  ak = (k-1)*log(Pk) - lW + qgamma(1-pbeta(Pk, k, L-k+1), shape=d)
  1 - pgamma(ak, shape=k+d-1)
}

ART.A <- function(p,ldmat)
{
  Z <- qchisq(p,df = 1,lower = F)
  M <- length(Z)
  ## transform Z-scores to chi-squares
  y <- Z^2
  ## calculate TQ-statistic
  yy <- y %*% rep(1, M)
  ## eigen decomposition of the LD matrix
  ee <- eigen(ldmat); eivec <- ee$vectors; eigva <- ee$values
  pc <- eivec %*% diag(sqrt(1/eigva)) %*% t(eivec)
  ## calculate decorreated statistics
  x <- (Z %*% pc)^2
  k <- round(0.5*M)
  #k <- M
  px <- 1-pchisq(x, df=1)
  P <- sort(px)
  wgt <- rep(1,k)
  z <- P
  z[1] <- ( 1 - P[1] )^L
  for(j in 2:L) z[j] <- ((1-P[j]) / (1-P[j-1]))^((L-(j-1)))
  p <- (1-z)[1:k]
  k = length(p)
  sumZ <- rep(0, k)
  y <- qnorm(p)
  z <- y
  gz <- z[1] * wgt[1]
  sumZ[1] <- gz
  for(i in 2:k) {
    gz <- p[ 1 : i ]
    for(j in 1:i) gz[j] <- z[j] * wgt[j]
    sumZ[i] <- sum(gz)
  }
  Lo = diag(k); Lo[lower.tri(Lo)] <- 1
  pSg <- Lo %*% diag(wgt[1:k]^2) %*% t(Lo)
  pCr <- cov2cor(pSg)
  sZ <- sumZ
  for(i in 1:k) {
    sZ[i] <- sumZ[i] / sqrt(diag(pSg))[i]
  }
  p.ARTA <- pmvnorm(lower = rep(-Inf,k), upper = rep(max(sZ), k), sigma = pCr)[1]
  return(p.ARTA)
}

ACAT<-function(p,Weights=NULL)
{
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(Weights)){
    Weights<-rep(1/length(p),length(p))
  }else{
    Weights<-Weights/sum(Weights)
  }
  
  #### check if there are very small non-zero p values
  is.small<-(p<1e-16)
  if (sum(is.small)==0){
    cct.stat<-sum(Weights*tan((0.5-p)*pi))
  }else{
    cct.stat<-sum((Weights[is.small]/p[is.small])/pi)
    cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-p[!is.small])*pi))
  }
  #### check if the test statistic is very large.
  if (cct.stat>1e+15){
    pval<-(1/cct.stat)/pi
  }else{
    pval<-1-pcauchy(cct.stat)
  }
  return(pval)
}

ACATO <- function(p)
{
  if (all(is.na(p))) return(NA)
  p <- p[!is.na(p)]
  #### check if there are very small non-zero p values
  is.small <- (p < 1e-16)
  if (sum(is.small) == 0) {
    cct.stat <- sum(tan((0.5 - p) * pi))/length(p)
  } else {
    cct.stat <- sum((1 / p[is.small]) / pi)
    cct.stat <- cct.stat + sum(tan((0.5 - p[!is.small]) * pi))
    cct.stat <- cct.stat/length(p)
  }
  #### check if the test statistic is very large.
  if (cct.stat > 1e+15){
    pval <- (1 / cct.stat) / pi
  } else {
    pval <- 1 - pcauchy(cct.stat)
  }
  pval
}

BJ <- function(p, ldmat=NULL) 
{
  library(Rcpp)
  sourceCpp("E:/ebb.cpp")
  
  Z <- qchisq(p,df = 1,lower = F)
  # Parse inputs, do some error checking.
  if (length(Z) > 2000) {
    stop("You have too many factors, please restrict to 2000")
  }
  t_vec <- sort(abs(Z), decreasing = TRUE)
  too_big <- which(t_vec > 8.2)
  if (length(too_big) > 0) {
    t_vec[too_big] <- 8.2
  }
  pairwise_cors <- ldmat[upper.tri(ldmat)]
  d <- length(t_vec)
  
  # Check for qualifying p-values under the null (the indicator part of the GBJ statistic)
  # and also that we are only considering 'first half' p-values
  p_values <- 1-pchisq(t_vec^2, df=1)
  BJ_indicator <- which( p_values < (1:d)/d )
  first_half <- 1:(ceiling(d/2))
  non_zero <- intersect(BJ_indicator, first_half)
  
  # If no indicies qualified, stop
  if (length(non_zero) == 0) {
    return ( p.BJ=1 )
  }
  
  #################################
  # Some indicies passed
  i_vec <- 1:d
  BJ_stats <- rep(0, d)
  BJ_stats[non_zero] <- i_vec[non_zero] * log(i_vec[non_zero]/(d*p_values[non_zero])) +
    (d-i_vec[non_zero]) * log((1-i_vec[non_zero]/d)/(1-p_values[non_zero]))
  BJ_stats[d] <- 0
  
  # Observed BJ statistic
  b <- max(BJ_stats[1:(ceiling(d/2))])
  
  # Calculate p-value
  if (b<=0) {
    return ( p.BJ=1 )
  }
  
  # BJ bounds
  BJ_p_bounds <- rep(NA, d)
  
  # Use uniroot to find the pvalue bounds.
  for ( jjj in 1:(ceiling(d/2)) ) {
    BJ_p_bounds[jjj] <- uniroot(f=function(x, k, d, b) {k*log(k/(d*x)) +
        (d-k)*log((1-k/d)/(1-x)) - b}, k=jjj, d=d,
        b=b, lower=0, upper=jjj/d, tol=(10^(-12)))$root
  }
  
  # The last half of the order statistic bounds
  BJ_p_bounds[(ceiling(d/2)+1):d] <- BJ_p_bounds[ceiling(d/2)]
  
  # Now put the bounds in terms of the Z statistics
  BJ_z_bounds <- qnorm(1 - BJ_p_bounds/2)
  BJ_z_bounds <- sort(BJ_z_bounds, decreasing=F)
  
  # qnorm can't handle more precision than 10^-16
  # Also crossprob_cor can only handle Z up to 8.2
  BJ_z_bounds[which(BJ_z_bounds > 8.2)]= 8.2
  
  # Send it to the C++.
  if (sum(abs(pairwise_cors)) == 0) {
    # For the independence flag in the c++, just have to send a number < -1.
    p.BJ <- ebb_crossprob_cor_R(d=d, bounds=BJ_z_bounds, correlations=rep(-999,2))
  } else {
    p.BJ <- ebb_crossprob_cor_R(d=d, bounds=BJ_z_bounds, correlations=pairwise_cors)
  }
  
  return (p.BJ)
}

GBJ <- function(p, ldmat=NULL)
{
  Z <- qchisq(p,df = 1,lower = F)
  # Parse inputs, do some error checking.
  if (length(Z) > 2000) {
    stop("You have too many factors, please restrict to 2000")
  }
  t_vec <- sort(abs(Z), decreasing = TRUE)
  too_big <- which(t_vec > 8.2)
  if (length(too_big) > 0) {
    t_vec[too_big] <- 8.2
  }
  pairwise_cors <- ldmat[upper.tri(ldmat)]
  d <- length(t_vec)
  
  # Move to BJ if no correlation at all
  if (sum(abs(pairwise_cors)) == 0) {
    p.BJ <- BJ(Z=t_vec, pairwise_cors=pairwise_cors)
    return (p.BJ)
  }
  
  # Calculate the observed GBJ statistic
  
  GBJ_stats <- GBJ_objective(t_vec=t_vec, d=d, pairwise_cors=pairwise_cors)
  gbj <- max(GBJ_stats)
  
  # Calculate p_value
  GBJ_p_list <- GBJ_pvalue(observed_gbj=gbj, d=d, pairwise_cors=pairwise_cors)
  p.GBJ=GBJ_p_list$GBJ_corp
  
  # If no NA in the p-value, then everything success and return GBJ output.
  return ( p.GBJ )
}

HC <- function(p, ldmat=NULL) 
{
  library(Rcpp)
  sourceCpp("E:/ebb.cpp")
  
  Z <- qchisq(p,df = 1,lower = F)
  # Parse inputs, do some error checking.
  if (length(Z) > 2000) {
    stop("You have too many factors, please restrict to 2000")
  }
  t_vec <- sort(abs(Z), decreasing = TRUE)
  too_big <- which(t_vec > 8.2)
  if (length(too_big) > 0) {
    t_vec[too_big] <- 8.2
  }
  pairwise_cors <- ldmat[upper.tri(ldmat)]
  d <- length(t_vec)
  
  # Calculate HC objectives
  p_values <- 1-pchisq(t_vec^2, df=1)
  i_vec <- 1:d
  HC_stats <- sqrt(d) * (i_vec/d - p_values) / sqrt(p_values*(1-p_values))
  
  # Observed HC statistic
  h <- max(HC_stats, na.rm=TRUE)
  
  # Calculate p-value
  if (h<=0) {
    return ( p.HC=1)
  }
  
  # BJ bounds
  HC_p_bounds <- rep(NA, d)
  
  # Explicit inverse of HC to find the p-value bounds
  HC_p_bounds <- ((2*i_vec+h^2)/d - sqrt((2*i_vec/d+h^2/d)^2 - 4*i_vec^2/d^2 - 4*i_vec^2*h^2/d^3))/(2*(1+h^2/d))
  HC_z_bounds <- qnorm(1-HC_p_bounds/2)
  HC_z_bounds <- sort(HC_z_bounds, decreasing=F)
  
  # qnorm can't handle more precision than 10^-16
  # Also crossprob_cor can only handle Z up to 8.2
  HC_z_bounds[which(HC_z_bounds > 8.2)]= 8.2
  
  # Send it to the C++.
  if (sum(abs(pairwise_cors)) == 0) {
    # For the independence flag in the c++, just have to send a number < -1.
    p.HC <- ebb_crossprob_cor_R(d=d, bounds=HC_z_bounds, correlations=rep(-999,2))
  } else {
    p.HC <- ebb_crossprob_cor_R(d=d, bounds=HC_z_bounds, correlations=pairwise_cors)
  }
  return ( p.HC)
}

GHC <- function(p, ldmat=NULL)
{
  Z <- qchisq(p,df = 1,lower = F)
  # Parse inputs, do some error checking.
  if (length(Z) > 2000) {
    stop("You have too many factors, please restrict to 2000")
  }
  t_vec <- sort(abs(Z), decreasing = TRUE)
  too_big <- which(t_vec > 8.2)
  if (length(too_big) > 0) {
    t_vec[too_big] <- 8.2
  }
  pairwise_cors <- ldmat[upper.tri(ldmat)]
  d <- length(t_vec)
  
  # Move to hC if no correlation at all
  if (sum(abs(pairwise_cors)) == 0) {
    p.GHC <- HC(test_stats=t_vec, pairwise_cors=pairwise_cors)
    return ( p.GHC )
  }
  
  # Calculate GHC objectives
  i_vec <- 1:d
  p_values <- 1-pchisq(t_vec^2, df=1)
  GHC_stats <- (i_vec - d*p_values) / sqrt(calc_var_nonzero_mu(d=d, t=t_vec, mu=0,
                                                               pairwise_cors=pairwise_cors))
  
  # Observed GHC statistic - sometimes a Z-statistic is 0 and so we get NA for variance
  ghc <- max(GHC_stats, na.rm=TRUE)
  
  # Calculate p-value
  if (ghc <= 0) {
    return (1)
  }
  
  # GHC bounds
  GHC_p_bounds <- rep(NA, d)
  
  # increase tolerance of uniroot for large ghc
  if(ghc>10) {
    my_tol <- (-100)
  } else {my_tol <- (-12)}
  
  # Use uniroot to find the pvalue bounds.
  GHC_lowerbound <- 10^(-20)
  for(kkk in 1:d) {
    # Sometimes run into precision errors, need to think about how
    # we can fix this so don't need the tryCatch
    temp_ghc <- tryCatch(uniroot(GHC_objective, k=kkk, d=d, offset=ghc,
                                 pairwise_cors=pairwise_cors, lower=GHC_lowerbound, upper=(1-10^(-12)),
                                 tol=(10^(my_tol))), error=function(e) e, warning=function(w) w)
    
    # If it doesn't work, just run the HC for them
    if(length(class(temp_ghc))>1) {
      if (ghc >= 200000) {
        HC_output <- HC(test_stats=t_vec, pairwise_cors=pairwise_cors)
        GHC_err_code <- '1: Pvalue likely less than 10^(-12), R/C++ not enough precision. Returning standard Higher Criticism test instead.'
        return ( HC_output$HC_pvalue )
      }
      
      # If evidence of underdispersion, again give them BJ p-value
      else if (sum(pairwise_cors) < 0) {
        HC_output <- HC(test_stats=t_vec, pairwise_cors=pairwise_cors)
        GHC_err_code <- '2: Error in numerical routines. Many apologies, please report to developer! Returning standard Higher Criticism test instead.'
        return ( HC_output$HC_pvalue )
      }
      
      # Any other errors, give them BJ p-value
      else {
        HC_output <- HC(test_stats=t_vec, pairwise_cors=pairwise_cors)
        GHC_err_code <- '3: Unknown error. Many apologies, please report to developer! Returning standard Higher Criticism test instead.'
        return ( HC_output$HC_pvalue )
      }
    }
    
    # It worked, keep going
    GHC_p_bounds[kkk] <- temp_ghc$root
    
    # small security measure to ensure that GHC bounds are increasing
    GHC_lowerbound <- GHC_p_bounds[kkk]
  }
  
  # now put the bounds in terms of the Z statistics
  GHC_z_bounds <- qnorm(1-GHC_p_bounds/2)
  GHC_z_bounds <- sort(GHC_z_bounds, decreasing=F)
  
  # qnorm can't handle more precision than 10^-16
  # Also crossprob_cor can only handle Z up to 8.2
  GHC_z_bounds[which(GHC_z_bounds > 8.2)]= 8.2
  
  # Send it to the C++.
  if (sum(abs(pairwise_cors)) == 0) {
    # For the independence flag in the c++, just have to send a number < -1.
    p.GHC <- ebb_crossprob_cor_R(d=d, bounds=GHC_z_bounds, correlations=rep(-999,2))
  } else {
    p.GHC <- ebb_crossprob_cor_R(d=d, bounds=GHC_z_bounds, correlations=pairwise_cors)
  }
  
  return ( p.GHC)
}
