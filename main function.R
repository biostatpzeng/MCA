MCA <- function(input.data, ref.data, weight.matrix = NULL, method = "BT",lambda = 0.99,
                t = 0.05,k = 0.5,a = 0.0383,p = c(0.1, 0.2, 0.3, 0.4, 1),
                path = NULL,o = c(1:8, Inf),f = 0.85,n = NULL)
{
  pkgs <- list("MASS","mvtnorm","Matrix", "CompQuadForm","Rcpp","corpcor",
               "mkatr","harmonicmeanp")
  checking<-unlist(lapply(pkgs, require, character.only = T))
  if(any(checking==F))
    stop("Please install the necessary packages first!")
  Z = as.numeric(paste(input.data$`t.value`))
  Z[Z == Inf] <- 60
  P = as.numeric(paste(input.data$P))
  if((min(P) < 0) & (max(P) > 1))
    stop("Please insure P value is between 0 and 1")
  M = length(Z)
  if (is.null(weight.matrix)) {W<-rep(1,M)} else {W<-weight.matrix}
  R = cor(as.matrix(ref.data), use = "p")
  R = lambda * R + diag(1-lambda, nrow=dim(R)[2], ncol=dim(R)[2])
  
  if (length(Z) <= 1) stop("less than 1 SNP.")
  if (!(length(Z) == dim(R)[1] & dim(R)[2] == dim(R)[1])) {
    stop("dimension do not match. Check dimension of Z and R.")
  }
  sumstat.GBM <- as.function(get(paste(method)))
  if(method == "TPM"){
    if((t < 0) & (t > 1)) stop("Please input a fixed value t between 0 and 1")
    return(sumstat.GBM(P,trunc = t))
  }else if(method %in% c("FCP","Simes")){
    return(sumstat.GBM(P))
  }else if(method %in% c("GATES","SimpleM")){
    return(sumstat.GBM(P,R))
  }else if(method %in% c("RTP","ART")){
    if((k < 0) & (k > 1)) stop("Please input a fixed value k between 0 and 1")
    k = min(k,1)*M
    return(sumstat.GBM(P,k = k))
  }else if(method == "ART.A"){
    if((k < 0) & (k > 1)) stop("Please input a fixed value k between 0 and 1")
    k = min(k,1)*M
    return(sumstat.GBM(P,W = W,k = k))
  }else if(method == "GM"){
    return(sumstat.GBM(P,a = a))
  }else if(method == "VEGAS"){
    return(sumstat.GBM(P,R,vegas.pct = p))
  }else if(method == "aSPUs"){
    return(sumstat.GBM(Z,R,pow = o))
  }else if(method == "PCA"){
    if(is.null(n))  stop("require the sample number n")
    if((f < 0) & (f > 1)) stop("Please input a fixed value fra between 0 and 1")
    return(sumstat.GBM(Z,R,W = W,fra = f,n=n))
  }else if(method %in% c("ACAT","HMP")){
    return(sumstat.GBM(P,W = W))
  }else if(method %in% c("BT","DOT","QT","SKAT","SKATO")){
    return(sumstat.GBM(Z,R,W = W))
  }else if(method %in% c("MLR","FLM")){
    sumstat.GBM <- as.function(get(paste(method)))
    if(is.null(n))  stop("require the sample number n")
    return(sumstat.GBM(Z,W = W,R,n = n))
  }else if(method %in% c("GBJ","BJ","HC","GHC","minP")){
    return(as.numeric(sumstat.GBM(Z,R)[2]))
  }else {
    stop("Please use the right method!")
  }
}

GBM_simple <- function(P,Z,R,W=NULL, method = "BT",
                       t = 0.05,k = 0.5,a = 0.0383,p = c(0.1, 0.2, 0.3, 0.4, 1),
                       path = NULL,o = c(1:8, Inf),f = 0.85,n = NULL){
  sumstat.GBM <- as.function(get(paste(method)))
  if(method == "TPM"){
    if((t < 0) & (t > 1)) stop("Please input a fixed value t between 0 and 1")
    return(sumstat.GBM(P,trunc = t))
  }else if(method %in% c("FCP","Simes")){
    return(sumstat.GBM(P))
  }else if(method %in% c("GATES","SimpleM")){
    return(sumstat.GBM(P,R))
  }else if(method %in% c("RTP","ART")){
    if((k < 0) & (k > 1)) stop("Please input a fixed value k between 0 and 1")
    k = min(k,1)*M
    return(sumstat.GBM(P,k = k))
  }else if(method == "ART.A"){
    if((k < 0) & (k > 1)) stop("Please input a fixed value k between 0 and 1")
    k = min(k,1)*M
    return(sumstat.GBM(P,W = W,k = k))
  }else if(method == "GM"){
    return(sumstat.GBM(P,a = a))
  }else if(method == "VEGAS"){
    return(sumstat.GBM(P,R,vegas.pct = p))
  }else if(method == "aSPUs"){
    return(sumstat.GBM(Z,R,pow = o))
  }else if(method == "PCA"){
    if(is.null(n))  stop("require the sample number n")
    if((f < 0) & (f > 1)) stop("Please input a fixed value fra between 0 and 1")
    return(sumstat.GBM(Z,R,W = W,fra = f,n=n))
  }else if(method %in% c("ACAT","HMP")){
    return(sumstat.GBM(P,W = W))
  }else if(method %in% c("BT","DOT","QT","SKAT","SKATO")){
    return(sumstat.GBM(Z,R,W = W))
  }else if(method %in% c("MLR","FLM")){
    sumstat.GBM <- as.function(get(paste(method)))
    if(is.null(n))  stop("require the sample number n")
    return(sumstat.GBM(Z,W = W,R,n = n))
  }else if(method %in% c("GBJ","BJ","HC","GHC","minP")){
    return(as.numeric(sumstat.GBM(Z,R)[2]))
  }else {
    stop("Please use the right method!")
  }
}


########################################################################

TPM <- function (P, trunc = 0.05) # add MCMC
{
  w <- prod(P^(P <= trunc))
  M <- length(P)
  if (w > trunc) {
    return(1)
  }
  else {
    pr <- 0
    for (k in 1:M) {
      s <- 0:(k - 1)
      term1 <- sum(w * (w <= (trunc^k)) * (((k * log(trunc)) - log(w))^s)/factorial(s))
      term2 <- (trunc^k) * (w > (trunc^k))
      pr <- pr + (choose(M, k) * ((1 - trunc)^(M - k)) * (term1 + term2))
    }
    return(pr)
  }
}

FCP <- function (P) 
{
  y = -2 * sum(log(P))
  pchisq(y, df = 2 * length(P), lower.tail = FALSE)
}

Simes <- function(P)
{
  P = P[order(P,decreasing = F)]
  M = length(P)
  P_Simes = NULL
  for(i in 1:M){
    P_Simes = c(P_Simes,(M*P[i]/i))
  }
  min(P_Simes)
}

GATES <- function (P, R) 
{
  P <- P[order(P,decreasing = F)]
  eff.snpcount.fun <- function(R) {
    snpcount.local <- dim(R)[1]
    ev <- eigen(R, only.values = TRUE)$values
    if (sum(ev < 0) != 0) {
      ev <- ev[ev > 0]
      ev <- ev/sum(ev) * snpcount.local
    }
    ev <- ev[ev > 1]
    snpcount.local - sum(ev - 1)
  }
  M <- length(P)
  eff.ld <- eff.snpcount.fun(R)
  candid <- sapply(1:M, function(i) {
    (eff.ld * P[i])/eff.snpcount.fun(R[1:i,1:i])
  })
  min(unlist(candid))
}

DOT_RTP <- function(Z,R, k=NULL,W=NULL) 
{
  M <- length(Z)
  if(is.null(k)) k <- round(0.5*M)
  y <- Z^2
  yy <- y %*% W
  ee <- eigen(R); eivec <- ee$vectors; eigva <- ee$values
  pc <- eivec %*% diag(sqrt(1/eigva)) %*% t(eivec)
  x <- (Z %*% pc)^2
  px <- 1-pchisq(x, df=1)
  w <- sum(-log(sort(px)[1:k]))
  integrate(function(x,y,m,n) 1-pgamma(log(qbeta(x,m+1,n-m))*m+y,m),
            lower=0,upper=1,w,k,M)$va
}

RTP <- function(P, k=NULL) 
{
  M <- length(P)
  if(is.null(k)) k <- round(0.5*M)
  w <- sum(-log(sort(P)[1:k]))
  integrate(function(x,y,m,n) 1-pgamma(log(qbeta(x,m+1,n-m))*m+y,m),
            lower=0,upper=1,w,k,M)$va
}

SimpleM <- function(P,R)
{
  M <- length(P)
  ev <- eigen(R, only.values = TRUE)$values
  ev_sort <- sort(ev, decreasing = TRUE)
  evs <- sum(ev_sort)
  M_eff_G <- 1
  for (k in 1:M) {
    temp <- sum(ev_sort[1:k])/evs
    if (temp >= 0.995) {
      M_eff_G <- k
      break
    }
  }
  1 - (1 - min(P))^M_eff_G
}

GM <- function(P,a=0.0383)
{
  M <- length(P)
  Y <- sum(qgamma(p=1-P,shape = a,scale = 1))
  1-pgamma(Y,shape = M*a,scale = 1)
}

DOT_ART <- function(Z,R,W=NULL,k=NULL) 
{
  M <- length(Z)
  y <- Z^2
  yy <- y %*% W
  ee <- eigen(R); eivec <- ee$vectors; eigva <- ee$values
  pc <- eivec %*% diag(sqrt(1/eigva)) %*% t(eivec)
  x <- (Z %*% pc)^2
  if(is.null(k)) k <- round(0.5*M)
  px <- 1-pchisq(x, df=1)
  P <- sort(px)
  Lw <- sum(log(P[1:(k-1)]))
  Pk <- P[k]
  d = (k-1)*(digamma(M+1) - digamma(k))
  ak = (k-1)*log(Pk) - Lw + qgamma(1-pbeta(Pk, k, M-k+1), shape=d)
  1 - pgamma(ak, shape=k+d-1)
}

ART <- function(P,k=NULL) 
{
  M <- length(P)
  P <- sort(P)
  Lw <- sum(log(P[1:(k-1)]))
  Pk <- P[k]
  d = (k-1)*(digamma(M+1) - digamma(k))
  ak = (k-1)*log(Pk) - Lw + qgamma(1-pbeta(Pk, k, M-k+1), shape=d)
  1 - pgamma(ak, shape=k+d-1)
}

DOT_ART.A <- function(Z,R,W=NULL,k=NULL)
{
  M <- length(Z)
  y <- Z^2
  yy <- y %*% W
  ee <- eigen(R); eivec <- ee$vectors; eigva <- ee$values
  pc <- eivec %*% diag(sqrt(1/eigva)) %*% t(eivec)
  x <- (Z %*% pc)^2
  if(is.null(k)) k <- round(0.5*M)
  px <- 1-pchisq(x, df=1)
  P <- sort(px)
  wgt <- rep(1,k)
  z <- P
  z[1] <- ( 1 - P[1] )^M
  for(j in 2:M) z[j] <- ((1-P[j]) / (1-P[j-1]))^((M-(j-1)))
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
  pmvnorm(lower = rep(-Inf,k), upper = rep(max(sZ), k), sigma = pCr)[1]
}

ART.A <- function(P,W=NULL,k=NULL)
{
  M <- length(P)
  if(is.null(k)) k <- round(0.5*M)
  P <- sort(P)
  if (is.null(W))   W <- rep(1, k)
  W <- W[1:k]
  z <- P
  z[1] <- ( 1 - P[1] )^M
  for(j in 2:M) z[j] <- ((1-P[j]) / (1-P[j-1]))^((M-(j-1)))
  p <- (1-z)[1:k]
  k = length(p)
  sumZ <- rep(0, k)
  y <- qnorm(p)
  z <- y
  gz <- z[1] * W[1]
  sumZ[1] <- gz
  for(i in 2:k) {
    gz <- p[ 1 : i ]
    for(j in 1:i) gz[j] <- z[j] * W[j]
    sumZ[i] <- sum(gz)
  }
  Lo = diag(k); Lo[lower.tri(Lo)] <- 1
  pSg <- Lo %*% diag(W[1:k]^2) %*% t(Lo)
  pCr <- cov2cor(pSg)
  sZ <- sumZ
  for(i in 1:k) {
    sZ[i] <- sumZ[i] / sqrt(diag(pSg))[i]
  }
  pmvnorm(lower = rep(-Inf,k), upper = rep(max(sZ), k), sigma = pCr)[1]
}

ACAT<-function(P,W=NULL)
{
  if (is.null(W)){
    W<-rep(1/length(P),length(P))
  }else{
    W<-W/sum(W)
  }
  
  if (sum(P<1e-16)==0){
    cct.stat<-sum(W*tan((0.5-P)*pi))
  }else{
    cct.stat<-sum((W[P<1e-16]/P[P<1e-16])/pi)
    cct.stat<-cct.stat+sum(W[!P<1e-16]*tan((0.5-P[!P<1e-16])*pi))
  }
  if (cct.stat>1e+15){
    P.ACAT<-(1/cct.stat)/pi
  }else{
    P.ACAT<-1-pcauchy(cct.stat)
  }
  return(P.ACAT)
}

HMP <- function(P,W=NULL){
  M <- length(P)
  W <- W/sum(W)
  P.HMP <- harmonicmeanp::p.hmp(p = P,w = W,L = M,
                                w.sum.tolerance = 1e-6,multilevel = TRUE)
  return(P.HMP)
}

# COMBAT <- function (P, R, vegas.pct = c(0.1, 0.2, 0.3, 0.4, 1), 
#                     pca_cut_perc = 0.995)
# {
#   set.seed(12345)
#   pval_gates <- GATES(P, R)
#   pval_vegas <- VEGAS(P, R)
#   pval_simpleM <- SimpleM(P, R)
#   gene_pvals <- c(GATES = pval_gates, pval_vegas, simpleM = pval_simpleM)
#   gene_pval_mat = apply(simul_pval_mat, 1, func1, R = R)
#   gene_pval_mat = t(gene_pval_mat)
#   method_cor <- cor(gene_pval_mat)
#   order_pvals <- order(gene_pvals)
#   sort_pvals <- gene_pvals[order_pvals]
#   method_cor <- method_cor[order_pvals, order_pvals]
#   P.COMBAT <- GATES(sort_pvals, method_cor)
#   return(P.COMBAT)
# }

COMBAT <- function (P, R, vegas.pct = c(0.1, 0.2, 0.3, 0.4, 1), pca_cut_perc = 0.995, 
                    nperm = 100, seed = 12345, ncores = 4) 
{
  pvalues <- as.numeric(P)
  n_snps <- length(pvalues)
  set.seed(seed)
  pval_gates <- GATES(P, R)
  pval_vegas <- VEGAS(P, R)
  pval_simpleM <- SimpleM(P, R)
  gene_pvals <- c(GATES = pval_gates, pval_vegas, simpleM = pval_simpleM)
  rd <- rmvnorm(nperm, mean = rep(0, n_snps), sigma = R)
  rd2 <- rd^2
  simul_pval_mat <- pchisq(rd2, 1, lower.tail = FALSE)
  func1 = function(x, R, vegas.pct, pca_cut_perc = 0.995) {
    p_gates <- GATES(P, R)
    p_vegas <- VEGAS(P, R)
    p_simpleM <- SimpleM(P, R)
    c(p_gates, p_vegas, p_simpleM)
  }
  if (ncores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    cl = parallel::makeCluster(ncores)
    parallel::clusterSetRNGStream(cl, .Random.seed)
    gene_pval_mat = parallel::parApply(cl, simul_pval_mat, 
                                       1, func1, R = R, vegas.pct = vegas.pct, pca_cut_perc = pca_cut_perc)
    parallel::stopCluster(cl)
  }
  else {
    gene_pval_mat = apply(simul_pval_mat, 1, func1, R = R, 
                          vegas.pct = vegas.pct, pca_cut_perc = pca_cut_perc)
  }
  gene_pval_mat = t(gene_pval_mat)
  method_cor <- cor(gene_pval_mat)
  order_pvals <- order(gene_pvals)
  sort_pvals <- gene_pvals[order_pvals]
  method_cor <- method_cor[order_pvals, order_pvals]
  p_combat_simes <- GATES(sort_pvals, method_cor)
  res <- c(COMBAT = p_combat_simes, gene_pvals)
  res
}

BT <- function (Z, R, W=NULL) 
{
  M = length(Z)
  if (is.null(W)) 
    W = rep(1, M)
  chi2 <- sum(W * Z)^2
  KKK <- sum(t(R * W) * W)
  chi2 <- chi2/KKK
  pchisq(chi2, 1, lower.tail = FALSE)
}

DOT <- function (Z, R,W=NULL) 
{
  ee <- eigen(R, TRUE)
  eivec <- ee$vectors
  eigva <- ee$values
  eigva <- sqrt(1/eigva)
  H <- eivec %*% (eigva * t(eivec))
  H <- 0.5 * (H + t(H))
  X <- H %*% Z
  ssq <- sum(X^2 * W)
  1 - pchisq(ssq, df=length(eigva))
}

QT <- function(Z, R,W=NULL, approx.method="chi-square")
{
  M <- length(Z)
  if (is.null(W)) {W = rep(1, M)}
  yy <- Z^2 %*% W
  rho.av2 <- sqrt(mean(R[lower.tri(R)]^2))
  if(approx.method=="chi-square"){
    P.QT <- 1-pchisq( (yy - (1 - rho.av2)*(M - 1) ) / (rho.av2*(M - 1) + 1), df=1)
  }else if(approx.method=="equicorrelation"){
    Sgm.e <- (1-rho.av2)*diag(1,M,M) + rho.av2 * matrix(1,M,M)
    P.QT <- CompQuadForm::imhof(yy, lambda = eigen(Sgm.e)$values, delta = rep(0, M), epsrel = 1e-11, limit = 2e5)$Qq
  }
  return(P.QT)
}

VEGAS <- function (P, R, vegas.pct = c(0.1, 0.2, 0.3, 0.4, 1), n_simul = 100000) 
{
  stopifnot(length(P) == ncol(R))
  vegas_vec <- ceiling(vegas.pct * ncol(R))
  vegas_vec <- sort(vegas_vec)
  if (vegas_vec[1] > 1) {
    vegas.pct <- c(0, vegas.pct)
    vegas_vec <- c(1, vegas_vec)
  }
  chisq_vec <- qchisq(P, 1, lower.tail = FALSE)
  chisq_vec[chisq_vec == Inf] <- 60
  n_snps <- length(P)
  n_tests <- length(vegas_vec)
  TS_obs <- rep(NA, n_tests)
  TS_obs[1] <- max(chisq_vec, na.rm = TRUE)
  chisq_vec <- sort(chisq_vec, decreasing = TRUE)
  for (j in 2:n_tests) TS_obs[j] <- sum(chisq_vec[1:vegas_vec[j]])
  rd <- rmvnorm(n_simul, mean = rep(0, n_snps), sigma = R)
  rd2 <- rd^2
  rd2 <- apply(rd2, 1, sort, decreasing = TRUE)
  pPerm0 <- rep(NA, n_tests)
  T0s <- apply(rd2, 2, max)
  pPerm0[1] <- (sum(T0s >= TS_obs[1]) + 1)/(length(T0s) + 1)
  for (j in 2:n_tests) {
    for (i in 1:n_simul) T0s[i] <- sum(rd2[1:vegas_vec[j], i])
    pPerm0[j] <- (sum(T0s >= TS_obs[j]) + 1)/(length(T0s) + 1)
  }
  v1 <- paste0("VEGAS.p", vegas.pct)
  v1[vegas_vec == ncol(R)] <- "VEGAS.all"
  v1[vegas_vec == 1] <- "VEGAS.max"
  names(pPerm0) <- v1
  pPerm0
}

aSPUs <- function(Z, R, pow = c(1:8, Inf), n.perm = 100000)
{
  M <- length(Z)
  eS <- eigen(R, symmetric = TRUE)
  ev <- eS$values
  CovSsqrt <- eS$vectors %*% diag(sqrt(pmax(ev, 0)), M)
  Ts = rep(NA, length(pow))
  for (j in 1:length(pow)) {
    if (pow[j] < Inf) 
      Ts[j] = sum(Z^pow[j])
    else Ts[j] = max(abs(Z))
  }
  pPerm0 = rep(NA, length(pow))
  T0s = numeric(n.perm)
  s <- sample(1:10^5, 1)
  for (j in 1:length(pow)) {
    set.seed(s)
    for (b in 1:n.perm) {
      U00 <- rnorm(M, 0, 1)
      U0 <- CovSsqrt %*% U00
      if (pow[j] < Inf) {
        T0s[b] = round(sum(U0^pow[j]), 8)
      }
      if (pow[j] == Inf) {
        T0s[b] = round(max(abs(U0)), 8)
      }
    }
    pPerm0[j] = sum(abs(Ts[j]) <= abs(T0s))/n.perm
    P0s = ((n.perm - rank(abs(T0s))) + 1)/(n.perm)
    if (j == 1) 
      minp0 = P0s
    else minp0[which(minp0 > P0s)] = P0s[which(minp0 > P0s)]
  }
  P.aSPUs <- (sum(minp0 <= min(pPerm0)) + 1)/(n.perm + 1)
  return(P.aSPUs)
}

MLR <- function(Z,R,W = NULL,n)
{
  M = length(Z)
  if (is.null(W))  W <- rep(1, M)
  Mat <- t(W) %*% R %*% W
  U05 = Mat %^% (-.5)
  U05Z = U05 %*% Z        # U^-.5 %*% Z
  R2 <- sum(U05Z^2) / n       # assimp M/n
  Fstat = ((n - M - 1) / M) * R2 / (1-R2)  # F-statistic
  P.MLR = as.double(pf(Fstat, M, n - M - 1, lower.tail = FALSE)) # F-test
  return(P.MLR)
}

FLM <- function(Z,R,W=NULL,n)
{
  M <- length(Z)
  if (is.null(W))  W <- rep(1, M)
  Mat <- t(W) %*% R %*% W
  BUB_1 <- Mat %^% (-1)
  
  BZstat <- as.vector(t(W) %*% Z)
  RSS <- n - sum(BZstat * (BUB_1 %*% BZstat))
  
  Fstat <- ((n - 25) / 25) * (n - RSS) / RSS   # F-statistic
  P.FLM <- pf(Fstat, 25, n - 25, lower.tail = FALSE)
  
  return(P.FLM)
}

PCA <- function(Z,R,W=NULL,fra = 0.85,n) 
{
  M <- length(Z)
  if (is.null(W))  W <- rep(1, M)
  WZ  <- as.vector(W * Z) * sqrt(n)
  WUW <- as.matrix(t(R * W) * W) * n
  ee <- eigen(WUW, TRUE)
  eivec <- ee$vectors
  eigva <- ee$values
  eigva[eigva < 0] <- 0
  prop.var <- eigva / sum(eigva)
  CPV <- cumsum(prop.var)
  M85 <- min(which(CPV >= fra))  # components for which Explained variance fraction is about 85%
  BBB <- as.matrix(eivec[,1:M85])
  GY <- as.vector(t(BBB) %*% WZ)
  CC <- as.matrix(t(BBB) %*% WUW %*% BBB)
  m <- qr(BBB)$rank
  if (m > 1) {
    RSS <- (n - sum(GY * as.vector((CC %^% (-1)) %*% GY)))
  } else { RSS <- (n - GY * GY / CC) }
  Fstat <- ((n - m - 1) / m) * (n - RSS) / RSS    # F-statistic
  p <- pf(Fstat, m, n - m - 1, lower.tail = FALSE)
  minP <- 100;
  minM <- 0
  if (p < minP) minM <- M85
  minP <- min(p, minP)
  c(minP)
}

SKAT <- function(Z, R, W = NULL) 
{
  M = length(Z)
  if (is.null(W))  W = rep(1, M)
  Zw = Z * W
  Rw = t(R * W) * W
  eR = eigen(Rw, sym = TRUE)
  lamR = abs(eR$val)
  Qv = sum(Zw^2)
  KATpval(Qv, lamR)
}

SKATO <- function(Z, R, W = NULL, rho = c((0:5/10)^2, 0.5, 1))
{
  M = length(Z)
  L = length(rho)
  if (L <= 2) return(NA)
  if (is.null(W))
    W = rep(1, M)
  Zw = Z * W
  Rw = t(R * W) * W
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
  for (k in 1:L1) qval[k] = Liu.qval.mod(minP, Lamk[[k]])
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
  P.SKATO = min(p.value, minP * L)
  return(P.SKATO)
  # list(p.value = c(A = p.value, S2 = pvalv, S = pvalb), 
  #      pval = pval, rho.est = rho[which.min(pval)])
}

if (T)
{
  prepare.GBJ <- function(path)
  {
    if(is.null(path)) stop("require the path of R package GBJ")
    f1 = list.files(path);f1 <- f1[-which(f1=="ebb.cpp")]
    for(i in 1:length(f1))  source(paste0(path,"/",f1[i]));rm(i)
    Rcpp::sourceCpp(paste0(path,"/ebb.cpp"))
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
      if ((tmp$ifault > 0) | (pval[i] <= 0) | (pval[i] >= 1)) 
        pval[i] = Sadd.pval(Q.all[i], lambda)
    }
    return(pval)
  }
  
  KAT.pval <- function (Q.all, lambda, acc = 1e-27, lim = 1e+06) 
  {
    pval = rep(0, length(Q.all))
    i1 = which(is.finite(Q.all))
    for (i in i1) {
      tmp = davies(Q.all[i], lambda, acc = acc, lim = lim)
      pval[i] = tmp$Qq
      if ((tmp$ifault > 0) | (pval[i] <= 0) | (pval[i] >= 1)) 
        pval[i] = Sadd.pval(Q.all[i], lambda)
    }
    return(pval)
  }
  
  Sadd.pval <- function (Q.all, lambda) 
  {
    sad = rep(1, length(Q.all))
    if (sum(Q.all > 0) > 0) {
      sad[Q.all > 0] = sapply(Q.all[Q.all > 0], saddle, lambda = lambda)
    }
    id = which(is.na(sad))
    if (length(id) > 0) {
      sad[id] = Liu.pval(Q.all[id], lambda)
    }
    return(sad)
  }
  
  saddle <- function (x, lambda) 
  {
    d = max(lambda)
    lambda = lambda/d
    x = x/d
    k0 = function(zeta) -sum(log(1 - 2 * zeta * lambda))/2
    kprime0 = function(zeta) sapply(zeta, function(zz) sum(lambda/(1 - 
                                                                     2 * zz * lambda)))
    kpprime0 = function(zeta) 2 * sum(lambda^2/(1 - 2 * zeta * 
                                                  lambda)^2)
    n = length(lambda)
    if (any(lambda < 0)) {
      lmin = max(1/(2 * lambda[lambda < 0])) * 0.99999
    }
    else if (x > sum(lambda)) {
      lmin = -0.01
    }
    else {
      lmin = -length(lambda)/(2 * x)
    }
    lmax = min(1/(2 * lambda[lambda > 0])) * 0.99999
    hatzeta = uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, 
                      upper = lmax, tol = 1e-08)$root
    w = sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
    v = hatzeta * sqrt(kpprime0(hatzeta))
    if (abs(hatzeta) < 1e-04) {
      return(NA)
    }
    else {
      return(pnorm(w + log(v/w)/w, lower.tail = FALSE))
    }
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
  
  Liu.qval.mod <- function (pval, lambda) 
  {
    c1 = rep(0, 4)
    c1[1] = sum(lambda)
    c1[2] = sum(lambda^2)
    c1[3] = sum(lambda^3)
    c1[4] = sum(lambda^4)
    muQ = c1[1]
    sigmaQ = sqrt(2 * c1[2])
    s1 = c1[3]/c1[2]^(3/2)
    s2 = c1[4]/c1[2]^2
    beta1 = sqrt(8) * s1
    beta2 = 12 * s2
    type1 = 0
    if (s1^2 > s2) {
      a = 1/(s1 - sqrt(s1^2 - s2))
      d = s1 * a^3 - a^2
      l = a^2 - 2 * d
    }
    else {
      type1 = 1
      l = 1/s2
      a = sqrt(l)
      d = 0
    }
    muX = l + d
    sigmaX = sqrt(2) * a
    df = l
    q.org = qchisq(pval, df = df, lower.tail = FALSE)
    (q.org - df)/sqrt(2 * df) * sigmaQ + muQ
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
  
  nsp <- function (R, eps = NULL, ...) 
  {
    if (is.null(eps)) 
      eps <- sqrt(.Machine$double.eps)
    CC <- eigen(R, TRUE)
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
    return((sum(T0s >= TS_obs[1]) + 1)/(length(T0s) + 1))
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
  
  ext_simes <- function (x, cor_r) 
  {
    eff.snpcount.fun <- function(R) {
      R <- as.matrix(R)
      snpcount.local <- dim(R)[1]
      if (snpcount.local <= 1) 
        return(1)
      ev <- eigen(R, only.values = TRUE)$values
      if (sum(ev < 0) != 0) {
        ev <- ev[ev > 0]
        ev <- ev/sum(ev) * snpcount.local
      }
      ev <- ev[ev > 1]
      snpcount.local - sum(ev - 1)
    }
    eff.snpcount.global <- eff.snpcount.fun(cor_r)
    n_values <- length(x)
    candid <- sapply(1:n_values, function(i) {
      (eff.snpcount.global * x[i])/eff.snpcount.fun(cor_r[1:i, 
                                                          1:i])
    })
    p_ext_simes <- min(candid)
    p_ext_simes
  }
}
