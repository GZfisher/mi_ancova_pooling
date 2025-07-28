# Fetch the character of the missing data
prelim.norm.new <- function (x) 
{
  #x <- all_mono
  if (is.vector(x)) 
    x <- matrix(x, length(x), 1)
  n <- nrow(x)
  p <- ncol(x)
  storage.mode(x) <- "double"
  r <- 1 * is.na(x) # add col of patient ID in x
  nmis <- as.integer(apply(r, 2, sum))
  names(nmis) <- dimnames(x)[[2]]
  mdp <- as.double((r %*% (2^((1:ncol(x)) - 1))) + 1)
  ro <- order(mdp)
  x <- matrix(x[ro, ], n, p)
  mdp <- mdp[ro]
  r <- matrix(r[ro, ], n, p)
  ro <- order(ro)
  mdpst <- as.integer(seq(along = mdp)[!duplicated(mdp)])
  mdp <- unique(mdp)
  npatt <- length(mdpst)
  r <- 1 - r
  r <- matrix(r[mdpst, ], npatt, p)
  if (npatt == 1) 
    tmp <- format(n)
  if (npatt > 1) 
    tmp <- format(c(mdpst[2:npatt], n + 1) - mdpst)
  dimnames(r) <- list(tmp, dimnames(x)[[2]])
  storage.mode(r) <- "integer"
  if (sum(is.na(x)) < length(x)) {
    mvcode <- as.double(max(x[!is.na(x)]) + 1000)
    x <- .na.to.snglcode(x, mvcode)
    tmp <- .Fortran("ctrsc", x, n, p, numeric(p), numeric(p), 
                    mvcode, PACKAGE = "norm")
    x <- tmp[[1]]
    xbar <- tmp[[4]]
    sdv <- tmp[[5]]
    x <- .code.to.na(x, mvcode)
  }
  if (sum(is.na(x)) == length(x)) {
    xbar <- rep(0, p)
    sdv <- rep(1, p)
  }
  d <- as.integer((2 + 3 * p + p^2)/2)
  psi <- .Fortran("mkpsi", p, matrix(as.integer(0), p + 1, 
                                     p + 1), PACKAGE = "norm")[[2]]
  if (npatt > 1) 
    nmdp <- as.integer(c(mdpst[-1], n + 1) - mdpst)
  if (npatt == 1) 
    nmdp <- n
  sj <- .Fortran("sjn", p, npatt, r, integer(p), PACKAGE = "norm")[[4]]
  nmon <- .Fortran("nmons", p, npatt, r, nmdp, sj, integer(p), 
                   PACKAGE = "norm")[[6]]
  last <- .Fortran("lasts", p, npatt, sj, integer(npatt), 
                   PACKAGE = "norm")[[4]]
  tmp <- .Fortran("layers", p, sj, integer(p), integer(1), 
                  PACKAGE = "norm")
  layer <- tmp[[3]]
  nlayer <- tmp[[4]] # number of missing variables
  list(x = x, n = n, p = p, r = r, nmis = nmis, ro = ro, mdpst = mdpst, 
       nmdp = nmdp, npatt = npatt, xbar = xbar, sdv = sdv, 
       d = d, psi = psi, sj = sj, nmon = nmon, last = last, 
       layer = layer, nlayer = nlayer)
}

# MCMC imputation function
step1 <- function(data, nimpute, emmaxits, maxits, seed) {
  s <- prelim.norm.new(data)     # Prepare data for norm package
  thetahat <- em.norm(s, maxits = emmaxits)   # Get initial parameters using EM algorithm
  rngseed(seed)                  # Set random seed
  
  all_mono_new <- data.frame(data)
  all_mono_new[,"impno"] <- 0   # Add imputation indicator for original data
  
  for (i in 1:nimpute) {
    theta <- mda.norm(s, thetahat, steps = maxits, showits = TRUE)   # Run MCMC to update parameters
    all_mono_one_time <- data.frame(imp.norm(s, theta, data))   # Impute missing values
    # Ensure monotone structure by setting post-missing values to NA for each row
    for (j in 1:nrow(data)) {
      last_num <- max(which(!is.na(data[j,])))
      if (last_num < ncol(data)) {
        all_mono_one_time[j, (last_num+1):ncol(data)] <- "is.na<-"(all_mono_one_time[j, (last_num+1):ncol(data)])
      }
    }
    all_mono_one_time$impno <- i   # Add imputation number
    all_mono_new <- rbind(all_mono_new, all_mono_one_time)  # Append result
    print(paste0("MCMC imputation: ", i, "..."))
  }
  return(all_mono_new)
}

# Monotone regression imputation using mice
step2 <- function(data, nimpute, method, formula_list, seed) {
  
  reg <- list()
  for (i in 1:nimpute) {
    seed = seed + 1
    # Filter data for current imputation
    reg[[paste0("x", i)]] <- data %>% filter(impno == i)
    # Sequentially apply regression formulas
    for (j in 1:length(formula_list)) {
      x <- mice::mice(reg[[paste0("x", i)]], m = 1, method = method, 
                      formulas = formula_list[j], seed = seed)
      reg[[paste0("x", i)]] <- mice::complete(x)
    }
  }
  
  complete <- do.call(rbind, reg)   # Combine imputed datasets
  return(complete)
}