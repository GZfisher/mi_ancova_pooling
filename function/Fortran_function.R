mda.norm <- function (s, theta, steps = 1, showits = FALSE) 
{
  s$x <- .na.to.snglcode(s$x, as.double(999))
  tobs <- .Fortran("tobsmn", s$p, s$psi, s$n, s$x, s$npatt, 
                   s$r, s$mdpst, s$nmdp, s$last, integer(s$p), s$sj, s$layer, 
                   s$nlayer, s$d, matrix(0, s$nlayer, s$d), PACKAGE = "norm")[[15]]
  if (showits) 
    cat(paste("Steps of Monotone Data Augmentation:", "\n"))
  for (i in 1:steps) {
    if (showits) 
      cat(paste(format(i), "...", sep = ""))
    s$x <- .Fortran("is2n", s$d, theta, s$p, s$psi, s$n, 
                    s$x, s$npatt, s$r, s$mdpst, s$nmdp, s$sj, s$last, 
                    integer(s$p), integer(s$p), double(s$p), theta, PACKAGE = "norm")[[6]]
    theta <- .Fortran("ps2n", s$p, s$psi, s$n, s$x, s$npatt, 
                      s$r, s$mdpst, s$nmdp, integer(s$p), integer(s$p), 
                      s$nmon, s$sj, s$nlayer, s$d, tobs, numeric(s$d), 
                      numeric(s$d), numeric(s$p + 1), numeric(s$d), PACKAGE = "norm")[[19]]
  }
  if (showits) 
    cat("\n")
  theta
}

# Finds the column numbers of the observed variables, 
# and stores them in the first noc elements of oc. Does not go beyond column=last.
gtoc <- function(p, npatt, r, patt, oc, last) {
  noc <- 0
  for (j in 1:last) {
    if (r[patt, j] == 1) {
      noc <- noc + 1
      oc[noc] <- j
    }
  }
  return(list(noc = noc, oc = oc))
}

# Finds the column numbers of the missing variables, and stores them
# in the first nmc elements of mc. Does not go beyond column=last.
gtmc <- function(p, npatt, r, patt, mc, last) {
  nmc <- 0
  for (j in 1:last) {
    if (r[patt, j] == 0) {
      nmc <- nmc + 1
      mc[nmc] <- j
    }
  }
  return(list(mc=mc, nmc=nmc))
}
# Translated version of the 'tobsmn' subroutine
# Calculates known parts of the sufficient statistics for monotone data augmentation
tobsmn <- function(p, psi, n, x, npatt, r, mdpst, nmdp, last, oc, sj, layer, nlayer, d) {
  tobs <- matrix(0, nrow = nlayer, ncol = d)
  storage.mode(tobs) <- "double"
  lpatt <- 0
  for (l in 1:nlayer) {
    # l=1
    # oc=integer(p)
    fpatt <- lpatt + 1
    j <- p + 1
    repeat {
      j <- j - 1
      if (layer[j] != l) next
      lpatt <- sj[j]
      for (patt in fpatt:lpatt) {
        res <- gtoc(p, npatt, r, patt, oc, last[patt])
        noc <- res$noc
        oc <- res$oc
        for (i in (mdpst[patt]:(mdpst[patt] + nmdp[patt] - as.double(1)))) {
          tobs[l, psi[1, 1]] <- tobs[l, psi[1, 1]] + 1
          for (j in 1:noc) {
            tobs[l, psi[1, oc[j]+1]] <- tobs[l, psi[1, oc[j]+1]] + x[i, oc[j]]
            for (k in j:noc) {
              tobs[l, psi[oc[j]+1, oc[k]+1]] <- tobs[l, psi[oc[j]+1, oc[k]+1]] + x[i, oc[j]] * x[i, oc[k]]
            }
          }
        }
      }
      break
    }
  }
  return(tobs)
}

# Performs sweep on a symmetric matrix in packed storage.
# Sweeps on pivot position. Sweeps only the (0:submat,0:submat)
# submatrix.
# If dir=1, performs ordinary sweep. If dir=-1, performs reverse sweep.
# need to check paper 
swp <- function(d, theta, pivot, p, psi, submat, dir) {
  a <- as.double(theta[psi[pivot+1, pivot+1]])
  theta[psi[pivot+1, pivot+1]] <- -1 / a
  for (j in 0:submat) {
    if (j != pivot) theta[psi[j+1, pivot+1]] <- theta[psi[j+1, pivot+1]] / a * dir
  }
  for (i in 0:submat) {
    for (j in i:submat) {
      if ((i != pivot) && (j != pivot)) {
        b <- as.double(theta[psi[i+1, pivot+1]])
        c <- as.double(theta[psi[j+1, pivot+1]])
        theta[psi[i+1, j+1]] <- theta[psi[i+1, j+1]] - a * b * c
      }
    }
  }
  return(theta)
}


# Sweeps theta to condition on the observed variables
swpobs <- function(d, theta, p, psi, npatt, r, patt) {
  # d=s$d
  # theta=thetahat
  # p=s$p
  # psi=s$psi
  # npatt=s$naptt
  # r=s$r
  # patt=10
  storage.mode(theta) <- "double"
  for (j in 1:p) {
    # print(c(j,patt))
    if (!is.na(theta[psi[j+1, j+1]])) {
      if ((r[patt, j] == 1) && (theta[psi[j+1, j+1]] > 0)) {
        theta <- swp(d, theta, j, p, psi, p, 1)
      } else if ((r[patt, j] == 0) && (theta[psi[j+1, j+1]] < 0)) {
        theta <- swp(d, theta, j, p, psi, p, -1)
      }
    }
  }
  return(theta)
}



# Define a function that extracts submatrix of theta corresponding to the columns of mc
sigex <- function(d, theta, extr, p, psi, mc, nmc) {
  # Loop over the columns of mc
  if (nmc>=1) {
    for (j in 1:nmc) {
      for (k in j:nmc) {
        # Assign the value of theta to extr based on the psi matrix
        extr[psi[mc[j]+1, mc[k]+1]] <- theta[psi[mc[j]+1, mc[k]+1]]
      }
    }
  }
  # Return the extr vector
  return(extr)
}

# Define a function that performs Cholesky decomposition on theta matrix
chols <- function(d, theta, p, psi, mc, nmc) {
  # Loop over the columns of mc
  if (nmc>=1) {
    for (i in 1:nmc) {
      tmp <- as.double(0)
      if (i>=2) {
        for (k in 1:(i-1)) {
          # Accumulate the squared values of theta
          tmp <- tmp + theta[psi[mc[k]+1, mc[i]+1]]^2
        }
        # Compute the diagonal element of theta
        theta[psi[mc[i]+1, mc[i]+1]] <- sqrt(theta[psi[mc[i]+1, mc[i]+1]] - tmp)
        if (nmc>=i+1) {
          for (j in (i+1):nmc) {
            tmp <- as.double(0)
            for (k in 1:(i-1)) {
              # Accumulate the product of theta values
              tmp <- tmp + theta[psi[mc[k]+1, mc[i]+1]] * theta[psi[mc[k]+1, mc[j]+1]]
            }
            # Compute the off-diagonal element of theta
            theta[psi[mc[i]+1, mc[j]+1]] <- (theta[psi[mc[i]+1, mc[j]+1]] - tmp) / theta[psi[mc[i]+1, mc[i]+1]]
          }
        }
      }
    }
  }
  # Return the theta matrix
  return(theta)
}

# Define a function that generates a random number from a standard normal distribution
gauss <- function() {
  # Save the values of alt and next in a persistent environment

  if (!exists("alt", envir = env)) {
    env$alt <- 0
  }
  if (!exists("next", envir = env)) {
    env$"next" <- 0
  }
  # Define the value of pi
  pi <- 3.141593
  # Check the value of alt
  if ((env$alt != 0) & (env$alt != 1)) {
    env$alt <- 0
  }
  if (env$alt == 0) {
    # Generate two random numbers from a uniform distribution
    u1 <- runif(1)
    u2 <- runif(1)
    # Use the Box-Muller method to generate a standard normal random number
    gauss <- sqrt(-2 * log(u1)) * cos(2 * pi * u2)
    # Save the next standard normal random number for later use
    env$"next" <- sqrt(-2 * log(u1)) * sin(2 * pi * u2)
    # Set the value of alt to 1
    env$alt <- 1
  } else {
    # Use the previously saved standard normal random number
    gauss <- env$"next"
    # Set the value of alt to 0
    env$alt <- 0
  }
  # Return the standard normal random number
  return(gauss)
}

# Define a function that performs I-step of monotone data augmentation
is2n <- function(d, theta, p, psi, n, x, npatt, r, mdpst, nmdp, sj, last, oc, mc, z, cc, norm_func) {
  # Generate a random number from a standard normal distribution
  junk <- norm_func
  # Loop over the patterns
  tt <- theta
  for (patt in 1:npatt) {
    # Swap the observations based on the pattern
    tt <- swpobs(d, tt, p, psi, npatt, r, patt)
    # Get the missing columns and the number of missing columns
    res_mc <- gtmc(p, npatt, r, patt, mc, last[patt])
    mc <- res_mc$mc
    nmc <- res_mc$nmc
    # Get the observed columns and the number of observed columns
    res_oc <- gtoc(p, npatt, r, patt, oc, last[patt])
    oc <- res_oc$oc
    noc <- res_oc$noc
    # Extract the submatrix of theta corresponding to the columns of mc
    cc <- sigex(d, tt, cc, p, psi, mc, nmc)
    # Perform Cholesky decomposition on the submatrix of theta
    cc <- chols(d,cc,p,psi,mc,nmc)
    # Loop over the rows of the pattern
    for (i in (mdpst[patt]:(mdpst[patt] + nmdp[patt] - as.double(1)))) {
      # Loop over the columns of mc
      if (nmc >=1){
        for (j in 1:nmc) {
          # Compute the mean of x[i, mc[j]] based on theta and the observed values
          x[i, mc[j]] <- tt[psi[1, mc[j]+1]]
          for (k in 1:noc) {
            x[i, mc[j]] <- x[i, mc[j]] + tt[psi[oc[k]+1, mc[j]+1]] * x[i, oc[k]]
          }
          # Generate a random number from a standard normal distribution
          z[mc[j]] <- norm_func
          # Add the random noise to x[i, mc[j]] based on the submatrix of theta
          for (k in 1:j) {
            x[i, mc[j]] <- x[i, mc[j]] + z[mc[k]] * cc[psi[mc[j]+1, mc[k]+1]]
          }
        }
      }

    }
  }
  # Return the x matrix
  # return(list(cc=cc,x=x,tt=tt))
  return(x)
}

chol2 <- function(d, theta, p, psi, last) {
  for (i in 0:last) {
    tmp <- 0
    if (i>=1) {
      for (k in 0:(i-1)) {
        tmp <- tmp + theta[psi[k+1, i+1]]^2
      }
      theta[psi[i+1, i+1]] <- sqrt(theta[psi[i+1, i+1]] - tmp)
    }
    
    if (last>=(i+1)) {
      for (j in (i+1):last) {
        tmp <- 0
        if (i>=1) {
          for (k in 0:(i-1)) {
            tmp <- tmp + theta[psi[k+1, i+1]] * theta[psi[k+1, j+1]]
          }
          theta[psi[i+1, j+1]] <- (theta[psi[i+1, j+1]] - tmp) / theta[psi[i+1, i+1]]
        }
        
      }
    }
    
  }
  return(theta)
}

ph2thn <- function(d, theta, p, psi) {
  if (p>=2) {
    for (j in 1:(p-1)) {
      theta <- swp(d, theta, j, p, psi, j, 1)
    }
    for (j in 1:(p-1)) {
      theta <- swp(d, theta, j, p, psi, p, -1)
    }
  }
  
  return(theta)
}

ps2n <- function(p, psi, n, x, npatt, r, mdpst, nmdp, oc, mc, nmon, sj, nlayer, d, tobs, t, c, v, theta) {
  # P-step of monotone data augmentation
  t <- rep(0, d)
  lsj <- 0
  l <- 0
  if (p>=1) {
    for (j in p:1) {
      if (sj[j] > lsj) {
        l <- l + 1
        t[psi[1,1]] <- t[psi[1,1]] + tobs[l, psi[1,1]]
        for (k in 1:j) {
          t[psi[1,k+1]] <- t[psi[1,k+1]] + tobs[l, psi[1,k+1]]
          for (m in k:j) {
            t[psi[m+1,k+1]] <- t[psi[m+1,k+1]] + tobs[l, psi[m+1,k+1]]
          }
        }
        for (patt in (lsj+1):sj[j]) {
          res_mc <- gtmc(p, npatt, r, patt, mc, j)
          mc <- res_mc$mc
          nmc <- res_mc$nmc
          res_oc <- gtoc(p, npatt, r, patt, oc, j)
          oc <- res_oc$oc
          noc <- res_oc$noc
          for (i in mdpst[patt]:(mdpst[patt] + nmdp[patt] - 1)) {
            if (nmc>=1) {
              for (k in 1:nmc) {
                t[psi[1,mc[k]+1]] <- t[psi[1,mc[k]+1]] + x[i, mc[k]]
                for (m in 1:noc) {
                  t[psi[mc[k]+1,oc[m]+1]] <- t[psi[mc[k]+1,oc[m]+1]] + x[i,mc[k]]*x[i,oc[m]]
                }
                for (m in 1:k) {
                  t[psi[mc[k]+1,mc[m]+1]] <- t[psi[mc[k]+1,mc[m]+1]] + x[i,mc[k]]*x[i,mc[m]]
                }
              }
            }
          }
        }
      }
      if (sj[j] > lsj) {
        for (k in 0:(j-1)) {
          t <- swp(d, t, k, p, psi, j, 1)
        }
      }
      df <- nmon[j] + 3*(p-j) - 1
      theta[psi[j+1,j+1]] <- t[psi[j+1,j+1]] / rchisq(1, df)
      for (k in 0:(j-1)) {
        for (m in k:(j-1)) {
          c[psi[k+1,m+1]] <- -theta[psi[j+1,j+1]] * t[psi[k+1,m+1]]
        }
      }
      c <- chol2(d, c, p, psi, j-1)
      for (k in 0:(j-1)) {
        v[k+1] <- rnorm(1)
        theta[psi[k+1,j+1]] <- t[psi[k+1,j+1]]
        for (m in 0:k) {
          theta[psi[k+1,j+1]] <- theta[psi[k+1,j+1]] + c[psi[m+1,k+1]] * v[m+1]
        }
      }
      if (j > 1) {
        if (sj[j-1] > sj[j]) {
          for (k in 0:(j-1)) {
            t <- swp(d, t, k, p, psi, j-1, -1)
          }
        } else if (sj[j-1] == sj[j]) {
          t <- swp(d, t, j-1, p, psi, j-1, -1)
        }
      }
      lsj <- sj[j]
    }
  }
  
  theta[psi[1,1]] <- -1
  theta <- ph2thn(d, theta, p, psi)
  return(theta)
}

# complete I-step and P-step
mda_r <- function (s, theta, steps = 1, showits = FALSE, criterion = 1e-04)  {
  s$x <- .na.to.snglcode(s$x, as.double(999))
  tobs <- tobsmn(s$p, s$psi, s$n, s$x, s$npatt, s$r, s$mdpst,  
                 s$nmdp, s$last, integer(s$p), s$sj, s$layer, s$nlayer, s$d)
  if (showits)
    cat(paste("Steps of Monotone Data Augmentation:", "\n"))
  it <- 0
  converged <- FALSE
  while ((!converged) & (it < steps)) {
    old <- theta
    # I-step: impute missing data given current parameters
    s$x <- is2n(s$d, theta, s$p, s$psi, s$n,
                s$x, s$npatt, s$r, s$mdpst, s$nmdp, s$sj, s$last,
                integer(s$p), integer(s$p), double(s$p), theta, rnorm(n=1))
    # P-step: update parameters given completed data
    theta <- ps2n(s$p, s$psi, s$n, s$x, s$npatt,
                  s$r, s$mdpst, s$nmdp, integer(s$p), integer(s$p),
                  s$nmon, s$sj, s$nlayer, s$d, tobs, numeric(s$d),
                  numeric(s$d), numeric(s$p + 1), numeric(s$d))
    it <- it + 1
    if (showits)
      cat(paste(format(it), "...", sep = ""))
    converged <- max(abs(old - theta)) <= criterion
  }
  if (showits)
    cat("\n")
  theta
  # for (i in 1:steps) {
  #   if (showits)
  #     cat(paste(format(i), "...", sep = ""))
  #   # I-step: impute missing data given current parameters
  #   s$x <- is2n(s$d, theta, s$p, s$psi, s$n,
  #               s$x, s$npatt, s$r, s$mdpst, s$nmdp, s$sj, s$last,
  #               integer(s$p), integer(s$p), double(s$p), theta, rnorm(n=1))
  #   # P-step: update parameters given completed data
  #   theta <- ps2n(s$p, s$psi, s$n, s$x, s$npatt,
  #                 s$r, s$mdpst, s$nmdp, integer(s$p), integer(s$p),
  #                 s$nmon, s$sj, s$nlayer, s$d, tobs, numeric(s$d),
  #                 numeric(s$d), numeric(s$p + 1), numeric(s$d))
  # }
  # if (showits)
  #   cat("\n")
  # theta
}

# s$x <- .na.to.snglcode(s$x, as.double(999))
# tobs1 <- .Fortran("tobsmn", s$p, s$psi, s$n, s$x, s$npatt, 
#                  s$r, s$mdpst, s$nmdp, s$last, integer(s$p), s$sj, s$layer, 
#                  s$nlayer, s$d, matrix(0, s$nlayer, s$d), PACKAGE = "norm")[[15]]
# 
# tobs2 <- tobsmn(s$p, s$psi, s$n, s$x, s$npatt, s$r, s$mdpst, 
#                 s$nmdp, s$last, integer(s$p), s$sj, s$layer, s$nlayer, s$d)
# 
# x_is2n <- .Fortran("is2n", s$d, thetahat, s$p, s$psi, s$n, 
#          s$x, s$npatt, s$r, s$mdpst, s$nmdp, s$sj, s$last, 
#          integer(s$p), integer(s$p), double(s$p), thetahat, PACKAGE = "norm")
# theta <- .Fortran("ps2n", s$p, s$psi, s$n, x_is2n[[6]], s$npatt, 
#                   s$r, s$mdpst, s$nmdp, integer(s$p), integer(s$p), 
#                   s$nmon, s$sj, s$nlayer, s$d, tobs1, numeric(s$d), 
#                   numeric(s$d), numeric(s$p + 1), numeric(s$d), PACKAGE = "norm")[[19]]
# env <- environment()
# x_is2n_r_gauss <- is2n(s$d, thetahat, s$p, s$psi, s$n, 
#                    s$x, s$npatt, s$r, s$mdpst, s$nmdp, s$sj, s$last, 
#                    integer(s$p), integer(s$p), double(s$p), thetahat, gauss())
# x_is2n_r_rnorm <- is2n(s$d, thetahat, s$p, s$psi, s$n, 
#                        s$x, s$npatt, s$r, s$mdpst, s$nmdp, s$sj, s$last, 
#                        integer(s$p), integer(s$p), double(s$p), thetahat, rnorm(n=1))
# theta_r <- ps2n(s$p, s$psi, s$n, x_is2n_r_rnorm$x, s$npatt, 
#                 s$r, s$mdpst, s$nmdp, integer(s$p), integer(s$p), 
#                 s$nmon, s$sj, s$nlayer, s$d, tobs1, numeric(s$d), 
#                 numeric(s$d), numeric(s$p + 1), numeric(s$d))
# sum(x_is2n[[16]]-x_is2n_r$cc)
# 
# for (patt in 1:71) {
#   print(patt)
#   res_mc <- gtmc(s$p, s$npatt, s$r, patt, integer(s$p), last[patt])
#   mc <- res_mc$mc
#   nmc <- res_mc$nmc
#   # th_n <- .Fortran("chols",s$d, thetahat, s$p, s$psi, mc, nmc, PACKAGE = "norm")[[2]]
#   th_nr <- swpobs(s$d, thetahat, s$p, s$psi, s$npatt, s$r, patt)
#   
#   c <- sigex(s$d,th_nr,thetahat,s$p,s$psi,mc,nmc)
# }
# 
# th_nr <- chols(s$d, th_nr, s$p, s$psi, mc, nmc)
