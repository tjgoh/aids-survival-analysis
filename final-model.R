library(msm)
library(expm)
library(ucminf)
library(ggplot2)
digits <- 3

# Using msm:
cat("\nUsing msm...\n")
# Define Q:
q <- 0.1
Q <- matrix( c(NA,q,0,0,q,
               q,NA,q,0,q,
               0,q,NA,q,q,
               0,0,q,NA,q,
               0,0,0,0,NA),5,5,byrow=TRUE)
# Fit model:
final_model <-msm(state~obstime, subject=id, data=dta_cov,
                           center=FALSE, 
                           qmatrix= Q, deathexact = TRUE, 
                           censor.states = 1:4,   censor = -2,
                           covariates = ~  obstime + drug + prevOI,
                           constraint = list(obstime=c(1,2,3,4,2,3,7,2,3,2), drug = c(1,666,70,50,666,71,50,666,72,666),
                                             prevOI = c(1,666,70,50,666,71,50,666,72,666)), #putting contraints on parameters
                           fixedpars = c(18,20,21, 24, 26, 27), 
                           control=list(fnscale=5000,maxit=500) 
) 
minus2LL   <- final_model$minus2loglik
p          <- final_model$estimates
p.se       <- sqrt(diag(final_model$covmat))
cat("-2loglik =", final_model$minus2loglik,"\n")
cat("AIC = ", final_model$minus2loglik + 2*length(final_model$opt$par))
cat("Parameter estimates:\n")
print(cbind(p=round(p,digits), 
            se=round(p.se,digits)),quote=FALSE) 

# Using manually written code:

# Prepare data for quick access:
dta_cov.split <- split(dta_cov, dta_cov$id)

# Prepare raw transition data:
o1 <- o2 <- dt <- dr <- ti <- pr <- rep(0, nrow(dta_cov))
n <- 0

for (i in 1:N) {
  # Extracting data for subject i:
  dta_cov.i <- dta_cov.split[[i]]
  O <- dta_cov.i$state 
  t <- dta_cov.i$obstime 
  drug <- dta_cov.i$drug[1]
  prev <- dta_cov.i$prevOI[1]
  # Loop over individual follow-up:
  for (j in 2:length(O)) { 
    n <- n + 1
    dt[n] <- t[j] - t[j - 1] # Time between checkups
    o1[n] <- O[j - 1] # State at start of interval
    o2[n] <- O[j] # State at end of interval
    dr[n] <- drug
    ti[n] <- t[j-1] # Time at start of interval
    pr[n] <- prev
  }
}    

# Function A for P matrix:
Pmatrix <- function(Q, t) { 
  decom <- eigen(Q * t)
  exp.D <- diag(exp(decom$values))
  U <- cbind(decom$vectors)
  U %*% exp.D %*% solve(U) 
}
# Defining function for P matrix:
Pmatrix <- function(Q, t) {
  expm(Q*t) #going from Q to P = exp[(t1-t2)Q]
}

# Loglikelihood:
p <- final_model$opt$par
loglikelihood <- function(p) {
  # Parameters:
  beta <- p
  if (mon) cat("beta = ", round(beta,digits),"\n")
  
  # Loop over intervals:
  loglik <- 0
  for (i in 1:n) {
    # Q matrix:
    Q <- matrix(0, D, D)
    # Off diagonal: 
    Q[1,2] <- exp(beta[1] + beta[11]*ti[i] + beta[16]*dr[i] + beta[19]*pr[i]) 
    Q[1,D] <- exp(beta[2] + beta[12]*ti[i] + beta[17]*dr[i] + beta[20]*pr[i]) 
    Q[2,1] <- exp(beta[3] + beta[13]*ti[i]) 
    Q[2,3] <- exp(beta[4] + beta[14]*ti[i] + beta[18]*dr[i] + beta[21]*pr[i])  
    Q[2,D] <- exp(beta[5] + beta[12]*ti[i] + beta[17]*dr[i] + beta[20]*pr[i]) 
    Q[3,2] <- exp(beta[6] + beta[13]*ti[i]) 
    Q[3,4] <- exp(beta[7] + beta[15]*ti[i] + beta[18]*dr[i] + beta[21]*pr[i]) 
    Q[3,D] <- exp(beta[8] + beta[12]*ti[i] + beta[17]*dr[i] + beta[20]*pr[i]) 
    Q[4,3] <- exp(beta[9] + beta[13]*ti[i]) 
    Q[4,D] <- exp(beta[10] + beta[12]*ti[i] + beta[17]*dr[i] + beta[20]*pr[i]) 
    
    # Diagonal:
    for (r in 1:(D - 1)){
      Q[r, r] <- -sum(Q[r, ]) #-sum of the diagonals
    }
    # Pmatrix:
    P <- Pmatrix(Q, dt[i]) #transition probability matrix - defined for time interval dt
    # Likelihood contribution:
    if(o2[i] == D){  #if C = died
      contribution <- P[o1[i], 1:(D - 1)] %*% Q[1:(D - 1), D] #multplied by q
    } 
    if(o2[i] == censored){ #B
      contribution <- P[o1[i], 1:(D - 1)]%*%rep(1,(D-1)) #
    }
    if(o2[i] != censored & o2[i]!=D){ #survived the study A
      contribution <- P[o1[i], o2[i]]
    }
    # Update likelihood:
    loglik <- loglik + log(contribution)
  }
  # Return:
  if (mon) cat("-2*Loglik = ", -2*loglik,"\n")
  -loglik
}

# Maximise with ucminf():
cat("\nUsing ucminf...\n")
mon <- FALSE #Monitoring
start_time <- Sys.time()
p.start <- rep(-3, 21)
control <- list(reltol = 1e-04, maxit = 500) 
max_final_model <- ucminf(par = p.start, fn = loglikelihood, hessian = 1) 
end_time <- Sys.time()
time.run <- end_time - start_time
print(time.run)
cat("Starting values =", p.start, "\n")
cat("Convergence code =", max_final_model$convergence, "\n")
cat("-2loglik =", 2 * max_final_model$value, "\n")
cat("AIC =", 2 * max_final_model$value + 2 * length(max_final_model$par), "\n")
p <- max_final_model$par
fisher <- solve(max_final_model$hessian)
p.se <- sqrt(diag(fisher))
print(cbind(p = round(p, digits), se = round(p.se, digits)), quote = FALSE)

