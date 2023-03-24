


# Function to compute the term exp(mx+v2x^2/2).
# @param x: a vector of predictors
# @param v2: posterior variance of b.
# @return a vector
compute_exp_quadratic <- function(x, m, v2){
  a = exp(m*x + v2*x^2/2)
  return(a)
}


# @param x: a vector of predictors
# @param y: a vector of outcomes
# @param d: a vector of status indicator. d == 1, outcome observed;
# d == 0, censored.
# @param h0: baseline hazard
# @param m: posterior mean of b.
# @param v2: posterior variance of b.
# @param s2: prior variance of b.
compute_elbo <- function(x, y, d, h0, m, v2, s2){
  a <- compute_exp_quadratic(x, m, v2)
  elbo = sum(d*(log(h0)+m*x)- h0*y*a)-(log(s2)-log(v2)+(m^2+v2)/s2)/2
  return(elbo)
}

# Function to calculate partial derivative w.r.t. m.
partial_dv_m <- function(m, x, y, d, s2, v2, h0){
  a <- compute_exp_quadratic(x, m, v2)
  res = sum(d*x - h0*x*y*a)-m/s2
  return(res)
}

# Function to calculate partial derivative w.r.t. v2.
partial_dv_v2 <- function(v2, x, y, s2, m, h0){
  a <- compute_exp_quadratic(x, m, v2)
  res = sum(h0*x^2*y*a)+1/s2-1/v2
  res = -res/2
  return(res)
}

# Function to minimize negative elbo w.r.t. v2 and m
# @param post: a vector contains (m, v2)
minF <- function(post, x, y, d, h0, s2){
  elbo = compute_elbo(x, y, d, h0, post[1], post[2], s2)
  return(-elbo)
}

# Function to minimize negative elbo w.r.t. m.
minF_m <- function(m, x, y, d, h0, v2, s2){
  elbo = compute_elbo(x, y, d, h0, m, v2, s2)
  return(-elbo)
}


# Note: the tol param is a place holder.
update_q <- function(x, y, d, h0, m, v2, s2, maxiter, tol, lower, upper){
  par.init = c(m, v2)
  progress <- matrix(NA, nrow = maxiter, ncol = 5)
  colnames(progress) = c("elbo", "h0", "m", "v2", "s2")

  for (iter in 1:maxiter){
    par = optim(par = par.init, fn = minF, x = x, y = y, d = d, h0 = h0, s2 = s2,
              method = "L-BFGS-B", lower = lower, upper = upper)$par
    m = par[1]
    v2 = par[2]
    s2 <- m^2+ v2
    h0 <- sum(d)/sum(y*exp(m*x+v2*x^2/2))
    elbo <- compute_elbo(x, y, d, h0, m, v2, s2)
    progress[iter, ] =  c(elbo, h0, m, v2, s2)
  }
  return(progress)
}



# v2 <- uniroot(function(v2) partial_dv_v2(v2, x, y, s2, m, h0),
#              c(1e-5,1e6))$root

#m <- uniroot(function(m) partial_dv_m(m, x, y, d, s2, h0),
#             c(-1e6,1e6))$root
