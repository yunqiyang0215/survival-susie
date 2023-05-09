
# Function to compute the term exp(mx+v2x^2/2).
# @param x: n by 1 vector for a predictor
# @param m: posterior mean of a single predictor, b|gamma
# @param logv2: log of the posterior variance of b.
# @return a vector
compute_exp_quadratic <- function(x, m, logv2){
  a = exp(m*x + exp(logv2)*x^2/2)
  return(a)
}


# @param h.vec: p by 1 vector, the elbo values for subproblem, q(b|gamma_j = 1), j =1,\cdots, p.
# @param alpha: p by 1 vector, posterior of gamma.
# @param pi: prior distribution of gamma.
compute_objective <- function(h.vec, alpha, pi){
  p = length(alpha)
  elbo = sum(alpha*h.vec) - sum(alpha*log(alpha/pi), na.rm = TRUE)
  return(elbo)
}


negative_h <- function(post, x, y, ss, sigma2){
  m = post[1]
  logv2 = post[2]
  a <- compute_exp_quadratic(x, m, logv2)
  h <- sum(y*x*m - ss*a) + (logv2 - log(sigma2) - (m^2+exp(logv2))/sigma2)/2
  return(-h)
}

safenormalize <- function (vec) {
  a <- max(vec)
  vec <- exp(vec - a)
  return(vec/sum(vec))
}


poisson_ser <- function(X, y, ss, m.vec, logv2.vec, sigma2, pi, alpha, maxiter, tol){
  p = length(alpha)
  progress <- list(elbo = rep(NA, maxiter), sigma2 = rep(NA, maxiter),
                   h.vec = matrix(NA, nrow = maxiter, ncol = p),
                   m.vec = matrix(NA, nrow = maxiter, ncol = p),
                   logv2.vec = matrix(NA, nrow = maxiter, ncol = p),
                   alpha = matrix(NA, nrow = maxiter, ncol = p))

  # Store initial values
  h.vec = rep(NA, p)
  for (j in 1:p){
    post = c(m.vec[j], logv2.vec[j])
    h.neg = negative_h(post, X[, j], y, ss, sigma2)
    h.vec[j] = -h.neg
  }
  progress$elbo[1] = compute_objective(h.vec, alpha, pi)
  progress$sigma2[1] = sigma2
  progress$h.vec[1, ] =  h.vec
  progress$m.vec[1, ] = m.vec
  progress$logv2.vec[1, ] = logv2.vec
  progress$alpha[1, ] = alpha


  for (iter in 2:maxiter){
    for (j in 1:p){
      par.init = c(m.vec[j], logv2.vec[j])
      res = optim(par = par.init, fn = negative_h, x = X[,j], y = y, ss = ss, sigma2 = sigma2,
                  method = "BFGS")

      m.vec[j] = res$par[1]
      logv2.vec[j] = res$par[2]
      h.vec[j] =  -res$value
    }

    h.weighted = h.vec + log(pi)
    alpha = safenormalize(h.weighted)
    sigma2 = sum(alpha*(m.vec^2 + exp(logv2.vec)))

    progress$elbo[iter] = compute_objective(h.vec, alpha, pi)
    progress$h.vec[iter, ] =  h.vec
    progress$m.vec[iter, ] = m.vec
    progress$logv2.vec[iter, ] = logv2.vec
    progress$alpha[iter, ] = alpha
    progress$sigma2[iter] = sigma2
  }
  progress$v2.vec = exp(progress$logv2.vec)
  return(progress)
}


# Function to compute the log-scalar for each observation for a specific l.
# @return length-n vector contains log-scalar
compute_log_smat <- function(X, m.mat, logv2.mat, l){
  n = nrow(X)
  log_scalar = rep(NA, n)
  for (i in 1:n){
    # Store the unweighted expectation exp(m_jx_{ij} + v_j^2x_{ij}^2/2) for observation i.
    exp.mat = matrix(NA, nrow = L, ncol = p)
    for (j in 1:p){
      exp.mat[, j] = sapply(1:L, function(l) compute_exp_quadratic(X[i,j], m.mat[l, j], logv2.mat[l, j]))
    }
    # log of the weighted expectation, weighted by alpha.
    logexp.vec = log(rowSums(exp.mat * alpha.mat))
    log_scalar[i] = sum(logexp.vec[-l])
  }
  return(log_scalar)
}

poisson_susie <- function(X, y, L, m.mat, logv2.mat, alpha.mat, pi.mat, sigma2.vec, maxiter.out){
  n = nrow(X)
  log_scalar_mat = matrix(NA, ncol = L, nrow = n)
  for (i in 1: maxiter.out){
    for (l in 1:L){
      log_scalar_mat[,l] = compute_log_smat(X, m.mat, logv2.mat, l)
      result <- poisson_ser(X, y, exp(log_scalar_mat[,l]), m.mat[l, ], logv2.mat[l, ], sigma2 = sigma2.vec[l], pi.mat[l, ], alpha.mat[l, ], maxiter = 20, tol)

      m.mat[l, ] = tail(result$m.vec, n =1)
      logv2.mat[l, ] = tail(result$logv2.vec, n = 1)
      alpha.mat[l, ] = tail(result$alpha, n = 1)
      sigma2.vec[l] = tail(result$sigma2, n =1)
    }
  }
  return(list(m.mat = m.mat, v2.mat = exp(logv2.mat), alpha.mat = alpha.mat, sigma2.vec = sigma2.vec))
}


