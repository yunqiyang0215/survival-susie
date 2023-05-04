
# Function to compute the term exp(mx+v2x^2/2).
# @param x: n by 1 vector for a predictor
# @param m: posterior mean of a single predictor, b|gamma
# @param v2: posterior variance of b.
# @return a vector
compute_exp_quadratic <- function(x, m, v2){
  a = exp(m*x + v2*x^2/2)
  return(a)
}


# @param h.vec: p by 1 vector, the elbo values for subproblem, q(b|gamma_j = 1), j =1,\cdots, p.
# @param alpha: p by 1 vector, posterior of gamma.
# @param pi: prior distribution of gamma.
compute_objective <- function(h.vec, alpha, pi){
  p = length(alpha)
  elbo = sum(alpha*h.vec) - sum(alpha*log(alpha/pi))
  return(elbo)
}


negative_h <- function(post, x, y, ss, sigma2){
  m = post[1]
  v2 = post[2]
  a <- compute_exp_quadratic(x, m, v2)
  h <- sum(y*x*m - ss*a) + (log(v2) - log(sigma2) - (m^2+v2)/sigma2)/2
  return(-h)
}

softmax <- function (vec) {
  a <- max(vec)
  vec <- exp(vec - a)
  return(vec/sum(vec))
}


update_q <- function(X, y, ss, m.vec, v2.vec, sigma2, pi, alpha, maxiter, tol, lower, upper){
  p = length(alpha)
  progress <- list(elbo = rep(NA, maxiter), sigma2 = rep(NA, maxiter),
                   h.vec = matrix(NA, nrow = maxiter, ncol = p),
                   m.vec = matrix(NA, nrow = maxiter, ncol = p),
                   v2.vec = matrix(NA, nrow = maxiter, ncol = p),
                   alpha = matrix(NA, nrow = maxiter, ncol = p))

  # Store initial values
  h.vec = rep(NA, p)
  for (j in 1:p){
    post = c(m.vec[j], v2.vec[j])
    h.neg = negative_h(post, X[, j], y, ss, sigma2)
    h.vec[j] = -h.neg
  }
  progress$elbo[1] = compute_objective(h.vec, alpha, pi)
  progress$sigma2[1] = sigma2
  progress$h.vec[1, ] =  h.vec
  progress$m.vec[1, ] = m.vec
  progress$v2.vec[1, ] = v2.vec
  progress$alpha[1, ] = alpha


  for (iter in 2:maxiter){
    for (j in 1:p){
      par.init = c(m.vec[j], v2.vec[j])
      res = optim(par = par.init, fn = negative_h, x = X[,j], y = y, ss = ss, sigma2 = sigma2,
                  method = "L-BFGS-B", lower = lower, upper = upper)

      m.vec[j] = res$par[1]
      v2.vec[j] = res$par[2]
      h.vec[j] =  -res$value
    }

    h.weighted = h.vec + log(pi)
    alpha = softmax(h.weighted)
    sigma2 = sum(alpha*(m.vec^2 + v2.vec))

    progress$elbo[iter] = compute_objective(h.vec, alpha, pi)
    progress$h.vec[iter, ] =  h.vec
    progress$m.vec[iter, ] = m.vec
    progress$v2.vec[iter, ] = v2.vec
    progress$alpha[iter, ] = alpha
  }
  return(progress)
}
