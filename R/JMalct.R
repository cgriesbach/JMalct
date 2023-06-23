#' Main function to carry out variable allocation for a joint model!
#'
#' @param y vector containing the longitudinal outcome
#' @param time vector containing the observed survival times
#' @param delta censoring indicator of the survival times
#' @param X design matrix of possibly high dimensional baseline covariates which underly the allocation process
#' @param t vector of longitudinal observation times
#' @param id cluster variable
#' @param mstop number of total boosting iterations, default is mstop = 500
#' @param nu learning rate of the single boosting updates, default is nu = 0.1

JMalct = function(y, time, delta, X, t, id, mstop = 500, nu = .1){
  ### lambda im exponenten!

  ### scale to unit variance
  # sdt = sd(t)
  # sdX = apply(X, 2, "sd")
  # X = scale(X, center = FALSE)
  # t = scale(t, center = FALSE)

  ### useful functions
  g.s = function(gamma, s) matrix(gamma, ncol = 2, byrow = T)[,s]

  step.risk = function(nu.l = 0, nu.s = 0, fit.l = 0, fit.s = 0, x = alpha){
    if(length(fit.l) == 1) fit.l = rep(0, N)
    etal = etal + nu.l*fit.l; etal.id = etal.id + nu.l*fit.l[first]; etas = etas + nu.s*fit.s

    long = -sum(dnorm(y, mean = etal, sd = sqrt(sigma2), log = TRUE))

    integral = exp(etas + x*etal.id) * ((exp(x*(betat + g.s(gamma, 2))*time) - 1)/(x*(betat + g.s(gamma, 2))))
    i.na = which(is.na(integral))
    integral[i.na] = exp(etas[i.na])*time[i.na]
    surv = -sum(delta*(etas + x*etal.id) - integral)

    return(long + surv)
  }

  alpha.risk = function(alpha) step.risk(x = alpha)

  ### basic definitions
  N = length(id)
  n = length(unique(id))
  id.t = as.numeric(table(id))
  first = rep(FALSE, N)
  for(i in 1:N) first[which.max(id==i)] = TRUE

  ### set up fixed effects
  Xl = X
  Xs = X[first,]
  po = ncol(X)

  ### set up random effects
  if(po > n){
    cormod = glmboost(y ~ cbind(t, X), control = boost_control(nu = 0.1, mstop = 100))
    # cv10f = cv(model.weights(cormod), type = "kfold")
    # cvm = cvrisk(cormod, folds = cv10f, papply = lapply)
    # cormod[mstop(cvm)]
    cc = sort(unique(cormod$xselect()))[-1] - 1
  }else{
    cc = 1:po
  }

  p1 = rep(seq(1, 2*n, 2), 2) + rep(0:1, each = n)
  p2 = rep(seq(1, 2*n, n), n) + rep(0:(n-1), each = 2)
  P1 = sparseMatrix(seq_along(p1), p1)
  P2 = sparseMatrix(seq_along(p2), p2)

  Xcor = cbind(1, Xs[,cc])
  Xcor = Xcor %*% chol2inv(chol(crossprod(Xcor))) %*% t(Xcor)
  Xcor = kronecker(diag(2), Xcor)
  Xcor = as.matrix(P2%*%(diag(2*n) - Xcor)%*%P1)

  Z = list()
  for(i in 1:n){
    Z[[i]] = cbind(1, t[id==i])
  }
  Z = as.matrix(bdiag(Z))
  lambda.df = mboost:::df2lambda(Z, 100, weights = 1)[2]
  S.ran = Xcor %*% chol2inv(chol(crossprod(Z) + lambda.df*diag(2*n))) %*% t(Z)

  ### initialize starting values
  int = mean(y)
  gamma = rep(0, 2*n)
  alpha = betat = 0
  betal = betas = rep(0, 2*po)
  lambda = log(sum(delta)/sum(time))
  sigma2 = var(y)

  offset = lme(y ~ 1, random = ~ t | id, control = lmeControl(opt = "optim", singular.ok = TRUE, returnObject = TRUE))
  int = offset$coefficients$fixed[1]
  gamma = Xcor %*% as.numeric(t(cbind(offset$coefficients$random$id)))
  sigma2 = offset$sigma^2

  INT = BETAT = BETAL = GAMMA = SIGMA2 = LAMBDA = ALPHA = BETAS = c()

  set.l = set.s = NULL

  ### initialize probing
  Xs = X[first,]
  Xs = cbind(Xs, apply(Xs, 2, sample))
  Xl = Xs[id,]
  # Xl = cbind(Xl, apply(Xl, 2, sample))
  Xs = Xl[first,]
  p = ncol(Xl)

  ### prepare matrices for component-wise fitting
  M.l = M.s = list()
  for(r in 1:p){
    x = cbind(1, t, Xl[,r])
    M.l[[r]] = solve(crossprod(x)) %*% t(x)
    x = cbind(1, Xs[,r])
    M.s[[r]] = solve(crossprod(x)) %*% t(x)
  }

  sb.l = sb.s = TRUE

  for(m in 1:mstop){
    if(!(sb.l | sb.s)) break

    etal = int + Xl%*%betal + t*betat + Z%*%gamma
    etal.id = int + Xl[first,]%*%betal
    etas = lambda + Xs%*%betas

    ###############################################################
    #### SA1 ######################################################
    ###############################################################
    if(length(set.l) < po & sb.l){
      u = y - etal - Xl%*%betas

      fits.l = matrix(0, 4, p)
      for(r in 1:p){
        fits.l[1:3, r] = M.l[[r]]%*%u
        fits.l[4, r] = sum((u - cbind(1, t, Xl[,r]) %*% fits.l[1:3, r])^2)
      }
      ranking = order(fits.l[4,])
      best.l = ranking[! ranking %in% set.l][1]

      opt.l = optimize(step.risk, interval = c(0, 100), fit.l = fits.l[3,best.l]*Xl[,best.l], nu.s = 0, fit.s = 0, x = alpha)

      df.l = ifelse(betal[best.l] == 0, 1, 0)

      if(best.l > po) sb.l = FALSE
    }else{
      opt.l = list()
      opt.l$objective = Inf
    }


    ###############################################################
    #### SB1 ######################################################
    ###############################################################
    if(length(set.s) < po & sb.s){
      u = delta - exp(etas + alpha*etal.id) * (exp(alpha*((betat + g.s(gamma, 2))*time)) - 1) / (alpha*(betat + g.s(gamma, 2)))
      u.na = which(is.na(u))
      u[u.na] = delta[u.na] - exp(etas[u.na])*time[u.na]

      fits.s = matrix(0, 3, p)
      for(r in 1:p){
        fits.s[1:2, r] = M.s[[r]]%*%u
        fits.s[3, r] = sum((u - cbind(1, Xs[,r]) %*% fits.s[1:2, r])^2)
      }
      ranking = order(fits.s[3,])
      best.s = ranking[! ranking %in% set.s][1]

      opt.s = optimize(step.risk, interval = c(0, 100), fit.s = fits.s[2,best.s]*Xs[,best.s], nu.l = 0, fit.l = 0, x = alpha)

      df.s = ifelse(betas[best.s] == 0, 1, 0)

      if(best.s > po) sb.s = FALSE
    }else{
      opt.s = list()
      opt.s$objective = Inf
    }

    long.best = opt.l$objective + df.l < opt.s$objective + df.s
    surv.best = opt.l$objective + df.l > opt.s$objective + df.s


    ###############################################################
    #### SA2 ######################################################
    ###############################################################
    if(long.best & sb.l){
      int = int + nu*opt.l$minimum*fits.l[1, best.l]
      betat = betat + nu*opt.l$minimum*fits.l[2, best.l]
      betal[best.l] = betal[best.l] + nu*opt.l$minimum*fits.l[3, best.l]

      set.s = unique(c(set.s, best.l))

      etal = int + Xl%*%betal + t*betat + Z%*%gamma
      gamma = gamma + nu*S.ran %*% (y - etal)

      etal = int + Xl%*%betal + t*betat + Z%*%gamma
      sigma2 = var(y - etal)

      if(m > 1){
        krit = sqrt(sum( (BETAL[m-1,] - betal)^2 )) / sqrt(sum( betal^2 ))
        if(krit < 1e-5) sb.l = FALSE
      }
    }


    ###############################################################
    #### SB2 ######################################################
    ###############################################################
    if(surv.best & sb.s){
      lambda = lambda + opt.s$minimum*fits.s[1,best.s]
      betas[best.s] = betas[best.s] + nu*opt.s$minimum*fits.s[2,best.s]

      set.l = unique(c(set.l, best.s))

      if(m > 1){
        krit = sqrt(sum( (BETAS[m-1,] - betas)^2 )) / sqrt(sum( betas^2 ))
        if(krit < 1e-5) sb.s = FALSE
      }
    }

    optim.int = c(alpha - nu, alpha + nu)
    alpha = optimize(alpha.risk, interval = optim.int)$minimum

    INT = c(INT, int)
    BETAT = c(BETAT, betat)
    BETAL = rbind(BETAL, betal)
    GAMMA = cbind(GAMMA, gamma)
    SIGMA2 = c(SIGMA2, sigma2)
    LAMBDA = c(LAMBDA, lambda)
    ALPHA = c(ALPHA, alpha)
    BETAS = rbind(BETAS, betas)

    if(m%%50 == 0) print(paste("Iteration", m, sb.l, sb.s))
  }

  structure(list(int = int, betat = betat, betal = betal[1:po], gamma = gamma,
                 betas = betas[1:po], alpha = alpha, lambda = lambda, sigma2 = sigma2,
                 INT = INT, BETAT = BETAT, BETAL = BETAL[,1:po], GAMMA = t(GAMMA),
                 BETAS = BETAS[,1:po], ALPHA = ALPHA, LAMBDA = LAMBDA, SIGMA2 = SIGMA2,
                 gamma0 = g.s(gamma, 1), gammat = g.s(gamma, 2)))
}
