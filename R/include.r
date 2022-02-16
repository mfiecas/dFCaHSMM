hsmm_viterbi <- function(res,Z,m){
  n <- dim(res$u.all)[1]
  nsub <- dim(res$u.all)[3]
  iv <- matrix(NA, nsub, n)
  for (Sub in 1:nsub){
    probs <- drop(res$u.all[,,Sub])
    xi <- matrix(NA,n,dim(res$u.all)[2])
    foo <- res$delta.hat[Sub,]*probs[1,]
    xi[1,] <- foo/sum(foo)
    lambda = get_lambdan(Z[Sub,],res$beta.hat)
    B = create_B_cpp(lambda,m,res$A.hat[,,Sub],tau=1,verbose=TRUE)

    for (i in 2:n){
      foo <- apply(xi[i-1,]*B, 2, max)*probs[i,]
      xi[i,]<- foo/sum(foo)
    }
    iv[Sub,n] <- which.max(xi[n,])
    for (i in (n-1):1){
      iv[Sub, i]<- which.max(B[,iv[Sub,i+1]]*xi[i,])
    }
  }
  # Convert states to state aggregates
  num.states = length(m)
  states = rep(NA,sum(m))
  states[1:m[1]] = 1
  for(s in 2:num.states){
    states[(sum(m[1:(s-1)])+1):sum(m[1:s])] = s
  }
  iv.hsmm = matrix(NA,nsub,n)
  for(Sub in 1:nsub){
    for(time in 1:n){
      iv.hsmm[Sub,time] = states[iv[Sub,time]]
    }
  }
  
  return(iv.hsmm)
}

pseudo_resid_hsmm = function (dat, res,m) {
  num_state = length(m)
  NSUB = dim(dat)[3]
  P = dim(dat)[2]
  TIME = dim(dat)[1]
  pSig = array(NA,dim=c(P,P,sum(m)))
  pmu = array(NA,dim=c(sum(m),P))
  resu1 <- matrix(NA, TIME, NSUB)
  for (sub in 1:NSUB){
    for (t in 1:TIME) {
      counter = 1
      for (i in 1:num_state){
        for(j in 1:m[i]){
          pSig[, , counter] = res$Sigma.hat[, , i] * res$u.all[t, counter, sub]  
          pmu[counter, ] = res$mu.hat[i, ] * res$u.all[t, counter, sub]
          counter = counter+1
        }
      }
      wmu = colSums(pmu)
      #varw = rowSums(pSig, dims = 2) - wmu %*% t(wmu)
      varw = rowSums(pSig, dims = 2)
      resu1[t, sub] = mahalanobis(dat[t, , sub], center = wmu, cov = varw)
    }
  }
  return(resu1)
}

get.metastate = function(state.sequence,metastates){
  # Convert states to metastates
  num.meta = length(metastates)
  metavec = rep(NA,num.meta)
  for(m in 1:num.meta){
    metavec[metastates[[m]]] = m
  }
  N = nrow(state.sequence)
  meta.hsmm = matrix(NA,nrow=N,ncol=ncol(state.sequence))
  for(n in 1:N){
    for(time in 1:ncol(state.sequence)){
      meta.hsmm[n,time] = metavec[state.sequence[n,time]]
    }
  }
  return(meta.hsmm)
}

hsmm_attempts = function(Xt,
                         Z,
                         params.init=NULL,
                         beta.init=NULL,
                         m,
                         ATTEMPTS=5,
                         tol=1e-4,
                         maxiter=250,
                         pars = NULL,
                         verbose=FALSE,
                         chol = FALSE,
                         refit=FALSE){
  if(is.null(pars)){
    pars = list(mu = 1,sd=.5)
  }

  result.attempts = list()
  P = dim(Xt)[2]
  num.states = nrow(params.init$A)
  for(attempt in 1:ATTEMPTS){
    if(verbose){
      print(sprintf('Working on HSMM, attempt = ...%g',attempt))
    }
    mu.init = params.init$mu
    Sigma.init = params.init$Sigma
    A.init = params.init$A

    result.attempts[[attempt]] = tryCatch({
      if(is.null(beta.init)){
        beta0 = matrix(rnorm(num.states*ncol(Z),
                                 mean=pars$mu,
                                 sd=pars$sd),
                       nrow=num.states,
                       ncol=ncol(Z)) 
      }
      else{
        beta0 = beta.init
      }

      hsmm_cpp(Xt,Z,
               A.init,
               mu.init,Sigma.init,
               beta_init = beta0,
               m = m,
               maxiter = maxiter,
               tol = tol,
               verbose = verbose,
               refit = refit)},
      error = function(e){
        print("Evaluation error. Restarting with new initial values.")
        })
  }
  return(result.attempts)
}

hsmm_permutation = function(results, Xt,Z,m,num.perms=100,ATTEMPTS=5,tol=1e-6,maxiter=250,penalty='none',chol=FALSE,verbose=FALSE){
  pars = list(mu = results$mu.hat,
              Sigma = results$Sigma.hat,
              A = results$A.hat)
  # beta.init = results$beta.hat
  # beta.init = NULL
  result.hsmm.perm = list()
  beta.perm = array(NA,dim=c(dim(results$beta.hat),num.perms))
  for(perm in 1:num.perms){
    print(sprintf("Permutation %g of %g",perm,num.perms))
    result.hsmm.perm = fit.ahsmm(Xt = Xt,
                                 Z = Z[sample(1:nrow(Z),nrow(Z),replace=FALSE),],
                                 m = m,
                                 mu.init = results$mu.hat,
                                 Sigma.init = results$Sigma.hat,
                                 A.init = results$A.hat,
                                 beta.init = NULL,
                                 ATTEMPTS = ATTEMPTS,
                                 tol=tol,
                                 maxiter = maxiter,
                                 pars = list(mu = 1, sd = 0.5),
                                 verbose=verbose,
                                 refit = TRUE)
    beta.perm[,,perm] = result.hsmm.perm$beta.hat
  }
  return(beta.perm)
}

em2_highdim_covar_perm = function(Xt,Z,PERM=500,result.hsmm,m,savefn,penalty='none',maxiter=250,tol=1e-5,verbose=TRUE){
  beta.hat.perm = array(NA,dim = c(dim(result.hsmm$beta.hat),PERM))
  N = nrow(Z)
  num.states = nrow(result.hsmm$A.hat)
  for(perm in 1:PERM){
    print(sprintf('Working on permutation %g of %g',perm,PERM))
    params.hmm.init = list(mu.expand = result.hsmm$mu.hat,
                           Sigma.expand = result.hsmm$Sigma.hat,
                           A.hat = result.hsmm$A.hat)
    hsmm.attempts.perm = em_hsmm_attempts(Xt,Z[sample(1:N),],params.hmm.init,m,ATTEMPTS=10,tol=tol,maxiter=maxiter,verbose=verbose,refit=TRUE)
    beta.hat.perm[,,perm] = hsmm.attempts.perm[[hsmm.attempts.perm %>%
                                                  lapply(extract.likelihood) %>%
                                                  unlist() %>%
                                                  which.max()]]$beta.hat
    save(beta.hat.perm, file = savefn)
  }
  return(beta.hat.perm)
}

swap = function(vec){
  return(c(vec[2],vec[1]))
}

norm.f = function(vec,square = TRUE){
  # Frobenius norm
  if(square){
    return((sum(vec^2)))
  }
  else{
    return(sqrt(sum(vec^2)))  
  }
}

logsumexp = function(vec){
  maxvec = max(vec)
  return(maxvec + log(sum(exp(vec-maxvec))))
}

covmat = function(Xt,mu,u){
  P = ncol(Xt)
  TIME = nrow(Xt)
  covest = matrix(0,nrow=P,ncol=P)
  for(time in 1:TIME){
    covest = covest + u[time]*t(t(Xt[time,]-mu))%*%((Xt[time,]-mu))
  }
  return(covest/sum(u))
}

standardize_list = function(a,time=NULL){
  num.states = length(a)
  if(is.null(time)){
    for(t in 1:length(a[[1]])){
      if(!is.na(a[[1]][t])){
        sum.a = 0
        for(j in 1:num.states){
          sum.a = sum.a + a[[j]][t]
        }
        for(j in 1:num.states){
          a[[j]][t] = a[[j]][t]/sum.a
        }
      }
    }
    return(a)
  }
  else{
    if(!is.na(a[[1]][time])){
      sum.a = 0
      for(j in 1:num.states){
        sum.a = sum.a + a[[j]][time]
      }
      for(j in 1:num.states){
        a[[j]][time] = a[[j]][time]/sum.a
      }
    }
    return(list(a=a,scale=sum.a))
  }
}

eye = function(P){
  return(diag(rep(1,P)))
}

pk = function(r,lambda,tau=0){
  # Shifted poisson
  # lambda - rate parameter
  # tau - shift parameter
  if(length(r)==1){
    return(dpois(r-tau,lambda))
  }
  else{
    vec = rep(NA,length(r))
    for(j in 1:length(vec)){
      vec[j] = dpois(r[j]-tau,lambda)
    }
    return(vec)
  }
}

Fk = function(r,lambda,tau=0){
  # Shifted poisson cdf
  # Shifted poisson
  # lambda - rate parameter
  # tau - shift parameter
  ppois(r-tau,lambda)
}

ck = function(r,lambda,tau=0){
  if(length(r)>1){
    result = rep(NA,length(r))
    for(i in 1:length(r)){
      result[i] = ck(r[i],lambda,tau)
    }
    return(result)
  }
  tryCatch({
    if(abs(Fk(r,lambda,tau) - 1) < 1e-10){
      return(1)
    }
    if(Fk(r,lambda,tau) < 1){
      if(pk(r,lambda,tau) < 1e-10){
        return(0)
      }
      return(pk(r,lambda,tau)/(1-Fk(r-1,lambda,tau)))
    }
  },
  error = function(cond){
    print("Error!")
    print(sprintf("r = %g, lambda = %g, tau = %g, Fk(r,lambda,tau) = %g",r,lambda,tau,Fk(r,lambda,tau)))
    return(NA)
  })
}

create_Bij = function(lambda,tau=0,aij,mi,mj){
  Bij = matrix(0,nrow=mi,ncol=mj)
  Bijvec = rep(0,mi)
  for(j in 1:mi){
    Bijvec[j] = aij*ck(j,lambda,tau)
  }
  Bij[,1] = Bijvec
  return(Bij)
}

create_Bii = function(lambda,tau=0,mi){
  diagB = rep(NA,mi-1)
  for(j in 1:(mi-1)){
    diagB[j] = 1-ck(j,lambda,tau)
  }
  Bii = rbind(cbind(rep(0,mi-1),diag(diagB)),
              c(0,rep(0,mi-2),1-ck(mi,lambda,tau))) 
  return(Bii)
}

create_B = function(lambda,tau=NULL,m,A){
  # lambda - vector of rate parameters for the Poisson distribution
  # tau - vector of location parameters for the shift of the Poisson distribution
  # m - vector of the degree of approximation. Higher numbers correspond to better approximation, but more computation
  # A - transition probability matrix
  
  if(is.null(tau)){
    tau = rep(0,num.states)
  }
  num.states = length(m)
  B = matrix(NA,nrow=sum(m),ncol=sum(m))
  
  i = 1
  Btemp = create_Bii(lambda[i],tau[i],m[i])
  for(j in 2:num.states){
    Btemp = cbind(Btemp,create_Bij(lambda[i],tau[i],A[i,j],m[i],m[j]))
  }
  B[1:m[1],] = Btemp
  
  for(i in 2:num.states){
    # "Lower triangle"
    Btemp = NULL
    for(j in 1:(i-1)){
      Btemp = cbind(Btemp,create_Bij(lambda[i],tau[i],A[i,j],m[i],m[j]))
    }
    Btemp = cbind(Btemp,create_Bii(lambda[i],tau[i],m[i]))
    # "Upper triangle
    if(i < num.states){
      for(j in (i+1):num.states){
        Btemp = cbind(Btemp,create_Bij(lambda[i],tau[i],A[i,j],m[i],m[j]))
      }
    }
    B[(sum(m[1:(i-1)])+1):(sum(m[1:i])),] = Btemp
  }
  return(B)
}

params_transform = function(params_vec,num.states,m){
  # Takes the vectorized parameter space and transform them to a list containing the vectors and matrices
  llambda.id = 1:num.states
  ltau.id = (num.states+1):(2*num.states)
  Avec = params_vec[-(1:(2*num.states))]
  if(length(Avec) != num.states^2){
    print("Error!")
    break;
  }
  temp.tau = round(exp(params_vec[ltau.id]))
  for(j in 1:num.states){
    if(temp.tau[j] > m[j]){
      temp.tau[j] = m[j]
    }
    temp.tau[j] = ifelse(temp.tau[j] < 0,0,temp.tau[j])
  }
  return(list(lambda = exp(params_vec[llambda.id]), # Make sure they are nonnegative
              tau = temp.tau, # Make sure they are nonnegative integers
              A = matrix(Avec,nrow=num.states,ncol=num.states)))
}

term2 = function(params_vec,num.states,m,sum_t_vjk){
  # Make sure to use the number of *aggregate* states (num.states, NOT num.Bstates)
  params = params_transform(params_vec,num.states,m)
  B = create_B(params$lambda,params$tau,m,params$A)
  B.0id = which(B==0)
  B0 = B
  B0[B.0id] = 1e-6
  # B.nonzero = which(B!=0)
  # return(-sum(log(B[B.nonzero])*sum_t_vjk[B.nonzero]))
  return(-sum(log(B0)*sum_t_vjk))
}

term2.covar = function(params_vec,num.states,m,sum_t_vjk,Z){
  # Make this more flexible in case we start to add more parameters.
  # For now, assume the model log(\lambda) = \beta0 + beta1*Z
  # params_vec = c(beta0(1),beta0(2),... 
  #               beta1(1),beta1(2),...
  #               tau(1),tau(2),.....
  #               c(A))
  # where theta(j) is the parameter for the j-th state
  
  # For now, Z is only a vector of covariates
  N = length(Z)
  # First, get the indices of the parameters in params_vec
  beta0.id = 1:num.states
  beta1.id = (num.states+1):(2*num.states)
  ltau.id = (2*num.states+1):(3*num.states)
  Avec = params_vec[-(1:(3*num.states))]
  if(length(Avec) != num.states^2){
    print("Error in params_vec!")
    break;
  }
  
  # Transform the parameters, then create the (expanded) transition probability matrix
  temp.tau = round(exp(params_vec[ltau.id]))
  # Make sure the shift parameter is not too large
  for(j in 1:num.states){
    if(temp.tau[j] > m[j]){
      temp.tau[j] = m[j] 
    }
    temp.tau[j] = ifelse(temp.tau[j] < 0, 0, temp.tau[j])
  }
  tau = temp.tau
  A = matrix(Avec,nrow=num.states,ncol=num.states)
  # Add up each subject's contributions
  # Make sure to transform the lambda parameter
  total = 0
  for(n in 1:N){
    lambdan = exp(params_vec[beta0.id] + params_vec[beta1.id]*Z[n])
    B = create_B(lambdan,tau,m,A)
    B.0id = which(B==0)
    B0 = B
    B0[B.0id] = 1e-6 # For numerical stability
    total = total + sum(log(B0)*sum_t_vjk[,,n])
  }
  return(-total) # Negative sign since nlm() minimizes, not maximizes
}

term2.nocovar = function(params_vec,num.states,m,sum_t_vjk){
  # Same as term2.covar, but no covariates
  # For now, assume the model log(\lambda) = \beta0
  # params_vec = c(beta0(1),beta0(2),... 
  #               c(A))
  # where theta(j) is the parameter for the j-th state
  
  N = dim(sum_t_vjk)[3]
  # First, get the indices of the parameters in params_vec
  beta0.id = 1:num.states
  Avec.temp = params_vec[-(beta0.id)]
  
  Amat = matrix(0,nrow=num.states,ncol=num.states)
  Amat[lower.tri(Amat)] = Avec.temp[1:(length(Avec.temp)/2)]
  Amat[upper.tri(Amat)] = Avec.temp[(length(Avec.temp)/2+1):length(Avec.temp)]
  # Amat = inv.logit.tpm(Amat)
  # Amat = Amat/apply(Amat,1,sum)
  
  # Create the (expanded) transition probability matrix
  # Make sure to transform the lambda parameter
  lambdan = exp(params_vec[beta0.id])
  # lambdan = params_vec[beta0.id] # Transformation should take place before the function
  # Find out where the 0s are so that you don't take the log and sum of these entries
  #Atemp = matrix(.1,nrow=num.states,ncol=num.states)
  #Atemp = Atemp/apply(Atemp,1,sum)
  # Btemp = create_B(lambdan,tau=rep(0,num.states),m,Amat)
  
  if(sum(!is.finite(lambdan))>0){
    return(NA)
  }
  
  B = create_B(lambdan,tau=NULL,m,Amat)
  B.dummy = create_B(rep(1,num.states),NULL,m,Amat) # Only purpose is to find the 0s
  B.0id = which(B.dummy==0)
  B0 = B[-B.0id]
  # B0[B.0id] = 1e-6 # For numerical stability
  total = 0
  for(n in 1:N){
    temp = sum_t_vjk[,,n]
    # print(sprintf("n=%g, sum_t_vjk=%g",n,sum(log(B0)*sum_t_vjk[,,n])))
    # total = total + sum(log(B0)*sum_t_vjk[,,n])
    total = total+sum(log(B0)*temp[-B.0id])
  }
  return(-total) # Negative sign since nlm() minimizes, not maximizes
}


term2.equality = function(params_vec,num.states,m,sum_t_vjk){
  # Equality constraints for term2.nocovar
  
  # First, get the indices of the parameters in params_vec
  beta0.id = 1:num.states
  Avec.temp = params_vec[-(beta0.id)]
  
  Amat = matrix(0,nrow=num.states,ncol=num.states)
  Amat[lower.tri(Amat)] = Avec.temp[1:(length(Avec.temp)/2)]
  Amat[upper.tri(Amat)] = Avec.temp[(length(Avec.temp)/2+1):length(Avec.temp)]
  
  return(rowSums(Amat))
}

expand_params = function(mu,Sigma,m){
  mu.expand = matrix(NA,nrow=sum(m),ncol=dim(mu)[2])
  Sigma.expand = array(NA,c(dim(Sigma)[1:2],sum(m)))
  for(j in 1:length(m)){
    for(k in 1:m[j]){
      if(j == 1){
        mu.expand[k,] = mu[j,]
        Sigma.expand[,,k] = Sigma[,,j]
      }
      else{
        mu.expand[sum(m[1:(j-1)])+k,] = mu[j,]
        Sigma.expand[,,sum(m[1:(j-1)])+k] = Sigma[,,j]
      }
    }
  }
  return(list(mu.expand = mu.expand,
              Sigma.expand = Sigma.expand))
}

em.hsmm = function(Xt,A,mu,Sigma,delta=NULL,
                   lambda,tau,m,
                   maxiter = 1000,tol=1e-4){
  # A - transition probability matrix
  num.states = nrow(A)
  num.Bstates = sum(m)
  TIME = nrow(Xt)
  params = params_transform(c(lambda,tau,c(A)),num.states,m)
  B = create_B(params$lambda,params$tau,m,params$A)
  mu.hat.next = vector(mode="list",length = num.Bstates)
  Sigma.hat.next = vector(mode="list",length = num.Bstates)
  llk = rep(NA,maxiter)
  
  for(iter in 1:maxiter){
    print(sprintf("Iteration %g of %g",iter,maxiter))
    mu_mat = t(matrix(unlist(mu),nrow=num.states,ncol=num.Bstates))
    Sigma_array = array(unlist(Sigma),dim=c(P,P,num.Bstates))
    fb = forward_backward_cpp(Xt,B,mu_mat,Sigma_array)
    # fb = forward_backward(Xt,B,mu,Sigma,delta)
    lalpha = fb$lalpha
    lbeta = fb$lbeta
    lallprobs = matrix(NA,nrow=TIME,ncol=num.Bstates)
    for(j in 1:num.Bstates){
      #allprobs[,j] = dmvnorm(Xt,mu[[j]],Sigma[[j]],log=FALSE)
      lallprobs[,j] = dmvnrm_arma(Xt,mu_mat[j,],Sigma_array[,,j],logd=TRUE)
    }
    
    # Get transition probabilities
    # This is the only thing different in the algorithm, since the tpm is a function of the duration distribution
    
    # Estimate the transition probability matrix
    # Get sum_t v_{jk}(t)
    llk[iter] = logsumexp(lalpha[,TIME])
    sum_t_vjk = B
    for(j in 1:num.Bstates){
      for(k in 1:num.Bstates){
        sum_t_vjk[j,k] = B[j,k]*sum(exp(lalpha[j,1:(TIME-1)] + lallprobs[2:TIME,k] + lbeta[k,2:TIME]-llk[iter]))
      }
    }
    params_vec = c(lambda,tau,c(A))
    # params.hat = params_transform(nlm(term2,params_vec,num.states,m,sum_t_vjk)$estimate,num.states,m)
    params.hat = nlminb(params_vec,term2,gradient=NULL,hessian=NULL,num.states,m,sum_t_vjk)$par
    B.next = create_B(params.hat$lambda,params.hat$tau,m,params.hat$A)
    
    # Get emission distribution parameters
    delta.next = solve(t(eye(nrow(B.next))-B.next+1),rep(1,nrow(B.next))) # Stationary distribution
    for(j in 1:num.states){
      # Average out the u's
      if(j == 1){
        u = colSums(exp(lalpha[1:m[1],] + lbeta[1:m[1],] - llk[iter]))
        for(k in 1:m[1]){
          mu.hat.next[[k]] = colSums(u*Xt)/sum(u)  
          Sigma.hat.next[[k]] = covmat_c(Xt,mu.hat.next[[k]],u)
        }
      }
      else{
        u = colSums(exp(lalpha[(m[j-1]+1):(m[j-1]+m[j]),] + lbeta[(m[j-1]+1):(m[j-1]+m[j]),] - llk[iter]))
        for(k in (m[j-1]+1):(m[j-1]+m[j])){
          mu.hat.next[[k]] = colSums(u*Xt)/sum(u)  
          Sigma.hat.next[[k]] = covmat_c(Xt,mu.hat.next[[k]],u)
        }
      }
    }
    
    # Determine convergence
    if(iter > 1){
      if(abs(llk[iter] - llk[iter-1]) < tol){
        print(sprintf("Converged after %g iterations",iter))
        return(list(mu.hat = mu.hat.next,
                    Sigma.hat = Sigma.hat.next,
                    lambda.hat = params.hat$lambda, # On original scale
                    tau.hat = params.hat$tau, # On original scale
                    A.hat = params.hat$A,
                    B.hat = B.next,
                    delta = delta.next,
                    llk = llk[1:iter]))
      }
    }
    mu = mu.hat.next
    Sigma = Sigma.hat.next
    B = B.next
    
    # Transform the variables
    lambda = log(params.hat$lambda)
    tau = log(params.hat$tau)
    A = params.hat$A
  }
  print(sprintf("No convergence after %g iterations",maxiter))
  return(NA)
}

# params.hmm.init = expand_params(result.em$mu.hat,result.em$Sigma.hat,m)
# mu = params.hmm.init$mu.expand
# Sigma = params.hmm.init$Sigma.expand
# beta = rbind(c(0,0),c(0,0))
# tau = c(1,1)

em.hsmm2 = function(Xt,Z,A,mu,Sigma,delta=NULL,
                    beta,tau,m,
                    maxiter = 1000,tol=1e-4,verbose=TRUE){
  # For multiple subjects
  # Xt is a list
  # mu and Sigma are common for all subjects, i.e., emission distributions are the same
  # The main difference is in the M-step, since maximizing the likelihood is more complicated
  # A - transition probability matrix
  N = length(Xt)
  num.states = nrow(A)
  num.Bstates = sum(m)
  TIME = nrow(Xt[[1]])
  mu.hat.next = vector(mode="list",length = num.Bstates)
  Sigma.hat.next = vector(mode="list",length = num.Bstates)
  llk = matrix(NA,nrow=N,ncol=maxiter)
  beta0 = beta[1,]
  beta1 = beta[2,]
  delta.next = list()
  
  for(iter in 1:maxiter){
    if(verbose){
      print(sprintf("Iteration %g of %g",iter,maxiter))
    }
    mu_mat = t(matrix(unlist(mu),ncol=num.Bstates))
    Sigma_array = array(unlist(Sigma),dim=c(P,P,num.Bstates))
    fb.all = vector("list",N)
    lalpha.all = vector("list",N)
    lbeta.all = vector("list",N)
    lprobs.all = vector("list",N)
    for(n in 1:N){
      B.next = create_B(exp(beta0 + beta1*Z[n]),tau,m,A)  
      fb.all[[n]] = forward_backward_cpp(Xt[,,n],B.next,mu_mat,Sigma_array)
      lalpha.all[[n]] = fb.all[[n]]$lalpha
      lbeta.all[[n]] = fb.all[[n]]$lbeta
      
      lprobs.all[[n]] = matrix(NA,nrow=TIME,ncol=num.Bstates)
      for(j in 1:num.Bstates){
        lprobs.all[[n]][,j] = dmvnrm_arma(Xt[,,n],mu[[j]],Sigma[[j]],logd=TRUE)
      }
    }
    
    # Get transition probabilities
    # This is the only thing different in the algorithm, since the tpm is a function of the duration distribution
    
    # Get sum_t v_{jk}(t)
    sum_t_vjk = array(NA,dim=c(num.Bstates,num.Bstates,N))
    for(n in 1:N){
      llk[n,iter] = logsumexp(lalpha.all[[n]][,TIME])
      B = create_B(exp(beta0 + beta1*Z[n]),tau,m,A)
      for(j in 1:num.Bstates){
        for(k in 1:num.Bstates){
          sum_t_vjk[j,k,n] = B[j,k]*sum(exp(lalpha.all[[n]][j,1:(TIME-1)] + lprobs.all[[n]][2:TIME,k] + lbeta.all[[n]][k,2:TIME]-llk[n,iter]))
        }
      }
    }
    params_vec = c(beta0,beta1,log(tau),c(A)) # Transform the variables if necessary
    # params.hat = nlm(term2.covar,params_vec,num.states,m,sum_t_vjk,Z)$estimate
    params.hat = nlminb(params_vec,term2.covar,gradient=NULL,hessian=NULL,num.states,m,sum_t_vjk,Z,
                        control = list(rel.tol=1e-6))$par
    beta0.hat = params.hat[1:num.states]
    beta1.hat = params.hat[(num.states+1):(2*num.states)]
    tau.hat = round(exp(params.hat[(2*num.states+1):(3*num.states)]))
    A.hat = matrix(params.hat[-(1:(3*num.states))],nrow=num.states,ncol=num.states)
    
    # Set up the updated emission distribution parameters
    for(k in 1:num.Bstates){
      mu.hat.next[[k]] = rep(0,P)
      Sigma.hat.next[[k]] = matrix(0,nrow=P,ncol=P)
    }
    # Get emission distribution parameters
    for(n in 1:N){
      B.next = create_B(exp(beta0.hat + beta1.hat*Z[n]),tau.hat,m,A.hat)  
      delta.next[[n]] = solve(t(eye(nrow(B.next))-B.next+1),rep(1,nrow(B.next))) # Stationary distribution
      
      for(j in 1:num.states){
        # Average out the u's
        if(j == 1){
          u = colSums(exp(lalpha.all[[n]][1:m[1],] + lbeta.all[[n]][1:m[1],] - llk[n,iter]))
          for(k in 1:m[1]){
            mu.hat.next[[k]] = mu.hat.next[[k]] + colSums(u*Xt[[n]])/sum(u)/N
            Sigma.hat.next[[k]] = Sigma.hat.next[[k]] + covmat_c(Xt[[n]],mu.hat.next[[k]],u)/N
          }
        }
        else{
          u = colSums(exp(lalpha.all[[n]][(m[j-1]+1):(m[j-1]+m[j]),] + lbeta.all[[n]][(m[j-1]+1):(m[j-1]+m[j]),] - llk[n,iter]))
          for(k in (m[j-1]+1):(m[j-1]+m[j])){
            mu.hat.next[[k]] =  mu.hat.next[[k]] + colSums(u*Xt[[n]])/sum(u)/N
            Sigma.hat.next[[k]] = Sigma.hat.next[[k]] + covmat_c(Xt[[n]],mu.hat.next[[k]],u)/N
          }
        }
      }
    }
    
    # Determine convergence
    if(iter > 1){
      if(abs(sum(llk[,iter]) - sum(llk[,iter-1])) < tol){
        print(sprintf("Converged after %g iterations",iter))
        return(list(mu.hat = mu.hat.next,
                    Sigma.hat = Sigma.hat.next,
                    beta0.hat = beta0.hat,
                    beta1.hat = beta1.hat,
                    tau.hat = tau.hat,
                    A.hat = A.hat,
                    delta = delta.next,
                    llk = llk[,1:iter]))
      }
    }
    
    # Update the variables
    mu = mu.hat.next
    Sigma = Sigma.hat.next
    beta0 = beta0.hat
    beta1 = beta1.hat
    tau = tau.hat
    A = A.hat
  }
  print(sprintf("No convergence after %g iterations",maxiter))
  return(NA)
}

HS.norm = function(A,B=NULL){
  if(is.null(B)){
    return(sqrt(sum(c(A^2))))
  }
  else{
    return(sqrt(sum(diag(t(A)%*%B))))
  }
}

get.shrinkage = function(Xt,S,u=NULL){
  # u = TIME x N array
  # Xt = list with N entries, each entry TIME x P matrix
  # S = P x P covariance matrix
  TIME = nrow(Xt[[1]])
  P = ncol(Xt[[1]])
  if(is.null(u)){
    u = rep(1,TIME)
  }
  
  mu = mean(diag(S))
  d2 = HS.norm(S-mu*eye(P))^2
  b2.temp = 0
  for(n in 1:N){
    for(time in 1:TIME){
      b2.temp = b2.temp + HS.norm(u[time,n]*t(t(Xt[[n]][time,]))%*%Xt[[n]][time,]-S)^2
    }
  }
  b2.temp = b2.temp/TIME/TIME/N
  b2 = min(b2.temp,d2)
  a2 = d2-b2
  weight = min(a2,1)
  return(list(S.hat = (1-weight)*mu*eye(P) + weight*S,
              weight = weight))
}

 #params.hmm.init = expand_params(result.em$mu.hat,result.em$Sigma.hat,m)
# mu = params.hmm.init$mu.expand
# Sigma = params.hmm.init$Sigma.expand
# beta = rbind(c(0,0),c(0,0))
# tau = c(1,1)
em.hsmm2.highdim = function(Xt,Z,A,mu,Sigma,delta=NULL,
                    beta,tau,m,
                    maxiter = 10000,tol=1e-4,penalty = 'glasso'){
  # For multiple subjects
  # Xt is a list
  # mu and Sigma are common for all subjects, i.e., emission distributions are the same
  # The main difference is in the M-step, since maximizing the likelihood is more complicated
  # A - transition probability matrix
  N = length(Xt)
  num.states = nrow(A)
  num.Bstates = sum(m)
  TIME = nrow(Xt[[1]])
  mu.hat.next = vector(mode="list",length = num.Bstates)
  Sigma.hat.next = vector(mode="list",length = num.Bstates)
  llk = matrix(NA,nrow=N,ncol=maxiter)
  beta0 = beta[1,]
  beta1 = beta[2,]
  delta.next = list()
  
  for(iter in 1:maxiter){
    print(sprintf("Iteration %g of %g",iter,maxiter))
    mu_mat = t(matrix(unlist(mu),ncol=num.Bstates))
    Sigma_array = array(unlist(Sigma),dim=c(P,P,num.Bstates))
    fb.all = vector("list",N)
    lalpha.all = vector("list",N)
    lbeta.all = vector("list",N)
    lprobs.all = vector("list",N)
    for(n in 1:N){
      B.next = create_B(exp(beta0 + beta1*Z[n]),tau,m,A)  
      fb.all[[n]] = forward_backward_cpp(Xt[[n]],B.next,mu_mat,Sigma_array)
      lalpha.all[[n]] = fb.all[[n]]$lalpha
      lbeta.all[[n]] = fb.all[[n]]$lbeta
      
      lprobs.all[[n]] = matrix(NA,nrow=TIME,ncol=num.Bstates)
      for(j in 1:num.Bstates){
        lprobs.all[[n]][,j] = dmvnrm_arma(Xt[[n]],mu[[j]],Sigma[[j]],logd=TRUE)
      }
    }
    
    # Get transition probabilities
    # This is the only thing different in the algorithm, since the tpm is a function of the duration distribution
    
    # Get sum_t v_{jk}(t)
    sum_t_vjk = array(NA,dim=c(num.Bstates,num.Bstates,N))
    for(n in 1:N){
      llk[n,iter] = logsumexp(lalpha.all[[n]][,TIME])
      for(j in 1:num.Bstates){
        for(k in 1:num.Bstates){
          B = create_B(exp(beta0 + beta1*Z[n]),tau,m,A)
          sum_t_vjk[j,k,n] = B[j,k]*sum(exp(lalpha.all[[n]][j,1:(TIME-1)] + lprobs.all[[n]][2:TIME,k] + lbeta.all[[n]][k,2:TIME]-llk[n,iter]))
        }
      }
    }
    params_vec = c(beta0,beta1,log(tau),c(A)) # Transform the variables if necessary
    # params.hat = nlm(term2.covar,params_vec,num.states,m,sum_t_vjk,Z)$estimate
    params.hat = nlminb(params_vec,term2.covar,gradient=NULL,hessian=NULL,num.states,m,sum_t_vjk,Z,
                        control = list(rel.tol=1e-6))$par
    beta0.hat = params.hat[1:num.states]
    beta1.hat = params.hat[(num.states+1):(2*num.states)]
    tau.hat = round(exp(params.hat[(2*num.states+1):(3*num.states)]))
    A.hat = matrix(params.hat[-(1:(3*num.states))],nrow=num.states,ncol=num.states)
    
    # Set up the updated emission distribution parameters
    for(k in 1:num.Bstates){
      mu.hat.next[[k]] = rep(0,P)
      Sigma.hat.next[[k]] = matrix(0,nrow=P,ncol=P)
    }
    # Get emission distribution parameters
    u.all = array(0,dim=c(TIME,num.states,N)) # probability per state
    Samp.Cov = Sigma.hat.next
    for(n in 1:N){
      B.next = create_B(exp(beta0.hat + beta1.hat*Z[n]),tau.hat,m,A.hat)  
      delta.next[[n]] = solve(t(eye(nrow(B.next))-B.next+1),rep(1,nrow(B.next))) # Stationary distribution
      for(j in 1:num.states){
        # Average out the u's
        if(j == 1){
          u = colSums(exp(lalpha.all[[n]][1:m[1],] + lbeta.all[[n]][1:m[1],] - llk[n,iter]))
          u.all[,j,n] = u
          for(k in 1:m[1]){
            mu.hat.next[[k]] = mu.hat.next[[k]] + colSums(u*Xt[[n]])/sum(u)/N
            Samp.Cov[[k]] = Samp.Cov[[k]] + covmat_c(Xt[[n]],mu.hat.next[[k]],u)/N
          }
        }
        else{
          u = colSums(exp(lalpha.all[[n]][(m[j-1]+1):(m[j-1]+m[j]),] + lbeta.all[[n]][(m[j-1]+1):(m[j-1]+m[j]),] - llk[n,iter]))
          u.all[,j,n] = u
          for(k in (m[j-1]+1):(m[j-1]+m[j])){
            mu.hat.next[[k]] =  mu.hat.next[[k]] + colSums(u*Xt[[n]])/sum(u)/N
            Samp.Cov[[k]] = Samp.Cov[[k]] + covmat_c(Xt[[n]],mu.hat.next[[k]],u)/N
          }
        }
      }
    }
    
    # Shrinkage/glasso on Covariance Matrix
    for(j in 1:num.states){
      if(j == 1){
        if(penalty == 'shrinkage'){
          S.hat = get.shrinkage(Xt,Samp.Cov[[1]],u.all[,1,])$S.hat  
        }
        if(penalty == 'glasso'){
          lambda = sqrt(2*sum(c(u.all[,j,]))*log(P))/2
          rho = 2*lambda/sum(c(u.all[,1,]))*sqrt(mean(c(u.all[,j,])))
          S.hat = glasso(Samp.Cov[[1]],rho=rho)$w  
        }
        for(k in 1:m[1]){
          Sigma.hat.next[[k]] = S.hat
        }
      }
      else{
        if(penalty == 'shrinkage'){
          S.hat = get.shrinkage(Xt,Samp.Cov[[sum(m[1:(j-1)])+1]],u.all[,j,])$S.hat 
        }
        if(penalty == 'glasso'){
          lambda = sqrt(2*sum(c(u.all[,j,]))*log(P))/2
          rho = 2*lambda/sum(c(u.all[,j,]))*sqrt(mean(c(u.all[,j,])))
          S.hat = glasso(Samp.Cov[[sum(m[1:(j-1)])+1]],rho=rho)$w  
        }
        
        for(k in (m[j-1]+1):(m[j-1]+m[j])){
          Sigma.hat.next[[k]] = S.hat
        }
      }
    }
    
    
    # Determine convergence
    if(iter > 1){
      if(abs(sum(llk[,iter]) - sum(llk[,iter-1])) < tol){
        print(sprintf("Converged after %g iterations",iter))
        return(list(mu.hat = mu.hat.next,
                    Sigma.hat = Sigma.hat.next,
                    beta0.hat = beta0.hat,
                    beta1.hat = beta1.hat,
                    tau.hat = tau.hat,
                    A.hat = A.hat,
                    delta = delta.next,
                    llk = llk[,1:iter]))
      }
    }
    
    # Update the variables
    mu = mu.hat.next
    Sigma = Sigma.hat.next
    beta0 = beta0.hat
    beta1 = beta1.hat
    tau = tau.hat
    A = A.hat
  }
  print(sprintf("No convergence after %g iterations",maxiter))
  return(NA)
}

em.hsmm2.highdim2 = function(Xt,A,mu,Sigma,delta=NULL,
                            beta,m,
                            maxiter = 10000,tol=1e-4,penalty = 'glasso',
                            verbose=FALSE){
  # For multiple subjects
  # Xt is a list
  # mu and Sigma are common for all subjects, i.e., emission distributions are the same
  # The main difference is in the M-step, since maximizing the likelihood is more complicated
  # A - transition probability matrix
  N = dim(Xt)[3]
  P = dim(Xt)[2]
  TIME = dim(Xt)[1]
  num.states = nrow(A)
  num.Bstates = sum(m)
  
  mu.hat.next = vector(mode="list",length = num.Bstates)
  Sigma.hat.next = vector(mode="list",length = num.Bstates)
  Omega.hat.next = vector(mode="list",length = num.Bstates)
  llk = matrix(NA,nrow=N,ncol=maxiter)
  beta0 = beta[1,]
  delta.next = list()
  
  for(iter in 1:maxiter){
    if(verbose){
      print(sprintf("Iteration %g of %g",iter,maxiter))
    }
    mu_mat = t(matrix(unlist(mu),ncol=num.Bstates))
    Sigma_array = array(unlist(Sigma),dim=c(P,P,num.Bstates))
    fb.all = vector("list",N)
    lalpha.all = vector("list",N)
    lbeta.all = vector("list",N)
    lprobs.all = vector("list",N)
    for(n in 1:N){
      B.next = create_B(exp(beta0),tau=NULL,m,A)  
      fb.all[[n]] = forward_backward_cpp(Xt[,,n],B.next,mu_mat,Sigma_array)
      lalpha.all[[n]] = fb.all[[n]]$lalpha
      lbeta.all[[n]] = fb.all[[n]]$lbeta
      
      lprobs.all[[n]] = matrix(NA,nrow=TIME,ncol=num.Bstates)
      for(j in 1:num.Bstates){
        lprobs.all[[n]][,j] = dmvnrm_arma(Xt[,,n],mu[[j]],Sigma[[j]],logd=TRUE)
      }
    }
    
    # Get transition probabilities
    # This is the only thing different in the algorithm, since the tpm is a function of the duration distribution
    
    # Get sum_t v_{jk}(t)
    sum_t_vjk = array(NA,dim=c(num.Bstates,num.Bstates,N))
    B = create_B(exp(beta0),tau = NULL,m,A) # Move this inside for loop if it will be a function of covariates
    for(n in 1:N){
      llk[n,iter] = logsumexp(lalpha.all[[n]][,TIME])
      for(j in 1:num.Bstates){
        for(k in 1:num.Bstates){
          sum_t_vjk[j,k,n] = B[j,k]*sum(exp(lalpha.all[[n]][j,1:(TIME-1)] + lprobs.all[[n]][2:TIME,k] + lbeta.all[[n]][k,2:TIME]-llk[n,iter]))
        }
      }
    }
    # params_vec = c(beta0,c(A)) # Transform the variables if necessary
    # Note that A is parameterized now so that the diagonals are exactly 0
    Atemp = A 
    if(num.states == 2){
      Atemp = matrix(1,nrow=2,ncol=2)
    }
    diag(Atemp) = 0
    
    Atemp = logit.tpm(Atemp/apply(Atemp,1,sum))
    params_vec = c(beta0,c(Atemp[upper.tri(Atemp)],Atemp[lower.tri(Atemp)]))
    # params.hat = nlm(term2.covar,params_vec,num.states,m,sum_t_vjk,Z)$estimate
    params.hat = nlminb(params_vec,term2.nocovar,gradient=NULL,hessian=NULL,num.states,m,sum_t_vjk,
                        control = list(rel.tol=1e-6))$par
    beta0.hat = params.hat[1:num.states]
    A.hat.vec = params.hat[-(1:num.states)]
    # A.hat = matrix(params.hat[-(1:(1*num.states))],nrow=num.states,ncol=num.states)
    A.hat = matrix(0,nrow=num.states,ncol=num.states)
    A.hat[lower.tri(A.hat)] = A.hat.vec[1:(length(A.hat.vec)/2)]
    A.hat[upper.tri(A.hat)] = A.hat.vec[(length(A.hat.vec)/2+1):length(A.hat.vec)]
    A.hat = inv.logit.tpm(A.hat)
    A.hat = A.hat/apply(A.hat,1,sum)
    
    # Set up the updated emission distribution parameters
    for(k in 1:num.Bstates){
      mu.hat.next[[k]] = rep(0,P)
      Sigma.hat.next[[k]] = matrix(0,nrow=P,ncol=P)
    }
    # Get emission distribution parameters
    u.all = array(0,dim=c(TIME,num.states,N)) # probability per state
    Samp.Cov = Sigma.hat.next
    for(n in 1:N){
      B.next = create_B(exp(beta0.hat),tau = NULL,m,A.hat)  
      delta.next[[n]] = solve(t(eye(nrow(B.next))-B.next+1),rep(1,nrow(B.next))) # Stationary distribution
      for(j in 1:num.states){
        # Average out the u's
        if(j == 1){
          u = colSums(exp(lalpha.all[[n]][1:m[1],] + lbeta.all[[n]][1:m[1],] - llk[n,iter]))
          u.all[,j,n] = u
          for(k in 1:m[1]){
            mu.hat.next[[k]] = mu.hat.next[[k]] + colSums(u*Xt[,,n])/sum(u)/N
            Samp.Cov[[k]] = Samp.Cov[[k]] + covmat_c(Xt[,,n],mu.hat.next[[k]],u)/N
          }
        }
        else{
          u = colSums(exp(lalpha.all[[n]][(sum(m[1:(j-1)])+1):sum(m[1:j]),] + lbeta.all[[n]][(sum(m[1:(j-1)])+1):sum(m[1:j]),] - llk[n,iter]))
          u.all[,j,n] = u
          for(k in (sum(m[1:(j-1)])+1):(sum(m[1:j]))){
            mu.hat.next[[k]] =  mu.hat.next[[k]] + colSums(u*Xt[,,n])/sum(u)/N
            Samp.Cov[[k]] = Samp.Cov[[k]] + covmat_c(Xt[,,n],mu.hat.next[[k]],u)/N
          }
        }
      }
    }
    
    # Shrinkage/glasso on Covariance Matrix
    for(j in 1:num.states){
      if(j == 1){
        if(penalty == 'none' | is.null(penalty)){
          S.hat = Samp.Cov[[j]]
        }
        if(penalty == 'shrinkage'){
          S.hat = get.shrinkage(Xt,Samp.Cov[[j]],u.all[,j,])$S.hat  
        }
        if(penalty == 'glasso'){
          lambda = sqrt(2*sum(c(u.all[,j,]))*log(P))/2
          rho = 2*lambda/sum(c(u.all[,j,]))*sqrt(mean(c(u.all[,j,])))
          S.hat = glasso(Samp.Cov[[j]],rho=rho)$w
          Omega.hat = glasso(Samp.Cov[[j]],rho=rho)$wi
        }
        for(k in 1:m[1]){
          Sigma.hat.next[[k]] = S.hat
          if(penalty=='glasso'){Omega.hat.next[[k]] = Omega.hat}
        }
      }
      else{
        if(penalty == 'none' | is.null(penalty)){
          S.hat = Samp.Cov[[sum(m[1:(j-1)])+1]]
        }
        if(penalty == 'shrinkage'){
          S.hat = get.shrinkage(Xt,Samp.Cov[[sum(m[1:(j-1)])+1]],u.all[,j,])$S.hat 
        }
        if(penalty == 'glasso'){
          lambda = sqrt(2*sum(c(u.all[,j,]))*log(P))/2
          rho = 2*lambda/sum(c(u.all[,j,]))*sqrt(mean(c(u.all[,j,])))
          S.hat = glasso(Samp.Cov[[sum(m[1:(j-1)])+1]],rho=rho)$w
          Omega.hat = glasso(Samp.Cov[[sum(m[1:(j-1)])+1]],rho=rho)$wi  
        }
        
        for(k in (sum(m[1:(j-1)])+1):(sum(m[1:j]))){
          Sigma.hat.next[[k]] = S.hat
          if(penalty=='glasso'){Omega.hat.next[[k]] = Omega.hat}
        }
      }
    }
    
    # Determine convergence
    if(iter > 1){
      if(abs(sum(llk[,iter]) - sum(llk[,iter-1])) < tol){
        print(sprintf("Converged after %g iterations",iter))
        if(penalty == 'shrinkage' | penalty == 'none' | is.null(penalty)){
          return(list(mu.hat = mu.hat.next,
                      Sigma.hat = Sigma.hat.next,
                      beta0.hat = beta0.hat,
                      A.hat = A.hat,
                      delta = delta.next,
                      llk = llk[,1:iter]))  
        }
        if(penalty == 'glasso'){
          return(list(mu.hat = mu.hat.next,
                      Sigma.hat = Sigma.hat.next,
                      Omega.hat = Omega.hat.next,
                      beta0.hat = beta0.hat,
                      A.hat = A.hat,
                      delta = delta.next,
                      llk = llk[,1:iter]))
        }
        
      }
    }
    
    # Update the variables
    mu = mu.hat.next
    Sigma = Sigma.hat.next
    beta0 = beta0.hat
    A = A.hat
  }
  print(sprintf("No convergence after %g iterations",maxiter))
  if(penalty == 'shrinkage'){
    return(list(mu.hat = mu.hat.next,
                Sigma.hat = Sigma.hat.next,
                beta0.hat = beta0.hat,
                A.hat = A.hat,
                delta = delta.next,
                llk = llk[,1:iter],
                fb.all = fb.all))  
  }
  if(penalty == 'glasso'){
    return(list(mu.hat = mu.hat.next,
                Sigma.hat = Sigma.hat.next,
                Omega.hat = Omega.hat.next,
                beta0.hat = beta0.hat,
                A.hat = A.hat,
                delta = delta.next,
                llk = llk[,1:iter],
                fb.all = fb.all))
  }
}

KL.divergence = function(Sigma1,Sigma2,Omega1=NULL,Omega2=NULL){
  if(is.null(Omega1)){
    Omega1 = solve(Sigma1)
  }
  if(is.null(Omega2)){
    Omega2 = solve(Sigma2)
  }
  return(sum(diag((Sigma1-Sigma2)%*%(Omega2-Omega1))))
}

get.distances = function(Sigma,Omega=NULL){
  # Sigma is a list of covariance matrices
  # Omega is a list of precision matrices
  K = dim(Sigma)[3]
  if(is.null(Omega)){
    Omega = vector("list",length = K)
    for(k in 1:K){
      Omega[[k]] = solve(Sigma[,,k])
    }
  }
  dist.mat = matrix(999,nrow=K,ncol=K)
  for(j in 1:(K-1)){
    for(k in (j+1):K){
      dist.mat[j,k] = KL.divergence(Sigma[,,j],Sigma[,,k],Omega[[j]],Omega[[k]])
      dist.mat[k,j] = dist.mat[j,k]
    }
  }
  return(dist.mat)
}

get.distances.hsmm.Sigma = function(result,m,Sigma){
  num.states = result$A.hat%>%nrow()
  dist.mat = matrix(NA,nrow=num.states,ncol=num.states)
  for(j in 1:num.states){
    for(k in 1:num.states){
      if(j == 1){
        dist.mat[j,k] = HS.norm(result$Sigma.hat[,,j]-Sigma[[k]])
      }
      else{
        dist.mat[j,k] = HS.norm(result$Sigma.hat[,,sum(m[1:(j-1)])+1]-Sigma[[k]])
      }
    }
  }
  return(dist.mat)
}

get.cost = function(result.em,penalty = 'BIC'){
  ll = colSums(result.em$llk)[length(colSums(result.em$llk))]  
  N = dim(result.em$u.all)[3]
  TIME = dim(result.em$u.all)[1]
  K = dim(result.em$u.all)[2]
  P = dim(result.em$Sigma.hat)[1]
  pi.hat = apply(apply(result.em$u.all,c(1,2),mean),2,mean)
  if(penalty == 'BIC'){
    return(-2*ll + log(N*TIME)*(K*(K-1) + K*choose(P,2)))
  }
  if(penalty == 'MMDL'){
    return(-2*ll + log(N*TIME)*K*(K-1) + sum(log(N*TIME*pi.hat)*choose(P,2)))
  }
}

IC.hmm = function(result){
  llk = extract.likelihood(result)
  num.states = nrow(result$A.hat)
  P = dim(result$Sigma.hat)[1]
  TIME = dim(result$u.all)[1]
  N = dim(result$u.all)[3]
  # num.params = num.states*(choose(P,2)+2*P) + num.states^2 # Common TPM
  num.params = num.states*(choose(P,2)+2*P) + N*num.states^2 # Subject-specific TPM
  ijhclss = matrix(NA,N,TIME) #alternatively
  iv = hmm_viterbi(result)
  for (i in 1:N) for (j in 1:TIME){
    z = iv[i,j]
    ijhclss[i,j]=-log(result$u.all[j,z,i])
  }
  hclss= sum(ijhclss)
  # hclss = -sum(c(result$u.all)*log(result$u.all),na.rm = T)
  return(list(AIC = -2*llk+2*num.params,
              BIC = -2*llk + num.params*log(TIME*N),
              ICL = -2*llk + num.params*log(TIME*N) + 2*hclss))
}

fit.ahsmm = function(Xt,
                   Z,
                   m,
                   mu.init,
                   Sigma.init,
                   A.init,
                   beta.init = NULL,
                   ATTEMPTS = 10,
                   tol = 1e-5,
                   maxiter = 250,
                   pars = list(mu = 1, sd = 0.5),
                   verbose = FALSE,
                   refit = FALSE)
  {
  params.init = list()
  params.init$A = A.init
  params.init$mu = mu.init
  params.init$Sigma = Sigma.init
  params.init$beta.init = beta.init
  
  result.hsmm.attempts = hsmm_attempts(Xt,
                                       Z,
                                       params.init = params.init,
                                       m = m,
                                       ATTEMPTS = ATTEMPTS,
                                       tol = tol,
                                       maxiter = maxiter,
                                       pars = pars,
                                       verbose = verbose,
                                       refit = refit)
  return(result.hsmm.attempts[[result.hsmm.attempts %>% 
                                 lapply(extract.likelihood) %>% 
                                 unlist() %>% 
                                 which.max()]])
}