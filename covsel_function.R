#Sparse inverse covariance selection function
require(CVTuningCov)

setwd('C:/Users/Emma Zohner/Documents/IndependentStudy')

covsel = function(D, lambda, rho, alpha){

  objective = function(S, X, Z, lambda){
    obj = sum(diag(S*Z)) - log(det(X)) + lambda*mynorm(Z) 
    return(obj)
  }
  
  shrinkage = function(a, kp){
    y = sign(a) * pmax(0, abs(a) - lambda)
  }
  
  mynorm = function(a){
    sum(abs(a))
  }

  #Global constants and defaults
  QUIET = 0; MAX_ITER = 1000; ABSTOL = 10^(-4); RELTOL = 10^(-2);
  history = data.frame(objval=c(1:k), r_norm=c(1:k), s_norm=c(1:k), eps_pri=c(1:k), eps_dual=c(1:k))
  
  #Data preprocessing
  S = cov(D)
  n = nrow(S)
  
  #ADMM solver
  X = matrix(0, n, n)
  Z = matrix(0, n, n)
  U = matrix(0, n, n)
  
  if (QUIET == 0) {write(c('iter', 'r norm', 'eps pri', 's norm', 'eps dual', 'objective'))}
  
  for (k in 1:MAX_ITER){
    #x-update
    E <- eigen(rho*(Z-U) - S)
    Q = Re(E$vectors)
    l = Re(E$values)
    L = diag(l)
    xi = (l + sqrt(l^2 + 4*rho))/(2*rho)
    X = Q%*%diag(xi)%*%t(Q)
    
    #z-update
    Zold = Z
    X_hat = alpha*X + (1 - alpha)*Zold
    Z = shrinkage(X_hat + U, lambda/rho)
    U = U + (X_hat - Z)
    
    #diagnostics, reporting, termination checks
    history$objval[k] = objective(S, X, Z, lambda)
    history$r_norm[k]  = norm(X - Z, "F");
    history$s_norm[k]  = norm(-rho*(Z - Zold),"F");
    
    history$eps_pri[k] = sqrt(n*n)*ABSTOL + RELTOL*max(norm(X,"F"), norm(Z,"F"));
    history$eps_dual[k]= sqrt(n*n)*ABSTOL + RELTOL*norm(rho*U,"F");
    
    
   
    
    if (history$r_norm[k] < history$eps_pri[k] && history$s_norm[k] < history$eps_dual[k]){break}

    
  }
  
    return(Z)
  
    #write.csv(history, 'covselhist.csv')

}


