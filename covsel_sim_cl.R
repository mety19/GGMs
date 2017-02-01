require(MASS)
require(pracma)
require(matrixcalc)
require(OpenMx)
require(GGMselect)

#1. CREATE MATRIX THETA - INV OF COV MATRIX
p = 200
n = 600
B = 100

for(b in 1:B){
  Theta = diag(p) #create Ipxp matrix
  
  for(i in 1:n){
    i = sample(1:p, 1) #select random row index
    j = sample(1:(i-1), 1) #select random col index in upper diag
    Theta[i, j] = Theta[j, i] = runif(1, 0, 0.5) #set upper and lower diag (symmetry) to same rand value from unif(0, 0.5) dist
  }
  
  E = eigen(Theta)
  Theta = E$vectors %*% diag(E$values -1.5 * min(E$values)) %*% #to keep evalues positive
    t(E$vectors)
  Theta = zapsmall(cov2cor(Theta))
}
is.positive.definite(Theta)

#2.  COVARIANCE MATRIX
Sigma = solve(Theta)

#3. SAMPLE FROM NORMAL DISTRIBUTION WITH COVARIANCE MATRIX AS IN 2.
spl = mvrnorm(n, rep(0,p), Sigma)

#4. ESTIMATION
estim_my = covsel(spl,0.1, 1, 1)

#5. COMPARISON

estim_gl = glasso(cov(spl), rho=0.09, thr=0.001, maxit=1000)$wi

norm(estim_gl-Theta)/norm(Theta)
norm(estim_my-Theta)/norm(Theta)

norm(estim_my-estim_gl)/norm(estim_gl)


