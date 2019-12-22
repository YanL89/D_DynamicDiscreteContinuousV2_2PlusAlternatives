#' compute utilities from design matrix and cx
#' 
#' @param X design matrix
#' @param b cx
#' @param n number of observations
#' @return matrix of utilities
#' @export
getUtilities = function(X, b, n){
	matrix( X %*% b, nrow = n, ncol = nrow(X) / n, byrow = TRUE)
	}

#' generates choices for probit
#' 
#' @param generic generic parameters
#' @param specific specific parameters
#' @param b cx
#' @param D the data set
#' @param errsDiff differences between errors
#' @param SDiff covariance of error differences
#' @return probit choices
#' @export
genChoiceProbit = function(generic, specific, b, D, errsDiff = NULL, SDiff = NULL, ASC = FALSE){
	nalt = length(specific)
	n = nrow(D)
	
	X = create_X(generic, specific, D, ASC)
	U = getUtilities(X, b, n)
	# Y1 = apply(U,1,which.max)
	# substract the first column (first alternative) to all other columns
	Udiff = U[,-1] - U[,1]
	
	if(is.null(errsDiff)){
		errsDiff = t(chol(SDiff)) %*% matrix(rnorm((nalt-1) * n), nalt -1, n)
		errsDiff = t(errsDiff)
		}
		
	#errsDiff = errsDiff / 0.01 
	Udiff = Udiff + errsDiff
	Udiff = cbind(rep(0,n), Udiff)
	Y2 = apply(Udiff,1,which.max)
	# mean(Y1 == Y2)
	Y2
	}

# dc4Params = function(x,nparPr,nparReg){
# 	betaPr = x[1:nparPr]
# 	betaReg = x[(nparPr + 1):(nparPr + nparReg)]
# 	L = vec2mat(c(sqrt(2),x[(nparPr + nparReg + 1):length(x)]))
# 	return(list(betaPr = betaPr, betaReg = betaReg, L = L))
# 	}

#' turns a lower triangular matrix into a vector of its elements
#' 
#' @param x the matrix
#' @return the vector of cx
#' @export
mat2vec = function(x){
	v = c()
	for(i in 1:nrow(x))
		v = c(v,x[i,1:i])
	v	
	}

#' turns a vector that contains the elements of a lower triangular
#' matrix into the matrix itself
#' 
#' @param x the vector
#' @return the matrix
#' @export
vec2mat = function(x){
	n = round((-1 + sqrt(1 + 8*length(x))) / 2 + 0.1)
	L = matrix(0,n,n)
	index = 1
	for(i in 1:n){
		for(j in 1:i){
			L[i,j] = x[index]
			index = index + 1
			}
		}
	L
	}
#' computes conditional covariance in multivariate normal
#' 
#' calculates the conditional covariance matrix of 
#' normal vector conditional on observing observation
#' INDEX posObs. Please note that the conditional
#' covariance does not depend on the observation itself
#' 
#' @param sigma the original covariance
#' @param posObs the position observed
#' @return the conditional covariance
#' @export
condCov = function(sigma, posObs){
	sigma11 = sigma[-posObs,-posObs]
	sigma22 = sigma[posObs,posObs]
	sigma12 = sigma[-posObs,posObs]
	sigma21 = sigma[posObs,-posObs]
	
	sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
	}

#' conditional expectation in narmal
#' 
#' calculate the conditional mean of a random vector
#' conditional on observing obs at index posObs
#' 
#' @param sigma original covariance
#' @param mu original expectation
#' @param posObs positios of observation
#' @param obs value of observatios
condMean = function(sigma, mu, posObs, obs){
	sigma12 = sigma[-posObs,posObs]
	sigma22 = sigma[posObs,posObs]
	mu1 = mu[-posObs]
	mu2 = mu[posObs]
	
	mu1 + sigma12 %*% solve(sigma22) %*% (obs - mu2)
	}

#' computes transformation matrix to get error differences
#' 
#' calculate the covariance of the differences between 
#' a random vector and one of its component (indexed by
#' baseUt) 
#' for example: returns the matrix that turn (e1, e2, e3, e4) into
#' (e2 - e1, e3 - e1, e4 - e1) for n = 4 and baseUt = 1
#' 
#' @param n dimension of vector
#' @param baseUt index to take differences against
#' @return tranformation matrix
#' @export
cov2diff = function(n, baseUt){
	ndiff = n - 1
	
	M = matrix(0, nrow = ndiff, ncol = n)
	M[,baseUt] = -1
	M[,-baseUt] = diag(ndiff)
		
	M
	}

#' compute tranformation matrix to change what index differences are taken against
#' 
#' generates a matrix M that transforms differences with
#' respect to error term index currentBase to differences
#' with respect to index newBase.
#' returns identity if the two bases are the same
#' @param n dimension of diff vector
#' @param currentBase index current differences are taken against
#' @param newBase index new differences should be taken against
#' @return the transforwation matrix
#' @export	
reparamErrors = function(n,currentBase, newBase){
	if(currentBase == newBase)
		return(diag(n))
	
	minus1_col = newBase - (newBase > currentBase)
	zero_row = currentBase - (currentBase > newBase)
	
	M = matrix(0,n,n)
	M[,minus1_col] = -1
	M[-zero_row,-minus1_col] = diag(n-1)
	
	M
	}

#' computes the likelihood of one probit observation with sim or genz
#' 
#' @param args arguments
#' @param genz true if genz should be used, otherwise sims need to be included in args
#' @return estimated likelihood
#' @export
probitL = function(args, genz = TRUE){
	if(genz)
		return(pmvProbitL(args))
	
	return(simProbitL(args))
	}

#' computes the likelihood of one probit observation with simulation
#' 
#' @param args arguments
#' @return estimated likelihood
#' @export
simProbitL = function(args){
	M = reparamErrors(args$nalt - 1, currentBase = 1, newBase = args$choice)
	muReparam = as.vector(M %*% args$mu)
	SReparam = M %*% args$S %*% t(M)
	
	#UDiff = args$U[- args$choice] - args$U[args$choice]
	UDiff = args$U[ args$choice ] - args$U[ - args$choice ]
  set.seed(1234)
	pmvnormSim(upper = UDiff,mean = muReparam, sigma = SReparam, draws = args$draws)
	}

#' MV normal proba with simulation p = P(N(mean,sigma) < upper)
#' 
#' @param upper upper bound on X
#' @param mean mean of X
#' @param sigma covariance of X
#' @param draws simulation draws matrix
#' @return estimated probability
#' @export
pmvnormSim = function(upper ,mean , sigma , draws){
	X = (t(chol(sigma)) %*% draws) + mean
	mean(colSums(X < upper) == nrow(sigma))
}

#' computes the likelihood of one probit observation with Genz algorithm
#' 
#' @param args arguments
#' @return estimated likelihood
#' @export
pmvProbitL = function(args){

	M = reparamErrors(args$nalt - 1, currentBase = 1, newBase = args$choice)
	muReparam = as.vector(M %*% args$mu)
	SReparam = M %*% args$S %*% t(M)
	
	#UDiff = args$U[- args$choice] - args$U[args$choice]
	UDiff = args$U[ args$choice ] - args$U[ - args$choice ]
  set.seed(1234)
	pmvnorm(upper = UDiff,mean = muReparam, sigma = SReparam)
	}

#' gets dc4 params
#' 
#' @param x cx all together in a vector
#' @param nparPr number of probit predictors
#' @param nparReg number of regression predictors
#' @return list of dc4 pamameters
#' @export
dc4Params = function(x,nparPr,nparReg){
minElem = 10e-6
  betaPr = x[1:nparPr]
  betaReg = x[(nparPr + 1):(nparPr + nparReg)]
  Lelem = x[(nparPr + nparReg + 1):length(x)]
#Lelem[Lelem < minElem] = minElem
  L = vec2mat(c(1,Lelem))
  return(list(betaPr = betaPr, betaReg = betaReg, L = L))
}
