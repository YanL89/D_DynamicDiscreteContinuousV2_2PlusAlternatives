#' generates a bunch of samples from the model specification and reports 
#' coverage properties of the MLE estimator
#' 
#' @param spec the model specification
#' @param simSpec specifications of the simulation
#' @export
#' @return a table with statistics
dynCheckAsympProperties = function(spec, simSpec){
	oldSpec = spec
  spec = checkFillDynSpec(oldSpec)
  B = simSpec$B
  npar = length(spec$b)
#  betaHatMat = matrix(0, npar, B)
 # sdMat = matrix(0, npar, B)
  
  dynList_dynGen = list()
	logitList_dynGen = list()

	dynList_logitGen = list()
	logitList_logitGen = list()

	logitB = spec$logitB
  for(b in 1:B){
    cat("iteration number: ")
    print(b)

		# dyn estimation
    C = simDynamic(spec)
    spec$C = C$Ch
    dynList_dynGen[[b]] = model(dyn, spec, NULL)
    
    # translate to logit
    toLogit = dyn2logitV2(spec)
    toLogit$specLogit$SD = "none"
    logitList_dynGen[[b]] = model(logit, toLogit$specLogit, toLogit$D)
    
    # gen from logit
    # dyn est
    C = genLogit(toLogit$specLogit$generic, toLogit$specLogit$specific, toLogit$D, logitB) -1
    newSpec = oldSpec
		newSpec$C = matrix(C, nrow = dim(spec$C) [1],  ncol = dim(spec$C) [2], byrow=FALSE)
  	spec = checkFillDynSpec(newSpec)
  	
  	modelFns = dyn
  	D = NULL		
		dynList_logitGen[[b]] = model(dyn, spec, NULL)
		
		# logit est
		toLogit$D$logitY = C
	  logitList_logitGen[[b]] = model(logit, toLogit$specLogit, toLogit$D)	
  }
  list(dynGenDynEst = dynList_dynGen
  	, dynGenLogitEst = logitList_dynGen
  	, logitGenDynEst =dynList_logitGen
  	, logitGenLogitEst = logitList_logitGen)
  #covers = rep(0, spec$b)
  #q = qnorm(1 - simSpec$alpha / 2)
  
  #for(i in 1:B){
  #  covers = covers + 1*((betaHatMat[,i] - q * sdMat[,i]) < spec$b &
  #                         (betaHatMat[,i] + q * sdMat[,i]) > spec$b)
  #}
  
  #par = c(spec$payoff_alt, spec$payoff_time, spec$payoff_global)
  #D = spec$modifyD(spec$D, spec$Time, spec$global, spec$First, 1)
  #par = c(par, getNames(spec, D))
  
  #data.frame(
  #  par = par,  
  #  b = spec$b,
  #  est_mean = rowMeans(betaHatMat)    
  #  )
}

compileDyn = function(L){
	B = length(L[[1]])
	nr = length(L[[1]][[1]]$results$beta_hat)
	
	MDynDyn = matrix(0,nr,B)
	MDynLogit = matrix(0,nr,B)
	MLogitDyn = matrix(0,nr,B)
	MLogitLogit = matrix(0,nr,B)
	
	for(b in 1:B){
		MDynDyn[,b] = L$dynGenDynEst[[b]]$results$beta_hat
		MDynLogit[,b] = L$dynGenLogitEst[[b]]$results$beta_hat
		MLogitDyn[,b] = L$logitGenDynEst[[b]]$results$beta_hat
		MLogitLogit[,b] = L$logitGenLogitEst[[b]]$results$beta_hat
	}
	list(MDynDyn = MDynDyn
			, MDynLogit = MDynLogit
			, MLogitDyn = MLogitDyn
			, MLogitLogit = MLogitLogit)
}

#' simulates choices of the dynamic model
#' 
#' @param spec the model specification. must include parameters
#' @export
#' @return a list with the choices, W, U and R
simDynamic = function(spec){
  nparH = length(spec$payoff_alt) + length(spec$payoff_time) 
  + length(spec$payoff_global)
  # parameter of payoff utilities
  alpha = spec$b[1:nparH]
  
  # parameter of (normal) utilities
  beta = spec$b[(nparH+1):length(spec$b)]
  
  Ch = matrix(0, spec$nObs, spec$nTime)
  
  Wlist = list()
  Rlist = list()
  Ulist = list()
  Z = spec$First
  for(t in 1:spec$nTime){
    M = list()
    H = list()
    X = list()
    U = list()
    V = list()
    for(l in 1:(spec$nLook+1)){
      #M[[l]] = getM(Z, spec$Time, spec$Global, t + l - 1, spec$payoff_alt
      #              , spec$payoff_time, spec$payoff_global)
      M[[l]] = getM(Z, t + l - 1, spec)

      H[[l]] = as.matrix(M[[l]]) %*% as.vector(alpha)
      if(0 == ncol(H[[l]]))
      	H[[l]] = rep(0,nrow(H[[l]]))
      
      Dt = spec$modifyD(spec$D, spec$Time, spec$Global, Z, t + l - 1)
      X[[l]] = create_X(spec$generic, spec$specific, Dt)
      U[[l]] = matrix(X[[l]] %*% beta, nrow = spec$nObs, ncol = spec$nAlt
                      , byrow = TRUE)
      V[[l]] = UtoV(U[[l]])
    }
    
    W = getW(H, V)
    R = getR(U[[1]])
    Wlist[[t]] = W
    Rlist[[t]] = R
    Ulist[[t]] = U[[1]]
    
    choices = getDynChoice(W,R,U[[1]])
    Ch[,t] = choices
    # now we also need to update Z for those that changed their alternative
    Z = updateZ(Z,spec$D[[t]],choices,spec$varNames)
  }
  list(Ch = Ch, W = Wlist, R = Rlist,U = Ulist)
}

#' updates Z for the new choices
#' 
#' @param Z the current Z to be updated
#' @param D alternatives attributes
#' @param choices the vector of choices at current time period
#' @param altVarNames a vector of attributes names for alternatives
#' @export
#' @return the new updated Z
updateZ = function(Z, D, choices, altVarNames){
  for(i in 1:nrow(Z))
    # if choice == 0 the person has not changed his item so his Z
    # does not need to be updated
    if(choices[i] != 0) 
      Z[i,altVarNames] = D[i,paste(altVarNames,choices[i],sep=".")]
  Z
}

#
#' select variables from the time data set
#' 
#' @param Time the time data set
#' @param vars the variables to get 
#' @param t the time period number
#' @param all whether all variables should be selected
#' @export
#' @return the corresponding columns in the data frame
Time_select_vars = function (Time, vars, t, all = FALSE){
	if(TRUE == all){
		names = unlist(strsplit(names(Time), ".", fixed = TRUE))
		names = names[seq(from = 1, by = 2, to = length(names))]
		vars = unique(names)
	}    #get names of variables in "time" data 
	if(0 == length(vars))
		return(data.frame(matrix(0,nrow(Time),0)))
  TimeSelect = data.frame(Time[,paste(vars,t,sep=".")])   #derive time data at t
  names(TimeSelect) = vars
  TimeSelect
}

#
#' computes M for the dynamic model
#' 
#' @param Z 
#' @param Time the time data set
#' @param Global the global data set
#' @param t the time period ID
#' @param payoff_alt the variables from the alternative
#' @param payoff_time the variables from time
#' @param payoff_global the glabal variables
#' @export
#' @return the vector M
#' 
#getM = function(Z,Time,Global, t, payoff_alt, payoff_time, payoff_global){
#  Timet = Time_select_vars(Time, payoff_time, t)
#  cbind(Z[,payoff_alt], Timet, Global[, payoff_global])
#}

#' computes M for the dynamic model
#' 
#' @param Z 
#' @param tPlusLook time period (including look ahead)
#' @param spec model specification
#' @export
#' @return the vector M
#' Yan: get data of attributes for pay-off 
getM = function(Z, tPlusLook, spec){
  Timet = Time_select_vars(spec$Time, spec$payoff_time, tPlusLook)   # attributes from Time data
  cbind(Z[,spec$payoff_alt], Timet, spec$Global[, spec$payoff_global])  # combine attributes from three data files
}


#' compute reservation utilities W
#' 
#' @param H the vector of one period payoffs
#' @param V the vector of max utilities (nu)
#' @export
#' @return W
getW = function(H, V){
  l = length(H)
  W = maxVec(V[[l]], H[[l]])
  
  l = l - 1
  while(l > 1){
    W = maxVec(V[[l]], H[[l]] + W) 
    l = l - 1
  }
  H[[1]] + W
}

#' computes tho mode of nu
#' 
#' @param V deterministic part of utilities
#' @export
#' @return the mode R
getR = function(V){
  log(rowSums(exp(V)))
}


#add by Yan
# compute final utility - fU of recursive logit
getfU = function(V, transition){
  l = length(V)
  expV = exp(V[[l]])
  R = log(expV %*% transition)
  
  l = l - 1
  while(l > 1){
    V[[l]] = V[[l]] + R
    expV = exp(V[[l]])
    R = log(expV %*% transition)
    l = l - 1
  }
  V[[1]] + R
}


#add by Yan
# compute final utility - fU of recursive probit
getProbitU = function(V, S, transition){
  l = length(V)
  n = nrow(V[[l]])
  m = ncol(V[[l]])
  R = matrix(0, n, m) 
  
  for(i in 1:n){
    mu1 = V[[l]][i, 1]
    mu2 = V[[l]][i, 2]
    q1 = (mu1 - mu2)/sqrt(S)
    q2 = (mu2 - mu1)/sqrt(S)
    expX = mu1*pnorm(q1) + mu2*pnorm(q2) + sqrt(S)*dnorm(q1)
    #refill value in R[i, 1]
    if(transition[1, 1]==1 && transition[1, 2]==1){
      R[i, 1] = expX
    } else if(transition[1, 1]==1 && transition[1, 2]==0){
      R[i, 1] = mu1  
    } else {
      R[i, 1] = mu2}
    #refill value in R[i, 2]
    if(transition[2, 1]==1 && transition[2, 2]==1){
      R[i, 2] = expX
    } else if(transition[2, 1]==1 && transition[2, 2]==0){
      R[i, 2] = mu1  
    } else {
      R[i, 2] = mu2}    
  }
  
  l = l - 1
  while(l > 1){
    V[[l]] = V[[l]] + R
    
    for(i in 1:n){
      mu1 = V[[l]][i, 1]
      mu2 = V[[l]][i, 2]
      q1 = (mu1 - mu2)/sqrt(S)
      q2 = (mu2 - mu1)/sqrt(S)
      expX = mu1*pnorm(q1) + mu2*pnorm(q2) + sqrt(S)*dnorm(q1)
      #refill value in R[i, 1]
      if(transition[1, 1]==1 && transition[1, 2]==1){
        R[i, 1] = expX
      } else if(transition[1, 1]==1 && transition[1, 2]==0){
        R[i, 1] = mu1  
      } else {
        R[i, 1] = mu2}
      #refill value in R[i, 2]
      if(transition[2, 1]==1 && transition[2, 2]==1){
        R[i, 2] = expX
      } else if(transition[2, 1]==1 && transition[2, 2]==0){
        R[i, 2] = mu1  
      } else {
        R[i, 2] = mu2}    
    }
    l = l - 1
  }
  V[[1]] + R
}


#' transforms U to nu = max(U)
#' 
#' @param U the utilities
#' @export
#' @return nu
UtoV = function(U){
  apply(U,1,max)
}

#' element-wise max of two vectors
#' 
#' @param v1 a vector
#' @param v2 another vector
#' @export
#' @return  v where v_i = max(v1_i, v2_i)
maxVec = function(v1,v2){
  apply(cbind(v1,v2),1,max)
}

#' element-wise min of two vectors
#' 
#' @param v1 a vector
#' @param v2 another vector
#' @export
#' @return  v where v_i = min(v1_i, v2_i)
minVec = function(v1,v2){
  apply(cbind(v1,v2),1,min)
}

#' calcunates the probability of keeping
#' 
#' @param W
#' @param R
#' @export
#' @return P(keep)
getPKeep = function(W,R){
  exp(-exp(-(W-R)))
}

#' generates a logit choice
#' 
#' @param U a matrix of utilities
#' @export 
#' @return a logit choice
getChoiceLogit = function(U){
  U = U + matrix(rgumbel(nrow(U) * ncol(U)), nrow(U), ncol(U))
  apply(U,1,which.max)
}

#' generate a dynamic choice for one time period
#' 
#' @param W
#' @param R
#' @param U
#' @export
#' @return a dynamic choice (for one time period)
getDynChoice = function(W, R, U){
  PK = getPKeep(W,R)
  C = rep(-1,length(W))
  
  CLogit = getChoiceLogit(U)
  for(i in 1:length(W)){
    draw = runif(1)
    if(draw <= PK[i])
      C[i] = 0
    if(draw > PK[i])
      C[i] = CLogit[i] 
  }
  C
}

#' calculates P(choire) for one time peried of the dynamic model
#' 
#' @param W
#' @param R
#' @param U
#' @param Choices
#' @export 
#' @return P(choices)
getDynProbas = function(W,R,U,Choices){
  P = rep(0, length(W))
  
  PK = getPKeep(W,R)
  PA = getPAlts(U)
  P[Choices == 0] = PK[Choices == 0]
  for(a in 1:ncol(U))
    P[Choices == a] = (1 - PK[Choices == a]) * PA[Choices == a, a]
  P
}


#add by Yan
getRLProbas = function(U,Choices){
  P = rep(0, nrow(U))
  
  PA = getPAlts(U)
  for(a in 0:ncol(U)-1)
    P[Choices == a] = PA[Choices == a, a+1]
  P
}


#add by Yan
#M is transformation matrix for variance-in-difference 
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


#add by Yan
pmvProbitL = function(args){
  M = reparamErrors(args$nalt - 1, currentBase = 1, newBase = args$choice)
  muReparam = as.vector(M %*% args$mu)
  SReparam = M %*% args$S %*% t(M)
  
  UDiff = args$U[ args$choice ] - args$U[ - args$choice ]
  #add for multivaraiate case
  for(i in 1:length(UDiff)){
    if(is.na(UDiff[i])==1){
      UDiff[i]=0
    }
  }

  set.seed(1234)
  if(ncol(SReparam)==1){
    pnorm(UDiff, mean = muReparam, sd = SReparam)
  } else {
    pmvnorm(upper = UDiff,mean = muReparam, sigma = SReparam)
  }
}


dc4Params = function(x,nvarPr,nvarReg, spec){
  betaPr = x[1:nvarPr]
  betaReg = x[(nvarPr + 1):(nvarPr + nvarReg)]
  
  ndim = spec$nAlt + 1
  numL = ndim*(ndim-1)/2-1
  L=list()
  for(t in 1:spec$nTime){
    L[[t]] = vec2mat(c(sqrt(2),x[(nvarPr + nvarReg + t*numL - 4):(nvarPr + nvarReg + t*numL)]))
  }
  return(list(betaPr = betaPr, betaReg = betaReg, L = L))
}
