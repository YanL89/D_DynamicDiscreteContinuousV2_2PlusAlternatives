source("dcmodels/R/dynUtils.R")
#' dynamic discrete-continuous model functions
dyn_dc = list(
  LLVec = function(b, args){
    params = dc4Params(b, nvarPr = ncol(args$X[[1]][[1]]),nvarReg = ncol(args$XReg[[1]]), args$spec)
    betaPr = params$betaPr
    betaReg = params$betaReg
    
    nReg = 1
    nalt = args$spec$nAlt
    ndim = nalt - 1 + nReg
    n = args$spec$nObs
    
    LL = 0
    PAllobs = rep(1, args$spec$nObs)
    PTime = matrix(1, args$spec$nObs, args$spec$nTime)
    
    U = list()

    for(t in 1:args$spec$nTime){
      L = params$L[[t]]
      S = L %*% t(L)
      
      LLReg=rep(1,args$spec$nObs)
      LLPr=rep(0,args$spec$nObs)
      
      #calculate marginal probability for continuous variable
      YReg = args$M[, t]   #observed Y
      YReg_hat = args$XReg[[t]] %*% as.matrix(betaReg)
      SReg = S[ndim, ndim]
      sd = sqrt(SReg)
      for(i in 1:args$spec$nObs){
        preg = -1
        preg = dnorm(YReg[i],YReg_hat[i],sd)
        #gurantee LLReg[i]<=1
        if(preg <= 1){
          LLReg[i] = preg
        }
      }
      
      #calculate conditional probability for discrete variable
      #step 1: calculate recursive utility at time t
      for(l in 1:(args$spec$nLook+1))
        U[[l]] = matrix(args$X[[t]][[l]] %*% betaPr, n, nalt, byrow = T)
      
      S_binary = S[1, 1]
      if(nalt > 2){
        fU = getfU(U, args$spec$transition)  #logsum of recursive logit
        } else if(nalt == 2) {
        fU = getProbitU(U, S_binary, args$spec$transition)  #obtain instant utility plus downstream utility at  time t
        }
      
      #step 2: load discrete choices at time t
      YPr = args$C[,t]
      
      #step 3: calculate mean and variance of conditional multivariate normal distribution
      SPr = S[-ndim, -ndim]
      S21 = as.matrix(S[ndim, -ndim])
      if(ncol(S21) == 1)
        S21 = t(S21)
      S12 = t(S21)
      res = as.matrix(YReg - YReg_hat)   #epsilon
      
      mu1 = as.matrix(rep(0,nalt -1))
      mu2 = as.matrix(rep(0,nReg))
      
      #step 4: calculate conditional probabilities for each individual
      for(i in 1:n){
        muCond = mu1 + S12 %*% solve(SReg) %*% (res[i] - mu2)
        Scond = SPr - S12 %*% solve(SReg) %*% S21
        LLPr[i] = pmvProbitL(list(U = fU[i,], choice = YPr[i], SCond = Scond
                            ,muCond = muCond, nalt = nalt, method = args$method))
      }
      LLPr = maxVec(LLPr,args$stop[,t])  #if out_of_market (alter=1), let the probability = 1
      PTime[, t] = LLReg * LLPr
      PAllobs = PAllobs * PTime[, t]
    }
    ##Yan add for loop and print(LL)
    for(i in 1:args$spec$nObs){
    LL = LL + log(PAllobs[i])
    }
    #print(LL)
    PAllobs
  },
  
 
  #Yan add
  computeArgs_Yan = function(spec, D){
    X = list()
    XReg = list()
    for(t in 1:spec$nTime){
      #data for regression
      DReg = spec$modifyD(spec$D, spec$Global, t)
      XReg[[t]] = as.matrix(DReg[, spec$reg])
      #data for probit
      X[[t]] = list()
      for(l in 1:(spec$nLook+1)){
        Dt = spec$modifyD(spec$D, spec$Global, t+l-1)   #get all attributes at t
        X[[t]][[l]] = create_X(spec$generic, spec$specific, Dt)
      }
    }
    #choices for probit
    C = spec$C
    C = C - min(C) + 1
    #vmt for regression
    M = spec$M
    # check if we reached stop point already
    stopMat_Yan = matrix(0, spec$nObs, spec$nTime) 
    delta = getElem(spec, name = "delta", default = 0.1)
    
    list(spec = spec, X = X, C = C, XReg = XReg, M = M, stop = stopMat_Yan, delta = delta, method = spec$method)    #Yan add
  },
  
  
  computeStart = function(spec, D){
    L = model(dyn, spec, D) 
    betaPr = L$results$beta_hat / (pi / sqrt(6))
    
    nalt = length(spec$specific)
    betaReg = rep(0, length(spec$reg))
    L = list()
    
    for(t in 1: spec$nTime){
    formReg = as.formula(paste(paste("vmt",t,sep=""), paste(spec$reg,collapse="+")
                               , sep = " ~ 0 +")) 
    D_Reg = cbind(spec$M[t], spec$D[[t]], spec$Global)
    R_Reg = lm(formReg, D_Reg)
    betaReg = betaReg + R_Reg$coefficients/spec$nTime  #average beta over time
    sigmaReg = mean(R_Reg$residuals ** 2)   #sigma22
    sigma2 = as.matrix(sigmaReg) 
    nReg = length(sigma2)

    ndim = nalt - 1 + nReg 
    S = diag(ndim)
    S[ndim, ndim] = sigmaReg
    # if errors have covariance identity
    # then the differences have covariance identity + 1
    S[-ndim, -ndim] = S[-ndim, -ndim] + 1     # different covariance matrix
    L[[t]] = t(chol(S))         #decholesky 
    
  }
    # the first element is set, we don't estimate it
    n = c(betaPr, betaReg)
    for(t in 1:spec$nTime){
    n = c(n, mat2vec(L[[t]])[-1])
  }
  n
  },
  
  computeMisc = function(spec, D){
    nalt = length(spec$specific)
    nReg = 1
    ntotaldim = nalt + nReg
    
    n = c()
    comIndex = 1
    for(com in spec$generic){
      n = c(n,paste("common",comIndex,sep="_"))
      comIndex = comIndex + 1
    }
    for(s in spec$specific)
      n = c(n,s)
    for(s in spec$reg)
      n = c(n,s)
    if(ntotaldim > 2){
      for(t in 1:spec$nTime){
        for(i in 2:(ntotaldim-1))
          for(j in 1:i){
            n = c(n, paste("L_",i,j,"_t",t,sep=""))
          }
      }
    }
    list(names = n)
  },
  
  computeSimData = function(args, D){
    list()
  }
) # dyn function list

dynApply = function(M,spec){
	print("YO NAYEL DO YOU SEE THAT")
  args = dyn$computeArgs(spec,NULL)
  b = M$results$beta_hat
  predCount = matrix(0, spec$nTime, spec$nAlt+1)
  
  # DIRTY copy from LL function :(
  nparH = ncol(args$M[[1]][[1]])
  alpha = b[1:nparH]
  beta = b[(nparH+1):length(b)]
  
  H = list()
  U = list()
  V = list()
  for(t in 1:args$spec$nTime){
    for(l in 1:(args$spec$nLook+1)){
      H[[l]] = as.matrix(args$M[[t]][[l]]) %*% alpha
      if(0 == ncol(H[[l]]))
      	H[[l]] = rep(0, args$spec$nObs)
      U[[l]] = matrix(args$X[[t]][[l]] %*% beta, args$spec$nObs, args$spec$nAlt, byrow = T)
      V[[l]] = UtoV(U[[l]])
    }
    W = getW(H,V)
    R = getR(U[[1]])
    #P = getDynProbas(W,R,U[[1]],args$C[,t])
    for(choice in 0:spec$nAlt){
    	P = getDynProbas(W,R,U[[1]],rep(choice,spec$nObs))
    	P = minVec(P,1 - args$stop[,t])
      predCount[t,choice+1] = sum(P)    
    	}

  }
  predCount
}


dynApply_Yan = function(M,spec){
  print("Hey Yan DO YOU SEE THAT")
  args = dyn$computeArgs_Yan(spec,NULL)
  b = M$results$beta_hat
  predCount = matrix(0, spec$nTime, spec$nAlt+1)
  
  # DIRTY copy from LL function :(
  nparH = ncol(args$M[[1]][[1]])
  alpha = b[1:nparH]
  beta = b[(nparH+1):length(b)]
  
  H = list()
  U = list()
  V = list()
  for(t in 1:args$spec$nTime){
    for(l in 1:(args$spec$nLook+1)){
      H[[l]] = as.matrix(args$M[[t]][[l]]) %*% alpha
      if(0 == ncol(H[[l]]))
        H[[l]] = rep(0, args$spec$nObs)
      U[[l]] = matrix(args$X[[t]][[l]] %*% beta, args$spec$nObs, args$spec$nAlt, byrow = T)
      V[[l]] = UtoV(U[[l]])
    }
    W = getW(H,V)
    R = getR(U[[1]])
    #P = getDynProbas(W,R,U[[1]],args$C[,t])
    for(choice in 0:spec$nAlt){
      P = getDynProbas(W,R,U[[1]],rep(choice,spec$nObs))
      P = minVec(P,1 - args$stop[,t])
      predCount[t,choice+1] = sum(P)    
    }
    
  }
  predCount
}


dynCount = function(spec){
  count = matrix(0, spec$nTime, spec$nAlt + 1)
  for(t in 1:spec$nTime)
    for(a in 0:spec$nAlt)
      count[t,a+1] = sum(spec$C[,t] == a)
  count
}
