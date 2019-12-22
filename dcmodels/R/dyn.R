source("dcmodels/R/dynUtils.R")
#' dyn model functions
dyn = list(
  LLVec = function(b, args){

    LL = 0
    
    U = list()
    PAllObs = rep(1,args$spec$nObs)

    for(t in 1:args$spec$nTime){
      for(l in 1:(args$spec$nLook+1))
        U[[l]] = matrix(args$X[[t]][[l]] %*% b, args$spec$nObs, args$spec$nAlt, byrow = T)
      
      fU = getfU(U, args$spec$transition)  #obtain instant utility plus downstream utility at  time t
      P = getRLProbas(fU, args$C[,t])
      #P = maxVec(P,args$stop[,t])  #if out_of_market (alter=1), let the probability = 1

      PAllObs = PAllObs * P
    }
    ##Yan add for loop and print(LL)
    for(i in 1:args$spec$nObs){
    LL = LL + log(PAllObs[i])
    }
    #print(LL)
    PAllObs
  },
  
  computeArgs = function(spec, D){
    X = list()
    
    for(t in 1:spec$nTime){
      X[[t]] = list()
      for(l in 1:(spec$nLook+1)){
        Dt = spec$modifyD(spec$D, spec$Global, t+l-1)   #get all attributes at t
        X[[t]][[l]] = create_X(spec$generic, spec$specific, Dt)
      }
    }
    
    # check if we reached stop point already
    stopMat = matrix(0, spec$nObs, spec$nTime)  
    for(t in 2:spec$nTime)    #Yan delete
    	stopMat[,t] = (spec$C[,t-1] %in% spec$stopAlt)  | (stopMat[,t-1] == 1)  
    
    list(spec = spec, X = X, C = spec$C, stop = stopMat)   
    #c(spec, M = M, X = X, stop = stopMat)
  },
  
  
  
  
  #Yan add
  computeArgs_Yan = function(spec, D){
    Z = spec$First
    X = list()
    
    for(t in 1:spec$nTime){
      X[[t]] = list()
      for(l in 1:(spec$nLook+1)){
        Dt = spec$modifyD(spec$D, spec$Global, t+l-1)   #get all attributes at t
        X[[t]][[l]] = create_X(spec$generic, spec$specific, Dt)
      }
      Z = updateZ(Z,spec$D[[t]],spec$C[,t],spec$varNames)   #update variables if make a purchase at time t
    }
    
    # check if we reached stop point already
    stopMat_Yan = matrix(0, spec$nObs, spec$nTime) 
    for(t in 2:spec$nTime){ 
      #Yan add: after t=13, if buy vehicle, out-of-market forever
      if(t > spec$nTime-4){
        stopMat_Yan[,t] = (spec$C[,t-1] %in% spec$stopAlt)  | (stopMat_Yan[,t-1] == 1)
        last_two_year_odd = stopMat_Yan[,t-1]*stopMat_Yan[,t-2]*stopMat_Yan[,t-3]*stopMat_Yan[,t-4]  
        for(i in 1:spec$nObs){
          if((last_two_year_odd[i]==1) & (t%%2 == 1)) 
            stopMat_Yan[i,t] = 0
        }
      }
      #Yan add: t is odd, (t-1) is even, buy vehicle at second time of a year, out-of-market for 4 time periods
      if((t < spec$nTime-3) & (t%%2 == 1)){
        a = (spec$C[,t-1] %in% spec$stopAlt)
        stopMat_Yan[,t] = stopMat_Yan[,t] + a
        stopMat_Yan[,t+1] = stopMat_Yan[,t+1] + a
        stopMat_Yan[,t+2] = stopMat_Yan[,t+2] + a
        stopMat_Yan[,t+3] = stopMat_Yan[,t+3] + a
      }
      #Yan add: t is even, (t-1) is odd, buy vehicle at first time of a year, out-of-market for 5 time periods
      if((t < spec$nTime-3) & (t%%2 == 0)){
        a = (spec$C[,t-1] %in% spec$stopAlt)
        stopMat_Yan[,t] = stopMat_Yan[,t] + a
        stopMat_Yan[,t+1] = stopMat_Yan[,t+1] + a
        stopMat_Yan[,t+2] = stopMat_Yan[,t+2] + a
        stopMat_Yan[,t+3] = stopMat_Yan[,t+3] + a
        stopMat_Yan[,t+4] = stopMat_Yan[,t+4] + a
      }      
    }
    
    
    #list(spec = spec, M = M, X = X, C = spec$C, stop = stopMat)   #Yan delete
    list(spec = spec, X = X, C = spec$C, stop = stopMat_Yan)    #Yan add
    #c(spec, M = M, X = X, stop = stopMat)
  },
  
  
  
  
  
  
  
  
  computeStart = function(spec, D){
    npar = length(spec$generic)
    for(s in spec$specific)
      npar = npar + length(s)
    rep(0, npar)
  },
  
  computeMisc = function(spec, D){
    n = c()
    comIndex = 1
    for(com in spec$generic){
      n = c(n,paste("common",comIndex,sep="_"))
      comIndex = comIndex + 1
    }
    for(s in spec$specific)
      n = c(n,s)
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
