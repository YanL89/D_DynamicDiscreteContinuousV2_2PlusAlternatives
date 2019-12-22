#' fits a model using maximum likelihood estimation
#' 
#' optimizes the LL function, finds standard deviations and reports this 
#' along with the value of the LL at the max
#' 
#' @param modelFns functions that specify the model
#' @param args arguments for the functions
#' @param start starting value for optimization
#' @param misc other stuff for optim
#' @param spec user specifications
#' @return list of results (param name, estimate, sd, t, p) and value of max LL
#' @export
modelFit = function(modelFns, args, start, misc, spec){
  # check for bounds
  big_num = 9999999
  lb = getLb(misc, start)
  ub = getUb(misc, start)
  
  #f = modelFns$LLVec
  delta = rep(spec$delta, length(start))

  O = optim(fn = LLWrapper, args = args, f = modelFns$LLVec
  	, method = "BFGS", hessian = F, par = start
      , control = list(fnscale = -1, ndeps = delta, abstol = spec$abstol, reltol = spec$reltol))
 
  betaHat = O$par
  SD = getSD(spec, betaHat, modelFns, args)
  t = betaHat / SD
  p = 2*pnorm(-abs(t))
  results = data.frame(name = misc[["names"]], beta_hat = betaHat, SD = SD, t = t, p = p)
  maxLL = O$value
  list(results = results, maxLL = maxLL)				
}

#' wrapper for \code{\link{modelFit}}
#' 
#' takes more high level user specifications and model functions to compute
#' the maximum log-likelihood estimates for the model
#' 
#' @param modelFns a list of model functions
#' @param spec list of user specification of the model and other params
#' @param D the data set to look for the data
#' @return a table of estimated coefficients and associated statistic + value of max LL
#' @export

model = function(modelFns, spec, D){
  D$const = 1
  
  # complete specification
  spec = addVal(spec,"delta", 0.0001)
  spec = addVal(spec,"reltol", 0.0001)
  spec = addVal(spec,"abstol",0.0001)
  spec = addVal(spec,"verbose",TRUE)
  spec = addVal(spec,"SD","hessian")
  
  # if users specifies starting values we take them, hence the special function
  start = getStart(modelFns, spec, D)  
  args = modelFns$computeArgs(spec, D)  
  args = c(args, spec)
  misc = modelFns$computeMisc(spec, D)   #add names of attributes
  
  M = modelFit(modelFns, args, start, misc, spec)
  
  class(M) = "dcModel"
  return(M)
}


model_Yan = function(modelFns, spec, D){
	D$const = 1

  # complete specification
	spec = addVal(spec,"delta", 0.01)
	spec = addVal(spec,"reltol", 0.0001)
	spec = addVal(spec,"abstol",0.0001)
	spec = addVal(spec,"verbose",TRUE)
	spec = addVal(spec,"SD","hessian")
	
  # if users specifies starting values we take them, hence the special function
  start = getStart(modelFns, spec, D)  
	args = modelFns$computeArgs_Yan(spec, D)
	args = c(args, spec)
  
  misc = modelFns$computeMisc(spec, D)   #add names of attributes
  
  M = modelFit(modelFns, args, start, misc, spec)

  class(M) = "dcModel"
  return(M)
}

#' wrapper for \code{\link{model}}
#' 
#' takes more high level user specifications and model functions to compute
#' the maximum log-likelihood estimates for the model
#' 
#' @param modelFns a list of model functions
#' @param spec list of user specification of the model and other params
#' @param D the data set to look for the data
#' @param fracKeep fraction of observations to keep in estimation
#' @return a table of estimated coefficients and associated statistic + value of max LL and model applied to validation sample
#' @export
modelAndApply = function(modelFns, spec, D, fracKeep){
  IdKeep = sample(nrow(D), nrow(D) * fracKeep)
  DKeep = D[IdKeep,]
  DVal = D[-IdKeep,]
  Val = model(modelFns, spec, DKeep)
  Apply = modelFns$apply(spec, Val, DVal)
  list(Val = Val, Apply = Apply)
}

#' calculate model elasticity with respect to specified variables
#' 
#' multiplie specified variables by 1% and applies the model to the
#' new data set
#' 
#' @param modelFns a list of model functions
#' @param spec list of user specification of the model and other params
#' @param D the data set to look for the data
#' @param varsElas variables to check elasticity
#' @param modelCx estimated model coefficients
#' @return tables of predicted value to compute elasticity
#' @export
modelElasticity = function(modelFns, spec, D, varsElas, modelCx = NULL){
  Val = modelCx
  if(is.null(Val))
    Val = model(modelFns, spec, D)
  Apply = modelFns$apply(Val, spec, D)
  cat("expected values for original data and model:\n")
  print(Apply$expected)
  cat("\n")
  for(v in varsElas){
    D2 = D
    D2[,v] = D[,v] * 1.01
    A2 = modelFns$apply(spec, Val, D2)
    cat("expected values for ",v," multiplied by 1.01 (101%):\n")
    print(A2$expected)
    cat("\n")    
  }
}
