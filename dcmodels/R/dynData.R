#' reads and veryfies data files
#' 
#' @param spec list of model specifications
#' @export
#' @return the original list plus data sets in it
checkFillDynSpec = function(spec){
  #  spec = pratt
  if("checked" %in% names(spec) && spec$checked)
    return(spec)

  spec$nTimeTotal = spec$nTime + spec$nLook
  spec$nAlt = length(spec$specific)  
  spec$C = getChoices(spec)
  spec$nObs = nrow(spec$C)
  spec$M = getMiles(spec)

  if( ! is.null(spec$Time)) {
  	spec$Time = read.table(spec$Time, header = TRUE)
  	} else {
  	spec$Time = data.frame(matrix(0, spec$nObs,0))
  	}
  checkTime(spec)
  #spec$timeNames = getTimeNames(spec)
  
  if( ! is.null(spec$Global)) {
  	spec$Global = getGlobal(spec)
  	} else {
  	spec$Global = data.frame(matrix(0, spec$nObs,0))
  	}

  spec$D = fillDynD(spec)
  spec$varNames = getVarNames(spec)  
  #spec$First = getFirst(spec)

  # needs to be after fillDynD
  checkCommonSpecific(spec)

  spec$checked = TRUE
  spec
}

#' get names of variables in Time
#' 
#' @param spec list of model specifications
#' @export
#' @return variables in Time
getTimeNames = function(spec){
	names = unlist(strsplit(names(spec$Time), ".", fixed = TRUE))
	nvar = length(names) / spec$nTimeTotal / 2
	names[seq(from = 1, by = 2*spec$nTimeTotal, length = nvar)]
}

#' verify if the time data set is consistent
#' 
#' @param spec list of model specifications
#' @export
checkTime = function(spec){
	if(0 == ncol(spec$Time))
		return()
  names = unlist(strsplit(names(spec$Time), ".", fixed = TRUE))
  nvar = length(names) / spec$nTimeTotal / 2
  for(i in 1:nvar){
    index = (i-1) * (spec$nTimeTotal * 2) + 1
    n = names[seq(from = index, length = spec$nTimeTotal, by = 2)]
    if(mean(n == n[1]) < 1){
      cat("Time data problem. the names of the variable # ",i
          , "in you time data set is not consistent. if you have"
          , "v1 and v2 with 3 times period you need to have the following"
          , "variables in the same order: v1.1, v1.2, v1.3, v2.1, v2.2"
          , "v2.3")
      stop("time data problem")
    }
    ids = strtoi(names[seq(from = index + 1
                           , length = spec$nTimeTotal, by = 2)])
    if(mean(ids == (1:spec$nTimeTotal)) < 1) {
      cat("time data problem. the time index of variable # ", i
          , "in you time data set is not ordered. if the variable is"
          , "v1 with 3 time perdiods, we expect v1.1, v1.2 and v1.3")
      stop("time data problem")
    }
  }
  
}

#' counts the number of time periods
#' 
#' @param spec list of model specifications
#' @export
#' @return the number of time periods
getNTime = function(Time){
  names = strsplit(names(Time), ".", fixed = TRUE)  
  ntime = 0
  for(i in 1:length(names)){
    if(2 != length(names[[i]])){
      cat("there is a problem with the ",i,"-th variable"
          ," in your data set. format should be <name>.<time>"
          ," for example \"somename.12\"")
      stop("illegal name in Time data frame")
    }
    newval = strtoi(names[[i]][2])
    if(newval > ntime)
      ntime = newval
  }
  ntime	
}

#' reads and checks consistency of attributes at each time peried
#' 
#' @param spec list of model specifications
#' @export
#' @return the number of time periods
fillDynD = function(spec){
  D = list()
  for(i in 1:spec$nTimeTotal){
    fname = paste(spec$D,i,sep="")
    fname = paste(fname,"txt",sep=".")
    # check if i-th data file exists
    if( ! file.exists(fname)){
      cat("data file # ", i, " does not exist\n")
      cat("I am expecting a file named", fname,"\n")
      stop("data file problem")
    }
    # read file and make sure the names are the same
    D[[i]] = read.table(fname, header = TRUE)
    if(mean(names(D[[i]]) == names(D[[1]])) < 1){
      cat("data file # ", i, " does not contain the same varibles"
          , "names than the others")
      stop("data file variable problem")
    }
    # make sure file has correct number of entries
    if(nrow(D[[i]]) != spec$nObs){
      cat("number of observations in ", i, "-th data set is not "
          , "consistent with the rest")
      stop("data size problem")
    }
    # make sure variable names follow required specifications
    nvar = spec$nvar_altspec
    for(k in 1:nvar){
      names = unlist(strsplit(names(D[[i]]), ".", fixed = TRUE))
      index = seq(from = 2*spec$nAlt*(k-1) + 1, length = spec$nAlt
                  , by = 2)
      n = names[index]
      id = strtoi(names[index+1])
      # check name
      if(mean(n == n[1]) < 1){
        cat("in data set ", i, ", varialbe # ", k,":\n")
        cat("the name is not consistent. names should be "
            , "var.1, var.2, ..., var.nalt\n")
        stop("var name problem")	
      }
      # check id
      if(mean(id == 0:(spec$nAlt-1)) < 1){
        cat("in data set ", i, ", varialbe # ", k,":\n")
        cat("the id part is not consistent. names should be "
            , "var.1, var.2, ..., var.nalt\n")
        stop("var name problem")	
      }
    }
  }
  D
}

#' get variable names for alternatives attributes
#'
#' @param spec list of model specifications
#' @export
#' @return a vector of all possible attributes for alternatives
getVarNames = function(spec){
  varNames = unlist(strsplit(names(spec$D[[1]]), ".", fixed = TRUE))
  varNames = varNames[seq(from = 1, length = length(varNames)/2/spec$nAlt
                          , by = 2*spec$nAlt)]
  varNames
}	

#' reads global variables
#'
#' @param spec list of model specifications
#' @export
#' @return a data frame of all global variables
getGlobal = function(spec){
  Global = read.table(spec$Global, header = TRUE)
  if(nrow(Global) != spec$nObs){
    stop("incorrect number of observation in Global data set")
  }
  Global$one = rep(1, spec$nObs)
  Global$zero = rep(0, spec$nObs)
  Global
}

#' reads attributes af alternative held at the beginning of the process
#'
#' @param spec list of model specifications
#' @export
#' @return a data frame of original items characteristics (Z at first time period)
getFirst = function(spec){
  First = read.table(spec$First, header = TRUE)
  if(nrow(First) != spec$nObs){
    stop("incorrect number of observation if first data file")
  }
  if(mean(colnames(First) %in% spec$varNames) < 1){
    cat("variables in First must be the same than Time and "
        , "other alternatives data")
    stop("problem with first data")
  }
  if(nrow(First) != spec$nObs){
    stop("incorrect number of observation in the first data set")
  }
  First
}

#' checks if variables specified in common and specific are available
#' 
#' @param spec list of model specifications
#' @export
checkCommonSpecific = function(spec){
  if(is.null(spec$modifyD)){
    car("you must specify a way to create the data set for computations"
        , "look up spec.R for an example")
    stop("must specify modifyD function")
  }
  for(gen in spec$generic)
    if(length(gen) != spec$nAlt){
      cat("problem with the specification of your common parameters"
          , " the length of each elements in common must be the"
          , " number of alternatives\n")
      stop("problem with common parameters")
    }
  
  D = spec$modifyD(spec$D, spec$Global, 1)
  
  names = c(unlist(spec$generic), unlist(spec$specific))
  names %in% colnames(D)
  if(mean(names %in% colnames(D)) < 1){
    cat("your modifyD function does not specify all variables in "
        , "your model specification (specific and common)\n")
    stop("modifyD / variable problem")
  }
}

#' reads and checks consistency af choices
#' 
#' @param spec list of model specifications
#' @export
#' @return a matrix af all choices
getChoices = function(spec){
  C = NULL
  if(is.character(spec$choices)){
    C = read.table(spec$choices, header = TRUE)
  }
  if(is.matrix(spec$choices)){
    C = spec$choices
  }
  
  #if(nrow(C) != spec$nObs){
  #  cat("incorrect number of observatiosn in choices\n")
  #  stop("choices problem")
  #}
  if(min(C) != 0 || max(C) != spec$nAlt - 1){
    print("problem with choices encoding. use choice = 0 for")
    print("no change and choice = 1,....nAlt for choice of")
    print("an alternative")
    stop("choice encoding problem")
  }
  C
}

# return a matrix af all miles
getMiles = function(spec){
  M = NULL
  if(is.character(spec$miles)){
    M = read.table(spec$miles, header = TRUE)
  }
  if(is.matrix(spec$miles)){
    M = spec$miles
  }
  M
}

