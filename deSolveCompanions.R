## Stuff from
#  https://rdrr.io/cran/deSolve/src/R/functions.R  
#  https://rdrr.io/cran/deSolve/src/R/checkevents.R 
#  https://rdrr.io/cran/deSolve/src/R/forcings.R 

## CONTENT OF https://rdrr.io/cran/deSolve/src/R/checkevents.R BELOW

### ============================================================================
### Check events data set
### Changes version 1.11: event can be an R-function, even if DLL model
###                       continueeroot: to continue even if a root is found
### ============================================================================

checkevents <- function (events, times, vars, dllname, root = FALSE) {
  
  if (is.null(events)) return(list())
  if (is.null(events$data) && is.null(events$func) &&
      is.null(events$terminalroot)) return(list())
  
  funevent <- events$func
  
  if (root) {  # check if root should trigger an event...
    Root <- events$root
    if (is.null(Root)) Root <- 0
    Root <- as.integer(Root)
  } else Root <- 0L
  
  maxroot <- events$maxroot
  if (is.null(maxroot)) maxroot <- 100  # number of roots to save.
  if (maxroot < 0)
    stop("events$maxroot should be > 0 in events")
  Terminalroot <- events$terminalroot
  
  if (! is.null(Terminalroot) && is.null(funevent))
    funevent <- function(t,y,p) return(y)  # dummy event function
  
  if (is.null(Terminalroot))
    Terminalroot <- 0  # at which roots simulation should continue
  
  
  ## ----------------------
  ## event in a function
  ## ----------------------
  if (!is.null(funevent)) {
    if (inherits (funevent, "CFunc")) {
      funevent <- body(funevent)[[2]]
      Type <- 3
    } else if (is.character(funevent)){
      if (is.null(dllname))
        stop("'dllname' should be given if 'events$func' is a string")
      if (is.loaded(funevent, PACKAGE = dllname, type = "") ||
          is.loaded(funevent, PACKAGE = dllname, type = "Fortran")) {
        funevent <- getNativeSymbolInfo(funevent, PACKAGE = dllname)$address
      } else
        stop(paste("'events$func' should be loaded ",funevent))
      Type <- 3
    } else {
      Type <- 2  # SHOULD ALSO CHECK THE FUNCTION if R-function....
      #      if (!is.null(dllname))      KARLINE: removed that 02/07/2011
      #       stop("'events$func' should be a string, events specified in compiled code if 'dllname' is not NULL")
    }
    if (Root == 0) {
      if (is.null(events$time))
        stop("either 'events$time' should be given and contain the times of the events, if 'events$func' is specified and no root function or your solver does not support root functions")
      eventtime <- sort(as.double(events$time)) # Karline: sorted that 4-01-2016
      
      if (any(!(eventtime %in% times))) {
        warning("Not all event times 'events$time' are in output 'times' so they are automatically included.")
        uniqueTimes <- cleanEventTimes(times, eventtime)
        if (length(uniqueTimes) < length(times))
          warning("Some time steps were very close to events - only the event times are used in these cases.")
        times <- sort(c(uniqueTimes, eventtime))
      }
    } else eventtime <- min(times) - 1  # never reached....
    return (list (Time = eventtime, SVar = NULL, Value = NULL,
                  Method = NULL, Type = as.integer(Type), func = funevent,
                  Rootsave = as.integer(maxroot), Root = Root,
                  Terminalroot = as.integer(Terminalroot), newTimes = times))    # added newTimes - Karline 4-01-2016
    
  }
  ## ----------------------
  ## event as a data series
  ## ----------------------
  eventdata <- events$data
  if (is.matrix(eventdata)) eventdata <- as.data.frame(eventdata)
  
  if (ncol(eventdata) < 3)
    stop("'event' should have at least 3 columns: state variable, time, value")
  
  if (!is.data.frame(eventdata))
    stop("'event' should be a data.frame with 3(4) columns: state variable, time, value, (method)")
  
  ## this should make check < 3 columns obsolete
  evtcols <-  c("var", "time", "value", "method")
  if (!all(evtcols %in% names(eventdata)))
    stop("structure of events does not match specification, see help('events')")
  
  ## make sure that event data frame has correct order
  eventdata <- eventdata[evtcols]
  
  ## variables, 1st column should be present
  if (is.factor(eventdata[,1]))
    eventdata[,1] <- as.character(eventdata[,1])
  
  if (is.character(eventdata[,1]))  {
    vv <- match(eventdata[,1], vars)
    if (is.character(eventdata[,1]))  {
      vv <- match(eventdata[,1],vars)
      if (any(is.na(vv)))
        stop("unknown state variable in 'event': ", paste(eventdata[,1][which(is.na(vv))], ","))
      eventdata[,1] <- vv
    } else if (max(eventdata[,1]) > length(vars))
      stop("unknown state variable in 'event': ", paste(eventdata[,1][which(is.na(vv))],","))
    eventdata[,1] <- vv
  } else if (max(eventdata[,1])>length(vars))
    stop("too many state variables in 'event'; should be < ", paste(length(vars)))
  
  ## 2nd and 3rd columns should be numeric
  if (!is.numeric(eventdata[,2]))
    stop("times in 'event', 2nd column should be numeric")
  
  if (!is.numeric(eventdata[,3]))
    stop("values in 'event', 3rd column should be numeric")
  
  ## Times in 'event' should be embraced by 'times'
  rt <- range(times)
  ii <- c(which(eventdata[,2] < rt[1]), which(eventdata[,2] > rt[2]))
  if (length(ii) > 0)
    eventdata <- eventdata [-ii,]
  if (any(!(eventdata[,2] %in% times))) {
    warning("Not all event times 'events$times' were in output 'times' so they are automatically included.")
    uniqueTimes <- cleanEventTimes(times, eventdata[,2])
    if (length(uniqueTimes) < length(times))
      warning("Some time steps were very close to events - only the event times are used in these cases.")
    times <- sort(c(uniqueTimes, eventdata[,2]))
  }
  
  if (any(!(eventdata[,2] %in% times))) {
    warning("Not all event times 'events$times' where in output 'times' so they are automatically included.")
    uniqueTimes <- cleanEventTimes(times, eventdata[,2])
    if (length(uniqueTimes) < length(times))
      warning("Some time steps were very close to events - only the event times are used in these cases.")
    times <- sort(c(uniqueTimes, eventdata[,2]))
  }
  ## check if times are ordered and if not, fix it
  if (any(diff(eventdata[,2]) < 0)) {
    warning("Time of events ('time' column of 'events') was not ordered.")
    ord <- order(eventdata[,2])
    eventdata <- eventdata[ord,]
  }
  
  ## 4th column: method; if not available: "replace" = method 1 - to date: 3 methods
  if (ncol(eventdata) ==3)
    eventdata$method <- rep(1,nrow(eventdata))
  else if (is.numeric(eventdata[,4])) {
    if (max(eventdata[,4]) > 3 | min(eventdata[,4]) < 1)
      stop("unknown method in 'event': should be >0 and < 4")
  } else {
    vv <- charmatch(eventdata[,4],c("replace","add","multiply"))
    if (any(is.na(vv)))
      stop("unknown method in 'event': ", paste(eventdata[,3][which(is.na(vv))],","),
           " should be one of 'replace', 'add', 'multiply'")
    eventdata$method <- vv
  }
  
  ## Check the other events elements (see optim code)
  con <- list(ties = "notordered", time = NULL, data = NULL, func = NULL, root = NULL)
  nmsC <- names(con)
  con[(namc <- names(events))] <- events
  if (length(noNms <- namc[!namc %in% nmsC]) > 0)
    warning("unknown names in events: ", paste(noNms, collapse = ", "))
  
  ## Check what needs to be done in case the time series is not "ordered"
  
  if (!identical(con$ties, "ordered")) { # see approx code
    
    ## first order with respect to time (2nd col), then to variable (1st col)
    if(length(x <- unique(eventdata[,1:2])) < nrow(eventdata)){
      ties <- mean
      if (missing(ties))
        warning("collapsing to unique 'x' values")
      eventdata <- aggregate(eventdata[,c(3, 4)], eventdata[,c(1, 2)], ties)
      ties <- mean
      if (missing(ties))
        warning("collapsing to unique 'x' values")
      eventdata <- aggregate(eventdata[,c(3,4)], eventdata[,c(1,2)], ties)
    }
  }
  
  return (list (Time = as.double(eventdata[,2]), SVar = as.integer(eventdata[,1]),
                Value = as.double(eventdata[,3]), Method = as.integer(eventdata[,4]),
                Rootsave = as.integer(maxroot),
                Type = 1L, Root = Root,
                Terminalroot = as.integer(Terminalroot),
                newTimes = times))
}

## CONTENT OF https://rdrr.io/cran/deSolve/src/R/functions.R  BELOW


## ========================================================================
## General functions of deSolve
## ========================================================================

timestep <- function (prev = TRUE) {
  out <- .Call("getTimestep", PACKAGE = "deSolve")
  if (prev)
    return(out[1])
  else
    return(out[2])
}


## ========================================================================
## Check solver input - livermore solvers and rk
## ========================================================================

checkInput <- function(y, times, func, rtol, atol,
  jacfunc, tcrit, hmin, hmax, hini, dllname, jacname = "jacfunc")
{
  if (!is.numeric(y))     stop("`y' must be numeric")
  n <- length(y)
  if (! is.null(times) && !is.numeric(times))
    stop("`times' must be NULL or numeric")
  if (!is.function(func) && !is.character(func))
    stop("`func' must be a function or character vector")
  if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
    stop("specify the name of the dll or shared library where func can be found (without extension)")
  if (!is.numeric(rtol))  stop("`rtol' must be numeric")
  if (!is.numeric(atol))  stop("`atol' must be numeric")
  if (!is.null(tcrit) & !is.numeric(tcrit)) stop("`tcrit' must be numeric")
  if (!is.null(jacfunc) && !(is.function(jacfunc) || is.character(jacfunc)))
    stop(paste(jacname," must be a function or character vector"))
  if (length(atol) > 1 && length(atol) != n)
    stop("`atol' must either be a scalar, or as long as `y'")
  if (length(rtol) > 1 && length(rtol) != n)
    stop("`rtol' must either be a scalar, or as long as `y'")
  if (!is.numeric(hmin))   stop("`hmin' must be numeric")
  if (hmin < 0)            stop("`hmin' must be a non-negative value")
  if (is.null(hmax))
    hmax <- if (is.null(times)) 0 else max(abs(diff(times)))
  if (!is.numeric(hmax))   stop("`hmax' must be numeric")
  if (hmax < 0)            stop("`hmax' must be a non-negative value")
  if (hmax == Inf)  hmax <- 0
  if (!is.null(hini))
   if(hini < 0)            stop("`hini' must be a non-negative value")
  return(hmax)
}

## ========================================================================
## Check solver input - euler and rk4
## ========================================================================

checkInputEuler <- function (y, times, func, dllname) {
    if (!is.numeric(y))  stop("`y' must be numeric")
    n <- length(y)
    if (! is.null(times) && !is.numeric(times))
        stop("`times' must be NULL or numeric")
    if (!is.function(func) && !is.character(func))
      stop("`func' must be a function or character vector")
    if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
      stop("You need to specify the name of the dll or shared library where func can be found (without extension)")
}

## ========================================================================
## Check ode function call - livermore solvers
## ========================================================================

checkFunc<- function (Func2, times, y, rho) {
    ## Call func once to figure out whether and how many "global"
    ## results it wants to return and some other safety checks
    tmp <- eval(Func2(times[1], y), rho)
    if (!is.list(tmp))
      stop("Model function must return a list\n")
    if (length(tmp[[1]]) != length(y))
      stop(paste("The number of derivatives returned by func() (",
                 length(tmp[[1]]),
                 ") must equal the length of the initial conditions vector (",
                 length(y), ")", sep = ""))

    ## use "unlist" here because some output variables are vectors/arrays
    Nglobal <- if (length(tmp) > 1)
      length(unlist(tmp[-1]))  else 0
    ## Karline: changed this:
    ##  Nmtot is now a list with names, dimensions,... for 1-D, 2-D vars
    Nmtot <- list()
    Nmtot$colnames <- attr(unlist(tmp[-1]), "names")

    Nmtot$lengthvar <- unlist(lapply(tmp, length))
    if (length(Nmtot$lengthvar) < Nglobal+1){
      Nmtot$dimvar <- lapply(tmp[-1], dim)
    }
   return(list(Nglobal = Nglobal, Nmtot = Nmtot))
}

## ========================================================================
## Check event function calls
## ========================================================================

checkEventFunc<- function (Func, times, y, rho) {
    ## Call func once
    tmp <- eval(Func(times[1], y), rho)
    if (length(tmp) != length(y))
      stop(paste("The number of values returned by events$func() (",
                 length(tmp),
                 ") must equal the length of the initial conditions vector (",
                 length(y), ")", sep = ""))
    if (!is.vector(tmp))
      stop("The event function 'events$func' must return a vector\n")
}

## ========================================================================
## Check ode function call - euler and rk solvers
## ========================================================================

checkFuncEuler<- function (Func, times, y, parms, rho, Nstates) {
      ## Call func once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
      tmp <- eval(Func(times[1], y, parms), rho)
      if (!is.list(tmp)) stop("Model function must return a list\n")
      if (length(tmp[[1]]) != Nstates)
        stop(paste("The number of derivatives returned by func() (",
                   length(tmp[[1]]),
                   "must equal the length of the initial conditions vector (",
                   Nstates, ")", sep=""))

      ## use "unlist" because output variables can be vectors/arrays
      Nglobal <- if (length(tmp) > 1)
        length(unlist(tmp[-1])) else 0
      Nmtot <- list()
      Nmtot$colnames <- attr(unlist(tmp[-1]), "names")
      Nmtot$lengthvar <- unlist(lapply(tmp, length))
      if (length(Nmtot$lengthvar) < Nglobal+1){
        Nmtot$dimvar <- lapply(tmp[-1], dim)
      }
   return(list(Nglobal = Nglobal, Nmtot = Nmtot))

}

## ========================================================================
## check ode DLL input
## ========================================================================

checkDLL <- function (func, jacfunc, dllname,
                      initfunc, verbose, nout, outnames, JT = 1) {

    if (sum(duplicated (c(func, initfunc, jacfunc))) > 0)
      stop("func, initfunc, or jacfunc cannot be the same")
    ModelInit <- NA
    if (! is.null(initfunc))  # to allow absence of initfunc
      if (inherits (initfunc, "CFunc"))
        ModelInit <- body(initfunc)[[2]]
      else if (is.loaded(initfunc, PACKAGE = dllname, type = "") ||
        is.loaded(initfunc, PACKAGE = dllname, type = "Fortran"))  {
        ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
      } else if (initfunc != dllname && ! is.null(initfunc))
        stop(paste("'initfunc' not loaded ", initfunc))

    ## Easier to deal with NA in C-code
    if (is.null(initfunc)) ModelInit <- NA

    ## copy value of func to funcname
    ## check to make sure it describes a function in a loaded dll

    funcname <- func
    ## get the pointer and put it in func

    if (inherits (func, "CFunc"))
        Func <- body(func)[[2]]
    else if(is.loaded(funcname, PACKAGE = dllname)) {
      Func <- getNativeSymbolInfo(funcname, PACKAGE = dllname)$address
    } else stop(paste("dyn function 'func' not loaded", funcname))

    ## Finally, is there a Jacobian?
    if (!is.null(jacfunc)) {
      if (!is.character(jacfunc))
        switch (JT,
          stop("If 'func' is dynloaded, so must 'jacfunc' be"),
          stop("If 'func' is dynloaded, so must 'jacvec' be")
        )
      jacfuncname <- jacfunc
      if (inherits (jacfunc, "CFunc"))
        JacFunc <- body(jacfunc)[[2]]

      else if(is.loaded(jacfuncname, PACKAGE = dllname))  {
        JacFunc <- getNativeSymbolInfo(jacfuncname, PACKAGE = dllname)$address
      } else stop(paste("cannot integrate: jac function not loaded ", jacfunc))
    } else JacFunc <- NULL
    Nglobal <- nout
    Nmtot <- list()
    if (is.null(outnames))
      { Nmtot$colnames   <- NULL} else
    if (length(outnames) == nout)
      { Nmtot$colnames   <- outnames} else
    if (length(outnames) > nout)
      Nmtot$colnames <- outnames[1:nout] else
    Nmtot$colnames <- c(outnames,(length(outnames)+1):nout)

    cnames <- outnames
    unames <- unique(outnames)
    if (length(cnames) > length(unames))
      Nmtot$lengthvar <- c(NA,
        sapply (unames, FUN = function(x) length(which(cnames == x))))

   return(list(ModelInit = ModelInit, Func = Func, JacFunc = JacFunc,
     Nglobal = Nglobal, Nmtot = Nmtot))
}

## =============================================================================
## print integration task
## =============================================================================
printtask <- function(itask, func, jacfunc) {
    printM("\n--------------------")
    printM("Time settings")
    printM("--------------------\n")
    if (itask==1) printM("  Normal computation of output values of y(t) at t = TOUT") else
    if (itask==2) printM("  Take one step only and return.")                          else
    if (itask==3) printM("  istop at the first internal mesh point at or beyond t = TOUT and return. ")  else
    if (itask==4) printM("  Normal computation of output values of y(t) at t = TOUT but without overshooting t = TCRIT.") else
    if (itask==5) printM("  Take one step, without passing TCRIT, and return.")
    printM("\n--------------------")
    printM("Integration settings")
    printM("--------------------\n")
    if (is.character(func)) printM(paste("  Model function a DLL: ", func)) else
    printM("  Model function an R-function: ")
    if (is.character(jacfunc)) printM(paste ("  Jacobian specified as a DLL: ", jacfunc)) else
    if (!is.null(jacfunc))     printM("  Jacobian specified as an R-function: ") else
    printM("  Jacobian not specified")
    cat("\n")
}

## =============================================================================
## Make Istate vector similar for all solvers.
## =============================================================================

setIstate <- function(istate, iin, iout)
{
  IstateOut <- rep(NA, 21)
  IstateOut[iout] <- istate[iin]
  IstateOut
}


## =============================================================================
## Output cleanup  - for the Livermore solvers
## =============================================================================

saveOut <- function (out, y, n, Nglobal, Nmtot, func, Func2,
  iin, iout, nr = 4) {
  troot  <- attr(out, "troot")
  istate <- attr(out, "istate")
  istate <- setIstate(istate,iin,iout)

  valroot  <- attr(out, "valroot")
  indroot <- attr(out, "indroot")

  Rstate <- attr(out, "rstate")
  rstate <- rep(NA,5)
  rstate[1:nr] <- Rstate[1:nr]

  nm <- c("time",
          if (!is.null(attr(y, "names"))) names(y) else as.character(1:n))
  if (Nglobal > 0) {
    nm <- c(nm,
            if (!is.null(Nmtot$colnames))
              Nmtot$colnames else as.character((n+1) : (n + Nglobal)))
  }
  attr(out,"istate") <- istate
  attr(out, "rstate") <- rstate
  if (! is.null(Nmtot$lengthvar))
    if (is.na(Nmtot$lengthvar[1]))Nmtot$lengthvar[1] <- length(y)
  attr(out, "lengthvar") <- Nmtot$lengthvar
  if (! is.null(troot)) attr(out, "troot") <-  troot
  if (! is.null(valroot)) attr(out, "valroot") <- matrix(nrow = n, valroot)
  if (! is.null(indroot)) attr(out, "indroot") <- indroot

  ii <- if (is.null(Nmtot$dimvar))
    NULL else !(unlist(lapply(Nmtot$dimvar, is.null))) # variables with dimension
  if (sum(ii) >0)
    attr(out, "dimvar") <- Nmtot$dimvar[ii]     # dimensions that are not null
  class(out) <- c("deSolve", "matrix")          # a differential equation
  dimnames(out) <- list(nm, NULL)
  return (t(out))
}

## =============================================================================
## Output cleanup  - for the Runge-Kutta solvers
## =============================================================================

saveOutrk <- function(out, y, n, Nglobal, Nmtot, iin, iout, transpose = FALSE)  {
  ## Names for the outputs
  nm <- c("time",
    if (!is.null(attr(y, "names"))) names(y) else as.character(1:n)
   )
  ## Global outputs
  if (Nglobal > 0) {
    nm  <- c(nm,
      if (!is.null(Nmtot$colnames))
        Nmtot$colnames else as.character((n + 1) : (n + Nglobal))
    )
  }
  ## Column names and state information
  dimnames(out) <- list(NULL, nm)
  istate <- attr(out, "istate")
  istate <- setIstate(istate, iin, iout)
  attr(out,"istate") <- istate
  if (! is.null(Nmtot$lengthvar))
    if (is.na(Nmtot$lengthvar[1])) Nmtot$lengthvar[1] <- length(y)
  attr(out, "lengthvar") <- Nmtot$lengthvar

  ii <- if (is.null(Nmtot$dimvar))
    NULL else !(unlist(lapply(Nmtot$dimvar, is.null))) # variables with dimension
  if (sum(ii) >0)
    attr(out, "dimvar") <- Nmtot$dimvar[ii] # only those which are not null
  class(out) <- c("deSolve", "matrix")      # output of a differential equation
  if (transpose)
    return(t(out))
  else
    return(out)
}



## CONTENT OF https://rdrr.io/cran/deSolve/src/R/forcings.R  BELOW

### ============================================================================
### Check forcing function data set, event inputs and time-lag input
### ============================================================================


checkforcings <- function (forcings, times, dllname,
                           initforc, verbose, fcontrol = list()) {
  
  
  ## Check the names of the initialiser function
  
  if (is.null(initforc))
    stop(paste("initforc should be loaded if there are forcing functions ",initforc))
  
  if (inherits (initforc, "CFunc")) {
    ModelForc <- body(initforc)[[2]]
  }  else if (is.loaded(initforc, PACKAGE = dllname, type = "") ||
              is.loaded(initforc, PACKAGE = dllname, type = "Fortran")) {
    ModelForc <- getNativeSymbolInfo(initforc, PACKAGE = dllname)$address
  } else
    stop(paste("initforc should be loaded if there are forcing functions ",initforc))
  
  ## Check the type of the forcing function data series
  
  if (is.data.frame(forcings)) forcings <- list(a=forcings)
  if (! is.list(forcings)) forcings <- list(a=forcings)
  nf <- length(forcings)
  #1 check if each forcing function consists of a 2-columned matrix
  for (i in 1:nf) {
    if (ncol(forcings[[i]]) != 2)
      stop("forcing function data sets should consist of two-colum matrix")
  }
  
  ## Check the control elements (see optim code)
  
  con <- list(method="linear", rule = 2, f = 0, ties = "ordered")
  nmsC <- names(con)
  con[(namc <- names(fcontrol))] <- fcontrol
  if (length(noNms <- namc[!namc %in% nmsC]) > 0)
    warning("unknown names in fcontrol: ", paste(noNms, collapse = ", "))
  
  method <- pmatch(con$method, c("linear", "constant"))
  if (is.na(method))
    stop("invalid interpolation method for forcing functions")
  # 1 if linear, 2 if constant...
  
  ## Check the timespan of the forcing function data series
  
  # time span of forcing function data sets should embrace simulation time...
  # although extrapolation is allowed if con$rule = 2 (the default)
  r_t <- range(times)
  
  for (i in 1:nf) {
    r_f <- range(forcings[[i]][,1])   # time range of this forcing function
    
    if (r_f[1] > r_t[1]) {
      if (con$rule == 2) {
        mint <- c(r_t[1],forcings[[i]][1,2] )
        forcings[[i]] <- rbind(mint,forcings[[i]])
        if(verbose)
          warning(paste("extrapolating forcing function data sets to first timepoint",i))
      } else
        stop(paste("extrapolating forcing function data sets to first timepoint",i))
    }
    
    nr   <- nrow(forcings[[i]])
    if (r_f[2] < r_t[2]) {
      if (con$rule == 2) {
        maxt <- c(r_t[2],forcings[[i]][nr,2] )
        forcings[[i]] <- rbind(forcings[[i]],maxt)
        if(verbose)
          warning(paste("extrapolating forcing function data sets to last timepoint",i))
      } else
        stop(paste("extrapolating forcing function data sets to last timepoint",i))
    }
    
  }
  
  ## Check what needs to be done in case the time series is not "ordered"
  
  if (!identical(con$ties, "ordered")) { # see approx code
    
    for (i in 1:nf) {
      
      x <- forcings[[i]][,1]
      nx <- length(x)
      if (length(ux <- unique(x)) < nx) {  # there are non-unique values
        y <- forcings[[i]][,2]
        ties <- con$tiesn
        if (missing(ties))
          warning("collapsing to unique 'x' values")
        y <- as.vector(tapply(y, x, ties))
        x <- sort(ux)
        forcings[[i]] <- cbind(x, y)
        
      } else {                             # values are unique, but need sorting
        y <- forcings[[i]][,2]
        o <- order(x)
        x <- x[o]
        y <- y[o]
        forcings[[i]] <-  cbind(x,y)
      }
    } # i
  }
  
  ## In case the interpolation is of type "constant" and f not equal to 0
  ## convert y-series, so that always the left value is taken
  if (method == 2 & con$f != 0) {
    for (i in 1:nf) {
      y <- forcings[[i]][,2]
      YY <- c(y,y[length(y)])[-1]
      forcings[[i]][,2] <- (1-con$f)*y + con$f*YY
    }
  }
  ## all forcings in one vector; adding index to start/end
  
  fmat <- tmat <- NULL
  imat <- rep(1,nf+1)
  
  for (i in 1:nf) {
    # Karline: check for NA in forcing series and remove those
    ii <- apply(forcings[[i]],1,function(x)any(is.na(x)))
    if (sum(ii) > 0) forcings[[i]] <- forcings[[i]][!ii,]
    tmat <- c(tmat, forcings[[i]][,1])
    fmat <- c(fmat, forcings[[i]][,2])
    imat[i+1]<-imat[i]+nrow(forcings[[i]])
  }
  
  storage.mode(tmat) <- storage.mode(fmat) <- "double"
  storage.mode(imat) <- "integer"
  
  # DIRTY trick not to inflate the number of arguments:
  # add method (linear/constant) to imat
  return(list(tmat = tmat, fmat = fmat, imat = c(imat, method),
              ModelForc = ModelForc))
}

### ============================================================================
### Check timelags data set - also passes "dllname" now  (not yet used)
### ============================================================================

checklags <- function (lags, dllname) {
  if (!is.null(lags)) {
    lags$islag = 1L
    if (is.null(lags$mxhist))
      lags$mxhist <- 1e4
    if (lags$mxhist <1)
      lags$mxhist <- 1e4
    lags$mxhist<-as.integer(lags$mxhist)
    if (is.null(lags$interpol))   # 1= hermitian, 2 = higher order interpolation
      lags$interpol <- 1
    lags$interpol<-as.integer(lags$interpol)
    lags$isfun <- 0L
  } else
    lags$islag <- 0L
  return(lags)
}