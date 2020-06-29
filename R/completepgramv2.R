#### regular DFT (in our definition)
#### n^{-1/2} standadized
#### omega_{1,n}~omega_{n,n}
#### can specify the frequencies. default: fundamental freq


fft2 <- function(x, freq = 2*pi*(1:length(x))/length(x)){
  n <- length(x);
  fand <- 2*pi*(1:length(x))/length(x) ##fundamental freq
  
  if(length(setdiff(freq,fand))+length(setdiff(fand,freq))==0){## fundamental freq
    newx <- x[c(n,(1:(n-1)))]
    dft <- fft(newx, inverse=TRUE)
    dft <- dft[c(2:n,1)]/sqrt(n)
    }
  
  if(length(setdiff(freq,fand))+length(setdiff(fand,freq))>0){ ##general freq
    newx <- x[c(n,(1:(n-1)))]
    dftmat <- exp(1i*freq%*%t(1:n))
    dft <- as.vector(dftmat%*%x)/sqrt(n)
    }
  
  return(dft)
  }

##### tapered DFT
#### n^{-1/2} standadized
#### regularization.type = "1": h_{t,n}, are sum to n
#### regularization.type = "2": h_{t,n}, are squared sum to n
#### p: amount of taper: default = 0.1

taperDFT <- function(x,freq = 2*pi*(1:length(x))/length(x), regularization.type = "1", p=0.1){
    xT <- spec.taper(x)
    n <- length(x)
    tap <- spec.taper(rep(1,n),p=p) #output: vector of weights
    
    if(regularization.type=="1"){
      taps <- sum(tap)
      xT <- xT*(n/taps)
      }
    
    if(regularization.type=="2"){
      taps <- sum(tap^2)
      xT <- xT*sqrt(n/taps)}
    
    return(fft2(xT,freq))
    }


#### predictive DFT
#### default: Using YW and order p is chosen by AIC. 
#### n^{-1/2} standadized
#### output: complex value

predictiveDFT <- function(x, freq = 2*pi*(1:length(x))/length(x),taper = FALSE, ar = NULL,...){
  n <- length(x);
  
  if(length(ar)!=0){ ### predetermined ar coefficients
    ord <- length(ar)
    phi <- ar
    }
  
  if((length(ar)==0)&(taper)){ ## no predermined ar coefficients & tapering
    xTap <- spec.taper(x)
    ar.fit <- ar(xTap, ...) ## fitting AR(p). Default is AIC and YW
    ord <- max(ar.fit$order) #order
    phi <- ar.fit$ar # Estimate of Coefficients
    }
  
  if((length(ar)==0)&(!taper)){ ## no predermined ar coefficients & NO tapering
    ar.fit <- ar(x, ...) ## fitting AR(p). Default is AIC and YW
    ord <- max(ar.fit$order) #order
    phi <- ar.fit$ar # Estimate of Coefficients
  }
  
  #freq = 2*pi*(1:n)/n # from 1~n
  m <- length(freq)
  phi.fcn <- rep(NA,m) # AR(p) characteristic function.
  
  if(ord==0){pred = 0}
  if(ord>0){
    for(w in 1:m){
      temp <- exp(-1i*(1:ord)*freq[w])
      phi.fcn[w] <- 1-sum(phi*temp)
      } ## characteristic function end
    
    JL <- rep(NA, m) ## predicted DFT on the left (JL) and right (JR) side.
    JR <- rep(NA, m)
    
    for(w in 1:m){
      tempL <- rep(NA, ord)
      tempR <- rep(NA, ord)
      for(ell in 1:ord){
        tempL[ell] <- sum(phi[ell:ord]*exp(-1i*(0:(ord-ell))*freq[w]))
        tempR[ell] <- sum(phi[ell:ord]*exp(1i*(1:(ord-ell+1))*freq[w]))
      }
      
      JL[w] <- sum(x[1:ord]*tempL)
      JR[w] <- sum(rev(x)[1:ord]*tempR)#*exp(1i*n*freq[w])
      
    }      

    JL <- (JL/phi.fcn)/sqrt(n) ## sqrt(n) standardization
    JR <- (JR/Conj(phi.fcn))/sqrt(n) ## sqrt(n) standardization
    
    pred <- JL +JR} #if(ord>0) end
  
  return(pred)
  }


#### complete DFT
#### default: Using YW and order p is chosen by AIC. 
#### n^{-1/2} standadized
#### output: complex value

completeDFT <- function(x, freq = 2*pi*(1:length(x))/length(x),...){
  temp <- fft2(x,freq) + predictiveDFT(x,freq, ...)

  return(temp)
  }

#### complete DFT
#### default: Using YW and order p is chosen by AIC. 
#### n^{-1/2} standadized
#### NO threshold

complete.pgram <- function(x, freq=2*pi*(1:length(x))/length(x), thres=NULL, ...){
  Le <- completeDFT(x,freq,...)
  Ri <- fft2(x,freq)
  
  rt <- Re(Le*Conj(Ri))
  if(length(thres)==1){
    rt <- ifelse(rt<thres,thres,rt)
    }
  
  return(rt)
}

#### complete DFT
#### default: Using YW and order p is chosen by AIC. 
#### n^{-1/2} standadized
#### NO threshold

tapered.complete.pgram = function(x, freq=2*pi*(1:length(x))/length(x), taperx = NULL, thres=NULL, ...){
  Le <- completeDFT(x,freq,...)
  
  if(length(taperx)==0){
    Ri <- taperDFT(x,freq)
    }
  
  if(length(taperx)!=0){
    Ri <- taperx}
  
  rt <- Re(Le*Conj(Ri))
  if(length(thres)==1){
    rt <- ifelse(rt<thres,thres,rt)
  }  
  return(rt)
}


