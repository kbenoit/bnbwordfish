


## load libraries
library(rstan)
#library(ca)
library(plyr)
library(coda)
library(mcmcplots)
#library(ggplot2)


compare2chainsplot <- function(model,coda=T,star=F,dim2=F,alpha=T){
  if(coda==T) stan.matrix <- as.matrix(model$coda)
  else stan.matrix <- as.matrix(model)
  stan.means1 <- colMeans(stan.matrix[1:(nrow(stan.matrix)/2),])
  stan.means2 <- colMeans(stan.matrix[(1+nrow(stan.matrix)/2):nrow(stan.matrix),])
  
  if(star==T)  s <- "_star\\["  
#   else if(star==F & coda==T) s <- "\\["
#   else s <- "."
  else s <- "\\["
  
  
  if(dim2 == F ) { 
    graphics.off()
    par(mar=rep(2,4))
    par(mfrow=c(2,2))
    
    plot(stan.means1[grep(paste0("^psi",s), names(stan.means1))],stan.means2[grep(paste0("^psi",s), names(stan.means2))], main="psi")
    abline(a=0,b=1)
    plot(stan.means1[grep(paste0("^beta",s), names(stan.means1))],stan.means2[grep(paste0("^beta",s), names(stan.means2))], main="beta")
    abline(a=0,b=1)
    abline(a=0,b=-1)
    plot(stan.means1[grep(paste0("^theta",s), names(stan.means1))],stan.means2[grep(paste0("^theta",s), names(stan.means2))],main="theta")
    abline(a=0,b=1)
    abline(a=0,b=-1)
    if(alpha==T){
      plot(stan.means1[grep(paste0("^alpha",s), names(stan.means1))],stan.means2[grep(paste0("^alpha",s), names(stan.means2))],main="alpha")
      abline(a=0,b=1)
    }
  } else { #  if(dim2 == T) 
    graphics.off()
    par(mar=rep(2,4))
    par(mfrow=c(3,2))
    plot(stan.means1[grep(paste0("^psi",s), names(stan.means1))],stan.means2[grep(paste0("^psi",s), names(stan.means2))],main="psi")
    abline(a=0,b=1)
    plot(stan.means1[grep(paste0("^beta1",s), names(stan.means1))],stan.means2[grep(paste0("^beta1",s), names(stan.means2))],main="beta1")
    abline(a=0,b=1)
    plot(stan.means1[grep(paste0("^theta1",s), names(stan.means1))],stan.means2[grep(paste0("^theta1",s), names(stan.means2))],main="theta1")
    abline(a=0,b=1)
    plot(stan.means1[grep(paste0("^beta2",s), names(stan.means1))],stan.means2[grep(paste0("^beta2",s), names(stan.means2))],main="beta2")
    abline(a=0,b=1)
    plot(stan.means1[grep(paste0("^theta2",s), names(stan.means1))],stan.means2[grep(paste0("^theta2",s), names(stan.means2))],main="theta2")
    abline(a=0,b=1)
    if(alpha==T){
      plot(stan.means1[grep(paste0("^alpha",s), names(stan.means1))],stan.means2[grep(paste0("^alpha",s), names(stan.means2))],main="alpha")
      abline(a=0,b=1)
    }
  } 
}



stantocoda <- function(x) {
  require(plyr)
  require(coda)
  stopifnot(class(x)=="stanfit")
  coda <- rstan::extract(x, inc_warmup=F, permuted=F) # extract array
  coda <- do.call(mcmc.list, alply(coda[,,-(length(coda[1,1,]))], 2, mcmc)) # remove lp_, convert to coda
}

### Stan init function
# generates inits based on corresp. analysis

as.stan.init <- function(data.s=data.s.filled, mod=c("1d","2d","2dnb","1dp","1dnb"),noisesd=inits.noisesd, chain_id=1,
                         n.bothzero=NULL, n.dim1only=NULL, n.dim2only=NULL, # for 2d model
                         group=NULL, init.model=NULL , phi=NULL # for 1d hier 
){ 
  # Note: the init function uses the data sorted as necessary by wordfishSTAN
  if(mod!="1dp") { # estimate corresp. analysis
    require(ca)
    myca <- ca(data.s) # nd=1 option currently not allowed, crashes
    alpha <- log(rowMeans(data.s))
    psi <- (log(myca$colmass)-mean(log(myca$colmass))) 
    psi <- psi-psi[ncol(data.s)]
    
    if(mod=="1d" | mod=="1dnb" ) {
      myca$colcoord <- myca$colcoord[,1]
      myca$rowcoord <- myca$rowcoord[,1]
      
      
      if (myca$colcoord[2] < myca$colcoord[1])  { # reflect and standardize if in reverse direction 
        theta <- as.numeric(scale(-myca$rowcoord)) 
        beta <- as.numeric(scale(-myca$colcoord)) 
      }  else { #  standardize as is  if in correct direction
        theta <- as.numeric(scale(myca$rowcoord))
        beta <-  as.numeric(scale(myca$colcoord))
        #beta <- as.numeric(-myca$colcoord*sd(myca$rowcoord)) # multiply accordingly to keep product constant
      }
      # optionally add random noise to beta and theta inits (for better check of convergence)
      if (noisesd > 0) { 
        beta <- beta+rnorm(ncol(data.s),0,sd=noisesd)
        theta <- theta+rnorm(nrow(data.s),0,sd=noisesd)
      } 
      
      # truncate extreme values of beta and theta
      beta[beta < -3] <- -3
      beta[beta > 3] <- 3
      theta[theta < -3] <- -3
      theta[theta > 3] <- 3
      
      ll <- list( mu_psi = mean(psi[-ncol(data.s)]),
                  sigma_psi = sd(psi[-ncol(data.s)]),
                  sigma_beta = sd(beta),
                  alpha=alpha,
                  psi_rest = psi[-ncol(data.s)],     
                  beta = beta,
                  theta= theta
      )      
      if(mod=="1dnb") {  ll$phi_inv <- runif(ncol(data.s),0,20) }
      
    } else if(mod=="2d" | mod=="2dnb"){
      if (myca$colcoord[2,1] > myca$colcoord[1,1])  { # reflect and standardize if in reverse direction 
        theta1 <- as.numeric(scale(-myca$rowcoord[,1])) 
        beta1 <- as.numeric(scale(-myca$colcoord[,1])) 
      }  else { #  standardize as is  if in correct direction
        theta1 <- as.numeric(scale(myca$rowcoord[,1]))
        beta1 <-  as.numeric(scale(myca$colcoord[,1]))
      }
      if (myca$colcoord[2,2] < myca$colcoord[1,2])  { # reflect and standardize if in reverse direction 
        theta2 <- as.numeric(scale(-myca$rowcoord[,2])) 
        beta2 <- as.numeric(scale(-myca$colcoord[,2])) 
      }  else { #  standardize as is  if in correct direction
        theta2 <- as.numeric(scale(myca$rowcoord[,2]))
        beta2 <-  as.numeric(scale(myca$colcoord[,2]))
      }
      beta1 <- beta1-beta1[2] # subtract ref cat
      beta2 <- beta2-beta2[1]
      
      # truncate extreme values of beta and theta
      beta1[beta1 < -3] <- -3
      beta1[beta1 > 3] <- 3
      theta1[theta1 < -3] <- -3
      theta1[theta1 > 3] <- 3
      beta2[beta2 < -3] <- -3
      beta2[beta2 > 3] <- 3
      theta2[theta2 < -3] <- -3
      theta2[theta2 > 3] <- 3
      
      # optionally add random noise to beta and theta inits to check convergence
      if (noisesd > 0) { 
        beta1 <- beta1+rnorm(ncol(data.s),0,sd=noisesd)
        theta1 <- theta1+rnorm(nrow(data.s),0,sd=noisesd)
        beta2 <- beta2+rnorm(ncol(data.s),0,sd=noisesd)
        theta2 <- theta2+rnorm(nrow(data.s),0,sd=noisesd)
      }
      
      # remove constrained items
      
      beta1 <- beta1[-c(1:(2+n.bothzero),(2+n.bothzero+n.dim1only+1):(2+n.bothzero+n.dim1only+n.dim2only))]
      beta2 <- beta2[-c(1:(2+n.bothzero),(2+n.bothzero+1):(2+n.bothzero+n.dim1only))]
      
      
      ll <- list( mu_psi = mean(psi[-1]),
                  sigma_psi = sd(psi[-1]),
                  mu_beta1 = rnorm(1,sd=2),
                  mu_beta2 = rnorm(1,sd=2),
                  sigma_beta1 = runif(1,0,4),
                  sigma_beta2 = runif(1,0,4),
                  alpha=alpha,
                  psi_rest = psi[-1],
                  beta1_rest = beta1, # note: these are not relative to constrained ones
                  beta2_rest = beta2,
                  theta1= theta1,
                  theta2=theta2,
                  sigma_theta1 = runif(1,0,4),
                  sigma_theta2 = runif(1,0,4)
      )   
      if(mod=="2dnb") {  
        ll$phi_inv <- runif(ncol(data.s),0,50)
      }
    }
  } else if(mod=="1dp") {
    
    model <- init.model
    theta <- model$theta
    beta <- model$beta
    if (noisesd > 0) {  # optionally add random noise to beta and theta inits to check convergence
      #       beta <- beta+rnorm(length(beta),0,sd=noisesd)
      theta <- theta+rnorm(length(theta),0,sd=noisesd)
    } 
    n_free <- length(beta)-1
    n_group <- length(unique(group))
    # the two lowest group values are the two reference groups
    # ref.groups <- c(unique(group)[order(unique(group))][1],unique(group)[order(unique(group))][2]) 
    # alpha_e <- model$alpha-mean(model$alpha) # note: is just distance from grand mean (not country-level)
    mu_theta <- aggregate(theta,by=list(group),mean)
    mu_theta <- mu_theta$x # partxx
    
    mu_psi <- model$psi[-ncol(data.s)]-model$psi[ncol(data.s)]
    mu_beta <- model$beta[-ncol(data.s)]-model$beta[ncol(data.s)] 
    # psi_rest <- matrix(rnorm(n_group*n_free,mu_psi,sd=1),ncol=n_group)   
    # beta_rest <- matrix(rnorm(n_group*n_free,mu_beta,sd=1),ncol=n_group) 
    
    ll <- list( alpha=model$alpha,
                # alpha_e = alpha_e, # rnorm(length(theta),mean=0,sd=noisesd),
                mu_psi = mu_psi,
                # psi_rest = psi_rest,
                psi_e = matrix(rnorm(n_group*n_free,mean=0,sd=noisesd),ncol=n_group),
                beta_e = matrix(rnorm(n_group*n_free,mean=0,sd=noisesd),ncol=n_group),
                # beta_rest = matrix(rnorm(n_group*n_free,mean=mu_beta,sd=noisesd),ncol=n_group),
                sigma_psi = runif(n_free,.5,3),               
                mu_beta = mu_beta,  
                # beta = beta,   
                sigma_beta= runif(n_free,.5,3),               
                mu_theta = mu_theta, #partxx
                theta= theta  ,
                cauchy_scale_psi = runif(1,0,10),       
                cauchy_scale_beta = runif(1,0,10),
                nu_psi = runif(1,0,100),       
                nu_beta = runif(1,0,100),
                epsilon = rep(0,nrow(data.s)*ncol(data.s)), 
                sigma_epsilon = runif(n_group,0,3)
                
                # tau_psi = runif(1,0,100),       
                # tau_beta = runif(1,0,100)                
    )   
    if(phi=="w") {
      ll$phi_inv <- runif(length(beta),0,300) 
    } else if(phi=="c") {
      ll$phi_inv <- runif(n_group,0,300) 
    }  else if(phi=="wc") {
      ll$phi_inv <- matrix(runif(n_group*ncol(data.s),0,300),ncol=n_group)
    } else if(phi=="add") {
      ll$phi_inv_t_rest <- runif(length(beta)-1,0,15) 
      ll$phi_inv_g <- runif(n_group,0,15) 
    } else if(phi=="d") {
      ll$phi_inv <- runif(nrow(data.s),0,300) 
    }
  }
  return(ll)
}



#######################################################
##### Stan implementation: using the above model files
#######################################################


wordfishSTAN1d <-  function(data, reuse=NULL, stan.code, 
                            dir.terms=c(1,2), dir.beta=TRUE,
                            psi.ref=ncol(data), 
                            iter=1000, warmup=500, n.thin=1, n.chains=2,
                            group.data=NULL, negbin=F, 
                            monitor.alpha=F, monitor.aux=F, return.raw.stan=F,
                            init.type=c("ca","random"),
                            inits.noisesd=0) {
  
  if(dir.beta==FALSE) stopifnot(is.null(group.data)) # direction must be set through beta if beta varies across groups
  mod <- "1d"
  if(negbin==T) { mod <- "1dnb"}
  
  # sort data
  # move one item into last col as reference item for psi
  index.col <- c((1:ncol(data))[-psi.ref], psi.ref)    
  data.s <- data[,index.col]
  # replace NA with 0 for ca for inits  
  data.s.filled <- data.s
  data.s.filled[is.na(data.s.filled)] <- 0
  
  # write data in long format (one doc after another)
  doc <- rep(1:nrow(data.s), each=ncol(data.s))
  term <- rep(1:ncol(data.s), nrow(data.s))
  y <- as.vector(t(as.matrix(data.s)))
  
  longdata <- list(D = nrow(data.s),
                   T = ncol(data.s),
                   N = nrow(data.s) * ncol(data.s),
                   doc = doc,
                   term = term,
                   y = y
  )
  if(!is.null(group.data)){
    longdata$grdoc <- group.data
    longdata$gr <- rep(group.data,each=ncol(data.s))
    longdata$G <- max(group.data)
  } 
  
  # identify NA cells
  
  if(any(is.na(y))) {
    del <- which(is.na(y))          
    longdata$N <- longdata$N-length(del)
    longdata$doc <- longdata$doc[-(del)]
    longdata$term <- longdata$term[-(del)]
    longdata$y <- longdata$y[-(del)]
  }
  
  
  # set monitor 
  pars <-  c("theta","psi", "beta")

  if (monitor.alpha== T ) { pars <- c(pars, "alpha")   }
  if (monitor.aux == T  ) { pars <- c(pars,  "sigma_beta",  "sigma_psi", "mu_psi")  } # 
  if (mod == "1dnb" ) { pars <- c(pars,"phi_inv") }

  require(rstan)
  #require(ca)
  require(coda)

  start.time <- proc.time()
  

    if(is.null(reuse)){
      if(init.type!="random") {
        
        
        
        stanmodel <- stan(model_code = stan.code, data = longdata,
                          init=lapply(1:n.chains, function(id) as.stan.init(data.s.filled,mod=mod,
                                                                            noisesd=inits.noisesd,
                                                                            chain_id= id)), 
                          iter=iter,warmup=warmup,thin=n.thin, chains = n.chains, pars = pars,
                          cores=min(n.chains,parallel::detectCores())
        )
      } else {
        stanmodel <- stan(model_code = stan.code, data = longdata,
                          init="random",
                          iter=iter,warmup=warmup,thin=n.thin, chains = n.chains, pars = pars,
                          cores=min(n.chains,parallel::detectCores()) )
      }
      
    } else {
      if(init.type!="random") {
        stanmodel <- stan(fit=reuse, data = longdata,
                          init=lapply(1:n.chains, function(id) as.stan.init(data.s.filled,mod=mod,
                                                                            noisesd=inits.noisesd,
                                                                            chain_id= id)),                         
                          iter=iter,warmup=warmup,thin=n.thin, chains = n.chains, pars = pars,
                          cores=min(n.chains,parallel::detectCores()))
      }  else {
        stanmodel <- stan(fit=reuse, data = longdata,
                          init="random",
                          iter=iter,warmup=warmup,thin=n.thin, chains = n.chains, pars = pars,
                          cores=min(n.chains,parallel::detectCores()))
      }
    }

  
  
  # convert output


duration <- (time.elapsed= proc.time() - start.time[3]) 
  s <- as.matrix(stanmodel) # draws in rows, parameters in columns 
  thetaind <- grep("^theta",colnames(s))
  sstar <- s
  
  if(is.null(group.data)){
    if(dir.beta==TRUE) reflect <- ifelse(s[,paste0("beta[",dir.terms[1],"]")] > s[,paste0("beta[",dir.terms[2],"]")] , -1,1) 
    else reflect <- ifelse(s[,paste0("theta[",dir.terms[1],"]")] > s[,paste0("theta[",dir.terms[2],"]")] , -1,1) 
    
    meantheta <- apply(s[,thetaind],1,mean) 
    sdtheta <- apply(s[,thetaind],1,sd)
    betaind <- grep("^beta",colnames(s))
    meanbeta <- apply(s[,betaind],1,mean)
    psiind <- grep("^psi",colnames(s))
    meanpsi <- apply(s[,psiind],1,mean)
    
    # transform beta
    sstar[,betaind] <- reflect*((sstar[,betaind]-meanbeta)*sdtheta)
    # transform theta
    sstar[,thetaind] <- reflect*((sstar[,thetaind]-meantheta)/sdtheta)
    # transform psi
    sstar[,psiind] <- sstar[,psiind] - meanpsi + meantheta*(s[,betaind]-meanbeta) # (note the use of raw beta from s, not sstar, here)
    
    if(monitor.alpha==TRUE) {
      alphaind <- grep("^alpha",colnames(s))
      sstar[,alphaind] <- sstar[,alphaind] + meanpsi + meanbeta*s[,thetaind] # (note the use of raw theta from s, not sstar, here)
    } 
    
    # re-order item parameters to original order (move ref cat from end back into order)
    sstar[,betaind] <- sstar[,betaind][,order(index.col)]
    sstar[,psiind] <- sstar[,psiind][,order(index.col)]
     
  } else if(!is.null(group.data)) {
    
    
      for(g in 1:longdata$G){
      # this will be -1 if the order in the group is reversed  
      reflect <- ifelse(s[,paste0("beta[",g,",",dir.terms[1],"]")] > s[,paste0("beta[",g,",",dir.terms[2],"]")] , -1,1) 
        
      meantheta <- apply(s[,thetaind[group.data==g]],1,mean) 
      sdtheta <- apply(s[,thetaind[group.data==g]],1,sd)
      betagind <- grep(paste0("^beta\\[",g,","),colnames(s))
      psigind <- grep(paste0("^psi.*",g,"\\]"),colnames(s))
      meanbeta <- apply(s[,betagind],1,mean)
      meanpsi <- apply(s[,psigind],1,mean)
      
      # transform beta
      sstar[,betagind] <- reflect*((sstar[,betagind]-meanbeta)*sdtheta)
      # transform theta
      sstar[,thetaind[group.data==g]] <- reflect*((sstar[,thetaind[group.data==g]]-meantheta)/sdtheta)
      # transform psi
      sstar[,psigind] <- sstar[,psigind] - meanpsi + meantheta*(s[,betagind]-meanbeta) # (note the use of raw beta from s, not sstar, here)
      
      if(monitor.alpha==TRUE) {
        alphaind <- grep("^alpha",colnames(s))
        sstar[,alphaind[group.data==g]] <- sstar[,alphaind[group.data==g]] + meanpsi + meanbeta*s[,thetaind[group.data==g]] # (note the use of raw theta from s, not sstar, here)
      }
  # re-order item parameters to original order (move ref cat from end back into order)
  # (after this, colnames will refer to original CMP order)
      sstar[,grep(paste0("^beta\\[",g,","),colnames(s))] <- sstar[,grep(paste0("^beta\\[",g,","),colnames(s))][,order(index.col)]
      sstar[,grep(paste0("^psi.*",g,"\\]"),colnames(s))] <- sstar[,grep(paste0("^psi.*",g,"\\]"),colnames(s))][,order(index.col)]
      }
  }
  
  # NOTE insert pearson and pr0 again?
  # reorder phi (which is not group-specific)
  if(negbin==TRUE) sstar[,grep("phi",colnames(s))] <- sstar[,grep("phi",colnames(s))][,order(index.col)]
    

  # turn sstar into a coda package object (mcmclist)
  # REM need to adapt for more than 2 chains, with laply?
  sscoda <- coda::as.mcmc.list(list(coda::mcmc(sstar[1:(nrow(sstar)/2),]),
                                    coda::mcmc(sstar[(nrow(sstar)/2+1):(nrow(sstar)),])
                                    ))
  # create a 3-D array from sstar for rstan::monitor
  # sstar <- as.matrix(sscoda)
  # REM need to make this more general for any n.chain, not just 2 like here:
  # create a vector (takes col by col), so this is chain 1 param 1 draws 1/N, chain 1 param 2 draws...
  tmp <- c(as.vector(sstar[1:(nrow(sstar)/2),]),as.vector(sstar[(nrow(sstar)/2+1):(nrow(sstar)),]) )
  # make an array (will be of dim:  iter/draws * param * chains)
  ssarr <- array(tmp, dim=c(nrow(sstar)/n.chains, ncol(sstar), n.chains ),
                 dimnames=list(draws=NULL,params=colnames(sstar),chains=1:n.chains)) 

  ssarr <- aperm(ssarr,c(1,3,2)) # this is (iterations/draws * chains * parameters) as needed for monitor
  tableout <- monitor(ssarr, warmup=0, print=FALSE )
  tableout <- data.frame(tableout)[, c(1,4,8:10)] # (the subsetting doesn't work without the data.frame - ?)

  retval <- list(coda=sscoda, summary=tableout,duration=duration)

  if(return.raw.stan==TRUE) retval$stanfit <- stanmodel
  return(retval)
}

#####################
### 2d R function ###
#####################

wordfishSTAN2d <- function(data, reuse=NA, # REM user needs to make sure that re-used object has same type of constraints/Stan code
                           base.terms=c(1,2), bothzero.terms=3,
                           dim1only.terms=c(4,5),  dim2only.terms=c(6,7), 
                           negbin=F,
                           iter=1000, warmup=500, n.thin=1, n.chains=2,
                           monitor.alpha=F, monitor.aux=F, return.raw.stan=F,
                           init.type=c("ca","random"),
                           inits.noisesd=0 ) {
  
  mod <- "2d"
  if(negbin==T) { mod <- "2dnb"}
  
  # determine pattern of constraints and sort data accordingly
  
  
  bothfree.terms <- (1:ncol(data))[-c(base.terms,bothzero.terms,dim1only.terms,dim2only.terms)]
  n.bothfree <- length(bothfree.terms)
  n.bothzero <- length(bothzero.terms)  
  n.dim1only <- length(dim1only.terms)
  n.dim2only <- length(dim2only.terms)
  
  # Type 1 beta pattern: loading on both present
  beta.both <- (n.bothfree > 0)
  stan.code <- stanPoisson2dBothYes
  if(beta.both==F & negbin==F) stan.code <- stanPoisson2dBothNo
  if(beta.both==T & negbin==T) stan.code <- stanNegbin2dWordBothYes
  if(beta.both==F & negbin==T) stan.code <- stanNegbin2dWordBothNo
  
  # sort data
  index.col  <- c(base.terms,bothzero.terms,dim1only.terms,dim2only.terms,bothfree.terms) # concatenation assumes that NULL if does not apply
  data.s <- data[,index.col]
  
  # write data in long format             
  longdata <- list(D = nrow(data.s),
                   T = ncol(data.s),
                   N = nrow(data.s) * ncol(data.s),
                   doc = rep(1:nrow(data.s), each=ncol(data.s)),
                   term = rep(1:ncol(data.s), nrow(data.s)),
                   y = as.vector(t(as.matrix(data.s))),
                   nbz= n.bothzero,
                   nd1 = n.dim1only,
                   nd2 = n.dim2only
                   #, mu_theta = c(0,0),
                   # Sigma_theta=matrix(c(5,0,0,5),nrow=2,byrow=T)
  )
  
  # set monitor 
  pars <- c("theta1","theta2",  "psi", "beta1", "beta2",
              "sigma_beta1", "sigma_beta2",   "sigma_theta1", "sigma_theta2"
            ) 
  
  if (monitor.aux== T ) { pars <- c(pars, "mu_psi","sigma_psi","mu_beta1","mu_beta2")   }
  if (monitor.alpha== T ) { pars <- c(pars, "alpha")   }
  if (negbin== T ) { pars <- c(pars, "phi_inv")   }
  
  
  require(rstan)
  #require(ca)
  require(coda)
  
  start.time <- proc.time()
  
  
    if(is.na(reuse)==T){
      if(init.type!="random") {
        stanmodel <- stan(model_code = stan.code, data = longdata,
                          init=lapply(1:n.chains, function(id) as.stan.init(data.s,mod=mod,
                                                                            n.bothzero=n.bothzero, n.dim1only=n.dim1only, 
                                                                            n.dim2only=n.dim2only,
                                                                            noisesd=inits.noisesd,
                                                                            chain_id= id)), 
                          iter=iter,warmup=warmup,thin=n.thin, chains = n.chains, pars = pars,
                          cores=min(n.chains,parallel::detectCores())
        )
      } else {
        stanmodel <- stan(model_code = stan.code, data = longdata,
                          init="random",
                          iter=iter,warmup=warmup,thin=n.thin, chains = n.chains, pars = pars,
                          cores=min(n.chains,parallel::detectCores()) )
      }
      
    } else {
      if(init.type!="random") {
        stanmodel <- stan(fit=reuse, data = longdata,
                          init=lapply(1:n.chains, function(id) as.stan.init(data.s,mod=mod,
                                                                            noisesd=inits.noisesd,
                                                                            n.bothzero=n.bothzero, n.dim1only=n.dim1only, 
                                                                            n.dim2only=n.dim2only,
                                                                            chain_id= id)),                         
                          iter=iter,warmup=warmup,thin=n.thin, chains = n.chains, pars = pars,
                          cores=min(n.chains,parallel::detectCores()) )
      }  else {
        stanmodel <- stan(fit=reuse, data = longdata,
                          init="random",
                          iter=iter,warmup=warmup,thin=n.thin, chains = n.chains, pars = pars,
                          cores=min(n.chains,parallel::detectCores()) )
      }
    }
    
  
  # convert output
  duration <- (time.elapsed= proc.time() - start.time[3]) 
  s <- as.matrix(stanmodel) # draws in rows, parameters in columns 
  
  theta1ind <- grep("^theta1",colnames(s))
  theta2ind <- grep("^theta2",colnames(s))
  beta1ind <- grep("^beta1",colnames(s))
  beta2ind <- grep("^beta2",colnames(s))
  psiind <- grep("^psi",colnames(s))
  
  sstar <- s
  
  meanpsi <- apply(s[,psiind],1,mean)
  sdtheta1 <- apply(s[,theta1ind],1,sd)
  sdtheta2 <- apply(s[,theta2ind],1,sd)
  
  # transform beta
    sstar[,beta1ind] <- sstar[,beta1ind] * sdtheta1
  sstar[,beta2ind] <- sstar[,beta2ind] * sdtheta2
    # transform theta
  sstar[,theta1ind] <- sstar[,theta1ind]/sdtheta1
  sstar[,theta2ind] <- sstar[,theta2ind]/sdtheta2
    # transform psi
    sstar[,psiind] <- sstar[,psiind] - meanpsi 
    
    if(monitor.alpha==TRUE) {
      alphaind <- grep("^alpha",colnames(s))
      sstar[,alphaind] <- sstar[,alphaind] + meanpsi 
    } 
    
    # re-order item parameters to original order 
    # (after this, colnames will refer to original CMP order)
    sstar[,beta1ind] <- sstar[,beta1ind][,order(index.col)]
    sstar[,beta2ind] <- sstar[,beta2ind][,order(index.col)]
    
    sstar[,psiind] <- sstar[,psiind][,order(index.col)]
    
  
  
  # NOTE insert pearson and pr0 again?
  # reorder phi (which is not group-specific)
    if(negbin==TRUE){
    sstar[,grep("phi",colnames(s))] <- sstar[,grep("phi",colnames(s))][,order(index.col)]
    }
  
  # turn sstar into a coda package object (mcmclist)
  # REM adapt chains, with laply?
  sscoda <- coda::as.mcmc.list(list(coda::mcmc(sstar[1:(nrow(sstar)/2),]),
                                    coda::mcmc(sstar[(nrow(sstar)/2+1):(nrow(sstar)),])
  ))
  # create a 3-D array from sstar for rstan::monitor
  # sstar <- as.matrix(sscoda)
  # REM need to make this more general for any n.chain, not just 2 like here:
  # create a vector (takes col by col), so this is chain 1 param 1 draws 1/N, chain 1 param 2 draws...
  tmp <- c(as.vector(sstar[1:(nrow(sstar)/2),]),as.vector(sstar[(nrow(sstar)/2+1):(nrow(sstar)),]) )
  
  # make an array (will be of dim:  iter/draws * param * chains)
  ssarr <- array(tmp, dim=c(nrow(sstar)/n.chains, ncol(sstar), n.chains ),
                 dimnames=list(draws=NULL,params=colnames(sstar),chains=1:n.chains)) 
  
  ssarr <- aperm(ssarr,c(1,3,2)) # this is (iterations/draws * chains * parameters) as needed for monitor
  tableout <- monitor(ssarr, warmup=0, print=FALSE )
  tableout <-  data.frame(tableout)[, c(1,4,8:10)]
  
  retval <- list(coda=sscoda, summary=tableout,duration=duration)
  
  if(return.raw.stan==TRUE) retval$stanfit <- stanmodel
  
 
  return(retval)
}



ggPlotLambdasbyRile <- function(data, 
                                 pdf=F, plotname=NULL) {
  
  stopifnot(is.null(plotname) | is.character(plotname))
  # stopifnot(((rile==T | eco == T) & is.null(group) )  | (rile==F & eco==F & !is.null(group))) 
  
  require(ggplot2)
  require(grid)
  
  data <- data[order(data$rile, data$lambda),]
  
  data <- rbind(data[1:13,],data[1:2,],data[14:43,],data[1:2,],data[44:56,],data[1,])
  data[c(14:15,46:47,61),c("lambda","cilo","cihi")] <- NA
  data$cats <- 1:nrow(data)
  data$catlab[c(14:15,46:47,61)] <- c("LEFT","","NEUTRAL","","RIGHT")
  data$catlab[-c(14:15,46:47,61)] <- paste("   ",data$catlab[-c(14:15,46:47,61)])
  
  ## caterpillar dotplot
  p <- ggplot(data, aes(lambda,cats)) +  #  y=reorder(cats,betas)
    
    xlab(expression(beta)) + ylab("Policy category")  +
    scale_x_continuous(breaks = seq(-3,3,by=1), limits=c(-3.3,3.3)  )  +
    scale_y_continuous(breaks=1:nrow(data),labels=data$catlab)  +
    geom_blank()  +
    theme(legend.position="none") +
      geom_errorbarh(aes(xmin=cilo, xmax=cihi), height=0) +
      geom_point(size=1.3) + 
    geom_vline(xintercept=0, linetype = "dashed", color="grey40") 
  
  
  p <- p  + theme_bw() +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dotted"))  +
    theme(panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank()) +
    theme(panel.grid.minor.x = element_blank()) +
    theme(plot.background = element_blank()) +
    theme(axis.text=element_text(size=12) ) +
    theme(axis.text.y=element_text(size=12, hjust=0) ) +
    theme(axis.title=element_text(size=14, face="bold") ) # axis text size to 7 from 14
  if(pdf == F) {
    return(p)
  } else {
    ggsave(filename=paste("figures/", plotname, ".pdf", sep=""), height=12, width=8)
  }
}

ggPlotPhi <- function(data, 
                                pdf=F, plotname=NULL) {
    
    stopifnot(is.null(plotname) | is.character(plotname))
    # stopifnot(((rile==T | eco == T) & is.null(group) )  | (rile==F & eco==F & !is.null(group))) 
    
    require(ggplot2)
    require(grid)
    

    data <- data[order(data$phi_inv),]
    

    data$cats <- 1:nrow(data)

    ## caterpillar dotplot
    p <- ggplot(data, aes(phi_inv,cats, shape=RILE, colour=RILE)) +  #  y=reorder(cats,betas)
        
        xlab(expression(phi^-1)) + ylab("Policy category")  +
        scale_x_continuous(breaks = seq(0,20,by=5), limits=c(-.3,21.5)  )  +
        scale_y_continuous(breaks=1:nrow(data),labels=data$catlab)  +
        scale_color_manual(values=c("red","black","blue")) +
        scale_shape_manual(values=c(15,1,17)) +
        geom_blank()  +
        geom_errorbarh(aes(xmin=cilo, xmax=cihi), height=0) +
        geom_point(size=2) +
        geom_text(aes(x=21, label=propzerolab), colour="black", size=3)
    
    p <- p  + theme_bw() +
        theme(panel.grid.major.y = element_line(color="grey", linetype="dotted"))  +
        theme(panel.grid.minor.y = element_blank()) +
        theme(panel.grid.major.x = element_blank()) +
        theme(panel.grid.minor.x = element_blank()) +
        theme(plot.background = element_blank()) +
        theme(axis.text=element_text(size=12) ) +
        theme(axis.text.y=element_text(size=12, hjust=0) ) +
        theme(axis.title=element_text(size=14, face="bold") ) + # axis text size to 7 from 14
        theme(legend.position="bottom") +
    if(pdf == F) {
        return(p)
    } else {
        ggsave(filename=paste("figures/", plotname, ".pdf", sep=""), height=12, width=8)
    }
}

ggPlotLambdas2D <- function(data, dim=1,
                                pdf=F, plotname=NULL) {
    
    stopifnot(is.null(plotname) | is.character(plotname))
    # stopifnot(((rile==T | eco == T) & is.null(group) )  | (rile==F & eco==F & !is.null(group))) 
    
    require(ggplot2)
    require(grid)
    
    data <- data[data[,paste0("pure",ifelse(dim==1,2,1))] != 1,]
    data <- data[order( data[,paste0("lambda",dim)]),]

    tmpind <- c((sum(data[,paste0("pure",dim)]==1)+1),(sum(data[,paste0("pure",dim)]==1)+2),nrow(data)+3)
 
 
    data <- rbind(data[data[,paste0("pure",dim)] == 1,],NA, NA,data[data[,paste0("pure",dim)] == 0,],NA )
    data$cats <- 1:nrow(data)
    
    data$catlab[tmpind] <- c(paste("PURE", ifelse(dim==1,"ECONOMIC","SOCIAL")),"","ECONOMIC AND SOCIAL")
    data$catlab[-tmpind] <- paste("   ",data$catlab[-tmpind])

    ## caterpillar dotplot
    p <- ggplot(data, aes_string(paste0("lambda",dim),"cats")) +  #  y=reorder(cats,betas)
        
        xlab(expression(beta)) + ylab("Policy category")  +
        scale_x_continuous(breaks = seq(-3,3,by=1), limits=c(-2.35,2.35)  )  +
        scale_y_continuous(breaks=1:nrow(data),labels=data$catlab)  +
        geom_blank()  +
        theme(legend.position="none") +
        geom_errorbarh(aes_string(xmin=paste0("cilo",dim), xmax=paste0("cihi",dim)), height=0) +
        geom_point(size=1.5) + 
        geom_vline(xintercept=0, linetype = "dashed", color="grey40") 
    
    
    p <- p  + theme_bw() +
        theme(panel.grid.major.y = element_line(color="grey", linetype="dotted"))  +
        theme(panel.grid.minor.y = element_blank()) +
        theme(panel.grid.major.x = element_blank()) +
        theme(panel.grid.minor.x = element_blank()) +
        theme(plot.background = element_blank()) +
        theme(axis.text=element_text(size=12) ) +
        theme(axis.text.y=element_text(size=12, hjust=0) ) +
        theme(axis.title=element_text(size=14, face="bold") ) # axis text size to 7 from 14
    if(pdf == F) {
        return(p)
    } else {
        ggsave(filename=paste("figures/", plotname, ".pdf", sep=""), height=12, width=8)
    }
}


ggPlotLambdasbyGroup <- function(data, dim=0,
                                 pdf=F, plotname=NULL) {
  
  stopifnot(is.null(plotname) | is.character(plotname))
  # stopifnot(((rile==T | eco == T) & is.null(group) )  | (rile==F & eco==F & !is.null(group))) 
  
  require(ggplot2)
  require(grid)
  
  target <- "lambda"
  if(dim != 0) target <- paste0(target,dim)
  data <- data[order(data[,target]),]
  data$cats <- 1:nrow(data)
  
  if(dim==0) {
    data$shadeymin <- ifelse(data$rile != "RILE neutral", data$cats-.5,NA)
    data$shadeymax <- ifelse(data$rile != "RILE neutral", data$cats+.5,NA)
    data$shadexmin <- ifelse(data$rile == "RILE right",0,NA)
    data$shadexmax <- ifelse(data$rile == "RILE right", Inf,NA)
    data$shadexmin[data$rile == "RILE left"] <- -Inf
    data$shadexmax[data$rile == "RILE left"] <- 0
  } 
  
  ## caterpillar dotplot
  p <- ggplot(data, aes(get(target),cats)) +  #  y=reorder(cats,betas)
    #  facet_grid( gr  ~ . , scales="free_y", space="free") + # 
    xlab(expression(beta)) + ylab("Policy category")  +
    scale_x_continuous(breaks = seq(-3,3,by=1), limits=c(-3,3)  )  +
    scale_y_continuous(breaks=1:nrow(data),labels=data$catlab)  +
    geom_blank()  +
    theme(legend.position="none") +
    geom_vline(xintercept=0, linetype = "dashed")
  
  if(dim == 0) {
    p <- p + geom_rect(aes(ymin=shadeymin,ymax=shadeymax, xmin=shadexmin, xmax=shadexmax),
                       fill="grey90") +
      geom_errorbarh(aes(xmin=cilo, xmax=cihi,color=gr), width=0) +
      scale_colour_manual(values = c("black","blue")) + 
      geom_point(aes(shape=gr,colour=gr),size=1.3) + scale_shape_manual(values = c(19,24,4))## # older: c(19,15,24)
  } else if (dim==1) {
    p <- p  + geom_errorbarh(aes(xmin=lam1cilo, xmax=lam1cihi,color=gr), width=0 ) +
      scale_colour_manual(values = c("black","blue")) +
      geom_point(aes(shape=gr,colour=gr),size=1.3) + scale_shape_manual(values = c(19,4))
  } else if (dim==2) { 
    p <- p  + geom_errorbarh(aes(xmin=lam2cilo, xmax=lam2cihi,color=gr), width=0 ) +
      scale_colour_manual(values = c("red","blue")) +
      geom_point(aes(shape=gr,colour=gr),size=1.3) + scale_shape_manual(values = c(24,4))
  }
  ##   specified shapes: solid ones for social and economic, hollow for others
  
  
  p <- p  +  theme(panel.grid.major.y = element_line(color="grey", linetype="dotted"))  +
    theme(panel.grid.minor.y = element_blank()) +
    theme(plot.background = element_blank()) +
    theme(axis.text=element_text(size=7) ) +
    theme(axis.title=element_text(size=9, face="bold") ) + theme_bw()# axis text size to 7 from 14
  if(pdf == F) {
    return(p)
  } else {
    ggsave(filename=paste("figures/", plotname, ".pdf", sep=""), height=12, width=8)
  }
}

ggPlotDiD <- function(data,target.varname,sort.varname, sig.varname,
                      xlim=c(-3,3),pdf=F, plotname=NULL) {
  
  stopifnot(is.null(plotname) | is.character(plotname))
  
  require(ggplot2)
  require(grid)
  
  # sort by target values in one group for plotting
  data$tmp <- ifelse(data$hiexpos==1 & data$post == 1, data[,sort.varname], NA)
  data$tmp <- ave(data$tmp, data$catnum, FUN=function(x) mean(x,na.rm=TRUE))
  data <- data[order(data$tmp, data$group),]
  data$cats <- rep(1:56, each=length(unique(data$group)))
  
  data$shadeymin <- ifelse(data$group == 1 & data[,sig.varname]==1, data$cats-.5, NA)
  data$shadeymax <- ifelse(data$group == 1 & data[,sig.varname]==1, data$cats+.5,NA)
  data$shadexmin <- ifelse(data$group == 1 & data[,sig.varname]==1, -Inf, NA)
  data$shadexmax <- ifelse(data$group == 1 & data[,sig.varname]==1, Inf, NA)
  
  ## caterpillar dotplot
  p <- ggplot(data, aes(y=cats)) +  # get(target.varname)
    xlab(expression(target.varname)) + ylab("Policy category")  +
    scale_x_continuous(breaks = seq(xlim[1],xlim[2],by=1), limits=xlim  )  +
    scale_y_continuous(breaks=1:56,labels=data$catlab[data$group==1])  +
    geom_blank()  +
    theme(legend.position="none") +
    geom_rect(aes(ymin=shadeymin,ymax=shadeymax, xmin=shadexmin, xmax=shadexmax),
              fill="grey90") +
    geom_vline(xintercept=0, linetype = "dashed") +
    geom_point(aes(x=get(target.varname),shape=hiexpos,colour=post),size=1.3) +  # + 
    scale_shape_manual(values = c(24,19)) + 
    scale_colour_manual(values = c("black","red"))
  ##   specified shapes: solid ones for social and economic, hollow for others
  
  
  p <- p  +  theme(panel.grid.major.y = element_line(color="grey", linetype="dotted"))  +
    theme(panel.grid.minor.y = element_blank()) +
    theme(plot.background = element_blank()) +
    theme(axis.text=element_text(size=7) ) +
    theme(axis.title=element_text(size=9, face="bold") ) + theme_bw() # axis text size to 7 from 14
  if(pdf == F) {
    return(p)
  } else {
    ggsave(filename=paste("figures/", plotname, ".pdf", sep=""), height=12, width=8)
  }
}

# Function for labeling items for CAP data
labelbetas <- function(beta,mydocitem,mydesc,myremoved){
  require(plyr)
  #tmp <- colnames(mydocitem)[-which(apply(mydocitem,2,sum)==0)] # remove all-zero categories from names
  tmp <- colnames(mydocitem)
  # remove chars from colnames of mydocitem
  for (i in 1:length(myremoved)){
    tmp <- sub(myremoved[i],"",tmp)
  }
  tmp2 <- data.frame(topic.str=tmp)
  itemdesc <- join(tmp2,mydesc,by="topic.str", type="left") # use join since it preserves the order!!!
  ret <- cbind(values=beta,desc=itemdesc)
  ret <- ret[order(ret$values),]
  return(ret[,c(3,4,1)])
}




plotBetas2d <- function(x,dim=1,select,groups,grouplabels, labels=cmpcatlab, spec="spec") {
  
  if (dim==1)  {
    betas <- x$beta1 
    betas.se <- x$beta1.se 
  } else { 
    betas <- x$beta2
    betas.se <- x$beta2.se 
  }
  
  zeros <- which(round(betas,5)==0)  
  groups <- factor(groups,labels=grouplabels) # use NA, 1 = pure, 2 = mixed
  
  par(mfrow=c(1,1), mar=c(4, 4, 2, 2)+.1)
  pdf(file=paste("figures/2dBetaXrile",spec,"dimension", dim, "pdf", sep="."), height=12, width=8)
  #pdf(file=paste("plotname", "pdf", sep="."), height=12, width=8)
  dotchart.ci(betas[-zeros], groups=groups[-zeros], labels=labels[-zeros],
              stderr=betas.se[-zeros],
              xlab=expression(hat(beta)),
              pch=20, xlim=c(-3,3),
              cex=1.3, cex.axis=1.3, cex.main=1.5, cex.lab= 1.3)
  abline(v=0, lty="dashed", col="red")
  dev.off()
}





corExpertsByCountry <- function(model,cmpData){

  tmp <- names(table(cmpData$countryname))
  thetas <- model$summary$mean[grepl("^theta",rownames(model$summary))]
  out <- data.frame(NULL)
  for(i in 1:length(tmp)) {
    out[i,1] <- tmp[i]
    out[i,2] <- round(cor(cmpData$rile[cmpData$countryname==tmp[i]],cmpData$experts[cmpData$countryname==tmp[i]],use="pairw"),2)
    out[i,3] <- round(cor(thetas[cmpData$countryname==tmp[i]],cmpData$experts[cmpData$countryname==tmp[i]],use="pairw"),2)
    out[i,4] <- sum(complete.cases(cmpData$experts[cmpData$countryname==tmp[i]]))
  }
  colnames(out) <- c("countryname","rile-experts", "theta-experts", "N")
  out <- merge(out, unique(cmpData[,c("countryname","country")]))
  return(out)
}


processCAP <- function(stanobj, data){
    
    print("Summary of convergence diagnostics Rhat:")
    print(summary(stanobj$summary[,"Rhat"]) )
    hirhat <- rownames(stanobj$summary)[stanobj$summary[,"Rhat"] > 1.03]
    
    dfl <- data.frame(lambda=stanobj$summary[grep("^beta",rownames(stanobj$summary)),1],
                      cilo=stanobj$summary[grep("^beta",rownames(stanobj$summary)),2],
                      cihi=stanobj$summary[grep("^beta",rownames(stanobj$summary)),3],
                      code=sub("count.","",colnames(data)[-(1:2)],fixed=T)
    )
    dfl <- merge(dfl, capcat, all.x=TRUE)
    print("'significant' discrimination params:")
    print(table(sign(dfl$cilo) == sign(dfl$cihi) ))
    print(dfl$topic[sign(dfl$cilo) == sign(dfl$cihi) ])
    
    dft <- data.frame(theta=stanobj$summary[grep("^theta",rownames(stanobj$summary)),1],
                      cilo=stanobj$summary[grep("^theta",rownames(stanobj$summary)),2],
                      cihi=stanobj$summary[grep("^theta",rownames(stanobj$summary)),3],
                      party = data[,1], year= data[,2]
    )
    return(list(lambda=dfl, theta=dft, hirhat=hirhat))
}
# prints summary and returns data.frames with lambdas and thetas


irfplotpoisson <- function(alpha,zeta,lambdarange,thetarange,
                           lambdastep=.2,thetastep=.1, fivelayout=TRUE,
                           main = "") {
  thetaseq <- seq(thetarange[1],thetarange[2],thetastep)
  lambdaseq <- seq(lambdarange[1],lambdarange[2],lambdastep)
  ymax <- exp(alpha + zeta + lambdarange[2]*thetarange[2])
  if (fivelayout==TRUE) {
    fivelty <- c(1,2,4,1,2)
    # fivecol <- c(rep("blue",2),"black",rep("red",2))
    fivecol <- c(rep("red",2),"black",rep("blue",2))
  }
  plot(NULL, type="n",  xlim=thetarange, ylim=c(0,ymax),
       xlab=bquote(paste("Latent position ",theta)), ylab="Expected number of items",
       main = main)
  for (i in 1:length(lambdaseq)) {
    y <- exp(alpha + zeta + lambdaseq[i]*thetaseq)
    par(new=T)
    plot(thetaseq,y, type="l", xlab="",ylab="", xlim=thetarange,ylim=c(0,ymax), axes=F, 
         lty=ifelse(fivelayout,fivelty[i],1), col=ifelse(fivelayout,fivecol[i],"black") )
  }
}

irfplotmultinomial <- function(zeta,lambdarange,thetarange,
                               lambdastep=.2,thetastep=.1, fivelayout=TRUE,
                               main = "") {
  thetaseq <- seq(thetarange[1],thetarange[2],thetastep)
  lambdaseq <- seq(lambdarange[1],lambdarange[2],lambdastep)
  exbtable <- matrix(nrow=length(lambdaseq),ncol = length(thetaseq))
  for (i in 1:length(lambdaseq)) {
    exbtable[i,] <- exp(zeta + lambdaseq[i]*thetaseq)
  }  
  if (fivelayout==TRUE) {
    fivelty <- c(1,2,4,1,2)
    # fivecol <- c(rep("blue",2),"black",rep("red",2))
    fivecol <- c(rep("red",2),"black",rep("blue",2))
  }
  plot(NULL, type="n",  xlim=thetarange, ylim=c(0,1),
       xlab=bquote(paste("Latent position ",theta)), ylab="Expected probability/proportion",
       main = main)
  for (i in 1:length(lambdaseq)) {
    y <- exbtable[i,]/colSums(exbtable)
    par(new=T)
    plot(thetaseq,y, type="l", xlab="",ylab="", xlim=thetarange,ylim=c(0,1), axes=F,
         lty=ifelse(fivelayout,fivelty[i],1), col=ifelse(fivelayout,fivecol[i],"black") )
  }
}
