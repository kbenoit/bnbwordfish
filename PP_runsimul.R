rm(list=ls())  ## clear memory of all objects
gc()

# set working dir to PP_replication folder



######################
### read functions ###
######################

source("HMC_Stan_code_PP.R")
source("HMC_R_functions_PP.R")

# library(quanteda)
# library(foreign)

######################
### read data      ###
######################

load("Data/PP_data.RData")


###############
### Models ####
###############

rstan_options(auto_write=TRUE)


### 1d, Western countries, long time period (post 1972)

m.1d.welo <- wordfishSTAN1d(data=cmpNwelo, stan.code=stanNegbin1dWord 
                               , negbin=T, psi.ref=56,  dir.terms=c(51,20),
                               init.type="random", 
                            iter=4000, warmup=1000, n.thin=3, n.chains=2, 
                               monitor.alpha=F, monitor.aux=F )

### 1d model: post 1990, democracies ###

m.1d.po90 <- wordfishSTAN1d(data=cmpNpo90, stan.code=stanNegbin1dWord 
                            , negbin=T, psi.ref=56,  dir.terms=c(51,20),
                            init.type="random", 
                            iter=4000, warmup=1000, n.thin=3, n.chains=2, 
                            monitor.alpha=T, monitor.aux=F )
    
### 2d model: post 1990, democracies ###

m.2d.po90 <- wordfishSTAN2d(data=cmpNpo90,  
                            base.terms=c(20,45), # free enterprise + // trad moral +
                            bothzero.terms=18, # pol. corruption
                            dim1only.terms=c(21:35,39:42,51,52,54),
                            dim2only.terms=c(11,12,46:48,55,56),
                            # implies loading on both: 43:44 nat way of life, 49:50 multicult and
                            # also: 1:10,13:17,19,36 environment 37 culture +,38 soc justice +, 53 farmers
                            init.type="random", 
                            iter=4000, warmup=1000, n.thin=3, n.chains=2,  
                            monitor.alpha=F, monitor.aux = T )

save(m.1d.welo, m.1d.po90, m.2d.po90, file="resultsResAndPol.RData")

###########
### CAP ###
###########

load("Data/PP_CAPdata.RData")

m.1d.cap.bel <- wordfishSTAN1d(data=dcbe[,-(1:2)] ,stan.code=stanNegbin1dWord # reuse=m.nb.word.cap.bel$stanfit 
                                    , negbin=T, psi.ref=1, 
                               dir.terms=c(grep("2103",colnames(dcbe))-2, # Public lands/natural resources
                                           grep("2011",colnames(dcbe))-2), # Branch relations
                                    init.type="random",
                               iter=10000, warmup=5000, n.thin=5, n.chains=2, 
                                    monitor.alpha=F, monitor.aux=T, return.raw.stan=T )

m.1d.cap.den <- wordfishSTAN1d(data=dcdk[,-(1:2)] , reuse=m.1d.cap.bel$stanfit ,
                               #stan.code=stanNegbin1dWord ,
                               negbin=T, psi.ref=1, 
                               dir.terms=c(26,29), 
                               dir.beta=FALSE,
                               init.type="random",
                               iter=10000, warmup=5000, n.thin=5, n.chains=2, 
                               monitor.alpha=F, monitor.aux=T, return.raw.stan=F )

save(m.1d.cap.bel, m.1d.cap.den, 
     file="resultsCAPPP.RData")


