rm(list=ls())  ## clear memory of all objects
gc()

# set working dir to PP_replication folder
if (Sys.getenv("USER")=="kbenoit") {
    setwd("~/Dropbox/Papers/manifesto_length/analysis/PPreplication")
} else if (Sys.getenv("USERNAME")=="Thomas") {
    setwd("C:/Users/Thomas/Dropbox/manifesto_length/analysis/PPreplication")
}

library(dplyr)

######################
### read functions ###
######################

source("HMC_R_functions_PP.R")

#####################
### read CMP data ###
#####################

load("Data/PP_data.RData")
source("cmpcatlab.R") # loads 56 category labels 
# category number perXXX:
cmpcatno <- as.numeric(substr(colnames(cmpNpo90),4,7))

## add more recent expert data 

dches14 <- readstata13::read.dta13("Data/2014_CHES_dataset_means.dta", convert.factors = FALSE)
dches14 <- dches14[,c("party_name","party_id","lrgen")] 
dches14$expyear <- 2014
dches17 <- readstata13::read.dta13("Data/CHES_means_2017.dta", convert.factors = FALSE) # (note that this hasn't got BEL or DEN)
dches17 <- dches17[,c("party_name","party_id","lrgen")] 
dches17$expyear <- 2017
dches1417 <- bind_rows(dches14, dches17)
rm(dches14, dches17)
dches1417 <- reshape(dches1417[,2:4], idvar="party_id", timevar="expyear", direction="wide")

# use info from party-facts
dpf <- readstata13::read.dta13("Data/partyfacts-external-parties.dta", convert.factors = FALSE)
dpf <- dpf[, c("partyfacts_id","dataset_key","dataset_party_id","year_first","year_last")]
dpf <- merge(dpf[dpf$dataset_key == "manifesto",],
             dpf[dpf$dataset_key == "ches",],
             by="partyfacts_id" )
dpf$dataset_party_id.x <- as.numeric(dpf$dataset_party_id.x) # manifesto
dpf$dataset_party_id.y <- as.numeric(dpf$dataset_party_id.y) # ches

dpf <- unique(dpf[,c(3:5,7:9)])
dpf <- dpf[dpf$year_last.y >= 2011,] # keep only cases of interest for CHES 2014,17 (three year window, see below)
dpf <- dpf[order(dpf$dataset_party_id.x),]
# fix multiply occurring cmp party_ids:
dpf <- add_count(dpf,dataset_party_id.x)
#View(dpf[dpf$n > 1,])
dpf <- dpf[!((dpf$dataset_party_id.x == 31421 & dpf$dataset_party_id.y == 603 ) |
               (dpf$dataset_party_id.x == 32230 & dpf$dataset_party_id.y == 850 ) |
               (dpf$dataset_party_id.x %in% c(32640, 32903, 32952) & dpf$dataset_party_id.y %in% c(851,852)) |
               dpf$dataset_party_id.x == 89031   ),]
dpf <- select(dpf,-n)
dpf <- add_count(dpf,dataset_party_id.x)
# drop 31421 603 ; 32230 850 (since 850 refers to 2017 which at least currently not included)
# c(32640, 32903, 32952) & 851 or 852, ches is 2017 too late to be relevant
# 89031 Macedonian party not of interest
stopifnot(dpf$n==1)


dpf <- merge(dpf,dches1417, by.x="dataset_party_id.y", by.y="party_id")
dpf <- dpf[,c("dataset_party_id.x","lrgen.2014","lrgen.2017")]
names(dpf)[1] <- "party"

# from older file
dex <- readstata13::read.dta13("Data/MPDataset_MPDS2014b_plusexperts.dta", convert.factors = FALSE) # this is based on: closest expert survey within 3 years pre/post the election
dex <- dex[,c("party","date","Posrile","ches_lrgen")]
x1 <- 1
y1 <- 0
x2 <- 20
y2 <- 10
a <- (y1 - y2)/(x1 - x2)
b <- y1 - x1*(y1 - y2)/(x1 - x2)
dex$experts <- a*dex$Posrile + b
dex$experts[which(is.na(dex$Posrile))] <- dex$ches_lrgen[which(is.na(dex$Posrile))]
dex <- dex[!is.na(dex$experts),]

d <- dplyr::left_join(d,dex,by=c("party","date"))

d <- dplyr::left_join(d,dpf,by=c("party"))
# note that up to here the lrgen. variables are just added, not incorporated into experts variable
# one way of doing this:
# # use CHES 2014 if within 3 years and closer than the 2010 value if that is not NA
# 
# 
d$experts_orig <- d$experts
# take 2014 value if closer than 2010 value
d$experts[!is.na(d$experts) & (abs(d$edate-as.Date("2014-07-02","%Y-%m-%d")) <=
                                 abs(d$edate-as.Date("2010-07-02","%Y-%m-%d")))] <-
  d$lrgen.2014[!is.na(d$experts) & (abs(d$edate-as.Date("2014-07-02","%Y-%m-%d")) <=
                                      abs(d$edate-as.Date("2010-07-02","%Y-%m-%d")))]
# take 2014 value if misssing beofre and withiin 3 years
d$experts[is.na(d$experts) & (abs(d$edate-as.Date("2014-07-02","%Y-%m-%d")) <= (3*365.25))] <-
  d$lrgen.2014[is.na(d$experts) & (abs(d$edate-as.Date("2014-07-02","%Y-%m-%d")) <= (3*365.25))]

# use CHES 2017 if within 3 years and closer than the 2014 value
d$experts[(abs(d$edate-as.Date("2014-07-02","%Y-%m-%d")) <= (3*365.25)) &
            (abs(d$edate-as.Date("2014-07-02","%Y-%m-%d")) <=
               abs(d$edate-as.Date("2010-07-02","%Y-%m-%d")))] <-
  d$lrgen.2014[(abs(d$edate-as.Date("2014-07-02","%Y-%m-%d")) <= (3*365.25)) &
                 (abs(d$edate-as.Date("2014-07-02","%Y-%m-%d")) <=
                    abs(d$edate-as.Date("2010-07-02","%Y-%m-%d"))) ]

################
### CAP data ###
################

load("Data/PP_CAPdata.RData")
load("resultsCAPPP.RData") 


cap.be <- processCAP(m.1d.cap.bel, dcbe)
cap.be$hirhat # R hat diag > 1.03

compare2chainsplot(m.1d.cap.bel,coda=T,alpha=F,star=F)



# CAP data use only current party names.
# need to add (partly time-varying) cmp code:
dbenam <- data.frame(
  pnam=c("CDH", "CDV", "ECOLO", "FDF", "FN", "GROEN",
         "MR", "NVA", "PS", "RW", "SPA", "SPIRIT", "VB", "VIVANT", "VLD"),
  cmpcode=c(21522,21521,21111,21912,NA,21112,
            21426,21916,21322,21911,21321,21330,21914,NA,21421)
  
)
cap.be$theta <- left_join(cap.be$theta,dbenam,by=c("party" = "pnam"))
cap.be$theta$cmpcode[cap.be$theta$party == "NVA" & cap.be$theta$year %in% c(1991,1995)] <- 21913 # VU
cap.be$theta$cmpcode[cap.be$theta$party == "NVA" & cap.be$theta$year ==1999] <- 21915 # VU-ID
cap.be$theta$cmpcode[cap.be$theta$party == "VB" & cap.be$theta$year >=2007] <- 21917 # VB
cap.be$theta$cmpcode[cap.be$theta$party == "SPA" & cap.be$theta$year %in% c(2003,2007)] <- 21221 # sp.a-spirit
cap.be$theta$cmpcode[cap.be$theta$party == "MR" & cap.be$theta$year ==1991] <- 21422 # prl
cap.be$theta$cmpcode[cap.be$theta$party == "MR" & cap.be$theta$year ==1995] <- 21423 # prl-fdf
cap.be$theta$cmpcode[cap.be$theta$party == "MR" & cap.be$theta$year ==1999] <- 21425 # prl-fdf-mcc
# cap fdf 03 no clear correspondence in cmp

cap.be$theta <- left_join(cap.be$theta, d[,c("party","year","rile","experts","lrgen.2014","lrgen.2017","parfam")],
                          by=c("cmpcode"="party", "year"="year"))

summary(cap.be$theta$experts)
# add missing values due to the fact that these cases are not in CMP and thus not in above expert data
cap.be$theta$experts[cap.be$theta$cmpcode == 21914 & cap.be$theta$year == 1999 ] <- 9.89 # CHES Vlaams Blok
cap.be$theta$experts[cap.be$theta$cmpcode == 21914 & cap.be$theta$year == 2003 ] <- a*18.86+b # PPMD Vlaams Blok
cap.be$theta$experts[cap.be$theta$party == "SPA" & cap.be$theta$year == 2007 ] <- 3.22 # CHES, 2006:"SPA", raw CAP 2007: "SPA"
summary(cap.be$theta$experts)

cor(cap.be$theta$theta,cap.be$theta$rile,use="pair" )
cor(cap.be$theta$theta,cap.be$theta$experts,use="pair" )
cor(cap.be$theta$experts,cap.be$theta$rile,use="pair" )



cap.dk <- processCAP(m.1d.cap.den, dcdk)
cap.dk$hirhat # R hat diag > 1.03

compare2chainsplot(m.1d.cap.den,coda=T,alpha=F,star=F)

ddknam <- data.frame(
  pnum=c(1:4,6,8,10,11,12,14,16,20),
  cmpcode=c(13320,13410,13620,13330, # soc dem, soc lib, cons, centre dem
            13230,13720, # social people, danish people
            13520,13420,13229, # chr-dem, lib, red-green
            13951,13001,NA # progr, lib all, minority
  ) )
cap.dk$theta <- left_join(cap.dk$theta,ddknam,by=c("party" = "pnum"))
cap.dk$theta <- left_join(cap.dk$theta, d[,c("party","year","rile","experts","lrgen.2014","lrgen.2017","parfam")],
                          by=c("cmpcode"="party", "year"="year"))
summary(cap.dk$theta$experts)

cap.dk$theta$experts[cap.dk$theta$party == 4 & cap.dk$theta$year == 2001 ] <- 5.57 # center dem CHES, 1999
# cap.dk$theta$experts[cap.dk$theta$party == 16 & cap.dk$theta$year == 2007 ] <- 8.2 # lib alliance CHES, 2010 don't add it since the 2010 CHES value is used for 2011
cap.dk$theta$experts[cap.dk$theta$party == 1 & cap.dk$theta$year == 2005 ] <- 4.11 # 05 early ele not in that cmp data CHES, 2006
cap.dk$theta$experts[cap.dk$theta$party == 2 & cap.dk$theta$year == 2005 ] <- 4.78 # 05 early ele not in that cmp data CHES, 2006
cap.dk$theta$experts[cap.dk$theta$party == 3 & cap.dk$theta$year == 2005 ] <- 7.11 # 05 early ele not in that cmp data CHES, 2006
# for parties  4 , 10, 20 for this ele no data in CHES
cap.dk$theta$experts[cap.dk$theta$party == 6 & cap.dk$theta$year == 2005 ] <- 2.33  # 05 early ele not in that cmp data CHES, 2006
cap.dk$theta$experts[cap.dk$theta$party == 8 & cap.dk$theta$year == 2005 ] <- 7.67  # 05 early ele not in that cmp data CHES, 2006
cap.dk$theta$experts[cap.dk$theta$party == 11 & cap.dk$theta$year == 2005 ] <- 7.22  # 05 early ele not in that cmp data CHES, 2006
cap.dk$theta$experts[cap.dk$theta$party == 12 & cap.dk$theta$year == 2005 ] <- 1  # 05 early ele not in that cmp data CHES, 2006
summary(cap.dk$theta$experts)

cor(cap.dk$theta$theta,cap.dk$theta$rile,use="pair" )
cor(cap.dk$theta$theta,cap.dk$theta$experts,use="pair" )
cor(cap.dk$theta$experts,cap.dk$theta$rile,use="pair" )



#####################
### load HMC runs ###
#####################

load("resultsPP.RData") 


##########################
### convergence checks ###
##########################

summary(m.1d.welo$summary[,"Rhat"])

compare2chainsplot(m.1d.welo,coda=T,alpha=F,star=F)


summary(m.1d.po90$summary[,"Rhat"])
compare2chainsplot(m.1d.po90,coda=T,alpha=T,star=F)

summary(m.2d.po90$summary[,"Rhat"]) # the NaN are the betas constrained to zero
View(m.2d.po90$summary[!is.na(m.2d.po90$summary[,"Rhat"]) & m.2d.po90$summary[,"Rhat"] > 1.05,])

compare2chainsplot(m.2d.po90,coda=T,alpha=F,star=F,dim2=T)

#######################################################
#### Figure 1: Expert Corr with Theta vs with Rile    ###
#######################################################


# Pooled correlations
sum(d$westlong & !is.na(d$experts))
cor(m.1d.welo$summary[grep("^theta",rownames(m.1d.welo$summary)),1], d$experts[d$westlong ==1], use="pair")
cor( d$rile[d$westlong ==1], d$experts[d$westlong ==1], use="pair")

sum(d$post90dem & !is.na(d$experts))
cor(m.1d.po90$summary[grep("^theta",rownames(m.1d.po90$summary)),1], d$experts[d$post90dem ==1], use="pair")
cor(d$rile[d$post90dem ==1], d$experts[d$post90dem ==1], use="pair")

# table of expert corrs by country 
# mytab <- corExpertsByCountry(m.1d.welo, d[d$westlong==1,])
mytab <- corExpertsByCountry(m.1d.po90, d[d$post90dem==1,])

mytab <- subset(mytab, N >= 4 , ) # & countryname != "Romania"

mytab$region <- ifelse(mytab$country < 80, 1, 3)
mytab$region[mytab$country > 59 & mytab$country < 80] <- 2
mytab$region <- factor(mytab$region)
levels(mytab$region) <- c("Western Europe", "Pacific/East", "Eastern Europe" )

mytab <- mytab[,c(1:4,6)]

names(mytab)[2:3] <- c("Rile", "theta")

#mytab$region <- relevel(mytab$region, 3)
mytab <- mytab[order(mytab$region, mytab$theta), ]

#pdf("figures/Fig1ValidVsRile.pdf", height = 8.5, width = 8)
setEPS(height = 8.5, width = 8) # how to set the resolution (min 800 dpi) here, or is it implied???
postscript("figures/Fig1ValidVsRile.eps")
dotchart(mytab[, 3], mytab$country, groups = factor(mytab$region), xlim = c(-.3, 1.1), 
         #col = alpha("black", 0.95),
         col="grey20",
         pch=16, xlab = "Correlation with expert placements")
mytab$ylocations <- c(19:37, 14:16, 1:11)
mycol <- ifelse(mytab$theta > mytab$Rile, "red", "green")
mysym <- ifelse(mytab$theta > mytab$Rile, 15,17)
segments(mytab$theta, mytab$ylocations, mytab$Rile, mytab$ylocations, col = "grey50")
points(mytab$Rile, mytab$ylocations, 
       col = mycol, # alpha(mycol, .9),
       pch = mysym)
legend(.83, 5, c(expression(theta), "Rile (worse)","Rile (better)"), pch = c(16, 15,17),
       #col = c(alpha("black", 0.95), alpha("darkred", .9),"green"), bty = "o", bg = "grey95")
       col = c("grey20","red","green"))    
dev.off()



###############################################################
####### Figure 2: Discriminiation Caterpillar by Rile     ########
###############################################################

rile <- c(3,3,1,2,1,1,1,3,3,3,2,1,2,3,3,3,3,3,2,2,2,1,1,3,1,2,3,3,3,3,1,1,2,3,3,3,3,3,1,2,1,3,2,3,2,3,2,2,3,3,1,3,3,3,3,3)
rile <- factor(rile,labels=c("RILE left","RILE right","RILE neutral"))
#rile <- factor(rile,levels(rile)[c(2,1,3)])
rile <- factor(rile,levels(rile)[c(1,3,2)])
# grouping (pure econ vs rest)
gr <- rep(2,56)
gr[c(21:35,39:42,51,52,54,20)] <- 1 
gr <- factor(gr,labels=c("Purely economic","Other"))

dfl.welo <- data.frame(lambda=m.1d.welo$summary[grep("^beta",rownames(m.1d.welo$summary)),1],
                       cilo=m.1d.welo$summary[grep("^beta",rownames(m.1d.welo$summary)),2],
                       cihi=m.1d.welo$summary[grep("^beta",rownames(m.1d.welo$summary)),3],
                       rile=rile,
                       gr=gr,
                       catnum=1:56,
                       catlab=cmpcatlab,
                       stringsAsFactors = FALSE
)

ggPlotLambdasbyRile(dfl.welo,
                    pdf=F, plotname = "Fig2DiscrimByRile")
ggsave("figures/Fig2DiscrimByRile.eps",
       height=12, width=8, dpi=800)
# the warnings refer to the empty lines added for better legibility



#################################################
### FIGURE 3 : CAP BEL/DEN theta expert plot  ###
#################################################

cap.be$theta$Country <- "Belgium"
cap.dk$theta$Country <- "Denmark"

cap.bedk.theta <- bind_rows(select(cap.be$theta,-party), select(cap.dk$theta,-party))

table(cap.bedk.theta$parfam[!is.na(cap.bedk.theta$experts)],useNA="ifany")

cap.bedk.theta <- cap.bedk.theta[!is.na(cap.bedk.theta$experts),]
cap.bedk.theta$parfam2 <- cap.bedk.theta$parfam
cap.bedk.theta$parfam2[cap.bedk.theta$cmpcode == 13330] <- 110 # CMP has Center Democrats as SocDem, but rather use Other
cap.bedk.theta$parfam2[cap.bedk.theta$parfam %in% c(90, 95)] <- 110 # ethnic/regional and special issue = Other
cap.bedk.theta$parfam2[cap.bedk.theta$cmpcode %in% c(13720, 21914, 21917)] <- 100 # Dansk Folkeparti and Vlaams Belang = Far-Right
# VB expert entries are not in above data, since CMP2014b hasn't got those in sample


table(cap.bedk.theta$parfam2[!is.na(cap.bedk.theta$experts)],useNA="ifany")


partycolours <- c("darkgreen","purple","red","blue","orange",
                  "black","brown","darkgrey")
partylabels <- c("Green","Socialist","Social Dem.","Liberal","Christian Dem.",
                 "Conservative","Far Right","Other")

require(ggplot2)
pcap <- ggplot(cap.bedk.theta, aes(x=theta,y=experts,colour=as.factor(parfam2))) +
    geom_point(shape=19,size=4) +
    scale_colour_manual(values = partycolours,name="Party family",
                        labels=partylabels) +
    stat_smooth(method="lm", se=FALSE, colour="darkgrey", linetype="dashed",fullrange=T) +
    geom_errorbarh(aes(xmin=cilo, xmax=cihi), size=1, height=0 ) + # 
    scale_x_continuous(name=expression(paste("Left-right positions (",theta,") from scaling")),
                                       limits=c(-3.25,2.8)) +
                           scale_y_continuous(name="Expert survey estimates of left-right",limits=c(0,10)) +
                           facet_wrap(~ Country, ncol=2) +
                           theme_bw() + 
                           theme(legend.position = "bottom") +
                           theme(title = element_text(size = rel(2.3))) +
                           theme(axis.text = element_text(size = rel(2.3))) +
                           theme(legend.text = element_text(size = rel(1.9))) + 
                           theme(strip.text.x = element_text(size = 26))
                       
                       pcap +   coord_flip()
                       #ggsave(filename="figures/Fig3CAPexperts.pdf", height=10, width=22)
                       ggsave(filename="figures/Fig3CAPexperts.eps", height=10, width=22, dpi=800)
                       # the warning message appears to come from fullrange = T, everything seems fine
                       
                       sum(sign(cap.be$lambda$cilo) == sign(cap.be$lambda$cihi))
                       nrow(cap.be$lambda)
                       sum(sign(cap.be$lambda$cilo) == sign(cap.be$lambda$cihi))/nrow(cap.be$lambda)
                       
                       sum(sign(cap.dk$lambda$cilo) == sign(cap.dk$lambda$cihi))
                       nrow(cap.dk$lambda)
                       sum(sign(cap.dk$lambda$cilo) == sign(cap.dk$lambda$cihi))/nrow(cap.dk$lambda)
                       
                       #View(cap.be$lambda[sign(cap.be$lambda$cilo) == sign(cap.be$lambda$cihi), ])
                       View(cap.be$lambda)
                       View(cap.dk$lambda)
                       

##################################################
####### Figure 4: 2d discrimination scatter ######
##################################################

# groups of categories as per simulation call                           
gr2d <- rep(3,56)
gr2d[c(20, 21:35,39:42,51,52,54)] <- 1 # pure econ
gr2d[c(45,11,12,46:48,55,56)] <- 2
gr2d[18] <- NA
gr2d <- factor(gr2d,labels=c("Economic","Non-economic","Mixed"))


dfl.2d <- data.frame(lambda1=m.2d.po90$summary[grep("^beta1",rownames(m.2d.po90$summary)),1],
                      lam1cilo=m.2d.po90$summary[grep("^beta1",rownames(m.2d.po90$summary)),2],
                       lam1cihi=m.2d.po90$summary[grep("^beta1",rownames(m.2d.po90$summary)),3],
                     lambda2=m.2d.po90$summary[grep("^beta2",rownames(m.2d.po90$summary)),1],
                     lam2cilo=m.2d.po90$summary[grep("^beta2",rownames(m.2d.po90$summary)),2],
                     lam2cihi=m.2d.po90$summary[grep("^beta2",rownames(m.2d.po90$summary)),3],
                     gr=gr2d,
                     rile=rile,
                       catno=cmpcatno, 
                       catlab=cmpcatlab )

# original submission:
dfl.2d$catlab[(abs(dfl.2d$lambda1) < .5 & abs(dfl.2d$lambda2) < .5) ] <- NA # | dfl.2d$lambda1 == 0 | dfl.2d$lambda2 == 0 

cor(dfl.2d$lambda1, dfl.2d$lambda2)
cor(dfl.2d$lambda1, dfl.2d$lambda2,method="spearman")


library(ggrepel)
p <- ggplot(dfl.2d, aes(lambda1,lambda2,colour=rile, shape=rile)) +
  geom_point(size=3) +
  geom_text_repel(aes(label=catlab), size=7.5, show.legend=FALSE, force=5) + # label=catno
  geom_hline(aes(yintercept = 0), linetype="dashed", colour="grey60") +
  geom_vline(aes(xintercept = 0), linetype="dashed",colour="grey60") +
  xlim(c(-2.2,2.2)) + ylim(c(-1.5,2)) +
  scale_color_manual(values = c("red","black","blue")) +
  labs(x=expression(paste(beta," on economic dimension")), 
       y=expression(paste(beta," on socio-cultural dimension"))) +
  theme_bw() + 
  # theme(aspect.ratio=1 )  +
  theme(legend.position="bottom") + 
    theme(axis.title = element_text(size=rel(2.3))) +
    theme(axis.text = element_text(size = rel(2.3))) +
    theme(legend.text = element_text(size = rel(2.3))) + 
    guides(shape=guide_legend(title=NULL)) + guides(colour=guide_legend(title=NULL))
p
ggsave(p, file="figures/Fig4Scatter2d.pdf", width=17.5, height=14) 
ggsave(p, file="figures/Fig4Scatter2d.eps", width=17.5, height=14, dpi=800) 
# width 15 and height 12 at first submission (added force spec above and increased measures to fit all labels [beause of ggrepel version changes?])


######################################################
### Figure 5 (Appendix) Item Characteristic curves ###
######################################################

graphics.off()
pdf(file="figures/Fig5irfunctions.pdf", height=6, width=12)
layout(matrix(c(1,2), ncol=2),
       height=c(4,4),width=c(8,8),respect=FALSE)  
#
irfplotpoisson(1,1,c(-1,1),c(-3,3),lambdastep=.5, 
               main = bquote(paste("Expected category counts (",alpha==1," and all ",psi==1,")"))   )
text(-2.75,130,bquote(beta==-1),pos=4,col="red")
text(-2.75,30,bquote(beta==-.5),pos=4,col="red")
text(-2.75,10,bquote(beta==0),pos=4,col="black")
text(1.8,30,bquote(beta==.5),pos=4,,col="blue")
text(2,130,bquote(beta==1),pos=4,col="blue")

irfplotmultinomial(1,c(-1,1),c(-3,3),lambdastep=.5, 
                   main = bquote(paste("Item response category functions (any ",alpha,", all ",psi==1,")"))   )
text(-2.75,.8,bquote(beta==-1),pos=4,col="red")
text(-2.75,.28,bquote(beta==-.5),pos=4,col="red")
text(-2.75,.11,bquote(beta==0),pos=4,col="black")
text(2,.28,bquote(beta==.5),pos=4,col="blue")
text(2,.8,bquote(beta==1),pos=4,col="blue")

dev.off()

##################################################
####### Figure 6 (appendix): Phi ######
##################################################

dfp.welo <- data.frame(phi_inv=m.1d.welo$summary[grep("^phi",rownames(m.1d.welo$summary)),1],
                       cilo=m.1d.welo$summary[grep("^phi",rownames(m.1d.welo$summary)),2],
                       cihi=m.1d.welo$summary[grep("^phi",rownames(m.1d.welo$summary)),3],
                       propzero=apply(cmpNwelo, 2, function(x) mean(x==0)),
                       RILE=rile,
                       gr=gr,
                       catnum=1:56,
                       catlab=cmpcatlab,
                       stringsAsFactors = FALSE
)
levels(dfp.welo$RILE) <- c("left","neutral","right")
dfp.welo$propzerolab <- as.character(round(dfp.welo$propzero,2))
dfp.welo$propzerolab[nchar(dfp.welo$propzerolab)==3] <- paste0(dfp.welo$propzerolab[nchar(dfp.welo$propzerolab)==3],"0")
dfp.welo$propzerolab <- paste0("(",dfp.welo$propzerolab,")")
ggPlotPhi(dfp.welo, pdf=T, plotname = "Fig6Phi")

######################################################
####### Figure 7 (Appendix) 2d lambdas caterpillar ###
######################################################


dfl.2d <- data.frame(lambda1=m.2d.po90$summary[grep("^beta1",rownames(m.2d.po90$summary)),1],
                       cilo1=m.2d.po90$summary[grep("^beta1",rownames(m.2d.po90$summary)),2],
                       cihi1=m.2d.po90$summary[grep("^beta1",rownames(m.2d.po90$summary)),3],
                     lambda2=m.2d.po90$summary[grep("^beta2",rownames(m.2d.po90$summary)),1],
                     cilo2=m.2d.po90$summary[grep("^beta2",rownames(m.2d.po90$summary)),2],
                     cihi2=m.2d.po90$summary[grep("^beta2",rownames(m.2d.po90$summary)),3],
                       catnum=1:56,
                     rile=rile,
                       catlab=cmpcatlab,
                       stringsAsFactors = FALSE
)
dfl.2d$pure1 <- ifelse(dfl.2d$lambda2 == 0, 1, 0)
dfl.2d$pure2 <- ifelse(dfl.2d$lambda1 == 0, 1, 0)
dfl.2d <- dfl.2d[!(dfl.2d$lambda1 == 0 & dfl.2d$lambda2 == 0),]

cor(dfl.2d$lambda1, dfl.2d$lambda2)
cor(dfl.2d$lambda1[dfl.2d$rile=="RILE right"], dfl.2d$lambda2[dfl.2d$rile=="RILE right"])
cor(dfl.2d$lambda1[dfl.2d$rile=="RILE left"], dfl.2d$lambda2[dfl.2d$rile=="RILE left"])


ggPlotLambdas2D(dfl.2d,dim=1,
                    pdf=TRUE, plotname = "Fig7aDiscrim2dEcon")
# removed 3 rows are labels, ok.

ggPlotLambdas2D(dfl.2d,dim=2,
                pdf=TRUE, plotname = "Fig7bDiscrim2dSocial")
# removed 3 rows are labels, ok.




