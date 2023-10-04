library(readstata13)
library(dplyr)
library(countrycode)

# set working dir to PP_replication folder



################
### CMP data ###
################
  


load(file="Data/CMP2018b.RData")
d <- d[d$date > 197212,]
d <- d[!(d$progtype %in% c(3,5,99)),]  # drop estimated, averaged and missing manifesto cases 
d <- d[!is.na(d$total),]  # drop the ones where total==0 
d$year <- floor(d$date/100)

dq <- read.dta13("Data/qog_tomerge_2019.dta", convert.factors=FALSE )
d <- merge(d,dq, by=c("countryname","year"), all.x=TRUE)
# View(d[is.na(d$ourdemocracy),c("countryname","year")])
# ok, missing only for Northern Ireland 1973 and post-communist cases in 1990 or early 1990s
d <- d[d$countryname != "Sri Lanka",] # since only two cases anyway
d <- d[d$countryname != "Malta",] # since only four cases in total anyway

## create variables for sub-setting cases for analysis
dagg <- summarize(group_by(d,countryname),
                  demshare=mean(ourdemocracy,na.rm=T), # share of democratic cases (note that this is affected by intra-country variation in no of parties)
                  mindemyear = min(ifelse(!is.na(ourdemocracy),year,NA),na.rm=T), # earliest known democratic election in this data
                  maxdemyear = max(ifelse(!is.na(ourdemocracy),year,NA),na.rm=T) # latest "
                  )
d <- merge(d,dagg)

# Western countries observed for a long period, from 1970s to now
d$westlong <- 0
d$westlong[(d$demshare > .9 & d$mindemyear < 1980) & d$ourdemocracy == 1] <- 1
# this includes the following three cases in general, which have one election in a non-democratic year each
# View(unique(d[d$countryname %in% c("Greece","Spain","Portugal"),c("countryname","date","ourdemocracy")]))
# (leaves out the one non-dem election in those)
table(d$westlong)

# post 1990 democracies
d$post90dem <- 0
d$post90dem[d$date > 199012 & d$demshare > .4 & !(d$countryname %in% c("Macedonia","Serbia")) & d$ourdemocracy == 1] <- 1
table(d$post90dem)
table(d$countryname[d$post90dem==1])
# Macedonia and Serbia not continuously democratic after onset of dem. -> take out
# includes South africa and Korea as more unusual cases
# View(unique(d[d$countryname %in% c("South Africa","South Korea","Croatia","Montenegro","Macedonia","Serbia")
# ,c("countryname","date","ourdemocracy")]))
d$groupto90 <- ifelse(d$date <= 199012, 1, 2)

d$cee <- ifelse(d$country > 74 & d$country < 100,1,0)

source("cmpcatlab.R") # loads 56 category labels 

### prepare counts ###

# NOTE duplicate documents linked to several parties may remain
corelo <- which(colnames(d)=="per101")
corehi <- which(colnames(d)=="per706")
stopifnot(corehi-corelo==55)

#  core 56 categories
stopifnot(colSums(d[,corelo:corehi])!=0)
#  post-communist 4-digit cats
fourlo <- which(colnames(d)=="per1011")
fourhi <- which(colnames(d)=="per7062")
stopifnot(fourhi-fourlo==53)
# there are a handful of NAs in those in two Russian and one Slovak party, replace with 0:
d[,fourlo:fourhi][is.na(d[,fourlo:fourhi])] <- 0
nonzeros <- colSums(d[,corelo:fourhi])!=0

# set the data to integer counts of category occurences - 
cmp.n.core <- round(d[, c(corelo:corehi)]/100 * d$total)   # core 56 categories
# core + extended categories (nonzero ones)
cmp.n.all  <- round(d[, c(corelo:fourhi)][,nonzeros]/100 * d$total)  # NOTE: we do not consider 'peruncod'

## combine post-comm cats into core
cmp.n.combined <- as.data.frame(t(cmp.n.all))
# create a variable with the "base" cmp category (without the 4th digit)
cmp.n.combined$cat <- substr(rownames(cmp.n.combined), 1, 6)
# aggregate the counts by the cat variable
cmp.n.combined <- aggregate(cmp.n.combined[,1:(ncol(cmp.n.combined)-1)], 
                            by=list(cmp.n.combined$cat), sum)
# transpose it back to the original dimensions, without the cat (aggregator) column
cmp.n.combined <- as.data.frame(t(cmp.n.combined[, 2:ncol(cmp.n.combined)]))
# rename the vars
names(cmp.n.combined) <- names(cmp.n.all)[c(1:56)]
cmpN <- cmp.n.combined
rm(cmp.n.combined)
cmpcatno <- as.numeric(substr(colnames(cmpN),4,7))

# inspect if shares add to 1
# NOTE there are some inconsistencies
#d$peruncod <- ifelse(is.na(d$peruncod),0,d$peruncod) # this sets missings in peruncod to 0
d$test <- (rowSums(cmpN) + d$total*d$peruncod/100)/d$total
d$testwouncod <- rowSums(cmpN)/d$total
summary(d$test)
# View(d[is.na(d$test) | (!is.na(d$test) & (d$test < .98 | d$test > 1.02)), 
#        c("countryname","partyname","date","testwouncod","test","peruncod")])

# counts for west-long set of cases:
cmpNwelo <- cmpN[d$westlong==1,]
stopifnot(colSums(cmpNwelo)!=0) # make sure no all zero categories
# counts for post 1990 set of cases:
cmpNpo90 <- cmpN[d$post90dem==1,]
stopifnot(colSums(cmpNpo90)!=0)

save(d, cmpNwelo, cmpNpo90,  file = "Data/PP_data.RData")

################
### CAP data ###
################

library(stringr)

### Master codebook ###
cap <- readLines("Data/CAP/CAPmasterCodebook.txt")
cap <- str_trim(cap)
cap <- cap[!grepl("^De",cap)] # Description entries (some have typos in that word)
cap <- str_split_fixed(cap, pattern="[\\p{P}\\p{S}]",2)
cap <- data.frame(cap, stringsAsFactors = FALSE)
names(cap) <- c("code","topic")
cap$code <- as.numeric(cap$code)
cap$topic <- str_trim(cap$topic)
cap$tmp <- ifelse(cap$topic %in% c("General","Other"), cap$code, 999999)
cap$tmp[cap$topic=="Other"] <- cap$tmp[cap$topic=="Other"]-99
cap$tmp <- cap$tmp/100
cap <- merge(cap,cap,by.x="tmp", by.y="code", all.x=TRUE)

cap$topic.x[!is.na(cap$topic.y)] <- paste(cap$topic.y[!is.na(cap$topic.y)], cap$topic.x[!is.na(cap$topic.y)])
cap <- cap[,c("code","topic.x")]
names(cap) <- c("code","topic")

### BEL ###
dcbe <- read.csv(file="Data/CAP/belgium_partymanifestos_3.csv")
summary(dcbe)
dcbe <- dcbe[dcbe$year > 1990,]
# add category labels:
dcbe <- merge(dcbe, cap, by.x="subtopic", by.y="code", all.x=TRUE)
# the 'topic' var has the label
table(dcbe$topic,useNA="ifany") # there are a few sentencs with codes lacking a correspondence in the master codebook
table(dcbe$subtopic[is.na(dcbe$topic)])
#View(dcbe[is.na(dcbe$topic),])
#leave them in

dcbe <- summarize(group_by(dcbe,party,year,subtopic),
                  count=n())

# dcbe <- mutate(group_by(dcbe,party,year),
#                 total=sum(count))
dcbe <- ungroup(dcbe)

dcbe <- reshape(as.data.frame(dcbe), idvar=c("party","year"), timevar="subtopic",
                  direction="wide")
# sort colnames 
dcbe <- dcbe[,sort(colnames(dcbe))]
dcbe <- dcbe[,c(ncol(dcbe)-1,ncol(dcbe),1:(ncol(dcbe)-2))]
dcbe[is.na(dcbe)] <- 0 # fill in missings due to reshape from long format
table(dcbe$year)


### DEN ###

dcdk1 <- read.csv(file="Data/CAP/denmark_party_manifestoes_1953-2011_part1.txt",
                 stringsAsFactors = FALSE)
dcdk2 <- readLines("Data/CAP/denmark_party_manifestoes_1953-2011_part2.txt")
# one solution: pull out the entries up to the 4th comma
tmp <- str_locate_all(dcdk2,",")
tmp2 <- plyr::laply(tmp,function(x) x[4,1])
tmp3 <- str_sub(dcdk2,1,tmp2-1)
writeLines(tmp3,"tmplines.txt")
dcdk2 <- read.csv("tmplines.txt", sep=",", header=FALSE)
file.remove("tmplines.txt")
names(dcdk2) <- names(dcdk1)[1:4]
dcdk <- rbind(dcdk1[,1:4], dcdk2)
summary(dcdk)
dcdk <- dcdk[dcdk$year > 1990,]
dcdk$party <- as.numeric(str_sub(dcdk$id,5,6))
table(dcdk$party[dcdk$year < 2011])
table(dcdk1$party[dcdk1$year > 1990]) # ok, same as prev. line
rm(dcdk1,dcdk2)

# add category labels:
dcdk <- merge(dcdk, cap, by.x="subtopic", by.y="code", all.x=TRUE)
# the 'topic' var has the ladkl
table(dcdk$topic,useNA="ifany") # there are ca. 6% of sentencs with codes lacking a correspondence in the master codebook
table(dcdk$subtopic[is.na(dcdk$topic)]) # all have code 2999
#View(dcdk[is.na(dcdk$topic),])
# leave them in

dcdk <- summarize(group_by(dcdk,party,year,subtopic),
                  count=n())

# dcdk <- mutate(group_by(dcdk,party,year),
#                 total=sum(count))
dcdk <- ungroup(dcdk)
dcdk <- reshape(as.data.frame(dcdk), idvar=c("party","year"), timevar="subtopic",
                direction="wide")
dcdk <- dcdk[,sort(colnames(dcdk))]
dcdk <- dcdk[,c(ncol(dcdk)-1,ncol(dcdk),1:(ncol(dcdk)-2))]

dcdk[is.na(dcdk)] <- 0 # fill in missings due to reshape from long format
table(dcdk$year)
summary(rowSums(dcdk[,3:ncol(dcdk)]))
dcdk <- dcdk[rowSums(dcdk[,3:ncol(dcdk)]) > 20,] # remove one manifesto with just 2 sentences 

capcat <- cap
save(capcat, dcbe,dcdk, file = "Data/PP_CAPdata.RData")

