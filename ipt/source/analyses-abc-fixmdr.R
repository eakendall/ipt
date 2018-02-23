## Script for running some of the paper analyses based on ABC

########################################################################
setwd("C:/Users/EAK/Google Drive/Latent TB models/khayelitsha")
source("source/khay-utility-fixmdr.R")
library(dplyr)
library(abind)
reload.source()

savelocation <- ""
date <- "20171213"
ntheta=20

## using output
set.seed(24543663)

draws <- lls.base <- steady.states <- lls.ipt <- thetas <- numeric()
for (x in c("a","bb"))#,"g","h")) # 50 each for 20171120
{
  out <- readRDS(paste0(savelocation,"generated_data/LHSdraws_and_lls.",date,x,".rds"))
  ## get the elements from the runs
  draws <- rbind(draws, out$draws[rep(1:nrow(out$draws), each=ntheta),])
  # lls.base <- append(lls.base, unlist(out$lls)[which(names(unlist(out$lls))=="ll")][rep(1:nrow(out$draws), each=ntheta)])
  lls.base <- rbind(lls.base, t(matrix(unlist(out$lls)[grep(names(unlist(out$lls)), pattern="^ll\\.")], nrow = 7))[rep(1:nrow(out$draws), each=ntheta),])
  steady.states <- rbind(steady.states, t(matrix(unlist(out$lls)[grep(names(unlist(out$lls)),pattern = "tail")], nrow = 79))[rep(1:nrow(out$draws), each=ntheta),])
  lls.ipt <- rbind(lls.ipt, t(matrix(unlist(out$lls)[grep(names(unlist(out$lls)), pattern="llthetas")], nrow=2)))
  thetas <- append(thetas, unlist(out$lls)[grep(names(unlist(out$lls)), pattern="^theta")])
}

allparams <- draws; allparams[,'theta'] <- thetas
colnames(steady.states) <- make.state.names()

sum.dist <- t(sapply(seq(ntheta, nrow(steady.states), by=ntheta),function(z) {
  get.summary.stats(
    rbind(cbind("time"=1,steady.states[z,,drop=F]),
          cbind("time"=2,steady.states[z,,drop=F])),
    draws[z,])[1,]}))
 


# ll1 <- lapply(X=1:nrow(sum.dist), function(x) log.lik(sum.dist[x,-c(1:2)]))
ll1b <- lapply(X=1:nrow(sum.dist), function(x) log.lik(sum.dist[x,-c(1:2)]))

# ll1t <- t(matrix(unlist(ll1), nrow = 7))
ll1tb <- t(matrix(unlist(ll1b), nrow = 7))
summary(sum.dist[,'tb.prev.on.art'])
# summary(rowSums(ll1t))
summary(rowSums(ll1tb))
# summary((ll1t))
summary((ll1tb))
#c(ll.not,  ll.prev,  ll.mort,  ll.mdr,  ll.hiv,  ll.ipt.elg.tb,  ll.cd4)


# ll1t20 <- ll1t[rep(1:nrow(ll1t), each=ntheta),]
ll1t20b <- ll1tb[rep(1:nrow(ll1tb), each=ntheta),]
rm(ll1tb)

# load(paste0(savelocation,"generated_data/lls.ipt.",date,".Rdata"))
# lls.ipt <- lls.ipt.vect
# # lapply(1:nrow(allparams), function(x) log.lik.ipt(run = rbind(append(0,steady.states[x,]), append(0,steady.states[x,])), params = allparams[x, ]))
# llmat <- cbind(ll1t20, lls.ipt)
llmatb <- cbind(ll1t20b, lls.ipt)
summary(llmatb)
gc()

# newlls <- rowSums(llmat)
newllsb <- rowSums(llmatb)
summary(newllsb)
rm(newllsb)

#ll.not,ll.prev,ll.mort,ll.mdr,ll.hiv,ll.ipt.elg.tb,ll.cd4; and from ipt, ll.ir and ll.irr

# newlls3b <- rowSums(llmatb[,c(1,3:7)]) + rowSums(llmatb[,8:9])
# summary(newlls3b)

newlls3c <- rowSums(llmatb[,c(1,3:7)]) + 6*llmatb[,9]
summary(newlls3c)



## sample with replacement with prob proportional to likelihood
param.id.draws <- sample(1:nrow(draws),100000,replace=T,prob=exp(newlls3c))
saveRDS(param.id.draws,paste0(savelocation,"generated_data/weighted-param-draws.",date,".rds"))
# rm(newlls3c)

length(unique(param.id.draws))
length(unique(rowSums(allparams[param.id.draws,-which(colnames(allparams)=='theta')])))
hist(param.id.draws, breaks = 100)

## now get steady states for weighted draws (i.e. our 'posterior')
ss.dist <- steady.states[param.id.draws,]
param.dist <- allparams[param.id.draws,]

gc()
un.ids <- param.id.draws[which(!duplicated(param.id.draws))]
re.map <- sapply(un.ids,function(x)
  which(param.id.draws==x))
  
weights <- numeric(); for (i in 1:length(un.ids)) weights[i] <- length(re.map[[i]])
summary(weights)
length(weights)

un.ss.dist <- ss.dist[which(!duplicated(param.id.draws)),]
un.param.dist <- param.dist[which(!duplicated(param.id.draws)),]

save(un.ss.dist,un.param.dist, un.ids, re.map, weights, 
        file=paste0(savelocation,"generated_data/ss_and_param_dist.",date,".rdata"))


# can run interventions here (separate source file set up to run in parallel)

# #re-create full ss.dist, param.dist:
# load(paste0("generated_data/ss_and_param_dist.",date,".rdata"), verbose = TRUE)
# idmap <- unlist(lapply(1:length(un.ids),function(x) {
#     rep(x,times=length(re.map[[x]]))
# }))
# ss.dist <- un.ss.dist[idmap,]
# param.dist <- un.param.dist[idmap,]

param.id.draws <- readRDS(paste0(savelocation,"generated_data/weighted-param-draws.",date,".rds"))

load(file=paste0(savelocation,"generated_data/ss_and_param_dist.",date,".rdata"))

# look at calibration
(table2 <- make.table.2(un.ss.dist, un.param.dist, weights))
table2
# make.table.2(steady.states[seq(1,10001,by=20),], draws[seq(1,10001,by=20),], rep(1,(10000/20+1)), q=c(0.5, seq(0,1,by=0.2)))
saveRDS(table2, file = paste0("table2.",date,".RDS"))
mean(thetas[param.id.draws]); quantile(thetas[param.id.draws], c(0.025,0.25,0.5,0.75,0.975))


# look at interventions
load(paste0("generated_data/sum.stat",date,".Rdata"), verbose = T)
load(paste0("generated_data/limited.sum.stat",date,".Rdata"), verbose = T)
load(paste0("generated_data/sensis.sum.stat",date,".Rdata"), verbose = T)

q <- c(0.025,0.5,0.975)

# projections
to.khay.pop(wtd.quantile(colSums(sum.stat.base[1:51,"tb.inc",]), weights = weights, q)*0.1)
to.khay.pop(wtd.quantile(colSums(sum.stat.base[1:51,"tb.inc",]), weights = weights, q)*0.1)/5
to.khay.pop(wtd.quantile(sum.stat.base[1,"tb.inc",], weights = weights, q))
to.khay.pop(wtd.quantile(sum.stat.base[1,"tb.mort",], weights = weights, q))
to.khay.pop(wtd.quantile(colSums(sum.stat.base[1:51,"tb.inc.hivpos",]*sum.stat.base[1:51,"hiv.prev",]), weights, q)*0.1)
to.khay.pop(wtd.quantile(colSums(sum.stat.base[1:51,"tb.inc.onart",]*sum.stat.base[1:51,"hiv.prev",]*sum.stat.base[1:51,"art.cov.among.all.hiv",]), weights, q)*0.1)
to.khay.pop(wtd.quantile(colSums(sum.stat.base[1:51,"tb.mort",]), weights, q)*0.1)
to.khay.pop(wtd.quantile(colSums(sum.stat.base[1:51,"tb.hiv.mort",]), weights, q)*0.1)

mainresults <- make.main.tables.and.figures(filetag=date,plot = T)
make.fig2.alt(sum.stat.int = sum.stat.int, sum.stat.base = sum.stat.base,weights = weights,times=seq(0,5,by=0.1))
make.fig2.sup(filetag = date)

#abstract
wtd.quantile(sum.stat.base[51,"tb.inc.onart",],weights,q)
wtd.quantile(sum.stat.int[51,"tb.inc.onart",],weights,q)
wtd.quantile(((sum.stat.base[51,"tb.inc.onart",]- sum.stat.int[51,"tb.inc.onart",])/sum.stat.base[51,"tb.inc.onart",]),weights,q)

wtd.quantile(sum.stat.base[51,"tb.inc",],weights,q)

#overview
mainresults[1:3,]

#proportion treated
mainresults[1:3,1]/392000
wtd.quantile(apply(sum.stat.int[,"onipt",],2,max), weights, q)
wtd.quantile(apply(sum.stat.int[,"onipt",]/sum.stat.int[,"hiv.prev",],2,max), weights, q)
wtd.quantile(apply(sum.stat.int[,"onipt",]/sum.stat.int[,"hiv.prev",]/sum.stat.int[,"art.cov.among.all.hiv",],2,max), weights, q)
wtd.quantile(apply(sum.stat.int[,"onipt",]/sum.stat.int[,"hiv.prev",],2,mean), weights, q)
wtd.quantile(apply(sum.stat.int[,"onipt",]/sum.stat.int[,"hiv.prev",]/sum.stat.int[,"art.cov.among.all.hiv",],2,mean), weights, q)

#follow up time
wtd.quantile(1e5*apply(sum.stat.int[1:51,"onipt",],2,sum)/(apply(sum.stat.int[1:51,"ipt.init",],2,sum)), weights, q)
wtd.quantile(1e5*apply(sum.stat.int[1:51,"postipt",],2,sum)/(apply(sum.stat.int[1:51,"ipt.init",],2,sum)), weights, q)


#averted
mainresults[1:3,1:12]

lateresults <- make.main.tables.and.figures(filetag=date, rows=1:151, plot = F)
lateresults

mainresults[4:6,]
mainresults[4:6,]/mainresults[1:3,]
lateresults[4:6,]/lateresults[1:3,]

wtd.quantile(
      apply((sum.stat.cont[1:51,"onipt",] - sum.stat.int[1:51,"onipt",])*100000, 2, sum)/
        apply(sum.stat.int[1:51,"tb.inc",] - sum.stat.cont[1:51,"tb.inc",], 2, sum),
  weights, q)

wtd.quantile(
  apply((sum.stat.cont[1:51,"onipt",] - sum.stat.int[1:51,"onipt",])*100000, 2, sum)/
    apply(sum.stat.int[1:51,"tb.mort",] - sum.stat.cont[1:51,"tb.mort",], 2, sum),
  weights, q)


(mainresults[4:6,2]/mainresults[4:6,1])/(mainresults[1:3,2]/mainresults[1:3,1])


#over 15 years: 
wtd.quantile(to.khay.pop(colSums((sum.stat.cont[,'onipt',]-sum.stat.int[,'onipt',])*1e5*0.1))/
           to.khay.pop(colSums((sum.stat.cont[,'tb.inc',]-sum.stat.int[,'tb.inc',])*0.1)),weights, q)

wtd.quantile(to.khay.pop(colSums((sum.stat.cont[,'onipt',]-sum.stat.int[,'onipt',])*1e5*0.1))/
         to.khay.pop(colSums((sum.stat.cont[,'tb.mort',]-sum.stat.int[,'tb.mort',])*0.1)),weights, q)




## now let's look at the composition of new infections over time
## how much is attributable to reduced infections within those who did not receive IPT

## using the same intervention but in a world where IPT has no effect
basediff <- colSums(sum.stat.base[1:51,"tb.inc",]-sum.stat.noeffect[1:51,"tb.inc",]); summary(basediff) #essentially zero - check.

## now get the number of cases averted by IPT effect while on IPT
#( with vs without ipt effect, sum over all 0.1-year time steps of tb.inc on, after, and without ipt)
inc.diff <- to.khay.pop(apply(sum.stat.noeffect[1:51,c("tb.inc"),]*.1,2,sum) - 
                              apply(sum.stat.int[1:51,c("tb.inc"),]*.1,2,sum))

inc.diff.post.ipt <- to.khay.pop(apply(sum.stat.noeffect[1:51,c("tb.inc.post.ipt"),]*.1,2,sum) - 
                                 apply(sum.stat.int[1:51,c("tb.inc.post.ipt"),]*.1,2,sum))

inc.diff.nonipt <- to.khay.pop(apply(sum.stat.noeffect[1:51,"tb.inc.nonipt",]*.1,2,sum)- 
                              apply(sum.stat.int[1:51,"tb.inc.nonipt",]*.1,2,sum))

inc.diff.on.ipt <- to.khay.pop(apply(sum.stat.noeffect[1:51,c("tb.inc.on.ipt"),]*.1,2,sum) - 
                                 apply(sum.stat.int[1:51,c("tb.inc.on.ipt"),]*.1,2,sum))


inc.diff.ipt <- inc.diff.on.ipt + inc.diff.post.ipt


summary(inc.diff.ipt)
summary(inc.diff.nonipt)

wtd.quantile(inc.diff.ipt, weights, q)
wtd.quantile(inc.diff.nonipt, weights, q)
#check:
wtd.quantile(inc.diff.ipt + inc.diff.nonipt, weights, q)  # should equal total cases averted with 12mo int in table 4. It's low because my get summary stats had too low a parameter for a latent fast on ipt, so the missing cases are there. 
mainresults[1:3,3]
wtd.quantile(inc.diff, weights, q)
wtd.quantile(inc.diff.ipt/inc.diff, weights, q)
wtd.quantile(inc.diff.on.ipt/inc.diff, weights, q)
wtd.quantile(inc.diff.post.ipt/inc.diff, weights, q)
wtd.quantile(inc.diff.nonipt/inc.diff, weights, q)

# magnitude of effect of IPT on progression
mean(thetas[param.id.draws]); quantile(thetas[param.id.draws], c(0.025,0.25,0.5,0.75,0.975))

# # add ART initiation with active TB diagnosis?
# apply(sum.stat.base[1:51,"art.init",]*.1,2,sum) #per 100k population
# apply(sum.stat.base[1:51,"tb.inc.hivpos",]*.1,2,sum) # per 100k hiv+
# apply(sum.stat.base[1:51,"tb.prev.hiv.pos.nonart",]*.1,1,sum) # per 100k hiv+ not on art
# apply(sum.stat.base[1:51,"art.cov.among.all.hiv",]*.1,2,sum) # proportion of hiv+

# # rate at which people who are HIV+ and art- start tb treatment, per year
# wtd.quantile(sum.stat.base[51,"tb.prev.hiv.pos.nonart",]*un.param.dist[,"tau_d0_h3"],weights,0.5) /
# # rate at which people who are HIV+ and art- start art (excluding initiations due to tb treatment)
# wtd.quantile(sum.stat.base[51,"art.init",]*10/1e5/(sum.stat.base[51,"hiv.prev",]*(1-sum.stat.base[51,"art.cov.among.all.hiv",])),weights,c(0.5, 0.1,0.9))
# 
# # about 14%, probably cose enough to the current 25% per Gary. We could add in this route of ART initiation 
# # (ie. when someone with HIV not on ART is diagnosed with active TB) and reduce the other ART initiation rates accordingly.
# 


#only prevents reinfection:
nothetaresults <- make.main.tables.and.figures(filetag = date, plot = F, notheta = T)
nothetaresults[4:6,]
nothetaresults[4:6,]/mainresults[1:3,]
table2["LTBI prevalence in HIV",]

wtd.quantile((apply(sum.stat.noeffect[2:51,c("tb.inc.post.ipt"),]*.1,2,sum)/apply(sum.stat.noeffect[2:51,c("tb.inc.on.ipt"),]*.1,2,sum)), weights,q)

# # with IPT only for ART initiators:
# ltdresults <- make.main.tables.and.figures(ltd = T, filetag=date, row=1:51, plot=F)
# ltdresults[4:6,1:5]/ltdresults[1:3, 1:5]

# # average CD4s on vs off
# lu <- make.look.up.list()
# wtd.quantile(rowSums(un.ss.dist[,lu$h1art.i])/rowSums(un.ss.dist[,lu$allart.i]), weights,q)
# wtd.quantile(rowSums(un.ss.dist[,lu$h2art.i])/rowSums(un.ss.dist[,lu$allart.i]), weights,q)

# # xpert only with IPT:
# load(paste0("generated_data/sum.statnoartscreen.norecal.",date,".Rdata"), verbose = T)
# to.khay.pop(wtd.quantile(colSums(sum.stat.base[1:51,"tb.inc",]), weights = weights, q)*0.1)
# to.khay.pop(wtd.quantile(colSums(sum.stat.base[1:51,"tb.inc.hivpos",]*sum.stat.base[1:51,"hiv.prev",]), weights, q)*0.1)
# to.khay.pop(wtd.quantile(colSums(sum.stat.base[1:51,"tb.inc.onart",]*sum.stat.base[1:51,"hiv.prev",]*sum.stat.base[1:51,"art.cov.among.all.hiv",]), weights, q)*0.1)
# to.khay.pop(wtd.quantile(colSums(sum.stat.base[1:51,"tb.mort",]), weights, q)*0.1)
# to.khay.pop(wtd.quantile(colSums(sum.stat.base[1:51,"tb.hiv.mort",]), weights, q)*0.1)
# 
# make.main.tables.and.figures(filetag = "noartscreen.norecal.20170626",plot = F)[1:3,]


wtd.quantile((sum.stat.art.base[51,"hiv.prev",])/(sum.stat.base[1,"hiv.prev",]), weights = weights, q)
wtd.quantile((sum.stat.art.base[51,"hiv.prev",]), weights = weights, q)
wtd.quantile((sum.stat.art.base[51,"art.cov",]), weights = weights, q)
wtd.quantile((sum.stat.base[51,"art.cov",]), weights = weights, q)
wtd.quantile((sum.stat.base[51,"art.cov.among.all.hiv",]), weights = weights, q)
wtd.quantile((sum.stat.art.base[51,"art.cov.among.all.hiv",]), weights = weights, q)
wtd.quantile((sum.stat.art.base[51,"tb.inc",]/sum.stat.base[51,"tb.inc",]), weights = weights, q)
# wtd.quantile((sum.stat.base[51,"art.init",]), weights = weights, q)
# wtd.quantile((sum.stat.art.base[51,"art.init",]), weights = weights, q)

artresults <- make.main.tables.and.figures(filetag = date, plot = F, art = T)
artresults[4:6,] # now includes new baseline with ART also
mainresults[1:3,]

# #clin trial sim
# load(paste0(savelocation,"generated_data/clintrialsim.Rdata"), verbose = T)
# m12 <- array(unlist(twelvemonth), dim=c(length(weights),33,51),
#              dimnames = list("run"=1:length(weights),"stat"=dimnames(sum.stat.base)[[2]], "time"=seq(0,5,by=0.1)))
# m36 <- array(unlist(thirtysixmonth), dim=c(51, 33,length(weights)),
#              list("time"=seq(0,5,by=0.1), "stat"=dimnames(sum.stat.base)[[2]], "run"=1:length(weights)))
# 
# 
# plot(seq(0,5,0.1), apply(m12[,"tb.inc",], 1, wtd.quantile, weights, 0.5), ylim=c(0,10000))
# points(seq(0,5,0.1), apply(m36[,"tb.inc",], 1, wtd.quantile, weights, 0.5))


## looking at prior vs. posteriors
load(paste0("generated_data/ss_and_param_dist.",date,".rdata"), verbose = TRUE)
idmap <- unlist(lapply(1:length(un.ids),function(x) {
    rep(x,times=length(re.map[[x]]))
}))
ss.dist <- un.ss.dist[idmap,]
param.dist <- un.param.dist[idmap,]

## notifications
# pdf("plots/tbnotified_prior_posterior.pdf")
par(mfrow=c(3,2), mar=c(5,4,1,3), oma=c(0,0,0,0))
plot(density(sum.stat.base[1,'tb.notified',],weights = weights/sum(weights)),xlab="",ylab="",main="")
title(xlab = "TB notifications", cex.lab = 1.2, line = 2)
curve(dlnorm(x, meanlog=log(1630), sdlog=log(1.1)), 600,2500,add=T,col=2,lty=2)
legend("topright",c("prior","posterior"),col=2:1,lty=c(2,1),bty="n")
# dev.off()

plot(density(sum.stat.base[1,'hiv.prev',],weights = weights/sum(weights)),xlab="",ylab="",main="")
title(xlab = "HIV prevalence", cex.lab = 1.2, line = 2)
ms <-  twCoefLogitnormCi(0.2224, 0.312)
curve(dlogitnorm(x, ms[1], ms[2]), 0,1,add=T,col=2,lty=2)
legend("topright",c("prior","posterior"),col=2:1,lty=c(2,1),bty="n")

plot(density(sum.stat.base[1,'tb.mort',],weights = weights/sum(weights)),xlab="",ylab="",main="")
title(xlab = "TB mortality rate", cex.lab = 1.2, line = 2)
curve(dlnorm(x, log(500), log(1.3)), add=T, col=2,lty=2)
legend("topright",c("prior","posterior"),col=2:1,lty=c(2,1),bty="n")

plot(density(sum.stat.base[1,'prop.MDR',],weights = weights/sum(weights)),xlab="",ylab="",main="")
title(xlab = "MDR proportion", cex.lab = 1.2, line = 2)
curve(dlogitnorm(x, logit(0.04), ((logit(0.055)-logit(0.03))/4)), add=T,col=2,lty=2)
legend("topright",c("prior","posterior"),col=2:1,lty=c(2,1),bty="n")

plot(density(sum.stat.base[1,'cd4.artinit',],weights = weights/sum(weights)),xlab="",ylab="",main="")
title(xlab = "Average CD4 at ART initiation", cex.lab = 1.2, line = 2)
curve(dlnorm(x, log(175), log(1.1)), add=T,col=2,lty=2)
legend("topright",c("prior","posterior"),col=2:1,lty=c(2,1),bty="n")

plot(density(sum.stat.base[1,'tb.prev.on.art',],weights = weights/sum(weights)),xlab="",ylab="",main="")
title(xlab = "Active TB prevalence on ART\n(prior to IPT intervention)", cex.lab = 1.2, line = 3)
ms <-  twCoefLogitnormCi(0.104, 0.131)
curve(dlogitnorm(x, ms[1], ms[2]), 0,1,add=T,col=2,lty=2)
legend("topright",c("prior","posterior"),col=2:1,lty=c(2,1),bty="n")


## input parameters
bounds.raw <- read.csv("parameters/parameters-working-fixmdr.csv", header=TRUE)[,c("param","value","min","max")]
rownames(bounds.raw) <- bounds.raw[,"param"]
ppplot <- function(paramname, altname)
{
  if(missing(altname)) altname <- paramname
  plot(density(un.param.dist[,paramname], weights=weights/sum(weights)),xlab="",#paramnames[which(colnames(un.param.dist)==paramname)],
       ylab="",main="")
  title(xlab = altname, cex.lab = 1.1, line = 2.5)
  curve(dunif(x, bounds.raw[paramname,"min"], bounds.raw[paramname,"max"]), add=T,col=2,lty=2)
 # legend("topright",c("prior","posterior"),col=2:1,lty=c(2,1),bty="n")
}
par(mfrow=c(3,4), mar=c(3,3,1,1))
for (name in colnames(un.param.dist)[!is.na(bounds.raw[,'min'])]) ppplot(name)
#skewed:
#Figure S3
par(mfrow=c(4,4), mar=c(4,3,1,1))
ppplot('beta', 'beta\n(transmission rate)')
legend("topright",c("prior","posterior"),col=2:1,lty=c(2,1),bty="n", cex=1.2)
ppplot('phi', 'phi   (relative TB\ntransmissibility in HIV)')
ppplot('prob.MDR', 'xi   (MDR proportion\nof new infections')
ppplot('arttheta', 'arttheta   (ART\nprotection against TB)')
ppplot('h1theta', 'h1theta   (reduced\nTB progression, CD4>200')
ppplot('h2theta', 'h2theta   (reduced\nTB progression, CD4>500')
ppplot('pi_h', 'pi_h3   (rapid\nprogression rate)')
ppplot('epsilon_h', 'epsilon_h3\n(slow progession rate)')
ppplot('tau_d0_h0', 'tau_d0_h0\n(TB treatment rate)')
ppplot('tau_d0_h3', 'tau_d0_h3\n(TB treatment rate)')
ppplot('tau_d1', 'tau_d1\n(MDR treatment rate)')
ppplot('mu_tb_h', 'mu_tb_h3 (TB-specific\nexcess mortality)')
ppplot('eta0', 'eta0 \n(HIV infection rate)')
ppplot('eta1', 'eta1   (early\nHIV progression rate)')
ppplot('chi2_h3', 'chi_h3   (ART initaition\nrate, advanced HIV)')
ppplot('rho', 'rho\n(LTBI prevalence at 15yo)')
# plot(1,lwd=0,axes=F,xlab="",ylab="", col='white');




# PRCC sensitivity analysis, among the unique sims resampled at least once
paramnames <- c(
  "Effective contacts per infectious person-year",
  "Relative transmissibility of TB in HIV+",
  "Proportion of TB infections that are MDR",
  "Rate of leaving of elevated risk for fast progression",
  "Relative risk of TB progression, ART vs no ART",
  "Relative risk of TB progression, middle versus low CD4",
  
  "Relative risk of TB progression, high vs low CD4",
  "Rapid progression rate after recent infection (HIV-)",
  "Rapid progression rate after recent infection (HIV+, CD4<200)",
  "Annual TB reactivation rate (HIV-)",
  "Annual TB reactivation rate (advanced HIV)",
  "Annual TB treatment rate (drug-susceptible, HIV-)",
  
  "TB treatment initiation rate (HIV+, CD4<200)",
  "Treatment initiation rate (active MDR TB)",
  "Relative risk of reinfection, latently-infected",
  "Background mortality rate",
  "TB specific excess mortality rate (HIV-)",
  "TB specific excess mortality rate  (HIV+)",
  
  "HIV specific excess mortality rate (CD4<200)",
  "HIV specific excess mortality rate (CD4 200-500)",
  "HIV specific excess mortality rate (CD4 >500)",
  "Annual rate of new HIV infection",
  "Rate of HIV progression from [CD4>500] to [200<CD4<500]",
  'Rate of HIV progression from [200<CD4<500] to [CD4<200]',

  'CD4+ reconstitution rate on ART, from CD4<200 to 200<CD4<500 ',
  'CD4+ reconstitution rate on ART, from 200<CD4<500 to CD4>500',
  '',
  'Annual rate of ART initiation, CD4 200-500',
  'Annual rate of ART initiation, CD4<200',
  'Percentage of 15 year-olds with latent TB infection',
  
  'Sensitivity of Xpert MTB/RIF screening for active TB',
  'Proportion of TB not initiating treatment after diagnosis',
  'Proportion of treated TB not initiating ART',
  "Relative risk of TB progression, IPT vs no IPT",
  "Rate of leaving IPT therapy",
  'iptwithart',
  'iptafterart',
  
  'artscreen',
  'iptscreen',
  'infonipt'
)

require(sensitivity)

load(paste0("generated_data/sum.stat",date,".Rdata"), verbose = T)
load(paste0("generated_data/limited.sum.stat",date,".Rdata"), verbose = T)
load(paste0("generated_data/sensis.sum.stat",date,".Rdata"), verbose = T)
load(file=paste0(savelocation,"generated_data/ss_and_param_dist.",date,".rdata"), verbose = T)

#% cases averted, % deaths averted, ynt
p.inc <- pcc(X = data.frame(un.param.dist[,]), y = (sum.stat.base[51,"tb.inc",]-sum.stat.int[51,"tb.inc",])/sum.stat.base[51,"tb.inc",],rank = T)
p.mort <- pcc(X = data.frame(un.param.dist[,]), y = (sum.stat.base[51,"tb.mort",]-sum.stat.int[51,"tb.mort",])/sum.stat.base[51,"tb.mort",],rank = T)
p.ynt <- pcc(X = data.frame(un.param.dist[,]), y = (colSums(sum.stat.int[1:51,"onipt",])/(colSums(sum.stat.base[1:51,"tb.inc",]-sum.stat.int[1:51,"tb.inc",]))),rank = T)
# p2 <- pcc(X = data.frame(t(sum.stat.base[1,,])), y = (sum.stat.base[51,"tb.inc",]-sum.stat.int[51,"tb.inc",])/sum.stat.base[51,"tb.inc",],rank = T)

n <- 5
plotthese <- union(union(rev(rev(order(abs(p.inc$PRCC$original)))[1:n]), rev(rev(order(abs(p.mort$PRCC$original)))[1:n])), rev(rev(order(abs(p.ynt$PRCC$original)))[1:n]))
length(plotthese)
plotthese <- plotthese[order(abs(p.inc$PRCC$original[plotthese]))]
# plotthese2 <- rev(rev(order(abs(p2$PRCC$original)))[1:15])


require('fields'); require('RColorBrewer')
par(mfrow=c(1,1), oma=c(8,0,2,0), mar=c(4,18,1,6))
icol <- c(rev(brewer.pal(9,"Oranges")),"white",brewer.pal(9,"Blues"))
image(z=t(cbind(p.inc$PRCC$original, p.mort$PRCC$original, p.ynt$PRCC$original)[plotthese,]), 
  zlim=c(-1,1), col=icol, xaxt='n', yaxt='n')
mtext("Correlation of parameters with key outcomes", side = 3, cex=1.2, outer=T, font=2)
axis(at = (0:2)/2,labels=FALSE, side=1)
text(x=(0:2)/2-0.25, y=par()$usr[3]-0.18*(par()$usr[4]-par()$usr[3]), labels=c("Incidence reduction", "Mortality reduction", "Years of IPT needed\nto prevent one case"), srt=45, xpd=NA) 
axis(at=seq(0,1, length=length(plotthese)), labels=paramnames[plotthese], cex.axis=0.8, side=2, xpd=NA, las=1)
axis(at=seq(0,1, length=length(plotthese)), 
 labels=paste0(signif(apply(un.param.dist[,plotthese],2,min),2),"-",signif(apply(un.param.dist[,plotthese],2,max),2)), cex.axis=0.7, side=4, xpd=NA, las=1)

text(x = 0, y=(0:8)/8,labels = round(p.inc$PRCC$original[plotthese],1), cex=0.8)
text(x =1/2, y=(0:8)/8,labels = round(p.mort$PRCC$original[plotthese],1), cex=0.8)
text(x =2/2, y=(0:8)/8,labels = round(p.ynt$PRCC$original[plotthese],1), cex=0.8)

box("plot")

par(oma=c( 0,0,0,1), mar=c(4,6,4,4))# reset margin to be much smaller.
image.plot(legend.cex=0.9, legend.only=TRUE, zlim=c(-1,1), col=icol, horizontal=TRUE, legend.lab ="Partial rank correlation coefficient (PRCC)") # image.plot tricked into  plotting in margin of old setting 
set.panel() # reset plotting device

# # not sure if i need this? 
# props.by.hiv.dist <- apply(ss.dist,1,props.by.hiv)
# rowMeans(props.by.hiv.dist)
# 
# no.hiv.death <- param.dist[,'mu_d'] + param.dist[,'mu_tb']
# mean(no.hiv.death)
# 
# weighted.death.rate <- (no.hiv.death*props.by.hiv.dist[1,] +
#                           (no.hiv.death+param.dist[,'mu_h1'])*props.by.hiv.dist[2,] +
#                           (no.hiv.death+param.dist[,'mu_h2'])*props.by.hiv.dist[3,] +
#                           (no.hiv.death+param.dist[,'mu_h3'])*props.by.hiv.dist[4,])
# mean(weighted.death.rate)
# 
# 1/(param.dist[,'tau_d0'] + weighted.death.rate) %>% mean
# 
# 



# 
#   table.s1 <- rbind(c(wtd.quantile((sum.stat.base[51,'tb.inc',]- sum.stat.int[51,'tb.inc',])/sum.stat.base[51,'tb.inc',], weights, 0.5),
# wtd.quantile((sum.stat.base[51,'tb.prev',]- sum.stat.int[51,'tb.prev',])/sum.stat.base[51,'tb.prev',], weights, 0.5),
#              wtd.quantile((sum.stat.base[51,'tb.mort',]- sum.stat.int[51,'tb.mort',])/sum.stat.base[51,'tb.mort',], weights, 0.5)),
#                     c(wtd.quantile((sum.stat.base[51,'tb.inc.hivpos',]- sum.stat.int[51,'tb.inc.hivpos',])/sum.stat.base[51,'tb.inc.hivpos',], weights, 0.5),
#                                    wtd.quantile((sum.stat.base[51,'tb.prev.hiv.pos',]- sum.stat.int[51,'tb.prev.hiv.pos',])/sum.stat.base[51,'tb.prev.hiv.pos',], weights, 0.5),
#                                                 wtd.quantile((sum.stat.base[51,'tb.hiv.mort',]- sum.stat.int[51,'tb.hiv.mort',])/sum.stat.base[51,'tb.hiv.mort',], weights, 0.5)),
#   c(wtd.quantile((sum.stat.noeffect[51,'tb.inc.eveript',]- sum.stat.int[51,'tb.inc.eveript',])/sum.stat.noeffect[51,'tb.inc.eveript',], weights, 0.5),
#     NA,#wtd.quantile((sum.stat.noeffect[51,'tb.prev.eveript',]- sum.stat.int[51,'tb.prev.eveript',])/sum.stat.noeffect[51,'tb.prev.eveript',], weights, 0.5),
#     NA)#wtd.quantile((sum.stat.noeffect[51,'tb.mort.eveript',]- sum.stat.int[51,'tb.mort.eveript',])/sum.stat.noeffect[51,'tb.mort.eveript',], weights, 0.5))
# )
# 
#   table.s2 <- rbind(c(wtd.quantile((sum.stat.base[51,'tb.inc',]- sum.stat.cont[51,'tb.inc',])/sum.stat.base[51,'tb.inc',], weights, 0.5),
#                       wtd.quantile((sum.stat.base[51,'tb.prev',]- sum.stat.cont[51,'tb.prev',])/sum.stat.base[51,'tb.prev',], weights, 0.5),
#                       wtd.quantile((sum.stat.base[51,'tb.mort',]- sum.stat.cont[51,'tb.mort',])/sum.stat.base[51,'tb.mort',], weights, 0.5)),
#                     c(wtd.quantile((sum.stat.base[51,'tb.inc.hivpos',]- sum.stat.cont[51,'tb.inc.hivpos',])/sum.stat.base[51,'tb.inc.hivpos',], weights, 0.5),
#                       wtd.quantile((sum.stat.base[51,'tb.prev.hiv.pos',]- sum.stat.cont[51,'tb.prev.hiv.pos',])/sum.stat.base[51,'tb.prev.hiv.pos',], weights, 0.5),
#                       wtd.quantile((sum.stat.base[51,'tb.hiv.mort',]- sum.stat.cont[51,'tb.hiv.mort',])/sum.stat.base[51,'tb.hiv.mort',], weights, 0.5)),
#                     c(NA,#wtd.quantile((sum.stat.noeffect[51,'tb.inc.eveript',]- sum.stat.cont[51,'tb.inc.eveript',])/sum.stat.noeffect[51,'tb.inc.eveript',], weights, 0.5),
#                       NA,#wtd.quantile((sum.stat.noeffect[51,'tb.prev.eveript',]- sum.stat.int[51,'tb.prev.eveript',])/sum.stat.noeffect[51,'tb.prev.eveript',], weights, 0.5),
#                       NA)#wtd.quantile((sum.stat.noeffect[51,'tb.mort.eveript',]- sum.stat.int[51,'tb.mort.eveript',])/sum.stat.noeffect[51,'tb.mort.eveript',], weights, 0.5))
#   )
  
## rownames(table.s1) <- c("no.hiv","hiv")
## print("Table S1")
## print(table.s1)

