## Script for running some of the paper analyses based on ABC

########################################################################
source("source/khay-utility-fixmdr.R")
library(dplyr)
library(abind)
library(parallel)
reload.source()

savelocation <- "../scratch/khaydata/"

## using output
set.seed(24543663)

run.basic <- as.logical(commandArgs(trailingOnly=TRUE)[1])
run.limited <- as.logical(commandArgs(trailingOnly=TRUE)[2])
run.sensis <- as.logical(commandArgs(trailingOnly=TRUE)[3])
date <- commandArgs(trailingOnly=TRUE)[4]


cores=detectCores()
print(cores)

load(file=paste0(savelocation,"generated_data/ss_and_param_dist.",date,".rdata"), verbose = T)

colnames(un.ss.dist) <- make.state.names()

lu <- make.look.up.list()

run.interventions <- function(un.ss.dist, un.param.dist, weights, filetag="")
{

  if(run.basic)
  {
    ## run basline epidemic
    params.base <- un.param.dist
    params.base[,'kappa'] <- 0
    params.base[,'kappa_onart'] <- 0
    params.base[,'gamma'] <- 1 #doesn't matter for baseline, right?
  
    ## get distribiton of base runs
    sum.stat.base <- mclapply(1:nrow(un.ss.dist),function(x) {
      print(x);
      tmp <- run.model(un.ss.dist[x,],0,15,
                       params=params.base[x,],
                       to.eq = FALSE) %>%
        get.summary.stats(.,params=params.base[x,])
      array(tmp, dim=c(dim(tmp),1), dimnames = list(NULL, dimnames(tmp)[[2]], x))
    }, mc.cores=cores) %>% abind(.,along=3)
  
    print("baseline projection completed")
    rm(params.base)
    
    #12mo IPT intervention
    params.int <- un.param.dist
    params.int[,'kappa'] <- 0.85
    params.int[,'kappa_onart'] <- .15
    params.int[,'gamma'] <- 1
  
    # 12 month intervention
    sum.stat.int <- mclapply(1:nrow(un.ss.dist),function(x) {
      print(x);
      tmp <- run.model(un.ss.dist[x,],0,15,
                       params=params.int[x,],
                       to.eq = FALSE) %>%
        get.summary.stats(.,params=params.int[x,])
      array(tmp, dim=c(dim(tmp),1), dimnames = list(NULL, dimnames(tmp)[[2]], x))
    }, mc.cores=cores) %>% abind(.,along=3)
  
    print("primary intervention completed")
    rm(params.int)
    
    # continuous intervention
    params.cont <- un.param.dist
    params.cont[,'kappa'] <- 0.85
    params.cont[,'kappa_onart'] <- .15
    params.cont[,'gamma'] <- 0 #no outflow from ipt
  
    sum.stat.cont <- mclapply(1:nrow(un.ss.dist),function(x) {
      print(x);
      tmp <- run.model(un.ss.dist[x,],0,15,
                       params=params.cont[x,],
                       to.eq = FALSE) %>%
        get.summary.stats(.,params=params.cont[x,])
      array(tmp, dim=c(dim(tmp),1), dimnames = list(NULL, dimnames(tmp)[[2]], x))
    }, mc.cores=cores) %>% abind(.,along=3)
    
    print("continuous intervention completed")
    rm(params.cont)
    
    saveRDS(apply(sum.stat.base[,c('tb.inc','tb.mort','onipt'),],c(1,2),median), file = paste0(savelocation,"generated_data/medianbaseline",filetag,".RDS"))
    saveRDS(apply(sum.stat.int[,c('tb.inc','tb.mort','onipt'),],c(1,2),median), file = paste0(savelocation,"generated_data/medianint_12month",filetag,".RDS"))
    saveRDS(apply(sum.stat.cont[,c('tb.inc','tb.mort','onipt'),],c(1,2),median), file = paste0(savelocation,"generated_data/medianint_continuous",filetag,".RDS"))
    
    save(sum.stat.base, sum.stat.int, sum.stat.cont, weights, file=paste0(savelocation,"generated_data/sum.stat",filetag,".Rdata"))
    rm(sum.stat.base, sum.stat.int, sum.stat.cont)
    gc()
  }
    
  
  if (run.limited)
  {
      # 12 month intervention, no effect from IPT (to classify population, and averted cases and deaths, by IPT history)
    params.noeffect <- un.param.dist
    params.noeffect[,'kappa'] <- 0.85
    params.noeffect[,'kappa_onart'] <- .15
    params.noeffect[,'gamma'] <- 1
    params.noeffect[,'inf_on_ipt'] <- 1 ## allow those on IPT to get infected
    params.noeffect[,'theta'] <- 1      ## no reduction in reactivation during or post ipt
    ## note, here ipt_screen is still 1, so IPT still has the benefit of screening people for active TB, but as long as all are screened anyway when starting ART, it shouldn't matter. (no screening is modeled with IPt initiation for those already on ART.)
  
    sum.stat.noeffect <- mclapply(1:nrow(un.ss.dist),function(x) {
      print(x);
      tmp <- run.model(un.ss.dist[x,],0,10,
                       params=params.noeffect[x,],
                       to.eq = FALSE) %>%
        get.summary.stats(.,params=params.noeffect[x,])
      array(tmp, dim=c(dim(tmp),1), dimnames = list(NULL, dimnames(tmp)[[2]], x))
    }, mc.cores=cores) %>% abind(.,along=3)
    
  
    print("no-effect intervention completed")
  
    # 12 month intervention, IPT only for ART initiators
    params.ltd <- un.param.dist
    params.ltd[,'kappa'] <- 0.85
    params.ltd[,'kappa_onart'] <- 0
    params.ltd[,'gamma'] <- 1
  
    sum.stat.ltd <- mclapply(1:nrow(un.ss.dist),function(x) {
      print(x);
      tmp <- run.model(un.ss.dist[x,],0,10,
                       params=params.ltd[x,],
                       to.eq = FALSE) %>%
        get.summary.stats(.,params=params.ltd[x,])
      array(tmp, dim=c(dim(tmp),1), dimnames = list(NULL, dimnames(tmp)[[2]], x))
    }, mc.cores=cores) %>% abind(.,along=3)
  
    print("initiators-only intervention completed")
  
    saveRDS(apply(sum.stat.noeffect[,c('tb.inc','tb.mort','onipt'),],c(1,2),median), file = paste0(savelocation,"generated_data/medianint_continuous",filetag,".RDS"))
    saveRDS(apply(sum.stat.ltd[,c('tb.inc','tb.mort','onipt'),],c(1,2),median), file = paste0(savelocation,"generated_data/medianint_ltd",filetag,".RDS"))
    save(sum.stat.noeffect, sum.stat.ltd, weights, file=paste0(savelocation,"generated_data/limited.sum.stat",filetag,".Rdata"))
    rm(sum.stat.noeffect, sum.stat.ltd, params.noeffect, params.ltd)
    gc()
  }
    
  if (run.sensis)
  {
    # increasing ART intervention
    params.art <- un.param.dist
    params.art[,'kappa'] <- 0.85
    params.art[,'kappa_onart'] <- .15
    params.art[,'gamma'] <- 1
    #setting projection=TRUE in dxdt below
    basearts <- rowSums(un.ss.dist[,c(lu$h1art.i,lu$h2art.i, lu$h3art.i)])/rowSums(un.ss.dist[,-c(lu$h0)])
  
    sum.stat.art <- mclapply(1:nrow(un.ss.dist),function(x) {
      print(x);
      tmp <- run.model(un.ss.dist[x,],0,10,
                       params=params.art[x,], projection=TRUE, baseart=basearts[x],
                       to.eq = FALSE) %>%
        get.summary.stats(.,params=params.art[x,])
      array(tmp, dim=c(dim(tmp),1), dimnames = list(NULL, dimnames(tmp)[[2]], x))
    }, mc.cores=cores) %>% abind(.,along=3)
  
    print("art-increase intervention completed")
  
    # increasing ART intervention, comparator
    params.art.base <- un.param.dist
    params.art.base[,'kappa'] <- 0
    params.art.base[,'kappa_onart'] <- 0
    params.art.base[,'gamma'] <- 1
    #setting projection=TRUE in dxdt below
    
    sum.stat.art.base <- mclapply(1:nrow(un.ss.dist),function(x) {
      print(x);
      tmp <- run.model(un.ss.dist[x,],0,10,
                       params=params.art.base[x,], projection=TRUE, 
                       to.eq = FALSE) %>%
        get.summary.stats(.,params=params.art.base[x,])
      array(tmp, dim=c(dim(tmp),1), dimnames = list(NULL, dimnames(tmp)[[2]], x))
    }, mc.cores=cores) %>% abind(.,along=3)
    
    print("art-increase baseline completed")
    
    # 12 month intervention, only prevents new infections
    params.notheta <- un.param.dist
    params.notheta[,'kappa'] <- 0.85
    params.notheta[,'kappa_onart'] <- 0.15
    params.notheta[,'gamma'] <- 1
    params.notheta[,'theta'] <- 1
  
    sum.stat.notheta <- mclapply(1:nrow(un.ss.dist),function(x) {
      print(x);
      tmp <- run.model(un.ss.dist[x,],0,10,
                       params=params.notheta[x,],
                       to.eq = FALSE) %>%
        get.summary.stats(.,params=params.notheta[x,])
      array(tmp, dim=c(dim(tmp),1), dimnames = list(NULL, dimnames(tmp)[[2]], x))
    }, mc.cores=cores) %>% abind(.,along=3)
  
    print("prevented reinfection only intervention completed")
  
    saveRDS(apply(sum.stat.art[,c('tb.inc','tb.mort','onipt'),],c(1,2),median), file = paste0(savelocation,"generated_data/medianint_art",filetag,".RDS"))
    saveRDS(apply(sum.stat.notheta[,c('tb.inc','tb.mort','onipt'),],c(1,2),median), file = paste0(savelocation,"generated_data/medianint_notheta",filetag,".RDS"))
  
    save(sum.stat.art, sum.stat.art.base, sum.stat.notheta, weights, file=paste0(savelocation,"generated_data/sensis.sum.stat",filetag,".Rdata"))
  }
  return()

}


run.interventions(un.ss.dist, un.param.dist, weights, filetag=date) #saves 'sum.stat's with date label
