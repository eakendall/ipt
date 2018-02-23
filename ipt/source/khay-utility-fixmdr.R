## -------------------------------------------- ##
## file for helper functions for khay ipt model ##
## -------------------------------------------- ##

# since we will source this file, let's load a few packages
library(deSolve)
library(RColorBrewer)
library(dplyr)
library(lhs)
library(logitnorm)
library(RgoogleMaps)
library(coda)
library('Hmisc')

##' reloads all the key files for the project
##' @return
##' @author asa
reload.source <- function(){
  print("sourcing all project files")
  source("source/khay-base-fixmdr.R")
  source("source/khay-utility-fixmdr.R")
}

# load.abc <- function(date){
#     rc <- readRDS(file=paste0("generated_data/ss_and_param_dist.",date,".rds"))
#     return(rc)
# }

##' loads parameters from csv file
##' @return
##' @author asa
load.params <- function(){
  rc <- read.csv("parameters/parameters-working-fixmdr.csv",as.is=T)[,1:2]
  params <-as.numeric(rc[,2])
  names(params) <- rc[,1]
  return(params)
}

##' gets column names for model
##' @param mod.out
##' @return
##' @author asa
get.colnames <- function(mod.out){
  state.names <- make.state.names()
  colnames(mod.out) <- c("time",state.names)
  return(mod.out)
}


##' makes state names vector
##' @title
##' @return
##' @author asa
make.state.names <- function(){
  
  s.states <-paste0("S-",c("0","1.1","2.1","3.1","1.2","2.2","3.2"))
  
  lf.states <-paste0("Lf-",c("0","1.1","2.1","3.1","1.2","2.2","3.2"))
  lf.states <-paste0(lf.states,"-",rep(0:1,each=length(lf.states)))
  
  ls.states <-paste0("Ls-",c("0","1.1","2.1","3.1","1.2","2.2","3.2"))
  ls.states <- paste0(ls.states,"-",rep(0:1,each=length(ls.states)))
  
  a.states <-paste0("A-",c("0","1.1","2.1","3.1","1.2","2.2","3.2"))
  a.states <- paste0(a.states,"-",rep(0:1,each=length(a.states)))
  
  sipt.states <-paste0("SIPT-",c("1.2","2.2","3.2"))
  #sipt.states <- paste0(sipt.states,"-",rep(0:1,each=length(sipt.states)))
  
  sipt2.states <-paste0("SIPT2-",c("1.2","2.2","3.2"))
  #sipt2.states <- paste0(sipt2.states,"-",rep(0:1,each=length(sipt2.states)))
  
  lipt.states <-paste0("LIPT-",c("1.2","2.2","3.2"))
  lipt.states <- paste0(lipt.states,"-",rep(0:1,each=length(lipt.states)))
  
  lipt2.states <-paste0("LIPT2-",c("1.2","2.2","3.2"))
  lipt2.states <- paste0(lipt2.states,"-",rep(0:1,each=length(lipt2.states)))
  
  liptf.states <-paste0("LIPTf-",c("1.2","2.2","3.2"))
  liptf.states <- paste0(liptf.states,"-",rep(0:1,each=length(liptf.states)))
  
  lipt2f.states <-paste0("LIPT2f-",c("1.2","2.2","3.2"))
  lipt2f.states <- paste0(lipt2f.states,"-",rep(0:1,each=length(lipt2f.states)))
  
  rc <- c(s.states,lf.states,ls.states,
          a.states,sipt.states,sipt2.states,
          lipt.states,lipt2.states, liptf.states,lipt2f.states)
  
  return(rc)
}


##' Gets key summary stats for a run
##' all stats are normalized to per 100k
##' @title
##' @param ode.lu$out
##' @param params expecting augmented params
##' @return
##' @author asa
get.summary.stats <- function(ode.out,params){
  
  params['iepsilon_0'] <- params['epsilon_0']
  params['iepsilon_h'] <- params['epsilon_h']
  params['ipi_0'] <- params['pi_0']
  params['ipi_h'] <- params['pi_h']
  
  if(colnames(ode.out)[1] != "time") {
    stop("expecting that the first column is time")
  }
  
  times <- ode.out[,1]
  ode.out <- ode.out[,-1]
  
  state.names <- make.state.names()
  
  # get lookup list
  lu <- make.look.up.list()
  
  
  ## ------------------ ##
  ## force of infection ##
  ## ------------------ ##
  inf.weight.d0 <-
    ode.out[,"A-0-0"] +
    params['phi']*(rowSums(ode.out[,intersect(lu$a.i,lu$d0.i)]) - ode.out[,"A-0-0"])
  
  inf.weight.d1 <- inf.weight.d0 * params['prob.MDR']/(1-params['prob.MDR'])
  
  N <- rowSums(ode.out)[1]

  foi <- params['beta'] / N * (inf.weight.d0 + inf.weight.d1)
  
  ## MDR cases are the only ones who get inected while on IPT
  ## except in sens. analyses where the flag inf_on_ipt is on
  inf.rate.on.ipt <-
    rowSums(ode.out[,c("SIPT-1.2","SIPT-2.2","SIPT-3.2")]*inf.weight.d1*params['beta']/N) +
    rowSums(ode.out[,c("LIPT-1.2-0","LIPT-2.2-0","LIPT-3.2-0")]*inf.weight.d1*params['psi']*params['beta']/N) +
    rowSums(ode.out[,c("LIPT-1.2-1","LIPT-2.2-1","LIPT-3.2-1")]*inf.weight.d1*params['psi']*params['beta']/N) +
    rowSums(ode.out[,c("LIPTf-1.2-0","LIPTf-2.2-0","LIPTf-3.2-0")]*inf.weight.d1*params['psi']*params['beta']/N) +
    rowSums(ode.out[,c("LIPTf-1.2-1","LIPTf-2.2-1","LIPTf-3.2-1")]*inf.weight.d1*params['psi']*params['beta']/N) +
    rowSums(ode.out[,c("SIPT-1.2","SIPT-2.2","SIPT-3.2")]*params['inf_on_ipt']*inf.weight.d0*params['beta']/N) + ## only in sens. analyis
    rowSums(ode.out[,c("LIPT-1.2-0","LIPT-2.2-0","LIPT-3.2-0")]*params['inf_on_ipt']*inf.weight.d0*params['psi']*params['beta']/N) +
    rowSums(ode.out[,c("LIPT-1.2-1","LIPT-2.2-1", "LIPT-3.2-1")]*params['inf_on_ipt']*inf.weight.d0*params['psi']*params['beta']/N) +
    rowSums(ode.out[,c("LIPTf-1.2-0","LIPTf-2.2-0", "LIPTf-3.2-0")]*params['inf_on_ipt']*inf.weight.d0*params['psi']*params['beta']/N) +
    rowSums(ode.out[,c("LIPTf-1.2-1","LIPTf-2.2-1","LIPTf-3.2-1")]*params['inf_on_ipt']*inf.weight.d0*params['psi']*params['beta']/N)
  
  inf.rate.post.ipt <-
    rowSums(ode.out[,c("SIPT2-1.2","SIPT2-2.2","SIPT2-3.2")]*inf.weight.d0*params['beta']/N) +
    rowSums(ode.out[,c("SIPT2-1.2","SIPT2-2.2","SIPT2-3.2")]*inf.weight.d1*params['beta']/N) +
    rowSums(ode.out[,c("LIPT2-1.2-0","LIPT2-2.2-0","LIPT2-3.2-0")]*inf.weight.d0*params["psi"]*params['beta']/N) +
    rowSums(ode.out[,c("LIPT2-1.2-0","LIPT2-2.2-0", "LIPT2-3.2-0")]*inf.weight.d1*params["psi"]*params['beta']/N) +
    rowSums(ode.out[,c("LIPT2-1.2-1","LIPT2-2.2-1","LIPT2-3.2-1")]*inf.weight.d0*params["psi"]*params['beta']/N) +
    rowSums(ode.out[,c("LIPT2-1.2-1","LIPT2-2.2-1","LIPT2-3.2-1")]*inf.weight.d1*params["psi"]*params['beta']/N) +
    rowSums(ode.out[,c("LIPT2f-1.2-0","LIPT2f-2.2-0","LIPT2f-3.2-0")]*inf.weight.d0*params["psi"]*params['beta']/N) +
    rowSums(ode.out[,c("LIPT2f-1.2-0","LIPT2f-2.2-0","LIPT2f-3.2-0")]*inf.weight.d1*params["psi"]*params['beta']/N) +
    rowSums(ode.out[,c("LIPT2f-1.2-1","LIPT2f-2.2-1","LIPT2f-3.2-0")]*inf.weight.d0*params["psi"]*params['beta']/N) +
    rowSums(ode.out[,c("LIPT2f-1.2-1","LIPT2f-2.2-1","LIPT2f-3.2-0")]*inf.weight.d1*params["psi"]*params['beta']/N)
  
  inf.rate.nonipt <-
    rowSums(ode.out[,lu$s.i]*inf.weight.d0*params['beta']/N) +
    rowSums(ode.out[,lu$s.i]*inf.weight.d1*params['beta']/N) +
    rowSums(ode.out[,intersect(lu$ls.i,lu$d0.i)]*inf.weight.d0*params['psi']*params['beta']/N) +
    rowSums(ode.out[,intersect(lu$ls.i,lu$d1.i)]*inf.weight.d0*params['psi']*params['beta']/N) +
    rowSums(ode.out[,intersect(lu$ls.i,lu$d0.i)]*inf.weight.d1*params['psi']*params['beta']/N) +
    rowSums(ode.out[,intersect(lu$ls.i,lu$d1.i)]*inf.weight.d1*params['psi']*params['beta']/N) +
    rowSums(ode.out[,intersect(lu$lf.i,lu$d0.i)]*inf.weight.d0*params['psi']*params['beta']/N) +
    rowSums(ode.out[,intersect(lu$lf.i,lu$d1.i)]*inf.weight.d0*params['psi']*params['beta']/N) +
    rowSums(ode.out[,intersect(lu$lf.i,lu$d0.i)]*inf.weight.d1*params['psi']*params['beta']/N) +
    rowSums(ode.out[,intersect(lu$lf.i,lu$d1.i)]*inf.weight.d1*params['psi']*params['beta']/N) # !! added lu$lf rows 8/19/16, although no change in state if same d (so shouldn't include in dxdt)
  
  ## TB prevalence overall
  tb.prev <-1e5*rowSums(ode.out[,lu$a.i])/rowSums(ode.out)
  
  ## TB prevalence, HIV+
  tb.prev.hiv.pos <-rowSums(ode.out[,setdiff(lu$a.i,lu$h0.i)])/rowSums(ode.out[,-c(lu$h0.i)])
  
  ## TB prevalance among those not on ARTs but HIV+
  tb.prev.hiv.pos.nonart <-rowSums(
    ode.out[,intersect(lu$a.i,lu$nonart.i)])/ 
    rowSums(ode.out[,lu$nonart.i])
  
  
  ## proportion of those going on "routine" ART (thus screened for IPT and not already known to have TB) with active TB
  ipt.elg.with.tb <-
    ( 
      ( rowSums(ode.out[,intersect(intersect(lu$a.i,lu$nonart.i), lu$h2.i)]) *params['chi2_h2'] + 
          rowSums(ode.out[,intersect(intersect(lu$a.i,lu$nonart.i),lu$h3.i)]) *params['chi2_h3'] + 
          rowSums(ode.out[,intersect(intersect(lu$a.i,lu$nonart.i),lu$h1.i)]) *params['chi2_h1'] ) 
    ) / 
    ( 
      ( rowSums(ode.out[,intersect(intersect(c(lu$s.i, lu$lf.i,lu$ls.i, lu$a.i), lu$nonart.i),lu$h2.i)]) *params['chi2_h2'] + 
          rowSums(ode.out[,intersect(intersect(c(lu$s.i, lu$lf.i,lu$ls.i, lu$a.i), lu$nonart.i),lu$h3.i)]) *params['chi2_h3'] + 
          rowSums(ode.out[,intersect(intersect(c(lu$s.i, lu$lf.i,lu$ls.i, lu$a.i), lu$nonart.i),lu$h1.i)]) *params['chi2_h1'] )
    )
  
  
  ## TB prevalance on ARTs (only applicable before IPT is introduced, corresponsonds to IPT initiators in Rankaga). 
  # Needs to include ART intiators also (28% of participants in rangaka)
  tb.prev.on.art <- 0.72* rowSums(
    ode.out[,intersect(lu$a.i,lu$allart.i)])/ 
    rowSums(ode.out[,lu$allart.i]) + 0.28*ipt.elg.with.tb
  
  
  
  
  ## LTBI prevalance (includes active TB too)
  ltbi.prob <-rowSums(ode.out[,c(lu$lf.i,lu$ls.i,lu$lipt.i,lu$lipt2.i,lu$liptf.i,lu$lipt2f.i,lu$a.i
  )])/rowSums(ode.out)
  
  ltbi.prob.in.hiv <- rowSums(ode.out[,setdiff(c(lu$lf.i,lu$ls.i,
                                                 lu$lipt.i,lu$lipt2.i,
                                                 lu$liptf.i,lu$lipt2f.i,
                                                 lu$a.i
  ),lu$h0.i)])/rowSums(ode.out[,-c(lu$h0.i)])
  
  ## TB Disease Incidence (includes fast progression from new infections and reinfections), all per 100k population
  ## Lf -> A , Ls -> A, LIPT f/s -> A, LITP2 f/s -> A
  tb.inc.nonipt <- 1e5*
    ( rowSums(ode.out[,intersect(lu$lf.i,lu$h0.i)])*params['pi_0']+
        rowSums(ode.out[,intersect(lu$lf.i,lu$h1art.i)])*params['pi_h']*params['arttheta']*params['h1theta'] +
        rowSums(ode.out[,intersect(lu$lf.i,lu$h1nonart.i)])*params['pi_h']*params['h1theta'] +
        rowSums(ode.out[,intersect(lu$lf.i,lu$h2art.i)])*params['pi_h']*params['arttheta']*params['h2theta'] +
        rowSums(ode.out[,intersect(lu$lf.i,lu$h2nonart.i)])*params['pi_h']*params['h2theta'] +
        rowSums(ode.out[,intersect(lu$lf.i,lu$h3art.i)])*params['pi_h']*params['arttheta'] +
        rowSums(ode.out[,intersect(lu$lf.i,lu$h3nonart.i)])*params['pi_h'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h0.i)])*params['epsilon_0'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h1art.i)])*params['epsilon_h']*params['arttheta']*params['h1theta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h1nonart.i)])*params['epsilon_h']*params['h1theta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h2art.i)])*params['epsilon_h']*params['arttheta']*params['h2theta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h2nonart.i)])*params['epsilon_h']*params['h2theta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h3art.i)])*params['epsilon_h']*params['arttheta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h3nonart.i)])*params['epsilon_h']) /rowSums(ode.out)
  
  tb.inc.on.ipt <- 1e5*
    (
      (ode.out[,intersect(intersect(lu$lipt.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h1art.i)])*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h2art.i)])*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h3art.i)])*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['ipi_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h1art.i)])*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h2art.i)])*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h3art.i)])*params['ipi_h']*params['arttheta']
    )/rowSums(ode.out)
  
  tb.inc.post.ipt <- 1e5*
    (
      (ode.out[,intersect(intersect(lu$lipt2.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h1art.i)])*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h2art.i)])*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h3art.i)])*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['ipi_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h1art.i)])*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h2art.i)])*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h3art.i)])*params['ipi_h']*params['arttheta']
    )/rowSums(ode.out)
  
  tb.inc <- tb.inc.nonipt + tb.inc.on.ipt + tb.inc.post.ipt
  
  ## mdr tb incidence
  tb.mdr.inc <- 1e5*(rowSums(ode.out[,"Lf-0-1",drop=F])*params['pi_0'] +
                       (ode.out[,intersect(intersect(lu$lf.i,lu$h1art.i),lu$d1.i)])*params['pi_h']*params['arttheta']*params['h1theta'] +
                       (ode.out[,intersect(intersect(lu$lf.i,lu$h1nonart.i),lu$d1.i)])*params['pi_h']*params['h1theta'] +
                       (ode.out[,intersect(intersect(lu$lf.i,lu$h2art.i),lu$d1.i)])*params['pi_h']*params['arttheta']*params['h2theta'] +
                       (ode.out[,intersect(intersect(lu$lf.i,lu$h2nonart.i),lu$d1.i)])*params['pi_h']*params['h2theta'] +
                       (ode.out[,intersect(intersect(lu$lf.i,lu$h3art.i),lu$d1.i)])*params['pi_h']*params['arttheta'] +
                       (ode.out[,intersect(intersect(lu$lf.i,lu$h3nonart.i),lu$d1.i)])*params['pi_h'] +
                       rowSums(ode.out[,"Ls-0-1",drop=F])*params['epsilon_0'] + 
                       (ode.out[,intersect(intersect(lu$ls.i,lu$h1art.i),lu$d1.i)])*params['epsilon_h']*params['arttheta']*params['h1theta'] +
                       (ode.out[,intersect(intersect(lu$ls.i,lu$h1nonart.i),lu$d1.i)])*params['epsilon_h']*params['h1theta'] +
                       (ode.out[,intersect(intersect(lu$ls.i,lu$h2art.i),lu$d1.i)])*params['epsilon_h']*params['arttheta']*params['h2theta'] +
                       (ode.out[,intersect(intersect(lu$ls.i,lu$h2nonart.i),lu$d1.i)])*params['epsilon_h']*params['h2theta'] +
                       (ode.out[,intersect(intersect(lu$ls.i,lu$h3art.i),lu$d1.i)])*params['epsilon_h']*params['arttheta'] +
                       (ode.out[,intersect(intersect(lu$ls.i,lu$h3nonart.i),lu$d1.i)])*params['epsilon_h'] + 
                       (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h1art.i)])*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
                       (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h2art.i)])*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
                       (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h3art.i)])*params['iepsilon_h']*params['arttheta'] +
                       (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h1art.i)])*params['ipi_h']*params['arttheta']*params['h1theta'] +
                       (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h2art.i)])*params['ipi_h']*params['arttheta']*params['h2theta'] +
                       (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h3art.i)])*params['ipi_h']*params['arttheta'] +
                       (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h1art.i)])*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
                       (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h2art.i)])*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
                       (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h3art.i)])*params['iepsilon_h']*params['arttheta'] +
                       (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h1art.i)])*params['ipi_h']*params['arttheta']*params['h1theta'] +
                       (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h2art.i)])*params['ipi_h']*params['arttheta']*params['h2theta'] +
                       (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h3art.i)])*params['ipi_h']*params['arttheta']
  )/rowSums(ode.out)
  
  tb.inc.hivpos <- 1e5*
    ( rowSums(ode.out[,intersect(lu$lf.i,lu$h1art.i)])*params['pi_h']*params['arttheta']*params['h1theta'] +
        rowSums(ode.out[,intersect(lu$lf.i,lu$h1nonart.i)])*params['pi_h']*params['h1theta'] +
        rowSums(ode.out[,intersect(lu$lf.i,lu$h2art.i)])*params['pi_h']*params['arttheta']*params['h2theta'] +
        rowSums(ode.out[,intersect(lu$lf.i,lu$h2nonart.i)])*params['pi_h']*params['h2theta'] +
        rowSums(ode.out[,intersect(lu$lf.i,lu$h3art.i)])*params['pi_h']*params['arttheta'] +
        rowSums(ode.out[,intersect(lu$lf.i,lu$h3nonart.i)])*params['pi_h'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h1art.i)])*params['epsilon_h']*params['arttheta']*params['h1theta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h1nonart.i)])*params['epsilon_h']*params['h1theta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h2art.i)])*params['epsilon_h']*params['arttheta']*params['h2theta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h2nonart.i)])*params['epsilon_h']*params['h2theta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h3art.i)])*params['epsilon_h']*params['arttheta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h3nonart.i)])*params['epsilon_h'] + 
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h1art.i)])*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h2art.i)])*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h3art.i)])*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['ipi_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h1art.i)])*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h2art.i)])*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h3art.i)])*params['ipi_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h1art.i)])*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h2art.i)])*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h3art.i)])*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['ipi_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h1art.i)])*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h2art.i)])*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h3art.i)])*params['ipi_h']*params['arttheta']
    )/rowSums(ode.out[,-c(lu$h0.i)])
  
  ## tb mortality not with hiv
  ## since WHO definintion for TB death is death for any reason with TB we will count this
  ## using the full populatin in the denominator here
  tb.mort.h0 <- 1e5*rowSums(ode.out[,intersect(lu$a.i,lu$h0.i)])*(params['mu_tb_0']+params['mu_d'])/rowSums(ode.out)
  
  tb.hiv.mort <-  1e5*
    (rowSums(ode.out[,intersect(lu$a.i,lu$h1nonart.i)]*(params['mu_tb_0']/2 + params['mu_tb_h']/2 + params['mu_h1'])) +
       rowSums(ode.out[,intersect(lu$a.i,lu$h2nonart.i)]*(params['mu_tb_0']/2 + params['mu_tb_h']/2 + params['mu_h2'])) +
       rowSums(ode.out[,intersect(lu$a.i,lu$h3nonart.i)]*(params['mu_tb_h'] + params['mu_h3'])) + 
       rowSums(ode.out[,intersect(lu$a.i,lu$h1art.i)]*(params['mu_tb_0']/2 + params['mu_tb_h']/2 + params['mu_h1']))*params['arttheta'] +
       rowSums(ode.out[,intersect(lu$a.i,lu$h2art.i)]*(params['mu_tb_0']/2 + params['mu_tb_h']/2 + params['mu_h2']))*params['arttheta'] +
       rowSums(ode.out[,intersect(lu$a.i,lu$h3art.i)]*(params['mu_tb_h'] + params['mu_h3']))*params['arttheta'])/rowSums(ode.out)
  
  tb.mort <- tb.mort.h0 + tb.hiv.mort
  
  ## Total population size (n in thousands) (all ages)
  pop.size <- rowSums(ode.out)
  
  ## HIV Prevalance
  hiv.prev <- rowSums(
    ode.out[,-lu$h0.i])/rowSums(ode.out)
  
  hiv.inc <- 1e5*rowSums(ode.out[,lu$h0.i]*params['eta0'])/rowSums(ode.out)
  
  ## HIV prevalence among notified TB cases
  hiv.prev.in.TB <- 1 - 
    (      (ode.out[,intersect(intersect(lu$a.i,lu$d0.i), lu$h0.i)])*params['tau_d0_h0'] + 
             (ode.out[,intersect(intersect(lu$a.i,lu$d1.i),lu$h0.i)])*params['tau_d1'] 
    )/
    (      (ode.out[,intersect(intersect(lu$a.i,lu$d0.i), lu$h0.i)])*params['tau_d0_h0'] + 
             rowSums(ode.out[,intersect(intersect(lu$a.i,lu$d0.i), lu$h1.i)])*(params['tau_d0_h0']+params['tau_d0_h3'])/2 + 
             rowSums(ode.out[,intersect(intersect(lu$a.i,lu$d0.i), lu$h2.i)])*(params['tau_d0_h0']+params['tau_d0_h3'])/2 + 
             rowSums(ode.out[,intersect(lu$a.i,lu$d1.i)])*params['tau_d1'] )    
  
  
  ## art coverage among those eligible (considering cd4<500 only, or >500 on art)
  art.cov <- rowSums(ode.out[,c(lu$h1art.i,lu$h2art.i,lu$h3art.i)])/rowSums(ode.out[,c(lu$h2.i,lu$h3.i,lu$h1art.i)])
  ## art coverage among all HIV+
  art.cov.among.all.hiv <- rowSums(ode.out[,c(lu$h1art.i,lu$h2art.i, lu$h3art.i)])/rowSums(ode.out[,-c(lu$h0)])
  
  ## art inititions at cd4 dependent rate, per 100k population per time step
  art.init.cd4 <-
    rowSums(ode.out[,c("S-1.1","Lf-1.1-0","Lf-1.1-1","Ls-1.1-0","Ls-1.1-1")]*params['chi2_h1']) +
    rowSums(ode.out[,c("S-2.1","Lf-2.1-0","Lf-2.1-1","Ls-2.1-0","Ls-2.1-1")]*params['chi2_h2']) +
    rowSums(ode.out[,c("S-3.1","Lf-3.1-0","Lf-3.1-1","Ls-3.1-0","Ls-3.1-1")]*params['chi2_h3']) +
    rowSums(ode.out[,c("A-1.1-0","A-1.1-1")]*params['chi2_h1']) +
    rowSums(ode.out[,c("A-2.1-0","A-2.1-1")]*params['chi2_h2']) +
    rowSums(ode.out[,c("A-3.1-0","A-3.1-1")]*params['chi2_h3'])
  
  ## art initiations due to acive TB treatment, pe 100k population per time step
  art.init.atb <-
    ode.out[,"A-1.1-0"]*(params['tau_d0_h0']+params['tau_d0_h3'])/2*(1-params['default_artatb']) +
    ode.out[,"A-2.1-0"]*(params['tau_d0_h0']+params['tau_d0_h3'])/2*(1-params['default_artatb']) +
    ode.out[,"A-3.1-0"]*(params['tau_d0_h3'])*(1-params['default_artatb']) +
    ode.out[,"A-1.1-1"]*(params['tau_d1'])*(1-params['default_artatb']) +
    ode.out[,"A-2.1-1"]*(params['tau_d1'])*(1-params['default_artatb']) +
    ode.out[,"A-3.1-1"]*(params['tau_d1'])*(1-params['default_artatb'])
  
  art.init <- art.init.cd4 + art.init.atb
  
  art.init.atbfrac <- art.init.atb/art.init
  
  # mean cd4 count at time of ART initiation (will be close to median given similar widths of cd4 compartments)
  cd4.artinit <- 
    (600*rowSums(ode.out[,c("S-1.1","Lf-1.1-0","Lf-1.1-1","Ls-1.1-0","Ls-1.1-1","A-1.1-0","A-1.1-1")])*params['chi2_h1'] + 
       316*rowSums(ode.out[,c("S-2.1","Lf-2.1-0","Lf-2.1-1","Ls-2.1-0","Ls-2.1-1", "A-2.1-0","A-2.1-1")])*params['chi2_h2'] + 
       50*rowSums(ode.out[,c("S-3.1","Lf-3.1-0","Lf-3.1-1","Ls-3.1-0","Ls-3.1-1", "A-3.1-0","A-3.1-1")])*params['chi2_h3'] + 
       600*ode.out[,"A-1.1-0"]*(params['tau_d0_h0']+params['tau_d0_h3'])/2*(1-params['default_artatb']) +
       316*ode.out[,"A-2.1-0"]*(params['tau_d0_h0']+params['tau_d0_h3'])/2*(1-params['default_artatb']) +
       50*ode.out[,"A-3.1-0"]*(params['tau_d0_h3'])*(1-params['default_artatb']) +
       600*ode.out[,"A-1.1-1"]*(params['tau_d1'])*(1-params['default_artatb']) +
       316*ode.out[,"A-2.1-1"]*(params['tau_d1'])*(1-params['default_artatb']) +
       50*ode.out[,"A-3.1-1"]*(params['tau_d1'])*(1-params['default_artatb']) )/art.init
  
  ## IPT initiations
  ipt.init <-
    rowSums(ode.out[,c("S-1.1","Lf-1.1-0","Lf-1.1-1","Ls-1.1-0","Ls-1.1-1")]*params['chi2_h1']*params['kappa']) +
    rowSums(ode.out[,c("S-2.1","Lf-2.1-0","Lf-2.1-1","Ls-2.1-0","Ls-2.1-1")]*params['chi2_h2']*params['kappa']) +
    rowSums(ode.out[,c("S-3.1","Lf-3.1-0","Lf-3.1-1","Ls-3.1-0","Ls-3.1-1")]*params['chi2_h3']*params['kappa']) +
    
    rowSums(ode.out[,c("A-1.1-0","A-1.1-1")])*
    ((1-params['art_screen'])*params['ipt_screen']*(1-params['default_acf'])*params['kappa']*params['chi2_h1'] + # IPT scrrening only
       params['art_screen']*(1-params['default_acf'])*params['kappa']*params['chi2_h1'] + # screening all starting IPT
       (1-params['ipt_screen'])*(1-params['art_screen'])*0) +
    
    rowSums(ode.out[,c("A-2.1-0","A-2.1-1")])*
    ((1-params['art_screen'])*params['ipt_screen']*(1-params['default_acf'])*params['kappa']*params['chi2_h2'] + # IPT scrrening only
       params['art_screen']*(1-params['default_acf'])*params['kappa']*params['chi2_h2'] + # screening all starting IPT
       (1-params['ipt_screen'])*(1-params['art_screen'])*0) +
    
    rowSums(ode.out[,c("A-3.1-0","A-3.1-1")])*
    ((1-params['art_screen'])*params['ipt_screen']*(1-params['default_acf'])*params['kappa']*params['chi2_h3'] + # those screened for IPT
       params['art_screen']*(1-params['default_acf'])*params['kappa']*params['chi2_h3'] + # with screening for all ART, same rate as with IPT screenign alone
       (1-params['ipt_screen'])*(1-params['art_screen'])*0) + # no actives receive effective IPT when not screened
    
    rowSums(ode.out[,c("S-1.2","S-2.2")]) * params['kappa_onart'] +
    rowSums(ode.out[,intersect(lu$lf.i,lu$h1art.i)]) *params['kappa_onart'] +
    rowSums(ode.out[,intersect(lu$lf.i,lu$h2art.i)]) *params['kappa_onart'] +
    rowSums(ode.out[,intersect(lu$ls.i,lu$h1art.i)]) *params['kappa_onart'] +
    rowSums(ode.out[,intersect(lu$ls.i,lu$h2art.i)]) *params['kappa_onart']
  
  # TB treatment initiations ~ notifications, per 100k
  # needs to include those also getting ART (i.e. simultaneously changing HIV-ART states)
  # and should also include those getting diagnosed through screening during ART initiation
  tb.notified <-
    (ode.out[,'A-0-0'])*params['tau_d0_h0'] + 
    rowSums(ode.out[,intersect(intersect(lu$a.i,lu$d0.i), lu$h1.i)])*(params['tau_d0_h0']+params['tau_d0_h3'])/2 + 
    rowSums(ode.out[,intersect(intersect(lu$a.i,lu$d0.i), lu$h2.i)])*(params['tau_d0_h0']+params['tau_d0_h3'])/2 + 
    rowSums(ode.out[,intersect(lu$a.i,lu$d1.i)])*params['tau_d1'] + 
    rowSums(ode.out[,intersect(lu$a.i, lu$h1nonart.i)])*params['chi2_h1']*params['sens_xpert']*
    (params['art_screen']+(1-params['art_screen'])*params['kappa']*params['ipt_screen']) +  
    rowSums(ode.out[,intersect(lu$a.i, lu$h2nonart.i)])*params['chi2_h2']*params['sens_xpert']*
    (params['art_screen']+(1-params['art_screen'])*params['kappa']*params['ipt_screen']) +  
    rowSums(ode.out[,intersect(lu$a.i, lu$h3nonart.i)])*params['chi2_h3']*params['sens_xpert']*
    (params['art_screen']+(1-params['art_screen'])*params['kappa']*params['ipt_screen']) 
  
  mdr.notified <-
    rowSums(ode.out[,intersect(lu$a.i,lu$d1.i)])*params['tau_d1'] + 
    (ode.out[,intersect(intersect(lu$a.i,lu$d1.i), lu$h1nonart.i)])*params['chi2_h1']*params['sens_xpert']*
    (params['art_screen']+(1-params['art_screen'])*params['kappa']*params['ipt_screen']) +  
    (ode.out[,intersect(intersect(lu$a.i,lu$d1.i), lu$h2nonart.i)])*params['chi2_h2']*params['sens_xpert']*
    (params['art_screen']+(1-params['art_screen'])*params['kappa']*params['ipt_screen']) +  
    (ode.out[,intersect(intersect(lu$a.i,lu$d1.i), lu$h3nonart.i)])*params['chi2_h3']*params['sens_xpert']*
    (params['art_screen']+(1-params['art_screen'])*params['kappa']*params['ipt_screen']) 
  
  prop.MDR <- mdr.notified/tb.notified
  
  ## TB incidence amoung those on ART, per 100k on ART
  tb.inc.onart <-
    1e5*(
      rowSums(ode.out[,intersect(lu$lf.i,lu$h1art.i)])*params['pi_h']*params['arttheta']*params['h1theta'] +
        rowSums(ode.out[,intersect(lu$lf.i,lu$h2art.i)])*params['pi_h']*params['arttheta']*params['h2theta'] +
        rowSums(ode.out[,intersect(lu$lf.i,lu$h3art.i)])*params['pi_h']*params['arttheta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h1art.i)])*params['epsilon_h']*params['arttheta']*params['h1theta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h2art.i)])*params['epsilon_h']*params['arttheta']*params['h2theta'] +
        rowSums(ode.out[,intersect(lu$ls.i,lu$h3art.i)])*params['epsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h1art.i)])*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h2art.i)])*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h3art.i)])*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['ipi_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h1art.i)])*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h2art.i)])*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h3art.i)])*params['ipi_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h1art.i)])*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h2art.i)])*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h3art.i)])*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['ipi_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h1art.i)])*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h2art.i)])*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h3art.i)])*params['ipi_h']*params['arttheta'])/
    rowSums(ode.out[,c(lu$h1art.i, lu$h2art.i, lu$h3art.i)])
  
  ## porportion of population on IPT
  onipt <-
    rowSums(ode.out[,c(lu$lipt.i,lu$liptf.i,lu$sipt.i)])/rowSums(ode.out)
  
  ## porportion of population after IPT
  postipt <-
    rowSums(ode.out[,lu$ipt2.i])/rowSums(ode.out)
  
  ## for table S1:
  tb.prev.eveript <-rowSums(
    ode.out[,intersect(lu$a.i,lu$allipt.i)])/ 
    rowSums(ode.out[,lu$allipt.i])
  
  tb.inc.eveript <- 1e5*
    (
      (ode.out[,intersect(intersect(lu$lipt.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h1art.i)])*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h2art.i)])*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt.i,lu$d1.i),lu$h3art.i)])*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['ipi_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h1art.i)])*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h2art.i)])*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$liptf.i,lu$d1.i),lu$h3art.i)])*params['ipi_h']*params['arttheta'] + 
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h1art.i)])*params['iepsilon_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h2art.i)])*params['iepsilon_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2.i,lu$d1.i),lu$h3art.i)])*params['iepsilon_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d0.i),lu$h1art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d0.i),lu$h2art.i)])*params['theta']*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d0.i),lu$h3art.i)])*params['theta']*params['ipi_h']*params['arttheta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h1art.i)])*params['ipi_h']*params['arttheta']*params['h1theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h2art.i)])*params['ipi_h']*params['arttheta']*params['h2theta'] +
        (ode.out[,intersect(intersect(lu$lipt2f.i,lu$d1.i),lu$h3art.i)])*params['ipi_h']*params['arttheta']
    )/rowSums(ode.out[,lu$allipt.i])
  
  
  ## does not include background mortality
  tb.mort.eveript <-  1e5*
    (rowSums(ode.out[,intersect(intersect(lu$a.i,lu$h1nonart.i),lu$allipt.i)]*(params['mu_tb_0']/2 + params['mu_tb_h']/2 + params['mu_h1'])) +
       rowSums(ode.out[,intersect(intersect(lu$a.i,lu$h2nonart.i),lu$allipt.i)]*(params['mu_tb_0']/2 + params['mu_tb_h']/2 + params['mu_h2'])) +
       rowSums(ode.out[,intersect(intersect(lu$a.i,lu$h3nonart.i),lu$allipt.i)]*(params['mu_tb_h'] + params['mu_h3'])) + 
       rowSums(ode.out[,intersect(intersect(lu$a.i,lu$h1art.i),lu$allipt.i)]*(params['mu_tb_0']/2 + params['mu_tb_h']/2 + params['mu_h1']))*params['arttheta'] +
       rowSums(ode.out[,intersect(intersect(lu$a.i,lu$h2art.i),lu$allipt.i)]*(params['mu_tb_0']/2 + params['mu_tb_h']/2 + params['mu_h2']))*params['arttheta'] +
       rowSums(ode.out[,intersect(intersect(lu$a.i,lu$h3art.i),lu$allipt.i)]*(params['mu_tb_h'] + params['mu_h3']))*params['arttheta'])/
    rowSums(ode.out[,lu$allipt.i])
  
  ##
  
  rc <- cbind(times,
              pop.size,
              tb.notified,
              tb.inc,
              tb.prev,
              hiv.inc,
              tb.mort.h0,
              tb.hiv.mort,
              tb.mort,
              tb.prev.hiv.pos,
              tb.prev.hiv.pos.nonart,
              tb.prev.on.art,
              tb.inc.hivpos,
              tb.mdr.inc,
              prop.MDR,
              ltbi.prob,
              ltbi.prob.in.hiv,
              hiv.prev,
              hiv.prev.in.TB,
              art.cov,
              art.cov.among.all.hiv,
              tb.inc.onart,
              ipt.elg.with.tb,
              onipt,
              postipt,
              ipt.init,
              art.init.cd4,
              art.init.atb,
              art.init,
              art.init.atbfrac,
              cd4.artinit,
              inf.rate.on.ipt,
              inf.rate.post.ipt,
              inf.rate.nonipt,
              tb.inc.on.ipt,
              tb.inc.post.ipt,
              tb.inc.nonipt,
              tb.prev.eveript,
              tb.inc.eveript,
              tb.mort.eveript
  )
  
  rc[is.na(rc)] <- 0
  
  return(rc)
}

##' Makes a list of indexes for different state groups
##' @return
##' @author asa
make.look.up.list <- function(state.names){
  
  if (missing(state.names)) state.names <- make.state.names()
  
  ## tb states
  s.i <- grep("S\\-",state.names,perl=T)
  lf.i <- grep("Lf\\-",state.names,perl=T)
  ls.i <- grep("Ls\\-",state.names,perl=T)
  a.i <- grep("A\\-",state.names,perl=T)
  sipt.i <- grep("SIPT\\-",state.names,perl=T)
  sipt2.i <- grep("SIPT2\\-",state.names,perl=T)
  lipt.i <- grep("LIPT\\-",state.names,perl=T)
  lipt2.i <- grep("LIPT2\\-",state.names,perl=T)
  liptf.i <- grep("LIPTf\\-",state.names,perl=T)
  lipt2f.i <- grep("LIPT2f\\-",state.names,perl=T)
  allipt.i <- c(sipt.i,sipt2.i,lipt.i,lipt2.i,liptf.i,lipt2f.i)
  ipt.i <- c(sipt.i,lipt.i,liptf.i)
  ipt2.i <- c(sipt2.i,lipt2.i,lipt2f.i)
  ## hiv compartment
  ## includes those on ART
  h1.i <- grep("\\-1\\.",state.names,perl=T) # hiv+  CD4 > 500
  h2.i <- grep("\\-2\\.",state.names,perl=T) # hiv+ 200 < CD4 < 500
  h3.i <- grep("\\-3\\.",state.names,perl=T) # hiv+ CD4 < 200
  h0.i <- setdiff(1:length(state.names),c(h1.i,h2.i,h3.i)) # hiv negative
  
  ## on art compartments
  h1art.i <- grep("\\-1\\.2",state.names,perl=T)
  h2art.i <- grep("\\-2\\.2",state.names,perl=T)
  h3art.i <- grep("\\-3\\.2",state.names,perl=T)
  allart.i <- c(h1art.i,h2art.i, h3art.i)
  h1nonart.i <- grep("\\-1\\.1",state.names,perl=T)
  h2nonart.i <- grep("\\-2\\.1",state.names,perl=T)
  h3nonart.i <- grep("\\-3\\.1",state.names,perl=T)
  nonart.i <- grep("\\.1",state.names,perl=T)
  
  ## drug resistant compartments
  d0.i <- grep("\\-0$",state.names,perl=T) # drug susceptible TB
  d0.i <- d0.i[-1] # remove S-0
  d1.i <- grep("\\-1$",state.names,perl=T) # MDR
  
  ret <- list(s.i=s.i,lf.i=lf.i,ls.i=ls.i,a.i=a.i,sipt.i=sipt.i,
              sipt2.i=sipt2.i,lipt.i=lipt.i,lipt2.i=lipt2.i, liptf.i=liptf.i,lipt2f.i=lipt2f.i,
              h0.i=h0.i,h1.i=h1.i,h2.i=h2.i,h3.i=h3.i,
              h1art.i=h1art.i,h2art.i=h2art.i,h3art.i=h3art.i,
              allart.i=allart.i,nonart.i=nonart.i,
              h1nonart.i=h1nonart.i,h2nonart.i=h2nonart.i,h3nonart.i=h3nonart.i,
              d0.i=d0.i,d1.i=d1.i,allipt.i=allipt.i,ipt.i=ipt.i,ipt2.i=ipt2.i)
  
  return(ret)
}

##' Takes the final state from a run and continue it
##' @param model ode model output
##' @param params params
##' @param extra.t how many more years
##' @param by time step
##' @return
##' @author asa
continue.run <- function(out,params,append=FALSE,extra.t=10,by=0.1,chg, print=T){
  lu <- make.look.up.list()
  if (missing(params))
    params <- load.params()
  
  params['iepsilon_0'] <- params['epsilon_0']
  params['iepsilon_h'] <- params['epsilon_h']
  params['ipi_0'] <- params['pi_0']
  params['ipi_h'] <- params['pi_h']
  
  if (missing(chg))
    chg <- make.chg.mat(params)
  
  new.state <- tail(out,1)[-1]
  rc <- ode(y=new.state,
            times=seq(tail(out,1)[1],tail(out,1)[1]+extra.t,by=.1),
            func=tb.dx.dt,
            parms=params,
            chg=chg,
            lu=lu)
  
  if (append)
    rc <- rbind(out,rc[-1,])
  
  state.names <- make.state.names()
  colnames(rc) <- c("time",state.names)
  if (print) print(final.stats(rc,params))
  
  return(rc)
}




plot.stats <- function(mod,params,sumstat){
  # generate summary stats if not supplied
  if (missing(sumstat))
    sumstat <- get.summary.stats(mod,params)
  
  times <- sumstat[,'times']
  par(mfrow=c(3,3),mar=c(2,2,2,2))
  plot(times,sumstat[,'tb.notified'],type="l")
  
  make.text.coords <- function(){
    c1 <- par("usr")[1] + (par("usr")[2] - par("usr")[1])*.5
    c2 <- par("usr")[3] + (par("usr")[4] - par("usr")[3])*.5
    c(c1,c2)
  }
  
  tmp <- make.text.coords()
  text(tmp[1],tmp[2],"TB Notifications")
  grid()
  plot(times,sumstat[,'tb.inc'],type="l")
  text(tmp[1],tmp[2],"TB Incidence")
  grid()
  plot(times,sumstat[,'tb.mdr.inc'],type="l")
  tmp <- make.text.coords()
  text(tmp[1],tmp[2],"MDR-TB Incidence")
  grid()
  plot(times,sumstat[,'tb.prev'],type="l")
  tmp <- make.text.coords()
  text(tmp[1],tmp[2],"TB Prevalance")
  grid()
  plot(times,sumstat[,'tb.prev.hiv.pos'],type="l")
  tmp <- make.text.coords()
  text(tmp[1],tmp[2],"%TB HIV+")
  grid()
  plot(times,sumstat[,'tb.mort'],type="l")
  tmp <- make.text.coords()
  text(tmp[1],tmp[2],"TB Mortality")
  grid()
  plot(times,sumstat[,'tb.hiv.mort'],type="l")
  tmp <- make.text.coords()
  text(tmp[1],tmp[2],"% TB Mortality HIV+")
  grid()
  #par(ask = TRUE)
  ## plot(times,sumstat[,'pop.size'],type="l")
  ## text(par("usr")[2]-3,par("usr")[3]+(par("usr")[4]-par("usr")[3])/1.2,"Population Size")
  ## grid()
  plot(times,sumstat[,'hiv.prev'],type="l")
  tmp <- make.text.coords()
  text(tmp[1],tmp[2],"HIV Prev.")
  grid()
  plot(times,sumstat[,'art.cov'],type="l")
  tmp <- make.text.coords()
  text(tmp[1],tmp[2],"ART Coverage \n (among all HIV+)")
  grid()
  #par(ask=FALSE)
}

##' Helper function to run the model with interventions
##' @param state
##' @param start.time
##' @param end.time
##' @param params
##' @param chg
##' @return
##' @author asa
run.model <- function(state,start.time,
                      end.time,
                      params,
                      chg,
                      to.eq,
                      by.time, 
                      projection=F, baseart=0.5){
  
  if(missing(by.time)) by.time <- 0.1
  if(missing(params)) params <- load.params()
  if(missing(chg)) chg <- make.chg.mat(params)
  
  lu <- make.look.up.list()
  
  rc <- ode(y=as.numeric(state),
            times=seq(start.time,end.time,by=by.time),
            func=tb.dx.dt,
            parms=params,
            chg=chg,
            lu=lu,
            projection=projection)
  
  if (to.eq){
    ## check the state.diffs
    state.diffs <- tail(apply(rc[,-1],2,diff),1)
    while(any(state.diffs > .1)){
      cat(sprintf("running to equilibrium, this may be a little while - new end.time = %s \n",end.time+20))
      rc <- ode(y=as.numeric(tail(rc[,-1],1)),
                times=seq(end.time,end.time+20,by=by.time),
                func=tb.dx.dt,
                parms=params,
                chg=chg,
                lu=lu,
                projection=projection, baseart=baseart)
      end.time <- end.time + 20
    }
  }
  
  rc <- get.colnames(rc)
  
  return(rc)
}

non.zero.cols <- function(ode.out){
  tmp <- which(colSums(ode.out)==0)
  if (length(tmp)==0) ode.out else ode.out[,-tmp]
  
}

##' simple objective function for fittting a single paramater
##' using squared difference loss
##' @title
##' @param param
##' @param param.name
##' @param params
##' @param stat.of.interest
##' @param value
##' @param lu
##' @param chg
##' @return
##' @author asa
obj.func <- function(param,param.name,params,stat.of.interest,value,
                     state.start,
                     lu,
                     chg){
  
  params[param.name] <- param
  if(missing(chg)) chg <- make.chg.mat(params)
  test <- ode(y=state.start,
              times=seq(0,1000,by=1),
              func=tb.dx.dt,
              parms=params,
              lu=lu,
              chg=chg)
  sum.stat <- tail(get.summary.stats(test,params),1)
  (sum.stat[,stat.of.interest] - value)^2
}

##' shows final summary stats for the last time step of a run
##' @param run output of ode()
##' @param params parameter vector
##' @param show.targets if TRUE displays (hardcoded) calibration targets for params
##' as an extra row below
##' @return
##' @author asa
final.stats <- function(run,params,show.targets=TRUE,targets){
  
  rc <-tail(get.summary.stats(run,params),1)[,-c(1:2)]
  if(show.targets){
    if (missing(targets)){
      targets <- rep(NA,length(rc))
      names(targets) <- names(rc)
      targets['tb.notified'] <- 1630 # MSF summary report and Lele
      targets['tb.prev'] <- 3200 # Claassens
      targets['prop.MDR'] <- .04   # weighted avergage of pretreated and new cases
      targets['tb.prev.on.art'] <- 0.1169 #0.026  #rangaka et al Lancet (under box 1)
      targets['hiv.prev'] <- 0.265 # adjusted antenatal prev survey data
      targets['tb.mort'] <- 500 # TB or HIV deaths in those with active TB, based on RSA 2013 WHO global report data
      targets['tb.mort.h0'] <- 50 # based on WHO 2011 RSA data in non-hivs
      targets['hiv.prev.in.TB'] <- 0.878
      targets['ltbi.prob.in.hiv'] <- 0.8 #0.612, Oni et al 2012 Plos One (IGRA) in HIV CLinics
      targets['ltbi.prob'] <- 0.8
      targets['art.cov'] <- 0.6
      targets['tb.inc'] <- 2000 # MSF report?
      targets['cd4.artinit'] <- 15
    }
    rc <- rbind(rc,targets,round((rc - targets)/targets,2))
    rownames(rc) <- c("modeled","target","pct.diff")
  }
  rc
}

##' quick and dirty function for running model from a two seeded cases
##' @title
##' @return
##' @author asa
run.and.print <- function(state.start,max.time=100,params){
  reload.source()
  ## starting state
  state.names <- make.state.names()
  if(missing(params)) params <- load.params()
  chg <- make.chg.mat(params)
  lu <- make.look.up.list()
  
  if(missing(state.start)){
    # one infected
    state.start <- rep(0,length(state.names))
    names(state.start) <- state.names
    
    N <- 100000
    state.start["A-0-0"] <- 1
    state.start["A-0-1"] <- 1
    state.start["S-0"] <- N- sum(state.start)
  }
  
  test <- ode(y=state.start,
              times=seq(0,max.time,by=.1),
              func=tb.dx.dt,
              parms=params,
              lu=lu,
              chg=chg)
  print(final.stats(test,params))
  plot.stats(test,params)
  return(test)
}


# ## gets effect size for progression to TB on IPT
# get.effect.size <- function(){
#     thetas <- seq(0,1,length=50)
#     sim.irrs <- sapply(thetas,function(x) get.irr(x,2, print=F))
#     approx(sim.irrs,thetas,.63)$y
# }
# 
##' function to make figure 2
##' @param sum.stat.int = array of states, dims=(time, state, run)
##' @param sum.stat.base = ""
##' @return
##' @author asa
make.fig2 <- function(sum.stat.int,
                      sum.stat.base, weights, gamma=1, filetag="", times = seq(0,5,by=.1)){
  
  if(missing(weights)) weights <- rep(1,dim(sum.stat.base)[1])
  
  cols.green <- brewer.pal(8,"Greens")
  cols.blue <- brewer.pal(8,"Blues")
  cols.orange <- brewer.pal(8,"Oranges")
  
  rows <- 0:length(times)
  ylim <- c(0,2800); ylim.m <- c(0,750); ylim.o <- c(0,0.03); times <- times
  layout(matrix(c(1,1,1,1,2,2,2,3,3),nrow=9))
  layout.show()
  par(mar=c(1,1,1,1),oma=c(2,6,1,1),mgp=c(1,.5,0),tck=-.01,las=1)
  plot(NA, xlim=c(0,5), ylim=ylim, xlab="",lty=2,axes=FALSE)
  polygon(c(times, rev(times)), c(apply(sum.stat.base[rows,'tb.inc',],1, wtd.quantile, weights, 0.975),rev(apply(sum.stat.base[rows,'tb.inc',],1, wtd.quantile, weights, 0.025))),
          col=AddAlpha(cols.green[3],.3),border=NA)
  polygon(c(times, rev(times)), c(apply(sum.stat.int[rows,'tb.inc',],1, wtd.quantile, weights, 0.975),rev(apply(sum.stat.int[rows,'tb.inc',],1, wtd.quantile, weights, 0.025))),
          col=AddAlpha(cols.green[5],0.3),border=NA)
  points(times, apply(sum.stat.base[rows,'tb.inc',],1, wtd.quantile, weights, 0.5), col=cols.green[4],type='l',lwd=2, lty=2)
  points(times, apply(sum.stat.int[rows,'tb.inc',],1, wtd.quantile, weights, 0.5), col=cols.green[6],type='l',lwd=2, lty=1)
  #     axis(3)
  axis(2,at=seq(ylim[1],ylim[2],by=300),cex.axis=0.95)
  axis(1,at=0:5,labels=rep("",6))
  #     grid()
  legend("bottomleft",c("",""),col=AddAlpha(cols.green[c(6,3)],0.3),lwd=10,lty=c(1,1),bty="n")
  legend("bottomleft",c("Intervention","Baseline"),col=cols.green[c(6,4)],lwd=2,lty=c(1,2),bty="n")
  mtext("TB Incidence\n(per 100,000 per year)", 2, cex=0.95, las=0, line=3)
  
  text(cex=1.05, x = 2.5, y=1900, paste0( round(wtd.quantile(colSums(to.khay.pop(sum.stat.base[1:51,"tb.inc",] - sum.stat.int[1:51,"tb.inc",])*0.1), weights=weights, probs = 0.5)),
                                          " TB cases\n(95% CrI ", round(wtd.quantile(colSums(to.khay.pop(sum.stat.base[1:51,"tb.inc",] - sum.stat.int[1:51,"tb.inc",])*0.1), weights=weights, probs = 0.025)),
                                          " - ", round(wtd.quantile(colSums(to.khay.pop(sum.stat.base[1:51,"tb.inc",] - sum.stat.int[1:51,"tb.inc",])*0.1), weights=weights, probs = 0.975)),
                                          ")\naverted over 5 years\nin population of 400,000"))
  text(cex=2, x=4.8, y=ylim[2]*0.9, "A")
  
  text(cex=1.05, x = 5, y=1800, paste0("Incidence reduction\nafter 5 years: ",sprintf("%.1f", round(100*wtd.quantile((sum.stat.base[51,'tb.inc',] - sum.stat.int[51,'tb.inc',])/sum.stat.base[51,'tb.inc',], weights, 0.5), 1)),"%\n(95% CrI ",
                                       sprintf("%.1f", round(100*wtd.quantile((sum.stat.base[51,'tb.inc',] - sum.stat.int[51,'tb.inc',])/sum.stat.base[51,'tb.inc',], weights, 0.025),1))," - ",
                                       sprintf("%.1f", round(100*wtd.quantile((sum.stat.base[51,'tb.inc',] - sum.stat.int[51,'tb.inc',])/sum.stat.base[51,'tb.inc',], weights, 0.975),1)),"%)"),
       pos=2)
  
  plot(NA, xlim=c(0,5), ylim=ylim.m, xlab="",lty=2,axes=FALSE)
  polygon(c(times, rev(times)), c(apply(sum.stat.base[rows,'tb.mort',],1, wtd.quantile, weights, 0.975),rev(apply(sum.stat.base[rows,'tb.mort',],1, wtd.quantile, weights, 0.025))),
          col=AddAlpha(cols.orange[3],0.3),border=NA)
  polygon(c(times, rev(times)), c(apply(sum.stat.int[rows,'tb.mort',],1, wtd.quantile, weights, 0.975),rev(apply(sum.stat.int[rows,'tb.mort',],1, wtd.quantile, weights, 0.025))),
          col=AddAlpha(cols.orange[6],0.3),border=NA)
  points(times, apply(sum.stat.base[rows,'tb.mort',],1, wtd.quantile, weights, 0.5), col=cols.orange[4],type='l',lwd=2, lty=2)
  points(times, apply(sum.stat.int[rows,'tb.mort',],1, wtd.quantile, weights, 0.5), col=cols.orange[6],type='l',lwd=2, lty=1)
  #     axis(3,at=0:5,labels=rep("",6))
  axis(2,at=seq(ylim.m[1],ylim.m[2],by=100),cex.axis=.95)
  axis(1,at=0:5,labels=rep("",6))
  text(cex=1.05, x = 2.5, y=300, paste0( round(wtd.quantile(colSums(to.khay.pop(sum.stat.base[1:51,"tb.mort",] - sum.stat.int[1:51,"tb.mort",])*0.1), weights=weights, probs = 0.5)),
                                         " TB deaths\n(95% CrI ", round(wtd.quantile(colSums(to.khay.pop(sum.stat.base[1:51,"tb.mort",] - sum.stat.int[1:51,"tb.mort",])*0.1), weights=weights, probs = 0.025)),
                                         " - ", round(wtd.quantile(colSums(to.khay.pop(sum.stat.base[1:51,"tb.mort",] - sum.stat.int[1:51,"tb.mort",])*0.1), weights=weights, probs = 0.975)),
                                         ")\naverted over 5 years\nin population of 400,000"))
  text(cex=2, x=4.8, y=ylim.m[2]*0.9, "B")
  text(cex=1.05, x = 5, y=270, paste0("Mortality reduction\nafter 5 years: ",sprintf("%.1f", round(100*wtd.quantile((sum.stat.base[51,'tb.mort',] - sum.stat.int[51,'tb.mort',])/sum.stat.base[51,'tb.mort',], weights, 0.5), 1)),"%\n(95% CrI ",
                                      sprintf("%.1f", round(100*wtd.quantile((sum.stat.base[51,'tb.mort',] - sum.stat.int[51,'tb.mort',])/sum.stat.base[51,'tb.mort',], weights, 0.025),1))," - ",
                                      sprintf("%.1f", round(100*wtd.quantile((sum.stat.base[51,'tb.mort',] - sum.stat.int[51,'tb.mort',])/sum.stat.base[51,'tb.mort',], weights, 0.975),1)),"%)"),
       pos=2)
  
  #     grid()s
  legend("bottomleft",c("",""),col=AddAlpha(cols.orange[c(6,3)],0.3),lwd=10,lty=c(1,1),bty="n")
  legend("bottomleft",c("Intervention","Baseline"),col=cols.orange[c(6,4)],lwd=2,lty=c(1,2),bty="n")
  mtext("TB Mortality\n(per 100,000 per year)", 2, cex=0.95, las=0, line=3)
  plot(NA, xlim=c(0,5), ylim=ylim.o, xlab="",lty=2,axes=FALSE)
  polygon(c(times, rev(times)), c(apply(sum.stat.base[rows,'onipt',],1, wtd.quantile, weights, 0.975),rev(apply(sum.stat.base[rows,'onipt',],1, wtd.quantile, weights, 0.025))),
          col=AddAlpha(cols.blue[3],0.3),border=NA)
  polygon(c(times, rev(times)), c(apply(sum.stat.int[rows,'onipt',],1, wtd.quantile, weights, 0.975),rev(apply(sum.stat.int[rows,'onipt',],1, wtd.quantile, weights, 0.025))),
          col=AddAlpha(cols.blue[6]),border=NA)
  points(times, apply(sum.stat.base[rows,'onipt',],1, wtd.quantile, weights, 0.5), col=cols.blue[4],type='l',lwd=2, lty=2)
  points(times, apply(sum.stat.int[rows,'onipt',],1, wtd.quantile, weights, 0.5), col=cols.blue[6],type='l',lwd=2, lty=1)
  #     axis(3,at=0:5,labels=rep("",6))
  axis(1,at=0:5)
  axis(2,cex.axis=.9, at=seq(0,0.03,by=0.005),labels = seq(0,3,by=0.5))
  #     grid()
  legend("topleft",c("",""),col=AddAlpha(cols.blue[c(6,3)],0.3),lwd=10,lty=c(1,1),bty="n")
  legend("topleft",c("Intervention","Baseline"),col=cols.blue[c(6,4)],lwd=2,lty=c(1,2),bty="n")
  mtext("Percent of\nPopulation\non IPT", 2, cex=0.95, las=0, line=3)
  text(cex=2, x=4.8, y=ylim.o[2]*0.9, "C")
  mtext("Years from Start of Intervention", 1, line=2, cex=0.9)
  
  #     mtext(ifelse(gamma==1, "12-month IPT","Continuous IPT"), outer=TRUE, line=-1, font=2)
}


make.fig2.alt <- function(sum.stat.int,
                          sum.stat.base, weights, gamma=1, filetag="", times = seq(0,5,by=.1)){
  
  if(missing(weights)) weights <- rep(1,dim(sum.stat.base)[1])
  
  cols.green <- brewer.pal(8,"Greens")
  cols.blue <- brewer.pal(8,"Blues")
  cols.orange <- brewer.pal(8,"Oranges")
  
  rows <- 0:length(times)
  ylim <- c(0,3000); ylim.m <- c(0,500); ylim.o <- c(0,0.03); times <- times
  layout(matrix(c(1,1,1,1,2,2,2,3,3,3),nrow=10))
  layout.show()
  par(mar=c(1,1,1,1),oma=c(2,4,1,1),mgp=c(1,.5,0),tck=-.01,las=1)
  plot(NA, xlim=c(0,5), ylim=ylim, xlab="",lty=2,axes=FALSE)
  
  #need cumsum fxn
  polygon(c(times, rev(times)), c(cumsum(apply((to.khay.pop(sum.stat.base[1:51,"tb.inc",] - sum.stat.int[1:51,"tb.inc",])*0.1), 1,wtd.quantile, weights=weights, probs = 0.975)), 
                                  rev(cumsum(apply((to.khay.pop(sum.stat.base[1:51,"tb.inc",] - sum.stat.int[1:51,"tb.inc",])*0.1), 1,wtd.quantile, weights=weights, probs = 0.025)))),
          col=AddAlpha(cols.green[5],.3),border=NA)
  points(times, cumsum(apply((to.khay.pop(sum.stat.base[1:51,"tb.inc",] - sum.stat.int[1:51,"tb.inc",])*0.1), 1,wtd.quantile, weights=weights, probs = 0.5)), col=cols.green[6],type='l',lwd=2, lty=1)
  points(times, rep(0,length(times)),col=cols.green[4],type='l',lwd=2, lty=2)
  #     axis(3)
  axis(2,at=seq(ylim[1],ylim[2],by=300),cex.axis=0.95)
  axis(1,at=0:5,labels=rep("",6))
  
  legend("topleft",c("",""),col=c(NA,AddAlpha(cols.green[5],0.3),NA),lwd=10,lty=1,bty="n")
  legend("topleft",c("Median projection", "95% credible interval", "Baseline"),col=c(cols.green[6],NA,cols.green[4]),lwd=2,lty=c(1,1,2),bty="n")
  mtext("TB cases averted", 2, cex=0.95, las=0, line=3)
  
  text(cex=2, x=0.2, y=ylim[2]*0.5, "A")
  
  text(cex=1.05, x = 5, y=ylim[2]*0.2, paste0("Median incidence rate reduction\nafter 5 years: ",sprintf("%.1f", round(100*wtd.quantile((sum.stat.base[51,'tb.inc',] - sum.stat.int[51,'tb.inc',])/sum.stat.base[51,'tb.inc',], weights, 0.5), 1)),"%\n(95% CrI ",
                                              sprintf("%.1f", round(100*wtd.quantile((sum.stat.base[51,'tb.inc',] - sum.stat.int[51,'tb.inc',])/sum.stat.base[51,'tb.inc',], weights, 0.025),1))," - ",
                                              sprintf("%.1f", round(100*wtd.quantile((sum.stat.base[51,'tb.inc',] - sum.stat.int[51,'tb.inc',])/sum.stat.base[51,'tb.inc',], weights, 0.975),1)),"%)"),
       pos=2)
  
  plot(NA, xlim=c(0,5), ylim=ylim.m, xlab="",lty=2,axes=FALSE)
  
  
  polygon(c(times, rev(times)), c(cumsum(apply((to.khay.pop(sum.stat.base[1:51,"tb.mort",] - sum.stat.int[1:51,"tb.mort",])*0.1), 1,wtd.quantile, weights=weights, probs = 0.975)), 
                                  rev(cumsum(apply((to.khay.pop(sum.stat.base[1:51,"tb.mort",] - sum.stat.int[1:51,"tb.mort",])*0.1), 1,wtd.quantile, weights=weights, probs = 0.025)))),
          col=AddAlpha(cols.orange[5],.3),border=NA)
  points(times, cumsum(apply((to.khay.pop(sum.stat.base[1:51,"tb.mort",] - sum.stat.int[1:51,"tb.mort",])*0.1), 1,wtd.quantile, weights=weights, probs = 0.5)), col=cols.orange[6],type='l',lwd=2, lty=1)
  points(times, rep(0,length(times)), col=cols.orange[4],type='l',lwd=2, lty=2)
  
  #     axis(3)
  axis(2,at=seq(ylim[1],ylim[2],by=100),cex.axis=0.95)
  axis(1,at=0:5,labels=rep("",6))
  
  legend("topleft",c("",""),col=c(NA,AddAlpha(cols.orange[5],0.3), NA),lwd=10,lty=1,bty="n")
  legend("topleft",c("Median projection", "95% credible interval", "Baseline"),col=c(cols.orange[6],NA,cols.orange[4]),lwd=2,lty=c(1,1,2),bty="n")
  mtext("TB deaths averted", 2, cex=0.95, las=0, line=3)
  
  text(cex=2, x=0.2, y=ylim.m[2]*0.5, "B")
  
  
  text(cex=1.05, x = 5, y=400, paste0("Median mortality rate reduction\nafter 5 years: ",sprintf("%.1f", round(100*wtd.quantile((sum.stat.base[51,'tb.mort',] - sum.stat.int[51,'tb.mort',])/sum.stat.base[51,'tb.mort',], weights, 0.5), 1)),"%\n(95% CrI ",
                                      sprintf("%.1f", round(100*wtd.quantile((sum.stat.base[51,'tb.mort',] - sum.stat.int[51,'tb.mort',])/sum.stat.base[51,'tb.mort',], weights, 0.025),1))," - ",
                                      sprintf("%.1f", round(100*wtd.quantile((sum.stat.base[51,'tb.mort',] - sum.stat.int[51,'tb.mort',])/sum.stat.base[51,'tb.mort',], weights, 0.975),1)),"%)"),
       pos=2)
  
  plot(NA, xlim=c(0,5), ylim=ylim.o, xlab="",lty=2,axes=FALSE)
  polygon(c(times, rev(times)), c(apply(sum.stat.base[rows,'onipt',],1, wtd.quantile, weights, 0.975),rev(apply(sum.stat.base[rows,'onipt',],1, wtd.quantile, weights, 0.025))),
          col=AddAlpha(cols.blue[3],0.3),border=NA)
  polygon(c(times, rev(times)), c(apply(sum.stat.int[rows,'onipt',],1, wtd.quantile, weights, 0.975),rev(apply(sum.stat.int[rows,'onipt',],1, wtd.quantile, weights, 0.025))),
          col=AddAlpha(cols.blue[6]),border=NA)
  points(times, apply(sum.stat.base[rows,'onipt',],1, wtd.quantile, weights, 0.5), col=cols.blue[4],type='l',lwd=2, lty=2)
  points(times, apply(sum.stat.int[rows,'onipt',],1, wtd.quantile, weights, 0.5), col=cols.blue[6],type='l',lwd=2, lty=1)
  #     axis(3,at=0:5,labels=rep("",6))
  axis(1,at=0:5)
  axis(2,cex.axis=.9, at=seq(0,0.03,by=0.005),labels = seq(0,3,by=0.5))
  #     grid()
  legend("topleft",c("",""),col=c(NA,AddAlpha(cols.blue[5],0.3), NA),lwd=10,lty=1,bty="n")
  legend("topleft",c("Median projection", "95% credible interval", "Baseline"),col=c(cols.blue[6],NA,cols.blue[4]),lwd=2,lty=c(1,1,2),bty="n")
  mtext("Percent of Population\non IPT", 2, cex=0.95, las=0, line=2)
  text(cex=2, x=0.2, y=ylim.o[2]*0.5, "C")
  mtext("Years from Start of Intervention", 1, line=2, cex=0.9)
  
  #     mtext(ifelse(gamma==1, "12-month IPT","Continuous IPT"), outer=TRUE, line=-1, font=2)
}



## figure 2 comparing continous to 12-month
make.fig2.sup <- function(filetag="", times=seq(0,5,by=0.1)){
  load(paste0(savelocation,"generated_data/sum.stat",filetag,".Rdata"))
  
  ## set up colors
  cols.green <- brewer.pal(8,"Greens")
  cols.blue <- brewer.pal(8,"Blues")
  cols.orange <- brewer.pal(8,"Oranges")
  
  ylim <- c(2000,2500); ylim.m <- c(400,600); ylim.o <- c(0,0.1); times <- seq(0,5,by=.1)
  ## plot
  layout(matrix(c(1,1,2,3),nrow=4))
  par(mar=c(2,1,1,1),oma=c(2,6,1,1),mgp=c(1,.5,0),tck=-.02,las=1)
  plot(times, apply(sum.stat.base[1:length(times),'tb.inc',],1, wtd.quantile, weights=weights, probs=0.5), ylim=ylim,type="l",col=cols.green[3],lwd=2,
       xlab="",ylab="",lty=2,axes=FALSE)
  axis(3)
  axis(2,at=seq(ylim[1],ylim[2],by=(ylim[2]-ylim[1])/10),cex.axis=.9)
  axis(1,at=0:5,labels=rep("",6))
  mtext("Median projected annual\nTB Incidence (per 100,000)", 2, cex=0.8, las=0, line=3)
  lines(times, apply(sum.stat.int[1:length(times),'tb.inc',],1, wtd.quantile, weights=weights, probs=0.5),col=cols.green[6],lwd=2)
  lines(times, apply(sum.stat.cont[1:length(times),'tb.inc',],1, wtd.quantile, weights=weights, probs=0.5),col=AddAlpha(cols.green[6],0.3),lwd=2)
  # grid()
  legend("bottomleft",c("Baseline", "12-month IPT","Continuous IPT"),
         col=c(cols.green[3],cols.green[6], AddAlpha(cols.green[6],.3)),lwd=2,lty=c(1,1,2),bty="n")
  text(cex=2, x=4.8, y=ylim[1] + (ylim[2]-ylim[1])*0.8, "A")
  
  plot(times, apply(sum.stat.base[1:length(times),'tb.mort',],1, wtd.quantile, weights=weights, probs=0.5),
       ylim=ylim.m,type="l",col=cols.orange[3],lwd=2,
       xlab="",ylab="",lty=2,axes=F)
  lines(times, apply(sum.stat.int[1:length(times),'tb.mort',],1, wtd.quantile, weights=weights, probs=0.5),
        col=cols.orange[6],lwd=2)
  lines(times, apply(sum.stat.cont[1:length(times),'tb.mort',],1, wtd.quantile, weights=weights, probs=0.5),
        col=AddAlpha(cols.orange[6],.3),lwd=2)
  #     axis(3,at=0:5,labels=rep("",6))
  axis(2,at=seq(ylim.m[1],ylim.m[2],by=(ylim.m[2]-ylim.m[1])/5),cex.axis=.9)
  axis(1,at=0:5,labels=rep("",6))
  mtext("Median projected annual\nTB Mortality (per 100,000)", 2, cex=0.8, las=0, line=3)
  # grid()
  legend("bottomleft",c("Baseline","12-month IPT","Continuous IPT"),
         col=c(cols.orange[3], cols.orange[6],AddAlpha(cols.orange[6],.3)),lwd=2,lty=c(2,1,1),bty="n")
  text(cex=2, x=4.8, y=ylim.m[1] + (ylim.m[2]-ylim.m[1])*0.8, "B")
  
  plot(times, apply(sum.stat.base[1:length(times),'onipt',],1, wtd.quantile, weights=weights, probs=0.5),
       ylim=ylim.o,type="h",col=cols.blue[3],lwd=2,
       xlab="",ylab="",lty=2,axes=FALSE)
  #     axis(3,at=0:5,labels=rep("",6))
  axis(1,at=0:5)
  axis(2,cex.axis=.9)
  lines(times, apply(sum.stat.int[1:length(times),'onipt',],1, wtd.quantile, weights=weights, probs=0.5),col=cols.blue[6],type="h",lwd=2)
  lines(times, apply(sum.stat.cont[1:length(times),'onipt',],1, wtd.quantile, weights=weights, probs=0.5),col=AddAlpha(cols.blue[6],.4),type="h",lwd=2)
  legend("topleft",c("12-month IPT", "Continuous IPT"),
         col=c(cols.blue[6],AddAlpha(cols.blue[6],.4)),lwd=2,lty=c(1,1),bty="n")
  text(cex=2, x=4.8, y=ylim.o[2]*0.9, "C")
  
  mtext("Median projected proportion\nof population on IPT", 2, cex=0.8, las=0, line=3)
  mtext("Years after start of intervention", 1, cex=0.9, las=0, line=2)
  
  # grid()
}


##' Takes output from the model and transforms it to a
##' raw number in the population of Khayelitsha
##' @title
##' @param data.point
##' @param pct.pop.over15 - pct populatino over 15 default value from 2001 census
##' @return
##' @author asa
to.khay.pop <- function(data.point,prop.pop.over15=(1-.296)){
  data.point*391749/1e5*prop.pop.over15
}




## get the proportion by HIV status
props.by.hiv <- function(start.state){
  lu <- make.look.up.list()
  
  c(sum(start.state[lu$h0.i]),
    sum(start.state[lu$h1.i]),
    sum(start.state[lu$h2.i]),
    sum(start.state[lu$h3.i]))/1e5
}


## function to make one of the summary tables for the paper
make.table.2 <- function(un.ss.dist,un.param.dist, weights, q=c(0.5,0.025,0.975)){
  states <- un.ss.dist
  colnames(states) <- make.state.names()
  
  sum.dist <- t(sapply(1:nrow(states),function(z) {
    get.summary.stats(
      rbind(cbind("time"=1,states[z,,drop=F]),
            cbind("time"=2,states[z,,drop=F])),
      un.param.dist[z,])[1,]}))
  
  ## now for the stats
  rc <- rbind(wtd.quantile(sum.dist[,'tb.notified'],weights,q),
              wtd.quantile(sum.dist[,'tb.prev'],weights,q),
              wtd.quantile(sum.dist[,'tb.mort'],weights,q),
              wtd.quantile(sum.dist[,'tb.hiv.mort'],weights,q),
              wtd.quantile(sum.dist[,'prop.MDR'],weights,q),
              wtd.quantile(sum.dist[,'hiv.prev'],weights,q),
              wtd.quantile(sum.dist[,'tb.inc'],weights,q),
              wtd.quantile(sum.dist[,'hiv.prev.in.TB'],weights,q),
              wtd.quantile(sum.dist[,'ltbi.prob.in.hiv'],weights,q),
              wtd.quantile(sum.dist[,'tb.prev.hiv.pos'],weights,q),
              wtd.quantile(sum.dist[,'art.cov'],weights,q),
              wtd.quantile(sum.dist[,'cd4.artinit'],weights,q),
              wtd.quantile(sum.dist[,'ipt.elg.with.tb'],weights,q),
              wtd.quantile(sum.dist[,'tb.inc.onart'],weights,q),
              wtd.quantile(sum.dist[,'tb.prev.on.art'],weights,q))
  
  rownames(rc) <- c("TB Notifications","TB prevalence", "TB mortality", "TB-HIV mortality",
                    "MDR proportion", "HIV prevalence", "TB Incidence", "HIV prevalence in TB", "LTBI prevalence in HIV",
                    "TB prevalence among HIV+", "ART coverage", "Median CD4 at ART initiation", "TB prevalence in IPT eligible", "TB incidence among those on ART", "TB prevalence among those on ART") 
  
  rc[c(5:6,8:11,13,15),] <- rc[c(5:6,8:11,13,15),]*100
  rc <- data.frame("estimate"=round(rc[,1],1),"CI-95%"=apply(apply(round(rc[,2:length(q)],2),1,paste0),1,paste0))
  # paste0("(",round(x[1],1),",",round(x[2],1),")")))
  
  return(rc)
  
}


make.main.tables.and.figures <- function(filetag="", plot=T,
                                         save.pdf=FALSE, ltd=FALSE, art=FALSE, notheta=FALSE,
                                         pdf.name, rows=1:51)
{
  
  load(paste0("generated_data/sum.stat",filetag,".Rdata"))
  if (ltd) 
  {
    load(paste0("generated_data/limited.sum.stat",filetag,".Rdata"))
    sum.stat.cont <- sum.stat.ltd
  }
  if (art)
  {
    load(paste0("generated_data/sensis.sum.stat",filetag,".Rdata"))
    sum.stat.cont <- sum.stat.art
    sum.stat.base <- sum.stat.art.base
  }
  if (notheta)
  {
    load(paste0("generated_data/sensis.sum.stat",filetag,".Rdata"))
    sum.stat.cont <- sum.stat.notheta
  }
  
  
  ## calculate cases and deaths averted, and make table
  q <- c(0.025,0.5,0.975)
  ## cases and deaths averted in Khay
  mort.base <-to.khay.pop(sum.stat.base[rows,"tb.mort",])*0.1 #0.1 is time interval
  mort.int <- to.khay.pop(sum.stat.int[rows,"tb.mort",])*0.1
  mort.cont <- to.khay.pop(sum.stat.cont[rows,"tb.mort",])*0.1
  cases.base <- to.khay.pop(sum.stat.base[rows,"tb.inc",])*0.1
  cases.int <- to.khay.pop(sum.stat.int[rows,"tb.inc",])*0.1
  cases.cont <- to.khay.pop(sum.stat.cont[rows,"tb.inc",])*0.1
  red.cases.int <- cases.base-cases.int
  red.cases.cont <- cases.base-cases.cont
  red.mort.int <- mort.base-mort.int
  red.mort.cont <- mort.base-mort.cont
  
  ## lets calculate the person time on IPT 
  ipt.int <- to.khay.pop(sum.stat.int[rows,"onipt",]*100000)*0.1 #0.1time, 1e5 population
  ipt.cont <- to.khay.pop(sum.stat.cont[rows,"onipt",]*100000)*0.1 #0.1time, 1e5 population
  
  table3.array <- array(c(
    wtd.quantile(to.khay.pop(colSums(sum.stat.int[rows,'ipt.init',]*.1)),weights,q),## number treated with IPT
    wtd.quantile(to.khay.pop(colSums(sum.stat.cont[rows,'ipt.init',]*.1)),weights,q),## number treated with IPT
    wtd.quantile(colSums(ipt.int),weights,q), ## person years on IPT
    wtd.quantile(colSums(ipt.cont),weights,q), ## person years on IPT
    wtd.quantile(colSums(red.cases.int),weights,q), ## tb cases averted
    wtd.quantile(colSums(red.cases.cont),weights,q), ## tb cases averted
    wtd.quantile((apply(cases.base,2,sum) - apply(cases.int,2,sum))/apply(cases.base,2,sum),weights,q), # % cum cases averted
    wtd.quantile((apply(cases.base,2,sum) - apply(cases.cont,2,sum))/apply(cases.base,2,sum),weights,q), # % cum cases averted
    wtd.quantile((sum.stat.base[max(rows),'tb.inc',] - sum.stat.int[max(rows),'tb.inc',])/sum.stat.base[max(rows),'tb.inc',],weights,q), ## % reduction in incidenceafter 5-years
    wtd.quantile((sum.stat.base[max(rows),'tb.inc',] - sum.stat.cont[max(rows),'tb.inc',])/sum.stat.base[max(rows),'tb.inc',],weights,q), ## % reduction in incidenceafter 5-years
    wtd.quantile(to.khay.pop(colSums(sum.stat.int[rows,'ipt.init',]*.1)/colSums(red.cases.int)),weights,q), ## nnt case
    wtd.quantile(to.khay.pop(colSums(sum.stat.cont[rows,'ipt.init',]*.1)/colSums(red.cases.cont)),weights,q), ## nnt case
    wtd.quantile(colSums(ipt.int)/colSums(red.cases.int),weights,q), ##person years per case averted
    wtd.quantile(colSums(ipt.cont)/colSums(red.cases.cont),weights,q), ##person years per case averted
    wtd.quantile(colSums(red.mort.int),weights,q), ## tb deaths averted
    wtd.quantile(colSums(red.mort.cont),weights,q), ## tb deaths averted
    wtd.quantile((apply(mort.base,2,sum) - apply(mort.int,2,sum))/apply(mort.base,2,sum),weights,q),
    wtd.quantile((apply(mort.base,2,sum) - apply(mort.cont,2,sum))/apply(mort.base,2,sum),weights,q),
    wtd.quantile((sum.stat.base[max(rows),'tb.mort',] - sum.stat.int[max(rows),'tb.mort',])/sum.stat.base[max(rows),'tb.mort',],weights,q), ##% reductions in mortalty rate after 5-years
    wtd.quantile((sum.stat.base[max(rows),'tb.mort',] - sum.stat.cont[max(rows),'tb.mort',])/sum.stat.base[max(rows),'tb.mort',],weights,q), ##% reductions in mortalty rate after 5-years
    wtd.quantile(to.khay.pop(colSums(sum.stat.int[rows,'ipt.init',]*.1)/colSums(red.mort.int)),weights,q), ## nnt death
    wtd.quantile(to.khay.pop(colSums(sum.stat.cont[rows,'ipt.init',]*.1)/colSums(red.mort.cont)),weights,q), ## nnt death
    wtd.quantile(colSums(ipt.int)/colSums(red.mort.int),weights,q), ##person years per death averted
    wtd.quantile(colSums(ipt.cont)/colSums(red.mort.cont),weights,q)), ##person years per death averted
    dim=c(6,12))
  
  tab3.names <- c("Number of people treated with IPT",
                  "Person-years on IPT",
                  "TB cases averted",
                  "% cumulative TB cases averted over 5 years", 
                  "Percent reduction in incidence rate after 5 years", 
                  "Number needed to treat to prevent one TB case",
                  "Number of person-years needed on IPT to prevent one TB case",
                  "TB deaths averted",
                  "% cumulative TB deaths averted over 5 years", 
                  "Percent reduction in mortality after 5 years",
                  "Number needed to treat to prevent one TB death",
                  "Number of person-years needed on IPT to prevent one TB death")
  
  colnames(table3.array) <- tab3.names; rownames(table3.array)<-rep(q,times=2) 
  
  if(plot) {
    ## make figure 2
    if(save.pdf){
      to.pdf(make.fig2(sum.stat.int,sum.stat.base,weights,filetag=filetag)
             ,paste0(ifelse(is.null(pdf.name),"plots/fig2",pdf.name),ifelse(gamma==0,"-sup",""),".pdf"),width=6,height=4)
      pdf(file="plots/fig2-sup.pdf",width=6,height=5); make.fig2.sup(filetag=filetag); dev.off()
    } else {
      make.fig2(sum.stat.int,sum.stat.base, weights, filetag=filetag)
    }}
  
  return(table3.array)
}

