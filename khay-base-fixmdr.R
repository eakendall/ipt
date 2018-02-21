## ---------------------------------- ##
## core functions for khay-IPT  model ##
## ---------------------------------- ##

##' This function makes a delta matrix for
##' @param params named list with parameters (assuming this is the augmented version with pseudo params)
##' @return array with transition rates between states
##' @author Andrew Azman
make.chg.mat <- function(params){
  
  params['iepsilon_0'] <- params['epsilon_0']
  params['iepsilon_h'] <- params['epsilon_h']
  params['ipi_0'] <- params['pi_0']
  params['ipi_h'] <- params['pi_h']
  
  # get state names
  state.names <- make.state.names()
  
  chg<-array(0,dim=c(length(state.names),length(state.names)))
  rownames(chg) <- colnames(chg) <- state.names
  
  ## -------------------- ##
  ## some helpful indexes ##
  ## -------------------- ##
  lu <- make.look.up.list(state.names)
  
  ## --------- ##
  ## mortality ##
  ## --------- ##
  
  ## general mortality:
  ## assign baseline mortality to diagonal elements of matrix
  diag(chg) <- diag(chg) + params['mu_d']
  
  ## hiv mortality
  diag(chg)[lu$h1art.i] <- diag(chg)[lu$h1art.i] + params['mu_h1']*params['arttheta']
  diag(chg)[lu$h2art.i] <- diag(chg)[lu$h2art.i] + params['mu_h2']*params['arttheta']
  diag(chg)[lu$h1.i] <- diag(chg)[lu$h1.i] + params['mu_h1']
  diag(chg)[lu$h2.i] <- diag(chg)[lu$h2.i] + params['mu_h2']
  diag(chg)[lu$h3.i] <- diag(chg)[lu$h3.i] + params['mu_h3']
  
  ## TB mortality, in those without HIV:
  diag(chg)[intersect(lu$a.i,lu$h0.i)] <- diag(chg)[intersect(lu$a.i,lu$h0.i)] + params['mu_tb_0']
  ## TB mortality, in those with HIV:
  diag(chg)[intersect(lu$a.i,lu$h1.i)] <- diag(chg)[intersect(lu$a.i,lu$h1.i)] + params['mu_tb_0']/2 + params['mu_tb_h']/2
  diag(chg)[intersect(lu$a.i,lu$h2.i)] <- diag(chg)[intersect(lu$a.i,lu$h2.i)] + params['mu_tb_0']/2 + params['mu_tb_h']/2
  diag(chg)[intersect(lu$a.i,lu$h3.i)] <- diag(chg)[intersect(lu$a.i,lu$h3.i)] + params['mu_tb_h']
  
  ## will take care of birth in dynamic updates
  
  ####################
  ## TB Transitions ##
  ####################
  
  ## new infections will be taken care of
  ## in dynamic updates
  
  ## ------------------------------------------------ ##
  ## fast latent to active disease (fast progression) ##
  ## ------------------------------------------------ ##
  
  ## hiv negatives
  chg[cbind(intersect(lu$lf.i,lu$h0.i),intersect(lu$a.i,lu$h0.i))] <-
    chg[cbind(intersect(lu$lf.i,lu$h0.i),intersect(lu$a.i,lu$h0.i))] +
    params['pi_0']
  
  ## CD4<200's
  chg[cbind(intersect(lu$lf.i,lu$h3nonart.i),intersect(lu$a.i,lu$h3nonart.i))] <-
    chg[cbind(intersect(lu$lf.i,lu$h3nonart.i),intersect(lu$a.i,lu$h3nonart.i))] +
    params['pi_h']

  chg[cbind(intersect(lu$lf.i,lu$h3art.i),intersect(lu$a.i,lu$h3art.i))] <-
    chg[cbind(intersect(lu$lf.i,lu$h3art.i),intersect(lu$a.i,lu$h3art.i))] +
    params['pi_h']*params['arttheta']
  
    
  chg[cbind(intersect(lu$lf.i,lu$h1nonart.i),intersect(lu$a.i,lu$h1nonart.i))] <-
    chg[cbind(intersect(lu$lf.i,lu$h1nonart.i),intersect(lu$a.i,lu$h1nonart.i))] +
    params['pi_h']*params['h1theta']
  
  chg[cbind(intersect(lu$lf.i,lu$h1art.i),intersect(lu$a.i,lu$h1art.i))] <-
    chg[cbind(intersect(lu$lf.i,lu$h1art.i),intersect(lu$a.i,lu$h1art.i))] +
    params['pi_h']*params['h1theta']*params['arttheta']
  
  chg[cbind(intersect(lu$lf.i,lu$h2nonart.i),intersect(lu$a.i,lu$h2nonart.i))] <-
    chg[cbind(intersect(lu$lf.i,lu$h2nonart.i),intersect(lu$a.i,lu$h2nonart.i))] +
    params['pi_h']*params['h2theta']
  
  chg[cbind(intersect(lu$lf.i,lu$h2art.i),intersect(lu$a.i,lu$h2art.i))] <-
    chg[cbind(intersect(lu$lf.i,lu$h2art.i),intersect(lu$a.i,lu$h2art.i))] +
    params['pi_h']*params['h2theta']*params['arttheta']
  
  ## -------------------------- ##
  ## fast latent to slow latent ##
  ## -------------------------- ##
  
  # chg[cbind(lu$lf.i,lu$ls.i)] <-
  #   chg[cbind(lu$lf.i,lu$ls.i)] + params['sigma'] # will reduce for reinfection in dxdt function (foi dependent)
  chg[cbind(lu$liptf.i,lu$lipt.i)] <-
    chg[cbind(lu$liptf.i,lu$lipt.i)] + params['sigma']
  chg[cbind(lu$lipt2f.i,lu$lipt2.i)] <-
    chg[cbind(lu$lipt2f.i,lu$lipt2.i)] + params['sigma']
  
  
  ## --------------------- ##
  ## slow latent to active ##
  ## --------------------- ##
  
  ## hiv negative
  chg[cbind(intersect(lu$ls.i,lu$h0.i),intersect(lu$a.i,lu$h0.i))] <-
    chg[cbind(intersect(lu$ls.i,lu$h0.i),intersect(lu$a.i,lu$h0.i))] +
    params['epsilon_0']
  
  ## CD4 <200
  chg[cbind(intersect(lu$ls.i,lu$h3nonart.i),intersect(lu$a.i,lu$h3nonart.i))] <-
    chg[cbind(intersect(lu$ls.i,lu$h3nonart.i),intersect(lu$a.i,lu$h3nonart.i))] +
    params['epsilon_h']
  
  chg[cbind(intersect(lu$ls.i,lu$h3art.i),intersect(lu$a.i,lu$h3art.i))] <-
    chg[cbind(intersect(lu$ls.i,lu$h3art.i),intersect(lu$a.i,lu$h3art.i))] +
    params['epsilon_h']*params['arttheta']
  
  chg[cbind(intersect(lu$ls.i,lu$h1nonart.i),intersect(lu$a.i,lu$h1nonart.i))] <-
    chg[cbind(intersect(lu$ls.i,lu$h1nonart.i),intersect(lu$a.i,lu$h1nonart.i))] +
    params['epsilon_h']*params['h1theta']

  chg[cbind(intersect(lu$ls.i,lu$h1art.i),intersect(lu$a.i,lu$h1art.i))] <-
    chg[cbind(intersect(lu$ls.i,lu$h1art.i),intersect(lu$a.i,lu$h1art.i))] +
    params['epsilon_h']*params['h1theta']*params['arttheta']
  
  chg[cbind(intersect(lu$ls.i,lu$h2nonart.i),intersect(lu$a.i,lu$h2nonart.i))] <-
    chg[cbind(intersect(lu$ls.i,lu$h2nonart.i),intersect(lu$a.i,lu$h2nonart.i))] +
    params['epsilon_h']*params['h2theta']

  chg[cbind(intersect(lu$ls.i,lu$h2art.i),intersect(lu$a.i,lu$h2art.i))] <-
    chg[cbind(intersect(lu$ls.i,lu$h2art.i),intersect(lu$a.i,lu$h2art.i))] +
    params['epsilon_h']*params['h2theta']*params['arttheta']
  
    
  ## --------------------- ##
  ## active to latent slow = TB treatment ##
  ## *now with potential for ART and IPT following TB treatment* ##
  ## (but not modeling time on treatment, so IPT will seem shorter) ##
  ## --------------------- ##
  
  ## drug susceptible
  #HIV-
  chg[cbind(intersect(intersect(lu$a.i,lu$d0.i),lu$h0.i),intersect(intersect(lu$ls.i,lu$d0.i),lu$h0.i))] <-
    chg[cbind(intersect(intersect(lu$a.i,lu$d0.i),lu$h0.i),intersect(intersect(lu$ls.i,lu$d0.i),lu$h0.i))]+params['tau_d0_h0']
   
  # HIV+
  # if were already on ART, just treat the TB
  chg[cbind(intersect(intersect(lu$a.i,lu$d0.i),lu$h3art.i),intersect(intersect(lu$ls.i,lu$d0.i),lu$h3art.i))] <-
    chg[cbind(intersect(intersect(lu$a.i,lu$d0.i),lu$h3art.i),intersect(intersect(lu$ls.i,lu$d0.i),lu$h3art.i))]+params['tau_d0_h3']

  chg[cbind(intersect(intersect(lu$a.i,lu$d0.i),lu$h1art.i),intersect(intersect(lu$ls.i,lu$d0.i),lu$h1art.i))] <-
    chg[cbind(intersect(intersect(lu$a.i,lu$d0.i),lu$h1art.i),intersect(intersect(lu$ls.i,lu$d0.i),lu$h1art.i))]+ mean(params[c('tau_d0_h3','tau_d0_h0')])
  
  chg[cbind(intersect(intersect(lu$a.i,lu$d0.i),lu$h2art.i),intersect(intersect(lu$ls.i,lu$d0.i),lu$h2art.i))] <-
    chg[cbind(intersect(intersect(lu$a.i,lu$d0.i),lu$h2art.i),intersect(intersect(lu$ls.i,lu$d0.i),lu$h2art.i))]+ mean(params[c('tau_d0_h3','tau_d0_h0')])

  # if not already on ART, treat the TB, start ART with some probability, and if starting ART, start IPT with probability kappa (e.g. 0.85)
  chg["A-3.1-0","Ls-3.1-0"] <- chg["A-3.1-0","Ls-3.1-0"]  + params['default_artatb']*params['tau_d0_h3']
  chg["A-3.1-0","LIPT-3.2-0"] <- chg["A-3.1-0","LIPT-3.2-0"]  + (1-params['default_artatb'])*params['kappa']*params['tau_d0_h3']
  chg["A-3.1-0","Ls-3.2-0"] <- chg["A-3.1-0","Ls-3.2-0"]  + (1-params['default_artatb'])*(1-params['kappa'])*params['tau_d0_h3']
  
  chg["A-2.1-0","Ls-2.1-0"] <- chg["A-2.1-0","Ls-2.1-0"]  + params['default_artatb']*mean(params[c('tau_d0_h3','tau_d0_h0')])
  chg["A-2.1-0","LIPT-2.2-0"] <- chg["A-2.1-0","LIPT-2.2-0"]  + (1-params['default_artatb'])*params['kappa']*mean(params[c('tau_d0_h3','tau_d0_h0')])
  chg["A-2.1-0","Ls-2.2-0"] <- chg["A-2.1-0","Ls-2.2-0"]  + (1-params['default_artatb'])*(1-params['kappa'])*mean(params[c('tau_d0_h3','tau_d0_h0')])
  
  chg["A-1.1-0","Ls-1.1-0"] <- chg["A-1.1-0","Ls-1.1-0"]  + params['default_artatb']*mean(params[c('tau_d0_h3','tau_d0_h0')])
  chg["A-1.1-0","LIPT-1.2-0"] <- chg["A-1.1-0","LIPT-1.2-0"]  + (1-params['default_artatb'])*params['kappa']*mean(params[c('tau_d0_h3','tau_d0_h0')])
  chg["A-1.1-0","Ls-1.2-0"] <- chg["A-1.1-0","Ls-1.2-0"]  + (1-params['default_artatb'])*(1-params['kappa'])*mean(params[c('tau_d0_h3','tau_d0_h0')])
  
  
  ## drug resistant
  chg[cbind(intersect(intersect(lu$a.i,lu$d1.i),lu$h0.i),intersect(intersect(lu$ls.i,lu$d1.i),lu$h0.i))] <-
    chg[cbind(intersect(intersect(lu$a.i,lu$d1.i),lu$h0.i),intersect(intersect(lu$ls.i,lu$d1.i),lu$h0.i))]+params['tau_d1']
  
  # if were already on ART, just treat the TB
  chg[cbind(intersect(intersect(lu$a.i,lu$d1.i),lu$allart.i),intersect(intersect(lu$ls.i,lu$d1.i),lu$allart.i))] <-
    chg[cbind(intersect(intersect(lu$a.i,lu$d1.i),lu$allart.i),intersect(intersect(lu$ls.i,lu$d1.i),lu$allart.i))]+params['tau_d1']
  
  # if not already on ART, treat the TB, start ART with some probability, and if starting ART, start IPT with probability kappa (e.g. 0.85)
  chg["A-3.1-1","Ls-3.1-1"] <- chg["A-3.1-1","Ls-3.1-1"]  + params['default_artatb']*params['tau_d1']
  chg["A-3.1-1","LIPT-3.2-1"] <- chg["A-3.1-1","LIPT-3.2-1"]  + (1-params['default_artatb'])*params['kappa']*params['tau_d1']
  chg["A-3.1-1","Ls-3.2-1"] <- chg["A-3.1-1","Ls-3.2-1"]  + (1-params['default_artatb'])*(1-params['kappa'])*params['tau_d1']
  
  chg["A-2.1-1","Ls-2.1-1"] <- chg["A-2.1-1","Ls-2.1-1"]  + params['default_artatb']*params['tau_d1']
  chg["A-2.1-1","LIPT-2.2-1"] <- chg["A-2.1-1","LIPT-2.2-1"]  + (1-params['default_artatb'])*params['kappa']*params['tau_d1']
  chg["A-2.1-1","Ls-2.2-1"] <- chg["A-2.1-1","Ls-2.2-1"]  + (1-params['default_artatb'])*(1-params['kappa'])*params['tau_d1']
  
  chg["A-1.1-1","Ls-1.1-1"] <- chg["A-1.1-1","Ls-1.1-1"]  + params['default_artatb']*params['tau_d1']
  chg["A-1.1-1","LIPT-1.2-1"] <- chg["A-1.1-1","LIPT-1.2-1"]  + (1-params['default_artatb'])*params['kappa']*params['tau_d1']
  chg["A-1.1-1","Ls-1.2-1"] <- chg["A-1.1-1","Ls-1.2-1"]  + (1-params['default_artatb'])*(1-params['kappa'])*params['tau_d1']
  
    
  ## ---------------------------------- ##
  ## on to IPT for those already on ART ##
  ## ---------------------------------- ##
  
  #recovery
  ## susceptibles (no infection) put on IPT
  chg[cbind(c("S-1.2","S-2.2","S-3.2"),c("SIPT-1.2","SIPT-2.2","SIPT-3.2"))] <-
    chg[cbind(c("S-1.2","S-2.2","S-3.2"),c("SIPT-1.2","SIPT-2.2","SIPT-3.2"))] + params['kappa_onart']
  
  ## Lf
  chg[cbind(intersect(lu$lf.i,lu$h1art.i),intersect(lu$liptf.i,lu$h1art.i))] <-
    chg[cbind(intersect(lu$lf.i,lu$h1art.i),intersect(lu$liptf.i,lu$h1art.i))] + params['kappa_onart']
  chg[cbind(intersect(lu$lf.i,lu$h2art.i),intersect(lu$liptf.i,lu$h2art.i))] <-
    chg[cbind(intersect(lu$lf.i,lu$h2art.i),intersect(lu$liptf.i,lu$h2art.i))] + params['kappa_onart']
  chg[cbind(intersect(lu$lf.i,lu$h3art.i),intersect(lu$liptf.i,lu$h3art.i))] <-
    chg[cbind(intersect(lu$lf.i,lu$h3art.i),intersect(lu$liptf.i,lu$h3art.i))] + params['kappa_onart']
  
  ## Ls
  chg[cbind(intersect(lu$ls.i,lu$h1art.i),intersect(lu$lipt.i,lu$h1art.i))] <-
    chg[cbind(intersect(lu$ls.i,lu$h1art.i),intersect(lu$lipt.i,lu$h1art.i))] + params['kappa_onart']
  chg[cbind(intersect(lu$ls.i,lu$h2art.i),intersect(lu$lipt.i,lu$h2art.i))] <-
    chg[cbind(intersect(lu$ls.i,lu$h2art.i),intersect(lu$lipt.i,lu$h2art.i))] + params['kappa_onart']
  chg[cbind(intersect(lu$ls.i,lu$h3art.i),intersect(lu$lipt.i,lu$h3art.i))] <-
    chg[cbind(intersect(lu$ls.i,lu$h3art.i),intersect(lu$lipt.i,lu$h3art.i))] + params['kappa_onart']
  
  ## --------------- ##
  ## IPT to POST-IPT ##
  ## --------------- ##
  chg[cbind(lu$sipt.i,lu$sipt2.i)] <-
    chg[cbind(lu$sipt.i,lu$sipt2.i)] + params['gamma']
  chg[cbind(lu$lipt.i,lu$lipt2.i)] <-
    chg[cbind(lu$lipt.i,lu$lipt2.i)] + params['gamma']
  chg[cbind(lu$liptf.i,lu$lipt2f.i)] <-
    chg[cbind(lu$liptf.i,lu$lipt2f.i)] + params['gamma']
  
  ## ------------------------ ##
  ## Progression while on IPT ##
  ## ------------------------ ##
  chg["LIPT-1.2-0","A-1.2-0"]<-chg["LIPT-1.2-0","A-1.2-0"] + params['theta']*params['h1theta']*params['iepsilon_h']*params['arttheta']
  chg["LIPT-2.2-0","A-2.2-0"]<-chg["LIPT-2.2-0","A-2.2-0"] + params['theta']*params['h2theta']*params['iepsilon_h']*params['arttheta']
  chg["LIPT-3.2-0","A-3.2-0"]<-chg["LIPT-3.2-0","A-3.2-0"] + params['theta']*params['iepsilon_h']*params['arttheta']

  chg["LIPTf-1.2-0","A-1.2-0"]<-chg["LIPTf-1.2-0","A-1.2-0"] + params['theta']*params['h1theta']*params['ipi_h']*params['arttheta']
  chg["LIPTf-2.2-0","A-2.2-0"]<-chg["LIPTf-2.2-0","A-2.2-0"] + params['theta']*params['h2theta']*params['ipi_h']*params['arttheta']
  chg["LIPTf-3.2-0","A-3.2-0"]<-chg["LIPTf-3.2-0","A-3.2-0"] + params['theta']*params['ipi_h']*params['arttheta']
  
  ## drug resistant receive no benifit from IPT
  chg["LIPT-1.2-1","A-1.2-1"]<-chg["LIPT-1.2-1","A-1.2-1"] + params['h1theta']*params['iepsilon_h']*params['arttheta']
  chg["LIPT-2.2-1","A-2.2-1"]<-chg["LIPT-2.2-1","A-2.2-1"] + params['h2theta']*params['iepsilon_h']*params['arttheta']
  chg["LIPT-3.2-1","A-3.2-1"]<-chg["LIPT-3.2-1","A-3.2-1"] + params['iepsilon_h']*params['arttheta']
  
  chg["LIPTf-1.2-1","A-1.2-1"]<-chg["LIPTf-1.2-1","A-1.2-1"] + params['h1theta']*params['ipi_h']*params['arttheta']
  chg["LIPTf-2.2-1","A-2.2-1"]<-chg["LIPTf-2.2-1","A-2.2-1"] + params['h2theta']*params['ipi_h']*params['arttheta']
  chg["LIPTf-3.2-1","A-3.2-1"]<-chg["LIPTf-3.2-1","A-3.2-1"] + params['ipi_h']*params['arttheta']
  
  ## -------------------- ##
  ## Progression Post IPT ##
  ## -------------------- ##
  chg["LIPT2-1.2-0","A-1.2-0"]<-chg["LIPT2-1.2-0","A-1.2-0"] + params['theta']*params['h1theta']*params['iepsilon_h']*params['arttheta']
  chg["LIPT2-2.2-0","A-2.2-0"]<-chg["LIPT2-2.2-0","A-2.2-0"] + params['theta']*params['h2theta']*params['iepsilon_h']*params['arttheta']
  chg["LIPT2-3.2-0","A-3.2-0"]<-chg["LIPT2-3.2-0","A-3.2-0"] + params['theta']*params['iepsilon_h']*params['arttheta']
  
  chg["LIPT2f-1.2-0","A-1.2-0"]<-chg["LIPT2f-1.2-0","A-1.2-0"] + params['theta']*params['h1theta']*params['ipi_h']*params['arttheta']
  chg["LIPT2f-2.2-0","A-2.2-0"]<-chg["LIPT2f-2.2-0","A-2.2-0"] + params['theta']*params['h2theta']*params['ipi_h']*params['arttheta']
  chg["LIPT2f-3.2-0","A-3.2-0"]<-chg["LIPT2f-3.2-0","A-3.2-0"] + params['theta']*params['ipi_h']*params['arttheta']
  
  chg["LIPT2-1.2-1","A-1.2-1"]<-chg["LIPT2-1.2-1","A-1.2-1"] + params['h1theta']*params['iepsilon_h']*params['arttheta']
  chg["LIPT2-2.2-1","A-2.2-1"]<-chg["LIPT2-2.2-1","A-2.2-1"] + params['h2theta']*params['iepsilon_h']*params['arttheta']
  chg["LIPT2-3.2-1","A-3.2-1"]<-chg["LIPT2-3.2-1","A-3.2-1"] + params['iepsilon_h']*params['arttheta']
  
  chg["LIPT2f-1.2-1","A-1.2-1"]<-chg["LIPT2f-1.2-1","A-1.2-1"] + params['h1theta']*params['ipi_h']*params['arttheta']
  chg["LIPT2f-2.2-1","A-2.2-1"]<-chg["LIPT2f-2.2-1","A-2.2-1"] + params['h2theta']*params['ipi_h']*params['arttheta']
  chg["LIPT2f-3.2-1","A-3.2-1"]<-chg["LIPT2f-3.2-1","A-3.2-1"] + params['ipi_h']*params['arttheta']
  
  
  #####################
  ## HIV Transitions ##
  #####################
  
  ## NEw HIV infections in dynamic updates (dxdt)
  
  ## ---------------------- ##
  ## Progression to CD4<500 ##
  ## ---------------------- ##
  
  chg[cbind(lu$h1nonart.i,lu$h2nonart.i)] <-chg[cbind(lu$h1nonart.i,lu$h2nonart.i)] + params['eta1']
  
  ## ---------------------- ##
  ## Progression to CD4<200 ##
  ## ---------------------- ##
  chg[cbind(lu$h2nonart.i,lu$h3nonart.i)] <- chg[cbind(lu$h2nonart.i,lu$h3nonart.i)] + params['eta2']
  

  ## ART initiation +- IPT now in dynamic updates
  

  ## ---------------------- ##
  ## CD4+ Rebound after ART ##
  ## ---------------------- ##
  
  chg[cbind(lu$h2art.i,lu$h1art.i)] <-
    chg[cbind(lu$h2art.i,lu$h1art.i)] + params['eta3']

  chg[cbind(lu$h3art.i,lu$h2art.i)] <-
    chg[cbind(lu$h3art.i,lu$h2art.i)] + params['eta4']
  
  return(chg)
}

##' delta function for model
##' @param t time
##' @param state named array with values for compartments
##' @param params named list of parameters
##' @param chg named array with flows between states
##' this can be made with make.chg.mat()
##' @param lu - lookup list of indices
##' @return updated state of system
##' @author andrew azman
tb.dx.dt <- function(t,state,params,chg,lu, projection=F, baseart=0.5) {
  
  ## a little set up
  dxdt <- numeric(nrow(chg))
  names(dxdt) <- names(state)
  names(state) <- rownames(chg)
  
  N <- sum(state)
  
  params['iepsilon_0'] <- params['epsilon_0']
  params['iepsilon_h'] <- params['epsilon_h']
  params['ipi_0'] <- params['pi_0']
  params['ipi_h'] <- params['pi_h']
  
  # time-varying elements for after calibration only
  # in future projections, double art initiation rate above 200 and add same for above 500
  # and consider linking hiv infection rate (eta0) to untreated-hiv prevalence, so that it starts low but falls over time
  
  if (projection)
  {
    params['chi2_h1'] <- params['chi2_h2'] <- params['chi2_h2']*(1+ 0.5*t)
    params['chi2_h3'] <- params['chi2_h3']*(1+ 0.1*t)
    params['eta0'] <- params['eta0']*sum(state[c(lu$h1art.i,lu$h2art.i,lu$h3art.i)])/sum(state[-c(lu$h0)])/(1-baseart)
  }
  
  ## the chg matrix doesn't include transitions with the force of infection
  ## (i.e., infection and re-infection) nor births (since we are keeping a
  ## constant population size and cannceling out deaths).
  
  ## --------------------- ##
  ## some name stuff first ##
  ## --------------------- ##
  
  state.names <- rownames(chg)
  
  ## ------------------ ##
  ## New HIV Infections ##
  ## ------------------ ##
  chg[cbind(lu$h0.i,intersect(lu$h1.i,lu$nonart.i))] <-
    chg[cbind(lu$h0.i,intersect(lu$h1.i,lu$nonart.i))] + params['eta0']
  
  
  ## ------------------ ##
  ## force of infection ##
  ## ------------------ ##
  
  inf.weight.d0 <-
    state["A-0-0"] +
    params['phi']*sum(state[intersect(lu$a.i,lu$d0.i)]) -
    params['phi']*state["A-0-0"]
  
  inf.weight.d1 <- inf.weight.d0 * params['prob.MDR']/(1-params['prob.MDR'])
  
  foi <- params['beta'] / N * (inf.weight.d0 + inf.weight.d1)
  
  if(is.integer(t/100)) {cat(paste0("foi ",foi,", prob.MDR ",params['prob.MDR'],"\n"))}
  ## -------------------------------- ##
  ## Dynamic transistions from S to L ##
  ## -------------------------------- ##
  
  ## first drug susceptible infections
  chg[cbind(lu$s.i,intersect(lu$lf.i,lu$d0.i))]<-inf.weight.d0*params['beta']/N
  
  ## MDR infections
  chg[cbind(lu$s.i,intersect(lu$lf.i,lu$d1.i))]<- inf.weight.d1*params['beta']/N
  
  ## --------------- ##
  ## SIPT --> Lf MDR ##
  ## --------------- ##
  chg[cbind(c("SIPT-1.2","SIPT-2.2","SIPT-3.2"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    inf.weight.d1*params['beta']/N
  
  
  ## note: adding this possibility of getting infected while on IPT for subanalyses but this param is 0 otherwise
  chg[cbind(c("SIPT-1.2","SIPT-2.2","SIPT-3.2" ),c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))] <-
    params['inf_on_ipt']*inf.weight.d0*params['beta']/N
  
  ## --------------- ##
  ## LIPT --> Lf MDR ##
  ## --------------- ##
  
  chg[cbind(c("LIPT-1.2-0","LIPT-2.2-0","LIPT-3.2-0"),
            c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    inf.weight.d1*params['psi']*params['beta']/N
  
  chg[cbind(c("LIPTf-1.2-0","LIPTf-2.2-0","LIPTf-3.2-0"),
            c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    inf.weight.d1*params['psi']*params['beta']/N
  
  chg[cbind(c("LIPT-1.2-1","LIPT-2.2-1","LIPT-3.2-1"),
            c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    inf.weight.d1*params['psi']*params['beta']/N
  
  chg[cbind(c("LIPTf-1.2-1","LIPTf-2.2-1","LIPTf-3.2-1"),
            c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    inf.weight.d1*params['psi']*params['beta']/N
  
  ## non-MDR transitions for some sensitivty analysis but flag should make this
  ## zero in main analyses
  chg[cbind(c("LIPT-1.2-0","LIPT-2.2-0","LIPT-3.2-0"),
            c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    params['inf_on_ipt']*inf.weight.d0*params['psi']*params['beta']/N
  
  chg[cbind(c("LIPTf-1.2-0","LIPTf-2.2-0","LIPTf-3.2-0"),
            c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    params['inf_on_ipt']*inf.weight.d0*params['psi']*params['beta']/N
  
  chg[cbind(c("LIPT-1.2-1","LIPT-2.2-1","LIPT-3.2-1"),
            c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    params['inf_on_ipt']*inf.weight.d0*params['psi']*params['beta']/N
  
  chg[cbind(c("LIPTf-1.2-1","LIPTf-2.2-1","LIPTf-3.2-1"),
            c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    params['inf_on_ipt']*inf.weight.d0*params['psi']*params['beta']/N
  
  ## ------------ ##
  ## SIPT2 --> Lf ##
  ## ------------ ##
  ## new infections post IPT
  chg[cbind(c("SIPT2-1.2","SIPT2-2.2","SIPT2-3.2"),c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    inf.weight.d0*params['beta']/N
  ## MDR infections
  chg[cbind(c("SIPT2-1.2","SIPT2-2.2","SIPT2-3.2"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    inf.weight.d1*params['beta']/N
  
  ## ------------ ##
  ## LIPT2 --> Lf ##
  ## ------------ ##
  ## new infections post-ipt
  chg[cbind(c("LIPT2-1.2-0","LIPT2-2.2-0","LIPT2-3.2-0"),c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    inf.weight.d0*params["psi"]*params['beta']/N
  
  chg[cbind(c("LIPT2f-1.2-0","LIPT2f-2.2-0","LIPT2f-3.2-0"),c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    inf.weight.d0*params["psi"]*params['beta']/N
  
  chg[cbind(c("LIPT2-1.2-0","LIPT2-2.2-0","LIPT2-3.2-0"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    inf.weight.d1*params["psi"]*params['beta']/N
  
  chg[cbind(c("LIPT2f-1.2-0","LIPT2f-2.2-0","LIPT2f-3.2-0"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    inf.weight.d1*params["psi"]*params['beta']/N
  
  ## MDR latents
  chg[cbind(c("LIPT2-1.2-1","LIPT2-2.2-1","LIPT2-3.2-1"),c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    inf.weight.d0*params['beta']/N*params["psi"]
  chg[cbind(c("LIPT2f-1.2-1","LIPT2f-2.2-1","LIPT2f-3.2-1"),c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    inf.weight.d0*params['beta']/N*params["psi"]
  
  chg[cbind(c("LIPT2-1.2-1","LIPT2-2.2-1","LIPT2-3.2-1"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    inf.weight.d1*params['beta']/N*params["psi"]
  chg[cbind(c("LIPT2f-1.2-1","LIPT2f-2.2-1","LIPT2-3.2-1"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    inf.weight.d1*params['beta']/N*params["psi"]
  
  ## ------------- ##``
  ## re-infections ##
  ## Ls --> Lf     ##
  ## ------------- ##
  
  ## both MDR and non-MDR can be infected with DS and for simplicity
  ## i am splitting this up here
  chg[cbind(intersect(lu$ls.i,lu$d0.i),intersect(lu$lf.i,lu$d0.i))]<-
    inf.weight.d0*params['beta']/N*params["psi"]
  chg[cbind(intersect(lu$ls.i,lu$d1.i),intersect(lu$lf.i,lu$d0.i))]<-
    inf.weight.d0*params['beta']/N*params["psi"]
  
  ## MDR infections
  chg[cbind(intersect(lu$ls.i,lu$d0.i),intersect(lu$lf.i,lu$d1.i))]<-
    inf.weight.d1*params['beta']/N*params["psi"]
  chg[cbind(intersect(lu$ls.i,lu$d1.i),intersect(lu$lf.i,lu$d1.i))]<-
    inf.weight.d1*params['beta']/N*params["psi"]
  
  ## ------------------------------- ##
  ## re-infections                   ##
  ## Lf --> Lf with other resistance ##
  ## ------------------------------- ##
  
  ## DS to MDR 
  chg[cbind(intersect(lu$lf.i,lu$d0.i),intersect(lu$lf.i,lu$d1.i))]<-
    inf.weight.d1*params['beta']/N*params["psi"]
  
  ## MDR to DS
  chg[cbind(intersect(lu$lf.i,lu$d1.i),intersect(lu$lf.i,lu$d0.i))]<-
    inf.weight.d0*params['beta']/N*params["psi"]
  
  
  ##-------------------------------- ##
  ## Progression minus reinfections: 
  ## adjustment of Lf->Ls to account for Lf->Lf ##
  ##-------------------------------- ##
  
  chg[cbind(intersect(lu$lf.i,lu$d0.i),intersect(lu$ls.i,lu$d0.i))]<-
    chg[cbind(intersect(lu$lf.i,lu$d0.i),intersect(lu$ls.i,lu$d0.i))] + 
    params['sigma'] * exp(-(inf.weight.d0+inf.weight.d1)*params['beta']/N*params["psi"]/params['sigma'])
  
  chg[cbind(intersect(lu$lf.i,lu$d1.i),intersect(lu$ls.i,lu$d1.i))]<-
    chg[cbind(intersect(lu$lf.i,lu$d1.i),intersect(lu$ls.i,lu$d1.i))] + 
    params['sigma'] * exp(-(inf.weight.d0+inf.weight.d1)*params['beta']/N*params["psi"]/params['sigma'])
  
  
  
  ## ----------------------- ##
  ## ART Initiation plus IPT ## 
  ##  (plus TB screen & ?Rx) ##
  ## ----------------------- ##
  
  ## S
  ## 1-\kappa proprtion of those going on ART don't go on IPT too
  chg["S-1.1","S-1.2"]<-chg["S-1.1","S-1.2"]+(1-params['kappa'])*params['chi2_h1']
  chg["S-2.1","S-2.2"]<-chg["S-2.1","S-2.2"]+(1-params['kappa'])*params['chi2_h2']
  chg["S-3.1","S-3.2"]<-chg["S-3.1","S-3.2"]+(1-params['kappa'])*params['chi2_h3']
  chg["S-1.1","SIPT-1.2"]<- chg["S-1.1","SIPT-1.2"]+params['kappa']*params['chi2_h1']
  chg["S-2.1","SIPT-2.2"]<- chg["S-2.1","SIPT-2.2"]+params['kappa']*params['chi2_h2']
  chg["S-3.1","SIPT-3.2"]<- chg["S-3.1","SIPT-3.2"]+params['kappa']*params['chi2_h3']
  
  ## Lf
  ## latents that don't get IPT stay as latent but now on ART
  chg[cbind(c("Lf-1.1-0","Lf-1.1-1","Ls-1.1-0","Ls-1.1-1"),
            c("Lf-1.2-0","Lf-1.2-1","Ls-1.2-0","Ls-1.2-1"))]<-
    chg[cbind(c("Lf-1.1-0","Lf-1.1-1","Ls-1.1-0","Ls-1.1-1"),
              c("Lf-1.2-0","Lf-1.2-1","Ls-1.2-0","Ls-1.2-1"))] +
    (1-params['kappa'])*params['chi2_h1']
  chg[cbind(c("Lf-2.1-0","Lf-2.1-1","Ls-2.1-0","Ls-2.1-1"),
            c("Lf-2.2-0","Lf-2.2-1","Ls-2.2-0","Ls-2.2-1"))]<-
    chg[cbind(c("Lf-2.1-0","Lf-2.1-1","Ls-2.1-0","Ls-2.1-1"),
              c("Lf-2.2-0","Lf-2.2-1","Ls-2.2-0","Ls-2.2-1"))] +
    (1-params['kappa'])*params['chi2_h2']
  chg[cbind(c("Lf-3.1-0","Lf-3.1-1","Ls-3.1-0","Ls-3.1-1"),
            c("Lf-3.2-0","Lf-3.2-1","Ls-3.2-0","Ls-3.2-1"))]<-
    chg[cbind(c("Lf-3.1-0","Lf-3.1-1","Ls-3.1-0","Ls-3.1-1"),
              c("Lf-3.2-0","Lf-3.2-1","Ls-3.2-0","Ls-3.2-1"))] +
    (1-params['kappa'])*params['chi2_h3']
  
  ## latents that get IPT
  chg[cbind(c("Lf-1.1-0","Lf-1.1-1","Ls-1.1-0","Ls-1.1-1"),
            c("LIPTf-1.2-0","LIPTf-1.2-1","LIPT-1.2-0","LIPT-1.2-1"))] <-
    chg[cbind(c("Lf-1.1-0","Lf-1.1-1","Ls-1.1-0","Ls-1.1-1"),
              c("LIPTf-1.2-0","LIPTf-1.2-1","LIPT-1.2-0","LIPT-1.2-1"))] +
    params['kappa']*params['chi2_h1']
  chg[cbind(c("Lf-2.1-0","Lf-2.1-1","Ls-2.1-0","Ls-2.1-1"),
            c("LIPTf-2.2-0","LIPTf-2.2-1","LIPT-2.2-0","LIPT-2.2-1"))] <-
    chg[cbind(c("Lf-2.1-0","Lf-2.1-1","Ls-2.1-0","Ls-2.1-1"),
              c("LIPTf-2.2-0","LIPTf-2.2-1","LIPT-2.2-0","LIPT-2.2-1"))] +
    params['kappa']*params['chi2_h2']
  chg[cbind(c("Lf-3.1-0","Lf-3.1-1","Ls-3.1-0","Ls-3.1-1"),
            c("LIPTf-3.2-0","LIPTf-3.2-1","LIPT-3.2-0","LIPT-3.2-1"))] <-
    chg[cbind(c("Lf-3.1-0","Lf-3.1-1","Ls-3.1-0","Ls-3.1-1"),
              c("LIPTf-3.2-0","LIPTf-3.2-1","LIPT-3.2-0","LIPT-3.2-1"))] +
    params['kappa']*params['chi2_h3']
  
  ## Actives
  # A -> A + ART
  ## A fraction of actives  are not screened (either in IPT or ART initiation-based screening)
  ## or are screened but default or are not detected
  chg[cbind(c("A-1.1-0","A-1.1-1"),c("A-1.2-0","A-1.2-1"))]<-
    chg[cbind(c("A-1.1-0","A-1.1-1"),c("A-1.2-0","A-1.2-1"))]+
    params['chi2_h1']* (
      (1-params['kappa'])*((1-params['art_screen']) + params['art_screen']*((1-params['sens_xpert']) + params['sens_xpert']*params['default_acf']) ) + 
        params['kappa']*((1-params['art_screen'])*((1-params['ipt_screen']) + params['ipt_screen']*((1-params['sens_xpert']) + params['sens_xpert']*params['default_acf'])) + 
                           params['art_screen']*((1-params['sens_xpert'])+params['sens_xpert']*params['default_acf'])) )
  
    chg[cbind(c("A-2.1-0","A-2.1-1"),c("A-2.2-0","A-2.2-1"))]<-
    chg[cbind(c("A-2.1-0","A-2.1-1"),c("A-2.2-0","A-2.2-1"))]+
    params['chi2_h2']* (
      (1-params['kappa'])*((1-params['art_screen']) + params['art_screen']*((1-params['sens_xpert']) + params['sens_xpert']*params['default_acf']) ) + 
        params['kappa']*((1-params['art_screen'])*((1-params['ipt_screen']) + params['ipt_screen']*((1-params['sens_xpert']) + params['sens_xpert']*params['default_acf'])) + 
                           params['art_screen']*((1-params['sens_xpert'])+params['sens_xpert']*params['default_acf'])) )
    
    chg[cbind(c("A-3.1-0","A-3.1-1"),c("A-3.2-0","A-3.2-1"))]<-
    chg[cbind(c("A-3.1-0","A-3.1-1"),c("A-3.2-0","A-3.2-1"))]+
    params['chi2_h3']* (
      (1-params['kappa'])*((1-params['art_screen']) + params['art_screen']*((1-params['sens_xpert']) + params['sens_xpert']*params['default_acf']) ) + 
        params['kappa']*((1-params['art_screen'])*((1-params['ipt_screen']) + params['ipt_screen']*((1-params['sens_xpert']) + params['sens_xpert']*params['default_acf'])) + 
                           params['art_screen']*((1-params['sens_xpert'])+params['sens_xpert']*params['default_acf'])) )
      

  # ART initiation, with TB screening (as part of the ART itself or through linked IPT screening), can lead to TB treatment
  # If IPT screening happens:
  ## The fraction who are screened for IPT and found to have active TB, will short circuit back to LIPT, 
  ## (since we didn't explcilty model time on treatment)
    ## They also still start ART (.1 -> .2).
  ## now adding IPT after active TB treatment for these ART initiators too, so that they become not Ls but L-IPT
  
  chg[cbind(c("A-1.1-0","A-1.1-1"),c("LIPT-1.2-0","LIPT-1.2-1"))]<-
    chg[cbind(c("A-1.1-0","A-1.1-1"),c("LIPT-1.2-0","LIPT-1.2-1"))]+
    params['chi2_h1'] * (
      params['kappa']*(params['art_screen']+ (1-params['art_screen'])*params['ipt_screen'])*params['sens_xpert']*(1-params['default_acf']))
  
  chg[cbind(c("A-1.1-0","A-1.1-1"),c("Ls-1.2-0","Ls-1.2-1"))]<-
    chg[cbind(c("A-1.1-0","A-1.1-1"),c("Ls-1.2-0","Ls-1.2-1"))]+
    params['chi2_h1'] * (
      (1-params['kappa'])*(params['art_screen'])*params['sens_xpert']*(1-params['default_acf']))
  
  chg[cbind(c("A-2.1-0","A-2.1-1"),c("LIPT-2.2-0","LIPT-2.2-1"))]<-
    chg[cbind(c("A-2.1-0","A-2.1-1"),c("LIPT-2.2-0","LIPT-2.2-1"))]+
    params['chi2_h2'] * (
      params['kappa']*(params['art_screen']+ (1-params['art_screen'])*params['ipt_screen'])*params['sens_xpert']*(1-params['default_acf']))
  
  chg[cbind(c("A-2.1-0","A-2.1-1"),c("Ls-2.2-0","Ls-2.2-1"))]<-
    chg[cbind(c("A-2.1-0","A-2.1-1"),c("Ls-2.2-0","Ls-2.2-1"))]+
    params['chi2_h2'] * (
      (1-params['kappa'])*(params['art_screen'])*params['sens_xpert']*(1-params['default_acf']))
  
  chg[cbind(c("A-3.1-0","A-3.1-1"),c("LIPT-3.2-0","LIPT-3.2-1"))]<-
    chg[cbind(c("A-3.1-0","A-3.1-1"),c("LIPT-3.2-0","LIPT-3.2-1"))]+
    params['chi2_h3'] * (
      params['kappa']*(params['art_screen']+ (1-params['art_screen'])*params['ipt_screen'])*params['sens_xpert']*(1-params['default_acf']))
  
  chg[cbind(c("A-3.1-0","A-3.1-1"),c("Ls-3.2-0","Ls-3.2-1"))]<-
    chg[cbind(c("A-3.1-0","A-3.1-1"),c("Ls-3.2-0","Ls-3.2-1"))]+
    params['chi2_h3'] * (
      (1-params['kappa'])*(params['art_screen'])*params['sens_xpert']*(1-params['default_acf']))
  
  


  ######################    
  ## in - outs - deaths
  ## remove deaths from chg (chg matrix after this will not have any non-zero diag. elements)
  tmp.deaths <- diag(chg)
  diag(chg) <- 0
  dxdt <- (t(chg) %*% state) - rowSums(chg * state) - tmp.deaths*state
  
  ## lets update the susceptible states to add in births into susceptibles
  dxdt['S-0',1] <-
    dxdt['S-0',1] + (1-params['rho'])*tmp.deaths%*%state
  ## and latents since this is an adult model
  
  ## don't really like the way we are apportioning MDR here
  ## since we are assuming the proportion is reflective of today
  ## when it is really about the past. Should be ok for this time frame
  ## though. ETA: and we're at equilibrium anyway so doesn't matter.
  dxdt['Ls-0-0',1] <-
    dxdt['Ls-0-0',1] + params['rho']*(1-params['prob.MDR'])*tmp.deaths %*% state
  ##  dxdt['Ls-0-0',1] + params['rho']*(1-params['prob.MDR'])*tmp.deaths %*% state
  
  dxdt['Ls-0-1',1] <-
    dxdt['Ls-0-1',1] + params['rho']*params['prob.MDR']*tmp.deaths %*% state
  ##  dxdt['Ls-0-1',1] + params['rho']*params['prob.MDR']*tmp.deaths %*% state
  
  return(list(dxdt))
}

##' delta function for model
##' attempting to replicate trial
##' Keeps fixed FOI for trial duration, and models deaths but not births
##' @param t time
##' @param state named array with values for compartments
##' @param params named list of parameters
##' @param chg named array with flows between states
##' this can be made with make.chg.mat()
##' @param lu - lookup list of indices
##' @param foi total force of infection (infectious weight sun * beta/N)
##' @return updated state of system
##' @author andrew azman
tb.dx.dt.trial <- function(t,
                           state,
                           params,
                           chg,
                           lu,
                           foi) {
  
  ## a little set up
  dxdt <- numeric(nrow(chg))
  names(dxdt) <- names(state)
  names(state) <- rownames(chg)
  
  N <- sum(state)
  
  params['iepsilon_0'] <- params['epsilon_0']
  params['iepsilon_h'] <- params['epsilon_h']
  params['ipi_0'] <- params['pi_0']
  params['ipi_h'] <- params['pi_h']
  
  ## the chg matrix doesn't include transitions with the force of infection
  ## (i.e., infection and re-infection) nor births (since we are keeping a
  ## constant population size and canceling out deaths).
  
  ## --------------------- ##
  ## some name stuff first ##
  ## --------------------- ##
  state.names <- rownames(chg)
  
  ## -------------------------------- ##
  ## Dynamic transistions from S to L ##
  ## -------------------------------- ##
  
  ## first drug susectople infections
  chg[cbind(lu$s.i,intersect(lu$lf.i,lu$d0.i))] <- foi*(1-params['prob.MDR'])
  
  ## MDR infections
  chg[cbind(lu$s.i,intersect(lu$lf.i,lu$d1.i))]<- foi*params['prob.MDR']
  
  # Infection on IPT, only possible for MDR 
  ## ----------- ##
  ## SIPT --> Lf ## 
  ## ----------- ##
  chg[cbind(c("SIPT-1.2","SIPT-2.2","SIPT-3.2"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    foi*(params['prob.MDR'])
  
  ## ----------- ##
  ## LIPT --> Lf ##
  ## ----------- ##
  chg[cbind(c("LIPT-1.2-0","LIPT-2.2-0","LIPT-3.2-0"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    foi*params['prob.MDR']*params["psi"]
  chg[cbind(c("LIPT-1.2-1","LIPT-2.2-1","LIPT-2.2-1"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    foi*params['prob.MDR']*params['psi']
  
  chg[cbind(c("LIPTf-1.2-0","LIPTf-2.2-0","LIPTf-3.2-0"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    foi*params['prob.MDR']*params["psi"]
  chg[cbind(c("LIPTf-1.2-1","LIPTf-2.2-1","LIPTf-3.2-1"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    foi*params['prob.MDR']*params['psi']
  
  # infection post IPT
  ## ------------ ##
  ## SIPT2 --> Lf ##
  ## ------------ ##
  chg[cbind(c("SIPT2-1.2","SIPT2-2.2","SIPT2-3.2"),c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    foi*(1-params['prob.MDR'])
  
  ## MDR infections
  chg[cbind(c("SIPT2-1.2","SIPT2-2.2","SIPT2-3.2"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    foi*(params['prob.MDR'])
  
  ## ------------ ##
  ## LIPT2 --> Lf ##
  ## ------------ ##
  chg[cbind(c("LIPT2-1.2-0","LIPT2-2.2-0","LIPT2-3.2-0"),c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    foi*(1-params['prob.MDR'])*params['psi']
  chg[cbind(c("LIPT2-1.2-0","LIPT2-2.2-0","LIPT2-3.2-0"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    foi*params['prob.MDR']*params['psi']
  
  chg[cbind(c("LIPT2f-1.2-0","LIPT2f-2.2-0","LIPT2f-3.2-0"),c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    foi*(1-params['prob.MDR'])*params['psi']
  chg[cbind(c("LIPT2f-1.2-0","LIPT2f-2.2-0","LIPT2f-3.2-0"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    foi*params['prob.MDR']*params['psi']
  
  ## infections of prior MDR
  chg[cbind(c("LIPT2-1.2-1","LIPT2-2.2-1","LIPT2-3.2-1"),c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    foi*(1-params['prob.MDR'])*params['psi']
  chg[cbind(c("LIPT2-1.2-1","LIPT2-2.2-1","LIPT2-3.2-1"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    foi*params['prob.MDR']*params['psi']
  
  chg[cbind(c("LIPT2f-1.2-1","LIPT2f-2.2-1","LIPT2f-3.2-1"),c("Lf-1.2-0","Lf-2.2-0","Lf-3.2-0"))]<-
    foi*(1-params['prob.MDR'])*params['psi']
  chg[cbind(c("LIPT2f-1.2-1","LIPT2f-2.2-1","LIPT2f-3.2-1"),c("Lf-1.2-1","Lf-2.2-1","Lf-3.2-1"))]<-
    foi*params['prob.MDR']*params['psi']
  
  ## -------------- ##
  ##  re-infections ##
  ## Ls --> Lf      ##
  ## -------------- ##
  
  ## both MDR and non-MDR can be infected with DS 
  chg[cbind(intersect(lu$ls.i,lu$d0.i),intersect(lu$lf.i,lu$d0.i))]<-
    foi*(1-params['prob.MDR'])*params["psi"]
  chg[cbind(intersect(lu$ls.i,lu$d1.i),intersect(lu$lf.i,lu$d0.i))]<-
    foi*(1-params['prob.MDR'])*params["psi"]
  
  ## MDR infections
  chg[cbind(intersect(lu$ls.i,lu$d0.i),intersect(lu$lf.i,lu$d1.i))]<-
    foi*params['prob.MDR']*params["psi"]
  chg[cbind(intersect(lu$ls.i,lu$d1.i),intersect(lu$lf.i,lu$d1.i))]<-
    foi*params['prob.MDR']*params["psi"]
  
  
  ## ------------------------------- ##
  ## re-infections                   ##
  ## Lf --> Lf with other resistance ##
  ## ------------------------------- ##
  
  ## MDR to DS
  chg[cbind(intersect(lu$lf.i,lu$d1.i),intersect(lu$lf.i,lu$d0.i))]<-
    foi*(1-params['prob.MDR'])*params["psi"]
  
  ## DS to MDR 
  chg[cbind(intersect(lu$lf.i,lu$d0.i),intersect(lu$lf.i,lu$d1.i))]<-
    foi*(params['prob.MDR'])*params["psi"]
  
  ##-------------------------------- ##
  ## Progression minus reinfections: 
  ## adjustment of Lf->Ls to account for Lf->Lf of same resistance ##
  ##-------------------------------- ##
  
  chg[cbind(intersect(lu$lf.i,lu$d0.i),intersect(lu$ls.i,lu$d0.i))]<-
    chg[cbind(intersect(lu$lf.i,lu$d0.i),intersect(lu$ls.i,lu$d0.i))] + 
    params['sigma'] * exp(-foi*params["psi"]/params['sigma'])
  
  chg[cbind(intersect(lu$lf.i,lu$d1.i),intersect(lu$ls.i,lu$d1.i))]<-
    chg[cbind(intersect(lu$lf.i,lu$d1.i),intersect(lu$ls.i,lu$d1.i))] + 
    params['sigma'] * exp(-foi*params["psi"]/params['sigma'])
  
  
  
  
  ## in - outs - deaths
  ## remove deaths from chg (WARNING: chg matrix after this will not have
  ## and non-zero diag. elements)
  tmp.deaths <- diag(chg)
  diag(chg) <- 0
  dxdt <- (t(chg) %*% state) - rowSums(chg * state) - tmp.deaths*state
  
  return(list(dxdt))
}



## create matrix of LHS draws
get.LHS.draws <- function(n.draws=1000, augment="", reuse=""){
  require(lhs)
  bounds.raw <- read.csv("parameters/parameters-working-fixmdr.csv", header=TRUE)[,c("param","value","min","max")]
  
  rem.params <- c('theta', 'gamma','kappa','kappa_onart','art_screen','ipt_screen','inf_on_ipt')
  bounds <- bounds.raw[-which(bounds.raw[,1] %in% rem.params),]
  
  ## draw LHS from unif(0,1)'s
  if (augment != "")   
  {
    load(paste0("lastLHS.",augment,".Rdata"), verbose = T)
    lastLHS <- unif.LHS
    unif.LHS <- augmentLHS(lastLHS, m=n.draws-nrow(lastLHS))
    save(unif.LHS, file = paste0("lastLHS.",date,".Rdata"))
    unif.LHS <- unif.LHS[(nrow(lastLHS)+1):n.draws,]
    n.draws <- n.draws-nrow(lastLHS)
  } else if (reuse != "")
  {
    load(paste0("lastLHS.",reuse,".Rdata"), verbose = T)
  } else 
  {
    unif.LHS <- randomLHS(n.draws,nrow(bounds))
    save(unif.LHS, file = paste0("lastLHS.",date,".Rdata"))
  }
  
  ## transform them to the appropriate bounds
  ## NOTE: for now not assuming anything structurally about the
  ## realtionship between parameterts (e.g., btw same param for HIV+ and negative)
  trans.param <- sapply(1:nrow(bounds),function(x) qunif(unif.LHS[,x], 
                                                         bounds[x,3],
                                                         bounds[x,4]))
  ## now add the static params
  trans.param <- cbind(trans.param,
                       matrix(rep(bounds.raw[which(bounds.raw[,1] %in% rem.params),2],each=n.draws),nrow=n.draws)
  )
  colnames(trans.param) <-c(as.character(bounds.raw[,1][-which(bounds.raw[,1] %in% rem.params)])
                            ,rem.params)
  return(trans.param)
}

## gets output for calibration
get.calib.outputs <- function(run,params){
  rc <- final.stats(run,params)[1,c('tb.notified','tb.inc','tb.prev','tb.mort',
                                    'prop.MDR','ltbi.prob','hiv.prev',
                                    'ipt.elg.with.tb','cd4.artinit','tb.inc.onart','tb.prev.on.art')]
  return(rc)
}


## log likelihood for 'fitting';
## change intervals to assymetrical
log.lik <- function(data) {#, initiators=T){
  ##incidence, assuming +-10% CI
  ll.not <- log(dnorm((log(data['tb.notified']) - log(1630))/log(1.05)))
  # check: exp(qnorm(c(0.025,0.5,0.975), log(1630), log(1.05)))
  
  ## prevalence assumed to be log-normal with 95% CIs to match the study
  ll.prev <- log(dnorm(log(data['tb.prev']/3200)/log(1.13)))
  # check: exp(qnorm(c(0.025,0.5,0.975), log(3200), log(1.13)))
  
  ## tb.mort
  ll.mort <- log(dnorm(log(data['tb.mort']/500)/log(1.3)))
  # check: exp(qnorm(c(0.025,0.5,0.975), log(500), log(1.3)))
  
  #logit for proportions:
  ll.mdr <- log(dnorm((logit(data['prop.MDR']) - logit(0.04))/((logit(0.055)-logit(0.03))/4)))
  # check: invlogit(qnorm(c(0.025,0.5,0.975), logit(0.04), ((logit(0.055)-logit(0.03))/4)))
                
  # ll.ltbi <- log(dtriangle(data['ltbi.prob'],a=2*0.811-1,b=1,c=0.811))
  
  # 0.8*c(33, 28, 39)
  ms <-  twCoefLogitnormCi(0.2224, 0.312)
  ll.hiv <- log(dnorm((logit(data['hiv.prev'])-ms[1])/ms[2]))
  # check: invlogit(qnorm(c(0.025,0.5,0.975), ms[1], ms[2]))
  
  
  ##mean CD4 in ART initiators, assuming +-20% CI
  ll.cd4 <- log(dnorm((log(data['cd4.artinit']) - log(175))/log(1.1)))
  # check: exp(qnorm(c(0.025,0.5,0.975), log(175), log(1.1)))
  
  
  # library(binom)
  # binom.confint(211+39,2138)
  ms <-  twCoefLogitnormCi(0.104, 0.131)
  # qlogitnorm(p = 0.5, mu = ms[1], sigma=ms[2])
  # if (initiators) ll.ipt.elg.tb <- log(dnorm((logit(data['ipt.elg.with.tb']) - ms[1])/ms[2])) else 
  ll.ipt.elg.tb <- log(dnorm((logit(data['tb.prev.on.art']) - ms[1])/ms[2]))
  
    
  ll.vec <- c(ll.not,ll.prev,ll.mort,ll.mdr,ll.hiv,ll.ipt.elg.tb,ll.cd4)
  ll.vec <- pmax(ll.vec,rep(-1000,length(ll.vec)))
  return(ll.vec)
}

## run model to steady state for each param set and
## return log lik. Flag if not reaching steady state.
one.run <- function(params,total.time=1000){
  init.run <- run.and.print(max.time=100,params=params)
  full.run <- continue.run(init.run,extra.t = total.time,by=1,params=params)
  
  while(max(apply(tail(full.run,10)[,-1],2,diff))>.01 & sum(tail(full.run[55:63]))/sum(tail(full.run[46:63]))<0.1)
  {
    print("adding more time")
    full.run <- continue.run(full.run,
                             extra.t = 500,
                             by=1,
                             params=params)
  }
  # plot.stats(full.run,params)
  
  # add sampling of multiple theta values for log.lik.ipt:
  thetas <- runif(n = 20, 0,1)
  llipts <- lapply(1:20, function(x) {params['theta'] <- thetas[x]; log.lik.ipt(full.run, params)})
  
  return(
    list(
      "ll"=log.lik(get.calib.outputs(full.run,params)), "tail"=tail(full.run,1)[-1],
      "thetas"=thetas, "llthetas"=llipts)
  )
}


##' psuedo log likelihood for the khay-trial ipt effect. Edited by EAK to use the FOI and pMDR for [the start of] each particular simulation.
##' @title
##' @param run
##' @param params
##' @param max.time
##' @return
##' @author asa
log.lik.ipt <- function(run,params,max.time=3){
  
  lu <- make.look.up.list()
  state.names <- make.state.names()
  
  params['iepsilon_0'] <- params['epsilon_0']
  params['iepsilon_h'] <- params['epsilon_h']
  params['ipi_0'] <- params['pi_0']
  params['ipi_h'] <- params['pi_h']
  
  steady.state <- continue.run(
    run,
    params=params,
    extra.t = 1,by=1,print = F) # turned printing on to see foi and probMDR

  DSprev <- sum(steady.state[nrow(steady.state), 1+ intersect(lu$a.i,lu$d0.i)])
    
  ## sample size
  N <- 622
  
  ## figure out distribtion of people in ART states
  inds.of.interest <- setdiff(lu$allart.i,lu$a.i)+1
  cols.of.interest <- steady.state[nrow(steady.state),inds.of.interest]
  prop.in.each <- cols.of.interest/sum(cols.of.interest)
  
  start.state.int <- numeric(ncol(steady.state)-1)
  names(start.state.int) <- make.state.names()
  start.state.int["SIPT-1.2"] <- N*prop.in.each["S-1.2"]
  start.state.int["SIPT-2.2"] <- N*prop.in.each["S-2.2"]
  start.state.int["SIPT-3.2"] <- N*prop.in.each["S-3.2"]
  start.state.int["LIPTf-1.2-0"] <- N*(prop.in.each["Lf-1.2-0"])
  start.state.int["LIPTf-1.2-1"] <- N*(prop.in.each["Lf-1.2-1"])
  start.state.int["LIPTf-2.2-0"] <- N*(prop.in.each["Lf-2.2-0"])
  start.state.int["LIPTf-2.2-1"] <- N*(prop.in.each["Lf-2.2-1"])
  start.state.int["LIPTf-3.2-0"] <- N*(prop.in.each["Lf-3.2-0"])
  start.state.int["LIPTf-3.2-1"] <- N*(prop.in.each["Lf-3.2-1"])
  start.state.int["LIPT-1.2-0"] <- N*(prop.in.each["Ls-1.2-0"])
  start.state.int["LIPT-1.2-1"] <- N*(prop.in.each["Ls-1.2-1"])
  start.state.int["LIPT-2.2-0"] <- N*(prop.in.each["Ls-2.2-0"])
  start.state.int["LIPT-2.2-1"] <- N*(prop.in.each["Ls-2.2-1"])
  start.state.int["LIPT-3.2-0"] <- N*(prop.in.each["Ls-3.2-0"])
  start.state.int["LIPT-3.2-1"] <- N*(prop.in.each["Ls-3.2-1"])
  
  start.state.con <- numeric(ncol(steady.state)-1)
  names(start.state.con) <- make.state.names()
  start.state.con[inds.of.interest-1]<- N*prop.in.each
  
  
  
  ## approx. steady state values
  # warning("check that foi and prob.MDR values are still sensible")
  foi <- params['beta']/100000 * DSprev
  prob.MDR <- params['prob.MDR']
  
  chg <- make.chg.mat(params)
  rc <- ode(y=as.numeric(start.state.int),
            times=seq(0,max.time,by=0.1),
            func=tb.dx.dt.trial,
            parms=params,
            chg=chg,
            lu=lu,
            foi=foi)
  
  colnames(rc) <- c("t",state.names)
  rc <- rc[,-1]
  
  ## now we will get the incidence and time at risk
  dt <- .1 # time step size
  
  ## LS -> A
  ii.11 <- sum(rc[,intersect(lu$ls.i,lu$h1art.i)]*params['epsilon_h']*params['h1theta']*params['arttheta']*dt)
  ii.12 <- sum(rc[,intersect(lu$ls.i,lu$h2art.i)]*params['epsilon_h']*params['h2theta']*params['arttheta']*dt)
  ii.13 <- sum(rc[,intersect(lu$ls.i,lu$h3art.i)]*params['epsilon_h']*params['arttheta']*dt)
  ii.21 <- sum(rc[,intersect(lu$lf.i,lu$h1art.i)]*params['pi_h']*params['h1theta']*params['arttheta']*dt)
  ii.22 <- sum(rc[,intersect(lu$lf.i,lu$h2art.i)]*params['pi_h']*params['h2theta']*params['arttheta']*dt)
  ii.23 <- sum(rc[,intersect(lu$lf.i,lu$h3art.i)]*params['pi_h']*params['arttheta']*dt)
  ii.31 <- sum(rc[,intersect(lu$lipt.i,lu$h1art.i)]*params['theta']*params['iepsilon_h']*params['h1theta']*params['arttheta']*dt)
  ii.32 <- sum(rc[,intersect(lu$lipt.i,lu$h2art.i)]*params['theta']*params['iepsilon_h']*params['h2theta']*params['arttheta']*dt)
  ii.33 <- sum(rc[,intersect(lu$lipt.i,lu$h3art.i)]*params['theta']*params['iepsilon_h']*params['arttheta']*dt)
  ii.41 <- sum(rc[,intersect(lu$lipt2.i,lu$h1art.i)]*params['theta']*params['iepsilon_h']*params['h1theta']*params['arttheta']*dt)
  ii.42 <- sum(rc[,intersect(lu$lipt2.i,lu$h2art.i)]*params['theta']*params['iepsilon_h']*params['h2theta']*params['arttheta']*dt)
  ii.43 <- sum(rc[,intersect(lu$lipt2.i,lu$h3art.i)]*params['theta']*params['iepsilon_h']*params['arttheta']*dt)
  ii.51 <- sum(rc[,intersect(lu$liptf.i,lu$h1art.i)]*params['theta']*params['ipi_h']*params['h1theta']*params['arttheta']*dt)
  ii.52 <- sum(rc[,intersect(lu$liptf.i,lu$h2art.i)]*params['theta']*params['ipi_h']*params['h2theta']*params['arttheta']*dt)
  ii.53 <- sum(rc[,intersect(lu$liptf.i,lu$h3art.i)]*params['theta']*params['ipi_h']*params['arttheta']*dt)
  ii.61 <- sum(rc[,intersect(lu$lipt2f.i,lu$h1art.i)]*params['theta']*params['ipi_h']*params['h1theta']*params['arttheta']*dt)
  ii.62 <- sum(rc[,intersect(lu$lipt2f.i,lu$h2art.i)]*params['theta']*params['ipi_h']*params['h2theta']*params['arttheta']*dt)
  ii.63 <- sum(rc[,intersect(lu$lipt2f.i,lu$h3art.i)]*params['theta']*params['ipi_h']*params['arttheta']*dt)
  
  new.cases.int <- ii.11+ii.12+ii.13+ii.21+ii.22+ii.23+ii.31+ii.32+ii.33+ii.41+ii.42+ii.43+ii.51+ii.52+ii.53+ii.61+ii.62+ii.63
  time.at.risk.int <- sum(rc[,-c(lu$a.i)]*dt)
  ir.int <- new.cases.int/time.at.risk.int
  
  ## now control gorup
  params['theta'] <- 1
  params['inf_on_ipt'] <- 1
  chg <- make.chg.mat(params)
  rc <- ode(y=as.numeric(start.state.con),
            times=seq(0,max.time,by=0.1),
            func=tb.dx.dt.trial,
            parms=params,
            chg=chg,
            lu=lu,
            foi=foi)
  
  colnames(rc) <- c("t",state.names)
  rc <- rc[,-1]
  ## now we will get the incidence and time at risk
  dt <- .1
  ## LS -> A
  ii.11 <- sum(rc[,intersect(lu$ls.i,lu$h1art.i)]*params['epsilon_h']*params['h1theta']*params['arttheta']*dt)
  ii.12 <- sum(rc[,intersect(lu$ls.i,lu$h2art.i)]*params['epsilon_h']*params['h2theta']*params['arttheta']*dt)
  ii.13 <- sum(rc[,intersect(lu$ls.i,lu$h3art.i)]*params['epsilon_h']*params['arttheta']*dt)
  ii.21 <- sum(rc[,intersect(lu$lf.i,lu$h1art.i)]*params['pi_h']*params['h1theta']*params['arttheta']*dt)
  ii.22 <- sum(rc[,intersect(lu$lf.i,lu$h2art.i)]*params['pi_h']*params['h2theta']*params['arttheta']*dt)
  ii.23 <- sum(rc[,intersect(lu$lf.i,lu$h3art.i)]*params['pi_h']*params['arttheta']*dt)
  ii.31 <- sum(rc[,intersect(lu$lipt.i,lu$h1art.i)]*params['theta']*params['iepsilon_h']*params['h1theta']*params['arttheta']*dt)
  ii.32 <- sum(rc[,intersect(lu$lipt.i,lu$h2art.i)]*params['theta']*params['iepsilon_h']*params['h2theta']*params['arttheta']*dt)
  ii.33 <- sum(rc[,intersect(lu$lipt.i,lu$h3art.i)]*params['theta']*params['iepsilon_h']*params['arttheta']*dt)
  ii.41 <- sum(rc[,intersect(lu$lipt2.i,lu$h1art.i)]*params['theta']*params['iepsilon_h']*params['h1theta']*params['arttheta']*dt)
  ii.42 <- sum(rc[,intersect(lu$lipt2.i,lu$h2art.i)]*params['theta']*params['iepsilon_h']*params['h2theta']*params['arttheta']*dt)
  ii.43 <- sum(rc[,intersect(lu$lipt2.i,lu$h3art.i)]*params['theta']*params['iepsilon_h']*params['arttheta']*dt)
  ii.51 <- sum(rc[,intersect(lu$liptf.i,lu$h1art.i)]*params['theta']*params['ipi_h']*params['h1theta']*params['arttheta']*dt)
  ii.52 <- sum(rc[,intersect(lu$liptf.i,lu$h2art.i)]*params['theta']*params['ipi_h']*params['h2theta']*params['arttheta']*dt)
  ii.53 <- sum(rc[,intersect(lu$liptf.i,lu$h3art.i)]*params['theta']*params['ipi_h']*params['arttheta']*dt)
  ii.61 <- sum(rc[,intersect(lu$lipt2f.i,lu$h1art.i)]*params['theta']*params['ipi_h']*params['h1theta']*params['arttheta']*dt)
  ii.62 <- sum(rc[,intersect(lu$lipt2f.i,lu$h2art.i)]*params['theta']*params['ipi_h']*params['h2theta']*params['arttheta']*dt)
  ii.63 <- sum(rc[,intersect(lu$lipt2f.i,lu$h3art.i)]*params['theta']*params['ipi_h']*params['arttheta']*dt)
  
  new.cases.cont <- ii.11+ii.12+ii.13+ii.21+ii.22+ii.23+ii.31+ii.32+ii.33+ii.41+ii.42+ii.43+ii.51+ii.52+ii.53+ii.61+ii.62+ii.63
  time.at.risk.cont <- sum(rc[,-c(lu$a.i)]*dt)
  ir.cont <- new.cases.cont/time.at.risk.cont
  
  ## now get the log liklihood
  # for control group incidence rate, target clin trial .036 (.028,.047)
  # sqrt(.028*.047)
  ll.ir <- log(dnorm((log(ir.int)-log(.036))/log(1.15)))
  # exp(qnorm(mean=log(0.036), sd=log(1.15), p=c(0.025,0.5,0.975)))
  # 94/63=63/42 = 1.5 = 1+2*(0.25)
  ll.irr <- log(dnorm((log(ir.int/ir.cont)-log(.63))/log(1.22)))
  # exp(qnorm(mean=log(0.63), sd=log(1.22), p=c(0.025,0.5,0.975)))
  
            
  return(list(ll.ir, ll.irr))
}


# recalibrate model to a (4x) lower-incidence setting (allowing wider CIs to improve fit), otherwise similar: 

## log likelihood for 'fitting';
## change intervals to assymetrical
log.lik.low <- function(data){
  ##incidence, assuming +-10% CI
  ll.not <- log(dnorm((log(data['tb.notified']) - log(1630/4))/log(1.1)))
  
  ll.prev <- log(dnorm(log(data['tb.prev']/(3200/4))/log(1.2)))
  
  ll.mort <- log(dnorm(log(data['tb.mort']/(500/4))/log(1.3)))
  
  ll.mdr <- log(dnorm((logit(data['prop.MDR']) - logit(0.04))/((logit(0.07)-logit(0.02))/4)))
  
  
  # 0.8*c(33, 28, 39)
  ms <-  twCoefLogitnormCi(0.2224, 0.312)
  ll.hiv <- log(dnorm((logit(data['hiv.prev'])-ms[1])/ms[2]))
  
  ##mean CD4 in ART initiators, assuming +-20% CI
  ll.cd4 <- log(dnorm((log(data['cd4.artinit']) - log(175))/log(1.1)))
  # check: exp(qnorm(c(0.025,0.5,0.975), log(175), log(1.1)))
  
  
  
  # library(binom)
  # binom.confint(211+39,2138)
  ms <-  twCoefLogitnormCi(0.104/4, 0.131/4)
  # qlogitnorm(p = 0.5, mu = ms[1], sigma=ms[2])
  ll.ipt.elg.tb <- log(dnorm((logit(data['tb.prev.on.art']) - ms[1])/ms[2]))
  
  ll.vec <- c(ll.not,ll.prev,ll.mort,ll.mdr,ll.hiv,ll.ipt.elg.tb, ll.cd4)
  ll.vec <- pmax(ll.vec,rep(-1000,length(ll.vec)))
  return(ll.vec)
}


log.lik.ipt.low <- function(run,params,max.time=3){
  
  lu <- make.look.up.list()
  state.names <- make.state.names()
  
  params['iepsilon_0'] <- params['epsilon_0']
  params['iepsilon_h'] <- params['epsilon_h']
  params['ipi_0'] <- params['pi_0']
  params['ipi_h'] <- params['pi_h']
  
  steady.state <- continue.run(
    run,
    params=params,
    extra.t = 1,by=1,print = F) # turned printing on to see foi and probMDR
  
  DSprev <- sum(steady.state[nrow(steady.state), 1+ intersect(lu$a.i,lu$d0.i)])
  
  ## sample size
  N <- 622
  
  ## figure out distribtion of people in ART states
  inds.of.interest <- setdiff(lu$allart.i,lu$a.i)+1
  cols.of.interest <- steady.state[nrow(steady.state),inds.of.interest]
  prop.in.each <- cols.of.interest/sum(cols.of.interest)
  
  start.state.int <- numeric(ncol(steady.state)-1)
  names(start.state.int) <- make.state.names()
  start.state.int["SIPT-1.2"] <- N*prop.in.each["S-1.2"]
  start.state.int["SIPT-2.2"] <- N*prop.in.each["S-2.2"]
  start.state.int["SIPT-3.2"] <- N*prop.in.each["S-3.2"]
  start.state.int["LIPTf-1.2-0"] <- N*(prop.in.each["Lf-1.2-0"])
  start.state.int["LIPTf-1.2-1"] <- N*(prop.in.each["Lf-1.2-1"])
  start.state.int["LIPTf-2.2-0"] <- N*(prop.in.each["Lf-2.2-0"])
  start.state.int["LIPTf-2.2-1"] <- N*(prop.in.each["Lf-2.2-1"])
  start.state.int["LIPTf-3.2-0"] <- N*(prop.in.each["Lf-3.2-0"])
  start.state.int["LIPTf-3.2-1"] <- N*(prop.in.each["Lf-3.2-1"])
  start.state.int["LIPT-1.2-0"] <- N*(prop.in.each["Ls-1.2-0"])
  start.state.int["LIPT-1.2-1"] <- N*(prop.in.each["Ls-1.2-1"])
  start.state.int["LIPT-2.2-0"] <- N*(prop.in.each["Ls-2.2-0"])
  start.state.int["LIPT-2.2-1"] <- N*(prop.in.each["Ls-2.2-1"])
  start.state.int["LIPT-3.2-0"] <- N*(prop.in.each["Ls-3.2-0"])
  start.state.int["LIPT-3.2-1"] <- N*(prop.in.each["Ls-3.2-1"])
  
  start.state.con <- numeric(ncol(steady.state)-1)
  names(start.state.con) <- make.state.names()
  start.state.con[inds.of.interest-1]<- N*prop.in.each
  
  
  
  ## approx. steady state values
  warning("check that foi and prob.MDR values are still sensible")
  foi <- params['beta']/100000 * DSprev
  prob.MDR <- params['prob.MDR']

  chg <- make.chg.mat(params)
  rc <- ode(y=as.numeric(start.state.int),
            times=seq(0,max.time,by=0.1),
            func=tb.dx.dt.trial,
            parms=params,
            chg=chg,
            lu=lu,
            foi=foi)
  
  colnames(rc) <- c("t",state.names)
  rc <- rc[,-1]
  
  ## now we will get the incidence and time at risk
  dt <- .1 # time step size
  
  ## LS -> A
  ii.11 <- sum(rc[,intersect(lu$ls.i,lu$h1art.i)]*params['epsilon_h']*params['h1theta']*params['arttheta']*dt)
  ii.12 <- sum(rc[,intersect(lu$ls.i,lu$h2art.i)]*params['epsilon_h']*params['h2theta']*params['arttheta']*dt)
  ii.13 <- sum(rc[,intersect(lu$ls.i,lu$h3art.i)]*params['epsilon_h']*params['arttheta']*dt)
  ii.21 <- sum(rc[,intersect(lu$lf.i,lu$h1art.i)]*params['pi_h']*params['h1theta']*params['arttheta']*dt)
  ii.22 <- sum(rc[,intersect(lu$lf.i,lu$h2art.i)]*params['pi_h']*params['h2theta']*params['arttheta']*dt)
  ii.23 <- sum(rc[,intersect(lu$lf.i,lu$h3art.i)]*params['pi_h']*params['arttheta']*dt)
  ii.31 <- sum(rc[,intersect(lu$lipt.i,lu$h1art.i)]*params['theta']*params['iepsilon_h']*params['h1theta']*params['arttheta']*dt)
  ii.32 <- sum(rc[,intersect(lu$lipt.i,lu$h2art.i)]*params['theta']*params['iepsilon_h']*params['h2theta']*params['arttheta']*dt)
  ii.33 <- sum(rc[,intersect(lu$lipt.i,lu$h3art.i)]*params['theta']*params['iepsilon_h']*params['arttheta']*dt)
  ii.41 <- sum(rc[,intersect(lu$lipt2.i,lu$h1art.i)]*params['theta']*params['iepsilon_h']*params['h1theta']*params['arttheta']*dt)
  ii.42 <- sum(rc[,intersect(lu$lipt2.i,lu$h2art.i)]*params['theta']*params['iepsilon_h']*params['h2theta']*params['arttheta']*dt)
  ii.43 <- sum(rc[,intersect(lu$lipt2.i,lu$h3art.i)]*params['theta']*params['iepsilon_h']*params['arttheta']*dt)
  ii.51 <- sum(rc[,intersect(lu$liptf.i,lu$h1art.i)]*params['theta']*params['ipi_h']*params['h1theta']*params['arttheta']*dt)
  ii.52 <- sum(rc[,intersect(lu$liptf.i,lu$h2art.i)]*params['theta']*params['ipi_h']*params['h2theta']*params['arttheta']*dt)
  ii.53 <- sum(rc[,intersect(lu$liptf.i,lu$h3art.i)]*params['theta']*params['ipi_h']*params['arttheta']*dt)
  ii.61 <- sum(rc[,intersect(lu$lipt2f.i,lu$h1art.i)]*params['theta']*params['ipi_h']*params['h1theta']*params['arttheta']*dt)
  ii.62 <- sum(rc[,intersect(lu$lipt2f.i,lu$h2art.i)]*params['theta']*params['ipi_h']*params['h2theta']*params['arttheta']*dt)
  ii.63 <- sum(rc[,intersect(lu$lipt2f.i,lu$h3art.i)]*params['theta']*params['ipi_h']*params['arttheta']*dt)
  
  new.cases.int <- ii.11+ii.12+ii.13+ii.21+ii.22+ii.23+ii.31+ii.32+ii.33+ii.41+ii.42+ii.43+ii.51+ii.52+ii.53+ii.61+ii.62+ii.63
  time.at.risk.int <- sum(rc[,-c(lu$a.i)]*dt)
  ir.int <- new.cases.int/time.at.risk.int
  
  ## now control gorup
  params['theta'] <- 1
  params['inf_on_ipt'] <- 1
  chg <- make.chg.mat(params)
  rc <- ode(y=as.numeric(start.state.con),
            times=seq(0,max.time,by=0.1),
            func=tb.dx.dt.trial,
            parms=params,
            chg=chg,
            lu=lu,
            foi=foi,
            prob.MDR=prob.MDR)
  
  colnames(rc) <- c("t",state.names)
  rc <- rc[,-1]
  ## now we will get the incidence and time at risk
  dt <- .1
  ## LS -> A
  ii.11 <- sum(rc[,intersect(lu$ls.i,lu$h1art.i)]*params['epsilon_h']*params['h1theta']*params['arttheta']*dt)
  ii.12 <- sum(rc[,intersect(lu$ls.i,lu$h2art.i)]*params['epsilon_h']*params['h2theta']*params['arttheta']*dt)
  ii.13 <- sum(rc[,intersect(lu$ls.i,lu$h3art.i)]*params['epsilon_h']*params['arttheta']*dt)
  ii.21 <- sum(rc[,intersect(lu$lf.i,lu$h1art.i)]*params['pi_h']*params['h1theta']*params['arttheta']*dt)
  ii.22 <- sum(rc[,intersect(lu$lf.i,lu$h2art.i)]*params['pi_h']*params['h2theta']*params['arttheta']*dt)
  ii.23 <- sum(rc[,intersect(lu$lf.i,lu$h3art.i)]*params['pi_h']*params['arttheta']*dt)
  ii.31 <- sum(rc[,intersect(lu$lipt.i,lu$h1art.i)]*params['theta']*params['iepsilon_h']*params['h1theta']*params['arttheta']*dt)
  ii.32 <- sum(rc[,intersect(lu$lipt.i,lu$h2art.i)]*params['theta']*params['iepsilon_h']*params['h2theta']*params['arttheta']*dt)
  ii.33 <- sum(rc[,intersect(lu$lipt.i,lu$h3art.i)]*params['theta']*params['iepsilon_h']*params['arttheta']*dt)
  ii.41 <- sum(rc[,intersect(lu$lipt2.i,lu$h1art.i)]*params['theta']*params['iepsilon_h']*params['h1theta']*params['arttheta']*dt)
  ii.42 <- sum(rc[,intersect(lu$lipt2.i,lu$h2art.i)]*params['theta']*params['iepsilon_h']*params['h2theta']*params['arttheta']*dt)
  ii.43 <- sum(rc[,intersect(lu$lipt2.i,lu$h3art.i)]*params['theta']*params['iepsilon_h']*params['arttheta']*dt)
  ii.51 <- sum(rc[,intersect(lu$liptf.i,lu$h1art.i)]*params['theta']*params['ipi_h']*params['h1theta']*params['arttheta']*dt)
  ii.52 <- sum(rc[,intersect(lu$liptf.i,lu$h2art.i)]*params['theta']*params['ipi_h']*params['h2theta']*params['arttheta']*dt)
  ii.53 <- sum(rc[,intersect(lu$liptf.i,lu$h3art.i)]*params['theta']*params['ipi_h']*params['arttheta']*dt)
  ii.61 <- sum(rc[,intersect(lu$lipt2f.i,lu$h1art.i)]*params['theta']*params['ipi_h']*params['h1theta']*params['arttheta']*dt)
  ii.62 <- sum(rc[,intersect(lu$lipt2f.i,lu$h2art.i)]*params['theta']*params['ipi_h']*params['h2theta']*params['arttheta']*dt)
  ii.63 <- sum(rc[,intersect(lu$lipt2f.i,lu$h3art.i)]*params['theta']*params['ipi_h']*params['arttheta']*dt)
  
  new.cases.cont <- ii.11+ii.12+ii.13+ii.21+ii.22+ii.23+ii.31+ii.32+ii.33+ii.41+ii.42+ii.43+ii.51+ii.52+ii.53+ii.61+ii.62+ii.63
  time.at.risk.cont <- sum(rc[,-c(lu$a.i)]*dt)
  ir.cont <- new.cases.cont/time.at.risk.cont
  
  ## now get the log liklihood
  # for control group incidence rate, target clin trial .036 (.028,.047)
  # sqrt(.028*.047)
  ll.ir <- 0#log(dnorm((log(ir.int)-log(.036))/log(1.15)))
  # exp(qnorm(mean=log(0.036), sd=log(1.15), p=c(0.025,0.5,0.975)))
  # 94/63=63/42 = 1.5 = 1+2*(0.25)
  ll.irr <- log(dnorm((log(ir.int/ir.cont)-log(.63))/log(1.22)))
  # exp(qnorm(mean=log(0.63), sd=log(1.22), p=c(0.025,0.5,0.975)))
  
  
  return(list(ll.ir, ll.irr)) # for low, no variation in ir target, only assume same irr
}

