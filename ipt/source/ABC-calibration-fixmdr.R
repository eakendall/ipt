## script to run LHS to get weights
library(parallel)
source("source/khay-utility-fixmdr.R")
reload.source()
# set.seed(3423235)
cores=detectCores()
savelocation <- "../scratch/khaydata/"
date <- "20171213" 
n <- as.numeric(commandArgs(trailingOnly=TRUE)[1])
newlabel <- commandArgs(trailingOnly=TRUE)[2]
toaugment <- commandArgs(trailingOnly=TRUE)[3]
# toreuse <- commandArgs(trailingOnly=TRUE)[3]

if (is.na(toaugment)) toaugment <- ""

print(paste0(c(n,".",toaugment,".",newlabel)))

print(cores)
## draw from LHS cube
draws <- get.LHS.draws(n.draws = n, augment=toaugment, reuse="")

## for each set of params run the model and calculate the log likelihood
## and steady-state
lls <- mclapply(1:nrow(draws),function(x) {print(x); one.run(draws[x,],300)}, mc.cores=cores)

saveRDS(list("draws"=draws,"lls"=lls),file=paste0(savelocation,"generated_data/LHSdraws_and_lls.",date,".rds"))



