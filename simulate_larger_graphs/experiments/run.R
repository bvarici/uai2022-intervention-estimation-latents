# Copyright (c) 2018-2020, Joris M. Mooij. All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# Use this script to run a method on simulated data

suppressMessages(library(rjson))
source('../../init.R',chdir=TRUE)

# parse command line arguments
args <- commandArgs(trailingOnly = TRUE)


basefilename<-args[1]

# read metadata
#metadatafile=paste(basefilename,'.json',sep='')
metadatafile=paste(basefilename,'-metadata.json',sep='')
print(metadatafile)
metadata<-fromJSON(file=metadatafile)
#pSysObs<-metadata$pSysObs
#stopifnot( pSysObs > 0 )
#pContext<-metadata$pContext
systemVars<-metadata$SystemVars
contextVars<-metadata$ContextVars

pSysObs<-length(systemVars)
stopifnot( pSysObs > 0 )
pContext<-length(contextVars)
p <- pSysObs + pContext
obsContext <- matrix(0,1,pContext)

#alg<-'fci'
alg<-args[2]
#mode<-'obs'
mode<-args[3]
miniter<-0
if( !is.na(args[4]) ) {
  miniter<-as.numeric(args[4])
}
maxiter<-0
if( !is.na(args[5]) ) {
  maxiter<-as.numeric(args[5])
}
alpha<-1e-2
if( !is.na(args[6]) )
  alpha<-as.numeric(args[6])
jcifci_test<-'gaussCIcontexttest'
if( !is.na(args[7]) )
  jcifci_test<-args[7]
doPDsep<-TRUE
if( !is.na(args[8]) )
  doPDsep<-(as.numeric(args[8]) != 0)
verbose<-0
if( !is.na(args[9]) )
  verbose<-(as.integer(args[9]))

cat('run.R: running with arguments',args,'\n')

# read and preprocess data
data<-read.csv(file=paste(basefilename,'-data.csv',sep=''),header=TRUE,sep=",")

for( iter in miniter:maxiter ) {
  cat('iter: ',iter,'\n')

set.seed(iter)
if( iter != 0 ) {
  subsamplefrac<-0.5
} else {
  subsamplefrac<-0.0
}

# start measuring time
start_time<-proc.time()[3]

# run causal-discovery
if( alg == 'fci' ) { # run FCI
  result<-fci_wrapper(data,systemVars,contextVars,alpha,verbose=verbose,subsamplefrac,test=jcifci_test,mode,obsContext,doPDsep)

  # save results
  #ifelse(!dir.exists(outdir), dir.create(outdir), FALSE)
  if( iter != 0 )
    outfile<-paste(basefilename,'-',alg,'-',mode,'-',iter,sep='')
  else
    outfile<-paste(basefilename,'-',alg,'-',mode,sep='')
  save(result,file=paste(outfile,'-result.Rdata',sep=''))
  write.csv(result$edge,file=paste(outfile,'-edge.csv',sep=''),row.names=FALSE)

} else
  stop('Unknown algorithm')

# stop measuring time
stop_time<-proc.time()[3]

if( exists('outfile') ) {
  # write time to file
  cat(file=paste(outfile,'-runtime.csv',sep=''),stop_time-start_time,'\n')
}

}
