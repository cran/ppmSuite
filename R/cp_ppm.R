## Fortran subroutines

## Compile: Fortran codes
## R CMD SHLIB File1.f  ... FileN.f  File1.c ... FileN.c
## -o FileName
## If ToolsLA.f is needed: add -llapack -lblas after FileName


## "MCMCLOGITP" subroutine
ccp_ppm = function(ydata,nu0,mu0,sigma0,mltypes,
                       thetas,devs,nburn,nskip,nsave, verbose=FALSE){

  ndata = dim(ydata)[2]
  nseries = dim(ydata)[1]
  invsigma0 = solve(sigma0)
  logdetsigma0 = log(det(sigma0))
  nthetas = 4
  thetas = thetas[,1:nthetas]
  mcmcc = matrix(2,nrow=nsave,ncol=((ndata-1)*nseries))
  mcmcp = matrix(0,nrow=nsave,ncol=((ndata-1)*nseries))
#  print(thetas)
  foo = .Fortran("mcmclogitp",
                 nburn=as.integer(nburn),
                 nskip=as.integer(nskip),
                 nsave=as.integer(nsave),
                 ndata=as.integer(ndata),
                 nseries=as.integer(nseries),
                 ydata=as.double(ydata),
                 nu0=as.double(nu0),
                 mu0=as.double(mu0),
                 invsigma0=as.double(invsigma0),
                 logdetsigma0=as.double(logdetsigma0),
                 mltypes=as.integer(mltypes),
                 nthetas=as.integer(nthetas),
                 thetas=as.double(thetas),
                 devs=as.double(devs),
                 verbose=as.integer(verbose),
                 mcmcc=as.integer(mcmcc),
                 mcmcp=as.double(mcmcp))

  C = matrix(foo$mcmcc,nrow=nsave,byrow=FALSE)
  P = matrix(foo$mcmcp,nrow=nsave,byrow=FALSE)

  return(list(C=C,P=P))
}



## "MCMCBETAS" subroutine
icp_ppm = function(ydata,a0,b0,mltypes,thetas,nburn,nskip,nsave,verbose=FALSE){

  ndata = dim(ydata)[2]
  nseries = dim(ydata)[1]
  nthetas = 4
  thetas = thetas[,1:nthetas]
  mcmcc = matrix(2,nrow=nsave,ncol=((ndata-1)*nseries))
  mcmcp = matrix(0,nrow=nsave,ncol=nseries)

  foo = .Fortran("mcmcbetas",
                 nburn=as.integer(nburn),
                 nskip=as.integer(nskip),
                 nsave=as.integer(nsave),
                 ndata=as.integer(ndata),
                 nseries=as.integer(nseries),
                 ydata=as.double(ydata),
                 a0=as.double(a0),
                 b0=as.double(b0),
                 mltypes=as.integer(mltypes),
                 nthetas=as.integer(nthetas),
                 thetas=as.double(thetas),
                 verbose=as.integer(verbose),
                 mcmcc=as.integer(mcmcc),
                 mcmcp=as.double(mcmcp))

  C = matrix(foo$mcmcc,nrow=nsave,byrow=FALSE)
  P = matrix(foo$mcmcp,nrow=nsave,byrow=FALSE)

  return(list(C=C,P=P))
}
