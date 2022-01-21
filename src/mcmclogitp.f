c=======================================================================
c=======================================================================
      subroutine mcmclogitp(nburn,nskip,nsave,ndata,nseries,ydata,
     & nu0,mu0,invsigma0,logdetsigma0,mltypes,nthetas,thetas,devs,
     & verbose,mcmcc,mcmcp)
c=======================================================================
c=======================================================================
c     BEGIN: MCMCLOGITP SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     "MCMCLOGITP" RETURNS POSTERIOR SAMPLES FROM THE FOLLOWING
c     BAYESIAN MODEL FOR CHANGE-POINT DETECTION:
c
c     A) PRIOR SPECIFICATION FOR CHANGE-POINT INDICATORS:
c
c        C=(C[I,T]:I=1,...,L AND T=1,...,N-1).
c
c        HERE, C[I,T]=1 IF TIME T+1 IN (Y[I,S]:S=1,...,N) IS A
c        CHANGE-POINT.
c
c     A.1) LET P=(P[I,T]:I=1,...,L AND T=1,...,N-1) BE A MATRIX WITH
c          ENTRIES ON [0,1]. THEN,
c
c          C[I,T]|P[I,T]~BERN(P[I,T]).
c
c          · BERN(·|P): BERNOULLI P.M.F. WITH PARAMETER P IN [0,1].
c
c     A.2) THE VECTORS (P[I,T]:I=1,...,L) ARE I.I.D. WITH COMMON
c          DISTRIBUTION LOGIT-T[L](NU[0],MU[0],SIGMA[0]), I.E.,
c
c          (LOGIT(P[I,T]):I=1,...,L)~T[L](NU[0],MU[0],SIGMA[0]).
c
c          HERE, T[L](NU[0],MU[0],SIGMA[0]) IS THE L-DIMENSIONAL
c          STUDENT DISTRIBUTION WITH DEGREES OF FREEDOM NU[0],
c          LOCATION VECTOR MU[0] AND SCALE MATRIX SIGMA[0].
c
c     B) DATA GENERATING MECHANISM: LET
c
c        E=(E[I,T]:I=1,...,L AND T=1,...,N).
c
c        HERE, E[I,1]=1 AND E[I,T+1]=E[I,T]+C[I,T], FOR ALL I=1,...,L
c        AND T=1,...,N-1. WITH THIS INFORMATION, (E[I,T]:T=1,...,N)
c        ARE THE CLUSTER LABELS FOR (Y[I,T]:T=1,...,N) AND
c        K[I]=1+SUM(C[I,T]:T=1,...,N-1) IS THE NUMBER OF CLUSTERS.
c
c     B.1) (Y[I,T]:T=1,...,N)|(E[I,T]:T=1,...,N)~
c          PRODUCT(ML[I](Y[I,T]:E[I,T]=J|THETA[I]):J=1,...,K[I]).
c
c          · ML[I](·|THETA[I]): MARGINAL LIKELIHOOD FOR
c            (Y(I,T):T=1,...,N), INDEXED BY A PARAMETER THETA[I].
c=======================================================================
c     INPUTS: MCMC
c=======================================================================
c     nburn: BURN-IN ITERATIONS
      integer nburn
c     nskip: SKIP BETWEEN SAVED ITERATIONS
      integer nskip
c     nsave: SAVED ITERATIONS
      integer nsave
c     verbose: PRINT TO SCREEN PROGRESSION OF MCMC
      integer verbose      
c=======================================================================
c     INPUTS: DATA
c=======================================================================
c     ndata: LENGTH OF THE TEMPORAL AXIS (N)
      integer ndata
c     nseries: NUMBER OF TIME SERIES (L)
      integer nseries
c     ydata: TIME SERIES (Y[I,T]:I=1,...,L AND T=1,...,N)
      real(8) ydata(nseries,ndata)
c=======================================================================
c     INPUTS: HYPER-PARAMETERS
c=======================================================================
c     nu0: DEGREES OF FREEDOM (NU[0])
      real(8) nu0
c     mu0: LOCATION VECTOR (MU[0])
      real(8) mu0(nseries)
c     invsigma0: INVERSE OF SCALE MATRIX (SIGMA[0])
      real(8) invsigma0(nseries,nseries)
c     logdetsigma0: LOG DETERMINANT OF SCALE MATRIX (SIGMA[0])
      real(8) logdetsigma0
c     mltypes: MARGINAL LIKELIHOOD TYPES (ML[I](·):T=1,...,N)
      integer mltypes(nseries)
c     nthetas: MAXIMUM OF {DIMENSION(THETA[I]):I=1,...,L}
      integer nthetas
c     thetas: PARAMETERS FOR EACH (ML[I](·):I=1,...,L)
      real(8) thetas(nseries,nthetas)
c     devs: STANDARD DEVIATIONS FOR RW-MH UPDATES
      real(8) devs(nseries,ndata-1)
      
      
c=======================================================================
c     OUTPUTS: MCMC SAMPLES
c=======================================================================
c     mcmcc: (C[I,T]:I=1,...,L AND T=1,...,N-1)
      integer mcmcc(nsave,(ndata-1)*nseries)
c     mcmcp: (P[I,T]:I=1,...,L AND T=1,...,N-1)
      real(8) mcmcp(nsave,(ndata-1)*nseries)
c=======================================================================
c     C++ FUNCTIONS
c=======================================================================
c     FUNCTIONS IN "TOOLSR2.C"
      real(8) binomd
      real(8) normr
      real(8) unifr
c=======================================================================
c     FORTRAN SUBROUTINES
c=======================================================================
c     SUBROUTINES IN "TOOLSGS2.F"
c     logml(···)
c     SUBROUTINES IN "TOOLSPD2.F"
c     logmvtd(···)
c=======================================================================
c     WORKING VARIABLES: MCMC
c=======================================================================
c     actual: COUNTER (ONGOING ITERATION)
      integer actual
c     stored: COUNTER (STORED ITERATION)
      integer stored
c=======================================================================
c     WORKING VARIABLES 1
c=======================================================================
c     EXCLUSIVE FOR STORING (Y[I,T]:T=1,...,N)
      real(8) y(ndata)
c     EXCLUSIVE FOR STORING (THETA[I]:I=1,...,L)
      real(8) theta(nthetas)
c     EXCLUSIVE FOR STORING LOG RATIOS
      real(8) llo
      real(8) lln
      real(8) llr
c     EXCLUSIVE FOR UPDATING (C[I,T]:I=1,...,L AND T=1,...,N-1)
      integer co(ndata-1)
      integer cn(ndata-1)
      integer eo(ndata)
      integer en(ndata)
      integer c(nseries,ndata-1)
c     EXCLUSIVE FOR UPDATING (P[I,T]:I=1,...,L AND T=1,...,N-1)
      real(8) po
      real(8) pn
      real(8) logitpo(nseries)
      real(8) logitpn(nseries)
      real(8) p(nseries,ndata-1)
c=======================================================================
c     WORKING VARIABLES 2
c=======================================================================
c     INDEXES
      integer i
      integer j
      integer k
      integer s
      integer t
c     OTHERS (REAL)
      real(8) cw
      real(8) rw
      real(8) sw
      real(8) uw
c=======================================================================
c     SETTING INITIAL VALUES
c=======================================================================
c     MCMC COUNTERS
      actual=1
      stored=1
c     CHANGE-POINT INDICATORS (C[I,T]=0:I=1,...,L AND T=1,...,N-1)
      do i=1,nseries
         do t=1,(ndata-1)
            c(i,t)=0
         end do
      end do
c     PROBABILITIES (P[I,T]:I=1,...,L AND T=1,...,N-1)
      do t=1,(ndata-1)
         do i=1,nseries
            p(i,t)=1.d0/dble(ndata)
c            p(i,t)=0.5
         end do
      end do
c=======================================================================
c     METROPOLIS-HASTINGS-WITHIN-GIBBS ALGORITHM
c=======================================================================
c     PRINT ON SCREEN: BEGINING OF MCMC ITERATIONS
      if(verbose.eq.1) then
c        print *,'Begining of MCMC iterations'
      end if
c=======================================================================
      call rndstart()
c     BEGIN: ITERATIONS
      do while (actual.le.(nburn+(nskip*nsave)))
c=======================================================================
c        UPDATING CHANGE-POINT INDICATORS
c=======================================================================
         do i=1,nseries
            k=mltypes(i)
            do s=1,ndata
               y(s)=ydata(i,s)
            end do
            do s=1,nthetas
               theta(s)=thetas(i,s)
            end do
            
            do t=1,(ndata-1)
               do s=1,(ndata-1)
                  co(s)=c(i,s)
                  cn(s)=c(i,s)
               end do
               co(t)=0
               cn(t)=1
               eo(1)=1
               en(1)=1
               do s=1,(ndata-1)
                  eo(s+1)=eo(s)+co(s)
                  en(s+1)=en(s)+cn(s)
               end do
               sw=0.d0
               j=en(t)
               call logml(k,ndata,y,nthetas,theta,en,j,rw)
c               print *, "ml_left", rw
               sw=sw+rw
               j=en(t+1)
               call logml(k,ndata,y,nthetas,theta,en,j,rw)
c               print *, "ml_right", rw
               sw=sw+rw
               j=eo(t)
c               print *, "ml_left + ml_right", sw
               call logml(k,ndata,y,nthetas,theta,eo,j,rw)
c               print *, "ml_all", rw
               sw=sw-rw
               llr=sw+(dlog(p(i,t))-dlog(1.d0-p(i,t)))
c               print *, "llr", llr
               uw=unifr(0.d0,1.d0)
               rw=dlog(uw)-dlog(1.d0-uw)
c               print *, "rw", rw
               if (llr.gt.rw) then
                  c(i,t)=1
               else
                  c(i,t)=0
               end if
c               print *, "c(i,t)", c(i,t)
            end do
         end do
         
c=======================================================================
c        UPDATING PROBABILITY PARAMETERS
c=======================================================================
         do t=1,(ndata-1)
c            print *, "t = ", t
            do i=1,nseries
c               print *, "i = ", i
c               print *, "devs(i,t) = ", devs(i,t)
               po=p(i,t)
               pn=normr(po,devs(i,t))
c               print *, "po", po
c               print *, "pn", pn
c               print *, "c(i,t)", c(i,t)
               if ((pn.gt.0.d0).and.(pn.lt.1.d0)) then
                  do s=1,nseries
                     logitpo(s)=dlog(p(s,t))-dlog(1.d0-p(s,t))
                     logitpn(s)=logitpo(s)
                  end do
                  logitpn(i)=dlog(pn)-dlog(1.d0-pn)
c                  print *, "logitpo", logitpo(1), logitpo(2)
                  cw=dble(c(i,t))
                  call logmvtd(nseries,logitpo,nu0,mu0,invsigma0,
     &                 logdetsigma0,rw)
c                  print *, "logmvtdo", rw
                  sw=binomd(cw,1.d0,po,1)
c                  print *, "binomdo", sw
                  llo=(rw+sw)-(dlog(po)+dlog(1.d0-po))
c                  llo=cw*dlog(po)+(1.d0-cw)*dlog(1.d0-po) + rw - 
c     &                    (dlog(po)+dlog(1.d0-po)) 
c                  print *, "llo", llo
c                  print *, "logitpn", logitpn(1), logitpn(2), logitpn(3)
c                  print *, "nu0",nu0
c                  print *, "mu0", mu0(1), mu0(2), mu0(3)
                  call logmvtd(nseries,logitpn,nu0,mu0,invsigma0,
     &                 logdetsigma0,rw)
c                  print *, "logmvtdn", rw
                  sw=binomd(cw,1.d0,pn,1)
c                  print *, "binomdn", sw
                  lln=(rw+sw)-(dlog(pn)+dlog(1.d0-pn))
c                  lln=cw*dlog(pn) + (1.d0-cw)*dlog(1.d0-pn) + 
c     &                 rw - (dlog(pn)+dlog(1.d0-pn))
c                  print *, "lln", lln
                  llr=lln-llo
c                  print *, "llr", llr
                  uw=unifr(0.d0,1.d0)
                  rw=dlog(uw)
c                  print *, "rw", rw
                  if (llr.gt.rw) then
                     p(i,t)=pn
                  end if
               end if
c               print *, "p(i,t)", p(i,t)
            end do
         end do
c=======================================================================
c        PRINT ON SCREEN: BURN-IN PHASE COMPLETED
         if(verbose.eq.1) then
           if (actual.eq.nburn) then
c             print *,'Burn-in phase completed'
           end if
         end if
c=======================================================================
c        PROCESSING: SAVED ITERATIONS
c        - AFTER NBURN ITERATIONS, NSAVE ITERATIONS BETWEEN NSKIP
c          STEPS
c=======================================================================
         if ((actual.gt.nburn).and.(mod(actual-nburn,nskip).eq.0)) then
c=======================================================================
c           STORING: MCMC SAMPLES
c=======================================================================
c           CHANGE-POINT INDICATORS AND PROBABILITIES
            s=0
            do i=1,nseries
               do t=1,(ndata-1)
                  s=s+1
                  mcmcc(stored,s)=c(i,t)
                  mcmcp(stored,s)=p(i,t)
               end do
            end do
c=======================================================================
c           PRINT ON SCREEN: STORED ITERATION (MULTIPLES OF 100)
            if(verbose.eq.1) then
              if (mod(stored,100).eq.0) then
c                 print *,stored,'stored MCMC iterations'
              end if
            end if
c=======================================================================
c           UPDATING: COUNTER (STORED ITERATION)
c=======================================================================
            stored=stored+1
         end if
c=======================================================================
c        UPDATING: COUNTER (ONGOING ITERATION)
c=======================================================================
         actual=actual+1
      end do
c     END: MCMC ITERATIONS
      call rndend()
c=======================================================================
c     PRINT ON SCREEN: END OF MCMC ITERATIONS
      if(verbose.eq.1) then
c        print *,'End of MCMC iterations'
      end if
c=======================================================================
      return
c     END: MCMCLOGITP SUBROUTINE
      end
c=======================================================================
c=======================================================================
