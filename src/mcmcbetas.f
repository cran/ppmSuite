c=======================================================================
c=======================================================================
      subroutine mcmcbetas(nburn,nskip,nsave,ndata,nseries,ydata,
     & a0,b0,mltypes,nthetas,thetas,verbose,mcmcc,mcmcp)
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
c     A.1) LET P=(P[I]:I=1,...,L) BE A VECTOR WITH ENTRIES ON [0,1].}
c          THEN,
c
c          C[I,T]|P[I]~BERN(P[I]).
c
c          · BERN(·|P): BERNOULLI P.M.F. WITH PARAMETER P IN [0,1].
c
c     A.2) THE SCALARS (P[I]:I=1,...,L) ARE INDEPENDENT, WHERE
c          P[I]~BETA(A[I,0],B[I,0]).
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
c     a0: SHAPE PARAMETERS (A[I,0]:I=1,...,L)
      real(8) a0(nseries)
c     b0: SHAPE PARAMETERS (B[I,0]:I=1,...,L)
      real(8) b0(nseries)
c     mltypes: MARGINAL LIKELIHOOD TYPES (ML[I](·):T=1,...,N)
      integer mltypes(nseries)
c     nthetas: MAXIMUM OF {DIMENSION(THETA[I]):I=1,...,L}
      integer nthetas
c     thetas: PARAMETERS FOR EACH (ML[I](·):I=1,...,L)
      real(8) thetas(nseries,nthetas)
c=======================================================================
c     OUTPUTS: MCMC SAMPLES
c=======================================================================
c     mcmcc: (C[I,T]:I=1,...,L AND T=1,...,N-1)
      integer mcmcc(nsave,(ndata-1)*nseries)
c     mcmcp: (P[I,T]:I=1,...,L AND T=1,...,N-1)
      real(8) mcmcp(nsave,nseries)
c=======================================================================
c     C++ FUNCTIONS
c=======================================================================
c     FUNCTIONS IN "TOOLSR2.C"
      real(8) betar
      real(8) unifr
c=======================================================================
c     FORTRAN SUBROUTINES
c=======================================================================
c     SUBROUTINES IN "TOOLSGS2.F"
c     logml(···)
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
      real(8) llr
c     EXCLUSIVE FOR UPDATING (C[I,T]:I=1,...,L AND T=1,...,N-1)
      integer co(ndata-1)
      integer cn(ndata-1)
      integer eo(ndata)
      integer en(ndata)
      integer c(nseries,ndata-1)
c     EXCLUSIVE FOR UPDATING (P[I]:I=1,...,L)
      real(8) pn
      real(8) p(nseries)
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
c      real(8) cw
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
c     PROBABILITIES (P[I]:I=1,...,L)
      do i=1,nseries
         p(i)=1.d0/dble(ndata)
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
            llo=dlog(p(i))-dlog(1.d0-p(i))
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
               sw=sw+rw
               j=en(t+1)
               call logml(k,ndata,y,nthetas,theta,en,j,rw)
               sw=sw+rw
               j=eo(t)
               call logml(k,ndata,y,nthetas,theta,eo,j,rw)
               sw=sw-rw
               llr=sw+llo
               uw=unifr(0.d0,1.d0)
               rw=dlog(uw)-dlog(1.d0-uw)
               if (llr.gt.rw) then
                  c(i,t)=1
               else
                  c(i,t)=0
               end if
            end do
         end do
c=======================================================================
c        UPDATING PROBABILITY PARAMETERS
c=======================================================================
         do i=1,nseries
            s=0
            do t=1,(ndata-1)
               s=s+c(i,t)
            end do
            rw=a0(i)+dble(s)
            sw=b0(i)+dble((ndata-1)-s)
            pn=betar(rw,sw)
            p(i)=pn
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
               end do
            end do
            do i=1,nseries
               mcmcp(stored,i)=p(i)
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
