c=======================================================================
c=======================================================================
c     TOOLS: (M)ARGINAL (L)IKELIHOODS
c
c     SUBROUTINE NAMES:
c     · LOGNORNIG
c     · LOGPOIGAM
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine lognornig(nobs,obs,npars,pars,labels,indj,val)
c=======================================================================
c=======================================================================
c     BEGIN: LOGNORNIG SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN A TIME SERIE (Y[T]:T=1,...,N), A SET OF CLUSTER LABELS
c     (E[T]:T=1,...,N) WITH 1=E[1]<=E[2]<=···<=E[N], AND
c     J IN {1,...,E[N]}, "LOGNORNIG" RETURNS THE LOG MARGINAL
c     LIKELIHOOD ML(Y[T]:T IN S[J]|THETA) INDUCED BY THE FOLLOWING
c     MODEL:
c
c     A) Y[T]|M,S^2~NORMAL(M,S^2):T IN S[J]={I IN {1,...,N}:E[I]=J}.
c
c        · GIVEN M AND S^2, OBSERVATIONS ARE INDEPENDENT.
c
c     B) M|S^2~NORMAL(MU,{(S^2)/KAPPA}).
c
c        · MEAN MU AND PRECISION KAPPA ARE FIXED.
c
c     C) S^2~INVERSE-GAMMA(ALPHA,BETA).
c
c        · SHAPE ALPHA AND SCALE BETA ARE FIXED.
c
c     IN THIS CASE, THETA=(MU,KAPPA,ALPHA,BETA).
c=======================================================================
c     INPUTS
c=======================================================================
c     nobs: LENGTH OF THE TIME SERIES (Y[T]:T=1,...,N)
      integer nobs
c     obs: TIME SERIES (Y[T]:T=1,...,N)
      real(8) obs(nobs)
c     npars: MAXIMUM NUMBER OF ADMISSIBLE PARAMETERS
      integer npars
c     pars: PARAMETER THETA=(MU,KAPPA,ALPHA,BETA)
      real(8) pars(npars)
c     labels: CLUSTER LABELS (E[T]:T=1,...,N)
      integer labels(nobs)
c     indj: INDEX J
      integer indj
c=======================================================================
c     OUTPUTS
c=======================================================================
c     val: LOG{ML(Y[T]:T IN S[J]|THETA)}
      real(8) val
c=======================================================================
c     C++ FUNCTIONS
c=======================================================================
c     FUNCTIONS IN "TOOLSR2.C"
      real(8) normd
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ss
c     EXCLUSIVE FOR STORING D
      real(8) dd
c     EXCLUSIVE FOR STORING LOG{ML(Y[T]:T IN S[J]|THETA)}
      real(8) muo
      real(8) kappao
      real(8) alphao
      real(8) betao
      real(8) mun
      real(8) kappan
      real(8) alphan
      real(8) betan
      real(8) sy
      real(8) my
      real(8) s2y
c     OTHERS
      real(8) logf1w
      real(8) logf2w
      real(8) logf3w
      real(8) xw(2)
c=======================================================================
c     SETTING VALUES FOR SOME WORKING VARIABLES
c=======================================================================
c      if (npars.lt.4) then
c         print *,'Error in lognornig subroutine: check npars'
c         stop
c      end if
      xw(1)=0.d0
      xw(2)=1.d0
      muo=pars(1)
      kappao=pars(2)
      alphao=pars(3)
      betao=pars(4)
c=======================================================================
c     ALGORITHM
c=======================================================================
      dd=0.d0
      sy=0.d0
      logf1w=0.d0
      do ss=1,nobs
         if (labels(ss).eq.indj) then
            dd=dd+1.d0
            sy=sy+obs(ss)
            logf1w=logf1w+normd(obs(ss),xw(1),dsqrt(xw(2)),1)
         end if
      end do
      my=sy/dd
      s2y=0.d0
      do ss=1,nobs
         if (labels(ss).eq.indj) then
            s2y=s2y+((obs(ss)-my)**2)
         end if
      end do
      mun=((muo*kappao)+sy)/(kappao+dd)
      kappan=kappao+dd
      alphan=alphao+(0.5d0*dd)
      betan=betao+(((0.5d0*(dd*kappao))*((my-muo)**2))/(dd+kappao))
      betan=betan+(0.5d0*s2y)
      call lognigd(xw,muo,kappao,alphao,betao,logf2w)
      call lognigd(xw,mun,kappan,alphan,betan,logf3w)
      val=logf1w+(logf2w-logf3w)
c=======================================================================
      return
c     END: LOGNORNIG SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine logpoigam(nobs,obs,npars,pars,labels,indj,val)
c=======================================================================
c=======================================================================
c     BEGIN: LOGPOIGAM SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN A TIME SERIE (Y[T]:T=1,...,N), A SET OF CLUSTER LABELS
c     (E[T]:T=1,...,N) WITH 1=E[1]<=E[2]<=···<=E[N], AND
c     J IN {1,...,E[N]}, "LOGNORNIG" RETURNS THE LOG MARGINAL
c     LIKELIHOOD ML(Y[T]:T IN S[J]|THETA) INDUCED BY THE FOLLOWING
c     MODEL:
c
c     A) Y[T]|L~POISSON(L):T S[J]={I IN {1,...,N}:E[I]=J}.
c
c        · GIVEN L, OBSERVATIONS ARE INDEPENDENT.
c
c     B) L~GAMMA(ALPHA,BETA).
c
c        · SHAPE ALPHA AND RATE BETA ARE FIXED.
c
c     IN THIS CASE, THETA=(ALPHA,BETA).
c=======================================================================
c     INPUTS
c=======================================================================
c     nobs: LENGTH OF THE TIME SERIES (Y[T]:T=1,...,N)
      integer nobs
c     obs: TIME SERIES (Y[T]:T=1,...,N)
      real(8) obs(nobs)
c     npars: MAXIMUM NUMBER OF ADMISSIBLE PARAMETERS
      integer npars
c     pars: PARAMETER THETA=(MU,KAPPA,ALPHA,BETA)
      real(8) pars(npars)
c     labels: CLUSTER LABELS (E[T]:T=1,...,N)
      integer labels(nobs)
c     indj: INDEX J
      integer indj
c=======================================================================
c     OUTPUTS
c=======================================================================
c     val: LOG{ML(Y[T]:T IN S[J]|THETA)}
      real(8) val
c=======================================================================
c     C++ FUNCTIONS
c=======================================================================
c     FUNCTIONS IN "TOOLSR2.C"
      real(8) gammad
      real(8) poisd
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ss
c     EXCLUSIVE FOR STORING D
      real(8) dd
c     EXCLUSIVE FOR STORING LOG{ML(Y[T]:T IN S[J]|THETA)}
      real(8) alphao
      real(8) betao
      real(8) alphan
      real(8) betan
      real(8) sy
c     OTHERS
      real(8) logf1w
      real(8) logf2w
      real(8) logf3w
c=======================================================================
c     SETTING VALUES FOR SOME WORKING VARIABLES
c=======================================================================
c      if (npars.lt.2) then
c         print *,'Error in logpoigam subroutine: check npars'
c         stop
c      end if
      alphao=pars(1)
      betao=pars(2)
c=======================================================================
c     ALGORITHM
c=======================================================================
      dd=0.d0
      sy=0.d0
      logf1w=0.d0
      do ss=1,nobs
         if (labels(ss).eq.indj) then
            dd=dd+1.d0
            sy=sy+obs(ss)
            logf1w=logf1w+poisd(obs(ss),1.d0,1)
         end if
      end do
      alphan=alphao+sy
      betan=betao+dd
      logf2w=gammad(1.d0,alphao,1.d0/betao,1)
      logf3w=gammad(1.d0,alphan,1.d0/betan,1)
      val=logf1w+(logf2w-logf3w)
c=======================================================================
      return
c     END: LOGPOIGAM SUBROUTINE
      end
c=======================================================================
c=======================================================================

