c=======================================================================
c=======================================================================
c     TOOLS: (G)IBBS (S)AMPLING
c
c     SUBROUTINE NAMES:
c     · LOGML
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine logml(imltype,inobs,iobs,inpars,ipars,ilabels,iindex,
     & oval)
c=======================================================================
c=======================================================================
c     BEGIN: LOGML SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN A TIME SERIE (Y[T]:T=1,...,N), A SET OF CLUSTER LABELS
c     (E[T]:T=1,...,N) WITH 1=E[1]<=E[2]<=···<=E[N], AND
c     J IN {1,...,E[N]}, "LOGML" RETURNS THE LOG MARGINAL LIKELIHOOD
c     ML(Y[T]:T IN S[J]|THETA) UNDER DIFFERENT STATISTICAL MODELS.
c     HERE, THETA IS A PARAMETER THAT CONTROLS ML(·) AND
c     S[J]={I IN {1,...,N}:E[I]=J}.
c=======================================================================
c     INPUTS
c=======================================================================
c     imltype: TYPE OF MARGINAL LIKELIHOOD TO BE EVALUATED
      integer imltype
c     inobs: LENGTH OF THE TIME SERIES (Y[T]:T=1,...,N)
      integer inobs
c     iobs: TIME SERIES (Y[T]:T=1,...,N)
      real(8) iobs(inobs)
c     inpars: MAXIMUM NUMBER OF ADMISSIBLE PARAMETERS
      integer inpars
c     ipars: PARAMETER THETA
      real(8) ipars(inpars)
c     ilabels: CLUSTER LABELS (E[T]:T=1,...,N)
      integer ilabels(inobs)
c     iindex: INDEX J
      integer iindex
c=======================================================================
c     OUTPUTS
c=======================================================================
c     oval: LOG{ML(Y[T]:T IN S[J]|THETA)}
      real(8) oval
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     OTHERS
      real(8) valw
c=======================================================================
c     ALGORITHM
c=======================================================================
      valw=0.d0
      if (imltype.eq.1) then
         call lognornig(inobs,iobs,inpars,ipars,ilabels,iindex,valw)
         oval=valw
      else if (imltype.eq.2) then
         call logpoigam(inobs,iobs,inpars,ipars,ilabels,iindex,valw)
         oval=valw
      else
c         print *,'Error in logml subroutine: check indtype'
c         stop
      end if
c=======================================================================
      return
c     END: LOGML SUBROUTINE
      end
c=======================================================================
c=======================================================================
