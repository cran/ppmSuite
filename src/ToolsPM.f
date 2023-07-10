c=======================================================================
c=======================================================================
c     TOOLS: (P)ROBABILITY (M)ODELS
c
c     SUBROUTINE NAMES:
c     · LOGPR2YCF
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine logpr2ycf(ntclusts,tclusts,assocg,shpsa,shpsb,logpr)
c=======================================================================
c=======================================================================
c     BEGIN: LOGPR2YCF SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN TWO TEMPORAL CLUSTERINGS C[1],C[2] OF {1,...,N}, WHERE
c     C[L]=(C[1,L],...,C[N-1,L]), "LOGPR2YCF" RETURNS THE LOG
c     PROBABILITY LOG{PR(C[1],C[2]|G)}, UNDER THE FOLLOWING
c     HIERARCHICAL MODEL:
c
c     A) PR(C[1],C[2]|P[1],P[2])=
c        PRODUCT(BERN(C[T,L]|P[L]):T=1,...,N-1 AND L=1,2).
c
c        · BERN(·|P[L]): BERNOULLI PMF WITH PARAMETER P[L].
c
c     B) (P[1],P[2])|G FOLLOWS A CONTINUOUS DISTRIBUTION SUPPORTED
c        ON [0,1]^2, WITH PDF GIVEN BY
c
c        PDF(P[1],P[2]|G)=
c        W[G](P[1],P[2])*PRODUCT(BET(P[L]|ALPHA[L],BETA[L]):L=1,2),
c
c        · W[G](P[1],P[2])=
c          1+{G*PRODUCT(P[L]-{ALPHA[L]/(ALPHA[L]+BETA[L])}:L=1,2)}.
c
c        · G IN [-L[B],U[B]]: CONTROLS THE CORRELATION BETWEEN P[1]
c          AND P[2].
c
c        · BET(·|ALPHA[L],BETA[L]): BETA PDF WITH SHAPE PARAMETERS
c          ALPHA[L] AND BETA[L].
c
c        · SHAPES (A[L],B[L]:L=1,2) ARE FIXED.
c
c     C) G~UNIFORM(-L[B],U[B]).
c
c        · L[B]=F[B]*[MAX{ALPHA[1]*ALPHA[2],BETA[1]*BETA[2]}]^(-1).
c
c        · U[B]=F[B]*[MAX{ALPHA[1]*BETA[2],ALPHA[2]*BETA[1]}]^(-1).
c
c        · F[B]=(ALPHA[1]+BETA[1])*(ALPHA[2]+BETA[2]).
c
c     WITH THIS INFORMATION,
c
c     LOG{PR(C[1],C[2]|G)}=LOG(F[1])+LOG(F[2]),
c
c     WHERE
c
c     LOG(F[1])=LOG{1+[G*PRODUCT(F(C[L])-{A[L]/(A[L]+B[L])}:L=1,2)]}
c
c     AND
c
c     LOG(F[2])=SUM(LOG{PR(C[L])}:L=1,2).
c
c     0) F(C[L])=(A[L]+SIGMA[L])/(A[L]+(N-1)+B[L]).
c
c        · SIGMA[L]=SUM(C[T,L]:T=1,...,N-1).
c
c     1) LOG{PR(C[L])}=
c        LOG{BE(A[L]+SIGMA[L],B[L]+(N-1)-SIGMA[L])}-LOG{BE(A[L],B[L])}.
c
c        · BE(·,·): BETA FUNCTION.
c=======================================================================
c     INPUTS
c=======================================================================
c     ntclusts: LENGTH OF THE TEMPORAL AXIS (N)
      integer ntclusts
c     tclusts: TEMPORAL CLUSTERINGS (C[L]:L=1,2)
      integer tclusts(ntclusts-1,2)
c     assocg: ASSOCIATION PARAMETER (G)
      real(8) assocg
c     shpsa: SHAPE PARAMETERS (A[L]:L=1,2)
      real(8) shpsa(2)
c     shpsb: SHAPE PARAMETERS (B[L]:L=1,2)
      real(8) shpsb(2)
c=======================================================================
c     OUTPUTS
c=======================================================================
c     logpr: LOG{PR(C[1],C[2]|G)}
      real(8) logpr
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ll,tt
c     EXCLUSIVE FOR STORING SIGMA[L]
      real(8) sumc
c     EXCLUSIVE FOR STORING A[L]
      real(8) shpa
c     EXCLUSIVE FOR STORING B[L]
      real(8) shpb
c     EXCLUSIVE FOR STORING LOG{PR(C[1],C[2]|G)}
      real(8) logf1
      real(8) logf2
c     OTHERS
      real(8) aw
      real(8) bw
      real(8) pw
      real(8) rw
      real(8) sw
c=======================================================================
c     ALGORITHM
c=======================================================================
c     LOGF2=LOG(F[2])
c     · LOG(F[2])=SUM(LOG{PR(C[L])}:L=1,2)
c     · LOG{PR(C[L])}=
c       LOG{BE(A[L]+SIGMA[L],B[L]+(N-1)-SIGMA[L])}-LOG{BE(A[L],B[L])}
      logf2=0.d0
c     PW=PRODUCT(F(C[L])-{A[L]/(A[L]+B[L])}:L=1,2)]
c     · F(C[L])=(A[L]+SIGMA[L])/(A[L]+(N-1)+B[L])
      pw=1.d0
      do ll=1,2
         sumc=0.d0
         do tt=1,(ntclusts-1)
            sumc=sumc+dble(tclusts(tt,ll))
         end do
         shpa=shpsa(ll)
         shpb=shpsb(ll)
c        AW=A[L]+SIGMA[L]
         aw=shpa+sumc
c        BW=B[L]+(N-1)-SIGMA[L]
         bw=(shpb+dble(ntclusts-1))-sumc
c        RW=LOG{BE(A[L],B[L])}
         rw=(log_gamma(shpa)+log_gamma(shpb))-log_gamma(shpa+shpb)
c        SW=LOG{BE(A[L]+SIGMA[L],B[L]+(N-1)-SIGMA[L])}
         sw=(log_gamma(aw)+log_gamma(bw))-log_gamma(aw+bw)
         logf2=logf2+(sw-rw)
c        RW=A[L]/(A[L]+B[L])
         rw=shpa/(shpa+shpb)
c        SW=(A[L]+SIGMA[L])/(A[L]+(N-1)+B[L])
         sw=(shpa+sumc)/(shpa+dble(ntclusts-1)+shpb)
         pw=pw*(sw-rw)
      end do
c     LOGF1=LOG(F[1])
c     · LOG(F[1])=
c       LOG{1+[G*PRODUCT(F(C[L])-{A[L]/(A[L]+B[L])}:L=1,2)]}
      logf1=dlog(1.d0+(assocg*pw))
c     LOGPR=LOG{PR(C[1],C[2]|G)}
      logpr=logf1+logf2
c=======================================================================
      return
c     END: LOGPR2YCF SUBROUTINE
      end
c=======================================================================
c=======================================================================

