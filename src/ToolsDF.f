c=======================================================================
c=======================================================================
c     TOOLS: (D)ATA (F)ACTORS
c
c     SUBROUTINE NAMES:
c     · LOGDFBINBET
c     · LOGDFNORIGA
c     · LOGDFNORNIG
c     · LOGDFNORNOR
c     · LOGDFPOIGAM
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine logdfbinbet(nobs,obs,indi1,indi2,ntrials,shpa,shpb,
     & logdf)
c=======================================================================
c=======================================================================
c     BEGIN: LOGDFBINBET SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN TWO INDEXES -1<I[1]<I[2]<N+1 AND A TIME SERIES
c     (Y[1],...,Y[N]), "LOGDFBINBET" RETURNS THE LOG DATA FACTOR
c     LOG{DF(Y[S]:S=I[1]+1,...,I[2])} INDUCED BY THE FOLLOWING MODEL:
c
c     A) Y[S]|P~BINOMIAL(M,P):S=I[1]+1,...,I[2].
c
c        · GIVEN P, OBSERVATIONS ARE INDEPENDENT.
c
c        · NUMBER OF TRIALS M IS FIXED.
c
c     B) P~BETA(A,B).
c
c        · SHAPES A AND B ARE FIXED.
c
c     IN THIS CASE,
c
c     LOG{DF(Y[S]:S=I[1]+1,...,I[2])}=LOG(F[1])+LOG(F[2])+LOG(F[3]),
c
c     WHERE
c
c     LOG(F[1])=LGAMMA(A+SIGMA)+LGAMMA(B+(D*M)-SIGMA)+LGAMMA(A+B),
c
c     LOG(F[2])=-{LGAMMA(A+(D*M)+B)+LGAMMA(A)+LGAMMA(B)},
c
c     AND
c
c     LOG(F[3])=SUM(LOG{BIN(M,Y[S])}:S=I[1]+1,...,I[2]).
c
c     0) LGAMMA(·): LOG-GAMMA FUNCTION.
c
c     1) SIGMA=SUM(Y[S]:S=I[1]+1,...,I[2]).
c
c     2) D=I[2]-I[1].
c
c     3) LOG{BIN(M,Y[S])}=LGAMMA(M+1)-LGAMMA(Y[S]+1)-LGAMMA(M-Y[S]+1).
c=======================================================================
c     INPUTS
c=======================================================================
c     nobs: LENGTH OF THE TEMPORAL AXIS (N)
      integer nobs
c     obs: TIME SERIES (Y[1],...,Y[N])
      real(8) obs(nobs)
c     indi1: INDEX (I[1])
      integer indi1
c     indi2: INDEX (I[2])
      integer indi2
c     ntrials: NUMBER OF TRIALS (M)
      real(8) ntrials
c     shpa: SHAPE PARAMETER (A)
      real(8) shpa
c     shpb: SHAPE PARAMETER (B)
      real(8) shpb
c=======================================================================
c     OUTPUTS
c=======================================================================
c     logdf: LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      real(8) logdf
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ii
c     EXCLUSIVE FOR STORING D
      integer dd
c     EXCLUSIVE FOR STORING LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      real(8) logf1
      real(8) logf2
c     OTHERS
      real(8) dw
      real(8) mw
      real(8) rw
      real(8) sw
      real(8) yw
      real(8) tt
c=======================================================================
c     SETTING VALUES FOR SOME WORKING VARIABLES
c=======================================================================
      dd=indi2-indi1
      dw=dble(dd)
      mw=ntrials
      rw=0.d0
c=======================================================================
c     ALGORITHM
c=======================================================================
c     SW=SIGMA
c     · SIGMA=SUM(Y[S]:S=I[1]+1,...,I[2])
      sw=0.d0
      do ii=1,dd
         sw=sw+obs(indi1+ii)
      end do
c     LOGF1=LOG(F[1])+LOG(F[2])
c     · LOG(F[1])=LGAMMA(A+SIGMA)+LGAMMA(B+(D*M)-SIGMA)+LGAMMA(A+B)
c     · LOG(F[2])=-{LGAMMA(A+(D*M)+B)+LGAMMA(A)+LGAMMA(B)}
      tt=shpb+(dw*mw)-sw
      rw=log_gamma(shpa+sw)+log_gamma(tt)+log_gamma(shpa+shpb)

      tt=shpa+(dw*mw)+shpb
      logf1=rw-(log_gamma(tt)+log_gamma(shpa)+log_gamma(shpb))

c     LOGF2=LOG(F[3])
c     · LOG(F[3])=SUM(LOG{BIN(M,Y[S])}:S=I[1]+1,...,I[2])
      logf2=0.d0
      do ii=1,dd
        yw=obs(indi1+ii)
        tt=log_gamma((mw-yw)+1.d0)
        rw=log_gamma(mw+1.d0)-(log_gamma(yw+1.d0)+tt)
        logf2=logf2+rw
      end do
c     LOGDF=LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      logdf=logf1+logf2
c=======================================================================
      return
c     END: LOGDFBINBET SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine logdfnoriga(nobs,obs,indi1,indi2,mmu,shpa,sclb,logdf)
c=======================================================================
c=======================================================================
c     BEGIN: LOGDFNORIGA SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN TWO INDEXES -1<I[1]<I[2]<N+1 AND A TIME SERIES
c     (Y[1],...,Y[N]), "LOGDFNORIGA" RETURNS THE LOG DATA FACTOR
c     LOG{DF(Y[S]:S=I[1]+1,...,I[2])} INDUCED BY THE FOLLOWING MODEL:
c
c     A) Y[S]|SIGMA^2~NORMAL(MU,SIGMA^2):S=I[1]+1,...,I[2].
c
c        · GIVEN SIGMA^2, OBSERVATIONS ARE INDEPENDENT.
c
c        · MEAN MU IS FIXED.
c
c     B) SIGMA^2~INVERSE-GAMMA(A,B).
c
c        · SHAPE A AND SCALE B ARE FIXED.
c
c     IN THIS CASE,
c
c     LOG{DF(Y[S]:S=I[1]+1,...,I[2])}=LOG(F[1])+LOG(F[2])+LOG(F[3]),
c
c     WHERE
c
c     LOG(F[1])=LGAMMA(A+(D/2))-LGAMMA(A),
c
c     LOG(F[2])=-[(D/2)*{LOG(PI)+LOG(2*B)}]
c
c     AND
c
c     LOG(F[3])={-A-(D/2)}*LOG{1+Q(X[D])}.
c
c     0) LGAMMA(·): LOG-GAMMA FUNCTION.
c
c     1) D=I[2]-I[1].
c
c     2) Q(X[D])=(X[D]-{MU*1[D]})'(A[D])(X[D]-{MU*1[D]}).
c
c        · X[D]=(Y[S]:S=I[1]+1,...,I[2])'.
c
c        · 1[D]: D-DIMENSIONAL VECTOR WITH ALL ENTRIES EQUAL 1.
c
c        · A[D]={1/(2*B)}*I[D].
c
c        · I[D]: D-DIMENSIONAL IDENTITY MATRIX.
c=======================================================================
c     INPUTS
c=======================================================================
c     nobs: LENGTH OF THE TEMPORAL AXIS (N)
      integer nobs
c     obs: TIME SERIES (Y[1],...,Y[N])
      real(8) obs(nobs)
c     indi1: INDEX (I[1])
      integer indi1
c     indi2: INDEX (I[2])
      integer indi2
c     mmu: MEAN PARAMETER (MU)
      real(8) mmu
c     shpa: SHAPE PARAMETER (A)
      real(8) shpa
c     sclb: SCALE PARAMETER (B)
      real(8) sclb
c=======================================================================
c     OUTPUTS
c=======================================================================
c     logdf: LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      real(8) logdf
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ii,jj
c     EXCLUSIVE FOR CONSTANT LOG(PI)
      real(8) logpi
c     EXCLUSIVE FOR STORING D
      integer dd
c     EXCLUSIVE FOR STORING LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      real(8) logf1
      real(8) logf2
c     OTHERS
      real(8) aw(indi2-indi1,indi2-indi1)
      real(8) dw
      real(8) rw
      real(8) sw
c=======================================================================
c     SETTING VALUES FOR SOME WORKING VARIABLES
c=======================================================================
      logpi=1.144729885849400174143427351353d0
      dd=indi2-indi1
c     AW=A[D]
c     · A[D]={1/(2*B)}*I[D]
      do ii=1,dd
         do jj=1,dd
            aw(ii,jj)=0.d0
         end do
         aw(ii,ii)=aw(ii,ii)+(0.5d0/sclb)
      end do
      dw=dble(dd)
c=======================================================================
c     ALGORITHM
c=======================================================================
c     LOGF1=LOG(F[1])+LOG(F[2])
c     · LOG(F[1])=LGAMMA(A+(D/2))-LGAMMA(A)
c     · LOG(F[2])=-[(D/2)*{LOG(PI)+LOG(2*B)}]
      rw=log_gamma(shpa+(0.5d0*dw))-log_gamma(shpa)
      logf1=rw-((0.5d0*dw)*(logpi+dlog(2.d0*sclb)))
c     LOGF2=LOG(F[3])
c     · {-A-(D/2)}*LOG{1+Q(X[D])}
c     · Q(X[D])=(X[D]-{MU*1[D]})'(A[D])(X[D]-{MU*1[D]})
      sw=0.d0
      do ii=1,dd
         do jj=1,dd
            rw=((obs(indi1+ii)-mmu)*aw(ii,jj))*(obs(indi1+jj)-mmu)
            sw=sw+rw
         end do
      end do
      logf2=(-shpa-(0.5d0*dw))*dlog(1.d0+sw)
c     LOGDF=LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      logdf=logf1+logf2
c=======================================================================
      return
c     END: LOGDFNORIGA SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine logdfnornig(nobs,obs,indi1,indi2,mmu0,pk0,shpa,sclb,
     & logdf)
c=======================================================================
c=======================================================================
c     BEGIN: LOGDFNORNIG SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN TWO INDEXES -1<I[1]<I[2]<N+1 AND A TIME SERIES
c     (Y[1],...,Y[N]), "LOGDFNORNIG" RETURNS THE LOG DATA FACTOR
c     LOG{DF(Y[S]:S=I[1]+1,...,I[2])} INDUCED BY THE FOLLOWING MODEL:
c
c     A) Y[S]|MU,SIGMA^2~NORMAL(MU,SIGMA^2):S=I[1]+1,...,I[2].
c
c        · GIVEN MU AND SIGMA^2, OBSERVATIONS ARE INDEPENDENT.
c
c     B) MU|SIGMA^2~NORMAL(MU[0],{(SIGMA^2)/K[0]}).
c
c        · MEAN MU[0] AND PRECISION K[0] ARE FIXED.
c
c     C) SIGMA^2~INVERSE-GAMMA(A,B).
c
c        · SHAPE A AND SCALE B ARE FIXED.
c
c     IN THIS CASE,
c
c     LOG{DF(Y[S]:S=I[1]+1,...,I[2])}=LOG(F[1])+LOG(F[2])+LOG(F[3]),
c
c     WHERE
c
c     LOG(F[1])=LGAMMA(A+(D/2))-LGAMMA(A),
c
c     LOG(F[2])=(1/2)*{LOG(KAPPA)-[D*{LOG(PI)-LOG(2*B)}]}
c
c     AND
c
c     LOG(F[3])={-A-(D/2)}*LOG{1+Q(X[D])}.
c
c     0) LGAMMA(·): LOG-GAMMA FUNCTION.
c
c     1) D=I[2]-I[1].
c
c     2) KAPPA=K[0]/(D+K[0]).
c
c     3) Q(X[D])=(X[D]-{MU[0]*1[D]})'(A[D])(X[D]-{MU[0]*1[D]}).
c
c        · X[D]=(Y[S]:S=I[1]+1,...,I[2])'.
c
c        · 1[D]: D-DIMENSIONAL VECTOR WITH ALL ENTRIES EQUAL 1.
c
c        · A[D]={1/(2*B)}*[I[D]+{(KAPPA-1)/D}*(1[D]1[D]')].
c
c        · I[D]: D-DIMENSIONAL IDENTITY MATRIX.
c=======================================================================
c     INPUTS
c=======================================================================
c     nobs: LENGTH OF THE TEMPORAL AXIS (N)
      integer nobs
c     obs: TIME SERIES (Y[1],...,Y[N])
      real(8) obs(nobs)
c     indi1: INDEX (I[1])
      integer indi1
c     indi2: INDEX (I[2])
      integer indi2
c     mmu0: MEAN PARAMETER (MU[0])
      real(8) mmu0
c     pk0: PRECISION PARAMETER (K[0])
      real(8) pk0
c     shpa: SHAPE PARAMETER (A)
      real(8) shpa
c     sclb: SCALE PARAMETER (B)
      real(8) sclb
c=======================================================================
c     OUTPUTS
c=======================================================================
c     logdf: LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      real(8) logdf
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ii,jj
c     EXCLUSIVE FOR CONSTANT LOG(PI)
      real(8) logpi
c     EXCLUSIVE FOR STORING D
      integer dd
c     EXCLUSIVE FOR STORING LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      real(8) kappa
      real(8) logf1
      real(8) logf2
c     OTHERS
      real(8) aw(indi2-indi1,indi2-indi1)
      real(8) dw
      real(8) rw
      real(8) sw
c=======================================================================
c     SETTING VALUES FOR SOME WORKING VARIABLES
c=======================================================================
      logpi=1.144729885849400174143427351353d0
      dd=indi2-indi1
      kappa=pk0/(dble(dd)+pk0)
c     AW=A[D]
c     · A[D]={1/(2*B)}*[I[D]+{(KAPPA-1)/D}*(1[D]1[D]')]
      do ii=1,dd
         do jj=1,dd
            aw(ii,jj)=(kappa-1.d0)/dble(dd)
         end do
         aw(ii,ii)=aw(ii,ii)+1.d0
      end do
      do ii=1,dd
         do jj=1,dd
            aw(ii,jj)=(0.5d0*aw(ii,jj))/sclb
         end do
      end do
      dw=dble(dd)
c=======================================================================
c     ALGORITHM
c=======================================================================
c     LOGF1=LOG(F[1])+LOG(F[2])
c     · LOG(F[1])=LGAMMA(A+(D/2))-LGAMMA(A)
c     · LOG(F[2])=(1/2)*{LOG(KAPPA)-[D*{LOG(PI)-LOG(2*B)}]}
      rw=log_gamma(shpa+(0.5d0*dw))-log_gamma(shpa)
      logf1=rw+(0.5d0*(dlog(kappa)-(dw*(logpi+dlog(2.d0*sclb)))))
c     LOGF2=LOG(F[3])
c     · LOG(F[3])={-A-(D/2)}*LOG{1+Q(X[D])}
c     · Q(X[D])=(X[D]-{MU[0]*1[D]})'(A[D])(X[D]-{MU[0]*1[D]})
      sw=0.d0
      do ii=1,dd
         do jj=1,dd
            rw=((obs(indi1+ii)-mmu0)*aw(ii,jj))*(obs(indi1+jj)-mmu0)
            sw=sw+rw
         end do
      end do
      logf2=(-shpa-(0.5d0*dw))*dlog(1.d0+sw)
c     LOGDF=LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      logdf=logf1+logf2
c=======================================================================
      return
c     END: LOGDFNORNIG SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine logdfnornor(nobs,obs,indi1,indi2,vsigma2,mmu0,
     & vsigma02,logdf)
c=======================================================================
c=======================================================================
c     BEGIN: LOGDFNORNOR SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN TWO INDEXES -1<I[1]<I[2]<N+1 AND A TIME SERIES
c     (Y[1],...,Y[N]), "LOGDFNORNOR" RETURNS THE LOG DATA FACTOR
c     LOG{DF(Y[S]:S=I[1]+1,...,I[2])} INDUCED BY THE FOLLOWING MODEL:
c
c     A) Y[S]|MU~NORMAL(MU,SIGMA^2):S=I[1]+1,...,I[2].
c
c        · GIVEN MU, OBSERVATIONS ARE INDEPENDENT.
c
c        · VARIANCE SIGMA^2 IS FIXED.
c
c     B) MU~NORMAL(MU[0],SIGMA[0]^2).
c
c        · MEAN MU[0] AND VARIANCE SIGMA[0]^2 ARE FIXED.
c
c     IN THIS CASE,
c
c     LOG{DF(Y[S]:S=I[1]+1,...,I[2])}=LOG(F[1])+LOG(F[2])+LOG(F[3]),
c
c     WHERE
c
c     LOG(F[1])=(-1/2)*{(D-1)*LOG(SIGMA^2)},
c
c     LOG(F[2])=(-1/2)*LOG[(SIGMA^2)+{D*(SIGMA[0]^2)}]
c
c     AND
c
c     LOG(F[3])=(-1/2)*[{D*LOG(2*PI)}+Q(X[D])].
c
c     0) D=I[2]-I[1].
c
c     1) Q(X[D])=(X[D]-{MU[0]*1[D]})'(A[D])(X[D]-{MU[0]*1[D]}).
c
c        · X[D]=(Y[S]:S=I[1]+1,...,I[2])'.
c
c        · 1[D]: D-DIMENSIONAL VECTOR WITH ALL ENTRIES EQUAL 1.
c
c        · A[D]={1/(SIGMA^2)}*[I[D]-{KAPPA*(1[D]1[D]')}].
c
c        · I[D]: D-DIMENSIONAL IDENTITY MATRIX.
c
c        · KAPPA=(SIGMA[0]^2)/[(SIGMA^2)+{D*(SIGMA[0]^2)}].
c=======================================================================
c     INPUTS
c=======================================================================
c     nobs: LENGTH OF THE TEMPORAL AXIS (N)
      integer nobs
c     obs: TIME SERIES (Y[1],...,Y[N])
      real(8) obs(nobs)
c     indi1: INDEX (I[1])
      integer indi1
c     indi2: INDEX (I[2])
      integer indi2
c     vsigma2: VARIANCE PARAMETER (SIGMA^2)
      real(8) vsigma2
c     mmu0: MEAN PARAMETER (MU[0])
      real(8) mmu0
c     vsigma02: VARIANCE PARAMETER (SIGMA[0]^2)
      real(8) vsigma02
c=======================================================================
c     OUTPUTS
c=======================================================================
c     logdf: LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      real(8) logdf
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ii,jj
c     EXCLUSIVE FOR CONSTANT LOG(2*PI)
      real(8) log2pi
c     EXCLUSIVE FOR STORING D
      integer dd
c     EXCLUSIVE FOR STORING LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      real(8) kappa
      real(8) logf1
      real(8) logf2
c     OTHERS
      real(8) aw(indi2-indi1,indi2-indi1)
      real(8) dw
      real(8) rw
      real(8) sw
c=======================================================================
c     SETTING VALUES FOR SOME WORKING VARIABLES
c=======================================================================
      log2pi=1.837877066409345483560659472811d0
      dd=indi2-indi1
      kappa=vsigma02/(vsigma2+(dble(dd)*vsigma02))
c     AW=A[D]
c     · A[D]={1/(SIGMA^2)}*[I[D]-{KAPPA*(1[D]1[D]')}]
      do ii=1,dd
         do jj=1,dd
            aw(ii,jj)=-kappa
         end do
         aw(ii,ii)=aw(ii,ii)+1.d0
      end do
      do ii=1,dd
         do jj=1,dd
            aw(ii,jj)=aw(ii,jj)/vsigma2
         end do
      end do
      dw=dble(dd)
c=======================================================================
c     ALGORITHM
c=======================================================================
c     LOGF1=LOG(F[1])+LOG(F[2])
c     · LOG(F[1])=(-1/2)*{(D-1)*LOG(SIGMA^2)}
c     · LOG(F[2])=(-1/2)*LOG[(SIGMA^2)+{D*(SIGMA[0]^2)}]
      rw=((dw-1.d0)*dlog(vsigma2))+dlog(vsigma2+(dw*vsigma02))
      logf1=-0.5d0*rw
c     LOGF2=LOG(F[3])
c     · LOG(F[3])=(-1/2)*[{D*LOG(2*PI)}+Q(X[D])]
c     · Q(X[D])=(X[D]-{MU[0]*1[D]})'(A[D])(X[D]-{MU[0]*1[D]})
      sw=0.d0
      do ii=1,dd
         do jj=1,dd
            rw=((obs(indi1+ii)-mmu0)*aw(ii,jj))*(obs(indi1+jj)-mmu0)
            sw=sw+rw
         end do
      end do
      logf2=-0.5d0*((dw*log2pi)+sw)
c     LOGDF=LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      logdf=logf1+logf2
c=======================================================================
      return
c     END: LOGDFNORNOR SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine logdfpoigam(nobs,obs,indi1,indi2,shpa,ratb,logdf)
c=======================================================================
c=======================================================================
c     BEGIN: LOGDFPOIGAM SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     GIVEN TWO INDEXES -1<I[1]<I[2]<N+1 AND A TIME SERIES
c     (Y[1],...,Y[N]), "LOGDFPOIGAM" RETURNS THE LOG DATA FACTOR
c     LOG{DF(Y[S]:S=I[1]+1,...,I[2])} INDUCED BY THE FOLLOWING MODEL:
c
c     A) Y[S]|LAMBDA~POISSON(LAMBDA):S=I[1]+1,...,I[2].
c
c        · GIVEN LAMBDA, OBSERVATIONS ARE INDEPENDENT.
c
c     B) LAMBDA~GAMMA(A,B).
c
c        · SHAPE A AND RATE B ARE FIXED.
c
c     IN THIS CASE,
c
c     LOG{DF(Y[S]:S=I[1]+1,...,I[2])}=LOG(F[1])+LOG(F[2])+LOG(F[3]),
c
c     WHERE
c
c     LOG(F[1])={A*LOG(B)}-LGAMMA(A),
c
c     LOG(F[2])=LGAMMA(A+SIGMA)-{(A+SIGMA)*LOG(B+D)}
c
c     AND
c
c     LOG(F[3])=SUM(-LOG(Y[S]!):S=I[1]+1,...,I[2]).
c
c     0) LGAMMA(·): LOG-GAMMA FUNCTION.
c
c     1) SIGMA=SUM(Y[S]:S=I[1]+1,...,I[2]).
c
c     2) D=I[2]-I[1].
c
c     3) LOG(Y[S]!)=LGAMMA(Y[S]+1).
c=======================================================================
c     INPUTS
c=======================================================================
c     nobs: LENGTH OF THE TEMPORAL AXIS (N)
      integer nobs
c     obs: TIME SERIES (Y[1],...,Y[N])
      real(8) obs(nobs)
c     indi1: INDEX (I[1])
      integer indi1
c     indi2: INDEX (I[2])
      integer indi2
c     shpa: SHAPE PARAMETER (A)
      real(8) shpa
c     ratb: RATE PARAMETER (B)
      real(8) ratb
c=======================================================================
c     OUTPUTS
c=======================================================================
c     logdf: LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      real(8) logdf
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ii
c     EXCLUSIVE FOR STORING D
      integer dd
c     EXCLUSIVE FOR STORING LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      real(8) logf1
      real(8) logf2
c     OTHERS
      real(8) dw
      real(8) rw
      real(8) sw
      real(8) yw
c=======================================================================
c     SETTING VALUES FOR SOME WORKING VARIABLES
c=======================================================================
      dd=indi2-indi1
      dw=dble(dd)
c=======================================================================
c     ALGORITHM
c=======================================================================
c     SW=SIGMA
c     · SIGMA=SUM(Y[S]:S=I[1]+1,...,I[2])
      sw=0.d0
      do ii=1,dd
         sw=sw+obs(indi1+ii)
      end do
c     LOGF1=LOG(F[1])+LOG(F[2])
c     · LOG(F[1])={A*LOG(B)}-LGAMMA(A)
c     · LOG(F[2])=LGAMMA(A+SIGMA)-{(A+SIGMA)*LOG(B+D)}
      rw=(shpa*dlog(ratb))-log_gamma(shpa)
      logf1=rw+(log_gamma(shpa+sw)-((shpa+sw)*dlog(ratb+dw)))
c     LOGF2=LOG(F[3])
c     ·LOG(F[3])=SUM(-LOG(Y[S]!):S=I[1]+1,...,I[2])
      logf2=0.d0
      do ii=1,dd
         yw=obs(indi1+ii)
         rw=log_gamma(yw+1.d0)
         logf2=logf2-rw
      end do
c     LOGDF=LOG{DF(Y[S]:S=I[1]+1,...,I[2])}
      logdf=logf1+logf2
c=======================================================================
      return
c     END: LOGDFPOIGAM SUBROUTINE
      end
c=======================================================================
c=======================================================================

