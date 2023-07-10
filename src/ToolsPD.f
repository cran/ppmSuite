c=======================================================================
c=======================================================================
c     TOOLS: (P)ROBABILITY (D)ISTRIBUTIONS
c
c     SUBROUTINE NAMES:
c     · LOGMVTD
c     · LOGNIGD
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine logmvtd(dmn,x,nu,mu,invsigma,logdetsigma,logdmvt)
c=======================================================================
c=======================================================================
c     BEGIN: LOGMVTD SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     "LOGMVTD" RETURNS THE LOG DENSITY, EVALUATED AT X IN (0,1)^D,
c     OF A D-DIMENSIONAL T DISTRIBUTION WITH DEGREES OF FREEDOM NU,
c     LOCATION VECTOR MU AND SCALE MATRIX SIGMA.
c=======================================================================
c     INPUTS
c=======================================================================
c     dmn: DIMENSION OF THE MULTIVARIATE T
      integer dmn
c     x: VECTOR WHERE THE DENSITY WILL BE EVALUATED
      real(8) x(dmn)
c     nu: DEGREES OF FREEDOM NU
      real(8) nu
c     mu: LOCATION VECTOR MU
      real(8) mu(dmn)
c     invsigma: INVERSE OF SCALE MATRIX SIGMA
      real(8) invsigma(dmn,dmn)
c     logdetsigma: LOG DETERMINANT OF SCALE MATRIX SIGMA
      real(8) logdetsigma
c=======================================================================
c     OUTPUTS
c=======================================================================
c     logdmvt: LOG{T[D](X|NU,MU,SIGMA)}
      real(8) logdmvt
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     INDEXES
      integer ii
      integer jj
c     EXCLUSIVE FOR CONSTANT LOG(PI)
      real(8) logpi
c     EXCLUSIVE FOR CALCULATING LOG{T[D](X|NU,MU,SIGMA)}
      real(8) qf
c     OTHERS
      real(8) aux1
      real(8) aux2
c=======================================================================
c     ALGORITHM
c=======================================================================
      logpi=dlog(4.d0*datan(1.d0))
      qf=0.d0
c      print *, "logpi", logpi
c      print *, "dmn", dmn
c      print *, "logdetsigma", logdetsigma
c      print *, "nu", nu
c      print *, "mu", mu(1), mu(2)
c      print *, "x", x(1), x(2)
      do ii=1,dmn
         do jj=1,dmn
            aux1=((x(ii)-mu(ii))*invsigma(ii,jj))*(x(jj)-mu(jj))
            qf=qf+aux1
         end do
      end do
      logdmvt=(-0.5d0*(nu+dble(dmn)))*dlog(1.d0+(qf/nu))
      aux1=log_gamma(0.5d0*(nu+dble(dmn)))-log_gamma(0.5d0*nu)
      aux2=((-0.5d0*dble(dmn))*(dlog(nu)+logpi))+(-0.5d0*logdetsigma)
      logdmvt=logdmvt+(aux1+aux2)
c=======================================================================
      return
c     END: LOGMVTD SUBROUTINE
      end
c=======================================================================
c=======================================================================
c
c
c=======================================================================
c=======================================================================
      subroutine lognigd(x,mu,kappa,alpha,beta,logdnig)
c=======================================================================
c=======================================================================
c     BEGIN: LOGNIGD SUBROUTINE
      implicit none
c=======================================================================
c     DESCRIPTION
c=======================================================================
c     "LOGNIGD" RETURNS THE LOG DENSITY, EVALUATED AT X=(M,S^2), OF A
c     NORMAL-INVERSE-GAMMA DISTRIBUTION GIVEN BY
c     LOG{NORMAL(M|MU,S^2/KAPPA)}+LOG{INV-GAMMA(S^2|ALPHA,BETA)}.
c=======================================================================
c     INPUTS
c=======================================================================
c     x: VECTOR WHERE THE DENSITY WILL BE EVALUATED
      real(8) x(2)
c     mu: MEAN (NORMAL COMPONENT)
      real(8) mu
c     kappa: PRECISION (NORMAL COMPONENT)
      real(8) kappa
c     alpha: SHAPE (INVERSE-GAMMA COMPONENT)
      real(8) alpha
c     beta: SCALE (INVERSE-GAMMA COMPONENT)
      real(8) beta
c=======================================================================
c     OUTPUTS
c=======================================================================
c     logdnig: LOG{NORMAL(MU|M,SIGMA^2/K)}+LOG{INV-GAMMA(SIGMA^2|A,B)}
      real(8) logdnig
c=======================================================================
c     C++ FUNCTIONS
c=======================================================================
c     FUNCTIONS IN "TOOLSR2.C"
      real(8) gammad
      real(8) normd
c=======================================================================
c     WORKING VARIABLES
c=======================================================================
c     OTHERS
      real(8) aux1
      real(8) aux2
c=======================================================================
c     ALGORITHM
c=======================================================================
      aux1=normd(x(1),mu,dsqrt(x(2)/kappa),1)
      aux2=gammad(1.d0/x(2),alpha,1.d0/beta,1)-(2.d0*dlog(x(2)))
      logdnig=aux1+aux2
c=======================================================================
      return
c     END: LOGNIGD SUBROUTINE
      end
c=======================================================================
c=======================================================================
