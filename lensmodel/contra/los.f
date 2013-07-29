c  Line-of-sight velocity dispersion of a tracer population
c  following Mamon & Lokas 2005, MNRAS, 363, 705 (eq. A16)
c
c  Written by Oleg Gnedin and Gary Mamon
c  Last modified: January 15, 2008

      subroutine losdispersion( ri, m, rhot, ANIS, ra, n, Rmax, los )
      implicit none
      integer ANIS, n, i, n_, ANIS_
      real*8 ri(n), m(n), rhot(n), los(n), ra, Rmax, R, sden1, los1,
     &     lnr_, lnm_, lnrhot_, ra_, surf_den, los_disp
      common lnr_(1000), lnm_(1000), lnrhot_(1000), n_, ANIS_, ra_
      external surf_den, los_disp, midpnt

      if(n .gt. 1000) then
         write(6,*) ': increase array size in los.f'
         stop
      else
         do i=1,n
            lnr_(i) = log(ri(i))
            lnm_(i) = log(m(i))
            lnrhot_(i) = log(rhot(i))
         enddo
         n_ = n
         ANIS_ = ANIS
         ra_ = ra
      endif

      ! integrate I(R) and sigma(R) from R to Rmax
      do i=1,n
         R = ri(i)
         los(i) = 0.0
         if(R .lt. Rmax) then
            call qromo(surf_den, R, Rmax, sden1, midpnt, R, 1.d-6)
            call qromo(los_disp, R, Rmax, los1, midpnt, R, 1.d-6)
            los(i) = dsqrt(los1/sden1)
         endif
      enddo
      return
      end


      function surf_den( x, R )
      implicit none
      integer n_, ANIS_
      real*8 surf_den, x, R, lnrho, lnr_, lnm_, lnrhot_, ra_
      common lnr_(1000), lnm_(1000), lnrhot_(1000), n_, ANIS_, ra_

      if(x .le. R) then
         surf_den = 0.0
      else
         call slinear(lnr_, lnrhot_, n_, log(x), lnrho)
         surf_den = 2.0*dexp(lnrho)*x/dsqrt(x*x-R*R)
      endif
      return
      end


      function los_disp( x, R )
      implicit none
      integer n_, ANIS_
      real*8 los_disp, x, R, kernel, lnrho, lnmr, lnr_, lnm_, lnrhot_,
     &     ra_, anis_kernel
      common lnr_(1000), lnm_(1000), lnrhot_(1000), n_, ANIS_, ra_
      external anis_kernel

      kernel = anis_kernel(x, R, ra_, ANIS_)
      call slinear(lnr_, lnrhot_, n_, log(x), lnrho)
      call slinear(lnr_, lnm_, n_, log(x), lnmr)
      los_disp = 2.0*kernel*dexp(lnrho+lnmr)/x
      return
      end


      function acosh( f )
      implicit none
      real*8 acosh, f
      acosh = dlog(f + dsqrt(f*f-1.d0))
      return
      end


      function anis_kernel( x, R, ra, ANIS )
      implicit none
      integer ANIS, k, n
      real*8 anis_kernel, x, R, ra, beta, u, ua,
     &     pi, f, c, s, y, z, tiny,
     &     gammln, betai, acosh, bico
      parameter (pi = 3.1415926535897932384626433832795028d0)
      external gammln, betai, acosh, bico

      anis_kernel = 0.d0
      u = x/R
      ua = ra/R
      y = 1.d0/(u*u)
      tiny = 1.d-10

c     Isotropic velocity distribution
      if(ANIS .eq. 0 .or. (ANIS .eq. 1 .and. dabs(ra) .lt. tiny)) then
         anis_kernel = dsqrt(1.d0-y)
         return
      endif

c     Radial velocity distribution
      if(ANIS .eq. 1 .and. dabs(ra-1.0) .lt. tiny) then
         anis_kernel = 0.25*pi*u - 0.5*dsqrt(1.0-y) - 0.5*u*asin(1.d0/u)
         return
      endif

c     Constant anisotropy parameter, beta=const
      if(ANIS .eq. 1) then
         beta = ra
c     if beta = 0.5 or -0.5
         if(dabs(beta-0.5) .lt. tiny .or. dabs(beta+0.5) .lt. tiny) then
            anis_kernel = u**(2.*beta-1.d0)*acosh(u)-beta*dsqrt(1.d0-y)
            return
         endif
c     if beta = -1, -2, ...
         if(beta .lt. -tiny .and. dabs(nint(beta)-beta) .lt. tiny) then
            n = nint(-beta)
            z = u*u-1.d0
            f = z**n/dfloat(2*n+1) + dfloat(n+1)
            do k=1,n-1
               f = f + z**k/dfloat(2*k+1)*dfloat(1+n-k)*bico(n,k)
            enddo
            anis_kernel = dsqrt(z)/u**(2*n+1)*f
            return
         endif

c         if(beta .gt. -0.5) then
c            g2 = sqrt(pi)*exp(gammln(beta+0.5)-gammln(beta+1.))
c            g1 = g2*beta/(beta-0.5)
c            g3 = g1
c            b2 = betai(beta+0.5, 0.5d0, 1./u**2)
c            b3 = b2 + (1./u**2)**(beta-0.5)*sqrt(1.-1./u**2)/beta/g2
c            f = (1.5-beta)*g1 + beta*b2*g2 - b3*g3
c            anis_kernel = 0.5d0*u**(2.*beta-1.)*f
c            if(f.ne.f) print *,'beta=',beta,' f=',f, ' u=',u
c            return
c         endif
c         print *,'no los dispersion for beta=',beta

c     for all other beta
         if(beta .lt. 1.0) then
            f = dsqrt(1.d0-y)/(1.d0-2.d0*beta)
     &           + 0.5d0*dsqrt(pi)
     &           *beta/(beta-0.5)*dexp(gammln(beta+0.5)-gammln(beta+1.))
     &           *(1.5d0-beta)*u**(2.d0*beta-1.d0)
     &           *(1.d0-betai(beta+0.5d0,0.5d0,y))
            anis_kernel = f
            if(f.ne.f) print *,'beta=',beta,' f=',f, ' u=',u
         else
            print *,'beta cannot be greater than 1: beta=',beta
         endif
      endif

c     Mamon-Lokas (2005) anisotropy model
      if(ANIS .eq. 2) then
         if(dabs(ua-1.0) .lt. tiny) then
            anis_kernel = (1.d0+1./u)*acosh(u) 
     &           - (8./u+7.)/6.*dsqrt((u-1.)/(u+1.))
         else
            if(ua .gt. 1.0) then
               c = acosh((ua*u+1.d0)/(u+ua))
               s = 1.d0
            else
               c = acos((ua*u+1.d0)/(u+ua))
               s = -1.d0
            endif
            f = ua*(ua**2-0.5)/(dabs(ua**2-1.))**1.5d0*(1.+ua/u)*c
            anis_kernel = 0.5d0/(ua**2-1.)*dsqrt(1.-1./u**2) 
     &           + (1.+ua/u)*acosh(u) - s*f
         endif
      endif

c     Osipkov-Merritt anisotropy model
      if(ANIS .eq. 3) then
         f = atan(dsqrt((u*u-1.d0)/(ua*ua+1.d0)))
         anis_kernel = (ua**2+0.5)/(ua**2+1.)**1.5d0*(u+ua**2/u)*f
     &        -0.5d0/(ua**2+1.)*dsqrt(1.-1./u**2)
      endif
      return
      end


      subroutine slinear( x_, y_, n, x, s )
      implicit none
      integer n, i, di
      real*8 x, x_(n), y_(n), s

      if(x_(n).gt.x_(1)) then
         di=1
         do i=2,n-1,di
            if(x .le. x_(i)) goto 999
         enddo
      else
         di=-1
         do i=n,3,di
            if(x .le. x_(i)) goto 999
         enddo
      endif
 999  continue
      s=y_(i)+(x-x_(i))*(y_(i)-y_(i-di))/(x_(i)-x_(i-di))
      return
      end


C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
      subroutine polint(xa,ya,n,x,y,dy)
      implicit none
      integer n,NMAX
      real*8 dy,x,y,xa(n),ya(n)
      parameter (NMAX=10)
      integer i,m,ns
      real*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
      subroutine qromo(func,a,b,ss,choose,par1,EPS)
      implicit none
      integer JMAX,JMAXP,K,KM
      real*8 func,a,b,ss,EPS,par1
      external func,choose
      parameter (JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
      integer j
      real*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.d0
      do 11 j=1,JMAX
        call choose(func,a,b,s(j),j,par1)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if(j.ge.19) WRITE(0,'(A,I3,A,1pe13.6)') 'j=',j,'  s=',ss
          if (dabs(dss).le.EPS*dabs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.d0
11    continue
      pause 'too many steps in qromo'
      end
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
      subroutine midpnt(func,a,b,s,n,par1)
      implicit none
      integer n
      real*8 a,b,s,func,par1
      external func
      integer it,j
      real*8 ddel,del,sum,tnm,x
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b),par1)
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x,par1)
          x=x+ddel
          sum=sum+func(x,par1)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
      function betai(a,b,x)
      real*8 betai,a,b,x
CU    USES betacf,gammln
      real*8 bt,betacf,gammln
      if(x.lt.0..or.x.gt.1.)pause 'bad argument x in betai'
      if(x.eq.0..or.x.eq.1.)then
        bt=0.
      else
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
      endif
      if(x.lt.(a+1.)/(a+b+2.))then
        betai=bt*betacf(a,b,x)/a
        return
      else
        betai=1.-bt*betacf(b,a,1.-x)/b
        return
      endif
      end
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
      function betacf(a,b,x)
      integer MAXIT
      real*8 betacf,a,b,x,EPS,FPMIN
      parameter (MAXIT=100,EPS=3.e-7,FPMIN=1.e-30)
      integer m,m2
      real*8 aa,c,d,del,h,qab,qam,qap
      qab=a+b
      qap=a+1.
      qam=a-1.
      c=1.
      d=1.-qab*x/qap
      if(abs(d).lt.FPMIN)d=FPMIN
      d=1./d
      h=d
      do 11 m=1,MAXIT
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        h=h*d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a or b too big, or MAXIT too small in betacf'
1     betacf=h
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
      function gammln(xx)
      real*8 gammln,xx
      integer j
      real*8 ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
      function gammp(a,x)
      real*8 a,gammp,x
CU    USES gcf,gser
      real*8 gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
      subroutine gcf(gammcf,a,x,gln)
      integer ITMAX
      real*8 a,gammcf,gln,x,EPS,FPMIN
      parameter (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammln
      integer i
      real*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
      subroutine gser(gamser,a,x,gln)
      integer ITMAX
      real*8 a,gamser,gln,x,EPS
      parameter (ITMAX=100,EPS=3.e-7)
CU    USES gammln
      integer n
      real*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
      function bico(n,k)
      integer k,n
      real*8 bico
CU    USES factln
      real*8 factln
      bico=nint(exp(factln(n)-factln(k)-factln(n-k)))
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
      function factln(n)
      integer n
      real*8 factln
CU    USES gammln
      real*8 a(100),gammln
      save a
      data a/100*-1./
      if (n.lt.0) pause 'negative factorial in factln'
      if (n.le.99) then
        if (a(n+1).lt.0.) a(n+1)=gammln(n+1.d0)
        factln=a(n+1)
      else
        factln=gammln(n+1.d0)
      endif
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
