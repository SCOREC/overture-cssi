      subroutine mgeosi (j,k)
c
c Initialize EOS for phase j and for material k
c
      implicit real*8 (a-h,o-z)
      common / matdatmg / kmat(2)
      common / romdatmg / aint0rom(2,2),r0rom(2,2),c0rom(2,2),
     *                    tolrom(2),lvrom(2)
      include 'PBX9501reactant-scaled.h'
      include 'HMXproducts-mod-scaled.h'
      include 'IdealGas.h'
      include 'StiffenedGas.h'
      include 'VirialGas.h'
      include 'IdealGas1p35.h'
c
      if (k.eq.1) then
        write(6,101)
  101   format('* PBX-9501 reactants (scaled)',/,
     *         '  rhoRef=.37778d0 g/cm^3',/,
     *         '  pRef  = 21.3d0 GPa')
        irhob=0
        nmgb=10000
        xminb=1.d-1
        xmaxb=1.d0
        irhoh=0
        nmgh=10000
        xminh=1.d-1
        xmaxh=1.d0
        nmgt=10000
        xmaxt=100.d0
c
      elseif (k.eq.2) then
        write(6,102)
  102   format('* HMX products (modified, scaled)',/,
     *         '  rhoRef=.37778d0 g/cm^3',/,
     *         '  pRef  = 21.3d0 GPa')
        irhob=1
        nmgb=10000
        xminb=1.d-5
        xmaxb=1.d1
        irhoh=1
        nmgh=10000
        xminh=1.d-5
        xmaxh=1.d1
        nmgt=2
        xmaxt=1.d0
c
      elseif (k.eq.3) then
        write(6,103)gam3
  103   format('* Ideal gas: gamma=',f7.4,/,
     *         '  (No look-up tables required)')
        irhob=0
        nmgb=2
        xminb=0.d0
        xmaxb=10.d0
        irhoh=0
        nmgh=2
        xminh=0.01d0
        xmaxh=10.d0
        nmgt=2
        xmaxt=100.d0
c
      elseif (k.eq.4) then
        write(6,104)gam4,pi4
  104   format('* Stiffened gas: gamma=',f7.4,/,
     *         '                    pi=',f7.4,/,
     *         '  (No look-up tables required)')
        irhob=0
        nmgb=2
        xminb=0.d0
        xmaxb=10.d0
        irhoh=0
        nmgh=2
        xminh=0.01d0
        xmaxh=10.d0
        nmgt=2
        xmaxt=100.d0
c
      elseif (k.eq.5) then
        write(6,105)gam5,b5
  105   format('* Virial gas: gamma=',f7.4,/,
     *         '                  b=',f7.4,/,
     *         '  (No look-up tables required)')
        irhob=1
        nmgb=2
        xminb=1.d-4
        xmaxb=10.d0
        irhoh=1
        nmgh=2
        xminh=1.d-4
        xmaxh=10.d0
        nmgt=2
        xmaxt=100.d0
c
      elseif (k.eq.6) then
        write(6,106)gam6
  106   format('* Ideal gas: gamma=',f7.4,/,
     *         '  (No look-up tables required)')
        irhob=0
        nmgb=2
        xminb=0.d0
        xmaxb=10.d0
        irhoh=0
        nmgh=2
        xminh=0.01d0
        xmaxh=10.d0
        nmgt=2
        xmaxt=100.d0
c
      else
        write(6,*)'Error (mginfo) : k not supported'
        stop
      end if
c
c set kmat flag for phase j
      kmat(j)=k
c
c tolerances for Romberg integration of MG isentrope
      tolrom(j)=1.d-12
      lvrom(j)=15
c
c set up look-up tables
      call getBmg (j,kmat(j),irhob,nmgb,xminb,xmaxb)
      call getHmg (j,kmat(j),irhoh,nmgh,xminh,xmaxh)
      call getTmg (j,kmat(j),nmgt,xmaxt)
c
      return
      end
c
c++++++++++++++=
c
      double precision function c2mg (k,V,P)
c
c compute the square of the sound speed for a given specific volume V
c and pressure P, for material k
c
      implicit real*8 (a-h,o-z)
      include 'PBX9501reactant-scaled.h'
      include 'HMXproducts-mod-scaled.h'
      include 'IdealGas.h'
      include 'StiffenedGas.h'
      include 'VirialGas.h'
      include 'IdealGas1p35.h'
c
      if (k.eq.1) then
c
c PBX-9501 reactant
        z=V/V0
        A=V/(a1s+a2s*z)
        Ap=a1s/(a1s+a2s*z)**2
        x=f*(z**e1-z**e2)
        y=1.d0+g*(z**e3-1.d0)
        xp=f*(e1*z**e1m1-e2*z**e2m1)/V0
        yp=g*(e3*z**e3m1)/V0
        P0=x*y
        P0p=x*yp+xp*y
        Bp=-Ap*P0-A*P0p-P0
        c2mg=(((Ap+1.d0)*P+Bp)*V**2)/A
c
      elseif (k.eq.2) then
c
c HMX products
        rho=1.d0/V
        Gam=gam0
        Gamp=0.d0                               ! deriv wrt rho
        if (rho.gt.rhoc) then
          z=rho-rhoc
          Gam=Gam+(d2g+d3g*z)*z**2
          Gamp=(d2gp+d3gp*z)*z
        end if
        Pref=Ag*dexp(akg*(rho-rhom))*rho**qq
        Prefp=(qq*V+akg)*Pref                   ! deriv wrt rho
        A=V/Gam
        Ap=(1.d0+rho*Gamp/Gam)/Gam              ! deriv wrt V
        Bp=-(Ap+1.d0)*Pref+A*Prefp*rho**2       ! deriv wrt V
        c2mg=(((Ap+1.d0)*P+Bp)*V**2)/A

      elseif (k.eq.3) then
c
c Ideal gas
        c2mg=gam3*V*P
c
      elseif (k.eq.4) then
c
c Stiffened gas
        c2mg=gam4*V*(P+pi4)
c
      elseif (k.eq.5) then
c
c Virial gas
        c2mg=(gam5*(V+b5)-b5**2/(V+b5))*P
c
      elseif (k.eq.6) then
c
c Ideal gas
        c2mg=gam6*V*P
c
      else
c
        write(6,*)'Error (c2mg) : k not supported, k=',k
        stop
c
      end if
c
      return
      end
c
c+++++++++++++++++++
c
      double precision function AmgV (k,V)
c
c compute A(V) for a given specific volume V, for material k
c
      implicit real*8 (a-h,o-z)
      include 'PBX9501reactant-scaled.h'
      include 'HMXproducts-mod-scaled.h'
      include 'IdealGas.h'
      include 'StiffenedGas.h'
      include 'VirialGas.h'
      include 'IdealGas1p35.h'
c
      if (k.eq.1) then
c
c PBX-9501 reactant
        z=V/V0
        AmgV=V/(a1s+a2s*z)
c
      elseif (k.eq.2) then
c
c HMX products
        rho=1.d0/V
        Gam=gam0
        if (rho.gt.rhoc) then
          z=rho-rhoc
          Gam=Gam+(d2g+d3g*z)*z**2
        end if
        AmgV=V/Gam
c
      elseif (k.eq.3) then
c
c Ideal gas
        AmgV=V/gm13
c
      elseif (k.eq.4) then
c
c Stiffened gas
        AmgV=V/gm14
c
      elseif (k.eq.5) then
c
c Virial gas
        AmgV=(V**2)/(gm15*(V+b5))
c
      elseif (k.eq.6) then
c
c Ideal gas
        AmgV=V/gm16
c
      else
c
        write(6,*)'Error (AmgV) : k not supported, k=',k
        stop
c
      end if
c
      return
      end
c
c+++++++++++++++++++
c
      double precision function amumg (k,V)
c
c compute mu(V)=A(V)*exp(I), I=int(1/A(V)) for material k
c at a given specific volume V.
c
c Note: for k=2, I=int(1/A(x), x=1/rhoc..V))
c
      implicit real*8 (a-h,o-z)
      include 'PBX9501reactant-scaled.h'
      include 'HMXproducts-mod-scaled.h'
      include 'IdealGas.h'
      include 'StiffenedGas.h'
      include 'VirialGas.h'
      include 'IdealGas1p35.h'
c
      if (k.eq.1) then
c
c PBX-9501 reactant
        z=V/V0
        A=V/(a1s+a2s*z)
        amumg=A*dexp(a2s*z)*z**a1s
c
      elseif (k.eq.2) then
c
c HMX products
        rho=1.d0/V
        if (rho.gt.rhoc) then
          z=rho-rhoc
          Gam=gam0+(d2g+d3g*z)*z**2
          arg=(e0g+rho*(e1g+rho*e2g))*z
          amumg=(V/Gam)*((rhoc/rho)**elog)*dexp(arg)
        else
          amumg=(V/gam0)*(rhoc/rho)**gam0
        end if
c
      elseif (k.eq.3) then
c
c Ideal gas
        amumg=(V**gam3)/gm13
c
      elseif (k.eq.4) then
c
c Stiffened gas
        amumg=(V**gam4)/gm14
c
      elseif (k.eq.5) then
c
c Virial gas
        amumg=(V**gp15)*dexp(-gm15*b5/V)/(gm15*(V+b5))
c
      elseif (k.eq.6) then
c
c Ideal gas
        amumg=(V**gam6)/gm16
c
      else
c
        write(6,*)'Error (amumg) : k not supported, k=',k
        stop
c
      end if
c
      return
      end
c
c+++++++++++++++++++
c
      subroutine getBmg (j,k,irho1,nmg1,xmin1,xmax1)
c
c compute the integral of Bprime(x) for x=xmin1:xn,
c where xn are nmg1 equally spaced nodes on xmin1:xmax1.
c
c If irho1.eq.0 , then x=specific volume
c         .ne.0 , then x=density
c
      implicit real*8 (a-h,o-z)
      parameter (nmgmax=200000)
      common / mgEOSb / irho(2),nmg(2),xmin(2),xmax(2),
     *                  Bmg(0:3,0:200000,2)
c
c certain EOS cases do not require a lookup table
      if (k.ge.3.and.k.le.6) return
c
      if (nmg1.gt.nmgmax) then
        write(6,*)'Error (getBmg) : nmg1 too big'
        stop
      end if
c
      irho(j)=irho1
c      write(6,*)'irho1=',irho1
      nmg(j)=nmg1
      xmin(j)=xmin1
      xmax(j)=xmax1
c
c get Bprime at the nodes
      dx=(xmax1-xmin1)/nmg1
      do n=0,nmg1
        xn=xmin1+n*dx
        Bmg(1,n,j)=Bprime(j,k,xn)
      end do
c
c get B at the nodes using Simpson's rule, and also compute
c divided differences for later Hermite interpolation
      Bmg(0,0,j)=0.d0
      do n=1,nmg1
        xm=xmin1+(n-.5d0)*dx
        Bpm=Bprime(j,k,xm)
        Bmg(0,n,j)=Bmg(0,n-1,j)
     *       +dx*(Bmg(1,n-1,j)+4.d0*Bpm+Bmg(1,n,j))/6.d0
        
        d1=(Bmg(0,n,j)-Bmg(0,n-1,j))/dx
        Bmg(2,n-1,j)=(d1-Bmg(1,n-1,j))/dx
        d2=(Bmg(1,n,j)-d1)/dx
        Bmg(3,n-1,j)=(d2-Bmg(2,n-1,j))/dx
      end do
c
      return
      end
c
c+++++++++++++++++++
c
      double precision function BmgV (j,V,ier)
c
c compute Bmg at specific volume V, for phase j material (specified in getBmg)
c
      implicit real*8 (a-h,o-z)
      common / mgEOSb / irho(2),nmg(2),xmin(2),xmax(2),
     *                  Bmg(0:3,0:200000,2)
      common / matdatmg / kmat(2)
      include 'StiffenedGas.h'
      data tiny / 1.d-8 /
c
      ier=0
c
c certain EOS cases can compute BmgV without a lookup table
      k=kmat(j)
      if (k.eq.3.or.k.eq.5.or.k.eq.6) then
        BmgV=0.d0
        return
      elseif (k.eq.4) then
        BmgV=gam4*pi4*V/gm14
        return
      end if
c
      x=V
      if (irho(j).ne.0) x=1.d0/x
c
c Check bounds
      if (x.le.xmin(j)) then
        ier=1
        write(6,*)'Error (BmgV) : x<xmin'
        return
      end if
c
      if (x.ge.xmax(j)) then
        ier=2
        write(6,*)'Error (BmgV) : x>xmax'
        return
      end if
c
c determine left index
      dx=(xmax(j)-xmin(j))/nmg(j)
      n=(x-xmin(j)-tiny)/dx
      n=max(min(n,nmg(j)-1),0)
      x1=xmin(j)+n*dx
c
c compute BmgV using Hermite interpolation
      z=x-x1
      BmgV=Bmg(0,n,j)+z*(Bmg(1,n,j)+z*(Bmg(2,n,j)+(z-dx)*Bmg(3,n,j)))
c
c      if (j.eq.1) then
c        write(6,*)V,BmgV
c        read(5,*)idum
c      end if
c
      return
      end
c
c+++++++++++++++++++
c
      double precision function Bprime (j,k,x)
c
c compute Bprime for phase j and material k.
c
c If irho(j).eq.0 , then x=specific volume
c           .ne.0 , then x=density
c
      implicit real*8 (a-h,o-z)
      common / mgEOSb / irho(2),nmg(2),xmin(2),xmax(2),
     *                  Bmg(0:3,0:200000,2)
      include 'PBX9501reactant-scaled.h'
      include 'HMXproducts-mod-scaled.h'
      include 'IdealGas.h'
      include 'StiffenedGas.h'
      include 'VirialGas.h'
      include 'IdealGas1p35.h'
c
      if (k.eq.1) then
c
c PBX-9501 reactant
        V=x
        if (irho(j).ne.0) then
          write(6,*)'Error (Bprime), irho.ne.0, k=',k
          stop
        end if
        z=V/V0
        A=V/(a1s+a2s*z)
        Ap=a1s/(a1s+a2s*z)**2
        xs=f*(z**e1-z**e2) 
        ys=1.d0+g*(z**e3-1.d0)
        xsp=f*(e1*z**e1m1-e2*z**e2m1)/V0
        ysp=g*(e3*z**e3m1)/V0
        P0=xs*ys
        P0p=xs*ysp+xsp*ys
        Bprime=-Ap*P0-A*P0p-P0
c
      elseif (k.eq.2) then
c
c HMX products
        rho=x
        if (irho(j).eq.0) then
          write(6,*)'Error (Bprime), irho.eq.0, k=',k
          stop
        end if
        Gam=gam0
        Gamp=0.d0                               ! deriv wrt rho
        if (rho.gt.rhoc) then
          z=rho-rhoc
          Gam=Gam+(d2g+d3g*z)*z**2
          Gamp=(d2gp+d3gp*z)*z
        end if
        Pref=Ag*dexp(akg*(rho-rhom))*rho**qq
        Prefp=(qq/rho+akg)*Pref                 ! deriv wrt rho
        A=1.d0/(rho*Gam)
        Ap=(1.d0+rho*Gamp/Gam)/Gam              ! deriv wrt V
        Bprime=-(Ap+1.d0)*Pref+A*Prefp*rho**2   ! deriv wrt V
        Bprime=-Bprime/rho**2                   ! deriv wrt rho
c
      elseif (k.eq.3.or.k.eq.6) then
c
c Ideal gas
        Bprime=0.d0
c
      elseif (k.eq.4) then
c
c Stiffened gas
        Bprime=gam4*pi4/gm14
c
      elseif (k.eq.5) then
c
c Virial gas
        Bprime=0.d0
c
      else
c
        write(6,*)'Error (Bprime) : k not supported, k=',k
        stop
c
      end if
c
      return
      end
c
c+++++++++++++++++++
c
      subroutine getHmg (j,k,irho1,nmg1,xmin1,xmax1)
c
c compute the integral of Hprime(x)=mu*Bprime/A for x=xmin1:xn,
c where xn are nmg1 equally spaced nodes on xmin1:xmax1.
c
c If irho1.eq.0 , then x=specific volume
c         .ne.0 , then x=density
c
      implicit real*8 (a-h,o-z)
      parameter (nmgmax=200000)
      common / mgEOSh / irho(2),nmg(2),xmin(2),xmax(2),
     *                  Hmg(0:3,0:200000,2)
c
c certain EOS cases do not require a lookup table
      if (k.ge.3.and.k.le.6) return
c
      if (nmg1.gt.nmgmax) then
        write(6,*)'Error (getHmg) : nmg1 too big'
        stop
      end if
c
      irho(j)=irho1
      nmg(j)=nmg1
      xmin(j)=xmin1
      xmax(j)=xmax1
c
c get Hprime at the nodes
      dx=(xmax1-xmin1)/nmg1
      do n=0,nmg1
        xn=xmin1+n*dx
        Hmg(1,n,j)=Hprime(j,k,xn)
      end do
c
c get H at the nodes using Simpson's rule, and also compute
c divided differences for later Hermite interpolation
      Hmg(0,0,j)=0.d0
      do n=1,nmg1
        xm=xmin1+(n-.5d0)*dx
        Hpm=Hprime(j,k,xm)
        Hmg(0,n,j)=Hmg(0,n-1,j)
     *               +dx*(Hmg(1,n-1,j)+4.d0*Hpm+Hmg(1,n,j))/6.d0
        d1=(Hmg(0,n,j)-Hmg(0,n-1,j))/dx
        Hmg(2,n-1,j)=(d1-Hmg(1,n-1,j))/dx
        d2=(Hmg(1,n,j)-d1)/dx
        Hmg(3,n-1,j)=(d2-Hmg(2,n-1,j))/dx
      end do
c
      return
      end
c
c+++++++++++++++++++
c
      double precision function HmgV (j,V,ier)
c
c compute Hmg at specific volume V, for phase j material (specified in getHmg)
c
      implicit real*8 (a-h,o-z)
      common / mgEOSh / irho(2),nmg(2),xmin(2),xmax(2),
     *                  Hmg(0:3,0:200000,2)
      common / matdatmg / kmat(2)
      include 'StiffenedGas.h'
      data tiny / 1.d-8 /
c
      ier=0
c
c certain EOS cases can compute HmgV without a lookup table
      k=kmat(j)
      if (k.eq.3.or.k.eq.5.or.k.eq.6) then
        HmgV=0.d0
        return
      elseif (k.eq.4) then
        HmgV=pi4*(V**gam4)/gm14
        return
      end if
c
      x=V
      if (irho(j).ne.0) x=1.d0/x
c
c Check bounds
      if (x.le.xmin(j)) then
        ier=1
        write(6,*)'Error (HmgV) : x<xmin'
        return
      end if
c
      if (x.ge.xmax(j)) then
        ier=2
        write(6,*)'Error (HmgV) : x>xmax'
        return
      end if
c
c determine left index
      dx=(xmax(j)-xmin(j))/nmg(j)
      n=(x-xmin(j)-tiny)/dx
      n=max(min(n,nmg(j)-1),0)
      x1=xmin(j)+n*dx
c
c compute HmgV using Hermite interpolation
      z=x-x1
      HmgV=Hmg(0,n,j)+z*(Hmg(1,n,j)+z*(Hmg(2,n,j)+(z-dx)*Hmg(3,n,j)))
c
      return
      end
c
c+++++++++++++++++++
c
      double precision function Hprime (j,k,x)
c
c compute Hprime for phase j and material k.
c
c If irho(j).eq.0 , then x=specific volume
c           .ne.0 , then x=density
c
c
      implicit real*8 (a-h,o-z)
      common / mgEOSh / irho(2),nmg(2),xmin(2),xmax(2),
     *                  Hmg(0:3,0:200000,2)
      include 'PBX9501reactant-scaled.h'
      include 'HMXproducts-mod-scaled.h'
      include 'IdealGas.h'
      include 'StiffenedGas.h'
      include 'VirialGas.h'
      include 'IdealGas1p35.h'
c
      if (k.eq.1) then
c
c PBX-9501 reactant
        V=x
        if (irho(j).ne.0) then
          write(6,*)'Error (Hprime), irho.ne.0, k=',k
          stop
        end if
        z=V/V0
        A=V/(a1s+a2s*z)
        Ap=a1s/(a1s+a2s*z)**2
        xs=f*(z**e1-z**e2)
        ys=1.d0+g*(z**e3-1.d0)
        xsp=f*(e1*z**e1m1-e2*z**e2m1)/V0
        ysp=g*(e3*z**e3m1)/V0
        P0=xs*ys
        P0p=xs*ysp+xsp*ys
        Bp=-Ap*P0-A*P0p-P0
        Hprime=Bp*dexp(a2s*z)*z**a1s
c
      elseif (k.eq.2) then
c
c HMX products
        rho=x
        V=1.d0/rho
        if (irho(j).eq.0) then
          write(6,*)'Error (Hprime), irho.eq.0, k=',k
          stop
        end if
        Gam=gam0
        Gamp=0.d0                               ! deriv wrt rho
        if (rho.gt.rhoc) then
          z=rho-rhoc
          Gam=gam0+(d2g+d3g*z)*z**2
          Gamp=(d2gp+d3gp*z)*z
          arg=(e0g+rho*(e1g+rho*e2g))*z
          amumg=(V/Gam)*((rhoc/rho)**elog)*dexp(arg)
        else
          amumg=(V/gam0)*(rhoc/rho)**gam0
        end if
        Pref=Ag*dexp(akg*(rho-rhom))*rho**qq
        Prefp=(qq*V+akg)*Pref                   ! deriv wrt rho
        A=V/Gam
        Ap=(1.d0+rho*Gamp/Gam)/Gam              ! deriv wrt V
        Bp=-(Ap+1.d0)*Pref+A*Prefp*rho**2       ! deriv wrt V
        Hprime=amumg*Bp/A                       ! deriv wrt V
        Hprime=-Hprime/rho**2                   ! deriv wrt rho
c
      elseif (k.eq.3.or.k.eq.6) then
c
c Ideal gas
        Hprime=0.d0
c
      elseif (k.eq.4) then
c
c Stiffened gas
        V=x
        if (irho(j).ne.0) then
          write(6,*)'Error (Hprime), irho.ne.0, k=',k
          stop
        end if
        Hprime=gam4*pi4*(V**gm14)/gm14
c
      elseif (k.eq.5) then
c
c Virial gas
        Hprime=0.d0
c
      else
c
        write(6,*)'Error (Hprime) : k not supported, k=',k
        stop
c
      end if
c
      return
      end
c
c+++++++++++++++++++
c
      double precision function Aprime (k,V)
c
c compute Aprime for a given specific volume V, for material k.
c
      implicit real*8 (a-h,o-z)
      include 'PBX9501reactant-scaled.h'
      include 'HMXproducts-mod-scaled.h'
      include 'IdealGas.h'
      include 'StiffenedGas.h'
      include 'VirialGas.h'
      include 'IdealGas1p35.h'
c
      if (k.eq.1) then
c
c PBX-9501 reactant
        z=V/V0
        Aprime=a1s/(a1s+a2s*z)**2
c
      elseif (k.eq.2) then
c
c HMX products
        rho=1.d0/V
        Gam=gam0
        Gamp=0.d0                               ! deriv wrt rho
        if (rho.gt.rhoc) then
          z=rho-rhoc
          Gam=Gam+(d2g+d3g*z)*z**2
          Gamp=(d2gp+d3gp*z)*z
        end if
        Aprime=(1.d0+rho*Gamp/Gam)/Gam          ! deriv wrt V
c
      elseif (k.eq.3) then
c
c Ideal gas
        Aprime=1.d0/gm13
c
      elseif (k.eq.4) then
c
c Stiffened gas
        Aprime=1.d0/gm14
c
      elseif (k.eq.5) then
c
c Virial gas
        Aprime=V*(V+2*b5)/(gm15*(V+b5)**2)
c
      elseif (k.eq.6) then
c
c Ideal gas
        Aprime=1.d0/gm16
c
      else
c
        write(6,*)'Error (Aprime) : k not supported, k=',k
        stop
c
      end if
c
      return
      end
c      
c+++++++++++++++++++
c
      subroutine energymg (j,rho,p,energy,er,ep,ideriv)
c
c Compute the internal energy per unit mass, energy=A(rho)*p+B(rho), for a given density (rho)
c and pressure (p).  If ideriv.ne.0, then return the derivatives of energy with respect to
c density (er) and pressure (ep).
c
      implicit real*8 (a-h,o-z)
      common / mgEOSb / irho(2),nmg(2),xmin(2),xmax(2),
     *                  Bmg(0:3,0:200000,2)
      common / matdatmg / kmat(2)
c
      k=kmat(j)
      Vol=1.d0/rho
      Ak=AmgV(k,Vol)
      Bk=BmgV(j,Vol,ier)
      if (ier.ne.0) then
        write(6,*)'Error (energymg) : error computing Bmg'
        stop
      end if
      energy=Ak*p+Bk
c
      er=0.d0
      ep=0.d0
      if (ideriv.eq.0) return
c
      Akp=-Aprime(k,Vol)/rho**2
      if (irho(j).eq.0) then
        Bkp=-Bprime(j,k,Vol)/rho**2   ! Bprime = dB/dV
      else
        Bkp= Bprime(j,k,rho)          ! Bprime = dB/drho
      end if
      er=Akp*p+Bkp
      ep=Ak
c
      return
      end
c      
c+++++++++++++++++++
c
      subroutine matlmg (k,rho,z,zp,ideriv)
c
c Compute A=z(1), B=z(2), mu=z(3) and H=z(4) for material k for a given density rho.
c If ideriv.ne.0, then return the derivatives of these quantities with respect to rho.
c
c This subroutine is needed for the exact calculation of the Jacobian matrix associated
c with the solid contact jump conditions.  Only the gas phase is needed for this.
c
      implicit real*8 (a-h,o-z)
      dimension z(4),zp(4)
      common / mgEOSb / irhob(2),nmgb(2),xminb(2),xmaxb(2),
     *                  Bmg(0:3,0:200000,2)
      common / mgEOSh / irhoh(2),nmgh(2),xminh(2),xmaxh(2),
     *                  Hmg(0:3,0:200000,2)
c
c the gas phase is assumed
      j=2
c
      Vol=1.d0/rho
      z(1)=AmgV(k,Vol)
      z(2)=BmgV(j,Vol,ier)
      if (ier.ne.0) then
        write(6,*)'Error (matlmg) : error computing Bmg'
        stop
      end if
      z(3)=amumg(k,Vol)
      z(4)=HmgV(j,Vol,ier)
      if (ier.ne.0) then
        write(6,*)'Error (matlmg) : error computing Hmg'
        stop
      end if
c
      do i=1,4
        zp(i)=0.d0
      end do
c
      if (ideriv.eq.0) return
c
      fact=-Vol**2
      Amgp=Aprime(k,Vol)
      zp(1)=Amgp*fact
      zp(3)=(z(3)*(1.d0+Amgp)/z(1))*fact
c
      if (irhob(j).eq.0) then
        x=Vol
        zp(2)=Bprime(j,k,x)*fact
      else
        x=rho
        zp(2)=Bprime(j,k,x)
      end if
c
      if (irhoh(j).eq.0) then
        x=Vol
        zp(4)=Hprime(j,k,x)*fact
      else
        x=rho
        zp(4)=Hprime(j,k,x)
      end if
c
      return
      end
c      
c+++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine getTmg (j,k,nmg1,xmax1)
c
c Set up a lookup table for the calculation of the temperature.  The basic
c problem involves an integral:
c
c x(y)=int(Cv(z),z=y0..y), where Cv is the specific heat, y0 is a reference
c temperature (or a scaled temperature) and y0 is a constant.
c
c Use a RK4 integrator to solve dy/dx=1/Cv(y), y(0)=y0, for 0<x<xmax1 using nmg1
c equally spaced point.
c
c Build a lookup table for a Hermite interpolation of the results.
c
      implicit real*8 (a-h,o-z)
      parameter (nmgmax=200000)
      common / mgEOST / nmg(2),xmax(2),Tmg(0:3,0:200000,2)
      include 'PBX9501reactant-scaled.h'
c
c certain EOS cases do not require a lookup table
      if (k.ne.1) return

      if (nmg1.gt.nmgmax) then
        write(6,*)'Error (getTmg) : nmg1 too big'
        stop
      end if
c
      nmg(j)=nmg1
      xmax(j)=xmax1
c
c Integrate dy/dx at nodes (using RK4)
      y=y0
      Tmg(0,0,j)=y
      dx=xmax1/nmg1
      do n=1,nmg1
        Tmg(1,n-1,j)=rCv(k,y)
        ak1=dx*Tmg(1,n-1,j)
        ak2=dx*rCv(k,y+.5d0*ak1)
        ak3=dx*rCv(k,y+.5d0*ak2)
        ak4=dx*rCv(k,y+ak3)
        y=y+(ak1+2.d0*(ak2+ak3)+ak4)/6.d0
        Tmg(0,n,j)=y
      end do
c
c construct the rest of the lookup table using divided differences
      do n=0,nmg1-1
        d1=(Tmg(0,n+1,j)-Tmg(0,n,j))/dx
        d2a=(d1-Tmg(1,n,j))/dx
        d2b=(Tmg(1,n+1,j)-d1)/dx
        d3=(d2b-d2a)/dx
        Tmg(2,n,j)=d2a
        Tmg(3,n,j)=d3
      end do
c             
      return
      end
c
c++++++++++++++++++++++++++++++++++
c
      double precision function rCv (k,y)
c
c Compute 1/Cv(y) for material k
c
      implicit real*8 (a-h,o-z)
      include 'PBX9501reactant-scaled.h'
c
      if (k.eq.1) then
        z=1.d0/y
        rCv=c3+z*(c2+z*(c1+z*c0))
      else
        write(6,*)'Error (rCv) : material not supported'
        stop
      end if
c
      return
      end
c
c+++++++++++++++++++++++++++++++++
c
      double precision function TmgV (j,k,V,p,ier)
c
c compute temperature at specific volume V and pressure P
c for phase j, material k
c
      implicit real*8 (a-h,o-z)
      common / mgEOST / nmg(2),xmax(2),Tmg(0:3,0:200000,2)
      data tiny / 1.d-8 /
      include 'PBX9501reactant-scaled.h'
      include 'HMXproducts-mod-scaled.h'
      include 'IdealGas.h'
      include 'StiffenedGas.h'
      include 'VirialGas.h'
      include 'IdealGas1p35.h'
c
      ier=0
c
      if (k.eq.1) then
c
c PBX-9501 reactants, compute x=(Pref-p)/thetap
        z=V/V0
        A=V/(a1s+a2s*z)
        xs=f*(z**e1-z**e2)
        ys=1.d0+g*(z**e3-1.d0)
        Pref=xs*ys
        theta=theta0*dexp(b*(1.d0-z))/z**a1s
        thetap=-theta/A
        x=(Pref-p)/thetap-cIntT
c
c check bounds
        if (x.le.0.d0) then
          ier=1
          write(6,*)'Error (TmgV) : x<0'
          return
        end if
c
        if (x.ge.xmax(j)) then
          ier=1
          write(6,*)'Error (TmgV) : x>xmax'
          return
        end if
c
c determine left index
        dx=xmax(j)/nmg(j)
        n=(x-tiny)/dx
        n=max(min(n,nmg(j)-1),0)
        x1=n*dx
c
c compute TmgV using Hermite interpolation
        z=x-x1
        TmgV=Tmg(0,n,j)+z*(Tmg(1,n,j)+z*(Tmg(2,n,j)+(z-dx)*Tmg(3,n,j)))
c
c scale the result by theta
        TmgV=theta*TmgV
c
      elseif (k.eq.2) then
c
c HMX products
        rho=1.d0/V
        Gam=gam0
        Tref=tcjgas*(rho/rhor)**gam0
        if (rho.gt.rhoc) then
          z=rho-rhoc
          Gam=Gam+(d2g+d3g*z)*z**2
          arg=(e0g+rho*(e1g+rho*e2g))*z
          Tref=Tref*((rho/rhoc)**d0gas)/dexp(arg)
        end if
        Pref=Ag*(rho**qq)*dexp(akg*(rho-rhom))
        TmgV=Tref+(p-Pref)/(rho*cvgas*Gam)
c
      elseif (k.eq.3) then
c
c Ideal Gas
        TmgV=(p*V/gm13)/Cv3
c
      elseif (k.eq.4) then
c
c Stiffened Gas
        TmgV=((p+pi4)*V/gm14)/Cv4
c
      elseif (k.eq.5) then
c
c Virial Gas
        fact=1.d0+b5/V
        TmgV=(p*V/(gm15*fact))/Cv5
c
      elseif (k.eq.6) then
c
c Ideal Gas with gamma=1.35
        TmgV=(p*V/gm16)/Cv6
c
      else
        write(6,*)'Error (TmgV) : Material k not supported'
        stop
      end if
c
      return
      end
c
c++++++++++++++++++++++++++++++++++
c
      double precision function Cvmg (k,V,T)
c
c Compute Cv=de/dT at constant volume for material k as a function
c of volume V and temperature T.
c
      implicit real*8 (a-h,o-z)
      include 'PBX9501reactant-scaled.h'
      include 'HMXproducts-mod-scaled.h'
      include 'IdealGas.h'
      include 'StiffenedGas.h'
      include 'VirialGas.h'
      include 'IdealGas1p35.h'
c
      if (k.eq.1) then
c
c PBX-9501 reactants
        z=V/V0
        theta=theta0*dexp(b*(1.d0-z))/z**a1s
        y=T/theta
        Cvmg=1.d0/rCv(k,y)
c
      elseif (k.eq.2) then
c
c HMX products
        Cvmg=cvgas
c
      elseif (k.eq.3) then
c
c Ideal Gas
        Cvmg=Cv3
c
      elseif (k.eq.4) then
c
c Stiffened Gas
        Cvmg=Cv4
c
      elseif (k.eq.5) then
c
c Virial Gas
        Cvmg=Cv5
c
      elseif (k.eq.6) then
c
c Ideal Gas with gamma=1.35
        Cvmg=Cv6
c
      else
        write(6,*)'Error (Cvmg) : Material k not supported'
        stop
      end if

      return
      end
