      subroutine cmpfluxMG (m,aj,a0,a1,a2,wl,wr,fl,fr,
     *                      speed,method,linsc)
c
c Compute Godunov fluxes for BN equations
c
c Input: wl,wr = left and right states in primitive variables
c
c Output: wl,wr = star-state and right-of-solid-contact state
c                 if solid contact velocity > 0
c               = left-of-solid-contact state and star-state
c                 if solid contact velocity < 0
c         fl,fr = fluxes including contribution from the solid
c                 contact layer
c         speed = largest wave speed (for time step calculation)
c
      implicit double precision (a-h,o-z)
      dimension wl(m),wr(m),fl(m),fr(m)
      dimension an(2),rm1(2,2),vm1(2,2),pm1(2,2),sm1(2,2)
      dimension rm(2,2),vm(2,2),pm(2,2),sm(2,2),fp(2,2),gp(2,2),dv(2,2)
      logical mgfxdbg
      common / mgflxerr / mgfxerr,mgfxdbg
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
      common / matdatmg / kmat(2)
c
c set flux debug flag
      mgfxdbg=.false.
      mgfxerr=0
c
c normalize metrics of the mapping
      rad=dsqrt(a1**2+a2**2)
c     an0=a0/rad            ! moving grids not implemented
      an(1)=a1/rad
      an(2)=a2/rad
c
c assign left state
      r(1,1)=wl(1)                    ! density, solid
      r(1,2)=wl(5)                    ! density, gas
      v(1,1)=an(1)*wl(2)+an(2)*wl(3)  ! velocity, solid
      v(1,2)=an(1)*wl(6)+an(2)*wl(7)  ! velocity, gas
      p(1,1)=wl(4)                    ! pressure, solid
      p(1,2)=wl(8)                    ! pressure, gas
      alpha(1,1)=wl(9)                ! vol, solid
      alpha(1,2)=1.d0-wl(9)           ! vol, gas
c
c assign right state
      r(2,1)=wr(1)                    ! density, solid
      r(2,2)=wr(5)                    ! density, gas
      v(2,1)=an(1)*wr(2)+an(2)*wr(3)  ! velocity, solid
      v(2,2)=an(1)*wr(6)+an(2)*wr(7)  ! velocity, gas
      p(2,1)=wr(4)                    ! pressure, solid
      p(2,2)=wr(8)                    ! pressure, gas
      alpha(2,1)=wr(9)                ! vol, solid
      alpha(2,2)=1.d0-wr(9)           ! vol, gas
c
c compute sound speeds, gamma and stiffening pressures
      do j=1,2
        kj=kmat(j)
        do i=1,2
          if (r(i,j).le.0.d0) then
            write(6,*)'Error (cmpfluxMG) : r(i,j).le.0, i,j=',i,j
            stop
          end if
          Volij=1.d0/r(i,j)
          cm2=c2mg(kj,Volij,p(i,j))
          if (cm2.lt.0.d0) then
            write(6,*)'Error (cmpfluxMG) : cm2.lt.0, i,j=',i,j
            stop
          end if
          c(i,j)=dsqrt(cm2)
          gamk(i,j)=Volij/AmgV(kj,Volij)+1.d0
          p0k(i,j)=max(r(i,j)*c(i,j)**2/gamk(i,j)-p(i,j),0.d0)
        end do
      end do
c
c method=0 => Full MG coupled middle state with MG flux
c method=1 => Approx stiffened middle state with MG flux
c method=2 => Approx stiffened gas with HLLC flux
      if (method.eq.0) then
c
c..compute decoupled middle state for each phase (full MG EOS)
        do j=1,2
          call cmpmiddlemg (j,rm(1,j),vm(1,j),pm(1,j),
     *                      fp(1,j),gp(1,j),sm(1,j))
        end do
c
        call cmpcouplemg (rm0,rm,vm,pm,fp,gp,sm,linsc)
c
c..compute fluxes
        aj1=rad*aj
        call cmpflxmg (m,aj1,an,rm0,rm,vm,pm,sm,wl,wr,fl,fr,method)
c
      elseif (method.eq.1) then
c
c..compute decoupled middle state for each phase (stiffened ideal EOS)
        do j=1,2
          call cmpmiddlesi (j,pm(1,j),vm(1,j),dv(1,j))
        end do
c
c..compute coupled middle states (approx stiffened), if necessary
        call cmpcouplesi (pm,vm,dv,rm0,rm,sm,linsc)
c
c..compute fluxes
        aj1=rad*aj
        call cmpflxmg (m,aj1,an,rm0,rm,vm,pm,sm,wl,wr,fl,fr,method)
c
      elseif (method.eq.2) then
c     
c..compute coupled middle pressures for approx stiffened and HLLC
        call cmpcouplehc (rm0,rm,vm,pm,sm,linsc)
c     
c..compute fluxes
        aj1=rad*aj
        call cmpflxhc (m,aj1,an,rm0,rm,vm,pm,sm,wl,wr,fl,fr)
c
      else
        write(6,*)'Error (gdfluxMG) : value of method not supported'
        stop
      end if
c
c..compute maximum wave speed
      speed=0.d0
      do i=1,2
        do j=1,2
          speed=max(sm(i,j)+dabs(vm(i,j)),c(i,j)+dabs(v(i,j)),speed)
        end do
      end do
      speed=rad*speed
c
c check debug flag
      if (mgfxerr.ne.0) then
        write(6,*)'Error (gdfluxMG) : mgfxerr.ne.0'
        write(6,*)'mgfxerr =',mgfxerr
        stop
      end if
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c HLLC flux subroutines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cmpcouplehc (rm0,rm,vm,pm,sm,linsc)
c
      implicit real*8 (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),sm(2,2)
      dimension z(2),dz(2),fp(2,2),gp(2),iq(2,2),a(4,5),dpm(2,2)
      dimension rmsave(2,2),vmsave(2,2),pmsave(2,2),smsave(2,2)
      dimension aa(4,5)
      logical icorr
      logical mgfxdbg
      common / mgflxerr / mgfxerr,mgfxdbg
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
c      common / flxdat / afxcnt(3),nfxcnt(3)
      data atol, tol, itmax, itcorr / 1.d-14, 1.d-8, 10, 2 /
      data sig, frac / 10.d0, 0.99d0 /
      data icorr / .true. /
      data idebug / 0 /
c
c      nfxcnt(1)=nfxcnt(1)+1
c
c loop over phases to compute decoupled middle pressures
      do j=1,2
c
c estimates for gamma, gm1 and p0 based on averages
        gam(j)=.5d0*(gamk(1,j)+gamk(2,j))
        gm1(j)=gam(j)-1.d0
        p0(j)=.5d0*(p0k(1,j)+p0k(2,j))
c
c recompute sound speeds
        do i=1,2
           c2mg=gam(j)*(p(i,j)+p0(j))/r(i,j)
           if (c2mg.le.0.d0) then
              write(6,*)'Error (cmpcouplehc) : c2mg.le.0, i,j=',i,j
              stop
           end if
           c(i,j)=dsqrt(c2mg)
        end do
c
c impedances and derivatives
        z(1)=r(1,j)*c(1,j)
        z(2)=r(2,j)*c(2,j)
        dz(1)=0.d0
        dz(2)=0.d0
c
c pstar
        pstar=(z(2)*p(1,j)+z(1)*p(2,j)
     *        -z(1)*z(2)*(v(2,j)-v(1,j)))/(z(1)+z(2))
c
c initialize the "q" flag
        iq(1,j)=0
        iq(2,j)=0
c
c possibly correct wave-speed estimates
        if (icorr) then
          ep=0.5d0*(gam(j)+1.d0)/gam(j)
          do i=1,2
            if (pstar.gt.p(i,j)) then
              q=dsqrt(1.d0+ep*((pstar+p0(j))/(p(i,j)+p0(j))-1.d0))
              dq=.5d0*ep/(q*(p(i,j)+p0(j)))
              dz(i)=z(i)*dq
              z(i)=z(i)*q
              iq(i,j)=1
            end if
          end do
          pstar=(z(2)*p(1,j)+z(1)*p(2,j)
     *           -z(1)*z(2)*(v(2,j)-v(1,j)))/(z(1)+z(2))
        end if
c
        pm(1,j)=pstar
        pm(2,j)=pstar
        vm(1,j)=v(1,j)-(pstar-p(1,j))/z(1)
        vm(2,j)=vm(1,j)
        fp(1,j)=(1.d0-(pm(1,j)-p(1,j))*dz(1)/z(1))/z(1)
        fp(2,j)=(1.d0-(pm(2,j)-p(2,j))*dz(2)/z(2))/z(2)
        sm(1,j)=v(1,j)-z(1)/r(1,j)
        sm(2,j)=v(2,j)+z(2)/r(2,j)
        do i=1,2
          rm(i,j)=r(i,j)*z(i)**2/(z(i)**2-r(i,j)*(pm(i,j)-p(i,j)))
        end do
        if (j.eq.2) then
          do i=1,2
            gp(i)=(1.d0-2.d0*(pm(i,j)-p(i,j))*dz(i)/z(i))
     *            *(rm(i,j)/z(i))**2
          end do
        end if
c
      end do
c
      call setrm0mg (rm0,rm,vm,iconf)
c
c use Newton's method to couple the phases (if necessary)
      if (dabs(alpha(1,1)-alpha(2,1)).lt.atol) return
c
c start with linearized solid contact, compute velocity difference
      dvm=vm(1,2)-vm(1,1)
c
c pressure scales (just the pstar's for each phase plus stiffening)
      pscale1=pm(1,1)+p0(1)
      pscale2=pm(1,2)+p0(2)
c
c compute alpha-bar_m and alpha_m
      if (alpha(1,1).gt.0.5d0) then
        if (alpha(2,1).lt.0.5d0) then
          asm=.5d0
        else
          asm=min(alpha(1,1),alpha(2,1))
        end if
      else
        if (alpha(2,1).gt.0.5d0) then
          asm=.5d0
        else
          asm=max(alpha(1,1),alpha(2,1))
        end if
      end if
      agm=1.d0-asm
c
c linear update for middle solid pressures
      eps=alpha(2,1)-alpha(1,1)
      del=(pm(1,2)-pm(1,1))*eps/(asm*(fp(1,1)+fp(2,1)))
      dpm1=-fp(2,1)*del
      dpm2= fp(1,1)*del
      pm(1,1)=pm(1,1)+dpm1
      pm(2,1)=pm(2,1)+dpm2
      err1=max(dabs(dpm1),dabs(dpm2))/pscale1
c
c limit velocity difference
      gampg=gam(2)*(pm(1,2)+p0(2))
      cm=frac*dsqrt(gampg/rm0)
      arg1=sig*(1.d0+dvm/cm)
      arg2=sig*(1.d0-dvm/cm)
      dvm=.5d0*cm*dlog(dcosh(arg1)/dcosh(arg2))/sig
c
c linear update for middle gas pressures
      del=gampg*dvm*eps/(agm*(gampg-rm0*dvm**2)*(fp(1,2)+fp(2,2)))
      dpm1= (rm0*dvm*fp(2,2)+1.d0)*del
      dpm2=-(rm0*dvm*fp(1,2)-1.d0)*del
      pm(1,2)=pm(1,2)+dpm1
      pm(2,2)=pm(2,2)+dpm2
      err2=max(dabs(dpm1),dabs(dpm2))/pscale2
c
c update middle velocities and densities
      call updatehc (rm,vm,pm,fp,gp,sm,iq)
c
      if (iconf.eq.1) then
c gas contact on the right
        rm0=rm(1,2)*((pm(2,2)+p0(2))
     *              /(pm(1,2)+p0(2)))**(1.d0/gam(2))
      else
c gas contact on the left
        rm0=rm(2,2)*((pm(1,2)+p0(2))
     *              /(pm(2,2)+p0(2)))**(1.d0/gam(2))
      end if
c
c check sonic condition on the middle states
      it=0
      isonic=0
      call chkvelmg (it,rm(1,2),vm,pm(1,2),isonic)
c
c check errors
      if (max(err1,err2).lt.tol.or.linsc.ne.0.or.isonic.ne.0) then
c        nfxcnt(2)=nfxcnt(2)+1
        return
      end if
c
      err1save=err1
      err2save=err2
c
c save the current state
      rm0save=rm0
      do i=1,2
         do j=1,2
            rmsave(i,j)=rm(i,j)
            vmsave(i,j)=vm(i,j)
            pmsave(i,j)=pm(i,j)
            smsave(i,j)=sm(i,j)
         end do
      end do
c
c Newton iteration to find middle pressures
      do it=1,itmax
c
c get Jacobian matrix (using finite differences)
c        call jacobhcfd (rm,vm,pm,aa,iconf,iq)
c
c get Jacobian matrix (exact)
        call jacobhc (rm,vm,pm,fp,gp,a,iconf)
c
c        write(77,*)' it,iconf=',it,iconf
c        do iii=1,4
c          write(77,777)(dabs(aa(iii,jjj)-a(iii,jjj)),jjj=1,5)
c  777     format(5(1x,1pe10.3))
c        end do
c        do iii=1,4
c          write(77,778)(a(iii,jjj),jjj=1,5)
c  778     format(5(1x,1pe15.8))
c        end do
c        do iii=1,4
c          write(77,778)(aa(iii,jjj),jjj=1,5)
c        end do
c
c
c compute corrections to the middle pressures
        call solvehc (a,dpm)
c
c update middle pressures
        pmmin=p(1,1)-dpm(1,1)
        do j=1,2
          do i=1,2
            pm(i,j)=pm(i,j)-dpm(i,j)
            pmmin=min(pm(i,j),pmmin)
          end do
        end do
c
c check for negative pressures
        if (pmmin.le.0.d0) then
          if (mgfxdbg) then
            write(6,*)'Warning (cmpcouplehc) : pm.le.0'
            write(6,*)'ul, ur ='
            write(6,*)alpha(1,1),alpha(2,1)
            write(6,*)r(1,1),r(2,1)
            write(6,*)v(1,1),v(2,1)
            write(6,*)p(1,1),p(2,1)
            write(6,*)r(1,2),r(2,2)
            write(6,*)v(1,2),v(2,2)
            write(6,*)p(1,2),p(2,2)
            write(6,*)'err1save, err2save=',err1save,err2save
            mgfxerr=1
          end if
          rm0=rm0save
          do i=1,2
            do j=1,2
              rm(i,j)=rmsave(i,j)
              vm(i,j)=vmsave(i,j)
              pm(i,j)=pmsave(i,j)
              sm(i,j)=smsave(i,j)
            end do
          end do
          return
        end if
c
c update iq, if desired
        if (it.le.itcorr.and.icorr) then
          do i=1,2
            do j=1,2
              if (pm(i,j).gt.p(i,j)) then
                iq(i,j)=1
              else
                iq(i,j)=0
              end if
            end do
          end do
        end if
c
c update middle velocities and densities
        call updatehc (rm,vm,pm,fp,gp,sm,iq)
c
c check for negative densities
        rmmin=rm(1,1)
        do i=1,2
          do j=1,2
            rmmin=min(rm(i,j),rmmin)
          end do
        end do
        if (rmmin.le.0.d0) then
          if (mgfxdbg) then
            write(6,*)'Warning (cmpcouplehc) : rm.le.0'
            write(6,*)'ul, ur ='
            write(6,*)alpha(1,1),alpha(2,1)
            write(6,*)r(1,1),r(2,1)
            write(6,*)v(1,1),v(2,1)
            write(6,*)p(1,1),p(2,1)
            write(6,*)r(1,2),r(2,2)
            write(6,*)v(1,2),v(2,2)
            write(6,*)p(1,2),p(2,2)
            write(6,*)'err1save, err2save=',err1save,err2save
            mgfxerr=2
          end if
          rm0=rm0save
          do i=1,2
            do j=1,2
              rm(i,j)=rmsave(i,j)
              vm(i,j)=vmsave(i,j)
              pm(i,j)=pmsave(i,j)
              sm(i,j)=smsave(i,j)
            end do
          end do
          return
        end if
c
c check sonic condition on the middle states
        isonic=0
        call chkvelmg (it,rm(1,2),vm,pm(1,2),isonic)
c
        if (isonic.ne.0) then
          if (mgfxdbg) then
            write(6,*)'Warning (cmpcouplehc) : error from chkvelmg'
            write(6,*)'ul, ur ='
            write(6,*)alpha(1,1),alpha(2,1)
            write(6,*)r(1,1),r(2,1)
            write(6,*)v(1,1),v(2,1)
            write(6,*)p(1,1),p(2,1)
            write(6,*)r(1,2),r(2,2)
            write(6,*)v(1,2),v(2,2)
            write(6,*)p(1,2),p(2,2)
            write(6,*)'err1save, err2save=',err1save,err2save
            mgfxerr=3
          end if
          rm0=rm0save
          do i=1,2
            do j=1,2
              rm(i,j)=rmsave(i,j)
              vm(i,j)=vmsave(i,j)
              pm(i,j)=pmsave(i,j)
              sm(i,j)=smsave(i,j)
            end do
          end do
          return
        end if
c
c error in solid and gas middle pressures
        err1=max(dabs(dpm(1,1)),dabs(dpm(2,1)))/pscale1
        err2=max(dabs(dpm(1,2)),dabs(dpm(2,2)))/pscale2
c
c print convergence information
        if (idebug.ne.0) then
          write(6,100)it,(pm(i,1),dpm(i,1),i=1,2),err1,
     *                   (pm(i,2),dpm(i,2),i=1,2),err2
  100     format(' it=',i2,2(2x,1pe15.8,1x,1pe10.3),2x,1pe10.3,/,
     *                  6x,2(2x,1pe15.8,1x,1pe10.3),2x,1pe10.3)
        end if
c
c check for convergence
        if (max(err1,err2).lt.tol) then
c        if (max(err1,err2).lt.tol.or.linsc.ne.0) then
          if (iconf.eq.1) then
c gas contact on the right
            rm0=rm(1,2)*((pm(2,2)+p0(2))
     *                  /(pm(1,2)+p0(2)))**(1.d0/gam(2))
          else
c gas contact on the left
            rm0=rm(2,2)*((pm(1,2)+p0(2))
     *                  /(pm(2,2)+p0(2)))**(1.d0/gam(2))
          end if
c          nfxcnt(3)=nfxcnt(3)+1
          return
        else
          if (it.eq.1) then
            if (iconf.eq.1) then
c gas contact on the right
              rm0save=rm(1,2)*((pm(2,2)+p0(2))
     *                  /(pm(1,2)+p0(2)))**(1.d0/gam(2))
            else
c gas contact on the left
              rm0save=rm(2,2)*((pm(1,2)+p0(2))
     *                  /(pm(2,2)+p0(2)))**(1.d0/gam(2))
            end if
            do i=1,2
              do j=1,2
                rmsave(i,j)=rm(i,j)
                vmsave(i,j)=vm(i,j)
                pmsave(i,j)=pm(i,j)
                smsave(i,j)=sm(i,j)
              end do
            end do
          end if
        end if
c
      end do
c
      if (mgfxdbg) then 
        write(6,*)'Warning (cmpcouplehc) : itmax exceeded'
        write(6,*)'ul, ur ='
        write(6,*)alpha(1,1),alpha(2,1)
        write(6,*)r(1,1),r(2,1)
        write(6,*)v(1,1),v(2,1)
        write(6,*)p(1,1),p(2,1)
        write(6,*)r(1,2),r(2,2)
        write(6,*)v(1,2),v(2,2)
        write(6,*)p(1,2),p(2,2)
        write(6,*)'err1save, err2save=',err1save,err2save
        mgfxerr=4
      end if
      rm0=rm0save
      do i=1,2
         do j=1,2
            rm(i,j)=rmsave(i,j)
            vm(i,j)=vmsave(i,j)
            pm(i,j)=pmsave(i,j)
            sm(i,j)=smsave(i,j)
         end do
      end do
c
      return
      end
c
c+++++++++++++++
c
      subroutine jacobhcfd (rm,vm,pm,a,iconf,iq)
c
c compute approximate Jacobian matrix of the solid contact jump conditions
c using finite differences
c
      implicit real*8 (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),a(4,5),iq(2,2)
      dimension fp(2,2),gp(2),sm(2,2),b(4),bp(4)
      data delta / 1.d-7 /
c
c get residual of the solid-contact jump conditions
      call updatehc (rm,vm,pm,fp,gp,sm,iq)
      call residhc (rm,vm,pm,b,iconf)
c      bmax=dabs(b(1))
c      do n=2,4
c        bmax=max(dabs(b(n)),bmax)
c      end do
c      write(6,*)'bmax=',bmax
c
      do j=1,2
        do i=1,2
          n=2*(j-1)+i
          dp=pm(i,j)*delta
          pm(i,j)=pm(i,j)+dp
          call updatehc (rm,vm,pm,fp,gp,sm,iq)
          call residhc (rm,vm,pm,bp,iconf)
          do m=1,4
            a(m,n)=(bp(m)-b(m))/dp
          end do
          pm(i,j)=pm(i,j)-dp
        end do
      end do
c
      do m=1,4
        a(m,5)=b(m)
      end do
c
      return
      end
c
c+++++++++++++++
c
      subroutine residhc (rm,vm,pm,b,iconf)
c
      implicit real*8 (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),b(4)
      dimension rg(2)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
c
      gi=1.d0/gam(2)
      rat=(pm(1,2)+p0(2))/(pm(2,2)+p0(2))
      if (iconf.eq.1) then
        rg(1)=rm(1,2)
        rg(2)=rm(1,2)/(rat**gi)
      else
        rg(1)=rm(2,2)*(rat**gi)
        rg(2)=rm(2,2)
      end if
c
      dv1=vm(1,2)-vm(1,1)      ! dv1=v(1,2)-v(1,1)-f(1,2)+f(1,1)
      dv2=vm(2,2)-vm(2,1)      ! dv2=v(2,2)-v(2,1)+f(2,2)-f(2,1)
      h1=gam(2)*(pm(1,2)+p0(2))/(gm1(2)*rg(1))
      h2=gam(2)*(pm(2,2)+p0(2))/(gm1(2)*rg(2))
c
      b(1)= vm(2,1)
     *     -vm(1,1)
      b(2)= alpha(2,2)*rg(2)*dv2
     *     -alpha(1,2)*rg(1)*dv1
      b(3)= alpha(2,1)*pm(2,1)+alpha(2,2)*(pm(2,2)+rg(2)*dv2**2)
     *     -alpha(1,1)*pm(1,1)-alpha(1,2)*(pm(1,2)+rg(1)*dv1**2)
      b(4)= h2+.5d0*dv2**2
     *     -h1-.5d0*dv1**2
c
      return
      end
c
c+++++++++++++++
c
      subroutine jacobhc (rm,vm,pm,fp,gp,a,iconf)
c
      implicit real*8 (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),fp(2,2),gp(2),a(4,5)
      dimension rg(2),drdp(2,2)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
c
      gi=1.d0/gam(2)
      rat=(pm(1,2)+p0(2))/(pm(2,2)+p0(2))
      if (iconf.eq.1) then
        rg(1)=rm(1,2)
        rg(2)=rm(1,2)/(rat**gi)
        drdp(1,1)=gp(1)                         ! drdp(i1,i2) = d(rg(i1))/d(pm(i2,2))
        drdp(1,2)=0.d0
        drdp(2,1)=rg(2)*(gp(1)/rm(1,2)-gi/(pm(1,2)+p0(2)))
        drdp(2,2)=rg(2)*gi/(pm(2,2)+p0(2))
      else
        rg(1)=rm(2,2)*(rat**gi)
        rg(2)=rm(2,2)
        drdp(1,1)=rg(1)*gi/(pm(1,2)+p0(2))      ! drdp(i1,i2) = d(rg(i1))/d(pm(i2,2))
        drdp(1,2)=rg(1)*(gp(2)/rm(2,2)-gi/(pm(2,2)+p0(2)))
        drdp(2,1)=0.d0
        drdp(2,2)=gp(2)
      end if
c
      dv1=vm(1,2)-vm(1,1)                       ! dv1=v(1,2)-v(1,1)-f(1,2)+f(1,1)
      dv2=vm(2,2)-vm(2,1)                       ! dv2=v(2,2)-v(2,1)+f(2,2)-f(2,1)
      h1=gam(2)*(pm(1,2)+p0(2))/(gm1(2)*rg(1))
      h2=gam(2)*(pm(2,2)+p0(2))/(gm1(2)*rg(2))
c
      a(1,1)=fp(1,1)
      a(1,2)=fp(2,1)
      a(1,3)=0.d0
      a(1,4)=0.d0
      a(1,5)= vm(2,1)
     *       -vm(1,1)
c
      a(2,1)=-alpha(1,2)*rg(1)*fp(1,1)
      a(2,2)=-alpha(2,2)*rg(2)*fp(2,1)
      a(2,3)= alpha(2,2)*(drdp(2,1)*dv2              )
     *       -alpha(1,2)*(drdp(1,1)*dv1-rg(1)*fp(1,2))
      a(2,4)= alpha(2,2)*(drdp(2,2)*dv2+rg(2)*fp(2,2))
     *       -alpha(1,2)*(drdp(1,2)*dv1              )
      a(2,5)= alpha(2,2)*rg(2)*dv2
     *       -alpha(1,2)*rg(1)*dv1
c
      a(3,1)=-alpha(1,1)-2.d0*alpha(1,2)*rg(1)*dv1*fp(1,1)
      a(3,2)= alpha(2,1)-2.d0*alpha(2,2)*rg(2)*dv2*fp(2,1)
      a(3,3)= alpha(2,2)*(     (drdp(2,1)*dv2                   )*dv2)
     *       -alpha(1,2)*(1.d0+(drdp(1,1)*dv1-2.d0*rg(1)*fp(1,2))*dv1)
      a(3,4)= alpha(2,2)*(1.d0+(drdp(2,2)*dv2+2.d0*rg(2)*fp(2,2))*dv2)
     *       -alpha(1,2)*(     (drdp(1,2)*dv1                   )*dv1)
      a(3,5)= alpha(2,1)*pm(2,1)+alpha(2,2)*(pm(2,2)+rg(2)*dv2**2)
     *       -alpha(1,1)*pm(1,1)-alpha(1,2)*(pm(1,2)+rg(1)*dv1**2)
c
      a(4,1)=-dv1*fp(1,1)
      a(4,2)=-dv2*fp(2,1)
      a(4,3)= (             -h2*drdp(2,1))/rg(2)
     *       -(gam(2)/gm1(2)-h1*drdp(1,1))/rg(1)+dv1*fp(1,2)
      a(4,4)= (gam(2)/gm1(2)-h2*drdp(2,2))/rg(2)+dv2*fp(2,2)
     *       -(             -h1*drdp(1,2))/rg(1)
      a(4,5)= h2+.5d0*dv2**2
     *       -h1-.5d0*dv1**2
c
c      dg(2,1)=-alpha(1,2)*(pm(1,2)**gi)*fp(1,1)
c      dg(2,2)=-alpha(2,2)*(pm(2,2)**gi)*fp(2,1)
c      dg(2,3)=alpha(1,2)*(pm(1,2)**gi)*(fp(1,2)-dv1*gi/pm(1,2))
c      dg(2,4)=alpha(2,2)*(pm(2,2)**gi)*(fp(2,2)+dv2*gi/pm(2,2))
c      dg(3,1)=-alpha(1,1)-2.d0*alpha(1,2)*rm(1)*dv1*fp(1,1)
c      dg(3,2)=alpha(2,1)-2.d0*alpha(2,2)*rm(2)*dv2*fp(2,1)
c      dg(3,3)=alpha(2,2)*drdp(2,1)*(dv2**2)-alpha(1,2)*
c     *                 (1.d0+drdp(1,1)*(dv1**2)-2.d0*rm(1)*dv1*fp(1,2))
c      dg(3,4)=-alpha(1,2)*drdp(1,2)*(dv1**2)+alpha(2,2)*
c     *                 (1.d0+drdp(2,2)*(dv2**2)+2.d0*rm(2)*dv2*fp(2,2))
c      dg(4,1)=-dv1*fp(1,1)
c      dg(4,2)=-dv2*fp(2,1)
c      dg(4,3)=dv1*fp(1,2)-gam(2)/gm1(2)*(1.0d0/rm(1)+
c     *                          drdp(2,1)*(pm(2,2)+p0(2))/(rm(2)**2)-
c     *                          drdp(1,1)*(pm(1,2)+p0(2))/(rm(1)**2))
c      dg(4,4)=dv2*fp(2,2)+gam(2)/gm1(2)*(1.0d0/rm(2)+
c     *                          drdp(1,2)*(pm(1,2)+p0(2))/(rm(1)**2)-
c     *                          drdp(2,2)*(pm(2,2)+p0(2))/(rm(2)**2))
c
      return
      end
c
c+++++++++++++++
c
      subroutine solvehc (a,dpm)
c
      implicit real*8 (a-h,o-z)
      dimension a(4,5),dpm(2,2)
      data tol / 1.d-12 /
c
c Gaussian elimination with partial pivoting
      do k=1,3
        kpiv=k
        apiv=dabs(a(k,k))
        do i=k+1,4
          if (dabs(a(i,k)).gt.apiv) then
            kpiv=i
            apiv=dabs(a(i,k))
          end if
        end do
        if (kpiv.ne.k) then
          do j=k,5
            atmp=a(k,j)
            a(k,j)=a(kpiv,j)
            a(kpiv,j)=atmp
          end do
        end if
        if (dabs(a(k,k)).lt.tol) then
          write(6,*)'Error (solvehc) : matrix singular??'
          stop
        end if
        do i=k+1,4
          fact=a(i,k)/a(k,k)
          do j=k+1,5
            a(i,j)=a(i,j)-fact*a(k,j)
          end do
        end do
      end do
c
c backward substitution
      a(4,5)=a(4,5)/a(4,4)
      do i=3,1,-1
        sum=a(i,5)
        do j=i+1,4
          sum=sum-a(i,j)*a(j,5)
        end do
        a(i,5)=sum/a(i,i)
      end do
c
c assign corrections
      dpm(1,1)=a(1,5)
      dpm(2,1)=a(2,5)
      dpm(1,2)=a(3,5)
      dpm(2,2)=a(4,5)
c
      return
      end
c
c+++++++++++++++
c
      subroutine updatehc (rm,vm,pm,fp,gp,sm,iq)
c
      implicit real*8 (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),fp(2,2),gp(2),sm(2,2),iq(2,2)
      dimension z(2),dz(2)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
c
c loop over phases
      do j=1,2
c
c impedances and derivatives
        z(1)=r(1,j)*c(1,j)
        z(2)=r(2,j)*c(2,j)
        dz(1)=0.d0
        dz(2)=0.d0
c
c possibly correct wave-speed estimates
        ep=0.5d0*(gam(j)+1.d0)/gam(j)
        do i=1,2
          if (iq(i,j).ne.0) then
            q=dsqrt(1.d0+ep*((pm(i,j)+p0(j))/(p(i,j)+p0(j))-1.d0))
            dq=.5d0*ep/(q*(p(i,j)+p0(j)))
            dz(i)=z(i)*dq
            z(i)=z(i)*q
          end if
        end do
c
        vm(1,j)=v(1,j)-(pm(1,j)-p(1,j))/z(1)
        vm(2,j)=v(2,j)+(pm(2,j)-p(2,j))/z(2)
        fp(1,j)=(1.d0-(pm(1,j)-p(1,j))*dz(1)/z(1))/z(1)
        fp(2,j)=(1.d0-(pm(2,j)-p(2,j))*dz(2)/z(2))/z(2)
        sm(1,j)=v(1,j)-z(1)/r(1,j)
        sm(2,j)=v(2,j)+z(2)/r(2,j)
        do i=1,2
          rm(i,j)=r(i,j)*z(i)**2/(z(i)**2-r(i,j)*(pm(i,j)-p(i,j)))
        end do
        if (j.eq.2) then
          do i=1,2
            gp(i)=(1.d0-2.d0*(pm(i,j)-p(i,j))*dz(i)/z(i))
     *            *(rm(i,j)/z(i))**2
          end do
        end if
c
      end do
c
      return
      end
c
c+++++++++++++++
c
      subroutine cmpflxhc (m,aj,an,rm0,rm,vm,pm,sm,wl,wr,fl,fr)
c
c compute fluxes based on coupled middle states
c
      implicit real*8 (a-h,o-z)
      dimension an(2),rm(2,2),vm(2,2),pm(2,2),sm(2,2),
     *          wl(m),wr(m),fl(m),fr(m)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
c
c solid contact source contributions
      source=pm(2,1)*alpha(2,1)-pm(1,1)*alpha(1,1)
      source2=vm(1,1)*source
c
c check if solid contact is to the right or left of x=0
      if (vm(1,1).gt.0.d0) then
c
c tangential component of the solid velocity from left state
        vt=an(1)*wl(3)-an(2)*wl(2)
c
c right-of-solid-contact state (only need solid velocity and gas pressure)
        wr(2)=an(1)*vm(1,1)-an(2)*vt
        wr(3)=an(2)*vm(1,1)+an(1)*vt
        wr(8)=pm(2,2)
c
c solid phase flux components
        call getfxhc (1,1,rm0,rm,vm,pm,sm,aj,an,vt,vt,fl(1))
c
c add on solid contact source contributions
        fr(1)=fl(1)
        fr(2)=fl(2)+aj*an(1)*source
        fr(3)=fl(3)+aj*an(2)*source
        fr(4)=fl(4)+aj*source2
c
c tangential component of the gas velocity from left state
        vt=an(1)*wl(7)-an(2)*wl(6)
c
c gas phase flux components
        if (vm(1,2).ge.0.d0) then
c gas contact (and solid contact) to the right
c compute tangential component of the gas velocity from left state
          call getfxhc (1,2,rm0,rm,vm,pm,sm,aj,an,vt,vt,fl(5))
        else
c gas contact to the left, solid contact to the right
c compute tangential component of the gas velocity from right state
          vt2=an(1)*wr(7)-an(2)*wr(6)
          call getfxhc (1,2,rm0,rm,vm,pm,sm,aj,an,vt,vt2,fl(5))
        end if
c
c add on solid contact source contributions
        fr(5)=fl(5)
        fr(6)=fl(6)-aj*an(1)*source
        fr(7)=fl(7)-aj*an(2)*source
        fr(8)=fl(8)-aj*source2
c
c solid volume fraction flux
        fl(9)= 0.d0
        fr(9)=-aj*vm(1,1)*(alpha(2,1)-alpha(1,1))
c
c flux for advected variables
        do i=10,m
          fl(i)= 0.d0
          fr(i)=-aj*vm(1,1)*(wr(i)-wl(i))
        end do
c
      else
c
c tangential component of the solid velocity from right state
        vt=an(1)*wr(3)-an(2)*wr(2)
c
c left-of-solid-contact state (only need solid velocity and gas pressure)
        wl(2)=an(1)*vm(1,1)-an(2)*vt
        wl(3)=an(2)*vm(1,1)+an(1)*vt
        wl(8)=pm(1,2)
c
c solid phase flux components
        call getfxhc (2,1,rm0,rm,vm,pm,sm,aj,an,vt,vt,fr(1))
c
c add on solid contact source contributions
        fl(1)=fr(1)
        fl(2)=fr(2)-aj*an(1)*source
        fl(3)=fr(3)-aj*an(2)*source
        fl(4)=fr(4)-aj*source2
c
c tangential component of the gas velocity from right state
        vt=an(1)*wr(7)-an(2)*wr(6)
c
c gas phase flux components
        if (vm(2,2).le.0.d0) then
c gas contact (and solid contact) to the left
c compute tangential component of the gas velocity from right state
          call getfxhc (2,2,rm0,rm,vm,pm,sm,aj,an,vt,vt,fr(5))
        else
c gas contact to the right, solid contact to the left
c compute tangential component of the gas velocity from left state
          vt2=an(1)*wl(7)-an(2)*wl(6)
          call getfxhc (2,2,rm0,rm,vm,pm,sm,aj,an,vt,vt2,fr(5))
        end if
c
c add on solid contact source contributions
        fl(5)=fr(5)
        fl(6)=fr(6)+aj*an(1)*source
        fl(7)=fr(7)+aj*an(2)*source
        fl(8)=fr(8)+aj*source2
c
c solid volume fraction flux
        fl(9)= aj*vm(1,1)*(alpha(2,1)-alpha(1,1))
        fr(9)= 0.d0
c
c flux for advected variables
        do i=10,m
          fl(i)= aj*vm(1,1)*(wr(i)-wl(i))
          fr(i)= 0.d0
        end do
c
      end if
c
      return
      end
c
c+++++++++++++++
c
      subroutine getfxhc (i,j,rm0,rm,vm,pm,sm,aj,an,vt,vt2,fx)
c
c compute solid or gas flux for i=side and j=phase
c
      implicit double precision (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),sm(2,2),an(2),fx(4)
      dimension uk(4),um(4),u0(4)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
      common / matdatmg / kmat(2)
c
c set isign=1 for i=1 and isign=-1 for i=2
      isign=3-2*i
c
c compute flux on the left (i=1) or right (i=2)
      Vols=1.d0/r(i,j)
      As=AmgV(kmat(j),Vols)
      Bs=BmgV(j,Vols,ier)
      if (ier.ne.0) then
        write(6,*)'Error (getfxhc) : cannot compute Bs'
        stop
      end if
      enk=As*p(i,j)+Bs
c
c x and y components of velocity in the left/right states
      v1=an(1)*v(i,j)-an(2)*vt
      v2=an(2)*v(i,j)+an(1)*vt
c
c add compaction potential for the solid phase case
      if (j.eq.1) then
        enk=enk+compac(alpha(i,j),0)
      end if
c
c add kinetic energy
      enk=r(i,j)*(enk+.5d0*(v1**2+v2**2))
c
c conserved variables
      uk(1)=alpha(i,j)*r(i,j)
      uk(2)=uk(1)*v1
      uk(3)=uk(1)*v2
      uk(4)=alpha(i,j)*enk
c
c compute flux
      fact=r(i,j)*v(i,j)
      fx(1)=alpha(i,j)*fact
      fx(2)=alpha(i,j)*(fact*v1+an(1)*p(i,j))
      fx(3)=alpha(i,j)*(fact*v2+an(2)*p(i,j))
      fx(4)=alpha(i,j)*(enk+p(i,j))*v(i,j)
c
      if (isign*sm(i,j).lt.0.d0) then
c middle state 1 (i=1) or 2 (i=2)
        um(1)=alpha(i,j)*rm(i,j)
        um(2)=um(1)*(an(1)*vm(i,j)-an(2)*vt)
        um(3)=um(1)*(an(2)*vm(i,j)+an(1)*vt)
        enm=rm(i,j)*(enk/r(i,j)+(vm(i,j)-v(i,j))
     *                *(vm(i,j)+p(i,j)/(r(i,j)*(sm(i,j)-v(i,j)))))
        um(4)=alpha(i,j)*enm
        do n=1,4
          fx(n)=fx(n)+sm(i,j)*(um(n)-uk(n))
        end do
        if (j.eq.2.and.isign*vm(i,j).lt.0.d0) then
c middle state 0 for the gas phase
          u0(1)=alpha(i,j)*rm0
          u0(2)=u0(1)*(an(1)*vm(i,j)-an(2)*vt2)
          u0(3)=u0(1)*(an(2)*vm(i,j)+an(1)*vt2)
          Vol0=1.d0/rm0
          A0=AmgV(kmat(j),Vol0)
          B0=BmgV(j,Vol0,ier)
          if (ier.ne.0) then
            write(6,*)'Error (getfxhc) : cannot compute Bs'
            stop
          end if
          en0=rm0*(A0*pm(i,j)+B0+.5d0*(vm(i,j)**2+vt2**2))
          u0(4)=alpha(i,j)*en0
          do n=1,4
            fx(n)=fx(n)+vm(i,j)*(u0(n)-um(n))
          end do
        end if
      end if
c
c scale by the jacobian
      do n=1,4
        fx(n)=aj*fx(n)
      end do
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Full Riemann flux subroutines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cmpmiddlemg (j,rm,vm,pm,fp,gp,sm)
c
c middle states of the Riemann problem for the Mie-Gruneisen EOS
c
      implicit real*8 (a-h,o-z)
      dimension rm(2),vm(2),pm(2),fp(2),gp(2),sm(2)
      dimension dv(2),drm(2),ishock(2)
      dimension rmsave(2),vmsave(2),pmsave(2),smsave(2)
      character*3 itype(2)
      logical mgfxdbg
      common / mgflxerr / mgfxerr,mgfxdbg
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / romdatmg / aint0(2,2),r0(2,2),c0(2,2),tolrom(2),lvrom(2)
      common / rardatmg / amumg1(2,2),Hmg1(2,2),initr(2,2)
      common / shkdatmg / Amg1(2,2),Bmg1(2,2),inits(2,2)
      common / matdatmg / kmat(2)
      data tol, itmax / 1.d-8, 10 /
      data ratio / 1.000001d0 /
      data idebug, iout / 0, 99 /
c
c print left and right states
      if (idebug.ne.0) then
        write(iout,*)'Msg (cmpmiddlemg) : left and right states'
        write(iout,40)r(1,j),r(2,j),v(1,j),v(2,j),p(1,j),p(2,j)
   40   format(' rho(L)=',1pe15.8,'   rho(R)=',1pe15.8,/,
     *         '   v(L)=',1pe15.8,'     v(R)=',1pe15.8,/,
     *         '   p(L)=',1pe15.8,'     p(R)=',1pe15.8)
      end if
c
c material type for phase j
      kj=kmat(j)
c
c get initial guess for middle densities
      call guessmg (j,kj,pstar,rm,iguess)
c      write(6,*)rm
c      read(5,*)idum
c
      if (idebug.ne.0) then
        if (iguess.eq.0) then
          write(iout,*)
     *        'Msg (cmpmiddlemg) : guess based on linearization'
        else
          if (iguess.gt.0) then
            write(iout,*)
     *        'Msg (cmpmiddlemg) : guess based on two-S approx'
          else
            write(iout,*)
     *        'Msg (cmpmiddlemg) : guess based on two-R approx'
          end if
        end if
      end if
c
c initialize Romberg integration
      do i=1,2
        aint0(i,j)=0.d0
        r0(i,j)=r(i,j)
        c0(i,j)=c(i,j)
      end do
c
c rarefaction/shock solution flags
      do i=1,2
        initr(i,j)=0
        inits(i,j)=0
      end do
c
c density scale
      rscale=.5d0*(r(1,j)+r(2,j))

c
c Newton iteration to find middle states
      do it=1,itmax
c
c evaluate rarefaction or shock from each side
        do i=1,2
c
          if (rm(i)/r(i,j).le.ratio) then
c
c rarefaction solution
            ishock(i)=0
            call raremg (i,j,kj,rm(i),dv(i),pm(i),fp(i),gp(i),sm(i),
     *                   ier)
            if (ier.ne.0) then
              if (it.eq.1) then
                write(6,*)'Error (cmpmiddlemg) : error from raremg'
                stop
              else
                if (mgfxdbg) then
                  write(6,*)'Warning (cmpmiddlemg) : error from raremg'
                end if
                do ii=1,2
                  rm(ii)=rmsave(ii)
                  vm(ii)=vmsave(ii)
                  pm(ii)=pmsave(ii)
                  sm(ii)=smsave(ii)
                end do
                mgfxerr=1
                return
              end if
            end if
c            write(6,*)'r=',i,j,kj,rm(i),dv(i),pm(i),fp(i),gp(i)
c            read(5,*)idum
c
c            delta=1.d-7
c            delta1=delta*rm(i)
c            rm1=rm(i)+delta1
c            call raremg (i,j,kj,rm1,dv1,pm1,fp1,gp1,sm1,ier)
c            fpi=(dv1-dv(i))/delta1
c            gpi=(pm1-pm(i))/delta1
c            fperr=dabs(fp(i)-fpi)/dabs(fpi)
c            gperr=dabs(gp(i)-gpi)/dabs(gpi)
c            write(6,*)'fp=',fp(i),fpi,fperr
c            write(6,*)'gp=',gp(i),gpi,gperr
c            fp(i)=fpi
c            gp(i)=gpi
c
          else
c
c shock solution
            ishock(i)=1
            call shockmg (i,j,kj,rm(i),dv(i),pm(i),fp(i),gp(i),sm(i),
     *                    ier)
            if (ier.ne.0) then
              if (it.eq.1) then
                write(6,*)'Error (cmpmiddlemg) : error from shockmg'
                stop
              else
                if (mgfxdbg) then
                 write(6,*)'Warning (cmpmiddlemg) : error from shockmg'
                end if
                do ii=1,2
                  rm(ii)=rmsave(ii)
                  vm(ii)=vmsave(ii)
                  pm(ii)=pmsave(ii)
                  sm(ii)=smsave(ii)
                end do
                mgfxerr=1
                return
              end if
            end if
c            write(6,*)'s=',i,j,kj,rm(i),dv(i),pm(i),fp(i),gp(i)
c            read(5,*)idum
c
c            delta=1.d-7
c            delta1=delta*rm(i)
c            rm1=rm(i)+delta1
c            call shockmg (i,j,kj,rm1,dv1,pm1,fp1,gp1,sm1,ier)
c            fpi=(dv1-dv(i))/delta1
c            gpi=(pm1-pm(i))/delta1
c            fperr=dabs(fp(i)-fpi)/dabs(fpi)
c            gperr=dabs(gp(i)-gpi)/dabs(gpi)
c            write(6,*)'fp=',fp(i),fpi,fperr
c            write(6,*)'gp=',gp(i),gpi,gperr
c            fp(i)=fpi
c            gp(i)=gpi
c
          end if
c
        end do
c
c compute middle velocities
        vm(1)=v(1,j)-dv(1)
        vm(2)=v(2,j)+dv(2)
c        write(6,*)'i=1',v(1,j),dv(1),vm(1)
c        write(6,*)'i=2',v(2,j),dv(2),vm(2)
c
        if (it.eq.1) then
          do i=1,2
            rmsave(i)=rm(i)
            smsave(i)=sm(i)
          end do
          pmsave(1)=.5d0*(pm(1)+pm(2))
          vmsave(1)=.5d0*(vm(1)+vm(2))
          pmsave(2)=pmsave(1)
          vmsave(2)=vmsave(1)
        end if
c
c compute density increments
        det=fp(1)*gp(2)+fp(2)*gp(1)
        drm(1)=(gp(2)*(vm(2)-vm(1))-fp(2)*(pm(2)-pm(1)))/det
        drm(2)=(gp(1)*(vm(2)-vm(1))+fp(1)*(pm(2)-pm(1)))/det
c
c update middle densities
        rm(1)=rm(1)-drm(1)
        rm(2)=rm(2)-drm(2)
c
        if (rm(1).le.0.d0.or.rm(2).le.0.d0) then
          if (mgfxdbg)
     *        write(6,*)'Warning (cmpmiddlemg) : negative densities'
          do i=1,2
            rm(i)=rmsave(i)
            vm(i)=vmsave(i)
            pm(i)=pmsave(i)
            sm(i)=smsave(i)
          end do
          mgfxerr=1
          return
        end if
c
        err=max(dabs(drm(1)),dabs(drm(2)))/rscale
c
c label transition type, rarefaction (R) or shock (S)
        do i=1,2
          if (ishock(i).eq.0) then
            write(itype(i),50)i
   50       format(i1,'-R')
          else
            write(itype(i),51)i
   51       format(i1,'-S')
          end if
        end do
c
c print convergence information
        if (idebug.ne.0) then
          write(iout,100)it,(itype(i),rm(i),drm(i),i=1,2),err
  100     format(' * it=',i2,' : ',a3,', rm(1)=',1pe15.8,
     *           ', drm(1)=',1pe9.2,/,
     *               '           ',a3,', rm(2)=',1pe15.8,
     *           ', drm(2)=',1pe9.2,', err=',1pe9.2)
        end if
c
        if (err.lt.tol) then
          pm(1)=.5d0*(pm(1)+pm(2))
          vm(1)=.5d0*(vm(1)+vm(2))
          pm(2)=pm(1)
          vm(2)=vm(1)
          if (idebug.ne.0) then
            write(iout,200)pm(1),vm(1)
  200       format(' * converged : pstar=',1pe15.8,', vstar=',1pe15.8)
          end if
          return
        end if
c
      end do
c
      if (mgfxdbg) then
       write(6,*)'Warning (cmpmiddlemg) : iteration failed to converge'
      end if
      do i=1,2
        rm(i)=rmsave(i)
        vm(i)=vmsave(i)
        pm(i)=pmsave(i)
        sm(i)=smsave(i)
      end do
      mgfxerr=1
c
      return
      end
c
c+++++++++++++++++
c
      subroutine guessmg (j,kj,pstar,rm,iguess)
c
c compute initial guess for middle densities based on a local
c stiffened gas approximation
c
      implicit real*8 (a-h,o-z)
      dimension rm(2)
      dimension a(2),b(2)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
      data pratio, pfact / 2.d0, 1.d-4 /
c
c estimates for gamma, gm1 and p0 based on averages
      gam(j)=.5d0*(gamk(1,j)+gamk(2,j))
      gm1(j)=gam(j)-1.d0
      p0(j)=.5d0*(p0k(1,j)+p0k(2,j))
c
c compute the average density and sound speed
      rma=.5d0*(r(1,j)+r(2,j))
      cma=.5d0*(c(1,j)+c(2,j))
c
c check for vacuum state (in the decoupled phases)
      if (2.d0*(c(1,j)+c(2,j))/gm1(j).le.v(2,j)-v(1,j)) then
        write(6,*)'Error (guessmg) : vacuum found, j=',j
        stop
      end if
c
c compute min/max pressures
      pmin=min(p(1,j)+p0(j),p(2,j)+p0(j))
      pmax=max(p(1,j)+p0(j),p(2,j)+p0(j))
c
c start with guess based on a linearization
      ppv=.5d0*(p(1,j)+p(2,j))
     *    -.125d0*(v(2,j)-v(1,j))*(r(1,j)+r(2,j))*(c(1,j)+c(2,j))
      ppv=max(ppv+p0(j),0.d0)
c
      if (pmax/pmin.le.pratio
     *   .and.pmin.le.ppv.and.pmax.ge.ppv) then
        pstar=ppv-p0(j)
        rm(1)=r(1,j)-(rma*cma*(v(2,j)-v(1,j))-p(2,j)+p(1,j))
     *               /(2.d0*cma**2)
        rm(2)=r(2,j)-(rma*cma*(v(2,j)-v(1,j))-p(1,j)+p(2,j))
     *               /(2.d0*cma**2)
        if (rm(1).le.0.d0.or.rm(2).le.0.d0) then
          write(6,*)'Error (guessmg) : initial density is negative'
          stop
        end if
        iguess=0
      else
        if (ppv.lt.pmin) then
c guess based on two rarefaction solution
          if (r(1,j).lt.r(2,j)) then
            gam(j)=gamk(1,j)
            p0(j)=p0k(1,j)
          else
            gam(j)=gamk(2,j)
            p0(j)=p0k(2,j)
          end if
          gm1(j)=gam(j)-1.d0
          em=0.5d0*gm1(j)/gam(j)
          arg1=c(1,j)/((p(1,j)+p0(j))**em)
          arg2=c(2,j)/((p(2,j)+p0(j))**em)
          arg3=(c(1,j)+c(2,j)-.5d0*gm1(j)*(v(2,j)-v(1,j)))/(arg1+arg2)
          pstar=arg3**(1.d0/em)-p0(j)
          rm(1)=r(1,j)*((pstar+p0(j))/(p(1,j)+p0(j)))**(1.d0/gam(j))
          rm(2)=r(2,j)*((pstar+p0(j))/(p(2,j)+p0(j)))**(1.d0/gam(j))
          iguess=-1
        else
c guess based on two shock approximate solution
          if (r(1,j).gt.r(2,j)) then
            gam(j)=gamk(1,j)
            p0(j)=p0k(1,j)
          else
            gam(j)=gamk(2,j)
            p0(j)=p0k(2,j)
          end if
          gm1(j)=gam(j)-1.d0
          gp1=gam(j)+1.d0
          do i=1,2
            a(i)=2.d0/(gp1*r(i,j))
            b(i)=gm1(j)*(p(i,j)+p0(j))/gp1+p0(j)
          end do
          gl=dsqrt(a(1)/(ppv+b(1)))
          gr=dsqrt(a(2)/(ppv+b(2)))
          pts=(gl*p(1,j)+gr*p(2,j)-v(2,j)+v(1,j))/(gl+gr)
          ptol=pfact*min(p(1,j),p(2,j))
          pstar=max(ptol,pts)
          do i=1,2
            fact=.5d0*(pstar+p(i,j))
            anum=fact+(pstar+gam(j)*p0(j))/gm1(j)
            denom=fact+(p(i,j)+gam(j)*p0(j))/gm1(j)
            rm(i)=r(i,j)*anum/denom
          end do
          iguess=1
        end if
      end if
c
      return
      end
c
c+++++++++++++++++++
c
      subroutine raremg (i,j,kj,rm,dv,pm,fp,gp,sm,ier)
c
c Rarefaction solution
c
      implicit real*8 (a-h,o-z)
      logical mgfxdbg
      common / mgflxerr / mgfxerr,mgfxdbg
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / rardatmg / amumg1(2,2),Hmg1(2,2),initr(2,2)
c      common / xtdat / xpos,tpos,lv
c
      ier=0
c      
c      if (rm.le.0.d0) then
c        if (mgfxdbg) write(6,*)'Error (raremg) : rm.le.0'
c        if (mgfxdbg) write(6,*)'x=',xpos,'t=',tpos,'lv=',lv
c        ier=1
c        return
c      end if
c
c set up and compute the pressure at rm along the isentrope passing through (r(i,j),p(i,j))
      if (initr(i,j).eq.0) then
        Vol1=1.d0/r(i,j)
        amumg1(i,j)=amumg(kj,Vol1)
        Hmg1(i,j)=HmgV(j,Vol1,ier)
        if (ier.ne.0) then
          if (mgfxdbg) write(6,*)'Error (raremg) : error computing Hmg1'
          return
        end if
        initr(i,j)=1
      end if
      Volm=1.d0/rm
      HmgVm=HmgV(j,Volm,ier)
      if (ier.ne.0) then
        if (mgfxdbg) write(6,*)'Error (raremg) : error computing HmgVm'
        return
      end if
c      write(6,*)Vol1,amumg1(i,j),Hmg1(i,j),HmgVm
      pm=(p(i,j)*amumg1(i,j)-(HmgVm-Hmg1(i,j)))/amumg(kj,Volm)
c
c compute the sound speed at rm,pm
      cm2=c2mg(kj,Volm,pm)
      if (cm2.lt.0.d0) then
        if (mgfxdbg) write(6,*)'Error (raremg) : cm2.lt.0, i,j=',i,j
        ier=2
        return
      end if
      cm=dsqrt(cm2)
c
c compute the change in velocity by integrating along a C_pm characteristic
      dv=getdvmg(i,j,kj,rm,cm,amumg1(i,j),Hmg1(i,j))
c
c derivatives for Newton's method
      gp=cm2
      fp=cm/rm
c
c wave speed is the sound speed for rarefaction solution
      sm=cm
c
      return
      end
c
c+++++++++++++++++++
c
      subroutine shockmg (i,j,kj,rm,dv,pm,fp,gp,sm,ier)
c
c Shock solution
c
      implicit real*8 (a-h,o-z)
      logical mgfxdbg
      common / mgflxerr / mgfxerr,mgfxdbg
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / shkdatmg / Amg1(2,2),Bmg1(2,2),inits(2,2)
c
      ier=0
c
      if (rm.le.0.d0) then
        if (mgfxdbg) write(6,*)'Error (shockmg) : rm.le.0'
        ier=1
        return
      end if
c
c set up and compute the pressure at rm along the Hugoniot passing through (r(i,j),p(i,j))
      if (inits(i,j).eq.0) then
        Vol1=1.d0/r(i,j)
        Amg1(i,j)=AmgV(kj,Vol1)
        Bmg1(i,j)=BmgV(j,Vol1,ier)
        if (ier.ne.0) then
          if (mgfxdbg)
     *        write(6,*)'Error (shockmg) : error computing Bmg1'
          return
        end if
        inits(i,j)=1
      end if
      Volm=1.d0/rm
      AmgVm=AmgV(kj,Volm)
      BmgVm=BmgV(j,Volm,ier)
      if (ier.ne.0) then
        if (mgfxdbg)
     *      write(6,*)'Error (shockmg) : error computing BmgVm'
        return
      end if
      z=r(i,j)*Volm
      fact=.5d0*(1.d0-z)
      denom=r(i,j)*AmgVm-fact
      if (denom.le.0.d0) then
        ier=2
        if (mgfxdbg)
     *      write(6,*)'Error (shockmg) : exceeded compression limit'
        return
      end if
      pm=(p(i,j)*(r(i,j)*Amg1(i,j)+fact)
     *         -r(i,j)*(BmgVm-Bmg1(i,j)))/denom
c
c compute the change in velocity using the Rayleigh line
      if (pm.le.p(i,j)) then
        ier=3
        if (mgfxdbg) then
          write(6,*)'Error (shockmg) : expansion shock'
          write(6,*)'rm,r(i,j) =',rm,r(i,j)
          write(6,*)'pm,p(i,j) =',pm,p(i,j)
        end if
        return
      end if
      fact=pm/p(i,j)-1.d0
      dv=dsqrt(p(i,j)*fact*(1.d0-z)/r(i,j))
c
c derivatives for Newton's method
      gp=(c2mg(kj,Volm,pm)*AmgVm*r(i,j)
     *       -.5d0*p(i,j)*fact*z**2/r(i,j))/denom
      fp=.5d0*dv*(gp/(pm-p(i,j))
     *       +Volm**2/(1.d0/r(i,j)-Volm))
c
c wave speed is -(U-vL) or (U-vR), where U is the shock velocity
c (Note: this divided difference is a valid calculation, but may be
c  inaccurate in the weak shock limit when z->1 where sm->cm.)
      sm=dv/(1.d0-z)
c
      return
      end
c
c+++++++++++++++++++
c
      double precision function getdvmg (i,j,kj,rm,cm,amumg1,Hmg1)
c
c evaluate FL(rho) or FR(rho) by integrating the C+ or C- characteristic equation
c numerically using a Romberg integration scheme (expensive).
c
      implicit real*8 (a-h,o-z)
      parameter (lvmax=20)
      dimension aint(2,lvmax)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / romdatmg / aint0(2,2),r0(2,2),c0(2,2),tolrom(2),lvrom(2)
c
c initialize Romberg integration
      dr=rm-r0(i,j)
      aint(1,1)=.5d0*(cm/rm+c0(i,j)/r0(i,j))*dr
c
c check whether to use a trapezoidal rule approximation
      if (lvrom(j).le.1) then
        aint0(i,j)=aint(1,1)+aint0(i,j)
        r0(i,j)=rm
        c0(i,j)=cm
        getdvmg=aint0(i,j)
        return
      end if
c
c compute integral using Romberg for n.le.lvrom, stop if the error is
c less than tolrom.  (lvrom.eq.2 is Simpson's rule, for example)
      nmax=min(lvrom(j),lvmax)
      do n=2,nmax
        nm1=n-1
        sum=0.d0
        kmax=2**(n-2)
        do k=1,kmax
          rk=r0(i,j)+(k-.5d0)*dr
          Vk=1.d0/rk
          HmgVk=HmgV(j,Vk,ier)
          if (ier.ne.0) then
            write(6,*)'Error (getdvmg) : error computing HmgVk'
            stop
          end if
          pk=(p(i,j)*amumg1-(HmgVk-Hmg1))/amumg(kj,Vk)
          ck2=c2mg(kj,Vk,pk)
          if (ck2.lt.0.d0) then
            write(6,*)'Error (getdvmg) : ck2.lt.0, i=',i
            stop
          end if
          sum=sum+dsqrt(ck2)/rk
        end do
        aint(2,1)=.5d0*(aint(1,1)+dr*sum)
        do l=1,nm1
          aint(2,l+1)=aint(2,l)+(aint(2,l)-aint(1,l))/(4**l-1.d0)
        end do
        err=dabs(aint(2,n)-aint(1,nm1))/c(i,j)
        if (err.lt.tolrom(j).or.n.ge.nmax) then
          aint0(i,j)=aint(2,n)+aint0(i,j)
          r0(i,j)=rm
          c0(i,j)=cm
          getdvmg=aint0(i,j)
          if (err.gt.tolrom(j)) then
            write(6,*)'Msg (getdvmg) : Rel_error > Romberg tolerance'
          end if
          return
        end if
        dr=.5d0*dr
        do l=1,n
          aint(1,l)=aint(2,l)
        end do
      end do
c
      write(6,*)'Error (getdvmg) : oops, bad logic'
c
      stop
      end
c
c+++++++++++++++
c
      subroutine cmpcouplemg (rm0,rm,vm,pm,fp,gp,sm,linsc)
c
c Use Newton's method to compute coupled middle states
c
      implicit real*8 (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),fp(2,2),gp(2,2),sm(2,2)
      dimension a(5,6),drm(2,2)
      dimension rmsave(2,2),vmsave(2,2),pmsave(2,2),smsave(2,2)
      character*3 itype(2,2)
      logical mgfxdbg
      common / mgflxerr / mgfxerr,mgfxdbg
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / matdatmg / kmat(2)
c      common / flxdat / afxcnt(3),nfxcnt(3)
c      common / xtdat / xpos,tpos,lv
      data atol, tol, itmax / 1.d-14, 1.d-8, 10 /
      data sig, frac / 10.d0, 0.99d0 /
      data c2min / 1.d-8 /
c
c print left and right volume fractions
c      write(6,*)'Msg (cmpcouplemg) : left and right volume fractions'
c      write(6,40)alpha(1,1),alpha(2,1)
c   40 format(' alpha_s(L)=',1pe15.8,'   alpha_s(R)=',1pe15.8)
c
c      nfxcnt(1)=nfxcnt(1)+1
c
      rmmin=rm(1,1)
      do i=1,2
        do j=1,2
          rmmin=min(rm(i,j),rmmin)
        end do
      end do
c
      if (rmmin.le.0.d0) then
        write(6,*)'Error (cmpcouplemg) : rm.le.0'
        write(6,101)rm(1,1),rm(2,1),rm(1,2),rm(2,2)
  101    format('rm=',4(1x,1pe15.8))
        stop
      end if
c
c initialize rm0
      call setrm0mg (rm0,rm,vm,iconf)
c
c check whether volume fractions on the left and right are nearly the same
      if (dabs(alpha(1,1)-alpha(2,1)).lt.atol) return
c
c start with linearized solid contact, compute velocity difference
      dvm=vm(1,2)-vm(1,1)
c
c density scales
      rscale1=max(rm(1,1),rm(2,1))
      rscale2=max(rm(1,2),rm(2,2))
c
c compute alpha-bar_m and alpha_m
      if (alpha(1,1).gt.0.5d0) then
        if (alpha(2,1).lt.0.5d0) then
          asm=.5d0
        else
          asm=min(alpha(1,1),alpha(2,1))
        end if
      else
        if (alpha(2,1).gt.0.5d0) then
          asm=.5d0
        else
          asm=max(alpha(1,1),alpha(2,1))
        end if
      end if
      agm=1.d0-asm
c
c linear update for middle solid densities
      eps=alpha(2,1)-alpha(1,1)
      del=(pm(1,2)-pm(1,1))*eps
     *       /(asm*(fp(1,1)*gp(2,1)+fp(2,1)*gp(1,1)))
      drm1=-fp(2,1)*del
      drm2= fp(1,1)*del
      rm(1,1)=rm(1,1)+drm1
      rm(2,1)=rm(2,1)+drm2
      err1=max(dabs(drm1),dabs(drm2))/rscale1
c
c limit velocity difference
      Vol0=1.d0/rm0
      cm2=max(c2mg(kmat(2),Vol0,pm(1,2)),c2min)
      gampg=rm0*cm2
      cm=frac*dsqrt(cm2)
      arg1=sig*(1.d0+dvm/cm)
      arg2=sig*(1.d0-dvm/cm)
      dvm=.5d0*cm*dlog(dcosh(arg1)/dcosh(arg2))/sig
c
c linear update for middle gas densities
      del=gampg*dvm*eps/(agm*(gampg-rm0*dvm**2)
     *           *(fp(1,2)*gp(2,2)+fp(2,2)*gp(1,2)))
      drm1=(gp(2,2)+rm0*dvm*fp(2,2))*del
      drm2=(gp(1,2)-rm0*dvm*fp(1,2))*del
      rm(1,2)=rm(1,2)+drm1
      rm(2,2)=rm(2,2)+drm2
      if (iconf.eq.1) then
        drm0=drm1+(gp(2,2)*drm2-gp(1,2)*drm1)/cm2
      else
        drm0=drm2-(gp(2,2)*drm2-gp(1,2)*drm1)/cm2
      end if
      rm0=rm0+drm0
      err2=max(dabs(drm1),dabs(drm2),dabs(drm0))/rscale2
c
c update middle velocities and pressures
      call updatemg (rm,vm,pm,fp,gp,sm,ier)
c
c check sonic condition on the middle states
      it=0
      isonic=0
      call chkvelmg (it,rm(1,2),vm,pm(1,2),isonic)
c
c      if (isonic.ne.0) then
c        if (mgfxdbg) then
c          write(6,*)'Warning (cmpcouplemg) : error from chkvelmg'
c          write(6,*)'x=',xpos,'t=',tpos,'lv=',lv
c        end if
c        mgfxerr=1
c        return
c      end if
c
c check errors
      if (max(err1,err2).lt.tol.or.linsc.ne.0.or.isonic.ne.0) then
c        nfxcnt(2)=nfxcnt(2)+1
        return
      end if
c
      err1save=err1
      err2save=err2
c
c reset flux error flag
      mgfxerr=0
c
      err1b=1.d10
      err2b=1.d10
      rm0save=rm0
      do i=1,2
        do j=1,2
          rmsave(i,j)=rm(i,j)
          vmsave(i,j)=vm(i,j)
          pmsave(i,j)=pm(i,j)
          smsave(i,j)=sm(i,j)
        end do
      end do
c
c Newton iteration to find rm0, rm(i,j)
      do it=1,itmax
c
c get Jacobian matrix (using finite differences)
c        call jacobmgfd (rm0,rm,vm,pm,a,iconf)
c
c get Jacobian matrix (exact)
        call jacobmg (rm0,rm,vm,pm,fp,gp,a,iconf)
c
c compute corrections to the middle densities
        call solvemg (a,drm,drm0)
c
c update middle densities
        rmmin=rm(1,1)-drm(1,1)
        do i=1,2
          do j=1,2
            rm(i,j)=rm(i,j)-drm(i,j)
            rmmin=min(rm(i,j),rmmin)
          end do
        end do
        rm0=rm0-drm0
        rmmin=min(rm0,rmmin)
c
        if (rmmin.le.0.d0) then
c          write(6,*)'Error (cmpcouplemg) : rm.le.0'
c          write(6,201)rm(1,1),rm(2,1),rm(1,2),rm(2,2)
c  201     format('rm=',4(1x,1pe15.8))
c          write(6,202)drm(1,1),drm(2,1),drm(1,2),drm(2,2)
c  202     format('drm=',4(1x,1pe15.8))
c          stop
          if (mgfxdbg) then
            write(6,*)'Warning (cmpcouplemg) : rm.le.0'
c            write(6,*)'x=',xpos,'t=',tpos,'lv=',lv
            write(6,*)'ul, ur ='
            write(6,*)alpha(1,1),alpha(2,1)
            write(6,*)r(1,1),r(2,1)
            write(6,*)v(1,1),v(2,1)
            write(6,*)p(1,1),p(2,1)
            write(6,*)r(1,2),r(2,2)
            write(6,*)v(1,2),v(2,2)
            write(6,*)p(1,2),p(2,2)
            write(6,*)'err1save, err2save=',err1save,err2save
          end if
          rm0=rm0save
          do i=1,2
            do j=1,2
              rm(i,j)=rmsave(i,j)
              vm(i,j)=vmsave(i,j)
              pm(i,j)=pmsave(i,j)
              sm(i,j)=smsave(i,j)
            end do
          end do
c          mgfxerr=1
          return
        end if
c
c update middle velocities and pressures
        call updatemg (rm,vm,pm,fp,gp,sm,ier)
c
        if (ier.ne.0) then
          if (mgfxdbg) then
            write(6,*)'Warning (cmpcouplemg) : error in update'
c            write(6,*)'x=',xpos,'t=',tpos,'lv=',lv
            write(6,*)'ul, ur ='
            write(6,*)alpha(1,1),alpha(2,1)
            write(6,*)r(1,1),r(2,1)
            write(6,*)v(1,1),v(2,1)
            write(6,*)p(1,1),p(2,1)
            write(6,*)r(1,2),r(2,2)
            write(6,*)v(1,2),v(2,2)
            write(6,*)p(1,2),p(2,2)
            write(6,*)'err1save, err2save=',err1save,err2save
          end if
          rm0=rm0save
          do i=1,2
            do j=1,2
              rm(i,j)=rmsave(i,j)
              vm(i,j)=vmsave(i,j)
              pm(i,j)=pmsave(i,j)
              sm(i,j)=smsave(i,j)
            end do
          end do
c          mgfxerr=1
          return
        end if
c
c check sonic condition on the middle states
        isonic=0
        call chkvelmg (it,rm(1,2),vm,pm(1,2),isonic)
c
        if (isonic.ne.0) then
          if (mgfxdbg) then
            write(6,*)'Warning (cmpcouplemg) : error from chkvelmg'
c            write(6,*)'x=',xpos,'t=',tpos,'lv=',lv
            write(6,*)'ul, ur ='
            write(6,*)alpha(1,1),alpha(2,1)
            write(6,*)r(1,1),r(2,1)
            write(6,*)v(1,1),v(2,1)
            write(6,*)p(1,1),p(2,1)
            write(6,*)r(1,2),r(2,2)
            write(6,*)v(1,2),v(2,2)
            write(6,*)p(1,2),p(2,2)
            write(6,*)'err1save, err2save=',err1save,err2save
          end if
          rm0=rm0save
          do i=1,2
            do j=1,2
              rm(i,j)=rmsave(i,j)
              vm(i,j)=vmsave(i,j)
              pm(i,j)=pmsave(i,j)
              sm(i,j)=smsave(i,j)
            end do
          end do
c          mgfxerr=1
          return
        end if
c
c label transition type, rarefaction (R) or shock (S)
c        do j=1,2
c          do i=1,2
c            if (rm(i,j).le.r(i,j)) then
c              write(itype(i,j),50)i
c   50         format(i1,'-R')
c            else
c              write(itype(i,j),51)i
c   51         format(i1,'-S')
c            end if
c          end do
c        end do
c
c error in solid middle densities
        rscale1=max(rm(1,1),rm(2,1))
        err1=max(dabs(drm(1,1)),dabs(drm(2,1)))/rscale1
c
c error in gas middle densities
        rscale2=max(rm(1,2),rm(2,2),rm0)
        err2=max(dabs(drm(1,2)),dabs(drm(2,2)),dabs(drm0))/rscale2
c
c print convergence information
c        write(6,100)it,itype(1,1),rm(1,1),drm(1,1),
c     *                 itype(2,1),rm(2,1),drm(2,1),err1,
c     *                 itype(1,2),rm(1,2),drm(1,2),
c     *                 itype(2,2),rm(2,2),drm(2,2),
c     *                            rm0,    drm0,    err2
c  100   format(' * it=',i2,', solid : ',
c     *             a3,', rm1=',1pe15.8,', drm1=',1pe9.2,/,
c     *         18x,a3,', rm2=',1pe15.8,', drm2=',1pe9.2,
c     *                                 ', err(rm_sol)=',1pe9.2,/,
c     *         12x,'gas : ',
c     *             a3,', rm1=',1pe15.8,', drm1=',1pe9.2,/,
c     *         18x,a3,', rm2=',1pe15.8,', drm2=',1pe9.2,/,
c     *            21x,'  rm0=',1pe15.8,', drm0=',1pe9.2,
c     *                                 ', err(rm_gas)=',1pe9.2)
c
c check for convergence
        if (max(err1,err2).lt.tol) then
c        if (max(err1,err2).lt.tol.or.linsc.ne.0) then
c          write(6,101)(vm(1,j),pm(1,j),vm(2,j),pm(2,j),j=1,2)
c  101     format(' * converged, solid : vm1=',1pe15.8,', pm1=',1pe15.8,/,
c     *           '                      vm2=',1pe15.8,', pm2=',1pe15.8,/,
c     *           '                gas : vm1=',1pe15.8,', pm1=',1pe15.8,/,
c     *           '                      vm2=',1pe15.8,', pm2=',1pe15.8)
c          nfxcnt(3)=nfxcnt(3)+1
          return
c        else
c           if (err1.lt.err1b.and.err2.lt.err2b.and.1.eq.2) then
c              err1b=err1
c              err2b=err2
c              rm0save=rm0
c              do i=1,2
c                 do j=1,2
c                    rmsave(i,j)=rm(i,j)
c                    vmsave(i,j)=vm(i,j)
c                    pmsave(i,j)=pm(i,j)
c                    smsave(i,j)=sm(i,j)
c                 end do
c              end do
c           end if
        end if
c
      end do
c
      if (mgfxdbg) then
        write(6,*)'Warning (cmpcouplemg) : itmax exceeded'
c        write(6,*)'x=',xpos,'t=',tpos,'lv=',lv
        write(6,*)'ul, ur ='
        write(6,*)alpha(1,1),alpha(2,1)
        write(6,*)r(1,1),r(2,1)
        write(6,*)v(1,1),v(2,1)
        write(6,*)p(1,1),p(2,1)
        write(6,*)r(1,2),r(2,2)
        write(6,*)v(1,2),v(2,2)
        write(6,*)p(1,2),p(2,2)
        write(6,*)'err1save, err2save=',err1save,err2save
      end if
      rm0=rm0save
      do i=1,2
        do j=1,2
          rm(i,j)=rmsave(i,j)
          vm(i,j)=vmsave(i,j)
          pm(i,j)=pmsave(i,j)
          sm(i,j)=smsave(i,j)
        end do
      end do
c      mgfxerr=1
c
      return
      end
c
c+++++++++++++++
c
      subroutine chkvelmg (it,r,v,p,isonic)
c
c Check sonic condition
c
      implicit real*8 (a-h,o-z)
      dimension r(2),v(2,2),p(2)
      common / matdatmg / kmat(2)
      data c2min, small, amax / 1.d-8, 1.d-12, 0.98d0 /
c
c loop over sides
      do i=1,2
        Vol=1.d0/r(i)
        c2=max(c2mg(kmat(2),Vol,p(i)),c2min)
        amach=dsqrt(max((v(i,2)-v(i,1))**2/c2,small))
        if (amach.gt.amax) then
c          write(6,*)'Warning (chkvelmg) : it,i,amach =',it,i,amach
          isonic=1
        end if
      end do
c
      return
      end
c
c+++++++++++++++
c
      subroutine setrm0mg (rm0,rm,vm,iconf)
c
c Initialize the value of the middle gas density rm0
c
      implicit real*8 (a-h,o-z)
      dimension rm(2,2),vm(2,2)
c
c check whether gas contact is left/right of solid contact
      if (vm(1,1).lt.vm(1,2)) then
c gas contact is right of solid contact (this is the case in figure 2 of the 2006 paper)
        iconf=1
        rm0=rm(1,2)
      else
c gas contact is right of solid contact
        iconf=2
        rm0=rm(2,2)
      end if
c
      return
      end
c
c+++++++++++++++
c
      subroutine residmg (rm0,rm,vm,pm,b,iconf)
c
c compute residual of the solid contact jump conditions
c
      implicit real*8 (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),b(5)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                  gamk(2,2),p0k(2,2)
      common / matdatmg / kmat(2)
c
c check whether gas contact is left/right of solid contact
      if (iconf.eq.1) then
c gas contact is right of solid contact
        rm1=rm(1,2)
        rm2=rm0
      else
c gas contact is left of solid contact
        rm1=rm0
        rm2=rm(2,2)
      end if
c
c MG quantities for the gas at rm1 and rm2
      k2=kmat(2)
      Vol1=1.d0/rm1
      Amg1=AmgV(k2,Vol1)
      Bmg1=BmgV(2,Vol1,ier)
      if (ier.ne.0) then
        write(6,*)'Error (residmg) : error computing Bmg1'
        stop
      end if
      amumg1=amumg(k2,Vol1)
      Hmg1=HmgV(2,Vol1,ier)
      if (ier.ne.0) then
        write(6,*)'Error (residmg) : error computing Hmg1'
        stop
      end if
      Vol2=1.d0/rm2
      Amg2=AmgV(k2,Vol2)
      Bmg2=BmgV(2,Vol2,ier)
      if (ier.ne.0) then
        write(6,*)'Error (residmg) : error computing Bmg2'
        stop
      end if
      amumg2=amumg(k2,Vol2)
      Hmg2=HmgV(2,Vol2,ier)
      if (ier.ne.0) then
        write(6,*)'Error (residmg) : error computing Hmg2'
        stop
      end if
c
c velocity differences
      dv1=vm(1,2)-vm(1,1)
      dv2=vm(2,2)-vm(2,1)
c
c residual of the solid contact jump conditions
      b(1)=vm(2,1)-vm(1,1)
      b(2)=alpha(2,2)*rm2*dv2-alpha(1,2)*rm1*dv1
      b(3)= alpha(2,1)*pm(2,1)+alpha(2,2)*(pm(2,2)+rm2*dv2**2)
     *     -alpha(1,1)*pm(1,1)-alpha(1,2)*(pm(1,2)+rm1*dv1**2)
      b(4)= (Amg2+Vol2)*pm(2,2)+Bmg2+.5d0*dv2**2
     *     -(Amg1+Vol1)*pm(1,2)-Bmg1-.5d0*dv1**2
      b(5)=amumg2*pm(2,2)+Hmg2-amumg1*pm(1,2)-Hmg1
c
      return
      end
c
c+++++++++++++++
c
      subroutine updatemg (rm,vm,pm,fp,gp,sm,ier)
c
c  update the velocity and pressure across acoustic fields
c
      implicit real*8 (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),fp(2,2),gp(2,2),sm(2,2),dv(2)
      logical mgfxdbg
      common / mgflxerr / mgfxerr,mgfxdbg
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                  gamk(2,2),p0k(2,2)
      common / matdatmg / kmat(2)
      data ratio / 1.0001d0 /
c
      ier=0
c
      do j=1,2
        kj=kmat(j)
c
        do i=1,2
          if (rm(i,j)/r(i,j).le.ratio) then
c rarefaction solution
            call raremg (i,j,kj,rm(i,j),dv(i),pm(i,j),
     *                   fp(i,j),gp(i,j),sm(i,j),ier)
            if (ier.ne.0) then
              if (mgfxdbg)
     *            write(6,*)'Error (updatemg) : error in raremg'
              return
            end if
          else
c shock solution
            call shockmg (i,j,kj,rm(i,j),dv(i),pm(i,j),
     *                    fp(i,j),gp(i,j),sm(i,j),ier)
            if (ier.ne.0) then
              if (mgfxdbg)
     *            write(6,*)'Error (updatemg) : error in shockmg'
              return
            end if
          end if
        end do
c
c compute middle velocities
        vm(1,j)=v(1,j)-dv(1)
        vm(2,j)=v(2,j)+dv(2)
c
      end do
c
      return
      end
c
c+++++++++++++++
c
      subroutine jacobmgfd (rm0,rm,vm,pm,a,iconf)
c
c compute approximate Jacobian matrix of the solid contact jump conditions
c using finite differences
c
      implicit real*8 (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),a(5,6)
      dimension fp(2,2),gp(2,2),sm(2,2),b(5),bp(5)
      data delta / 1.d-7 /
c
c get residual of the solid-contact jump conditions
      call residmg (rm0,rm,vm,pm,b,iconf)
c      bmax=dabs(b(1))
c      do n=2,5
c        bmax=max(dabs(b(n)),bmax)
c      end do
c      write(6,*)'bmax=',bmax
c
      do i=1,2
        do j=1,2
          n=2*(i-1)+j
          dr=rm(i,j)*delta
          rm(i,j)=rm(i,j)+dr
          call updatemg (rm,vm,pm,fp,gp,sm,ier)
          if (ier.ne.0) then
            write(6,*)'Error (jacobmgfd) : error from updatemg'
            stop
          end if
          call residmg (rm0,rm,vm,pm,bp,iconf)
          do m=1,5
            a(m,n)=(bp(m)-b(m))/dr
          end do
          rm(i,j)=rm(i,j)-dr
        end do
      end do
c
      dr=rm0*delta
      rm0=rm0+dr
      call updatemg (rm,vm,pm,fp,gp,sm,ier)
      if (ier.ne.0) then
        write(6,*)'Error (jacobmgfd) : error from updatemg'
        stop
      end if
      call residmg (rm0,rm,vm,pm,bp,iconf)
      do m=1,5
        a(m,5)=(bp(m)-b(m))/dr
      end do
      rm0=rm0-dr
c
      do m=1,5
        a(m,6)=b(m)
      end do
c
      return
      end
c
c+++++++++++++++
c
      subroutine jacobmg (rm0,rm,vm,pm,fp,gp,a,iconf)
c
c compute Jacobian matrix of the solid contact jump conditions
c (Requires the subroutine matlmg in materials.f)
c
      implicit real*8 (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),fp(2,2),gp(2,2),a(5,6)
      dimension z1(4),z1p(4),z2(4),z2p(4)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / matdatmg / kmat(2)
c
      k2=kmat(2)
c
c velocity differences
      dv1=vm(1,2)-vm(1,1)
      dv2=vm(2,2)-vm(2,1)
c
c build the Jacobian
      a(1,1)=fp(1,1)
      a(1,3)=fp(2,1)
      a(1,2)=0.d0
      a(1,4)=0.d0
      a(1,5)=0.d0
c
      if (iconf.eq.1) then
c gas contact is right of solid contact
        ideriv=1
        rm1=rm(1,2)
        call matlmg (k2,rm1,z1,z1p,ideriv)
        rm2=rm0
        call matlmg (k2,rm2,z2,z2p,ideriv)
c
        a(2,1)=-alpha(1,2)*rm1*fp(1,1)
        a(2,3)=-alpha(2,2)*rm2*fp(2,1)
        a(2,2)= alpha(1,2)*(rm1*fp(1,2)-dv1)
        a(2,4)= alpha(2,2)*rm2*fp(2,2)
        a(2,5)= alpha(2,2)*dv2
c
        a(3,1)=-alpha(1,1)*gp(1,1)-2.d0*alpha(1,2)*rm1*dv1*fp(1,1)
        a(3,3)= alpha(2,1)*gp(2,1)-2.d0*alpha(2,2)*rm2*dv2*fp(2,1)
        a(3,2)=-alpha(1,2)*(gp(1,2)-2.d0*rm1*dv1*fp(1,2))
        a(3,4)= alpha(2,2)*(gp(2,2)+2.d0*rm2*dv2*fp(2,2))
        a(3,2)= a(3,2)-alpha(1,2)*dv1**2
        a(3,5)=        alpha(2,2)*dv2**2
c
        a(4,1)=-dv1*fp(1,1)
        a(4,3)=-dv2*fp(2,1)
        a(4,2)= dv1*fp(1,2)-(z1(1)+1.d0/rm1)*gp(1,2)
     *                     -(z1p(1)-1.d0/rm1**2)*pm(1,2)-z1p(2)
        a(4,4)= dv2*fp(2,2)+(z2(1)+1.d0/rm2)*gp(2,2)
        a(4,5)=            +(z2p(1)-1.d0/rm2**2)*pm(2,2)+z2p(2)
c
        a(5,1)=0.d0
        a(5,3)=0.d0
        a(5,2)=-z1(3)*gp(1,2)-z1p(3)*pm(1,2)-z1p(4)
        a(5,4)= z2(3)*gp(2,2)
        a(5,5)=              +z2p(3)*pm(2,2)+z2p(4)
c
      else
c gas contact is left of solid contact
        ideriv=1
        rm1=rm0
        call matlmg (k2,rm1,z1,z1p,ideriv)
        rm2=rm(2,2)
        call matlmg (k2,rm2,z2,z2p,ideriv)
c
        a(2,1)=-alpha(1,2)*rm1*fp(1,1)
        a(2,3)=-alpha(2,2)*rm2*fp(2,1)
        a(2,2)= alpha(1,2)*rm1*fp(1,2)
        a(2,4)= alpha(2,2)*(rm2*fp(2,2)+dv2)
        a(2,5)=-alpha(1,2)*dv1
c
        a(3,1)=-alpha(1,1)*gp(1,1)-2.d0*alpha(1,2)*rm1*dv1*fp(1,1)
        a(3,3)= alpha(2,1)*gp(2,1)-2.d0*alpha(2,2)*rm2*dv2*fp(2,1)
        a(3,2)=-alpha(1,2)*(gp(1,2)-2.d0*rm1*dv1*fp(1,2))
        a(3,4)= alpha(2,2)*(gp(2,2)+2.d0*rm2*dv2*fp(2,2))
        a(3,4)= a(3,4)+alpha(2,2)*dv2**2
        a(3,5)=       -alpha(1,2)*dv1**2
c
        a(4,1)=-dv1*fp(1,1)
        a(4,3)=-dv2*fp(2,1)
        a(4,2)= dv1*fp(1,2)-(z1(1)+1.d0/rm1)*gp(1,2)
        a(4,4)= dv2*fp(2,2)+(z2(1)+1.d0/rm2)*gp(2,2)
     *                     +(z2p(1)-1.d0/rm2**2)*pm(2,2)+z2p(2)
        a(4,5)=            -(z1p(1)-1.d0/rm1**2)*pm(1,2)-z1p(2)
c
        a(5,1)=0.d0
        a(5,3)=0.d0
        a(5,2)=-z1(3)*gp(1,2)
        a(5,4)= z2(3)*gp(2,2)+z2p(3)*pm(2,2)+z2p(4)
        a(5,5)=              -z1p(3)*pm(1,2)-z1p(4)
c
      end if
c
c residual of the solid contact jump conditions
      a(1,6)= vm(2,1)-vm(1,1)
      a(2,6)= alpha(2,2)*rm2*dv2-alpha(1,2)*rm1*dv1
      a(3,6)= alpha(2,1)*pm(2,1)+alpha(2,2)*(pm(2,2)+rm2*dv2**2)
     *       -alpha(1,1)*pm(1,1)-alpha(1,2)*(pm(1,2)+rm1*dv1**2)
      a(4,6)= (z2(1)+1.d0/rm2)*pm(2,2)+z2(2)+.5d0*dv2**2
     *       -(z1(1)+1.d0/rm1)*pm(1,2)-z1(2)-.5d0*dv1**2
      a(5,6)= z2(3)*pm(2,2)+z2(4)-z1(3)*pm(1,2)-z1(4)
c
      return
      end
c
c+++++++++++++++
c
      subroutine solvemg (a,drm,drm0)
c
      implicit real*8 (a-h,o-z)
      dimension a(5,6),drm(2,2)
      data tol / 1.d-12 /
c
c Gaussian elimination with partial pivoting
      do k=1,4
        kpiv=k
        apiv=dabs(a(k,k))
        do i=k+1,5
          if (dabs(a(i,k)).gt.apiv) then
            kpiv=i
            apiv=dabs(a(i,k))
          end if
        end do
        if (kpiv.ne.k) then
          do j=k,6
            atmp=a(k,j)
            a(k,j)=a(kpiv,j)
            a(kpiv,j)=atmp
          end do
        end if
        if (dabs(a(k,k)).lt.tol) then
          write(6,*)'Error (solvemg) : matrix singular??'
          stop
        end if
        do i=k+1,5
          fact=a(i,k)/a(k,k)
          do j=k+1,6
            a(i,j)=a(i,j)-fact*a(k,j)
          end do
        end do
      end do
c
c backward substitution
      a(5,6)=a(5,6)/a(5,5)
      do i=4,1,-1
        sum=a(i,6)
        do j=i+1,5
          sum=sum-a(i,j)*a(j,6)
        end do
        a(i,6)=sum/a(i,i)
      end do
c
c assign corrections
      drm(1,1)=a(1,6)
      drm(1,2)=a(2,6)
      drm(2,1)=a(3,6)
      drm(2,2)=a(4,6)
      drm0=a(5,6)
c
      return
      end
c
c+++++++++++++++
c
      subroutine cmpflxmg (m,aj,an,rm0,rm,vm,pm,sm,wl,wr,fl,fr,method)
c
c compute fluxes based on coupled middle states
c
      implicit real*8 (a-h,o-z)
      dimension an(2),rm(2,2),vm(2,2),pm(2,2),sm(2,2),
     *          wl(m),wr(m),fl(m),fr(m)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
c
c gas-state "zero" flag, default value
      izero=0
c
c solid contact source contributions
      source=pm(2,1)*alpha(2,1)-pm(1,1)*alpha(1,1)
      source2=vm(1,1)*source
c
c check if solid contact is to the right or left of x=0
      if (vm(1,1).gt.0.d0) then
c
c tangential component of the solid velocity from left state
        vt=an(1)*wl(3)-an(2)*wl(2)
c
c right-of-solid-contact state (only need solid velocity and gas pressure)
        wr(2)=an(1)*vm(1,1)-an(2)*vt
        wr(3)=an(2)*vm(1,1)+an(1)*vt
        wr(8)=pm(2,2)
c
c solid phase flux components
        call getfxmg (1,1,rm0,rm,vm,pm,sm,
     *                aj,an,vt,wl(1),fl(1),izero,method)
c
c add on solid contact source contributions
        fr(1)=fl(1)
        fr(2)=fl(2)+aj*an(1)*source
        fr(3)=fl(3)+aj*an(2)*source
        fr(4)=fl(4)+aj*source2
c
c gas phase flux components
        if (vm(1,2).ge.0.d0) then
c gas contact (and solid contact) to the right
c compute tangential component of the gas velocity from left state
          vt=an(1)*wl(7)-an(2)*wl(6)
          call getfxmg (1,2,rm0,rm,vm,pm,sm,
     *                  aj,an,vt,wl(5),fl(5),izero,method)
        else
c gas contact to the left, solid contact to the right
c compute tangential component of the gas velocity from right state
          izero=1
          vt=an(1)*wr(7)-an(2)*wr(6)
          call getfxmg (1,2,rm0,rm,vm,pm,sm,
     *                  aj,an,vt,wl(5),fl(5),izero,method)
        end if
c
c add on solid contact source contributions
        fr(5)=fl(5)
        fr(6)=fl(6)-aj*an(1)*source
        fr(7)=fl(7)-aj*an(2)*source
        fr(8)=fl(8)-aj*source2
c
c solid volume fraction flux
        fl(9)= 0.d0
        fr(9)=-aj*vm(1,1)*(alpha(2,1)-alpha(1,1))
c
c flux for advected variables
        do i=10,m
          fl(i)= 0.d0
          fr(i)=-aj*vm(1,1)*(wr(i)-wl(i))
        end do
c
      else
c
c tangential component of the solid velocity from right state
        vt=an(1)*wr(3)-an(2)*wr(2)
c
c left-of-solid-contact state (only need solid velocity and gas pressure)
        wl(2)=an(1)*vm(1,1)-an(2)*vt
        wl(3)=an(2)*vm(1,1)+an(1)*vt
        wl(8)=pm(1,2)
c
c solid phase flux components
        call getfxmg (2,1,rm0,rm,vm,pm,sm,
     *                aj,an,vt,wr(1),fr(1),izero,method)
c
c add on solid contact source contributions
        fl(1)=fr(1)
        fl(2)=fr(2)-aj*an(1)*source
        fl(3)=fr(3)-aj*an(2)*source
        fl(4)=fr(4)-aj*source2
c
c gas phase flux components
        if (vm(2,2).le.0.d0) then
c gas contact (and solid contact) to the left
c compute tangential component of the gas velocity from right state
          vt=an(1)*wr(7)-an(2)*wr(6)
          call getfxmg (2,2,rm0,rm,vm,pm,sm,
     *                  aj,an,vt,wr(5),fr(5),izero,method)
        else
c gas contact to the right, solid contact to the left
c compute tangential component of the gas velocity from left state
          izero=1
          vt=an(1)*wl(7)-an(2)*wl(6)
          call getfxmg (2,2,rm0,rm,vm,pm,sm,
     *                  aj,an,vt,wr(5),fr(5),izero,method)
        end if
c
c add on solid contact source contributions
        fl(5)=fr(5)
        fl(6)=fr(6)+aj*an(1)*source
        fl(7)=fr(7)+aj*an(2)*source
        fl(8)=fr(8)+aj*source2
c
c solid volume fraction flux
        fl(9)= aj*vm(1,1)*(alpha(2,1)-alpha(1,1))
        fr(9)= 0.d0
c
c flux for advected variables
        do i=10,m
          fl(i)= aj*vm(1,1)*(wr(i)-wl(i))
          fr(i)= 0.d0
        end do
c
      end if
c
      return
      end
c
c+++++++++++++++
c
      subroutine getfxmg (i,j,rm0,rm,vm,pm,sm,aj,an,vt,wstar,fx,
     *                    izero,method)
c
c compute solid or gas flux for i=side and j=phase
c
      implicit double precision (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),sm(2,2),an(2),wstar(4),fx(4)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
      common / sigdatmg / gp1(2),em(2),ep(2),a(2,2),b(2,2)
      common / matdatmg / kmat(2)
c
c      if (j.eq.2) then
c        write(6,*)'an(1),an(2),vm,vt=',an(1),an(2),vm(i,j),vt
c      end if
c
      if (izero.ne.0) then
c special case in which the star state of the gas lies between
c the contacts of the gas and solid phases
        wstar(1)=rm0
        wstar(2)=an(1)*vm(i,j)-an(2)*vt
        wstar(3)=an(2)*vm(i,j)+an(1)*vt
        wstar(4)=pm(i,j)
      else
c
c remaining cases: set isign=1 for i=1 and isign=-1 for i=2
        isign=3-2*i
c
        if (pm(i,j).gt.p(i,j)) then
c shock cases
          sp=isign*v(i,j)-sm(i,j)
          if (sp.ge.0.d0) then
c left (i=1) or right (i=2)
c            write(6,*)'shock, left/right',i,j
            wstar(1)=r(i,j)
            wstar(2)=an(1)*v(i,j)-an(2)*vt
            wstar(3)=an(2)*v(i,j)+an(1)*vt
            wstar(4)=p(i,j)
          else
c middle left (i=1) or middle right (i=2)
c            write(6,*)'shock, middle',i,j
            wstar(1)=rm(i,j)
            wstar(2)=an(1)*vm(i,j)-an(2)*vt
            wstar(3)=an(2)*vm(i,j)+an(1)*vt
            wstar(4)=pm(i,j)
          end if
        else
c rarefaction cases
          if (isign*v(i,j)-c(i,j).ge.0.d0) then
c left (i=1) or right (i=2)
c            write(6,*)'rarefaction, left/right',i,j
            wstar(1)=r(i,j)
            wstar(2)=an(1)*v(i,j)-an(2)*vt
            wstar(3)=an(2)*v(i,j)+an(1)*vt
            wstar(4)=p(i,j)
          else
c left middle or sonic (i=1) or right middle or sonic (i=2)
            if (isign*vm(i,j)-sm(i,j).gt.0.d0) then
c sonic left (i=1) or right (i=2)
              if (method.eq.0) then
c                write(6,*)'rarefaction, sonic',i,j
                call sonicmg (i,j,rm,vm,pm,sm,rsonic,vsonic,psonic)
                wstar(1)=rsonic
                wstar(2)=an(1)*vsonic-an(2)*vt
                wstar(3)=an(2)*vsonic+an(1)*vt
                wstar(4)=psonic
              elseif (method.eq.1) then
                arg=(2.d0+isign*gm1(j)*v(i,j)/c(i,j))/gp1(j)
                vn=c(i,j)*arg*isign
                wstar(1)=r(i,j)*arg**(2.d0/gm1(j))
                wstar(2)=an(1)*vn-an(2)*vt
                wstar(3)=an(2)*vn+an(1)*vt
                wstar(4)=(p(i,j)+p0(j))*arg**(2.d0*gam(j)/gm1(j))-p0(j)
              else
                write(6,*)'Error (getfxmg) : method not supported'
                stop
              end if
            else
c middle left (i=1) or right (i=2)
c              write(6,*)'rarefaction, middle',i,j
              wstar(1)=rm(i,j)
              wstar(2)=an(1)*vm(i,j)-an(2)*vt
              wstar(3)=an(2)*vm(i,j)+an(1)*vt
              wstar(4)=pm(i,j)
            end if
          end if
        end if
      end if
c
c square of velocity and normal component
      q2=wstar(2)**2+wstar(3)**2
      vn=an(1)*wstar(2)+an(2)*wstar(3)

c compute internal energy at star state
      Vols=1.d0/wstar(1)
      As=AmgV(kmat(j),Vols)
      Bs=BmgV(j,Vols,ier)
      if (ier.ne.0) then
        write(6,*)'Error (getfxmg) : cannot compute Bs'
        stop
      end if
      energy=As*wstar(4)+Bs
c
c add compaction potential for the solid phase case
      if (j.eq.1) then
        energy=energy+compac(alpha(i,j),0)
      end if
c
c add kinetic energy
      energy=wstar(1)*(energy+.5d0*q2)
c
c compute flux
      fact1=alpha(i,j)*aj
      fact2=wstar(1)*vn
      fx(1)=fact1*fact2
      fx(2)=fact1*(fact2*wstar(2)+an(1)*wstar(4))
      fx(3)=fact1*(fact2*wstar(3)+an(2)*wstar(4))
      fx(4)=fact1*(energy+wstar(4))*vn
c
      return
      end
c
c+++++++++++++++
c
      subroutine sonicmg (i,j,rm,vm,pm,sm,rsonic,vsonic,psonic)
c
c compute solid or gas flux for i=side and j=phase
c
      implicit double precision (a-h,o-z)
      dimension rm(2,2),vm(2,2),pm(2,2),sm(2,2)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / matdatmg / kmat(2)
      data tol, itmax / 1.d-13, 20 /
c
c set isign=1 for i=1 and isign=-1 for i=2
      isign=3-2*i
c
c material type for phase j
      kj=kmat(j)
c
c initialize the iteration
      g0=isign*vm(i,j)-sm(i,j)
      g1=isign*v(i,j)-c(i,j)
      r0=rm(i,j)
      r1=r(i,j)
      r2=r1-g1*(r1-r0)/(g1-g0)
c
c density scale
      rscale=.5d0*(r(1,j)+r(2,j))
c
c secant iteration
      do it=1,itmax
        call raremg (i,j,kj,r2,dv,p2,fp2,gp2,s2,ier)
        if (ier.ne.0) then
          write(6,*)'Error (sonicmg) : error from raremg'
          stop
        end if
        v2=v(i,j)-isign*dv
c
c residual and change in density
        g2=isign*v2-s2
        dr=g2*(r2-r1)/(g2-g1)
        err=dabs(dr)/rscale
c
c check for convergence
        if (err.lt.tol) then
          rsonic=r2
          vsonic=v2
          psonic=p2
          return
        end if
c reset
        r1=r2
        g1=g2
        r2=r2-dr

      end do
c
      write(6,*)'Error (sonicmg) : iteration failed to converge'
c
      stop
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Stiffened gas flux subroutines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cmpmiddlesi (j,pm,vm,dv)
c
c middle states of the Riemann problem for a stiffened ideal gas EOS
c
      implicit real*8 (a-h,o-z)
      dimension pm(2),vm(2),dv(2)
      dimension rm(2),vdif(2,2),ishock(2)
      dimension pmsave(2),vmsave(2)
      character*3 itype(2)
      logical mgfxdbg
      common / mgflxerr / mgfxerr,mgfxdbg
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
      common / sigdatmg / gp1(2),em(2),ep(2),a(2,2),b(2,2)
      common / matdatmg / kmat(2)
      data tol, itmax / 1.d-8, 10 /
      data idebug, iout / 0, 99 /
c
c print left and right states
      if (idebug.ne.0) then
        write(iout,*)'Msg (cmpmiddlesi) : left and right states'
        write(iout,40)r(1,j),r(2,j),v(1,j),v(2,j),p(1,j),p(2,j)
   40   format(' rho(L)=',1pe15.8,'   rho(R)=',1pe15.8,/,
     *         '   v(L)=',1pe15.8,'     v(R)=',1pe15.8,/,
     *         '   p(L)=',1pe15.8,'     p(R)=',1pe15.8)
      end if
c
c material type for phase j
      kj=kmat(j)
c
c get initial guess for middle densities
      call guessmg (j,kj,pstar,rm,iguess)
c      write(6,*)pstar
c      read(5,*)idum
c
      if (idebug.ne.0) then
        if (iguess.eq.0) then
          write(iout,*)
     *        'Msg (cmpmiddlesi) : guess based on linearization'
        else
          if (iguess.gt.0) then
            write(iout,*)
     *        'Msg (cmpmiddlesi) : guess based on two-S approx'
          else
            write(iout,*)
     *        'Msg (cmpmiddlesi) : guess based on two-R approx'
          end if
        end if
      end if
c
c some constants and re-compute sound speeds
c (Note: gam(j), gm1(j) and p0(j) are set in guessmg)
      gp1(j)=gam(j)+1.d0
      em(j)=0.5d0*gm1(j)/gam(j)
      ep(j)=0.5d0*gp1(j)/gam(j)
      do i=1,2
        a(i,j)=2.d0/(gp1(j)*r(i,j))
        b(i,j)=gm1(j)*(p(i,j)+p0(j))/gp1(j)+p0(j)
        c2=gam(j)*(p(i,j)+p0(j))/r(i,j)
        if (c2.le.0.d0) then
          write(6,*)'Error (cmpmiddlesi) : c2.le.0, i,j =',i,j
          stop
        end if
        c(i,j)=dsqrt(c2)
      end do
c
c Newton iteration to find middle states
      do it=1,itmax
c
c evaluate rarefaction or shock from each side
        do i=1,2
c
          if (pstar.le.p(i,j)) then
c
c rarefaction solution
            ishock(i)=0
            arg=(pstar+p0(j))/(p(i,j)+p0(j))
            fact=2.d0*c(i,j)/gm1(j)
            vdif(1,i)=fact*(arg**em(j)-1.d0)
            vdif(2,i)=1.d0/(r(i,j)*c(i,j)*arg**ep(j))
c
          else
c
c shock solution
            ishock(i)=1
            arg=pstar+b(i,j)
            fact=dsqrt(a(i,j)/arg)
            diff=pstar-p(i,j)
            vdif(1,i)=fact*diff
            vdif(2,i)=fact*(1.d0-0.5d0*diff/arg)
c
          end if
c
        end do
c
c save guess
        if (it.eq.1) then
          pmsave(1)=pstar
          vmsave(1)=.5d0*(v(1,j)-vdif(1,1)+v(2,j)+vdif(1,2))
          pmsave(2)=pmsave(1)
          vmsave(2)=vmsave(1)
        end if
c
c determine change to pressure in the middle state
        dp=(vdif(1,1)+vdif(1,2)+v(2,j)-v(1,j))/(vdif(2,1)+vdif(2,2))
        pstar=pstar-dp
c
        if (pstar.le.0.d0) then
          if (mgfxdbg) then
           write(6,*)'Warning (cmpmiddlesi) : negative middle pressure'
          end if
          do i=1,2
            pm(i)=pmsave(i)
            vm(i)=vmsave(i)
          end do
          mgfxerr=1
          return
        end if
c
        err=dabs(dp)/(pstar+p0(j))
c
c label transition type, rarefaction (R) or shock (S)
        do i=1,2
          if (ishock(i).eq.0) then
            write(itype(i),50)i
   50       format(i1,'-R')
          else
            write(itype(i),51)i
   51       format(i1,'-S')
          end if
        end do
c
c print convergence information
        if (idebug.ne.0) then
          write(iout,100)it,itype(1),itype(2),pstar,dp,err
  100     format(' * it=',i2,' : ',a3,',',a3,',  pstar=',1pe15.8,
     *           ', dp=',1pe10.3,', err=',1pe10.3)
        end if
c
        if (err.lt.tol) then
          pm(1)=pstar
          vm(1)=.5d0*(v(1,j)-vdif(1,1)+v(2,j)+vdif(1,2))
          pm(2)=pm(1)
          vm(2)=vm(1)
          if (idebug.ne.0) then
            write(iout,200)pm(1),vm(1)
  200       format(' * converged : pstar=',1pe15.8,', vstar=',1pe15.8)
          end if
          dv(1)=vdif(2,1)
          dv(2)=vdif(2,2)
          return
        end if
c
      end do
c
      if (mgfxdbg) then
       write(6,*)'Warning (cmpmiddlesi) : iteration failed to converge'
      end if
      do i=1,2
        pm(i)=pmsave(i)
        vm(i)=vmsave(i)
      end do
      mgfxerr=1
c
      return
      end
c
c+++++++++++++++
c
      subroutine cmpcouplesi (pm,vm,dv,rm0,rm,sm,linsc)
c
c compute coupled middle states assuming a stiffenend ideal EOS
c for each phase
c
      implicit real*8 (a-h,o-z)
      dimension pm(2,2),vm(2,2),dv(2,2),rm(2,2),sm(2,2)
      dimension pm0(2,2),vm0(2,2),dummy(2),gdum(4)
      logical ctest
      logical mgfxdbg
      common / mgflxerr / mgfxerr,mgfxdbg
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
      common / sigdatmg / gp1(2),em(2),ep(2),a(2,2),b(2,2)
      common / scales / pscale,vscale,gfact(4)
c      common / flxdat / afxcnt(3),nfxcnt(3)
      data sig, frac, tol / 10.d0, .99d0, 1.d-8 /
c
c      nfxcnt(1)=nfxcnt(1)+1
c
c check whether the phases are coupled
      atol=1.d-14
      if (dabs(alpha(2,1)-alpha(1,1)).gt.atol) then
c
c reset flux error flag
        mgfxerr=0
c
c compute gas density
        dvm=vm(1,2)-vm(1,1)
        if (dvm.ge.0.d0) then
          if (pm(1,2).gt.p(1,2)) then
            r0=r(1,2)*(gp1(2)*(pm(1,2)+p0(2))/(p(1,2)+p0(2))+gm1(2))
     *               /(gm1(2)*(pm(1,2)+p0(2))/(p(1,2)+p0(2))+gp1(2))
          else
            r0=r(1,2)*((pm(1,2)+p0(2))/(p(1,2)+p0(2)))**(1.d0/gam(2))
          end if
        else
          if (pm(2,2).gt.p(2,2)) then
            r0=r(2,2)*(gp1(2)*(pm(2,2)+p0(2))/(p(2,2)+p0(2))+gm1(2))
     *               /(gm1(2)*(pm(2,2)+p0(2))/(p(2,2)+p0(2))+gp1(2))
          else
            r0=r(2,2)*((pm(2,2)+p0(2))/(p(2,2)+p0(2)))**(1.d0/gam(2))
          end if
        end if
c
c compute pressure and velocity scales
        pscale=max(pm(1,1),pm(1,2))
        vscale=dsqrt(gam(2)*(pm(1,2)+p0(2))/r0)
        gfact(1)=1.d0/vscale
        gfact(2)=gfact(1)/pscale**(1.d0/gam(2))
        gfact(3)=1.d0/pscale
        gfact(4)=gfact(1)**2
c
c "middle" value for alpha to linearize about
c       asm=.5d0*(alpha(1,1)+alpha(2,1))
        if (alpha(1,1).gt.0.5d0) then
          if (alpha(2,1).lt.0.5d0) then
            asm=.5d0
          else
            asm=min(alpha(1,1),alpha(2,1))
          end if
        else
          if (alpha(2,1).gt.0.5d0) then
            asm=.5d0
          else
            asm=max(alpha(1,1),alpha(2,1))
          end if
        end if
        agm=1.d0-asm
c
c pressure scales (just the pstar's for each phase plus stiffening)
        pscale1=pm(1,1)+p0(1)
        pscale2=pm(1,2)+p0(2)
c
c linear update for middle solid pressures
        eps=alpha(2,1)-alpha(1,1)
        del=(pm(1,2)-pm(1,1))*eps/(asm*(dv(1,1)+dv(2,1)))
        dpm1=-dv(2,1)*del
        dpm2= dv(1,1)*del
        pm(1,1)=pm(1,1)+dpm1
        pm(2,1)=pm(2,1)+dpm2
        err1=max(dabs(dpm1),dabs(dpm2))/pscale1
c
c limit velocity difference
        gampg=gam(2)*(pm(1,2)+p0(2))
        cm=frac*dsqrt(gampg/r0)
c        if (dabs(dvm/cm).ge.1.d0) then
c          write(6,*)dvm,cm
c          stop
c        end if
        arg1=sig*(1.d0+dvm/cm)
        arg2=sig*(1.d0-dvm/cm)
        dvm=.5d0*cm*dlog(dcosh(arg1)/dcosh(arg2))/sig
c
c linear update for middle gas pressures
        del=gampg*dvm*eps/(agm*(gampg-r0*dvm**2)*(dv(1,2)+dv(2,2)))
        dpm1= (r0*dvm*dv(2,2)+1.d0)*del
        dpm2=-(r0*dvm*dv(1,2)-1.d0)*del
        pm(1,2)=pm(1,2)+dpm1
        pm(2,2)=pm(2,2)+dpm2
        err2=max(dabs(dpm1),dabs(dpm2))/pscale2
c
c compute residual of the jump conditions
        call getgsi (pm,vm,dummy,gdum,resid)
        resid=max(err1,err2)
c
c if residual (error) is too big, then do Newton
        if (resid.gt.tol.and.linsc.eq.0) then
c
c save current middle pressures and velocities
          do i=1,2
            do j=1,2
              pm0(i,j)=pm(i,j)
              vm0(i,j)=vm(i,j)
            end do
          end do
c
          call newtonsi (pm,vm,isonic,ctest)
c
c if Newton fails, then reset middle pressures to their linearized values
          if (ctest.and.isonic.eq.0) then
c            nfxcnt(3)=nfxcnt(3)+1
          else
            do i=1,2
              do j=1,2
                pm(i,j)=pm0(i,j)
                vm(i,j)=vm0(i,j)
              end do
            end do
c            nfxcnt(2)=nfxcnt(2)+1
c could set flux error flag to be nonzero, but this is not done at the moment
c            mgfxerr=1
          end if
c
        else
c          nfxcnt(2)=nfxcnt(2)+1
        end if
c
      end if
c
c compute middle densities across shocks/rarefactions
      do j=1,2
        do i=1,2
          if (pm(i,j).gt.p(i,j)) then
            rm(i,j)=r(i,j)*(gp1(j)*(pm(i,j)+p0(j))
     *                             /(p(i,j)+p0(j))+gm1(j))
     *                    /(gm1(j)*(pm(i,j)+p0(j))
     *                             /(p(i,j)+p0(j))+gp1(j))
          else
            rm(i,j)=r(i,j)*((pm(i,j)+p0(j))
     *                      /(p(i,j)+p0(j)))**(1.d0/gam(j))
          end if
        end do
      end do
c
c compute the density between the solid and gas contacts
      dvm=vm(1,2)-vm(1,1)
      if (dvm.ge.0.d0) then
        rm0=rm(1,2)*((pm(2,2)+p0(2))/(pm(1,2)+p0(2)))**(1.d0/gam(2))
      else
        rm0=rm(2,2)*((pm(1,2)+p0(2))/(pm(2,2)+p0(2)))**(1.d0/gam(2))
      end if
c
c compute wave velocities
      do j=1,2
        do i=1,2
          if (pm(i,j).gt.p(i,j)) then
            sm(i,j)=c(i,j)*dsqrt(ep(j)*(pm(i,j)+p0(j))
     *                                 /(p(i,j)+p0(j))+em(j))
          else
            sm(i,j)=dsqrt(gam(j)*(pm(i,j)+p0(j))/rm(i,j))
          end if
        end do
      end do
c
      return
      end
c
c+++++++++++++++
c
      subroutine newtonsi (pm,vm,isonic,ctest)
c
      implicit double precision (a-h,o-z)
      dimension pm(2,2),vm(2,2),rm(2),g(4),dg(4,4),dga(4,4)
      logical ctest
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
      common / sigdatmg / gp1(2),em(2),ep(2),a(2,2),b(2,2)
      common / scales / pscale,vscale,gfact(4)
      data pfact, tol, itmax / 1.d-4, 1.d-8, 10 /
      data idebug / 0 /
c
c solid pressure is used to compute source
      source0=pm(2,1)*alpha(2,1)-pm(1,1)*alpha(1,1)
c
c Newton iteration to find pm
      it=0
    1 it=it+1
c
c get right-hand-side vector
        call getgsi (pm,vm,rm,g,resid)
c
c check sonic condition on the middle states
        isonic=0
        call chkvelmg (it,rm,vm,pm(1,2),isonic)
        if (isonic.ne.0) return
c
c get Jacobian matrix (analytic)
        call getdgasi (dg,pm,vm,rm)
c
c get Jacobian matrix (finite difference)
c        call getdgfsi (pm,dg,g)
c        if (idebug.gt.1.and.it.eq.1) then
c          call getdgasi (dga,pm,vm,rm)
c          write(88,*)'** diffs in Jacobian'
c          do ii=1,4
c            write(88,200)ii,(dga(ii,jj)-dg(ii,jj),jj=1,4)
c  200       format(5x,i2,4(1x,1pe10.3))
c          end do
c        end if
c
c compute correction, overwrite
        call solvesi (dg,g,ier)
        if (ier.ne.0) then
          ctest=.false.
          return
        end if
c
c update pm
        pm(1,1)=dmax1(pm(1,1)-g(1),pfact*pm(1,1))
        pm(2,1)=dmax1(pm(2,1)-g(2),pfact*pm(2,1))
        pm(1,2)=dmax1(pm(1,2)-g(3),pfact*pm(1,2))
        pm(2,2)=dmax1(pm(2,2)-g(4),pfact*pm(2,2))
c
c relative error in the source calculation
        source=pm(2,1)*alpha(2,1)-pm(1,1)*alpha(1,1)
        err1=dabs(source-source0)/min(pm(1,1)+p0(1),pm(2,1)+p0(1))
c
c relative error in the solid pressures
        err2=max(dabs(g(1))/(pm(1,1)+p0(1)),dabs(g(2))/(pm(2,1)+p0(1)))
c
c relative error in the gas pressures
        err3=max(dabs(g(3))/(pm(1,2)+p0(2)),dabs(g(4))/(pm(2,2)+p0(2)))
c
        if (idebug.ne.0) then
          if (it.eq.1) write(88,*)'** Newton iteration'
          write(88,100)it,err1,err2,err3
  100     format('    it=',i2,' err=',3(1x,1pe10.3))
        end if
c
c check for convergence
      if (max(err1,err2,err3).gt.tol) then
        if (it.lt.itmax) then
          source0=source
          goto 1
        else
          ctest=.false.
          return
        end if
      end if
c
c if converged, then check solution
      ctest=.true.
      call soltestsi (pm,vm,ctest)
c
      return
      end
c
c+++++++++++++++
c
      subroutine solvesi (dg,g,ier)
c
      implicit double precision (a-h,o-z)
      dimension dg(4,4),g(4),a(4,5)
      data tol / 1.d-12 /
c
      ier=0
c
c set up augmented matrix a=[dg,g] and compute Frobenius norm of dg
      anorm=0.d0
      do i=1,4
        do j=1,4
          a(i,j)=dg(i,j)
          anorm=anorm+dg(i,j)**2
        end do
        a(i,5)=g(i)
      end do
      anorm=dsqrt(anorm)
c
c tolerance for singular system
      atol=anorm*tol
c
c Gaussian elimination with partial pivoting
      do k=1,3
        kpiv=k
        apiv=dabs(a(k,k))
        do i=k+1,4
          if (dabs(a(i,k)).gt.apiv) then
            kpiv=i
            apiv=dabs(a(i,k))
          end if
        end do
        if (kpiv.ne.k) then
          do j=k,5
            atmp=a(k,j)
            a(k,j)=a(kpiv,j)
            a(kpiv,j)=atmp
          end do
        end if
        if (dabs(a(k,k)).lt.atol) then
          ier=1
          return
        end if
        do i=k+1,4
          fact=a(i,k)/a(k,k)
          do j=k+1,5
            a(i,j)=a(i,j)-fact*a(k,j)
          end do
        end do
      end do
c
c backward substitution
      g(4)=a(4,5)/a(4,4)
      do i=3,1,-1
        sum=a(i,5)
        do j=i+1,4
          sum=sum-a(i,j)*g(j)
        end do
        g(i)=sum/a(i,i)
      end do
c
      return
      end
c
c+++++++++++++++
c
      subroutine soltestsi (pm,vm,ctest)
c
      implicit double precision (a-h,o-z)
      dimension pm(2,2),vm(2,2),sp(3,2)
      logical ctest
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
      common / sigdatmg / gp1(2),em(2),ep(2),a(2,2),b(2,2)
      data tiny / 1.d-8 /
c
      if (vm(1,1).lt.vm(1,2)) then
        sp(2,2)=vm(2,2)
      else
        sp(2,2)=vm(1,2)
      end if
      sp(2,1)=vm(1,1)
c
      j=1
      if (pm(1,j).gt.p(1,j)) then
        sp(1,j)=v(1,j)-c(1,j)*dsqrt(ep(j)*(pm(1,j)+p0(j))
     *                                    /(p(1,j)+p0(j))+em(j))
        if (sp(1,j)-tiny.gt.v(1,j)-c(1,j)) then
          ctest=.false.
          return
        end if
      else
        rm=r(1,j)*((pm(1,j)+p0(j))/(p(1,j)+p0(j)))**(1.d0/gam(j))
        sp(1,j)=vm(1,j)-dsqrt(gam(j)*(pm(1,j)+p0(j))/rm)
        if (sp(1,j)+tiny.lt.v(1,j)-c(1,j)) then
          ctest=.false.
          return
        end if
      end if
      if (pm(2,j).gt.p(2,j)) then
        sp(3,j)=v(2,j)+c(2,j)*dsqrt(ep(j)*(pm(2,j)+p0(j))
     *                                    /(p(2,j)+p0(j))+em(j))
        if (sp(3,j)+tiny.lt.v(2,j)+c(2,j)) then
          ctest=.false.
          return
        end if
      else
        rm=r(2,j)*((pm(2,j)+p0(j))/(p(2,j)+p0(j)))**(1.d0/gam(j))
        sp(3,j)=vm(2,j)+dsqrt(gam(j)*(pm(2,j)+p0(j))/rm)
        if (sp(3,j)-tiny.gt.v(2,j)+c(2,j)) then
          ctest=.false.
          return
        end if
      end if
c
      j=2
      if (pm(1,j).gt.p(1,j)) then
        sp(1,j)=v(1,j)-c(1,j)*dsqrt(ep(j)*(pm(1,j)+p0(j))
     *                                    /(p(1,j)+p0(j))+em(j))
        if (sp(1,j)-tiny.gt.v(1,j)-c(1,j)) then
          ctest=.false.
          return
        end if
      else
        rm=r(1,j)*((pm(1,j)+p0(j))/(p(1,j)+p0(j)))**(1.d0/gam(j))
        sp(1,j)=vm(1,j)-dsqrt(gam(j)*(pm(1,j)+p0(j))/rm)
        if (sp(1,j)+tiny.lt.v(1,j)-c(1,j)) then
          ctest=.false.
          return
        end if
      end if
      if (pm(2,j).gt.p(2,j)) then
        sp(3,j)=v(2,j)+c(2,j)*dsqrt(ep(j)*(pm(2,j)+p0(j))
     *                                    /(p(2,j)+p0(j))+em(j))
        if (sp(3,j)+tiny.lt.v(2,j)+c(2,j)) then
          ctest=.false.
          return
        end if
      else
        rm=r(2,j)*((pm(2,j)+p0(j))/(p(2,j)+p0(j)))**(1.d0/gam(j))
        sp(3,j)=vm(2,j)+dsqrt(gam(j)*(pm(2,j)+p0(j))/rm)
        if (sp(3,j)-tiny.gt.v(2,j)+c(2,j)) then
          ctest=.false.
          return
        end if
      end if
c
      do j=1,2
        if (sp(2,1)+tiny.lt.sp(1,j).or.
     *      sp(2,1)-tiny.gt.sp(3,j)) then
          ctest=.false.
          return
        end if
        if (sp(2,j)+tiny.lt.sp(1,j).or.
     *      sp(2,j)-tiny.gt.sp(3,j)) then
          ctest=.false.
          return
        end if
      end do
c
      return
      end
c
c+++++++++++++++
c
      subroutine getgsi (pm,vm,rm,g,resid)
c
      implicit double precision (a-h,o-z)
      dimension pm(2,2),vm(2,2),rm(2),g(4)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
      common / sigdatmg / gp1(2),em(2),ep(2),a(2,2),b(2,2)
      common / scales / pscale,vscale,gfact(4)
c
c compute velocities across acoustic fields
      do j=1,2
        do i=1,2
          isign=2*i-3
          if (pm(i,j).gt.p(i,j)) then
            arg=pm(i,j)+b(i,j)
            fact=dsqrt(a(i,j)/arg)
            vm(i,j)=v(i,j)+isign*fact*(pm(i,j)-p(i,j))
          else
            arg=(pm(i,j)+p0(j))/(p(i,j)+p0(j))
            fact=2.d0*c(i,j)/gm1(j)
            vm(i,j)=v(i,j)+isign*fact*(arg**em(j)-1.d0)
          end if
        end do
      enddo
c
c check whether gas contact is to the left or
c to the right of solid contact, compute gas densities
      if (vm(1,1).gt.vm(1,2)) then
c gas contact on the left
        if (pm(2,2).gt.p(2,2)) then
          rm(2)=r(2,2)*(gp1(2)*(pm(2,2)+p0(2))/(p(2,2)+p0(2))+gm1(2))
     *                /(gm1(2)*(pm(2,2)+p0(2))/(p(2,2)+p0(2))+gp1(2))
        else
          rm(2)=r(2,2)*((pm(2,2)+p0(2))/(p(2,2)+p0(2)))**(1.d0/gam(2))
        end if
        rm(1)=rm(2)*((pm(1,2)+p0(2))/(pm(2,2)+p0(2)))**(1.d0/gam(2))
      else
c gas contact on the right
        if (pm(1,2).gt.p(1,2)) then
          rm(1)=r(1,2)*(gp1(2)*(pm(1,2)+p0(2))/(p(1,2)+p0(2))+gm1(2))
     *                /(gm1(2)*(pm(1,2)+p0(2))/(p(1,2)+p0(2))+gp1(2))
        else
          rm(1)=r(1,2)*((pm(1,2)+p0(2))/(p(1,2)+p0(2)))**(1.d0/gam(2))
        end if
        rm(2)=rm(1)*((pm(2,2)+p0(2))/(pm(1,2)+p0(2)))**(1.d0/gam(2))
      end if
c
c compute constraints across solid particle path
      v1=vm(1,2)-vm(1,1)
      v2=vm(2,2)-vm(2,1)
      h1=gam(2)*(pm(1,2)+p0(2))/(gm1(2)*rm(1))
      h2=gam(2)*(pm(2,2)+p0(2))/(gm1(2)*rm(2))
      g(1)=vm(2,1)-vm(1,1)
      g(2)= alpha(2,2)*((pm(2,2)+p0(2))**(1.d0/gam(2)))*v2
     *     -alpha(1,2)*((pm(1,2)+p0(2))**(1.d0/gam(2)))*v1
      g(3)= alpha(2,1)*pm(2,1)+alpha(2,2)*(pm(2,2)+rm(2)*v2**2)
     *     -alpha(1,1)*pm(1,1)-alpha(1,2)*(pm(1,2)+rm(1)*v1**2)
      g(4)=.5d0*v2**2+h2-.5d0*v1**2-h1
c
c scale residual
      do i=1,4
        g(i)=gfact(i)*g(i)
      end do
c
c compute residual
      resid=max(dabs(g(1)),dabs(g(2)),dabs(g(3)),dabs(g(4)))
c
      return
      end
c
c+++++++++++++++
c
      subroutine getdgasi (dg,pm,vm,rm)
c
c exact Jacobian
c
      implicit double precision (a-h,o-z)
      dimension dg(4,4),pm(2,2),vm(2,2),rm(2)
      dimension fp(2,2),gk(2),gp(2),drdp(2,2)
      common / prmdatmg / alpha(2,2),r(2,2),v(2,2),p(2,2),c(2,2),
     *                    gamk(2,2),p0k(2,2)
      common / gasdatmg / gam(2),gm1(2),p0(2)
      common / sigdatmg / gp1(2),em(2),ep(2),a(2,2),b(2,2)
      common / scales / pscale,vscale,gfact(4)
c
      j=1
      do i=1,2
        if (pm(i,j).gt.p(i,j)) then
          fp(i,j)=(gp1(j)*(pm(i,j)+p0(j))+(3.d0*gam(j)-1.d0)*
     *             (p(i,j)+p0(j)))/dsqrt(2.d0*r(i,j)*(gp1(j)*
     *             (pm(i,j)+p0(j))+gm1(j)*(p(i,j)+p0(j)))**3)
        else
          fp(i,j)=(((pm(i,j)+p0(j))/(p(i,j)+p0(j)))**(-ep(j)))/
     *                                (r(i,j)*c(i,j))
        end if
      end do
c
      j=2
      do i=1,2
        if (pm(i,j).gt.p(i,j)) then
          fp(i,j)=(gp1(j)*(pm(i,j)+p0(j))+(3.d0*gam(j)-1.d0)*
     *             (p(i,j)+p0(j)))/dsqrt(2.d0*r(i,j)*(gp1(j)*
     *             (pm(i,j)+p0(j))+gm1(j)*(p(i,j)+p0(j)))**3)
          gk(i)=r(i,2)*(gp1(2)*(pm(i,2)+p0(2))+gm1(2)*(p(i,2)+p0(2)))/
     *                 (gm1(2)*(pm(i,2)+p0(2))+gp1(2)*(p(i,2)+p0(2)))
          gp(i)=4.0d0*gam(2)*r(i,2)*(p(i,2)+p0(2))/
     *          ((gm1(2)*(pm(i,2)+p0(2))+gp1(2)*(p(i,2)+p0(2)))**2)
        else
          fp(i,j)=((pm(i,j)/p(i,j))**(-ep(j)))/(r(i,j)*c(i,j))
          gk(i)=r(i,2)*((pm(i,2)+p0(2))/(p(i,2)+p0(2)))**(1/gam(2))
          gp(i)=gk(i)/(gam(2)*(pm(i,2)+p0(2)))
        end if
      end do
c
      gi=1.d0/gam(2)
      rat=((pm(1,2)+p0(2))/(pm(2,2)+p0(2)))**gi
      if (vm(1,1).gt.vm(1,2)) then
        drdp(1,1)=gk(2)*gi*rat/(pm(1,2)+p0(2))
        drdp(1,2)=(gp(2)-gk(2)*gi/(pm(2,2)+p0(2)))*rat
        drdp(2,1)=0.d0
        drdp(2,2)=gp(2)
      else
        drdp(1,1)=gp(1)
        drdp(1,2)=0.d0
        drdp(2,1)=(gp(1)-gk(1)*gi/(pm(1,2)+p0(2)))/rat
        drdp(2,2)=gk(1)*gi/(pm(2,2)+p0(2))/rat
      end if
c
      dv1=vm(1,2)-vm(1,1)
      dv2=vm(2,2)-vm(2,1)
c
      dg(1,1)=fp(1,1)
      dg(1,2)=fp(2,1)
      dg(1,3)=0.d0
      dg(1,4)=0.d0
      dg(2,1)=-alpha(1,2)*((pm(1,2)+p0(2))**gi)*fp(1,1)
      dg(2,2)=-alpha(2,2)*((pm(2,2)+p0(2))**gi)*fp(2,1)
      dg(2,3)= alpha(1,2)*((pm(1,2)+p0(2))**gi)*(fp(1,2)
     *                                         -dv1*gi/(pm(1,2)+p0(2)))
      dg(2,4)= alpha(2,2)*((pm(2,2)+p0(2))**gi)*(fp(2,2)
     *                                         +dv2*gi/(pm(2,2)+p0(2)))
      dg(3,1)=-alpha(1,1)-2.d0*alpha(1,2)*rm(1)*dv1*fp(1,1)
      dg(3,2)=alpha(2,1)-2.d0*alpha(2,2)*rm(2)*dv2*fp(2,1)
      dg(3,3)=alpha(2,2)*drdp(2,1)*(dv2**2)-alpha(1,2)*
     *                 (1.d0+drdp(1,1)*(dv1**2)-2.d0*rm(1)*dv1*fp(1,2))
      dg(3,4)=-alpha(1,2)*drdp(1,2)*(dv1**2)+alpha(2,2)*
     *                 (1.d0+drdp(2,2)*(dv2**2)+2.d0*rm(2)*dv2*fp(2,2))
      dg(4,1)=-dv1*fp(1,1)
      dg(4,2)=-dv2*fp(2,1)
      dg(4,3)= dv1*fp(1,2)
     *        -gam(2)/gm1(2)*( 1.d0/rm(1)
     *                        +drdp(2,1)*(pm(2,2)+p0(2))/(rm(2)**2)
     *                        -drdp(1,1)*(pm(1,2)+p0(2))/(rm(1)**2))
      dg(4,4)= dv2*fp(2,2)
     *        +gam(2)/gm1(2)*( 1.d0/rm(2)
     *                        +drdp(1,2)*(pm(1,2)+p0(2))/(rm(1)**2)
     *                        -drdp(2,2)*(pm(2,2)+p0(2))/(rm(2)**2))
c
      do i=1,4
        do j=1,4
          dg(i,j)=gfact(i)*dg(i,j)
        end do
      end do
c
      return
      end
c
c+++++++++++++++
c
      subroutine getdgfsi (pm,dg,g0)
c
      implicit double precision (a-h,o-z)
      dimension pm(2,2),vm(2,2),rm(2),dg(4,4),g0(4),g(4)
      data delta / 1.d-6 /
c
c finite difference approximation of the jacobian
      do i=1,2
        do j=1,2
          n=2*(j-1)+i
          dpm=delta*pm(i,j)
          pm(i,j)=pm(i,j)+dpm
          call getgsi (pm,vm,rm,g,resid)
          do k=1,4
            dg(k,n)=(g(k)-g0(k))/dpm
          end do
          pm(i,j)=pm(i,j)-dpm
        end do
      end do
c
      return
      end
