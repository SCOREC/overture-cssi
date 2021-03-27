! included in interfaceOpt3dOrder4.bf90

! -------------------------------------------------------------------------
! Macro: Evaluate DISPERSIVE forcing terms, FOURTH-ORDER AND 3D
!   This macro can be usedto eval values in either domain 1 or domain 2
!
! Input:
!   FACE : LEFT or : RIGHT
!   fev(n) : forcing on E equation: E_{tt} = c^2 Delta(E) + ... + fev
!   fpv(n,jv) : forcing on equation for P_{n,jv} 
! Output
!   fp(n) : 
!   c2PttLEsum   : coeff of L(E) in P.tt     (second-order)
!   c4PttLEsum   : coeff of L(E) in P.tt     (fourth-order)
!   c4PttLLEsum  : coeff of L*L(E) in P.tt 
! ------------------------------------------------------------------------
#beginMacro getUnifiedDispersiveForcing3dOrder4(FACE,k1,k2,k3,fp,fpv,fev,p,pn,pm,q,qn,qm,u,un,um, \
            dispersionModel,nonlinearModel,numberOfPolarizationVectors,numberOfAtomicLevels,alphaP,c,\
            c2PttEsum,c2PttLEsum,c4PttLEsum,c4PttLLEsum,c2PttttLEsum,c2PttttLLEsum,a0v,a1v,\
            pnec,prc,peptc,b0v,b1v,LE,LLE,LEm,LfE,LfP,fEt,fEtt,fPt,fPtt,fnv,fntv,\
            pevtt,pevttx,pevtty,pevttz,pevttL,pevtttt,evx,evy,evz,evnx,evny,evnz,\
            fevx,fevy,fevz,fpvx,fpvy,fpvz,LEx,LEy,LEz,fPttx,fPtty,fPttz,fLPtt,fPtttt)

  do n=0,nd-1
    fp(n)=0.
    fPttx(n) =0.
    fPtty(n) =0.
    fPttz(n) =0.
    fLPtt(n) =0.
    fPtttt(n)=0.
  end do

  !------------------
  ! MLA
  !------------------
  if( dispersionModel.ne.noDispersion .and. nonlinearModel.ne.noNonlinearModel) then

    nce = pxc+nd*numberOfPolarizationVectors

    dtsq = dt**2

    ! -----------------------------------------
    ! order 2 (E, P, N) at the interface (fictitious step)
    !------------------------------------------
 
    do n=0,nd-1
      ec = ex +n
      ev0  =  un(k1,k2,k3,ec)
      ev   =  u(k1,k2,k3,ec) ! time where we need to fill in ghost points

      pSum = 0.
      do jv=0,numberOfPolarizationVectors-1

        ! a0=a0v(jv)
        ! a1=a1v(jv)
        ! b0=b0v(jv)
        ! b1=b1v(jv)
        ! alpha=alphaP

        pc = n + jv*nd  

        pv0 =  pn(k1,k2,k3,pc)
        pv  =  p(k1,k2,k3,pc)
 
        pvn = 2.*pv-pn(k1,k2,k3,pc) + 0.5*dt*b1v(jv)*pn(k1,k2,k3,pc) - dtsq*b0v(jv)*pv + dtsq*fpv(n,jv)

        do na = 0,numberOfAtomicLevels-1 ! \Delta N^n*E^n
          pvn = pvn + dtsq*pnec(jv,na)*q(k1,k2,k3,na)*ev
        enddo ! na

        pvec(n,jv)= pvn/( 1.+.5*dt*b1v(jv) ) ! time + dt

        pSum = pSum + pvec(n,jv) -2.*pv + pv0 ! keep sum

      enddo ! jv

      evec(n) = (2.*ev-ev0) + dtsq*LE(n) - alphaP*pSum + dtsq*fev(n) ! cSq is already in LE

    enddo ! n

    ! N outside of space loop
    ! 1st derivative
    do na=0,numberOfAtomicLevels-1
      qt(na) = fnv(na)
      do jv = 0,numberOfAtomicLevels-1
        qt(na) = qt(na)+prc(na,jv)*q(k1,k2,k3,jv)
      enddo
      do n=0,nd-1
        do jv=0,numberOfPolarizationVectors-1
          qt(na) = qt(na) + peptc(na,jv)*u(k1,k2,k3,ex+n)*(pvec(n,jv)-pn(k1,k2,k3,n+jv*nd))/(2.*dt)
        enddo
      enddo
    enddo

    ! 2nd derivative
    do na=0,numberOfAtomicLevels-1
      qtt(na) = fntv(na)
      do jv = 0,numberOfAtomicLevels-1
        qtt(na) = qtt(na)+prc(na,jv)*qt(jv)
      enddo
      do n=0,nd-1
        do jv=0,numberOfPolarizationVectors-1
          qtt(na) = qtt(na) + peptc(na,jv)*(evec(n)-un(k1,k2,k3,ex+n))/(2.*dt)*(pvec(n,jv)-pn(k1,k2,k3,n+jv*nd))/(2.*dt)\
                            + peptc(na,jv)*u(k1,k2,k3,ex+n)*(pvec(n,jv)-2.*p(k1,k2,k3,n+jv*nd)+pn(k1,k2,k3,n+jv*nd))/(dtsq)
        enddo
      enddo
    enddo

    ! taylor expansion
    do na=0,numberOfAtomicLevels-1
      qv(na) = q(k1,k2,k3,na) + dt*qt(na) + dtsq/2.*qtt(na)
    enddo
    
    !----------------------------------------
    ! order 4 update of P at interface (fictitious step)
    !----------------------------------------

    ! second order accurate terms
    do n=0,nd-1
      do jv = 0,numberOfPolarizationVectors-1
        ! #If #p eq "p1"
        ! call ogderiv(ep, 0,0,0,0, xy1(k1,k2,k3,0),xy1(k1,k2,k3,1),0.,t,pxc+jv*nd+n, pe(n)   )
        ! call ogderiv(ep, 1,0,0,0, xy1(k1,k2,k3,0),xy1(k1,k2,k3,1),0.,t,pxc+jv*nd+n, pet(n)   )
        ! call ogderiv(ep, 2,0,0,0, xy1(k1,k2,k3,0),xy1(k1,k2,k3,1),0.,t,pxc+jv*nd+n, pett(n)   )
        ! #Else
        ! call ogderiv(ep, 0,0,0,0, xy2(k1,k2,k3,0),xy2(k1,k2,k3,1),0.,t,pxc+jv*nd+n, pe(n)   )
        ! call ogderiv(ep, 1,0,0,0, xy2(k1,k2,k3,0),xy2(k1,k2,k3,1),0.,t,pxc+jv*nd+n, pet(n)   )
        ! call ogderiv(ep, 2,0,0,0, xy2(k1,k2,k3,0),xy2(k1,k2,k3,1),0.,t,pxc+jv*nd+n, pett(n)   )
        ! #End
        ! pvec(n,jv) = pe(n)
        ! ptv(n,jv) = pet(n)
        ! pttv(n,jv) = pett(n)

        ptv(n,jv) = (pvec(n,jv)-pn(k1,k2,k3,n+jv*nd))/(2.*dt)

        pttv(n,jv) = (pvec(n,jv)-2.*p(k1,k2,k3,n+jv*nd)+pn(k1,k2,k3,n+jv*nd))/dtsq


        ptttv(n,jv) = -b1v(jv)*pttv(n,jv)-b0v(jv)*ptv(n,jv)+fPt(n,jv)
        do na = 0,numberOfAtomicLevels-1 ! update using ODE
          ptttv(n,jv) = ptttv(n,jv) + pnec(jv,na)*qt(na)*u(k1,k2,k3,ex+n) \
                                    + pnec(jv,na)*q(k1,k2,k3,na)*(evec(n)-un(k1,k2,k3,ex+n))/(2.*dt)
        enddo

        pttttv(n,jv) = -b1v(jv)*ptttv(n,jv)-b0v(jv)*pttv(n,jv)+fPtt(n,jv) ! update using ODE
        do na = 0,numberOfAtomicLevels-1
          pttttv(n,jv) = pttttv(n,jv) + pnec(jv,na)*qtt(na)*u(k1,k2,k3,ex+n) \
                                   + 2.*pnec(jv,na)*qt(na)*(evec(n)-un(k1,k2,k3,ex+n))/(2.*dt) \
                                      + pnec(jv,na)*q(k1,k2,k3,na)*(evec(n)-2.*u(k1,k2,k3,ex+n)+un(k1,k2,k3,ex+n))/dtsq
          ! print *, '++++++Dispersive forcing+++++++++++++++'
          ! print *, 'check P derivatives: ', n,jv,ptv(n,jv),pttv(n,jv),ptttv(n,jv),pttttv(n,jv)
          ! print *, na,pnec(jv,na),q(k1,k2,k3,na),qt(na),qtt(na)
        enddo
        ! print *, '++++++Dispersive forcing+++++++++++++++'
        ! print *, 'check P time derivatives: ', n,jv,ptv(n,jv),pttv(n,jv),ptttv(n,jv),pttttv(n,jv)
        ! print *, fPt(n,jv),fPtt(n,jv)
      enddo
      ! print *, evec(n),u(k1,k2,k3,ex+n),un(k1,k2,k3,ex+n)
    enddo

    do n=0,nd-1

      ec = ex +n 
      ev   =  u(k1,k2,k3,ec) ! time where we need to fill in ghost points

      ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
      #perl $ORDER=4;
      ! #perl $ORDER=2; ! should use order 2 since E is only filled at first ghost lines at t
      #If #FACE eq "LEFT"
        opEvalJacobianDerivatives(rsxy1,k1,k2,k3,aj1,1)
        ! uu1 in the next statement defines names of intermediate values
        evalSecondDerivs3d(rsxy1,aj1,u1,k1,k2,k3,ec,uu1,e1)
        ! write(*,'("FACE: u1x,u1y,u1xx,u1yy,u1Lap=",5(1pe12.4))') u1x,u1y,u1xx,u1yy,u1Lap
   
        evx0  = e1x
        evy0  = e1y
        evz0  = e1z
        evLap = e1xx+e1yy+e1zz ! these values use the second order predicted values in the first ghost lines
   
      #Elif #FACE eq "RIGHT"
        opEvalJacobianDerivatives(rsxy2,k1,k2,k3,aj2,1)
        ! uu2 in the next statement defines names of intermediate values
        evalSecondDerivs3d(rsxy2,aj2,u2,k1,k2,k3,ec,uu2,e2)
        ! write(*,'("FACE: u2x,u2y,u2xx,u2yy,u2Lap=",5(1pe12.4))') u2x,u2y,u2xx,u2yy,u2Lap
   
        evx0  = e2x
        evy0  = e2y
        evz0  = e2z
        evLap = e2xx+e2yy+e2zz
   
      #Else
         write(*,'(" interface3d:ERROR: unknown FACE for E")')
         stop 7777
      #End

      do jv=0,numberOfPolarizationVectors-1

        pc = n + jv*nd 
   
        ! Left: u1x,u1y, u1xx, u1yy, u1Lap (ex)
        !       v1x,v1y, v1xx, v1yy, v1Lap (ey) 

        ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
        #perl $ORDER=4;
        ! #perl $ORDER=2;
        #If #FACE eq "LEFT"
          opEvalJacobianDerivatives(rsxy1,k1,k2,k3,aj1,1)
          ! uu1 in the next statement defines names of intermediate values
          evalSecondDerivs3d(rsxy1,aj1,p1,k1,k2,k3,pc,uu1,p1)
          ! write(*,'("FACE: p1x,p1y,p1xx,p1yy,p1Lap=",5(1pe12.4))') p1x,p1y,p1xx,p1yy,p1Lap
   
          LP  = p1Lap ! removed c^2
          evalSecondDerivs3d(rsxy1,aj1,p1n,k1,k2,k3,pc,uu1,p1n)
          ! write(*,'("FACE: p1nxx,p1nyy,p1nLap=",3e12.4)') p1nxx,p1nyy,p1nLap
          LPm = p1nLap ! removed c^2
   
          pvx  = p1x
          pvy  = p1y
          pvz  = p1z

          pvnx = p1nx
          pvny = p1ny
          pvnz = p1nz
   
        #Elif #FACE eq "RIGHT"
          opEvalJacobianDerivatives(rsxy2,k1,k2,k3,aj2,1)
          ! uu1 in the next statement defines names of intermediate values
          evalSecondDerivs3d(rsxy2,aj2,p2,k1,k2,k3,pc,uu2,p2)
          LP  = p2Lap ! removed c^2
          evalSecondDerivs3d(rsxy2,aj2,p2n,k1,k2,k3,pc,uu2,p2n)
          LPm = p2nLap ! removed c^2
   
          pvx  = p2x
          pvy  = p2y
          pvz  = p2z

          pvnx = p2nx
          pvny = p2ny
          pvnz = p2nz
   
        #Else
          write(*,'(" interface3d:ERROR: unknown FACE")')
          stop 7777
        #End
        
        pv0 =  pn(k1,k2,k3,pc)
        pv  =  p(k1,k2,k3,pc)
 
        pvn = 2.*pv-pn(k1,k2,k3,pc) + 0.5*dt*b1v(jv)*pn(k1,k2,k3,pc) + dt**4/12.*pttttv(n,jv) + dt**4/6.*b1v(jv)*ptttv(n,jv) - dtsq*b0v(jv)*pv + dtsq*fpv(n,jv)
        ! pvn = 2.*pv-pn(k1,k2,k3,pc) + 0.5*dt*b1v(jv)*pn(k1,k2,k3,pc) - dtsq*b0v(jv)*pv + dtsq*fpv(n,jv)

        do na = 0,numberOfAtomicLevels-1 ! \Delta N^n*E^n
          pvn = pvn + dtsq*pnec(jv,na)*q(k1,k2,k3,na)*ev
        enddo ! na

        pvec(n,jv)= pvn/( 1.+.5*dt*b1v(jv) ) ! time + dt

        ! 4th order accurate term
        fp(n) = fp(n) + (pvec(n,jv)-2.*p(k1,k2,k3,pc)+pn(k1,k2,k3,pc))/dtsq  - dt**2/12.*pttttv(n,jv)
        
        ! 2nd order accurate terms
        fPtttt(n) = fPtttt(n) + pttttv(n,jv)

        #perl $ORDER=4; ! N is extrapolated in the ghost points
        ! #perl $ORDER=2;
        do na=0,numberOfAtomicLevels-1
          #If #FACE eq "LEFT"
            opEvalJacobianDerivatives(rsxy1,k1,k2,k3,aj1,1)
            ! uu1 in the next statement defines names of intermediate values
            evalSecondDerivs3d(rsxy1,aj1,q1,k1,k2,k3,na,qq1,q1)
            ! write(*,'("FACE: p1x,p1y,p1xx,p1yy,p1Lap=",5(1pe12.4))') p1x,p1y,p1xx,p1yy,p1Lap
       
            qvx  = q1x
            qvy  = q1y
            qvz  = q1z
            qvLap  = q1Lap
       
          #Elif #FACE eq "RIGHT"
            opEvalJacobianDerivatives(rsxy2,k1,k2,k3,aj2,1)
            ! uu2 in the next statement defines names of intermediate values
            evalSecondDerivs3d(rsxy2,aj2,q2,k1,k2,k3,na,qq2,q2)
       
            qvx  = q2x
            qvy  = q2y
            qvz  = q2z
            qvLap  = q2Lap
       
          #Else
             write(*,'(" interface3d:ERROR: unknown FACE for N")')
             stop 7777
          #End

          qex(na) = evx0*q(k1,k2,k3,na)+qvx*ev
          qey(na) = evy0*q(k1,k2,k3,na)+qvy*ev
          qez(na) = evz0*q(k1,k2,k3,na)+qvz*ev
          qeLap(na) = ev*qvLap+q(k1,k2,k3,na)*evLap+2.*evx0*qvx+2.*evy0*qvy+2.*evz0*qvz
        enddo

        ! laplacian
        LPn = 2.*LP-LPm + 0.5*dt*b1v(jv)*LPm - dtsq*b0v(jv)*LP + dtsq*LfP(n,jv)

        do na = 0,numberOfAtomicLevels-1 ! \Delta N^n*E^n
          LPn = LPn + dtsq*pnec(jv,na)*qeLap(na)
        enddo
        ! time derivatives
        fLPtt(n) = fLPtt(n) + c**2*(LPn/(1.+.5*dt*b1v(jv)) - 2.*LP + LPm)/dtsq ! added c^2 to be consistent with GDM codes

        ! x
        pttxa = 2.*pvx-pvnx + 0.5*dt*b1v(jv)*pvnx - dtsq*b0v(jv)*pvx + dtsq*fpvx(n,jv)

        do na = 0,numberOfAtomicLevels-1 ! \Delta N^n*E^n
          pttxa = pttxa + dtsq*pnec(jv,na)*qex(na)
        enddo
        ! time derivatives
        fPttx(n) = fPttx(n) + (pttxa/(1.+.5*dt*b1v(jv)) - 2.*pvx + pvnx)/dtsq

        ! y
        pttya = 2.*pvy-pvny + 0.5*dt*b1v(jv)*pvny - dtsq*b0v(jv)*pvy + dtsq*fpvy(n,jv)

        do na = 0,numberOfAtomicLevels-1 ! \Delta N^n*E^n
          pttya = pttya + dtsq*pnec(jv,na)*qey(na)
        enddo
        ! time derivatives
        fPtty(n) = fPtty(n) + (pttya/(1.+.5*dt*b1v(jv)) - 2.*pvy + pvny)/dtsq

        ! z
        pttza = 2.*pvz-pvnz + 0.5*dt*b1v(jv)*pvnz - dtsq*b0v(jv)*pvz + dtsq*fpvz(n,jv)

        do na = 0,numberOfAtomicLevels-1 ! \Delta N^n*E^n
          pttza = pttza + dtsq*pnec(jv,na)*qez(na)
        enddo
        ! time derivatives
        fPttz(n) = fPttz(n) + (pttza/(1.+.5*dt*b1v(jv)) - 2.*pvz + pvnz)/dtsq
   
        if( .false. .and. twilightZone.eq.1 )then
          write(*,'("")')
          write(*,'("DI4:FACE: k1,k2=",2i3," jv=",i2," n=",i2," ptta,ptte=",2e12.4)') k1,k2,jv,n,(pvec(n,jv)-2.*p(k1,k2,k3,pc)+pn(k1,k2,k3,pc))/dtsq  - dt**2/12.*pttttv(n,jv),pevtt(n,jv)
          write(*,'("        : pttxa,pttxe=",2(1pe12.4)," pttya,pttye=",2(1pe12.4)," pttza,pttze=",2(1pe12.4))') (pttxa/(1.+.5*dt*b1v(jv)) - 2.*pvx + pvnx)/dtsq,pevttx(n,jv),(pttya/(1.+.5*dt*b1v(jv)) - 2.*pvy + pvny)/dtsq,pevtty(n,jv),(pttza/(1.+.5*dt*b1v(jv)) - 2.*pvz + pvnz)/dtsq,pevttz(n,jv)
          write(*,'("        : ptttta,ptttte=",2(1pe12.4))') pttttv(n,jv),pevtttt(n,jv)
          write(*,'("        : pttLa,pttLe=",2(1pe12.4))') (LPn/(1.+.5*dt*b1v(jv)) - 2.*LP + LPm)/dtsq,pevttL(n,jv)
   
          ! write(*,'(" c2PttLE,c2PttE,c2PttEm,c2PttP,c2PttPm,c2PttfE,c2PttfP=",7(1pe12.4))') c2PtLE,c2PtE,c2PtEm,c2PtP,c2PtPm,c2PtfE,c2PtfP
          ! write(*,'(" LEx,evx,evnx,pvx,pvnx,fevx,fpvx=",7(1pe12.4))') LEx(n),evx(n),evnx(n),pvx,pvnx,fevx(n),fpvx(n,jv)
          ! write(*,'(" LEy,evy,evny,pvy,pvny,fevy,fpvy=",7(1pe12.4))') LEy(n),evy(n),evny(n),pvy,pvny,fevy(n),fpvy(n,jv)
   
          ! write(*,'(" LE,LLE,LEm,LP,LPm=",5e12.4)') LE(n),LLE(n),LEm(n),LP,LPm
          ! write(*,'(" LfE,LfP,fEt,fEtt,fPt,fPtt=",6e12.4)') LfE(n),LfP(n,jv),fEt(n),fEtt(n),fPt(n,jv),fPtt(n,jv)
          ! write(*,'(" ev,evn,pv,pvn,fev,fpv=",6e12.4)')ev,evn,pv,pvn,fev(n),fpv(n,jv)
        end if
   
      end do 
    end do ! end do n 
  !-----------------------
  ! GDM
  !----------------------
  elseif( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.noNonlinearModel) then

    c2PttEsum=0.
    c2PttLEsum=0.
    c4PttLEsum=0.
    c4PttLLEsum=0.
   
    c2PttttLEsum=0.
    c2PttttLLEsum=0.
   
    do jv=0,numberOfPolarizationVectors-1
      a0=a0v(jv)
      a1=a1v(jv)
      b0=b0v(jv)
      b1=b1v(jv)
      alpha=alphaP
   
      ! Second-order coefficients: 
      ! Ptt = c2PttLE*LE + c2PttE*E + c2PttEm*Em + c2PttP*P + c2PttPm*Pm + c2PttfE*fE + c2PttfP*fP
      #Include interfaceAdeGdmOrder2.h 
   
      ! Fourth-order coefficients
      ! Ptt = c4PttLE*LE + c4PttE*E + c4PttEm*Em + c4PttP*P + c4PttPm*Pm + c4PttfE*fE + c4PttfP*fP
      !    + c4PttLLE*LLE + c4PttLP*LP + c4PttLEm*LEm + c4PttLPm*LPm+ c4PttLfE*LfE + c4PttLfP*LfP
      !    + c4PttfEt*fEt + c4PttfEtt*fEtt + c4PttfPt*fPt+ c4PttfPtt*fPtt
      #Include interfaceAdeGdmOrder4.h 
   
      if( .false. .and. twilightZone.eq.1 )then
        write(*,'(" FACE: alpha,dt,d4,d1=",4e12.4)') alpha,dt,d4,d1
        write(*,'(" a0,a1,b0,b1=",4e12.4)') a0,a1,b0,b1
        write(*,'(" c2EtE,c2PtE,c2PttfP=",3e12.4)') c2EtE,c2PtE,c2PttfP
        write(*,'(" c4PttLE,c4PttE,c4PttEm,c4PttP,c4PttPm,c4PttfE=",6e12.4)') c4PttLE,c4PttE,c4PttEm,c4PttP,c4PttPm,c4PttfE
        write(*,'(" c4PttfP,c4PttLLE,c4PttLP,c4PttLEm,c4PttLPm,c4PttLfE,c4PttLfP=",7e12.4)') c4PttfP,c4PttLLE,c4PttLP,c4PttLEm,c4PttLPm,c4PttLfE,c4PttLfP
        write(*,'(" c4PttfEt,c4PttfEtt,c4PttfPt,c4PttfPtt=",4e12.4)') c4PttfEt,c4PttfEtt,c4PttfPt,c4PttfPtt
      end if
   
   
      ! Coeff of E in P.tt (4th order)
      c2PttEsum = c2PttEsum + c2PttE

      ! Coeff of LE in P.tt (4th order)
      c2PttLEsum  = c2PttLEsum + c2PttLE
   
      ! Coeff of LE in P.tt (4th order)
      c4PttLEsum  = c4PttLEsum + c4PttLE
   
      ! Coeff of LLE in P.tt
      c4PttLLEsum = c4PttLLEsum + c4PttLLE
   
      ! Coeff of LE and LLE in P.tttt
      c2PttttLEsum =c2PttttLEsum +c2PttttLE
      c2PttttLLEsum=c2PttttLLEsum+c2PttttLLE
   
      do n=0,nd-1
        pc = n + jv*nd 
        ec = ex +n
   
        pv   =  p(k1,k2,k3,pc)
        pvn  =  pn(k1,k2,k3,pc)
   
        ev    =  u(k1,k2,k3,ec)
        evn   =  un(k1,k2,k3,ec)
   
        ! Left: u1x,u1y, u1xx, u1yy, u1Lap (ex)
        !       v1x,v1y, v1xx, v1yy, v1Lap (ey) 
   
   
        ! These next derivatives may only be needed to order2, (use order 4 for testing TZ polynomials)
        #perl $ORDER=4;
        ! perl $ORDER=2;
        #If #FACE eq "LEFT"
          opEvalJacobianDerivatives(rsxy1,k1,k2,k3,aj1,1)
          ! uu1 in the next statement defines names of intermediate values
          evalSecondDerivs3d(rsxy1,aj1,p1,k1,k2,k3,pc,uu1,p1)
          ! write(*,'("FACE: p1x,p1y,p1xx,p1yy,p1Lap=",5(1pe12.4))') p1x,p1y,p1xx,p1yy,p1Lap
   
          LP  = (c1**2)*p1Lap
          evalSecondDerivs3d(rsxy1,aj1,p1n,k1,k2,k3,pc,uu1,p1n)
          ! write(*,'("FACE: p1nxx,p1nyy,p1nLap=",3e12.4)') p1nxx,p1nyy,p1nLap
          LPm = (c1**2)*p1nLap
   
          pvx  = p1x
          pvy  = p1y
          pvz  = p1z

          pvnx = p1nx
          pvny = p1ny
          pvnz = p1nz
   
        #Elif #FACE eq "RIGHT"
          opEvalJacobianDerivatives(rsxy2,k1,k2,k3,aj2,1)
          ! uu1 in the next statement defines names of intermediate values
          evalSecondDerivs3d(rsxy2,aj2,p2,k1,k2,k3,pc,uu2,p2)
          LP  = (c2**2)*p2Lap
          evalSecondDerivs3d(rsxy2,aj2,p2n,k1,k2,k3,pc,uu2,p2n)
          LPm = (c2**2)*p2nLap
   
          pvx  = p2x
          pvy  = p2y
          pvz  = p2z

          pvnx = p2nx
          pvny = p2ny
          pvnz = p2nz
   
        #Else
          write(*,'(" interface3d:ERROR: unknown FACE")')
          stop 7777
        #End
       
   
        ! Accumulate: SUM_m Pm,tt
        ptta = c4PttLE*LE(n) + c4PttE*ev + c4PttEm*evn + c4PttP*pv + c4PttPm*pvn + c4PttfE*fev(n) + c4PttfP*fpv(n,jv) \
           + c4PttLLE*LLE(n) + c4PttLP*LP + c4PttLEm*LEm(n) + c4PttLPm*LPm+ c4PttLfE*LfE(n) + c4PttLfP*LfP(n,jv) \
           + c4PttfEt*fEt(n) + c4PttfEtt*fEtt(n) + c4PttfPt*fPt(n,jv)+ c4PttfPtt*fPtt(n,jv)
   
        ! ---- Compute fp = P.tt 
        fp(n) = fp(n) + ptta
   
        ! ----- Compute fPttx = (P.tt).x , fPtty = (P.tt).y  (second order)
        pttxa = c2PttLE*LEx(n) + c2PttE*evx(n) + c2PttEm*evnx(n) + c2PttP*pvx + c2PttPm*pvnx + c2PttfE*fevx(n) + c2PttfP*fpvx(n,jv)
        pttya = c2PttLE*LEy(n) + c2PttE*evy(n) + c2PttEm*evny(n) + c2PttP*pvy + c2PttPm*pvny + c2PttfE*fevy(n) + c2PttfP*fpvy(n,jv)
        pttza = c2PttLE*LEz(n) + c2PttE*evz(n) + c2PttEm*evnz(n) + c2PttP*pvz + c2PttPm*pvnz + c2PttfE*fevz(n) + c2PttfP*fpvz(n,jv)
   
        ! write(debugFile,'(" pttxa,fevx(n),fpvx(n,jv)=",3(1pe10.2))') pttxa,fevx(n),fpvx(n,jv)
        ! write(debugFile,'(" LEx(n),evx(n),evnx(n),pvx,pvnx=",5(1pe10.2))') LEx(n),evx(n),evnx(n),pvx,pvnx

        fPttx(n) = fPttx(n) + pttxa
        fPtty(n) = fPtty(n) + pttya
        fPttz(n) = fPttz(n) + pttza
   
        ! ----- Compute fLPtt = L(P.tt) (second order)
        Lptta = c2PttLE*LLE(n) + c2PttE*LE(n) + c2PttEm*LEm(n) + c2PttP*LP + c2PttPm*LPm + c2PttfE*LfE(n) + c2PttfP*LfP(n,jv)
   
        fLPtt(n) = fLPtt(n) + Lptta
   
        ! ----- Compute fPtttt = P.tttt
        ptttta= c2PttttLE*LE(n) + c2PttttE*ev + c2PttttEm*evn + c2PttttP*pv + c2PttttPm*pvn + c2PttttfE*fev(n) + c2PttttfP*fpv(n,jv) \
           + c2PttttLLE*LLE(n) + c2PttttLP*LP + c2PttttLEm*LEm(n) + c2PttttLPm*LPm+ c2PttttLfE*LfE(n) + c2PttttLfP*LfP(n,jv) \
           + c2PttttfEt*fEt(n) + c2PttttfEtt*fEtt(n) + c2PttttfPt*fPt(n,jv)+ c2PttttfPtt*fPtt(n,jv)
   
        fPtttt(n) = fPtttt(n) + ptttta
   
        if( .false. .and. twilightZone.eq.1 )then
          write(*,'("")')
          write(*,'("DI4:FACE: k1,k2=",2i3," jv=",i2," n=",i2," ptta,ptte=",2e12.4)') k1,k2,jv,n,ptta,pevtt(n,jv)
          write(*,'("        : pttxa,pttxe=",2(1pe12.4)," pttya,pttye=",2(1pe12.4))') pttxa,pevttx(n,jv),pttya,pevtty(n,jv)
          write(*,'("        : ptttta,ptttte=",2(1pe12.4),/)') ptttta,pevtttt(n,jv)
   
          write(*,'(" c2PttLE,c2PttE,c2PttEm,c2PttP,c2PttPm,c2PttfE,c2PttfP=",7(1pe12.4))') c2PtLE,c2PtE,c2PtEm,c2PtP,c2PtPm,c2PtfE,c2PtfP
          write(*,'(" LEx,evx,evnx,pvx,pvnx,fevx,fpvx=",7(1pe12.4))') LEx(n),evx(n),evnx(n),pvx,pvnx,fevx(n),fpvx(n,jv)
          write(*,'(" LEy,evy,evny,pvy,pvny,fevy,fpvy=",7(1pe12.4))') LEy(n),evy(n),evny(n),pvy,pvny,fevy(n),fpvy(n,jv)
   
          ! write(*,'(" LE,LLE,LEm,LP,LPm=",5e12.4)') LE(n),LLE(n),LEm(n),LP,LPm
          ! write(*,'(" LfE,LfP,fEt,fEtt,fPt,fPtt=",6e12.4)') LfE(n),LfP(n,jv),fEt(n),fEtt(n),fPt(n,jv),fPtt(n,jv)
          ! write(*,'(" ev,evn,pv,pvn,fev,fpv=",6e12.4)')ev,evn,pv,pvn,fev(n),fpv(n,jv)
   
        end if
   
      end do ! end do n 

    end do ! jv
  !---------------------
  ! no dispersion
  !---------------------
  elseif( dispersionModel.eq.noDispersion .and. nonlinearModel.eq.noNonlinearModel) then
    ! do nothing
  end if
#endMacro

! -------------------------------------------------------------------------------
! Macro: Evaluate the TZ forcings GDM FOURTH-ORDER 3D
! -------------------------------------------------------------------------------
#beginMacro evalTZForcingUnified3dOrder4(xy,i1,i2,i3,dispersionModel,nonlinearModel,numberOfPolarizationVectors,numberOfAtomicLevels,\
                                         c,alphaP,pnec,prc,peptc,a0v,a1v,b0v,b1v,fpv,fpSum,fev,\
                                         LfE,fEt,fEtt,LfP,fPt,fPtt,fnv,fntv,fnttv,fntttv,\
                                         pevtt,pevttx,pevtty,pevttz,pevttL,pevtttt,fevx,fevy,fevz,fpvx,fpvy,fpvz,\
                                         pevttSum,pevttxSum,pevttySum,pevttzSum,pevttLSum,pevttttSum)

!----------------
! MLA
!----------------
if( dispersionModel.ne.noDispersion .and. nonlinearModel.ne.noNonlinearModel) then

  nce = pxc+nd*numberOfPolarizationVectors

  do n=0,nd-1

    fpSum(n)=0.
    pevttSum(n)=0.

    pevttxSum(n)=0.
    pevttySum(n)=0.
    pevttzSum(n)=0.

    pevttLSum(n)=0.
    pevttttSum(n)=0.

    petttSum=0.

    call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, es(n)   ) 
    call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, est(n)  )
    call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estt(n) )

    call ogderiv(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esx(n) )
    call ogderiv(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esy(n) )
    call ogderiv(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esz(n) )

    call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxx(n) )
    call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyy(n) )
    call ogderiv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, eszz(n) )

    call ogderiv(ep, 1,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estx(n) )
    call ogderiv(ep, 1,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esty(n) )
    call ogderiv(ep, 1,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estz(n) )

    call ogderiv(ep, 2,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttx(n) )
    call ogderiv(ep, 2,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estty(n) )
    call ogderiv(ep, 2,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttz(n) )

    call ogderiv(ep, 0,3,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxx(n) )
    call ogderiv(ep, 0,2,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxy(n) )
    call ogderiv(ep, 0,1,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxyy(n) )
    call ogderiv(ep, 0,0,3,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyyy(n) )

    call ogderiv(ep, 0,2,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxz(n) )
    call ogderiv(ep, 0,1,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxzz(n) )
    call ogderiv(ep, 0,0,2,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyyz(n) )
    call ogderiv(ep, 0,0,1,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyzz(n) )
    call ogderiv(ep, 0,0,0,3, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, eszzz(n) )

    call ogderiv(ep, 3,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttt(n) )
    call ogderiv(ep, 4,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estttt(n) )

    call ogderiv(ep, 1,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estxx(n) )
    call ogderiv(ep, 1,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estyy(n) )
    call ogderiv(ep, 1,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estzz(n) )


    call ogderiv(ep, 2,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttxx(n) )
    call ogderiv(ep, 2,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttyy(n) )
    call ogderiv(ep, 2,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttzz(n) )

    call ogderiv(ep, 0,4,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxxx(n) )
    call ogderiv(ep, 0,2,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxyy(n) )
    call ogderiv(ep, 0,0,4,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyyyy(n) )

    call ogderiv(ep, 0,0,0,4, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, eszzzz(n) )
    call ogderiv(ep, 0,2,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxzz(n) )
    call ogderiv(ep, 0,0,2,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyyzz(n) )

    ! L = Delta (removed c^2)
    esL  = ( esxx(n)   + esyy(n)   + eszz(n)  )
    estL = ( estxx(n)  + estyy(n)  + estzz(n) )
    esttL= ( esttxx(n) + esttyy(n) + esttzz(n) )

    esLx  = ( esxxx(n) + esxyy(n) + esxzz(n))
    esLy  = ( esxxy(n) + esyyy(n) + esyzz(n) )
    esLz  = ( esxxz(n) + esyyz(n) + eszzz(n) )

    ! L^2 : (xx + yy + zz)*( xx + yy + zz )
    esLL = ( esxxxx(n) + esyyyy(n) + eszzzz(n) + 2.*( esxxyy(n) + esxxzz(n) + esyyzz(n)) )

    do jv=0,numberOfPolarizationVectors-1
      ! The TZ component is offset by pxc
      pc = pxc + jv*nd
      call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pe(n)   )
      call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pet(n)  )
      call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pett(n) )
      call ogderiv(ep, 3,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettt(n) )
      call ogderiv(ep, 4,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petttt(n) )

      call ogderiv(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pex(n) )
      call ogderiv(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pey(n) )
      call ogderiv(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pez(n) )

      call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pexx(n) )
      call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, peyy(n) )
      call ogderiv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pezz(n) )

      call ogderiv(ep, 1,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petx(n) )
      call ogderiv(ep, 1,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pety(n) )
      call ogderiv(ep, 1,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petz(n) )

      call ogderiv(ep, 1,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petxx(n) )
      call ogderiv(ep, 1,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petyy(n) )
      call ogderiv(ep, 1,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petzz(n) )

      call ogderiv(ep, 2,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettx(n) )
      call ogderiv(ep, 2,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petty(n) )
      call ogderiv(ep, 2,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettz(n) )

      call ogderiv(ep, 2,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettxx(n) )
      call ogderiv(ep, 2,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettyy(n) )
      call ogderiv(ep, 2,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettzz(n) )

      peL  = ( pexx(n)   + peyy(n)   + pezz(n)  ) ! deleted c^2
      petL = ( petxx(n)  + petyy(n)  + petzz(n) )
      pettL= ( pettxx(n) + pettyy(n) + pettzz(n) )

      ! Normal TZ forcing for P_{n,jv} equation: 
      fpv(n,jv) = pett(n)   + b1v(jv)*pet(n)   + b0v(jv)*pe(n)
      fPt(n,jv) = pettt(n)  + b1v(jv)*pett(n)  + b0v(jv)*pet(n)
      fPtt(n,jv)= petttt(n) + b1v(jv)*pettt(n) + b0v(jv)*pett(n)
      LfP(n,jv) = pettL     + b1v(jv)*petL     + b0v(jv)*peL

      fpvx(n,jv)= pettx(n)  + b1v(jv)*petx(n)  + b0v(jv)*pex(n)
      fpvy(n,jv)= petty(n)  + b1v(jv)*pety(n)  + b0v(jv)*pey(n)
      fpvz(n,jv)= pettz(n)  + b1v(jv)*petz(n)  + b0v(jv)*pez(n)

      do na = 0,numberOfAtomicLevels-1
        call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0  )
        call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0t  )
        call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0tt  )

        call ogderiv(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0x  )
        call ogderiv(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0y  )
        call ogderiv(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0z  )

        call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0xx  )
        call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0yy  )
        call ogderiv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0zz  )

        fpv(n,jv) = fpv(n,jv) - pnec(jv,na)*q0*es(n) ! adding \Delta N*E
        fPt(n,jv) = fPt(n,jv) - pnec(jv,na)*q0t*es(n) - pnec(jv,na)*q0*est(n)
        fPtt(n,jv) = fPtt(n,jv) - pnec(jv,na)*q0tt*es(n)- 2.0*pnec(jv,na)*q0t*est(n) - pnec(jv,na)*q0*estt(n)
        LfP(n,jv) = LfP(n,jv) - pnec(jv,na)*(q0xx*es(n)+2.*q0x*esx(n)+q0*esxx(n) \
                                           + q0yy*es(n)+2.*q0y*esy(n)+q0*esyy(n) \
                                           + q0zz*es(n)+2.*q0z*esz(n)+q0*eszz(n))
        fpvx(n,jv) = fpvx(n,jv) - pnec(jv,na)*q0x*es(n) - pnec(jv,na)*q0*esx(n)
        fpvy(n,jv) = fpvy(n,jv) - pnec(jv,na)*q0y*es(n) - pnec(jv,na)*q0*esy(n)
        fpvz(n,jv) = fpvz(n,jv) - pnec(jv,na)*q0z*es(n) - pnec(jv,na)*q0*esz(n)
      enddo

      ! write(*,'(" n=",i4," LfP=",e10.4," pettL,petL,peL,esL,estL=",5e12.4)') n,LfP(n,jv),pettL,petL,peL,esL,estL
      ! write(*,'(" pe,pet,pett,pettt,petttt=",5e12.4)') pe(n),pet(n),pett(n),pettt(n),petttt(n)

      ! write(*,'("TZ: n,jv=",2i4," pex,pey,pexx,peyy=",4(1pe12.4))') n,jv,pex(n),pey(n),pexx(n),peyy(n)

      ! Save ptt for checking later
      pevtt(n,jv)=pett(n)
      pevttx(n,jv)=pettx(n)
      pevtty(n,jv)=petty(n)
      pevttz(n,jv)=pettz(n)

      pevtttt(n,jv)=petttt(n)
      pevttL(n,jv) = pettL
      pevttLSum(n) = pevttLSum(n)  + c**2*pettL ! added c^2 to be consistence with GDM
      pevttttSum(n)= pevttttSum(n) + petttt(n) 

      ! Keep some sums: 
      fpSum(n)   = fpSum(n)  + fpv(n,jv)
      pevttSum(n)  = pevttSum(n)  + pett(n) 
      pevttxSum(n) = pevttxSum(n) + pettx(n)
      pevttySum(n) = pevttySum(n) + petty(n)
      pevttzSum(n) = pevttzSum(n) + pettz(n)

      petttSum  = petttSum  + pettt(n) 


    end do 

    ! TZ forcing for E_{n} equation:
    ! E_tt - c^2 Delta E + alphaP*Ptt  = 
    fev(n) = estt(n)   - c**2*esL   + alphaP*pevttSum(n)
    fEt(n) = esttt(n)  - c**2*estL  + alphaP*petttSum
    fEtt(n)= estttt(n) - c**2*esttL + alphaP*pevttttSum(n)

    fevx(n) = esttx(n) - c**2*esLx   + alphaP*pevttxSum(n)
    fevy(n) = estty(n) - c**2*esLy   + alphaP*pevttySum(n)
    fevz(n) = esttz(n) - c**2*esLz   + alphaP*pevttzSum(n)
    

    ! write(*,'("--> fEtt=",e10.2," estttt,esttL,pettttSum=",3e10.2)')  fEtt(n),estttt(n),esttL,pettttSum
    LfE(n) = esttL     - c**2*esLL  + alphaP*pevttLSum(n)/c**2 ! added c**2
 end do !n

 !--------------------------------
  ! outside of dimension loop for N
  !--------------------------------
  do na=0,numberOfAtomicLevels-1
    ! na-th level
    call ogderiv(ep, 1,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0t )
    call ogderiv(ep, 2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0tt)
    call ogderiv(ep, 3,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0ttt)
    call ogderiv(ep, 4,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+na, q0tttt)
    ! initialize
    fnv(na)  = q0t ! forcing for \partial_tN_\ell = alpha_{\ell,k}N_k+\beta_{\ell,m}E\cdot\partial_tP_k
    fntv(na) = q0tt ! next derivative
    fnttv(na) = q0ttt
    fntttv(na) = q0tttt

    ! relaxation (alpha_{\ell,m})
    do jv=0,numberOfAtomicLevels-1
      call ogderiv(ep, 0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+jv, q0 )
      call ogderiv(ep, 1,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+jv, q0t)
      call ogderiv(ep, 2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+jv, q0tt)
      call ogderiv(ep, 3,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, nce+jv, q0ttt)
      fnv(na)  = fnv(na)  - prc(na,jv)*q0
      fntv(na) = fntv(na) - prc(na,jv)*q0t
      fnttv(na) = fnttv(na) - prc(na,jv)*q0tt
      fntttv(na) = fntttv(na) - prc(na,jv)*q0ttt
    enddo

    ! dot product (\beta_{\ell,k})
    do n=0,nd-1 ! loop over dim
      call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, es(n)   ) 
      call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, est(n)  )
      call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estt(n) )
      call ogderiv(ep, 3,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttt(n) )
      ! corresponding polarization vector
      do jv=0,numberOfPolarizationVectors-1 
        pc = pxc + jv*nd
        call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pe(n)   )
        call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pet(n)  )
        call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pett(n) )
        call ogderiv(ep, 3,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettt(n) )
        call ogderiv(ep, 4,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petttt(n) ) 
        fnv(na)  = fnv(na) - peptc(na,jv)*es(n)*pet(n)
        fntv(na) = fntv(na) - peptc(na,jv)*est(n)*pet(n) - peptc(na,jv)*es(n)*pett(n)
        fnttv(na) = fnttv(na) - peptc(na,jv)*estt(n)*pet(n) - 2.d0*peptc(na,jv)*est(n)*pett(n) - peptc(na,jv)*es(n)*pettt(n)
        fntttv(na) = fntttv(na) - peptc(na,jv)*esttt(n)*pet(n) \
                           - 3.d0*peptc(na,jv)*estt(n)*pett(n) \
                           - 3.d0*peptc(na,jv)*est(n)*pettt(n) \
                                - peptc(na,jv)*es(n)*petttt(n)
      enddo
    enddo

  enddo ! na

!----------------
! GDM
!----------------
elseif( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.noNonlinearModel) then
  do n=0,nd-1

    fpSum(n)=0.
    pevttSum(n)=0.

    pevttxSum(n)=0.
    pevttySum(n)=0.
    pevttzSum(n)=0.

    pevttLSum(n)=0.
    pevttttSum(n)=0.

    petttSum=0.

    call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, es(n)   ) 
    call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, est(n)  )
    call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estt(n) )

    call ogderiv(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esx(n) )
    call ogderiv(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esy(n) )
    call ogderiv(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esz(n) )

    call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxx(n) )
    call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyy(n) )
    call ogderiv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, eszz(n) )

    call ogderiv(ep, 1,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estx(n) )
    call ogderiv(ep, 1,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esty(n) )
    call ogderiv(ep, 1,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estz(n) )

    call ogderiv(ep, 2,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttx(n) )
    call ogderiv(ep, 2,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estty(n) )
    call ogderiv(ep, 2,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttz(n) )

    call ogderiv(ep, 0,3,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxx(n) )
    call ogderiv(ep, 0,2,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxy(n) )
    call ogderiv(ep, 0,1,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxyy(n) )
    call ogderiv(ep, 0,0,3,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyyy(n) )

    call ogderiv(ep, 0,2,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxz(n) )
    call ogderiv(ep, 0,1,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxzz(n) )
    call ogderiv(ep, 0,0,2,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyyz(n) )
    call ogderiv(ep, 0,0,1,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyzz(n) )
    call ogderiv(ep, 0,0,0,3, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, eszzz(n) )

    call ogderiv(ep, 3,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttt(n) )
    call ogderiv(ep, 4,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estttt(n) )

    call ogderiv(ep, 1,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estxx(n) )
    call ogderiv(ep, 1,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estyy(n) )
    call ogderiv(ep, 1,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, estzz(n) )


    call ogderiv(ep, 2,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttxx(n) )
    call ogderiv(ep, 2,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttyy(n) )
    call ogderiv(ep, 2,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esttzz(n) )

    call ogderiv(ep, 0,4,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxxx(n) )
    call ogderiv(ep, 0,2,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxyy(n) )
    call ogderiv(ep, 0,0,4,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyyyy(n) )

    call ogderiv(ep, 0,0,0,4, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, eszzzz(n) )
    call ogderiv(ep, 0,2,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esxxzz(n) )
    call ogderiv(ep, 0,0,2,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,ex+n, esyyzz(n) )

    ! L = c^2*Delta
    esL  = c**2*( esxx(n)   + esyy(n)   + eszz(n)  )
    estL = c**2*( estxx(n)  + estyy(n)  + estzz(n) )
    esttL= c**2*( esttxx(n) + esttyy(n) + esttzz(n) )

    esLx  = c**2*( esxxx(n) + esxyy(n) + esxzz(n))
    esLy  = c**2*( esxxy(n) + esyyy(n) + esyzz(n) )
    esLz  = c**2*( esxxz(n) + esyyz(n) + eszzz(n) )

    ! L^2 : (xx + yy + zz)*( xx + yy + zz )
    esLL = c**4*( esxxxx(n) + esyyyy(n) + eszzzz(n) + 2.*( esxxyy(n) + esxxzz(n) + esyyzz(n)) )

    do jv=0,numberOfPolarizationVectors-1
      ! The TZ component is offset by pxc
      pc = pxc + jv*nd
      call ogderiv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pe(n)   )
      call ogderiv(ep, 1,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pet(n)  )
      call ogderiv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pett(n) )
      call ogderiv(ep, 3,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettt(n) )
      call ogderiv(ep, 4,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petttt(n) )

      call ogderiv(ep, 0,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pex(n) )
      call ogderiv(ep, 0,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pey(n) )
      call ogderiv(ep, 0,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pez(n) )

      call ogderiv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pexx(n) )
      call ogderiv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, peyy(n) )
      call ogderiv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pezz(n) )

      call ogderiv(ep, 1,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petx(n) )
      call ogderiv(ep, 1,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pety(n) )
      call ogderiv(ep, 1,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petz(n) )

      call ogderiv(ep, 1,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petxx(n) )
      call ogderiv(ep, 1,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petyy(n) )
      call ogderiv(ep, 1,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petzz(n) )

      call ogderiv(ep, 2,1,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettx(n) )
      call ogderiv(ep, 2,0,1,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, petty(n) )
      call ogderiv(ep, 2,0,0,1, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettz(n) )

      call ogderiv(ep, 2,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettxx(n) )
      call ogderiv(ep, 2,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettyy(n) )
      call ogderiv(ep, 2,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,pc+n, pettzz(n) )

      peL  = c**2*( pexx(n)   + peyy(n)   + pezz(n)  )
      petL = c**2*( petxx(n)  + petyy(n)  + petzz(n) )
      pettL= c**2*( pettxx(n) + pettyy(n) + pettzz(n) )

      ! Normal TZ forcing for P_{n,jv} equation: 
      fpv(n,jv) = pett(n)   + b1v(jv)*pet(n)   + b0v(jv)*pe(n)   - a0v(jv)*es(n)   - a1v(jv)*est(n)
      fPt(n,jv) = pettt(n)  + b1v(jv)*pett(n)  + b0v(jv)*pet(n)  - a0v(jv)*est(n)  - a1v(jv)*estt(n)
      fPtt(n,jv)= petttt(n) + b1v(jv)*pettt(n) + b0v(jv)*pett(n) - a0v(jv)*estt(n) - a1v(jv)*esttt(n)
      LfP(n,jv) = pettL     + b1v(jv)*petL     + b0v(jv)*peL     - a0v(jv)*esL     - a1v(jv)*estL

      fpvx(n,jv)= pettx(n)  + b1v(jv)*petx(n)  + b0v(jv)*pex(n)  - a0v(jv)*esx(n)  - a1v(jv)*estx(n)
      fpvy(n,jv)= petty(n)  + b1v(jv)*pety(n)  + b0v(jv)*pey(n)  - a0v(jv)*esy(n)  - a1v(jv)*esty(n)
      fpvz(n,jv)= pettz(n)  + b1v(jv)*petz(n)  + b0v(jv)*pez(n)  - a0v(jv)*esz(n)  - a1v(jv)*estz(n)

      ! write(*,'(" n=",i4," LfP=",e10.4," pettL,petL,peL,esL,estL=",5e12.4)') n,LfP(n,jv),pettL,petL,peL,esL,estL
      ! write(*,'(" pe,pet,pett,pettt,petttt=",5e12.4)') pe(n),pet(n),pett(n),pettt(n),petttt(n)

      ! write(*,'("TZ: n,jv=",2i4," pex,pey,pexx,peyy=",4(1pe12.4))') n,jv,pex(n),pey(n),pexx(n),peyy(n)

      ! Save ptt for checking later
      pevtt(n,jv)=pett(n)
      pevttx(n,jv)=pettx(n)
      pevtty(n,jv)=petty(n)
      pevttz(n,jv)=pettz(n)

      pevtttt(n,jv)=petttt(n)
      pevttLSum(n) = pevttLSum(n)  + pettL
      pevttttSum(n)= pevttttSum(n) + petttt(n) 

      ! Keep some sums: 
      fpSum(n)   = fpSum(n)  + fpv(n,jv)
      pevttSum(n)  = pevttSum(n)  + pett(n) 
      pevttxSum(n) = pevttxSum(n) + pettx(n)
      pevttySum(n) = pevttySum(n) + petty(n)
      pevttzSum(n) = pevttzSum(n) + pettz(n)

      petttSum  = petttSum  + pettt(n) 


    end do 

    ! TZ forcing for E_{n} equation:
    ! E_tt - c^2 Delta E + alphaP*Ptt  = 
    fev(n) = estt(n)   - esL   + alphaP*pevttSum(n)
    fEt(n) = esttt(n)  - estL  + alphaP*petttSum
    fEtt(n)= estttt(n) - esttL + alphaP*pevttttSum(n)

    fevx(n) = esttx(n) - esLx   + alphaP*pevttxSum(n)
    fevy(n) = estty(n) - esLy   + alphaP*pevttySum(n)
    fevz(n) = esttz(n) - esLz   + alphaP*pevttzSum(n)
    

    ! write(*,'("--> fEtt=",e10.2," estttt,esttL,pettttSum=",3e10.2)')  fEtt(n),estttt(n),esttL,pettttSum
    LfE(n) = esttL     - esLL  + alphaP*pevttLSum(n)
    
 end do

!----------------
! no dispersion
!----------------
elseif( dispersionModel.eq.noDispersion .and. nonlinearModel.eq.noNonlinearModel) then
  do n=0,nd-1

    fpSum(n)=0.
    pevttSum(n)=0.

    pevttxSum(n)=0.
    pevttySum(n)=0.
    pevttzSum(n)=0.

    pevttLSum(n)=0.
    pevttttSum(n)=0.
  enddo
end if

#endMacro 

!-------------------------------------------------------------------------------------------
! Macro: Evaluate TZ forcing for dispersive equations in 3D 
!
! Output
!    fpv1(n,jv) : RHS To Pv_{n,jv} equation on domain 1
!    fpv2(n,jv) : RHS To Pv_{n,jv} equation on domain 2
!    fev1(n)    : RHS to E_{n} equation on domain 1
!    fev2(n)    : RHS to E_{n} equation on domain 2
!-------------------------------------------------------------------------------------------
#beginMacro getUnifiedDispersiveTZForcing3dOrder4(fpv1,fpv2,fev1,fev2,fnv1,fntv1,fnv2,fntv2)

  if( twilightZone.eq.1 )then
    evalTZForcingUnified3dOrder4(xy1,i1,i2,i3,dispersionModel1,nonlinearModel1,numberOfPolarizationVectors1,numberOfAtomicLevels1,\
                             c1,alphaP1,pnec1,prc1,peptc1,a0v1,a1v1,b0v1,b1v1,fpv1,fpSum1,fev1,\
                             LfE1,fEt1,fEtt1,LfP1,fPt1,fPtt1,fnv1,fntv1,fnttv1,fntttv1,\
                             pevtt1,pevttx1,pevtty1,pevttz1,pevttL1,pevtttt1,fevx1,fevy1,fevz1,fpvx1,fpvy1,fpvz1,\
                             pevttSum1,pevttxSum1,pevttySum1,pevttzSum1,pevttLSum1,pevttttSum1)

    evalTZForcingUnified3dOrder4(xy2,j1,j2,j3,dispersionModel2,nonlinearModel2,numberOfPolarizationVectors2,numberOfAtomicLevels2,\
                             c2,alphaP2,pnec2,prc2,peptc2,a0v2,a1v2,b0v2,b1v2,fpv2,fpSum2,fev2,\
                             LfE2,fEt2,fEtt2,LfP2,fPt2,fPtt2,fnv2,fntv2,fnttv2,fntttv2,\
                             pevtt2,pevttx2,pevtty2,pevttz2,pevttL2,pevtttt2,fevx2,fevy2,fevz2,fpvx2,fpvy2,fpvz2,\
                             pevttSum2,pevttxSum2,pevttySum2,pevttzSum2,pevttLSum2,pevttttSum2)
  end if

#endMacro

! --------------------------------------------------------------------------
! Macro: Evaluate the MLA jump conditions in 3D, order=4
! --------------------------------------------------------------------------
#beginMacro eval3dJumpUnifiedDispersiveOrder4()

! primary IC
 divE1 = u1x+v1y+w1z
 curlE1x = w1y-v1z
 curlE1y = u1z-w1x
 curlE1z = v1x-u1y
 nDotCurlE1=an1*curlE1x+an2*curlE1y+an3*curlE1z
 nDotLapE1 = an1*u1Lap + an2*v1Lap + an3*w1Lap

 divE2 = u2x+v2y+w2z
 curlE2x = w2y-v2z
 curlE2y = u2z-w2x
 curlE2z = v2x-u2y
 nDotCurlE2=an1*curlE2x+an2*curlE2y+an3*curlE2z
 nDotLapE2 = an1*u2Lap + an2*v2Lap + an3*w2Lap


 ! [ div(E) n + (I-nn^T)( curl(E)/mu ] = 0 
 f(0)= ( divE1*an1 + (curlE1x- nDotCurlE1*an1)/mu1 ) - ( divE2*an1 + (curlE2x- nDotCurlE2*an1)/mu2 )
 f(1)= ( divE1*an2 + (curlE1y- nDotCurlE1*an2)/mu1 ) - ( divE2*an2 + (curlE2y- nDotCurlE2*an2)/mu2 )
 f(2)= ( divE1*an3 + (curlE1z- nDotCurlE1*an3)/mu1 ) - ( divE2*an3 + (curlE2z- nDotCurlE2*an3)/mu2 )

 ! [  n n^T( Delta(E) /mu + (I-n n^T)( c^2 Delta(E) - alphaP*Ptt ] = 0 
 ! ---> [  c^2 Delta(E) + (1/mu -c^2) n n^T Delta(E) - (I-n n^T)( alphaP*Ptt ] = 0 
 cem1=(1.-1./eps1)/mu1
 cem2=(1.-1./eps2)/mu2

 nDotFp1 = an1*fp1(0) + an2*fp1(1) + an3*fp1(2)
 nDotFp2 = an1*fp2(0) + an2*fp2(1) + an3*fp2(2)

 f(3)= ( u1Lap/(epsmu1) + cem1*nDotLapE1*an1 - alphaP1*( fp1(0)-an1*nDotFp1) ) - \
       ( u2Lap/(epsmu2) + cem2*nDotLapE2*an1 - alphaP2*( fp2(0)-an1*nDotFp2) )
 f(4)= ( v1Lap/(epsmu1) + cem1*nDotLapE1*an2 - alphaP1*( fp1(1)-an2*nDotFp1) ) - \
       ( v2Lap/(epsmu2) + cem2*nDotLapE2*an2 - alphaP2*( fp2(1)-an2*nDotFp2) )
 f(5)= ( w1Lap/(epsmu1) + cem1*nDotLapE1*an3 - alphaP1*( fp1(2)-an3*nDotFp1) ) - \
       ( w2Lap/(epsmu2) + cem2*nDotLapE2*an3 - alphaP2*( fp2(2)-an3*nDotFp2) )

! compatibility conditions
 divLapE1 = u1xxx+v1xxy+w1xxz + u1xyy+v1yyy+w1yyz + u1xzz+v1yzz+w1zzz
 curlLapE1x = w1xxy-v1xxz + w1yyy-v1yyz + w1yzz-v1zzz
 curlLapE1y = u1xxz-w1xxx + u1yyz-w1xyy + u1zzz-w1xzz
 curlLapE1z = v1xxx-u1xxy + v1xyy-u1yyy + v1xzz-u1yzz
 nDotCurlLapE1=an1*curlLapE1x+an2*curlLapE1y+an3*curlLapE1z
 nDotLapSqE1 = an1*u1LapSq + an2*v1LapSq + an3*w1LapSq

 divLapE2 = u2xxx+v2xxy+w2xxz + u2xyy+v2yyy+w2yyz + u2xzz+v2yzz+w2zzz
 curlLapE2x = w2xxy-v2xxz + w2yyy-v2yyz + w2yzz-v2zzz
 curlLapE2y = u2xxz-w2xxx + u2yyz-w2xyy + u2zzz-w2xzz
 curlLapE2z = v2xxx-u2xxy + v2xyy-u2yyy + v2xzz-u2yzz
 nDotCurlLapE2=an1*curlLapE2x+an2*curlLapE2y+an3*curlLapE2z
 nDotLapSqE2 = an1*u2LapSq + an2*v2LapSq + an3*w2LapSq

 ! Compute curl(fPtt) and n^T curl(fPtt)
 curlfPttx1=fPtty1(2)-fPttz1(1) 
 curlfPtty1=fPttz1(0)-fPttx1(2)
 curlfPttz1=fPttx1(1)-fPtty1(0)
 nDotCurlfPtt1 = an1*curlfPttx1 + an2*curlfPtty1 + an3*curlfPttz1

 curlfPttx2=fPtty2(2)-fPttz2(1) 
 curlfPtty2=fPttz2(0)-fPttx2(2)
 curlfPttz2=fPttx2(1)-fPtty2(0)
 nDotCurlfPtt2 = an1*curlfPttx2 + an2*curlfPtty2 + an3*curlfPttz2

 ! Compute n^T( fLPtt) and n^T( fPtttt )
 nDotfLPtt1 = an1*fLPtt1(0)  + an2*fLPtt1(1)  + an3*fLPtt1(2)
 nDotfPtttt1= an1*fPtttt1(0) + an2*fPtttt1(1) + an3*fPtttt1(2) 

 nDotfLPtt2 = an1*fLPtt2(0)  + an2*fLPtt2(1)  + an3*fLPtt2(2)
 nDotfPtttt2= an1*fPtttt2(0) + an2*fPtttt2(1) + an3*fPtttt2(2) 

 ! [ c^2 div(Delta(E)) n + (1/mu)*( I - n n^T )( curl( c^2 Delta^2(E) -alphaP*P_tt) ] = 0
 !
 f(6)= ( ( divLapE1*an1 + (curlLapE1x- nDotCurlLapE1*an1)/mu1 )/(epsmu1) - (alphaP1/mu1)*( curlfPttx1 - an1*nDotCurlfPtt1 ) ) - \
       ( ( divLapE2*an1 + (curlLapE2x- nDotCurlLapE2*an1)/mu2 )/(epsmu2) - (alphaP2/mu2)*( curlfPttx2 - an1*nDotCurlfPtt2 ) )

 f(7)= ( ( divLapE1*an2 + (curlLapE1y- nDotCurlLapE1*an2)/mu1 )/(epsmu1) - (alphaP1/mu1)*( curlfPtty1 - an2*nDotCurlfPtt1 ) ) - \
       ( ( divLapE2*an2 + (curlLapE2y- nDotCurlLapE2*an2)/mu2 )/(epsmu2) - (alphaP2/mu2)*( curlfPtty2 - an2*nDotCurlfPtt2 ) )

 f(8)= ( ( divLapE1*an3 + (curlLapE1z- nDotCurlLapE1*an3)/mu1 )/(epsmu1) - (alphaP1/mu1)*( curlfPttz1 - an3*nDotCurlfPtt1 ) ) - \
       ( ( divLapE2*an3 + (curlLapE2z- nDotCurlLapE2*an3)/mu2 )/(epsmu2) - (alphaP2/mu2)*( curlfPttz2 - an3*nDotCurlfPtt2 ) )

 ! f(6)= ( ( divLapE1*an1 + (curlLapE1x- nDotCurlLapE1*an1)/mu1 )/(epsmu1) - (alphaP1/mu1)*( curlfPttx1 - an1*nDotCurlfPtt1 ) ) - \
 !       ( ( divLapE2*an1 + (curlLapE2x- nDotCurlLapE2*an1)/mu2 )/(epsmu2) - (alphaP2/mu2)*( curlfPttx2 - an1*nDotCurlfPtt2 ) )

 ! f(7)= ( ( divLapE1*an2 + (curlLapE1y- nDotCurlLapE1*an2)/mu1 )/(epsmu1) - (alphaP1/mu1)*( curlfPtty1 - an2*nDotCurlfPtt1 ) ) - \
 !       ( ( divLapE2*an2 + (curlLapE2y- nDotCurlLapE2*an2)/mu2 )/(epsmu2) - (alphaP2/mu2)*( curlfPtty2 - an2*nDotCurlfPtt2 ) )

 ! f(8)= ( ( divLapE1*an3 + (curlLapE1z- nDotCurlLapE1*an3)/mu1 )/(epsmu1) - (alphaP1/mu1)*( curlfPttz1 - an3*nDotCurlfPtt1 ) ) - \
 !       ( ( divLapE2*an3 + (curlLapE2z- nDotCurlLapE2*an3)/mu2 )/(epsmu2) - (alphaP2/mu2)*( curlfPttz2 - an3*nDotCurlfPtt2 ) )

 ! write(debugFile,'(" f(6)=",1pe10.2," curlfPttx1,nDotCurlfPtt1=",2(1pe10.2))') f(6),curlfPttx1,nDotCurlfPtt1
 ! write(debugFile,'(" f(7)=",1pe10.2," curlfPttx2,nDotCurlfPtt2=",2(1pe10.2))') f(7),curlfPttx2,nDotCurlfPtt2
 ! write(debugFile,'(" fPttx2=",3(1pe10.2))') fPttx2(0),fPttx2(1),fPttx2(2)
 ! write(debugFile,'(" fPtty2=",3(1pe10.2))') fPtty2(0),fPtty2(1),fPtty2(2)
 ! write(debugFile,'(" fPttz2=",3(1pe10.2))') fPttz2(0),fPttz2(1),fPttz2(2)
 ! write(debugFile,'(" divLapE1,curlLapE1x,nDotCurlLapE1=",3(1pe10.2))') divLapE1,curlLapE1x,nDotCurlLapE1
 ! write(debugFile,'(" divLapE2,curlLapE2x,nDotCurlLapE2=",3(1pe10.2))') divLapE2,curlLapE2x,nDotCurlLapE2


 ! [ (1/mu)* n n^T( c^2*Delta^2(E) - alphaP*Delta(P_tt) ) + (I-n n^T)( c^4 Delta^2(E) - c^2 alphaP*Delta(P_tt) - alphaP*P_tttt ]=0
 !  Note: fLptt = c^2*Delta( Ptt ) 
 ! f(9)= ( ( u1LapSq/(epsmu1) + cem1*nDotLapSqE1*an1 )/(epsmu1) - alphaP1*( (fLPtt1(0)-an1*nDotfLPtt1)/epsmu1 + fPtttt1(0)-an1*nDotfPtttt1  ) ) - \
 !       ( ( u2LapSq/(epsmu2) + cem2*nDotLapSqE2*an1 )/(epsmu2) - alphaP2*( (fLPtt2(0)-an1*nDotfLPtt2)/epsmu2 + fPtttt2(0)-an1*nDotfPtttt2  ) )
 !
 ! f(10)=( ( v1LapSq/(epsmu1) + cem1*nDotLapSqE1*an2 )/(epsmu1) - alphaP1*( (fLPtt1(1)-an2*nDotfLPtt1)/epsmu1 + fPtttt1(1)-an2*nDotfPtttt1  ) ) - \
 !       ( ( v2LapSq/(epsmu2) + cem2*nDotLapSqE2*an2 )/(epsmu2) - alphaP2*( (fLPtt2(1)-an2*nDotfLPtt2)/epsmu2 + fPtttt2(1)-an2*nDotfPtttt2  ) )
 !
 ! f(11)=( ( w1LapSq/(epsmu1) + cem1*nDotLapSqE1*an3 )/(epsmu1) - alphaP1*( (fLPtt1(2)-an3*nDotfLPtt1)/epsmu1 + fPtttt1(2)-an3*nDotfPtttt1  ) ) - \
 !       ( ( w2LapSq/(epsmu2) + cem2*nDotLapSqE2*an3 )/(epsmu2) - alphaP2*( (fLPtt2(2)-an3*nDotfLPtt2)/epsmu2 + fPtttt2(2)-an3*nDotfPtttt2  ) )

 ! AGAIN -- Oct 24, 2018:
 ! [ (1/mu)* n n^T( c^2 div(Delta(E) - alphaP* Delta(Ptt) )  + (I-n n^T)( c^4 Delta^2(E) - alphaP*{ c^2 *Delta(P_tt) + P_tttt} ]=0
 !  Note: fLptt = Delta( Ptt ) ! removed c^2
 ! nDotEP1 = (nDotLapSqE1/epsmu1 -alphaP1*nDotfLPtt1)/mu1
 ! nDotEP2 = (nDotLapSqE2/epsmu2 -alphaP2*nDotfLPtt2)/mu2

 ! f(9)= ( nDotEP1*an1 + (u1LapSq-nDotLapSqE1*an1)/epsmu1**2 -alphaP1*( (fLPtt1(0)-an1*nDotfLPtt1)/epsmu1 + fPtttt1(0)-an1*nDotfPtttt1 ) ) - \
 !       ( nDotEP2*an1 + (u2LapSq-nDotLapSqE2*an1)/epsmu2**2 -alphaP2*( (fLPtt2(0)-an1*nDotfLPtt2)/epsmu2 + fPtttt2(0)-an1*nDotfPtttt2 ) )
                                                                                             
 ! f(10)=( nDotEP1*an2 + (v1LapSq-nDotLapSqE1*an2)/epsmu1**2 -alphaP1*( (fLPtt1(1)-an2*nDotfLPtt1)/epsmu1 + fPtttt1(1)-an2*nDotfPtttt1 ) ) - \
 !       ( nDotEP2*an2 + (v2LapSq-nDotLapSqE2*an2)/epsmu2**2 -alphaP2*( (fLPtt2(1)-an2*nDotfLPtt2)/epsmu2 + fPtttt2(1)-an2*nDotfPtttt2 ) )
                                                                                             
 ! f(11)=( nDotEP1*an3 + (w1LapSq-nDotLapSqE1*an3)/epsmu1**2 -alphaP1*( (fLPtt1(2)-an3*nDotfLPtt1)/epsmu1 + fPtttt1(2)-an3*nDotfPtttt1 ) ) - \
 !       ( nDotEP2*an3 + (w2LapSq-nDotLapSqE2*an3)/epsmu2**2 -alphaP2*( (fLPtt2(2)-an3*nDotfLPtt2)/epsmu2 + fPtttt2(2)-an3*nDotfPtttt2 ) )

 !  Note: fLptt = c^2*Delta( Ptt ) 
 nDotEP1 = (nDotLapSqE1/epsmu1 -alphaP1*nDotfLPtt1*epsmu1)/mu1
 nDotEP2 = (nDotLapSqE2/epsmu2 -alphaP2*nDotfLPtt2*epsmu2)/mu2

 f(9)= ( nDotEP1*an1 + (u1LapSq-nDotLapSqE1*an1)/epsmu1**2 -alphaP1*( (fLPtt1(0)-an1*nDotfLPtt1) + fPtttt1(0)-an1*nDotfPtttt1 ) ) - \
       ( nDotEP2*an1 + (u2LapSq-nDotLapSqE2*an1)/epsmu2**2 -alphaP2*( (fLPtt2(0)-an1*nDotfLPtt2) + fPtttt2(0)-an1*nDotfPtttt2 ) )
                                                                                             
 f(10)=( nDotEP1*an2 + (v1LapSq-nDotLapSqE1*an2)/epsmu1**2 -alphaP1*( (fLPtt1(1)-an2*nDotfLPtt1) + fPtttt1(1)-an2*nDotfPtttt1 ) ) - \
       ( nDotEP2*an2 + (v2LapSq-nDotLapSqE2*an2)/epsmu2**2 -alphaP2*( (fLPtt2(1)-an2*nDotfLPtt2) + fPtttt2(1)-an2*nDotfPtttt2 ) )
                                                                                             
 f(11)=( nDotEP1*an3 + (w1LapSq-nDotLapSqE1*an3)/epsmu1**2 -alphaP1*( (fLPtt1(2)-an3*nDotfLPtt1) + fPtttt1(2)-an3*nDotfPtttt1 ) ) - \
       ( nDotEP2*an3 + (w2LapSq-nDotLapSqE2*an3)/epsmu2**2 -alphaP2*( (fLPtt2(2)-an3*nDotfLPtt2) + fPtttt2(2)-an3*nDotfPtttt2 ) )

 ! For testing extrap the 2nd ghost line :
 ! *e678*
 ! f(6)= u1(i1-2*is1,i2-2*is2,i3-2*is3,ex) - extrap5(u1,i1-is1,i2-is2,i3-is3,ex,is1,is2,is3)
 ! f(7)= u1(i1-2*is1,i2-2*is2,i3-2*is3,ey) - extrap5(u1,i1-is1,i2-is2,i3-is3,ey,is1,is2,is3)
 ! f(8)= u1(i1-2*is1,i2-2*is2,i3-2*is3,ez) - extrap5(u1,i1-is1,i2-is2,i3-is3,ez,is1,is2,is3)   

 ! f(9)= u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)-extrap5(u2,j1-js1,j2-js2,j3-js3,ex,js1,js2,js3)
 ! f(10)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)-extrap5(u2,j1-js1,j2-js2,j3-js3,ey,js1,js2,js3)
 ! f(11)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)-extrap5(u2,j1-js1,j2-js2,j3-js3,ez,js1,js2,js3)

 ! write(*,'("++++++++++++ --> before tz 3d:order4-c: f(solve)=",12e10.2)') (f(n),n=0,11)

 if( twilightZone.eq.1 )then

   x1=xy1(i1,i2,i3,0)
   y1=xy1(i1,i2,i3,1)
   z1=xy1(i1,i2,i3,2)
   call ogderiv(ep, 0,1,0,0, x1,y1,z1,t, ex, uex  )
   call ogderiv(ep, 0,0,1,0, x1,y1,z1,t, ex, uey  )
   call ogderiv(ep, 0,0,0,1, x1,y1,z1,t, ex, uez  )
   call ogderiv(ep, 0,2,0,0, x1,y1,z1,t, ex, uexx )
   call ogderiv(ep, 0,0,2,0, x1,y1,z1,t, ex, ueyy )
   call ogderiv(ep, 0,0,0,2, x1,y1,z1,t, ex, uezz )

   call ogderiv(ep, 0,1,0,0, x1,y1,z1,t, ey, vex  )
   call ogderiv(ep, 0,0,1,0, x1,y1,z1,t, ey, vey  )
   call ogderiv(ep, 0,0,0,1, x1,y1,z1,t, ey, vez  )
   call ogderiv(ep, 0,2,0,0, x1,y1,z1,t, ey, vexx )
   call ogderiv(ep, 0,0,2,0, x1,y1,z1,t, ey, veyy )
   call ogderiv(ep, 0,0,0,2, x1,y1,z1,t, ey, vezz )

   call ogderiv(ep, 0,1,0,0, x1,y1,z1,t, ez, wex  )
   call ogderiv(ep, 0,0,1,0, x1,y1,z1,t, ez, wey  )
   call ogderiv(ep, 0,0,0,1, x1,y1,z1,t, ez, wez  )
   call ogderiv(ep, 0,2,0,0, x1,y1,z1,t, ez, wexx )
   call ogderiv(ep, 0,0,2,0, x1,y1,z1,t, ez, weyy )
   call ogderiv(ep, 0,0,0,2, x1,y1,z1,t, ez, wezz )

   call ogderiv(ep, 0,3,0,0, x1,y1,z1,t, ex, uexxx )
   call ogderiv(ep, 0,0,3,0, x1,y1,z1,t, ex, ueyyy )
   call ogderiv(ep, 0,0,0,3, x1,y1,z1,t, ex, uezzz )
   call ogderiv(ep, 0,2,1,0, x1,y1,z1,t, ex, uexxy )
   call ogderiv(ep, 0,2,0,1, x1,y1,z1,t, ex, uexxz )
   call ogderiv(ep, 0,0,2,1, x1,y1,z1,t, ex, ueyyz )

   call ogderiv(ep, 0,1,2,0, x1,y1,z1,t, ex, uexyy ) 
   call ogderiv(ep, 0,1,0,2, x1,y1,z1,t, ex, uexzz )
   call ogderiv(ep, 0,0,1,2, x1,y1,z1,t, ex, ueyzz )


   call ogderiv(ep, 0,4,0,0, x1,y1,z1,t, ex, uexxxx )
   call ogderiv(ep, 0,2,2,0, x1,y1,z1,t, ex, uexxyy )
   call ogderiv(ep, 0,2,0,2, x1,y1,z1,t, ex, uexxzz )
   call ogderiv(ep, 0,0,4,0, x1,y1,z1,t, ex, ueyyyy )
   call ogderiv(ep, 0,0,2,2, x1,y1,z1,t, ex, ueyyzz )
   call ogderiv(ep, 0,0,0,4, x1,y1,z1,t, ex, uezzzz )

   call ogderiv(ep, 0,3,0,0, x1,y1,z1,t, ey, vexxx )
   call ogderiv(ep, 0,0,3,0, x1,y1,z1,t, ey, veyyy )
   call ogderiv(ep, 0,0,0,3, x1,y1,z1,t, ey, vezzz )
   call ogderiv(ep, 0,2,1,0, x1,y1,z1,t, ey, vexxy )
   call ogderiv(ep, 0,2,0,1, x1,y1,z1,t, ey, vexxz )
   call ogderiv(ep, 0,0,2,1, x1,y1,z1,t, ey, veyyz )

   call ogderiv(ep, 0,1,2,0, x1,y1,z1,t, ey, vexyy )
   call ogderiv(ep, 0,1,0,2, x1,y1,z1,t, ey, vexzz )
   call ogderiv(ep, 0,0,1,2, x1,y1,z1,t, ey, veyzz )

   call ogderiv(ep, 0,4,0,0, x1,y1,z1,t, ey, vexxxx )
   call ogderiv(ep, 0,2,2,0, x1,y1,z1,t, ey, vexxyy )
   call ogderiv(ep, 0,2,0,2, x1,y1,z1,t, ey, vexxzz )
   call ogderiv(ep, 0,0,4,0, x1,y1,z1,t, ey, veyyyy )
   call ogderiv(ep, 0,0,2,2, x1,y1,z1,t, ey, veyyzz )
   call ogderiv(ep, 0,0,0,4, x1,y1,z1,t, ey, vezzzz )

   call ogderiv(ep, 0,3,0,0, x1,y1,z1,t, ez, wexxx )
   call ogderiv(ep, 0,0,3,0, x1,y1,z1,t, ez, weyyy )
   call ogderiv(ep, 0,0,0,3, x1,y1,z1,t, ez, wezzz )
   call ogderiv(ep, 0,2,1,0, x1,y1,z1,t, ez, wexxy )
   call ogderiv(ep, 0,2,0,1, x1,y1,z1,t, ez, wexxz )
   call ogderiv(ep, 0,0,2,1, x1,y1,z1,t, ez, weyyz )

   call ogderiv(ep, 0,1,2,0, x1,y1,z1,t, ez, wexyy ) 
   call ogderiv(ep, 0,1,0,2, x1,y1,z1,t, ez, wexzz )
   call ogderiv(ep, 0,0,1,2, x1,y1,z1,t, ez, weyzz )

   call ogderiv(ep, 0,4,0,0, x1,y1,z1,t, ez, wexxxx )
   call ogderiv(ep, 0,2,2,0, x1,y1,z1,t, ez, wexxyy )
   call ogderiv(ep, 0,2,0,2, x1,y1,z1,t, ez, wexxzz )
   call ogderiv(ep, 0,0,4,0, x1,y1,z1,t, ez, weyyyy )
   call ogderiv(ep, 0,0,2,2, x1,y1,z1,t, ez, weyyzz )
   call ogderiv(ep, 0,0,0,4, x1,y1,z1,t, ez, wezzzz )


   ueLap = uexx+ueyy+uezz
   veLap = vexx+veyy+vezz
   weLap = wexx+weyy+wezz

   curlEex = wey-vez
   curlEey = uez-wex
   curlEez = vex-uey
   nDotCurlEe=an1*curlEex+an2*curlEey+an3*curlEez
   nDotLapEe=an1*ueLap+an2*veLap+an3*weLap

   f(0)= f(0) - ( (curlEex- nDotCurlEe*an1)*(1./mu1-1./mu2) )
   f(1)= f(1) - ( (curlEey- nDotCurlEe*an2)*(1./mu1-1./mu2) )
   f(2)= f(2) - ( (curlEez- nDotCurlEe*an3)*(1./mu1-1./mu2) )

   nDotPevttSum1 = an1*pevttSum1(0) + an2*pevttSum1(1) + an3*pevttSum1(2)
   nDotPevttSum2 = an1*pevttSum2(0) + an2*pevttSum2(1) + an3*pevttSum2(2)

   f(3)= f(3) - ( ueLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an1*(cem1-cem2) ) \
              + alphaP1*( pevttSum1(0) -an1*nDotPevttSum1 ) \
              - alphaP2*( pevttSum2(0) -an1*nDotPevttSum2 ) 

   f(4)= f(4) - ( veLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an2*(cem1-cem2) ) \
              + alphaP1*( pevttSum1(1) -an2*nDotPevttSum1 ) \
              - alphaP2*( pevttSum2(1) -an2*nDotPevttSum2 )

   f(5)= f(5) - ( weLap*(1./epsmu1-1./epsmu2) + nDotLapEe*an3*(cem1-cem2) )  \
              + alphaP1*( pevttSum1(2) -an3*nDotPevttSum1 ) \
              - alphaP2*( pevttSum2(2) -an3*nDotPevttSum2 )


   ueLapSq = uexxxx+ueyyyy+uezzzz +2.*( uexxyy+uexxzz+ueyyzz )
   veLapSq = vexxxx+veyyyy+vezzzz +2.*( vexxyy+vexxzz+veyyzz )
   weLapSq = wexxxx+weyyyy+wezzzz +2.*( wexxyy+wexxzz+weyyzz )

   curlLapEex = wexxy-vexxz + weyyy-veyyz + weyzz-vezzz
   curlLapEey = uexxz-wexxx + ueyyz-wexyy + uezzz-wexzz
   curlLapEez = vexxx-uexxy + vexyy-ueyyy + vexzz-ueyzz
   nDotCurlLapEe=an1*curlLapEex+an2*curlLapEey+an3*curlLapEez
   nDotLapSqEe=an1*ueLapSq+an2*veLapSq+an3*weLapSq

   curlPevttxSum1 = pevttySum1(2) - pevttzSum1(1)
   curlPevttySum1 = pevttzSum1(0) - pevttxSum1(2)
   curlPevttzSum1 = pevttxSum1(1) - pevttySum1(0)
   nDotCurlPevttSum1 = an1*curlPevttxSum1 + an2*curlPevttySum1 + an3*curlPevttzSum1

   ! print *, '1 pevtt x,y,z diff 0',fPttx1(0)-pevttxSum1(0),fPtty1(0)-pevttySum1(0),fPttz1(0)-pevttzSum1(0)
   ! print *, '1 pevtt x,y,z diff 1',fPttx1(1)-pevttxSum1(1),fPtty1(1)-pevttySum1(1),fPttz1(1)-pevttzSum1(1)
   ! print *, '1 pevtt x,y,z diff 2',fPttx1(2)-pevttxSum1(2),fPtty1(2)-pevttySum1(2),fPttz1(2)-pevttzSum1(2)

   ! print *, '2 pevtt x,y,z diff 0',fPttx2(0)-pevttxSum2(0),fPtty2(0)-pevttySum2(0),fPttz2(0)-pevttzSum2(0)
   ! print *, '2 pevtt x,y,z diff 1',fPttx2(1)-pevttxSum2(1),fPtty2(1)-pevttySum2(1),fPttz2(1)-pevttzSum2(1)
   ! print *, '2 pevtt x,y,z diff 2',fPttx2(2)-pevttxSum2(2),fPtty2(2)-pevttySum2(2),fPttz2(2)-pevttzSum2(2)

   curlPevttxSum2 = pevttySum2(2) - pevttzSum2(1)
   curlPevttySum2 = pevttzSum2(0) - pevttxSum2(2)
   curlPevttzSum2 = pevttxSum2(1) - pevttySum2(0)
   nDotCurlPevttSum2 = an1*curlPevttxSum2 + an2*curlPevttySum2 + an3*curlPevttzSum2

   divLapEe = uexxx+vexxy+wexxz + uexyy+veyyy+weyyz + uexzz+veyzz+wezzz

   ! f(6)= ( ( divLapE1*an1 + (curlLapE1x- nDotCurlLapE1*an1)/mu1 )/(epsmu1) - (alphaP1/mu1)*( curlfPttx1 - an1*nDotCurlfPtt1 ) ) - \
   !     ( ( divLapE2*an1 + (curlLapE2x- nDotCurlLapE2*an1)/mu2 )/(epsmu2) - (alphaP2/mu2)*( curlfPttx2 - an1*nDotCurlfPtt2 ) )

   ! f(7)= ( ( divLapE1*an2 + (curlLapE1y- nDotCurlLapE1*an2)/mu1 )/(epsmu1) - (alphaP1/mu1)*( curlfPtty1 - an2*nDotCurlfPtt1 ) ) - \
   !     ( ( divLapE2*an2 + (curlLapE2y- nDotCurlLapE2*an2)/mu2 )/(epsmu2) - (alphaP2/mu2)*( curlfPtty2 - an2*nDotCurlfPtt2 ) )

   ! f(8)= ( ( divLapE1*an3 + (curlLapE1z- nDotCurlLapE1*an3)/mu1 )/(epsmu1) - (alphaP1/mu1)*( curlfPttz1 - an3*nDotCurlfPtt1 ) ) - \
   !     ( ( divLapE2*an3 + (curlLapE2z- nDotCurlLapE2*an3)/mu2 )/(epsmu2) - (alphaP2/mu2)*( curlfPttz2 - an3*nDotCurlfPtt2 ) )


   f(6) = f(6) - ( divLapEe*an1*(1./epsmu1-1./epsmu2) + (curlLapEex- nDotCurlLapEe*an1)*(1./(mu1*epsmu1)-1./(mu2*epsmu2)) )  \
               + (alphaP1/mu1)*( curlPevttxSum1 - an1*nDotCurlPevttSum1 ) \
               - (alphaP2/mu2)*( curlPevttxSum2 - an1*nDotCurlPevttSum2 ) 

   ! print *, 'diff in eqn 6',divLapE1-divLapEe,curlLapE1x-curlLapEex,nDotCurlLapE1-nDotCurlLapEe,curlfPttx1-curlPevttxSum1,nDotCurlfPtt1-nDotCurlPevttSum1

   f(7) = f(7) - ( divLapEe*an2*(1./epsmu1-1./epsmu2) + (curlLapEey- nDotCurlLapEe*an2)*(1./(mu1*epsmu1)-1./(mu2*epsmu2)) ) \
               + (alphaP1/mu1)*( curlPevttySum1 - an2*nDotCurlPevttSum1 ) \
               - (alphaP2/mu2)*( curlPevttySum2 - an2*nDotCurlPevttSum2 ) 
   ! print *, 'diff in eqn 7',divLapE1-divLapEe,curlLapE1y-curlLapEey,nDotCurlLapE1-nDotCurlLapEe,curlfPtty1-curlPevttySum1,nDotCurlfPtt1-nDotCurlPevttSum1

   f(8) = f(8) - ( divLapEe*an3*(1./epsmu1-1./epsmu2) + (curlLapEez- nDotCurlLapEe*an3)*(1./(mu1*epsmu1)-1./(mu2*epsmu2)) ) \
               + (alphaP1/mu1)*( curlPevttzSum1 - an3*nDotCurlPevttSum1 ) \
               - (alphaP2/mu2)*( curlPevttzSum2 - an3*nDotCurlPevttSum2 ) 

   ! print *, 'diff in eqn 8',divLapE1-divLapEe,curlLapE1z-curlLapEez,nDotCurlLapE1-nDotCurlLapEe,curlfPttz1-curlPevttzSum1,nDotCurlfPtt1-nDotCurlPevttSum1

   nDotPevttLSum1 = an1*pevttLSum1(0)  + an2*pevttLSum1(1)  + an3*pevttLSum1(2)
   nDotPevttttSum1= an1*pevttttSum1(0) + an2*pevttttSum1(1) + an3*pevttttSum1(2)

   nDotPevttLSum2 = an1*pevttLSum2(0)  + an2*pevttLSum2(1)  + an3*pevttLSum2(2)
   nDotPevttttSum2= an1*pevttttSum2(0) + an2*pevttttSum2(1) + an3*pevttttSum2(2)
   
   ! 
   ! (nDotLapSqE1*an1/epsmu1 -alphaP1*nDotfLPtt1*an1)/mu1
   ! + (u1LapSq-nDotLapSqE1*an1)/epsmu1**2 
   !  -alphaP1*( (fLPtt1(0)-an1*nDotfLPtt1)/epsmu1 + fPtttt1(0)-an1*nDotfPtttt1 ) )
   ! without c^2 in fLPtt
   ! f(9) = f(9) - ( ueLapSq*(1./epsmu1**2-1./epsmu2**2) + nDotLapSqEe*an1*(cem1/epsmu1-cem2/epsmu2) ) \
   !             + alphaP1*(  (pevttLSum1(0)/epsmu1+cem1*an1*nDotPevttLSum1) + (pevttttSum1(0)-an1*nDotPevttttSum1) ) \
   !             - alphaP2*(  (pevttLSum2(0)/epsmu2+cem2*an1*nDotPevttLSum2) + (pevttttSum2(0)-an1*nDotPevttttSum2) )

   ! f(10)= f(10)- ( veLapSq*(1./epsmu1**2-1./epsmu2**2) + nDotLapSqEe*an2*(cem1/epsmu1-cem2/epsmu2) ) \
   !             + alphaP1*(  (pevttLSum1(1)/epsmu1+cem1*an2*nDotPevttLSum1) + (pevttttSum1(1)-an2*nDotPevttttSum1) ) \
   !             - alphaP2*(  (pevttLSum2(1)/epsmu2+cem2*an2*nDotPevttLSum2) + (pevttttSum2(1)-an2*nDotPevttttSum2) )

   ! f(11)= f(11)- ( weLapSq*(1./epsmu1**2-1./epsmu2**2) + nDotLapSqEe*an3*(cem1/epsmu1-cem2/epsmu2) )  \
   !             + alphaP1*(  (pevttLSum1(2)/epsmu1+cem1*an3*nDotPevttLSum1) + (pevttttSum1(2)-an3*nDotPevttttSum1) ) \
   !             - alphaP2*(  (pevttLSum2(2)/epsmu2+cem2*an3*nDotPevttLSum2) + (pevttttSum2(2)-an3*nDotPevttttSum2) )
    
   ! with c^2 in fLPtt
   f(9) = f(9) - ( ueLapSq*(1./epsmu1**2-1./epsmu2**2) + nDotLapSqEe*an1*(cem1/epsmu1-cem2/epsmu2) ) \
               + alphaP1*(  (pevttLSum1(0)+cem1*an1*nDotPevttLSum1*epsmu1) + (pevttttSum1(0)-an1*nDotPevttttSum1) ) \
               - alphaP2*(  (pevttLSum2(0)+cem2*an1*nDotPevttLSum2*epsmu2) + (pevttttSum2(0)-an1*nDotPevttttSum2) )

   f(10)= f(10)- ( veLapSq*(1./epsmu1**2-1./epsmu2**2) + nDotLapSqEe*an2*(cem1/epsmu1-cem2/epsmu2) ) \
               + alphaP1*(  (pevttLSum1(1)+cem1*an2*nDotPevttLSum1*epsmu1) + (pevttttSum1(1)-an2*nDotPevttttSum1) ) \
               - alphaP2*(  (pevttLSum2(1)+cem2*an2*nDotPevttLSum2*epsmu2) + (pevttttSum2(1)-an2*nDotPevttttSum2) )

   f(11)= f(11)- ( weLapSq*(1./epsmu1**2-1./epsmu2**2) + nDotLapSqEe*an3*(cem1/epsmu1-cem2/epsmu2) )  \
               + alphaP1*(  (pevttLSum1(2)+cem1*an3*nDotPevttLSum1*epsmu1) + (pevttttSum1(2)-an3*nDotPevttttSum1) ) \
               - alphaP2*(  (pevttLSum2(2)+cem2*an3*nDotPevttLSum2*epsmu2) + (pevttttSum2(2)-an3*nDotPevttttSum2) )
 end if

 ! write(*,'("++++++++++++ --> after tz 3d:order4-c: f(solve)=",12e10.2)') (f(n),n=0,11)

#endMacro


! ==========================================================================================
!   Evaluate the jump conditions (including compatibility) for the MLA interface equations 
! ==========================================================================================
#beginMacro evaluateUnifiedInterfaceEquations3dOrder4()

 ! Evaluate TZ forcing for dispersive equations in 3D 
 getUnifiedDispersiveTZForcing3dOrder4(fpv1,fpv2,fev1,fev2,fnv1,fntv1,fnv2,fntv2)

 evalDerivs3dOrder4()

 ! Store c^2*Delta(E) in a vector 
 LE1(0)=(c1**2)*u1Lap
 LE1(1)=(c1**2)*v1Lap
 LE1(2)=(c1**2)*w1Lap
 
 LE2(0)=(c2**2)*u2Lap
 LE2(1)=(c2**2)*v2Lap
 LE2(2)=(c2**2)*w2Lap
 
 ! Store L^2(E) 
 LLE1(0)=(c1**4)*u1LapSq
 LLE1(1)=(c1**4)*v1LapSq
 LLE1(2)=(c1**4)*w1LapSq
 
 LLE2(0)=(c2**4)*u2LapSq
 LLE2(1)=(c2**4)*v2LapSq
 LLE2(2)=(c2**4)*w2LapSq
 
 ! Store (LE).x an (LE).y 
 LEx1(0) = (c1**2)*( u1xxx + u1xyy + u1xzz )
 LEx1(1) = (c1**2)*( v1xxx + v1xyy + v1xzz )
 LEx1(2) = (c1**2)*( w1xxx + w1xyy + w1xzz )

 LEy1(0) = (c1**2)*( u1xxy + u1yyy + u1yzz )
 LEy1(1) = (c1**2)*( v1xxy + v1yyy + v1yzz )
 LEy1(2) = (c1**2)*( w1xxy + w1yyy + w1yzz )

 LEz1(0) = (c1**2)*( u1xxz + u1yyz + u1zzz )
 LEz1(1) = (c1**2)*( v1xxz + v1yyz + v1zzz )
 LEz1(2) = (c1**2)*( w1xxz + w1yyz + w1zzz )

 LEx2(0) = (c2**2)*( u2xxx + u2xyy + u2xzz )
 LEx2(1) = (c2**2)*( v2xxx + v2xyy + v2xzz )
 LEx2(2) = (c2**2)*( w2xxx + w2xyy + w2xzz )

 LEy2(0) = (c2**2)*( u2xxy + u2yyy + u2yzz )
 LEy2(1) = (c2**2)*( v2xxy + v2yyy + v2yzz )
 LEy2(2) = (c2**2)*( w2xxy + w2yyy + w2yzz )

 LEz2(0) = (c2**2)*( u2xxz + u2yyz + u2zzz )
 LEz2(1) = (c2**2)*( v2xxz + v2yyz + v2zzz )
 LEz2(2) = (c2**2)*( w2xxz + w2yyz + w2zzz )

 ! We also need derivatives at the old time:
 ! These next derivatives may only be needed to order2, but use order 4 for now so exact for degree 4
 #perl $ORDER=4;
 ! perl $ORDER=2;
 opEvalJacobianDerivatives(rsxy1,i1,i2,i3,aj1,1)
 evalSecondDerivs3d(rsxy1,aj1,u1n,i1,i2,i3,ex,uu1,u1n)
 evalSecondDerivs3d(rsxy1,aj1,u1n,i1,i2,i3,ey,vv1,v1n)
 evalSecondDerivs3d(rsxy1,aj1,u1n,i1,i2,i3,ez,ww1,w1n)
 ! Here are c^2*Delta(E) at the old time: 
 LE1m(0) = (c1**2)*u1nLap
 LE1m(1) = (c1**2)*v1nLap
 LE1m(2) = (c1**2)*w1nLap

 opEvalJacobianDerivatives(rsxy2,j1,j2,j3,aj2,1)
 evalSecondDerivs3d(rsxy2,aj2,u2n,j1,j2,j3,ex,uu2,u2n)
 evalSecondDerivs3d(rsxy2,aj2,u2n,j1,j2,j3,ey,vv2,v2n)
 evalSecondDerivs3d(rsxy2,aj2,u2n,j1,j2,j3,ez,ww2,w2n)
 LE2m(0) = (c2**2)*u2nLap
 LE2m(1) = (c2**2)*v2nLap
 LE2m(2) = (c2**2)*w2nLap

 evx1(0) = u1x
 evx1(1) = v1x
 evx1(2) = w1x

 evy1(0) = u1y
 evy1(1) = v1y 
 evy1(2) = w1y 

 evz1(0) = u1z
 evz1(1) = v1z 
 evz1(2) = w1z 

 evnx1(0) = u1nx
 evnx1(1) = v1nx
 evnx1(2) = w1nx

 evny1(0) = u1ny
 evny1(1) = v1ny 
 evny1(2) = w1ny 

 evnz1(0) = u1nz
 evnz1(1) = v1nz 
 evnz1(2) = w1nz 

 evx2(0) = u2x
 evx2(1) = v2x
 evx2(2) = w2x

 evy2(0) = u2y
 evy2(1) = v2y 
 evy2(2) = w2y 

 evz2(0) = u2z
 evz2(1) = v2z 
 evz2(2) = w2z 

 evnx2(0) = u2nx
 evnx2(1) = v2nx
 evnx2(2) = w2nx

 evny2(0) = u2ny
 evny2(1) = v2ny 
 evny2(2) = w2ny 

 evnz2(0) = u2nz
 evnz2(1) = v2nz 
 evnz2(2) = w2nz 

 ! eval dispersive forcings for domain 1
 getUnifiedDispersiveForcing3dOrder4(LEFT,i1,i2,i3, fp1, fpv1,fev1,p1,p1n,p1m,q1,q1n,q1m,u1,u1n,u1m, \
              dispersionModel1,nonlinearModel1,numberOfPolarizationVectors1,numberOfAtomicLevels1,alphaP1,c1,\
              c2PttEsum1,c2PttLEsum1,c4PttLEsum1,c4PttLLEsum1,c2PttttLEsum1,c2PttttLLEsum1,a0v1,a1v1,\
              pnec1,prc1,peptc1,b0v1,b1v1,LE1,LLE1,LE1m,LfE1,LfP1,fEt1,fEtt1,fPt1,fPtt1,fnv1,fntv1,\
              pevtt1,pevttx1,pevtty1,pevttz1,pevttL1,pevtttt1,evx1,evy1,evz1,evnx1,evny1,evnz1,\
              fevx1,fevy1,fevz1,fpvx1,fpvy1,fpvz1,LEx1,LEy1,LEz1,fPttx1,fPtty1,fPttz1,fLPtt1,fPtttt1) 

 ! eval dispersive forcings for domain 2
 getUnifiedDispersiveForcing3dOrder4(RIGHT,j1,j2,j3, fp2, fpv2,fev2,p2,p2n,p2m, q2,q2n,q2m,u2,u2n,u2m, \
              dispersionModel2,nonlinearModel2,numberOfPolarizationVectors2,numberOfAtomicLevels2,alphaP2,c2,\
              c2PttEsum2,c2PttLEsum2,c4PttLEsum2,c4PttLLEsum2,c2PttttLEsum2,c2PttttLLEsum2,a0v2,a1v2,\
              pnec2,prc2,peptc2,b0v2,b1v2,LE2,LLE2,LE2m,LfE2,LfP2,fEt2,fEtt2,fPt2,fPtt2,fnv2,fntv2,\
              pevtt2,pevttx2,pevtty2,pevttz2,pevttL2,pevtttt2,evx2,evy2,evz2,evnx2,evny2,evnz2,\
              fevx2,fevy2,fevz2,fpvx2,fpvy2,fpvz2,LEx2,LEy2,LEz2,fPttx2,fPtty2,fPttz2,fLPtt2,fPtttt2) 


 eval3dJumpUnifiedDispersiveOrder4()

 ! eval3dJumpOrder4New()
#endMacro

 ! ==========================================================================================
!         Fill in FIRST SET of 6 interface equations -- GDM/MLA -- .
! 
! Input:
!  e0,e1,e2,e3,e4,e5 : equation numbers (0,1,2,3,4,5) or (6,7,8,9,10,11)
! ==========================================================================================
#beginMacro fillUnifiedDispersiveEquations3dOrder4a(e0,e1,e2,e3,e4,e5)

 ! Equations 0,1,2:
 ! [ div(E) n  + (1/mu)*(I-n n^T)( curl(E) ] = 0 

 ! Equation 0: 
 ! (u.x+v.y+w.z)*an1 + ( w1y-v1z - nDotCurlE1*an1)/mu1

 ! n^T curl(E) = an1*((Ez)_y - (Ey)_z+ + an2*((Ex)_z - (Ez)_x) + an3*((Ey)_x - (Ex)_y )
 ! (I - n n^T) curl(E) = [ (Ez)_y - (Ey)_z - an1*( n^T curl(E) )]   [ a11*Ex + a12*Ey + a13*Ez ]
 !                       [ (Ex)_z - (Ez)_x - an2*(             )] = [ a21*Ex + a22*Ey + a23*Ez ]
 !                       [ (Ey)_x - (Ex)_y - an3*(             )]   [ a31*Ex + a32*Et + a33*Ez ]
 !       a11 =     -an1*( an2*Dz -an3*Dy )
 !       a12 = -Dz -an1*( an3*Dx -an1*Dz )
 !       a13 =  Dy -an1*( an1*Dy -an2*Dz )

 ! ------- Equation 0 ------
 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e0,0) = ( c1m1x*an1 + (           - (c1m1z*an2-c1m1y*an3)*an1 )/mu1 ) ! coeff of u1(-1)
 a12(e0,1) = ( c1m1y*an1 + (    -c1m1z - (c1m1x*an3-c1m1z*an1)*an1 )/mu1 ) ! coeff of v1(-1)
 a12(e0,2) = ( c1m1z*an1 + ( c1m1y     - (c1m1y*an1-c1m1x*an2)*an1 )/mu1 ) ! coeff of w1(-1)

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e0,3) = ( c1m2x*an1 + (           - (c1m2z*an2-c1m2y*an3)*an1 )/mu1 ) ! coeff of u1(-1)
 a12(e0,4) = ( c1m2y*an1 + (    -c1m2z - (c1m2x*an3-c1m2z*an1)*an1 )/mu1 ) ! coeff of v1(-1)
 a12(e0,5) = ( c1m2z*an1 + ( c1m2y     - (c1m2y*an1-c1m2x*an2)*an1 )/mu1 ) ! coeff of w1(-1)

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e0,6) =-( c2m1x*an1 + (           - (c2m1z*an2-c2m1y*an3)*an1 )/mu2 ) ! coeff of u2(-1)
 a12(e0,7) =-( c2m1y*an1 + (    -c2m1z - (c2m1x*an3-c2m1z*an1)*an1 )/mu2 ) ! coeff of v2(-1)
 a12(e0,8) =-( c2m1z*an1 + ( c2m1y     - (c2m1y*an1-c2m1x*an2)*an1 )/mu2 ) ! coeff of w2(-1)

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e0,9) =-( c2m2x*an1 + (           - (c2m2z*an2-c2m2y*an3)*an1 )/mu2 ) ! coeff of u2(-1)
 a12(e0,10)=-( c2m2y*an1 + (    -c2m2z - (c2m2x*an3-c2m2z*an1)*an1 )/mu2 ) ! coeff of v2(-1)
 a12(e0,11)=-( c2m2z*an1 + ( c2m2y     - (c2m2y*an1-c2m2x*an2)*an1 )/mu2 ) ! coeff of w2(-1)

 ! *bug* fixed wdh July 6, 2019 "+" to minus "-" 8 places 
 ! Equation 1:
 ! (u.x+v.y+w.z)*an2 + ( u1z-w1x - nDotCurlE1*an2)/mu1
 a12(e1,0) = ( c1m1x*an2 + ( c1m1z     - (c1m1z*an2-c1m1y*an3)*an2 )/mu1 ) ! coeff of u1(-1)
 a12(e1,1) = ( c1m1y*an2 + (           - (c1m1x*an3-c1m1z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
 a12(e1,2) = ( c1m1z*an2 + (    -c1m1x - (c1m1y*an1-c1m1x*an2)*an2 )/mu1 ) ! coeff of w1(-1)

 a12(e1,3) = ( c1m2x*an2 + ( c1m2z     - (c1m2z*an2-c1m2y*an3)*an2 )/mu1 ) ! coeff of u1(-1)
 a12(e1,4) = ( c1m2y*an2 + (           - (c1m2x*an3-c1m2z*an1)*an2 )/mu1 ) ! coeff of v1(-1)
 a12(e1,5) = ( c1m2z*an2 + (    -c1m2x - (c1m2y*an1-c1m2x*an2)*an2 )/mu1 ) ! coeff of w1(-1)

 a12(e1,6) =-( c2m1x*an2 + ( c2m1z     - (c2m1z*an2-c2m1y*an3)*an2 )/mu2 ) ! coeff of u2(-1)
 a12(e1,7) =-( c2m1y*an2 + (           - (c2m1x*an3-c2m1z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
 a12(e1,8) =-( c2m1z*an2 + (    -c2m1x - (c2m1y*an1-c2m1x*an2)*an2 )/mu2 ) ! coeff of w2(-1)

 a12(e1,9) =-( c2m2x*an2 + ( c2m2z     - (c2m2z*an2-c2m2y*an3)*an2 )/mu2 ) ! coeff of u2(-1)
 a12(e1,10)=-( c2m2y*an2 + (           - (c2m2x*an3-c2m2z*an1)*an2 )/mu2 ) ! coeff of v2(-1)
 a12(e1,11)=-( c2m2z*an2 + (    -c2m2x - (c2m2y*an1-c2m2x*an2)*an2 )/mu2 ) ! coeff of w2(-1)

 ! Equation 2:
 ! (u.x+v.y+w.z)*an3 + ( v1x-u1y - nDotCurlE1*an3)/mu1
 a12(e2,0) = ( c1m1x*an3 + (    -c1m1y - (c1m1z*an2-c1m1y*an3)*an3 )/mu1 ) ! coeff of u1(-1)
 a12(e2,1) = ( c1m1y*an3 + ( c1m1x     - (c1m1x*an3-c1m1z*an1)*an3 )/mu1 ) ! coeff of v1(-1)
 a12(e2,2) = ( c1m1z*an3 + (           - (c1m1y*an1-c1m1x*an2)*an3 )/mu1 ) ! coeff of w1(-1)

 a12(e2,3) = ( c1m2x*an3 + (    -c1m2y - (c1m2z*an2-c1m2y*an3)*an3 )/mu1 ) ! coeff of u1(-1)
 a12(e2,4) = ( c1m2y*an3 + ( c1m2x     - (c1m2x*an3-c1m2z*an1)*an3 )/mu1 ) ! coeff of v1(-1)
 a12(e2,5) = ( c1m2z*an3 + (           - (c1m2y*an1-c1m2x*an2)*an3 )/mu1 ) ! coeff of w1(-1)

 a12(e2,6) =-( c2m1x*an3 + (    -c2m1y - (c2m1z*an2-c2m1y*an3)*an3 )/mu2 ) ! coeff of u2(-1)
 a12(e2,7) =-( c2m1y*an3 + ( c2m1x     - (c2m1x*an3-c2m1z*an1)*an3 )/mu2 ) ! coeff of v2(-1)
 a12(e2,8) =-( c2m1z*an3 + (           - (c2m1y*an1-c2m1x*an2)*an3 )/mu2 ) ! coeff of w2(-1)

 a12(e2,9) =-( c2m2x*an3 + (    -c2m2y - (c2m2z*an2-c2m2y*an3)*an3 )/mu2 ) ! coeff of u2(-1)
 a12(e2,10)=-( c2m2y*an3 + ( c2m2x     - (c2m2x*an3-c2m2z*an1)*an3 )/mu2 ) ! coeff of v2(-1)
 a12(e2,11)=-( c2m2z*an3 + (           - (c2m2y*an1-c2m2x*an2)*an3 )/mu2 ) ! coeff of w2(-1)

 !  (1/mu1)* u1Lap*ani*anj + (delta_ij -ani*anj)*( u1Lap/epsmu1 - alphaP1*( c4PttLEsum1*u1Lap + c4PttLLEsum1*u1LapSq )= ...

 ! ------------- Equation 3 -------------
 !  u1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an1 + (delta_ij -an1*anj)*( - alphaP1*( c4PttLEsum1*u1Lap + c4PttLLEsum1*u1LapSq ) )

 ! coeffs not used in MLA
 if( nonlinearModel1.ne.noNonlinearModel )then
   c4PttLEsum1 = 0.
   c4PttLLEsum1 = 0. 
 endif
 if( nonlinearModel2.ne.noNonlinearModel )then
   c4PttLEsum2 = 0.
   c4PttLLEsum2 = 0.
 endif

 eqn1Coeffm1 = -alphaP1*( (c4PttLEsum1/epsmu1)*clap1m1+(c4PttLLEsum1/epsmu1**2)*clapSq1m1 )
 eqn1Coeffm2 = -alphaP1*( (c4PttLEsum1/epsmu1)*clap1m2+(c4PttLLEsum1/epsmu1**2)*clapSq1m2 )

 eqn2Coeffm1 = -alphaP2*( (c4PttLEsum2/epsmu2)*clap2m1+(c4PttLLEsum2/epsmu2**2)*clapSq2m1 )
 eqn2Coeffm2 = -alphaP2*( (c4PttLEsum2/epsmu2)*clap2m2+(c4PttLLEsum2/epsmu2**2)*clapSq2m2 )

 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e3,0) = ( clap1m1/(epsmu1) + cem1*( an1*clap1m1                         )*an1 + (1.-an1*an1)*( eqn1Coeffm1 ) )
 a12(e3,1) = (                    cem1*(             an2*clap1m1             )*an1 + (  -an1*an2)*( eqn1Coeffm1 ) )
 a12(e3,2) = (                    cem1*(                         an3*clap1m1 )*an1 + (  -an1*an3)*( eqn1Coeffm1 ) )
                                                                                                                   
 ! coeffs of u1(-2), v1(-2), w1(-2)                                                                                
 a12(e3,3) = ( clap1m2/(epsmu1) + cem1*( an1*clap1m2                         )*an1 + (1.-an1*an1)*( eqn1Coeffm2 ) )
 a12(e3,4) = (                    cem1*(             an2*clap1m2             )*an1 + (  -an1*an2)*( eqn1Coeffm2 ) )
 a12(e3,5) = (                    cem1*(                         an3*clap1m2 )*an1 + (  -an1*an3)*( eqn1Coeffm2 ) )
                                                                                                                   
 ! coeffs of u2(-1), v2(-1), w2(-1)                                                                                
 a12(e3,6) =-( clap2m1/(epsmu2) + cem2*( an1*clap2m1                         )*an1 + (1.-an1*an1)*( eqn2Coeffm1 ) )
 a12(e3,7) =-(                    cem2*(             an2*clap2m1             )*an1 + (  -an1*an2)*( eqn2Coeffm1 ) )
 a12(e3,8) =-(                    cem2*(                         an3*clap2m1 )*an1 + (  -an1*an3)*( eqn2Coeffm1 ) )
                                                                                                                   
 ! coeffs of u2(-2), v2(-2), w2(-2)                                                                                
 a12(e3,9) =-( clap2m2/(epsmu2) + cem2*( an1*clap2m2                         )*an1 + (1.-an1*an1)*( eqn2Coeffm2 ) ) 
 a12(e3,10)=-(                    cem2*(             an2*clap2m2             )*an1 + (  -an1*an2)*( eqn2Coeffm2 ) )
 a12(e3,11)=-(                    cem2*(                         an3*clap2m2 )*an1 + (  -an1*an3)*( eqn2Coeffm2 ) )

 ! Equations 4,5,6
 !  [ (1/mu)* Lap(E) n n^T + (I-n n^T)* c^2 Lap(E) - alphaP*Ptt ) ] = 0 

 ! -------------- Equation 4 -----------------
 !  v1Lap/(epsmu1) + cem1*( an1*u1Lap + an2*v1Lap + an3*w1Lap )*an2+ (delta_ij -an2*anj)*( - alphaP1*( c4PttLEsum1*u1Lap + c4PttLLEsum1*u1LapSq ) )

 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e4,0) = (                    cem1*( an1*clap1m1                         )*an2 + (  -an2*an1)*( eqn1Coeffm1 ) )  
 a12(e4,1) = ( clap1m1/(epsmu1) + cem1*(             an2*clap1m1             )*an2 + (1.-an2*an2)*( eqn1Coeffm1 ) )
 a12(e4,2) = (                    cem1*(                         an3*clap1m1 )*an2 + (  -an2*an3)*( eqn1Coeffm1 ) )
                                                                                                                   
 ! coeffs of u1(-2), v1(-2), w1(-2)                                                                                
 a12(e4,3) = (                    cem1*( an1*clap1m2                         )*an2 + (  -an2*an1)*( eqn1Coeffm2 ) ) 
 a12(e4,4) = ( clap1m2/(epsmu1) + cem1*(             an2*clap1m2             )*an2 + (1.-an2*an2)*( eqn1Coeffm2 ) )
 a12(e4,5) = (                    cem1*(                         an3*clap1m2 )*an2 + (  -an2*an3)*( eqn1Coeffm2 ) )
                                                                                                                   
 ! coeffs of u2(-1), v2(-1), w2(-1)                                                                                
 a12(e4,6) =-(                    cem2*( an1*clap2m1                         )*an2 + (  -an2*an1)*( eqn2Coeffm1 ) ) 
 a12(e4,7) =-( clap2m1/(epsmu2) + cem2*(             an2*clap2m1             )*an2 + (1.-an2*an2)*( eqn2Coeffm1 ) )
 a12(e4,8) =-(                    cem2*(                         an3*clap2m1 )*an2 + (  -an2*an3)*( eqn2Coeffm1 ) )
                                                                                                                   
 ! coeffs of u2(-2), v2(-2), w2(-2)                                                                                
 a12(e4,9) =-(                    cem2*( an1*clap2m2                         )*an2 + (  -an2*an1)*( eqn2Coeffm2 ) ) 
 a12(e4,10)=-( clap2m2/(epsmu2) + cem2*(             an2*clap2m2             )*an2 + (1.-an2*an2)*( eqn2Coeffm2 ) )
 a12(e4,11)=-(                    cem2*(                         an3*clap2m2 )*an2 + (  -an2*an3)*( eqn2Coeffm2 ) )

 ! ----- Equation 5 -----------

 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e5,0) = (                    cem1*( an1*clap1m1                         )*an3 + (  -an3*an1)*( eqn1Coeffm1 ) )
 a12(e5,1) = (                    cem1*(             an2*clap1m1             )*an3 + (  -an3*an2)*( eqn1Coeffm1 ) )
 a12(e5,2) = ( clap1m1/(epsmu1) + cem1*(                         an3*clap1m1 )*an3 + (1.-an3*an3)*( eqn1Coeffm1 ) )
                                                                                                                   
 ! coeffs of u1(-2), v1(-2), w1(-2)                                                                                
 a12(e5,3) = (                    cem1*( an1*clap1m2                         )*an3 + (  -an3*an1)*( eqn1Coeffm2 ) )
 a12(e5,4) = (                    cem1*(             an2*clap1m2             )*an3 + (  -an3*an2)*( eqn1Coeffm2 ) )
 a12(e5,5) = ( clap1m2/(epsmu1) + cem1*(                         an3*clap1m2 )*an3 + (1.-an3*an3)*( eqn1Coeffm2 ) )
                                                                                                                   
 ! coeffs of u2(-1), v2(-1), w2(-1)                                                                                
 a12(e5,6) =-(                    cem2*( an1*clap2m1                         )*an3 + (  -an3*an1)*( eqn2Coeffm1 ) )
 a12(e5,7) =-(                    cem2*(             an2*clap2m1             )*an3 + (  -an3*an2)*( eqn2Coeffm1 ) )
 a12(e5,8) =-( clap2m1/(epsmu2) + cem2*(                         an3*clap2m1 )*an3 + (1.-an3*an3)*( eqn2Coeffm1 ) )
                                                                                                                   
 ! coeffs of u2(-2), v2(-2), w2(-2)                                                                                
 a12(e5,9) =-(                    cem2*( an1*clap2m2                         )*an3 + (  -an3*an1)*( eqn2Coeffm2 ) )
 a12(e5,10)=-(                    cem2*(             an2*clap2m2             )*an3 + (  -an3*an2)*( eqn2Coeffm2 ) )
 a12(e5,11)=-( clap2m2/(epsmu2) + cem2*(                         an3*clap2m2 )*an3 + (1.-an3*an3)*( eqn2Coeffm2 ) )
#endMacro


! ==========================================================================================
!         Fill in SECOND SET of 6 interface equations -- GDM/MLA -- .
! 
! Input:
!  e0,e1,e2,e3,e4,e5 : equation numbers  (6,7,8,9,10,11)
! ==========================================================================================
#beginMacro fillUnifiedDispersiveEquations3dOrder4b(e0,e1,e2,e3,e4,e5)

 ! Equations 0,1,2:
 ! [ c^2 Lap(div(E)) n + (1/mu)*(I-n n^T)( c^2 Lap(curl(E)) - alphaP*curl(Ptt) ) ] = 0 
 ! curl(E) = [ (Ez)_y - (Ey)_z,  (Ex)_z - (Ez)_x , (Ey)_x - (Ex)_y ]

 ! n^T curl(E) = an1*((Ez)_y - (Ey)_z+ + an2*((Ex)_z - (Ez)_x) + an3*((Ey)_x - (Ex)_y )
 ! (I - n n^T) curl(E) = [ (Ez)_y - (Ey)_z - an1*( n^T curl(E) )]   [ a11*Ex + a12*Ey + a13*Ez ]
 !                       [ (Ex)_z - (Ez)_x - an2*(             )] = [ a21*Ex + a22*Ey + a23*Ez ]
 !                       [ (Ey)_x - (Ex)_y - an3*(             )]   [ a31*Ex + a32*Et + a33*Ez ]
 !       a11 =     -an1*( an2*Dz -an3*Dy )
 !       a12 = -Dz -an1*( an3*Dx -an1*Dz )
 !       a13 =  Dy -an1*( an1*Dy -an2*Dz )
 if( nonlinearModel1.ne.noNonlinearModel )then
   c2PttLEsum1 = 0.
   c2PttEsum1 = 0.
 endif
 if( nonlinearModel2.ne.noNonlinearModel )then
   c2PttLEsum2 = 0.
   c2PttEsum2 = 0.
 endif

 eqn1Coeff  = (1./epsmu1 -alphaP1*c2PttLEsum1/epsmu1)/mu1
 eqn1Coeffb = -alphaP1*c2PttEsum1/mu1 

 eqn2Coeff  = (1./epsmu2 -alphaP2*c2PttLEsum2/epsmu2)/mu2 
 eqn2Coeffb = -alphaP2*c2PttEsum2/mu2 

 ! ------- Equation 6 ------
 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e0,0) = ( clap1m1x*an1/epsmu1 + (              - (clap1m1z*an2-clap1m1y*an3)*an1 )*eqn1Coeff + (           - (c1m1z*an2-c1m1y*an3)*an1 )*eqn1Coeffb )
 a12(e0,1) = ( clap1m1y*an1/epsmu1 + (    -clap1m1z - (clap1m1x*an3-clap1m1z*an1)*an1 )*eqn1Coeff + (    -c1m1z - (c1m1x*an3-c1m1z*an1)*an1 )*eqn1Coeffb )
 a12(e0,2) = ( clap1m1z*an1/epsmu1 + ( clap1m1y     - (clap1m1y*an1-clap1m1x*an2)*an1 )*eqn1Coeff + ( c1m1y     - (c1m1y*an1-c1m1x*an2)*an1 )*eqn1Coeffb )

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e0,3) = ( clap1m2x*an1/epsmu1 + (              - (clap1m2z*an2-clap1m2y*an3)*an1 )*eqn1Coeff + (           - (c1m2z*an2-c1m2y*an3)*an1 )*eqn1Coeffb ) 
 a12(e0,4) = ( clap1m2y*an1/epsmu1 + (    -clap1m2z - (clap1m2x*an3-clap1m2z*an1)*an1 )*eqn1Coeff + (    -c1m2z - (c1m2x*an3-c1m2z*an1)*an1 )*eqn1Coeffb ) 
 a12(e0,5) = ( clap1m2z*an1/epsmu1 + ( clap1m2y     - (clap1m2y*an1-clap1m2x*an2)*an1 )*eqn1Coeff + ( c1m2y     - (c1m2y*an1-c1m2x*an2)*an1 )*eqn1Coeffb ) 

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e0,6) =-( clap2m1x*an1/epsmu2 + (              - (clap2m1z*an2-clap2m1y*an3)*an1 )*eqn2Coeff + (           - (c2m1z*an2-c2m1y*an3)*an1 )*eqn2Coeffb )
 a12(e0,7) =-( clap2m1y*an1/epsmu2 + (    -clap2m1z - (clap2m1x*an3-clap2m1z*an1)*an1 )*eqn2Coeff + (    -c2m1z - (c2m1x*an3-c2m1z*an1)*an1 )*eqn2Coeffb )
 a12(e0,8) =-( clap2m1z*an1/epsmu2 + ( clap2m1y     - (clap2m1y*an1-clap2m1x*an2)*an1 )*eqn2Coeff + ( c2m1y     - (c2m1y*an1-c2m1x*an2)*an1 )*eqn2Coeffb )

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e0,9) =-( clap2m2x*an1/epsmu2 + (              - (clap2m2z*an2-clap2m2y*an3)*an1 )*eqn2Coeff + (           - (c2m2z*an2-c2m2y*an3)*an1 )*eqn2Coeffb ) 
 a12(e0,10)=-( clap2m2y*an1/epsmu2 + (    -clap2m2z - (clap2m2x*an3-clap2m2z*an1)*an1 )*eqn2Coeff + (    -c2m2z - (c2m2x*an3-c2m2z*an1)*an1 )*eqn2Coeffb ) 
 a12(e0,11)=-( clap2m2z*an1/epsmu2 + ( clap2m2y     - (clap2m2y*an1-clap2m2x*an2)*an1 )*eqn2Coeff + ( c2m2y     - (c2m2y*an1-c2m2x*an2)*an1 )*eqn2Coeffb ) 

 ! ----- Equation 7 -----
 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e1,0) = ( clap1m1x*an2/epsmu1 + ( clap1m1z     - (clap1m1z*an2-clap1m1y*an3)*an2 )*eqn1Coeff + ( c1m1z     - (c1m1z*an2-c1m1y*an3)*an2 )*eqn1Coeffb ) 
 a12(e1,1) = ( clap1m1y*an2/epsmu1 + (              - (clap1m1x*an3-clap1m1z*an1)*an2 )*eqn1Coeff + (           - (c1m1x*an3-c1m1z*an1)*an2 )*eqn1Coeffb ) 
 a12(e1,2) = ( clap1m1z*an2/epsmu1 + (    -clap1m1x - (clap1m1y*an1-clap1m1x*an2)*an2 )*eqn1Coeff + (    -c1m1x - (c1m1y*an1-c1m1x*an2)*an2 )*eqn1Coeffb ) 

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e1,3) = ( clap1m2x*an2/epsmu1 + ( clap1m2z     - (clap1m2z*an2-clap1m2y*an3)*an2 )*eqn1Coeff + ( c1m2z     - (c1m2z*an2-c1m2y*an3)*an2 )*eqn1Coeffb ) 
 a12(e1,4) = ( clap1m2y*an2/epsmu1 + (              - (clap1m2x*an3-clap1m2z*an1)*an2 )*eqn1Coeff + (           - (c1m2x*an3-c1m2z*an1)*an2 )*eqn1Coeffb ) 
 a12(e1,5) = ( clap1m2z*an2/epsmu1 + (    -clap1m2x - (clap1m2y*an1-clap1m2x*an2)*an2 )*eqn1Coeff + (    -c1m2x - (c1m2y*an1-c1m2x*an2)*an2 )*eqn1Coeffb ) 

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e1,6) =-( clap2m1x*an2/epsmu2 + ( clap2m1z     - (clap2m1z*an2-clap2m1y*an3)*an2 )*eqn2Coeff + ( c2m1z     - (c2m1z*an2-c2m1y*an3)*an2 )*eqn2Coeffb ) 
 a12(e1,7) =-( clap2m1y*an2/epsmu2 + (              - (clap2m1x*an3-clap2m1z*an1)*an2 )*eqn2Coeff + (           - (c2m1x*an3-c2m1z*an1)*an2 )*eqn2Coeffb ) 
 a12(e1,8) =-( clap2m1z*an2/epsmu2 + (    -clap2m1x - (clap2m1y*an1-clap2m1x*an2)*an2 )*eqn2Coeff + (    -c2m1x - (c2m1y*an1-c2m1x*an2)*an2 )*eqn2Coeffb ) 

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e1,9) =-( clap2m2x*an2/epsmu2 + ( clap2m2z     - (clap2m2z*an2-clap2m2y*an3)*an2 )*eqn2Coeff + ( c2m2z     - (c2m2z*an2-c2m2y*an3)*an2 )*eqn2Coeffb ) 
 a12(e1,10)=-( clap2m2y*an2/epsmu2 + (              - (clap2m2x*an3-clap2m2z*an1)*an2 )*eqn2Coeff + (           - (c2m2x*an3-c2m2z*an1)*an2 )*eqn2Coeffb ) 
 a12(e1,11)=-( clap2m2z*an2/epsmu2 + (    -clap2m2x - (clap2m2y*an1-clap2m2x*an2)*an2 )*eqn2Coeff + (    -c2m2x - (c2m2y*an1-c2m2x*an2)*an2 )*eqn2Coeffb ) 

 ! ----- Equation 8 -----

 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e2,0) = ( clap1m1x*an3/epsmu1 + (    -clap1m1y - (clap1m1z*an2-clap1m1y*an3)*an3 )*eqn1Coeff + (    -c1m1y - (c1m1z*an2-c1m1y*an3)*an3 )*eqn1Coeffb ) 
 a12(e2,1) = ( clap1m1y*an3/epsmu1 + ( clap1m1x     - (clap1m1x*an3-clap1m1z*an1)*an3 )*eqn1Coeff + ( c1m1x     - (c1m1x*an3-c1m1z*an1)*an3 )*eqn1Coeffb ) 
 a12(e2,2) = ( clap1m1z*an3/epsmu1 + (              - (clap1m1y*an1-clap1m1x*an2)*an3 )*eqn1Coeff + (           - (c1m1y*an1-c1m1x*an2)*an3 )*eqn1Coeffb ) 

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e2,3) = ( clap1m2x*an3/epsmu1 + (    -clap1m2y - (clap1m2z*an2-clap1m2y*an3)*an3 )*eqn1Coeff + (    -c1m2y - (c1m2z*an2-c1m2y*an3)*an3 )*eqn1Coeffb ) 
 a12(e2,4) = ( clap1m2y*an3/epsmu1 + ( clap1m2x     - (clap1m2x*an3-clap1m2z*an1)*an3 )*eqn1Coeff + ( c1m2x     - (c1m2x*an3-c1m2z*an1)*an3 )*eqn1Coeffb ) 
 a12(e2,5) = ( clap1m2z*an3/epsmu1 + (              - (clap1m2y*an1-clap1m2x*an2)*an3 )*eqn1Coeff + (           - (c1m2y*an1-c1m2x*an2)*an3 )*eqn1Coeffb ) 

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e2,6) =-( clap2m1x*an3/epsmu2 + (    -clap2m1y - (clap2m1z*an2-clap2m1y*an3)*an3 )*eqn2Coeff + (    -c2m1y - (c2m1z*an2-c2m1y*an3)*an3 )*eqn2Coeffb ) 
 a12(e2,7) =-( clap2m1y*an3/epsmu2 + ( clap2m1x     - (clap2m1x*an3-clap2m1z*an1)*an3 )*eqn2Coeff + ( c2m1x     - (c2m1x*an3-c2m1z*an1)*an3 )*eqn2Coeffb ) 
 a12(e2,8) =-( clap2m1z*an3/epsmu2 + (              - (clap2m1y*an1-clap2m1x*an2)*an3 )*eqn2Coeff + (           - (c2m1y*an1-c2m1x*an2)*an3 )*eqn2Coeffb ) 

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e2,9) =-( clap2m2x*an3/epsmu2 + (    -clap2m2y - (clap2m2z*an2-clap2m2y*an3)*an3 )*eqn2Coeff + (    -c2m2y - (c2m2z*an2-c2m2y*an3)*an3 )*eqn2Coeffb ) 
 a12(e2,10)=-( clap2m2y*an3/epsmu2 + ( clap2m2x     - (clap2m2x*an3-clap2m2z*an1)*an3 )*eqn2Coeff + ( c2m2x     - (c2m2x*an3-c2m2z*an1)*an3 )*eqn2Coeffb ) 
 a12(e2,11)=-( clap2m2z*an3/epsmu2 + (              - (clap2m2y*an1-clap2m2x*an2)*an3 )*eqn2Coeff + (           - (c2m2y*an1-c2m2x*an2)*an3 )*eqn2Coeffb ) 

 ! Equations 9,10,11
 !  [ (1/mu)*n^T( c^2*LapSqE - alphaP*Lap(Ptt) )n  + (I-n n^T)*( c^4LapSq(E) - c^2*alphaP*Lap(Ptt) -alphaP*Ptttt ]
 !  Note: LE = c^2*Delta(E)

 ! first term:  (c^2*LapSqE - alphaP*Lap(Ptt) )/mu 

 ! not used in MLA
 if( nonlinearModel1.ne.noNonlinearModel )then
   c2PttEsum1 = 0.
   c2PttLEsum1 = 0.
   c2PttttLEsum1 = 0.
   c2PttttLLEsum1 = 0.
 endif

 coeffLap1a   =               -alphaP1*( c2PttEsum1         )/mu1
 coeffLapSq1a = 1./epsmu1/mu1 -alphaP1*( c2PttLEsum1/epsmu1 )/mu1  
 ! second term:  c^4LapSq(E) - alphaP*( c^2*Lap(Ptt) + Ptttt )
 coeffLap1b   =               -alphaP1*( c2PttEsum1  + c2PttttLEsum1  )/epsmu1
 coeffLapSq1b = 1./epsmu1**2  -alphaP1*( c2PttLEsum1 + c2PttttLLEsum1 )/epsmu1**2
 
 eqn1Coeff1a = clap1m1*coeffLap1a + clapSq1m1*coeffLapSq1a
 eqn1Coeff1b = clap1m1*coeffLap1b + cLapSq1m1*coeffLapSq1b 

 eqn1Coeff2a = clap1m2*coeffLap1a + clapSq1m2*coeffLapSq1a
 eqn1Coeff2b = clap1m2*coeffLap1b + cLapSq1m2*coeffLapSq1b 

 ! not used in MLA
 if( nonlinearModel2.ne.noNonlinearModel )then
   c2PttEsum2 = 0.
   c2PttLEsum2 = 0.
   c2PttttLEsum2 = 0.
   c2PttttLLEsum2 = 0.
 endif


 ! first term:  c^2*LapSqE - alphaP*Delta(Ptt)
 coeffLap2a   =               -alphaP2*( c2PttEsum2         )/mu2
 coeffLapSq2a = 1./epsmu2/mu2 -alphaP2*( c2PttLEsum2/epsmu2 )/mu2   
 ! second term:  c^4LapSq(E) - alphaP*( c^2*Lap(Ptt) + Ptttt )
 coeffLap2b   =               -alphaP2*( c2PttEsum2  + c2PttttLEsum2  )/epsmu2
 coeffLapSq2b = 1./epsmu2**2  -alphaP2*( c2PttLEsum2 + c2PttttLLEsum2 )/epsmu2**2

 eqn2Coeff1a = clap2m1*coeffLap2a + clapSq2m1*coeffLapSq2a 
 eqn2Coeff1b = clap2m1*coeffLap2b + cLapSq2m1*coeffLapSq2b 

 eqn2Coeff2a = clap2m2*coeffLap2a + clapSq2m2*coeffLapSq2a 
 eqn2Coeff2b = clap2m2*coeffLap2b + cLapSq2m2*coeffLapSq2b

 ! -----Equation 9-----
 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e3,0) = ( an1*an1*eqn1Coeff1a + (1.-an1*an1)*eqn1Coeff1b )
 a12(e3,1) = ( an1*an2*eqn1Coeff1a + (  -an1*an2)*eqn1Coeff1b )
 a12(e3,2) = ( an1*an3*eqn1Coeff1a + (  -an1*an3)*eqn1Coeff1b )

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e3,3) = ( an1*an1*eqn1Coeff2a + (1.-an1*an1)*eqn1Coeff2b )
 a12(e3,4) = ( an1*an2*eqn1Coeff2a + (  -an1*an2)*eqn1Coeff2b )
 a12(e3,5) = ( an1*an3*eqn1Coeff2a + (  -an1*an3)*eqn1Coeff2b )

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e3,6) =-( an1*an1*eqn2Coeff1a + (1.-an1*an1)*eqn2Coeff1b )
 a12(e3,7) =-( an1*an2*eqn2Coeff1a + (  -an1*an2)*eqn2Coeff1b )
 a12(e3,8) =-( an1*an3*eqn2Coeff1a + (  -an1*an3)*eqn2Coeff1b )

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e3,9) =-( an1*an1*eqn2Coeff2a + (1.-an1*an1)*eqn2Coeff2b )
 a12(e3,10)=-( an1*an2*eqn2Coeff2a + (  -an1*an2)*eqn2Coeff2b )
 a12(e3,11)=-( an1*an3*eqn2Coeff2a + (  -an1*an3)*eqn2Coeff2b )

 ! -----Equation 10 -----
 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e4,0) = ( an2*an1*eqn1Coeff1a + (  -an2*an1)*eqn1Coeff1b )
 a12(e4,1) = ( an2*an2*eqn1Coeff1a + (1.-an2*an2)*eqn1Coeff1b )
 a12(e4,2) = ( an2*an3*eqn1Coeff1a + (  -an2*an3)*eqn1Coeff1b )

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e4,3) = ( an2*an1*eqn1Coeff2a + (  -an2*an1)*eqn1Coeff2b )
 a12(e4,4) = ( an2*an2*eqn1Coeff2a + (1.-an2*an2)*eqn1Coeff2b )
 a12(e4,5) = ( an2*an3*eqn1Coeff2a + (  -an2*an3)*eqn1Coeff2b )

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e4,6) =-( an2*an1*eqn2Coeff1a + (  -an2*an1)*eqn2Coeff1b )
 a12(e4,7) =-( an2*an2*eqn2Coeff1a + (1.-an2*an2)*eqn2Coeff1b )
 a12(e4,8) =-( an2*an3*eqn2Coeff1a + (  -an2*an3)*eqn2Coeff1b )

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e4,9) =-( an2*an1*eqn2Coeff2a + (  -an2*an1)*eqn2Coeff2b )
 a12(e4,10)=-( an2*an2*eqn2Coeff2a + (1.-an2*an2)*eqn2Coeff2b )
 a12(e4,11)=-( an2*an3*eqn2Coeff2a + (  -an2*an3)*eqn2Coeff2b )

 ! -----Equation 11 -----
 ! coeffs of u1(-1), v1(-1), w1(-1)
 a12(e5,0) = ( an3*an1*eqn1Coeff1a + (  -an3*an1)*eqn1Coeff1b )
 a12(e5,1) = ( an3*an2*eqn1Coeff1a + (  -an3*an2)*eqn1Coeff1b )
 a12(e5,2) = ( an3*an3*eqn1Coeff1a + (1.-an3*an3)*eqn1Coeff1b )

 ! coeffs of u1(-2), v1(-2), w1(-2)
 a12(e5,3) = ( an3*an1*eqn1Coeff2a + (  -an3*an1)*eqn1Coeff2b )
 a12(e5,4) = ( an3*an2*eqn1Coeff2a + (  -an3*an2)*eqn1Coeff2b )
 a12(e5,5) = ( an3*an3*eqn1Coeff2a + (1.-an3*an3)*eqn1Coeff2b )

 ! coeffs of u2(-1), v2(-1), w2(-1) 
 a12(e5,6) =-( an3*an1*eqn2Coeff1a + (  -an3*an1)*eqn2Coeff1b )
 a12(e5,7) =-( an3*an2*eqn2Coeff1a + (  -an3*an2)*eqn2Coeff1b )
 a12(e5,8) =-( an3*an3*eqn2Coeff1a + (1.-an3*an3)*eqn2Coeff1b )

 ! coeffs of u2(-2), v2(-2), w2(-2)
 a12(e5,9) =-( an3*an1*eqn2Coeff2a + (  -an3*an1)*eqn2Coeff2b )
 a12(e5,10)=-( an3*an2*eqn2Coeff2a + (  -an3*an2)*eqn2Coeff2b )
 a12(e5,11)=-( an3*an3*eqn2Coeff2a + (1.-an3*an3)*eqn2Coeff2b )

#endMacro


#beginMacro assignUnifiedInterfaceGhost3dOrder4()
  INFO("34c-Unified-Dispersive")

! --- initialize some forcing functions ---
do n=0,nd-1
  fev1(n)=0.
  LfE1(n)=0.
  fEt1(n)=0.
  fEtt1(n)=0.

  fev2(n)=0.
  LfE2(n)=0.
  fEt2(n)=0.
  fEtt2(n)=0.

  fevx1(n)=0.
  fevy1(n)=0.
  fevz1(n)=0.

  fevx2(n)=0.
  fevy2(n)=0.
  fevz2(n)=0.
  
  if (dispersionModel1 .ne. 0) then
    do jv=0,numberOfPolarizationVectors1-1
      fpv1(n,jv)=0.
      LfP1(n,jv)=0.
      fPt1(n,jv)=0.
      fPtt1(n,jv)=0.

      fpvx1(n,jv)=0.
      fpvy1(n,jv)=0.
      fpvz1(n,jv)=0.
    end do
  endif
  if (dispersionModel2 .ne. 0) then
    do jv=0,numberOfPolarizationVectors2-1
      fpv2(n,jv)=0.
      LfP2(n,jv)=0.
      fPt2(n,jv)=0.
      fPtt2(n,jv)=0.

      fpvx2(n,jv)=0.
      fpvy2(n,jv)=0.
      fpvz2(n,jv)=0.
    end do
  endif
end do
! forcing functions for N
if (nonlinearModel1 .ne. 0) then
  do jv = 0,numberOfAtomicLevels1-1
      fnv1(jv) = 0.
      fntv1(jv) = 0.
  enddo
endif
if (nonlinearModel2 .ne. 0) then
  do jv = 0,numberOfAtomicLevels2-1
      fnv2(jv) = 0.
      fntv2(jv) = 0.
  enddo
endif

  err2=0.
  count=0
  beginLoopsMask3d()

    ! here is the normal (assumed to be the same on both sides)
    an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2,an3)
    an2=rsxy1(i1,i2,i3,axis1,1)
    an3=rsxy1(i1,i2,i3,axis1,2)
    aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
    an1=an1/aNorm
    an2=an2/aNorm
    an3=an3/aNorm

    cem1=(1.-1./eps1)/mu1
    cem2=(1.-1./eps2)/mu2

    ! ---- Evaluate the jump conditions using the wrong values at the ghost points  ------
    evaluateUnifiedInterfaceEquations3dOrder4()

   !  write(debugFile,'(" f(6)=",1pe10.2," curlfPttx1,nDotCurlfPtt1=",2(1pe10.2))') f(6),curlfPttx1,nDotCurlfPtt1
   !  write(debugFile,'(" f(7)=",1pe10.2," curlfPttx2,nDotCurlfPtt2=",2(1pe10.2))') f(7),curlfPttx2,nDotCurlfPtt2


    if( debug.gt.4 )then
     write(debugFile,'(/," --> 3d-order4-curv: i1,i2,i3=",3i4," an1,an2,an3=",3e11.3)') i1,i2,i3,an1,an2,an3
     write(debugFile,'(" --> 3d-order4-curv: i1,i2,i3=",3i4," f(start)=",12e10.2)') i1,i2,i3,(f(n),n=0,11)
     ! '
    end if
    if( debug.gt.8 )then
     write(debugFile,'(" --> u1x,u1y,u1z,v1x,v1y,v1z=",6f8.4)') u1x,u1y,u1z,v1x,v1y,v1z
     write(debugFile,'(" --> u2x,u2y,u2z,v2x,v2y,v2z=",6f8.4)') u2x,u2y,u2z,v2x,v2y,v2z

     write(debugFile,'(" --> vv1r,vv1s,vv1t         =",3e9.2)') vv1r,vv1s,vv1t
     do k3=-1,1
     do k2=-1,1
     write(debugFile,'(" --> v1: =",3f8.4)') u1(i1-1,i2+k2,i3+k3,ey),u1(i1,i2+k2,i3+k3,ey),u1(i1+1,i2+k2,i3+k3,ey)
     end do
     end do
     do k3=-1,1
     do k2=-1,1
     write(debugFile,'(" --> v2: =",3f8.4)') u2(j1-1,j2+k2,j3+k3,ey),u2(j1,j2+k2,j3+k3,ey),u2(j1+1,j2+k2,j3+k3,ey)
     end do
     end do
     ! '
    end if

    ! here is the matrix of coefficients for the unknowns u1(-1),v1(-1),w1(-1),  u2(-1),v2(-1),w2(-1)
    ! Solve:
    !     
    !       A [ U ] = A [ U(old) ] - [ f ]


    c1m1x = -is*8.*rsxy1(i1,i2,i3,axis1,0)/(12.*dr1(axis1))    ! coeff of u1(-1) from D.x
    c1m1y = -is*8.*rsxy1(i1,i2,i3,axis1,1)/(12.*dr1(axis1))    ! coeff of u1(-1) from D.y 
    c1m1z = -is*8.*rsxy1(i1,i2,i3,axis1,2)/(12.*dr1(axis1))    ! coeff of u1(-1) from D.z

    c1m2x =  is   *rsxy1(i1,i2,i3,axis1,0)/(12.*dr1(axis1))    ! coeff of u1(-2) from D.x
    c1m2y =  is   *rsxy1(i1,i2,i3,axis1,1)/(12.*dr1(axis1))    ! coeff of u1(-2) from D.y 
    c1m2z =  is   *rsxy1(i1,i2,i3,axis1,2)/(12.*dr1(axis1))    ! coeff of u1(-2) from D.z

    c2m1x = -js*8.*rsxy2(j1,j2,j3,axis2,0)/(12.*dr2(axis2))    ! coeff of u2(-1) from D.x
    c2m1y = -js*8.*rsxy2(j1,j2,j3,axis2,1)/(12.*dr2(axis2))
    c2m1z = -js*8.*rsxy2(j1,j2,j3,axis2,2)/(12.*dr2(axis2))

    c2m2x =  js   *rsxy2(j1,j2,j3,axis2,0)/(12.*dr2(axis2))
    c2m2y =  js   *rsxy2(j1,j2,j3,axis2,1)/(12.*dr2(axis2))
    c2m2z =  js   *rsxy2(j1,j2,j3,axis2,2)/(12.*dr2(axis2))

    rxx1(0,0,0)=aj1rxx
    rxx1(0,0,1)=aj1rxy
    rxx1(0,0,2)=aj1rxz
    rxx1(0,1,1)=aj1ryy
    rxx1(0,1,2)=aj1ryz
    rxx1(0,2,2)=aj1rzz

    rxx1(1,0,0)=aj1sxx
    rxx1(1,0,1)=aj1sxy
    rxx1(1,0,2)=aj1sxz
    rxx1(1,1,1)=aj1syy
    rxx1(1,1,2)=aj1syz
    rxx1(1,2,2)=aj1szz

    rxx1(2,0,0)=aj1txx
    rxx1(2,0,1)=aj1txy
    rxx1(2,0,2)=aj1txz
    rxx1(2,1,1)=aj1tyy
    rxx1(2,1,2)=aj1tyz
    rxx1(2,2,2)=aj1tzz


    rxx2(0,0,0)=aj2rxx
    rxx2(0,0,1)=aj2rxy
    rxx2(0,0,2)=aj2rxz
    rxx2(0,1,1)=aj2ryy
    rxx2(0,1,2)=aj2ryz
    rxx2(0,2,2)=aj2rzz

    rxx2(1,0,0)=aj2sxx
    rxx2(1,0,1)=aj2sxy
    rxx2(1,0,2)=aj2sxz
    rxx2(1,1,1)=aj2syy
    rxx2(1,1,2)=aj2syz
    rxx2(1,2,2)=aj2szz

    rxx2(2,0,0)=aj2txx
    rxx2(2,0,1)=aj2txy
    rxx2(2,0,2)=aj2txz
    rxx2(2,1,1)=aj2tyy
    rxx2(2,1,2)=aj2tyz
    rxx2(2,2,2)=aj2tzz


    ! clap1m1 : coeff of u(-1) from lap1 = u1.xx + u1.yy + u1.zz
    ! clap1m2 : coeff of u(-2) from lap1 = u1.xx + u1.yy + u1.zz

    
    clap1m1=4./3.*(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2+rsxy1(i1,i2,i3,axis1,2)**2)/(dr1(axis1)**2) \
              -is*2./3.*(rxx1(axis1,0,0)+rxx1(axis1,1,1)+rxx1(axis1,2,2))/(2.*dr1(axis1))
    clap1m2=-1./12.*(rsxy1(i1,i2,i3,axis1,0)**2+rsxy1(i1,i2,i3,axis1,1)**2+rsxy1(i1,i2,i3,axis1,2)**2)/(dr1(axis1)**2) \
              +is*1./12.*(rxx1(axis1,0,0)+rxx1(axis1,1,1)+rxx1(axis1,2,2))/(2.*dr1(axis1)) 

    clap2m1=4/3.*(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2+rsxy2(j1,j2,j3,axis2,2)**2)/(dr2(axis2)**2) \
              -js*2./3.*(rxx2(axis2,0,0)+rxx2(axis2,1,1)+rxx2(axis2,2,2))/(2.*dr2(axis2)) 
    clap2m2=-1./12.*(rsxy2(j1,j2,j3,axis2,0)**2+rsxy2(j1,j2,j3,axis2,1)**2+rsxy2(j1,j2,j3,axis2,2)**2)/(dr2(axis2)**2) \
              +js*1./12.*(rxx2(axis2,0,0)+rxx2(axis2,1,1)+rxx2(axis2,2,2))/(2.*dr2(axis2))


    ! rx1(m,n) = D r_m / D x_n
    ! rxx1(m,n1,n2) = D^2 r /( D x_n1 D X_n2 )
    ! rxxx1(m,n1,n2,n3) 

    ! first derivatives
    rx1(0,0)=rsxy1(i1,i2,i3,0,0)
    rx1(0,1)=rsxy1(i1,i2,i3,0,1)
    rx1(0,2)=rsxy1(i1,i2,i3,0,2)
    rx1(1,0)=rsxy1(i1,i2,i3,1,0)
    rx1(1,1)=rsxy1(i1,i2,i3,1,1)
    rx1(1,2)=rsxy1(i1,i2,i3,1,2)
    rx1(2,0)=rsxy1(i1,i2,i3,2,0)
    rx1(2,1)=rsxy1(i1,i2,i3,2,1)
    rx1(2,2)=rsxy1(i1,i2,i3,2,2)

    rx2(0,0)=rsxy2(j1,j2,j3,0,0)
    rx2(0,1)=rsxy2(j1,j2,j3,0,1)
    rx2(0,2)=rsxy2(j1,j2,j3,0,2)
    rx2(1,0)=rsxy2(j1,j2,j3,1,0)
    rx2(1,1)=rsxy2(j1,j2,j3,1,1)
    rx2(1,2)=rsxy2(j1,j2,j3,1,2)
    rx2(2,0)=rsxy2(j1,j2,j3,2,0)
    rx2(2,1)=rsxy2(j1,j2,j3,2,1)
    rx2(2,2)=rsxy2(j1,j2,j3,2,2)


    ! 3rd derivatives: (only some are needed)
    ! note for last 3 entries - entries must increase or stay the same
    rxxx1(axis1,0,0,0) = aj1rxxx
    rxxx1(axis1,1,1,1) = aj1ryyy
    rxxx1(axis1,2,2,2) = aj1rzzz
    rxxx1(axis1,0,0,1) = aj1rxxy
    rxxx1(axis1,0,0,2) = aj1rxxz
    rxxx1(axis1,0,1,1) = aj1rxyy
    rxxx1(axis1,0,2,2) = aj1rxzz
    rxxx1(axis1,1,1,2) = aj1ryyz
    rxxx1(axis1,1,2,2) = aj1ryzz

    rxxx2(axis2,0,0,0) = aj2rxxx
    rxxx2(axis2,1,1,1) = aj2ryyy
    rxxx2(axis2,2,2,2) = aj2rzzz
    rxxx2(axis2,0,0,1) = aj2rxxy
    rxxx2(axis2,0,0,2) = aj2rxxz
    rxxx2(axis2,0,1,1) = aj2rxyy
    rxxx2(axis2,0,2,2) = aj2rxzz
    rxxx2(axis2,1,1,2) = aj2ryyz
    rxxx2(axis2,1,2,2) = aj2ryzz



    ! Some 4th derivatives are needed by LapSq: 
    rxxxx1(axis1,0,0,0,0) = aj1rxxxx
    rxxxx1(axis1,1,1,1,1) = aj1ryyyy
    rxxxx1(axis1,2,2,2,2) = aj1rzzzz

    rxxxx1(axis1,0,0,1,1) = aj1rxxyy
    rxxxx1(axis1,0,0,2,2) = aj1rxxzz
    rxxxx1(axis1,1,1,2,2) = aj1ryyzz

    rxxxx2(axis2,0,0,0,0) = aj2rxxxx
    rxxxx2(axis2,1,1,1,1) = aj2ryyyy
    rxxxx2(axis2,2,2,2,2) = aj2rzzzz

    rxxxx2(axis2,0,0,1,1) = aj2rxxyy
    rxxxx2(axis2,0,0,2,2) = aj2rxxzz
    rxxxx2(axis2,1,1,2,2) = aj2ryyzz


 
    ! coeff of u1(-1) and u1(-2) from Lap^2
    ! dr1a, dr2a : used for avoidInterfaceIterations - tangential spacings are dsBig to eliminate mixed derivatives
    clapSq1m1=lapSqCoeff3DOrder2a(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1,rxxxx1)
    clapSq1m2=lapSqCoeff3DOrder2b(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1,rxxxx1)


    ! coeff of u2(-1) and u2(-2) from Lap^2
    clapSq2m1=lapSqCoeff3DOrder2a(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2,rxxxx2)
    clapSq2m2=lapSqCoeff3DOrder2b(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2,rxxxx2)


    ! coeff of u1(-1) from D.x(Lap), D.y(Lap) and D.z(Lap): (divideb by eps*mu)
    clap1m1x = xLapCoeff3DOrder2a(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1) 
    clap1m1y = yLapCoeff3DOrder2a(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)  
    clap1m1z = zLapCoeff3DOrder2a(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)  

    ! coeff of u1(-2) from D.x(Lap), D.y(Lap) and D.z(Lap): 
    clap1m2x = xLapCoeff3DOrder2b(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)  
    clap1m2y = yLapCoeff3DOrder2b(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)  
    clap1m2z = zLapCoeff3DOrder2b(is,dr1a(axis1),dr1a(axis1p1),dr1a(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1)  


    ! coeff of u2(-1) from D.x(Lap), D.y(Lap) and D.z(Lap): 
    clap2m1x = xLapCoeff3DOrder2a(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)  
    clap2m1y = yLapCoeff3DOrder2a(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2) 
    clap2m1z = zLapCoeff3DOrder2a(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)  

    ! coeff of u2(-2) from D.x(Lap), D.y(Lap) and D.z(Lap): 
    clap2m2x = xLapCoeff3DOrder2b(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)  
    clap2m2y = yLapCoeff3DOrder2b(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)  
    clap2m2z = zLapCoeff3DOrder2b(js,dr2a(axis2),dr2a(axis2p1),dr2a(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2)  


    ! cdivE1 =  u.c1x + v.c1y + w.c1z
    ! nDotCurlE1 = (w1y-v1z)*an1 + (u1z-w1x)*an2 + (v1x-u1y)*an3

    ! 12 Unknowns:
    !   u1(-1), v1(-1), w1(-1), u1(-2), v1(-2), w1(-2), 
    !   u2(-1), v2(-1), w2(-1), u2(-2), v2(-2), w2(-2)  
    ! 12 Equations: 
    !    a12(eqn,unknown) 


    ! fill equations 0,..,5
    fillUnifiedDispersiveEquations3dOrder4a(0,1,2,3,4,5)



    ! coeff of u1(-1) and u1(-2) from Lap^2
    ! clapSq1m1=lapSqCoeff3DOrder2a(is,dr1(axis1),dr1(axis1p1),dr1(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1,rxxxx1)
    ! clapSq1m2=lapSqCoeff3DOrder2b(is,dr1(axis1),dr1(axis1p1),dr1(axis1p2),axis1,axis1p1,axis1p2,rx1,rxx1,rxxx1,rxxxx1)


    ! coeff of u2(-1) and u2(-2) from Lap^2
    ! clapSq2m1=lapSqCoeff3DOrder2a(js,dr2(axis2),dr2(axis2p1),dr2(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2,rxxxx2)
    ! clapSq2m2=lapSqCoeff3DOrder2b(js,dr2(axis2),dr2(axis2p1),dr2(axis2p2),axis2,axis2p1,axis2p2,rx2,rxx2,rxxx2,rxxxx2)

    ! write(debugFile,'(" --> 3d-order4-c: i1,i2,i3=",3i4," c1m1x,c1m1y,c1m1z,clap1m1,clap1m2=",12e10.2)') i1,i2,i3,c1m1x,c1m1y,c1m1z,clap1m1,clap1m2

    ! fill equations 6...11
    fillUnifiedDispersiveEquations3dOrder4b(6,7,8,9,10,11)

    ! Equation 6..11 :  extrapolate 2nd ghost point 
    cex1=1.
    cex2=-5. ! ** fix me ** orderOfExtrapolation for 2nd ghost point 
    ! *e678*
    ! do ii=6,11
    !   do jj=0,11
    !     a12(ii,jj)=0.
    !   end do
    ! end do
    ! a12(6,0)  = cex2   ! u1(-1)
    ! a12(6,3)  = cex1   ! u1(-2)

    ! a12(7,1)  = cex2   ! v1(-1)
    ! a12(7,4)  = cex1   ! v1(-2)

    ! a12(8,2)  = cex2   ! w1(-1)
    ! a12(8,5)  = cex1   ! w1(-2)

    ! a12(9,6)  = cex2   ! u2(-1)
    ! a12(9,9)  = cex1   ! u2(-2)
    ! a12(10,7) = cex2   ! v2(-1)
    ! a12(10,10)= cex1   ! v2(-2)

    ! a12(11,8) = cex2   ! w2(-1)
    ! a12(11,11)= cex1   ! w2(-2)


    ! --- check matrix coefficients by delta function approach ----
    if( checkCoeff.eq.1 .and. it.le.1 )then
      numberOfEquations=12
      checkCoefficients(i1,i2,i3, j1,j2,j3,numberOfEquations,a12,evaluateUnifiedInterfaceEquations3dOrder4 )
    end if

    ! fill in the current values for the unknowns: 
    q(0) = u1(i1-is1,i2-is2,i3-is3,ex)
    q(1) = u1(i1-is1,i2-is2,i3-is3,ey)
    q(2) = u1(i1-is1,i2-is2,i3-is3,ez)
    q(3) = u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)
    q(4) = u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)
    q(5) = u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)

    q(6) = u2(j1-js1,j2-js2,j3-js3,ex)
    q(7) = u2(j1-js1,j2-js2,j3-js3,ey)
    q(8) = u2(j1-js1,j2-js2,j3-js3,ez)
    q(9) = u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)
    q(10)= u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)
    q(11)= u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)

    ! subtract off the contributions from the wrong values at the ghost points:
    numberOfEquations=12
    do n=0,numberOfEquations-1
      f(n) = (a12(n,0)*q(0)+a12(n,1)*q(1)+a12(n,2)*q(2)+a12(n,3)*q(3)+a12(n,4)*q(4)+a12(n,5)*q(5)+\
              a12(n,6)*q(6)+a12(n,7)*q(7)+a12(n,8)*q(8)+a12(n,9)*q(9)+a12(n,10)*q(10)+a12(n,11)*q(11) ) - f(n)
    end do
    if( debug.gt.3 )then
      write(debugFile,'(" --> 3d:order4-c: f(subtract)=",12e10.2)') (f(n),n=0,11)
    end if
    if( .false. )then
      do n=0,numberOfEquations-1
        write(debugFile,'("a(i,j)=",12(1pe10.2))') (a12(n,m),m=0,11)
      end do
    end if

    ! solve A Q = F
    ! factor the matrix
    call dgeco( a12(0,0), numberOfEquations, numberOfEquations, ipvt(0),rcond,work(0))


    if( debug.gt.3 )then
      write(debugFile,'(" --> 3d:order4-c: rcond=",e10.2)') rcond
      write(debugFile,'(" eqn1Coeff1a,eqn1Coeff1b,eqn1Coeff2a,eqn1Coeff2b=",4(1pe10.2))') eqn1Coeff1a,eqn1Coeff1b,eqn1Coeff2a,eqn1Coeff2b
      write(debugFile,'(" eqn2Coeff1a,eqn2Coeff1b,eqn2Coeff2a,eqn2Coeff2b=",4(1pe10.2))') eqn2Coeff1a,eqn2Coeff1b,eqn2Coeff2a,eqn2Coeff2b
    end if

    ! solve
    job=0
    call dgesl( a12(0,0), numberOfEquations, numberOfEquations, ipvt(0), f(0), job)
    if( debug.gt.3 )then
      write(debugFile,'(" --> 3d:order4-c: f(solve)=",12e10.2)') (f(n),n=0,11)
      write(debugFile,'(" --> 3d:order4-c:        q=",12e10.2)') (q(n),n=0,11)
      write(debugFile,'(" --> 3d:order4-c:      f-q=",12e10.2)') (f(n)-q(n),n=0,11)
    end if
    ! fill in the answer:
    if( useJacobiUpdate.eq.0 )then
      u1(i1-is1,i2-is2,i3-is3,ex)      =f(0 )*omega+(1.-omega)*q(0)
      u1(i1-is1,i2-is2,i3-is3,ey)      =f(1 )*omega+(1.-omega)*q(1)
      u1(i1-is1,i2-is2,i3-is3,ez)      =f(2 )*omega+(1.-omega)*q(2)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=f(3 )*omega+(1.-omega)*q(3)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=f(4 )*omega+(1.-omega)*q(4)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=f(5 )*omega+(1.-omega)*q(5)
  
      u2(j1-js1,j2-js2,j3-js3,ex)      =f(6 )*omega+(1.-omega)*q(6)
      u2(j1-js1,j2-js2,j3-js3,ey)      =f(7 )*omega+(1.-omega)*q(7)
      u2(j1-js1,j2-js2,j3-js3,ez)      =f(8 )*omega+(1.-omega)*q(8)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=f(9 )*omega+(1.-omega)*q(9)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=f(10)*omega+(1.-omega)*q(10)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=f(11)*omega+(1.-omega)*q(11)
    else
      ! Jacobi-update -- save results in the work space

      wk1(i1-is1,i2-is2,i3-is3,ex)      =f(0 )*omega+(1.-omega)*q(0)
      wk1(i1-is1,i2-is2,i3-is3,ey)      =f(1 )*omega+(1.-omega)*q(1)
      wk1(i1-is1,i2-is2,i3-is3,ez)      =f(2 )*omega+(1.-omega)*q(2)
      wk1(i1-2*is1,i2-2*is2,i3-2*is3,ex)=f(3 )*omega+(1.-omega)*q(3)
      wk1(i1-2*is1,i2-2*is2,i3-2*is3,ey)=f(4 )*omega+(1.-omega)*q(4)
      wk1(i1-2*is1,i2-2*is2,i3-2*is3,ez)=f(5 )*omega+(1.-omega)*q(5)
  
      wk2(j1-js1,j2-js2,j3-js3,ex)      =f(6 )*omega+(1.-omega)*q(6)
      wk2(j1-js1,j2-js2,j3-js3,ey)      =f(7 )*omega+(1.-omega)*q(7)
      wk2(j1-js1,j2-js2,j3-js3,ez)      =f(8 )*omega+(1.-omega)*q(8)
      wk2(j1-2*js1,j2-2*js2,j3-2*js3,ex)=f(9 )*omega+(1.-omega)*q(9)
      wk2(j1-2*js1,j2-2*js2,j3-2*js3,ey)=f(10)*omega+(1.-omega)*q(10)
      wk2(j1-2*js1,j2-2*js2,j3-2*js3,ez)=f(11)*omega+(1.-omega)*q(11)

    end if

    ! compute the maximum change in the solution for this iteration
    do n=0,11
      err=max(err,abs(q(n)-f(n)))
      err2 = err2 + (q(n)-f(n))**2
      count = count + 1
    end do

    if( .true. .or. debug.gt.3 )then ! re-evaluate

     do n=0,11
       errv(n)=abs(q(n)-f(n))
     end do
     ! evalDerivs3dOrder4()
     ! eval3dJumpOrder4New()

     evaluateUnifiedInterfaceEquations3dOrder4()

     res=0.
     do n=0,11
       res=max(res,abs(f(n)))
     end do
     if( .false. .and. res.gt.1.e-9 )then
       write(debugFile,'(" --> ERR: 3d-GDM-order4-c: it=",i3," i1,i2,i3=",3i4," f(re-eval)=",12e10.2)') it,i1,i2,i3,f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),f(11)
       write(debugFile,'("     it=",i3," err in ghost=",12e10.2)') it,(errv(n),n=0,11)
     end if
    end if

  endLoopsMask3d()

  err2 = sqrt( err2/count )

  if( checkCoeff.eq.1 .and. it.le.1 )then
    write(*,'("+++++ IGDM34c: check coeff in interface: max(diff) = ",1pe8.2)') coeffDiff
  end if

  if( useJacobiUpdate.ne.0 )then
    ! Jacobi update -- copy work-space to solution arrays
    beginLoopsMask3d()
      u1(i1-is1,i2-is2,i3-is3,ex)      = wk1(i1-is1,i2-is2,i3-is3,ex)      
      u1(i1-is1,i2-is2,i3-is3,ey)      = wk1(i1-is1,i2-is2,i3-is3,ey)      
      u1(i1-is1,i2-is2,i3-is3,ez)      = wk1(i1-is1,i2-is2,i3-is3,ez)      
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ex)= wk1(i1-2*is1,i2-2*is2,i3-2*is3,ex)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ey)= wk1(i1-2*is1,i2-2*is2,i3-2*is3,ey)
      u1(i1-2*is1,i2-2*is2,i3-2*is3,ez)= wk1(i1-2*is1,i2-2*is2,i3-2*is3,ez)
                                                                           
      u2(j1-js1,j2-js2,j3-js3,ex)      = wk2(j1-js1,j2-js2,j3-js3,ex)      
      u2(j1-js1,j2-js2,j3-js3,ey)      = wk2(j1-js1,j2-js2,j3-js3,ey)      
      u2(j1-js1,j2-js2,j3-js3,ez)      = wk2(j1-js1,j2-js2,j3-js3,ez)      
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ex)= wk2(j1-2*js1,j2-2*js2,j3-2*js3,ex)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ey)= wk2(j1-2*js1,j2-2*js2,j3-2*js3,ey)
      u2(j1-2*js1,j2-2*js2,j3-2*js3,ez)= wk2(j1-2*js1,j2-2*js2,j3-2*js3,ez)
    endLoopsMask3d()
  end if

#endMacro