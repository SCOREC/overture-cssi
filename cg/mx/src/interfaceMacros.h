! -*- mode: F90; -*-
! ==================================================================================================
! Macro: compute the normal vector (an1,an2) in 2D
! GRIDTYPE (input) : curvilinear or rectangular
! ==================================================================================================
#beginMacro getNormal2d(an1,an2,GRIDTYPE)
  #If #GRIDTYPE eq "curvilinear"
   an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
   an2=rsxy1(i1,i2,i3,axis1,1)
   aNorm=max(epsx,sqrt(an1**2+an2**2))
   an1=an1/aNorm
   an2=an2/aNorm
  #Elif #GRIDTYPE eq "rectangular"
   an1=an1Cartesian
   an2=an2Cartesian
  #Else
     stop 1111
  #End
#endMacro

! ==================================================================================================
! Macro: compute the normal vector (an1,an2,an3) in 3D
! GRIDTYPE (input) : curvilinear or rectangular
! ==================================================================================================
#beginMacro getNormal3d(an1,an2,an3,GRIDTYPE)
  #If #GRIDTYPE eq "curvilinear"
   an1=rsxy1(i1,i2,i3,axis1,0)  
   an2=rsxy1(i1,i2,i3,axis1,1)
   an3=rsxy1(i1,i2,i3,axis1,2)
   aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
   an1=an1/aNorm
   an2=an2/aNorm
   an3=an3/aNorm
  #Elif #GRIDTYPE eq "rectangular"
   an1=an1Cartesian
   an2=an2Cartesian
   an3=an3Cartesian
  #Else
     stop 3333
  #End
#endMacro

! ==================================================================================================
! This macro will assign the jump conditions on the boundary
! DIM (input): number of dimensions (2 or 3)
! GRIDTYPE (input) : curvilinear or rectangular
! ==================================================================================================
#beginMacro boundaryJumpConditions(DIM,GRIDTYPE)
if( dispersive.gt.0 )then

 ! ****** DISPERSIVE CASE *******
 #If #DIM eq "2"
  if( .true. .and. useImpedanceInterfaceProjection.eq.1 )then
    ! --------------------------------------------------------------
    ! ------ Use impedance weighting to project the interface ------
    ! --------------------------------------------------------------
    !  (see maxwell.pdf)
    !

    if( t.le.3*dt .and. debug.gt.0  )then
      write(*,'("cgmx:interface3d: PROJECT INTERFACE DISPERSIVE in 2D")')
    end if


    ! --- init polarization vectors in case one side is non-dispersive --- 
    do n=0,nd-1   
      p1v(n)=0.
      p2v(n)=0.
    end do
    g1=0.
    g2=0.

    ! --- LOOP over the interface ---
    beginGhostLoopsMask2d()

      getNormal2d(an1,an2,GRIDTYPE)

      ! left state:
      ex1=u1(i1,i2,i3,ex)
      ey1=u1(i1,i2,i3,ey)
      hz1=u1(i1,i2,i3,hz)

      if( dispersionModel1.ne.noDispersion )then
        ! eval sum of Pv vectors 
        do n=0,nd-1
          p1v(n)=0.
          do jv=0,numberOfPolarizationVectors1-1
            pc = jv*nd
            p1v(n)=p1v(n) + p1(i1,i2,i3,pc+n)
          end do
        end do
      end if

      ! displacement vector 
      D1v(0) = eps1*(ex1 + alphaP1*p1v(0))
      D1v(1) = eps1*(ey1 + alphaP1*p1v(1))

      ! right state:
      ex2=u2(j1,j2,j3,ex)
      ey2=u2(j1,j2,j3,ey)
      hz2=u2(j1,j2,j3,hz)

      if( dispersionModel2.ne.noDispersion )then
        ! eval sum of Pv vectors 
        do n=0,nd-1
          p2v(n)=0.
          do jv=0,numberOfPolarizationVectors2-1
            pc = jv*nd
            p2v(n)=p2v(n) + p2(j1,j2,j3,pc+n)
          end do
        end do
      end if

      ! displacement vector 
      D2v(0) = eps2*(ex2 + alphaP2*p2v(0))
      D2v(1) = eps2*(ey2 + alphaP2*p2v(1))

      ! normal components 
      ! nDotE1 = an1*ex1+an2*ey1
      ! nDotE2 = an1*ex2+an2*ey2
   
      nDotD1 = an1*D1v(0)+an2*D1v(1)
      nDotD2 = an1*D2v(0)+an2*D2v(1)
   
      nDotP1 = an1*p1v(0)+an2*p1v(1)
      nDotP2 = an1*p2v(0)+an2*p2v(1)
   

      ! Interface displacement-vector is an impedance weighted average
      nDotDI = ( eta1*nDotD1 + eta2*nDotD2 )/(eta1+eta2)

      ! Normal components of n.E on the left and right:
      nDotE1I = nDotDI/eps1 -alphaP1*nDotP1
      nDotE2I = nDotDI/eps2 -alphaP2*nDotP2

      ! nDotE1I = nDotE1 + (1./(1.*eps1))*(nDotDI-nDotD1)
      ! nDotE2I = nDotE2 + (1./(1.*eps2))*(nDotDI-nDotD2)


      if( twilightZone.eq.1)then
       ! --- adjust for TZ forcing ---
       ! *NOTE* Here we assume the exact solution for E is the same on both sides
       ! *NOTE* exact soltuion for P may be different 

       ! See notes...
       !    gm = Ee - D_e^I /epsm + alphaP*Pe
       !    D_e^I = (eta1*D1e + eta2*D2e )/( eta1+eta2 )


       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ue )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, ve )
       nDotEe = an1*ue+an2*ve
 
       ! Compute sum of exact polarization vectors: 
       do n=0,nd-1
        p1ev(n)=0. 
        do jv=0,numberOfPolarizationVectors1-1
          ! NOTE: the TZ component number is offset by pxc
          pc = pxc + jv*nd
          call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,pc+n, pe(n)  )
          p1ev(n) = p1ev(n) + pe(n)
        end do
        p2ev(n)=0. 
        do jv=0,numberOfPolarizationVectors2-1
          ! NOTE: the TZ component number is offset by pxc
          pc = pxc + jv*nd
          call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,pc+n, pe(n)  )
          p2ev(n) = p2ev(n) + pe(n)
        end do
       end do

       ! normal components 
       nDotP1e = an1*p1ev(0) + an2*p1ev(1)
       nDotP2e = an1*p2ev(0) + an2*p2ev(1)
       nDotDe = ( eta1*eps1*(nDotEe+alphaP1*nDotP1e) + eta2*eps2*(nDotEe+alphaP2*nDotP2e) )/(eta1+eta2)

       ! epsNDotEI + g1 = eps1*nDotE 
       ! epsNDotEI + g2 = eps2*nDotE 
       g1= nDotEe - nDotDe/eps1 + alphaP1*nDotP1e 
       g2= nDotEe - nDotDe/eps2 + alphaP2*nDotP2e 

       ! Adjust n.E 
       nDotE1I = nDotE1I+g1  ! nDotE for interface on left 
       nDotE2I = nDotE2I+g2  ! nDotE for interface on right 

      end if


      ! inverse impedance average of tangential components  (set values for full vector and then correct below)

      exI = ( eta1i*ex1 + eta2i*ex2 )/( eta1i + eta2i )
      eyI = ( eta1i*ey1 + eta2i*ey2 )/( eta1i + eta2i )

      ! hz : impedance weighted average 
      hzI = ( eta1*hz1 + eta2*hz2 )/(eta1+eta2) 

      nDotEI= an1*exI+an2*eyI  ! we need to subtract off normal component of (exI,eyI) 

      u1(i1,i2,i3,ex) = exI + (nDotE1I - nDotEI)*an1
      u1(i1,i2,i3,ey) = eyI + (nDotE1I - nDotEI)*an2
      u1(i1,i2,i3,hz) = hzI 

      u2(j1,j2,j3,ex) = exI + (nDotE2I - nDotEI)*an1
      u2(j1,j2,j3,ey) = eyI + (nDotE2I - nDotEI)*an2
      u2(j1,j2,j3,hz) = hzI 



      if( .false. .and. twilightZone.eq.1 )then
          ! check errors on the boundary
          do n=0,nd-1
            call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t,ex+n, es(n)   ) 
          end do
          f(0) =  u1(i1,i2,i3,ex) -es(0)
          f(1) =  u1(i1,i2,i3,ey) -es(1)
          do n=0,nd-1
            call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),0.,t,ex+n, est(n)   ) 
          end do
          f(2) =  u2(j1,j2,j3,ex) -est(0)
          f(3) =  u2(j1,j2,j3,ey) -est(1)
          write(*,'(" PROJECT-INTERFACE ERRORS =",4e10.2)') f(0),f(1),f(2),f(3) 
  
      end if

    endLoopsMask2d()

  else 

    write(*,'("cgmx:interface3d: FINISH ME -- TWO-D dispersive bndry jump conditions")')

  end if
 #Else


  ! ******************** 3D PROJECT INTERFACE DISPERSIVE  ******************************

  if( .true. .and. useImpedanceInterfaceProjection.eq.1 )then
    ! --------------------------------------------------------------
    ! ------ Use impedance weighting to project the interface ------
    ! --------------------------------------------------------------
    !  (see maxwell.pdf)
    !

    if( t.le.3*dt .and. debug.gt.0 )then
      write(*,'("cgmx:interface3d: PROJECT INTERFACE DISPERSIVE in 3D")')
    end if

    ! --- init polarization vectors in case one side is non-dispersive --- 
    do n=0,nd-1   
      p1v(n)=0.
      p2v(n)=0.
    end do
    g1=0.
    g2=0.

    ! --- LOOP over the interface ---
    beginGhostLoopsMask3d()

      getNormal3d(an1,an2,an3,GRIDTYPE)

      ! left state:
      ex1=u1(i1,i2,i3,ex)
      ey1=u1(i1,i2,i3,ey)
      ez1=u1(i1,i2,i3,ez)

      if( dispersionModel1.ne.noDispersion )then
        ! eval sum of Pv vectors 
        do n=0,nd-1
          p1v(n)=0.
          do jv=0,numberOfPolarizationVectors1-1
            pc = jv*nd
            p1v(n)=p1v(n) + p1(i1,i2,i3,pc+n)
          end do
        end do
      end if

      ! left displacement vector 
      D1v(0) = eps1*(ex1 + alphaP1*p1v(0))
      D1v(1) = eps1*(ey1 + alphaP1*p1v(1))
      D1v(2) = eps1*(ez1 + alphaP1*p1v(2))

      ! right state:
      ex2=u2(j1,j2,j3,ex)
      ey2=u2(j1,j2,j3,ey)
      ez2=u2(j1,j2,j3,ez)

      if( dispersionModel2.ne.noDispersion )then
        ! eval sum of Pv vectors 
        do n=0,nd-1
          p2v(n)=0.
          do jv=0,numberOfPolarizationVectors2-1
            pc = jv*nd
            p2v(n)=p2v(n) + p2(j1,j2,j3,pc+n)
          end do
        end do
      end if

      ! right displacement vector 
      D2v(0) = eps2*(ex2 + alphaP2*p2v(0))
      D2v(1) = eps2*(ey2 + alphaP2*p2v(1))
      D2v(2) = eps2*(ez2 + alphaP2*p2v(2))

      ! normal components of D1 and D2
      nDotD1 = an1*D1v(0)+an2*D1v(1)+an3*D1v(2)
      nDotD2 = an1*D2v(0)+an2*D2v(1)+an3*D2v(2)
   
      ! normal components of P1 and P2
      nDotP1 = an1*p1v(0)+an2*p1v(1)+an3*p1v(2)
      nDotP2 = an1*p2v(0)+an2*p2v(1)+an3*p2v(2)
   

      ! Interface displacement-vector is an impedance weighted average
      nDotDI = ( eta1*nDotD1 + eta2*nDotD2 )/(eta1+eta2)

      ! Normal components of n.E on the left and right:
      nDotE1I = nDotDI/eps1 -alphaP1*nDotP1
      nDotE2I = nDotDI/eps2 -alphaP2*nDotP2


      if( twilightZone.eq.1)then
       ! --- adjust for TZ forcing ---
       ! *NOTE* Here we assume the exact solution for E is the same on both sides
       ! *NOTE* exact soltuion for P may be different 

       ! See notes...
       !    gm = Ee - D_e^I /epsm + alphaP*Pe
       !    D_e^I = (eta1*D1e + eta2*D2e )/( eta1+eta2 )


       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, ue )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, ve )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, we )
       nDotEe = an1*ue+an2*ve+an3*we
 
       ! Compute sum of exact polarization vectors: 
       do n=0,nd-1
        p1ev(n)=0. 
        do jv=0,numberOfPolarizationVectors1-1
          ! NOTE: the TZ component number is offset by pxc
          pc = pxc + jv*nd
          call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t,pc+n, pe(n)  )
          p1ev(n) = p1ev(n) + pe(n)
        end do
        p2ev(n)=0. 
        do jv=0,numberOfPolarizationVectors2-1
          ! NOTE: the TZ component number is offset by pxc
          pc = pxc + jv*nd
          call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),xy1(i1,i2,i3,2),t,pc+n, pe(n)  )
          p2ev(n) = p2ev(n) + pe(n)
        end do
       end do

       ! normal components 
       nDotP1e = an1*p1ev(0) + an2*p1ev(1) + an3*p1ev(2)
       nDotP2e = an1*p2ev(0) + an2*p2ev(1) + an3*p2ev(2)
       nDotDe = ( eta1*eps1*(nDotEe+alphaP1*nDotP1e) + eta2*eps2*(nDotEe+alphaP2*nDotP2e) )/(eta1+eta2)

       ! epsNDotEI + g1 = eps1*nDotE 
       ! epsNDotEI + g2 = eps2*nDotE 
       g1= nDotEe - nDotDe/eps1 + alphaP1*nDotP1e 
       g2= nDotEe - nDotDe/eps2 + alphaP2*nDotP2e 

       ! Adjust n.E 
       nDotE1I = nDotE1I+g1  ! nDotE for interface on left 
       nDotE2I = nDotE2I+g2  ! nDotE for interface on right 

      end if

      ! inverse impedance average of tangential components  (set values for full vector and then correct below)

      exI = ( eta1i*ex1 + eta2i*ex2 )/( eta1i + eta2i )
      eyI = ( eta1i*ey1 + eta2i*ey2 )/( eta1i + eta2i )
      ezI = ( eta1i*ez1 + eta2i*ez2 )/( eta1i + eta2i )

      nDotEI= an1*exI+an2*eyI+an3*ezI  ! we need to subtract off normal component of (exI,eyI) 

      u1(i1,i2,i3,ex) = exI + (nDotE1I - nDotEI)*an1
      u1(i1,i2,i3,ey) = eyI + (nDotE1I - nDotEI)*an2
      u1(i1,i2,i3,ez) = ezI + (nDotE1I - nDotEI)*an3

      u2(j1,j2,j3,ex) = exI + (nDotE2I - nDotEI)*an1
      u2(j1,j2,j3,ey) = eyI + (nDotE2I - nDotEI)*an2
      u2(j1,j2,j3,ez) = ezI + (nDotE2I - nDotEI)*an3


      if( .false. .and. twilightZone.eq.1 )then
          ! check errors on the boundary
          do n=0,nd-1
            call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t,ex+n, es(n)   ) 
          end do
          f(0) =  u1(i1,i2,i3,ex) -es(0)
          f(1) =  u1(i1,i2,i3,ey) -es(1)
          f(2) =  u1(i1,i2,i3,ez) -es(2)
          do n=0,nd-1
            call ogderiv(ep, 0,0,0,0, xy2(j1,j2,j3,0),xy2(j1,j2,j3,1),xy1(i1,i2,i3,2),t,ex+n, est(n)   ) 
          end do
          f(3) =  u2(j1,j2,j3,ex) -est(0)
          f(4) =  u2(j1,j2,j3,ey) -est(1)
          f(5) =  u2(j1,j2,j3,ez) -est(2)
          write(*,'(" PROJECT-INTERFACE ERRORS =",6e10.2)') f(0),f(1),f(2),f(3),f(4),f(5) 
  
      end if

    endLoopsMask3d()

  else 

    stop 1234

  end if

 #End

else

 ! ****** NON-DISPERSIVE CASE *******

 #If #DIM eq "2"
  if( useImpedanceInterfaceProjection.eq.1 )then
    ! --------------------------------------------------------------
    ! ------ Use impedance weighting to project the interface ------
    ! --------------------------------------------------------------
    !  (see maxwell.pdf)
    !
    beginGhostLoopsMask2d()
      ! eps2 n.u2 = eps1 n.u1
      !     tau.u2 = tau.u1

      getNormal2d(an1,an2,GRIDTYPE)

      ! left state:
      ex1=u1(i1,i2,i3,ex)
      ey1=u1(i1,i2,i3,ey)
      hz1=u1(i1,i2,i3,hz)

      ! right state:
      ex2=u2(j1,j2,j3,ex)
      ey2=u2(j1,j2,j3,ey)
      hz2=u2(j1,j2,j3,hz)

      ! normal components 
      nDotE1 = an1*ex1+an2*ey1
      nDotE2 = an1*ex2+an2*ey2
   
      ! The interface value of (eps n.E) is an impedance average of (eps* n.E)
      !
      !      (eps n.E)_I = [ eta1*( eps1*n.E1 ) +  eta2*( eps2*n.E2 ) ]/[ eta1 + eta2 ]
      !
      ! We then set 
      !   eps1*nDotE1_I = (eps n.E)_I 
      !   eps2*nDotE2_I = (eps n.E)_I 
      epsNDotEI = ( eta1*(eps1*nDotE1) + eta2*(eps2*nDotE2) )/(eta1+eta2) ! (eps * n. E)_I 

      if( twilightZone.eq.1)then
       ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ue )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, ve )
       nDotEe = an1*ue+an2*ve

       ! epsNDotEI + g1 = eps1*nDotE 
       ! epsNDotEI + g2 = eps2*nDotE 
       g1= (eps1-eps2) * nDotEe *eta2/(eta1+eta2)
       g2= (eps2-eps1) * nDotEe *eta1/(eta1+eta2)

      else
        g1=0.
        g2=0. 
      end if

      nDotE1I = (epsNDotEI+g1)/eps1  ! nDotE for interface on left 
      nDotE2I = (epsNDotEI+g2)/eps2  ! nDotE for interface on right 

      ! inverse impedance average of tangential components  (do full vector and correct below)

      exI = ( eta1i*ex1 + eta2i*ex2 )/( eta1i + eta2i )
      eyI = ( eta1i*ey1 + eta2i*ey2 )/( eta1i + eta2i )

      ! hz : impedance weighted average 
      hzI = ( eta1*hz1 + eta2*hz2 )/(eta1+eta2) 

      nDotEI= an1*exI+an2*eyI  ! we need to subtract off normal component of (exI,eyI) 


      u1(i1,i2,i3,ex) = exI + (nDotE1I - nDotEI)*an1
      u1(i1,i2,i3,ey) = eyI + (nDotE1I - nDotEI)*an2
      u1(i1,i2,i3,hz) = hzI 

      u2(j1,j2,j3,ex) = exI + (nDotE2I - nDotEI)*an1
      u2(j1,j2,j3,ey) = eyI + (nDotE2I - nDotEI)*an2
      u2(j1,j2,j3,hz) = hzI 


    endLoopsMask2d()

  
  else if( eps1.lt.eps2 )then
    epsRatio=eps1/eps2
    beginGhostLoopsMask2d()
      ! eps2 n.u2 = eps1 n.u1
      !     tau.u2 = tau.u1

      #If #GRIDTYPE eq "curvilinear"
       an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
       an2=rsxy1(i1,i2,i3,axis1,1)
       aNorm=max(epsx,sqrt(an1**2+an2**2))
       an1=an1/aNorm
       an2=an2/aNorm
      #Elif #GRIDTYPE eq "rectangular"
       an1=an1Cartesian
       an2=an2Cartesian
      #Else
         stop 1111
      #End
      ua=u1(i1,i2,i3,ex)
      ub=u1(i1,i2,i3,ey)
      nDotU = an1*ua+an2*ub
      if( twilightZone.eq.1 )then
       ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ue )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, ve )
       nDotU = nDotU - (an1*ue+an2*ve)
      end if

      ! u2 equals u1 but with normal component = eps1/eps2*(n.u1)
      u2(j1,j2,j3,ex) = ua + (nDotU*epsRatio - nDotU)*an1
      u2(j1,j2,j3,ey) = ub + (nDotU*epsRatio - nDotU)*an2
      u2(j1,j2,j3,hz) = u1(i1,i2,i3,hz)


    endLoopsMask2d()
  else
    epsRatio=eps2/eps1
    beginGhostLoopsMask2d()
      ! eps2 n.u2 = eps1 n.u1
      !     tau.u2 = tau.u1

      #If #GRIDTYPE eq "curvilinear"
       an1=rsxy1(i1,i2,i3,axis1,0)
       an2=rsxy1(i1,i2,i3,axis1,1)
       aNorm=max(epsx,sqrt(an1**2+an2**2))
       an1=an1/aNorm
       an2=an2/aNorm
      #Elif #GRIDTYPE eq "rectangular"
       an1=an1Cartesian
       an2=an2Cartesian
      #Else
        stop 1112
      #End
      ua=u2(j1,j2,j3,ex)
      ub=u2(j1,j2,j3,ey)

      nDotU = an1*ua+an2*ub
      if( twilightZone.eq.1 )then
       ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ue )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, ve )
! write(*,'(" jump: x,y=",2e10.2," ua,ue=",2e10.2," ub,ve=",2e10.2)') xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),ua,ue,ub,ve
       nDotU = nDotU - (an1*ue+an2*ve)
      end if

      u1(i1,i2,i3,ex) = ua + (nDotU*epsRatio - nDotU)*an1
      u1(i1,i2,i3,ey) = ub + (nDotU*epsRatio - nDotU)*an2
      u1(i1,i2,i3,hz) = u2(j1,j2,j3,hz)
    endLoopsMask2d()
  end if
 #Else
  ! ******************** 3D PROJECT INTERFACE ******************************
  if( useImpedanceInterfaceProjection.eq.1 )then
    ! --------------------------------------------------------------
    ! ------ Use impedance weighting to project the interface ------
    ! --------------------------------------------------------------
    !  (see maxwell.pdf)
    !
    beginGhostLoopsMask3d()
      #If #GRIDTYPE eq "curvilinear"
       an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
       an2=rsxy1(i1,i2,i3,axis1,1)
       an3=rsxy1(i1,i2,i3,axis1,2)
       aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
       an1=an1/aNorm
       an2=an2/aNorm
       an3=an3/aNorm
      #Elif #GRIDTYPE eq "rectangular"
       an1=an1Cartesian
       an2=an2Cartesian
       an3=an3Cartesian
      #Else
         stop 1111
      #End
      ! left state:
      ex1=u1(i1,i2,i3,ex)
      ey1=u1(i1,i2,i3,ey)
      ez1=u1(i1,i2,i3,ez)

      ! right state:
      ex2=u2(j1,j2,j3,ex)
      ey2=u2(j1,j2,j3,ey)
      ez2=u2(j1,j2,j3,ez)

      ! normal components 
      nDotE1 = an1*ex1+an2*ey1+an3*ez1
      nDotE2 = an1*ex2+an2*ey2+an3*ez2

      ! The interface value of (eps n.E) is an impedance average of (eps* n.E)
      !
      !      (eps n.E)_I = [ eta1*( eps1*n.E1 ) +  eta2*( eps2*n.E2 ) ]/[ eta1 + eta2 ]
      !
      ! We then set 
      !   eps1*nDotE1_I = (eps n.E)_I 
      !   eps2*nDotE2_I = (eps n.E)_I 
      epsNDotEI = ( eta1*(eps1*nDotE1) + eta2*(eps2*nDotE2) )/(eta1+eta2) ! (eps * n. E)_I 

      if( twilightZone.eq.1)then
       ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, ue )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, ve )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, we )
       nDotEe = an1*ue+an2*ve+an3*we

       ! epsNDotEI + g1 = eps1*nDotE 
       ! epsNDotEI + g2 = eps2*nDotE 
       g1= (eps1-eps2) * nDotEe *eta2/(eta1+eta2)
       g2= (eps2-eps1) * nDotEe *eta1/(eta1+eta2)

      else
        g1=0.
        g2=0. 
      end if

      nDotE1I = (epsNDotEI+g1)/eps1  ! nDotE for interface on left 
      nDotE2I = (epsNDotEI+g2)/eps2  ! nDotE for interface on right 

      ! inverse impedance average of tangential components  (do full vector and correct below)
      exI = ( eta1i*ex1 + eta2i*ex2 )/( eta1i + eta2i )
      eyI = ( eta1i*ey1 + eta2i*ey2 )/( eta1i + eta2i )
      ezI = ( eta1i*ez1 + eta2i*ez2 )/( eta1i + eta2i )

      nDotEI= an1*exI+an2*eyI+an3*ezI  ! we need to subtract off normal component of (exI,eyI,ezI) 

      u1(i1,i2,i3,ex) = exI + (nDotE1I - nDotEI)*an1
      u1(i1,i2,i3,ey) = eyI + (nDotE1I - nDotEI)*an2
      u1(i1,i2,i3,ez) = ezI + (nDotE1I - nDotEI)*an3

      u2(j1,j2,j3,ex) = exI + (nDotE2I - nDotEI)*an1
      u2(j1,j2,j3,ey) = eyI + (nDotE2I - nDotEI)*an2
      u2(j1,j2,j3,ez) = ezI + (nDotE2I - nDotEI)*an3

    endLoopsMask3d()

  else if( eps1.lt.eps2 )then
    epsRatio=eps1/eps2
    beginGhostLoopsMask3d()
      ! eps2 n.u2 = eps1 n.u1
      !     tau.u2 = tau.u1

      #If #GRIDTYPE eq "curvilinear"
       an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
       an2=rsxy1(i1,i2,i3,axis1,1)
       an3=rsxy1(i1,i2,i3,axis1,2)
       aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
       an1=an1/aNorm
       an2=an2/aNorm
       an3=an3/aNorm
      #Elif #GRIDTYPE eq "rectangular"
       an1=an1Cartesian
       an2=an2Cartesian
       an3=an3Cartesian
      #Else
         stop 1111
      #End
      ua=u1(i1,i2,i3,ex)
      ub=u1(i1,i2,i3,ey)
      uc=u1(i1,i2,i3,ez)
      nDotU = an1*ua+an2*ub+an3*uc
      if( twilightZone.eq.1 )then
       ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, ue )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, ve )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, we )
       nDotU = nDotU - (an1*ue+an2*ve+an3*we)
      end if
      ! u2 equals u1 but with normal component = eps1/eps2*(n.u1)
      u2(j1,j2,j3,ex) = ua + (nDotU*epsRatio - nDotU)*an1
      u2(j1,j2,j3,ey) = ub + (nDotU*epsRatio - nDotU)*an2
      u2(j1,j2,j3,ez) = uc + (nDotU*epsRatio - nDotU)*an3

!   write(*,'(" jump(1): (i1,i2,i3)=",3i3," j1,j2,j3=",3i3)') i1,i2,i3,j1,j2,j3
!   write(*,'(" jump(1): x,y,z=",3e10.2," ua,ue=",2e10.2," ub,ve=",2e10.2," uc,we=",2e10.2)') xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),ua,ue,ub,ve,uc,we

    endLoopsMask3d()
  else
    epsRatio=eps2/eps1
    beginGhostLoopsMask3d()
      ! eps2 n.u2 = eps1 n.u1
      !     tau.u2 = tau.u1

      #If #GRIDTYPE eq "curvilinear"
       an1=rsxy1(i1,i2,i3,axis1,0)   ! normal (an1,an2)
       an2=rsxy1(i1,i2,i3,axis1,1)
       an3=rsxy1(i1,i2,i3,axis1,2)
       aNorm=max(epsx,sqrt(an1**2+an2**2+an3**2))
       an1=an1/aNorm
       an2=an2/aNorm
       an3=an3/aNorm
      #Elif #GRIDTYPE eq "rectangular"
       an1=an1Cartesian
       an2=an2Cartesian
       an3=an3Cartesian
      #Else
         stop 1111
      #End
      ua=u2(j1,j2,j3,ex)
      ub=u2(j1,j2,j3,ey)
      uc=u2(j1,j2,j3,ez)

      nDotU = an1*ua+an2*ub+an3*uc
      if( twilightZone.eq.1 )then
       ! adjust for TZ forcing (here we assume the exact solution is the same on both sides)
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ex, ue )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ey, ve )
       call ogderiv(ep, 0,0,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),t, ez, we )
       nDotU = nDotU - (an1*ue+an2*ve+an3*we)

!   write(*,'(" jump(2): (i1,i2,i3)=",3i3," j1,j2,j3=",3i3)') i1,i2,i3,j1,j2,j3
!   write(*,'(" jump(2): x,y,z=",3e10.2," ua,ue=",2e10.2," ub,ve=",2e10.2," uc,we=",2e10.2)') xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),xy1(i1,i2,i3,2),ua,ue,ub,ve,uc,we

      end if

      u1(i1,i2,i3,ex) = ua + (nDotU*epsRatio - nDotU)*an1
      u1(i1,i2,i3,ey) = ub + (nDotU*epsRatio - nDotU)*an2
      u1(i1,i2,i3,ez) = uc + (nDotU*epsRatio - nDotU)*an3
    endLoopsMask3d()
  end if
 #End
end if
#endMacro

! ===========================================================================================
!  Smooth P on the interface 
! ===========================================================================================
#beginMacro smoothInterfaceP()
 if( .true. .and. orderOfAccuracy.eq.4  )then

   ! ----- Add some smoothing on the boundary for P ------
   ! Fourth-order filter: 
   !    p <- p + cd*[ -1 p(i-2) + 4*p(i-1) -6 p(i) + 4 p(i+1) - p(i+2) ]
   !           + cd*[ -1 p(j-2) + 4*p(j-1) -6 p(j) + 4 p(j+1) - p(j+2) ]
   ! cd = 1.(sum-of-coefficients)
   ! In 2d: cd=16*2
   ! In 3d: cd=16*3
   write(*,'(" interface-project: add smoothing to P")') 
   if( nd.eq.2 )then
    ! --- LOOP over the interface ---
    ! Question: should we use a Jacobi iteration?
    beginGhostLoopsMask2d()
     if( dispersionModel1.ne.noDispersion )then
       do n=0,nd-1
        do jv=0,numberOfPolarizationVectors1-1
          pc = jv*nd+n
          ! smooth in the normal direction: 
          p1(i1,i2,i3,pc) = ( -( p1(i1-2*is1,i2-2*is2,i3,pc)+p1(i1+2*is1,i2+2*is2,i3,pc) ) \
                           +4.*( p1(i1-  is1,i2-  is2,i3,pc)+p1(i1+  is1,i2+  is2,i3,pc) ) \
                           +10.* p1(i1,i2,i3,pc) )/16.
          ! p1(i1,i2,i3,pc) = ( -( p1(i1-2,i2,i3,pc)+p1(i1+2,i2,i3,pc)+p1(i1,i2-2,i3,pc)+p1(i1,i2+2,i3,pc) ) \
          !                 +4.*( p1(i1-1,i2,i3,pc)+p1(i1+1,i2,i3,pc)+p1(i1,i2-1,i3,pc)+p1(i1,i2+1,i3,pc) ) \
          !                 +20.* p1(i1,i2,i3,pc) )/32.
        end do
       end do
     end if
     if( dispersionModel2.ne.noDispersion )then
       do n=0,nd-1
         do jv=0,numberOfPolarizationVectors2-1
          pc = jv*nd+n 
          p2(j1,j2,j3,pc) = ( -( p2(j1-2*js1,j2-2*js2,j3,pc)+p2(j1+2*js1,j2+2*js2,j3,pc) ) \
                           +4.*( p2(j1-  js1,j2-  js2,j3,pc)+p2(j1+  js1,j2+  js2,j3,pc) ) \
                           +10.* p2(j1,j2,j3,pc) )/16.

           ! p2(j1,j2,j3,pc) = ( -( p2(j1-2,j2,j3,pc)+p2(j1+2,j2,j3,pc)+p2(j1,j2-2,j3,pc)+p2(j1,j2+2,j3,pc) ) \
           !                 +4.*( p2(j1-1,j2,j3,pc)+p2(j1+1,j2,j3,pc)+p2(j1,j2-1,j3,pc)+p2(j1,j2+1,j3,pc) ) \
           !                 +20.* p2(j1,j2,j3,pc) )/32.
         end do
       end do
     end if
    endLoopsMask2d()
   else
     write(*,'("cgmx:interface3d:smoothInterfaceP: finish me 3D")') 
     stop 2975
   end if

 end if
#endMacro

! ----------------------------------------------------------------------------------
!  Macro:
!    --- perturb a component of E by delta
! ----------------------------------------------------------------------------------
#beginMacro perturbComponent(n2,delta)
 if( nd.eq.2 .or. orderOfAccuracy.eq.2 )then
  if( n2.lt.nd )then
    u1(i1-is1,i2-is2,i3-is3,ex+n2   )=u1(i1-is1,i2-is2,i3-is3,ex+n2   )+(delta)
  else if( n2.lt.2*nd )then
    u2(j1-js1,j2-js2,j3-js3,ex+n2-nd)= u2(j1-js1,j2-js2,j3-js3,ex+n2-nd)+(delta)
  else if( n2.lt.3*nd )then
    u1(i1-2*is1,i2-2*is2,i3-2*is3,ex+n2-2*nd)=u1(i1-2*is1,i2-2*is2,i3-2*is3,ex+n2-2*nd)+(delta)
  else if( n2.lt.4*nd )then
    u2(j1-2*js1,j2-2*js2,j3-2*js3,ex+n2-3*nd)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ex+n2-3*nd)+(delta)
  else
    ! this should not happen
    stop 6363
  end if
 else
  ! 3D order 4 has a different ordering from 3D order 2
  !  (Ex1,Ey1,Ez1)(-1) , (Ex1,Ey1,Ez1)(-2), (Ex2,Ey2,Ez2)(-1) , (Ex2,Ey2,Ez2)(-2), 
  if( n2.lt.nd )then
    u1(i1-is1,i2-is2,i3-is3,ex+n2   )=u1(i1-is1,i2-is2,i3-is3,ex+n2   )+(delta)
  else if( n2.lt.2*nd )then
    u1(i1-2*is1,i2-2*is2,i3-2*is3,ex+n2-nd)=u1(i1-2*is1,i2-2*is2,i3-2*is3,ex+n2-nd)+(delta)
  else if( n2.lt.3*nd )then
    u2(j1-js1,j2-js2,j3-js3,ex+n2-2*nd)= u2(j1-js1,j2-js2,j3-js3,ex+n2-2*nd)+(delta)
  else if( n2.lt.4*nd )then
    u2(j1-2*js1,j2-2*js2,j3-2*js3,ex+n2-3*nd)=u2(j1-2*js1,j2-2*js2,j3-2*js3,ex+n2-3*nd)+(delta)
  else
    ! this should not happen
    stop 6363
  end if
 end if
#endMacro

! ----------------------------------------------------------------------------------
!  Macro:
!    --- check matrix coefficients by delta function approach ----
! Input:
!   am : matrix of coefficients is stored here 
!   evalInterfaceEquations : macro that evaluates the interface equations
! ----------------------------------------------------------------------------------
#beginMacro checkCoefficients(i1,i2,i3, j1,j2,j3, numberOfEquations,am,evalInterfaceEquations )

  ! hw1 = half stencil width
  hw1=orderOfAccuracy/2
  hw2=hw1
  if( nd.eq.2 )then
    hw3=0
  else
    hw3=hw1
  end if

  write(*,'("CHECK-COEFF: i1,i2,i3=",3i3," hw1,hw2,hw3=",3i2)') i1,i2,i3,hw1,hw2,hw3

  ! First eval equations with no pertutbation --> save in f0 
  evalInterfaceEquations()
  do n1=0,numberOfEquations-1
   f0(n1)=f(n1)
  end do

  delta=1.  ! perturb E by this amount 
  do n2=0,numberOfEquations-1

    ! pertub one component: 
    perturbComponent(n2,delta)

    evalInterfaceEquations()
    
    ! compute the difference
    do n1=0,numberOfEquations-1
     f(n1)=f(n1)-f0(n1)
    end do

    ! write(*,'(" u1x,v1y,u2x,v2y=",4(1pe10.2))') u1x,v1y, u2x,v2y
    if( .true. )then
     if( numberOfEquations.eq.4 )then
      write(*,'(" am(*,",i1,")=",4(1pe10.2)," diff(*)=",4(1pe8.1) )') n2,(am(n1,n2),n1=0,numberOfEquations-1),(am(n1,n2)-f(n1),n1=0,numberOfEquations-1)
     else if( numberOfEquations.eq.6 )then
      write(*,'(" am(*,",i1,")=",6(1pe10.2)," diff(*)=",6(1pe8.1) )') n2,(am(n1,n2),n1=0,numberOfEquations-1),(am(n1,n2)-f(n1),n1=0,numberOfEquations-1)
     else if( numberOfEquations.eq.8 )then
      write(*,'(" am(*,",i1,")=",8(1pe10.2)," diff(*)=",8(1pe8.1) )') n2,(am(n1,n2),n1=0,numberOfEquations-1),(am(n1,n2)-f(n1),n1=0,numberOfEquations-1)
     else if( numberOfEquations.eq.12 )then
      write(*,'(" am(*,",i1,")=",12(1pe10.2)," diff(*)=",12(1pe8.1) )') n2,(am(n1,n2),n1=0,numberOfEquations-1),(am(n1,n2)-f(n1),n1=0,numberOfEquations-1)
     else 
       stop 8181
     end if
    else
     if( numberOfEquations.eq.4 )then
      write(*,'(" am(*,",i1,")=",4(1pe10.2)," f(*)=",4(1pe10.2)," diff(*)=",4(1pe8.1) )') n2,(am(n1,n2),n1=0,numberOfEquations-1),(f(n1),n1=0,numberOfEquations-1),(am(n1,n2)-f(n1),n1=0,numberOfEquations-1)
     else if( numberOfEquations.eq.6 )then
      write(*,'(" am(*,",i1,")=",6(1pe10.2)," f(*)=",6(1pe10.2)," diff(*)=",6(1pe8.1) )') n2,(am(n1,n2),n1=0,numberOfEquations-1),(f(n1),n1=0,numberOfEquations-1),(am(n1,n2)-f(n1),n1=0,numberOfEquations-1)
     else if( numberOfEquations.eq.8 )then
      write(*,'(" am(*,",i1,")=",8(1pe10.2)," f(*)=",8(1pe10.2)," diff(*)=",8(1pe8.1) )') n2,(am(n1,n2),n1=0,numberOfEquations-1),(f(n1),n1=0,numberOfEquations-1),(am(n1,n2)-f(n1),n1=0,numberOfEquations-1)
     else if( numberOfEquations.eq.12 )then
      write(*,'(" am(*,",i1,")=",12(1pe10.2)," f(*)=",12(1pe10.2)," diff(*)=",12(1pe8.1) )') n2,(am(n1,n2),n1=0,numberOfEquations-1),(f(n1),n1=0,numberOfEquations-1),(am(n1,n2)-f(n1),n1=0,numberOfEquations-1)
     else 
       stop 8181
     end if
    end if

    do n1=0,numberOfEquations-1
      coeffDiff = max(coeffDiff,abs(am(n1,n2)-f(n1)))
    end do

    ! reset pertubation
    perturbComponent(n2,-delta)

  end do 

  ! restore 
  evalInterfaceEquations()

#endMacro

! ----------------------------------------------------------------------------------
!  Macro:
!    --- eval matrix coefficients by delta function approach ----
! Input:
!   am : matrix of coefficients is stored here 
!   evalInterfaceEquations : macro that evaluates the interface equations
! ----------------------------------------------------------------------------------
#beginMacro evalCoefficients(i1,i2,i3, j1,j2,j3, numberOfEquations,am,evalInterfaceEquations )

  ! hw1 = half stencil width
  hw1=orderOfAccuracy/2
  hw2=hw1
  if( nd.eq.2 )then
    hw3=0
  else
    hw3=hw1
  end if

  !  write(*,'("EVAL-COEFF: i1,i2,i3=",3i3," hw1,hw2,hw3=",3i2)') i1,i2,i3,hw1,hw2,hw3

  ! First eval equartions with no pertutbation --> save in f0 
  evalInterfaceEquations()
  do n1=0,numberOfEquations-1
   f0(n1)=f(n1)
  end do

  delta=1.  ! perturb E by this amount 
  do n2=0,numberOfEquations-1

    ! pertub one component: 
    perturbComponent(n2,delta)

    evalInterfaceEquations()
    
    ! compute the difference
    do n1=0,numberOfEquations-1
     f(n1)=f(n1)-f0(n1)
     am(n1,n2) = f(n1)
    end do

    ! reset pertubation
    perturbComponent(n2,-delta)

  end do 

  ! restore 
  evalInterfaceEquations()

#endMacro



! ----------------------------------------------------------------------------------
!  Macro:
!    --- save matrix coefficients by delta function approach ----
! Input:
!   am : matrix of coefficients is stored here 
!   evalInterfaceEquations : macro that evaluates the interface equations
! ----------------------------------------------------------------------------------
#beginMacro saveCoefficients(i1,i2,i3, j1,j2,j3, numberOfEquations,am,evalInterfaceEquations )


  ieqn = ieqn +1 ! counts equations (initialize to -1)

  ! hw1 = half stencil width
  hw1=orderOfAccuracy/2
  hw2=hw1
  if( nd.eq.2 )then
    hw3=0
  else
    hw3=hw1
  end if

  !  write(*,'("EVAL-COEFF: i1,i2,i3=",3i3," hw1,hw2,hw3=",3i2)') i1,i2,i3,hw1,hw2,hw3

  ! First eval equations with no perturbation --> save in f0 
  evalInterfaceEquations()
  do n1=0,numberOfEquations-1
   f0(n1)=f(n1)
  end do

  delta=1.  ! perturb E by this amount 
  ja=-1      ! counts nonzeros in each equation 
  ! ------ LEFT SIDE ---------
  if( is1.ne.0 )then
    i1a=i1-hw1*is1
    i1b=i1-is1
    i1c=is1
  else
    i1a=i1-hw1
    i1b=i1+hw1
    i1c=1
  end if
  if( is2.ne.0 )then
    i2a=i2-hw2*is2
    i2b=i2-is2
    i2c=is2
  else
    i2a=i2-hw2
    i2b=i2+hw2
    i2c=1
  end if
  if( is3.ne.0 )then
    i3a=i3-hw3*is3
    i3b=i3-is3
    i3c=is3
  else
    i3a=i3-hw3
    i3b=i3+hw3
    i3c=1
  end if


  do i3p=i3a,i3b,i3c
  do i2p=i2a,i2b,i2c
  do i1p=i1a,i1b,i1c
  do n2=0,nd-1
    ja = ja+1 ! next non-zero

    ! pertub one component: 
    u1(i1p,i2p,i3p,ex+n2)=u1(i1p,i2p,i3p,ex+n2)+(delta)

    evalInterfaceEquations()
    
    ! compute the difference
    do n1=0,numberOfEquations-1
     f(n1)=f(n1)-f0(n1)
     am(n1,ja) = f(n1)
    end do

    ! reset pertubation
    u1(i1p,i2p,i3p,ex+n2)=u1(i1p,i2p,i3p,ex+n2)-(delta)

  end do 
  end do 
  end do 
  end do 

  ! -------- RIGHT SIDE --------
  if( js1.ne.0 )then
    j1a=j1-hw1*js1
    j1b=j1-js1
    j1c=js1
  else
    j1a=j1-hw1
    j1b=j1+hw1
    j1c=1
  end if
  if( js2.ne.0 )then
    j2a=j2-hw2*js2
    j2b=j2-js2
    j2c=js2
  else
    j2a=j2-hw2
    j2b=j2+hw2
    j2c=1
  end if
  if( js3.ne.0 )then
    j3a=j3-hw3*js3
    j3b=j3-js3
    j3c=js3
  else
    j3a=j3-hw3
    j3b=j3+hw3
    j3c=1
  end if

  do j3p=j3a,j3b,j3c
  do j2p=j2a,j2b,j2c
  do j1p=j1a,j1b,j1c
  do n2=0,nd-1
    ja = ja+1 ! next non-zero

    ! pertub one component: 
    u2(j1p,j2p,j3p,ex+n2)=u2(j1p,j2p,j3p,ex+n2)+(delta)

    evalInterfaceEquations()
    
    ! compute the difference
    do n1=0,numberOfEquations-1
     f(n1)=f(n1)-f0(n1)
     am(n1,ja) = f(n1)
    end do

    ! reset pertubation
    u2(j1p,j2p,j3p,ex+n2)=u2(j1p,j2p,j3p,ex+n2)-(delta)

  end do 
  end do 
  end do 
  end do 

  ! output equations
  write(coeffFile,'("ieqn, i1,i2,i3, j1,j2,j3")')
  write(coeffFile,'(i10,2x,6(i10,1x))') ieqn, i1,i2,i3, j1,j2,j3
  write(coeffFile,'(" nnz=",i6," (num. non-zeros)")') ja+1
  do n1=0,numberOfEquations-1
    write(coeffFile,'(500(e24.16,1x))') (am(n1,k),k=0,ja)
  end do


  ! restore 
  evalInterfaceEquations()

#endMacro


  
