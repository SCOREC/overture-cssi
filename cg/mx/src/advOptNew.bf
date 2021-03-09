!
! Optimized advance routines for cgmx
!   This version includes routines for dispersion
!
! These next include files will define the macros that will define the difference approximations
! The actual macro is called below
#Include "defineDiffOrder2f.h"
#Include "defineDiffOrder4f.h"


! ---------------------------------------------------------------------------
! Macro : beginLoopsMask
! ---------------------------------------------------------------------------
#beginMacro beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  do i3=n3a,n3b
  do i2=n2a,n2b
  do i1=n1a,n1b
    if( mask(i1,i2,i3).gt.0 )then
#endMacro

! ---------------------------------------------------------------------------
! Macro : endLoopsMask
! ---------------------------------------------------------------------------
#beginMacro endLoopsMask()
    end if
  end do
  end do
  end do
#endMacro


#beginMacro loopse6(e1,e2,e3,e4,e5,e6)
if( useWhereMask.ne.0 )then
  do i3=n3a,n3b
  do i2=n2a,n2b
  do i1=n1a,n1b
    if( mask(i1,i2,i3).gt.0 )then
      e1
      e2
      e3
      e4
      e5
      e6
    end if
  end do
  end do
  end do
else
  do i3=n3a,n3b
  do i2=n2a,n2b
  do i1=n1a,n1b
      e1
      e2
      e3
      e4
      e5
      e6
  end do
  end do
  end do
end if
#endMacro


#beginMacro loopse9(e1,e2,e3,e4,e5,e6,e7,e8,e9)
if( useWhereMask.ne.0 )then
 do i3=n3a,n3b
 do i2=n2a,n2b
 do i1=n1a,n1b
  if( mask(i1,i2,i3).gt.0 )then
   e1
   e2
   e3
   e4
   e5
   e6
   e7
   e8
   e9
  end if
 end do
 end do
 end do
else
 do i3=n3a,n3b
 do i2=n2a,n2b
 do i1=n1a,n1b
  e1
  e2
  e3
  e4
  e5
  e6
  e7
  e8
  e9
 end do
 end do
 end do
end if
#endMacro

#beginMacro loopse18(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18)
if( useWhereMask.ne.0 )then
  do i3=n3a,n3b
  do i2=n2a,n2b
  do i1=n1a,n1b
    if( mask(i1,i2,i3).gt.0 )then
      e1
      e2
      e3
      e4
      e5
      e6
      e7
      e8
      e9
      e10
      e11
      e12
      e13
      e14
      e15
      e16
      e17
      e18
    end if
  end do
  end do
  end do
else
  do i3=n3a,n3b
  do i2=n2a,n2b
  do i1=n1a,n1b
      e1
      e2
      e3
      e4
      e5
      e6
      e7
      e8
      e9
      e10
      e11
      e12
      e13
      e14
      e15
      e16
      e17
      e18
  end do
  end do
  end do
end if
#endMacro


! This macro is used for variable dissipation in 2D
#beginMacro loopse6VarDis(e1,e2,e3,e4,e5,e6)
if( useWhereMask.ne.0 )then
  do i3=n3a,n3b
  do i2=n2a,n2b
  do i1=n1a,n1b
    if( varDis(i1,i2,i3).gt.0. .and. mask(i1,i2,i3).gt.0 )then
      e1
      e2
      e3
!     write(*,'(" i=",3i3," varDis=",e10.2," diss=",3e10.2)') i1,i2,i3,varDis(i1,i2,i3),dis(i1,i2,i3,ex),\
!         dis(i1,i2,i3,ey),dis(i1,i2,i3,ez)
    else
      e4
      e5
      e6
    end if
  end do
  end do
  end do
else
  do i3=n3a,n3b
  do i2=n2a,n2b
  do i1=n1a,n1b
    if( varDis(i1,i2,i3).gt.0. )then
      e1
      e2
      e3
    else
      e4
      e5
      e6
    end if
  end do
  end do
  end do
end if
#endMacro

! This macro is used for variable dissipation in 3D
#beginMacro loopse12VarDis(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12)
if( useWhereMask.ne.0 )then
  do i3=n3a,n3b
  do i2=n2a,n2b
  do i1=n1a,n1b
    if( varDis(i1,i2,i3).gt.0. .and. mask(i1,i2,i3).gt.0 )then
      e1
      e2
      e3
      e7
      e8
      e9
    else
      e4
      e5
      e6
      e10
      e11
      e12
    end if
  end do
  end do
  end do
else
  do i3=n3a,n3b
  do i2=n2a,n2b
  do i1=n1a,n1b
    if( varDis(i1,i2,i3).gt.0. )then
      e1
      e2
      e3
      e7
      e8
      e9
    else
      e4
      e5
      e6
      e10
      e11
      e12
    end if
  end do
  end do
  end do
end if
#endMacro

! This macro is used for variable dissipation in 3D
#beginMacro loopsVarDis3D(e1,e2,e3,e4,e5,e6,h1,h2,h3,h4,h5,h6)

 if( solveForE.ne.0 .and. solveForH.ne.0 )then
   loopse12VarDis(e1,e2,e3,e4,e5,e6,h1,h2,h3,h4,h5,h6)
 else if( solveForE.ne.0 ) then
   loopse6VarDis(e1,e2,e3,e4,e5,e6)
 else
   loopse6VarDis(h1,h2,h3,h4,h5,h6)
 end if

#endMacro



! Optionally add the forcing terms
#beginMacro loopsF2D(f1,f2,f3,e1,e2,e3,e4,e5,e6,e7,e8,e9)
if( addForcing.eq.0 )then
  loopse9(e1,e2,e3,e4,e5,e6,e7,e8,e9)
else
! add forcing to the first 3 equations
  loopse9(e1+f1,e2+f2,e3+f3,e4,e5,e6,e7,e8,e9)
end if
#endMacro

! Optionally add the forcing terms
! Optionally solve for E or H or both
#beginMacro loopsF3D(fe1,fe2,fe3,e1,e2,e3,e4,e5,e6,e7,e8,e9,fh1,fh2,fh3,h1,h2,h3,h4,h5,h6,h7,h8,h9)
if( addForcing.eq.0 )then

  if( solveForE.ne.0 .and. solveForH.ne.0 )then
    loopse18(e1,e2,e3,e4,e5,e6,e7,e8,e9,h1,h2,h3,h4,h5,h6,h7,h8,h9)
  else if( solveForE.ne.0 ) then
    loopse9(e1,e2,e3,e4,e5,e6,e7,e8,e9)
  else
    loopse9(h1,h2,h3,h4,h5,h6,h7,h8,h9)
  end if

else
! add forcing to the equations

  if( solveForE.ne.0 .and. solveForH.ne.0 )then
    loopse18(e1+fe1,e2+fe2,e3+fe3,e4,e5,e6,e7,e8,e9,h1+fh1,h2+fh2,h3+fh3,h4,h5,h6,h7,h8,h9)
  else if( solveForE.ne.0 ) then
    loopse9(e1+fe1,e2+fe2,e3+fe3,e4,e5,e6,e7,e8,e9)
  else
    loopse9(h1+fh1,h2+fh2,h3+fh3,h4,h5,h6,h7,h8,h9)
  end if

end if
#endMacro


! Optionally add the dissipation and or forcing terms
#beginMacro loopsF2DD(f1,f2,f3,e1,e2,e3,e4,e5,e6,e7,e8,e9)
if( addForcing.eq.0 .and. .not.addDissipation )then
  loopse9(e1,e2,e3,e4,e5,e6,e7,e8,e9)
else if( addForcing.ne.0 .and. .not.addDissipation )then
! add forcing to the first 3 equations
  loopse9(e1+f1,e2+f2,e3+f3,e4,e5,e6,e7,e8,e9)
else if( addForcing.eq.0 .and. addDissipation )then
! add dissipation to the first 3 equations
  loopse9(e1+dis(i1,i2,i3,ex),e2+dis(i1,i2,i3,ey),e3+dis(i1,i2,i3,hz),e4,e5,e6,e7,e8,e9)
else
!  add forcing and dissipation
  loopse9(e1+f1+dis(i1,i2,i3,ex),e2+f2+dis(i1,i2,i3,ey),e3+f3+dis(i1,i2,i3,hz),e4,e5,e6,e7,e8,e9)  
end if
#endMacro


! Optionally add the dissipation and or forcing terms
! Optionally solve for E or H or both
#beginMacro loopsF3DD(fe1,fe2,fe3,e1,e2,e3,e4,e5,e6,e7,e8,e9,fh1,fh2,fh3,h1,h2,h3,h4,h5,h6,h7,h8,h9)
if( addForcing.eq.0 .and. .not.addDissipation )then

  if( solveForE.ne.0 .and. solveForH.ne.0 )then
    ! stop 6654
    loopse18(e1,e2,e3,e4,e5,e6,e7,e8,e9,h1,h2,h3,h4,h5,h6,h7,h8,h9)
  else if( solveForE.ne.0 ) then
    loopse9(e1,e2,e3,e4,e5,e6,e7,e8,e9)
  else
    stop 9987
!***    loopse9(h1,h2,h3,h4,h5,h6,h7,h8,h9)
  end if

else if( addForcing.ne.0 .and. .not.addDissipation )then
! add forcing to the equations

  if( solveForE.ne.0 .and. solveForH.ne.0 )then
   ! stop 6654
    loopse18(e1+fe1,e2+fe2,e3+fe3,e4,e5,e6,e7,e8,e9,h1+fh1,h2+fh2,h3+fh3,h4,h5,h6,h7,h8,h9)
  else if( solveForE.ne.0 ) then
    loopse9(e1+fe1,e2+fe2,e3+fe3,e4,e5,e6,e7,e8,e9)
  else
    stop 9987
!***    loopse9(h1+fh1,h2+fh2,h3+fh3,h4,h5,h6,h7,h8,h9)
  end if

else if( addForcing.eq.0 .and. addDissipation )then
! add dissipation to the equations

  if( solveForE.ne.0 .and. solveForH.ne.0 )then
    ! stop 6654
    loopse18(e1+dis(i1,i2,i3,ex),e2+dis(i1,i2,i3,ey),e3+dis(i1,i2,i3,ez),e4,e5,e6,e7,e8,e9,h1+dis(i1,i2,i3,hx),h2+dis(i1,i2,i3,hy),h3+dis(i1,i2,i3,hz),h4,h5,h6,h7,h8,h9)
  else if( solveForE.ne.0 ) then
    loopse9(e1+dis(i1,i2,i3,ex),e2+dis(i1,i2,i3,ey),e3+dis(i1,i2,i3,ez),e4,e5,e6,e7,e8,e9)
  else
    stop 6654
!***    loopse9(h1+dis(i1,i2,i3,hx),h2+dis(i1,i2,i3,hy),h3+dis(i1,i2,i3,hz),h4,h5,h6,h7,h8,h9)
  end if

else
! add dissipation and forcing to the equations

  if( solveForE.ne.0 .and. solveForH.ne.0 )then
    ! stop 6654
    loopse18(e1+fe1+dis(i1,i2,i3,ex),e2+fe2+dis(i1,i2,i3,ey),e3+fe3+dis(i1,i2,i3,ez),e4,e5,e6,e7,e8,e9,h1+fh1+dis(i1,i2,i3,hx),h2+fh2+dis(i1,i2,i3,hy),h3+fh3+dis(i1,i2,i3,hz),h4,h5,h6,h7,h8,h9)
  else if( solveForE.ne.0 ) then
    loopse9(e1+fe1+dis(i1,i2,i3,ex),e2+fe2+dis(i1,i2,i3,ey),e3+fe3+dis(i1,i2,i3,ez),e4,e5,e6,e7,e8,e9)
  else
    stop 6654
!****    loopse9(h1+fh1+dis(i1,i2,i3,hx),h2+fh2+dis(i1,i2,i3,hy),h3+fh3+dis(i1,i2,i3,hz),h4,h5,h6,h7,h8,h9)
  end if

end if
#endMacro

! The next macro is used for curvilinear girds where the Laplacian term is precomputed.
#beginMacro loopsFC(e1,e2,e3,e4,e5,e6,h1,h2,h3,h4,h5,h6)

if( nd.eq.2 )then
  ! This next line assumes we solve for ex,ey and hz
  loopse9(e1,e2,e4,e5,h3,h6,,,)

else

  if( solveForE.ne.0 .and. solveForH.ne.0 )then
    loopse18(e1,e2,e3,e4,e5,e6,h1,h2,h3,h4,h5,h6,,,,,,)
  else if( solveForE.ne.0 ) then
    loopse9(e1,e2,e3,e4,e5,e6,,,)
  else
    loopse9(h1,h2,h3,h4,h5,h6,,,)
  end if

end if
#endMacro

! -------------------------------------------------------------------------------------------
! The next macro is used for curvilinear girds where the Laplacian term is precomputed.
! Optionally add dissipation too
! -------------------------------------------------------------------------------------------
#beginMacro loopsFCD(e1,e2,e3,e4,e5,e6,h1,h2,h3,h4,h5,h6)

if( .not.addDissipation )then
 if( nd.eq.2 )then
  ! This next line assumes we solve for ex,ey and hz
  loopse9(e1,e2,e4,e5,h3,h6,,,)
 else
  if( solveForE.ne.0 .and. solveForH.ne.0 )then
    loopse18(e1,e2,e3,e4,e5,e6,h1,h2,h3,h4,h5,h6,,,,,,)
  else if( solveForE.ne.0 ) then
    loopse9(e1,e2,e3,e4,e5,e6,,,)
  else
    loopse9(h1,h2,h3,h4,h5,h6,,,)
  end if
 end if
else ! add dissipation too
 if( nd.eq.2 )then
  ! This next line assumes we solve for ex,ey and hz
  loopse9(e1+dis(i1,i2,i3,ex),e2+dis(i1,i2,i3,ey),e4,e5,h3+dis(i1,i2,i3,hz),h6,,,)
 else
  if( solveForE.ne.0 .and. solveForH.ne.0 )then
    loopse18(e1+dis(i1,i2,i3,ex),e2+dis(i1,i2,i3,ey),e3+dis(i1,i2,i3,ez),e4,e5,e6,h1+dis(i1,i2,i3,hx),h2+dis(i1,i2,i3,hy),h3+dis(i1,i2,i3,hz),h4,h5,h6,,,,,,)
  else if( solveForE.ne.0 ) then
    loopse9(e1+dis(i1,i2,i3,ex),e2+dis(i1,i2,i3,ey),e3+dis(i1,i2,i3,ez),e4,e5,e6,,,)
  else
    loopse9(h1+dis(i1,i2,i3,hx),h2+dis(i1,i2,i3,hy),h3+dis(i1,i2,i3,hz),h4,h5,h6,,,)
  end if
 end if
end if
#endMacro


! -------------------------------------------------------------------------------------------
! The next macro is used for curvilinear girds where the Laplacian term is precomputed.
! Optionally add dissipation too
! -------------------------------------------------------------------------------------------
#beginMacro loopsFCD2DF(e1,e2,h3)
if( .not.addDissipation )then
  ! This next line assumes we solve for ex,ey and hz
  loopse9(e1,e2,h3,,,,,,)
else ! add dissipation too
  ! This next line assumes we solve for ex,ey and hz
  loopse9(e1+dis(i1,i2,i3,ex),e2+dis(i1,i2,i3,ey),h3+dis(i1,i2,i3,hz),,,,,,)
end if
#endMacro

! -------------------------------------------------------------------------------------------
! The next macro is used for curvilinear girds where the Laplacian term is precomputed.
! Optionally add dissipation too
! -------------------------------------------------------------------------------------------
#beginMacro loopsFCD3DF(e1,e2,e3,e4,e5,e6,h1,h2,h3,h4,h5,h6)

if( .not.addDissipation )then
  if( solveForE.ne.0 .and. solveForH.ne.0 )then
    loopse18(e1,e2,e3,e4,e5,e6,h1,h2,h3,h4,h5,h6,,,,,,)
  else if( solveForE.ne.0 ) then
    loopse9(e1,e2,e3,e4,e5,e6,,,)
  else
    loopse9(h1,h2,h3,h4,h5,h6,,,)
  end if
else ! add dissipation too
  if( solveForE.ne.0 .and. solveForH.ne.0 )then
    loopse18(e1+dis(i1,i2,i3,ex),e2+dis(i1,i2,i3,ey),e3+dis(i1,i2,i3,ez),e4,e5,e6,h1+dis(i1,i2,i3,hx),h2+dis(i1,i2,i3,hy),h3+dis(i1,i2,i3,hz),h4,h5,h6,,,,,,)
  else if( solveForE.ne.0 ) then
    loopse9(e1+dis(i1,i2,i3,ex),e2+dis(i1,i2,i3,ey),e3+dis(i1,i2,i3,ez),e4,e5,e6,,,)
  else
    loopse9(h1+dis(i1,i2,i3,hx),h2+dis(i1,i2,i3,hy),h3+dis(i1,i2,i3,hz),h4,h5,h6,,,)
  end if
end if
#endMacro

! ======================================================================================
!   Evaluate the TZ exact solution in 2D
! ======================================================================================
#beginMacro OGDERIV2D( ntd,nxd,nyd,nzd,i1,i2,i3,t, n,val)
  call ogDeriv(ep, ntd,nxd,nyd,nzd, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, n,val )
#endMacro

! ======================================================================================
!   Evaluate the TZ exact solution in 3D
! ======================================================================================
#beginMacro OGDERIV3D( ntd,nxd,nyd,nzd,i1,i2,i3,t, n,val)
  call ogDeriv(ep, ntd,nxd,nyd,nzd, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, n,val )
#endMacro

#defineMacro LAP2D2(U,i1,i2,i3,c) \
                       (U(i1+1,i2,i3,c)-2.*U(i1,i2,i3,c)+U(i1-1,i2,i3,c))*dxsqi\
                      +(U(i1,i2+1,i3,c)-2.*U(i1,i2,i3,c)+U(i1,i2-1,i3,c))*dysqi
#defineMacro LAP3D2(U,i1,i2,i3,c) \
                       (U(i1+1,i2,i3,c)-2.*U(i1,i2,i3,c)+U(i1-1,i2,i3,c))*dxsqi\
                      +(U(i1,i2+1,i3,c)-2.*U(i1,i2,i3,c)+U(i1,i2-1,i3,c))*dysqi\
                      +(U(i1,i2,i3+1,c)-2.*U(i1,i2,i3,c)+U(i1,i2,i3-1,c))*dzsqi

#defineMacro LAP2D2POW2(U,i1,i2,i3,c) ( 6.*U(i1,i2,i3,c)   \
                      - 4.*(U(i1+1,i2,i3,c)+U(i1-1,i2,i3,c))    \
                      +(U(i1+2,i2,i3,c)+U(i1-2,i2,i3,c)) )*dxi4 \
                      +( 6.*U(i1,i2,i3,c)    \
                      -4.*(U(i1,i2+1,i3,c)+U(i1,i2-1,i3,c))    \
                      +(U(i1,i2+2,i3,c)+U(i1,i2-2,i3,c)) )*dyi4  \
                      +( 8.*U(i1,i2,i3,c)     \
                      -4.*(U(i1+1,i2,i3,c)+U(i1-1,i2,i3,c)+U(i1,i2+1,i3,c)+U(i1,i2-1,i3,c))   \
                      +2.*(U(i1+1,i2+1,i3,c)+U(i1-1,i2+1,i3,c)+U(i1+1,i2-1,i3,c)+U(i1-1,i2-1,i3,c)) )*dxdyi2

#defineMacro LAP3D2POW2(U,i1,i2,i3,c) ( 6.*U(i1,i2,i3,c)   \
        - 4.*(U(i1+1,i2,i3,c)+U(i1-1,i2,i3,c))    \
            +(U(i1+2,i2,i3,c)+U(i1-2,i2,i3,c)) )*dxi4 \
       +(  +6.*U(i1,i2,i3,c)    \
         -4.*(U(i1,i2+1,i3,c)+U(i1,i2-1,i3,c))    \
            +(U(i1,i2+2,i3,c)+U(i1,i2-2,i3,c)) )*dyi4\
       +(  +6.*U(i1,i2,i3,c)    \
         -4.*(U(i1,i2,i3+1,c)+U(i1,i2,i3-1,c))    \
            +(U(i1,i2,i3+2,c)+U(i1,i2,i3-2,c)) )*dzi4\
        +(8.*U(i1,i2,i3,c)     \
         -4.*(U(i1+1,i2,i3,c)+U(i1-1,i2,i3,c)+U(i1,i2+1,i3,c)+U(i1,i2-1,i3,c))   \
         +2.*(U(i1+1,i2+1,i3,c)+U(i1-1,i2+1,i3,c)+U(i1+1,i2-1,i3,c)+U(i1-1,i2-1,i3,c)) )*dxdyi2 \
        +(8.*U(i1,i2,i3,c)     \
         -4.*(U(i1+1,i2,i3,c)+U(i1-1,i2,i3,c)+U(i1,i2,i3+1,c)+U(i1,i2,i3-1,c))   \
         +2.*(U(i1+1,i2,i3+1,c)+U(i1-1,i2,i3+1,c)+U(i1+1,i2,i3-1,c)+U(i1-1,i2,i3-1,c)) )*dxdzi2 \
        +(8.*U(i1,i2,i3,c)     \
         -4.*(U(i1,i2+1,i3,c)+U(i1,i2-1,i3,c)+U(i1,i2,i3+1,c)+U(i1,i2,i3-1,c))   \
         +2.*(U(i1,i2+1,i3+1,c)+U(i1,i2-1,i3+1,c)+U(i1,i2+1,i3-1,c)+U(i1,i2-1,i3-1,c)) )*dydzi2 

#defineMacro LAP2D4(U,i1,i2,i3,c) ( -30.*U(i1,i2,i3,c)     \
        +16.*(U(i1+1,i2,i3,c)+U(i1-1,i2,i3,c))     \
            -(U(i1+2,i2,i3,c)+U(i1-2,i2,i3,c)) )*dxsq12i + \
       ( -30.*U(i1,i2,i3,c)     \
        +16.*(U(i1,i2+1,i3,c)+U(i1,i2-1,i3,c))     \
            -(U(i1,i2+2,i3,c)+U(i1,i2-2,i3,c)) )*dysq12i

#defineMacro LAP3D4(U,i1,i2,i3,c) ( -30.*U(i1,i2,i3,c)     \
        +16.*(U(i1+1,i2,i3,c)+U(i1-1,i2,i3,c))     \
            -(U(i1+2,i2,i3,c)+U(i1-2,i2,i3,c)) )*dxsq12i + \
       ( -30.*U(i1,i2,i3,c)     \
        +16.*(U(i1,i2+1,i3,c)+U(i1,i2-1,i3,c))     \
            -(U(i1,i2+2,i3,c)+U(i1,i2-2,i3,c)) )*dysq12i+ \
       ( -30.*U(i1,i2,i3,c)      \
        +16.*(U(i1,i2,i3+1,c)+U(i1,i2,i3-1,c))      \
            -(U(i1,i2,i3+2,c)+U(i1,i2,i3-2,c)) )*dzsq12i

#defineMacro LAP2D6(U,i1,i2,i3,c) \
               c00lap2d6*U(i1,i2,i3,c)     \
              +c10lap2d6*(U(i1+1,i2,i3,c)+U(i1-1,i2,i3,c)) \
              +c01lap2d6*(U(i1,i2+1,i3,c)+U(i1,i2-1,i3,c)) \
              +c20lap2d6*(U(i1+2,i2,i3,c)+U(i1-2,i2,i3,c)) \
              +c02lap2d6*(U(i1,i2+2,i3,c)+U(i1,i2-2,i3,c)) \
              +c30lap2d6*(U(i1+3,i2,i3,c)+U(i1-3,i2,i3,c)) \
              +c03lap2d6*(U(i1,i2+3,i3,c)+U(i1,i2-3,i3,c))

#defineMacro LAP3D6(U,i1,i2,i3,c) \
               c000lap3d6*U(i1,i2,i3,c) \
              +c100lap3d6*(U(i1+1,i2,i3,c)+U(i1-1,i2,i3,c)) \
              +c010lap3d6*(U(i1,i2+1,i3,c)+U(i1,i2-1,i3,c)) \
              +c001lap3d6*(U(i1,i2,i3+1,c)+U(i1,i2,i3-1,c)) \
              +c200lap3d6*(U(i1+2,i2,i3,c)+U(i1-2,i2,i3,c)) \
              +c020lap3d6*(U(i1,i2+2,i3,c)+U(i1,i2-2,i3,c)) \
              +c002lap3d6*(U(i1,i2,i3+2,c)+U(i1,i2,i3-2,c)) \
              +c300lap3d6*(U(i1+3,i2,i3,c)+U(i1-3,i2,i3,c)) \
              +c030lap3d6*(U(i1,i2+3,i3,c)+U(i1,i2-3,i3,c)) \
              +c003lap3d6*(U(i1,i2,i3+3,c)+U(i1,i2,i3-3,c))


! ** evaluate the laplacian on the 9 points centred at (i1,i2,i3)
#beginMacro getLapValues2dOrder2(n)
 uLap(-1,-1,n) = uLaplacian22(i1-1,i2-1,i3,n)
 uLap( 0,-1,n) = uLaplacian22(i1  ,i2-1,i3,n)
 uLap(+1,-1,n) = uLaplacian22(i1+1,i2-1,i3,n)

 uLap(-1, 0,n) = uLaplacian22(i1-1,i2  ,i3,n)
 uLap( 0, 0,n) = uLaplacian22(i1  ,i2  ,i3,n)
 uLap(+1, 0,n) = uLaplacian22(i1+1,i2  ,i3,n)

 uLap(-1,+1,n) = uLaplacian22(i1-1,i2+1,i3,n)
 uLap( 0,+1,n) = uLaplacian22(i1  ,i2+1,i3,n)
 uLap(+1,+1,n) = uLaplacian22(i1+1,i2+1,i3,n)
#endMacro


! ** evaluate the square of the Laplacian for a component ****
#beginMacro evalLapSq2dOrder2(n)
 getLapValues2dOrder2(n)
 uLaprr2 = (uLap(+1, 0,n)-2.*uLap( 0, 0,n)+uLap(-1, 0,n))/(dr(0)**2)
 uLapss2 = (uLap( 0,+1,n)-2.*uLap( 0, 0,n)+uLap( 0,-1,n))/(dr(1)**2)
 uLaprs2 = (uLap(+1,+1,n)-uLap(-1,+1,n)-uLap(+1,-1,n)+uLap(-1,-1,n))/(4.*dr(0)*dr(1))
 uLapr2  = (uLap(+1, 0,n)-uLap(-1, 0,n))/(2.*dr(0))
 uLaps2  = (uLap( 0,+1,n)-uLap( 0,-1,n))/(2.*dr(1))

 uLapSq(n) =(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*uLaprr2\
        +2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*uLaprs2\
        +(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*uLapss2\
        +(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*uLapr2\
        +(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*uLaps2
 ! write(*,'(" n : uLaprr2,uLapss2,uLaprs2,uLapr2,uLaps2=",5f6.2)') uLaprr2,uLapss2,uLaprs2,uLapr2,uLaps2
#endMacro

! ** evaluate the square of the Laplacian for [ex,ey,hz] ****
#beginMacro getLapSq2dOrder2()
 evalLapSq2dOrder2(ex)
 evalLapSq2dOrder2(ey)
 evalLapSq2dOrder2(hz)
 ! write(*,'("addForcing,adc=",i2,f5.2,", uLapSq(n)=",3e9.2)') addForcing,adc,uLapSq(ex),uLapSq(ey),uLapSq(hz)
#endMacro

! ==== loops for curvilinear, with forcing, dissipation in 2D
! Optionally add the dissipation and or forcing terms
#beginMacro loopsFCD2D(expr0,f1,f2,f3,expr1,expr2,expr3)
if( addForcing.eq.0 .and. .not.addDissipation )then
 loopse9(expr0,expr1,expr2,expr3,,,,,)
else if( addForcing.ne.0 .and. .not.addDissipation )then
! add forcing to the first 3 equations
 loopse9(expr0,expr1+f1,expr2+f2,expr3+f3,,,,,)
else if( addForcing.eq.0 .and. addDissipation )then
! add dissipation to the first 3 equations
 loopse9(expr0,expr1+dis(i1,i2,i3,ex),expr2+dis(i1,i2,i3,ey),expr3+dis(i1,i2,i3,hz),,,,,)
else
!  add forcing and dissipation
 loopse9(expr0,expr1+f1+dis(i1,i2,i3,ex),expr2+f2+dis(i1,i2,i3,ey)+dis(i1,i2,i3,hz),expr3+f3,,,,,)  
end if
#endMacro

! ==== loops for curvilinear, with forcing, dissipation in 2D
! Optionally add the dissipation and or forcing terms
#beginMacro loopsFCD2DA(expr0,f1,expr1)
if( addForcing.eq.0 .and. .not.addDissipation )then
 loopse9(expr0,expr1,,,,,,,)
else if( addForcing.ne.0 .and. .not.addDissipation )then
! add forcing to the first 3 equations
 loopse9(expr0,expr1+f1,,,,,,,)
else if( addForcing.eq.0 .and. addDissipation )then
! add dissipation to the first 3 equations
 loopse9(expr0,expr1+dis(i1,i2,i3,ex),,,,,,,)
else
!  add forcing and dissipation
 loopse9(expr0,expr1+f1+dis(i1,i2,i3,ex),,,,,,,)
end if
#endMacro

! ===========================================================================================
! Macro: compute the coefficients in the sosup dissipation for curvilinear grids
! ===========================================================================================
#beginMacro getSosupDissipationCoeff2d(adxSosup)
 do dir=0,1
   ! diss-coeff ~= 1/(change in x along direction r(dir) )
   ! Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
   adxSosup(dir) = adSosup*uDotFactor*sqrt( rsxy(i1,i2,i3,dir,0)**2 + rsxy(i1,i2,i3,dir,1)**2 )/dr(dir) 
 end do
#endMacro

! ===========================================================================================
! Macro: Output some debug info for the first few time-steps
! ===========================================================================================
#beginMacro INFO(string)
if( t.le.3.*dt .and. debug.gt.1 )then
  write(*,'("advOPT>>>",string)')
end if
#endMacro

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (4th-order difference used with 2nd-order scheme) 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#defineMacro sosupDiss2d4(uDot,i1,i2,i3,n) \
              ( -6.*uDot(i1,i2,i3,n)     \
                +4.*(uDot(i1+1,i2,i3,n)+uDot(i1-1,i2,i3,n))     \
                   -(uDot(i1+2,i2,i3,n)+uDot(i1-2,i2,i3,n)) )*adxSosup(0) + \
              ( -6.*uDot(i1,i2,i3,n)     \
                +4.*(uDot(i1,i2+1,i3,n)+uDot(i1,i2-1,i3,n))     \
                   -(uDot(i1,i2+2,i3,n)+uDot(i1,i2-2,i3,n)) )*adxSosup(1)

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (6th-order difference used with 4th-order scheme)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#defineMacro sosupDiss2d6(uDot,i1,i2,i3,n) \
             (-20.*uDot(i1,i2,i3,n)     \
               +15.*(uDot(i1+1,i2,i3,n)+uDot(i1-1,i2,i3,n))     \
                -6.*(uDot(i1+2,i2,i3,n)+uDot(i1-2,i2,i3,n))     \
                   +(uDot(i1+3,i2,i3,n)+uDot(i1-3,i2,i3,n))  )*adxSosup(0) + \
              (-20.*uDot(i1,i2,i3,n)     \
               +15.*(uDot(i1,i2+1,i3,n)+uDot(i1,i2-1,i3,n))     \
                -6.*(uDot(i1,i2+2,i3,n)+uDot(i1,i2-2,i3,n))     \
                   +(uDot(i1,i2+3,i3,n)+uDot(i1,i2-3,i3,n))  )*adxSosup(1)

! ===========================================================================================
! Macro:     UPWIND DISSIPATION, RECTANGULAR, 2D, ORDER 2
! ===========================================================================================
#beginMacro updateUpwindDissipationRectangular2dOrder2()
 adxSosup(0)=cdSosupx*uDotFactor
 adxSosup(1)=cdSosupy*uDotFactor
 if( updateSolution.eq.1 .and. updateDissipation.eq.1 )then
  ! advance + sosup dissipation:
  INFO("FD22r-UP...update-solution-and-dissipation")
  adxSosup(0)=cdSosupx ! for D-minus-t do not scale by .5
  adxSosup(1)=cdSosupy
  loopsF2DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
           un(i1,i2,i3,ex)=maxwell2dr(i1,i2,i3,ex)+sosupDiss2d4(DmtU,i1,i2,i3,ex),\
           un(i1,i2,i3,ey)=maxwell2dr(i1,i2,i3,ey)+sosupDiss2d4(DmtU,i1,i2,i3,ey),\
           un(i1,i2,i3,hz)=maxwell2dr(i1,i2,i3,hz)+sosupDiss2d4(DmtU,i1,i2,i3,hz),,,,,,)
 else if( updateSolution.eq.1 )then
    ! advance to time n+1
  INFO("FD22r-UP...update-solution")
  loopsF2DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
          un(i1,i2,i3,ex)=maxwell2dr(i1,i2,i3,ex),\
          un(i1,i2,i3,ey)=maxwell2dr(i1,i2,i3,ey),\
          un(i1,i2,i3,hz)=maxwell2dr(i1,i2,i3,hz),,,,,,)
 else if( updateDissipation.eq.1 )then
  if( sosupDissipationOption.eq.0 .and. computeUt.eq.1 )then
   ! apply sosup dissipation to time n+1 (use precomputed v=uDot)
   INFO("FD22r-UP...update-un-with-dissipation-using-v")
   loopse6(un(i1,i2,i3,ex)=un(i1,i2,i3,ex)+sosupDiss2d4(v,i1,i2,i3,ex),\
           un(i1,i2,i3,ey)=un(i1,i2,i3,ey)+sosupDiss2d4(v,i1,i2,i3,ey),\
           un(i1,i2,i3,hz)=un(i1,i2,i3,hz)+sosupDiss2d4(v,i1,i2,i3,hz),,,) 
  else if( sosupDissipationOption.eq.0 )then
   ! apply sosup dissipation to time n+1
   INFO("FD22r-UP...update-un-with-dissipation")
   loopse6(un(i1,i2,i3,ex)=un(i1,i2,i3,ex)+sosupDiss2d4(DztU,i1,i2,i3,ex),\
           un(i1,i2,i3,ey)=un(i1,i2,i3,ey)+sosupDiss2d4(DztU,i1,i2,i3,ey),\
           un(i1,i2,i3,hz)=un(i1,i2,i3,hz)+sosupDiss2d4(DztU,i1,i2,i3,hz),,,) 
   else
   ! apply sosup dissipation to time n using times n-1 and n-2
   ! assume un holds u(t-2*dt) on input
   ! NOTE: the dissipation is added to u in a Gauss-Siedel fashion
   INFO("FD22r-UP...update-u-with-dissipation")
   loopse6(u(i1,i2,i3,ex)=u(i1,i2,i3,ex)+sosupDiss2d4(DzstU,i1,i2,i3,ex),\
           u(i1,i2,i3,ey)=u(i1,i2,i3,ey)+sosupDiss2d4(DzstU,i1,i2,i3,ey),\
           u(i1,i2,i3,hz)=u(i1,i2,i3,hz)+sosupDiss2d4(DzstU,i1,i2,i3,hz),,,) 
  end if
 else
   write(*,'("advOpt:FD22r-UP ERROR: unexpected option? sosupDissipationOption=",i2)') sosupDissipationOption
   stop 1010
 end if
#endMacro


! ===========================================================================================
! Macro:     UPWIND DISSIPATION, RECTANGULAR, 2D, ORDER 4
! ===========================================================================================
#beginMacro updateUpwindDissipationRectangular2dOrder4()
  if( useNewForcingMethod.ne.0 )then
   write(*,'(" finish me: useSosupDissipation && useNewForcingMethod")')
   stop 7733
  end if

 adxSosup(0)=cdSosupx*uDotFactor
 adxSosup(1)=cdSosupy*uDotFactor
 if( updateSolution.eq.1 .and. updateDissipation.eq.1 )then
  ! advance + sosup dissipation:
  INFO("FD44r-UP...update-solution-and-dissipation")
  adxSosup(0)=cdSosupx ! for D-minus-t do not scale by .5
  adxSosup(1)=cdSosupy
  loopsF2DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
         un(i1,i2,i3,ex)=maxwell2dr44me(i1,i2,i3,ex)+sosupDiss2d6(DmtU,i1,i2,i3,ex),\
         un(i1,i2,i3,ey)=maxwell2dr44me(i1,i2,i3,ey)+sosupDiss2d6(DmtU,i1,i2,i3,ey),\
         un(i1,i2,i3,hz)=maxwell2dr44me(i1,i2,i3,hz)+sosupDiss2d6(DmtU,i1,i2,i3,hz),,,,,,) 

 else if( updateSolution.eq.1 )then
   ! advance to time n+1
  INFO("FD44r-UP...update-solution")
  loopsF2DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
       un(i1,i2,i3,ex)=maxwell2dr44me(i1,i2,i3,ex),\
       un(i1,i2,i3,ey)=maxwell2dr44me(i1,i2,i3,ey),\
       un(i1,i2,i3,hz)=maxwell2dr44me(i1,i2,i3,hz),,,,,,)

 else if( updateDissipation.eq.1 )then

  if( sosupDissipationOption.eq.0 .and. computeUt.eq.1 )then
   ! apply sosup dissipation to time n+1 (use precomputed v=uDot)
   INFO("FD44r-UP...update-un-with-dissipation-using-v")
   loopse6(un(i1,i2,i3,ex)=un(i1,i2,i3,ex)+sosupDiss2d6(v,i1,i2,i3,ex),\
           un(i1,i2,i3,ey)=un(i1,i2,i3,ey)+sosupDiss2d6(v,i1,i2,i3,ey),\
           un(i1,i2,i3,hz)=un(i1,i2,i3,hz)+sosupDiss2d6(v,i1,i2,i3,hz),,,) 

  else if( sosupDissipationOption.eq.0 )then
   ! apply sosup dissipation to time n+1
   INFO("FD44r-UP...update-un-with-dissipation")
   loopse6(un(i1,i2,i3,ex)=un(i1,i2,i3,ex)+sosupDiss2d6(DztU,i1,i2,i3,ex),\
           un(i1,i2,i3,ey)=un(i1,i2,i3,ey)+sosupDiss2d6(DztU,i1,i2,i3,ey),\
           un(i1,i2,i3,hz)=un(i1,i2,i3,hz)+sosupDiss2d6(DztU,i1,i2,i3,hz),,,) 
   else
   ! apply sosup dissipation to time n using times n-1 and n-2
   ! assume un holds u(t-2*dt) on input
   ! NOTE: the dissipation is added to u in a Gauss-Siedel fashion
   INFO("FD44r-UP...update-u-with-dissipation")
   loopse6(u(i1,i2,i3,ex)=u(i1,i2,i3,ex)+sosupDiss2d6(DzstU,i1,i2,i3,ex),\
           u(i1,i2,i3,ey)=u(i1,i2,i3,ey)+sosupDiss2d6(DzstU,i1,i2,i3,ey),\
           u(i1,i2,i3,hz)=u(i1,i2,i3,hz)+sosupDiss2d6(DzstU,i1,i2,i3,hz),,,) 
  end if
 else
   write(*,'("advOpt:FD44r-UP ERROR: unexpected option? sosupDissipationOption=",i2)') sosupDissipationOption
   stop 1010
 end if
#endMacro



! ===========================================================================================
! Macro:     UPWIND DISSIPATION, CURVILINEAR, 2D, ORDER 2
! ===========================================================================================
#beginMacro updateUpwindDissipationCurvilinear2dOrder2()
 if( t.le.2.*dt )then
   write(*,'(" advOpt: FD22 + sosup-dissipation for curvilinear")')
 end if

 if( useNewForcingMethod.ne.0 )then
  write(*,'(" finish me: useSosupDissipation && useNewForcingMethod")')
  stop 7739
 end if

 ! FD22 (curvilinear grid) with Sosup (wide stencil dissiption)
 if( updateSolution.eq.1 .and. updateDissipation.eq.1 )then
  ! advance + sosup dissipation:
  ! note: forcing is already added to the rhs.
  INFO("FD22c-UP...update-solution-and-dissipation")
  uDotFactor=1. ! for D-minus-t do not scale by .5
  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
   getSosupDissipationCoeff2d(adxSosup)
   do m=0,2 ! ex, ey, hz
     ec=ex+m
     un(i1,i2,i3,ec)=maxwellc22(i1,i2,i3,ec)+sosupDiss2d4(DmtU,i1,i2,i3,ec)
   end do
  endLoopsMask()
  uDotFactor=.5 ! reset

else if( updateSolution.eq.1 )then
   ! advance to time n+1
  INFO("FD22c-UP...update-solution")
  ! note: forcing is already added to the rhs.
  if( updateSolution.eq.1 )then
   beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
    do m=0,2 ! ex, ey, hz
      ec=ex+m
      un(i1,i2,i3,ec)=maxwellc22(i1,i2,i3,ec)
    end do
   endLoopsMask()
  end if
 else if( updateDissipation.eq.1 )then
  ! --- add dissipation only ----
  if( sosupDissipationOption.eq.0 .and. computeUt.eq.1 )then
   ! apply sosup dissipation to time n+1 (use precomputed v=uDot)
   INFO("FD22c-UP...update-un-with-dissipation-using-v")
   beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
    getSosupDissipationCoeff2d(adxSosup)
    do m=0,2 ! ex, ey, hz
      ec=ex+m
      un(i1,i2,i3,ec)=un(i1,i2,i3,ec)+sosupDiss2d4(v,i1,i2,i3,ec)
    end do
   endLoopsMask()

  else if( sosupDissipationOption.eq.0 .and. computeUt.eq.0 )then

   ! apply sosup dissipation to time n+1
   INFO("FD22c-UP...update-un-with-dissipation")
   beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
    getSosupDissipationCoeff2d(adxSosup)
    do m=0,2 ! ex, ey, hz
      ec=ex+m
      un(i1,i2,i3,ec)=un(i1,i2,i3,ec)+sosupDiss2d4(DztU,i1,i2,i3,ec)
    end do
   endLoopsMask()
  else

   ! apply sosup dissipation to time n using times n-1 and n-2
   ! assume un holds u(t-2*dt) on input
   ! NOTE: the dissipation is added to u in a Gauss-Siedel fashion
   INFO("FD22c-UP...update-u-with-dissipation")
   beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
    getSosupDissipationCoeff2d(adxSosup)
    do m=0,2 ! ex, ey, hz
      ec=ex+m
      u(i1,i2,i3,ec)=u(i1,i2,i3,ec)+sosupDiss2d4(DzstU,i1,i2,i3,ec)
    end do
   endLoopsMask()
  end if
 else
   write(*,'("advOpt:FD22c-UP ERROR: unexpected option? sosupDissipationOption=",i2)') sosupDissipationOption
   stop 2020
 end if



#endMacro

! ===========================================================================================
! Macro:     UPWIND DISSIPATION, CURVILINEAR, 2D, ORDER 4
! ===========================================================================================
#beginMacro updateUpwindDissipationCurvilinear2dOrder4()
  if( t.le.2.*dt )then
    write(*,'(" advOpt: FD44 + upwind-dissipation for curvilinear")')
  end if

  if( useNewForcingMethod.ne.0 )then
   write(*,'(" finish me: FD44 + sosup-dissipation && useNewForcingMethod")')
   stop 4487
  end if

 ! FD44 (curvilinear grid) with Sosup (wide stencil dissiption)
 if( updateSolution.eq.1 .and. updateDissipation.eq.1 )then
  ! advance + sosup dissipation:
  ! note: forcing is already added to the rhs.
  INFO("FD44c-UP...update-solution-and-dissipation")
  uDotFactor=1. ! for D-minus-t do not scale by .5
  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
   getSosupDissipationCoeff2d(adxSosup)
   do m=0,2 ! ex, ey, hz
     ec=ex+m
     un(i1,i2,i3,ec)=maxwellc44me(i1,i2,i3,ec)+sosupDiss2d6(DmtU,i1,i2,i3,ec)
   end do
  endLoopsMask()
  uDotFactor=.5 ! reset
 else if( updateSolution.eq.1 )then
   ! advance to time n+1
  INFO("FD44c-UP...update-solution")
  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
   getSosupDissipationCoeff2d(adxSosup)
   do m=0,2 ! ex, ey, hz
     ec=ex+m
     un(i1,i2,i3,ec)=maxwellc44me(i1,i2,i3,ec)
   end do
  endLoopsMask()

 else if( updateDissipation.eq.1 )then
  ! ----- add dissipation only  ----

  if( sosupDissipationOption.eq.0 .and. computeUt.eq.1 )then
   ! apply sosup dissipation to time n+1 (use precomputed v=uDot)
   INFO("FD44c-UP...update-un-with-dissipation-using-v")
   beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
    getSosupDissipationCoeff2d(adxSosup)
    do m=0,2 ! ex, ey, hz
      ec=ex+m
      un(i1,i2,i3,ec)=un(i1,i2,i3,ec) + sosupDiss2d6(v,i1,i2,i3,ec)
    end do
   endLoopsMask()
  else if( sosupDissipationOption.eq.0 .and. computeUt.eq.0 )then
   ! apply sosup dissipation to time n+1
   INFO("FD44c-UP...update-un-with-dissipation")
   beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
    getSosupDissipationCoeff2d(adxSosup)
    do m=0,2 ! ex, ey, hz
      ec=ex+m
      ! Use D-zero = (un-um)
      un(i1,i2,i3,ec)=un(i1,i2,i3,ec) + sosupDiss2d6(DztU,i1,i2,i3,ec)
    end do
   endLoopsMask()
   else
   ! apply sosup dissipation to time n using times n-1 and n-2
   ! assume un holds u(t-2*dt) on input
   ! NOTE: the dissipation is added to u in a Gauss-Siedel fashion
   INFO("FD44c-UP...update-u-with-dissipation")
   beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
    getSosupDissipationCoeff2d(adxSosup)
    do m=0,2 ! ex, ey, hz
      ec=ex+m
      ! Use special D-zero = (u - un) , u=u(t), um=u(t-dt), un=u(t-2*dt)
      u(i1,i2,i3,ec)=u(i1,i2,i3,ec) + sosupDiss2d6(DzstU,i1,i2,i3,ec)
    end do
   endLoopsMask()
  end if
 else
   write(*,'("advOpt:FD44c-UP ERROR: unexpected option? sosupDissipationOption=",i2)') sosupDissipationOption
   stop 5050
 end if


#endMacro

! ========================================================================
! Macro: Getting forcing for GDM
!  Input:
!    ec : E component
!    pc : P component
!  Output:
!    fe : forcing for E is updated
!    fpv(iv) : forcing for polarization vector iv=0,1,2,...
! ========================================================================
#beginMacro getGDMForcing(ec,pc)
 if( addForcing.ne.0 )then
   ! fp = dtsq*f(i1,i2,i3,pc)
   if( forcingOption.eq.twilightZoneForcing )then

     if( nd.eq.2 )then
       OGDERIV2D( 0,0,0,0,i1,i2,i3,t, ec, e0  )
       OGDERIV2D( 1,0,0,0,i1,i2,i3,t, ec, e0t )
     else
       OGDERIV3D( 0,0,0,0,i1,i2,i3,t, ec, e0  )
       OGDERIV3D( 1,0,0,0,i1,i2,i3,t, ec, e0t )
     end if

     do iv=0,numberOfPolarizationVectors-1
       pce = pc+iv*nd
       if( nd.eq.2 )then
         OGDERIV2D( 0,0,0,0,i1,i2,i3,t, pce, p0  )
         OGDERIV2D( 1,0,0,0,i1,i2,i3,t, pce, p0t )
         OGDERIV2D( 2,0,0,0,i1,i2,i3,t, pce, p0tt)
       else
         OGDERIV3D( 0,0,0,0,i1,i2,i3,t, pce, p0  )
         OGDERIV3D( 1,0,0,0,i1,i2,i3,t, pce, p0t )
         OGDERIV3D( 2,0,0,0,i1,i2,i3,t, pce, p0tt)
       end if
       fe = fe + dtsq*alphaP*p0tt
       ! write(*,'(" fe,p0tt=",2e12.4)') fe,p0tt
       fpv(iv) = dtsq*( p0tt + b1v(iv)*p0t + b0v(iv)*p0 - a0v(iv)*e0 - a1v(iv)*e0t )
     end do

     !if( abs(fp-fp2).gt. 1.e-14 )then
     !  write(*,'(" (i1,i2)=",2i6," fp,fp2,diff=",3e12.4)') i1,i2,fp,fp2,fp-fp2
     !else
     !  fp=fp2
     !end if

   else
     do iv=0,numberOfPolarizationVectors-1
       fpv(iv)=0.
     end do
   end if
 end if
#endMacro

! ========================================================================
! Macro: Getting forcing for E and P in multilevel atomic (MLA) system (2nd order)
!  Input:
!    ec : E component
!    pc : P component
!  Output:
!    fe : forcing for E is updated
!    fpv(iv) : forcing for polarization vector iv=0,1,2,...
! ========================================================================
#beginMacro getEPMLAForcing(ec,pc)
  if( addForcing.ne.0 )then
    if( forcingOption.eq.twilightZoneForcing )then

      if( nd.eq.2 )then
        OGDERIV2D( 0,0,0,0,i1,i2,i3,t, ec, e0  )
        OGDERIV2D( 1,0,0,0,i1,i2,i3,t, ec, e0t )
        OGDERIV2D( 2,0,0,0,i1,i2,i3,t, ec, e0tt )
        OGDERIV2D( 0,2,0,0,i1,i2,i3,t, ec, e0xx )
        OGDERIV2D( 0,0,2,0,i1,i2,i3,t, ec, e0yy )
      else
        OGDERIV3D( 0,0,0,0,i1,i2,i3,t, ec, e0  )
        OGDERIV3D( 1,0,0,0,i1,i2,i3,t, ec, e0t )
        OGDERIV3D( 2,0,0,0,i1,i2,i3,t, ec, e0tt )
        OGDERIV3D( 0,2,0,0,i1,i2,i3,t, ec, e0xx )
        OGDERIV3D( 0,0,2,0,i1,i2,i3,t, ec, e0yy )
        OGDERIV3D( 0,0,0,2,i1,i2,i3,t, ec, e0zz )
      end if

      if( nd.eq.2 )then
        fe =  e0tt-csq * (e0xx + e0yy)
      else
        fe =  e0tt-csq * (e0xx + e0yy + e0zz)
      endif

      nce = pxc+nd*numberOfPolarizationVectors

      do iv=0,numberOfPolarizationVectors-1
        pce = pc+iv*nd
        if( nd.eq.2 )then
          OGDERIV2D( 0,0,0,0,i1,i2,i3,t, pce, p0  )
          OGDERIV2D( 1,0,0,0,i1,i2,i3,t, pce, p0t )
          OGDERIV2D( 2,0,0,0,i1,i2,i3,t, pce, p0tt)
        else
          OGDERIV3D( 0,0,0,0,i1,i2,i3,t, pce, p0  )
          OGDERIV3D( 1,0,0,0,i1,i2,i3,t, pce, p0t )
          OGDERIV3D( 2,0,0,0,i1,i2,i3,t, pce, p0tt)
        end if
        ! adding polarizations
        fe =  fe + alphaP*p0tt ! sum over P

        ! left hand side of gdm equations
        fpv(iv) = p0tt + b1v(iv)*p0t + b0v(iv)*p0
        do na = 0,numberOfAtomicLevels-1
          if( nd.eq.2 )then
            OGDERIV2D( 0,0,0,0,i1,i2,i3,t, nce+na, q0  )
          else
            OGDERIV3D( 0,0,0,0,i1,i2,i3,t, nce+na, q0  )
          end if
          ! adding \Delta N*E
          fpv(iv) = fpv(iv) - pnec(iv,na)*q0*e0
        enddo
      end do
    end if
  else
      fe = 0.
      do iv=0,numberOfPolarizationVectors-1
        fpv(iv)=0.
      end do
  end if
#endMacro


! ========================================================================
! Macro: Getting forcing for E and P in multilevel atomic (MLA) system (4th order)
!  Input:
!    ec : E component
!    pc : P component
!  Output:
!    fe : forcing for E is updated
!    fpv(iv) : forcing for polarization vector iv=0,1,2,...
! ========================================================================
#beginMacro getEPMLAForcing44(ec,pc)
  if( addForcing.ne.0 )then
    if( forcingOption.eq.twilightZoneForcing )then

      if( nd.eq.2 )then
        OGDERIV2D( 0,0,0,0,i1,i2,i3,t, ec, e0  )
        OGDERIV2D( 1,0,0,0,i1,i2,i3,t, ec, e0t )
        OGDERIV2D( 2,0,0,0,i1,i2,i3,t, ec, e0tt )
        OGDERIV2D( 3,0,0,0,i1,i2,i3,t, ec, e0ttt )
        OGDERIV2D( 4,0,0,0,i1,i2,i3,t, ec, e0tttt )
        OGDERIV2D( 0,1,0,0,i1,i2,i3,t, ec, e0x )
        OGDERIV2D( 0,2,0,0,i1,i2,i3,t, ec, e0xx )
        OGDERIV2D( 1,2,0,0,i1,i2,i3,t, ec, e0xxt )
        OGDERIV2D( 0,0,1,0,i1,i2,i3,t, ec, e0y )
        OGDERIV2D( 0,0,2,0,i1,i2,i3,t, ec, e0yy )
        OGDERIV2D( 1,0,2,0,i1,i2,i3,t, ec, e0yyt )
        OGDERIV2D( 2,2,0,0,i1,i2,i3,t, ec, e0xxtt )
        OGDERIV2D( 2,0,2,0,i1,i2,i3,t, ec, e0yytt )
        OGDERIV2D( 0,4,0,0,i1,i2,i3,t, ec, e0xxxx )
        OGDERIV2D( 0,2,2,0,i1,i2,i3,t, ec, e0xxyy )
        OGDERIV2D( 0,0,4,0,i1,i2,i3,t, ec, e0yyyy )
      else
        OGDERIV3D( 0,0,0,0,i1,i2,i3,t, ec, e0  )
        OGDERIV3D( 1,0,0,0,i1,i2,i3,t, ec, e0t )
        OGDERIV3D( 2,0,0,0,i1,i2,i3,t, ec, e0tt )
        OGDERIV3D( 3,0,0,0,i1,i2,i3,t, ec, e0ttt )
        OGDERIV3D( 4,0,0,0,i1,i2,i3,t, ec, e0tttt )
        OGDERIV3D( 0,1,0,0,i1,i2,i3,t, ec, e0x )
        OGDERIV3D( 0,2,0,0,i1,i2,i3,t, ec, e0xx )
        OGDERIV3D( 1,2,0,0,i1,i2,i3,t, ec, e0xxt )
        OGDERIV3D( 0,0,1,0,i1,i2,i3,t, ec, e0y )
        OGDERIV3D( 0,0,2,0,i1,i2,i3,t, ec, e0yy )
        OGDERIV3D( 1,0,2,0,i1,i2,i3,t, ec, e0yyt )
        OGDERIV3D( 0,0,0,1,i1,i2,i3,t, ec, e0z )
        OGDERIV3D( 0,0,0,2,i1,i2,i3,t, ec, e0zz )
        OGDERIV3D( 1,0,0,2,i1,i2,i3,t, ec, e0zzt )
        OGDERIV3D( 2,2,0,0,i1,i2,i3,t, ec, e0xxtt )
        OGDERIV3D( 2,0,2,0,i1,i2,i3,t, ec, e0yytt )
        OGDERIV3D( 2,0,0,2,i1,i2,i3,t, ec, e0zztt )
        OGDERIV3D( 0,4,0,0,i1,i2,i3,t, ec, e0xxxx )
        OGDERIV3D( 0,0,4,0,i1,i2,i3,t, ec, e0yyyy )
        OGDERIV3D( 0,0,0,4,i1,i2,i3,t, ec, e0zzzz )
        OGDERIV3D( 0,2,2,0,i1,i2,i3,t, ec, e0xxyy )
        OGDERIV3D( 0,2,0,2,i1,i2,i3,t, ec, e0xxzz )
        OGDERIV3D( 0,0,2,2,i1,i2,i3,t, ec, e0yyzz )
      end if

      if( nd.eq.2 )then
        fe =  e0tt-csq * (e0xx + e0yy)
        lapfe = e0xxtt+e0yytt-csq*(e0xxxx+2.*e0xxyy+e0yyyy) ! laplacian of fe
        fet =  e0ttt-csq * (e0xxt + e0yyt) ! time derivatives of fe
        fett = e0tttt-csq * (e0xxtt + e0yytt)
      else
        fe =  e0tt-csq * (e0xx + e0yy + e0zz)
        lapfe = e0xxtt+e0yytt+e0zztt-csq*(e0xxxx+e0yyyy+e0zzzz+2.*(e0xxyy+e0xxzz+e0yyzz))
        fet =  e0ttt-csq * (e0xxt + e0yyt + e0zzt)
        fett =  e0tttt-csq * (e0xxtt + e0yytt + e0zztt)
      endif

      nce = pxc+nd*numberOfPolarizationVectors

      do iv=0,numberOfPolarizationVectors-1
        pce = pc+iv*nd
        if( nd.eq.2 )then
          OGDERIV2D( 0,0,0,0,i1,i2,i3,t, pce, p0  )
          OGDERIV2D( 0,2,0,0,i1,i2,i3,t, pce, p0xx  )
          OGDERIV2D( 0,0,2,0,i1,i2,i3,t, pce, p0yy  )
          OGDERIV2D( 1,0,0,0,i1,i2,i3,t, pce, p0t )
          OGDERIV2D( 1,2,0,0,i1,i2,i3,t, pce, p0xxt )
          OGDERIV2D( 1,0,2,0,i1,i2,i3,t, pce, p0yyt )
          OGDERIV2D( 2,0,0,0,i1,i2,i3,t, pce, p0tt)
          OGDERIV2D( 3,0,0,0,i1,i2,i3,t, pce, p0ttt)
          OGDERIV2D( 4,0,0,0,i1,i2,i3,t, pce, p0tttt)
          OGDERIV2D( 2,2,0,0,i1,i2,i3,t, pce, p0xxtt)
          OGDERIV2D( 2,0,2,0,i1,i2,i3,t, pce, p0yytt)
          ! laplacians
          fp02x     = p0xxtt + b1v(iv)*p0xxt + b0v(iv)*p0xx
          fp02y     = p0yytt + b1v(iv)*p0yyt + b0v(iv)*p0yy
          fp02      = fp02x  + fp02y
        else
          OGDERIV3D( 0,0,0,0,i1,i2,i3,t, pce, p0  )
          OGDERIV3D( 0,2,0,0,i1,i2,i3,t, pce, p0xx  )
          OGDERIV3D( 0,0,2,0,i1,i2,i3,t, pce, p0yy  )
          OGDERIV3D( 0,0,0,2,i1,i2,i3,t, pce, p0zz  )
          OGDERIV3D( 1,0,0,0,i1,i2,i3,t, pce, p0t )
          OGDERIV3D( 1,2,0,0,i1,i2,i3,t, pce, p0xxt )
          OGDERIV3D( 1,0,2,0,i1,i2,i3,t, pce, p0yyt )
          OGDERIV3D( 1,0,0,2,i1,i2,i3,t, pce, p0zzt )
          OGDERIV3D( 2,0,0,0,i1,i2,i3,t, pce, p0tt)
          OGDERIV3D( 3,0,0,0,i1,i2,i3,t, pce, p0ttt)
          OGDERIV3D( 4,0,0,0,i1,i2,i3,t, pce, p0tttt)
          OGDERIV3D( 2,2,0,0,i1,i2,i3,t, pce, p0xxtt)
          OGDERIV3D( 2,0,2,0,i1,i2,i3,t, pce, p0yytt)
          OGDERIV3D( 2,0,0,2,i1,i2,i3,t, pce, p0zztt)
          ! laplacians
          fp02x     = p0xxtt + b1v(iv)*p0xxt + b0v(iv)*p0xx
          fp02y     = p0yytt + b1v(iv)*p0yyt + b0v(iv)*p0yy
          fp02z     = p0zztt + b1v(iv)*p0zzt + b0v(iv)*p0zz
          fp02      = fp02x  + fp02y + fp02z
        end if

        ! adding polarizations
        fe =  fe + alphaP*p0tt ! sum over P
        fet =  fet + alphaP*p0ttt
        fett =  fett + alphaP*p0tttt
        if( nd.eq.2 )then
            lapfe = lapfe + alphaP*p0xxtt + alphaP*p0yytt
          else
            lapfe = lapfe + alphaP*p0xxtt + alphaP*p0yytt + alphaP*p0zztt
        endif

        ! initialize laplacian of forcing functions of P
        lapfpv(iv) = fp02
        ! 
        fpv(iv)   = p0tt   + b1v(iv)*p0t   + b0v(iv)*p0
        fptv(iv)  = p0ttt  + b1v(iv)*p0tt  + b0v(iv)*p0t
        fpttv(iv) = p0tttt + b1v(iv)*p0ttt + b0v(iv)*p0tt
        do na = 0,numberOfAtomicLevels-1
          if( nd.eq.2 )then
            OGDERIV2D( 0,0,0,0,i1,i2,i3,t, nce+na, q0  )
            OGDERIV2D( 0,1,0,0,i1,i2,i3,t, nce+na, q0x  )
            OGDERIV2D( 0,2,0,0,i1,i2,i3,t, nce+na, q0xx  )
            OGDERIV2D( 0,0,1,0,i1,i2,i3,t, nce+na, q0y  )
            OGDERIV2D( 0,0,2,0,i1,i2,i3,t, nce+na, q0yy  )
            OGDERIV2D( 1,0,0,0,i1,i2,i3,t, nce+na, q0t  )
            OGDERIV2D( 2,0,0,0,i1,i2,i3,t, nce+na, q0tt  )
            fp02 = q0xx*e0+2.*q0x*e0x+q0*e0xx \
                 + q0yy*e0+2.*q0y*e0y+q0*e0yy
          else
            OGDERIV3D( 0,0,0,0,i1,i2,i3,t, nce+na, q0  )
            OGDERIV3D( 0,1,0,0,i1,i2,i3,t, nce+na, q0x  )
            OGDERIV3D( 0,2,0,0,i1,i2,i3,t, nce+na, q0xx  )
            OGDERIV3D( 0,0,1,0,i1,i2,i3,t, nce+na, q0y  )
            OGDERIV3D( 0,0,2,0,i1,i2,i3,t, nce+na, q0yy  )
            OGDERIV3D( 0,0,0,1,i1,i2,i3,t, nce+na, q0z  )
            OGDERIV3D( 0,0,0,2,i1,i2,i3,t, nce+na, q0zz  )
            OGDERIV3D( 1,0,0,0,i1,i2,i3,t, nce+na, q0t  )
            OGDERIV3D( 2,0,0,0,i1,i2,i3,t, nce+na, q0tt  )
            fp02 = q0xx*e0+2.*q0x*e0x+q0*e0xx \
                 + q0yy*e0+2.*q0y*e0y+q0*e0yy \
                 + q0zz*e0+2.*q0z*e0z+q0*e0zz
          end if
          
          ! adding \Delta N*E terms
          lapfpv(iv) = lapfpv(iv) - pnec(iv,na)*fp02
          fpv(iv) = fpv(iv) - pnec(iv,na)*q0*e0
          fptv(iv) = fptv(iv) - pnec(iv,na)*q0t*e0 \
                              - pnec(iv,na)*q0*e0t
          fpttv(iv) = fpttv(iv) - pnec(iv,na)*q0tt*e0 \
                             - 2.*pnec(iv,na)*q0t*e0t \
                                - pnec(iv,na)*q0*e0tt
        enddo
      end do
    end if
  else
      fe = 0.
      lapfe = 0.
      fet  = 0.
      fett = 0.
      do iv=0,numberOfPolarizationVectors-1
        fpv(iv)=0.
        fptv(iv)=0.
        fpttv(iv)=0.
        lapfpv(iv)=0.
      end do
  end if
#endMacro

! ========================================================================
! Macro: Getting forcing for N in multilevel atomic (MLA) system (2nd order)
!  Input:
!    na : N component
!  Output:
!    fnv(iv) : forcing for population density vector iv=0,1,2,...
! ========================================================================
#beginMacro getMLAForcing(na)
  if( addForcing.ne.0 )then

    if( forcingOption.eq.twilightZoneForcing ) then

      !
      ! for carrier population density
      !
      ! first place for nonlinear model
      nce = pxc+nd*numberOfPolarizationVectors

      !
      ! na-th level
      if( nd.eq.2 )then
        OGDERIV2D( 1,0,0,0,i1,i2,i3,t, nce+na, q0t )
        OGDERIV2D( 2,0,0,0,i1,i2,i3,t, nce+na, q0tt)
      else
        OGDERIV3D( 1,0,0,0,i1,i2,i3,t, nce+na, q0t )
        OGDERIV3D( 2,0,0,0,i1,i2,i3,t, nce+na, q0tt)
      end if
      ! initialize
      fnv(na)  = q0t ! forcing for \partial_tN_\ell = alpha_{\ell,k}N_k+\beta_{\ell,m}E\cdot\partial_tP_k
      fntv(na) = q0tt ! next derivative

      ! relaxation (alpha_{\ell,m})
      do iv=0,numberOfAtomicLevels-1
        if( nd.eq.2 )then
          OGDERIV2D( 0,0,0,0,i1,i2,i3,t, nce+iv, q0  )
          OGDERIV2D( 1,0,0,0,i1,i2,i3,t, nce+iv, q0t )
        else
          OGDERIV3D( 0,0,0,0,i1,i2,i3,t, nce+iv, q0  )
          OGDERIV3D( 1,0,0,0,i1,i2,i3,t, nce+iv, q0t )
        end if
        fnv(na)  = fnv(na)  - prc(na,iv)*q0
        fntv(na) = fntv(na) - prc(na,iv)*q0t
      enddo

      ! dot product (\beta_{\ell,k})
      do m=0,nd-1 ! loop over dim
        ! electric field
        if ( nd.eq.2 ) then
          OGDERIV2D( 0,0,0,0,i1,i2,i3,t, m, e0  )
          OGDERIV2D( 1,0,0,0,i1,i2,i3,t, m, e0t )
        else
          OGDERIV3D( 0,0,0,0,i1,i2,i3,t, m, e0  )
          OGDERIV3D( 1,0,0,0,i1,i2,i3,t, m, e0t )
        endif
        ! corresponding polarization vector
        do iv=0,numberOfPolarizationVectors-1  
          if( nd.eq.2 )then
            OGDERIV2D( 1,0,0,0,i1,i2,i3,t, pxc+m+iv*nd, p0t )
            OGDERIV2D( 2,0,0,0,i1,i2,i3,t, pxc+m+iv*nd, p0tt)
          else
            OGDERIV3D( 1,0,0,0,i1,i2,i3,t, pxc+m+iv*nd, p0t )
            OGDERIV3D( 2,0,0,0,i1,i2,i3,t, pxc+m+iv*nd, p0tt)
          end if
          fnv(na)  = fnv(na) - peptc(na,iv)*e0*p0t
          fntv(na) = fntv(na) - peptc(na,iv)*e0t*p0t - peptc(na,iv)*e0*p0tt
        enddo
      enddo
    end if
  ! no forcing
  else
      fnv(na)  = 0.
      fntv(na) = 0.
  end if
#endMacro

! ========================================================================
! Macro: Getting forcing for N in multilevel atomic (MLA) system (4th order)
!  Input:
!    na : N component
!  Output:
!    fnv(iv) : forcing for population density vector iv=0,1,2,...
! ========================================================================
#beginMacro getMLAForcing44(na)
  if( addForcing.ne.0 )then

    if( forcingOption.eq.twilightZoneForcing ) then

      !
      ! for carrier population density
      !
      ! first place for nonlinear model
      nce = pxc+nd*numberOfPolarizationVectors

      !
      ! na-th level
      if( nd.eq.2 )then
        OGDERIV2D( 1,0,0,0,i1,i2,i3,t, nce+na, q0t )
        OGDERIV2D( 2,0,0,0,i1,i2,i3,t, nce+na, q0tt)
        OGDERIV2D( 3,0,0,0,i1,i2,i3,t, nce+na, q0ttt)
        OGDERIV2D( 4,0,0,0,i1,i2,i3,t, nce+na, q0tttt)
      else
        OGDERIV3D( 1,0,0,0,i1,i2,i3,t, nce+na, q0t )
        OGDERIV3D( 2,0,0,0,i1,i2,i3,t, nce+na, q0tt)
        OGDERIV3D( 3,0,0,0,i1,i2,i3,t, nce+na, q0ttt)
        OGDERIV3D( 4,0,0,0,i1,i2,i3,t, nce+na, q0tttt)
      end if
      ! initialize
      fnv(na)    = q0t ! forcing for \partial_tN_\ell = alpha_{\ell,k}N_k+\beta_{\ell,m}E\cdot\partial_tP_k
      fntv(na)   = q0tt ! next derivative
      fnttv(na)  = q0ttt
      fntttv(na) = q0tttt

      ! relaxation (alpha_{\ell,m})
      do iv=0,numberOfAtomicLevels-1
        if( nd.eq.2 )then
          OGDERIV2D( 0,0,0,0,i1,i2,i3,t, nce+iv, q0  )
          OGDERIV2D( 1,0,0,0,i1,i2,i3,t, nce+iv, q0t )
          OGDERIV2D( 2,0,0,0,i1,i2,i3,t, nce+iv, q0tt)
          OGDERIV2D( 3,0,0,0,i1,i2,i3,t, nce+iv, q0ttt)
        else
          OGDERIV3D( 0,0,0,0,i1,i2,i3,t, nce+iv, q0  )
          OGDERIV3D( 1,0,0,0,i1,i2,i3,t, nce+iv, q0t )
          OGDERIV3D( 2,0,0,0,i1,i2,i3,t, nce+iv, q0tt)
          OGDERIV3D( 3,0,0,0,i1,i2,i3,t, nce+iv, q0ttt)
        end if
        fnv(na)    = fnv(na)    - prc(na,iv)*q0
        fntv(na)   = fntv(na)   - prc(na,iv)*q0t
        fnttv(na)  = fnttv(na)  - prc(na,iv)*q0tt
        fntttv(na) = fntttv(na) - prc(na,iv)*q0ttt
      enddo

      ! dot product (\beta_{\ell,k})
      do m=0,nd-1 ! loop over dim
        ! electric field
        if ( nd.eq.2 ) then
          OGDERIV2D( 0,0,0,0,i1,i2,i3,t, m, e0  )
          OGDERIV2D( 1,0,0,0,i1,i2,i3,t, m, e0t )
          OGDERIV2D( 2,0,0,0,i1,i2,i3,t, m, e0tt )
          OGDERIV2D( 3,0,0,0,i1,i2,i3,t, m, e0ttt )
        else
          OGDERIV3D( 0,0,0,0,i1,i2,i3,t, m, e0  )
          OGDERIV3D( 1,0,0,0,i1,i2,i3,t, m, e0t )
          OGDERIV3D( 2,0,0,0,i1,i2,i3,t, m, e0tt )
          OGDERIV3D( 3,0,0,0,i1,i2,i3,t, m, e0ttt )
        endif
        ! corresponding polarization vector
        do iv=0,numberOfPolarizationVectors-1  
          if( nd.eq.2 )then
            OGDERIV2D( 1,0,0,0,i1,i2,i3,t, pxc+m+iv*nd, p0t )
            OGDERIV2D( 2,0,0,0,i1,i2,i3,t, pxc+m+iv*nd, p0tt)
            OGDERIV2D( 3,0,0,0,i1,i2,i3,t, pxc+m+iv*nd, p0ttt)
            OGDERIV2D( 4,0,0,0,i1,i2,i3,t, pxc+m+iv*nd, p0tttt)
          else
            OGDERIV3D( 1,0,0,0,i1,i2,i3,t, pxc+m+iv*nd, p0t )
            OGDERIV3D( 2,0,0,0,i1,i2,i3,t, pxc+m+iv*nd, p0tt)
            OGDERIV3D( 3,0,0,0,i1,i2,i3,t, pxc+m+iv*nd, p0ttt )
            OGDERIV3D( 4,0,0,0,i1,i2,i3,t, pxc+m+iv*nd, p0tttt)
          end if
          fnv(na)  = fnv(na) - peptc(na,iv)*e0*p0t
          fntv(na) = fntv(na) - peptc(na,iv)*e0t*p0t - peptc(na,iv)*e0*p0tt
          fnttv(na) = fnttv(na) - peptc(na,iv)*e0tt*p0t \
                             - 2.*peptc(na,iv)*e0t*p0tt \
                                - peptc(na,iv)*e0*p0ttt
          fntttv(na) = fntttv(na) - peptc(na,iv)*e0ttt*p0t \
                               - 3.*peptc(na,iv)*e0tt*p0tt \
                               - 3.*peptc(na,iv)*e0t*p0ttt \
                                  - peptc(na,iv)*e0*p0tttt
        enddo
      enddo
   
   
    end if
  ! no forcing
  else
      fnv(na)    = 0.
      fntv(na)   = 0.
      fnttv(na)  = 0.
      fntttv(na) = 0.
  end if
#endMacro


! ========================================================================
! Macro: Getting forcing for fourth-order accurate GDM
!  Input:
!    ec : E component
!    pc : P component
!  Output:
!    fpv(iv) : forcing for polarization vector iv=0,1,2,...
!    More forcing vectors, etc.
! ========================================================================
#beginMacro getGDMForcing44(ec,pc,OGDERIV)

 f1 = 0.
 f4 = 0.
 f5 = 0.

 ! *wdh* March 4, 2018: initialize this:
 fe00=0.
 do iv=0,numberOfPolarizationVectors-1
   f2v(iv) = 0.
   f3v(iv) = 0.
   f6v(iv) = 0.
   ! *wdh* March 4, 2018: initialize this:
   fp00v(iv)=0. 
   ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f6=",3e16.8)') i1,i2,m,f2v(iv),f3v(iv),f6v(iv)
 end do

 if( addForcing.ne.0 )then


       OGDERIV( 0,0,0,0,i1,i2,i3,t, ec, e0    )
       OGDERIV( 1,0,0,0,i1,i2,i3,t, ec, e0t   )
       OGDERIV( 2,0,0,0,i1,i2,i3,t, ec, e0tt  )
       OGDERIV( 3,0,0,0,i1,i2,i3,t, ec, e0ttt )
       OGDERIV( 4,0,0,0,i1,i2,i3,t, ec, e0tttt)
       OGDERIV( 0,2,0,0,i1,i2,i3,t, ec, e0xx  )
       OGDERIV( 0,0,2,0,i1,i2,i3,t, ec, e0yy  )
       OGDERIV( 0,0,0,2,i1,i2,i3,t, ec, e0zz  )
       OGDERIV( 1,2,0,0,i1,i2,i3,t, ec, e0xxt )
       OGDERIV( 1,0,2,0,i1,i2,i3,t, ec, e0yyt )
       OGDERIV( 1,0,0,2,i1,i2,i3,t, ec, e0zzt )
       OGDERIV( 2,2,0,0,i1,i2,i3,t, ec, e0xxtt)
       OGDERIV( 2,0,2,0,i1,i2,i3,t, ec, e0yytt)
       OGDERIV( 2,0,0,2,i1,i2,i3,t, ec, e0zztt)
       OGDERIV( 0,4,0,0,i1,i2,i3,t, ec, e0xxxx)
       OGDERIV( 0,2,2,0,i1,i2,i3,t, ec, e0xxyy)
       OGDERIV( 0,2,0,2,i1,i2,i3,t, ec, e0xxzz)
       OGDERIV( 0,0,4,0,i1,i2,i3,t, ec, e0yyyy)
       OGDERIV( 0,0,2,2,i1,i2,i3,t, ec, e0yyzz)
       OGDERIV( 0,0,0,4,i1,i2,i3,t, ec, e0zzzz)

       pSum0tt    = 0
       pSum0ttt   = 0
       pSum0tttt  = 0
       pSum0xxtt  = 0
       pSum0yytt  = 0
       pSum0zztt  = 0

       do iv=0,numberOfPolarizationVectors-1
         pce = pc+iv*nd

         OGDERIV( 0,0,0,0,i1,i2,i3,t, pc, p0    )
         OGDERIV( 1,0,0,0,i1,i2,i3,t, pc, p0t   )
         OGDERIV( 2,0,0,0,i1,i2,i3,t, pc, p0tt  )
         OGDERIV( 3,0,0,0,i1,i2,i3,t, pc, p0ttt )
         OGDERIV( 4,0,0,0,i1,i2,i3,t, pc, p0tttt)
         OGDERIV( 0,2,0,0,i1,i2,i3,t, pc, p0xx  )
         OGDERIV( 0,0,2,0,i1,i2,i3,t, pc, p0yy  )
         OGDERIV( 0,0,0,2,i1,i2,i3,t, pc, p0zz  )
         OGDERIV( 1,2,0,0,i1,i2,i3,t, pc, p0xxt )
         OGDERIV( 1,0,2,0,i1,i2,i3,t, pc, p0yyt )
         OGDERIV( 1,0,0,2,i1,i2,i3,t, pc, p0zzt )
         OGDERIV( 2,2,0,0,i1,i2,i3,t, pc, p0xxtt)
         OGDERIV( 2,0,2,0,i1,i2,i3,t, pc, p0yytt)
         OGDERIV( 2,0,0,2,i1,i2,i3,t, pc, p0zztt)

         ! Derivatives of OG individual p eqn forcing terms fp
         fp00      = p0tt   + b1v(iv)*p0t   + b0v(iv)*p0   - a1v(iv)*e0t   - a0v(iv)*e0
         fp10      = p0ttt  + b1v(iv)*p0tt  + b0v(iv)*p0t  - a1v(iv)*e0tt  - a0v(iv)*e0t
         fp20      = p0tttt + b1v(iv)*p0ttt + b0v(iv)*p0tt - a1v(iv)*e0ttt - a0v(iv)*e0tt
         fp02x     = p0xxtt + b1v(iv)*p0xxt + b0v(iv)*p0xx - a1v(iv)*e0xxt - a0v(iv)*e0xx
         fp02y     = p0yytt + b1v(iv)*p0yyt + b0v(iv)*p0yy - a1v(iv)*e0yyt - a0v(iv)*e0yy
         fp02z     = p0zztt + b1v(iv)*p0zzt + b0v(iv)*p0zz - a1v(iv)*e0zzt - a0v(iv)*e0zz
         fp02      = fp02x  + fp02y + fp02z
         fp00v(iv) = fp00

         ! Building derivatives of full P summation terms
         pSum0tt   = pSum0tt    + p0tt
         pSum0ttt  = pSum0ttt   + p0ttt
         pSum0tttt = pSum0tttt  + p0tttt
         pSum0xxtt = pSum0xxtt  + p0xxtt
         pSum0yytt = pSum0yytt  + p0yytt
         pSum0zztt = pSum0zztt  + p0zztt

         ! Forcing on EACH individual p_ttt is (fp)_t - b1*fp :
         f2v(iv) = fp10 - b1v(iv)*fp00

         ! Forcing on EACH individual p eqn:
         f3v(iv) = dtsq * (fp00 + (dtsq/12.) * fp20)

         ! Forcing for EACH individual pxx equation at second order predictor is (xx derivative of fp)
         f6v(iv) = dtsq * fp02

       end do

       ! Build derivatives of E eqn forcing term
       fe00  = e0tt   + alphaP*pSum0tt   - csq*(e0xx   + e0yy    + e0zz  )
       fe10  = e0ttt  + alphaP*pSum0ttt  - csq*(e0xxt  + e0yyt   + e0zzt )
       fe20  = e0tttt + alphaP*pSum0tttt - csq*(e0xxtt + e0yytt  + e0zztt)
       fe02x = e0xxtt + alphaP*pSum0xxtt - csq*(e0xxxx + e0xxyy  + e0xxzz)
       fe02y = e0yytt + alphaP*pSum0yytt - csq*(e0xxyy + e0yyyy  + e0yyzz)
       fe02z = e0zztt + alphaP*pSum0zztt - csq*(e0xxzz + e0yyzz  + e0zzzz)
       fe02  = fe02x  + fe02y + fe02z
       fe    = fe00

       ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") e0xxxx,e0xxyy=",2e16.8)') i1,i2,m,e0xxxx,e0xxyy

       ! Forcing on Ettt = c^2 Etxx - alphaP Pttt is fet
       f1 = fe10

       ! Forcing on E equation is
       f4 = dtsq * (fe00 + (dtsq/12.) * (fe20 + csq * fe02))

       ! Forcing for Exx equation at second order predictor is (xx derivative of fe)
       f5 = dtsq * fe02

 end if

#endMacro

! ========================================================================
! Macro: Add forcing to Hz for GDM in 2D
! ========================================================================
#beginMacro addtForcingHz()
 if( addForcing.ne.0 )then
   if( forcingOption.eq.twilightZoneForcing )then
     ! TESTING ...
     ! OGDERIV2D( 0,0,0,0,i1,i2,i3,t+dt, hz, hz0)
     ! un(i1,i2,i3,hz)=hz0


     ! OGDERIV2D( ntd,nxd,nyd,nzd,i1,i2,i3,t, n,val)
     ! OGDERIV2D( 1,0,0,0,i1,i2,i3,t+.5*dt, hz, hz0t)
     ! OGDERIV2D( 0,1,0,0,i1,i2,i3,t+.5*dt, ey, ey0x)
     ! OGDERIV2D( 0,0,1,0,i1,i2,i3,t+.5*dt, ex, ex0y)
     ! fhz = ( ey0x - ex0y )/mu + hz0t

     OGDERIV2D( 1,0,0,0,i1,i2,i3,t-dt, hz, hz0t)
     OGDERIV2D( 0,1,0,0,i1,i2,i3,t-dt, ey, ey0x)
     OGDERIV2D( 0,0,1,0,i1,i2,i3,t-dt, ex, ex0y)
     fhz = -.5*( ( ey0x - ex0y )/mu + hz0t )

     OGDERIV2D( 1,0,0,0,i1,i2,i3,t   , hz, hz0t)
     OGDERIV2D( 0,1,0,0,i1,i2,i3,t   , ey, ey0x)
     OGDERIV2D( 0,0,1,0,i1,i2,i3,t   , ex, ex0y)
     fhz = fhz + 1.5*( ( ey0x - ex0y )/mu + hz0t )

     un(i1,i2,i3,hz)=un(i1,i2,i3,hz) + dt*fhz

   else
     un(i1,i2,i3,hz) = un(i1,i2,i3,hz) + dt*f(i1,i2,i3,hz) ! first order only **FIX ME**
   end if
 end if
#endMacro

! ===========================================================================================
! Macro:     DISPERSIVE: RECTANGULAR, 2D, ORDER 2
! ===========================================================================================
#beginMacro updateRectangular2dOrder2Dispersive()

 if( .true. .and. numberOfPolarizationVectors.eq.1 )then
  INFO("FD22r-dispersive");
  fp=0
  fe=0.

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! Advance Hz first:
    ! For now solve H_t = -(1/mu)*(  (E_y)_x - (E_x)_y )
    !   USE AB2 -- note: this is just a quadrature so stability is not an issue
    un(i1,i2,i3,hz) = u(i1,i2,i3,hz) -(dt/mu)*( 1.5*ux22r(i1,i2,i3,ey) -.5*umx22r(i1,i2,i3,ey) \
                                               -1.5*uy22r(i1,i2,i3,ex) +.5*umy22r(i1,i2,i3,ex) )

    addtForcingHz()

    do m=0,1
     pc=pxc+m
     ec=ex+m

     if( addForcing.ne.0 )then
       fe = dtsq*f(i1,i2,i3,ec)
       getGDMForcing(ec,pc)
       fp=fpv(0)
     end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     pv0 = p(i1,i2,i3,m)
     pvm0=pm(i1,i2,i3,m)

     rhsE = maxwell2dr(i1,i2,i3,ec) + alphaP*(2.*pv0-pvm0) + fe
     rhsP = 2.*pv0-pvm0 + .5*dt*( b1*pvm0 -a1*evm ) + dt*dt*( -b0*pv0 + a0*ev ) + fp

     deti = 1./(1.+.5*dt*(b1+alphaP*a1))
     un(i1,i2,i3,ec) = ((1.+.5*dt*b1)*rhsE -alphaP*rhsP)*deti
     pn(i1,i2,i3,m)  = (.5*a1*dt*rhsE            + rhsP)*deti

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") pvm0,pv0,,pn=",3e16.8)') i1,i2,m,pvm0,pv0,pn(i1,i2,i3,m)

    end do
  endLoopsMask()

 else

  ! ------- 2D DISPERSIVE RECTANGULAR MULTIPLE PV -------

  INFO("FD22r-dispersive-MULTI-PV");

  fe=0.
  ! -- first compute some coefficients ---
  beta=0.
  do iv=0,numberOfPolarizationVectors-1
    betav(iv) = 1./( 1.+.5*dt*b1v(iv) )
    beta = beta + .5*dt*a1v(iv)*betav(iv)
    fpv(iv)=0.  ! initialize if not used
  end do

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! Advance Hz first:
    ! For now solve H_t = -(1/mu)*(  (E_y)_x - (E_x)_y )
    !   USE AB2 -- note: this is just a quadrature so stability is not an issue
    un(i1,i2,i3,hz) = u(i1,i2,i3,hz) -(dt/mu)*( 1.5*ux22r(i1,i2,i3,ey) -.5*umx22r(i1,i2,i3,ey) \
                                               -1.5*uy22r(i1,i2,i3,ex) +.5*umy22r(i1,i2,i3,ex) )

    addtForcingHz()

    ! -- loop over components of the vector --
    do m=0,1
     pc=pxc+m
     ec=ex+m

     if( addForcing.ne.0 )then
       fe = dtsq*f(i1,i2,i3,ec)
       ! Compute fpv(iv) :
       getGDMForcing(ec,pc)
     end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     rhsP = 0.
     pSum=0.
     do iv=0,numberOfPolarizationVectors-1
       pv(iv) = p(i1,i2,i3,m+iv*nd)
       pvm(iv)=pm(i1,i2,i3,m+iv*nd)

       rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + fpv(iv)
       rhsP = rhsP + betav(iv)*rhspv(iv)
       pSum = pSum + 2.*pv(iv) - pvm(iv)


       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pv,pvm,rhspv,rhsP,pSum,fpv=",6e16.8)') i1,i2,m,pv(iv),pvm(iv),rhspv(iv),rhsP,pSum,fpv(iv) 


       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") ev,evm,rhspv,rhsP,pSum,fpv=",6e16.8)') i1,i2,m,ev,evm,a1v(iv),a0v(iv),b1v(iv),b0v(iv)


       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") first term=",1e16.8)') i1,i2,m, 2.*pv(iv)-pvm(iv)

       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") two terms",2e16.8)') i1,i2,m,.5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm),dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev )


    end do

     rhsE = maxwell2dr(i1,i2,i3,ec) + alphaP*( pSum - rhsP ) + fe

     evn = rhsE / (1.+ alphaP*beta)


     write(*,'(" (i2,i2,m)=(",i3,i3,i2,") rhsE,evn,fe=",3e16.8)') i1,i2,m,rhsE,evn,fe

     ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") rhsE,rhsP,fe00,evn=",4e16.8)') i1,i2,m,rhsE,rhsP,fe,evn    

     un(i1,i2,i3,ec) = evn
     do iv=0,numberOfPolarizationVectors-1
       pn(i1,i2,i3,m+iv*nd)  = betav(iv)*( .5*dt*a1v(iv)*evn + rhspv(iv) )

       ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pvn(iv)=",4e16.8)') i1,i2,m,pn(i1,i2,i3,m+iv*nd)   

     end do



   end do ! m=0,1
  endLoopsMask()

 end if
#endMacro

! ===========================================================================================
! Macro:     DISPERSIVE: CURVILINEAR, 2D, ORDER 2
! ===========================================================================================
#beginMacro updateCurvilinear2dOrder2Dispersive()
 if( addDissipation )then
   write(*,'(" -- finish me : dispersion and AD")')
   stop 8256
 end if
 if( useNewForcingMethod.ne.0 )then
  write(*,'(" finish me: dispersion && useNewForcingMethod")')
  stop 7733
 end if

 fp=0.
 fe=0.

 if( .true. .and. numberOfPolarizationVectors.eq.1 )then
  ! **** PROBABLY NO NEED FOR THIS SPECIAL CASE ****

  ! ------- 2D DISPERSIVE CURVILINEAR NP=1 ------
  INFO("FD22c-dispersive")

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! Advance Hz first:
    ! For now solve H_t = -(1/mu)*(  (E_y)_x - (E_x)_y )
    !   USE AB2 -- note: this is just a quadrature so stability is not an inssue
    un(i1,i2,i3,hz) = u(i1,i2,i3,hz) -(dt/mu)*( 1.5*ux22(i1,i2,i3,ey) -.5*umx22(i1,i2,i3,ey) \
                                               -1.5*uy22(i1,i2,i3,ex) +.5*umy22(i1,i2,i3,ex) )
    addtForcingHz()

    !  --- advance E and P ---
    do m=0,1
     pc=pxc+m
     ec=ex+m

     if( addForcing.ne.0 )then ! forcing in E equation already added to f 
       ! fe = dtsq*f(i1,i2,i3,ec)  ! this term is already included
       fe=0.
       getGDMForcing(ec,pc)
       fp=fpv(0)
     end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     pv0 = p(i1,i2,i3,m)
     pvm0=pm(i1,i2,i3,m)

     ! pv = u(i1,i2,i3,pc)
     ! pvm=um(i1,i2,i3,pc)

     rhsE = maxwellc22(i1,i2,i3,ec) + alphaP*(2.*pv0-pvm0) + fe
     rhsP = 2.*pv0-pvm0 + .5*dt*( b1*pvm0 -a1*evm ) + dt*dt*( -b0*pv0 + a0*ev ) + fp
     deti = 1./(1.+.5*dt*(b1+alphaP*a1))
     un(i1,i2,i3,ec) = ((1.+.5*dt*b1)*rhsE -alphaP*rhsP)*deti
     pn(i1,i2,i3,m)  = (.5*a1*dt*rhsE            + rhsP)*deti

     ! pn(i1,i2,i3,m) = un(i1,i2,i3,pc)
     ! un(i1,i2,i3,pc) = pn(i1,i2,i3,m)

    end do
  endLoopsMask()

 else

  ! ------- 2D DISPERSIVE CURVILINEAR MULTIPLE PV -------

  INFO("FD22c-dispersive-MULTI-PV");

  fe=0.
  ! -- first compute some coefficients ---
  beta=0.
  do iv=0,numberOfPolarizationVectors-1
    betav(iv) = 1./( 1.+.5*dt*b1v(iv) )
    beta = beta + .5*dt*a1v(iv)*betav(iv)
    fpv(iv)=0.  ! initialize if not used
  end do

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! Advance Hz first:
    ! For now solve H_t = -(1/mu)*(  (E_y)_x - (E_x)_y )
    !   USE AB2 -- note: this is just a quadrature so stability is not an issue
    un(i1,i2,i3,hz) = u(i1,i2,i3,hz) -(dt/mu)*( 1.5*ux22(i1,i2,i3,ey) -.5*umx22(i1,i2,i3,ey) \
                                               -1.5*uy22(i1,i2,i3,ex) +.5*umy22(i1,i2,i3,ex) )

    addtForcingHz()

    ! -- loop over components of the vector --
    do m=0,1
     pc=pxc+m
     ec=ex+m

     if( addForcing.ne.0 )then
       ! fe = dtsq*f(i1,i2,i3,ec)  ! this term is already include
       fe=0.
       ! Compute fpv(iv) :
       getGDMForcing(ec,pc)
     end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     rhsP = 0.
     pSum=0.
     do iv=0,numberOfPolarizationVectors-1
       pv(iv) = p(i1,i2,i3,m+iv*nd)
       pvm(iv)=pm(i1,i2,i3,m+iv*nd)

       rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + fpv(iv)
       rhsP = rhsP + betav(iv)*rhspv(iv)
       pSum = pSum + 2.*pv(iv) - pvm(iv)
     end do

     rhsE = maxwellc22(i1,i2,i3,ec) + alphaP*( pSum - rhsP ) + fe

     evn = rhsE / (1.+ alphaP*beta)
     un(i1,i2,i3,ec) = evn
     do iv=0,numberOfPolarizationVectors-1
       pn(i1,i2,i3,m+iv*nd)  = betav(iv)*( .5*dt*a1v(iv)*evn + rhspv(iv) )
    end do

   end do ! m=0,1
  endLoopsMask()

 end if
#endMacro




! ===========================================================================================
! Macro:     DISPERSIVE: RECTANGULAR, 3D, ORDER 2
!          *** THREE DIMENSIONS ***
! ===========================================================================================
#beginMacro updateRectangular3dOrder2Dispersive()

 if( .true. .and. numberOfPolarizationVectors.eq.1 )then
  INFO("FD22r-3D-dispersive");
  fp=0
  fe=0.

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    do m=0,nd-1
     pc=pxc+m
     ec=ex+m

     if( addForcing.ne.0 )then
       fe = dtsq*f(i1,i2,i3,ec)
       getGDMForcing(ec,pc)
       fp=fpv(0)
     end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     pv0 = p(i1,i2,i3,m)
     pvm0=pm(i1,i2,i3,m)

     rhsE = maxwell3dr(i1,i2,i3,ec) + alphaP*(2.*pv0-pvm0) + fe
     rhsP = 2.*pv0-pvm0 + .5*dt*( b1*pvm0 -a1*evm ) + dt*dt*( -b0*pv0 + a0*ev ) + fp

     deti = 1./(1.+.5*dt*(b1+alphaP*a1))
     un(i1,i2,i3,ec) = ((1.+.5*dt*b1)*rhsE -alphaP*rhsP)*deti
     pn(i1,i2,i3,m)  = (.5*a1*dt*rhsE            + rhsP)*deti

     if( .false. )then
       OGDERIV3D( 0,0,0,0,i1,i2,i3,t+dt, pc, p0  )
       write(*,'(" (i1,i2,m)=(",i3,i3,i2,") p,p0=",2e16.8)') i1,i2,m, pn(i1,i2,i3,m),p0
       pn(i1,i2,i3,m)=p0
     end if
     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") pvm0,pv0,,pn=",3e16.8)') i1,i2,m,pvm0,pv0,pn(i1,i2,i3,m)

    end do
  endLoopsMask()

 else

  ! ------- 3D DISPERSIVE RECTANGULAR MULTIPLE PV (ORDER 2) -------

  INFO("FD22r-3D-dispersive-MULTI-PV");

  fe=0.
  ! -- first compute some coefficients ---
  beta=0.
  do iv=0,numberOfPolarizationVectors-1
    betav(iv) = 1./( 1.+.5*dt*b1v(iv) )
    beta = beta + .5*dt*a1v(iv)*betav(iv)
    fpv(iv)=0.  ! initialize if not used
  end do

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! -- loop over components of the vector --
    do m=0,nd-1
     pc=pxc+m
     ec=ex+m

     if( addForcing.ne.0 )then
       fe = dtsq*f(i1,i2,i3,ec)
       ! Compute fpv(iv) :
       getGDMForcing(ec,pc)
     end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     rhsP = 0.
     pSum=0.
     do iv=0,numberOfPolarizationVectors-1
       pv(iv) = p(i1,i2,i3,m+iv*nd)
       pvm(iv)=pm(i1,i2,i3,m+iv*nd)

       rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + fpv(iv)
       rhsP = rhsP + betav(iv)*rhspv(iv)
       pSum = pSum + 2.*pv(iv) - pvm(iv)
     end do

     rhsE = maxwell3dr(i1,i2,i3,ec) + alphaP*( pSum - rhsP ) + fe

     evn = rhsE / (1.+ alphaP*beta)
     un(i1,i2,i3,ec) = evn
     do iv=0,numberOfPolarizationVectors-1
       pn(i1,i2,i3,m+iv*nd)  = betav(iv)*( .5*dt*a1v(iv)*evn + rhspv(iv) )
    end do

   end do ! m=0,1
  endLoopsMask()

 end if
#endMacro


! ===========================================================================================
! Macro:     DISPERSIVE: CURVILINEAR, 3D, ORDER 2
!          *** THREE DIMENSIONS ***
! ===========================================================================================
#beginMacro updateCurvilinear3dOrder2Dispersive()
 if( addDissipation )then
   write(*,'(" -- finish me : dispersion and AD")')
   stop 8256
 end if
 if( useNewForcingMethod.ne.0 )then
  write(*,'(" finish me: dispersion && useNewForcingMethod")')
  stop 7733
 end if

 fp=0.
 fe=0.

 if( .true. .and. numberOfPolarizationVectors.eq.1 )then
  ! **** PROBABLY NO NEED FOR THIS SPECIAL CASE ****

  ! ------- 3D DISPERSIVE CURVILINEAR NP=1 (ORDER 2) ------
  INFO("FD22c-3D-dispersive")

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    !  --- advance E and P ---
    do m=0,nd-1
     pc=pxc+m
     ec=ex+m

     if( addForcing.ne.0 )then ! forcing in E equation already added to f 
       ! fe = dtsq*f(i1,i2,i3,ec)  ! this term is already included
       fe=0.
       getGDMForcing(ec,pc)
       fp=fpv(0)
     end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     pv0 = p(i1,i2,i3,m)
     pvm0=pm(i1,i2,i3,m)

     ! pv = u(i1,i2,i3,pc)
     ! pvm=um(i1,i2,i3,pc)

     rhsE = maxwellc23(i1,i2,i3,ec) + alphaP*(2.*pv0-pvm0) + fe
     rhsP = 2.*pv0-pvm0 + .5*dt*( b1*pvm0 -a1*evm ) + dt*dt*( -b0*pv0 + a0*ev ) + fp
     deti = 1./(1.+.5*dt*(b1+alphaP*a1))
     un(i1,i2,i3,ec) = ((1.+.5*dt*b1)*rhsE -alphaP*rhsP)*deti
     pn(i1,i2,i3,m)  = (.5*a1*dt*rhsE            + rhsP)*deti

    end do
  endLoopsMask()

 else

  ! ------- 3D DISPERSIVE CURVILINEAR MULTIPLE PV -------

  INFO("FD22c-3D-dispersive-MULTI-PV");

  fe=0.
  ! -- first compute some coefficients ---
  beta=0.
  do iv=0,numberOfPolarizationVectors-1
    betav(iv) = 1./( 1.+.5*dt*b1v(iv) )
    beta = beta + .5*dt*a1v(iv)*betav(iv)
    fpv(iv)=0.  ! initialize if not used
  end do

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

     ! -- loop over components of the vector --
    do m=0,nd-1
     pc=pxc+m
     ec=ex+m

     if( addForcing.ne.0 )then
       ! fe = dtsq*f(i1,i2,i3,ec)  ! this term is already include
       fe=0.
       ! Compute fpv(iv) :
       getGDMForcing(ec,pc)
     end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     rhsP = 0.
     pSum=0.
     do iv=0,numberOfPolarizationVectors-1
       pv(iv) = p(i1,i2,i3,m+iv*nd)
       pvm(iv)=pm(i1,i2,i3,m+iv*nd)

       rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + fpv(iv)
       rhsP = rhsP + betav(iv)*rhspv(iv)
       pSum = pSum + 2.*pv(iv) - pvm(iv)
     end do

     rhsE = maxwellc23(i1,i2,i3,ec) + alphaP*( pSum - rhsP ) + fe

     evn = rhsE / (1.+ alphaP*beta)
     un(i1,i2,i3,ec) = evn
     do iv=0,numberOfPolarizationVectors-1
       pn(i1,i2,i3,m+iv*nd)  = betav(iv)*( .5*dt*a1v(iv)*evn + rhspv(iv) )
    end do

   end do ! m=0,..,2
  endLoopsMask()

 end if
#endMacro



! **********************************************************************************
! ******************** FOURTH ORDER DISPERSIVE MACROS ******************************
! **********************************************************************************

! ===========================================================================================
! Macro:     DISPERSIVE: RECTANGULAR, 2D, ORDER 4
! ===========================================================================================
#beginMacro updateRectangular2dOrder4Dispersive()

 if( .true. .and. numberOfPolarizationVectors.eq.1 )then
  INFO("FD44r-dispersive");

  write(*,'(" DISPERSIVE: RECTANGULAR, 2D, ORDER 4 (ONE PV) ")')

  fp=0
  fe=0.

  ! Coefficients in equation for Pttt
  ! Pttt = b1ttt*P_t + b0ttt*P + a0ttt*E + a1ttt*E_t + a2ttt*Ett
  b1ttt=b1*b1-b0
  b0ttt=b1*b0
  a0ttt=-a0*b1
  a1ttt=a0-a1*b1
  a2ttt=a1

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! Advance Hz first:
    ! For now solve H_t = -(1/mu)*(  (E_y)_x - (E_x)_y )
    !   USE AB2 -- note: this is just a quadrature so stability is not an issue
    un(i1,i2,i3,hz) = u(i1,i2,i3,hz) -(dt/mu)*( 1.5*ux42r(i1,i2,i3,ey) -.5*umx42r(i1,i2,i3,ey) \
                                               -1.5*uy42r(i1,i2,i3,ex) +.5*umy42r(i1,i2,i3,ex) )

    addtForcingHz()

    do m=0,1
     pc=pxc+m
     ec=ex+m

     !if( addForcing.ne.0 )then
     !  fe = dtsq*f(i1,i2,i3,ec)
     !  getGDMForcing(ec,pc)
     !  fp=fpv(0)
     !end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     pv0 = p(i1,i2,i3,m)
     pvm0=pm(i1,i2,i3,m)

     !rhsE = maxwell2dr44me(i1,i2,i3,ec) + alphaP*(2.*pv0-pvm0) + fe
     !rhsP = 2.*pv0-pvm0 + .5*dt*( b1*pvm0 -a1*evm ) + dt*dt*( -b0*pv0 + a0*ev ) + fp

     !deti = 1./(1.+.5*dt*(b1+alphaP*a1))
     !un(i1,i2,i3,ec) = ((1.+.5*dt*b1)*rhsE -alphaP*rhsP)*deti
     !pn(i1,i2,i3,m)  = (.5*a1*dt*rhsE            + rhsP)*deti

     ! New 4th Order code starts here
     ! pvx = px22r(i1,i2,i3,m)


     f1 = 0
     f2 = 0
     f3 = 0
     f4 = 0
     f5 = 0
     f6 = 0


     if( addForcing.ne.0 )then
       OGDERIV2D( 0,0,0,0,i1,i2,i3,t, pc, p0    )
       OGDERIV2D( 1,0,0,0,i1,i2,i3,t, pc, p0t   )
       OGDERIV2D( 2,0,0,0,i1,i2,i3,t, pc, p0tt  )
       OGDERIV2D( 3,0,0,0,i1,i2,i3,t, pc, p0ttt )
       OGDERIV2D( 4,0,0,0,i1,i2,i3,t, pc, p0tttt)
       OGDERIV2D( 0,2,0,0,i1,i2,i3,t, pc, p0xx  )
       OGDERIV2D( 0,0,2,0,i1,i2,i3,t, pc, p0yy  )
       OGDERIV2D( 1,2,0,0,i1,i2,i3,t, pc, p0xxt )
       OGDERIV2D( 1,0,2,0,i1,i2,i3,t, pc, p0yyt )
       OGDERIV2D( 2,2,0,0,i1,i2,i3,t, pc, p0xxtt)
       OGDERIV2D( 2,0,2,0,i1,i2,i3,t, pc, p0yytt)


       OGDERIV2D( 0,0,0,0,i1,i2,i3,t, ec, e0    )
       OGDERIV2D( 1,0,0,0,i1,i2,i3,t, ec, e0t   )
       OGDERIV2D( 2,0,0,0,i1,i2,i3,t, ec, e0tt  )
       OGDERIV2D( 3,0,0,0,i1,i2,i3,t, ec, e0ttt )
       OGDERIV2D( 4,0,0,0,i1,i2,i3,t, ec, e0tttt)
       OGDERIV2D( 0,2,0,0,i1,i2,i3,t, ec, e0xx  )
       OGDERIV2D( 0,0,2,0,i1,i2,i3,t, ec, e0yy  )
       OGDERIV2D( 1,2,0,0,i1,i2,i3,t, ec, e0xxt )
       OGDERIV2D( 1,0,2,0,i1,i2,i3,t, ec, e0yyt )
       OGDERIV2D( 2,2,0,0,i1,i2,i3,t, ec, e0xxtt)
       OGDERIV2D( 2,0,2,0,i1,i2,i3,t, ec, e0yytt)
       OGDERIV2D( 0,4,0,0,i1,i2,i3,t, ec, e0xxxx)
       OGDERIV2D( 0,2,2,0,i1,i2,i3,t, ec, e0xxyy)
       OGDERIV2D( 0,0,4,0,i1,i2,i3,t, ec, e0yyyy)


       fp00  = p0tt   + b1*p0t   + b0*p0   - a1*e0t   - a0*e0
       fp10  = p0ttt  + b1*p0tt  + b0*p0t  - a1*e0tt  - a0*e0t
       fp20  = p0tttt + b1*p0ttt + b0*p0tt - a1*e0ttt - a0*e0tt
       fp02x = p0xxtt + b1*p0xxt + b0*p0xx - a1*e0xxt - a0*e0xx
       fp02y = p0yytt + b1*p0yyt + b0*p0yy - a1*e0yyt - a0*e0yy
       fp02  = fp02x  + fp02y

       fe00  = e0tt   + alphaP*p0tt   - csq*(e0xx   + e0yy)
       fe10  = e0ttt  + alphaP*p0ttt  - csq*(e0xxt  + e0yyt)
       fe20  = e0tttt + alphaP*p0tttt - csq*(e0xxtt + e0yytt)
       fe02x = e0xxtt + alphaP*p0xxtt - csq*(e0xxxx + e0xxyy)
       fe02y = e0yytt + alphaP*p0yytt - csq*(e0xxyy + e0yyyy)
       fe02  = fe02x  + fe02y

       ! f1 does not appear in this version of the 4th order scheme
       f1 = 0

       ! Forcing on P_ttt is (fp)_t - b1*fp :
       ! f2 = fp(1,0,x) - b1 * fp(0,0,x);
       ! f2 = p0t - b1 * p0 (incorrect)
       f2 = fp10 - b1*fp00

       ! Forcing on P eqn:
       ! fp + (dt^2/12) * fp_tt
       ! f3 = dt^2 * (fp(0,0,x) + (dt^2/12) * (- a1 * fe(1,0,x) + fp(2,0,x)));  % FIX ME *wdh* Jan 8, 2018
       ! f3 = dtsq * (p0 + (dtsq/12.) * (-a1 * e0t + p0tt )) incorrect
       f3 = dtsq * (fp00 + (dtsq/12.) * (-a1 * fe10 + fp20))

       ! Forcing on E equation is
       ! fe + (dt^2)/12 * ( fe_tt + c^2*Delta(fe) )
       ! f4 = dt^2 * (fe(0,0,x) + (dt^2/12) * (fe(2,0,x) + c^2 * fe(0,2,x)));
       ! f4 = dtsq * (e0 + (dtsq/12.) * (e0tt + csq * (e0xx + e0yy))) incorrect
       f4 = dtsq * (fe00 + (dtsq/12.) * (fe20 + csq * fe02))


       ! Forcing for Exx equation at second order predictor is (xx derivative
       ! of fe)
       ! f5 = dt^2 * fe(0,2,x);
       ! f5 = dtsq * (e0xx + e0yy)
       f5 = dtsq * fe02

       ! Forcing for Pxx equation at second order predictor is (xx derivative
       ! of fp)
       ! f6 = dt^2 * fp(0,2,x);
       ! f6 = dtsq * (p0xx + p0yy)
       f6 = dtsq * fp02

     end if

     elap4   = lap2d4(i1,i2,i3,ec)
     elap4m  = lap2d4m(i1,i2,i3,ec)
     elapsq2 = lap2d2Pow2(i1,i2,i3,ec)
     elap2m  = lap2d2m(i1,i2,i3,ec)
     plap2   = plap2d2(i1,i2,i3,m)
     plap2m  = plap2d2m(i1,i2,i3,m)
     plap4   = plap2d4(i1,i2,i3,m)
     plap4m  = plap2d4m(i1,i2,i3,m)
     ! Matlab code starts here

     write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f6=",3e16.8)') i1,i2,m,f2,f3,f6    
     write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f1,f4,f5=",3e16.8)') i1,i2,m,f1,f4,f5

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") elap4,elap4m,elapsq2,elap2m,plap2,plap2m,plap4,plap4m=",4e16.8)') i1,i2,m,elap4,elap4m,elapsq2,elap2m
     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") plap2,plap2m,plap4,plap4m=",4e16.8)') i1,i2,m,plap2,plap2m,plap4,plap4m

     Exxxx   = elapsq2

     En   = ev
     Enm1 = evm

     ! Pn   = pv0
     Pnm1 = pvm0

     ! compute c^2*Delta( E_t ) using
     ! Et = (1.5*En - 2*Enm1 +.5*Enm2)/dt;

     Exxn   = elap4
     Exxnm1 = elap4m


     betaP = 1/( 1+ (b1+alphaP*a1)*dt/2 )
     Pxxn   = plap4
     Pxxnm1 = plap4m



     ! Approximation for P_ttt
     ! Pttt^*    = b1ttt*(Pnp1 -Pnm1)/(2.*dt) + b0ttt*pv0 + a0ttt*En + a1ttt*(Enp1-Enm1)/(2.*dt) + a2ttt*(Enp1-2*En+Enm1)/dt^2; 
     cStar=1  !  cStar=1 : use Pttt^*, 0: use old

     PtttStarRhs = b1ttt*(0.   -Pnm1)/(2.*dt) + b0ttt*pv0 + a0ttt*En + a1ttt*(0.  -Enm1)/(2.*dt) + a2ttt*(0.  -2*En+Enm1)/dtsq + f2

     ! Use second order update to predict new E and P for future space derivatives
     !   fexx = fe;
     !   fpxx = fp;

     rhsExxnp1Predict = 2.*Exxn - Exxnm1 + cdtsq*Exxxx + alphaP*(2.*Pxxn - Pxxnm1) + f5
     rhsPxxnp1Predict = 2.*Pxxn - Pxxnm1 + .5*dt*( b1*Pxxnm1 - a1*Exxnm1 ) + dtsq*( -b0*Pxxn + a0*Exxn ) + f6
     deti = 1./(1.+.5*dt*(b1+alphaP*a1))
     Exxnp1Predict = ((1.+.5*dt*b1)*rhsExxnp1Predict - alphaP*rhsPxxnp1Predict)*deti
     Pxxnp1Predict = (.5*a1*dt*rhsExxnp1Predict            + rhsPxxnp1Predict)*deti


     ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") Exxnp1Predict,Pxxnp1Predict=",2e16.8)') i1,i2,m,Exxnp1Predict,Pxxnp1Predict





     QxxStar = (Pxxnp1Predict - 2.*Pxxn + Pxxnm1)/dtsq
     Etxx = (Exxnp1Predict - Exxnm1)/(2.*dt)

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") QxxStar,Exxn,Pxxn=",3e16.8)') i1,i2,m,QxxStar,Exxnp1Predict,Pxxnp1Predict

     ! A = zeros(2,2);
     ! b = zeros(2,1);

     !------------------------------- P EQUATION ----------------------------------------
     ! D+tD-t P^n = -b1*P_t -b0*P + a0*E + a1*E_t + (dt^2/12)*( D+tD-t Q^n )
     !            = -b1*[D0t(P^n) ]  - b0*P^n + a0*E^n  + a1*[ D0t(E^n) - (dt^2/12)*E_tttStar ] ...
     !              + (dt^2/12)*( -b0*D+tD-t(P^n) + a0*D+tD-t(E^n) )
     !              + (b1+alphaP*a1)*(dt^2)/12*P_ttt^*


     !
     ! We can either use:
     !       Qt = (Q^{n+1}-Q^{n-1})/(2*dt)   <-- adds a term to A(1,3)
     !       Qt = ( 1.5*Qn-2.*Qnm1 +.5*Qnm2 )/dt;  <-- remove term from A(1,3) 

     A(1,1) =   -a1*dt/2.  -a0*dtsq/12.   + (b1+alphaP*a1)*(dtsq/12.)*dtsq*( -a1ttt/(2.*dt) - a2ttt/dtsq ) 
     A(1,2) = 1+ b1*dt/2.  +b0*dtsq/12.   + (b1+alphaP*a1)*(dtsq/12.)*dtsq*( -b1ttt/(2.*dt) )   

     b(1) = (2.*pv0-Pnm1)\
         + b1*dt*Pnm1/(2.)\
         - dtsq*b0*pv0\
         + dtsq*a0*En\
         - a1*(dt/2.)*Enm1\
         - a1*(dt**4/12.)*( csq*Etxx )\
         + (dtsq/12.)*( -b0*(0. -2.*pv0+Pnm1) + a0*(0. -2.*En+Enm1) )\
         + (b1+alphaP*a1)*(dtsq/12.)*dtsq*( PtttStarRhs )\
         + f3


     write(*,'(" (i1,i2,m)=(",i3,i3,i2,") rhspv=",1e16.8)') i1,i2,m,b(1)




      ! ------------------------------- E EQUATION ----------------------------------------
      ! D+tD-t E^n = c^2*Delta(E^n) -alphaP*D+tD-t P^n + (dt^2/(12)*( [c^2*Delta]^2 E^n  -alphaP*c^2*Qxx^*  )

     A(2,1) = 1.        ! coeff of E^{n+1}
     A(2,2) = alphaP    ! coeff of P^{n+1}

     b(2) = (2.*En-Enm1)\
            +csq*dtsq*Exxn\
            +(2.*pv0-Pnm1)*alphaP\
            +csq**2*dt**4*Exxxx/(12.)\
            -csq*dt**4*alphaP*QxxStar/(12.)\
            + f4

     write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b2,f4=",2e16.8)') i1,i2,m,b(2),f4


     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") Exxn,Exxxx,QxxStar=",6e16.8)') i1,i2,m,Exxn,Exxxx,QxxStar

     deti = 1./(A(1,1)*A(2,2) - A(1,2)*A(2,1))

     y(1) = deti *( A(2,2) * b(1) - A(1,2) * b(2))
     y(2) = deti *(-A(2,1) * b(1) + A(1,1) * b(2))

     write(*,'(" A(1,1),A(1,2),A(2,1),A(2,2),b(1),b(2)=",6e16.8)') A(1,1),A(1,2),A(2,1),A(2,2),b(1),b(2)


     !u(j,1) = y(1);
     !u(j,2) = y(2);


     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f4,f5,f6=",5e16.8)') i1,i2,m,f2,f3,f4,f5,f6

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") e2,e4,p2,p4=",4e16.8)') i1, i2, m, un(i1,i2,i3,ec), y(1), pn(i1,i2,i3,m), y(2)

     ! write(*,'(" En,Enm1,Exxn,pv0,Pnm1,alphaP,csq=",7e16.8)') En,Enm1,Exxn,pv0,Pnm1,alphaP,csq

     ! write(*,'(" dt,Exxxx,QxxStar=",7e16.8)') dt,Exxxx,QxxStar

     ! write(*,'(" Pxxnp1Predict,Pxxn,Pxxnm1=",7e16.8)') Pxxnp1Predict,Pxxn,Pxxnm1

     ! write(*,'(" rhsExxnp1Predict, rhsPxxnp1Predict, f5, f6=",4e16.8)') rhsExxnp1Predict, rhsPxxnp1Predict, f5, f6

     ! write(*,'(" Exxn,Exxnm1,cdtsq,Exxxx,Pxxn,Pxxnm1=",6e16.8)') Exxn,Exxnm1,cdtsq,Exxxx,Pxxn,Pxxnm1

     un(i1,i2,i3,ec) = y(1)
     pn(i1,i2,i3,m)  = y(2)

     end do
  endLoopsMask()

 else

  ! ------- 2D DISPERSIVE RECTANGULAR MULTIPLE PV (ORDER 4)------
  ! XYZ
  INFO("FD44r-dispersive-MULTI-PV");
  write(*,'(" DISPERSIVE: RECTANGULAR, 2D, ORDER 4, MULTIPLE PV:FINISH ME (YOU ARE HERE)  ")')

  fe=0.
  ! -- first compute some coefficients ---
  beta=0.
  do iv=0,numberOfPolarizationVectors-1
    betav(iv) = 1./( 1.+.5*dt*b1v(iv) )
    beta = beta + .5*dt*a1v(iv)*betav(iv)
    fpv(iv)=0.  ! initialize if not used
    b1tttv(iv)=b1v(iv)*b1v(iv)-b0v(iv)
    b0tttv(iv)=b1v(iv)*b0v(iv)
    a0tttv(iv)=-a0v(iv)*b1v(iv)
    a1tttv(iv)=a0v(iv)-a1v(iv)*b1v(iv)
    a2tttv(iv)=a1v(iv)
  end do

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! Advance Hz first:
    ! For now solve H_t = -(1/mu)*(  (E_y)_x - (E_x)_y )
    !   USE AB2 -- note: this is just a quadrature so stability is not an issue
    un(i1,i2,i3,hz) = u(i1,i2,i3,hz) -(dt/mu)*( 1.5*ux42r(i1,i2,i3,ey) -.5*umx42r(i1,i2,i3,ey) \
                                               -1.5*uy42r(i1,i2,i3,ex) +.5*umy42r(i1,i2,i3,ex) )

    addtForcingHz()

    ! -- loop over components of the vector --
    do m=0,1
     pc=pxc+m
     ec=ex+m

     f1 = 0
     f4 = 0
     f5 = 0

     do iv=0,numberOfPolarizationVectors-1
       f2v(iv) = 0
       f3v(iv) = 0
       f6v(iv) = 0
       ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f6=",3e16.8)') i1,i2,m,f2v(iv),f3v(iv),f6v(iv)
     end do

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f1,f4,f5=",3e16.8)') i1,i2,m,f1,f4,f5

     if( addForcing.ne.0 )then
       ! 2nd order forcing used in 2nd order prediction step
       ! fe = dtsq*f(i1,i2,i3,ec)
       ! Compute fpv(iv) :
       ! getGDMForcing(ec,pc)

       ! 4th order forcing used in 4th order scheme

       OGDERIV2D( 0,0,0,0,i1,i2,i3,t, ec, e0    )
       OGDERIV2D( 1,0,0,0,i1,i2,i3,t, ec, e0t   )
       OGDERIV2D( 2,0,0,0,i1,i2,i3,t, ec, e0tt  )
       OGDERIV2D( 3,0,0,0,i1,i2,i3,t, ec, e0ttt )
       OGDERIV2D( 4,0,0,0,i1,i2,i3,t, ec, e0tttt)
       OGDERIV2D( 0,2,0,0,i1,i2,i3,t, ec, e0xx  )
       OGDERIV2D( 0,0,2,0,i1,i2,i3,t, ec, e0yy  )
       OGDERIV2D( 1,2,0,0,i1,i2,i3,t, ec, e0xxt )
       OGDERIV2D( 1,0,2,0,i1,i2,i3,t, ec, e0yyt )
       OGDERIV2D( 2,2,0,0,i1,i2,i3,t, ec, e0xxtt)
       OGDERIV2D( 2,0,2,0,i1,i2,i3,t, ec, e0yytt)
       OGDERIV2D( 0,4,0,0,i1,i2,i3,t, ec, e0xxxx)
       OGDERIV2D( 0,2,2,0,i1,i2,i3,t, ec, e0xxyy)
       OGDERIV2D( 0,0,4,0,i1,i2,i3,t, ec, e0yyyy)

       pSum0tt    = 0
       pSum0ttt   = 0
       pSum0tttt  = 0
       pSum0xxtt  = 0
       pSum0yytt  = 0

       do iv=0,numberOfPolarizationVectors-1
         pce = pc+iv*nd

         OGDERIV2D( 0,0,0,0,i1,i2,i3,t, pce, p0    )
         OGDERIV2D( 1,0,0,0,i1,i2,i3,t, pce, p0t   )
         OGDERIV2D( 2,0,0,0,i1,i2,i3,t, pce, p0tt  )
         OGDERIV2D( 3,0,0,0,i1,i2,i3,t, pce, p0ttt )
         OGDERIV2D( 4,0,0,0,i1,i2,i3,t, pce, p0tttt)
         OGDERIV2D( 0,2,0,0,i1,i2,i3,t, pce, p0xx  )
         OGDERIV2D( 0,0,2,0,i1,i2,i3,t, pce, p0yy  )
         OGDERIV2D( 1,2,0,0,i1,i2,i3,t, pce, p0xxt )
         OGDERIV2D( 1,0,2,0,i1,i2,i3,t, pce, p0yyt )
         OGDERIV2D( 2,2,0,0,i1,i2,i3,t, pce, p0xxtt)
         OGDERIV2D( 2,0,2,0,i1,i2,i3,t, pce, p0yytt)

         ! Derivatives of OG individual p eqn forcing terms fp           
         fp00      = p0tt   + b1v(iv)*p0t   + b0v(iv)*p0   - a1v(iv)*e0t   - a0v(iv)*e0
         fp10      = p0ttt  + b1v(iv)*p0tt  + b0v(iv)*p0t  - a1v(iv)*e0tt  - a0v(iv)*e0t
         fp20      = p0tttt + b1v(iv)*p0ttt + b0v(iv)*p0tt - a1v(iv)*e0ttt - a0v(iv)*e0tt
         fp02x     = p0xxtt + b1v(iv)*p0xxt + b0v(iv)*p0xx - a1v(iv)*e0xxt - a0v(iv)*e0xx
         fp02y     = p0yytt + b1v(iv)*p0yyt + b0v(iv)*p0yy - a1v(iv)*e0yyt - a0v(iv)*e0yy
         fp02      = fp02x  + fp02y
         fp00v(iv) = fp00

         ! Building derivatives of full P summation terms
         pSum0tt   = pSum0tt    + p0tt
         pSum0ttt  = pSum0ttt   + p0ttt
         pSum0tttt = pSum0tttt  + p0tttt
         pSum0xxtt = pSum0xxtt  + p0xxtt
         pSum0yytt = pSum0yytt  + p0yytt

         ! Forcing on EACH individual p_ttt is (fp)_t - b1*fp :
         f2v(iv) = fp10 - b1v(iv)*fp00

         ! Forcing on EACH individual p eqn:
         f3v(iv) = dtsq * (fp00 + (dtsq/12.) * fp20)

         ! Forcing for EACH individual pxx equation at second order predictor is (xx derivative of fp)
         f6v(iv) = dtsq * fp02

         if( iv.eq.0 )then
           write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f6=",3e16.8)') i1,i2,m,f2v(iv),f3v(iv),f6v(iv)
         end if
       end do

       ! Build derivatives of E eqn forcing term
       fe00  = e0tt   + alphaP*pSum0tt   - csq*(e0xx   + e0yy)
       fe10  = e0ttt  + alphaP*pSum0ttt  - csq*(e0xxt  + e0yyt)
       fe20  = e0tttt + alphaP*pSum0tttt - csq*(e0xxtt + e0yytt)
       fe02x = e0xxtt + alphaP*pSum0xxtt - csq*(e0xxxx + e0xxyy)
       fe02y = e0yytt + alphaP*pSum0yytt - csq*(e0xxyy + e0yyyy)
       fe02  = fe02x  + fe02y
       fe    = fe00

       ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") e0xxxx,e0xxyy=",2e16.8)') i1,i2,m,e0xxxx,e0xxyy

       ! Forcing on Ettt = c^2 Etxx - alphaP Pttt is fet
       f1 = fe10

       ! Forcing on E equation is
       f4 = dtsq * (fe00 + (dtsq/12.) * (fe20 + csq * fe02))

       ! Forcing for Exx equation at second order predictor is (xx derivative of fe)
       f5 = dtsq * fe02

       write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f1,f4,f5=",3e16.8)') i1,i2,m,f1,f4,f5
     end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     rhsP = 0.
     pSum=0.

     ! Fourth order below

     elap4   = lap2d4(i1,i2,i3,ec)
     elap4m  = lap2d4m(i1,i2,i3,ec)
     elapsq2 = lap2d2Pow2(i1,i2,i3,ec)
     elap2m  = lap2d2m(i1,i2,i3,ec)

     rhsPxx = 0.
     pxxSum = 0.

     exxv  = elap4
     exxvm = elap4m

     ! First we do the second order prediction on (i) E, (ii) Exx, (iii) Pxx, (iv) Individual pk
     do iv=0,numberOfPolarizationVectors-1

       pv(iv)    = p(i1,i2,i3,m+iv*nd)
       pvm(iv)   = pm(i1,i2,i3,m+iv*nd)
       pxxv(iv)  = plap2d4(i1,i2,i3,m+iv*nd)
       pxxvm(iv) = plap2d4m(i1,i2,i3,m+iv*nd)

       ! rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + fp00v(iv)
       rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + dtsq*fp00v(iv)
       rhsP = rhsP + betav(iv)*rhspv(iv)
       pSum = pSum + 2.*pv(iv) - pvm(iv)


       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pv,pvm,rhspv,rhsP,pSum,fpv=",6e16.8)') i1,i2,m,pv(iv),pvm(iv),rhspv(iv),rhsP,pSum,fp00v(iv)

       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") ev,evm,rhspv,rhsP,pSum,fpv=",6e16.8)') i1,i2,m,ev,evm,a1v(iv),a0v(iv),b1v(iv),b0v(iv)


       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") first term=",1e16.8)') i1,i2,m,2.*pv(iv)-pvm(iv)

       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") two terms",2e16.8)') i1,i2,m,.5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ),dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev )


       rhspxxv(iv) = 2.*pxxv(iv)-pxxvm(iv) + .5*dt*( b1v(iv)*pxxvm(iv) -a1v(iv)*exxvm ) + dtSq*( -b0v(iv)*pxxv(iv) + a0v(iv)*exxv ) + f6v(iv)
       rhsPxx = rhsPxx + betav(iv)*rhspxxv(iv)
       pxxSum = pxxSum + 2.*pxxv(iv) - pxxvm(iv)




     end do

     rhsE   = maxwell2dr44me(i1,i2,i3,ec)     + alphaP*( pSum   - rhsP   ) + dtsq * fe00 
     rhsE   = 2.*ev   - evm   + cdtsq*elap4   + alphaP*( pSum   - rhsP   ) + dtsq * fe00

     rhsExx = 2.*exxv - exxvm + cdtsq*elapsq2 + alphaP*( pxxSum - rhsPxx ) + f5



     evn   = rhsE   / (1.+ alphaP*beta)
     exxvn = rhsExx / (1.+ alphaP*beta)



     write(*,'(" (i2,i2,m)=(",i3,i3,i2,") rhsE,evn,fe00=",3e16.8)') i1,i2,m,rhsE,evn,fe00


     ! Update x derivative of P
     Pxxn  = beta * exxvn + rhsPxx

     ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") exxvn,Pxxn=",2e16.8)') i1,i2,m,exxvn,Pxxn


     ! Now we predict individual pk (pvn) and their third time derivative to second order
     PtttStar = 0

     do iv=0,numberOfPolarizationVectors-1
       pvn(iv) = betav(iv)*( .5*dt*a1v(iv)*evn + rhspv(iv) )
       ! Now update PtttStar for correction terms and use in Ettt
       ptttStarv(iv) = b1tttv(iv)*(pvn(iv)-pvm(iv))/(2.*dt) + b0tttv(iv)*pv(iv) + a0tttv(iv)*ev + a1tttv(iv)*(evn-evm)/(2.*dt) + a2tttv(iv)*(evn-2*ev+evm)/dtsq + f2v(iv)
       PtttStar = PtttStar + ptttStarv(iv)


       write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pvm,pv,pvn,evm,ev,evn=",6e16.8)') i1,i2,m,pvm(iv),pv(iv),pvn(iv),evm,ev,evn      

       write(*,'(" (i2,i2,m)=(",i3,i3,i2,") evn-evm, evn-2*ev+evm=",2e16.8)') i1,i2,m,evn-evm,evn-2*ev+evm 

       write(*,'(" (i2,i2,m)=(",i3,i3,i2,") 1diff,2diff,3diff,f2v(iv)=",4e16.8)') i1,i2,m,(pvn(iv)-pvm(iv))/(2.*dt),(evn-evm)/(2.*dt),(evn-2*ev+evm)/dtsq,f2v(iv)
       write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pvn(iv),ptttStarv,dtsq=",3e16.8)') i1,i2,m,pvn(iv),ptttStarv(iv),dtsq
     end do

     ! Second Order Updates Complete, now we construct necessary terms
     ! LapPtt using prediction
     QxxStar = (Pxxn  - pxxSum)/dtsq
     EtxxStar    = (exxvn -  exxvm)/(2.*dt)
     EtttStar    = csq*EtxxStar - alphaP*PtttStar + f1

     write(*,'(" (i2,i2,m)=(",i3,i3,i2,") QxxStar,exxvn,Pxxn,EtxxStar,PtttStar,EtttStar=",7e16.8)') i1,i2,m,QxxStar,exxvn,Pxxn,EtxxStar,PtttStar,EtttStar,f1

     rhsP4 = 0
     LHSev = 0

     do iv=0,numberOfPolarizationVectors-1
         ! Coeff for pk(n+1) for  LHS of invidiual p equation
         LHSpv(iv) = 1+ b1v(iv)*dt/2.  + b0v(iv)*dtsq/12.

         ! Build coeff for ev for left hand side of full P equation
         LHSev     = LHSev + ((-a1v(iv)*dt/2.  -a0v(iv)*dtsq/12.)/LHSpv(iv))

         ! LHS for pk for pk equation
         rhspv(iv) = ((2.*pv(iv)-pvm(iv))\
         + b1v(iv)*dt*pvm(iv)/(2.)\
         - dtsq*b0v(iv)*pv(iv)\
         + dtsq*a0v(iv)*ev\
         - a1v(iv)*(dt/2.)*evm\
         - a1v(iv)*(dt**4/12.)*( EtttStar )\
         + b1v(iv)*(dt**4/12.)*( ptttStarv(iv) )\
         + (dtsq/12.)*( -b0v(iv)*(0. -2.*pv(iv)+pvm(iv)) + a0v(iv)*(0. -2.*ev+evm) )\
         + f3v(iv))

        write(*,'(" (i1,i2,m)=(",i3,i3,i2,") LHSpv,rhspv=",2e16.8)') i1,i2,m,LHSpv(iv),rhspv(iv)


        rhsP4 = rhsP4 + ( rhspv(iv) )/LHSpv(iv)

     end do

     ! We have now built the equation for P
     A(1,1) = LHSev
     A(1,2) = 1.

     b(1) = rhsP4

     ! Now we build the equation for E in terms of E and P

     A(2,1) = 1.        ! coeff of E^{n+1}
     A(2,2) = alphaP    ! coeff of P^{n+1}

     ! Note that Psum - rhsP here is same as for second order code
     !  (but using 4th order accurate version of p to compute them)
     b(2) = (2.*ev-evm)\
            +csq*dtsq*elap4\
            +alphaP*( Psum )\
            +csq**2*dt**4*elapsq2/(12.)\
            -csq*dt**4*alphaP*QxxStar/(12.)\
            + f4


     write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b2,f4=",2e16.8)') i1,i2,m,b(2),f4


     deti = 1./(A(1,1)*A(2,2) - A(1,2)*A(2,1))

     y(1) = deti *( A(2,2) * b(1) - A(1,2) * b(2))
     y(2) = deti *(-A(2,1) * b(1) + A(1,1) * b(2))

     write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b1,b2,y1,y2,deti=",5e16.8)') i1,i2,m,b(1),b(2),y(1),y(2),deti

     ! Update E^{n+1}
     un(i1,i2,i3,ec) = y(1)
     evn             = y(1)

     ! Update pk using new E^{n+1} = evn
     do iv=0,numberOfPolarizationVectors-1
       rhspv(iv) = rhspv(iv) + (a1v(iv) * dt/(2.) + a0v(iv)*dtsq/(12.))*evn
       pn(i1,i2,i3,m+iv*nd)  = (1/LHSpv(iv)) * rhspv(iv)
     end do



   end do ! m=0,1
  endLoopsMask()

 end if
#endMacro



! ===========================================================================================
! Macro:     DISPERSIVE: CURVILINEAR (NEW), 2D, ORDER 4
! ===========================================================================================
#beginMacro updateCurvilinear2dOrder4Dispersive()

 if( .true. .and. numberOfPolarizationVectors.eq.1 )then
  INFO("FD44r-dispersive");

  write(*,'(" DISPERSIVE: CURVILINEAR, 2D, ORDER 4 (ONE PV) ")')

  fp=0
  fe=0.

  ! Coefficients in equation for Pttt
  ! Pttt = b1ttt*P_t + b0ttt*P + a0ttt*E + a1ttt*E_t + a2ttt*Ett
  b1ttt=b1*b1-b0
  b0ttt=b1*b0
  a0ttt=-a0*b1
  a1ttt=a0-a1*b1
  a2ttt=a1

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! Advance Hz first:
    ! For now solve H_t = -(1/mu)*(  (E_y)_x - (E_x)_y )
    !   USE AB2 -- note: this is just a quadrature so stability is not an issue
    un(i1,i2,i3,hz) = u(i1,i2,i3,hz) -(dt/mu)*( 1.5*ux42r(i1,i2,i3,ey) -.5*umx42r(i1,i2,i3,ey) \
                                               -1.5*uy42r(i1,i2,i3,ex) +.5*umy42r(i1,i2,i3,ex) )

    addtForcingHz()

    do m=0,1
     pc=pxc+m
     ec=ex+m

     !if( addForcing.ne.0 )then
     !  fe = dtsq*f(i1,i2,i3,ec)
     !  getGDMForcing(ec,pc)
     !  fp=fpv(0)
     !end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     pv0 = p(i1,i2,i3,m)
     pvm0=pm(i1,i2,i3,m)

     !rhsE = maxwell2dr44me(i1,i2,i3,ec) + alphaP*(2.*pv0-pvm0) + fe
     !rhsP = 2.*pv0-pvm0 + .5*dt*( b1*pvm0 -a1*evm ) + dt*dt*( -b0*pv0 + a0*ev ) + fp

     !deti = 1./(1.+.5*dt*(b1+alphaP*a1))
     !un(i1,i2,i3,ec) = ((1.+.5*dt*b1)*rhsE -alphaP*rhsP)*deti
     !pn(i1,i2,i3,m)  = (.5*a1*dt*rhsE            + rhsP)*deti

     ! New 4th Order code starts here
     ! pvx = px22r(i1,i2,i3,m)


     f1 = 0
     f2 = 0
     f3 = 0
     f4 = 0
     f5 = 0
     f6 = 0


     if( addForcing.ne.0 )then
       OGDERIV2D( 0,0,0,0,i1,i2,i3,t, pc, p0    )
       OGDERIV2D( 1,0,0,0,i1,i2,i3,t, pc, p0t   )
       OGDERIV2D( 2,0,0,0,i1,i2,i3,t, pc, p0tt  )
       OGDERIV2D( 3,0,0,0,i1,i2,i3,t, pc, p0ttt )
       OGDERIV2D( 4,0,0,0,i1,i2,i3,t, pc, p0tttt)
       OGDERIV2D( 0,2,0,0,i1,i2,i3,t, pc, p0xx  )
       OGDERIV2D( 0,0,2,0,i1,i2,i3,t, pc, p0yy  )
       OGDERIV2D( 1,2,0,0,i1,i2,i3,t, pc, p0xxt )
       OGDERIV2D( 1,0,2,0,i1,i2,i3,t, pc, p0yyt )
       OGDERIV2D( 2,2,0,0,i1,i2,i3,t, pc, p0xxtt)
       OGDERIV2D( 2,0,2,0,i1,i2,i3,t, pc, p0yytt)


       OGDERIV2D( 0,0,0,0,i1,i2,i3,t, ec, e0    )
       OGDERIV2D( 1,0,0,0,i1,i2,i3,t, ec, e0t   )
       OGDERIV2D( 2,0,0,0,i1,i2,i3,t, ec, e0tt  )
       OGDERIV2D( 3,0,0,0,i1,i2,i3,t, ec, e0ttt )
       OGDERIV2D( 4,0,0,0,i1,i2,i3,t, ec, e0tttt)
       OGDERIV2D( 0,2,0,0,i1,i2,i3,t, ec, e0xx  )
       OGDERIV2D( 0,0,2,0,i1,i2,i3,t, ec, e0yy  )
       OGDERIV2D( 1,2,0,0,i1,i2,i3,t, ec, e0xxt )
       OGDERIV2D( 1,0,2,0,i1,i2,i3,t, ec, e0yyt )
       OGDERIV2D( 2,2,0,0,i1,i2,i3,t, ec, e0xxtt)
       OGDERIV2D( 2,0,2,0,i1,i2,i3,t, ec, e0yytt)
       OGDERIV2D( 0,4,0,0,i1,i2,i3,t, ec, e0xxxx)
       OGDERIV2D( 0,2,2,0,i1,i2,i3,t, ec, e0xxyy)
       OGDERIV2D( 0,0,4,0,i1,i2,i3,t, ec, e0yyyy)


       fp00  = p0tt   + b1*p0t   + b0*p0   - a1*e0t   - a0*e0
       fp10  = p0ttt  + b1*p0tt  + b0*p0t  - a1*e0tt  - a0*e0t
       fp20  = p0tttt + b1*p0ttt + b0*p0tt - a1*e0ttt - a0*e0tt
       fp02x = p0xxtt + b1*p0xxt + b0*p0xx - a1*e0xxt - a0*e0xx
       fp02y = p0yytt + b1*p0yyt + b0*p0yy - a1*e0yyt - a0*e0yy
       fp02  = fp02x  + fp02y

       fe00  = e0tt   + alphaP*p0tt   - csq*(e0xx   + e0yy)
       fe10  = e0ttt  + alphaP*p0ttt  - csq*(e0xxt  + e0yyt)
       fe20  = e0tttt + alphaP*p0tttt - csq*(e0xxtt + e0yytt)
       fe02x = e0xxtt + alphaP*p0xxtt - csq*(e0xxxx + e0xxyy)
       fe02y = e0yytt + alphaP*p0yytt - csq*(e0xxyy + e0yyyy)
       fe02  = fe02x  + fe02y

       ! f1 does not appear in this version of the 4th order scheme
       f1 = 0

       ! Forcing on P_ttt is (fp)_t - b1*fp :
       ! f2 = fp(1,0,x) - b1 * fp(0,0,x);
       ! f2 = p0t - b1 * p0 (incorrect)
       f2 = fp10 - b1*fp00

       ! Forcing on P eqn:
       ! fp + (dt^2/12) * fp_tt
       ! f3 = dt^2 * (fp(0,0,x) + (dt^2/12) * (- a1 * fe(1,0,x) + fp(2,0,x)));  % FIX ME *wdh* Jan 8, 2018
       ! f3 = dtsq * (p0 + (dtsq/12.) * (-a1 * e0t + p0tt )) incorrect
       f3 = dtsq * (fp00 + (dtsq/12.) * (-a1 * fe10 + fp20))

       ! Forcing on E equation is
       ! fe + (dt^2)/12 * ( fe_tt + c^2*Delta(fe) )
       ! f4 = dt^2 * (fe(0,0,x) + (dt^2/12) * (fe(2,0,x) + c^2 * fe(0,2,x)));
       ! f4 = dtsq * (e0 + (dtsq/12.) * (e0tt + csq * (e0xx + e0yy))) incorrect
       f4 = dtsq * (fe00 + (dtsq/12.) * (fe20 + csq * fe02))


       ! Forcing for Exx equation at second order predictor is (xx derivative
       ! of fe)
       ! f5 = dt^2 * fe(0,2,x);
       ! f5 = dtsq * (e0xx + e0yy)
       f5 = dtsq * fe02

       ! Forcing for Pxx equation at second order predictor is (xx derivative
       ! of fp)
       ! f6 = dt^2 * fp(0,2,x);
       ! f6 = dtsq * (p0xx + p0yy)
       f6 = dtsq * fp02

     end if

     elap4   = lap2d4(i1,i2,i3,ec)
     elap4m  = lap2d4m(i1,i2,i3,ec)
     elapsq2 = lap2d2Pow2(i1,i2,i3,ec)
     elap2m  = lap2d2m(i1,i2,i3,ec)
     plap2   = plap2d2(i1,i2,i3,m)
     plap2m  = plap2d2m(i1,i2,i3,m)
     plap4   = plap2d4(i1,i2,i3,m)
     plap4m  = plap2d4m(i1,i2,i3,m)
     ! Matlab code starts here

     elap4   = ulaplacian42(i1,i2,i3,ec)
     elap4m  = umlaplacian42(i1,i2,i3,ec)
     elapsq2 = vlaplacian22(i1,i2,i3,ec)
     elap2m  = umlaplacian22(i1,i2,i3,ec)
     plap2   = plaplacian22(i1,i2,i3,m)
     plap2m  = pmlaplacian22(i1,i2,i3,m)
     plap4   = plaplacian42(i1,i2,i3,m)
     plap4m  = pmlaplacian42(i1,i2,i3,m)

     !write(*,'(" (i1,i2,m)=(",i3,i3,i2,") elap4,elap4m,elapsq2,elap2m,plap2,plap2m,plap4,plap4m=",4e16.8)') i1,i2,m,elap4,elap4m,elapsq2,elap2m
     !write(*,'(" (i1,i2,m)=(",i3,i3,i2,") you are here  plap2,plap2m,plap4,plap4m=",4e16.8)') i1,i2,m,plap2,plap2m,plap4,plap4m

     Exxxx   = elapsq2

     En   = ev
     Enm1 = evm

     ! Pn   = pv0
     Pnm1 = pvm0

     ! compute c^2*Delta( E_t ) using
     ! Et = (1.5*En - 2*Enm1 +.5*Enm2)/dt;

     Exxn   = elap4
     Exxnm1 = elap4m


     betaP = 1/( 1+ (b1+alphaP*a1)*dt/2 )
     Pxxn   = plap4
     Pxxnm1 = plap4m



     ! Approximation for P_ttt
     ! Pttt^*    = b1ttt*(Pnp1 -Pnm1)/(2.*dt) + b0ttt*pv0 + a0ttt*En + a1ttt*(Enp1-Enm1)/(2.*dt) + a2ttt*(Enp1-2*En+Enm1)/dt^2; 
     cStar=1  !  cStar=1 : use Pttt^*, 0: use old

     PtttStarRhs = b1ttt*(0.   -Pnm1)/(2.*dt) + b0ttt*pv0 + a0ttt*En + a1ttt*(0.  -Enm1)/(2.*dt) + a2ttt*(0.  -2*En+Enm1)/dtsq + f2

     ! Use second order update to predict new E and P for future space derivatives
     !   fexx = fe;
     !   fpxx = fp;

     rhsExxnp1Predict = 2.*Exxn - Exxnm1 + cdtsq*Exxxx + alphaP*(2.*Pxxn - Pxxnm1) + f5
     rhsPxxnp1Predict = 2.*Pxxn - Pxxnm1 + .5*dt*( b1*Pxxnm1 - a1*Exxnm1 ) + dtsq*( -b0*Pxxn + a0*Exxn ) + f6
     deti = 1./(1.+.5*dt*(b1+alphaP*a1))
     Exxnp1Predict = ((1.+.5*dt*b1)*rhsExxnp1Predict - alphaP*rhsPxxnp1Predict)*deti
     Pxxnp1Predict = (.5*a1*dt*rhsExxnp1Predict            + rhsPxxnp1Predict)*deti

     QxxStar = (Pxxnp1Predict - 2.*Pxxn + Pxxnm1)/dtsq
     Etxx = (Exxnp1Predict - Exxnm1)/(2.*dt)



     ! A = zeros(2,2);
     ! b = zeros(2,1);

     !------------------------------- P EQUATION ----------------------------------------
     ! D+tD-t P^n = -b1*P_t -b0*P + a0*E + a1*E_t + (dt^2/12)*( D+tD-t Q^n )
     !            = -b1*[D0t(P^n) ]  - b0*P^n + a0*E^n  + a1*[ D0t(E^n) - (dt^2/12)*E_tttStar ] ...
     !              + (dt^2/12)*( -b0*D+tD-t(P^n) + a0*D+tD-t(E^n) )
     !              + (b1+alphaP*a1)*(dt^2)/12*P_ttt^*


     !
     ! We can either use:
     !       Qt = (Q^{n+1}-Q^{n-1})/(2*dt)   <-- adds a term to A(1,3)
     !       Qt = ( 1.5*Qn-2.*Qnm1 +.5*Qnm2 )/dt;  <-- remove term from A(1,3) 

     A(1,1) =   -a1*dt/2.  -a0*dtsq/12.   + (b1+alphaP*a1)*(dtsq/12.)*dtsq*( -a1ttt/(2.*dt) - a2ttt/dtsq ) 
     A(1,2) = 1+ b1*dt/2.  +b0*dtsq/12.   + (b1+alphaP*a1)*(dtsq/12.)*dtsq*( -b1ttt/(2.*dt) )   

     b(1) = (2.*pv0-Pnm1)\
         + b1*dt*Pnm1/(2.)\
         - dtsq*b0*pv0\
         + dtsq*a0*En\
         - a1*(dt/2.)*Enm1\
         - a1*(dt**4/12.)*( csq*Etxx )\
         + (dtsq/12.)*( -b0*(0. -2.*pv0+Pnm1) + a0*(0. -2.*En+Enm1) )\
         + (b1+alphaP*a1)*(dtsq/12.)*dtsq*( PtttStarRhs )\
         + f3

      ! ------------------------------- E EQUATION ----------------------------------------
      ! D+tD-t E^n = c^2*Delta(E^n) -alphaP*D+tD-t P^n + (dt^2/(12)*( [c^2*Delta]^2 E^n  -alphaP*c^2*Qxx^*  )

     A(2,1) = 1.        ! coeff of E^{n+1}
     A(2,2) = alphaP    ! coeff of P^{n+1}

     b(2) = (2.*En-Enm1)\
            +csq*dtsq*Exxn\
            +(2.*pv0-Pnm1)*alphaP\
            +csq**2*dt**4*Exxxx/(12.)\
            -csq*dt**4*alphaP*QxxStar/(12.)\
            + f4

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") Exxn,Exxxx,QxxStar=",6e16.8)') i1,i2,m,Exxn,Exxxx,QxxStar

     deti = 1./(A(1,1)*A(2,2) - A(1,2)*A(2,1))

     y(1) = deti *( A(2,2) * b(1) - A(1,2) * b(2))
     y(2) = deti *(-A(2,1) * b(1) + A(1,1) * b(2))

     ! write(*,'(" A(1,1),A(1,2),A(2,1),A(2,2),b(1),b(2)=",6e16.8)') A(1,1),A(1,2),A(2,1),A(2,2),b(1),b(2)


     !u(j,1) = y(1);
     !u(j,2) = y(2);


     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f4,f5,f6=",5e16.8)') i1,i2,m,f2,f3,f4,f5,f6

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") e2,e4,p2,p4=",4e16.8)') i1, i2, m, un(i1,i2,i3,ec), y(1), pn(i1,i2,i3,m), y(2)

     ! write(*,'(" En,Enm1,Exxn,pv0,Pnm1,alphaP,csq=",7e16.8)') En,Enm1,Exxn,pv0,Pnm1,alphaP,csq

     ! write(*,'(" dt,Exxxx,QxxStar=",7e16.8)') dt,Exxxx,QxxStar

     ! write(*,'(" Pxxnp1Predict,Pxxn,Pxxnm1=",7e16.8)') Pxxnp1Predict,Pxxn,Pxxnm1

     ! write(*,'(" rhsExxnp1Predict, rhsPxxnp1Predict, f5, f6=",4e16.8)') rhsExxnp1Predict, rhsPxxnp1Predict, f5, f6

     ! write(*,'(" Exxn,Exxnm1,cdtsq,Exxxx,Pxxn,Pxxnm1=",6e16.8)') Exxn,Exxnm1,cdtsq,Exxxx,Pxxn,Pxxnm1

     un(i1,i2,i3,ec) = y(1)
     pn(i1,i2,i3,m)  = y(2)

     end do
  endLoopsMask()

 else

  ! ------- 2D DISPERSIVE CURVILINEAR MULTIPLE PV (ORDER 4)-------
  ! XYZ2

  INFO("FD44r-dispersive-MULTI-PV");
  write(*,'(" DISPERSIVE: CURVILINEAR, 2D, ORDER 4, MULTIPLE PV:FINISH ME (YOU ARE HERE)  ")')

  fe=0.
  ! -- first compute some coefficients ---
  beta=0.
  do iv=0,numberOfPolarizationVectors-1
    betav(iv) = 1./( 1.+.5*dt*b1v(iv) )
    beta = beta + .5*dt*a1v(iv)*betav(iv)
    fpv(iv)=0.  ! initialize if not used
    b1tttv(iv)=b1v(iv)*b1v(iv)-b0v(iv)
    b0tttv(iv)=b1v(iv)*b0v(iv)
    a0tttv(iv)=-a0v(iv)*b1v(iv)
    a1tttv(iv)=a0v(iv)-a1v(iv)*b1v(iv)
    a2tttv(iv)=a1v(iv)
  end do

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! Advance Hz first:
    ! For now solve H_t = -(1/mu)*(  (E_y)_x - (E_x)_y )
    !   USE AB2 -- note: this is just a quadrature so stability is not an issue
    un(i1,i2,i3,hz) = u(i1,i2,i3,hz) -(dt/mu)*( 1.5*ux42r(i1,i2,i3,ey) -.5*umx42r(i1,i2,i3,ey) \
                                               -1.5*uy42r(i1,i2,i3,ex) +.5*umy42r(i1,i2,i3,ex) )

    addtForcingHz()

    ! -- loop over components of the vector --
    do m=0,1
     pc=pxc+m
     ec=ex+m

     f1 = 0
     f4 = 0
     f5 = 0

     do iv=0,numberOfPolarizationVectors-1
       f2v(iv) = 0
       f3v(iv) = 0
       f6v(iv) = 0
       ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f6=",3e16.8)') i1,i2,m,f2v(iv),f3v(iv),f6v(iv)
     end do

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f1,f4,f5=",3e16.8)') i1,i2,m,f1,f4,f5

     if( addForcing.ne.0 )then
       ! 2nd order forcing used in 2nd order prediction step
       ! fe = dtsq*f(i1,i2,i3,ec)
       ! Compute fpv(iv) :
       ! getGDMForcing(ec,pc)

       ! 4th order forcing used in 4th order scheme

       OGDERIV2D( 0,0,0,0,i1,i2,i3,t, ec, e0    )
       OGDERIV2D( 1,0,0,0,i1,i2,i3,t, ec, e0t   )
       OGDERIV2D( 2,0,0,0,i1,i2,i3,t, ec, e0tt  )
       OGDERIV2D( 3,0,0,0,i1,i2,i3,t, ec, e0ttt )
       OGDERIV2D( 4,0,0,0,i1,i2,i3,t, ec, e0tttt)
       OGDERIV2D( 0,2,0,0,i1,i2,i3,t, ec, e0xx  )
       OGDERIV2D( 0,0,2,0,i1,i2,i3,t, ec, e0yy  )
       OGDERIV2D( 1,2,0,0,i1,i2,i3,t, ec, e0xxt )
       OGDERIV2D( 1,0,2,0,i1,i2,i3,t, ec, e0yyt )
       OGDERIV2D( 2,2,0,0,i1,i2,i3,t, ec, e0xxtt)
       OGDERIV2D( 2,0,2,0,i1,i2,i3,t, ec, e0yytt)
       OGDERIV2D( 0,4,0,0,i1,i2,i3,t, ec, e0xxxx)
       OGDERIV2D( 0,2,2,0,i1,i2,i3,t, ec, e0xxyy)
       OGDERIV2D( 0,0,4,0,i1,i2,i3,t, ec, e0yyyy)

       pSum0tt    = 0
       pSum0ttt   = 0
       pSum0tttt  = 0
       pSum0xxtt  = 0
       pSum0yytt  = 0

       do iv=0,numberOfPolarizationVectors-1
         pce = pc+iv*nd

         OGDERIV2D( 0,0,0,0,i1,i2,i3,t, pce, p0    )
         OGDERIV2D( 1,0,0,0,i1,i2,i3,t, pce, p0t   )
         OGDERIV2D( 2,0,0,0,i1,i2,i3,t, pce, p0tt  )
         OGDERIV2D( 3,0,0,0,i1,i2,i3,t, pce, p0ttt )
         OGDERIV2D( 4,0,0,0,i1,i2,i3,t, pce, p0tttt)
         OGDERIV2D( 0,2,0,0,i1,i2,i3,t, pce, p0xx  )
         OGDERIV2D( 0,0,2,0,i1,i2,i3,t, pce, p0yy  )
         OGDERIV2D( 1,2,0,0,i1,i2,i3,t, pce, p0xxt )
         OGDERIV2D( 1,0,2,0,i1,i2,i3,t, pce, p0yyt )
         OGDERIV2D( 2,2,0,0,i1,i2,i3,t, pce, p0xxtt)
         OGDERIV2D( 2,0,2,0,i1,i2,i3,t, pce, p0yytt)

         ! Derivatives of OG individual p eqn forcing terms fp           
         fp00      = p0tt   + b1v(iv)*p0t   + b0v(iv)*p0   - a1v(iv)*e0t   - a0v(iv)*e0
         fp10      = p0ttt  + b1v(iv)*p0tt  + b0v(iv)*p0t  - a1v(iv)*e0tt  - a0v(iv)*e0t
         fp20      = p0tttt + b1v(iv)*p0ttt + b0v(iv)*p0tt - a1v(iv)*e0ttt - a0v(iv)*e0tt
         fp02x     = p0xxtt + b1v(iv)*p0xxt + b0v(iv)*p0xx - a1v(iv)*e0xxt - a0v(iv)*e0xx
         fp02y     = p0yytt + b1v(iv)*p0yyt + b0v(iv)*p0yy - a1v(iv)*e0yyt - a0v(iv)*e0yy
         fp02      = fp02x  + fp02y
         fp00v(iv) = fp00

         ! Building derivatives of full P summation terms
         pSum0tt   = pSum0tt    + p0tt
         pSum0ttt  = pSum0ttt   + p0ttt
         pSum0tttt = pSum0tttt  + p0tttt
         pSum0xxtt = pSum0xxtt  + p0xxtt
         pSum0yytt = pSum0yytt  + p0yytt

         ! Forcing on EACH individual p_ttt is (fp)_t - b1*fp :
         f2v(iv) = fp10 - b1v(iv)*fp00

         ! Forcing on EACH individual p eqn:
         f3v(iv) = dtsq * (fp00 + (dtsq/12.) * fp20)

         ! Forcing for EACH individual pxx equation at second order predictor is (xx derivative of fp)
         f6v(iv) = dtsq * fp02

         ! if( iv.eq.0 )then
         !   write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f6=",3e16.8)') i1,i2,m,f2v(iv),f3v(iv),f6v(iv)
         ! end if
       end do

       ! Build derivatives of E eqn forcing term
       fe00  = e0tt   + alphaP*pSum0tt   - csq*(e0xx   + e0yy)
       fe10  = e0ttt  + alphaP*pSum0ttt  - csq*(e0xxt  + e0yyt)
       fe20  = e0tttt + alphaP*pSum0tttt - csq*(e0xxtt + e0yytt)
       fe02x = e0xxtt + alphaP*pSum0xxtt - csq*(e0xxxx + e0xxyy)
       fe02y = e0yytt + alphaP*pSum0yytt - csq*(e0xxyy + e0yyyy)
       fe02  = fe02x  + fe02y
       fe    = fe00

       ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") e0xxxx,e0xxyy=",2e16.8)') i1,i2,m,e0xxxx,e0xxyy

       ! Forcing on Ettt = c^2 Etxx - alphaP Pttt is fet
       f1 = fe10

       ! Forcing on E equation is
       f4 = dtsq * (fe00 + (dtsq/12.) * (fe20 + csq * fe02))

       ! Forcing for Exx equation at second order predictor is (xx derivative of fe)
       f5 = dtsq * fe02

       write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f1,f4,f5=",3e16.8)') i1,i2,m,f1,f4,f5
     end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     rhsP = 0.
     pSum=0.

     ! Fourth order below

     elap4   = lap2d4(i1,i2,i3,ec)
     elap4m  = lap2d4m(i1,i2,i3,ec)
     elapsq2 = lap2d2Pow2(i1,i2,i3,ec)
     elap2m  = lap2d2m(i1,i2,i3,ec)


     elap4   = ulaplacian42(i1,i2,i3,ec)
     elap4m  = umlaplacian42(i1,i2,i3,ec)
     elapsq2 = vlaplacian22(i1,i2,i3,ec)
     elap2m  = umlaplacian22(i1,i2,i3,ec)

     rhsPxx = 0.
     pxxSum = 0.

     exxv  = elap4
     exxvm = elap4m

     ! First we do the second order prediction on (i) E, (ii) Exx, (iii) Pxx, (iv) Individual pk
     do iv=0,numberOfPolarizationVectors-1

       ! We use these in all versions
       pv(iv)    = p(i1,i2,i3,m+iv*nd)
       pvm(iv)   = pm(i1,i2,i3,m+iv*nd)

       ! We used these terms for the rectangular code
       ! pxxv(iv)  = plap2d4(i1,i2,i3,m+iv*nd)
       ! pxxvm(iv) = plap2d4m(i1,i2,i3,m+iv*nd)

       ! We used these for the single PV curvilinear code (probably not 2nd order ones)
       ! plap2   = plaplacian22(i1,i2,i3,m)
       ! plap2m  = pmlaplacian22(i1,i2,i3,m)
       ! plap4   = plaplacian42(i1,i2,i3,m)
       ! plap4m  = pmlaplacian42(i1,i2,i3,m)

       ! We use these for the multiple PV curvilinear code
       pxxv(iv)   = plaplacian42(i1,i2,i3,m+iv*nd)
       pxxvm(iv)  = pmlaplacian42(i1,i2,i3,m+iv*nd)

       ! rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + fp00v(iv)
       rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + dtsq*fp00v(iv)
       rhsP = rhsP + betav(iv)*rhspv(iv)
       pSum = pSum + 2.*pv(iv) - pvm(iv)


       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pv,pvm,rhspv,rhsP,pSum,fpv=",6e16.8)') i1,i2,m,pv(iv),pvm(iv),rhspv(iv),rhsP,pSum,fp00v(iv)

       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") ev,evm,rhspv,rhsP,pSum,fpv=",6e16.8)') i1,i2,m,ev,evm,a1v(iv),a0v(iv),b1v(iv),b0v(iv)


       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") first term=",1e16.8)') i1,i2,m,2.*pv(iv)-pvm(iv)

       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") two terms",2e16.8)') i1,i2,m,.5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ),dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev )


       rhspxxv(iv) = 2.*pxxv(iv)-pxxvm(iv) + .5*dt*( b1v(iv)*pxxvm(iv) -a1v(iv)*exxvm ) + dtSq*( -b0v(iv)*pxxv(iv) + a0v(iv)*exxv ) + f6v(iv)
       rhsPxx = rhsPxx + betav(iv)*rhspxxv(iv)
       pxxSum = pxxSum + 2.*pxxv(iv) - pxxvm(iv)




     end do

     ! rhsE   = maxwell2dr44me(i1,i2,i3,ec)     + alphaP*( pSum   - rhsP   ) + dtsq * fe00 
     rhsE   = 2.*ev   - evm   + cdtsq*elap4   + alphaP*( pSum   - rhsP   ) + dtsq * fe00
     rhsExx = 2.*exxv - exxvm + cdtsq*elapsq2 + alphaP*( pxxSum - rhsPxx ) + f5



     evn   = rhsE   / (1.+ alphaP*beta)
     exxvn = rhsExx / (1.+ alphaP*beta)



     write(*,'(" (i2,i2,m)=(",i3,i3,i2,") rhsE,evn,fe00=",3e16.8)') i1,i2,m,rhsE,evn,fe00


     ! Update x derivative of P
     Pxxn  = beta * exxvn + rhsPxx

     ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") exxvn,Pxxn=",2e16.8)') i1,i2,m,exxvn,Pxxn


     ! Now we predict individual pk (pvn) and their third time derivative to second order
     PtttStar = 0

     do iv=0,numberOfPolarizationVectors-1
       pvn(iv) = betav(iv)*( .5*dt*a1v(iv)*evn + rhspv(iv) )
       ! Now update PtttStar for correction terms and use in Ettt
       ptttStarv(iv) = b1tttv(iv)*(pvn(iv)-pvm(iv))/(2.*dt) + b0tttv(iv)*pv(iv) + a0tttv(iv)*ev + a1tttv(iv)*(evn-evm)/(2.*dt) + a2tttv(iv)*(evn-2*ev+evm)/dtsq + f2v(iv)
       PtttStar = PtttStar + ptttStarv(iv)


       write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pvm,pv,pvn,evm,ev,evn=",6e16.8)') i1,i2,m,pvm(iv),pv(iv),pvn(iv),evm,ev,evn      

       write(*,'(" (i2,i2,m)=(",i3,i3,i2,") evn-evm, evn-2*ev+evm=",2e16.8)') i1,i2,m,evn-evm,evn-2*ev+evm 

       write(*,'(" (i2,i2,m)=(",i3,i3,i2,") 1diff,2diff,3diff,f2v(iv)=",4e16.8)') i1,i2,m,(pvn(iv)-pvm(iv))/(2.*dt),(evn-evm)/(2.*dt),(evn-2*ev+evm)/dtsq,f2v(iv)
       write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pvn(iv),ptttStarv,dtsq=",3e16.8)') i1,i2,m,pvn(iv),ptttStarv(iv),dtsq
     end do

     ! Second Order Updates Complete, now we construct necessary terms
     ! LapPtt using prediction
     QxxStar = (Pxxn  - pxxSum)/dtsq
     EtxxStar    = (exxvn -  exxvm)/(2.*dt)
     EtttStar    = csq*EtxxStar - alphaP*PtttStar + f1

     write(*,'(" (i2,i2,m)=(",i3,i3,i2,") QxxStar,exxvn,Pxxn,EtxxStar,PtttStar,EtttStar=",7e16.8)') i1,i2,m,QxxStar,exxvn,Pxxn,EtxxStar,PtttStar,EtttStar,f1

     rhsP4 = 0
     LHSev = 0

     do iv=0,numberOfPolarizationVectors-1
         ! Coeff for pk(n+1) for  LHS of invidiual p equation
         LHSpv(iv) = 1+ b1v(iv)*dt/2.  + b0v(iv)*dtsq/12.

         ! Build coeff for ev for left hand side of full P equation
         LHSev     = LHSev + ((-a1v(iv)*dt/2.  -a0v(iv)*dtsq/12.)/LHSpv(iv))

         ! LHS for pk for pk equation
         rhspv(iv) = ((2.*pv(iv)-pvm(iv))\
         + b1v(iv)*dt*pvm(iv)/(2.)\
         - dtsq*b0v(iv)*pv(iv)\
         + dtsq*a0v(iv)*ev\
         - a1v(iv)*(dt/2.)*evm\
         - a1v(iv)*(dt**4/12.)*( EtttStar )\
         + b1v(iv)*(dt**4/12.)*( ptttStarv(iv) )\
         + (dtsq/12.)*( -b0v(iv)*(0. -2.*pv(iv)+pvm(iv)) + a0v(iv)*(0. -2.*ev+evm) )\
         + f3v(iv))

        write(*,'(" (i1,i2,m)=(",i3,i3,i2,") LHSpv,rhspv=",2e16.8)') i1,i2,m,LHSpv(iv),rhspv(iv)


        rhsP4 = rhsP4 + ( rhspv(iv) )/LHSpv(iv)

     end do

     ! We have now built the equation for P
     A(1,1) = LHSev
     A(1,2) = 1.

     b(1) = rhsP4

     ! Now we build the equation for E in terms of E and P

     A(2,1) = 1.        ! coeff of E^{n+1}
     A(2,2) = alphaP    ! coeff of P^{n+1}

     ! Note that Psum here is same as for second order code
     !  (but using 4th order accurate version of p to compute them)
     b(2) = (2.*ev-evm)\
            +csq*dtsq*elap4\
            +alphaP*( Psum )\
            +csq**2*dt**4*elapsq2/(12.)\
            -csq*dt**4*alphaP*QxxStar/(12.)\
            + f4


     write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b2,f4=",2e16.8)') i1,i2,m,b(2),f4


     deti = 1./(A(1,1)*A(2,2) - A(1,2)*A(2,1))

     y(1) = deti *( A(2,2) * b(1) - A(1,2) * b(2))
     y(2) = deti *(-A(2,1) * b(1) + A(1,1) * b(2))

     write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b1,b2,y1,y2,deti=",5e16.8)') i1,i2,m,b(1),b(2),y(1),y(2),deti

     ! Update E^{n+1}
     un(i1,i2,i3,ec) = y(1)
     evn             = y(1)

     ! Update pk using new E^{n+1} = evn
     do iv=0,numberOfPolarizationVectors-1
       rhspv(iv) = rhspv(iv) + (a1v(iv) * dt/(2.) + a0v(iv)*dtsq/(12.))*evn
       pn(i1,i2,i3,m+iv*nd)  = (1/LHSpv(iv)) * rhspv(iv)
     end do



   end do ! m=0,1
  endLoopsMask()

 end if
#endMacro

! XYZ3


! ===========================================================================================
! Macro:     DISPERSIVE: RECTANGULAR, 3D, ORDER 4
! ===========================================================================================
#beginMacro updateRectangular3dOrder4Dispersive()

 if( .true. .and. numberOfPolarizationVectors.eq.1 )then
  INFO("FD44r-dispersive");

  write(*,'(" DISPERSIVE: RECTANGULAR, 3D, ORDER 4 (ONE PV) ")')

  fp=0
  fe=0.

  ! Coefficients in equation for Pttt
  ! Pttt = b1ttt*P_t + b0ttt*P + a0ttt*E + a1ttt*E_t + a2ttt*Ett
  b1ttt=b1*b1-b0
  b0ttt=b1*b0
  a0ttt=-a0*b1
  a1ttt=a0-a1*b1
  a2ttt=a1

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! Advance Hz first:
    ! For now solve H_t = -(1/mu)*(  (E_y)_x - (E_x)_y )
    !   USE AB2 -- note: this is just a quadrature so stability is not an issue
    un(i1,i2,i3,hz) = u(i1,i2,i3,hz) -(dt/mu)*( 1.5*ux42r(i1,i2,i3,ey) -.5*umx42r(i1,i2,i3,ey) \
                                               -1.5*uy42r(i1,i2,i3,ex) +.5*umy42r(i1,i2,i3,ex) )

    addtForcingHz()

    do m=0,nd-1
     pc=pxc+m
     ec=ex+m

     !if( addForcing.ne.0 )then
     !  fe = dtsq*f(i1,i2,i3,ec)
     !  getGDMForcing(ec,pc)
     !  fp=fpv(0)
     !end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     pv0 = p(i1,i2,i3,m)
     pvm0=pm(i1,i2,i3,m)

     !rhsE = maxwell3dr44me(i1,i2,i3,ec) + alphaP*(2.*pv0-pvm0) + fe
     !rhsP = 2.*pv0-pvm0 + .5*dt*( b1*pvm0 -a1*evm ) + dt*dt*( -b0*pv0 + a0*ev ) + fp

     !deti = 1./(1.+.5*dt*(b1+alphaP*a1))
     !un(i1,i2,i3,ec) = ((1.+.5*dt*b1)*rhsE -alphaP*rhsP)*deti
     !pn(i1,i2,i3,m)  = (.5*a1*dt*rhsE            + rhsP)*deti

     ! New 4th Order code starts here
     ! pvx = px22r(i1,i2,i3,m)


     f1 = 0
     f2 = 0
     f3 = 0
     f4 = 0
     f5 = 0
     f6 = 0


     if( addForcing.ne.0 )then
       OGDERIV3D( 0,0,0,0,i1,i2,i3,t, pc, p0    )
       OGDERIV3D( 1,0,0,0,i1,i2,i3,t, pc, p0t   )
       OGDERIV3D( 2,0,0,0,i1,i2,i3,t, pc, p0tt  )
       OGDERIV3D( 3,0,0,0,i1,i2,i3,t, pc, p0ttt )
       OGDERIV3D( 4,0,0,0,i1,i2,i3,t, pc, p0tttt)
       OGDERIV3D( 0,2,0,0,i1,i2,i3,t, pc, p0xx  )
       OGDERIV3D( 0,0,2,0,i1,i2,i3,t, pc, p0yy  )
       OGDERIV3D( 0,0,0,2,i1,i2,i3,t, pc, p0zz  )
       OGDERIV3D( 1,2,0,0,i1,i2,i3,t, pc, p0xxt )
       OGDERIV3D( 1,0,2,0,i1,i2,i3,t, pc, p0yyt )
       OGDERIV3D( 1,0,0,2,i1,i2,i3,t, pc, p0zzt )
       OGDERIV3D( 2,2,0,0,i1,i2,i3,t, pc, p0xxtt)
       OGDERIV3D( 2,0,2,0,i1,i2,i3,t, pc, p0yytt)
       OGDERIV3D( 2,0,0,2,i1,i2,i3,t, pc, p0zztt)



       OGDERIV3D( 0,0,0,0,i1,i2,i3,t, ec, e0    )
       OGDERIV3D( 1,0,0,0,i1,i2,i3,t, ec, e0t   )
       OGDERIV3D( 2,0,0,0,i1,i2,i3,t, ec, e0tt  )
       OGDERIV3D( 3,0,0,0,i1,i2,i3,t, ec, e0ttt )
       OGDERIV3D( 4,0,0,0,i1,i2,i3,t, ec, e0tttt)
       OGDERIV3D( 0,2,0,0,i1,i2,i3,t, ec, e0xx  )
       OGDERIV3D( 0,0,2,0,i1,i2,i3,t, ec, e0yy  )
       OGDERIV3D( 0,0,0,2,i1,i2,i3,t, ec, e0zz  )
       OGDERIV3D( 1,2,0,0,i1,i2,i3,t, ec, e0xxt )
       OGDERIV3D( 1,0,2,0,i1,i2,i3,t, ec, e0yyt )
       OGDERIV3D( 1,0,0,2,i1,i2,i3,t, ec, e0zzt )
       OGDERIV3D( 2,2,0,0,i1,i2,i3,t, ec, e0xxtt)
       OGDERIV3D( 2,0,2,0,i1,i2,i3,t, ec, e0yytt)
       OGDERIV3D( 2,0,0,2,i1,i2,i3,t, ec, e0zztt)
       OGDERIV3D( 0,4,0,0,i1,i2,i3,t, ec, e0xxxx)
       OGDERIV3D( 0,2,2,0,i1,i2,i3,t, ec, e0xxyy)
       OGDERIV3D( 0,2,0,2,i1,i2,i3,t, ec, e0xxzz)
       OGDERIV3D( 0,0,4,0,i1,i2,i3,t, ec, e0yyyy)
       OGDERIV3D( 0,0,2,2,i1,i2,i3,t, ec, e0yyzz)
       OGDERIV3D( 0,0,0,4,i1,i2,i3,t, ec, e0zzzz)


       fp00  = p0tt   + b1*p0t   + b0*p0   - a1*e0t   - a0*e0
       fp10  = p0ttt  + b1*p0tt  + b0*p0t  - a1*e0tt  - a0*e0t
       fp20  = p0tttt + b1*p0ttt + b0*p0tt - a1*e0ttt - a0*e0tt
       fp02x = p0xxtt + b1*p0xxt + b0*p0xx - a1*e0xxt - a0*e0xx
       fp02y = p0yytt + b1*p0yyt + b0*p0yy - a1*e0yyt - a0*e0yy
       fp02z = p0yytt + b1*p0zzt + b0*p0zz - a1*e0zzt - a0*e0zz
       fp02  = fp02x  + fp02y + fp02z

       fe00  = e0tt   + alphaP*p0tt   - csq*(e0xx   + e0yy   + e0zz  )
       fe10  = e0ttt  + alphaP*p0ttt  - csq*(e0xxt  + e0yyt  + e0zzt )
       fe20  = e0tttt + alphaP*p0tttt - csq*(e0xxtt + e0yytt + e0zztt)
       fe02x = e0xxtt + alphaP*p0xxtt - csq*(e0xxxx + e0xxyy + e0xxzz)
       fe02y = e0yytt + alphaP*p0yytt - csq*(e0xxyy + e0yyyy + e0yyzz)
       fe02z = e0zztt + alphaP*p0zztt - csq*(e0zzzz + e0xxzz + e0yyzz)
       fe02  = fe02x  + fe02y + fe02z

       ! f1 does not appear in this version of the 4th order scheme
       f1 = 0

       ! Forcing on P_ttt is (fp)_t - b1*fp :
       ! f2 = fp(1,0,x) - b1 * fp(0,0,x);
       ! f2 = p0t - b1 * p0 (incorrect)
       f2 = fp10 - b1*fp00

       ! Forcing on P eqn:
       ! fp + (dt^2/12) * fp_tt
       ! f3 = dt^2 * (fp(0,0,x) + (dt^2/12) * (- a1 * fe(1,0,x) + fp(2,0,x)));  % FIX ME *wdh* Jan 8, 2018
       ! f3 = dtsq * (p0 + (dtsq/12.) * (-a1 * e0t + p0tt )) incorrect
       f3 = dtsq * (fp00 + (dtsq/12.) * (-a1 * fe10 + fp20))

       ! Forcing on E equation is
       ! fe + (dt^2)/12 * ( fe_tt + c^2*Delta(fe) )
       ! f4 = dt^2 * (fe(0,0,x) + (dt^2/12) * (fe(2,0,x) + c^2 * fe(0,2,x)));
       ! f4 = dtsq * (e0 + (dtsq/12.) * (e0tt + csq * (e0xx + e0yy))) incorrect
       f4 = dtsq * (fe00 + (dtsq/12.) * (fe20 + csq * fe02))


       ! Forcing for Exx equation at second order predictor is (xx derivative
       ! of fe)
       ! f5 = dt^2 * fe(0,2,x);
       ! f5 = dtsq * (e0xx + e0yy)
       f5 = dtsq * fe02

       ! Forcing for Pxx equation at second order predictor is (xx derivative
       ! of fp)
       ! f6 = dt^2 * fp(0,2,x);
       ! f6 = dtsq * (p0xx + p0yy)
       f6 = dtsq * fp02

     end if

     ! YOU ARE HERE 3D

     elap4   = lap3d4(i1,i2,i3,ec)
     elap4m  = lap3d4m(i1,i2,i3,ec)
     elapsq2 = lap3d2Pow2(i1,i2,i3,ec)
     elap2m  = lap3d2m(i1,i2,i3,ec)
     plap2   = plap3d2(i1,i2,i3,m)
     plap2m  = plap3d2m(i1,i2,i3,m)
     plap4   = plap3d4(i1,i2,i3,m)
     plap4m  = plap3d4m(i1,i2,i3,m)
     ! Matlab code starts here

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f6=",3e16.8)') i1,i2,m,f2,f3,f6
     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f1,f4,f5=",3e16.8)') i1,i2,m,f1,f4,f5

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") elap4,elap4m,elapsq2,elap2m,plap2,plap2m,plap4,plap4m=",4e16.8)') i1,i2,m,elap4,elap4m,elapsq2,elap2m
     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") plap2,plap2m,plap4,plap4m=",4e16.8)') i1,i2,m,plap2,plap2m,plap4,plap4m

     Exxxx   = elapsq2

     En   = ev
     Enm1 = evm

     ! Pn   = pv0
     Pnm1 = pvm0

     ! compute c^2*Delta( E_t ) using
     ! Et = (1.5*En - 2*Enm1 +.5*Enm2)/dt;

     Exxn   = elap4
     Exxnm1 = elap4m


     betaP = 1/( 1+ (b1+alphaP*a1)*dt/2 )
     Pxxn   = plap4
     Pxxnm1 = plap4m



     ! Approximation for P_ttt
     ! Pttt^*    = b1ttt*(Pnp1 -Pnm1)/(2.*dt) + b0ttt*pv0 + a0ttt*En + a1ttt*(Enp1-Enm1)/(2.*dt) + a2ttt*(Enp1-2*En+Enm1)/dt^2;
     cStar=1  !  cStar=1 : use Pttt^*, 0: use old

     PtttStarRhs = b1ttt*(0.   -Pnm1)/(2.*dt) + b0ttt*pv0 + a0ttt*En + a1ttt*(0.  -Enm1)/(2.*dt) + a2ttt*(0.  -2*En+Enm1)/dtsq + f2

     ! Use second order update to predict new E and P for future space derivatives
     !   fexx = fe;
     !   fpxx = fp;

     ! These are actually Laplacians
     rhsExxnp1Predict = 2.*Exxn - Exxnm1 + cdtsq*Exxxx + alphaP*(2.*Pxxn - Pxxnm1) + f5
     rhsPxxnp1Predict = 2.*Pxxn - Pxxnm1 + .5*dt*( b1*Pxxnm1 - a1*Exxnm1 ) + dtsq*( -b0*Pxxn + a0*Exxn ) + f6
     deti = 1./(1.+.5*dt*(b1+alphaP*a1))
     Exxnp1Predict = ((1.+.5*dt*b1)*rhsExxnp1Predict - alphaP*rhsPxxnp1Predict)*deti
     Pxxnp1Predict = (.5*a1*dt*rhsExxnp1Predict            + rhsPxxnp1Predict)*deti


     write(*,'(" (i2,i2,m)=(",i3,i3,i2,") Exxnp1Predict,Pxxnp1Predict=",2e16.8)') i1,i2,m,Exxnp1Predict,Pxxnp1Predict





     QxxStar = (Pxxnp1Predict - 2.*Pxxn + Pxxnm1)/dtsq
     Etxx = (Exxnp1Predict - Exxnm1)/(2.*dt)

     write(*,'(" (i1,i2,m)=(",i3,i3,i2,") QxxStar,Exxn,Pxxn=",3e16.8)') i1,i2,m,QxxStar,Exxnp1Predict,Pxxnp1Predict

     ! A = zeros(2,2);
     ! b = zeros(2,1);

     !------------------------------- P EQUATION ----------------------------------------
     ! D+tD-t P^n = -b1*P_t -b0*P + a0*E + a1*E_t + (dt^2/12)*( D+tD-t Q^n )
     !            = -b1*[D0t(P^n) ]  - b0*P^n + a0*E^n  + a1*[ D0t(E^n) - (dt^2/12)*E_tttStar ] ...
     !              + (dt^2/12)*( -b0*D+tD-t(P^n) + a0*D+tD-t(E^n) )
     !              + (b1+alphaP*a1)*(dt^2)/12*P_ttt^*


     !
     ! We can either use:
     !       Qt = (Q^{n+1}-Q^{n-1})/(2*dt)   <-- adds a term to A(1,3)
     !       Qt = ( 1.5*Qn-2.*Qnm1 +.5*Qnm2 )/dt;  <-- remove term from A(1,3)

     A(1,1) =   -a1*dt/2.  -a0*dtsq/12.   + (b1+alphaP*a1)*(dtsq/12.)*dtsq*( -a1ttt/(2.*dt) - a2ttt/dtsq )
     A(1,2) = 1+ b1*dt/2.  +b0*dtsq/12.   + (b1+alphaP*a1)*(dtsq/12.)*dtsq*( -b1ttt/(2.*dt) )

     b(1) = (2.*pv0-Pnm1)\
         + b1*dt*Pnm1/(2.)\
         - dtsq*b0*pv0\
         + dtsq*a0*En\
         - a1*(dt/2.)*Enm1\
         - a1*(dt**4/12.)*( csq*Etxx )\
         + (dtsq/12.)*( -b0*(0. -2.*pv0+Pnm1) + a0*(0. -2.*En+Enm1) )\
         + (b1+alphaP*a1)*(dtsq/12.)*dtsq*( PtttStarRhs )\
         + f3


     !write(*,'(" (i1,i2,m)=(",i3,i3,i2,") rhspv=",1e16.8)') i1,i2,m,b(1)




      ! ------------------------------- E EQUATION ----------------------------------------
      ! D+tD-t E^n = c^2*Delta(E^n) -alphaP*D+tD-t P^n + (dt^2/(12)*( [c^2*Delta]^2 E^n  -alphaP*c^2*Qxx^*  )

     A(2,1) = 1.        ! coeff of E^{n+1}
     A(2,2) = alphaP    ! coeff of P^{n+1}

     b(2) = (2.*En-Enm1)\
            +csq*dtsq*Exxn\
            +(2.*pv0-Pnm1)*alphaP\
            +csq**2*dt**4*Exxxx/(12.)\
            -csq*dt**4*alphaP*QxxStar/(12.)\
            + f4

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b2,f4=",2e16.8)') i1,i2,m,b(2),f4


     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") Exxn,Exxxx,QxxStar=",6e16.8)') i1,i2,m,Exxn,Exxxx,QxxStar

     deti = 1./(A(1,1)*A(2,2) - A(1,2)*A(2,1))

     y(1) = deti *( A(2,2) * b(1) - A(1,2) * b(2))
     y(2) = deti *(-A(2,1) * b(1) + A(1,1) * b(2))

     ! write(*,'(" A(1,1),A(1,2),A(2,1),A(2,2),b(1),b(2)=",6e16.8)') A(1,1),A(1,2),A(2,1),A(2,2),b(1),b(2)


     !u(j,1) = y(1);
     !u(j,2) = y(2);


     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f4,f5,f6=",5e16.8)') i1,i2,m,f2,f3,f4,f5,f6

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") e2,e4,p2,p4=",4e16.8)') i1, i2, m, un(i1,i2,i3,ec), y(1), pn(i1,i2,i3,m), y(2)

     ! write(*,'(" En,Enm1,Exxn,pv0,Pnm1,alphaP,csq=",7e16.8)') En,Enm1,Exxn,pv0,Pnm1,alphaP,csq

     ! write(*,'(" dt,Exxxx,QxxStar=",7e16.8)') dt,Exxxx,QxxStar

     ! write(*,'(" Pxxnp1Predict,Pxxn,Pxxnm1=",7e16.8)') Pxxnp1Predict,Pxxn,Pxxnm1

     ! write(*,'(" rhsExxnp1Predict, rhsPxxnp1Predict, f5, f6=",4e16.8)') rhsExxnp1Predict, rhsPxxnp1Predict, f5, f6

     ! write(*,'(" Exxn,Exxnm1,cdtsq,Exxxx,Pxxn,Pxxnm1=",6e16.8)') Exxn,Exxnm1,cdtsq,Exxxx,Pxxn,Pxxnm1

     un(i1,i2,i3,ec) = y(1)
     pn(i1,i2,i3,m)  = y(2)

     end do
  endLoopsMask()

 else

  ! ------- 3D DISPERSIVE RECTANGULAR MULTIPLE PV (ORDER 4)------
  ! XYZ
  INFO("FD44r-dispersive-MULTI-PV");
  write(*,'(" DISPERSIVE: RECTANGULAR, 3D, ORDER 4, MULTIPLE PV  ")')

  fe=0.
  ! -- first compute some coefficients ---
  beta=0.
  do iv=0,numberOfPolarizationVectors-1
    betav(iv) = 1./( 1.+.5*dt*b1v(iv) )
    beta = beta + .5*dt*a1v(iv)*betav(iv)
    fpv(iv)=0.  ! initialize if not used
    b1tttv(iv)=b1v(iv)*b1v(iv)-b0v(iv)
    b0tttv(iv)=b1v(iv)*b0v(iv)
    a0tttv(iv)=-a0v(iv)*b1v(iv)
    a1tttv(iv)=a0v(iv)-a1v(iv)*b1v(iv)
    a2tttv(iv)=a1v(iv)
  end do

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! Advance Hz first:
    ! For now solve H_t = -(1/mu)*(  (E_y)_x - (E_x)_y )
    !   USE AB2 -- note: this is just a quadrature so stability is not an issue
    un(i1,i2,i3,hz) = u(i1,i2,i3,hz) -(dt/mu)*( 1.5*ux42r(i1,i2,i3,ey) -.5*umx42r(i1,i2,i3,ey) \
                                               -1.5*uy42r(i1,i2,i3,ex) +.5*umy42r(i1,i2,i3,ex) )

    addtForcingHz()

    ! -- loop over components of the vector --
    do m=0,nd-1
     pc=pxc+m
     ec=ex+m

     f1 = 0
     f4 = 0
     f5 = 0

     do iv=0,numberOfPolarizationVectors-1
       f2v(iv) = 0
       f3v(iv) = 0
       f6v(iv) = 0
       ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f6=",3e16.8)') i1,i2,m,f2v(iv),f3v(iv),f6v(iv)
     end do

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f1,f4,f5=",3e16.8)') i1,i2,m,f1,f4,f5

     if( addForcing.ne.0 )then
       ! 2nd order forcing used in 2nd order prediction step
       ! fe = dtsq*f(i1,i2,i3,ec)
       ! Compute fpv(iv) :
       ! getGDMForcing(ec,pc)

       ! 4th order forcing used in 4th order scheme

       OGDERIV3D( 0,0,0,0,i1,i2,i3,t, ec, e0    )
       OGDERIV3D( 1,0,0,0,i1,i2,i3,t, ec, e0t   )
       OGDERIV3D( 2,0,0,0,i1,i2,i3,t, ec, e0tt  )
       OGDERIV3D( 3,0,0,0,i1,i2,i3,t, ec, e0ttt )
       OGDERIV3D( 4,0,0,0,i1,i2,i3,t, ec, e0tttt)
       OGDERIV3D( 0,2,0,0,i1,i2,i3,t, ec, e0xx  )
       OGDERIV3D( 0,0,2,0,i1,i2,i3,t, ec, e0yy  )
       OGDERIV3D( 0,0,0,2,i1,i2,i3,t, ec, e0zz  )
       OGDERIV3D( 1,2,0,0,i1,i2,i3,t, ec, e0xxt )
       OGDERIV3D( 1,0,2,0,i1,i2,i3,t, ec, e0yyt )
       OGDERIV3D( 1,0,0,2,i1,i2,i3,t, ec, e0zzt )
       OGDERIV3D( 2,2,0,0,i1,i2,i3,t, ec, e0xxtt)
       OGDERIV3D( 2,0,2,0,i1,i2,i3,t, ec, e0yytt)
       OGDERIV3D( 2,0,0,2,i1,i2,i3,t, ec, e0zztt)
       OGDERIV3D( 0,4,0,0,i1,i2,i3,t, ec, e0xxxx)
       OGDERIV3D( 0,2,2,0,i1,i2,i3,t, ec, e0xxyy)
       OGDERIV3D( 0,2,0,2,i1,i2,i3,t, ec, e0xxzz)
       OGDERIV3D( 0,0,4,0,i1,i2,i3,t, ec, e0yyyy)
       OGDERIV3D( 0,0,2,2,i1,i2,i3,t, ec, e0yyzz)
       OGDERIV3D( 0,0,0,4,i1,i2,i3,t, ec, e0zzzz)

       pSum0tt    = 0
       pSum0ttt   = 0
       pSum0tttt  = 0
       pSum0xxtt  = 0
       pSum0yytt  = 0
       pSum0zztt  = 0

       do iv=0,numberOfPolarizationVectors-1
         pce = pc+iv*nd

         OGDERIV3D( 0,0,0,0,i1,i2,i3,t, pc, p0    )
         OGDERIV3D( 1,0,0,0,i1,i2,i3,t, pc, p0t   )
         OGDERIV3D( 2,0,0,0,i1,i2,i3,t, pc, p0tt  )
         OGDERIV3D( 3,0,0,0,i1,i2,i3,t, pc, p0ttt )
         OGDERIV3D( 4,0,0,0,i1,i2,i3,t, pc, p0tttt)
         OGDERIV3D( 0,2,0,0,i1,i2,i3,t, pc, p0xx  )
         OGDERIV3D( 0,0,2,0,i1,i2,i3,t, pc, p0yy  )
         OGDERIV3D( 0,0,0,2,i1,i2,i3,t, pc, p0zz  )
         OGDERIV3D( 1,2,0,0,i1,i2,i3,t, pc, p0xxt )
         OGDERIV3D( 1,0,2,0,i1,i2,i3,t, pc, p0yyt )
         OGDERIV3D( 1,0,0,2,i1,i2,i3,t, pc, p0zzt )
         OGDERIV3D( 2,2,0,0,i1,i2,i3,t, pc, p0xxtt)
         OGDERIV3D( 2,0,2,0,i1,i2,i3,t, pc, p0yytt)
         OGDERIV3D( 2,0,0,2,i1,i2,i3,t, pc, p0zztt)


         ! Derivatives of OG individual p eqn forcing terms fp
         fp00      = p0tt   + b1v(iv)*p0t   + b0v(iv)*p0   - a1v(iv)*e0t   - a0v(iv)*e0
         fp10      = p0ttt  + b1v(iv)*p0tt  + b0v(iv)*p0t  - a1v(iv)*e0tt  - a0v(iv)*e0t
         fp20      = p0tttt + b1v(iv)*p0ttt + b0v(iv)*p0tt - a1v(iv)*e0ttt - a0v(iv)*e0tt
         fp02x     = p0xxtt + b1v(iv)*p0xxt + b0v(iv)*p0xx - a1v(iv)*e0xxt - a0v(iv)*e0xx
         fp02y     = p0yytt + b1v(iv)*p0yyt + b0v(iv)*p0yy - a1v(iv)*e0yyt - a0v(iv)*e0yy
         fp02z     = p0zztt + b1v(iv)*p0zzt + b0v(iv)*p0zz - a1v(iv)*e0zzt - a0v(iv)*e0zz
         fp02      = fp02x  + fp02y + fp02z
         fp00v(iv) = fp00

         ! Building derivatives of full P summation terms
         pSum0tt   = pSum0tt    + p0tt
         pSum0ttt  = pSum0ttt   + p0ttt
         pSum0tttt = pSum0tttt  + p0tttt
         pSum0xxtt = pSum0xxtt  + p0xxtt
         pSum0yytt = pSum0yytt  + p0yytt
         pSum0zztt = pSum0zztt  + p0zztt

         ! Forcing on EACH individual p_ttt is (fp)_t - b1*fp :
         f2v(iv) = fp10 - b1v(iv)*fp00

         ! Forcing on EACH individual p eqn:
         f3v(iv) = dtsq * (fp00 + (dtsq/12.) * fp20)

         ! Forcing for EACH individual pxx equation at second order predictor is (xx derivative of fp)
         f6v(iv) = dtsq * fp02

         if( iv.eq.0 )then
       !    write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f6=",3e16.8)') i1,i2,m,f2v(iv),f3v(iv),f6v(iv)
         end if
       end do

       ! Build derivatives of E eqn forcing term
       fe00  = e0tt   + alphaP*pSum0tt   - csq*(e0xx   + e0yy    + e0zz  )
       fe10  = e0ttt  + alphaP*pSum0ttt  - csq*(e0xxt  + e0yyt   + e0zzt )
       fe20  = e0tttt + alphaP*pSum0tttt - csq*(e0xxtt + e0yytt  + e0zztt)
       fe02x = e0xxtt + alphaP*pSum0xxtt - csq*(e0xxxx + e0xxyy  + e0xxzz)
       fe02y = e0yytt + alphaP*pSum0yytt - csq*(e0xxyy + e0yyyy  + e0yyzz)
       fe02z = e0zztt + alphaP*pSum0zztt - csq*(e0xxzz + e0yyzz  + e0zzzz)
       fe02  = fe02x  + fe02y + fe02z
       fe    = fe00

       ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") e0xxxx,e0xxyy=",2e16.8)') i1,i2,m,e0xxxx,e0xxyy

       ! Forcing on Ettt = c^2 Etxx - alphaP Pttt is fet
       f1 = fe10

       ! Forcing on E equation is
       f4 = dtsq * (fe00 + (dtsq/12.) * (fe20 + csq * fe02))

       ! Forcing for Exx equation at second order predictor is (xx derivative of fe)
       f5 = dtsq * fe02

      ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f1,f4,f5=",3e16.8)') i1,i2,m,f1,f4,f5
     end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     rhsP = 0.
     pSum=0.

     ! Fourth order below

     elap4   = lap3d4(i1,i2,i3,ec)
     elap4m  = lap3d4m(i1,i2,i3,ec)
     elapsq2 = lap3d2Pow2(i1,i2,i3,ec)
     elap2m  = lap3d2m(i1,i2,i3,ec)

     !elap4   = ulaplacian42(i1,i2,i3,ec)
     !elap4m  = umlaplacian42(i1,i2,i3,ec)
     !elapsq2 = vlaplacian22(i1,i2,i3,ec)
     !elap2m  = umlaplacian22(i1,i2,i3,ec)

     rhsPxx = 0.
     pxxSum = 0.

     exxv  = elap4
     exxvm = elap4m

     ! First we do the second order prediction on (i) E, (ii) Exx, (iii) Pxx, (iv) Individual pk
     do iv=0,numberOfPolarizationVectors-1

       pv(iv)    = p(i1,i2,i3,m+iv*nd)
       pvm(iv)   = pm(i1,i2,i3,m+iv*nd)
       pxxv(iv)  = plap3d4(i1,i2,i3,m+iv*nd)
       pxxvm(iv) = plap3d4m(i1,i2,i3,m+iv*nd)

       ! rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + fp00v(iv)
       rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + dtsq*fp00v(iv)
       rhsP = rhsP + betav(iv)*rhspv(iv)
       pSum = pSum + 2.*pv(iv) - pvm(iv)


       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pv,pvm,rhspv,rhsP,pSum,fpv=",6e16.8)') i1,i2,m,pv(iv),pvm(iv),rhspv(iv),rhsP,pSum,fp00v(iv)

       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") ev,evm,rhspv,rhsP,pSum,fpv=",6e16.8)') i1,i2,m,ev,evm,a1v(iv),a0v(iv),b1v(iv),b0v(iv)


       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") first term=",1e16.8)') i1,i2,m,2.*pv(iv)-pvm(iv)

       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") two terms",2e16.8)') i1,i2,m,.5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ),dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev )


       rhspxxv(iv) = 2.*pxxv(iv)-pxxvm(iv) + .5*dt*( b1v(iv)*pxxvm(iv) -a1v(iv)*exxvm ) + dtSq*( -b0v(iv)*pxxv(iv) + a0v(iv)*exxv ) + f6v(iv)
       rhsPxx = rhsPxx + betav(iv)*rhspxxv(iv)
       pxxSum = pxxSum + 2.*pxxv(iv) - pxxvm(iv)




     end do

     rhsE   = maxwell2dr44me(i1,i2,i3,ec)     + alphaP*( pSum   - rhsP   ) + dtsq * fe00
     rhsE   = 2.*ev   - evm   + cdtsq*elap4   + alphaP*( pSum   - rhsP   ) + dtsq * fe00

     rhsExx = 2.*exxv - exxvm + cdtsq*elapsq2 + alphaP*( pxxSum - rhsPxx ) + f5



     evn   = rhsE   / (1.+ alphaP*beta)
     exxvn = rhsExx / (1.+ alphaP*beta)



     ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") rhsE,evn,fe00=",3e16.8)') i1,i2,m,rhsE,evn,fe00


     ! Update x derivative of P
     Pxxn  = beta * exxvn + rhsPxx

     write(*,'(" (i2,i2,m)=(",i3,i3,i2,") exxvn,Pxxn=",2e16.8)') i1,i2,m,exxvn,Pxxn


     ! Now we predict individual pk (pvn) and their third time derivative to second order
     PtttStar = 0

     do iv=0,numberOfPolarizationVectors-1
       pvn(iv) = betav(iv)*( .5*dt*a1v(iv)*evn + rhspv(iv) )
       ! Now update PtttStar for correction terms and use in Ettt
       ptttStarv(iv) = b1tttv(iv)*(pvn(iv)-pvm(iv))/(2.*dt) + b0tttv(iv)*pv(iv) + a0tttv(iv)*ev + a1tttv(iv)*(evn-evm)/(2.*dt) + a2tttv(iv)*(evn-2*ev+evm)/dtsq + f2v(iv)
       PtttStar = PtttStar + ptttStarv(iv)


      ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pvm,pv,pvn,evm,ev,evn=",6e16.8)') i1,i2,m,pvm(iv),pv(iv),pvn(iv),evm,ev,evn

      ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") evn-evm, evn-2*ev+evm=",2e16.8)') i1,i2,m,evn-evm,evn-2*ev+evm

      !  write(*,'(" (i2,i2,m)=(",i3,i3,i2,") 1diff,2diff,3diff,f2v(iv)=",4e16.8)') i1,i2,m,(pvn(iv)-pvm(iv))/(2.*dt),(evn-evm)/(2.*dt),(evn-2*ev+evm)/dtsq,f2v(iv)
      ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pvn(iv),ptttStarv,dtsq=",3e16.8)') i1,i2,m,pvn(iv),ptttStarv(iv),dtsq
     end do

     ! Second Order Updates Complete, now we construct necessary terms
     ! LapPtt using prediction
     QxxStar = (Pxxn  - pxxSum)/dtsq
     EtxxStar    = (exxvn -  exxvm)/(2.*dt)
     EtttStar    = csq*EtxxStar - alphaP*PtttStar + f1

     write(*,'(" (i2,i2,m)=(",i3,i3,i2,") QxxStar,exxvn,Pxxn,EtxxStar,PtttStar,EtttStar=",7e16.8)') i1,i2,m,QxxStar,exxvn,Pxxn,EtxxStar,PtttStar,EtttStar,f1

     rhsP4 = 0
     LHSev = 0

     do iv=0,numberOfPolarizationVectors-1
         ! Coeff for pk(n+1) for  LHS of invidiual p equation
         LHSpv(iv) = 1+ b1v(iv)*dt/2.  + b0v(iv)*dtsq/12.

         ! Build coeff for ev for left hand side of full P equation
         LHSev     = LHSev + ((-a1v(iv)*dt/2.  -a0v(iv)*dtsq/12.)/LHSpv(iv))

         ! LHS for pk for pk equation
         rhspv(iv) = ((2.*pv(iv)-pvm(iv))\
         + b1v(iv)*dt*pvm(iv)/(2.)\
         - dtsq*b0v(iv)*pv(iv)\
         + dtsq*a0v(iv)*ev\
         - a1v(iv)*(dt/2.)*evm\
         - a1v(iv)*(dt**4/12.)*( EtttStar )\
         + b1v(iv)*(dt**4/12.)*( ptttStarv(iv) )\
         + (dtsq/12.)*( -b0v(iv)*(0. -2.*pv(iv)+pvm(iv)) + a0v(iv)*(0. -2.*ev+evm) )\
         + f3v(iv))

       ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") LHSpv,rhspv=",2e16.8)') i1,i2,m,LHSpv(iv),rhspv(iv)


        rhsP4 = rhsP4 + ( rhspv(iv) )/LHSpv(iv)

     end do

     ! We have now built the equation for P
     A(1,1) = LHSev
     A(1,2) = 1.

     b(1) = rhsP4

     ! Now we build the equation for E in terms of E and P

     A(2,1) = 1.        ! coeff of E^{n+1}
     A(2,2) = alphaP    ! coeff of P^{n+1}

     ! Note that Psum - rhsP here is same as for second order code
     !  (but using 4th order accurate version of p to compute them)
     b(2) = (2.*ev-evm)\
            +csq*dtsq*elap4\
            +alphaP*( Psum )\
            +csq**2*dt**4*elapsq2/(12.)\
            -csq*dt**4*alphaP*QxxStar/(12.)\
            + f4


     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b2,f4=",2e16.8)') i1,i2,m,b(2),f4


     deti = 1./(A(1,1)*A(2,2) - A(1,2)*A(2,1))

     y(1) = deti *( A(2,2) * b(1) - A(1,2) * b(2))
     y(2) = deti *(-A(2,1) * b(1) + A(1,1) * b(2))

     write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b1,b2,y1,y2,deti=",5e16.8)') i1,i2,m,b(1),b(2),y(1),y(2),deti

     ! Update E^{n+1}
     un(i1,i2,i3,ec) = y(1)
     evn             = y(1)

     ! Update pk using new E^{n+1} = evn
     do iv=0,numberOfPolarizationVectors-1
       rhspv(iv) = rhspv(iv) + (a1v(iv) * dt/(2.) + a0v(iv)*dtsq/(12.))*evn
       pn(i1,i2,i3,m+iv*nd)  = (1/LHSpv(iv)) * rhspv(iv)
     end do



   end do ! m=0,1
  endLoopsMask()

 end if
#endMacro


! ===========================================================================================
! Macro:     DISPERSIVE: CURVILINEAR, 3D, ORDER 4
! ===========================================================================================
#beginMacro updateCurvilinear3dOrder4Dispersive()

 if( .true. .and. numberOfPolarizationVectors.eq.1 )then
  INFO("FD44-dispersive");

  write(*,'(" DISPERSIVE: CURVILINEAR, 3D, ORDER 4 (ONE PV) ")')

  fp=0
  fe=0.

  ! Coefficients in equation for Pttt
  ! Pttt = b1ttt*P_t + b0ttt*P + a0ttt*E + a1ttt*E_t + a2ttt*Ett
  b1ttt=b1*b1-b0
  b0ttt=b1*b0
  a0ttt=-a0*b1
  a1ttt=a0-a1*b1
  a2ttt=a1

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! Advance Hz first:
    ! For now solve H_t = -(1/mu)*(  (E_y)_x - (E_x)_y )
    !   USE AB2 -- note: this is just a quadrature so stability is not an issue
    un(i1,i2,i3,hz) = u(i1,i2,i3,hz) -(dt/mu)*( 1.5*ux42r(i1,i2,i3,ey) -.5*umx42r(i1,i2,i3,ey) \
                                               -1.5*uy42r(i1,i2,i3,ex) +.5*umy42r(i1,i2,i3,ex) )

    addtForcingHz()

    do m=0,nd-1
     pc=pxc+m
     ec=ex+m

     !if( addForcing.ne.0 )then
     !  fe = dtsq*f(i1,i2,i3,ec)
     !  getGDMForcing(ec,pc)
     !  fp=fpv(0)
     !end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     pv0 = p(i1,i2,i3,m)
     pvm0=pm(i1,i2,i3,m)

     !rhsE = maxwell3dr44me(i1,i2,i3,ec) + alphaP*(2.*pv0-pvm0) + fe
     !rhsP = 2.*pv0-pvm0 + .5*dt*( b1*pvm0 -a1*evm ) + dt*dt*( -b0*pv0 + a0*ev ) + fp

     !deti = 1./(1.+.5*dt*(b1+alphaP*a1))
     !un(i1,i2,i3,ec) = ((1.+.5*dt*b1)*rhsE -alphaP*rhsP)*deti
     !pn(i1,i2,i3,m)  = (.5*a1*dt*rhsE            + rhsP)*deti

     ! New 4th Order code starts here
     ! pvx = px22r(i1,i2,i3,m)


     f1 = 0
     f2 = 0
     f3 = 0
     f4 = 0
     f5 = 0
     f6 = 0


     if( addForcing.ne.0 )then
       OGDERIV3D( 0,0,0,0,i1,i2,i3,t, pc, p0    )
       OGDERIV3D( 1,0,0,0,i1,i2,i3,t, pc, p0t   )
       OGDERIV3D( 2,0,0,0,i1,i2,i3,t, pc, p0tt  )
       OGDERIV3D( 3,0,0,0,i1,i2,i3,t, pc, p0ttt )
       OGDERIV3D( 4,0,0,0,i1,i2,i3,t, pc, p0tttt)
       OGDERIV3D( 0,2,0,0,i1,i2,i3,t, pc, p0xx  )
       OGDERIV3D( 0,0,2,0,i1,i2,i3,t, pc, p0yy  )
       OGDERIV3D( 0,0,0,2,i1,i2,i3,t, pc, p0zz  )
       OGDERIV3D( 1,2,0,0,i1,i2,i3,t, pc, p0xxt )
       OGDERIV3D( 1,0,2,0,i1,i2,i3,t, pc, p0yyt )
       OGDERIV3D( 1,0,0,2,i1,i2,i3,t, pc, p0zzt )
       OGDERIV3D( 2,2,0,0,i1,i2,i3,t, pc, p0xxtt)
       OGDERIV3D( 2,0,2,0,i1,i2,i3,t, pc, p0yytt)
       OGDERIV3D( 2,0,0,2,i1,i2,i3,t, pc, p0zztt)



       OGDERIV3D( 0,0,0,0,i1,i2,i3,t, ec, e0    )
       OGDERIV3D( 1,0,0,0,i1,i2,i3,t, ec, e0t   )
       OGDERIV3D( 2,0,0,0,i1,i2,i3,t, ec, e0tt  )
       OGDERIV3D( 3,0,0,0,i1,i2,i3,t, ec, e0ttt )
       OGDERIV3D( 4,0,0,0,i1,i2,i3,t, ec, e0tttt)
       OGDERIV3D( 0,2,0,0,i1,i2,i3,t, ec, e0xx  )
       OGDERIV3D( 0,0,2,0,i1,i2,i3,t, ec, e0yy  )
       OGDERIV3D( 0,0,0,2,i1,i2,i3,t, ec, e0zz  )
       OGDERIV3D( 1,2,0,0,i1,i2,i3,t, ec, e0xxt )
       OGDERIV3D( 1,0,2,0,i1,i2,i3,t, ec, e0yyt )
       OGDERIV3D( 1,0,0,2,i1,i2,i3,t, ec, e0zzt )
       OGDERIV3D( 2,2,0,0,i1,i2,i3,t, ec, e0xxtt)
       OGDERIV3D( 2,0,2,0,i1,i2,i3,t, ec, e0yytt)
       OGDERIV3D( 2,0,0,2,i1,i2,i3,t, ec, e0zztt)
       OGDERIV3D( 0,4,0,0,i1,i2,i3,t, ec, e0xxxx)
       OGDERIV3D( 0,2,2,0,i1,i2,i3,t, ec, e0xxyy)
       OGDERIV3D( 0,2,0,2,i1,i2,i3,t, ec, e0xxzz)
       OGDERIV3D( 0,0,4,0,i1,i2,i3,t, ec, e0yyyy)
       OGDERIV3D( 0,0,2,2,i1,i2,i3,t, ec, e0yyzz)
       OGDERIV3D( 0,0,0,4,i1,i2,i3,t, ec, e0zzzz)


       fp00  = p0tt   + b1*p0t   + b0*p0   - a1*e0t   - a0*e0
       fp10  = p0ttt  + b1*p0tt  + b0*p0t  - a1*e0tt  - a0*e0t
       fp20  = p0tttt + b1*p0ttt + b0*p0tt - a1*e0ttt - a0*e0tt
       fp02x = p0xxtt + b1*p0xxt + b0*p0xx - a1*e0xxt - a0*e0xx
       fp02y = p0yytt + b1*p0yyt + b0*p0yy - a1*e0yyt - a0*e0yy
       fp02z = p0yytt + b1*p0zzt + b0*p0zz - a1*e0zzt - a0*e0zz
       fp02  = fp02x  + fp02y + fp02z

       fe00  = e0tt   + alphaP*p0tt   - csq*(e0xx   + e0yy   + e0zz  )
       fe10  = e0ttt  + alphaP*p0ttt  - csq*(e0xxt  + e0yyt  + e0zzt )
       fe20  = e0tttt + alphaP*p0tttt - csq*(e0xxtt + e0yytt + e0zztt)
       fe02x = e0xxtt + alphaP*p0xxtt - csq*(e0xxxx + e0xxyy + e0xxzz)
       fe02y = e0yytt + alphaP*p0yytt - csq*(e0xxyy + e0yyyy + e0yyzz)
       fe02z = e0zztt + alphaP*p0zztt - csq*(e0zzzz + e0xxzz + e0yyzz)
       fe02  = fe02x  + fe02y + fe02z

       ! f1 does not appear in this version of the 4th order scheme
       f1 = 0

       ! Forcing on P_ttt is (fp)_t - b1*fp :
       ! f2 = fp(1,0,x) - b1 * fp(0,0,x);
       ! f2 = p0t - b1 * p0 (incorrect)
       f2 = fp10 - b1*fp00

       ! Forcing on P eqn:
       ! fp + (dt^2/12) * fp_tt
       ! f3 = dt^2 * (fp(0,0,x) + (dt^2/12) * (- a1 * fe(1,0,x) + fp(2,0,x)));  % FIX ME *wdh* Jan 8, 2018
       ! f3 = dtsq * (p0 + (dtsq/12.) * (-a1 * e0t + p0tt )) incorrect
       f3 = dtsq * (fp00 + (dtsq/12.) * (-a1 * fe10 + fp20))

       ! Forcing on E equation is
       ! fe + (dt^2)/12 * ( fe_tt + c^2*Delta(fe) )
       ! f4 = dt^2 * (fe(0,0,x) + (dt^2/12) * (fe(2,0,x) + c^2 * fe(0,2,x)));
       ! f4 = dtsq * (e0 + (dtsq/12.) * (e0tt + csq * (e0xx + e0yy))) incorrect
       f4 = dtsq * (fe00 + (dtsq/12.) * (fe20 + csq * fe02))


       ! Forcing for Exx equation at second order predictor is (xx derivative
       ! of fe)
       ! f5 = dt^2 * fe(0,2,x);
       ! f5 = dtsq * (e0xx + e0yy)
       f5 = dtsq * fe02

       ! Forcing for Pxx equation at second order predictor is (xx derivative
       ! of fp)
       ! f6 = dt^2 * fp(0,2,x);
       ! f6 = dtsq * (p0xx + p0yy)
       f6 = dtsq * fp02

     end if

     ! YOU ARE HERE 3D

     elap4   = lap3d4(i1,i2,i3,ec)
     elap4m  = lap3d4m(i1,i2,i3,ec)
     elapsq2 = lap3d2Pow2(i1,i2,i3,ec)
     elap2m  = lap3d2m(i1,i2,i3,ec)
     plap2   = plap3d2(i1,i2,i3,m)
     plap2m  = plap3d2m(i1,i2,i3,m)
     plap4   = plap3d4(i1,i2,i3,m)
     plap4m  = plap3d4m(i1,i2,i3,m)
     ! Matlab code starts here

     ! Curvilinear
     elap4   = ulaplacian43(i1,i2,i3,ec)
     elap4m  = umlaplacian43(i1,i2,i3,ec)
     elapsq2 = vlaplacian23(i1,i2,i3,ec)
     elap2m  = umlaplacian23(i1,i2,i3,ec)
     plap2   = plaplacian23(i1,i2,i3,m)
     plap2m  = pmlaplacian23(i1,i2,i3,m)
     plap4   = plaplacian43(i1,i2,i3,m)
     plap4m  = pmlaplacian43(i1,i2,i3,m)

     !write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f6=",3e16.8)') i1,i2,m,f2,f3,f6
     !write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f1,f4,f5=",3e16.8)') i1,i2,m,f1,f4,f5

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") elap4,elap4m,elapsq2,elap2m,plap2,plap2m,plap4,plap4m=",4e16.8)') i1,i2,m,elap4,elap4m,elapsq2,elap2m
     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") plap2,plap2m,plap4,plap4m=",4e16.8)') i1,i2,m,plap2,plap2m,plap4,plap4m

     Exxxx   = elapsq2

     En   = ev
     Enm1 = evm

     ! Pn   = pv0
     Pnm1 = pvm0

     ! compute c^2*Delta( E_t ) using
     ! Et = (1.5*En - 2*Enm1 +.5*Enm2)/dt;

     Exxn   = elap4
     Exxnm1 = elap4m


     betaP = 1/( 1+ (b1+alphaP*a1)*dt/2 )
     Pxxn   = plap4
     Pxxnm1 = plap4m



     ! Approximation for P_ttt
     ! Pttt^*    = b1ttt*(Pnp1 -Pnm1)/(2.*dt) + b0ttt*pv0 + a0ttt*En + a1ttt*(Enp1-Enm1)/(2.*dt) + a2ttt*(Enp1-2*En+Enm1)/dt^2;
     cStar=1  !  cStar=1 : use Pttt^*, 0: use old

     PtttStarRhs = b1ttt*(0.   -Pnm1)/(2.*dt) + b0ttt*pv0 + a0ttt*En + a1ttt*(0.  -Enm1)/(2.*dt) + a2ttt*(0.  -2*En+Enm1)/dtsq + f2

     ! Use second order update to predict new E and P for future space derivatives
     !   fexx = fe;
     !   fpxx = fp;

     ! These are actually Laplacians
     rhsExxnp1Predict = 2.*Exxn - Exxnm1 + cdtsq*Exxxx + alphaP*(2.*Pxxn - Pxxnm1) + f5
     rhsPxxnp1Predict = 2.*Pxxn - Pxxnm1 + .5*dt*( b1*Pxxnm1 - a1*Exxnm1 ) + dtsq*( -b0*Pxxn + a0*Exxn ) + f6
     deti = 1./(1.+.5*dt*(b1+alphaP*a1))
     Exxnp1Predict = ((1.+.5*dt*b1)*rhsExxnp1Predict - alphaP*rhsPxxnp1Predict)*deti
     Pxxnp1Predict = (.5*a1*dt*rhsExxnp1Predict            + rhsPxxnp1Predict)*deti


     ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") Exxnp1Predict,Pxxnp1Predict=",2e16.8)') i1,i2,m,Exxnp1Predict,Pxxnp1Predict





     QxxStar = (Pxxnp1Predict - 2.*Pxxn + Pxxnm1)/dtsq
     Etxx = (Exxnp1Predict - Exxnm1)/(2.*dt)

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") QxxStar,Exxn,Pxxn=",3e16.8)') i1,i2,m,QxxStar,Exxnp1Predict,Pxxnp1Predict

     ! A = zeros(2,2);
     ! b = zeros(2,1);

     !------------------------------- P EQUATION ----------------------------------------
     ! D+tD-t P^n = -b1*P_t -b0*P + a0*E + a1*E_t + (dt^2/12)*( D+tD-t Q^n )
     !            = -b1*[D0t(P^n) ]  - b0*P^n + a0*E^n  + a1*[ D0t(E^n) - (dt^2/12)*E_tttStar ] ...
     !              + (dt^2/12)*( -b0*D+tD-t(P^n) + a0*D+tD-t(E^n) )
     !              + (b1+alphaP*a1)*(dt^2)/12*P_ttt^*


     !
     ! We can either use:
     !       Qt = (Q^{n+1}-Q^{n-1})/(2*dt)   <-- adds a term to A(1,3)
     !       Qt = ( 1.5*Qn-2.*Qnm1 +.5*Qnm2 )/dt;  <-- remove term from A(1,3)

     A(1,1) =   -a1*dt/2.  -a0*dtsq/12.   + (b1+alphaP*a1)*(dtsq/12.)*dtsq*( -a1ttt/(2.*dt) - a2ttt/dtsq )
     A(1,2) = 1+ b1*dt/2.  +b0*dtsq/12.   + (b1+alphaP*a1)*(dtsq/12.)*dtsq*( -b1ttt/(2.*dt) )

     b(1) = (2.*pv0-Pnm1)\
         + b1*dt*Pnm1/(2.)\
         - dtsq*b0*pv0\
         + dtsq*a0*En\
         - a1*(dt/2.)*Enm1\
         - a1*(dt**4/12.)*( csq*Etxx )\
         + (dtsq/12.)*( -b0*(0. -2.*pv0+Pnm1) + a0*(0. -2.*En+Enm1) )\
         + (b1+alphaP*a1)*(dtsq/12.)*dtsq*( PtttStarRhs )\
         + f3


     !write(*,'(" (i1,i2,m)=(",i3,i3,i2,") rhspv=",1e16.8)') i1,i2,m,b(1)




      ! ------------------------------- E EQUATION ----------------------------------------
      ! D+tD-t E^n = c^2*Delta(E^n) -alphaP*D+tD-t P^n + (dt^2/(12)*( [c^2*Delta]^2 E^n  -alphaP*c^2*Qxx^*  )

     A(2,1) = 1.        ! coeff of E^{n+1}
     A(2,2) = alphaP    ! coeff of P^{n+1}

     b(2) = (2.*En-Enm1)\
            +csq*dtsq*Exxn\
            +(2.*pv0-Pnm1)*alphaP\
            +csq**2*dt**4*Exxxx/(12.)\
            -csq*dt**4*alphaP*QxxStar/(12.)\
            + f4

     !write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b2,f4=",2e16.8)') i1,i2,m,b(2),f4


     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") Exxn,Exxxx,QxxStar=",6e16.8)') i1,i2,m,Exxn,Exxxx,QxxStar

     deti = 1./(A(1,1)*A(2,2) - A(1,2)*A(2,1))

     y(1) = deti *( A(2,2) * b(1) - A(1,2) * b(2))
     y(2) = deti *(-A(2,1) * b(1) + A(1,1) * b(2))

     !write(*,'(" A(1,1),A(1,2),A(2,1),A(2,2),b(1),b(2)=",6e16.8)') A(1,1),A(1,2),A(2,1),A(2,2),b(1),b(2)


     !u(j,1) = y(1);
     !u(j,2) = y(2);


     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f4,f5,f6=",5e16.8)') i1,i2,m,f2,f3,f4,f5,f6

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") e2,e4,p2,p4=",4e16.8)') i1, i2, m, un(i1,i2,i3,ec), y(1), pn(i1,i2,i3,m), y(2)

     ! write(*,'(" En,Enm1,Exxn,pv0,Pnm1,alphaP,csq=",7e16.8)') En,Enm1,Exxn,pv0,Pnm1,alphaP,csq

     ! write(*,'(" dt,Exxxx,QxxStar=",7e16.8)') dt,Exxxx,QxxStar

     ! write(*,'(" Pxxnp1Predict,Pxxn,Pxxnm1=",7e16.8)') Pxxnp1Predict,Pxxn,Pxxnm1

     ! write(*,'(" rhsExxnp1Predict, rhsPxxnp1Predict, f5, f6=",4e16.8)') rhsExxnp1Predict, rhsPxxnp1Predict, f5, f6

     ! write(*,'(" Exxn,Exxnm1,cdtsq,Exxxx,Pxxn,Pxxnm1=",6e16.8)') Exxn,Exxnm1,cdtsq,Exxxx,Pxxn,Pxxnm1

     un(i1,i2,i3,ec) = y(1)
     pn(i1,i2,i3,m)  = y(2)

     end do
  endLoopsMask()

 else

  ! ------- 3D DISPERSIVE CURVILINEAR MULTIPLE PV (ORDER 4)------
  ! XYZ
  INFO("FD44-dispersive-MULTI-PV");
  write(*,'(" DISPERSIVE: CURVILINEAR, 3D, ORDER 4, MULTIPLE PV  ")')

  fe=0.
  ! -- first compute some coefficients ---
  beta=0.
  do iv=0,numberOfPolarizationVectors-1
    betav(iv) = 1./( 1.+.5*dt*b1v(iv) )
    beta = beta + .5*dt*a1v(iv)*betav(iv)
    fpv(iv)=0.  ! initialize if not used
    b1tttv(iv)=b1v(iv)*b1v(iv)-b0v(iv)
    b0tttv(iv)=b1v(iv)*b0v(iv)
    a0tttv(iv)=-a0v(iv)*b1v(iv)
    a1tttv(iv)=a0v(iv)-a1v(iv)*b1v(iv)
    a2tttv(iv)=a1v(iv)
  end do

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

    ! Advance Hz first:
    ! For now solve H_t = -(1/mu)*(  (E_y)_x - (E_x)_y )
    !   USE AB2 -- note: this is just a quadrature so stability is not an issue
    un(i1,i2,i3,hz) = u(i1,i2,i3,hz) -(dt/mu)*( 1.5*ux42r(i1,i2,i3,ey) -.5*umx42r(i1,i2,i3,ey) \
                                               -1.5*uy42r(i1,i2,i3,ex) +.5*umy42r(i1,i2,i3,ex) )

    addtForcingHz()

    ! -- loop over components of the vector --
    do m=0,nd-1
     pc=pxc+m
     ec=ex+m

     f1 = 0
     f4 = 0
     f5 = 0

     do iv=0,numberOfPolarizationVectors-1
       f2v(iv) = 0
       f3v(iv) = 0
       f6v(iv) = 0
       ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f6=",3e16.8)') i1,i2,m,f2v(iv),f3v(iv),f6v(iv)
     end do

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f1,f4,f5=",3e16.8)') i1,i2,m,f1,f4,f5

     if( addForcing.ne.0 )then
       ! 2nd order forcing used in 2nd order prediction step
       ! fe = dtsq*f(i1,i2,i3,ec)
       ! Compute fpv(iv) :
       ! getGDMForcing(ec,pc)

       ! 4th order forcing used in 4th order scheme

       OGDERIV3D( 0,0,0,0,i1,i2,i3,t, ec, e0    )
       OGDERIV3D( 1,0,0,0,i1,i2,i3,t, ec, e0t   )
       OGDERIV3D( 2,0,0,0,i1,i2,i3,t, ec, e0tt  )
       OGDERIV3D( 3,0,0,0,i1,i2,i3,t, ec, e0ttt )
       OGDERIV3D( 4,0,0,0,i1,i2,i3,t, ec, e0tttt)
       OGDERIV3D( 0,2,0,0,i1,i2,i3,t, ec, e0xx  )
       OGDERIV3D( 0,0,2,0,i1,i2,i3,t, ec, e0yy  )
       OGDERIV3D( 0,0,0,2,i1,i2,i3,t, ec, e0zz  )
       OGDERIV3D( 1,2,0,0,i1,i2,i3,t, ec, e0xxt )
       OGDERIV3D( 1,0,2,0,i1,i2,i3,t, ec, e0yyt )
       OGDERIV3D( 1,0,0,2,i1,i2,i3,t, ec, e0zzt )
       OGDERIV3D( 2,2,0,0,i1,i2,i3,t, ec, e0xxtt)
       OGDERIV3D( 2,0,2,0,i1,i2,i3,t, ec, e0yytt)
       OGDERIV3D( 2,0,0,2,i1,i2,i3,t, ec, e0zztt)
       OGDERIV3D( 0,4,0,0,i1,i2,i3,t, ec, e0xxxx)
       OGDERIV3D( 0,2,2,0,i1,i2,i3,t, ec, e0xxyy)
       OGDERIV3D( 0,2,0,2,i1,i2,i3,t, ec, e0xxzz)
       OGDERIV3D( 0,0,4,0,i1,i2,i3,t, ec, e0yyyy)
       OGDERIV3D( 0,0,2,2,i1,i2,i3,t, ec, e0yyzz)
       OGDERIV3D( 0,0,0,4,i1,i2,i3,t, ec, e0zzzz)

       pSum0tt    = 0
       pSum0ttt   = 0
       pSum0tttt  = 0
       pSum0xxtt  = 0
       pSum0yytt  = 0
       pSum0zztt  = 0

       do iv=0,numberOfPolarizationVectors-1
         pce = pc+iv*nd

         OGDERIV3D( 0,0,0,0,i1,i2,i3,t, pc, p0    )
         OGDERIV3D( 1,0,0,0,i1,i2,i3,t, pc, p0t   )
         OGDERIV3D( 2,0,0,0,i1,i2,i3,t, pc, p0tt  )
         OGDERIV3D( 3,0,0,0,i1,i2,i3,t, pc, p0ttt )
         OGDERIV3D( 4,0,0,0,i1,i2,i3,t, pc, p0tttt)
         OGDERIV3D( 0,2,0,0,i1,i2,i3,t, pc, p0xx  )
         OGDERIV3D( 0,0,2,0,i1,i2,i3,t, pc, p0yy  )
         OGDERIV3D( 0,0,0,2,i1,i2,i3,t, pc, p0zz  )
         OGDERIV3D( 1,2,0,0,i1,i2,i3,t, pc, p0xxt )
         OGDERIV3D( 1,0,2,0,i1,i2,i3,t, pc, p0yyt )
         OGDERIV3D( 1,0,0,2,i1,i2,i3,t, pc, p0zzt )
         OGDERIV3D( 2,2,0,0,i1,i2,i3,t, pc, p0xxtt)
         OGDERIV3D( 2,0,2,0,i1,i2,i3,t, pc, p0yytt)
         OGDERIV3D( 2,0,0,2,i1,i2,i3,t, pc, p0zztt)


         ! Derivatives of OG individual p eqn forcing terms fp
         fp00      = p0tt   + b1v(iv)*p0t   + b0v(iv)*p0   - a1v(iv)*e0t   - a0v(iv)*e0
         fp10      = p0ttt  + b1v(iv)*p0tt  + b0v(iv)*p0t  - a1v(iv)*e0tt  - a0v(iv)*e0t
         fp20      = p0tttt + b1v(iv)*p0ttt + b0v(iv)*p0tt - a1v(iv)*e0ttt - a0v(iv)*e0tt
         fp02x     = p0xxtt + b1v(iv)*p0xxt + b0v(iv)*p0xx - a1v(iv)*e0xxt - a0v(iv)*e0xx
         fp02y     = p0yytt + b1v(iv)*p0yyt + b0v(iv)*p0yy - a1v(iv)*e0yyt - a0v(iv)*e0yy
         fp02z     = p0zztt + b1v(iv)*p0zzt + b0v(iv)*p0zz - a1v(iv)*e0zzt - a0v(iv)*e0zz
         fp02      = fp02x  + fp02y + fp02z
         fp00v(iv) = fp00

         ! Building derivatives of full P summation terms
         pSum0tt   = pSum0tt    + p0tt
         pSum0ttt  = pSum0ttt   + p0ttt
         pSum0tttt = pSum0tttt  + p0tttt
         pSum0xxtt = pSum0xxtt  + p0xxtt
         pSum0yytt = pSum0yytt  + p0yytt
         pSum0zztt = pSum0zztt  + p0zztt

         ! Forcing on EACH individual p_ttt is (fp)_t - b1*fp :
         f2v(iv) = fp10 - b1v(iv)*fp00

         ! Forcing on EACH individual p eqn:
         f3v(iv) = dtsq * (fp00 + (dtsq/12.) * fp20)

         ! Forcing for EACH individual pxx equation at second order predictor is (xx derivative of fp)
         f6v(iv) = dtsq * fp02

         if( iv.eq.0 )then
           !write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f2,f3,f6=",3e16.8)') i1,i2,m,f2v(iv),f3v(iv),f6v(iv)
         end if
       end do

       ! Build derivatives of E eqn forcing term
       fe00  = e0tt   + alphaP*pSum0tt   - csq*(e0xx   + e0yy    + e0zz  )
       fe10  = e0ttt  + alphaP*pSum0ttt  - csq*(e0xxt  + e0yyt   + e0zzt )
       fe20  = e0tttt + alphaP*pSum0tttt - csq*(e0xxtt + e0yytt  + e0zztt)
       fe02x = e0xxtt + alphaP*pSum0xxtt - csq*(e0xxxx + e0xxyy  + e0xxzz)
       fe02y = e0yytt + alphaP*pSum0yytt - csq*(e0xxyy + e0yyyy  + e0yyzz)
       fe02z = e0zztt + alphaP*pSum0zztt - csq*(e0xxzz + e0yyzz  + e0zzzz)
       fe02  = fe02x  + fe02y + fe02z
       fe    = fe00

       ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") e0xxxx,e0xxyy=",2e16.8)') i1,i2,m,e0xxxx,e0xxyy

       ! Forcing on Ettt = c^2 Etxx - alphaP Pttt is fet
       f1 = fe10

       ! Forcing on E equation is
       f4 = dtsq * (fe00 + (dtsq/12.) * (fe20 + csq * fe02))

       ! Forcing for Exx equation at second order predictor is (xx derivative of fe)
       f5 = dtsq * fe02

       ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") f1,f4,f5=",3e16.8)') i1,i2,m,f1,f4,f5
     end if

     ! GDM:
     !   (E^{n+1} -2 E^n + E^{n-1})/dt^2 = c^2*Delta(E) -alphaP*(P^{n+1} -2 P^n + P^{n-1})/dt^2
     !   (P^{n+1} -2 P^n + P^{n-1})/dt^2 + b1* (P^{n+1} - P^{n-1})/(2*dt) + b0*P^n =
     !                             a0*E^n + a1*(E^{n+1} - E^{n-1})/(2*dt)
     ! =>
     !            E^{n+1} +       alphaP*P^{n+1} = rhsE
     !  -.5*a1*dt*E^{n+1} + (1+.5*b1*dt)*P^{n+1} = rhsP
     !
     ev = u(i1,i2,i3,ec)
     evm=um(i1,i2,i3,ec)

     rhsP = 0.
     pSum=0.

     ! Fourth order below

     elap4   = lap3d4(i1,i2,i3,ec)
     elap4m  = lap3d4m(i1,i2,i3,ec)
     elapsq2 = lap3d2Pow2(i1,i2,i3,ec)
     elap2m  = lap3d2m(i1,i2,i3,ec)

     elap4   = ulaplacian43(i1,i2,i3,ec)
     elap4m  = umlaplacian43(i1,i2,i3,ec)
     elapsq2 = vlaplacian23(i1,i2,i3,ec)
     elap2m  = umlaplacian23(i1,i2,i3,ec)

     rhsPxx = 0.
     pxxSum = 0.

     exxv  = elap4
     exxvm = elap4m

     ! First we do the second order prediction on (i) E, (ii) Exx, (iii) Pxx, (iv) Individual pk
     do iv=0,numberOfPolarizationVectors-1

       pv(iv)    = p(i1,i2,i3,m+iv*nd)
       pvm(iv)   = pm(i1,i2,i3,m+iv*nd)
       pxxv(iv)  = plap3d4(i1,i2,i3,m+iv*nd)
       pxxvm(iv) = plap3d4m(i1,i2,i3,m+iv*nd)

       pxxv(iv)  = plaplacian43(i1,i2,i3,m+iv*nd)
       pxxvm(iv) = pmlaplacian43(i1,i2,i3,m+iv*nd)

       ! rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + fp00v(iv)
       rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + dtsq*fp00v(iv)
       rhsP = rhsP + betav(iv)*rhspv(iv)
       pSum = pSum + 2.*pv(iv) - pvm(iv)


       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pv,pvm,rhspv,rhsP,pSum,fpv=",6e16.8)') i1,i2,m,pv(iv),pvm(iv),rhspv(iv),rhsP,pSum,fp00v(iv)

       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") ev,evm,rhspv,rhsP,pSum,fpv=",6e16.8)') i1,i2,m,ev,evm,a1v(iv),a0v(iv),b1v(iv),b0v(iv)


       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") first term=",1e16.8)') i1,i2,m,2.*pv(iv)-pvm(iv)

       !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") two terms",2e16.8)') i1,i2,m,.5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ),dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev )


       rhspxxv(iv) = 2.*pxxv(iv)-pxxvm(iv) + .5*dt*( b1v(iv)*pxxvm(iv) -a1v(iv)*exxvm ) + dtSq*( -b0v(iv)*pxxv(iv) + a0v(iv)*exxv ) + f6v(iv)
       rhsPxx = rhsPxx + betav(iv)*rhspxxv(iv)
       pxxSum = pxxSum + 2.*pxxv(iv) - pxxvm(iv)




     end do

     ! rhsE   = maxwell2d44me(i1,i2,i3,ec)     + alphaP*( pSum   - rhsP   ) + dtsq * fe00
     rhsE   = 2.*ev   - evm   + cdtsq*elap4   + alphaP*( pSum   - rhsP   ) + dtsq * fe00

     rhsExx = 2.*exxv - exxvm + cdtsq*elapsq2 + alphaP*( pxxSum - rhsPxx ) + f5



     evn   = rhsE   / (1.+ alphaP*beta)
     exxvn = rhsExx / (1.+ alphaP*beta)



     ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") rhsE,evn,fe00=",3e16.8)') i1,i2,m,rhsE,evn,fe00


     ! Update x derivative of P
     Pxxn  = beta * exxvn + rhsPxx

     ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") exxvn,Pxxn=",2e16.8)') i1,i2,m,exxvn,Pxxn


     ! Now we predict individual pk (pvn) and their third time derivative to second order
     PtttStar = 0

     do iv=0,numberOfPolarizationVectors-1
       pvn(iv) = betav(iv)*( .5*dt*a1v(iv)*evn + rhspv(iv) )
       ! Now update PtttStar for correction terms and use in Ettt
       ptttStarv(iv) = b1tttv(iv)*(pvn(iv)-pvm(iv))/(2.*dt) + b0tttv(iv)*pv(iv) + a0tttv(iv)*ev + a1tttv(iv)*(evn-evm)/(2.*dt) + a2tttv(iv)*(evn-2*ev+evm)/dtsq + f2v(iv)
       PtttStar = PtttStar + ptttStarv(iv)


       ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pvm,pv,pvn,evm,ev,evn=",6e16.8)') i1,i2,m,pvm(iv),pv(iv),pvn(iv),evm,ev,evn

       ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") evn-evm, evn-2*ev+evm=",2e16.8)') i1,i2,m,evn-evm,evn-2*ev+evm

       ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") 1diff,2diff,3diff,f2v(iv)=",4e16.8)') i1,i2,m,(pvn(iv)-pvm(iv))/(2.*dt),(evn-evm)/(2.*dt),(evn-2*ev+evm)/dtsq,f2v(iv)
       ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pvn(iv),ptttStarv,dtsq=",3e16.8)') i1,i2,m,pvn(iv),ptttStarv(iv),dtsq
     end do

     ! Second Order Updates Complete, now we construct necessary terms
     ! LapPtt using prediction
     QxxStar = (Pxxn  - pxxSum)/dtsq
     EtxxStar    = (exxvn -  exxvm)/(2.*dt)
     EtttStar    = csq*EtxxStar - alphaP*PtttStar + f1

     ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") QxxStar,exxvn,Pxxn,EtxxStar,PtttStar,EtttStar=",7e16.8)') i1,i2,m,QxxStar,exxvn,Pxxn,EtxxStar,PtttStar,EtttStar,f1

     rhsP4 = 0
     LHSev = 0

     do iv=0,numberOfPolarizationVectors-1
         ! Coeff for pk(n+1) for  LHS of invidiual p equation
         LHSpv(iv) = 1+ b1v(iv)*dt/2.  + b0v(iv)*dtsq/12.

         ! Build coeff for ev for left hand side of full P equation
         LHSev     = LHSev + ((-a1v(iv)*dt/2.  -a0v(iv)*dtsq/12.)/LHSpv(iv))

         ! LHS for pk for pk equation
         rhspv(iv) = ((2.*pv(iv)-pvm(iv))\
         + b1v(iv)*dt*pvm(iv)/(2.)\
         - dtsq*b0v(iv)*pv(iv)\
         + dtsq*a0v(iv)*ev\
         - a1v(iv)*(dt/2.)*evm\
         - a1v(iv)*(dt**4/12.)*( EtttStar )\
         + b1v(iv)*(dt**4/12.)*( ptttStarv(iv) )\
         + (dtsq/12.)*( -b0v(iv)*(0. -2.*pv(iv)+pvm(iv)) + a0v(iv)*(0. -2.*ev+evm) )\
         + f3v(iv))

        ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") LHSpv,rhspv=",2e16.8)') i1,i2,m,LHSpv(iv),rhspv(iv)


        rhsP4 = rhsP4 + ( rhspv(iv) )/LHSpv(iv)

     end do

     ! We have now built the equation for P
     A(1,1) = LHSev
     A(1,2) = 1.

     b(1) = rhsP4

     ! Now we build the equation for E in terms of E and P

     A(2,1) = 1.        ! coeff of E^{n+1}
     A(2,2) = alphaP    ! coeff of P^{n+1}

     ! Note that Psum - rhsP here is same as for second order code
     !  (but using 4th order accurate version of p to compute them)
     b(2) = (2.*ev-evm)\
            +csq*dtsq*elap4\
            +alphaP*( Psum )\
            +csq**2*dt**4*elapsq2/(12.)\
            -csq*dt**4*alphaP*QxxStar/(12.)\
            + f4


     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b2,f4=",2e16.8)') i1,i2,m,b(2),f4


     deti = 1./(A(1,1)*A(2,2) - A(1,2)*A(2,1))

     y(1) = deti *( A(2,2) * b(1) - A(1,2) * b(2))
     y(2) = deti *(-A(2,1) * b(1) + A(1,1) * b(2))

     ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b1,b2,y1,y2,deti=",5e16.8)') i1,i2,m,b(1),b(2),y(1),y(2),deti

     ! Update E^{n+1}
     un(i1,i2,i3,ec) = y(1)
     evn             = y(1)

     ! Update pk using new E^{n+1} = evn
     do iv=0,numberOfPolarizationVectors-1
       rhspv(iv) = rhspv(iv) + (a1v(iv) * dt/(2.) + a0v(iv)*dtsq/(12.))*evn
       pn(i1,i2,i3,m+iv*nd)  = (1/LHSpv(iv)) * rhspv(iv)
     end do



   end do ! m=0,1
  endLoopsMask()

 end if
#endMacro



! **********************************************************************************
! Macro updateDispersive 
!
! Initial version from Michael Jenkinson 2018
!
! **********************************************************************************
#beginMacro updateDispersive(DIM,ORDER,GRIDTYPE)

  if( t.le.3*dt )then
    INFO("update-dispersive_dim=DIM _order=ORDER _gridType=GRIDTYPE");
  end if

  fe=0.
  ! -- first compute some coefficients ---
  beta=0.
  do iv=0,numberOfPolarizationVectors-1
    betav(iv) = 1./( 1.+.5*dt*b1v(iv) )
    beta = beta + .5*dt*a1v(iv)*betav(iv)
    fpv(iv)=0.  ! initialize if not used
    b1tttv(iv)=b1v(iv)*b1v(iv)-b0v(iv)
    b0tttv(iv)=b1v(iv)*b0v(iv)
    a0tttv(iv)=-a0v(iv)*b1v(iv)
    a1tttv(iv)=a0v(iv)-a1v(iv)*b1v(iv)
    a2tttv(iv)=a1v(iv)
  end do

  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

  do m=0,nd-1
    pc=pxc+m
    ec=ex+m

    ! This is only needed for the second order code
    if( addForcing.ne.0 )then ! forcing in E equation already added to f
      ! wdh: Keep this term now: NOTE : fe is replaced for fourth-order below
      fe = dtsq*f(i1,i2,i3,ec)
      ! this next function will adjust fe by affing -alphaP*Ptt
      getGDMForcing(ec,pc)

      fp=fpv(0)
    end if

    ev = u(i1,i2,i3,ec)
    evm=um(i1,i2,i3,ec)

    do iv=0,numberOfPolarizationVectors-1
      pv(iv) = p(i1,i2,i3,m+iv*nd)
      pvm(iv)=pm(i1,i2,i3,m+iv*nd)
    end do

    rhsP = 0.
    pSum = 0.

    ! write(*,*) 'Inside updateDispersive'

    #If #ORDER eq "2"

        ! write(*,*) 'Inside updateDispersive order=2'

        #If #DIM eq "2"
            #If #GRIDTYPE eq "rectangular"

              ! INFO("FD22r-2D-dispersive-Any-PV");
              ! write(*,*) 'Inside updateDispersive rectangular order=2'

              elap2 = lap2d2(i1,i2,i3,ec)

            #Elif #GRIDTYPE eq "curvilinear"

              ! INFO("FD22c-2D-dispersive-any-PV");

              ! write(*,*) 'Inside updateDispersive curvilinear order=2'

              elap2 = ulaplacian22(i1,i2,i3,ec)

            #End


        #Elif #DIM eq "3"
            #If #GRIDTYPE eq "rectangular"

              ! INFO("FD22r-3D-dispersive-Any-PV");

              elap2 = lap3d2(i1,i2,i3,ec)

            #Elif #GRIDTYPE eq "curvilinear"

              ! INFO("FD22c-3D-dispersive-any-PV");

              elap2 = ulaplacian23(i1,i2,i3,ec)

            #End


        #Else
            ! Stop message
            stop 123
        #End

        do iv=0,numberOfPolarizationVectors-1

          rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + fpv(iv)
          rhsP = rhsP + betav(iv)*rhspv(iv)
          pSum = pSum + 2.*pv(iv) - pvm(iv)

          !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") two terms",2e16.8)') i1,i2,m,.5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm),dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev )

        end do

        rhsE = (2.*ev-evm) + csq*dtsq*elap2 + alphaP*( pSum - rhsP ) + fe

        evn = rhsE / (1.+ alphaP*beta)

        un(i1,i2,i3,ec) = evn

        do iv=0,numberOfPolarizationVectors-1
          pn(i1,i2,i3,m+iv*nd)  = betav(iv)*( .5*dt*a1v(iv)*evn + rhspv(iv) )
          ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pvn(iv)=",4e16.8)') i1,i2,m,pn(i1,i2,i3,m+iv*nd)
        end do

        ! End of second order code

    #Elif #ORDER eq "4"

        ! write(*,*) 'Inside updateDispersive order=4'

        #If #DIM eq "2"
            #If #GRIDTYPE eq "rectangular"

              ! INFO("FD44r-2D-dispersive-Any-PV");
              ! write(*,*) 'Inside updateDispersive 2D rectangular order=4'

              elap4   = lap2d4(i1,i2,i3,ec)
              elap4m  = lap2d4m(i1,i2,i3,ec)
              elapsq2 = lap2d2Pow2(i1,i2,i3,ec)
              elap2m  = lap2d2m(i1,i2,i3,ec)

              do iv=0,numberOfPolarizationVectors-1
                pxxv(iv)  = plap2d4(i1,i2,i3,m+iv*nd)
                pxxvm(iv) = plap2d4m(i1,i2,i3,m+iv*nd)
              end do

            #Elif #GRIDTYPE eq "curvilinear"

              ! INFO("FD44c-2D-dispersive-any-PV");

              ! write(*,*) 'Inside updateDispersive 2D curvilinear order=4'

              elap4   = ulaplacian42(i1,i2,i3,ec)
              elap4m  = umlaplacian42(i1,i2,i3,ec)
              elapsq2 = vlaplacian22(i1,i2,i3,ec)
              elap2m  = umlaplacian22(i1,i2,i3,ec)

              do iv=0,numberOfPolarizationVectors-1
                pxxv(iv)  = plaplacian42(i1,i2,i3,m+iv*nd)
                pxxvm(iv) = pmlaplacian42(i1,i2,i3,m+iv*nd)
              end do

            #End


        #Elif #DIM eq "3"
            #If #GRIDTYPE eq "rectangular"

              ! INFO("FD44r-3D-dispersive-Any-PV");
              ! write(*,*) 'Inside updateDispersive 3D rectangular order=4'

              elap4   = lap3d4(i1,i2,i3,ec)
              elap4m  = lap3d4m(i1,i2,i3,ec)
              elapsq2 = lap3d2Pow2(i1,i2,i3,ec)
              elap2m  = lap3d2m(i1,i2,i3,ec)

              do iv=0,numberOfPolarizationVectors-1
                pxxv(iv)  = plap3d4(i1,i2,i3,m+iv*nd)
                pxxvm(iv) = plap3d4m(i1,i2,i3,m+iv*nd)
              end do

            #Elif #GRIDTYPE eq "curvilinear"

              ! INFO("FD44c-3D-dispersive-any-PV");
              ! write(*,*) 'Inside updateDispersive 3D curvilinear order=4'

              elap4   = ulaplacian43(i1,i2,i3,ec)
              elap4m  = umlaplacian43(i1,i2,i3,ec)
              elapsq2 = vlaplacian23(i1,i2,i3,ec)
              elap2m  = umlaplacian23(i1,i2,i3,ec)

              do iv=0,numberOfPolarizationVectors-1
                pxxv(iv)  = plaplacian43(i1,i2,i3,m+iv*nd)
                pxxvm(iv) = pmlaplacian43(i1,i2,i3,m+iv*nd)
              end do

            #End
         #End

         ! Bug fixed, May 28, 2018 -- use 2D or 3D versions of ogderiv *wdh* 
         #If #DIM eq "2" 
           getGDMForcing44(ec,pc,OGDERIV2D)
         #Else
           getGDMForcing44(ec,pc,OGDERIV3D)
         #End

         rhsPxx = 0.
         pxxSum = 0.

         exxv  = elap4
         exxvm = elap4m

         ! First we do the second order prediction on (i) E, (ii) Exx, (iii) Pxx, (iv) Individual pk
         do iv=0,numberOfPolarizationVectors-1

           rhspv(iv) = 2.*pv(iv)-pvm(iv) + .5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ) + dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev ) + dtsq*fp00v(iv)
           rhsP = rhsP + betav(iv)*rhspv(iv)
           pSum = pSum + 2.*pv(iv) - pvm(iv)

           !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pv,pvm,rhspv,rhsP,pSum,fpv=",6e16.8)') i1,i2,m,pv(iv),pvm(iv),rhspv(iv),rhsP,pSum,fp00v(iv)
           !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") ev,evm,rhspv,rhsP,pSum,fpv=",6e16.8)') i1,i2,m,ev,evm,a1v(iv),a0v(iv),b1v(iv),b0v(iv)
           !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") first term=",1e16.8)') i1,i2,m,2.*pv(iv)-pvm(iv)
           !write(*,'(" (i2,i2,m)=(",i3,i3,i2,") two terms",2e16.8)') i1,i2,m,.5*dt*( b1v(iv)*pvm(iv) -a1v(iv)*evm ),dtSq*( -b0v(iv)*pv(iv) + a0v(iv)*ev )

           rhspxxv(iv) = 2.*pxxv(iv)-pxxvm(iv) + .5*dt*( b1v(iv)*pxxvm(iv) -a1v(iv)*exxvm ) + dtSq*( -b0v(iv)*pxxv(iv) + a0v(iv)*exxv ) + f6v(iv)
           rhsPxx = rhsPxx + betav(iv)*rhspxxv(iv)
           pxxSum = pxxSum + 2.*pxxv(iv) - pxxvm(iv)

         end do

         rhsE   = 2.*ev   - evm   + cdtsq*elap4   + alphaP*( pSum   - rhsP   ) + dtsq * fe00
         rhsExx = 2.*exxv - exxvm + cdtsq*elapsq2 + alphaP*( pxxSum - rhsPxx ) + f5

         evn   = rhsE   / (1.+ alphaP*beta)
         exxvn = rhsExx / (1.+ alphaP*beta)

         ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") rhsE,evn,fe00=",3e16.8)') i1,i2,m,rhsE,evn,fe00

         ! Update x derivative of P
         Pxxn  = beta * exxvn + rhsPxx

         ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") exxvn,Pxxn=",2e16.8)') i1,i2,m,exxvn,Pxxn

         ! Now we predict individual pk (pvn) and their third time derivative to second order
         PtttStar = 0

         do iv=0,numberOfPolarizationVectors-1
           pvn(iv) = betav(iv)*( .5*dt*a1v(iv)*evn + rhspv(iv) )
           ! Now update PtttStar for correction terms and use in Ettt
           ptttStarv(iv) = b1tttv(iv)*(pvn(iv)-pvm(iv))/(2.*dt) + b0tttv(iv)*pv(iv) + a0tttv(iv)*ev + a1tttv(iv)*(evn-evm)/(2.*dt) + a2tttv(iv)*(evn-2*ev+evm)/dtsq + f2v(iv)
           PtttStar = PtttStar + ptttStarv(iv)

           ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pvm,pv,pvn,evm,ev,evn=",6e16.8)') i1,i2,m,pvm(iv),pv(iv),pvn(iv),evm,ev,evn
           ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") evn-evm, evn-2*ev+evm=",2e16.8)') i1,i2,m,evn-evm,evn-2*ev+evm
           ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") 1diff,2diff,3diff,f2v(iv)=",4e16.8)') i1,i2,m,(pvn(iv)-pvm(iv))/(2.*dt),(evn-evm)/(2.*dt),(evn-2*ev+evm)/dtsq,f2v(iv)   
           ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") pvn(iv),ptttStarv,dtsq=",3e16.8)') i1,i2,m,pvn(iv),ptttStarv(iv),dtsq
         end do

         ! Second Order Updates Complete, now we construct necessary terms
         ! LapPtt using prediction
         QxxStar = (Pxxn  - pxxSum)/dtsq
         EtxxStar    = (exxvn -  exxvm)/(2.*dt)
         EtttStar    = csq*EtxxStar - alphaP*PtttStar + f1

         ! write(*,'(" (i2,i2,m)=(",i3,i3,i2,") QxxStar,exxvn,Pxxn,EtxxStar,PtttStar,EtttStar=",7e16.8)') i1,i2,m,QxxStar,exxvn,Pxxn,EtxxStar,PtttStar,EtttStar,f1

         rhsP4 = 0
         LHSev = 0

         do iv=0,numberOfPolarizationVectors-1
           ! Coeff for pk(n+1) for  LHS of invidiual p equation
           LHSpv(iv) = 1+ b1v(iv)*dt/2.  + b0v(iv)*dtsq/12.

           ! Build coeff for ev for left hand side of full P equation
           LHSev     = LHSev + ((-a1v(iv)*dt/2.  -a0v(iv)*dtsq/12.)/LHSpv(iv))

           ! LHS for pk for pk equation
           rhspv(iv) = ((2.*pv(iv)-pvm(iv))\
                       + b1v(iv)*dt*pvm(iv)/(2.)\
                       - dtsq*b0v(iv)*pv(iv)\
                       + dtsq*a0v(iv)*ev\
                       - a1v(iv)*(dt/2.)*evm\
                       - a1v(iv)*(dt**4/12.)*( EtttStar )\
                       + b1v(iv)*(dt**4/12.)*( ptttStarv(iv) )\
                       + (dtsq/12.)*( -b0v(iv)*(0. -2.*pv(iv)+pvm(iv)) + a0v(iv)*(0. -2.*ev+evm) )\
                       + f3v(iv))

           ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") LHSpv,rhspv=",2e16.8)') i1,i2,m,LHSpv(iv),rhspv(iv)

           rhsP4 = rhsP4 + ( rhspv(iv) )/LHSpv(iv)

         end do

         ! We have now built the equation for P
         A(1,1) = LHSev
         A(1,2) = 1.

         b(1) = rhsP4

         ! Now we build the equation for E in terms of E and P

         A(2,1) = 1.        ! coeff of E^{n+1}
         A(2,2) = alphaP    ! coeff of P^{n+1}

         ! Note that Psum - rhsP here is same as for second order code
         !  (but using 4th order accurate version of p to compute them)
         b(2) = (2.*ev-evm)\
                +csq*dtsq*elap4\
                +alphaP*( Psum )\
                +csq**2*dt**4*elapsq2/(12.)\
                -csq*dt**4*alphaP*QxxStar/(12.)\
                + f4


         ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b2,f4=",2e16.8)') i1,i2,m,b(2),f4


         deti = 1./(A(1,1)*A(2,2) - A(1,2)*A(2,1))

         y(1) = deti *( A(2,2) * b(1) - A(1,2) * b(2))
         y(2) = deti *(-A(2,1) * b(1) + A(1,1) * b(2))

         ! write(*,'(" (i1,i2,m)=(",i3,i3,i2,") b1,b2,y1,y2,deti=",5e16.8)') i1,i2,m,b(1),b(2),y(1),y(2),deti

         ! Update E^{n+1}
         un(i1,i2,i3,ec) = y(1)
         evn             = y(1)

         ! Update pk using new E^{n+1} = evn
         do iv=0,numberOfPolarizationVectors-1
           rhspv(iv) = rhspv(iv) + (a1v(iv) * dt/(2.) + a0v(iv)*dtsq/(12.))*evn
           pn(i1,i2,i3,m+iv*nd)  = (1/LHSpv(iv)) * rhspv(iv)

           ! write(*,'("advOpt: i1,i2=",2i3," f,fe,fp,pn=",4e12.4)') i1,i2,f(i1,i2,i3,ec),fe,fp,pn(i1,i2,i3,m+iv*nd)
         end do

         ! End of fourth order code

    #Else
            ! Stop message if order not 2 or 4
            stop 123
    #End


    ! End of fourth order code

  end do !m=0,nd-1

  endLoopsMask()

#endMacro

! **********************************************************************************
! Macro updateMultilevelAtomic (MLA)
!
! Initial version: July 2020
!
! **********************************************************************************
#beginMacro updateMultilevelAtomic(DIM,ORDER,GRIDTYPE)

  if( t.le.3*dt )then
    INFO("update-MULTI-LEVEL-ATOMIC_dim=DIM _order=ORDER _gridType=GRIDTYPE");
  end if
  
  ! initialize forcing functions
  fe=0.
  fet = 0.
  fett = 0.
  lapfe = 0.
  do iv=0,numberOfPolarizationVectors-1
    betav(iv) = 1./( 1.+.5*dt*b1v(iv) ) ! coefficients
    fpv(iv)=0.  ! initialize if not used
    fptv(iv)=0. 
    fpttv(iv)=0. 
    lapfpv(iv)=0.
  end do

  ! index location for first TZ nonlinear variable: 
  nce = pxc+nd*numberOfPolarizationVectors
  ! write(*,'(" *** UpadateMLA: pxc=",i2," numberOfPolarizationVectors=",i4," nce=",i4)') pxc,numberOfPolarizationVectors,nce 

  ! loop over space
  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  #If #ORDER eq "2"
    do m=0,nd-1
      pc=pxc+m
      ec=ex+m

      if( addForcing.ne.0 )then
        getEPMLAForcing(ec,pc) ! forcing functions for E and P
      end if

      ! ping old values
      ev = u(i1,i2,i3,ec)
      evm=um(i1,i2,i3,ec)
      do iv=0,numberOfPolarizationVectors-1
        pv(iv) = p(i1,i2,i3,m+iv*nd)
        pvm(iv)=pm(i1,i2,i3,m+iv*nd)
      end do

      ! write(*,*) 'Inside updateDispersive'

      ! write(*,*) 'Inside updateDispersive order=2'
      
      ! get laplacians
      #If #DIM eq "2"
          #If #GRIDTYPE eq "rectangular"

            ! INFO("FD22r-2D-dispersive-Any-PV");
            ! write(*,*) 'Inside updateDispersive rectangular order=2'

            ! elap2 = ulaplacian22r(i1,i2,i3,ec)
            elap2 = lap2d2(i1,i2,i3,ec)

          #Elif #GRIDTYPE eq "curvilinear"

            ! INFO("FD22c-2D-dispersive-any-PV");

            ! write(*,*) 'Inside updateDispersive curvilinear order=2'

            elap2 = ulaplacian22(i1,i2,i3,ec)

          #End


      #Elif #DIM eq "3"
          #If #GRIDTYPE eq "rectangular"

            ! INFO("FD22r-3D-dispersive-Any-PV");

            elap2 = lap3d2(i1,i2,i3,ec)

          #Elif #GRIDTYPE eq "curvilinear"

            ! INFO("FD22c-3D-dispersive-any-PV");

            elap2 = ulaplacian23(i1,i2,i3,ec)

          #End

      #Else
          ! Stop message
          stop 123
      #End

      ! second order update of P_m
      pSum = 0.
      do iv=0,numberOfPolarizationVectors-1

        pvn(iv) = 2.*pv(iv)-pvm(iv) + 0.5*dt*b1v(iv)*pvm(iv) - dtsq*b0v(iv)*pv(iv) + dtsq*fpv(iv)

        do na = 0,numberOfAtomicLevels-1 ! \Delta N^n*E^n
          pvn(iv) = pvn(iv) + dtsq*pnec(iv,na)*q(i1,i2,i3,na)*ev
        enddo

        pn(i1,i2,i3,m+iv*nd) = betav(iv)*pvn(iv)

        pSum = pSum + betav(iv)*pvn(iv) -2.*pv(iv) + pvm(iv)

      end do
      
      ! second order update of E
      evn = (2.*ev-evm) + csq*dtsq*elap2 - alphaP*pSum + dtsq*fe
      un(i1,i2,i3,ec) = evn

      ! End of second order code
    end do !m=0,nd-1 over space dim

    ! outside of dimension loop
    ! E(t+dt) and P(t+dt) are now known 
    ! --- second order update of N ---
    ! MLA
    do na=0,numberOfAtomicLevels-1
      ! forcing function
      getMLAForcing(na)
    enddo
    ! N_t
    do na=0,numberOfAtomicLevels-1
      ! first time derivative
      qt(na) = fnv(na)
      do iv=0,numberOfAtomicLevels-1 ! relaxation
        qt(na) = qt(na)+prc(na,iv)*q(i1,i2,i3,iv) 
      enddo
      do m=0,nd-1 ! dot product
        do iv = 0,numberOfPolarizationVectors-1 ! loop over pc
          qt(na) = qt(na)+peptc(na,iv)*u(i1,i2,i3,m)*(pn(i1,i2,i3,m+iv*nd)-pm(i1,i2,i3,m+iv*nd))/(2.*dt)
        enddo
      enddo
    enddo
    ! N_tt
    do na=0,numberOfAtomicLevels-1
      qtt(na) = fntv(na)
      do iv=0,numberOfAtomicLevels-1 ! relaxation
        qtt(na) = qtt(na)+prc(na,iv)*qt(iv)
      enddo
      do m=0,nd-1 ! dot product
        do iv = 0,numberOfPolarizationVectors-1 ! loop over pc
          qtt(na) = qtt(na)+peptc(na,iv)*(un(i1,i2,i3,m)-um(i1,i2,i3,m))/(2.*dt)*(pn(i1,i2,i3,m+iv*nd)-pm(i1,i2,i3,m+iv*nd))/(2.0*dt) \
                           +peptc(na,iv)*u(i1,i2,i3,m)*(pn(i1,i2,i3,m+iv*nd)-2.0*p(i1,i2,i3,m+iv*nd)+pm(i1,i2,i3,m+iv*nd))/dtsq
        enddo
      enddo
    enddo

    ! fill in the population densities "N"
    do na=0,numberOfAtomicLevels-1
      qn(i1,i2,i3,na) = q(i1,i2,i3,na) + dt*qt(na) + dt**2/2.*qtt(na)
    end do

  #Elif #ORDER eq '4'
    !
    ! second order predictions
    !
    do m=0,nd-1
      pc=pxc+m
      ec=ex+m

      if ( addForcing.ne.0) then 
        getEPMLAForcing(ec,pc) ! forcing functions for E and P
      end if

      ! ping old values
      ev = u(i1,i2,i3,ec)
      evm=um(i1,i2,i3,ec)
      do iv=0,numberOfPolarizationVectors-1
        pv(iv) = p(i1,i2,i3,m+iv*nd)
        pvm(iv)=pm(i1,i2,i3,m+iv*nd)
      end do

      ! get laplacians
      #If #DIM eq "2"
          #If #GRIDTYPE eq "rectangular"

            ! INFO("FD22r-2D-dispersive-Any-PV");
            ! write(*,*) 'Inside updateDispersive rectangular order=2'

            elap2 = lap2d2(i1,i2,i3,ec)

          #Elif #GRIDTYPE eq "curvilinear"

            ! INFO("FD22c-2D-dispersive-any-PV");
            ! write(*,*) 'Inside updateDispersive curvilinear order=2'

            elap2 = ulaplacian22(i1,i2,i3,ec)

          #End

      #Elif #DIM eq "3"
          #If #GRIDTYPE eq "rectangular"

            ! INFO("FD22r-3D-dispersive-Any-PV");

            elap2 = lap3d2(i1,i2,i3,ec)

          #Elif #GRIDTYPE eq "curvilinear"

            ! INFO("FD22c-3D-dispersive-any-PV");

            elap2 = ulaplacian23(i1,i2,i3,ec)

          #End

      #Else
          ! Stop message
          stop 123
      #End

      ! second order update of P_m
      pSum = 0.
      do iv=0,numberOfPolarizationVectors-1

        pvn(iv) = 2.*pv(iv)-pvm(iv) + 0.5*dt*b1v(iv)*pvm(iv) - dtsq*b0v(iv)*pv(iv) + dtsq*fpv(iv)

        do na = 0,numberOfAtomicLevels-1 ! \Delta N^n*E^n
          pvn(iv) = pvn(iv) + dtsq*pnec(iv,na)*q(i1,i2,i3,na)*ev
        enddo

        pn(i1,i2,i3,m+iv*nd) = betav(iv)*pvn(iv)

        pSum = pSum + betav(iv)*pvn(iv) -2.*pv(iv) + pvm(iv)

      end do
      
      ! second order update of E
      evn = (2.*ev-evm) + csq*dtsq*elap2 - alphaP*pSum + dtsq*fe
      un(i1,i2,i3,ec) = evn

      ! End of second order code
    end do ! m=0,nd-1 over space dim

    ! outside of dimension loop
    ! --- second order update of N ---
    ! MLA
    do na=0,numberOfAtomicLevels-1
      ! forcing functions for both the second order codes and fourth order codes
      getMLAForcing44(na)
    enddo
    ! N_t
    do na=0,numberOfAtomicLevels-1
      qt(na) = fnv(na)
      do iv=0,numberOfAtomicLevels-1 ! relaxation
        qt(na) = qt(na)+prc(na,iv)*q(i1,i2,i3,iv) 
      enddo
      do m=0,nd-1 ! dot product
        do iv = 0,numberOfPolarizationVectors-1 ! loop over pc
          qt(na) = qt(na)+peptc(na,iv)*u(i1,i2,i3,m)*(pn(i1,i2,i3,m+iv*nd)-pm(i1,i2,i3,m+iv*nd))/(2.*dt)
        enddo
      enddo
    enddo
    ! N_tt
    do na=0,numberOfAtomicLevels-1
      qtt(na) = fntv(na)
      do iv=0,numberOfAtomicLevels-1 ! relaxation
        qtt(na) = qtt(na)+prc(na,iv)*qt(iv)
      enddo
      do m=0,nd-1 ! dot product
        do iv = 0,numberOfPolarizationVectors-1 ! loop over pc
          qtt(na) = qtt(na)+peptc(na,iv)*(un(i1,i2,i3,m)-um(i1,i2,i3,m))/(2.*dt)*(pn(i1,i2,i3,m+iv*nd)-pm(i1,i2,i3,m+iv*nd))/(2.0*dt) \
                           +peptc(na,iv)*u(i1,i2,i3,m)*(pn(i1,i2,i3,m+iv*nd)-2.0*p(i1,i2,i3,m+iv*nd)+pm(i1,i2,i3,m+iv*nd))/dtsq
        enddo
      enddo
    enddo
    ! fill in the population densities "N"
    do na=0,numberOfAtomicLevels-1
      qn(i1,i2,i3,na) = q(i1,i2,i3,na) + dt*qt(na) + dt**2/2.*qtt(na)
    end do

    !----------------------
    ! fourth order update
    !----------------------
    ! write(*,*) 'Inside updateMultilevelAtomic 2D rectangular order=4'
    do m=0,nd-1 ! loop over dimension
      pc=pxc+m
      ec=ex+m

      if( addForcing.ne.0 )then
          getEPMLAForcing44(ec,pc)
      end if

      ! ping the current indexed values of E and P
      evn=un(i1,i2,i3,ec)
      ev = u(i1,i2,i3,ec)
      evm=um(i1,i2,i3,ec)

      do iv=0,numberOfPolarizationVectors-1
        pvn(iv)=pn(i1,i2,i3,m+iv*nd)
        pv(iv) = p(i1,i2,i3,m+iv*nd)
        pvm(iv)=pm(i1,i2,i3,m+iv*nd)
      end do

      ! write(*,*) 'Inside updateDispersive order=4'

      ! space derivatives of E, P, N
      #If #DIM eq "2"
        #If #GRIDTYPE eq "rectangular"

          ! INFO("FD44r-2D-dispersive-Any-PV");
          ! write(*,*) 'Inside updateDispersive 2D rectangular order=4'

          ! laplacians of E
          elap4   = lap2d4(i1,i2,i3,ec)
          elapsq2 = lap2d2Pow2(i1,i2,i3,ec)
          elap2   = lap2d2(i1,i2,i3,ec)
          elap2m  = lap2d2m(i1,i2,i3,ec)

          ! first order derivatives of E
          ex2v = ex2(i1,i2,i3,ec)
          ey2v = ey2(i1,i2,i3,ec)

          ! laplacians of P
          do iv=0,numberOfPolarizationVectors-1
            pxxv(iv)  = plap2d4(i1,i2,i3,m+iv*nd)
            pxxvm(iv) = plap2d4m(i1,i2,i3,m+iv*nd)
          end do

          
          do na=0,numberOfAtomicLevels-1
            qlap2(na) = qlap2d2(i1,i2,i3,na) ! laplacian of N
            qx2v(na) =  qx2(i1,i2,i3,na) ! first order derivatives of N
            qy2v(na) =  qy2(i1,i2,i3,na)
            ! chain rule for \Delta (E*N)
            qelap2(na) = ev*qlap2(na)+q(i1,i2,i3,na)*elap2+2.*ex2v*qx2v(na)+2.*ey2v*qy2v(na)
          enddo

        #Elif #GRIDTYPE eq "curvilinear"

          ! INFO("FD44c-2D-dispersive-any-PV");

          ! write(*,*) 'Inside updateDispersive 2D curvilinear order=4'

          ! laplacians of E
          elap4   = ulaplacian42(i1,i2,i3,ec)
          elapsq2 = vlaplacian22(i1,i2,i3,ec)
          elap2   = ulaplacian22(i1,i2,i3,ec)
          elap2m  = umlaplacian22(i1,i2,i3,ec)

          ! first order derivatives of E
          ex2v = ux22(i1,i2,i3,ec)
          ey2v = uy22(i1,i2,i3,ec)

          ! laplacians of P
          do iv=0,numberOfPolarizationVectors-1
            pxxv(iv)  = plaplacian42(i1,i2,i3,m+iv*nd)
            pxxvm(iv) = pmlaplacian42(i1,i2,i3,m+iv*nd)
          end do

          do na=0,numberOfAtomicLevels-1
            qlap2(na) = qlaplacian22(i1,i2,i3,na)
            qx2v(na) =  qx22(i1,i2,i3,na)
            qy2v(na) =  qy22(i1,i2,i3,na)
            ! chain rule for \Delta (E*N)
            qelap2(na) = ev*qlap2(na)+q(i1,i2,i3,na)*elap2+2.*ex2v*qx2v(na)+2.*ey2v*qy2v(na)
          enddo

        #End


      #Elif #DIM eq "3"
        #If #GRIDTYPE eq "rectangular"

          ! INFO("FD44r-3D-dispersive-Any-PV");
          ! write(*,*) 'Inside updateDispersive 3D rectangular order=4'
          
          ! laplacians of E
          elap4   = lap3d4(i1,i2,i3,ec)
          elapsq2 = lap3d2Pow2(i1,i2,i3,ec)
          elap2   = lap3d2(i1,i2,i3,ec)
          elap2m  = lap3d2m(i1,i2,i3,ec)

          ! first order derivatives of E
          ex2v = ex2(i1,i2,i3,ec)
          ey2v = ey2(i1,i2,i3,ec)
          ez2v = ez2(i1,i2,i3,ec)

          ! laplacians of P
          do iv=0,numberOfPolarizationVectors-1
            pxxv(iv)  = plap3d4(i1,i2,i3,m+iv*nd)
            pxxvm(iv) = plap3d4m(i1,i2,i3,m+iv*nd)
          end do

          do na=0,numberOfAtomicLevels-1
            qlap2(na) = qlap3d2(i1,i2,i3,na)
            qx2v(na) =  qx2(i1,i2,i3,na)
            qy2v(na) =  qy2(i1,i2,i3,na)
            qz2v(na) =  qz2(i1,i2,i3,na)
            ! chain rule for \Delta (E*N)
            qelap2(na) = ev*qlap2(na)+q(i1,i2,i3,na)*elap2+2.*ex2v*qx2v(na)+2.*ey2v*qy2v(na)+2.*ez2v*qz2v(na)
          enddo

        #Elif #GRIDTYPE eq "curvilinear"

          ! INFO("FD44c-3D-dispersive-any-PV");
          ! write(*,*) 'Inside updateDispersive 3D curvilinear order=4'

          ! laplacians of E
          elap4   = ulaplacian43(i1,i2,i3,ec)
          elapsq2 = vlaplacian23(i1,i2,i3,ec)
          elap2   = ulaplacian23(i1,i2,i3,ec)
          elap2m  = umlaplacian23(i1,i2,i3,ec)

          ! first order derivatives of E
          ex2v = ux23(i1,i2,i3,ec)
          ey2v = uy23(i1,i2,i3,ec)
          ez2v = uz23(i1,i2,i3,ec)

          ! laplacians of P
          do iv=0,numberOfPolarizationVectors-1
            pxxv(iv)   = plaplacian43(i1,i2,i3,m+iv*nd)
            pxxvm(iv) = pmlaplacian43(i1,i2,i3,m+iv*nd)
          end do

          do na=0,numberOfAtomicLevels-1
            qlap2(na) = qlaplacian23(i1,i2,i3,na)
            qx2v(na) =  qx23(i1,i2,i3,na)
            qy2v(na) =  qy23(i1,i2,i3,na)
            qz2v(na) =  qz23(i1,i2,i3,na)
            ! chain rule for \Delta (E*N)
            qelap2(na) = ev*qlap2(na)+q(i1,i2,i3,na)*elap2+2.*ex2v*qx2v(na)+2.*ey2v*qy2v(na)+2.*ez2v*qz2v(na)
          enddo

        #End
      #End

      !--------------------
      ! fourth order P
      !--------------------
      ! second order accurate terms
      et = (evn-evm)/(2.0*dt)
      ett = (evn-2.*ev+evm)/dtsq
      ptttSum = 0.
      do iv = 0,numberOfPolarizationVectors-1
        ptv(m,iv) = (pvn(iv)-pvm(iv))/(2.0*dt)
        pttv(m,iv) = (pvn(iv)-2.0*pv(iv)+pvm(iv))/dtsq
        ptttv(m,iv) = -b1v(iv)*pttv(m,iv)-b0v(iv)*ptv(m,iv)+fptv(iv)
        do na = 0,numberOfAtomicLevels-1 ! update using ODE
          ptttv(m,iv) = ptttv(m,iv) + pnec(iv,na)*qt(na)*ev \
                                    + pnec(iv,na)*q(i1,i2,i3,na)*et
        enddo
        pttttv(m,iv) = -b1v(iv)*ptttv(m,iv)-b0v(iv)*pttv(m,iv)+fpttv(iv) ! update using ODE
        do na = 0,numberOfAtomicLevels-1
          pttttv(m,iv) = pttttv(m,iv) + pnec(iv,na)*qtt(na)*ev \
                                   + 2.*pnec(iv,na)*qt(na)*et \
                                      + pnec(iv,na)*q(i1,i2,i3,na)*ett
        enddo
        ptttSum = ptttSum + ptttv(m,iv)
      enddo

      ! update 4th order accurate P
      pSum = 0.
      pxxSum = 0.
      do iv=0,numberOfPolarizationVectors-1

        pvn(iv) = 2.*pv(iv)-pvm(iv) + 0.5*dt*b1v(iv)*pvm(iv) + dt**4/12.*pttttv(m,iv) + dt**4/6.*b1v(iv)*ptttv(m,iv) - dtsq*b0v(iv)*pv(iv) + dtsq*fpv(iv)

        do na = 0,numberOfAtomicLevels-1 ! \Delta N^n*E^n
          pvn(iv) = pvn(iv) + dtsq*pnec(iv,na)*q(i1,i2,i3,na)*ev
        enddo

        pvn(iv) = betav(iv)*pvn(iv)

        pn(i1,i2,i3,m+iv*nd) = pvn(iv) ! store updated P

        pSum = pSum + pvn(iv) - 2.*pv(iv) + pvm(iv) ! sum P_tt

        ! laplacians of P^{n+1}
        pxxvn(iv) = 2.*pxxv(iv)-pxxvm(iv) + 0.5*dt*b1v(iv)*pxxvm(iv) - dtsq*b0v(iv)*pxxv(iv) + dtsq*lapfpv(iv)

        do na = 0,numberOfAtomicLevels-1 ! \Delta N^n*E^n
          pxxvn(iv) = pxxvn(iv) + dtsq*pnec(iv,na)*qelap2(na)
        enddo
        ! time derivatives of laplacian of P
        pxxSum = pxxSum + betav(iv)*pxxvn(iv) - 2.*pxxv(iv) + pxxvm(iv)

      end do

      !-------------------
      ! fourth order E
      !-------------------
      
      evn =   (2.*ev-evm) + csq*dtsq*elap4 - alphaP*pSum + dtsq*fe \
            + dt**4/12.*(csq**2*elapsq2 - alphaP*csq*pxxSum/dtsq + csq*lapfe + fett)
      
      ! store 4th order accurate E
      un(i1,i2,i3,ec) = evn

      ! laplacian of E at new time
      elap2n = 2.*elap2-elap2m+dtsq*csq*elapsq2-alphaP*pxxSum + dtsq*lapfe
      ettt = (csq*elap2n-csq*elap2m)/(2.0*dt)-alphaP*ptttSum + fet ! E_ttt
      
      ! store in (Ex,Ey,E_z)
      ettv(m) = (evn-2.*ev+evm)/dtsq
      etttv(m) = ettt

      ! fourth order accurate terms
      etv(m) = (evn-evm)/(2.*dt)-dtsq/6.*ettt
      do iv = 0,numberOfPolarizationVectors-1
        ptv(m,iv) = (pvn(iv)-pvm(iv))/(2.*dt)-dtsq/6.*ptttv(m,iv)
        pttv(m,iv) = -b1v(iv)*ptv(m,iv)-b0v(iv)*pv(iv) + fpv(iv)
        do na = 0,numberOfAtomicLevels-1
          pttv(m,iv) = pttv(m,iv) + pnec(iv,na)*q(i1,i2,i3,na)*ev
        enddo
        ! second order accurate terms
        ptttv(m,iv) = -b1v(iv)*pttv(m,iv)-b0v(iv)*ptv(m,iv)+fptv(iv)
        do na = 0,numberOfAtomicLevels-1
          ptttv(m,iv) = ptttv(m,iv) + pnec(iv,na)*qt(na)*ev \
                                     +pnec(iv,na)*q(i1,i2,i3,na)*etv(m)
        enddo
        pttttv(m,iv) = -b1v(iv)*ptttv(m,iv)-b0v(iv)*pttv(m,iv)+fpttv(iv)
        do na = 0,numberOfAtomicLevels-1
          pttttv(m,iv) = pttttv(m,iv) + pnec(iv,na)*qtt(na)*ev \
                                   +2.0*pnec(iv,na)*qt(na)*etv(m) \
                                      + pnec(iv,na)*q(i1,i2,i3,na)*ettv(m)
        enddo
      enddo 

    enddo ! m=0,nd-1
    
    !-------------------
    ! fourth order N
    !-------------------
    ! fourth order accurate terms
    ! N_t
    do na=0,numberOfAtomicLevels-1
      qt(na) = fnv(na)
      do iv=0,numberOfAtomicLevels-1 ! relaxation
        qt(na) = qt(na)+prc(na,iv)*q(i1,i2,i3,iv) 
      enddo
      do m=0,nd-1 ! dot product
        do iv = 0,numberOfPolarizationVectors-1 ! loop over pc
          qt(na) = qt(na)+peptc(na,iv)*u(i1,i2,i3,m)*ptv(m,iv)
        enddo
      enddo
    enddo
    ! N_tt
    do na=0,numberOfAtomicLevels-1
      qtt(na) = fntv(na)
      do iv=0,numberOfAtomicLevels-1 ! relaxation
        qtt(na) = qtt(na)+prc(na,iv)*qt(iv)
      enddo
      do m=0,nd-1 ! dot product
        do iv = 0,numberOfPolarizationVectors-1 ! loop over pc
          qtt(na) = qtt(na)+peptc(na,iv)*etv(m)*ptv(m,iv) \
                           +peptc(na,iv)*u(i1,i2,i3,m)*pttv(m,iv)
        enddo
      enddo
    enddo

    ! N_ttt
    do na=0,numberOfAtomicLevels-1
      qttt(na) = fnttv(na)
      do iv=0,numberOfAtomicLevels-1 ! relaxation
        qttt(na) = qttt(na)+prc(na,iv)*qtt(iv)
      enddo
      do m=0,nd-1 ! dot product
        do iv = 0,numberOfPolarizationVectors-1 ! loop over pc
          qttt(na) = qttt(na)+peptc(na,iv)*ettv(m)*ptv(m,iv) \
                          +2.*peptc(na,iv)*etv(m)*pttv(m,iv) \
                             +peptc(na,iv)*u(i1,i2,i3,m)*ptttv(m,iv)
        enddo
      enddo
    enddo

    ! N_tttt
    do na=0,numberOfAtomicLevels-1
      qtttt(na) = fntttv(na)
      do iv=0,numberOfAtomicLevels-1 ! relaxation
        qtttt(na) = qtttt(na)+prc(na,iv)*qttt(iv)
      enddo
      do m=0,nd-1 ! dot product
        do iv = 0,numberOfPolarizationVectors-1 ! loop over pc
          qtttt(na) = qtttt(na)+peptc(na,iv)*etttv(m)*ptv(m,iv) \
                            +3.*peptc(na,iv)*ettv(m)*pttv(m,iv) \
                            +3.*peptc(na,iv)*etv(m)*ptttv(m,iv) \
                               +peptc(na,iv)*u(i1,i2,i3,m)*pttttv(m,iv)
        enddo
      enddo
    enddo

    do na=0,numberOfAtomicLevels-1
      qn(i1,i2,i3,na) = q(i1,i2,i3,na) + dt*qt(na) + dt**2/2.*qtt(na) + dt**3/6.*qttt(na) + dt**4/24.*qtttt(na)
    end do

  #Else
    ! Stop message if order not 2 or 4
    stop 123
  #End

  endLoopsMask()

#endMacro

! **********************************************************************************
! Macro ADV_MAXWELL:
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
!  GRIDTYPE : rectangular, curvilinear
! **********************************************************************************
#beginMacro ADV_MAXWELL(NAME,DIM,ORDER,GRIDTYPE)
 subroutine NAME(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                 mask,xy,rsxy,  um,u,un,f,fa, v, pm,p,pn, qm, q,qn, bc, dis, varDis, ipar, rpar, ierr )
!======================================================================
!   Advance a time step for Maxwells equations
!     OPTIMIZED version for rectangular grids.
! nd : number of space dimensions
!
! ipar(0)  = option : option=0 - Maxwell+Artificial diffusion
!                           =1 - AD only
!
!  dis(i1,i2,i3) : temp space to hold artificial dissipation
!  varDis(i1,i2,i3) : coefficient of the variable artificial dissipation
!======================================================================
 implicit none
 integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b

 real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
 real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 ! Polarization vectors
 real pm(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
 real p(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
 real pn(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
 ! Vectors for the nonlinear multilevel Atomic model 

 real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)

 ! for nonlinear model
 real qm(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 real q(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 real qn(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 real qe(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,nd4a:nd4b)

 real dis(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 real varDis(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
 real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)

 integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
 integer bc(0:1,0:2),ierr

 integer ipar(0:*)
 real rpar(0:*)

!     ---- local variables -----
 integer m1a,m1b,m2a,m2b,m3a,m3b,numGhost,nStart,nEnd,debug

 integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderInTime,axis,dir
 integer addForcing,orderOfDissipation,option
 integer useWhereMask,useWhereMaskSave,solveForE,solveForH,grid,useVariableDissipation
 integer useCurvilinearOpt,useConservative,combineDissipationWithAdvance,useDivergenceCleaning
 integer useNewForcingMethod,numberOfForcingFunctions,fcur,fnext,fprev
 integer ex,ey,ez, hx,hy,hz
 real t,cc,dt,dy,dz,cdt,cdtdx,cdtdy,cdtdz,adc,adcdt,add,adddt
 real dt4by12
 real eps,mu,sigmaE,sigmaH,kx,ky,kz,divergenceCleaningCoefficient
 logical addDissipation

 real ep ! holds the pointer to the TZ function

 real dx(0:2),dr(0:2)

 real dx2i,dy2i,dz2i,dxsqi,dysqi,dzsqi,dxi,dyi,dzi
 real dx12i,dy12i,dz12i,dxsq12i,dysq12i,dzsq12i,dxy4i,dxz4i,dyz4,time0,time1

 real dxi4,dyi4,dzi4,dxdyi2,dxdzi2,dydzi2

 real uLap(-1:1,-1:1,0:5),uLapSq(0:5)
 real uLaprr2,uLapss2,uLaprs2,uLapr2,uLaps2

 real c0,c1,csq,dtsq,cdtsq,cdtsq12,lap(0:20),cdtSqBy12
 real c40,c41,c42,c43
 real c60,c61,c62,c63,c64,c65
 real c80,c81,c82,c83,c84,c85,c86,c87

 real c00lap2d6,c10lap2d6,c01lap2d6,c20lap2d6,c02lap2d6,c30lap2d6,c03lap2d6
 real c00lap2d8,c10lap2d8,c01lap2d8,c20lap2d8,c02lap2d8,c30lap2d8,c03lap2d8,c40lap2d8,c04lap2d8
 real c000lap3d6,c100lap3d6,c010lap3d6,c001lap3d6,\
                 c200lap3d6,c020lap3d6,c002lap3d6,\
                 c300lap3d6,c030lap3d6,c003lap3d6
 real c000lap3d8,c100lap3d8,c010lap3d8,c001lap3d8,\
                 c200lap3d8,c020lap3d8,c002lap3d8,\
                 c300lap3d8,c030lap3d8,c003lap3d8,\
                 c400lap3d8,c040lap3d8,c004lap3d8

 integer rectangular,curvilinear
 parameter( rectangular=0, curvilinear=1 )

 integer timeSteppingMethod
 integer defaultTimeStepping,adamsSymmetricOrder3,rungeKuttaFourthOrder,\
         stoermerTimeStepping,modifiedEquationTimeStepping
 parameter(defaultTimeStepping=0,adamsSymmetricOrder3=1,\
           rungeKuttaFourthOrder=2,stoermerTimeStepping=3,modifiedEquationTimeStepping=4)

 ! Dispersion models
      #Include "dispersionModelsFortranInclude.h"

 ! Nonlinear models
      #Include "nonlinearModelsFortranInclude.h"
      integer nonlinearModel
           
 ! forcing options
      #Include "forcingDefineFortranInclude.h"

!...........start statement function
 integer kd,m
 real rx,ry,rz,sx,sy,sz,tx,ty,tz

 declareDifferenceOrder2(u,RX)
 declareDifferenceOrder2(un,none)
 declareDifferenceOrder2(v,none)
 declareDifferenceOrder2(um,none)
 declareDifferenceOrder2(ff,none)
 declareDifferenceOrder2(p,none)
 declareDifferenceOrder2(pm,none)
 ! MLA
 declareDifferenceOrder2(q,none)

 declareDifferenceOrder4(u,RX)
 declareDifferenceOrder4(un,none)
 declareDifferenceOrder4(um,none)
 declareDifferenceOrder4(v,none)
 declareDifferenceOrder4(p,none)
 declareDifferenceOrder4(pm,none)

 real maxwell2dr,maxwell3dr,maxwellr44,maxwellr66,maxwellr88
 real maxwellc22,maxwellc44,maxwellc66,maxwellc88, maxwellc23
 real maxwell2dr44me,maxwell2dr66me,maxwell2dr88me
 real maxwell3dr44me,maxwell3dr66me,maxwell3dr88me
 real maxwellc44me,maxwellc66me,maxwellc88me
 real max2dc44me,max2dc44me2,max3dc44me
 real mxdc2d2Ex,mxdc2d2Ey,mxdc2d4Ex,mxdc2d4Ey, mxdc2d4cEx,mxdc2d4cEy
 real mxdc2d2cEx,mxdc2d2cEy
 real mxdc3d2Ex,mxdc3d2Ey,mxdc3d2Ez,mxdc3d2Hx,mxdc3d2Hy,mxdc3d2Hz
 real mxdc3d2cEx,mxdc3d2cEy,mxdc3d2cEz,mxdc3d2cHx,mxdc3d2cHy,mxdc3d2cHz
 real mxdc2d4cConsEx,mxdc2d4cConsEy,mxdc2d4cConsEz
 real mxdc3d4Ex,mxdc3d4Ey,mxdc3d4Ez,mxdc3d4Hx,mxdc3d4Hy,mxdc3d4Hz

 real DptU,DmtU,DztU, DzstU

 real fhz
 real hz0t,hz0x,hz0y
 real ex0,ex0t,ex0x,ex0y,ex0z
 real ey0,ey0t,ey0x,ey0y,ey0z
 real ez0,ez0t,ez0x,ez0y,ez0z

 real p0,p0t,p0tt
 real e0,e0t,e0tt

 real cdt4by360,cdt6by20160

 real lap2d2,lap3d2,lap2d4,lap3d4,lap2d6,lap3d6,lap2d8,lap3d8,lap2d2Pow2,lap3d2Pow2,lap2d2Pow3,lap3d2Pow3,\
      lap2d2Pow4,lap3d2Pow4,lap2d4Pow2,lap3d4Pow2,lap2d4Pow3,lap3d4Pow3,lap2d6Pow2,lap3d6Pow2
 real lap2d2m,lap3d2m

 real ex2,ey2,ez2,qx2,qy2,qz2,qlap2d2,qlap3d2,ex2v,ey2v,ez2v

 real du,fd22d,fd23d,fd42d,fd43d,fd62d,fd63d,fd82d,fd83d
 real elap4, elap4m, lap2d4m, lap3d4m, elapsq2, elap2n, elap2, elap2m, plap2, plap2m, plap4, plap4m, plap2d2, plap2d2m, plap2d4, plap2d4m

 real a0ttt, a1ttt, a2ttt, b0ttt, b1ttt, Exx, Exxxx, En, Enm1, Pnm1, Exxn, Exxnm1, betaP, Pxxn, Pxxnm1, Etxx,\
          cStar, PtttStarRhs, rhsExxnp1Predict, rhsPxxnp1Predict, Exxnp1Predict, Pxxnp1Predict, A(2,2), b(2), y(2),\
          f1, f2, f3, f4, f5, f6

 real p0ttt, p0tttt, p0xx, p0yy, p0xxt, p0yyt, p0xxtt, p0yytt, e0ttt, e0tttt,\
            e0x,e0y,e0z, \
            e0xx, e0yy, e0xxt, e0yyt, e0xxtt, e0yytt, e0xxxx, e0xxyy, e0yyyy
 real fp00, fp10, fp20, fp02x, fp02y, fp02, fe00, fe10, fe20, fe02x, fe02y, fe02
 real p0zz, p0zzt, p0zztt, e0zz, e0zzt, e0zztt, e0xxzz, e0yyzz, e0zzzz, fp02z, fe02z
 real plap3d2, plap3d2m, plap3d4, plap3d4m
 real q0x,q0y,q0z,q0xx,q0yy,q0zz

 ! forcing correction functions:
 real lap2d2f,f2drme44, lap3d2f, f3drme44, f2dcme44, f3dcme44, ff

 real cdSosupx,cdSosupy,cdSosupz, adSosup,sosupParameter, uDotFactor, adxSosup(0:2)
 integer useSosupDissipation,sosupDissipationOption
 integer updateSolution,updateDissipation,computeUt,forcingOption

 ! div cleaning:
 real dc,dcp,cdc0,cdc1,cdcxx,cdcyy,cdczz,cdcEdx,cdcEdy,cdcEdz,cdcHdx,cdcHdy,cdcHdz,cdcf
 real cdcE,cdcELap,cdcELapsq,cdcELapm,cdcHzxLap,cdcHzyLap
 real cdcH,cdcHLap,cdcHLapsq,cdcHLapm

 ! dispersion
 integer dispersionModel,numberOfPolarizationVectors,pxc,pyc,pzc,iv
 integer ec,pc,pce
 real gamma,omegap
 real gammaDt,omegapDtSq,ptt, fe,fp,fp2
 ! Generalized dispersion model parameters
 real alphaP, a0,a1,b0,b1
 real ev,evm,evn,pv0,pvm0,pvx,pvxE,deti,rhsE,rhsP

 integer gdmParOption
 integer maxNumberOfParameters,maxNumberOfPolarizationVectors
 parameter( maxNumberOfParameters=4, maxNumberOfPolarizationVectors=20 )
 real gdmPar(0:maxNumberOfParameters-1,0:maxNumberOfPolarizationVectors-1)
 real a0v,a1v,b0v,b1v
 real beta, pSum
 real pv(0:maxNumberOfPolarizationVectors-1)
 real pvm(0:maxNumberOfPolarizationVectors-1)
 real rhspv(0:maxNumberOfPolarizationVectors-1)
 real betav(0:maxNumberOfPolarizationVectors-1)
 real fpv(0:maxNumberOfPolarizationVectors-1)
 real fptv(0:maxNumberOfPolarizationVectors-1)
 real fpttv(0:maxNumberOfPolarizationVectors-1)
 real lapfpv(0:maxNumberOfPolarizationVectors-1)

! More Generalized dispersion model parameters
 real pSum0tt, pSum0ttt, pSum0tttt, pSum0xxtt, pSum0yytt, pSum0zztt, rhsPxx, pxxSum, PtttStar, QxxStar, EtxxStar, EtttStar, rhsP4, LHSev, exxv, exxvm, exxvn, rhsExx
 real b1tttv(0:maxNumberOfPolarizationVectors-1)
 real b0tttv(0:maxNumberOfPolarizationVectors-1)
 real a1tttv(0:maxNumberOfPolarizationVectors-1)
 real a0tttv(0:maxNumberOfPolarizationVectors-1)
 real a2tttv(0:maxNumberOfPolarizationVectors-1)
 real f2v(0:maxNumberOfPolarizationVectors-1)
 real f3v(0:maxNumberOfPolarizationVectors-1)
 real f6v(0:maxNumberOfPolarizationVectors-1)
 real pxxvn(0:maxNumberOfPolarizationVectors-1)
 real pxxv(0:maxNumberOfPolarizationVectors-1)
 real pxxvm(0:maxNumberOfPolarizationVectors-1)
 real rhspxxv(0:maxNumberOfPolarizationVectors-1)
 real LHSpv(0:maxNumberOfPolarizationVectors-1)
 real fp00v(0:maxNumberOfPolarizationVectors-1)
 real ptttStarv(0:maxNumberOfPolarizationVectors-1)
 real pvn(0:maxNumberOfPolarizationVectors-1)

 ! ----- multilevel atomic model -----
 integer numberOfAtomicLevels,maxPar,m1,m2,na,nce
 real q0,q0t,q0tt,q0ttt,q0tttt
 real pnec,prc,peptc
 parameter( maxPar=20 )
 real nlPar(0:maxPar-1,0:maxPar-1,0:2)
 real fnv(0:maxPar-1)
 real fntv(0:maxPar-1)
 real fnttv(0:maxPar-1)
 real fntttv(0:maxPar-1)
 real qvec(0:maxPar-1)
 real qt(0:maxPar-1)
 real qtt(0:maxPar-1)
 real qttt(0:maxPar-1)
 real qtttt(0:maxPar-1)
 real qelap2(0:maxPar-1)
 real qx2v(0:maxPar-1),qy2v(0:maxPar-1),qz2v(0:maxPar-1)
 real qlap2(0:maxPar-1)
 real et,ett,ettt
 real etv(0:2),ettv(0:2),etttv(0:2)
 real ptv(0:2,0:maxPar-1),pttv(0:2,0:maxPar-1),ptttv(0:2,0:maxPar-1),pttttv(0:2,0:maxPar-1)
 real ptttSum,lapfe,fet,fett
 

! .......statement functions for GDM parameters
 a0v(iv) = gdmPar(0,iv)
 a1v(iv) = gdmPar(1,iv)
 b0v(iv) = gdmPar(2,iv)
 b1v(iv) = gdmPar(3,iv)

 ! ..... statement functions for multilevel atomic model
 ! pnec  = polarizationNECoefficients
 ! prc   = populationRelaxationCoefficients
 ! peptc = populationEPtCoefficients
 pnec(m1,m2)  = nlPar(m1,m2,0)
 prc(m1,m2)  = nlPar(m1,m2,1)
 peptc(m1,m2) = nlPar(m1,m2,2)
 

!.......statement functions for jacobian
 rx(i1,i2,i3)=rsxy(i1,i2,i3,0,0)
 ry(i1,i2,i3)=rsxy(i1,i2,i3,0,1)
 rz(i1,i2,i3)=rsxy(i1,i2,i3,0,2)
 sx(i1,i2,i3)=rsxy(i1,i2,i3,1,0)
 sy(i1,i2,i3)=rsxy(i1,i2,i3,1,1)
 sz(i1,i2,i3)=rsxy(i1,i2,i3,1,2)
 tx(i1,i2,i3)=rsxy(i1,i2,i3,2,0)
 ty(i1,i2,i3)=rsxy(i1,i2,i3,2,1)
 tz(i1,i2,i3)=rsxy(i1,i2,i3,2,2)

!     The next macro will define the difference approximation statement functions
 defineDifferenceOrder2Components1(u,RX)
 defineDifferenceOrder4Components1(u,RX)

 defineDifferenceOrder2Components1(un,none)
 defineDifferenceOrder4Components1(un,none)

 defineDifferenceOrder2Components1(v,none)
 defineDifferenceOrder4Components1(v,none)

 defineDifferenceOrder2Components1(um,none)
 defineDifferenceOrder4Components1(um,none)

 defineDifferenceOrder2Components1(ff,none)

 defineDifferenceOrder2Components1(p,none)
 defineDifferenceOrder4Components1(p,none)

 defineDifferenceOrder2Components1(pm,none)
 defineDifferenceOrder4Components1(pm,none)

 ! MLA
 defineDifferenceOrder2Components1(q,none)

 ! 2nd-order in space and time
 maxwell2dr(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+\
            cdtdx*(u(i1-1,i2,i3,n)+u(i1+1,i2,i3,n)-2.*u(i1,i2,i3,n))+\
            cdtdy*(u(i1,i2-1,i3,n)+u(i1,i2+1,i3,n)-2.*u(i1,i2,i3,n))


 du(i1,i2,i3,c)=u(i1,i2,i3,c)-um(i1,i2,i3,c)

 ! D-zero in time (really undivided)
 DztU(i1,i2,i3,n) = (un(i1,i2,i3,n)-um(i1,i2,i3,n))

 ! D-plus in time (really undivided) (add factor of 2 below since formula assumes D0)
 DptU(i1,i2,i3,n) = (un(i1,i2,i3,n)-u(i1,i2,i3,n))
 ! D-minus in time (add factor of 2 below since formula assumes D0)
 DmtU(i1,i2,i3,n) = (u(i1,i2,i3,n)-um(i1,i2,i3,n))*2.

 ! special D-zero in time : assume u=u(t), um=u(t-dt),  un=u(t-2*dt)
 DzstU(i1,i2,i3,n) = (u(i1,i2,i3,n)-un(i1,i2,i3,n))


 maxwell3dr(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+\
            cdtdx*(u(i1-1,i2,i3,n)+u(i1+1,i2,i3,n)-2.*u(i1,i2,i3,n))+\
            cdtdy*(u(i1,i2-1,i3,n)+u(i1,i2+1,i3,n)-2.*u(i1,i2,i3,n))+\
            cdtdz*(u(i1,i2,i3-1,n)+u(i1,i2,i3+1,n)-2.*u(i1,i2,i3,n))

 ! these use pre-computed RHS in f
 maxwellc22(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+dtsq*f(i1,i2,i3,n)
 maxwellc23(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+dtsq*f(i1,i2,i3,n)

 ! 2D, 2nd-order, div cleaning:
 !    D+tD-t( E ) + alpha*( D0t E ) = c^2 Delta(E) + alpha*( (1/eps) Curl ( H ) )

 ! - rectangular:
 mxdc2d2Ex(i1,i2,i3) = cdc0*u(i1,i2,i3,ex)+cdc1*um(i1,i2,i3,ex)\
            +cdcxx*(u(i1-1,i2,i3,ex)+u(i1+1,i2,i3,ex)-2.*u(i1,i2,i3,ex))\
            +cdcyy*(u(i1,i2-1,i3,ex)+u(i1,i2+1,i3,ex)-2.*u(i1,i2,i3,ex))\
            +cdcEdy*( u(i1,i2+1,i3,hz)-u(i1,i2-1,i3,hz) )

 mxdc2d2Ey(i1,i2,i3) = cdc0*u(i1,i2,i3,ey)+cdc1*um(i1,i2,i3,ey)\
            +cdcxx*(u(i1-1,i2,i3,ey)+u(i1+1,i2,i3,ey)-2.*u(i1,i2,i3,ey))\
            +cdcyy*(u(i1,i2-1,i3,ey)+u(i1,i2+1,i3,ey)-2.*u(i1,i2,i3,ey))\
            -cdcEdx*( u(i1+1,i2,i3,hz)-u(i1-1,i2,i3,hz) )

 ! - 2D curvilinear:  (assumes f contains Delta u )
 mxdc2d2cEx(i1,i2,i3) = cdc0*u(i1,i2,i3,ex)+cdc1*um(i1,i2,i3,ex)\
                       + cdcf*f(i1,i2,i3,ex)\
                       + cdcE*uy22(i1,i2,i3,hz)

 mxdc2d2cEy(i1,i2,i3) = cdc0*u(i1,i2,i3,ey)+cdc1*um(i1,i2,i3,ey)\
                       + cdcf*f(i1,i2,i3,ey)\
                       - cdcE*ux22(i1,i2,i3,hz)

#If #DIM eq "3"
 ! 3D, 2nd-order, div cleaning:
 !    D+tD-t( E ) + alpha*( D0t E ) = c^2 Delta(E) + alpha*( (1/eps) Curl ( H ) )

 ! - 3D rectangular:
 mxdc3d2Ex(i1,i2,i3) = cdc0*u(i1,i2,i3,ex)+cdc1*um(i1,i2,i3,ex)\
            +cdcxx*(u(i1-1,i2,i3,ex)+u(i1+1,i2,i3,ex)-2.*u(i1,i2,i3,ex))\
            +cdcyy*(u(i1,i2-1,i3,ex)+u(i1,i2+1,i3,ex)-2.*u(i1,i2,i3,ex))\
            +cdczz*(u(i1,i2,i3-1,ex)+u(i1,i2,i3+1,ex)-2.*u(i1,i2,i3,ex))\
            +cdcEdy*( u(i1,i2+1,i3,hz)-u(i1,i2-1,i3,hz) )\
            -cdcEdz*( u(i1,i2,i3+1,hy)-u(i1,i2,i3-1,hy) )

 mxdc3d2Ey(i1,i2,i3) = cdc0*u(i1,i2,i3,ey)+cdc1*um(i1,i2,i3,ey)\
            +cdcxx*(u(i1-1,i2,i3,ey)+u(i1+1,i2,i3,ey)-2.*u(i1,i2,i3,ey))\
            +cdcyy*(u(i1,i2-1,i3,ey)+u(i1,i2+1,i3,ey)-2.*u(i1,i2,i3,ey))\
            +cdczz*(u(i1,i2,i3-1,ey)+u(i1,i2,i3+1,ey)-2.*u(i1,i2,i3,ey))\
            +cdcEdz*( u(i1,i2,i3+1,hx)-u(i1,i2,i3-1,hx) )\
            -cdcEdx*( u(i1+1,i2,i3,hz)-u(i1-1,i2,i3,hz) )

 mxdc3d2Ez(i1,i2,i3) = cdc0*u(i1,i2,i3,ez)+cdc1*um(i1,i2,i3,ez)\
            +cdcxx*(u(i1-1,i2,i3,ez)+u(i1+1,i2,i3,ez)-2.*u(i1,i2,i3,ez))\
            +cdcyy*(u(i1,i2-1,i3,ez)+u(i1,i2+1,i3,ez)-2.*u(i1,i2,i3,ez))\
            +cdczz*(u(i1,i2,i3-1,ez)+u(i1,i2,i3+1,ez)-2.*u(i1,i2,i3,ez))\
            +cdcEdx*( u(i1+1,i2,i3,hy)-u(i1-1,i2,i3,hy) )\
            -cdcEdy*( u(i1,i2+1,i3,hx)-u(i1,i2-1,i3,hx) )

 mxdc3d2Hx(i1,i2,i3) = cdc0*u(i1,i2,i3,hx)+cdc1*um(i1,i2,i3,hx)\
            +cdcxx*(u(i1-1,i2,i3,hx)+u(i1+1,i2,i3,hx)-2.*u(i1,i2,i3,hx))\
            +cdcyy*(u(i1,i2-1,i3,hx)+u(i1,i2+1,i3,hx)-2.*u(i1,i2,i3,hx))\
            +cdczz*(u(i1,i2,i3-1,hx)+u(i1,i2,i3+1,hx)-2.*u(i1,i2,i3,hx))\
            -cdcHdy*( u(i1,i2+1,i3,ez)-u(i1,i2-1,i3,ez) )\
            +cdcHdz*( u(i1,i2,i3+1,ey)-u(i1,i2,i3-1,ey) )

 mxdc3d2Hy(i1,i2,i3) = cdc0*u(i1,i2,i3,hy)+cdc1*um(i1,i2,i3,hy)\
            +cdcxx*(u(i1-1,i2,i3,hy)+u(i1+1,i2,i3,hy)-2.*u(i1,i2,i3,hy))\
            +cdcyy*(u(i1,i2-1,i3,hy)+u(i1,i2+1,i3,hy)-2.*u(i1,i2,i3,hy))\
            +cdczz*(u(i1,i2,i3-1,hy)+u(i1,i2,i3+1,hy)-2.*u(i1,i2,i3,hy))\
            -cdcHdz*( u(i1,i2,i3+1,ex)-u(i1,i2,i3-1,ex) )\
            +cdcHdx*( u(i1+1,i2,i3,ez)-u(i1-1,i2,i3,ez) )

 mxdc3d2Hz(i1,i2,i3) = cdc0*u(i1,i2,i3,hz)+cdc1*um(i1,i2,i3,hz)\
            +cdcxx*(u(i1-1,i2,i3,hz)+u(i1+1,i2,i3,hz)-2.*u(i1,i2,i3,hz))\
            +cdcyy*(u(i1,i2-1,i3,hz)+u(i1,i2+1,i3,hz)-2.*u(i1,i2,i3,hz))\
            +cdczz*(u(i1,i2,i3-1,hz)+u(i1,i2,i3+1,hz)-2.*u(i1,i2,i3,hz))\
            -cdcHdx*( u(i1+1,i2,i3,ey)-u(i1-1,i2,i3,ey) )\
            +cdcHdy*( u(i1,i2+1,i3,ex)-u(i1,i2-1,i3,ex) )

 ! - 3D curvilinear:  (assumes f contains Delta u )
 mxdc3d2cEx(i1,i2,i3) = cdc0*u(i1,i2,i3,ex)+cdc1*um(i1,i2,i3,ex)\
                       + cdcf*f(i1,i2,i3,ex)\
                       + cdcE*( uy23(i1,i2,i3,hz)\
                               -uz23(i1,i2,i3,hy))
 mxdc3d2cEy(i1,i2,i3) = cdc0*u(i1,i2,i3,ey)+cdc1*um(i1,i2,i3,ey)\
                       + cdcf*f(i1,i2,i3,ey)\
                       + cdcE*( uz23(i1,i2,i3,hx)\
                               -ux23(i1,i2,i3,hz))
 mxdc3d2cEz(i1,i2,i3) = cdc0*u(i1,i2,i3,ez)+cdc1*um(i1,i2,i3,ez)\
                       + cdcf*f(i1,i2,i3,ez)\
                       + cdcE*( ux23(i1,i2,i3,hy)\
                               -uy23(i1,i2,i3,hx))

 mxdc3d2cHx(i1,i2,i3) = cdc0*u(i1,i2,i3,hx)+cdc1*um(i1,i2,i3,hx)\
                       + cdcf*f(i1,i2,i3,hx)\
                       + cdcH*(-uy23(i1,i2,i3,ez)\
                               +uz23(i1,i2,i3,ey))
 mxdc3d2cHy(i1,i2,i3) = cdc0*u(i1,i2,i3,hy)+cdc1*um(i1,i2,i3,hy)\
                       + cdcf*f(i1,i2,i3,hy)\
                       + cdcH*(-uz23(i1,i2,i3,ex)\
                               +ux23(i1,i2,i3,ez))
 mxdc3d2cHz(i1,i2,i3) = cdc0*u(i1,i2,i3,hz)+cdc1*um(i1,i2,i3,hz)\
                       + cdcf*f(i1,i2,i3,hz)\
                       + cdcH*(-ux23(i1,i2,i3,ey)\
                               +uy23(i1,i2,i3,ex))


#End



!- ! Stoermer: 4th order in space and 4th order in time:
!- maxwellr44(i1,i2,i3,n)=2.*u(I1,I2,I3,n)-um(I1,I2,I3,n)+\
!-    c40*lap(n)+c41*v(I1,I2,I3,n)+c42*vvt2(I1,I2,I3,n)+c43*ut3(I1,I2,I3,n)
!-
!- maxwellc44(i1,i2,i3,n)=2.*u(I1,I2,I3,n)-um(I1,I2,I3,n)+\
!-    c40*f(i1,i2,i3,n)+c41*v(I1,I2,I3,n)+c42*vvt2(I1,I2,I3,n)+c43*ut3(I1,I2,I3,n)
!-
!- ! Stoermer: 6th order in space and 6th order in time:
!- maxwellr66(i1,i2,i3,n)=2.*u(I1,I2,I3,n)-um(I1,I2,I3,n)+\
!-    c60*lap(n)+c61*v(I1,I2,I3,n)+c62*vvt2(I1,I2,I3,n)+c63*ut3(I1,I2,I3,n)+\
!-    c64*vvt4(I1,I2,I3,n)+c65*ut5(I1,I2,I3,n)
!-
!- maxwellc66(i1,i2,i3,n)=2.*u(I1,I2,I3,n)-um(I1,I2,I3,n)+\
!-    c60*f(i1,i2,i3,n)+c61*v(I1,I2,I3,n)+c62*vvt2(I1,I2,I3,n)+c63*ut3(I1,I2,I3,n)+\
!-    c64*vvt4(I1,I2,I3,n)+c65*ut5(I1,I2,I3,n)
!-
!- ! Stoermer: 8th order in space and 8th order in time:
!- maxwellr88(i1,i2,i3,n)=2.*u(I1,I2,I3,n)-um(I1,I2,I3,n)+\
!-    c80*lap(n)+c81*v(I1,I2,I3,n)+c82*vvt2(I1,I2,I3,n)+c83*ut3(I1,I2,I3,n)+\
!-    c84*vvt4(I1,I2,I3,n)+c85*ut5(I1,I2,I3,n)+c86*ut6(I1,I2,I3,n)+c87*ut7(I1,I2,I3,n)
!-
!- maxwellc88(i1,i2,i3,n)=2.*u(I1,I2,I3,n)-um(I1,I2,I3,n)+\
!-    c80*f(i1,i2,i3,n)+c81*v(I1,I2,I3,n)+c82*vvt2(I1,I2,I3,n)+c83*ut3(I1,I2,I3,n)+\
!-    c84*vvt4(I1,I2,I3,n)+c85*ut5(I1,I2,I3,n)+c86*ut6(I1,I2,I3,n)+c87*ut7(I1,I2,I3,n)


 !    *** 2nd order ***
 lap2d2(i1,i2,i3,c)=(u(i1+1,i2,i3,c)-2.*u(i1,i2,i3,c)+u(i1-1,i2,i3,c))*dxsqi\
                   +(u(i1,i2+1,i3,c)-2.*u(i1,i2,i3,c)+u(i1,i2-1,i3,c))*dysqi

 lap2d2m(i1,i2,i3,c)=(um(i1+1,i2,i3,c)-2.*um(i1,i2,i3,c)+um(i1-1,i2,i3,c))*dxsqi\
                    +(um(i1,i2+1,i3,c)-2.*um(i1,i2,i3,c)+um(i1,i2-1,i3,c))*dysqi

 plap2d2(i1,i2,i3,c) =(p(i1+1,i2,i3,c)-2.*p(i1,i2,i3,c)+p(i1-1,i2,i3,c))*dxsqi\
                     +(p(i1,i2+1,i3,c)-2.*p(i1,i2,i3,c)+p(i1,i2-1,i3,c))*dysqi

 plap2d2m(i1,i2,i3,c)=(pm(i1+1,i2,i3,c)-2.*pm(i1,i2,i3,c)+pm(i1-1,i2,i3,c))*dxsqi\
                     +(pm(i1,i2+1,i3,c)-2.*pm(i1,i2,i3,c)+pm(i1,i2-1,i3,c))*dysqi

 lap3d2(i1,i2,i3,c)=(u(i1+1,i2,i3,c)-2.*u(i1,i2,i3,c)+u(i1-1,i2,i3,c))*dxsqi\
                   +(u(i1,i2+1,i3,c)-2.*u(i1,i2,i3,c)+u(i1,i2-1,i3,c))*dysqi\
                   +(u(i1,i2,i3+1,c)-2.*u(i1,i2,i3,c)+u(i1,i2,i3-1,c))*dzsqi

 lap3d2m(i1,i2,i3,c)=(um(i1+1,i2,i3,c)-2.*um(i1,i2,i3,c)+um(i1-1,i2,i3,c))*dxsqi\
                    +(um(i1,i2+1,i3,c)-2.*um(i1,i2,i3,c)+um(i1,i2-1,i3,c))*dysqi\
                    +(um(i1,i2,i3+1,c)-2.*um(i1,i2,i3,c)+um(i1,i2,i3-1,c))*dzsqi


 plap3d2(i1,i2,i3,c)=(p(i1+1,i2,i3,c)-2.*p(i1,i2,i3,c)+p(i1-1,i2,i3,c))*dxsqi\
                    +(p(i1,i2+1,i3,c)-2.*p(i1,i2,i3,c)+p(i1,i2-1,i3,c))*dysqi\
                    +(p(i1,i2,i3+1,c)-2.*p(i1,i2,i3,c)+p(i1,i2,i3-1,c))*dzsqi

 plap3d2m(i1,i2,i3,c)=(pm(i1+1,i2,i3,c)-2.*pm(i1,i2,i3,c)+pm(i1-1,i2,i3,c))*dxsqi\
                    +(pm(i1,i2+1,i3,c)-2.*pm(i1,i2,i3,c)+pm(i1,i2-1,i3,c))*dysqi\
                    +(pm(i1,i2,i3+1,c)-2.*pm(i1,i2,i3,c)+pm(i1,i2,i3-1,c))*dzsqi

 lap2d2f(i1,i2,i3,c,m)=(fa(i1+1,i2,i3,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1-1,i2,i3,c,m))*dxsqi\
                      +(fa(i1,i2+1,i3,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1,i2-1,i3,c,m))*dysqi

 lap3d2f(i1,i2,i3,c,m)=(fa(i1+1,i2,i3,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1-1,i2,i3,c,m))*dxsqi\
                      +(fa(i1,i2+1,i3,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1,i2-1,i3,c,m))*dysqi\
                      +(fa(i1,i2,i3+1,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1,i2,i3-1,c,m))*dzsqi

 ! MLA
 ! first derivatives of E
 ex2(i1,i2,i3,c)=(u(i1+1,i2,i3,c)-u(i1-1,i2,i3,c))/(2.*dx(0))
 ey2(i1,i2,i3,c)=(u(i1,i2+1,i3,c)-u(i1,i2-1,i3,c))/(2.*dx(1))
 ez2(i1,i2,i3,c)=(u(i1,i2,i3+1,c)-u(i1,i2,i3-1,c))/(2.*dx(2))

 ! first derivatives of N
 qx2(i1,i2,i3,c)=(q(i1+1,i2,i3,c)-q(i1-1,i2,i3,c))/(2.*dx(0))
 qy2(i1,i2,i3,c)=(q(i1,i2+1,i3,c)-q(i1,i2-1,i3,c))/(2.*dx(1))
 qz2(i1,i2,i3,c)=(q(i1,i2,i3+1,c)-q(i1,i2,i3-1,c))/(2.*dx(2))

 ! laplacians of N
 qlap2d2(i1,i2,i3,c) = (q(i1+1,i2,i3,c)-2.*q(i1,i2,i3,c)+q(i1-1,i2,i3,c))*dxsqi\
                       +(q(i1,i2+1,i3,c)-2.*q(i1,i2,i3,c)+q(i1,i2-1,i3,c))*dysqi
 qlap3d2(i1,i2,i3,c) = (q(i1+1,i2,i3,c)-2.*q(i1,i2,i3,c)+q(i1-1,i2,i3,c))*dxsqi\
                       +(q(i1,i2+1,i3,c)-2.*q(i1,i2,i3,c)+q(i1,i2-1,i3,c))*dysqi\
                       +(q(i1,i2,i3+1,c)-2.*q(i1,i2,i3,c)+q(i1,i2,i3-1,c))*dzsqi


 ! 2D laplacian squared = u.xxxx + 2 u.xxyy + u.yyyy
 lap2d2Pow2(i1,i2,i3,c)= ( 6.*u(i1,i2,i3,c)   \
   - 4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))    \
       +(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) )*dxi4 \
   +( 6.*u(i1,i2,i3,c)    \
    -4.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))    \
       +(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) )*dyi4  \
   +( 8.*u(i1,i2,i3,c)     \
    -4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)+u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))   \
    +2.*(u(i1+1,i2+1,i3,c)+u(i1-1,i2+1,i3,c)+u(i1+1,i2-1,i3,c)+u(i1-1,i2-1,i3,c)) )*dxdyi2

 ! 3D laplacian squared = u.xxxx + u.yyyy + u.zzzz + 2 (u.xxyy + u.xxzz + u.yyzz )
 lap3d2Pow2(i1,i2,i3,c)= ( 6.*u(i1,i2,i3,c)   \
   - 4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))    \
       +(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) )*dxi4 \
  +(  +6.*u(i1,i2,i3,c)    \
    -4.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))    \
       +(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) )*dyi4\
  +(  +6.*u(i1,i2,i3,c)    \
    -4.*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))    \
       +(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) )*dzi4\
   +(8.*u(i1,i2,i3,c)     \
    -4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)+u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))   \
    +2.*(u(i1+1,i2+1,i3,c)+u(i1-1,i2+1,i3,c)+u(i1+1,i2-1,i3,c)+u(i1-1,i2-1,i3,c)) )*dxdyi2 \
   +(8.*u(i1,i2,i3,c)     \
    -4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)+u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))   \
    +2.*(u(i1+1,i2,i3+1,c)+u(i1-1,i2,i3+1,c)+u(i1+1,i2,i3-1,c)+u(i1-1,i2,i3-1,c)) )*dxdzi2 \
   +(8.*u(i1,i2,i3,c)     \
    -4.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)+u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))   \
    +2.*(u(i1,i2+1,i3+1,c)+u(i1,i2-1,i3+1,c)+u(i1,i2+1,i3-1,c)+u(i1,i2-1,i3-1,c)) )*dydzi2

 lap2d2Pow3(i1,i2,i3,c)=LAP2D2(lap2d2Pow2,i1,i2,i3,c)

 lap3d2Pow3(i1,i2,i3,c)=LAP3D2(lap3d2Pow2,i1,i2,i3,c)

 lap2d2Pow4(i1,i2,i3,c)=LAP2D2POW2(lap2d2Pow2,i1,i2,i3,c)
 lap3d2Pow4(i1,i2,i3,c)=LAP3D2POW2(lap3d2Pow2,i1,i2,i3,c)

!    ** 4th order ****

 lap2d4(i1,i2,i3,c)=( -30.*u(i1,i2,i3,c)     \
   +16.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))     \
       -(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) )*dxsq12i + \
  ( -30.*u(i1,i2,i3,c)     \
   +16.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))     \
       -(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) )*dysq12i

 lap2d4m(i1,i2,i3,c)=( -30.*um(i1,i2,i3,c)     \
   +16.*(um(i1+1,i2,i3,c)+um(i1-1,i2,i3,c))     \
       -(um(i1+2,i2,i3,c)+um(i1-2,i2,i3,c)) )*dxsq12i + \
  ( -30.*um(i1,i2,i3,c)     \
   +16.*(um(i1,i2+1,i3,c)+um(i1,i2-1,i3,c))     \
       -(um(i1,i2+2,i3,c)+um(i1,i2-2,i3,c)) )*dysq12i

 plap2d4(i1,i2,i3,c)=( -30.*p(i1,i2,i3,c)     \
   +16.*(p(i1+1,i2,i3,c)+p(i1-1,i2,i3,c))     \
       -(p(i1+2,i2,i3,c)+p(i1-2,i2,i3,c)) )*dxsq12i + \
  ( -30.*p(i1,i2,i3,c)     \
   +16.*(p(i1,i2+1,i3,c)+p(i1,i2-1,i3,c))     \
       -(p(i1,i2+2,i3,c)+p(i1,i2-2,i3,c)) )*dysq12i

 plap2d4m(i1,i2,i3,c)=( -30.*pm(i1,i2,i3,c)     \
   +16.*(pm(i1+1,i2,i3,c)+pm(i1-1,i2,i3,c))     \
       -(pm(i1+2,i2,i3,c)+pm(i1-2,i2,i3,c)) )*dxsq12i + \
  ( -30.*pm(i1,i2,i3,c)     \
   +16.*(pm(i1,i2+1,i3,c)+pm(i1,i2-1,i3,c))     \
       -(pm(i1,i2+2,i3,c)+pm(i1,i2-2,i3,c)) )*dysq12i


 lap3d4(i1,i2,i3,c)=lap2d4(i1,i2,i3,c)+ \
  ( -30.*u(i1,i2,i3,c)      \
   +16.*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))      \
       -(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) )*dzsq12i

 lap3d4m(i1,i2,i3,c)=lap2d4m(i1,i2,i3,c)+ \
  ( -30.*um(i1,i2,i3,c)      \
   +16.*(um(i1,i2,i3+1,c)+um(i1,i2,i3-1,c))      \
       -(um(i1,i2,i3+2,c)+um(i1,i2,i3-2,c)) )*dzsq12i

 plap3d4(i1,i2,i3,c)=plap2d4(i1,i2,i3,c)+ \
  ( -30.*p(i1,i2,i3,c)      \
   +16.*(p(i1,i2,i3+1,c)+p(i1,i2,i3-1,c))      \
       -(p(i1,i2,i3+2,c)+p(i1,i2,i3-2,c)) )*dzsq12i

  plap3d4m(i1,i2,i3,c)=plap2d4m(i1,i2,i3,c)+ \
  ( -30.*pm(i1,i2,i3,c)      \
   +16.*(pm(i1,i2,i3+1,c)+pm(i1,i2,i3-1,c))      \
       -(pm(i1,i2,i3+2,c)+pm(i1,i2,i3-2,c)) )*dzsq12i


 lap2d4Pow2(i1,i2,i3,c)=LAP2D4(lap2d4,i1,i2,i3,c)
 lap3d4Pow2(i1,i2,i3,c)=LAP3D4(lap3d4,i1,i2,i3,c)

 lap2d4Pow3(i1,i2,i3,c)=LAP2D4(lap2d4Pow2,i1,i2,i3,c)
 lap3d4Pow3(i1,i2,i3,c)=LAP3D4(lap3d4Pow2,i1,i2,i3,c)

!     *** 6th order ***

 lap2d6(i1,i2,i3,c)= \
          c00lap2d6*u(i1,i2,i3,c)     \
         +c10lap2d6*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)) \
         +c01lap2d6*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)) \
         +c20lap2d6*(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) \
         +c02lap2d6*(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) \
         +c30lap2d6*(u(i1+3,i2,i3,c)+u(i1-3,i2,i3,c)) \
         +c03lap2d6*(u(i1,i2+3,i3,c)+u(i1,i2-3,i3,c))

 lap3d6(i1,i2,i3,c)=\
          c000lap3d6*u(i1,i2,i3,c) \
         +c100lap3d6*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)) \
         +c010lap3d6*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)) \
         +c001lap3d6*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c)) \
         +c200lap3d6*(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) \
         +c020lap3d6*(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) \
         +c002lap3d6*(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) \
         +c300lap3d6*(u(i1+3,i2,i3,c)+u(i1-3,i2,i3,c)) \
         +c030lap3d6*(u(i1,i2+3,i3,c)+u(i1,i2-3,i3,c)) \
         +c003lap3d6*(u(i1,i2,i3+3,c)+u(i1,i2,i3-3,c))

 lap2d6Pow2(i1,i2,i3,c)=LAP2D6(lap2d6,i1,i2,i3,c)
 lap3d6Pow2(i1,i2,i3,c)=LAP3D6(lap3d6,i1,i2,i3,c)


!     *** 8th order ***

 lap2d8(i1,i2,i3,c)=c00lap2d8*u(i1,i2,i3,c)      \
          +c10lap2d8*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))     \
          +c01lap2d8*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)) \
          +c20lap2d8*(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c))  \
          +c02lap2d8*(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) \
          +c30lap2d8*(u(i1+3,i2,i3,c)+u(i1-3,i2,i3,c))  \
          +c03lap2d8*(u(i1,i2+3,i3,c)+u(i1,i2-3,i3,c)) \
          +c40lap2d8*(u(i1+4,i2,i3,c)+u(i1-4,i2,i3,c))  \
          +c04lap2d8*(u(i1,i2+4,i3,c)+u(i1,i2-4,i3,c))

 lap3d8(i1,i2,i3,c)=c000lap3d8*u(i1,i2,i3,c)      \
          +c100lap3d8*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))     \
          +c010lap3d8*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)) \
          +c001lap3d8*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c)) \
          +c200lap3d8*(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c))  \
          +c020lap3d8*(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) \
          +c002lap3d8*(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) \
          +c300lap3d8*(u(i1+3,i2,i3,c)+u(i1-3,i2,i3,c))  \
          +c030lap3d8*(u(i1,i2+3,i3,c)+u(i1,i2-3,i3,c)) \
          +c003lap3d8*(u(i1,i2,i3+3,c)+u(i1,i2,i3-3,c)) \
          +c400lap3d8*(u(i1+4,i2,i3,c)+u(i1-4,i2,i3,c))  \
          +c040lap3d8*(u(i1,i2+4,i3,c)+u(i1,i2-4,i3,c)) \
          +c004lap3d8*(u(i1,i2,i3+4,c)+u(i1,i2,i3-4,c))

! ******* artificial dissipation ******

!      (2nd difference)
 fd22d(i1,i2,i3,c)= \
 (     ( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c) ) \
  -4.*du(i1,i2,i3,c) )
!
 fd23d(i1,i2,i3,c)=\
 (     ( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c)+du(i1,i2,i3-1,c)+du(i1,i2,i3+1,c) ) \
   -6.*du(i1,i2,i3,c) )

!     -(fourth difference)
 fd42d(i1,i2,i3,c)= \
 (    -( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c) ) \
   +4.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c) ) \
  -12.*du(i1,i2,i3,c) )
!
 fd43d(i1,i2,i3,c)=\
 (    -( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c)+du(i1,i2,i3-2,c)+du(i1,i2,i3+2,c) ) \
   +4.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c)+du(i1,i2,i3-1,c)+du(i1,i2,i3+1,c) ) \
  -18.*du(i1,i2,i3,c) )

 ! (sixth  difference)
 fd62d(i1,i2,i3,c)= \
 (     ( du(i1-3,i2,i3,c)+du(i1+3,i2,i3,c)+du(i1,i2-3,i3,c)+du(i1,i2+3,i3,c) ) \
   -6.*( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c) ) \
  +15.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c) ) \
  -40.*du(i1,i2,i3,c) )

 fd63d(i1,i2,i3,c)=\
 (     ( du(i1-3,i2,i3,c)+du(i1+3,i2,i3,c)+du(i1,i2-3,i3,c)+du(i1,i2+3,i3,c)+du(i1,i2,i3-3,c)+du(i1,i2,i3+3,c) ) \
   -6.*( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c)+du(i1,i2,i3-2,c)+du(i1,i2,i3+2,c) ) \
  +15.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c)+du(i1,i2,i3-1,c)+du(i1,i2,i3+1,c) ) \
  -60.*du(i1,i2,i3,c) )

 ! -(eighth  difference)
 fd82d(i1,i2,i3,c)= \
 (    -( du(i1-4,i2,i3,c)+du(i1+4,i2,i3,c)+du(i1,i2-4,i3,c)+du(i1,i2+4,i3,c) ) \
   +8.*( du(i1-3,i2,i3,c)+du(i1+3,i2,i3,c)+du(i1,i2-3,i3,c)+du(i1,i2+3,i3,c) ) \
  -28.*( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c) ) \
  +56.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c) ) \
 -140.*du(i1,i2,i3,c) )

 fd83d(i1,i2,i3,c)=\
 (    -( du(i1-4,i2,i3,c)+du(i1+4,i2,i3,c)+du(i1,i2-4,i3,c)+du(i1,i2+4,i3,c)+du(i1,i2,i3-4,c)+du(i1,i2,i3+4,c) ) \
   +8.*( du(i1-3,i2,i3,c)+du(i1+3,i2,i3,c)+du(i1,i2-3,i3,c)+du(i1,i2+3,i3,c)+du(i1,i2,i3-3,c)+du(i1,i2,i3+3,c) ) \
  -28.*( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c)+du(i1,i2,i3-2,c)+du(i1,i2,i3+2,c) ) \
  +56.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c)+du(i1,i2,i3-1,c)+du(i1,i2,i3+1,c) ) \
 -210.*du(i1,i2,i3,c) )


!     **** Modified equation method: ****

 maxwell2dr44me(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+cdtsq*lap2d4(i1,i2,i3,n)\
                            +cdtsq12*lap2d2Pow2(i1,i2,i3,n)
 maxwell3dr44me(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+cdtsq*lap3d4(i1,i2,i3,n)\
                            +cdtsq12*lap3d2Pow2(i1,i2,i3,n)

  maxwell2dr66me(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+cdtsq*lap2d6(i1,i2,i3,n)\
                            +cdtsq12  *lap2d4Pow2(i1,i2,i3,n)\
                            +cdt4by360*lap2d2Pow3(i1,i2,i3,n)
 maxwell3dr66me(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+cdtsq*lap3d6(i1,i2,i3,n)\
                            +cdtsq12*  lap3d4Pow2(i1,i2,i3,n)\
                            +cdt4by360*lap3d2Pow3(i1,i2,i3,n)

 maxwell2dr88me(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+cdtsq*lap2d8(i1,i2,i3,n)\
                            +cdtsq12*lap2d6Pow2(i1,i2,i3,n)\
                            +cdt4by360*lap2d4Pow3(i1,i2,i3,n)+cdt6by20160*lap2d2Pow4(i1,i2,i3,n)
 maxwell3dr88me(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+cdtsq*lap3d8(i1,i2,i3,n)\
                            +cdtsq12*lap3d6Pow2(i1,i2,i3,n)\
                            +cdt4by360*lap3d4Pow3(i1,i2,i3,n)+cdt6by20160*lap3d2Pow4(i1,i2,i3,n)

! *********NEW forcing method (for user defined forcing)**********
!    -- forcing correction for modified equation method ---
!        RHS = f + (dt^2/12)*( c^2 * Delta f + f_tt )
!  Approximate the term in brackets to 2nd-order
! ---- Cartesian grids:
! f2drme44(i1,i2,i3,n) = fa(i1,i2,i3,n,fcur)
! f2drme44(i1,i2,i3,n) = f(i1,i2,i3,n)
f2drme44(i1,i2,i3,n) = fa(i1,i2,i3,n,fcur)+cdtSqBy12*lap2d2f(i1,i2,i3,n,fcur) \
                       +(fa(i1,i2,i3,n,fnext)-2.*fa(i1,i2,i3,n,fcur)+fa(i1,i2,i3,n,fprev))/(12.)

f3drme44(i1,i2,i3,n) = fa(i1,i2,i3,n,fcur)+cdtSqBy12*lap3d2f(i1,i2,i3,n,fcur) \
                       +(fa(i1,i2,i3,n,fnext)-2.*fa(i1,i2,i3,n,fcur)+fa(i1,i2,i3,n,fprev))/(12.)

! ---- Curvilinear grids
ff(i1,i2,i3,n) = fa(i1,i2,i3,n,fcur)
f2dcme44(i1,i2,i3,n) = fa(i1,i2,i3,n,fcur)+cdtSqBy12*ffLaplacian22(i1,i2,i3,n) \
                       +(fa(i1,i2,i3,n,fnext)-2.*fa(i1,i2,i3,n,fcur)+fa(i1,i2,i3,n,fprev))/(12.)

f3dcme44(i1,i2,i3,n) = fa(i1,i2,i3,n,fcur)+cdtSqBy12*ffLaplacian23(i1,i2,i3,n) \
                       +(fa(i1,i2,i3,n,fnext)-2.*fa(i1,i2,i3,n,fcur)+fa(i1,i2,i3,n,fprev))/(12.)


 ! f  = csq*Lap4(u)+f,  v= (csq*Lap2)**2
 maxwellc44me(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+dtsq*f(i1,i2,i3,n)+dt4by12*v(I1,I2,I3,n)
 ! these next are only valid for second order accuracy in time:
 maxwellc66me(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+dtsq*f(i1,i2,i3,n)
 maxwellc88me(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+dtsq*f(i1,i2,i3,n)

 ! for non-conservative modified-equation:
 max2dc44me(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+cdtsq*uLaplacian42(i1,i2,i3,n)+cdtsq12*uLapSq(n)

 ! This version for the 2-stage computation:

!$$$ vr2(i1,i2,i3,kd)=(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))*d12(0)
!$$$ vs2(i1,i2,i3,kd)=(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))*d12(1)
!$$$
!$$$ vrr2(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)) )*d22(0)
!$$$ vss2(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd)) )*d22(1)
!$$$ vrs2(i1,i2,i3,kd)=(vr2(i1,i2+1,i3,kd)-vr2(i1,i2-1,i3,kd))*d12(1)
!$$$
!$$$ vlaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*\
!$$$      vrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*\
!$$$      sy(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**\
!$$$      2)*vss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*vr2(i1,\
!$$$      i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*vs2(i1,i2,i3,kd)

 max2dc44me2(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+cdtsq*uLaplacian42(i1,i2,i3,n)\
                                                       +cdtsq12*vLaplacian22(i1,i2,i3,n)

 max3dc44me(i1,i2,i3,n)=2.*u(i1,i2,i3,n)-um(i1,i2,i3,n)+cdtsq*uLaplacian43(i1,i2,i3,n)\
                                                      +cdtsq12*vLaplacian23(i1,i2,i3,n)

 ! 2D, 4th-order, div cleaning:
 ! We could further optimize the D0x(lap2d2) ...
!!$ mxdc2d4Ex(i1,i2,i3) = cdc0*u(i1,i2,i3,ex)+cdc1*um(i1,i2,i3,ex)\
!!$            +cdcLap*lap2d4(i1,i2,i3,ex)\
!!$            +cdcLapsq*lap2d2Pow2(i1,i2,i3,ex)\
!!$            +cdcHz*uy42r(i1,i2,i3,hz)\
!!$            +cdcHzyLap*( lap2d2(i1,i2+1,i3,hz)-lap2d2(i1,i2-1,i3,hz) )
!!$
!!$ mxdc2d4Ey(i1,i2,i3) = cdc0*u(i1,i2,i3,ey)+cdc1*um(i1,i2,i3,ey)\
!!$            +cdcLap*lap2d4(i1,i2,i3,ey)\
!!$            +cdcLapsq*lap2d2Pow2(i1,i2,i3,ey)\
!!$            -cdcHz*ux42r(i1,i2,i3,hz)\
!!$            -cdcHzxLap*( lap2d2(i1+1,i2,i3,hz)-lap2d2(i1-1,i2,i3,hz) )
#If #DIM eq "2"
 ! 2D, 4th-order, rectangular, div cleaning:
 ! new version : here we replace curl( Delta H ) by an approx. to E_ttt that uses Delta E_t
 mxdc2d4Ex(i1,i2,i3) = cdc0*u(i1,i2,i3,ex)+cdc1*um(i1,i2,i3,ex)\
            +cdcELap*lap2d4(i1,i2,i3,ex)\
            +cdcELapsq*lap2d2Pow2(i1,i2,i3,ex)\
            +cdcE*uy42r(i1,i2,i3,hz)\
            +cdcELapm*( lap2d2m(i1,i2,i3,ex) )

 mxdc2d4Ey(i1,i2,i3) = cdc0*u(i1,i2,i3,ey)+cdc1*um(i1,i2,i3,ey)\
            +cdcELap*lap2d4(i1,i2,i3,ey)\
            +cdcELapsq*lap2d2Pow2(i1,i2,i3,ey)\
            -cdcE*ux42r(i1,i2,i3,hz)\
            +cdcELapm*( lap2d2m(i1,i2,i3,ey) )

 ! 2d, 4th order, curvilinear (conservative), div cleaning (f=Lap(E), v=Lapsq(E))
 mxdc2d4cConsEx(i1,i2,i3) = cdc0*u(i1,i2,i3,ex)+cdc1*um(i1,i2,i3,ex)\
            +cdcELap*f(i1,i2,i3,ex)\
            +cdcELapsq*v(i1,i2,i3,ex)\
            +cdcE*uy42(i1,i2,i3,hz)\
            +cdcELapm*( umLaplacian22(i1,i2,i3,ex) )

 mxdc2d4cConsEy(i1,i2,i3) = cdc0*u(i1,i2,i3,ey)+cdc1*um(i1,i2,i3,ey)\
            +cdcELap*f(i1,i2,i3,ey)\
            +cdcELapsq*v(i1,i2,i3,ey)\
            -cdcE*ux42(i1,i2,i3,hz)\
            +cdcELapm*( umLaplacian22(i1,i2,i3,ey) )


 ! 2D, 4th-order, curvilinear, div cleaning: **check me**
 mxdc2d4cEx(i1,i2,i3) = cdc0*u(i1,i2,i3,ex)+cdc1*um(i1,i2,i3,ex)\
            +cdcELap*uLaplacian42(i1,i2,i3,ex)\
            +cdcELapsq*vLaplacian22(i1,i2,i3,ex)\
            +cdcE*uy42(i1,i2,i3,hz)\
            +cdcELapm*( umLaplacian22(i1,i2,i3,ex) )

 mxdc2d4cEy(i1,i2,i3) = cdc0*u(i1,i2,i3,ey)+cdc1*um(i1,i2,i3,ey)\
            +cdcELap*uLaplacian42(i1,i2,i3,ey)\
            +cdcELapsq*vLaplacian22(i1,i2,i3,ey)\
            -cdcE*ux42(i1,i2,i3,hz)\
            +cdcELapm*( umLaplacian22(i1,i2,i3,ey) )
#End
#If #DIM eq "3"
 ! 3D, 4th-order, rectangular, div cleaning:
 mxdc3d4Ex(i1,i2,i3) = cdc0*u(i1,i2,i3,ex)+cdc1*um(i1,i2,i3,ex)\
            +cdcELap*lap3d4(i1,i2,i3,ex)\
            +cdcELapsq*lap3d2Pow2(i1,i2,i3,ex)\
            +cdcE*( uy43r(i1,i2,i3,hz) \
                   -uz43r(i1,i2,i3,hy) )\
            +cdcELapm*( lap3d2m(i1,i2,i3,ex) )

 mxdc3d4Ey(i1,i2,i3) = cdc0*u(i1,i2,i3,ey)+cdc1*um(i1,i2,i3,ey)\
            +cdcELap*lap3d4(i1,i2,i3,ey)\
            +cdcELapsq*lap3d2Pow2(i1,i2,i3,ey)\
            +cdcE*( uz43r(i1,i2,i3,hx) \
                   -ux43r(i1,i2,i3,hz) )\
            +cdcELapm*( lap3d2m(i1,i2,i3,ey) )

 mxdc3d4Ez(i1,i2,i3) = cdc0*u(i1,i2,i3,ez)+cdc1*um(i1,i2,i3,ez)\
            +cdcELap*lap3d4(i1,i2,i3,ez)\
            +cdcELapsq*lap3d2Pow2(i1,i2,i3,ez)\
            +cdcE*( ux43r(i1,i2,i3,hy) \
                   -uy43r(i1,i2,i3,hx) )\
            +cdcELapm*( lap3d2m(i1,i2,i3,ez) )

 mxdc3d4Hx(i1,i2,i3) = cdc0*u(i1,i2,i3,hx)+cdc1*um(i1,i2,i3,hx)\
            +cdcHLap*lap3d4(i1,i2,i3,hx)\
            +cdcHLapsq*lap3d2Pow2(i1,i2,i3,hx)\
            +cdcH*(-uy43r(i1,i2,i3,ez) \
                   +uz43r(i1,i2,i3,ey) )\
            +cdcHLapm*( lap3d2m(i1,i2,i3,hx) )

 mxdc3d4Hy(i1,i2,i3) = cdc0*u(i1,i2,i3,hy)+cdc1*um(i1,i2,i3,hy)\
            +cdcHLap*lap3d4(i1,i2,i3,hy)\
            +cdcHLapsq*lap3d2Pow2(i1,i2,i3,hy)\
            +cdcH*(-uz43r(i1,i2,i3,ex) \
                   +ux43r(i1,i2,i3,ez) )\
            +cdcHLapm*( lap3d2m(i1,i2,i3,hy) )

 mxdc3d4Hz(i1,i2,i3) = cdc0*u(i1,i2,i3,hz)+cdc1*um(i1,i2,i3,hz)\
            +cdcHLap*lap3d4(i1,i2,i3,hz)\
            +cdcHLapsq*lap3d2Pow2(i1,i2,i3,hz)\
            +cdcH*(-ux43r(i1,i2,i3,ey) \
                   +uy43r(i1,i2,i3,ex) )\
            +cdcHLapm*( lap3d2m(i1,i2,i3,ez) )

#End
!...........end   statement functions



 cc    =rpar(0)  ! this is c
 dt    =rpar(1)
 dx(0) =rpar(2)
 dx(1) =rpar(3)
 dx(2) =rpar(4)
 adc   =rpar(5)  ! coefficient of artificial dissipation
 add   =rpar(6)  ! coefficient of divergence damping
 dr(0) =rpar(7)
 dr(1) =rpar(8)
 dr(2) =rpar(9)
 eps   =rpar(10)
 mu    =rpar(11)
 kx    =rpar(12)
 ky    =rpar(13)
 kz    =rpar(14)
 sigmaE=rpar(15)  ! electric conductivity (for lossy materials, complex index of refraction)
 sigmaH=rpar(16)  ! magnetic conductivity
 divergenceCleaningCoefficient=rpar(17)
 t     =rpar(18)
 ep    =rpar(19)  ! for TZ forcing  -- new *wdh* Sept 2, 2017
 rpar(20)=0.  ! return the time used for adding dissipation

 ! dx(0) = dr(0)
 ! dx(1) = dr(1)
 ! dx(2) = dr(2)

 ! Drude-Lorentz dispersion model:
 gamma= rpar(21)
 omegap=rpar(22)

 sosupParameter=rpar(23)

 ! No need to pass these anymore:
 alphaP=rpar(24)
 a0    =rpar(25)
 a1    =rpar(26)
 b0    =rpar(27)
 b1    =rpar(28)


 dy=dx(1)  ! Are these needed?
 dz=dx(2)

 ! timeForArtificialDissipation=rpar(6) ! return value

 option             =ipar(0)
 gridType           =ipar(1)
 orderOfAccuracy    =ipar(2)
 orderInTime        =ipar(3)
 addForcing         =ipar(4)
 orderOfDissipation =ipar(5)
 ex                 =ipar(6)
 ey                 =ipar(7)
 ez                 =ipar(8)
 hx                 =ipar(9)
 hy                 =ipar(10)
 hz                 =ipar(11)
 solveForE          =ipar(12)
 solveForH          =ipar(13)
 useWhereMask       =ipar(14)
 timeSteppingMethod =ipar(15)
 useVariableDissipation        =ipar(16)
 useCurvilinearOpt             =ipar(17)
 useConservative               =ipar(18)
 combineDissipationWithAdvance =ipar(19)
 useDivergenceCleaning         =ipar(20)
 useNewForcingMethod           =ipar(21)
 numberOfForcingFunctions      =ipar(22)
 fcur                          =ipar(23)
 dispersionModel     =ipar(24)
 pxc                 =ipar(25)
 pyc                 =ipar(26)
 pzc                 =ipar(27)
 numberOfPolarizationVectors =ipar(28)
 grid                =ipar(29)
 nonlinearModel      =ipar(30)
 debug               =ipar(31)
 
 ! rxc                 =ipar(31) ! for future use
 ! ryc                 =ipar(32) ! for future use
 ! rzc                 =ipar(33) ! for future use

 useSosupDissipation   =ipar(34)
 sosupDissipationOption=ipar(35)
 updateSolution        =ipar(36)
 updateDissipation     =ipar(37)
 computeUt             =ipar(38)
 forcingOption         =ipar(39) ! new *wdh* Sept 2, 2017

 fprev = mod(fcur-1+numberOfForcingFunctions,max(1,numberOfForcingFunctions))
 fnext = mod(fcur+1                         ,max(1,numberOfForcingFunctions))

 if( t.le.3*dt .and. debug.gt.1 )then
   write(*,'(/,">>>> Inside advOptNew... t=",e10.3," grid=",i3)') t,grid
 end if

 ! addDissipation=.true. if we add the dissipation in the dis(i1,i2,i3,c) array
 !  if combineDissipationWithAdvance.ne.0 we compute the dissipation on the fly in the time step
 !  rather than pre-computing it in diss(i1,i2,i3,c)
 addDissipation = adc.gt.0. .and. combineDissipationWithAdvance.eq.0
 adcdt=adc*dt

 csq=cc**2
 dtsq=dt**2

 cdt=cc*dt

 cdtsq=(cc**2)*(dt**2)
 cdtsq12=cdtsq*cdtsq/12.  ! c^4 dt^4 /14
 cdt4by360=(cdt)**4/360.
 cdt6by20160=cdt**6/(8.*7.*6.*5.*4.*3.)

 cdtSqBy12= cdtsq/12.   ! c^2*dt*2/12

 dt4by12=dtsq*dtsq/12.

 cdtdx = (cc*dt/dx(0))**2
 cdtdy = (cc*dt/dy)**2
 cdtdz = (cc*dt/dz)**2

 dxsqi=1./(dx(0)**2)
 dysqi=1./(dy**2)
 dzsqi=1./(dz**2)

 dxsq12i=1./(12.*dx(0)**2)
 dysq12i=1./(12.*dy**2)
 dzsq12i=1./(12.*dz**2)

 dxi4=1./(dx(0)**4)
 dyi4=1./(dy**4)
 dxdyi2=1./(dx(0)*dx(0)*dy*dy)

 dzi4=1./(dz**4)
 dxdzi2=1./(dx(0)*dx(0)*dz*dz)
 dydzi2=1./(dy*dy*dz*dz)

 gammaDt=gamma*dt
 omegapDtSq=(omegap*dt)**2

 ! write(*,*) 'Before calling getGDMparameters...'

 gdmParOption=1 ! scale a0 and a1 by eps 
 if( dispersionModel.ne.noDispersion )then
   ! get the gdm parameters
   !   gdmPar(0:3,iv) = (a0,a1,b0,b1)
   ! This routine returns numberOfPolarizationVectors (no need to pass)
   call getGDMParameters( grid,alphaP,gdmPar,numberOfPolarizationVectors, maxNumberOfParameters,maxNumberOfPolarizationVectors,gdmParOption )

   if( alphaP.ne.0 .and. abs(eps*alphaP-1.) .gt. 1.e-10 )then
     write(*,'(" advOptNew: ERROR alphaP != 1/eps ")') 
     stop 2288
   end if

   if( t.eq.0. .and. dispersionModel.ne.noDispersion .and. debug.gt.1 )then
     ! ---- Dispersive Maxwell ----
     write(*,'("--advOptNew-- dispersionModel =",i4," px,py,pz =",3i3)') dispersionModel,pxc,pyc,pzc
     write(*,'("--advOptNew-- GDM: numberOfPolarizationVectors =",i4," alphaP =",e8.2)') numberOfPolarizationVectors,alphaP
     write(*,'("--advOptNew-- GDM: alphaP,a0,a1,b0,b1 =",5(1p,e10.2))') alphaP,a0,a1,b0,b1

     do iv=0,numberOfPolarizationVectors-1
       write(*,'("--advOptNew-- GDM: eqn=",i3," a0,a1,b0,b1 =",4(1p,e10.2))') iv,a0v(iv),a1v(iv),b0v(iv),b1v(iv)
     end do
  end if
 end if

 
 if( nonlinearModel .ne. noNonlinearModel )then
   write(*,'("--advOptNew-- nonlinearModel =",i4,"(1=multilevelAtomic)")') nonlinearModel
   call getMultilevelAtomicParameters( grid, nlPar, maxPar, maxPar, numberOfPolarizationVectors, numberOfAtomicLevels )

   write(*,'("multilevelAtomic: numberOfPolarizationVectors =",i4,"  numberOfAtomicLevels =",i4)') numberOfPolarizationVectors, numberOfAtomicLevels
   write(*,'("polarizationNECoefficients:")')
   do m1=0,numberOfPolarizationVectors-1
     write(*,'( 10(e12.3,1x) )') (pnec(m1,m2),m2=0,numberOfAtomicLevels-1)
   end do 

   write(*,'("populationRelaxationCoefficients:")')
   do m1=0,numberOfAtomicLevels-1
     write(*,'( 10(e12.3,1x) )') (prc(m1,m2),m2=0,numberOfAtomicLevels-1)
   end do 

   write(*,'("populationEPtCoefficients:")')
   do m1=0,numberOfAtomicLevels-1
     write(*,'( 10(e12.3,1x) )') (peptc(m1,m2),m2=0,numberOfPolarizationVectors-1)
   end do 

   
 end if 

 if( useSosupDissipation.ne.0 )then

  ! Coefficients in the sosup dissipation from Jordan Angel
  if( orderOfAccuracy.eq.2 )then
   adSosup=cc*dt*1./8.
  else if( orderOfAccuracy.eq.4 )then
    adSosup=cc*dt*5./288.
  else
    stop 1005
  end if

  uDotFactor=.5  ! By default uDot is D-zero and so we scale (un-um) by .5 --> .5*(un-um)/(dt)

  ! sosupParameter=gamma in sosup scheme  0<= gamma <=1   0=centered scheme
  adSosup=sosupParameter*adSosup

  if( t.le.2*dt .and. debug.gt.1 )then
    write(*,'("advOPT: useSosup dissipation, t,dt,adSosup=",3e10.2)') t,dt,adSosup
    write(*,'("advOPT: sosupDissipationOption=",i2)') sosupDissipationOption
    write(*,'("advOPT: updateDissipation=",i2)') updateDissipation
    write(*,'("advOPT: updateSolution=",i2)') updateSolution
    write(*,'("advOPT: useNewForcingMethod=",i2)') useNewForcingMethod
  end if
  ! Coefficients of the sosup dissipation with Cartesian grids:
  cdSosupx= adSosup/dx(0)
  cdSosupy= adSosup/dx(1)
  cdSosupz= adSosup/dx(2)


 end if

 if( useDivergenceCleaning.eq.1 )then
   ! Here are the coefficients that define the div cleaning formulae
   !    D+tD-t( E ) + alpha*( D0t E ) = c^2 Delta(E) + alpha*( (1/eps) Curl ( H ) )
   if( orderOfAccuracy.eq.2 )then
     ! 2D, 2nd-order, div cleaning:
     dc = divergenceCleaningCoefficient
     dcp = 1. + dc*dt*.5
     cdc0 = 2./dcp
     cdc1 = -(1.-dc*dt*.5)/dcp

     cdcxx = (cc*dt)**2/(dx(0)**2)/dcp
     cdcyy = (cc*dt)**2/(dx(1)**2)/dcp
     cdczz = (cc*dt)**2/(dx(2)**2)/dcp
     ! for div(H) damping in E eqn:
     cdcEdx= dc*dt**2/(eps*2.*dx(0))/dcp
     cdcEdy= dc*dt**2/(eps*2.*dx(1))/dcp
     cdcEdz= dc*dt**2/(eps*2.*dx(2))/dcp
     ! for div(E) damping in H eqn:
     cdcHdx= dc*dt**2/(mu*2.*dx(0))/dcp
     cdcHdy= dc*dt**2/(mu*2.*dx(1))/dcp
     cdcHdz= dc*dt**2/(mu*2.*dx(2))/dcp

     ! These next two are for the curvilinear case:
     cdcf = (cc*dt)**2/dcp
     cdcE = dc*dt**2/(eps)/dcp
     cdcH = dc*dt**2/(mu )/dcp

     if( t.eq.0. )then
       write(*,'(" advOpt: order=2 : div clean: dc,cc,dt,eps,mu=",5e10.2)') dc,cc,dt,eps,mu
       write(*,'(" advOpt: div clean: cdc0,cdc1,cdcxx,cdcyy,cdcHdy,cdcHdx=",6e10.2)') cdc0,cdc1,cdcxx,cdcyy,cdcHdy,cdcHdx
     end if
   else if( orderOfAccuracy.eq.4 )then

     dc = divergenceCleaningCoefficient
     dcp = 1. + dc*dt*.5
     cdc0 = 2./dcp
     cdc1 = -(1.-dc*dt*.5)/dcp

     cdcE= dc*dt**2/(eps)/dcp
     cdcELap= ((cc*dt)**2/dcp)*( 1. + dc*dt/(6.*eps) )
     cdcELapsq = ((cc*dt)**4/12./dcp)*( 1. + dc*dt/eps )
     cdcELapm = ((cc*dt)**2/dcp)*( - dc*dt/(6.*eps) )

     cdcH= dc*dt**2/(mu )/dcp
     cdcHLap= ((cc*dt)**2/dcp)*( 1. + dc*dt/(6.*mu) )
     cdcHLapsq = ((cc*dt)**4/12./dcp)*( 1. + dc*dt/mu )
     cdcHLapm = ((cc*dt)**2/dcp)*( - dc*dt/(6.*mu ) )

     if( t.eq.0. )then
       write(*,'(" advOpt: order=4 :  div clean: dc,cc,dt,eps,mu=",5e10.2)') dc,cc,dt,eps,mu
       write(*,'(" advOpt: div clean: cdc0,cdc1,cdcELap,cdcELapsq,cdcE,cdcELapm=",8e10.2)') cdc0,cdc1,cdcELap,cdcELapsq,cdcE,cdcELapm
     end if




   else
    write(*,'(" advOpt.bf: un-implemented orderOfAccuracy for div-cleaning")')
    stop 2277
   end if
 end if

 if( orderOfAccuracy.eq.6 )then
   if( nd.eq.2 )then
     c00lap2d6=csq*(-49./18.)*(1./dx(0)**2+1./dy**2)
     c10lap2d6=csq*(1.5     )*(1./dx(0)**2)
     c01lap2d6=csq*(1.5     )*(1./dy**2)
     c20lap2d6=csq*(-3./20. )*(1./dx(0)**2)
     c02lap2d6=csq*(-3./20. )*(1./dy**2)
     c30lap2d6=csq*(1./90.  )*(1./dx(0)**2)
     c03lap2d6=csq*(1./90.  )*(1./dy**2)
   else
     c000lap3d6=csq*(-49./18.)*(1./dx(0)**2+1./dy**2+1./dz**2)
     c100lap3d6=csq*(1.5     )*(1./dx(0)**2)
     c010lap3d6=csq*(1.5     )*(1./dy**2)
     c001lap3d6=csq*(1.5     )*(1./dz**2)
     c200lap3d6=csq*(-3./20. )*(1./dx(0)**2)
     c020lap3d6=csq*(-3./20. )*(1./dy**2)
     c002lap3d6=csq*(-3./20. )*(1./dz**2)
     c300lap3d6=csq*(1./90.  )*(1./dx(0)**2)
     c030lap3d6=csq*(1./90.  )*(1./dy**2)
     c003lap3d6=csq*(1./90.  )*(1./dz**2)
   end if
 end if
 if( orderOfAccuracy.eq.8 )then
   if( nd.eq.2 )then
     c00lap2d8=csq*(-205./72.)*(1./dx(0)**2+1./dy**2)
     c10lap2d8=csq*(8./5.    )*(1./dx(0)**2)
     c01lap2d8=csq*(8./5.    )*(1./dy**2)
     c20lap2d8=csq*(-1./5.   )*(1./dx(0)**2)
     c02lap2d8=csq*(-1./5.   )*(1./dy**2)
     c30lap2d8=csq*(8./315.  )*(1./dx(0)**2)
     c03lap2d8=csq*(8./315.  )*(1./dy**2)
     c40lap2d8=csq*(-1./560. )*(1./dx(0)**2)
     c04lap2d8=csq*(-1./560. )*(1./dy**2)
   else
     c000lap3d8=csq*(-205./72.)*(1./dx(0)**2+1./dy**2+1./dz**2)
     c100lap3d8=csq*(8./5.    )*(1./dx(0)**2)
     c010lap3d8=csq*(8./5.    )*(1./dy**2)
     c001lap3d8=csq*(8./5.    )*(1./dz**2)
     c200lap3d8=csq*(-1./5.   )*(1./dx(0)**2)
     c020lap3d8=csq*(-1./5.   )*(1./dy**2)
     c002lap3d8=csq*(-1./5.   )*(1./dz**2)
     c300lap3d8=csq*(8./315.  )*(1./dx(0)**2)
     c030lap3d8=csq*(8./315.  )*(1./dy**2)
     c003lap3d8=csq*(8./315.  )*(1./dz**2)
     c400lap3d8=csq*(-1./560. )*(1./dx(0)**2)
     c040lap3d8=csq*(-1./560. )*(1./dy**2)
     c004lap3d8=csq*(-1./560. )*(1./dz**2)
   end if
 end if

! ! For stoermer: -- no longer used
! if( orderInTime.eq.4 )then
!   c40=( 7./6. )*dtsq
!   c41=(-5./12.)*dtsq
!   c42=( 1./3. )*dtsq
!   c43=(-1./12.)*dtsq
! else if( orderInTime.eq.6 )then
!   c60=( 317./240.)*dtsq    ! from stoermer.maple
!   c61=(-266./240.)*dtsq
!   c62=( 374./240.)*dtsq
!   c63=(-276./240.)*dtsq
!   c64=( 109./240.)*dtsq
!   c65=( -18./240.)*dtsq
! else if( orderInTime.eq.8 )then
!
!!     g := 1/60480 (236568 fv[4] + 88324 fv[0] - 121797 fv[1] + 245598 fv[2]
!!     + 33190 fv[6] - 4125 fv[7] - 300227 fv[3] - 117051 fv[5])
!
!   c80=(  88324./60480.)*dtsq ! from stoermer.maple
!   c81=(-121797./60480.)*dtsq
!   c82=( 245598./60480.)*dtsq
!   c83=(-300227./60480.)*dtsq
!   c84=( 236568./60480.)*dtsq
!   c85=(-117051./60480.)*dtsq
!   c86=(  33190./60480.)*dtsq
!   c87=(  -4125./60480.)*dtsq
! end if


 if( computeUt.eq.1 .and. updateDissipation.eq.1 )then
   ! precompute "uDot" = dt*du/dt used in the dissipation and store in v
   ! we uDot at enough ghost points for the dissipation operator
   if( t.le.3.*dt )then
     write(*,'(" advOPT>>> Eval uDot...")')
   end if
   numGhost=orderOfAccuracy/2
   if( useSosupDissipation.eq.1 )then
     numGhost=numGhost+1
   end if
   m1a=n1a-numGhost
   m1b=n1b+numGhost
   m2a=n2a-numGhost
   m2b=n2b+numGhost
   if( nd.eq.2 )then
    m3a=n3a
    m3b=n3b
   else
     m3a=n3a-numGhost
     m3b=n3b+numGhost
   end if
   nStart=ex
   if( nd.eq.2 )then
      nEnd=hz
   else
      nEnd=ez
   end if
   ! Use Dot( un )
   do n=nStart,nEnd
     beginLoopsMask(i1,i2,i3,m1a,m1b,m2a,m2b,m3a,m3b)
       v(i1,i2,i3,n)=un(i1,i2,i3,n)-um(i1,i2,i3,n)
     endLoopsMask()
   end do
 endif

  ! This next function will:
  !   (1) optionally compute the dissipation and fill in the diss array
  !            if: (adc.gt.0. .and. combineDissipationWithAdvance.eq.0
  !   (2) add the divergence damping
  !         if( add.gt.0. )
 if( nd.eq.2 .and. orderOfAccuracy.eq.2 )then
   call advMxDiss2dOrder2(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,\
     nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,  um,u,un,f, v,\
     pm,p,pn,qm,q,qn, bc, dis, varDis, ipar, rpar, ierr )
 else if(  nd.eq.2 .and. orderOfAccuracy.eq.4 )then
   call advMxDiss2dOrder4(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,\
     nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,  um,u,un,f, v,\
     pm,p,pn,qm,q,qn, bc, dis, varDis, ipar, rpar, ierr )
 else if( nd.eq.3 .and. orderOfAccuracy.eq.2 )then
   call advMxDiss3dOrder2(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,\
     nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,  um,u,un,f, v,\
     pm,p,pn,qm,q,qn, bc, dis, varDis, ipar, rpar, ierr )
 else if(  nd.eq.3 .and. orderOfAccuracy.eq.4 )then
   call advMxDiss3dOrder4(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,\
     nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,  um,u,un,f, v,\
     pm,p,pn,qm,q,qn, bc, dis, varDis, ipar, rpar, ierr )
 else
   if( (adc.gt.0. .and. combineDissipationWithAdvance.eq.0) .or. add.gt.0. )then
     stop 1116
   end if
 end if


 if( option.eq.1 ) then
   return
 end if


! write(*,'(" advMaxwell: timeSteppingMethod=",i2)') timeSteppingMethod
 if( timeSteppingMethod.eq.defaultTimeStepping )then
  write(*,'(" advMaxwell:ERROR: timeSteppingMethod=defaultTimeStepping -- this should be set")')
    ! '
  stop 83322
 end if

 if( dispersionModel.ne.noDispersion .and. useConservative.eq.1 )then
   write(*,'("advOpt:ERROR: useConservative not implemented for dispersion model")') 
   stop 1213
 end if

!      *********************************************************
!      *************** Dispersive Update Here ******************
!      *********************************************************
! if( dispersionModel.ne.noDispersion )then

!   updateDispersive(DIM,ORDER,GRIDTYPE)

! end if


 if( gridType.eq.rectangular )then
  ! write(*,*) 'Inside advMaxwell rectangular marker 1...'

 #If #GRIDTYPE eq "rectangular"

!       **********************************************
!       *************** rectangular ******************
!       **********************************************
   ! write(*,*) 'Inside advMaxwell rectangular marker 2...'

 #If #ORDER eq "2"

   if( t.le.3*dt .and. debug>3 )then
     write(*,*) 'Inside advOptNew rectangular order=2...'
   end if
   #If #DIM eq "2"
    if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.noNonlinearModel )then
      ! --- linear dispersion model --- 

      updateDispersive(2,2,rectangular)

    else if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.multilevelAtomic )then

      ! --- Multilevel Atomic (Maxwell-Bloch) nonlinear model --- 

      updateMultilevelAtomic(2,2,rectangular)


    else if( useDivergenceCleaning.eq.0 )then

     ! FD22 with no dissipation
     loopsF2DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
              un(i1,i2,i3,ex)=maxwell2dr(i1,i2,i3,ex),\
              un(i1,i2,i3,ey)=maxwell2dr(i1,i2,i3,ey),\
              un(i1,i2,i3,hz)=maxwell2dr(i1,i2,i3,hz),,,,,,)
    else
     ! 2D, 2nd-order, div cleaning:
     !    D+tD-t( E ) + alpha*( D0t E ) = c^2 Delta(E) + alpha*( (1/eps) Curl ( H ) )
     write(*,'("advMaxwell: advance 2D, 2nd-order, rectangular, div cleaning... t=",e10.2)') t
     loopsF2DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
              un(i1,i2,i3,ex)=mxdc2d2Ex(i1,i2,i3),\
              un(i1,i2,i3,ey)=mxdc2d2Ey(i1,i2,i3),\
              un(i1,i2,i3,hz)=maxwell2dr(i1,i2,i3,hz),,,,,,)
    endif

   #Else
    ! ****** RECTANGULAR THREE DIMENSIONS SECOND-ORDER ****

    if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.noNonlinearModel )then
      ! --dispersion model --

      ! updateRectangular3dOrder2Dispersive()

      updateDispersive(3,2,rectangular)

    else if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.multilevelAtomic )then

      ! --- Multilevel Atomic (Maxwell-Bloch) nonlinear model --- 

      updateMultilevelAtomic(3,2,rectangular)

    else if( useDivergenceCleaning.eq.0 )then
     loopsF3DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,ez),\
              un(i1,i2,i3,ex)=maxwell3dr(i1,i2,i3,ex),\
              un(i1,i2,i3,ey)=maxwell3dr(i1,i2,i3,ey),\
              un(i1,i2,i3,ez)=maxwell3dr(i1,i2,i3,ez),,,,,,,\
              dtsq*f(i1,i2,i3,hx),dtsq*f(i1,i2,i3,hy),dtsq*f(i1,i2,i3,hz),\
              un(i1,i2,i3,hx)=maxwell3dr(i1,i2,i3,hx),\
              un(i1,i2,i3,hy)=maxwell3dr(i1,i2,i3,hy),\
              un(i1,i2,i3,hz)=maxwell3dr(i1,i2,i3,hz),,,,,,)
    else
     ! 3D, 2nd-order, div cleaning:
     write(*,'("advMaxwell: advance 3D, 2nd-order, rectangular, div cleaning... t=",e10.2)') t
     loopsF3DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,ez),\
              un(i1,i2,i3,ex)=mxdc3d2Ex(i1,i2,i3),\
              un(i1,i2,i3,ey)=mxdc3d2Ey(i1,i2,i3),\
              un(i1,i2,i3,ez)=mxdc3d2Ez(i1,i2,i3),,,,,,,\
              dtsq*f(i1,i2,i3,hx),dtsq*f(i1,i2,i3,hy),dtsq*f(i1,i2,i3,hz),\
              un(i1,i2,i3,hx)=mxdc3d2Hx(i1,i2,i3),\
              un(i1,i2,i3,hy)=mxdc3d2Hy(i1,i2,i3),\
              un(i1,i2,i3,hz)=mxdc3d2Hz(i1,i2,i3),,,,,,)
    end if
   #End

   ! write(*,*) 'Finished w/ order 2...'

 #Elif #ORDER eq "4"

   ! ======================================================================================
   ! ==================   4th order in space and 4th order in time: =======================
   ! ======================================================================================
   if( t.le.3*dt .and. debug.gt.3 )then
     write(*,*) 'Inside advMaxwell order=4 YOU ARE HERE'
   end if

   if( timeSteppingMethod.eq.modifiedEquationTimeStepping )then

     #If #DIM eq "2"

      ! ------------------------------------------------------------------------------
      !    2D : 4th order modified equation (rectangular)
      ! ------------------------------------------------------------------------------

      if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.noNonlinearModel )then
        ! --dispersion model --

        ! updateRectangular2dOrder4Dispersive()

        updateDispersive(2,4,rectangular)

      else if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.multilevelAtomic )then

        ! --- Multilevel Atomic (Maxwell-Bloch) nonlinear model --- 

        updateMultilevelAtomic(2,4,rectangular)

      else if( useDivergenceCleaning.eq.0 )then

!-       if( useSosupDissipation.ne.0 )then
!-
!-         ! FD44 (rectangular grid) with Sosup dissipation (wide stencil dissiption)
!-         updateUpwindDissipationRectangular2dOrder4()
!-
!-       else 
         if( useNewForcingMethod.eq.1 ) then

         ! fix forcing for ME scheme to be 4th-order
         if( combineDissipationWithAdvance.eq.0 )then
          write(*,*) 'advOpt: 2d, rect, 4th-order fix-force modified equation'
          loopsF2DD(dtsq*f2drme44(i1,i2,i3,ex),dtsq*f2drme44(i1,i2,i3,ey),dtsq*f2drme44(i1,i2,i3,hz),\
                 un(i1,i2,i3,ex)=maxwell2dr44me(i1,i2,i3,ex),\
                 un(i1,i2,i3,ey)=maxwell2dr44me(i1,i2,i3,ey),\
                 un(i1,i2,i3,hz)=maxwell2dr44me(i1,i2,i3,hz),,,,,,)

         else
          ! modified equation and dissipation in one loop
          write(*,*) 'advOpt: 2d, rect, 4th-order, diss, fix-force modified equation'
          loopsF2D(dtsq*f2drme44(i1,i2,i3,ex),dtsq*f2drme44(i1,i2,i3,ey),dtsq*f2drme44(i1,i2,i3,hz),\
                   un(i1,i2,i3,ex)=maxwell2dr44me(i1,i2,i3,ex)+adcdt*fd42d(i1,i2,i3,ex),\
                   un(i1,i2,i3,ey)=maxwell2dr44me(i1,i2,i3,ey)+adcdt*fd42d(i1,i2,i3,ey),\
                   un(i1,i2,i3,hz)=maxwell2dr44me(i1,i2,i3,hz)+adcdt*fd42d(i1,i2,i3,hz),,,,,,)
         end if

       else if( combineDissipationWithAdvance.eq.0 )then
        ! write(*,*) 'advOpt: 2d, rect, modified equation'
        loopsF2DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
               un(i1,i2,i3,ex)=maxwell2dr44me(i1,i2,i3,ex),\
               un(i1,i2,i3,ey)=maxwell2dr44me(i1,i2,i3,ey),\
               un(i1,i2,i3,hz)=maxwell2dr44me(i1,i2,i3,hz),,,,,,)

       else
        ! modified equation and dissipation in one loop
        loopsF2D(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
                 un(i1,i2,i3,ex)=maxwell2dr44me(i1,i2,i3,ex)+adcdt*fd42d(i1,i2,i3,ex),\
                 un(i1,i2,i3,ey)=maxwell2dr44me(i1,i2,i3,ey)+adcdt*fd42d(i1,i2,i3,ey),\
                 un(i1,i2,i3,hz)=maxwell2dr44me(i1,i2,i3,hz)+adcdt*fd42d(i1,i2,i3,hz),,,,,,)
       end if
      else
       ! 2D, 4th-order, div cleaning:
       !    D+tD-t( E ) + alpha*( D0t E ) = c^2 Delta(E) + alpha*( (1/eps) Curl ( H ) )
       write(*,'("advMaxwell: advance 2D, 4th-order, rect, div cleaning... t=",e10.2,", adcdt=",e10.2 )') t,adcdt

       if( useNewForcingMethod.eq.1 ) then
         write(*,'(" advOpt: FINISH ME")')
         stop 10123
       end if

       if( combineDissipationWithAdvance.eq.0 )then
         loopsF2DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
                un(i1,i2,i3,ex)=mxdc2d4Ex(i1,i2,i3),\
                un(i1,i2,i3,ey)=mxdc2d4Ey(i1,i2,i3),\
                un(i1,i2,i3,hz)=maxwell2dr44me(i1,i2,i3,hz),,,,,,)
       else
         loopsF2DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
                un(i1,i2,i3,ex)=mxdc2d4Ex(i1,i2,i3)+adcdt*fd42d(i1,i2,i3,ex),\
                un(i1,i2,i3,ey)=mxdc2d4Ey(i1,i2,i3)+adcdt*fd42d(i1,i2,i3,ey),\
                un(i1,i2,i3,hz)=maxwell2dr44me(i1,i2,i3,hz)+adcdt*fd42d(i1,i2,i3,hz),,,,,,)
       end if

      end if

     #Else

      ! ------------------------------------------------------------------------------
      !    3D : 4th order modified equation (rectangular)
      ! ------------------------------------------------------------------------------
      ! write(*,*) 'Inside advMaxwell order=4, 3D NOW YOU ARE HERE'

      if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.noNonlinearModel )then
      ! --dispersion model --
        ! ZZZ
        ! updateRectangular3dOrder4Dispersive()

        updateDispersive(3,4,rectangular)

      else if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.multilevelAtomic )then

       ! --- Multilevel Atomic (Maxwell-Bloch) nonlinear model --- 

        updateMultilevelAtomic(3,4,rectangular)

      else if( useDivergenceCleaning.eq.0 )then

        if( useNewForcingMethod.eq.1 ) then

         if( combineDissipationWithAdvance.eq.0 )then
          ! 4th order modified equation
          loopsF3DD(dtsq*f3drme44(i1,i2,i3,ex),dtsq*f3drme44(i1,i2,i3,ey),dtsq*f3drme44(i1,i2,i3,ez),\
                un(i1,i2,i3,ex)=maxwell3dr44me(i1,i2,i3,ex),\
                un(i1,i2,i3,ey)=maxwell3dr44me(i1,i2,i3,ey),\
                un(i1,i2,i3,ez)=maxwell3dr44me(i1,i2,i3,ez),,,,,,,\
                dtsq*f(i1,i2,i3,hx),dtsq*f(i1,i2,i3,hy),dtsq*f(i1,i2,i3,hz),\
                un(i1,i2,i3,hx)=maxwell3dr44me(i1,i2,i3,hx),\
                un(i1,i2,i3,hy)=maxwell3dr44me(i1,i2,i3,hy),\
                un(i1,i2,i3,hz)=maxwell3dr44me(i1,i2,i3,hz),,,,,,)
          else
           ! 4th order modified equation and dissipation in one loop
           loopsF3D(dtsq*f3drme44(i1,i2,i3,ex),dtsq*f3drme44(i1,i2,i3,ey),dtsq*f3drme44(i1,i2,i3,ez),\
                un(i1,i2,i3,ex)=maxwell3dr44me(i1,i2,i3,ex)+adcdt*fd43d(i1,i2,i3,ex),\
                un(i1,i2,i3,ey)=maxwell3dr44me(i1,i2,i3,ey)+adcdt*fd43d(i1,i2,i3,ey),\
                un(i1,i2,i3,ez)=maxwell3dr44me(i1,i2,i3,ez)+adcdt*fd43d(i1,i2,i3,ez),,,,,,,\
                dtsq*f(i1,i2,i3,hx),dtsq*f(i1,i2,i3,hy),dtsq*f(i1,i2,i3,hz),\
                un(i1,i2,i3,hx)=maxwell3dr44me(i1,i2,i3,hx)+adcdt*fd43d(i1,i2,i3,hx),\
                un(i1,i2,i3,hy)=maxwell3dr44me(i1,i2,i3,hy)+adcdt*fd43d(i1,i2,i3,hy),\
                un(i1,i2,i3,hz)=maxwell3dr44me(i1,i2,i3,hz)+adcdt*fd43d(i1,i2,i3,hz),,,,,,)
          end if

        else

         if( combineDissipationWithAdvance.eq.0 )then
          ! 4th order modified equation
          loopsF3DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,ez),\
                un(i1,i2,i3,ex)=maxwell3dr44me(i1,i2,i3,ex),\
                un(i1,i2,i3,ey)=maxwell3dr44me(i1,i2,i3,ey),\
                un(i1,i2,i3,ez)=maxwell3dr44me(i1,i2,i3,ez),,,,,,,\
                dtsq*f(i1,i2,i3,hx),dtsq*f(i1,i2,i3,hy),dtsq*f(i1,i2,i3,hz),\
                un(i1,i2,i3,hx)=maxwell3dr44me(i1,i2,i3,hx),\
                un(i1,i2,i3,hy)=maxwell3dr44me(i1,i2,i3,hy),\
                un(i1,i2,i3,hz)=maxwell3dr44me(i1,i2,i3,hz),,,,,,)
          else
           ! 4th order modified equation and dissipation in one loop
           loopsF3D(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,ez),\
                un(i1,i2,i3,ex)=maxwell3dr44me(i1,i2,i3,ex)+adcdt*fd43d(i1,i2,i3,ex),\
                un(i1,i2,i3,ey)=maxwell3dr44me(i1,i2,i3,ey)+adcdt*fd43d(i1,i2,i3,ey),\
                un(i1,i2,i3,ez)=maxwell3dr44me(i1,i2,i3,ez)+adcdt*fd43d(i1,i2,i3,ez),,,,,,,\
                dtsq*f(i1,i2,i3,hx),dtsq*f(i1,i2,i3,hy),dtsq*f(i1,i2,i3,hz),\
                un(i1,i2,i3,hx)=maxwell3dr44me(i1,i2,i3,hx)+adcdt*fd43d(i1,i2,i3,hx),\
                un(i1,i2,i3,hy)=maxwell3dr44me(i1,i2,i3,hy)+adcdt*fd43d(i1,i2,i3,hy),\
                un(i1,i2,i3,hz)=maxwell3dr44me(i1,i2,i3,hz)+adcdt*fd43d(i1,i2,i3,hz),,,,,,)
          end if
         end if

      else
         ! -- div clean
        write(*,'("advMaxwell: advance 3D, 4th-order, rect, div cleaning... t=",e10.2,", adcdt=",e10.2 )') t,adcdt

        if( useNewForcingMethod.eq.1 ) then
          write(*,'(" advOpt: FINISH ME")')
          stop 10124
        end if

        if( combineDissipationWithAdvance.eq.0 )then
         ! 4th order modified equation
         loopsF3DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,ez),\
               un(i1,i2,i3,ex)=mxdc3d4Ex(i1,i2,i3),\
               un(i1,i2,i3,ey)=mxdc3d4Ey(i1,i2,i3),\
               un(i1,i2,i3,ez)=mxdc3d4Ez(i1,i2,i3),,,,,,,\
               dtsq*f(i1,i2,i3,hx),dtsq*f(i1,i2,i3,hy),dtsq*f(i1,i2,i3,hz),\
               un(i1,i2,i3,hx)=mxdc3d4Hx(i1,i2,i3),\
               un(i1,i2,i3,hy)=mxdc3d4Hy(i1,i2,i3),\
               un(i1,i2,i3,hz)=mxdc3d4Hz(i1,i2,i3),,,,,,)
         else
!         ! 4th order modified equation and dissipation in one loop
          loopsF3D(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,ez),\
               un(i1,i2,i3,ex)=mxdc3d4Ex(i1,i2,i3)+adcdt*fd43d(i1,i2,i3,ex),\
               un(i1,i2,i3,ey)=mxdc3d4Ey(i1,i2,i3)+adcdt*fd43d(i1,i2,i3,ey),\
               un(i1,i2,i3,ez)=mxdc3d4Ez(i1,i2,i3)+adcdt*fd43d(i1,i2,i3,ez),,,,,,,\
               dtsq*f(i1,i2,i3,hx),dtsq*f(i1,i2,i3,hy),dtsq*f(i1,i2,i3,hz),\
               un(i1,i2,i3,hx)=mxdc3d4Hx(i1,i2,i3)+adcdt*fd43d(i1,i2,i3,hx),\
               un(i1,i2,i3,hy)=mxdc3d4Hy(i1,i2,i3)+adcdt*fd43d(i1,i2,i3,hy),\
               un(i1,i2,i3,hz)=mxdc3d4Hz(i1,i2,i3)+adcdt*fd43d(i1,i2,i3,hz),,,,,,)
         end if
       end if

     #End

   else  ! not modified equation

     ! We no longer support Stoermer
     stop 4444

   end if

 #Elif #ORDER eq "6"
   ! *** else if( orderOfAccuracy.eq.6 .and. orderInTime.eq.6 )then

   ! 6th order in space and 6th order in time:
   ! write(*,*) 'Inside advMaxwell order=6...'
   if( timeSteppingMethod.eq.modifiedEquationTimeStepping )then

     #If #DIM eq "2"
       ! 6th order modified equation


       ! write(*,*) 'advOpt: 2d, rect, modified equation'
       loopsF2D(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
              un(i1,i2,i3,ex)=maxwell2dr66me(i1,i2,i3,ex),\
              un(i1,i2,i3,ey)=maxwell2dr66me(i1,i2,i3,ey),\
              un(i1,i2,i3,hz)=maxwell2dr66me(i1,i2,i3,hz),,,,,,)

     #Else
       ! 6th order modified equation
       loopsF3D(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,ez),\
              un(i1,i2,i3,ex)=maxwell3dr66me(i1,i2,i3,ex),\
              un(i1,i2,i3,ey)=maxwell3dr66me(i1,i2,i3,ey),\
              un(i1,i2,i3,ez)=maxwell3dr66me(i1,i2,i3,ez),,,,,,,\
              dtsq*f(i1,i2,i3,hx),dtsq*f(i1,i2,i3,hy),dtsq*f(i1,i2,i3,hz),\
              un(i1,i2,i3,hx)=maxwell3dr66me(i1,i2,i3,hx),\
              un(i1,i2,i3,hy)=maxwell3dr66me(i1,i2,i3,hy),\
              un(i1,i2,i3,hz)=maxwell3dr66me(i1,i2,i3,hz),,,,,,)

     #End

   else

    ! We no longer support Stoermer
     stop 4444

   end if

 #Elif #ORDER eq "8"
   ! *** else if( orderOfAccuracy.eq.8 .and. orderInTime.eq.8 )then

   ! 8th order in space and 8th order in time:
   ! write(*,*) 'Inside advMaxwell order=8...'
   if( timeSteppingMethod.eq.modifiedEquationTimeStepping )then

     #If #DIM eq "2"
       ! 8th order modified equation


       ! write(*,*) 'advOpt: 2d, rect, modified equation'
       loopsF2D(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
              un(i1,i2,i3,ex)=maxwell2dr88me(i1,i2,i3,ex),\
              un(i1,i2,i3,ey)=maxwell2dr88me(i1,i2,i3,ey),\
              un(i1,i2,i3,hz)=maxwell2dr88me(i1,i2,i3,hz),,,,,,)

     #Else
       ! 8th order modified equation
       loopsF3D(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,ez),\
              un(i1,i2,i3,ex)=maxwell3dr88me(i1,i2,i3,ex),\
              un(i1,i2,i3,ey)=maxwell3dr88me(i1,i2,i3,ey),\
              un(i1,i2,i3,ez)=maxwell3dr88me(i1,i2,i3,ez),,,,,,,\
              dtsq*f(i1,i2,i3,hx),dtsq*f(i1,i2,i3,hy),dtsq*f(i1,i2,i3,hz),\
              un(i1,i2,i3,hx)=maxwell3dr88me(i1,i2,i3,hx),\
              un(i1,i2,i3,hy)=maxwell3dr88me(i1,i2,i3,hy),\
              un(i1,i2,i3,hz)=maxwell3dr88me(i1,i2,i3,hz),,,,,,)

     #End

   else

    ! We no longer support Stoermer
     stop 4444

   end if

 #Else
   write(*,*) 'advMaxwell:ERROR orderOfAccuracy,orderInTime=',orderOfAccuracy,orderInTime
   stop 1

 #End

 #End

 else

 #If #GRIDTYPE eq "curvilinear"

!       **********************************************
!       *************** curvilinear ******************
!       **********************************************

   if( (useCurvilinearOpt.eq.1 .or. dispersionModel.ne.noDispersion) .and. useConservative.eq.0 )then

    ! ****************************************************************************
    ! *************** OPTIMIZED-CURVILINEAR AND NON-CONSERVATIVE *****************
    ! ****************************************************************************

    #If #ORDER eq "2"

     if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.noNonlinearModel )then

      #If #DIM eq "2"
       !updateCurvilinear2dOrder2Dispersive()

       updateDispersive(2,2,curvilinear)

      #Else
       ! updateCurvilinear3dOrder2Dispersive()

       updateDispersive(3,2,curvilinear)
      #End

     else if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.multilevelAtomic )then

      ! --- Multilevel Atomic (Maxwell-Bloch) nonlinear model --- 

      #If #DIM eq "2"
       
       updateMultilevelAtomic(2,2,curvilinear)

      #Else
       
       updateMultilevelAtomic(3,2,curvilinear)
      
      #End


     else

      ! --- Todo: non-conservative operators could be inlined here ---
      !   -- these might be faster than precomputing

        !$$$     loopsFCD(un(i1,i2,i3,ex)=maxwellc22(i1,i2,i3,ex),\
        !$$$              un(i1,i2,i3,ey)=maxwellc22(i1,i2,i3,ey),\
        !$$$              un(i1,i2,i3,ez)=maxwellc22(i1,i2,i3,ez),\
        !$$$               ,,,\
        !$$$              un(i1,i2,i3,hx)=maxwellc22(i1,i2,i3,hx),\
        !$$$              un(i1,i2,i3,hy)=maxwellc22(i1,i2,i3,hy),\
        !$$$              un(i1,i2,i3,hz)=maxwellc22(i1,i2,i3,hz),\
        !$$$              ,,)
        
        write(*,'(" advOptNew:ERROR: Curvilinear-optimized: useCurvilinearOpt=",i3," dispersionModel=",i2," useConservative=",i2)') useCurvilinearOpt,dispersionModel,useConservative
        write(*,'(" advOptNew:ERROR: FINISH ME")') 
        stop 88044
     end if

   #Elif #ORDER eq "4"

!-     if( useSosupDissipation.ne.0 .and. updateDissipation.eq.1 )then
!-
!-       ! ---- use sosup dissipation (wider stencil) ---
!-
!-      write(*,'(" finish me: FD44 non-cons && useSosupDissipation")')
!-      stop 4486
!-
!-
!-       if( useNewForcingMethod.ne.0 )then
!-        write(*,'(" finish me: dispersion && useNewForcingMethod")')
!-        stop 4487
!-       end if
!-
!-     else 
     if( timeSteppingMethod.eq.modifiedEquationTimeStepping )then

       ! ----------------------------------------------------------
       ! ----- 4th order in space and 4th order in time ------------
       ! ---- Modified equation, NON-CONSERVATIVE difference ----
       ! ----------------------------------------------------------


       !   cdtsq*uLaplacian42(i1,i2,i3,n)+cdtsq12*uLapSq(n)
       ! write(*,*) 'advOpt: 2d, curv, FULL modified equation'

      #If #DIM eq "2"

!$$$ loopsFCD2D($$getLapSq2dOrder2(),\
!$$$                  dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
!$$$                  un(i1,i2,i3,ex)=max2dc44me(i1,i2,i3,ex),\
!$$$                  un(i1,i2,i3,ey)=max2dc44me(i1,i2,i3,ey),\
!$$$                  un(i1,i2,i3,hz)=max2dc44me(i1,i2,i3,hz))

       ! do one at a time -- is this faster ? NO

!$$$ loopsFCD2DA($$evalLapSq2dOrder2(ex),\
!$$$             dtsq*f(i1,i2,i3,ex),\
!$$$             un(i1,i2,i3,ex)=max2dc44me(i1,i2,i3,ex))
!$$$
!$$$ loopsFCD2DA($$evalLapSq2dOrder2(ey),\
!$$$             dtsq*f(i1,i2,i3,ey),\
!$$$             un(i1,i2,i3,ey)=max2dc44me(i1,i2,i3,ey))
!$$$
!$$$ loopsFCD2DA($$evalLapSq2dOrder2(hz),\
!$$$             dtsq*f(i1,i2,i3,hz),\
!$$$             un(i1,i2,i3,hz)=max2dc44me(i1,i2,i3,hz))


        ! first evaluate Laplacian to 2nd-order
       ! *** need to evaluate on one additional line ***
       n1a=n1a-1
       n1b=n1b+1
       n2a=n2a-1
       n2b=n2b+1
       ! **** for this first loop we cannot use the mask --
       useWhereMaskSave=useWhereMask
       useWhereMask=0
       loopse9(v(i1,i2,i3,ex)=uLaplacian22(i1,i2,i3,ex),\
               v(i1,i2,i3,ey)=uLaplacian22(i1,i2,i3,ey),\
               v(i1,i2,i3,hz)=uLaplacian22(i1,i2,i3,hz),,,,,,)

       ! write(*,*) 'advOpt: 2d, rect, modified equation'
       n1a=n1a+1
       n1b=n1b-1
       n2a=n2a+1
       n2b=n2b-1
       useWhereMask=useWhereMaskSave

       if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.noNonlinearModel )then

         ! --dispersive model --
         ! updateCurvilinear2dOrder4Dispersive()

         updateDispersive(2,4,curvilinear)

       else if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.multilevelAtomic )then

        ! --- Multilevel Atomic (Maxwell-Bloch) nonlinear model --- 
       
        updateMultilevelAtomic(2,4,curvilinear)

       else if( useDivergenceCleaning.eq.0 )then

        if( useNewForcingMethod.eq.1 ) then
          ! fix forcing for ME scheme to be 4th-order
          loopsF2DD(dtsq*f2dcme44(i1,i2,i3,ex),dtsq*f2dcme44(i1,i2,i3,ey),dtsq*f2dcme44(i1,i2,i3,hz),\
                un(i1,i2,i3,ex)=max2dc44me2(i1,i2,i3,ex),\
                un(i1,i2,i3,ey)=max2dc44me2(i1,i2,i3,ey),\
                un(i1,i2,i3,hz)=max2dc44me2(i1,i2,i3,hz),,,,,,)
        else
          loopsF2DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
                un(i1,i2,i3,ex)=max2dc44me2(i1,i2,i3,ex),\
                un(i1,i2,i3,ey)=max2dc44me2(i1,i2,i3,ey),\
                un(i1,i2,i3,hz)=max2dc44me2(i1,i2,i3,hz),,,,,,)
        end if
       else
        if( useNewForcingMethod.eq.1 ) then
          write(*,'(" advOpt: FINISH ME")')
          stop 10134
        end if
        loopsF2DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,hz),\
              un(i1,i2,i3,ex)=mxdc2d4cEx(i1,i2,i3),\
              un(i1,i2,i3,ey)=mxdc2d4cEy(i1,i2,i3),\
              un(i1,i2,i3,hz)=max2dc44me2(i1,i2,i3,hz),,,,,,)
       end if


      #Elif #DIM == "3"

       ! *** need to evaluate on one additional line ***
       n1a=n1a-1
       n1b=n1b+1
       n2a=n2a-1
       n2b=n2b+1
       n3a=n3a-1
       n3b=n3b+1
       ! **** for this first loop we cannot use the mask --
       useWhereMaskSave=useWhereMask
       useWhereMask=0
       if( solveForE.ne.0 .and. solveForH.ne.0 )then
         stop 6666
!$$$        loopse9(v(i1,i2,i3,ex)=uLaplacian23(i1,i2,i3,ex),\
!$$$                v(i1,i2,i3,ey)=uLaplacian23(i1,i2,i3,ey),\
!$$$                v(i1,i2,i3,ez)=uLaplacian23(i1,i2,i3,ez),\
!$$$                v(i1,i2,i3,hx)=uLaplacian23(i1,i2,i3,hx),\
!$$$                v(i1,i2,i3,hy)=uLaplacian23(i1,i2,i3,hy),\
!$$$                v(i1,i2,i3,hz)=uLaplacian23(i1,i2,i3,hz),,,)
       else if( solveForE.ne.0 )then
        loopse9(v(i1,i2,i3,ex)=uLaplacian23(i1,i2,i3,ex),\
                v(i1,i2,i3,ey)=uLaplacian23(i1,i2,i3,ey),\
                v(i1,i2,i3,ez)=uLaplacian23(i1,i2,i3,ez),,,,,,)
       else
!$$$        loopse9(v(i1,i2,i3,hx)=uLaplacian23(i1,i2,i3,hx),\
!$$$                v(i1,i2,i3,hy)=uLaplacian23(i1,i2,i3,hy),\
!$$$                v(i1,i2,i3,hz)=uLaplacian23(i1,i2,i3,hz),,,,,,)
       end if

       ! write(*,*) 'advOpt: 2d, rect, modified equation'
       n1a=n1a+1
       n1b=n1b-1
       n2a=n2a+1
       n2b=n2b-1
       n3a=n3a+1
       n3b=n3b-1
       useWhereMask=useWhereMaskSave
       ! 4th order modified equation

       if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.noNonlinearModel )then
       ! --dispersion model --

         ! updateCurvilinear3dOrder4Dispersive()

         updateDispersive(3,4,curvilinear)

       else if( dispersionModel.ne.noDispersion .and. nonlinearModel.eq.multilevelAtomic )then

        ! --- Multilevel Atomic (Maxwell-Bloch) nonlinear model --- 
       
        updateMultilevelAtomic(3,4,curvilinear)

       else if( useDivergenceCleaning.eq.0 )then
         loopsF3DD(dtsq*f(i1,i2,i3,ex),dtsq*f(i1,i2,i3,ey),dtsq*f(i1,i2,i3,ez),\
              un(i1,i2,i3,ex)=max3dc44me(i1,i2,i3,ex),\
              un(i1,i2,i3,ey)=max3dc44me(i1,i2,i3,ey),\
              un(i1,i2,i3,ez)=max3dc44me(i1,i2,i3,ez),,,,,,,\
              dtsq*f(i1,i2,i3,hx),dtsq*f(i1,i2,i3,hy),dtsq*f(i1,i2,i3,hz),\
              un(i1,i2,i3,hx)=max3dc44me(i1,i2,i3,hx),\
              un(i1,i2,i3,hy)=max3dc44me(i1,i2,i3,hy),\
              un(i1,i2,i3,hz)=max3dc44me(i1,i2,i3,hz),,,,,,)
       else
         ! finish me
         stop 1005
       end if
      #End

     else

      stop 22743

     end if

   #Else
     write(*,'("MX: advOpt: curv-grid order=ORDER not implemeted")')
     stop 11155
   #End


   else if( useCurvilinearOpt.eq.1 .and. useConservative.eq.1 )then

    ! *************** conservative *****************

    stop 94422



   else

     ! **********************************************************************************
     ! **************** USE PRE-COMPUTED SPATIAL OPERATORS ******************************
     ! **********************************************************************************

     !  --> The Laplacian and Laplacian squared have already been computed by the calling program
     !  --> For example, mainly when using conservative operators


   #If #ORDER eq "2"

!-    if( useSosupDissipation.ne.0 )then
!-
!-      ! ---- use sosup dissipation (wider stencil) ---
!-      #If #DIM eq "2"
!-        updateUpwindDissipationCurvilinear2dOrder2()
!-      #Else
!-        write(*,'(" PRE-COMPUTED SPATIAL OPERATORS + GDM Not available")' )
!-        stop 3313
!-      #End
!-
!-    else 
    if( dispersionModel.ne.noDispersion )then

      stop 8843

      ! -- dispersive model --
      #If #DIM eq "2"
       ! updateCurvilinear2dOrder2Dispersive()

       ! updateDispersive(2,2,curvilinear)

      #Else
       ! updateCurvilinear3dOrder2Dispersive()

       ! updateDispersive(3,2,curvilinear)
      #End

    else if( useDivergenceCleaning.eq.0 )then

     ! --- currently 2nd-order conservative and non-conservative opertaors are done here ---
     ! --- non-dispersive ---
     ! ---- THIS IS FOR 2D OR 3D -----
     loopsFCD(un(i1,i2,i3,ex)=maxwellc22(i1,i2,i3,ex),\
              un(i1,i2,i3,ey)=maxwellc22(i1,i2,i3,ey),\
              un(i1,i2,i3,ez)=maxwellc22(i1,i2,i3,ez),\
              ,,,\
              un(i1,i2,i3,hx)=maxwellc22(i1,i2,i3,hx),\
              un(i1,i2,i3,hy)=maxwellc22(i1,i2,i3,hy),\
              un(i1,i2,i3,hz)=maxwellc22(i1,i2,i3,hz),\
              ,,)
    else
       ! 2D, 2nd-order, curvilinear, div cleaning:
       !    D+tD-t( E ) + alpha*( D0t E ) = c^2 Delta(E) + alpha*( (1/eps) Curl ( H ) )
     #If #DIM eq "2"
      write(*,'("advMaxwell: 2D, 2nd-order, curv, div cleaning... t=",e10.2,", adcdt=",e10.2 )') t,adcdt
      loopsFCD(un(i1,i2,i3,ex)=mxdc2d2cEx(i1,i2,i3),\
               un(i1,i2,i3,ey)=mxdc2d2cEy(i1,i2,i3),\
               un(i1,i2,i3,ez)=maxwellc22(i1,i2,i3,ez),\
               ,,,\
               un(i1,i2,i3,hx)=maxwellc22(i1,i2,i3,hx),\
               un(i1,i2,i3,hy)=maxwellc22(i1,i2,i3,hy),\
               un(i1,i2,i3,hz)=maxwellc22(i1,i2,i3,hz),\
               ,,)

     #Else
      write(*,'("advMaxwell: 3D, 2nd-order, curv, div cleaning... t=",e10.2,", adcdt=",e10.2 )') t,adcdt
      loopsFCD(un(i1,i2,i3,ex)=mxdc3d2cEx(i1,i2,i3),\
               un(i1,i2,i3,ey)=mxdc3d2cEy(i1,i2,i3),\
               un(i1,i2,i3,ez)=mxdc3d2cEz(i1,i2,i3),\
               ,,,\
               un(i1,i2,i3,hx)=mxdc3d2cHx(i1,i2,i3),\
               un(i1,i2,i3,hy)=mxdc3d2cHy(i1,i2,i3),\
               un(i1,i2,i3,hz)=mxdc3d2cHz(i1,i2,i3),\
               ,,)
     #End

    end if

   #Elif #ORDER eq "4"

!-     if( useSosupDissipation.ne.0 )then
!-
!-       ! ---- use sosup dissipation (wider stencil) ---
!-       updateUpwindDissipationCurvilinear2dOrder4()
!-
!-
!-     else 
     if( dispersionModel.ne.noDispersion )then

      ! --dispersive model --
      stop 7777

     else if( timeSteppingMethod.eq.modifiedEquationTimeStepping )then

       ! ******* 4th order in space and 4th order in time ********
       ! **************** CURVILINEAR *****************************
       ! **************** CONSERVATIVE DIFFERENCE *****************

       ! write(*,*) 'advOpt: 2d, curv, modified equation'

       if( useDivergenceCleaning.eq.0 )then
        if( useNewForcingMethod.eq.1 ) then

         if( addForcing.eq.0 )then
          ! NOTE: in 2D the ez,hx and hy entries are ignored
          loopsFCD(un(i1,i2,i3,ex)=maxwellc44me(i1,i2,i3,ex),\
                   un(i1,i2,i3,ey)=maxwellc44me(i1,i2,i3,ey),\
                   un(i1,i2,i3,ez)=maxwellc44me(i1,i2,i3,ez),,,,\
                   un(i1,i2,i3,hx)=maxwellc44me(i1,i2,i3,hx),\
                   un(i1,i2,i3,hy)=maxwellc44me(i1,i2,i3,hy),\
                   un(i1,i2,i3,hz)=maxwellc44me(i1,i2,i3,hz),,,)
         else
          ! Add forcing to ME scheme
          if( nd.eq.2 )then
           loopsFCD2DF(un(i1,i2,i3,ex)=maxwellc44me(i1,i2,i3,ex)+dtsq*f2dcme44(i1,i2,i3,ex),\
                       un(i1,i2,i3,ey)=maxwellc44me(i1,i2,i3,ey)+dtsq*f2dcme44(i1,i2,i3,ey),\
                       un(i1,i2,i3,hz)=maxwellc44me(i1,i2,i3,hz)+dtsq*f2dcme44(i1,i2,i3,hz))
          else
           loopsFCD3DF(un(i1,i2,i3,ex)=maxwellc44me(i1,i2,i3,ex)+dtsq*f3dcme44(i1,i2,i3,ex),\
                       un(i1,i2,i3,ey)=maxwellc44me(i1,i2,i3,ey)+dtsq*f3dcme44(i1,i2,i3,ey),\
                       un(i1,i2,i3,ez)=maxwellc44me(i1,i2,i3,ez)+dtsq*f3dcme44(i1,i2,i3,ez),,,,\
                       un(i1,i2,i3,hx)=maxwellc44me(i1,i2,i3,hx)+dtsq*f3dcme44(i1,i2,i3,hx),\
                       un(i1,i2,i3,hy)=maxwellc44me(i1,i2,i3,hy)+dtsq*f3dcme44(i1,i2,i3,hy),\
                       un(i1,i2,i3,hz)=maxwellc44me(i1,i2,i3,hz)+dtsq*f3dcme44(i1,i2,i3,hz),,,)
          end if

         end if
        else
         loopsFCD(un(i1,i2,i3,ex)=maxwellc44me(i1,i2,i3,ex),\
                  un(i1,i2,i3,ey)=maxwellc44me(i1,i2,i3,ey),\
                  un(i1,i2,i3,ez)=maxwellc44me(i1,i2,i3,ez),,,,\
                  un(i1,i2,i3,hx)=maxwellc44me(i1,i2,i3,hx),\
                  un(i1,i2,i3,hy)=maxwellc44me(i1,i2,i3,hy),\
                  un(i1,i2,i3,hz)=maxwellc44me(i1,i2,i3,hz),,,)
        end if
       else
        ! 2D, 4th-order, curvilinear, div cleaning:
        write(*,'("advMaxwell: 2D, 4th-order, curv, cons, div cleaning... t=",e10.2 )') t

        if( useNewForcingMethod.eq.1 ) then
          write(*,'(" advOpt: FINISH ME")')
          stop 10138
        end if

        #If #DIM eq "2"
         loopsFCD(un(i1,i2,i3,ex)=mxdc2d4cConsEx(i1,i2,i3),\
                  un(i1,i2,i3,ey)=mxdc2d4cConsEy(i1,i2,i3),\
                  un(i1,i2,i3,ez)=maxwellc44me(i1,i2,i3,ez),,,,\
                  un(i1,i2,i3,hx)=maxwellc44me(i1,i2,i3,hx),\
                  un(i1,i2,i3,hy)=maxwellc44me(i1,i2,i3,hy),\
                  un(i1,i2,i3,hz)=maxwellc44me(i1,i2,i3,hz),,,)
        #Else
          stop 4481
        #End
       end if
     else
       ! write(*,*) 'Inside advMaxwell curv, order=4...'

      ! We no longer support Stoermer
      stop 4444


     end if

   #Elif #ORDER eq "6"

     ! 6th order in space and 6th order in time:
     ! write(*,*) 'Inside advMaxwell order=6...'

     if( timeSteppingMethod.eq.modifiedEquationTimeStepping )then

       if( orderInTime.ne.2 )then
          write(*,'("MX: advOpt:ERROR curv-grid conservative orderInTime.ne.2")')
          stop 77155
       end if


       loopsFCD(un(i1,i2,i3,ex)=maxwellc66me(i1,i2,i3,ex),\
               un(i1,i2,i3,ey)=maxwellc66me(i1,i2,i3,ey),\
               un(i1,i2,i3,ez)=maxwellc66me(i1,i2,i3,ez),,,,\
               un(i1,i2,i3,hx)=maxwellc66me(i1,i2,i3,hx),\
               un(i1,i2,i3,hy)=maxwellc66me(i1,i2,i3,hy),\
               un(i1,i2,i3,hz)=maxwellc66me(i1,i2,i3,hz),,,)

     else
      ! We no longer support Stoermer
      stop 4444


     end if

   #Elif #ORDER eq "8"

     ! 8th order in space and 8th order in time:
     ! write(*,*) 'Inside advMaxwell order=8...'

     if( timeSteppingMethod.eq.modifiedEquationTimeStepping )then

       ! ** for now we just do 2nd-order in time **
       if( orderInTime.ne.2 )then
          stop 88188
       end if


       loopsFCD(un(i1,i2,i3,ex)=maxwellc88me(i1,i2,i3,ex),\
               un(i1,i2,i3,ey)=maxwellc88me(i1,i2,i3,ey),\
               un(i1,i2,i3,ez)=maxwellc88me(i1,i2,i3,ez),,,,\
               un(i1,i2,i3,hx)=maxwellc88me(i1,i2,i3,hx),\
               un(i1,i2,i3,hy)=maxwellc88me(i1,i2,i3,hy),\
               un(i1,i2,i3,hz)=maxwellc88me(i1,i2,i3,hz),,,)

     else

      ! We no longer support Stoermer
      stop 4444

     end if
  #Else
     write(*,*) 'advMaxwell:ERROR orderOfAccuracy,orderInTime=',orderOfAccuracy,orderInTime
     stop 2
  #End

  end if

 #End
 end if

 return
 end

#endMacro



#beginMacro buildFile(NAME,DIM,ORDER,GRIDTYPE)
#beginFile NAME.f
 ADV_MAXWELL(NAME,DIM,ORDER,GRIDTYPE)
#endFile
#endMacro

! ----- Some of these are still  built with advOpt.bf -------

      buildFile(advMx2dOrder2r,2,2,rectangular)
      buildFile(advMx3dOrder2r,3,2,rectangular)
!
      buildFile(advMx2dOrder2c,2,2,curvilinear)
      buildFile(advMx3dOrder2c,3,2,curvilinear)
!
      buildFile(advMx2dOrder4r,2,4,rectangular)
      buildFile(advMx3dOrder4r,3,4,rectangular)
!
      buildFile(advMx2dOrder4c,2,4,curvilinear)
      buildFile(advMx3dOrder4c,3,4,curvilinear)
!
!      buildFile(advMx2dOrder6r,2,6,rectangular)
!      buildFile(advMx3dOrder6r,3,6,rectangular)
!
!       ! build these for testing symmetric operators -- BC's not implemented yet
!      buildFile(advMx2dOrder6c,2,6,curvilinear)
!      buildFile(advMx3dOrder6c,3,6,curvilinear)
!
!      buildFile(advMx2dOrder8r,2,8,rectangular)
!      buildFile(advMx3dOrder8r,3,8,rectangular)
!
!       ! build these for testing symmetric operators -- BC's not implemented yet
!      buildFile(advMx2dOrder8c,2,8,curvilinear)
!      buildFile(advMx3dOrder8c,3,8,curvilinear)



! build an empty version of high order files so we do not have to compile the full version
#beginMacro ADV_MAXWELL_NULL(NAME,DIM,ORDER,GRIDTYPE)
 subroutine NAME(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                 mask,xy,rsxy,  um,u,un,f,fa, v, pm,p,pn, qm,q,qn, bc, dis, varDis, ipar, rpar, ierr )
!======================================================================
!   Advance a time step for Maxwells eqution
!     OPTIMIZED version for rectangular grids.
! nd : number of space dimensions
!
! ipar(0)  = option : option=0 - Maxwell+Artificial diffusion
!                           =1 - AD only
!
!  dis(i1,i2,i3) : temp space to hold artificial dissipation
!  varDis(i1,i2,i3) : coefficient of the variable artificial dissipation
!======================================================================
  write(*,'("ERROR: null version of NAME called")')
  stop 9922
  return
  end
#endMacro


#beginMacro buildFileNull(NAME,DIM,ORDER,GRIDTYPE)
#beginFile NAME ## Null.f
 ADV_MAXWELL_NULL(NAME,DIM,ORDER,GRIDTYPE)
#endFile
#endMacro

!      buildFileNull(advMx2dOrder6r,2,6,rectangular)
!      buildFileNull(advMx3dOrder6r,3,6,rectangular)
!
!      buildFileNull(advMx2dOrder6c,2,6,curvilinear)
!      buildFileNull(advMx3dOrder6c,3,6,curvilinear)
!
!      buildFileNull(advMx2dOrder8r,2,8,rectangular)
!      buildFileNull(advMx3dOrder8r,3,8,rectangular)
!
!      buildFileNull(advMx2dOrder8c,2,8,curvilinear)
!      buildFileNull(advMx3dOrder8c,3,8,curvilinear)



! ******* THIS IS NOT CURRENTLY USED -- see version in advOpt.bf *******************
      subroutine advMaxwellNew(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                            mask,xy,rx,  um,u,un,f,fa, v, vvt2,ut3,vvt4,ut5,ut6,ut7, bc, dis, varDis, ipar, rpar, ierr )
!======================================================================
!   Advance a time step for Maxwells eqution
!     OPTIMIZED version for rectangular grids.
! nd : number of space dimensions
!
! ipar(0)  = option : option=0 - Maxwell+Artificial diffusion
!                           =1 - AD only
!======================================================================
      implicit none
      integer nd, n1a,n1b,n2a,n2b,n3a,n3b,
     & nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b

      real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
      real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real vvt2(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real ut3(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real vvt4(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real ut5(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real ut6(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real ut7(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real dis(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      real varDis(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
      real rx(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)

      integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      integer bc(0:1,0:2),ierr

      integer ipar(0:*)
      real rpar(0:*)

!     ---- local variables -----
      integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderInTime
      integer addForcing,orderOfDissipation,option
      integer useWhereMask,solveForE,solveForH,grid
      integer ex,ey,ez, hx,hy,hz

      integer rectangular,curvilinear
      parameter( rectangular=0, curvilinear=1 )
!...........end   statement functions


      ! write(*,*) 'Inside advMaxwell...'

      orderOfAccuracy    =ipar(2)
      gridType           =ipar(1)

      if( orderOfAccuracy.eq.2 )then

        if( nd.eq.2 .and. gridType.eq.rectangular ) then
          call advMx2dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
          call advMx2dOrder2c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if( nd.eq.3 .and. gridType.eq.rectangular ) then
          call advMx3dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                             mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
          call advMx3dOrder2c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                             mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else
          stop 2271
        end if

      else if( orderOfAccuracy.eq.4 ) then
        if( nd.eq.2 .and. gridType.eq.rectangular )then
          call advMx2dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(nd.eq.2 .and. gridType.eq.curvilinear )then
          call advMx2dOrder4c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.rectangular )then
          call advMx3dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
          call advMx3dOrder4c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
       else
         stop 8843
       end if

!
      else if( orderOfAccuracy.eq.6 ) then
        if( nd.eq.2 .and. gridType.eq.rectangular )then
          call advMx2dOrder6r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(nd.eq.2 .and. gridType.eq.curvilinear )then
          call advMx2dOrder6c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.rectangular )then
          call advMx3dOrder6r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
          call advMx3dOrder6c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
       else
         stop 8843
       end if

      else if( orderOfAccuracy.eq.8 ) then

        if( nd.eq.2 .and. gridType.eq.rectangular )then
          call advMx2dOrder8r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(nd.eq.2 .and. gridType.eq.curvilinear )then
          call advMx2dOrder8c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.rectangular )then
          call advMx3dOrder8r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
          call advMx3dOrder8c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,xy,rx, um,u,un,f,fa, v,vvt2,ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
       else
         stop 8843
       end if

      else
        write(*,'(" advMaxwell:ERROR: un-implemented order of accuracy =",i6)') orderOfAccuracy
          ! '
        stop 11122
      end if

      return
      end
