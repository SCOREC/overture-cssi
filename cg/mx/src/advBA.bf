!
! Optimized routines for BA Maxwell
!

! These next include files will define the macros that will define the difference approximations
! The actual macro is called below
#Include "defineDiffOrder2f.h"
#Include "defineDiffOrder4f.h"


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

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (4th-order difference used with 2nd-order scheme) 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#defineMacro sosupDiss2d4(u,i1,i2,i3,n) \
              ( -6.*u(i1,i2,i3,n)     \
                +4.*(u(i1+1,i2,i3,n)+u(i1-1,i2,i3,n))     \
                   -(u(i1+2,i2,i3,n)+u(i1-2,i2,i3,n)) )*adxSosup(0) + \
              ( -6.*u(i1,i2,i3,n)     \
                +4.*(u(i1,i2+1,i3,n)+u(i1,i2-1,i3,n))     \
                   -(u(i1,i2+2,i3,n)+u(i1,i2-2,i3,n)) )*adxSosup(1)

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (3D - 4th-order difference used with 2nd-order scheme) 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#defineMacro sosupDiss3d4(u,i1,i2,i3,n) \
              ( -6.*u(i1,i2,i3,n)     \
                +4.*(u(i1+1,i2,i3,n)+u(i1-1,i2,i3,n))     \
                   -(u(i1+2,i2,i3,n)+u(i1-2,i2,i3,n)) )*adxSosup(0) + \
              ( -6.*u(i1,i2,i3,n)     \
                +4.*(u(i1,i2+1,i3,n)+u(i1,i2-1,i3,n))     \
                   -(u(i1,i2+2,i3,n)+u(i1,i2-2,i3,n)) )*adxSosup(1) + \
              ( -6.*u(i1,i2,i3,n)     \
                +4.*(u(i1,i2,i3+1,n)+u(i1,i2,i3-1,n))     \
                   -(u(i1,i2,i3+2,n)+u(i1,i2,i3-2,n)) )*adxSosup(2)


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (6th-order difference used with 4th-order scheme)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#defineMacro sosupDiss2d6(u,i1,i2,i3,n) \
             (-20.*u(i1,i2,i3,n)     \
               +15.*(u(i1+1,i2,i3,n)+u(i1-1,i2,i3,n))     \
                -6.*(u(i1+2,i2,i3,n)+u(i1-2,i2,i3,n))     \
                   +(u(i1+3,i2,i3,n)+u(i1-3,i2,i3,n))  )*adxSosup(0) + \
              (-20.*u(i1,i2,i3,n)     \
               +15.*(u(i1,i2+1,i3,n)+u(i1,i2-1,i3,n))     \
                -6.*(u(i1,i2+2,i3,n)+u(i1,i2-2,i3,n))     \
                   +(u(i1,i2+3,i3,n)+u(i1,i2-3,i3,n))  )*adxSosup(1)

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (3D 6th-order difference used with 4th-order scheme)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#defineMacro sosupDiss3d6(u,i1,i2,i3,n) \
             (-20.*u(i1,i2,i3,n)     \
               +15.*(u(i1+1,i2,i3,n)+u(i1-1,i2,i3,n))     \
                -6.*(u(i1+2,i2,i3,n)+u(i1-2,i2,i3,n))     \
                   +(u(i1+3,i2,i3,n)+u(i1-3,i2,i3,n))  )*adxSosup(0) + \
              (-20.*u(i1,i2,i3,n)     \
               +15.*(u(i1,i2+1,i3,n)+u(i1,i2-1,i3,n))     \
                -6.*(u(i1,i2+2,i3,n)+u(i1,i2-2,i3,n))     \
                   +(u(i1,i2+3,i3,n)+u(i1,i2-3,i3,n))  )*adxSosup(1) + \
              (-20.*u(i1,i2,i3,n)     \
               +15.*(u(i1,i2,i3+1,n)+u(i1,i2,i3-1,n))     \
                -6.*(u(i1,i2,i3+2,n)+u(i1,i2,i3-2,n))     \
                   +(u(i1,i2,i3+3,n)+u(i1,i2,i3-3,n))  )*adxSosup(2)

                   
! =========================================================================================
! **OLD***
! Macro: Update isotropic equations 2D, Order 2
! ========================================================================================
#beginMacro updateMx2dOrder2()
  ! eps* u_t = w_y
  ! eps* v_t = - w_x
  ! mu* w_t = u_y - v_x 

  write(*,'("advBA: advance ISOTROPIC 2D, 2nd-order, rectangular... t=",e10.2)') t

  !  Taylor2 time-stepping: 
  !    u(n+1) = u(n) + dt*ut + (dt^2/2)*utt 
  if( solveForAllFields.ne.0 )then

    if( .not.methodOfLines )then

      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  
        if( addForcing.ne.0 )then  ! do this for now *fix me*
          do m=0,5
            fv(m)=f(i1,i2,i3,m)
          end do
        end if 
  
        un(i1,i2,i3,ex)=u(i1,i2,i3,ex) + (dt/eps)*uy22r(i1,i2,i3,hz) +.5*cdtsq*lap2d2(i1,i2,i3,ex) + dt*fv(ex)
        un(i1,i2,i3,ey)=u(i1,i2,i3,ey) - (dt/eps)*ux22r(i1,i2,i3,hz) +.5*cdtsq*lap2d2(i1,i2,i3,ey) + dt*fv(ey)
        un(i1,i2,i3,ez)=u(i1,i2,i3,ez) + (dt/eps)*(ux22r(i1,i2,i3,hy)-uy22r(i1,i2,i3,hx)) \
                                                                     +.5*cdtsq*lap2d2(i1,i2,i3,ez) + dt*fv(ez)
  
        un(i1,i2,i3,hx)=u(i1,i2,i3,hx) + (dt/mu)*( -uy22r(i1,i2,i3,ez) ) +.5*cdtsq*lap2d2(i1,i2,i3,hx) + dt*fv(hx)
        un(i1,i2,i3,hy)=u(i1,i2,i3,hy) + (dt/mu)*(  ux22r(i1,i2,i3,ez) ) +.5*cdtsq*lap2d2(i1,i2,i3,hy) + dt*fv(hy)
  
        un(i1,i2,i3,hz)=u(i1,i2,i3,hz) + (dt/mu)*(uy22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ey)) \
                                                                   +.5*cdtsq*lap2d2(i1,i2,i3,hz) + dt*fv(hz)
      endLoopsMask()

    else
      ! METHOD OF LINES
      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  
        if( addForcing.ne.0 )then  ! do this for now *fix me*
          do m=0,5
            fv(m)=f(i1,i2,i3,m)
          end do
        end if 
  
        un(i1,i2,i3,ex)= + (1./eps)*uy22r(i1,i2,i3,hz) + fv(ex)
        un(i1,i2,i3,ey)= - (1./eps)*ux22r(i1,i2,i3,hz) + fv(ey)
        un(i1,i2,i3,ez)= + (1./eps)*(ux22r(i1,i2,i3,hy)-uy22r(i1,i2,i3,hx)) + fv(ez)
  
        un(i1,i2,i3,hx)= + (1./mu)*( -uy22r(i1,i2,i3,ez) ) + fv(hx)
        un(i1,i2,i3,hy)= + (1./mu)*(  ux22r(i1,i2,i3,ez) ) + fv(hy)
        un(i1,i2,i3,hz)= + (1./mu)*(uy22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ey))  + fv(hz)
  
      endLoopsMask()

    end if 

  else
     ! TEz polarization
    if( .not.methodOfLines )then
      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  
       if( addForcing.ne.0 )then  ! do this for now *fix me*
         fv(ex)=f(i1,i2,i3,ex)
         fv(ey)=f(i1,i2,i3,ey)
         fv(hz)=f(i1,i2,i3,hz)
       else
         do c=ex,hz
           fv(c)= 0.
         end do
       end if 
  
       un(i1,i2,i3,ex)=u(i1,i2,i3,ex) + (dt/eps)*uy22r(i1,i2,i3,hz) +.5*cdtsq*lap2d2(i1,i2,i3,ex) + dt*fv(ex)
       un(i1,i2,i3,ey)=u(i1,i2,i3,ey) - (dt/eps)*ux22r(i1,i2,i3,hz) +.5*cdtsq*lap2d2(i1,i2,i3,ey) + dt*fv(ey)
       un(i1,i2,i3,hz)=u(i1,i2,i3,hz) + (dt/mu)*(uy22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ey)) \
                                                                   +.5*cdtsq*lap2d2(i1,i2,i3,hz) + dt*fv(hz)
      endLoopsMask()
    else
      ! METHOD OF LINES
      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  
       if( addForcing.ne.0 )then  ! do this for now *fix me*
         fv(ex)=f(i1,i2,i3,ex)
         fv(ey)=f(i1,i2,i3,ey)
         fv(hz)=f(i1,i2,i3,hz)
       else
         do c=ex,hz
           fv(c)= 0.
         end do
       end if 
  
       un(i1,i2,i3,ex)= + (1./eps)*uy22r(i1,i2,i3,hz) + fv(ex)
       un(i1,i2,i3,ey)= - (1./eps)*ux22r(i1,i2,i3,hz) + fv(ey)
       un(i1,i2,i3,hz)= + (1./mu)*(uy22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ey)) + fv(hz)
      endLoopsMask()

    end if

  end if

#endMacro

! =========================================================================================
! **OLD***
! Macro: Update BI-ANISTROPIC equations 2D, Order 2
! ========================================================================================
#beginMacro updateBA2dOrder2()
  ! eps* u_t = w_y
  ! eps* v_t = - w_x
  ! mu* w_t = u_y - v_x 

  if( t.lt.2*dt )then
    write(*,'("advBA: advance BI-ANISTROPIC 2D, 2nd-order, rectangular... t=",e10.2)') t
  end if
  mr=0
  
  !  Taylor2 time-stepping: 
  !    u(n+1) = u(n) + dt*ut + (dt^2/2)*utt 
  if( solveForAllFields.ne.0 )then

    if( .not.methodOfLines )then

      ! Finish method of lines for general K0
      !     dt*K0i*curl + (dt^2/2) * K0i^2*( curl(curl ))
       
      stop 111
      eps = 1./K0i(0,0,0)
      mu = 1./K0i(5,5,0)
       
      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  
        if( addForcing.ne.0 )then  ! do this for now *fix me*
          do m=0,5
            fv(m)=f(i1,i2,i3,m)
          end do
        end if 
  
        un(i1,i2,i3,ex)=u(i1,i2,i3,ex) + (dt/eps)*uy22r(i1,i2,i3,hz) +.5*cdtsq*lap2d2(i1,i2,i3,ex) + dt*fv(ex)
        un(i1,i2,i3,ey)=u(i1,i2,i3,ey) - (dt/eps)*ux22r(i1,i2,i3,hz) +.5*cdtsq*lap2d2(i1,i2,i3,ey) + dt*fv(ey)
        un(i1,i2,i3,ez)=u(i1,i2,i3,ez) + (dt/eps)*(ux22r(i1,i2,i3,hy)-uy22r(i1,i2,i3,hx)) \
                                                                     +.5*cdtsq*lap2d2(i1,i2,i3,ez) + dt*fv(ez)
  
        un(i1,i2,i3,hx)=u(i1,i2,i3,hx) + (dt/mu)*( -uy22r(i1,i2,i3,ez) ) +.5*cdtsq*lap2d2(i1,i2,i3,hx) + dt*fv(hx)
        un(i1,i2,i3,hy)=u(i1,i2,i3,hy) + (dt/mu)*(  ux22r(i1,i2,i3,ez) ) +.5*cdtsq*lap2d2(i1,i2,i3,hy) + dt*fv(hy)
  
        un(i1,i2,i3,hz)=u(i1,i2,i3,hz) + (dt/mu)*(uy22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ey)) \
                                                                   +.5*cdtsq*lap2d2(i1,i2,i3,hz) + dt*fv(hz)
      endLoopsMask()

    else
      ! METHOD OF LINES
      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  
        if( addForcing.ne.0 )then  ! do this for now *fix me*
          do m=0,5
            fv(m)=f(i1,i2,i3,m)
          end do
          ! write(*,'(" advBA: i1,i2=",2i4," force(Hz)=",e10.2)') i1,i2,fv(5) 
        end if 
  
        ! compute components of the curl(H) and -curl(E)
        curl(0) =  uy22r(i1,i2,i3,hz)
        curl(1) = -ux22r(i1,i2,i3,hz)
        curl(2) = (ux22r(i1,i2,i3,hy)-uy22r(i1,i2,i3,hx))

        curl(3) = -uy22r(i1,i2,i3,ez)
        curl(4) =  ux22r(i1,i2,i3,ez)
        curl(5) =-(ux22r(i1,i2,i3,ey)-uy22r(i1,i2,i3,ex))

        if( numberOfMaterialRegions.gt.1 )then
          mr = matMask(i1,i2,i3)
          if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
             stop 9999
          end if            
        end if       
        do m=0,5 
          un(i1,i2,i3,m)= K0i(m,0,mr)*curl(0) + K0i(m,1,mr)*curl(1) + K0i(m,2,mr)*curl(2) + \
                          K0i(m,3,mr)*curl(3) + K0i(m,4,mr)*curl(4) + K0i(m,5,mr)*curl(5) + fv(m)
        end do
        !  write(*,'(" i1,i2=",2i3," curl=",6(1pe10.3,1x)," hz=",e14.6)') i1,i2,(curl(m),m=0,5),u(i1,i2,i3,hz)
        !  write(*,'("  un=",6(1pe10.3,1x))') (un(i1,i2,i3,m),m=0,5)
        
!        un(i1,i2,i3,ex)= + (1./eps)*uy22r(i1,i2,i3,hz) + fv(ex)
!        un(i1,i2,i3,ey)= - (1./eps)*ux22r(i1,i2,i3,hz) + fv(ey)
!        un(i1,i2,i3,ez)= + (1./eps)*(ux22r(i1,i2,i3,hy)-uy22r(i1,i2,i3,hx)) + fv(ez)
!  
!        un(i1,i2,i3,hx)= + (1./mu)*( -uy22r(i1,i2,i3,ez) ) + fv(hx)
!        un(i1,i2,i3,hy)= + (1./mu)*(  ux22r(i1,i2,i3,ez) ) + fv(hy)
!        un(i1,i2,i3,hz)= + (1./mu)*(uy22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ey))  + fv(hz)

      endLoopsMask()

    end if 

  else
     ! TEz polarization
    if( .not.methodOfLines )then
      stop 111

      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  
       if( addForcing.ne.0 )then  ! do this for now *fix me*
         fv(ex)=f(i1,i2,i3,ex)
         fv(ey)=f(i1,i2,i3,ey)
         fv(hz)=f(i1,i2,i3,hz)
       else
         do c=ex,hz
           fv(c)= 0.
         end do
       end if 
  
       un(i1,i2,i3,ex)=u(i1,i2,i3,ex) + (dt/eps)*uy22r(i1,i2,i3,hz) +.5*cdtsq*lap2d2(i1,i2,i3,ex) + dt*fv(ex)
       un(i1,i2,i3,ey)=u(i1,i2,i3,ey) - (dt/eps)*ux22r(i1,i2,i3,hz) +.5*cdtsq*lap2d2(i1,i2,i3,ey) + dt*fv(ey)
       un(i1,i2,i3,hz)=u(i1,i2,i3,hz) + (dt/mu)*(uy22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ey)) \
                                                                   +.5*cdtsq*lap2d2(i1,i2,i3,hz) + dt*fv(hz)
      endLoopsMask()
    else
      ! METHOD OF LINES -- TEz 

      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  
       if( addForcing.ne.0 )then  ! do this for now *fix me*
         fv(ex)=f(i1,i2,i3,ex)
         fv(ey)=f(i1,i2,i3,ey)
         fv(hz)=f(i1,i2,i3,hz)
       else
         do c=ex,hz
           fv(c)= 0.
         end do
       end if 
  
        ! compute components of the curl(H) and -curl(E)
        curl(0) =  uy22r(i1,i2,i3,hz)
        curl(1) = -ux22r(i1,i2,i3,hz)
        ! curl(2) = (ux22r(i1,i2,i3,hy)-uy22r(i1,i2,i3,hx))

        ! curl(3) = -uy22r(i1,i2,i3,ez)
        ! curl(4) =  ux22r(i1,i2,i3,ez)
        curl(2) =-(ux22r(i1,i2,i3,ey)-uy22r(i1,i2,i3,ex))

        if( numberOfMaterialRegions.gt.1 )then
          mr = matMask(i1,i2,i3)
          if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
             stop 9999
          end if            
        end if 
        do m=0,2
          un(i1,i2,i3,m)= Ki(m,0,mr)*curl(0) + Ki(m,1,mr)*curl(1) + Ki(m,2,mr)*curl(2) + fv(m)
        end do 

        ! un(i1,i2,i3,ex)= + (1./eps)*uy22r(i1,i2,i3,hz) + fv(ex)
        ! un(i1,i2,i3,ey)= - (1./eps)*ux22r(i1,i2,i3,hz) + fv(ey)
        ! un(i1,i2,i3,hz)= + (1./mu)*(uy22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ey)) + fv(hz)
      endLoopsMask()

    end if

  end if

#endMacro


! =========================================================================================
! **OLD***
! Macro: Update isotropic equations 2D, Order 4
! ========================================================================================
#beginMacro updateMx2dOrder4()
  if( solveForAllFields.ne.0 ) then

    if( .false. )then
      write(*,'("advBA: FOS 2D, Order 4, rect... t=",e10.2," addForcing=",i3)') t,addForcing 
    end if       


    if( .not.methodOfLines )then

     !  Taylor4 time-stepping: 
     !    u(n+1) = u(n) + dt*ut + (dt^2/2)*utt + (dt*^3/6)*uttt + (dt^4/24)*utttt
     ! FINISH ME 
     beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

       if( addForcing.ne.0 )then ! do this for now *fix me*
         do m=0,5
           fv(m)=f(i1,i2,i3,m)
         end do 
       end if

       un(i1,i2,i3,ex)=u(i1,i2,i3,ex) + (dt/eps)*uy42r(i1,i2,i3,hz) +.5*cdtsq*lap2d4(i1,i2,i3,ex) \
                  + (cdtsq*dt/(6.*eps))*( uxxy22r(i1,i2,i3,hz) + uyyy22r(i1,i2,i3,hz) ) +  \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,ex) \
                  + dt*fv(ex)

       un(i1,i2,i3,ey)=u(i1,i2,i3,ey) - (dt/eps)*ux42r(i1,i2,i3,hz) +.5*cdtsq*lap2d4(i1,i2,i3,ey) \
                  - (cdtsq*dt/(6.*eps))*( uxxx22r(i1,i2,i3,hz) + uxyy22r(i1,i2,i3,hz) ) + \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,ey) \
                  + dt*fv(ey)

       un(i1,i2,i3,ez)=u(i1,i2,i3,ez) - (dt/eps)*(uy42r(i1,i2,i3,hx)-ux42r(i1,i2,i3,hy)) +.5*cdtsq*lap2d4(i1,i2,i3,ez) \
             - (cdtsq*dt/(6.*eps))*( (uxxy22r(i1,i2,i3,hx) + uyyy22r(i1,i2,i3,hx)) \
                                    -(uxxx22r(i1,i2,i3,hy) + uxyy22r(i1,i2,i3,hy)) ) \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,ez) \
                             + dt*fv(ez)

       un(i1,i2,i3,hx)=u(i1,i2,i3,hx) - (dt/mu )*uy42r(i1,i2,i3,ez) +.5*cdtsq*lap2d4(i1,i2,i3,hx) \
                  - (cdtsq*dt/(6.*mu ))*( uxxy22r(i1,i2,i3,ez) + uyyy22r(i1,i2,i3,ez) ) +  \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,hx) \
                  + dt*fv(hx)

       un(i1,i2,i3,hy)=u(i1,i2,i3,hy) + (dt/mu )*ux42r(i1,i2,i3,ez) +.5*cdtsq*lap2d4(i1,i2,i3,hy) \
                  + (cdtsq*dt/(6.*mu ))*( uxxx22r(i1,i2,i3,ez) + uxyy22r(i1,i2,i3,ez) ) + \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,hy) \
                  + dt*fv(hy)

       un(i1,i2,i3,hz)=u(i1,i2,i3,hz) + (dt/mu)*(uy42r(i1,i2,i3,ex)-ux42r(i1,i2,i3,ey)) +.5*cdtsq*lap2d4(i1,i2,i3,hz) \
             + (cdtsq*dt/(6.*mu))*( (uxxy22r(i1,i2,i3,ex) + uyyy22r(i1,i2,i3,ex)) \
                                   -(uxxx22r(i1,i2,i3,ey) + uxyy22r(i1,i2,i3,ey)) ) \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,hz) \
                             + dt*fv(hz)
     endLoopsMask()
   else
     ! METHOD OF LINES
     beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

       if( addForcing.ne.0 )then  ! do this for now *fix me*
         do m=0,5
           fv(m)=f(i1,i2,i3,m)
         end do
       end if 

       un(i1,i2,i3,ex)= + (1./eps)*uy42r(i1,i2,i3,hz) + fv(ex)
       un(i1,i2,i3,ey)= - (1./eps)*ux42r(i1,i2,i3,hz) + fv(ey)
       un(i1,i2,i3,ez)= + (1./eps)*(ux42r(i1,i2,i3,hy)-uy42r(i1,i2,i3,hx)) + fv(ez)

       un(i1,i2,i3,hx)= + (1./mu)*( -uy42r(i1,i2,i3,ez) ) + fv(hx)
       un(i1,i2,i3,hy)= + (1./mu)*(  ux42r(i1,i2,i3,ez) ) + fv(hy)
       un(i1,i2,i3,hz)= + (1./mu)*(uy42r(i1,i2,i3,ex)-ux42r(i1,i2,i3,ey))  + fv(hz)

     endLoopsMask()

   end if

  else if( .true. )then 

    ! TEz polarization

    if( .false. )then
      write(*,'("advBA: FOS 2D, Order 4, rect... t=",e10.2," addForcing=",i3)') t,addForcing 
    end if       

   !  Taylor4 time-stepping: 
   !    u(n+1) = u(n) + dt*ut + (dt^2/2)*utt + (dt*^3/6)*uttt + (dt^4/24)*utttt

   if( .not.methodOfLines )then

     beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

       if( addForcing.ne.0 )then ! do this for now *fix me*
         fv(ex)=f(i1,i2,i3,ex)
         fv(ey)=f(i1,i2,i3,ey)
         fv(hz)=f(i1,i2,i3,hz)
       end if

       un(i1,i2,i3,ex)=u(i1,i2,i3,ex) + (dt/eps)*uy42r(i1,i2,i3,hz) +.5*cdtsq*lap2d4(i1,i2,i3,ex) \
                  + (cdtsq*dt/(6.*eps))*( uxxy22r(i1,i2,i3,hz) + uyyy22r(i1,i2,i3,hz) ) +  \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,ex) \
                  + dt*fv(ex)

       un(i1,i2,i3,ey)=u(i1,i2,i3,ey) - (dt/eps)*ux42r(i1,i2,i3,hz) +.5*cdtsq*lap2d4(i1,i2,i3,ey) \
                  - (cdtsq*dt/(6.*eps))*( uxxx22r(i1,i2,i3,hz) + uxyy22r(i1,i2,i3,hz) ) + \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,ey) \
                  + dt*fv(ey)

       un(i1,i2,i3,hz)=u(i1,i2,i3,hz) + (dt/mu)*(uy42r(i1,i2,i3,ex)-ux42r(i1,i2,i3,ey)) +.5*cdtsq*lap2d4(i1,i2,i3,hz) \
             + (cdtsq*dt/(6.*mu))*( (uxxy22r(i1,i2,i3,ex) + uyyy22r(i1,i2,i3,ex)) \
                                   -(uxxx22r(i1,i2,i3,ey) + uxyy22r(i1,i2,i3,ey)) ) \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,hz) \
                             + dt*fv(hz)
     endLoopsMask()
   else
     ! METHOD OF LINES
     beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

      if( addForcing.ne.0 )then  ! do this for now *fix me*
        fv(ex)=f(i1,i2,i3,ex)
        fv(ey)=f(i1,i2,i3,ey)
        fv(hz)=f(i1,i2,i3,hz)
      else
        do c=ex,hz
          fv(c)= 0.
        end do
      end if 

      un(i1,i2,i3,ex)= + (1./eps)*uy42r(i1,i2,i3,hz) + fv(ex)
      un(i1,i2,i3,ey)= - (1./eps)*ux42r(i1,i2,i3,hz) + fv(ey)
      un(i1,i2,i3,hz)= + (1./mu)*(uy42r(i1,i2,i3,ex)-ux42r(i1,i2,i3,ey)) + fv(hz)
     endLoopsMask()
   end if
  end if
#endMacro

        
! =========================================================================================
! **OLD***
! Macro: Update BI-ANISTROPIC equations 2D, Order 4
! ========================================================================================
#beginMacro updateBA2dOrder4()

  mr=0 
  if( solveForAllFields.ne.0 ) then

    if( .false. )then
      write(*,'("advBA: FOS 2D, Order 4, rect... t=",e10.2," addForcing=",i3)') t,addForcing 
    end if       


    if( .not.methodOfLines )then

      stop 111
      eps = 1./K0i(0,0,0)
      mu = 1./K0i(5,5,0)

     !  Taylor4 time-stepping: 
     !    u(n+1) = u(n) + dt*ut + (dt^2/2)*utt + (dt*^3/6)*uttt + (dt^4/24)*utttt
     ! FINISH ME 
     beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

       if( addForcing.ne.0 )then ! do this for now *fix me*
         do m=0,5
           fv(m)=f(i1,i2,i3,m)
         end do 
       end if

       un(i1,i2,i3,ex)=u(i1,i2,i3,ex) + (dt/eps)*uy42r(i1,i2,i3,hz) +.5*cdtsq*lap2d4(i1,i2,i3,ex) \
                  + (cdtsq*dt/(6.*eps))*( uxxy22r(i1,i2,i3,hz) + uyyy22r(i1,i2,i3,hz) ) +  \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,ex) \
                  + dt*fv(ex)

       un(i1,i2,i3,ey)=u(i1,i2,i3,ey) - (dt/eps)*ux42r(i1,i2,i3,hz) +.5*cdtsq*lap2d4(i1,i2,i3,ey) \
                  - (cdtsq*dt/(6.*eps))*( uxxx22r(i1,i2,i3,hz) + uxyy22r(i1,i2,i3,hz) ) + \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,ey) \
                  + dt*fv(ey)

       un(i1,i2,i3,ez)=u(i1,i2,i3,ez) - (dt/eps)*(uy42r(i1,i2,i3,hx)-ux42r(i1,i2,i3,hy)) +.5*cdtsq*lap2d4(i1,i2,i3,ez) \
             - (cdtsq*dt/(6.*eps))*( (uxxy22r(i1,i2,i3,hx) + uyyy22r(i1,i2,i3,hx)) \
                                    -(uxxx22r(i1,i2,i3,hy) + uxyy22r(i1,i2,i3,hy)) ) \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,ez) \
                             + dt*fv(ez)

       un(i1,i2,i3,hx)=u(i1,i2,i3,hx) - (dt/mu )*uy42r(i1,i2,i3,ez) +.5*cdtsq*lap2d4(i1,i2,i3,hx) \
                  - (cdtsq*dt/(6.*mu ))*( uxxy22r(i1,i2,i3,ez) + uyyy22r(i1,i2,i3,ez) ) +  \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,hx) \
                  + dt*fv(hx)

       un(i1,i2,i3,hy)=u(i1,i2,i3,hy) + (dt/mu )*ux42r(i1,i2,i3,ez) +.5*cdtsq*lap2d4(i1,i2,i3,hy) \
                  + (cdtsq*dt/(6.*mu ))*( uxxx22r(i1,i2,i3,ez) + uxyy22r(i1,i2,i3,ez) ) + \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,hy) \
                  + dt*fv(hy)

       un(i1,i2,i3,hz)=u(i1,i2,i3,hz) + (dt/mu)*(uy42r(i1,i2,i3,ex)-ux42r(i1,i2,i3,ey)) +.5*cdtsq*lap2d4(i1,i2,i3,hz) \
             + (cdtsq*dt/(6.*mu))*( (uxxy22r(i1,i2,i3,ex) + uyyy22r(i1,i2,i3,ex)) \
                                   -(uxxx22r(i1,i2,i3,ey) + uxyy22r(i1,i2,i3,ey)) ) \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,hz) \
                             + dt*fv(hz)
     endLoopsMask()
   else
     ! METHOD OF LINES
     beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

       if( addForcing.ne.0 )then  ! do this for now *fix me*
         do m=0,5
           fv(m)=f(i1,i2,i3,m)
         end do
       end if 

        ! compute components of the curl(H) and -curl(E)
        curl(0) =  uy42r(i1,i2,i3,hz)
        curl(1) = -ux42r(i1,i2,i3,hz)
        curl(2) = (ux42r(i1,i2,i3,hy)-uy42r(i1,i2,i3,hx))

        curl(3) = -uy42r(i1,i2,i3,ez)
        curl(4) =  ux42r(i1,i2,i3,ez)
        curl(5) =-(ux42r(i1,i2,i3,ey)-uy42r(i1,i2,i3,ex))

        if( numberOfMaterialRegions.gt.1 )then
          mr = matMask(i1,i2,i3)
          if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
             stop 9999
          end if            
        end if            
        do m=0,5 
          un(i1,i2,i3,m)= K0i(m,0,mr)*curl(0) + K0i(m,1,mr)*curl(1) + K0i(m,2,mr)*curl(2) + \
                          K0i(m,3,mr)*curl(3) + K0i(m,4,mr)*curl(4) + K0i(m,5,mr)*curl(5) + fv(m)
        end do

       ! un(i1,i2,i3,ex)= + (1./eps)*uy42r(i1,i2,i3,hz) + fv(ex)
       ! un(i1,i2,i3,ey)= - (1./eps)*ux42r(i1,i2,i3,hz) + fv(ey)
       ! un(i1,i2,i3,ez)= + (1./eps)*(ux42r(i1,i2,i3,hy)-uy42r(i1,i2,i3,hx)) + fv(ez)

       ! un(i1,i2,i3,hx)= + (1./mu)*( -uy42r(i1,i2,i3,ez) ) + fv(hx)
       ! un(i1,i2,i3,hy)= + (1./mu)*(  ux42r(i1,i2,i3,ez) ) + fv(hy)
       ! un(i1,i2,i3,hz)= + (1./mu)*(uy42r(i1,i2,i3,ex)-ux42r(i1,i2,i3,ey))  + fv(hz)

     endLoopsMask()

   end if

  else if( .true. )then 

    ! TEz polarization

    if( .false. )then
      write(*,'("advBA: FOS 2D, Order 4, rect... t=",e10.2," addForcing=",i3)') t,addForcing 
    end if       

   !  Taylor4 time-stepping: 
   !    u(n+1) = u(n) + dt*ut + (dt^2/2)*utt + (dt*^3/6)*uttt + (dt^4/24)*utttt

   if( .not.methodOfLines )then

     stop 111

     beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

       if( addForcing.ne.0 )then ! do this for now *fix me*
         fv(ex)=f(i1,i2,i3,ex)
         fv(ey)=f(i1,i2,i3,ey)
         fv(hz)=f(i1,i2,i3,hz)
       end if

       un(i1,i2,i3,ex)=u(i1,i2,i3,ex) + (dt/eps)*uy42r(i1,i2,i3,hz) +.5*cdtsq*lap2d4(i1,i2,i3,ex) \
                  + (cdtsq*dt/(6.*eps))*( uxxy22r(i1,i2,i3,hz) + uyyy22r(i1,i2,i3,hz) ) +  \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,ex) \
                  + dt*fv(ex)

       un(i1,i2,i3,ey)=u(i1,i2,i3,ey) - (dt/eps)*ux42r(i1,i2,i3,hz) +.5*cdtsq*lap2d4(i1,i2,i3,ey) \
                  - (cdtsq*dt/(6.*eps))*( uxxx22r(i1,i2,i3,hz) + uxyy22r(i1,i2,i3,hz) ) + \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,ey) \
                  + dt*fv(ey)

       un(i1,i2,i3,hz)=u(i1,i2,i3,hz) + (dt/mu)*(uy42r(i1,i2,i3,ex)-ux42r(i1,i2,i3,ey)) +.5*cdtsq*lap2d4(i1,i2,i3,hz) \
             + (cdtsq*dt/(6.*mu))*( (uxxy22r(i1,i2,i3,ex) + uyyy22r(i1,i2,i3,ex)) \
                                   -(uxxx22r(i1,i2,i3,ey) + uxyy22r(i1,i2,i3,ey)) ) \
                  + (cdtsq*cdtsq/24.)*lap2d2Pow2(i1,i2,i3,hz) \
                             + dt*fv(hz)
     endLoopsMask()
   else
     ! METHOD OF LINES
     beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

       if( addForcing.ne.0 )then  ! do this for now *fix me*
         fv(ex)=f(i1,i2,i3,ex)
         fv(ey)=f(i1,i2,i3,ey)
         fv(hz)=f(i1,i2,i3,hz)
       else
         do c=ex,hz
           fv(c)= 0.
         end do
       end if 

       ! compute components of the curl(H) and -curl(E)
       curl(0) =  uy42r(i1,i2,i3,hz)
       curl(1) = -ux42r(i1,i2,i3,hz)
       ! curl(2) = (ux42r(i1,i2,i3,hy)-uy42r(i1,i2,i3,hx))

       ! curl(3) = -uy42r(i1,i2,i3,ez)
       ! curl(4) =  ux42r(i1,i2,i3,ez)
       curl(2) =-(ux42r(i1,i2,i3,ey)-uy42r(i1,i2,i3,ex))

       
       if( numberOfMaterialRegions.gt.1 )then
         mr = matMask(i1,i2,i3)
          if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
             stop 9999
          end if            
       end if            
       do m=0,2
         un(i1,i2,i3,m)= Ki(m,0,mr)*curl(0) + Ki(m,1,mr)*curl(1) + Ki(m,2,mr)*curl(2) + fv(m)
       end do 

       ! un(i1,i2,i3,ex)= + (1./eps)*uy42r(i1,i2,i3,hz) + fv(ex)
       ! un(i1,i2,i3,ey)= - (1./eps)*ux42r(i1,i2,i3,hz) + fv(ey)
       ! un(i1,i2,i3,hz)= + (1./mu)*(uy42r(i1,i2,i3,ex)-ux42r(i1,i2,i3,ey)) + fv(hz)
     endLoopsMask()
   end if
  end if
#endMacro

        
! =========================================================================================
! **OLD***
! Macro: Update BI-ANISTROPIC equations 3D, Order 2
! ========================================================================================
#beginMacro updateBA3dOrder2()

  if( t.lt.2*dt )then
    write(*,'("advBA: advance BI-ANISTROPIC 3D, 2nd-order, rectangular... t=",e10.2)') t
  end if
  mr=0
  
  if( methodOfLines )then

    ! METHOD OF LINES -- TEz 

    beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

      if( addForcing.ne.0 )then  ! do this for now *fix me*
        do m=ex,hz
          fv(m)=f(i1,i2,i3,m)
        end do
      end if 

      ! compute components of the curl(H) and -curl(E)
      curl(0) =  uy23r(i1,i2,i3,hz)-uz23r(i1,i2,i3,hy)
      curl(1) =  uz23r(i1,i2,i3,hx)-ux23r(i1,i2,i3,hz)
      curl(2) =  ux23r(i1,i2,i3,hy)-uy23r(i1,i2,i3,hx)

      curl(3) =-(uy23r(i1,i2,i3,ez)-uz23r(i1,i2,i3,ey))
      curl(4) =-(uz23r(i1,i2,i3,ex)-ux23r(i1,i2,i3,ez))
      curl(5) =-(ux23r(i1,i2,i3,ey)-uy23r(i1,i2,i3,ex))

      if( numberOfMaterialRegions.gt.1 )then
        mr = matMask(i1,i2,i3)
        if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
           stop 9999
        end if            
      end if       
      do m=ex,hz
        un(i1,i2,i3,m)= K0i(m,0,mr)*curl(0) + K0i(m,1,mr)*curl(1) + K0i(m,2,mr)*curl(2) + \
                        K0i(m,3,mr)*curl(3) + K0i(m,4,mr)*curl(4) + K0i(m,5,mr)*curl(5) + fv(m)
      end do

    endLoopsMask()

  else
    write(*,*) 'Finish me: BAMX 3D modified equation Order 2'
    stop 9367
  end if

#endMacro

! =========================================================================================
! **OLD***
! Macro: Update BI-ANISTROPIC equations 3D, Order 4
! ========================================================================================
#beginMacro updateBA3dOrder4()

  if( t.lt.2*dt )then
    write(*,'("advBA: advance BI-ANISTROPIC 3D, 4th-order, rectangular... t=",e10.2)') t
  end if
  mr=0
  
  if( methodOfLines )then

    ! METHOD OF LINES 

    beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

      if( addForcing.ne.0 )then  ! do this for now *fix me*
        do m=ex,hz
          fv(m)=f(i1,i2,i3,m)
        end do
      else
        do m=ex,hz
          fv(m)=0.
        end do
      end if 

      ! compute components of the curl(H) and -curl(E)
      curl(0) =  uy43r(i1,i2,i3,hz)-uz43r(i1,i2,i3,hy)
      curl(1) =  uz43r(i1,i2,i3,hx)-ux43r(i1,i2,i3,hz)
      curl(2) =  ux43r(i1,i2,i3,hy)-uy43r(i1,i2,i3,hx)

      curl(3) =-(uy43r(i1,i2,i3,ez)-uz43r(i1,i2,i3,ey))
      curl(4) =-(uz43r(i1,i2,i3,ex)-ux43r(i1,i2,i3,ez))
      curl(5) =-(ux43r(i1,i2,i3,ey)-uy43r(i1,i2,i3,ex))

      if( numberOfMaterialRegions.gt.1 )then
        mr = matMask(i1,i2,i3)
        if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
           stop 9999
        end if            
      end if       
      do m=0,5
        un(i1,i2,i3,m)= K0i(m,0,mr)*curl(0) + K0i(m,1,mr)*curl(1) + K0i(m,2,mr)*curl(2) + \
                        K0i(m,3,mr)*curl(3) + K0i(m,4,mr)*curl(4) + K0i(m,5,mr)*curl(5) + fv(m)
      end do

    endLoopsMask()

  else
    write(*,*) 'Finish me: BAMX 3D modified equation Order 4'
    stop 9367
  end if

#endMacro



! =========================================================================================
! **OLD***
! Macro: Update BI-ANISTROPIC GDM equations
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================
#beginMacro updateBAGDM(DIM,ORDER,GRIDTYPE)

  if( t.lt.2*dt )then
    write(*,'("advBA: advance BA GDM, dim=DIM, order=ORDER grid=GRIDTYPE... t=",e10.2)') t
  end if
  mr=0
  
  if( solveForAllFields.ne.0 )then

    ! --- SOLVE FOR ALL FIELDS ---

    if( .not.methodOfLines )then

      ! --- TAYLOR TIME-STEPPING --- 
      stop 111
       
    else

      ! --- METHOD OF LINES (RK) ---
      ! zero out some forcing terms 
      do m=0,numPolarizationTerms-1
         pv(m)=0.
      end do
 
      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  
        if( addForcing.ne.0 )then  ! do this for now *fix me*
          do m=0,5
            fv(m)=f(i1,i2,i3,m)
          end do
        end if 
  
        if( numberOfMaterialRegions.gt.1 )then
          mr = matMask(i1,i2,i3)
          if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
             stop 9999
          end if            
        end if       

        if( forcingOption.eq.twilightZoneForcing )then

          if( nd.eq.2 )then
            do m=0,5
              ec = m 
              OGDERIV2D( 0,0,0,0,i1,i2,i3,t, ec, ev(m)  )
              OGDERIV2D( 1,0,0,0,i1,i2,i3,t, ec, evt(m) )
            end do 
            ! eval the polarization terms and time derivatives 
            do m=0,numPolarizationTerms-1
              pc = m+6  ! TZ index 
              OGDERIV2D( 0,0,0,0,i1,i2,i3,t, pc, pv(m)  )
              OGDERIV2D( 1,0,0,0,i1,i2,i3,t, pc, pvt(m) )
            end do    
          else
            do m=0,5
              ec = m 
              OGDERIV3D( 0,0,0,0,i1,i2,i3,t, ec, ev(m)  )
              OGDERIV3D( 1,0,0,0,i1,i2,i3,t, ec, evt(m) )
            end do 
            ! eval the polarization terms and time derivatives 
            do m=0,numPolarizationTerms-1
              pc = m+6  ! TZ index 
              OGDERIV3D( 0,0,0,0,i1,i2,i3,t, pc, pv(m)  )
              OGDERIV3D( 1,0,0,0,i1,i2,i3,t, pc, pvt(m) )
            end do    

         end if

         !write(*,'(" i1,i2=",2i3," xy=",2(f5.2,1x))') i1,i2,xy(i1,i2,i3,0),xy(i1,i2,i3,1)
         !write(*,'(" i1,i2=",2i3," ev=",6(f6.3,1x))') i1,i2,(ev(m),m=0,5)
         !write(*,'(" i1,i2=",2i3," pv=",10(f6.3,1x))') i1,i2,(pv(m),m=0,numPolarizationTerms-1)

          pc=0  ! P
          qc=1  ! Q = P.t 
          do k1=1,6
            do k2=1,6
              ec=k2-1 ! E or H component 
              do n=1,Np(k1,k2,mr)
                a0 = gdmPar(1,n,k1,k2,mr)
                a1 = gdmPar(2,n,k1,k2,mr)
                b0 = gdmPar(3,n,k1,k2,mr)
                b1 = gdmPar(4,n,k1,k2,mr)

                ! TZ forcing for polarization equations: 
                fp(pc) = pvt(pc) - pv(qc)
                fp(qc) = pvt(qc) - (  a0*ev(ec) + a1*evt(ec) - b0*pv(pc)- b1*pv(qc) )

                pc=pc+2
                qc=qc+2
              end do
            end do
          end do 
       else
          ! no TZ : set polarization forcing to zero 
          do m=0,numPolarizationTerms-1
            fp(m)=0.
          end do 
       end if


        ! compute components of the curl(H) and -curl(E)
       #If #DIM eq "2"
         ! --- 2D -----
        #If #ORDER eq "2" 
         curl(0) =  uy22r(i1,i2,i3,hz)
         curl(1) = -ux22r(i1,i2,i3,hz)
         curl(2) = (ux22r(i1,i2,i3,hy)-uy22r(i1,i2,i3,hx))
 
         curl(3) = -uy22r(i1,i2,i3,ez)
         curl(4) =  ux22r(i1,i2,i3,ez)
         curl(5) =-(ux22r(i1,i2,i3,ey)-uy22r(i1,i2,i3,ex))
        #Else
         curl(0) =  uy42r(i1,i2,i3,hz)
         curl(1) = -ux42r(i1,i2,i3,hz)
         curl(2) = (ux42r(i1,i2,i3,hy)-uy42r(i1,i2,i3,hx))
 
         curl(3) = -uy42r(i1,i2,i3,ez)
         curl(4) =  ux42r(i1,i2,i3,ez)
         curl(5) =-(ux42r(i1,i2,i3,ey)-uy42r(i1,i2,i3,ex))
        #End
       #Else
         ! --- 3D -----
        #If #ORDER eq "2" 
          curl(0) =  uy23r(i1,i2,i3,hz)-uz23r(i1,i2,i3,hy)
          curl(1) =  uz23r(i1,i2,i3,hx)-ux23r(i1,i2,i3,hz)
          curl(2) =  ux23r(i1,i2,i3,hy)-uy23r(i1,i2,i3,hx)
    
          curl(3) =-(uy23r(i1,i2,i3,ez)-uz23r(i1,i2,i3,ey))
          curl(4) =-(uz23r(i1,i2,i3,ex)-ux23r(i1,i2,i3,ez))
          curl(5) =-(ux23r(i1,i2,i3,ey)-uy23r(i1,i2,i3,ex))
        #Else
          curl(0) =  uy43r(i1,i2,i3,hz)-uz43r(i1,i2,i3,hy)
          curl(1) =  uz43r(i1,i2,i3,hx)-ux43r(i1,i2,i3,hz)
          curl(2) =  ux43r(i1,i2,i3,hy)-uy43r(i1,i2,i3,hx)
    
          curl(3) =-(uy43r(i1,i2,i3,ez)-uz43r(i1,i2,i3,ey))
          curl(4) =-(uz43r(i1,i2,i3,ex)-ux43r(i1,i2,i3,ez))
          curl(5) =-(ux43r(i1,i2,i3,ey)-uy43r(i1,i2,i3,ex))
        #End
 
       #End  

        if( debug.gt.3 )then
          write(*,'("----- i1,i2=",2i3)') i1,i2
        end if

        ! ---- Compute q = p.t = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, qc )
        pc=0  ! P
        qc=1  ! Q = P.t 
        do m=0,5
          ptSum(m)=0.   ! local total Pt or Mt 
        end do 
        do k1=1,6
          ec=k1-1 ! E or H component 
          do k2=1,6
            do n=1,Np(k1,k2,mr)
              ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)

        ! if( debug.gt.3 )then
        !   write(*,'("  q,qv,ptSum=",3(e10.2,1x))') p(i1,i2,i3,qc),pv(qc),ptSum(ec)
        ! end if

              qc=qc+2
            end do
          end do
          ! subtract off P.t = sum Q_m
          curl(ec) = curl(ec) - ptSum(ec)
        end do 

        ! if( debug.gt.3 )then
        !   write(*,'("   ptSum=",6(f6.3,1x))') (ptSum(m),m=0,5)
        ! end if

        do m=0,5 
          un(i1,i2,i3,m)= K0i(m,0,mr)*curl(0) + K0i(m,1,mr)*curl(1) + K0i(m,2,mr)*curl(2) + \
                          K0i(m,3,mr)*curl(3) + K0i(m,4,mr)*curl(4) + K0i(m,5,mr)*curl(5) + fv(m)
        end do

        ! if( .false. )then
        !   write(*,'(" i1,i2=",2i3," ut=",6(f6.3,1x))') i1,i2,(un(i1,i2,i3,m),m=0,5)
        ! end if

        ! --- compute time derivatives of P and Q
        ! p.t = q
        ! q.t = a0*E + a1*Et - b0*p - b1*q   
        ! or 
        ! q.t = a0*H + a1*Ht - b0*p - b1*q   
        pc=0  ! P
        qc=1  ! Q = P.t 
        do k1=1,6
          do k2=1,6
            ec=k2-1 ! This GDM term involves this E or H component 

            do n=1,Np(k1,k2,mr)
              a0 = gdmPar(1,n,k1,k2,mr)
              a1 = gdmPar(2,n,k1,k2,mr)
              b0 = gdmPar(3,n,k1,k2,mr)
              b1 = gdmPar(4,n,k1,k2,mr)

              pct=pc+6  ! p.t is stored in un here 
              qct=qc+6  ! q.t is stored in un here

              un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
              un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
              ! test: 
              ! un(i1,i2,i3,qct) = a0*ev(ec)         + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 

              ! if( debug.gt.3 )then
              !   write(*,'(" i1,i2=",2i4," k1,k2=",2i2," Np=",i3," ec=",i1," t=",e9.2)') i1,i2,k1,k2,Np(k1,k2,mr),ec,t
              !   write(*,'("   ... pc,qc,ec=",3i3," a0,a1,b0,b1=",4e10.2)') pc,qc,ec,a0,a1,b0,b1
              !   write(*,'("   ... p ,q =",2f7.3)') p(i1,i2,i3,pc),p(i1,i2,i3,qc)
              ! end if 
              ! write(*,'("   ... p ,pe =",2f7.3," q ,qe =",2f7.3)') p(i1,i2,i3,pc),pv(pc),p(i1,i2,i3,qc),pv(qc)
              ! write(*,'("   ... E ,Ee =",2f7.3," Et,Ete =",2f7.3)') u(i1,i2,i3,ec),ev(ec),un(i1,i2,i3,ec),evt(ec)
              ! write(*,'("   ... pt,pte=",2f7.3," qt,qte=",2f7.3)') un(i1,i2,i3,pct),pvt(pc),un(i1,i2,i3,qct),pvt(qc)

              ! test: -- set to exact time-derivative 
              ! un(i1,i2,i3,pct) = pvt(pc)
              ! un(i1,i2,i3,qct) = pvt(qc)
             
              pc=pc+2
              qc=qc+2
            end do
          end do
        end do 

      endLoopsMask()

      if( .false. .or. debug > 15 )then
        stop  4444
     end if
     
    end if 

  else

    ! TEz polarization
     
    stop 5555

    if( .not.methodOfLines )then
      stop 111

      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  
       if( addForcing.ne.0 )then  ! do this for now *fix me*
         fv(ex)=f(i1,i2,i3,ex)
         fv(ey)=f(i1,i2,i3,ey)
         fv(hz)=f(i1,i2,i3,hz)
       else
         do c=ex,hz
           fv(c)= 0.
         end do
       end if 
  
       un(i1,i2,i3,ex)=u(i1,i2,i3,ex) + (dt/eps)*uy22r(i1,i2,i3,hz) +.5*cdtsq*lap2d2(i1,i2,i3,ex) + dt*fv(ex)
       un(i1,i2,i3,ey)=u(i1,i2,i3,ey) - (dt/eps)*ux22r(i1,i2,i3,hz) +.5*cdtsq*lap2d2(i1,i2,i3,ey) + dt*fv(ey)
       un(i1,i2,i3,hz)=u(i1,i2,i3,hz) + (dt/mu)*(uy22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ey)) \
                                                                   +.5*cdtsq*lap2d2(i1,i2,i3,hz) + dt*fv(hz)
      endLoopsMask()
    else
      ! METHOD OF LINES -- TEz 

      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
  
       if( addForcing.ne.0 )then  ! do this for now *fix me*
         fv(ex)=f(i1,i2,i3,ex)
         fv(ey)=f(i1,i2,i3,ey)
         fv(hz)=f(i1,i2,i3,hz)
       else
         do c=ex,hz
           fv(c)= 0.
         end do
       end if 
  
        ! compute components of the curl(H) and -curl(E)
        curl(0) =  uy22r(i1,i2,i3,hz)
        curl(1) = -ux22r(i1,i2,i3,hz)
        ! curl(2) = (ux22r(i1,i2,i3,hy)-uy22r(i1,i2,i3,hx))

        ! curl(3) = -uy22r(i1,i2,i3,ez)
        ! curl(4) =  ux22r(i1,i2,i3,ez)
        curl(2) =-(ux22r(i1,i2,i3,ey)-uy22r(i1,i2,i3,ex))

        if( numberOfMaterialRegions.gt.1 )then
          mr = matMask(i1,i2,i3)
          if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
             stop 9999
          end if            
        end if 
        do m=0,2
          un(i1,i2,i3,m)= Ki(m,0,mr)*curl(0) + Ki(m,1,mr)*curl(1) + Ki(m,2,mr)*curl(2) + fv(m)
        end do 

        ! un(i1,i2,i3,ex)= + (1./eps)*uy22r(i1,i2,i3,hz) + fv(ex)
        ! un(i1,i2,i3,ey)= - (1./eps)*ux22r(i1,i2,i3,hz) + fv(ey)
        ! un(i1,i2,i3,hz)= + (1./mu)*(uy22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ey)) + fv(hz)
      endLoopsMask()

    end if

  end if

#endMacro


! =========================================================================================
! ** OPTIMIZED VERSION**
! Macro: Update BI-ANISTROPIC GDM equations
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
!   ETAXOP: for supergrid layer in x-direction set to etxa(i1)*, otherwise leave null
! ========================================================================================
#beginMacro updateBAGDMOpt(DIM,ORDER,GRIDTYPE,POLAR,ETAXOP,ETAYOP,ETAZOP)

  if( t.lt.2*dt )then
    write(*,'("advBA: advance BA GDM dim=DIM order=ORDER grid=GRIDTYPE polar=POLAR... t=",e10.2)') t
  end if
  
  ! ---- Precompute some indirection arrays to make the dispersion loops go faster ----
  do mr=0,numberOfMaterialRegions-1
    m=0
    pc=0
    qc=1
    do k1=1,6
      ec=k1-1 ! E or H component 
      do k2=1,6
        do n=1,Np(k1,k2,mr)
          m=m+1
          ecIndex1(m,mr)=ec
          qcIndex1(m,mr)=qc
          ! ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)
          qc=qc+2
        end do
      end do
      ! subtract off P.t = sum Q_m
      ! curl(ec) = curl(ec) - ptSum(ec)
    end do 
    numTerms1(mr)=m
    if( numTerms1(mr).gt.maxNumPolarizationTerms )then
      write(*,'(" ERROR numTerms1=",i6," too big")') numTerms1(mr) 
      stop 1616
    end if 
  
    pc=0  ! P
    qc=1  ! Q = P.t 
    m=0 
    do k1=1,6
      do k2=1,6
        ec=k2-1 ! This GDM term involves this E or H component 
  
        do n=1,Np(k1,k2,mr)
          a0 = gdmPar(1,n,k1,k2,mr)
          a1 = gdmPar(2,n,k1,k2,mr)
          b0 = gdmPar(3,n,k1,k2,mr)
          b1 = gdmPar(4,n,k1,k2,mr)
  
          pct=pc+numComp  ! p.t is stored in un here 
          qct=qc+numComp  ! q.t is stored in un here
  
          m=m+1
          ecIndex2(m,mr)=ec
          pcIndex2(m,mr)=pc
          a0v(m,mr)=a0
          a1v(m,mr)=a1
          b0v(m,mr)=b0
          b1v(m,mr)=b1
          
          ! un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
          ! un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
  
          pc=pc+2
          qc=qc+2
        end do
      end do
    end do 
    numTerms2(mr)=m
    if( numTerms2(mr).gt.maxNumPolarizationTerms )then
      write(*,'(" ERROR numTerms2=",i6," too big")') numTerms2(mr)
      stop 1616
    end if 
  end do ! end do mr
  mr=0

  if( .not.methodOfLines )then

    ! --- TAYLOR TIME-STEPPING --- 
    stop 111
     
  else

    ! --- METHOD OF LINES (RK) ---
    ! zero out some forcing terms 
    do m=0,numPolarizationTerms-1
       pv(m)=0.
      fp(m)=0.
    end do

 
    beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

      if( addForcing.ne.0 )then  ! do this for now *fix me*
        #If #POLAR eq "TEZ" 
          do m=0,2
           fv(m)=f(i1,i2,i3,m)
          end do
        #Else
          do m=0,5
           fv(m)=f(i1,i2,i3,m)
          end do
        #End 
      end if 

      if( numberOfMaterialRegions.gt.1 )then
        mr = matMask(i1,i2,i3)
        if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
           stop 9999
        end if            
      end if       

      if( forcingOption.eq.twilightZoneForcing )then

        if( nd.eq.2 )then
          do m=0,numComp-1
            ec = m 
            OGDERIV2D( 0,0,0,0,i1,i2,i3,t, ec, ev(m)  )
            OGDERIV2D( 1,0,0,0,i1,i2,i3,t, ec, evt(m) )
          end do 
          ! eval the polarization terms and time derivatives 
          do m=0,numPolarizationTerms-1
            pc = m+numComp  ! TZ index 
            OGDERIV2D( 0,0,0,0,i1,i2,i3,t, pc, pv(m)  )
            OGDERIV2D( 1,0,0,0,i1,i2,i3,t, pc, pvt(m) )
          end do    
        else
          do m=0,numComp-1
            ec = m 
            OGDERIV3D( 0,0,0,0,i1,i2,i3,t, ec, ev(m)  )
            OGDERIV3D( 1,0,0,0,i1,i2,i3,t, ec, evt(m) )
          end do 
          ! eval the polarization terms and time derivatives 
          do m=0,numPolarizationTerms-1
            pc = m+numComp  ! TZ index 
            OGDERIV3D( 0,0,0,0,i1,i2,i3,t, pc, pv(m)  )
            OGDERIV3D( 1,0,0,0,i1,i2,i3,t, pc, pvt(m) )
          end do    

       end if

       !write(*,'(" i1,i2=",2i3," xy=",2(f5.2,1x))') i1,i2,xy(i1,i2,i3,0),xy(i1,i2,i3,1)
       !write(*,'(" i1,i2=",2i3," ev=",6(f6.3,1x))') i1,i2,(ev(m),m=0,5)
       !write(*,'(" i1,i2=",2i3," pv=",10(f6.3,1x))') i1,i2,(pv(m),m=0,numPolarizationTerms-1)

        ! FIX ME -- use opt version here too
        pc=0  ! P
        qc=1  ! Q = P.t 
        do k1=1,6
          do k2=1,6
            ec=k2-1 ! E or H component 
            do n=1,Np(k1,k2,mr)
              a0 = gdmPar(1,n,k1,k2,mr)
              a1 = gdmPar(2,n,k1,k2,mr)
              b0 = gdmPar(3,n,k1,k2,mr)
              b1 = gdmPar(4,n,k1,k2,mr)

              ! TZ forcing for polarization equations: 
              fp(pc) = pvt(pc) - pv(qc)
              fp(qc) = pvt(qc) - (  a0*ev(ec) + a1*evt(ec) - b0*pv(pc)- b1*pv(qc) )

              pc=pc+2
              qc=qc+2
            end do
          end do
        end do 
     end if


      ! compute components of the curl(H) and -curl(E)
     #If #DIM eq "2"
       ! --- 2D -----
      #If #ORDER eq "2" 
       #If #POLAR eq "TEZ" 
         curl(0) =                            ETAYOP uy22r(i1,i2,i3,hz)
         curl(1) = -ETAXOP ux22r(i1,i2,i3,hz)
         curl(2) =-(ETAXOP ux22r(i1,i2,i3,ey)-ETAYOP uy22r(i1,i2,i3,ex))
       #Else
         curl(0) =                            ETAYOP uy22r(i1,i2,i3,hz)
         curl(1) = -ETAXOP ux22r(i1,i2,i3,hz)
         curl(2) = (ETAXOP ux22r(i1,i2,i3,hy)-ETAYOP uy22r(i1,i2,i3,hx))
 
         curl(3) =                           -ETAYOP uy22r(i1,i2,i3,ez)
         curl(4) =  ETAXOP ux22r(i1,i2,i3,ez)
         curl(5) =-(ETAXOP ux22r(i1,i2,i3,ey)-ETAYOP uy22r(i1,i2,i3,ex))
       #End   
      #Else
       #If #POLAR eq "TEZ" 
         curl(0) =                            ETAYOP uy42r(i1,i2,i3,hz)
         curl(1) = -ETAXOP ux42r(i1,i2,i3,hz)
         curl(2) =-(ETAXOP ux42r(i1,i2,i3,ey)-ETAYOP uy42r(i1,i2,i3,ex))
       #End
         curl(0) =                            ETAYOP uy42r(i1,i2,i3,hz)
         curl(1) = -ETAXOP ux42r(i1,i2,i3,hz)
         curl(2) = (ETAXOP ux42r(i1,i2,i3,hy)-ETAYOP uy42r(i1,i2,i3,hx))
 
         curl(3) =                           -ETAYOP uy42r(i1,i2,i3,ez)
         curl(4) =  ETAXOP ux42r(i1,i2,i3,ez)
         curl(5) =-(ETAXOP ux42r(i1,i2,i3,ey)-ETAYOP uy42r(i1,i2,i3,ex))
      #End
     #Else
       ! --- 3D -----
      #If #ORDER eq "2" 
        curl(0) =  ETAYOP uy23r(i1,i2,i3,hz)-ETAZOP uz23r(i1,i2,i3,hy)
        curl(1) =  ETAZOP uz23r(i1,i2,i3,hx)-ETAXOP ux23r(i1,i2,i3,hz)
        curl(2) =  ETAXOP ux23r(i1,i2,i3,hy)-ETAYOP uy23r(i1,i2,i3,hx)
  
        curl(3) =-(ETAYOP uy23r(i1,i2,i3,ez)-ETAZOP uz23r(i1,i2,i3,ey))
        curl(4) =-(ETAZOP uz23r(i1,i2,i3,ex)-ETAXOP ux23r(i1,i2,i3,ez))
        curl(5) =-(ETAXOP ux23r(i1,i2,i3,ey)-ETAYOP uy23r(i1,i2,i3,ex))
      #Else
        curl(0) =  ETAYOP uy43r(i1,i2,i3,hz)-ETAZOP uz43r(i1,i2,i3,hy)
        curl(1) =  ETAZOP uz43r(i1,i2,i3,hx)-ETAXOP ux43r(i1,i2,i3,hz)
        curl(2) =  ETAXOP ux43r(i1,i2,i3,hy)-ETAYOP uy43r(i1,i2,i3,hx)
  
        curl(3) =-(ETAYOP uy43r(i1,i2,i3,ez)-ETAZOP uz43r(i1,i2,i3,ey))
        curl(4) =-(ETAZOP uz43r(i1,i2,i3,ex)-ETAXOP ux43r(i1,i2,i3,ez))
        curl(5) =-(ETAXOP ux43r(i1,i2,i3,ey)-ETAYOP uy43r(i1,i2,i3,ex))
      #End
 
     #End  

      if( debug.gt.3 )then
        write(*,'("----- i1,i2=",2i3)') i1,i2
      end if

      ! ---- Compute q = p.t = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, qc )

      ! opt version 
      #If #POLAR eq "TEZ" 
        do m=0,2
          ptSum(m)=0
        end do 
      #Else
        do m=0,5
          ptSum(m)=0
        end do 
      #End 
      do m=1,numTerms1(mr)
         ec = ecIndex1(m,mr)
         qc = qcIndex1(m,mr)
         ! upiv(m) = p(i1,i2,i3,qc-1)  ! didn't seem to help
         ! uqiv(m) = p(i1,i2,i3,qc )
         ! uptSum(ec) = ptSum(ec) + qiv(m) - pv(qc)
         ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)
      end do
      #If #POLAR eq "TEZ" 
        do m=0,2
          curl(m) = curl(m) - ptSum(m)
        end do 
      #Else
        do m=0,5
          curl(m) = curl(m) - ptSum(m)
        end do 
      #End

      ! if( debug.gt.3 )then
      !   write(*,'("   ptSum=",6(f6.3,1x))') (ptSum(m),m=0,5)
      ! end if

      #If #POLAR eq "TEZ" 
        do m=0,2 
          un(i1,i2,i3,m)= Ki(m,0,mr)*(curl(0)+fv(0)) + Ki(m,1,mr)*(curl(1)+fv(1)) + Ki(m,2,mr)*(curl(2)+fv(2))
        end do
      #Else
        do m=0,5 
          un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(2)) + \
                          K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4)) + K0i(m,5,mr)*(curl(5)+fv(5))
        end do
      #End

      ! if( .false. )then
      !   write(*,'(" i1,i2=",2i3," ut=",6(f6.3,1x))') i1,i2,(un(i1,i2,i3,m),m=0,5)
      ! end if

      ! --- compute time derivatives of P and Q
      ! p.t = q
      ! q.t = a0*E + a1*Et - b0*p - b1*q   
      ! or 
      ! q.t = a0*H + a1*Ht - b0*p - b1*q   

      ! optimized version
      do m=1,numTerms2(mr)

        ec = ecIndex2(m,mr)
        pc = pcIndex2(m,mr)
        qc = pc+1
        pct=pc+numComp  ! p.t is stored in un here 
        qct=qc+numComp  ! q.t is stored in un here
        
        a0 = a0v(m,mr)
        a1 = a1v(m,mr)
        b0 = b0v(m,mr)
        b1 = b1v(m,mr)
          
        un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
        un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 

        ! uun(i1,i2,i3,pct) = qiv(m) + fp(pc) 
        ! uun(i1,i2,i3,qct) = a0*uv(ec) + a1*unv(ec) - b0*piv(m)- b1*qiv(m) + fp(qc) 

      end do 

    endLoopsMask()

    if( .false. .or. debug > 15 )then
      stop  4444
   end if
   
  end if 

#endMacro



  
! =========================================================================================
! Optimized version -- non-dispersive
! Macro: Update BI-ANISTROPIC equations
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
!   POLAR : polarization, NONE or TEZ
!   ETAXOP: for supergrid layer in x-direction set to etxa(i1)*, otherwise leave null
! ========================================================================================
#beginMacro updateBAOpt(DIM,ORDER,GRIDTYPE,POLAR,ETAXOP,ETAYOP,ETAZOP)

  if( t.lt.2*dt )then
    write(*,'("advBA: advance BA dim=DIM, order=ORDER grid=GRIDTYPE polar=POLAR... t=",e10.2)') t
  end if
  
  mr=0

  if( .not.methodOfLines )then

    ! --- TAYLOR TIME-STEPPING --- 
    stop 111
     
  else

    ! --- METHOD OF LINES (RK) ---
    beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

      if( addForcing.ne.0 )then  ! do this for now *fix me*
        #If #POLAR eq "TEZ"
          do m=0,2
           fv(m)=f(i1,i2,i3,m)
          end do
        #Else
          do m=0,5
           fv(m)=f(i1,i2,i3,m)
          end do
        #End
      end if 

      if( numberOfMaterialRegions.gt.1 )then
        mr = matMask(i1,i2,i3)
        if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
           stop 9999
        end if            
      end if       

     #If #DIM eq "2"
       ! --- 2D -----
      #If #ORDER eq "2" 
       #If #POLAR eq "TEZ"
         curl(0) =                            ETAYOP uy22r(i1,i2,i3,hz)
         curl(1) = -ETAXOP ux22r(i1,i2,i3,hz)
         curl(2) =-(ETAXOP ux22r(i1,i2,i3,ey)-ETAYOP uy22r(i1,i2,i3,ex))
       #Else
         curl(0) =                            ETAYOP uy22r(i1,i2,i3,hz)
         curl(1) = -ETAXOP ux22r(i1,i2,i3,hz)
         curl(2) = (ETAXOP ux22r(i1,i2,i3,hy)-ETAYOP uy22r(i1,i2,i3,hx))
 
         curl(3) =                           -ETAYOP uy22r(i1,i2,i3,ez)
         curl(4) =  ETAXOP ux22r(i1,i2,i3,ez)
         curl(5) =-(ETAXOP ux22r(i1,i2,i3,ey)-ETAYOP uy22r(i1,i2,i3,ex))
       #End
      #Else
       #If #POLAR eq "TEZ"
         curl(0) =                            ETAYOP uy42r(i1,i2,i3,hz)
         curl(1) = -ETAXOP ux42r(i1,i2,i3,hz)
         curl(2) =-(ETAXOP ux42r(i1,i2,i3,ey)-ETAYOP uy42r(i1,i2,i3,ex))
       #Else
         curl(0) =                            ETAYOP uy42r(i1,i2,i3,hz)
         curl(1) = -ETAXOP ux42r(i1,i2,i3,hz)
         curl(2) = (ETAXOP ux42r(i1,i2,i3,hy)-ETAYOP uy42r(i1,i2,i3,hx))
 
         curl(3) =                           -ETAYOP uy42r(i1,i2,i3,ez)
         curl(4) =  ETAXOP ux42r(i1,i2,i3,ez)
         curl(5) =-(ETAXOP ux42r(i1,i2,i3,ey)-ETAYOP uy42r(i1,i2,i3,ex))
       #End
      #End
     #Else
       ! --- 3D -----
      #If #ORDER eq "2" 
        curl(0) =  ETAYOP uy23r(i1,i2,i3,hz)-ETAZOP uz23r(i1,i2,i3,hy)
        curl(1) =  ETAZOP uz23r(i1,i2,i3,hx)-ETAXOP ux23r(i1,i2,i3,hz)
        curl(2) =  ETAXOP ux23r(i1,i2,i3,hy)-ETAYOP uy23r(i1,i2,i3,hx)
  
        curl(3) =-(ETAYOP uy23r(i1,i2,i3,ez)-ETAZOP uz23r(i1,i2,i3,ey))
        curl(4) =-(ETAZOP uz23r(i1,i2,i3,ex)-ETAXOP ux23r(i1,i2,i3,ez))
        curl(5) =-(ETAXOP ux23r(i1,i2,i3,ey)-ETAYOP uy23r(i1,i2,i3,ex))
      #Else
        curl(0) =  ETAYOP uy43r(i1,i2,i3,hz)-ETAZOP uz43r(i1,i2,i3,hy)
        curl(1) =  ETAZOP uz43r(i1,i2,i3,hx)-ETAXOP ux43r(i1,i2,i3,hz)
        curl(2) =  ETAXOP ux43r(i1,i2,i3,hy)-ETAYOP uy43r(i1,i2,i3,hx)
  
        curl(3) =-(ETAYOP uy43r(i1,i2,i3,ez)-ETAZOP uz43r(i1,i2,i3,ey))
        curl(4) =-(ETAZOP uz43r(i1,i2,i3,ex)-ETAXOP ux43r(i1,i2,i3,ez))
        curl(5) =-(ETAXOP ux43r(i1,i2,i3,ey)-ETAYOP uy43r(i1,i2,i3,ex))
      #End
 
     #End  

      #If #POLAR eq "TEZ"
        do m=0,2 
          un(i1,i2,i3,m)= Ki(m,0,mr)*(curl(0)+fv(0)) + Ki(m,1,mr)*(curl(1)+fv(1)) + Ki(m,2,mr)*(curl(2)+fv(2))
        end do
      #Else
        do m=0,5 
          un(i1,i2,i3,m)= K0i(m,0,mr)*(curl(0)+fv(0)) + K0i(m,1,mr)*(curl(1)+fv(1)) + K0i(m,2,mr)*(curl(2)+fv(2)) + \
                          K0i(m,3,mr)*(curl(3)+fv(3)) + K0i(m,4,mr)*(curl(4)+fv(4)) + K0i(m,5,mr)*(curl(5)+fv(5))
        end do
     #End 
 
    endLoopsMask()

    if( .false. .or. debug > 15 )then
      stop  4747
   end if
   
  end if ! end MOL

#endMacro


! =========================================================================
!  Macro to call super-grid for difference cases, depending on which axis
!  has an absorning layer   
! =========================================================================
#beginMacro updateSuperGrid(updateRoutine,DIM,ORDER,GRIDTYPE,POLAR)
  #If #DIM eq "2" 
    !  --- TWO-DIMENSIONS ---
    if( useAbsorbingLayer(0).eq.1 .and. useAbsorbingLayer(1).eq.1 )then
      updateRoutine(DIM,ORDER,GRIDTYPE,POLAR,etax(i1)*,etay(i2)*,)

    else if( useAbsorbingLayer(0).eq.1 .and. useAbsorbingLayer(1).eq.0 )then
      updateRoutine(DIM,ORDER,GRIDTYPE,POLAR,etax(i1)*,,)

    else if( useAbsorbingLayer(0).eq.0 .and. useAbsorbingLayer(1).eq.1 )then
      updateRoutine(DIM,ORDER,GRIDTYPE,POLAR,,etay(i2)*,)

    else
      updateRoutine(DIM,ORDER,GRIDTYPE,POLAR,,,)
    end if 
  #Else
    !  --- THREE-DIMENSIONS ---
    if( useAbsorbingLayer(0).eq.1 .and. useAbsorbingLayer(1).eq.1 .and. useAbsorbingLayer(2).eq.1 )then
      updateRoutine(DIM,ORDER,GRIDTYPE,POLAR,etax(i1)*,etay(i2)*,etaz(i3)*)

    else if( useAbsorbingLayer(0).eq.1 .and. useAbsorbingLayer(1).eq.1 .and. useAbsorbingLayer(2).eq.0 )then
      updateRoutine(DIM,ORDER,GRIDTYPE,POLAR,etax(i1)*,etay(i2)*,)

    else if( useAbsorbingLayer(0).eq.1 .and. useAbsorbingLayer(1).eq.0 .and. useAbsorbingLayer(2).eq.1 )then
      updateRoutine(DIM,ORDER,GRIDTYPE,POLAR,etax(i1)*,,etaz(i3)*)

    else if( useAbsorbingLayer(0).eq.0 .and. useAbsorbingLayer(1).eq.1 .and. useAbsorbingLayer(2).eq.1 )then
      updateRoutine(DIM,ORDER,GRIDTYPE,POLAR,,etay(i2)*,etaz(i3)*)

    else if( useAbsorbingLayer(0).eq.1 .and. useAbsorbingLayer(1).eq.0 .and. useAbsorbingLayer(2).eq.0 )then
      updateRoutine(DIM,ORDER,GRIDTYPE,POLAR,etax(i1)*,,)

    else if( useAbsorbingLayer(0).eq.0 .and. useAbsorbingLayer(1).eq.1 .and. useAbsorbingLayer(2).eq.0 )then
      updateRoutine(DIM,ORDER,GRIDTYPE,POLAR,,etay(i2)*,)

    else if( useAbsorbingLayer(0).eq.0 .and. useAbsorbingLayer(1).eq.0 .and. useAbsorbingLayer(2).eq.1 )then
      updateRoutine(DIM,ORDER,GRIDTYPE,POLAR,,,etaz(i3)*)

    else
      updateRoutine(DIM,ORDER,GRIDTYPE,POLAR,,,)
    end if 
  #End
#endMacro   


! =========================================================================
!  Macro to call super-grid or non-supergrid version of BA update 
! =========================================================================
#beginMacro updateBA(DIM,ORDER,GRIDTYPE,POLAR)
  if( dispersionModel.eq.noDispersion )then

    if( useSuperGrid.eq.0  ) then
      updateBAOpt(DIM,ORDER,GRIDTYPE,POLAR,,,)
    else
      ! --- SUPERGRID ---
      if( t.le.3*dt )then
        write(*,'(" USE SUPERGRID...")' )
      end if

      updateSuperGrid(updateBAOpt,DIM,ORDER,GRIDTYPE,POLAR)

    end if

  else

    if( useSuperGrid.eq.0  ) then
      updateBAGDMOpt(DIM,ORDER,GRIDTYPE,POLAR,,,)
    else
      ! --- SUPERGRID ---
      if( t.le.3*dt )then
        write(*,'(" USE SUPERGRID...")' )
      end if

      updateSuperGrid(updateBAGDMOpt,DIM,ORDER,GRIDTYPE,POLAR)

    end if
  end if   
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
   mask,rsxy,  um,u,un,f,fa, K0i, matMask, \
   pm,p,pt,xy, etax,etay,etaz, bc, dis, varDis, ipar, rpar, ierr )
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

 ! Polarization vectors
 real pm(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
 real p(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
 real pt(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)

 real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)

 ! Super-grid layer functions 
 real etax(nd1a:nd1b), etay(nd2a:nd2b), etaz(nd3a:nd3b)
 
! real vvt2(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
! real ut3(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
! real vvt4(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
!  real ut5(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
!  real ut6(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
!  real ut7(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 real dis(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 real varDis(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
 real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)

 integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)

 real K0i(0:5,0:5,0:*)  ! material matrix 
 integer matMask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)

 integer bc(0:1,0:2),ierr

 integer ipar(0:*)
 real rpar(0:*)
      
!     ---- local variables -----
 integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderInTime,k1,k2,debug,numComp
 integer addForcing,orderOfDissipation,option
 integer useWhereMask,useWhereMaskSave,solveForE,solveForH,grid,useVariableDissipation
 integer useCurvilinearOpt,useConservative,combineDissipationWithAdvance,useDivergenceCleaning
 integer useNewForcingMethod,numberOfForcingFunctions,fcur,fnext,fprev
 integer ex,ey,ez, hx,hy,hz, solveForAllFields, useSuperGrid, useAbsorbingLayer(0:2)
 real t,cc,dt,dy,dz,cdt,cdtdx,cdtdy,cdtdz,adc,adcdt,add,adddt
 real dt4by12
 real eps,mu,sigmaE,sigmaH,kx,ky,kz,divergenceCleaningCoefficient
 logical addDissipation, updateInterior, methodOfLines
 real fv(0:10) , curl(0:5)
 real ep 

 integer maxRegions,NpMax
 parameter( maxRegions=100,NpMax=10 )  ! FIX ME 

 integer numberOfMaterialRegions, mr
 real Ki(0:2,0:2,0:maxRegions) ! 3x3 material matrix for TEz polarization

 integer Np(6,6,0:maxRegions-1)
 real gdmPar(4,NpMax,6,6,0:maxRegions-1), ptSum(0:5) 
 real a0,a1,b0,b1
 integer ec,pc,qc, pct,qct

 real e0,e0t,p0,p0t,q0,q0t
 real ev(0:5),evt(0:5)
 integer numPolarizationTerms
 integer maxNumPolarizationTerms
 parameter( maxNumPolarizationTerms=200 )
 real pv(0:maxNumPolarizationTerms), pvt(0:maxNumPolarizationTerms) ! 
 real fp(0:maxNumPolarizationTerms)

 real dx(0:2),dr(0:2)
 real adxSosup(0:2), sigma1

 real dx2i,dy2i,dz2i,dxsqi,dysqi,dzsqi,dxi,dyi,dzi
 real dx12i,dy12i,dz12i,dxsq12i,dysq12i,dzsq12i,dxy4i,dxz4i,dyz4,time0,time1

 real dxi4,dyi4,dzi4,dxdyi2,dxdzi2,dydzi2

 real uLap(-1:1,-1:1,0:5),uLapSq(0:5)
 real uLaprr2,uLapss2,uLaprs2,uLapr2,uLaps2

 real c0,c1,csq,dtsq,cdtsq,cdtsq12,cdtSqBy12, csqdt
!  real c0,c1,csq,dtsq,cdtsq,cdtsq12,lap(0:20),cdtSqBy12, csqdt
!  real c40,c41,c42,c43
!  real c60,c61,c62,c63,c64,c65
!  real c80,c81,c82,c83,c84,c85,c86,c87


 
! real c00lap2d6,c10lap2d6,c01lap2d6,c20lap2d6,c02lap2d6,c30lap2d6,c03lap2d6
!  real c00lap2d8,c10lap2d8,c01lap2d8,c20lap2d8,c02lap2d8,c30lap2d8,c03lap2d8,c40lap2d8,c04lap2d8
!  real c000lap3d6,c100lap3d6,c010lap3d6,c001lap3d6,\
!                  c200lap3d6,c020lap3d6,c002lap3d6,\
!                  c300lap3d6,c030lap3d6,c003lap3d6
!  real c000lap3d8,c100lap3d8,c010lap3d8,c001lap3d8,\
!                  c200lap3d8,c020lap3d8,c002lap3d8,\
!                  c300lap3d8,c030lap3d8,c003lap3d8,\
!                  c400lap3d8,c040lap3d8,c004lap3d8

 integer rectangular,curvilinear
 parameter( rectangular=0, curvilinear=1 )

 integer timeSteppingMethod
 integer defaultTimeStepping,adamsSymmetricOrder3,rungeKutta,\
         stoermerTimeStepping,modifiedEquationTimeStepping
 parameter(defaultTimeStepping=0,adamsSymmetricOrder3=1,\
           rungeKutta=2,stoermerTimeStepping=3,modifiedEquationTimeStepping=4)

 integer materialType
 integer isotropic,bianisotropic
 parameter( isotropic=0, bianisotropic=1 )

!...........start statement function
 integer kd,m
 real rx,ry,rz,sx,sy,sz,tx,ty,tz

 declareDifferenceOrder2(u,RX)
 declareDifferenceOrder2(un,none)
 declareDifferenceOrder2(v,none)
 declareDifferenceOrder2(um,none)
 declareDifferenceOrder2(ff,none)

 declareDifferenceOrder4(u,RX)
 declareDifferenceOrder4(un,none)
 declareDifferenceOrder4(v,none)

 real maxwell2dr,maxwell3dr,maxwellr44,maxwellr66,maxwellr88
 real maxwellc22,maxwellc44,maxwellc66,maxwellc88
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

! real vr2,vs2,vrr2,vss2,vrs2,vLaplacian22

 real cdt4by360,cdt6by20160

 real lap2d2,lap3d2,lap2d4,lap3d4,lap2d6,lap3d6,lap2d8,lap3d8,lap2d2Pow2,lap3d2Pow2,lap2d2Pow3,lap3d2Pow3,\
      lap2d2Pow4,lap3d2Pow4,lap2d4Pow2,lap3d4Pow2,lap2d4Pow3,lap3d4Pow3,lap2d6Pow2,lap3d6Pow2
 real lap2d2m,lap3d2m
 real du,fd22d,fd23d,fd42d,fd43d,fd62d,fd63d,fd82d,fd83d

 ! forcing correction functions: 
 real lap2d2f,f2drme44, lap3d2f, f3drme44, f2dcme44, f3dcme44, ff

 ! div cleaning: 
 real dc,dcp,cdc0,cdc1,cdcxx,cdcyy,cdczz,cdcEdx,cdcEdy,cdcEdz,cdcHdx,cdcHdy,cdcHdz,cdcf
 real cdcE,cdcELap,cdcELapsq,cdcELapm,cdcHzxLap,cdcHzyLap
 real cdcH,cdcHLap,cdcHLapsq,cdcHLapm

 ! dispersion
 integer dispersionModel,pxc,pyc,pzc,qxc,qyc,qzc,rxc,ryc,rzc

 integer numTerms1(0:maxRegions),ecIndex1(maxNumPolarizationTerms,0:maxRegions),qcIndex1(maxNumPolarizationTerms,0:maxRegions)
 integer numTerms2(0:maxRegions),ecIndex2(maxNumPolarizationTerms,0:maxRegions),pcIndex2(maxNumPolarizationTerms,0:maxRegions)
 real a0v(maxNumPolarizationTerms,0:maxRegions),a1v(maxNumPolarizationTerms,0:maxRegions)
 real b0v(maxNumPolarizationTerms,0:maxRegions),b1v(maxNumPolarizationTerms,0:maxRegions)
 ! real uv(0:5), unv(0:5)
 logical useOpt

! Dispersion models
 #Include "dispersionModelsFortranInclude.h"

 integer forcingOption   
 ! forcing options
 #Include "forcingDefineFortranInclude.h"

 ! boundary conditions parameters
 ! #Include "bcDefineFortranInclude.h"

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

 ! efineDifferenceOrder2Components1(v,none)
 ! defineDifferenceOrder4Components1(v,none)

 ! defineDifferenceOrder2Components1(um,none)

 ! defineDifferenceOrder2Components1(ff,none)

 !    *** 2nd order ***
 lap2d2(i1,i2,i3,c)=(u(i1+1,i2,i3,c)-2.*u(i1,i2,i3,c)+u(i1-1,i2,i3,c))*dxsqi\
                   +(u(i1,i2+1,i3,c)-2.*u(i1,i2,i3,c)+u(i1,i2-1,i3,c))*dysqi
 
!-  lap2d2m(i1,i2,i3,c)=(um(i1+1,i2,i3,c)-2.*um(i1,i2,i3,c)+um(i1-1,i2,i3,c))*dxsqi\
!-                     +(um(i1,i2+1,i3,c)-2.*um(i1,i2,i3,c)+um(i1,i2-1,i3,c))*dysqi
!- 
!-  lap3d2(i1,i2,i3,c)=(u(i1+1,i2,i3,c)-2.*u(i1,i2,i3,c)+u(i1-1,i2,i3,c))*dxsqi\
!-                    +(u(i1,i2+1,i3,c)-2.*u(i1,i2,i3,c)+u(i1,i2-1,i3,c))*dysqi\
!-                    +(u(i1,i2,i3+1,c)-2.*u(i1,i2,i3,c)+u(i1,i2,i3-1,c))*dzsqi
!- 
!-  lap3d2m(i1,i2,i3,c)=(um(i1+1,i2,i3,c)-2.*um(i1,i2,i3,c)+um(i1-1,i2,i3,c))*dxsqi\
!-                     +(um(i1,i2+1,i3,c)-2.*um(i1,i2,i3,c)+um(i1,i2-1,i3,c))*dysqi\
!-                     +(um(i1,i2,i3+1,c)-2.*um(i1,i2,i3,c)+um(i1,i2,i3-1,c))*dzsqi
!- 
!-  lap2d2f(i1,i2,i3,c,m)=(fa(i1+1,i2,i3,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1-1,i2,i3,c,m))*dxsqi\
!-                       +(fa(i1,i2+1,i3,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1,i2-1,i3,c,m))*dysqi
!- 
!-  lap3d2f(i1,i2,i3,c,m)=(fa(i1+1,i2,i3,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1-1,i2,i3,c,m))*dxsqi\
!-                       +(fa(i1,i2+1,i3,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1,i2-1,i3,c,m))*dysqi\
!-                       +(fa(i1,i2,i3+1,c,m)-2.*fa(i1,i2,i3,c,m)+fa(i1,i2,i3-1,c,m))*dzsqi
!- 
!- 
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
!- 
!-  ! 3D laplacian squared = u.xxxx + u.yyyy + u.zzzz + 2 (u.xxyy + u.xxzz + u.yyzz )
!-  lap3d2Pow2(i1,i2,i3,c)= ( 6.*u(i1,i2,i3,c)   \
!-    - 4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))    \
!-        +(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) )*dxi4 \
!-   +(  +6.*u(i1,i2,i3,c)    \
!-     -4.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))    \
!-        +(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) )*dyi4\
!-   +(  +6.*u(i1,i2,i3,c)    \
!-     -4.*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))    \
!-        +(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) )*dzi4\
!-    +(8.*u(i1,i2,i3,c)     \
!-     -4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)+u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))   \
!-     +2.*(u(i1+1,i2+1,i3,c)+u(i1-1,i2+1,i3,c)+u(i1+1,i2-1,i3,c)+u(i1-1,i2-1,i3,c)) )*dxdyi2 \
!-    +(8.*u(i1,i2,i3,c)     \
!-     -4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)+u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))   \
!-     +2.*(u(i1+1,i2,i3+1,c)+u(i1-1,i2,i3+1,c)+u(i1+1,i2,i3-1,c)+u(i1-1,i2,i3-1,c)) )*dxdzi2 \
!-    +(8.*u(i1,i2,i3,c)     \
!-     -4.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)+u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))   \
!-     +2.*(u(i1,i2+1,i3+1,c)+u(i1,i2-1,i3+1,c)+u(i1,i2+1,i3-1,c)+u(i1,i2-1,i3-1,c)) )*dydzi2 
!- 
!-  lap2d2Pow3(i1,i2,i3,c)=LAP2D2(lap2d2Pow2,i1,i2,i3,c)
!- 
!-  lap3d2Pow3(i1,i2,i3,c)=LAP3D2(lap3d2Pow2,i1,i2,i3,c)
!- 
!-  lap2d2Pow4(i1,i2,i3,c)=LAP2D2POW2(lap2d2Pow2,i1,i2,i3,c)
!-  lap3d2Pow4(i1,i2,i3,c)=LAP3D2POW2(lap3d2Pow2,i1,i2,i3,c)
!-  

!    ** 4th order ****

lap2d4(i1,i2,i3,c)=( -30.*u(i1,i2,i3,c)     \
  +16.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))     \
      -(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) )*dxsq12i + \
 ( -30.*u(i1,i2,i3,c)     \
  +16.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))     \
      -(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) )*dysq12i 
!- 
!-  lap3d4(i1,i2,i3,c)=lap2d4(i1,i2,i3,c)+ \
!-   ( -30.*u(i1,i2,i3,c)      \
!-    +16.*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))      \
!-        -(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) )*dzsq12i 
!- 
!-  lap2d4Pow2(i1,i2,i3,c)=LAP2D4(lap2d4,i1,i2,i3,c)
!-  lap3d4Pow2(i1,i2,i3,c)=LAP3D4(lap3d4,i1,i2,i3,c)
!- 
!-  lap2d4Pow3(i1,i2,i3,c)=LAP2D4(lap2d4Pow2,i1,i2,i3,c)
!-  lap3d4Pow3(i1,i2,i3,c)=LAP3D4(lap3d4Pow2,i1,i2,i3,c)
!- 
!- !     *** 6th order ***
!- 
!-  lap2d6(i1,i2,i3,c)= \
!-           c00lap2d6*u(i1,i2,i3,c)     \
!-          +c10lap2d6*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)) \
!-          +c01lap2d6*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)) \
!-          +c20lap2d6*(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) \
!-          +c02lap2d6*(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) \
!-          +c30lap2d6*(u(i1+3,i2,i3,c)+u(i1-3,i2,i3,c)) \
!-          +c03lap2d6*(u(i1,i2+3,i3,c)+u(i1,i2-3,i3,c)) 
!- 
!-  lap3d6(i1,i2,i3,c)=\
!-           c000lap3d6*u(i1,i2,i3,c) \
!-          +c100lap3d6*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)) \
!-          +c010lap3d6*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)) \
!-          +c001lap3d6*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c)) \
!-          +c200lap3d6*(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) \
!-          +c020lap3d6*(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) \
!-          +c002lap3d6*(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) \
!-          +c300lap3d6*(u(i1+3,i2,i3,c)+u(i1-3,i2,i3,c)) \
!-          +c030lap3d6*(u(i1,i2+3,i3,c)+u(i1,i2-3,i3,c)) \
!-          +c003lap3d6*(u(i1,i2,i3+3,c)+u(i1,i2,i3-3,c))
!- 
!-  lap2d6Pow2(i1,i2,i3,c)=LAP2D6(lap2d6,i1,i2,i3,c)
!-  lap3d6Pow2(i1,i2,i3,c)=LAP3D6(lap3d6,i1,i2,i3,c)
!- 
!- 
!- !     *** 8th order ***
!- 
!-  lap2d8(i1,i2,i3,c)=c00lap2d8*u(i1,i2,i3,c)      \
!-           +c10lap2d8*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))     \
!-           +c01lap2d8*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)) \
!-           +c20lap2d8*(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c))  \
!-           +c02lap2d8*(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) \
!-           +c30lap2d8*(u(i1+3,i2,i3,c)+u(i1-3,i2,i3,c))  \
!-           +c03lap2d8*(u(i1,i2+3,i3,c)+u(i1,i2-3,i3,c)) \
!-           +c40lap2d8*(u(i1+4,i2,i3,c)+u(i1-4,i2,i3,c))  \
!-           +c04lap2d8*(u(i1,i2+4,i3,c)+u(i1,i2-4,i3,c))
!- 
!-  lap3d8(i1,i2,i3,c)=c000lap3d8*u(i1,i2,i3,c)      \
!-           +c100lap3d8*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))     \
!-           +c010lap3d8*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)) \
!-           +c001lap3d8*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c)) \
!-           +c200lap3d8*(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c))  \
!-           +c020lap3d8*(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) \
!-           +c002lap3d8*(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) \
!-           +c300lap3d8*(u(i1+3,i2,i3,c)+u(i1-3,i2,i3,c))  \
!-           +c030lap3d8*(u(i1,i2+3,i3,c)+u(i1,i2-3,i3,c)) \
!-           +c003lap3d8*(u(i1,i2,i3+3,c)+u(i1,i2,i3-3,c)) \
!-           +c400lap3d8*(u(i1+4,i2,i3,c)+u(i1-4,i2,i3,c))  \
!-           +c040lap3d8*(u(i1,i2+4,i3,c)+u(i1,i2-4,i3,c)) \
!-           +c004lap3d8*(u(i1,i2,i3+4,c)+u(i1,i2,i3-4,c))
!- 
!- ! ******* artificial dissipation ******
!-  du(i1,i2,i3,c)=u(i1,i2,i3,c)-um(i1,i2,i3,c)
!- 
!- !      (2nd difference)
!-  fd22d(i1,i2,i3,c)= \
!-  (     ( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c) ) \
!-   -4.*du(i1,i2,i3,c) )
!- !
!-  fd23d(i1,i2,i3,c)=\
!-  (     ( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c)+du(i1,i2,i3-1,c)+du(i1,i2,i3+1,c) ) \
!-    -6.*du(i1,i2,i3,c) )
!- 
!- !     -(fourth difference)
!-  fd42d(i1,i2,i3,c)= \
!-  (    -( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c) ) \
!-    +4.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c) ) \
!-   -12.*du(i1,i2,i3,c) )
!- !
!-  fd43d(i1,i2,i3,c)=\
!-  (    -( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c)+du(i1,i2,i3-2,c)+du(i1,i2,i3+2,c) ) \
!-    +4.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c)+du(i1,i2,i3-1,c)+du(i1,i2,i3+1,c) ) \
!-   -18.*du(i1,i2,i3,c) )
!- 
!-  ! (sixth  difference)
!-  fd62d(i1,i2,i3,c)= \
!-  (     ( du(i1-3,i2,i3,c)+du(i1+3,i2,i3,c)+du(i1,i2-3,i3,c)+du(i1,i2+3,i3,c) ) \
!-    -6.*( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c) ) \
!-   +15.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c) ) \
!-   -40.*du(i1,i2,i3,c) )
!- 
!-  fd63d(i1,i2,i3,c)=\
!-  (     ( du(i1-3,i2,i3,c)+du(i1+3,i2,i3,c)+du(i1,i2-3,i3,c)+du(i1,i2+3,i3,c)+du(i1,i2,i3-3,c)+du(i1,i2,i3+3,c) ) \
!-    -6.*( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c)+du(i1,i2,i3-2,c)+du(i1,i2,i3+2,c) ) \
!-   +15.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c)+du(i1,i2,i3-1,c)+du(i1,i2,i3+1,c) ) \
!-   -60.*du(i1,i2,i3,c) )
!- 
!-  ! -(eighth  difference)
!-  fd82d(i1,i2,i3,c)= \
!-  (    -( du(i1-4,i2,i3,c)+du(i1+4,i2,i3,c)+du(i1,i2-4,i3,c)+du(i1,i2+4,i3,c) ) \
!-    +8.*( du(i1-3,i2,i3,c)+du(i1+3,i2,i3,c)+du(i1,i2-3,i3,c)+du(i1,i2+3,i3,c) ) \
!-   -28.*( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c) ) \
!-   +56.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c) ) \
!-  -140.*du(i1,i2,i3,c) )
!- 
!-  fd83d(i1,i2,i3,c)=\
!-  (    -( du(i1-4,i2,i3,c)+du(i1+4,i2,i3,c)+du(i1,i2-4,i3,c)+du(i1,i2+4,i3,c)+du(i1,i2,i3-4,c)+du(i1,i2,i3+4,c) ) \
!-    +8.*( du(i1-3,i2,i3,c)+du(i1+3,i2,i3,c)+du(i1,i2-3,i3,c)+du(i1,i2+3,i3,c)+du(i1,i2,i3-3,c)+du(i1,i2,i3+3,c) ) \
!-   -28.*( du(i1-2,i2,i3,c)+du(i1+2,i2,i3,c)+du(i1,i2-2,i3,c)+du(i1,i2+2,i3,c)+du(i1,i2,i3-2,c)+du(i1,i2,i3+2,c) ) \
!-   +56.*( du(i1-1,i2,i3,c)+du(i1+1,i2,i3,c)+du(i1,i2-1,i3,c)+du(i1,i2+1,i3,c)+du(i1,i2,i3-1,c)+du(i1,i2,i3+1,c) ) \
!-  -210.*du(i1,i2,i3,c) )


!...........end   statement functions


 useOpt=.true. ! is true, use new optimized version

 debug=0

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
 ep    =rpar(19)  ! for TZ forcing  

 rpar(20)=0.  ! return the time used for adding dissipation

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
 ! qxc                 =ipar(28)
 grid                =ipar(29)
 ! qzc                 =ipar(30)
 ! rxc                 =ipar(31)
 ! ryc                 =ipar(32)
 ! rzc                 =ipar(33)

 ! useSosupDissipation,
 ! sosupDissipationOption,
 updateInterior        =ipar(36).eq.1
 addDissipation        =ipar(37).eq.1

 forcingOption         =ipar(39) 

 ! addDissipation=.true. if we add the dissipation in the dis(i1,i2,i3,c) array
 !  if combineDissipationWithAdvance.ne.0 we compute the dissipation on the fly in the time step
 !  rather than pre-computing it in diss(i1,i2,i3,c)
 ! OLD addDissipation = adc.gt.0. .and. combineDissipationWithAdvance.eq.0
 ! adcdt=adc*dt


 solveForAllFields       = ipar(40)
 materialType            = ipar(41) 
 numberOfMaterialRegions = ipar(42)
 useSuperGrid            = ipar(43) 
 useAbsorbingLayer(0)    = ipar(44)
 useAbsorbingLayer(1)    = ipar(45)
 useAbsorbingLayer(2)    = ipar(46)
 
 if( nd.eq.2 .and. solveForAllFields.eq.0 )then
   numComp=3  ! TEZ has 3 components 
 else
   numComp=6
 end if
    

 if( t.le.2*dt )then
   write(*,*) 'Inside NAME...'
   write(*,'("addForcing=",i2, " solveForAllFields=",i2," useSuperGrid=",i2)') addForcing,solveForAllFields
   write(*,'(" useSuperGrid=",i2," useAbsorbingLayer(0:2)=",3i2)') useSuperGrid,(useAbsorbingLayer(n),n=0,2)
   write(*,'("addDissipation=",l2, " updateInterior=",l2)') addDissipation,updateInterior
   write(*,'("dispersionModel=",i2)') dispersionModel

   if( timeSteppingMethod.eq.modifiedEquationTimeStepping )then
      write(*,*) 'timeSteppingMethod =  modifiedEquationTimeStepping'
   else if( timeSteppingMethod.eq.rungeKutta ) then
      write(*,*) 'timeSteppingMethod = rungeKutta'
   end if

   if( .false. .and. useSuperGrid.ne.0 )then
     write(*,'(" etax=",100(f6.3,1x))') (etax(i1),i1=nd1a,nd1b)
  end if
  
 end if


 if( numberOfMaterialRegions.gt.maxRegions )then
    stop 1002
 end if
 ! 3x3 Material matrix for TEz polarization
 ! We use ex=0,ey=1 and hz=5 entries in K0i(0:5,0:5) 
 do mr=0,numberOfMaterialRegions-1
   Ki(0,0,mr) = K0i(0,0,mr)
   Ki(0,1,mr) = K0i(0,1,mr)
   Ki(0,2,mr) = K0i(0,5,mr)

   Ki(1,0,mr) = K0i(1,0,mr)
   Ki(1,1,mr) = K0i(1,1,mr)
   Ki(1,2,mr) = K0i(1,5,mr)

   Ki(2,0,mr) = K0i(5,0,mr)
   Ki(2,1,mr) = K0i(5,1,mr)
   Ki(2,2,mr) = K0i(5,5,mr)
 end do

 if( t.lt.dt )then
  
   write(*,*) 'materialType=',materialType
   write(*,*) 'numberOfMaterialRegions=',numberOfMaterialRegions
   do mr=0,numberOfMaterialRegions-1
     write(*,*) 'Material region=',mr 
     write(*,'("K0i=",6("[",6(f6.3,1x),"]",/,4x))') ((K0i(i1,i2,mr),i1=0,5),i2=0,5)
     if( solveForAllFields .eq. 0 )then
       write(*,'("Ki =",3("[",3(f6.3,1x),"]",/,4x))') ((Ki(i1,i2,mr),i1=0,2),i2=0,2)
     end if
   end do
   if( numberOfMaterialRegions>1 .and. n1b-n1a.lt.20 )then
     do i3=n3a,n3b
     do i2=n2a,n2b
       write(*,*) 'matMask=',(matMask(i1,i2,i3),i1=n1a,n1b)
     end do
     end do
  
   end if 
   
 end if

 if( dispersionModel.ne.noDispersion )then
   ! get the BA gdm parameters
   !  gdmPar(4,NpMax,6,6,0:maxRegions-1)


   call getBAGDMParameters( grid, gdmPar,Np, NpMax,maxRegions )

   ! count the total number of polarization terms: 
   numPolarizationTerms=0 
   do mr=0,numberOfMaterialRegions-1 
     do k1=1,6
       do k2=1,6
         numPolarizationTerms = numPolarizationTerms + Np(k1,k2,mr)
       end do
     end do
   end do
   numPolarizationTerms = numPolarizationTerms*2  ! we store p and p.t
   if( numPolarizationTerms > maxNumPolarizationTerms )then
      write(*,'("advBA: ERROR: numPolarizationTerms > maxNumPolarizationTerms")');
      stop 1234
   end if      

   if( t.eq.0. .and. dispersionModel.ne.noDispersion )then
     ! ---- Dispersive Maxwell ----
     write(*,'("--advOpt-- dispersionModel=",i4," numPolarizationTerms=",i6)') dispersionModel,numPolarizationTerms
     !write(*,'("--advOpt-- GDM: numberOfPolarizationVectors=",i4," alphaP=",e8.2)') numberOfPolarizationVectors,alphaP
     !write(*,'("--advOpt-- GDM: alphaP,a0,a1,b0,b1=",5(1p,e10.2))') alphaP,a0,a1,b0,b1

     !do iv=0,numberOfPolarizationVectors-1
     !  write(*,'("--advOpt-- GDM: eqn=",i3," a0,a1,b0,b1=",4(1p,e10.2))') iv,a0v(iv),a1v(iv),b0v(iv),b1v(iv)
     !end do

     if( .false. )then
      do mr=0,numberOfMaterialRegions-1 
        write(*,'("BA-GDM: material region mr=",i2)') mr
        ! write(*,'((5x,6i3))') ((Np(k1,k2,mr),k1=1,6),k2=1,6)
        do k1=1,6
        do k2=1,6
           if( Np(k1,k2,mr) > 0 )then
             write(*,'("  K(",i1,",",i1,") : Np=",i3," :")') k1,k2,Np(k1,k2,mr)
             do n=1,Np(k1,k2,mr)
               write(*,'("    n=",i3," [a0,a1,b0,b1]=",4(e9.3,1x))') n,(gdmPar(m,n,k1,k2,mr),m=1,4)
             end do 
           end if
        end do           
        end do
      end do
     end if
    
    ! stop 3333

  end if
 end if


 methodOfLines=.false.
 if( timeSteppingMethod.ne.modifiedEquationTimeStepping )then
    methodOfLines=.true.
 end if    



 fprev = mod(fcur-1+numberOfForcingFunctions,max(1,numberOfForcingFunctions))
 fnext = mod(fcur+1                         ,max(1,numberOfForcingFunctions))
 


 csq=cc**2
 dtsq=dt**2

 cdt=cc*dt

 cdtsq=(cc**2)*(dt**2)
 cdtsq12=cdtsq*cdtsq/12.  ! c^4 dt^4 /14 
 cdt4by360=(cdt)**4/360.
 cdt6by20160=cdt**6/(8.*7.*6.*5.*4.*3.)

 csqdt = cc**2*dt

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

 do m=0,10
    fv(m)=0.  ! temp for forcing
 end do    

 if( t.eq.0. .and. dispersionModel.ne.noDispersion )then
    write(*,'("--advOpt-- dispersionModel=",i4," px,py,pz=",3i2)') dispersionModel,pxc,pyc,pzc
 end if


! write(*,'(" advMaxwell: timeSteppingMethod=",i2)') timeSteppingMethod
 if( timeSteppingMethod.eq.defaultTimeStepping )then
  write(*,'(" advMaxwell:ERROR: timeSteppingMethod=defaultTimeStepping -- this should be set")')
    ! '
  stop 83322
 end if


 if( .not. addDissipation )then
   ! **new way ** Dec 17, 2019 
   ! All methods are done here 
   #If #DIM eq "2"
     if( solveForAllFields.eq.0 )then
       ! TEZ polarization
       updateBA(DIM,ORDER,GRIDTYPE,TEZ)       
     else
       updateBA(DIM,ORDER,GRIDTYPE,NONE)       
     end if 
   #Else
     updateBA(DIM,ORDER,GRIDTYPE,NONE)       
   #End
 end if 


 if( gridType.eq.rectangular )then

 #If #GRIDTYPE eq "rectangular"

!       **********************************************
!       *************** rectangular ******************
!       **********************************************


  #If #ORDER eq "2" 

   #If #DIM eq "2"

    ! ------------------------------------------------------------------------------
    !    2D : 2nd order (rectangular)
    ! ------------------------------------------------------------------------------
    
    if( addDissipation )then

      if( t.le.3*dt )then
        write(*,'("advBA: order=2: addDissipation=",l2," adc=",e10.2)') addDissipation,adc
      end if
     
      sigma1=1.
      adxSosup(0)=sigma1/(2*16.)
      adxSosup(1)=sigma1/(2*16.)
      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
        do m=ex,hz
          un(i1,i2,i3,m)=u(i1,i2,i3,m) + sosupDiss2d4(u,i1,i2,i3,m)
        end do         
      endLoopsMask()

!     else if( materialType.eq.isotropic )then
!
!       updateMx2dOrder2()
!
!     else if( dispersionModel.eq.noDispersion )then
!       updateBA2dOrder2()
!     else
!       ! updateBAGDM2dOrder2()
!       ! updateBAGDM(DIM,ORDER,GRIDTYPE)
!      if( useOpt )then
!        if( useSuperGrid.eq.0  ) then
!          updateBAGDMOpt(2,2,rectangular,,,)
!         else
!          ! --- SUPERGRID ---
!          write(*,'(" USE SUPERGRID...")' )
!          updateBAGDMOpt(2,2,rectangular,etax(i1)*,etay(i2)*,etaz(i3)*)
!         end if
!      else 
!        updateBAGDM(2,2,rectangular)
!      end if 
    end if 

   #Else

    ! ------------------------------------------------------------------------------
    !    3D : Order 2 (rectangular)
    ! ------------------------------------------------------------------------------

    if( addDissipation )then

      if( t.le.3*dt )then
        write(*,'("advBA: order=2: addDissipation=",l2," adc=",e10.2)') addDissipation,adc
      end if        
      sigma1=1.
      adxSosup(0)=sigma1/(2*16.)
      adxSosup(1)=sigma1/(2*16.)
      adxSosup(2)=sigma1/(2*16.)
      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
        do m=ex,hz
          un(i1,i2,i3,m)=u(i1,i2,i3,m) + sosupDiss3d4(u,i1,i2,i3,m)
        end do         
      endLoopsMask()
    end if 

!     else if( materialType.eq.isotropic )then
!       ! updateMx2dOrder2()
!       stop 3333
!     else if( dispersionModel.eq.noDispersion )then
!       updateBA3dOrder2()
!     else
!      if( useOpt )then
!        if( useSuperGrid.eq.0  ) then
!          updateBAGDMOpt(3,2,rectangular,,,)
!        else
!          updateBAGDMOpt(3,2,rectangular,etax(i1)*,etay(i2)*,etaz(i3)*)
!        end if 
!      else 
!        updateBAGDM(3,2,rectangular)
!      end if 
!     end if 

   #End

  #Elif #ORDER eq "4" 

   #If #DIM eq "2"

    ! ------------------------------------------------------------------------------
    !    2D : 4th order (rectangular)
    ! ------------------------------------------------------------------------------


    if( addDissipation )then

      if( t.le.3*dt )then
        write(*,'("advBA: order=4: addDissipation=",l2," adc=",e10.2)') addDissipation,adc
      end if           
      sigma1=adc
      adxSosup(0)=sigma1/(2*64.)
      adxSosup(1)=sigma1/(2*64.)
      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
        do m=ex,hz
          un(i1,i2,i3,m)=u(i1,i2,i3,m) + sosupDiss2d6(u,i1,i2,i3,m)
        end do         
      endLoopsMask()
    end if
      
!    else if( materialType.eq.isotropic )then
!      updateMx2dOrder4()
!
!    else if( dispersionModel.eq.noDispersion )then
!      updateBA2dOrder4()
!    else
!      ! updateBAGDM(DIM,ORDER,GRIDTYPE)
!      if( useOpt )then
!        if( useSuperGrid.eq.0  ) then
!          updateBAGDMOpt(2,4,rectangular,,,)
!        else
!          updateBAGDMOpt(2,4,rectangular,etax(i1)*,etay(i2)*,etaz(i3)*)
!        end if
!      else 
!        updateBAGDM(2,4,rectangular)
!      end if 
!    end if 

   #Else

    ! ------------------------------------------------------------------------------
    !    3D : 4th order  (rectangular)
    ! ------------------------------------------------------------------------------

    if( addDissipation )then

      if( t.le.3*dt )then
        write(*,'("advBA: order=4: addDissipation=",l2," adc=",e10.2)') addDissipation,adc
      end if
            
      sigma1=adc
      adxSosup(0)=sigma1/(2*64.)
      adxSosup(1)=sigma1/(2*64.)
      adxSosup(2)=sigma1/(2*64.)
      beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
        do m=ex,hz
          un(i1,i2,i3,m)=u(i1,i2,i3,m) + sosupDiss3d6(u,i1,i2,i3,m)
        end do         
      endLoopsMask()
    end if
!    else if( materialType.eq.isotropic )then
!      ! updateMx3dOrder4()
!      stop 7676
!    else if( dispersionModel.eq.noDispersion )then
!      updateBA3dOrder4()
!    else
!      if( useOpt )then
!        if( useSuperGrid.eq.0  ) then
!          updateBAGDMOpt(3,4,rectangular,,,)
!        else
!          updateBAGDMOpt(3,4,rectangular,etax(i1)*,etay(i2)*,etaz(i3)*)
!        end if 
!      else 
!        updateBAGDM(3,4,rectangular)
!      end if 
!    end if 

   #End

      
  #Elif #ORDER eq "6" 

   stop 6666
   
  #Elif #ORDER eq "8"

  stop 8888

  #Else
     write(*,*) 'advBA:ERROR orderOfAccuracy,orderInTime=',orderOfAccuracy,orderInTime
    stop 1

  #End

  ! End if rectnagular 
#End

 else               

 #If #GRIDTYPE eq "curvilinear"

!       **********************************************
!       *************** curvilinear ******************
!       **********************************************


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


     buildFile(advBA2dOrder2r,2,2,rectangular)
     buildFile(advBA3dOrder2r,3,2,rectangular)

!      buildFile(advMx2dOrder2c,2,2,curvilinear)
!      buildFile(advMx3dOrder2c,3,2,curvilinear)

     buildFile(advBA2dOrder4r,2,4,rectangular)
     buildFile(advBA3dOrder4r,3,4,rectangular)

!      buildFile(advMx2dOrder4c,2,4,curvilinear)
!      buildFile(advMx3dOrder4c,3,4,curvilinear)

!       buildFile(advMx2dOrder6r,2,6,rectangular)
!       buildFile(advMx3dOrder6r,3,6,rectangular)
! 
!        ! build these for testing symmetric operators -- BC's not implemented yet
!       buildFile(advMx2dOrder6c,2,6,curvilinear)
!       buildFile(advMx3dOrder6c,3,6,curvilinear)
! 
!       buildFile(advMx2dOrder8r,2,8,rectangular)
!       buildFile(advMx3dOrder8r,3,8,rectangular)
! 
!        ! build these for testing symmetric operators -- BC's not implemented yet
!       buildFile(advMx2dOrder8c,2,8,curvilinear)
!       buildFile(advMx3dOrder8c,3,8,curvilinear)



! build an empty version of high order files so we do not have to compile the full version
#beginMacro ADV_MAXWELL_NULL(NAME,DIM,ORDER,GRIDTYPE)
 subroutine NAME(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                 mask,rsxy,  um,u,un,f,fa, K0i,matMask, pm,p,pn,xy, etax,etay,etaz, bc, dis, varDis, ipar, rpar, ierr )
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

!       buildFileNull(advMx2dOrder6r,2,6,rectangular)
!       buildFileNull(advMx3dOrder6r,3,6,rectangular)
! 
!       buildFileNull(advMx2dOrder6c,2,6,curvilinear)
!       buildFileNull(advMx3dOrder6c,3,6,curvilinear)
! 
!       buildFileNull(advMx2dOrder8r,2,8,rectangular)
!       buildFileNull(advMx3dOrder8r,3,8,rectangular)
! 
!       buildFileNull(advMx2dOrder8c,2,8,curvilinear)
!       buildFileNull(advMx3dOrder8c,3,8,curvilinear)



      subroutine advBA(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                       mask,rx,  um,u,un,f,fa, K0i, matMask, pm,p,pn,xy,etax,etay,etaz, bc, dis, varDis, ipar, rpar, ierr )
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

      ! Polarization vectors
      real pm(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real p(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
      real pn(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)

      real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)

      ! Super-grid layer functions 
      real etax(nd1a:nd1b), etay(nd2a:nd2b), etaz(nd3a:nd3b)

      ! real ut5(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      ! real ut6(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      ! real ut7(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real dis(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      real varDis(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      real rx(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)

      integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)

      real K0i(0:5,0:5,0:*)
      integer matMask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)

      integer bc(0:1,0:2),ierr

      integer ipar(0:*)
      real rpar(0:*)
      
!     ---- local variables -----
      integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderInTime
      integer addForcing,orderOfDissipation,option
      integer useSosupDissipation,sosupDissipationOption
      integer updateSolution,updateDissipation,computeUt
      integer useWhereMask,solveForE,solveForH,grid
      integer ex,ey,ez, hx,hy,hz

      integer rectangular,curvilinear
      parameter( rectangular=0, curvilinear=1 )
!...........end   statement functions

      !  write(*,*) 'Inside advBA...'

      orderOfAccuracy    =ipar(2)
      gridType           =ipar(1)
      useSosupDissipation=ipar(34)
      sosupDissipationOption=ipar(35)
      updateSolution        =ipar(36)
      updateDissipation     =ipar(37)
      computeUt             =ipar(38)

      if( orderOfAccuracy.eq.2 )then
 
        if( nd.eq.2 .and. gridType.eq.rectangular ) then
          ! BA version 
          call advBA2dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,rx, um,u,un,f,fa,K0i,matMask, pm,p,pn,xy,etax,etay,etaz,bc, dis,varDis, ipar, rpar, ierr )
        else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
          write(*,*) 'advBA -- unimplemented option'
          stop 1234
        else if( nd.eq.3 .and. gridType.eq.rectangular ) then
          ! BA version 
          call advBA3dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,rx, um,u,un,f,fa,K0i,matMask, pm,p,pn,xy,etax,etay,etaz,bc, dis,varDis, ipar, rpar, ierr )
        else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
          write(*,*) 'advBA -- unimplemented option'
          stop 1234
        else
          stop 2271
        end if
  
       else if( orderOfAccuracy.eq.4 ) then
  
        if( nd.eq.2 .and. gridType.eq.rectangular )then
          ! BA version 
          ! write(*,'(" call advBA2dOrder4r...")')
          call advBA2dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,rx, um,u,un,f,fa,K0i,matMask, pm,p,pn,xy,etax,etay,etaz,bc, dis,varDis, ipar, rpar, ierr )
        else if(nd.eq.2 .and. gridType.eq.curvilinear )then
          write(*,*) 'advBA -- unimplemented option'
          stop 1234
        else if(  nd.eq.3 .and. gridType.eq.rectangular )then
          ! BA version 
          call advBA3dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                              mask,rx, um,u,un,f,fa,K0i,matMask, pm,p,pn,xy,etax,etay,etaz,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
          write(*,*) 'advBA -- unimplemented option'
          stop 1234
  
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








