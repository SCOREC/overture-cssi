! -*-f90-*-
! ----------------------------------------------------------------------------------
!  Macro:
!    --- Compute stencil coefficients by delta function approach ----
!
! The initial comptability equations take the form
!      AS * UGhost - AS1 * U1 - AS2 * U2 - ASF * F  = 0   
!
! The stencil is found by multiplying by AS^(-1)
!     UGhost = AS^(-1) AS1 * U1 + AS^(-1) AS2 * U2 + AS^(-1) * ASF = 0   
!
!  which we redefine AS1, AS2 and ASF to be 
!     UGhost = AS1 * U1 + AS2 * U2 + ASF * F  = 0   
! Input:
!   am : matrix of coefficients is stored here 
!   evalInterfaceEquations : macro that evaluates the interface equations
! ----------------------------------------------------------------------------------
#beginMacro saveStencilCoefficients(i1,i2,i3, j1,j2,j3, numberOfEquations,maxStencilSize,am,as,as1,as2,evalInterfaceEquations )

  ! hw1 = half stencil width
  hw1=orderOfAccuracy/2
  hw2=hw1
  if( nd.eq.2 )then
    hw3=0
  else
    hw3=hw1
  end if

  ! write(*,'("SAVE-STENCIL-COEFF: i1,i2,i3=",3i3," hw1,hw2,hw3=",3i2)') i1,i2,i3,hw1,hw2,hw3

  ! First eval equations with no perturbation --> save in f0 
  evalInterfaceEquations()
  do n1=0,numberOfEquations-1
   f0(n1)=f(n1)
  end do

  maxStencilDiff=0.

  ! ------- STEP 1: Compute the matrix that multiplies the unknown ghost points -------

  delta=1.  ! perturb E by this amount 
  do n2=0,numberOfEquations-1

    ! pertub one component: 
    perturbComponent(n2,delta)

    evalInterfaceEquations()
    
    ! compute the difference
    do n1=0,numberOfEquations-1
     f(n1)=f(n1)-f0(n1)
     as(n1,n2)=f(n1)
     maxStencilDiff=max(maxStencilDiff,am(n1,n2)-as(n1,n2))
    end do

    write(*,'("SAVE-STENCIL-COEFF: i1,i2,i3=",3i3," maxDiff=",e12.4)') i1,i2,i3,maxStencilDiff

    ! reset pertubation
    perturbComponent(n2,-delta)

  end do 
  coeffDiff = max(coeffDiff,maxStencilDiff)


  ! ------- STEP 2: Compute the matrices that multiply the boundary and interior points  -------

  ! ------ "LEFT" SIDE OF INTERFACE : u1 ---------
  if( is1.ne.0 )then
    i1a=i1
    i1b=i1+hw1*is1
    i1c=is1
  else
    i1a=i1-hw1
    i1b=i1+hw1
    i1c=1
  end if
  if( is2.ne.0 )then
    i2a=i2
    i2b=i2+hw2*is2
    i2c=is2
  else
    i2a=i2-hw2
    i2b=i2+hw2
    i2c=1
  end if
  if( is3.ne.0 )then
    i3a=i3
    i3b=i3+hw3*is3
    i3c=is3
  else
    i3a=i3-hw3
    i3b=i3+hw3
    i3c=1
  end if


  numCol=-1
  do i3p=i3a,i3b,i3c
  do i2p=i2a,i2b,i2c
  do i1p=i1a,i1b,i1c
  do n2=0,nd-1
    numCol=numCol+1 ! column in matrix 

    ! pertub one component: 
    u1(i1p,i2p,i3p,ex+n2)=u1(i1p,i2p,i3p,ex+n2)+(delta)

    evalInterfaceEquations()
    
    ! compute the difference
    do n1=0,numberOfEquations-1
      f(n1)=f(n1)-f0(n1)
      ! save the matrix (rectangular) (NOTE: Save minus)
      as1(n1,numCol) = - f(n1)
    end do

    ! reset pertubation
    u1(i1p,i2p,i3p,ex+n2)=u1(i1p,i2p,i3p,ex+n2)-(delta)

  end do 
  end do 
  end do 
  end do 

  if( numCol+1 .ne. maxStencilSize )then
    write(*,'("saveStencil:ERROR: numCol+1=",i3," is NOT equal to maxStencilSize=",i3)') numCol+1,maxStencilSize
    stop 7222
  end if
  
  write(*,'("as1: multiplies interior points on left (u1)")') 
  do n1=0,numberOfEquations-1
    write(*,'("as1(",i2,",:)=[",200(e12.2,1x),"]")') n1,(as1(n1,n2),n2=0,numCol)
  end do 

  ! -------- "RIGHT" SIDE OF THE INTERFACE : u2 --------
  if( js1.ne.0 )then
    j1a=j1
    j1b=j1+hw1*js1
    j1c=js1
  else
    j1a=j1-hw1
    j1b=j1+hw1
    j1c=1
  end if
  if( js2.ne.0 )then
    j2a=j2
    j2b=j2+hw2*js2
    j2c=js2
  else
    j2a=j2-hw2
    j2b=j2+hw2
    j2c=1
  end if
  if( js3.ne.0 )then
    j3a=j3
    j3b=j3+hw3*js3
    j3c=js3
  else
    j3a=j3-hw3
    j3b=j3+hw3
    j3c=1
  end if

  numCol=-1
  do j3p=j3a,j3b,j3c
  do j2p=j2a,j2b,j2c
  do j1p=j1a,j1b,j1c
  do n2=0,nd-1
    numCol=numCol+1 ! column in matrix 

    ! pertub one component: 
    u2(j1p,j2p,j3p,ex+n2)=u2(j1p,j2p,j3p,ex+n2)+(delta)

    evalInterfaceEquations()
    
    ! compute the difference
    do n1=0,numberOfEquations-1
     f(n1)=f(n1)-f0(n1)
     ! save the matrix (rectangular) (NOTE: save minus) 
     as2(n1,numCol) = - f(n1)
    end do

    ! reset pertubation
    u2(j1p,j2p,j3p,ex+n2)=u2(j1p,j2p,j3p,ex+n2)-(delta)

  end do 
  end do 
  end do 
  end do 

  write(*,'("as2: multiplies interior points on right (u2)")') 
  do n1=0,numberOfEquations-1
    write(*,'("as2(",i2,",:)=[",200(e12.2,1x),"]")') n1,(as2(n1,n2),n2=0,numCol)
  end do 

  ! restore 
  evalInterfaceEquations()

#endMacro

! ----------------------------------------------------------------------------------
!  Macro:
!    --- check the stencil coefficients after stage I
! Input:
! ----------------------------------------------------------------------------------
#beginMacro checkStencilCoefficients(i1,i2,i3, j1,j2,j3, numberOfEquations,as,as1,as2 )

  ! -- first compute the forcing ---
  do n1=0,numberOfEquations
   res(n1)=0.
  end do 

  ! ** FIX ME IN GENERAL **
  if( nd.ne.2 .or. orderOfAccuracy.ne.2 )then
    stop 8888
  end if

  if( twilightZone.eq.1 )then
    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )
 
    ueLap = uexx + ueyy
    veLap = vexx + veyy
    res(3) = - ( tau1*ueLap +tau2*veLap )*(1./eps1-1./eps2)
  end if

  write(*,'(" checkStencilCoeff: Stage 1: res=",10(e12.2,1x))') (res(n1),n1=0,numberOfEquations-1)

  ! add contributions from interior pts of u1 
  numCol=-1
  do i3p=i3a,i3b,i3c
  do i2p=i2a,i2b,i2c
  do i1p=i1a,i1b,i1c
  do n2=0,nd-1
    numCol=numCol+1 ! column in matrix 

    do n1=0,numberOfEquations-1
      res(n1) = res(n1) - as1(n1,numCol)*u1(i1p,i2p,i3p,ex+n2)
    end do

  end do 
  end do 
  end do 
  end do 

  write(*,'(" checkStencilCoeff: Stage 2: res=",10(e12.2,1x))') (res(n1),n1=0,numberOfEquations-1)

  ! add contributions from interior pts of u2 
  numCol=-1
  do j3p=j3a,j3b,j3c
  do j2p=j2a,j2b,j2c
  do j1p=j1a,j1b,j1c
  do n2=0,nd-1
    numCol=numCol+1 ! column in matrix 

    do n1=0,numberOfEquations-1
      res(n1) = res(n1) - as2(n1,numCol)*u2(j1p,j2p,j3p,ex+n2)
    end do

  end do 
  end do 
  end do 
  end do 

  write(*,'(" checkStencilCoeff: Stage 3: res=",10(e12.2,1x))') (res(n1),n1=0,numberOfEquations-1)

  ! This is for order=2 only !
  do n1=0,numberOfEquations-1
   res(n1) = res(n1) + as(n1,0)*u1(i1-is1,i2-is2,i3,ex) + as(n1,1)*u1(i1-is1,i2-is2,i3,ey)
   res(n1) = res(n1) + as(n1,2)*u2(j1-js1,j2-js2,j3,ex) + as(n1,3)*u2(j1-js1,j2-js2,j3,ey)
  end do


  do n1=0,numberOfEquations-1
    maxStencilRes = max(maxStencilRes,res(n1))
  end do
  write(*,'(" checkStencil: res=",50(e12.2,1x))') (res(n1),n1=0,numberOfEquations-1)

#endMacro


! ----------------------------------------------------------------------------------
!  Macro:
!    --- Compute actual stencil coefficients ----
!
! The initial comptability equations take the form
!      AS * UGhost - AS1 * U1 - AS2 * U2 - ASF * F  = 0   
!
! The stencil is found by multiplying by AS^(-1)
!     UGhost = AS^(-1) AS1 * U1 + AS^(-1) AS2 * U2 + AS^(-1) * ASF 
!
!  which we redefine AS1, AS2 and ASF to be 
!     UGhost = AS1 * U1 + AS2 * U2 + ASF * F 
! 
! Input:
! ----------------------------------------------------------------------------------
#beginMacro computeStencilCoefficients( numberOfEquations,maxStencilSize,as,as1,as2,asf )

  ! factor matrix that multiplies ghost point unknowns 
  call dgetrf( numberOfEquations,numberOfEquations,as(0,0),numberOfEquations,ipvt(0),info )
  if( info .ne.0 )then
    write(*,'("computeStencil: error return from dgetrf")')
    stop 9921
  end if
  
  nrhs=maxStencilSize
  call dgetrs( 'N', numberOfEquations, nrhs, as(0,0), numberOfEquations, ipvt(0), as1(0,0), numberOfEquations, info )
  if( info .ne.0 )then
    write(*,'("computeStencil: error return from dgetrs")')
    stop 9921
  end if

  call dgetrs( 'N', numberOfEquations, nrhs, as(0,0), numberOfEquations, ipvt(0), as2(0,0), numberOfEquations, info )
  if( info .ne.0 )then
    write(*,'("computeStencil: error return from dgetrs")')
    stop 9921
  end if


  ! initialize asf = - identity
  do n2=0,numberOfEquations-1
    do n1=0,numberOfEquations-1
      asf(n1,n2)=0.
    end do
    asf(n2,n2)= -1.
  end do
  call dgetrs( 'N', numberOfEquations, numberOfEquations, as(0,0), numberOfEquations, ipvt(0), asf(0,0), numberOfEquations, info )
  if( info .ne.0 )then
    write(*,'("computeStencil: error return from dgetrs")')
    stop 9921
  end if

  write(*,'("STAGE II: as1: multiplies interior points on left (u1)")') 
  do n1=0,numberOfEquations-1
    write(*,'("as1(",i2,",:)=[",200(e12.2,1x),"]")') n1,(as1(n1,n2),n2=0,nrhs-1)
  end do 

  write(*,'("STAGE II: as2: multiplies interior points on right (u2)")') 
  do n1=0,numberOfEquations-1
    write(*,'("as2(",i2,",:)=[",200(e12.2,1x),"]")') n1,(as2(n1,n2),n2=0,nrhs-1)
  end do 
    
  write(*,'("STAGE II: asf: multiplies RHS f ")') 
  do n1=0,numberOfEquations-1
    write(*,'("asf(",i2,",:)=[",200(e12.2,1x),"]")') n1,(asf(n1,n2),n2=0,numberOfEquations-1)
  end do 
    
#endMacro 



! ----------------------------------------------------------------------------------
!  Macro:
!    ---  Evaluate the stencil and check the error ---
! Input:
! ----------------------------------------------------------------------------------
#beginMacro checkStencil(i1,i2,i3, j1,j2,j3, numberOfEquations,aSave,as,as1,as2,asf,ue )

  ! -- first compute the forcing ---
  do n1=0,numberOfEquations-1
   f0(n1)=0.   ! holds the forcing 
   res(n1)=0.  ! holds solution
  end do 

  ! ** FIX ME IN GENERAL **
  if( nd.ne.2 .or. orderOfAccuracy.ne.2 )then
    stop 8888
  end if

  if( twilightZone.eq.1 )then
    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, uexx )
    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ex, ueyy )
    call ogderiv(ep, 0,2,0,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, vexx )
    call ogderiv(ep, 0,0,2,0, xy1(i1,i2,i3,0),xy1(i1,i2,i3,1),0.,t, ey, veyy )
 
    ueLap = uexx + ueyy
    veLap = vexx + veyy
    f0(3) = ( tau1*ueLap +tau2*veLap )*(1./eps1-1./eps2)
  end if



  ! Compute contributions from the forcing: ASF * F 
  do n1=0,numberOfEquations-1
    do n2=0,numberOfEquations-1
     res(n1) = res(n1) + asf(n1,n2)*f0(n2);
    end do
  end do

  write(*,'(" checkStencil: Stage 1: f0=",10(e12.2,1x))') (f0(n1),n1=0,numberOfEquations-1)
  write(*,'(" checkStencil: Stage 1: res=",10(e12.2,1x))') (res(n1),n1=0,numberOfEquations-1)


  ! add contributions from interior pts of u1 
  numCol=-1
  do i3p=i3a,i3b,i3c
  do i2p=i2a,i2b,i2c
  do i1p=i1a,i1b,i1c
  do n2=0,nd-1
    numCol=numCol+1 ! column in matrix 

    do n1=0,numberOfEquations-1
      res(n1) = res(n1) + as1(n1,numCol)*u1(i1p,i2p,i3p,ex+n2)
    end do

  end do 
  end do 
  end do 
  end do 

  write(*,'(" checkStencil: Stage 2: res=",10(e12.2,1x))') (res(n1),n1=0,numberOfEquations-1)

  ! compare to stage 1 res
  do n1=0,numberOfEquations-1
    f0(n1)=0.
    do n2=0,numberOfEquations-1
     f0(n1) = f0(n1) + a4Save(n1,n2)*res(n2);
    end do
  end do
  write(*,'(" checkStencil: Stage 2: aSave*res=",10(e12.2,1x))') (f0(n1),n1=0,numberOfEquations-1)


  ! add contributions from interior pts of u2 
  numCol=-1
  do j3p=j3a,j3b,j3c
  do j2p=j2a,j2b,j2c
  do j1p=j1a,j1b,j1c
  do n2=0,nd-1
    numCol=numCol+1 ! column in matrix 

    do n1=0,numberOfEquations-1
      res(n1) = res(n1) + as2(n1,numCol)*u2(j1p,j2p,j3p,ex+n2)
    end do

  end do 
  end do 
  end do 
  end do 

  write(*,'(" checkStencil: Stage 3: res=",10(e12.2,1x))') (res(n1),n1=0,numberOfEquations-1)

  ! compare to stage 1 res
  do n1=0,numberOfEquations-1
    f0(n1)=0.
    do n2=0,numberOfEquations-1
     f0(n1) = f0(n1) + a4Save(n1,n2)*res(n2);
    end do
  end do
  write(*,'(" checkStencil: Stage 3: aSave*res=",10(e12.2,1x))') (f0(n1),n1=0,numberOfEquations-1)

  write(*,'(" checkStencil: [true,u-stencil]=",10("[",e12.2,",",e12.2,"]",2x))') (ue(n1),res(n1),n1=0,numberOfEquations-1)

  do n1=0,numberOfEquations-1
    ! ue(:) holds solution for u on ghost 
    res(n1) = res(n1) - ue(n1) 
    maxStencilErr = max(maxStencilErr,res(n1));
  end do
  write(*,'(" checkStencil: Err in u =",50(e12.2,1x))') (res(n1),n1=0,numberOfEquations-1)

  if( .true. )then
    write(*,*) 'Stop here for now'
    stop 1234
  end if
  

#endMacro
  
