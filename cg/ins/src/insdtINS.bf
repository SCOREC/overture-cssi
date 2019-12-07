! ===========================================================================
!   Incompressible Navier-Stokes : explicit discretization 
! ===========================================================================


#Include "insdt.h"

#Include "advectionMacros.h"

!====================================================================
! This macro will build the statements that form the body of the loop
!
! IMPEXP: EXPLICIT, EXPLICIT_ONLY, BOTH
! SCALAR: NONE
!         PASSIVE - include equations for a passive scalar
! AXISYMMETRIC : YES or NO
! UPWIND : CENTERED, UPWIND or BWENO
!====================================================================
#beginMacro buildEquationsUpwind(IMPEXP,SCALAR,DIM,ORDER,GRIDTYPE,AXISYMMETRIC,UPWIND)


#If #IMPEXP == "EXPLICIT"

 ! INS, no AD

 #If #UPWIND == "CENTERED"
   ! --- centered approximations for advection ---
   ! .. this section is currently not used.
   #If #DIM == "2"
      ut(i1,i2,i3,uc)= -UU(uc)*UX(uc)-UU(vc)*UY(uc)-UX(pc)+nu*ULAP(uc)
      ut(i1,i2,i3,vc)= -UU(uc)*UX(vc)-UU(vc)*UY(vc)-UY(pc)+nu*ULAP(vc)
   #Else
    ut(i1,i2,i3,uc)= -UU(uc)*UX(uc)-UU(vc)*UY(uc)-UU(wc)*UZ(uc)-UX(pc)+nu*ULAP(uc)
    ut(i1,i2,i3,vc)= -UU(uc)*UX(vc)-UU(vc)*UY(vc)-UU(wc)*UZ(vc)-UY(pc)+nu*ULAP(vc)
    ut(i1,i2,i3,wc)= -UU(uc)*UX(wc)-UU(vc)*UY(wc)-UU(wc)*UZ(wc)-UZ(pc)+nu*ULAP(wc)
   #End
 #Else
   ! --- upwind approximations ---
   getAdvection(u,uu,i1,i2,i3,SCALAR,DIM,ORDER,GRIDTYPE,UPWIND, agu) 
  
   #If #DIM == "2"
      ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-UX(pc)+nu*ULAP(uc) 
      ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-UY(pc)+nu*ULAP(vc)
!       u(i1,i2,i3,3)= agu(uc,3) ! au in r direction 
!       u(i1,i2,i3,4)= agu(uc,8) ! au in s direction
!       u(i1,i2,i3,5)= agu(uc,4)  ! wpl uur
!       u(i1,i2,i3,6)= agu(uc,5)  ! wpr uur
!       u(i1,i2,i3,7)= agu(uc,6)  ! wml uur
!       u(i1,i2,i3,8)= agu(uc,7)  ! wmr uur
!       u(i1,i2,i3,9)= agu(uc,9)  ! wpl uus
!       u(i1,i2,i3,10)= agu(uc,10)  ! wpr uus
!       u(i1,i2,i3,11)= agu(uc,11)  ! wml uus
!       u(i1,i2,i3,12)= agu(uc,12)  ! wmr uus
! !
!       u(i1,i2,i3,13)= agu(vc,4)  ! wpl uvr
!       u(i1,i2,i3,14)= agu(vc,5)  ! wpr uvr
!       u(i1,i2,i3,15)= agu(vc,6)  ! wml uvr
!       u(i1,i2,i3,16)= agu(vc,7)  ! wmr uvr
!       u(i1,i2,i3,17)= agu(vc,9)  ! wpl uvs
!       u(i1,i2,i3,18)= agu(vc,10)  ! wpr uvs
!       u(i1,i2,i3,19)= agu(vc,11)  ! wml uvs
!       u(i1,i2,i3,20)= agu(vc,12)  ! wmr uvs
!       u(i1,i2,i3,21)= agu(uc,13)  ! au at where fix is turned on 
!       u(i1,i2,i3,22)= agu(vc,13)  ! au at where fix is turned on 

   #Else
      ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc)+agu(wc,uc))-UX(pc)+nu*ULAP(uc)
      ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc)+agu(wc,vc))-UY(pc)+nu*ULAP(vc)
      ut(i1,i2,i3,wc)= -(agu(uc,wc)+agu(vc,wc)+agu(wc,wc))-UZ(pc)+nu*ULAP(wc)
    ! ut(i1,i2,i3,uc)= -UU(uc)*UX(uc)-UU(vc)*UY(uc)-UU(wc)*UZ(uc)-UX(pc)+nu*ULAP(uc)
    ! ut(i1,i2,i3,vc)= -UU(uc)*UX(vc)-UU(vc)*UY(vc)-UU(wc)*UZ(vc)-UY(pc)+nu*ULAP(vc)
    ! ut(i1,i2,i3,wc)= -UU(uc)*UX(wc)-UU(vc)*UY(wc)-UU(wc)*UZ(wc)-UZ(pc)+nu*ULAP(wc)
   #End
 #End
 
  #If #AXISYMMETRIC == "YES"
   ! -- add on axisymmetric corrections ---
   ri=radiusInverse(i1,i2,i3)
   if( ri.ne.0. )then
     ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( UY(uc)*ri )
     ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( (UY(vc)-UU(vc)*ri)*ri )
   else
     ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( UYY(uc) )
     ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( .5*UYY(vc) )
   end if
  #End
 
 #Elif #IMPEXP == "EXPLICIT_ONLY" || #IMPEXP == "BOTH"
  ! explicit terms only, no diffusion
  #If #UPWIND == "CENTERED"
    ! centered approximation
    #If #DIM == "2"
      ut(i1,i2,i3,uc)= -UU(uc)*UX(uc)-UU(vc)*UY(uc)-UX(pc)
      ut(i1,i2,i3,vc)= -UU(uc)*UX(vc)-UU(vc)*UY(vc)-UY(pc)
    #Else
      ut(i1,i2,i3,uc)= -UU(uc)*UX(uc)-UU(vc)*UY(uc)-UU(wc)*UZ(uc)-UX(pc)
      ut(i1,i2,i3,vc)= -UU(uc)*UX(vc)-UU(vc)*UY(vc)-UU(wc)*UZ(vc)-UY(pc)
      ut(i1,i2,i3,wc)= -UU(uc)*UX(wc)-UU(vc)*UY(wc)-UU(wc)*UZ(wc)-UZ(pc)
    #End

  #Else
    ! upwind approximation/bweno
    getAdvection(u,uu,i1,i2,i3,SCALAR,DIM,ORDER,GRIDTYPE,UPWIND, agu) 
    #If #DIM == "2"
      ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc))-UX(pc)
      ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc))-UY(pc)
    #Else
      ut(i1,i2,i3,uc)= -(agu(uc,uc)+agu(vc,uc)+agu(wc,uc))-UX(pc)
      ut(i1,i2,i3,vc)= -(agu(uc,vc)+agu(vc,vc)+agu(wc,vc))-UY(pc)
      ut(i1,i2,i3,wc)= -(agu(uc,wc)+agu(vc,wc)+agu(wc,wc))-UZ(pc)
    #End
  #End

 #Else
   stop 688
 
 #End
 
 #If #IMPEXP == "BOTH"
  ! include implicit terms - diffusion
  #If #DIM == "2"
   uti(i1,i2,i3,uc)= nu*ULAP(uc)
   uti(i1,i2,i3,vc)= nu*ULAP(vc)
  #Elif #DIM == "3"
   uti(i1,i2,i3,uc)= nu*ULAP(uc)
   uti(i1,i2,i3,vc)= nu*ULAP(vc)
   uti(i1,i2,i3,wc)= nu*ULAP(wc)
  #Else
    stop 11
  #End
 
  #If #AXISYMMETRIC == "YES"
   ri=radiusInverse(i1,i2,i3)
   if( ri.ne.0. )then
     uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( UY(uc)*ri )
     uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( (UY(vc)-UU(vc)*ri)*ri )
   else
     uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( UYY(uc) )
     uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( .5*UYY(vc) )
   end if
  #End
 
 #End
 #endMacro 
 
 !====================================================================
 ! *OLD*
 ! This macro will build the statements that form the body of the loop
 !
 ! IMPEXP: EXPLICIT, EXPLICIT_ONLY, BOTH
 ! SCALAR: NONE
 !         PASSIVE - include equations for a passive scalar
 ! AXISYMMETRIC : YES or NO
 !====================================================================
 #beginMacro buildEquations(IMPEXP,SCALAR,DIM,ORDER,GRIDTYPE,AXISYMMETRIC)
 
 #If #IMPEXP == "EXPLICIT"
 
  ! INS, no AD
  #If #DIM == "2"
   ut(i1,i2,i3,uc)= -UU(uc)*UX(uc)-UU(vc)*UY(uc)-UX(pc)+nu*ULAP(uc)
   ut(i1,i2,i3,vc)= -UU(uc)*UX(vc)-UU(vc)*UY(vc)-UY(pc)+nu*ULAP(vc)
  #Else
   ut(i1,i2,i3,uc)= -UU(uc)*UX(uc)-UU(vc)*UY(uc)-UU(wc)*UZ(uc)-UX(pc)+nu*ULAP(uc)
   ut(i1,i2,i3,vc)= -UU(uc)*UX(vc)-UU(vc)*UY(vc)-UU(wc)*UZ(vc)-UY(pc)+nu*ULAP(vc)
   ut(i1,i2,i3,wc)= -UU(uc)*UX(wc)-UU(vc)*UY(wc)-UU(wc)*UZ(wc)-UZ(pc)+nu*ULAP(wc)
  #End
 
  #If #AXISYMMETRIC == "YES"
   ! -- add on axisymmetric corrections ---
   ri=radiusInverse(i1,i2,i3)
   if( ri.ne.0. )then
     ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( UY(uc)*ri )
     ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( (UY(vc)-UU(vc)*ri)*ri )
   else
     ut(i1,i2,i3,uc)=ut(i1,i2,i3,uc)+nu*( UYY(uc) )
     ut(i1,i2,i3,vc)=ut(i1,i2,i3,vc)+nu*( .5*UYY(vc) )
   end if
  #End
 
 #Elif #IMPEXP == "EXPLICIT_ONLY" || #IMPEXP == "BOTH"
  ! explicit terms only, no diffusion
  #If #DIM == "2"
   ut(i1,i2,i3,uc)= -UU(uc)*UX(uc)-UU(vc)*UY(uc)-UX(pc)
   ut(i1,i2,i3,vc)= -UU(uc)*UX(vc)-UU(vc)*UY(vc)-UY(pc)
  #Else
   ut(i1,i2,i3,uc)= -UU(uc)*UX(uc)-UU(vc)*UY(uc)-UU(wc)*UZ(uc)-UX(pc)
   ut(i1,i2,i3,vc)= -UU(uc)*UX(vc)-UU(vc)*UY(vc)-UU(wc)*UZ(vc)-UY(pc)
   ut(i1,i2,i3,wc)= -UU(uc)*UX(wc)-UU(vc)*UY(wc)-UU(wc)*UZ(wc)-UZ(pc)
  #End
 #Else
   stop 688
 
 #End
 
 #If #IMPEXP == "BOTH"
  ! include implicit terms - diffusion
  #If #DIM == "2"
   uti(i1,i2,i3,uc)= nu*ULAP(uc)
   uti(i1,i2,i3,vc)= nu*ULAP(vc)
  #Elif #DIM == "3"
   uti(i1,i2,i3,uc)= nu*ULAP(uc)
   uti(i1,i2,i3,vc)= nu*ULAP(vc)
   uti(i1,i2,i3,wc)= nu*ULAP(wc)
  #Else
    stop 11
  #End
 
  #If #AXISYMMETRIC == "YES"
   ri=radiusInverse(i1,i2,i3)
   if( ri.ne.0. )then
     uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( UY(uc)*ri )
     uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( (UY(vc)-UU(vc)*ri)*ri )
   else
     uti(i1,i2,i3,uc)=uti(i1,i2,i3,uc)+nu*( UYY(uc) )
     uti(i1,i2,i3,vc)=uti(i1,i2,i3,vc)+nu*( .5*UYY(vc) )
   end if
  #End
 
 #End
 
 #endMacro 
 
 
       buildFile(insdtINS2dOrder2,2,2)
       buildFile(insdtINS2dOrder4,2,4)
       buildFile(insdtINS3dOrder2,3,2)
       buildFile(insdtINS3dOrder4,3,4)
 
