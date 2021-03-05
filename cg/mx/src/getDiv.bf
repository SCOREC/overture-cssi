!
! Optimized routine to compute the divergence 
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

  

! =========================================================================================
! Macro: Compute induced fields D and B for the BI-ANISTROPIC equations
!
!   POLAR -- polarization: TEZ or NONE
! ========================================================================================
#beginMacro computeInducedFields(POLAR)

  if( t.lt.2*dt )then
    write(*,'("computeInducedFields polar=POLAR... t=",e10.2)') t
  end if
  
  ! loop bounds -- include ghost points 
  mr=0  ! default value for a single material
  beginLoopsMask(i1,i2,i3,m1a,m1b,m2a,m2b,m3a,m3b)

    if( numberOfMaterialRegions.gt.1 )then
      mr = matMask(i1,i2,i3)
      if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
         stop 9999
      end if            
    end if       

    ! -- Compute [D,B] = K0*U + P 
    #If #POLAR eq "TEZ" 
      do m=0,2 
        v(i1,i2,i3,m)= K03(m,0,mr)*u(i1,i2,i3,0) + K03(m,1,mr)*u(i1,i2,i3,1) + K03(m,2,mr)*u(i1,i2,i3,2) 
      end do
    #Else
      do m=0,5 
        v(i1,i2,i3,m)= K0(m,0,mr)*u(i1,i2,i3,0) + K0(m,1,mr)*u(i1,i2,i3,1) + K0(m,2,mr)*u(i1,i2,i3,2) + \
                       K0(m,3,mr)*u(i1,i2,i3,3) + K0(m,4,mr)*u(i1,i2,i3,4) + K0(m,5,mr)*u(i1,i2,i3,5) 
      end do
    #End

      
  endLoopsMask()

#endMacro

! =========================================================================================
! Macro: Compute induced fields D and B for the BI-ANISTROPIC GDM equations
!
!   POLAR -- polarization: TEZ or NONE
! ========================================================================================
#beginMacro computeInducedFieldsGDM(POLAR)

  if( t.lt.2*dt )then
    write(*,'("computeInducedFieldsGDM polar=POLAR... t=",e10.2)') t
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
          pcIndex1(m,mr)=pc ! ** check me 
          pc=pc+2
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
  
  end do ! end do mr
  mr=0

  ! loop bounds -- include ghost points **finish me**
  ! include ghost 

  mr=0  ! default value for a single material
  beginLoopsMask(i1,i2,i3,m1a,m1b,m2a,m2b,m3a,m3b)

    if( numberOfMaterialRegions.gt.1 )then
      mr = matMask(i1,i2,i3)
      if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
         stop 9999
      end if            
    end if       

    ! ---- Compute p = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, pc )
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
       pc = pcIndex1(m,mr)
       ptSum(ec) = ptSum(ec) + p(i1,i2,i3,pc) 
    end do

    ! -- Compute [D,B] = K0*U + P 
    #If #POLAR eq "TEZ" 
      do m=0,2 
        v(i1,i2,i3,m)= K03(m,0,mr)*u(i1,i2,i3,0) + K03(m,1,mr)*u(i1,i2,i3,1) + K03(m,2,mr)*u(i1,i2,i3,2) + ptSum(m) 
      end do
    #Else
      do m=0,5 
        v(i1,i2,i3,m)= K0(m,0,mr)*u(i1,i2,i3,0) + K0(m,1,mr)*u(i1,i2,i3,1) + K0(m,2,mr)*u(i1,i2,i3,2) + \
                       K0(m,3,mr)*u(i1,i2,i3,3) + K0(m,4,mr)*u(i1,i2,i3,4) + K0(m,5,mr)*u(i1,i2,i3,5) + ptSum(m) 
      end do
    #End

      
  endLoopsMask()

#endMacro



! =========================================================================================
!
!   Compute the divergence of D and B for the BA MAXWELL
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================
#beginMacro computeDivBA(DIM,ORDER,GRIDTYPE,POLAR)

  if( t.lt.2*dt )then
    write(*,'("computeDiv (BA) dim=DIM order=ORDER grid=GRIDTYPE polar=POLAR... t=",e10.2)') t
  end if
  
  ! Loop bounds -- exclude super-grid layer ? 
  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

   #If #DIM eq "2"

     ! --- 2D -----

     #If #POLAR eq "TEZ" 
       ! ---- TEZ polarization Ex,Ey,Hz  -----
       #If #GRIDTYPE eq "rectangular"
         #If #ORDER eq "2" 
           d1x = vx22r(i1,i2,i3,ex)
           d2x = vx22r(i1,i2,i3,ey)
           d1y = vy22r(i1,i2,i3,ex)
           d2y = vy22r(i1,i2,i3,ey)

           b3x =vx22r(i1,i2,i3,hz)
           b3y =vy22r(i1,i2,i3,hz)
         #Else
           d1x = vx42r(i1,i2,i3,ex)
           d2x = vx42r(i1,i2,i3,ey)
           d1y = vy42r(i1,i2,i3,ex)
           d2y = vy42r(i1,i2,i3,ey)
         #End
       #Else
         #If #ORDER eq "2" 
           d1x = vx22(i1,i2,i3,ex)
           d2x = vx22(i1,i2,i3,ey)
           d1y = vy22(i1,i2,i3,ex)
           d2y = vy22(i1,i2,i3,ey)
         #End
           d1x = vx42(i1,i2,i3,ex)
           d2x = vx42(i1,i2,i3,ey)
           d1y = vy42(i1,i2,i3,ex)
           d2y = vy42(i1,i2,i3,ey)
         #Else
       #End
  
       dDiv = d1x + d2y
       if( saveDivergence.eq.1 )then
         divD(i1,i2,i3) = dDiv
       end if 
       divDmax=max(divDmax,abs(dDiv))
       gradDmax = max( gradDmax,abs(d1x),abs(d2x),abs(d1y),abs(d2y))
       gradBmax = max( gradBmax,abs(b3x),abs(b3y) )

     #Else
       ! -- solve for all fields 

       #If #GRIDTYPE eq "rectangular"
         #If #ORDER eq "2" 
           d1x = vx22r(i1,i2,i3,ex)
           d2x = vx22r(i1,i2,i3,ey)
           d3x = vx22r(i1,i2,i3,ez)
  
           d1y = vy22r(i1,i2,i3,ex)
           d2y = vy22r(i1,i2,i3,ey)
           d3y = vy22r(i1,i2,i3,ez)
  
           b1x = vx22r(i1,i2,i3,hx)
           b2x = vx22r(i1,i2,i3,hy)
           b3x = vx22r(i1,i2,i3,hz)
  
           b1y = vy22r(i1,i2,i3,hx)
           b2y = vy22r(i1,i2,i3,hy)
           b3y = vy22r(i1,i2,i3,hz)
  
         #Else
           d1x = vx42r(i1,i2,i3,ex)
           d2x = vx42r(i1,i2,i3,ey)
           d3x = vx42r(i1,i2,i3,ez)
  
           d1y = vy42r(i1,i2,i3,ex)
           d2y = vy42r(i1,i2,i3,ey)
           d3y = vy42r(i1,i2,i3,ez)
  
           b1x = vx42r(i1,i2,i3,hx)
           b2x = vx42r(i1,i2,i3,hy)
           b3x = vx42r(i1,i2,i3,hz)
  
           b1y = vy42r(i1,i2,i3,hx)
           b2y = vy42r(i1,i2,i3,hy)
           b3y = vy42r(i1,i2,i3,hz)
  
         #End
       #Else
         ! -- curvilinear grid --
         #If #ORDER eq "2" 
           d1x = vx22(i1,i2,i3,ex)
           d2x = vx22(i1,i2,i3,ey)
           d3x = vx22(i1,i2,i3,ez)
  
           d1y = vy22(i1,i2,i3,ex)
           d2y = vy22(i1,i2,i3,ey)
           d3y = vy22(i1,i2,i3,ez)
  
           b1x = vx22(i1,i2,i3,hx)
           b2x = vx22(i1,i2,i3,hy)
           b3x = vx22(i1,i2,i3,hz)
  
           b1y = vy22(i1,i2,i3,hx)
           b2y = vy22(i1,i2,i3,hy)
           b3y = vy22(i1,i2,i3,hz)
  
         #Else
           d1x = vx42(i1,i2,i3,ex)
           d2x = vx42(i1,i2,i3,ey)
           d3x = vx42(i1,i2,i3,ez)
  
           d1y = vy42(i1,i2,i3,ex)
           d2y = vy42(i1,i2,i3,ey)
           d3y = vy42(i1,i2,i3,ez)
  
           b1x = vx42(i1,i2,i3,hx)
           b2x = vx42(i1,i2,i3,hy)
           b3x = vx42(i1,i2,i3,hz)
  
           b1y = vy42(i1,i2,i3,hx)
           b2y = vy42(i1,i2,i3,hy)
           b3y = vy42(i1,i2,i3,hz)
         #End
       #End


       dDiv = d1x + d2y
       if( saveDivergence.eq.1 )then
         divD(i1,i2,i3) = dDiv
       end if 
       divDmax=max(divDmax,abs(dDiv))

       gradDmax = max( gradDmax,abs(d1x),abs(d2x),abs(d3x),abs(d1y),abs(d2y),abs(d3y))

       bDiv = b1x + b2y 
       if( saveDivergence.eq.1 )then
         divB(i1,i2,i3)=bDiv
       end if 
       divBmax = max(divBmax,abs(bDiv))
       gradBmax =  max( gradBmax,abs(b1x),abs(b2x),abs(b3x),abs(b1y),abs(b2y),abs(b3y))

     #End   

   #Else
     ! --- 3D -----

     #If #GRIDTYPE eq "rectangular"
       #If #ORDER eq "2" 
         d1x = vx23r(i1,i2,i3,ex)
         d2x = vx23r(i1,i2,i3,ey)
         d3x = vx23r(i1,i2,i3,ez)

         d1y = vy23r(i1,i2,i3,ex)
         d2y = vy23r(i1,i2,i3,ey)
         d3y = vy23r(i1,i2,i3,ez)

         d1z = vz23r(i1,i2,i3,ex)
         d2z = vz23r(i1,i2,i3,ey)
         d3z = vz23r(i1,i2,i3,ez)

         b1x = vx23r(i1,i2,i3,hx)
         b2x = vx23r(i1,i2,i3,hy)
         b3x = vx23r(i1,i2,i3,hz)

         b1y = vy23r(i1,i2,i3,hx)
         b2y = vy23r(i1,i2,i3,hy)
         b3y = vy23r(i1,i2,i3,hz)

         b1z = vz23r(i1,i2,i3,hx)
         b2z = vz23r(i1,i2,i3,hy)
         b3z = vz23r(i1,i2,i3,hz)
       #Else
         d1x = vx43r(i1,i2,i3,ex)
         d2x = vx43r(i1,i2,i3,ey)
         d3x = vx43r(i1,i2,i3,ez)

         d1y = vy43r(i1,i2,i3,ex)
         d2y = vy43r(i1,i2,i3,ey)
         d3y = vy43r(i1,i2,i3,ez)

         d1z = vz43r(i1,i2,i3,ex)
         d2z = vz43r(i1,i2,i3,ey)
         d3z = vz43r(i1,i2,i3,ez)

         b1x = vx43r(i1,i2,i3,hx)
         b2x = vx43r(i1,i2,i3,hy)
         b3x = vx43r(i1,i2,i3,hz)

         b1y = vy43r(i1,i2,i3,hx)
         b2y = vy43r(i1,i2,i3,hy)
         b3y = vy43r(i1,i2,i3,hz)

         b1z = vz43r(i1,i2,i3,hx)
         b2z = vz43r(i1,i2,i3,hy)
         b3z = vz43r(i1,i2,i3,hz)
       #End
     #Else
       ! -- curvilinear grid --
       #If #ORDER eq "2" 
         d1x = vx23(i1,i2,i3,ex)
         d2x = vx23(i1,i2,i3,ey)
         d3x = vx23(i1,i2,i3,ez)

         d1y = vy23(i1,i2,i3,ex)
         d2y = vy23(i1,i2,i3,ey)
         d3y = vy23(i1,i2,i3,ez)

         d1z = vz23(i1,i2,i3,ex)
         d2z = vz23(i1,i2,i3,ey)
         d3z = vz23(i1,i2,i3,ez)

         b1x = vx23(i1,i2,i3,hx)
         b2x = vx23(i1,i2,i3,hy)
         b3x = vx23(i1,i2,i3,hz)

         b1y = vy23(i1,i2,i3,hx)
         b2y = vy23(i1,i2,i3,hy)
         b3y = vy23(i1,i2,i3,hz)

         b1z = vz23(i1,i2,i3,hx)
         b2z = vz23(i1,i2,i3,hy)
         b3z = vz23(i1,i2,i3,hz)
       #Else
         d1x = vx43(i1,i2,i3,ex)
         d2x = vx43(i1,i2,i3,ey)
         d3x = vx43(i1,i2,i3,ez)

         d1y = vy43(i1,i2,i3,ex)
         d2y = vy43(i1,i2,i3,ey)
         d3y = vy43(i1,i2,i3,ez)

         d1z = vz43(i1,i2,i3,ex)
         d2z = vz43(i1,i2,i3,ey)
         d3z = vz43(i1,i2,i3,ez)

         b1x = vx43(i1,i2,i3,hx)
         b2x = vx43(i1,i2,i3,hy)
         b3x = vx43(i1,i2,i3,hz)

         b1y = vy43(i1,i2,i3,hx)
         b2y = vy43(i1,i2,i3,hy)
         b3y = vy43(i1,i2,i3,hz)

         b1z = vz43(i1,i2,i3,hx)
         b2z = vz43(i1,i2,i3,hy)
         b3z = vz43(i1,i2,i3,hz)
       #End
     #End

     dDiv = d1x + d2y + d3z 
     if( saveDivergence.eq.1 )then
       divD(i1,i2,i3) = dDiv
     end if 
     divDmax=max(divDmax,abs(dDiv))

     bDiv = b1x + b2y + b3z 
     if( saveDivergence.eq.1 )then
       divB(i1,i2,i3)=bDiv
     end if 
     divBmax = max(divBmax,abs(bDiv))

     gradDmax = max( gradDmax,abs(d1x),abs(d2x),abs(d3x),abs(d1y),abs(d2y),abs(d3y),abs(d1z),abs(d2z),abs(d3z))
     gradBmax = max( gradBmax,abs(b1x),abs(b2x),abs(b3x),abs(b1y),abs(b2y),abs(b3y),abs(b1z),abs(b2z),abs(b3z))

   #End
        
  endLoopsMask()



#endMacro


! =========================================================================================
!
!   Compute the divergence of E for the isotropic MAXWELL
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================
#beginMacro computeDiv(DIM,ORDER,GRIDTYPE)

  if( t.lt.2*dt )then
    write(*,'("computeDiv dim=DIM order=ORDER grid=GRIDTYPE ... t=",e10.2)') t
  end if
  
  ! FINISH ME -- divMax, gradmax ...
  stop 1234

  ! Loop bounds -- exclude absorbing layers ? 
  beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

   #If #DIM eq "2"
     ! --- 2D -----

     #If #ORDER eq "2" 

       #If #GRIDTYPE eq "rectangular"
         divD(i1,i2,i3) = vx22r(i1,i2,i3,ex) + vy22r(i1,i2,i3,ey)
       #Else
         divD(i1,i2,i3) = vx22(i1,i2,i3,ex) + vy22(i1,i2,i3,ey)
       #End

     #Else

       #If #GRIDTYPE eq "rectangular"
         divD(i1,i2,i3) = vx42r(i1,i2,i3,ex) + vy42r(i1,i2,i3,ey)
       #Else
         divD(i1,i2,i3) = vx42(i1,i2,i3,ex) + vy42(i1,i2,i3,ey)
       #End


     #End

   #Else
     ! --- 3D -----

    #If #ORDER eq "2" 

      #If #GRIDTYPE eq "rectangular"
        divD(i1,i2,i3) = vx23r(i1,i2,i3,ex) + vy23r(i1,i2,i3,ey) + vz23r(i1,i2,i3,ez)
      #Else
        divD(i1,i2,i3) = vx23(i1,i2,i3,ex) + vy23(i1,i2,i3,ey) + vz23(i1,i2,i3,ez)
      #End

    #Else

      #If #GRIDTYPE eq "rectangular"
        divD(i1,i2,i3) = vx43r(i1,i2,i3,ex) + vy43r(i1,i2,i3,ey) + vz43r(i1,i2,i3,ez)
      #Else
        divD(i1,i2,i3) = vx43(i1,i2,i3,ex) + vy43(i1,i2,i3,ey) + vz43(i1,i2,i3,ez)
      #End
    #End
 
   #End  

   divDmax=max(divDmax,divD(i1,i2,i3))

  endLoopsMask()



#endMacro



  
! **********************************************************************************
! Macro GET_DIVERGENCE
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
!  GRIDTYPE : rectangular, curvilinear
! **********************************************************************************
#beginMacro GET_DIVERGENCE(NAME,DIM,ORDER,GRIDTYPE)
 subroutine NAME(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
   mask,rsxy,  u,p,v, divD, divB, K0, matMask, \
   ipar, rpar, ierr )
!======================================================================
!   Compute the divergence 
!
!   v : tempoary storage to hold D and B for the BA Maxwell
!   divD(i1,i2,i3,0:1) : return div(E) ( or div(D) for BA) here 
!   divB(i1,i2,i3,0:1) : return div(H) ( or Div(B) for BA) here 
!======================================================================
 implicit none
 integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b

 real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
 real divD(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
 real divB(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)

 ! Polarization vectors
 real p(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
 real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)

 integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)

 real K0(0:5,0:5,0:*)  ! material matrix 
 integer matMask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)

 integer ierr

 integer ipar(0:*)
 real rpar(0:*)
      
 !     ---- local variables -----
 integer c,i1,i2,i3,n,gridType,orderOfAccuracy,debug
 integer m1a,m1b,m2a,m2b,m3a,m3b,numGhost,k1,k2 
 integer saveDivergence
 real t,dt
 real d1x,d2x,d3x, d1y,d2y,d3y, d1z,d2z,d3z
 real b1x,b2x,b3x, b1y,b2y,b3y, b1z,b2z,b3z
 real dDiv, bDiv 
 real divDmax,divBmax, gradDmax, gradBmax

 integer ex,ey,ez, hx,hy,hz, solveForAllFields, grid 

 integer maxRegions,NpMax
 parameter( maxRegions=100,NpMax=10 )  ! FIX ME 

 integer numberOfMaterialRegions, mr
 real K03(0:2,0:2,0:maxRegions) ! 3x3 material matrix for TEz polarization

 integer Np(6,6,0:maxRegions-1)
 real gdmPar(4,NpMax,6,6,0:maxRegions-1), ptSum(0:5) 
 integer ec,pc,qc, pct,qct

 integer numPolarizationTerms
 integer maxNumPolarizationTerms
 parameter( maxNumPolarizationTerms=200 )

 real dx(0:2),dr(0:2)

 integer rectangular,curvilinear
 parameter( rectangular=0, curvilinear=1 )

 !...........start statement function
 integer kd,m
 real rx,ry,rz,sx,sy,sz,tx,ty,tz

 declareDifferenceOrder2(v,RX)
 declareDifferenceOrder4(v,RX)

 ! dispersion
 integer dispersionModel

 integer numTerms1(0:maxRegions),ecIndex1(maxNumPolarizationTerms,0:maxRegions),pcIndex1(maxNumPolarizationTerms,0:maxRegions)

! Dispersion models
 #Include "dispersionModelsFortranInclude.h"

!  integer forcingOption   
!  ! forcing options
!  #Include "forcingDefineFortranInclude.h"

 integer method,nfdtd,bamx
 parameter( nfdtd=5, bamx=7 )

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
 defineDifferenceOrder2Components1(v,RX)
 defineDifferenceOrder4Components1(v,RX)


!...........end   statement functions


 dx(0) =rpar(0)
 dx(1) =rpar(1)
 dx(2) =rpar(2)
 dr(0) =rpar(3)
 dr(1) =rpar(4)
 dr(2) =rpar(5)
 t     =rpar(6)
 ! This are provided on input and output (max is taken over all grids)
 divDmax = rpar(7)  ! max(div(D)) returned here 
 gradDmax= rpar(8)  ! max(grad(D)) returned here 
 divBmax = rpar(9)  ! max(div(B)) returned here 
 gradBmax= rpar(10) ! max(grad(B)) returned here 

 saveDivergence     =ipar(0)
 method             =ipar(1)
 gridType           =ipar(2)
 orderOfAccuracy    =ipar(3)
 ex                 =ipar(4)
 ey                 =ipar(5)
 ez                 =ipar(6)
 hx                 =ipar(7)
 hy                 =ipar(8)
 hz                 =ipar(9)
 dispersionModel    =ipar(10)
 grid               =ipar(11)

 solveForAllFields  =ipar(12)
 numberOfMaterialRegions = ipar(13)
 debug              =ipar(14)
 

 if( t.le.2*dt .and. debug.gt.0 )then
   write(*,*) 'Inside NAME...'
   write(*,'("solveForAllFields=",i2)') solveForAllFields
   write(*,'("dispersionModel=",i2)') dispersionModel

 end if


 if( method.eq.bamx )then
   if( numberOfMaterialRegions.gt.maxRegions )then
     write(*,*) 'getDiv: Error: numberOfMaterialRegions=',numberOfMaterialRegions,' is bigger than maxRegions=',maxRegions
     write(*,*) 'FIX ME BILL!'
     stop 1002
   end if
   ! 3x3 Material matrix for TEz polarization
   ! We use ex=0,ey=1 and hz=5 entries in K0i(0:5,0:5) 
   do mr=0,numberOfMaterialRegions-1
     K03(0,0,mr) = K0(0,0,mr)
     K03(0,1,mr) = K0(0,1,mr)
     K03(0,2,mr) = K0(0,5,mr)
  
     K03(1,0,mr) = K0(1,0,mr)
     K03(1,1,mr) = K0(1,1,mr)
     K03(1,2,mr) = K0(1,5,mr)
  
     K03(2,0,mr) = K0(5,0,mr)
     K03(2,1,mr) = K0(5,1,mr)
     K03(2,2,mr) = K0(5,5,mr)
   end do
 end if

 if( method.eq.bamx .and. t.lt.dt )then
  
   write(*,*) 'numberOfMaterialRegions=',numberOfMaterialRegions
   do mr=0,numberOfMaterialRegions-1
     write(*,*) 'Material region=',mr 
     write(*,'("K0=",6("[",6(f6.3,1x),"]",/,4x))') ((K0(i1,i2,mr),i1=0,5),i2=0,5)
     if( solveForAllFields .eq. 0 )then
       write(*,'("K03 =",3("[",3(f6.3,1x),"]",/,4x))') ((K03(i1,i2,mr),i1=0,2),i2=0,2)
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

 if( method.eq.bamx .and. dispersionModel.ne.noDispersion )then
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
     if( debug.gt.1 )then
       write(*,'("--getDiv-- dispersionModel=",i4," numPolarizationTerms=",i6)') dispersionModel,numPolarizationTerms
     end if
     
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


 if( t.eq.0. .and. dispersionModel.ne.noDispersion .and. debug>1 )then
    write(*,'("--getDiv-- dispersionModel=",i4)') dispersionModel
 end if

  ! We need to evaluate D and B at extra points so that we can take the divergence
  numGhost = orderOfAccuracy/2

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

  ! divDmax=0.
  ! divBmax=0.
  ! gradDmax=0.
  ! gradBmax=0.
  if( method.ne.bamx )then
   ! --- Isotropic Maxwell ----

    computeDiv(DIM,ORDER,GRIDTYPE)

  else
    ! --- BA EQUATIONS ---
    if( dispersionModel.eq.noDispersion )then
      if( nd.eq.2 .and. solveForAllFields.eq.0 )then
        computeInducedFields(TEZ)
        computeDivBA(DIM,ORDER,GRIDTYPE,TEZ)
      else
        computeInducedFields(NONE)
        computeDivBA(DIM,ORDER,GRIDTYPE,NONE)
      end if 
    else
      if( nd.eq.2 .and. solveForAllFields.eq.0 )then
        computeInducedFieldsGDM(TEZ)
        computeDivBA(DIM,ORDER,GRIDTYPE,TEZ)
      else
        computeInducedFieldsGDM(NONE)
        computeDivBA(DIM,ORDER,GRIDTYPE,NONE)
      end if 
    end if 

  end if


  ! return max values in rpar
  rpar(7) =divDmax
  rpar(8) =gradDMax
  rpar(9) =divBmax
  rpar(10)=gradBMax

  return
  end 

#endMacro


 


#beginMacro buildFile(NAME,DIM,ORDER,GRIDTYPE)
#beginFile NAME.f
 GET_DIVERGENCE(NAME,DIM,ORDER,GRIDTYPE)
#endFile
#endMacro


     buildFile(getDiv2dOrder2r,2,2,rectangular)
     buildFile(getDiv3dOrder2r,3,2,rectangular)

     buildFile(getDiv2dOrder2c,2,2,curvilinear)
     buildFile(getDiv3dOrder2c,3,2,curvilinear)

     buildFile(getDiv2dOrder4r,2,4,rectangular)
     buildFile(getDiv3dOrder4r,3,4,rectangular)

     buildFile(getDiv2dOrder4c,2,4,curvilinear)
     buildFile(getDiv3dOrder4c,3,4,curvilinear)


      subroutine getDiv(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                        mask,rsxy,  u,p,v, divD,divB, K0, matMask,ipar, rpar, ierr )
!======================================================================
!   Advance a time step for Maxwells eqution
!     OPTIMIZED version for rectangular grids.
! nd : number of space dimensions
!
!======================================================================
      implicit none
      integer nd, n1a,n1b,n2a,n2b,n3a,n3b,
     & nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b

      real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real p(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)

      real divD(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      real divB(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)

      integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      integer matMask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)

      real K0(0:5,0:5,0:*)

      integer ierr

      integer ipar(0:*)
      real rpar(0:*)
      
!     ---- local variables -----
      integer saveDivergence,method,gridType,orderOfAccuracy
 
      integer nfdtd,bamx
      parameter( nfdtd=5, bamx=7 )

      integer rectangular,curvilinear
      parameter( rectangular=0, curvilinear=1 )
!...........end   statement functions

      ! write(*,*) 'Inside getDiv...'

      saveDivergence     =ipar(0)
      method             =ipar(1)
      gridType           =ipar(2)
      orderOfAccuracy    =ipar(3)

      if( orderOfAccuracy.eq.2 )then
 
        if( nd.eq.2 .and. gridType.eq.rectangular ) then
          call getDiv2dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                        mask,rsxy,  u,p,v, divD,divB, K0, matMask,ipar, rpar, ierr )
        else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
          write(*,*) 'getDiv -- unimplemented option'
          stop 1234
        else if( nd.eq.3 .and. gridType.eq.rectangular ) then
          call getDiv3dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                        mask,rsxy,  u,p,v, divD,divB, K0, matMask,ipar, rpar, ierr )
        else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
          write(*,*) 'getDiv -- unimplemented option'
          stop 1234
        else
          stop 2271
        end if
  
       else if( orderOfAccuracy.eq.4 ) then
  
        if( nd.eq.2 .and. gridType.eq.rectangular )then
          ! write(*,'(" call getDiv2dOrder4r...")')
          call getDiv2dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                        mask,rsxy,  u,p,v, divD,divB, K0, matMask,ipar, rpar, ierr )
        else if(nd.eq.2 .and. gridType.eq.curvilinear )then
          write(*,*) 'getDiv -- unimplemented option'
          stop 1234
        else if(  nd.eq.3 .and. gridType.eq.rectangular )then
          call getDiv3dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,\
                        mask,rsxy,  u,p,v, divD,divB, K0, matMask,ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
          write(*,*) 'getDiv -- unimplemented option'
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








