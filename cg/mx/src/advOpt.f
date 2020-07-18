! This file automatically generated from advOpt.bf with bpp.
!
! Optimized advance routines for cgmx
!

! These next include files will define the macros that will define the difference approximations
! The actual macro is called below
c Use this next macro to declare the statement functions that are defined below
c To include derivatives of rx use OPTION=RX


c Define statement functions for difference approximations of order 2 
c To include derivatives of rx use OPTION=RX
c To include derivatives of rx use OPTION=RX



c Use this next macro to declare the statement functions that are defined below
c To include derivatives of rx use OPTION=RX


c Define statement functions for difference approximations of order 4 
c To include derivatives of rx use OPTION=RX
c To include derivatives of rx use OPTION=RX


! ---------------------------------------------------------------------------
! Macro : beginLoopsMask
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : endLoopsMask
! ---------------------------------------------------------------------------







! This macro is used for variable dissipation in 2D

! This macro is used for variable dissipation in 3D

! This macro is used for variable dissipation in 3D



! Optionally add the forcing terms

! Optionally add the forcing terms
! Optionally solve for E or H or both


! Optionally add the dissipation and or forcing terms


! Optionally add the dissipation and or forcing terms
! Optionally solve for E or H or both

! The next macro is used for curvilinear girds where the Laplacian term is precomputed.

! -------------------------------------------------------------------------------------------
! The next macro is used for curvilinear girds where the Laplacian term is precomputed.
! Optionally add dissipation too
! -------------------------------------------------------------------------------------------


! -------------------------------------------------------------------------------------------
! The next macro is used for curvilinear girds where the Laplacian term is precomputed.
! Optionally add dissipation too
! -------------------------------------------------------------------------------------------

! -------------------------------------------------------------------------------------------
! The next macro is used for curvilinear girds where the Laplacian term is precomputed.
! Optionally add dissipation too
! -------------------------------------------------------------------------------------------











! ** evaluate the laplacian on the 9 points centred at (i1,i2,i3)


! ** evaluate the square of the Laplacian for a component ****

! ** evaluate the square of the Laplacian for [ex,ey,hz] ****

! ==== loops for curvilinear, with forcing, dissipation in 2D
! Optionally add the dissipation and or forcing terms

! ==== loops for curvilinear, with forcing, dissipation in 2D
! Optionally add the dissipation and or forcing terms

! **********************************************************************************
! Macro ADV_MAXWELL:
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
!  GRIDTYPE : rectangular, curvilinear
! **********************************************************************************






! some newer version of these are created in advOptNew.bf

!      buildFile(advMx2dOrder2r,2,2,rectangular)
!      buildFile(advMx3dOrder2r,3,2,rectangular)

!      buildFile(advMx2dOrder2c,2,2,curvilinear)
!      buildFile(advMx3dOrder2c,3,2,curvilinear)

!      buildFile(advMx2dOrder4r,2,4,rectangular)
!      buildFile(advMx3dOrder4r,3,4,rectangular)

!      buildFile(advMx2dOrder4c,2,4,curvilinear)
!      buildFile(advMx3dOrder4c,3,4,curvilinear)


       ! build these for testing symmetric operators -- BC's not implemented yet


       ! build these for testing symmetric operators -- BC's not implemented yet



! build an empty version of high order files so we do not have to compile the full version









      subroutine advMaxwell(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,
     & nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx,  um,u,un,f,fa, v,vvt2,ut3,
     & vvt4,ut5,ut6,ut7, bc, dis, varDis, ipar, rpar, ierr )
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
      integer useSosupDissipation,sosupDissipationOption
      integer updateSolution,updateDissipation,computeUt
      integer useWhereMask,solveForE,solveForH,grid
      integer ex,ey,ez, hx,hy,hz

      integer rectangular,curvilinear
      parameter( rectangular=0, curvilinear=1 )
!...........end   statement functions


      ! write(*,*) 'Inside advMaxwell...'

      orderOfAccuracy    =ipar(2)
      gridType           =ipar(1)
      useSosupDissipation=ipar(34)
      sosupDissipationOption=ipar(35)
      updateSolution        =ipar(36)
      updateDissipation     =ipar(37)
      computeUt             =ipar(38)

      if( useSosupDissipation.eq.0 .or. (updateSolution.eq.1 .and. 
     & updateDissipation.eq.0 .and.computeUt.eq.0 ) )then

       ! -- FD schemes : dispersive and non-dispersive -- 
       !   These are also used when the sosup dissipation is applied in a separate stage 

       if( orderOfAccuracy.eq.2 )then

        if( nd.eq.2 .and. gridType.eq.rectangular ) then
          call advMx2dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
          call advMx2dOrder2c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if( nd.eq.3 .and. gridType.eq.rectangular ) then
          call advMx3dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
          call advMx3dOrder2c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else
          stop 2271
        end if

       else if( orderOfAccuracy.eq.4 ) then
        if( nd.eq.2 .and. gridType.eq.rectangular )then
          call advMx2dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(nd.eq.2 .and. gridType.eq.curvilinear )then
          call advMx2dOrder4c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.rectangular )then
          call advMx3dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
          call advMx3dOrder4c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
       else
         stop 8843
       end if

!
       else if( orderOfAccuracy.eq.6 ) then
        if( nd.eq.2 .and. gridType.eq.rectangular )then
          call advMx2dOrder6r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(nd.eq.2 .and. gridType.eq.curvilinear )then
          call advMx2dOrder6c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.rectangular )then
          call advMx3dOrder6r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
          call advMx3dOrder6c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
       else
         stop 8843
       end if

       else if( orderOfAccuracy.eq.8 ) then

        if( nd.eq.2 .and. gridType.eq.rectangular )then
          call advMx2dOrder8r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(nd.eq.2 .and. gridType.eq.curvilinear )then
          call advMx2dOrder8c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.rectangular )then
          call advMx3dOrder8r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
          call advMx3dOrder8c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
       else
         stop 8843
       end if

       else
        write(*,'(" advMaxwell:ERROR: un-implemented order of accuracy 
     & =",i6)') orderOfAccuracy
          ! '
        stop 11122
       end if

      else
        ! ========= FINITE DIFFERENCE WITH UPWIND DISSIPATION ========
        ! *NEW July 23, 2017

       if( orderOfAccuracy.eq.2 )then

        if( nd.eq.2 .and. gridType.eq.rectangular ) then
          call advMxUp2dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
          call advMxUp2dOrder2c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if( nd.eq.3 .and. gridType.eq.rectangular ) then
          call advMxUp3dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
          call advMxUp3dOrder2c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else
          stop 6629
        end if

       else if( orderOfAccuracy.eq.4 ) then
        if( nd.eq.2 .and. gridType.eq.rectangular )then
          call advMxUp2dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(nd.eq.2 .and. gridType.eq.curvilinear )then
          call advMxUp2dOrder4c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.rectangular )then
          call advMxUp3dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
          call advMxUp3dOrder4c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rx, um,u,un,f,fa, v,vvt2,
     & ut3,vvt4,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, ierr )
       else
         stop 9255
       end if

      else
        write(*,'(" advMaxwell:ERROR: un-implemented order of accuracy 
     & =",i6)') orderOfAccuracy
          ! '
        stop 20208
       end if


      end if

      return
      end








