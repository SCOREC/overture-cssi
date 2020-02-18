! This file automatically generated from getDiv.bf with bpp.
!
! Optimized routine to compute the divergence 
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



! =========================================================================================
! Macro: Compute induced fields D and B for the BI-ANISTROPIC equations
!
!   POLAR -- polarization: TEZ or NONE
! ========================================================================================

! =========================================================================================
! Macro: Compute induced fields D and B for the BI-ANISTROPIC GDM equations
!
!   POLAR -- polarization: TEZ or NONE
! ========================================================================================



! =========================================================================================
!
!   Compute the divergence of D and B for the BA MAXWELL
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================


! =========================================================================================
!
!   Compute the divergence of E for the isotropic MAXWELL
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================




! **********************************************************************************
! Macro GET_DIVERGENCE
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
!  GRIDTYPE : rectangular, curvilinear
! **********************************************************************************












      subroutine getDiv(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,
     & nd3a,nd3b,nd4a,nd4b,mask,rsxy,  u,p,v, divD,divB, K0, matMask,
     & ipar, rpar, ierr )
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
          call getDiv2dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,  u,p,v, divD,divB, K0,
     &  matMask,ipar, rpar, ierr )
        else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
          write(*,*) 'getDiv -- unimplemented option'
          stop 1234
        else if( nd.eq.3 .and. gridType.eq.rectangular ) then
          call getDiv3dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,  u,p,v, divD,divB, K0,
     &  matMask,ipar, rpar, ierr )
        else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
          write(*,*) 'getDiv -- unimplemented option'
          stop 1234
        else
          stop 2271
        end if

       else if( orderOfAccuracy.eq.4 ) then

        if( nd.eq.2 .and. gridType.eq.rectangular )then
          ! write(*,'(" call getDiv2dOrder4r...")')
          call getDiv2dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,  u,p,v, divD,divB, K0,
     &  matMask,ipar, rpar, ierr )
        else if(nd.eq.2 .and. gridType.eq.curvilinear )then
          write(*,*) 'getDiv -- unimplemented option'
          stop 1234
        else if(  nd.eq.3 .and. gridType.eq.rectangular )then
          call getDiv3dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rsxy,  u,p,v, divD,divB, K0,
     &  matMask,ipar, rpar, ierr )
        else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
          write(*,*) 'getDiv -- unimplemented option'
          stop 1234

        else
          stop 8843
        end if

       else

        write(*,'(" advMaxwell:ERROR: un-implemented order of accuracy 
     & =",i6)') orderOfAccuracy
          ! '
        stop 11122

      end if


      return
      end








