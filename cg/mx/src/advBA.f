! This file automatically generated from advBA.bf with bpp.
!
! Optimized routines for BA Maxwell
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


! ======================================================================================
!   Evaluate the TZ exact solution in 2D
! ======================================================================================

! ======================================================================================
!   Evaluate the TZ exact solution in 3D
! ======================================================================================


! ---------------------------------------------------------------------------
! Macro : beginLoopsMask
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : endLoopsMask
! ---------------------------------------------------------------------------







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











! ** evaluate the laplacian on the 9 points centred at (i1,i2,i3)


! ** evaluate the square of the Laplacian for a component ****

! ** evaluate the square of the Laplacian for [ex,ey,hz] ****

! ==== loops for curvilinear, with forcing, dissipation in 2D
! Optionally add the dissipation and or forcing terms

! ==== loops for curvilinear, with forcing, dissipation in 2D
! Optionally add the dissipation and or forcing terms

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (4th-order difference used with 2nd-order scheme) 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (3D - 4th-order difference used with 2nd-order scheme) 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (6th-order difference used with 4th-order scheme)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (3D 6th-order difference used with 4th-order scheme)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! =========================================================================================
! Macro: Update isotropic equations 2D, Order 2
! ========================================================================================

! =========================================================================================
! Macro: Update BI-ANISTROPIC equations 2D, Order 2
! ========================================================================================


! =========================================================================================
! Macro: Update isotropic equations 2D, Order 4
! ========================================================================================


! =========================================================================================
! Macro: Update BI-ANISTROPIC equations 2D, Order 4
! ========================================================================================


! =========================================================================================
! Macro: Update BI-ANISTROPIC equations 3D, Order 2
! ========================================================================================

! =========================================================================================
! Macro: Update BI-ANISTROPIC equations 3D, Order 4
! ========================================================================================



!- ! =========================================================================================
!- ! Macro: Update BI-ANISTROPIC GDM equations 2D, Order 2
!- ! ========================================================================================
!- #beginMacro updateBAGDM2dOrder2Old()
!-   ! eps* u_t = w_y
!-   ! eps* v_t = - w_x
!-   ! mu* w_t = u_y - v_x 
!- 
!-   if( t.lt.2*dt )then
!-     write(*,'("advBA: *** advance BA GDM, 2nd-order, rectangular... t=",e10.2)') t
!-     write(*,'("  ++++ solveForAllFields,methodOfLines=",2i2)') solveForAllFields,methodOfLines
!-   end if
!-   mr=0
!-   
!-   if( solveForAllFields.ne.0 )then
!- 
!-     ! --- SOLVE FOR ALL FIELDS ---
!- 
!-     if( .not.methodOfLines )then
!- 
!-       ! --- TAYLOR TIME-STEPPING --- 
!-       stop 111
!-        
!-     else
!- 
!-   if( t.lt.2*dt )then
!-     write(*,'("    +++++++ (2) ")') 
!-   end if
!- 
!-       ! --- METHOD OF LINES (RK) ---
!-       if( forcingOption.eq.twilightZoneForcing )then
!-          ! zero out some forcing terms 
!-          do m=0,numPolarizationTerms-1
!-            pv(m)=0.
!-          end do
!-       end if    
!-       beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
!-   
!-         if( addForcing.ne.0 )then  ! do this for now *fix me*
!-           do m=0,5
!-             fv(m)=f(i1,i2,i3,m)
!-           end do
!-         end if 
!-   
!-         if( numberOfMaterialRegions.gt.1 )then
!-           mr = matMask(i1,i2,i3)
!-           if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
!-              stop 9999
!-           end if            
!-         end if       
!- 
!-         if( forcingOption.eq.twilightZoneForcing )then
!- 
!-           if( nd.eq.2 )then
!-             do m=0,5
!-               ec = m 
!-               OGDERIV2D( 0,0,0,0,i1,i2,i3,t, ec, ev(m)  )
!-               OGDERIV2D( 1,0,0,0,i1,i2,i3,t, ec, evt(m) )
!-             end do 
!-             ! eval the polarization terms and time derivatives 
!-             do m=0,numPolarizationTerms-1
!-               pc = m+6  ! TZ index 
!-               OGDERIV2D( 0,0,0,0,i1,i2,i3,t, pc, pv(m)  )
!-               OGDERIV2D( 1,0,0,0,i1,i2,i3,t, pc, pvt(m) )
!-             end do    
!-           else
!-             stop 3333
!-             ! OGDERIV3D( 0,0,0,0,i1,i2,i3,t, ec, e0  )
!-             ! OGDERIV3D( 1,0,0,0,i1,i2,i3,t, ec, e0t )
!-          end if
!- 
!-          !write(*,'(" i1,i2=",2i3," xy=",2(f5.2,1x))') i1,i2,xy(i1,i2,i3,0),xy(i1,i2,i3,1)
!-          !write(*,'(" i1,i2=",2i3," ev=",6(f6.3,1x))') i1,i2,(ev(m),m=0,5)
!-          !write(*,'(" i1,i2=",2i3," pv=",10(f6.3,1x))') i1,i2,(pv(m),m=0,numPolarizationTerms-1)
!- 
!-           pc=0  ! P
!-           qc=1  ! Q = P.t 
!-           do k1=1,6
!-             ec=k1-1 ! E or H component 
!-             do k2=1,6
!-               do n=1,Np(k1,k2,mr)
!-                 a0 = gdmPar(1,n,k1,k2,mr)
!-                 a1 = gdmPar(2,n,k1,k2,mr)
!-                 b0 = gdmPar(3,n,k1,k2,mr)
!-                 b1 = gdmPar(4,n,k1,k2,mr)
!- 
!-                 ! TZ forcing for polarization equations: 
!-                 fp(pc) = pvt(pc) - pv(qc)
!-                 fp(qc) = pvt(qc) - (  a0*ev(ec) + a1*evt(ec) - b0*pv(pc)- b1*pv(qc) )
!- 
!-                 pc=pc+2
!-                 qc=qc+2
!-               end do
!-             end do
!-           end do 
!-        else
!-           ! no TZ : set polarization forcing to zero 
!-           do m=0,numPolarizationTerms-1
!-             fp(m)=0.
!-           end do 
!-        end if
!- 
!- 
!-         ! compute components of the curl(H) and -curl(E)
!-         curl(0) =  uy22r(i1,i2,i3,hz)
!-         curl(1) = -ux22r(i1,i2,i3,hz)
!-         curl(2) = (ux22r(i1,i2,i3,hy)-uy22r(i1,i2,i3,hx))
!- 
!-         curl(3) = -uy22r(i1,i2,i3,ez)
!-         curl(4) =  ux22r(i1,i2,i3,ez)
!-         curl(5) =-(ux22r(i1,i2,i3,ey)-uy22r(i1,i2,i3,ex))
!- 
!-         ! ---- Compute q = p.t = SUM_k2 SUM_n  p(i1,i2,i3, n,k1,k2, qc )
!-         pc=0  ! P
!-         qc=1  ! Q = P.t 
!-         do m=0,5
!-           ptSum(m)=0.   ! local total Pt or Mt 
!-         end do 
!-         do k1=1,6
!-           ec=k1-1 ! E or H component 
!-           do k2=1,6
!-             do n=1,Np(k1,k2,mr)
!-               ptSum(ec) = ptSum(ec) + p(i1,i2,i3,qc) - pv(qc)
!-               qc=qc+2
!-             end do
!-           end do
!-           ! subtract off P.t = sum Q_m
!-           curl(ec) = curl(ec) - ptSum(ec)
!-         end do 
!- 
!-         if( debug.gt.3 )then
!-           write(*,'(" i1,i2=",2i3," ptSum=",6(f6.3,1x))') i1,i2,(ptSum(m),m=0,5)
!-         end if
!- 
!-         do m=0,5 
!-           un(i1,i2,i3,m)= K0i(m,0,mr)*curl(0) + K0i(m,1,mr)*curl(1) + K0i(m,2,mr)*curl(2) + !-                           K0i(m,3,mr)*curl(3) + K0i(m,4,mr)*curl(4) + K0i(m,5,mr)*curl(5) + fv(m)
!-         end do
!- 
!-         if( .true. )then
!-           write(*,'(" i1,i2=",2i3," ut=",6(f6.3,1x))') i1,i2,(un(i1,i2,i3,m),m=0,5)
!-         end if
!- 
!-         ! --- compute time derivatives of P and Q
!-         ! p.t = q
!-         ! q.t = a0*E + a1*Et - b0*p - b1*q   
!-         ! or 
!-         ! q.t = a0*H + a1*Ht - b0*p - b1*q   
!-         pc=0  ! P
!-         qc=1  ! Q = P.t 
!-         do k1=1,6
!-           ec=k1-1 ! E or H component 
!-           do k2=1,6
!-             do n=1,Np(k1,k2,mr)
!-               a0 = gdmPar(1,n,k1,k2,mr)
!-               a1 = gdmPar(2,n,k1,k2,mr)
!-               b0 = gdmPar(3,n,k1,k2,mr)
!-               b1 = gdmPar(4,n,k1,k2,mr)
!- 
!-               pct=pc+6  ! p.t is stored in un here 
!-               qct=qc+6  ! q.t is stored in un here
!- 
!-               un(i1,i2,i3,pct) = p(i1,i2,i3,qc) + fp(pc) 
!-               un(i1,i2,i3,qct) = a0*u(i1,i2,i3,ec) + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
!-               ! test: 
!-               ! un(i1,i2,i3,qct) = a0*ev(ec)         + a1*un(i1,i2,i3,ec) - b0*p(i1,i2,i3,pc)- b1*p(i1,i2,i3,qc) + fp(qc) 
!- 
!-               if( debug.gt.3 )then
!-                 write(*,'(" i1,i2=",2i4," k1,k2=",2i2," Np=",i3," ec=",i1," t=",e9.2)') i1,i2,k1,k2,Np(k1,k2,mr),ec,t
!-                 write(*,'("   ... pc,qc,ec=",3i3," a0,a1,b0,b1=",4e10.2)') pc,qc,ec,a0,a1,b0,b1
!-                 write(*,'("   ... p ,q =",2f7.3)') p(i1,i2,i3,pc),p(i1,i2,i3,qc)
!-               end if 
!-               ! write(*,'("   ... p ,pe =",2f7.3," q ,qe =",2f7.3)') p(i1,i2,i3,pc),pv(pc),p(i1,i2,i3,qc),pv(qc)
!-               ! write(*,'("   ... E ,Ee =",2f7.3," Et,Ete =",2f7.3)') u(i1,i2,i3,ec),ev(ec),un(i1,i2,i3,ec),evt(ec)
!-               ! write(*,'("   ... pt,pte=",2f7.3," qt,qte=",2f7.3)') un(i1,i2,i3,pct),pvt(pc),un(i1,i2,i3,qct),pvt(qc)
!- 
!-               ! test: -- set to exact time-derivative 
!-               ! un(i1,i2,i3,pct) = pvt(pc)
!-               ! un(i1,i2,i3,qct) = pvt(qc)
!-              
!-               pc=pc+2
!-               qc=qc+2
!-             end do
!-           end do
!-         end do 
!- 
!-       endLoopsMask()
!- 
!-       stop  4444
!- 
!-     end if 
!- 
!-   else
!- 
!- 
!-      stop 5555
!- 
!-      ! TEz polarization
!-     if( .not.methodOfLines )then
!-       stop 111
!- 
!-       beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
!-   
!-        if( addForcing.ne.0 )then  ! do this for now *fix me*
!-          fv(ex)=f(i1,i2,i3,ex)
!-          fv(ey)=f(i1,i2,i3,ey)
!-          fv(hz)=f(i1,i2,i3,hz)
!-        else
!-          do c=ex,hz
!-            fv(c)= 0.
!-          end do
!-        end if 
!-   
!-        un(i1,i2,i3,ex)=u(i1,i2,i3,ex) + (dt/eps)*uy22r(i1,i2,i3,hz) +.5*cdtsq*lap2d2(i1,i2,i3,ex) + dt*fv(ex)
!-        un(i1,i2,i3,ey)=u(i1,i2,i3,ey) - (dt/eps)*ux22r(i1,i2,i3,hz) +.5*cdtsq*lap2d2(i1,i2,i3,ey) + dt*fv(ey)
!-        un(i1,i2,i3,hz)=u(i1,i2,i3,hz) + (dt/mu)*(uy22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ey)) !-                                                                    +.5*cdtsq*lap2d2(i1,i2,i3,hz) + dt*fv(hz)
!-       endLoopsMask()
!-     else
!-       ! METHOD OF LINES -- TEz 
!- 
!-       beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
!-   
!-        if( addForcing.ne.0 )then  ! do this for now *fix me*
!-          fv(ex)=f(i1,i2,i3,ex)
!-          fv(ey)=f(i1,i2,i3,ey)
!-          fv(hz)=f(i1,i2,i3,hz)
!-        else
!-          do c=ex,hz
!-            fv(c)= 0.
!-          end do
!-        end if 
!-   
!-         ! compute components of the curl(H) and -curl(E)
!-         curl(0) =  uy22r(i1,i2,i3,hz)
!-         curl(1) = -ux22r(i1,i2,i3,hz)
!-         ! curl(2) = (ux22r(i1,i2,i3,hy)-uy22r(i1,i2,i3,hx))
!- 
!-         ! curl(3) = -uy22r(i1,i2,i3,ez)
!-         ! curl(4) =  ux22r(i1,i2,i3,ez)
!-         curl(2) =-(ux22r(i1,i2,i3,ey)-uy22r(i1,i2,i3,ex))
!- 
!-         if( numberOfMaterialRegions.gt.1 )then
!-           mr = matMask(i1,i2,i3)
!-           if( mr.lt.0 .or. mr.ge.numberOfMaterialRegions )then  ! do this for now 
!-              stop 9999
!-           end if            
!-         end if 
!-         do m=0,2
!-           un(i1,i2,i3,m)= Ki(m,0,mr)*curl(0) + Ki(m,1,mr)*curl(1) + Ki(m,2,mr)*curl(2) + fv(m)
!-         end do 
!- 
!-         ! un(i1,i2,i3,ex)= + (1./eps)*uy22r(i1,i2,i3,hz) + fv(ex)
!-         ! un(i1,i2,i3,ey)= - (1./eps)*ux22r(i1,i2,i3,hz) + fv(ey)
!-         ! un(i1,i2,i3,hz)= + (1./mu)*(uy22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ey)) + fv(hz)
!-       endLoopsMask()
!- 
!-     end if
!- 
!-   end if
!- 
!- #endMacro


! =========================================================================================
! Macro: Update BI-ANISTROPIC GDM equations
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================


! =========================================================================================
! **NEW VERSION** Try to optimized
! Macro: Update BI-ANISTROPIC GDM equations
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================




! **********************************************************************************
! Macro ADV_MAXWELL:
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
!  GRIDTYPE : rectangular, curvilinear
! **********************************************************************************








!      buildFile(advMx2dOrder2c,2,2,curvilinear)
!      buildFile(advMx3dOrder2c,3,2,curvilinear)


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



      subroutine advBA(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,
     & nd3a,nd3b,nd4a,nd4b,mask,rx,  um,u,un,f,fa, K0i, matMask, pm,p,
     & pn,xy,ut5,ut6,ut7, bc, dis, varDis, ipar, rpar, ierr )
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

      real ut5(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real ut6(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      real ut7(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
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
          call advBA2dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rx, um,u,un,f,fa,K0i,
     & matMask, pm,p,pn,xy,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, 
     & ierr )
        else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
          write(*,*) 'advBA -- unimplemented option'
          stop 1234
        else if( nd.eq.3 .and. gridType.eq.rectangular ) then
          ! BA version 
          call advBA3dOrder2r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rx, um,u,un,f,fa,K0i,
     & matMask, pm,p,pn,xy,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, 
     & ierr )
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
          call advBA2dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rx, um,u,un,f,fa,K0i,
     & matMask, pm,p,pn,xy,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, 
     & ierr )
        else if(nd.eq.2 .and. gridType.eq.curvilinear )then
          write(*,*) 'advBA -- unimplemented option'
          stop 1234
        else if(  nd.eq.3 .and. gridType.eq.rectangular )then
          ! BA version 
          call advBA3dOrder4r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,
     & nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,rx, um,u,un,f,fa,K0i,
     & matMask, pm,p,pn,xy,ut5,ut6,ut7,bc, dis,varDis, ipar, rpar, 
     & ierr )
        else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
          write(*,*) 'advBA -- unimplemented option'
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








