! This file automatically generated from evalDispersionRelation.bf with bpp.
! =====================================================================================================
! 
! -----------------------------------------------------------------------------------
! Evaluate the dispersion relation for the generalized dispersion model (GDM)
!  With multiple polarization vectors 
! -----------------------------------------------------------------------------------
!
!       E_tt - c^2 Delta(E) = (1/eps)* P_tt
!       P_tt + b1 P_1 + b0 = eps*( a0*E + a1*E_t )
!
!  Input:
!      mode : mode to choose, i.e. which root to choose. If mode=-1 then the default root is chosen.
!                The default root is the one with largest NEGATIVE imaginary part.
!      Np : number of polarization vectors  
!      c,k, 
!      a0(0:Np-1), a1(0:Np-1), b0(0:Np-1), b1(0:Np-1) : GDM parameters
! Output:
!      sr(0:nEig-1),si(0:nEig-1) : real and imaginary part of eigenvalues: nEig = 2*NpP+2 
!      srm, sim : eigenvalue with largest imaginary part 
!      chir(0:Np-1),chii(0:Np-1) : real and imaginary parts of the electric susceptibility 
!                    chi(j) = ( a0(j) + a1(s)*s)/(s^2 + b1(j)*s + b0(j) )
!      chiSumr, chiSumi : real and imaginary parts of chiSum = sum_j chi(j) 
!
! =====================================================================================================
      subroutine evalEigGDM( mode, Np, c,k, a0,a1,b0,b1, sr,si, srm,
     & sim,chir,chii, chiSumr,chiSumi )

      implicit none

      integer mode, Np
      real c,k, srm,sim
      real a0(0:*),a1(0:*),b0(0:*),b1(0:*), chir(0:*),chii(0:*), 
     & chiSumr, chiSumi
      real sr(0:*),si(0:*)

      ! local variables 
      real ck,ck2, eps
      integer nd,lda, lwork, i, iMode, j
      ! lwork>= 3*lda : for good performance lwork must generally be larger
      parameter( lda=10, lwork=10*lda )
      real a(0:lda-1,0:lda-1), work(lwork), vr(1), vl(1)
      integer info

      complex*16 s, chi, chiSum

      ck=c*k
      ck2=ck**2
      ! Companion matrix for GDM model Np File written by CG/DMX/matlab/gdm.maple
      if( Np .eq. 0 )then
! Companion matrix for GDM model Np=0, File written by CG/DMX/matlab/gdm.maple
      a(0,0) = 0
      a(0,1) = -ck2
      a(1,0) = 1
      a(1,1) = 0
      else if( Np .eq. 1 )then
! Companion matrix for GDM model Np=1, File written by CG/DMX/matlab/gdm.maple
      a(0,0) = 0
      a(0,1) = 0
      a(0,2) = 0
      a(0,3) = -ck2*b0(0)
      a(1,0) = 1
      a(1,1) = 0
      a(1,2) = 0
      a(1,3) = -ck2*b1(0)
      a(2,0) = 0
      a(2,1) = 1
      a(2,2) = 0
      a(2,3) = -ck2-b0(0)-a0(0)
      a(3,0) = 0
      a(3,1) = 0
      a(3,2) = 1
      a(3,3) = -b1(0)-a1(0)
      else if( Np .eq. 2 )then
! Companion matrix for GDM model Np=2, File written by CG/DMX/matlab/gdm.maple
      a(0,0) = 0
      a(0,1) = 0
      a(0,2) = 0
      a(0,3) = 0
      a(0,4) = 0
      a(0,5) = -ck2*b0(0)*b0(1)
      a(1,0) = 1
      a(1,1) = 0
      a(1,2) = 0
      a(1,3) = 0
      a(1,4) = 0
      a(1,5) = -ck2*b0(0)*b1(1)-ck2*b1(0)*b0(1)
      a(2,0) = 0
      a(2,1) = 1
      a(2,2) = 0
      a(2,3) = 0
      a(2,4) = 0
      a(2,5) = -ck2*b0(0)-ck2*b1(0)*b1(1)-(ck2+b0(0))*b0(1)-a0(0)*b0(1)
     & -a0(1)*b0(0)
      a(3,0) = 0
      a(3,1) = 0
      a(3,2) = 1
      a(3,3) = 0
      a(3,4) = 0
      a(3,5) = -ck2*b1(0)-(ck2+b0(0))*b1(1)-b1(0)*b0(1)-a0(0)*b1(1)-a1(
     & 0)*b0(1)-a0(1)*b1(0)-a1(1)*b0(0)
      a(4,0) = 0
      a(4,1) = 0
      a(4,2) = 0
      a(4,3) = 1
      a(4,4) = 0
      a(4,5) = -b1(0)*b1(1)-a1(0)*b1(1)-a1(1)*b1(0)-b0(1)-a0(1)-b0(0)-
     & a0(0)-ck2
      a(5,0) = 0
      a(5,1) = 0
      a(5,2) = 0
      a(5,3) = 0
      a(5,4) = 1
      a(5,5) = -b1(0)-b1(1)-a1(0)-a1(1)
      else if( Np .eq. 3 )then
! Companion matrix for GDM model Np=3, File written by CG/DMX/matlab/gdm.maple
      a(0,0) = 0
      a(0,1) = 0
      a(0,2) = 0
      a(0,3) = 0
      a(0,4) = 0
      a(0,5) = 0
      a(0,6) = 0
      a(0,7) = -ck2*b0(0)*b0(1)*b0(2)
      a(1,0) = 1
      a(1,1) = 0
      a(1,2) = 0
      a(1,3) = 0
      a(1,4) = 0
      a(1,5) = 0
      a(1,6) = 0
      a(1,7) = -ck2*b0(0)*b0(1)*b1(2)-(ck2*b0(0)*b1(1)+ck2*b1(0)*b0(1))
     & *b0(2)
      a(2,0) = 0
      a(2,1) = 1
      a(2,2) = 0
      a(2,3) = 0
      a(2,4) = 0
      a(2,5) = 0
      a(2,6) = 0
      a(2,7) = -ck2*b0(0)*b0(1)-(ck2*b0(0)*b1(1)+ck2*b1(0)*b0(1))*b1(2)
     & -(ck2*b0(0)+ck2*b1(0)*b1(1)+(ck2+b0(0))*b0(1))*b0(2)-a0(0)*b0(
     & 1)*b0(2)-a0(1)*b0(0)*b0(2)-a0(2)*b0(0)*b0(1)
      a(3,0) = 0
      a(3,1) = 0
      a(3,2) = 1
      a(3,3) = 0
      a(3,4) = 0
      a(3,5) = 0
      a(3,6) = 0
      a(3,7) = -ck2*b0(0)*b1(1)-ck2*b1(0)*b0(1)-(ck2*b0(0)+ck2*b1(0)*
     & b1(1)+(ck2+b0(0))*b0(1))*b1(2)-(ck2*b1(0)+(ck2+b0(0))*b1(1)+b1(
     & 0)*b0(1))*b0(2)-a0(0)*b0(1)*b1(2)-(a0(0)*b1(1)+a1(0)*b0(1))*b0(
     & 2)-a0(1)*b0(0)*b1(2)-(a0(1)*b1(0)+a1(1)*b0(0))*b0(2)-a0(2)*b0(
     & 0)*b1(1)-(a0(2)*b1(0)+a1(2)*b0(0))*b0(1)
      a(4,0) = 0
      a(4,1) = 0
      a(4,2) = 0
      a(4,3) = 1
      a(4,4) = 0
      a(4,5) = 0
      a(4,6) = 0
      a(4,7) = -ck2*b0(0)-ck2*b1(0)*b1(1)-(ck2+b0(0))*b0(1)-(ck2*b1(0)+
     & (ck2+b0(0))*b1(1)+b1(0)*b0(1))*b1(2)-(b1(0)*b1(1)+b0(1)+b0(0)+
     & ck2)*b0(2)-a0(0)*b0(1)-(a0(0)*b1(1)+a1(0)*b0(1))*b1(2)-(a1(0)*
     & b1(1)+a0(0))*b0(2)-a0(1)*b0(0)-(a0(1)*b1(0)+a1(1)*b0(0))*b1(2)-
     & (a1(1)*b1(0)+a0(1))*b0(2)-a0(2)*b0(0)-(a0(2)*b1(0)+a1(2)*b0(0))
     & *b1(1)-(a1(2)*b1(0)+a0(2))*b0(1)
      a(5,0) = 0
      a(5,1) = 0
      a(5,2) = 0
      a(5,3) = 0
      a(5,4) = 1
      a(5,5) = 0
      a(5,6) = 0
      a(5,7) = -ck2*b1(0)-(ck2+b0(0))*b1(1)-b1(0)*b0(1)-(b1(0)*b1(1)+
     & b0(1)+b0(0)+ck2)*b1(2)-(b1(0)+b1(1))*b0(2)-a0(0)*b1(1)-a1(0)*
     & b0(1)-(a1(0)*b1(1)+a0(0))*b1(2)-a1(0)*b0(2)-a0(1)*b1(0)-a1(1)*
     & b0(0)-(a1(1)*b1(0)+a0(1))*b1(2)-a1(1)*b0(2)-a0(2)*b1(0)-a1(2)*
     & b0(0)-(a1(2)*b1(0)+a0(2))*b1(1)-a1(2)*b0(1)
      a(6,0) = 0
      a(6,1) = 0
      a(6,2) = 0
      a(6,3) = 0
      a(6,4) = 0
      a(6,5) = 1
      a(6,6) = 0
      a(6,7) = -b1(0)*b1(1)-b0(1)-b0(0)-ck2-(b1(0)+b1(1))*b1(2)-b0(2)-
     & a1(1)*b1(2)-a1(0)*b1(2)-a1(2)*b1(1)-a1(2)*b1(0)-a1(0)*b1(1)-a1(
     & 1)*b1(0)-a0(2)-a0(1)-a0(0)
      a(7,0) = 0
      a(7,1) = 0
      a(7,2) = 0
      a(7,3) = 0
      a(7,4) = 0
      a(7,5) = 0
      a(7,6) = 1
      a(7,7) = -b1(0)-b1(1)-b1(2)-a1(0)-a1(1)-a1(2)
      else if( Np .eq. 4 )then
! Companion matrix for GDM model Np=4, File written by CG/DMX/matlab/gdm.maple
      a(0,0) = 0
      a(0,1) = 0
      a(0,2) = 0
      a(0,3) = 0
      a(0,4) = 0
      a(0,5) = 0
      a(0,6) = 0
      a(0,7) = 0
      a(0,8) = 0
      a(0,9) = -ck2*b0(0)*b0(1)*b0(2)*b0(3)
      a(1,0) = 1
      a(1,1) = 0
      a(1,2) = 0
      a(1,3) = 0
      a(1,4) = 0
      a(1,5) = 0
      a(1,6) = 0
      a(1,7) = 0
      a(1,8) = 0
      a(1,9) = -ck2*b0(0)*b0(1)*b0(2)*b1(3)-(ck2*b0(0)*b0(1)*b1(2)+(
     & ck2*b0(0)*b1(1)+ck2*b1(0)*b0(1))*b0(2))*b0(3)
      a(2,0) = 0
      a(2,1) = 1
      a(2,2) = 0
      a(2,3) = 0
      a(2,4) = 0
      a(2,5) = 0
      a(2,6) = 0
      a(2,7) = 0
      a(2,8) = 0
      a(2,9) = -ck2*b0(0)*b0(1)*b0(2)-(ck2*b0(0)*b0(1)*b1(2)+(ck2*b0(0)
     & *b1(1)+ck2*b1(0)*b0(1))*b0(2))*b1(3)-(ck2*b0(0)*b0(1)+(ck2*b0(
     & 0)*b1(1)+ck2*b1(0)*b0(1))*b1(2)+(ck2*b0(0)+ck2*b1(0)*b1(1)+(
     & ck2+b0(0))*b0(1))*b0(2))*b0(3)-a0(0)*b0(1)*b0(2)*b0(3)-a0(1)*
     & b0(0)*b0(2)*b0(3)-a0(2)*b0(0)*b0(1)*b0(3)-a0(3)*b0(0)*b0(1)*b0(
     & 2)
      a(3,0) = 0
      a(3,1) = 0
      a(3,2) = 1
      a(3,3) = 0
      a(3,4) = 0
      a(3,5) = 0
      a(3,6) = 0
      a(3,7) = 0
      a(3,8) = 0
      a(3,9) = -ck2*b0(0)*b0(1)*b1(2)-(ck2*b0(0)*b1(1)+ck2*b1(0)*b0(1))
     & *b0(2)-(ck2*b0(0)*b0(1)+(ck2*b0(0)*b1(1)+ck2*b1(0)*b0(1))*b1(2)
     & +(ck2*b0(0)+ck2*b1(0)*b1(1)+(ck2+b0(0))*b0(1))*b0(2))*b1(3)-(
     & ck2*b0(0)*b1(1)+ck2*b1(0)*b0(1)+(ck2*b0(0)+ck2*b1(0)*b1(1)+(
     & ck2+b0(0))*b0(1))*b1(2)+(ck2*b1(0)+(ck2+b0(0))*b1(1)+b1(0)*b0(
     & 1))*b0(2))*b0(3)-a0(0)*b0(1)*b0(2)*b1(3)-(a0(0)*b0(1)*b1(2)+(
     & a0(0)*b1(1)+a1(0)*b0(1))*b0(2))*b0(3)-a0(1)*b0(0)*b0(2)*b1(3)-(
     & a0(1)*b0(0)*b1(2)+(a0(1)*b1(0)+a1(1)*b0(0))*b0(2))*b0(3)-a0(2)*
     & b0(0)*b0(1)*b1(3)-(a0(2)*b0(0)*b1(1)+(a0(2)*b1(0)+a1(2)*b0(0))*
     & b0(1))*b0(3)-a0(3)*b0(0)*b0(1)*b1(2)-(a0(3)*b0(0)*b1(1)+(a0(3)*
     & b1(0)+a1(3)*b0(0))*b0(1))*b0(2)
      a(4,0) = 0
      a(4,1) = 0
      a(4,2) = 0
      a(4,3) = 1
      a(4,4) = 0
      a(4,5) = 0
      a(4,6) = 0
      a(4,7) = 0
      a(4,8) = 0
      a(4,9) = -ck2*b0(0)*b0(1)-(ck2*b0(0)*b1(1)+ck2*b1(0)*b0(1))*b1(2)
     & -(ck2*b0(0)+ck2*b1(0)*b1(1)+(ck2+b0(0))*b0(1))*b0(2)-(ck2*b0(0)
     & *b1(1)+ck2*b1(0)*b0(1)+(ck2*b0(0)+ck2*b1(0)*b1(1)+(ck2+b0(0))*
     & b0(1))*b1(2)+(ck2*b1(0)+(ck2+b0(0))*b1(1)+b1(0)*b0(1))*b0(2))*
     & b1(3)-(ck2*b0(0)+ck2*b1(0)*b1(1)+(ck2+b0(0))*b0(1)+(ck2*b1(0)+(
     & ck2+b0(0))*b1(1)+b1(0)*b0(1))*b1(2)+(b1(0)*b1(1)+b0(1)+b0(0)+
     & ck2)*b0(2))*b0(3)-a0(0)*b0(1)*b0(2)-(a0(0)*b0(1)*b1(2)+(a0(0)*
     & b1(1)+a1(0)*b0(1))*b0(2))*b1(3)-(a0(0)*b0(1)+(a0(0)*b1(1)+a1(0)
     & *b0(1))*b1(2)+(a1(0)*b1(1)+a0(0))*b0(2))*b0(3)-a0(1)*b0(0)*b0(
     & 2)-(a0(1)*b0(0)*b1(2)+(a0(1)*b1(0)+a1(1)*b0(0))*b0(2))*b1(3)-(
     & a0(1)*b0(0)+(a0(1)*b1(0)+a1(1)*b0(0))*b1(2)+(a1(1)*b1(0)+a0(1))
     & *b0(2))*b0(3)-a0(2)*b0(0)*b0(1)-(a0(2)*b0(0)*b1(1)+(a0(2)*b1(0)
     & +a1(2)*b0(0))*b0(1))*b1(3)-(a0(2)*b0(0)+(a0(2)*b1(0)+a1(2)*b0(
     & 0))*b1(1)+(a1(2)*b1(0)+a0(2))*b0(1))*b0(3)-a0(3)*b0(0)*b0(1)-(
     & a0(3)*b0(0)*b1(1)+(a0(3)*b1(0)+a1(3)*b0(0))*b0(1))*b1(2)-(a0(3)
     & *b0(0)+(a0(3)*b1(0)+a1(3)*b0(0))*b1(1)+(a1(3)*b1(0)+a0(3))*b0(
     & 1))*b0(2)
      a(5,0) = 0
      a(5,1) = 0
      a(5,2) = 0
      a(5,3) = 0
      a(5,4) = 1
      a(5,5) = 0
      a(5,6) = 0
      a(5,7) = 0
      a(5,8) = 0
      a(5,9) = -ck2*b0(0)*b1(1)-ck2*b1(0)*b0(1)-(ck2*b0(0)+ck2*b1(0)*
     & b1(1)+(ck2+b0(0))*b0(1))*b1(2)-(ck2*b1(0)+(ck2+b0(0))*b1(1)+b1(
     & 0)*b0(1))*b0(2)-(ck2*b0(0)+ck2*b1(0)*b1(1)+(ck2+b0(0))*b0(1)+(
     & ck2*b1(0)+(ck2+b0(0))*b1(1)+b1(0)*b0(1))*b1(2)+(b1(0)*b1(1)+b0(
     & 1)+b0(0)+ck2)*b0(2))*b1(3)-(ck2*b1(0)+(ck2+b0(0))*b1(1)+b1(0)*
     & b0(1)+(b1(0)*b1(1)+b0(1)+b0(0)+ck2)*b1(2)+(b1(0)+b1(1))*b0(2))*
     & b0(3)-a0(0)*b0(1)*b1(2)-(a0(0)*b1(1)+a1(0)*b0(1))*b0(2)-(a0(0)*
     & b0(1)+(a0(0)*b1(1)+a1(0)*b0(1))*b1(2)+(a1(0)*b1(1)+a0(0))*b0(2)
     & )*b1(3)-(a0(0)*b1(1)+a1(0)*b0(1)+(a1(0)*b1(1)+a0(0))*b1(2)+a1(
     & 0)*b0(2))*b0(3)-a0(1)*b0(0)*b1(2)-(a0(1)*b1(0)+a1(1)*b0(0))*b0(
     & 2)-(a0(1)*b0(0)+(a0(1)*b1(0)+a1(1)*b0(0))*b1(2)+(a1(1)*b1(0)+
     & a0(1))*b0(2))*b1(3)-(a0(1)*b1(0)+a1(1)*b0(0)+(a1(1)*b1(0)+a0(1)
     & )*b1(2)+a1(1)*b0(2))*b0(3)-a0(2)*b0(0)*b1(1)-(a0(2)*b1(0)+a1(2)
     & *b0(0))*b0(1)-(a0(2)*b0(0)+(a0(2)*b1(0)+a1(2)*b0(0))*b1(1)+(a1(
     & 2)*b1(0)+a0(2))*b0(1))*b1(3)-(a0(2)*b1(0)+a1(2)*b0(0)+(a1(2)*
     & b1(0)+a0(2))*b1(1)+a1(2)*b0(1))*b0(3)-a0(3)*b0(0)*b1(1)-(a0(3)*
     & b1(0)+a1(3)*b0(0))*b0(1)-(a0(3)*b0(0)+(a0(3)*b1(0)+a1(3)*b0(0))
     & *b1(1)+(a1(3)*b1(0)+a0(3))*b0(1))*b1(2)-(a0(3)*b1(0)+a1(3)*b0(
     & 0)+(a1(3)*b1(0)+a0(3))*b1(1)+a1(3)*b0(1))*b0(2)
      a(6,0) = 0
      a(6,1) = 0
      a(6,2) = 0
      a(6,3) = 0
      a(6,4) = 0
      a(6,5) = 1
      a(6,6) = 0
      a(6,7) = 0
      a(6,8) = 0
      a(6,9) = -ck2*b1(0)*b1(1)-(a0(2)*b1(0)+a1(2)*b0(0)+(a1(2)*b1(0)+
     & a0(2))*b1(1)+a1(2)*b0(1))*b1(3)-(a1(2)*b1(1)+a1(2)*b1(0)+a0(2))
     & *b0(3)-a0(3)*b0(0)-(a0(3)*b1(0)+a1(3)*b0(0))*b1(1)-(a1(3)*b1(0)
     & +a0(3))*b0(1)-(a0(3)*b1(0)+a1(3)*b0(0)+(a1(3)*b1(0)+a0(3))*b1(
     & 1)+a1(3)*b0(1))*b1(2)-(a1(3)*b1(1)+a1(3)*b1(0)+a0(3))*b0(2)-(
     & ck2*b1(0)+(ck2+b0(0))*b1(1)+b1(0)*b0(1)+(b1(0)*b1(1)+b0(1)+b0(
     & 0)+ck2)*b1(2)+(b1(0)+b1(1))*b0(2))*b1(3)-(b1(0)*b1(1)+b0(1)+b0(
     & 0)+ck2+(b1(0)+b1(1))*b1(2)+b0(2))*b0(3)-(a0(0)*b1(1)+a1(0)*b0(
     & 1)+(a1(0)*b1(1)+a0(0))*b1(2)+a1(0)*b0(2))*b1(3)-(a1(0)*b1(2)+
     & a1(0)*b1(1)+a0(0))*b0(3)-(a0(1)*b1(0)+a1(1)*b0(0)+(a1(1)*b1(0)+
     & a0(1))*b1(2)+a1(1)*b0(2))*b1(3)-(a1(1)*b1(2)+a1(1)*b1(0)+a0(1))
     & *b0(3)-(ck2*b1(0)+(ck2+b0(0))*b1(1)+b1(0)*b0(1))*b1(2)-(b1(0)*
     & b1(1)+b0(1)+b0(0)+ck2)*b0(2)-(a0(0)*b1(1)+a1(0)*b0(1))*b1(2)-(
     & a1(0)*b1(1)+a0(0))*b0(2)-(a0(1)*b1(0)+a1(1)*b0(0))*b1(2)-(a1(1)
     & *b1(0)+a0(1))*b0(2)-a0(2)*b0(0)-(a0(2)*b1(0)+a1(2)*b0(0))*b1(1)
     & -(a1(2)*b1(0)+a0(2))*b0(1)-(ck2+b0(0))*b0(1)-a0(0)*b0(1)-a0(1)*
     & b0(0)-ck2*b0(0)
      a(7,0) = 0
      a(7,1) = 0
      a(7,2) = 0
      a(7,3) = 0
      a(7,4) = 0
      a(7,5) = 0
      a(7,6) = 1
      a(7,7) = 0
      a(7,8) = 0
      a(7,9) = -(a1(0)*b1(2)+a1(0)*b1(1)+a0(0))*b1(3)-a1(0)*b0(3)-(a1(
     & 1)*b1(2)+a1(1)*b1(0)+a0(1))*b1(3)-a1(1)*b0(3)-(a1(2)*b1(1)+a1(
     & 2)*b1(0)+a0(2))*b1(3)-a1(2)*b0(3)-a0(3)*b1(0)-a1(3)*b0(0)-(a1(
     & 3)*b1(0)+a0(3))*b1(1)-a1(3)*b0(1)-(a1(3)*b1(1)+a1(3)*b1(0)+a0(
     & 3))*b1(2)-a1(3)*b0(2)-(b1(0)*b1(1)+b0(1)+b0(0)+ck2+(b1(0)+b1(1)
     & )*b1(2)+b0(2))*b1(3)-(b1(0)+b1(1)+b1(2))*b0(3)-(b1(0)*b1(1)+b0(
     & 1)+b0(0)+ck2)*b1(2)-(b1(0)+b1(1))*b0(2)-(a1(0)*b1(1)+a0(0))*b1(
     & 2)-a1(0)*b0(2)-(a1(1)*b1(0)+a0(1))*b1(2)-a1(1)*b0(2)-a0(2)*b1(
     & 0)-a1(2)*b0(0)-(a1(2)*b1(0)+a0(2))*b1(1)-a1(2)*b0(1)-(ck2+b0(0)
     & )*b1(1)-b1(0)*b0(1)-a0(0)*b1(1)-a1(0)*b0(1)-a0(1)*b1(0)-a1(1)*
     & b0(0)-ck2*b1(0)
      a(8,0) = 0
      a(8,1) = 0
      a(8,2) = 0
      a(8,3) = 0
      a(8,4) = 0
      a(8,5) = 0
      a(8,6) = 0
      a(8,7) = 1
      a(8,8) = 0
      a(8,9) = -b1(0)*b1(1)-b0(1)-b0(0)-ck2-(b1(0)+b1(1))*b1(2)-b0(2)-(
     & b1(0)+b1(1)+b1(2))*b1(3)-b0(3)-a1(2)*b1(3)-a1(1)*b1(3)-a1(0)*
     & b1(3)-a1(3)*b1(2)-a1(3)*b1(1)-a1(3)*b1(0)-a1(1)*b1(2)-a1(0)*b1(
     & 2)-a1(2)*b1(1)-a1(2)*b1(0)-a1(0)*b1(1)-a1(1)*b1(0)-a0(3)-a0(2)-
     & a0(1)-a0(0)
      a(9,0) = 0
      a(9,1) = 0
      a(9,2) = 0
      a(9,3) = 0
      a(9,4) = 0
      a(9,5) = 0
      a(9,6) = 0
      a(9,7) = 0
      a(9,8) = 1
      a(9,9) = -b1(0)-b1(1)-b1(2)-b1(3)-a1(0)-a1(1)-a1(2)-a1(3)
      else
        write(*,'(" evalEigGDM:ERROR: dispersion relation not 
     & implemented for Np=",i6)') Np
        stop 1234
      end if

      nd = 2*Np+2 ! order of the matrix

      ! Compute eigenvalues of a general non-symmetric matrix 
      call dgeev( 'N','N',nd,a,lda,sr,si,vl,1,vr,1,work,lwork,info )
      write(*,'("evalEigGDM: return from dgeev: info=",i8)') info

      if( .true. )then
          write(*,'(" evalEigGDM: input: mode=",i6)') mode
        do i=0,nd-1
          write(*,'(" evalEigGDM: i=",i3," s=(",1P,e20.12,",", 1P,
     & e20.12,")")') i,sr(i),si(i)
       end do
      end if
      if( mode.ge.0 )then
       ! Take user supplied mode: 
       iMode=min(mode,nd-1)
      else
        ! choose s with largest NEGATIVE imaginary part
        iMode=0
        sim=si(iMode)
        do i=0,nd-1
          if( si(i) .lt. sim )then
            iMode=i
           sim=si(i)
          end if
        end do
      end if
      sim=si(iMode)
      srm=sr(iMode)


      ! chi(j) = ( a0(j) + a1(s)*s)/(s^2 + b1(j)*s + b0(j) )
      ! P = eps*chi(s)*E
      eps=1.e-14 ! FIX ME
      chiSum=0.
      do j=0,Np-1
        s = cmplx(srm,sim)
        if( abs(real(s)) + abs(imag(s)) .gt.eps )then
          chi = (a0(j)+a1(j)*s)/(s**2+b1(j)*s+b0(j))
        else
          chi=0.
        end if

        chiSum= chiSum + chi

        chir(j) = real(chi)
        chii(j) = imag(chi)
      end do

      chiSumr = real(chiSum)
      chiSumi = imag(chiSum)
      call flush(6)

      return
      end




! =====================================================================================================
! 
! -----------------------------------------------------------------------------------
! Evaluate the "INVERSE" dispersion relation (compute k=*(kr,ki) given s=(sr,si) 
! for the generalized dispersion model (GDM) With multiple polarization vectors 
! -----------------------------------------------------------------------------------
!
!       E_tt - c^2 Delta(E) = - (1/eps) P_tt
!       P_tt + b1 P_1 + b0 = eps( a0*E + a1*E_t )
!
!  Input:
!      sr,si : s=(sr,si)
!      Np : number of polarization vectors  
!      c  : 
!      a0(0:Np-1), a1(0:Np-1), b0(0:Np-1), b1(0:Np-1) : GDM parameters
! Output:
!      kr,ki : k=(kr,ki) = complex wave number 
!      chir(0:Np-1),chii(0:Np-1) : real and imaginary parts of the electric susceptibility 
!                    chi(j) = ( a0(j) + a1(s)*s)/(s^2 + b1(j)*s + b0(j) )
!      chiSumr, chiSumi : real and imaginary parts of chiSum = sum_j chi(j) 
!
! =====================================================================================================
      subroutine evalInverseGDM( c, sr,si, Np, a0,a1,b0,b1, kr,ki,chir,
     & chii, chiSumr, chiSumi )

      implicit none

      real c,sr,si
      integer Np
      real a0(0:*),a1(0:*),b0(0:*),b1(0:*), chir(0:*),chii(0:*), 
     & chiSumr, chiSumi
      real kr,ki

      ! local variables 
      real eps
      integer j
      complex*16 s, chi, chiSum , k

      s=cmplx(sr,si)

      chiSum = 0.
      do j=0,Np-1
        ! Phat = eps * chi(s) * Ehat
        chi = (a0(j)+a1(j)*s)/(s**2+b1(j)*s+b0(j))
        chir(j) = real(chi)
        chii(j) = imag(chi)

        chiSum =chiSum + chi
      end do

      ! (ck)^2 = -s^2 - alphaP* s^2 *Sum_k{ chi_k(s) }
      k = csqrt(-s*s - s*s*( chiSum ) )/c

      kr = real(k)
      ki = imag(k)

      chiSumr = real(chiSum)
      chiSumi = imag(chiSum)

      return
      end










!- ! -----------------------------------------------------------------------------------
!- ! Evaluate the dispersion relation for the generalized dispersion model (GDM)
!- ! -----------------------------------------------------------------------------------
!- !
!- !       E_tt - c^2 Delta(E) = -alphaP P_tt
!- !       P_tt + b1 P_1 + b0 = a0*E + a1*E_t 
!- !
!- !  Input:
!- !      c,k, a0,a1,b0,b1,alphaP
!- ! Output:
!- !      sr,si : real and imaginary part of s
!- !      psir,psii : real and imaginary parts of psi : P = psi*E 
!- !
!-       subroutine evalGeneralizedDispersionRelation( c,k, a0,a1,b0,b1,alphaP, sr,si, psir,psii )
!- 
!- 
!-       ! implicit none
!-       implicit complex*16 (t)
!- 
!-       real c,k, a0,a1,b0,b1,alphaP, psir,psii
!-       real ck,ck2,  f ,ap, cki,ck2i 
!-       real sr,si
!-       complex*16 ai, ss,s,cComplex, psi
!- 
!-       ! OLD -- should not be used 
!-       stop 1111 
!- 
!-       ! cComplex = dcmplx(c,0.) ! convert c to complex to force complex arithmetic below
!-       ck=c*k
!-       ck2=ck*ck
!-       cki=1/ck
!-       ck2i=1./ck2 
!- 
!-       ap=alphaP
!-  
!- ! File generated by overtureFramework/cg/mx/codes/gdm.maple
!- 
!- #Include "generalizedDispersionRelation1.h"
!- 
!-       ! The valid root will have an imaginary part close to I*ck 
!-       if( abs(dimag(ss)) < .01*ck )then
!- #Include "generalizedDispersionRelation3.h"     
!-         write(*,'("evalGDM: Use root 3")')
!-         if( abs(dimag(ss)) < .01*ck )then
!-           write(*,'("evalGDM: INVALID root found!?")')
!-           stop 6666
!-         end if
!-       else
!-         write(*,'("evalGDM: Use root 1")')
!-       end if
!- 
!-       s =ss*ck
!- 
!-       if( dimag(s)>0. )then
!-         s=dconjg(s) ! choose right moving wave, imag(s) < 0 
!-       end if 
!- 
!-       sr= dreal(s)
!-       si= dimag(s)
!- 
!-       ! P = psi(s)*E
!-       psi = (a0+a1*s)/(s**2+b1*s+b0)
!-       psir = dreal(psi)
!-       psii = dimag(psi)  ! NOTE: CHOOSE RIGHT MOVING WAVE
!- 
!-       ! check root:
!-       f = cabs((ss**2 + 1)*(ss**2 + cki*b1*ss+ck2i*b0) + ap*ss**2*cki*( a1*ss+ cki*a0 ) )
!- 
!-       write(*,'("evalGDM: c,k,a0,a1,b0,b1=",6(1P,e20.12)," alphaP=",1Pe20.12)') c,k,a0,a1,b0,b1,alphaP
!-       write(*,'("evalGDM: sr,si=",2(1P,e24.15)," |f|=",e12.4)') sr,si,f
!-       write(*,'("evalGDM: psir,psii=",2e12.4)') psir,psii
!- 
!- 
!-       return
!-       end   
!- 
!- 
!- ! Evaluate the dispersion relation for the Drude model 
!- !
!- !    E_tt + c^2 k^2 E = -alphaP*P_tt
!- !
!- !  Input:
!- !      c0,eps,gam,omegap,k : 
!- ! Output:
!- !      reS,imS
!- !
!-       subroutine evalDispersionRelation( c0,eps,gam,omegap,k, reS,imS )
!- 
!- 
!-       implicit none
!- 
!-       real c0,eps,gam,omegap,k, ck2, epsi, om2, det
!-       real reS,imS
!-       complex*16 ai, c, s
!- 
!-       ! OLD -- should not be used 
!-       stop 2222
!- 
!-       ck2=(c0*k)**2
!-       epsi=1./eps
!-       om2=omegap**2
!- 
!-        ! ai=cmplx(0.,1.)  ! i 
!-        c = cmplx(c0,0.) ! convert c to complex to force complex arithmetic below
!- ! File generated by overtureFramework/cg/mx/codes/dispersion.maple
!- ! Here is root 3 from the dispersion relation exp( i*k*x + s*t) .
!- ! s = -1/12*(36*epsi*gam*om2-72*ck2*gam-8*gam^3+12*(12*epsi^3*om2^3-3*epsi^2*gam^2*om2^2+36*ck2*epsi^2*om2^2-60*ck2*epsi*gam^2*om2+12*ck2*gam^4+36*ck2^2*epsi*om2+24*ck2^2*gam^2+12*ck2^3)^(1/2))^(1/3)+3*(1/3*epsi*om2+1/3*ck2-1/9*gam^2)/(36*epsi*gam*om2-72*ck2*gam-8*gam^3+12*(12*epsi^3*om2^3-3*epsi^2*gam^2*om2^2+36*ck2*epsi^2*om2^2-60*ck2*epsi*gam^2*om2+12*ck2*gam^4+36*ck2^2*epsi*om2+24*ck2^2*gam^2+12*ck2^3)^(1/2))^(1/3)-1/3*gam+1/2*I*3^(1/2)*(1/6*(36*epsi*gam*om2-72*ck2*gam-8*gam^3+12*(12*epsi^3*om2^3-3*epsi^2*gam^2*om2^2+36*ck2*epsi^2*om2^2-60*ck2*epsi*gam^2*om2+12*ck2*gam^4+36*ck2^2*epsi*om2+24*ck2^2*gam^2+12*ck2^3)^(1/2))^(1/3)+6*(1/3*epsi*om2+1/3*ck2-1/9*gam^2)/(36*epsi*gam*om2-72*ck2*gam-8*gam^3+12*(12*epsi^3*om2^3-3*epsi^2*gam^2*om2^2+36*ck2*epsi^2*om2^2-60*ck2*epsi*gam^2*om2+12*ck2*gam^4+36*ck2^2*epsi*om2+24*ck2^2*gam^2+12*ck2^3)^(1/2))^(1/3))
!-       s = -(((36*epsi*gam*om2)-(72*ck2*gam)-(8*gam ** 3)+12.*sqrt((12*epsi ** 3*om2 ** 3-3*epsi ** 2*om2 ** 2*gam ** 2+36*epsi ** 2*om2 ** 2*ck2-60*epsi*om2*ck2*gam ** 2+12*ck2*gam ** 4+36*epsi*om2*ck2 ** 2+24*ck2 ** 2*gam ** 2+12*ck2 ** 3))) ** (1./3.)/12.)+(3.*((epsi*om2)/3.+(ck2)/3.-(gam ** 2)/0.9E1)*((36*epsi*gam*om2)-(72*ck2*gam)-(8*gam ** 3)+12.*sqrt((12*epsi ** 3*om2 ** 3-3*epsi ** 2*om2 ** 2*gam ** 2+36*epsi ** 2*om2 ** 2*ck2-60*epsi*om2*ck2*gam ** 2+12*ck2*gam ** 4+36*epsi*om2*ck2 ** 2+24*ck2 ** 2*gam ** 2+12*ck2 ** 3))) ** (-1./3.))-((gam)/3.)+cmplx(0, 1./2.)*sqrt(3.)*(((36*epsi*gam*om2)-(72*ck2*gam)-(8*gam ** 3)+12.*sqrt((12*epsi ** 3*om2 ** 3-3*epsi ** 2*om2 ** 2*gam ** 2+36*epsi ** 2*om2 ** 2*ck2-60*epsi*om2*ck2*gam ** 2+12*ck2*gam ** 4+36*epsi*om2*ck2 ** 2+24*ck2 ** 2*gam ** 2+12*ck2 ** 3))) ** (1./3.)/6.+6.*((epsi*om2)/3.+(ck2)/3.-(gam ** 2)/0.9E1)*((36*epsi*gam*om2)-(72*ck2*gam)-(8*gam ** 3)+12.*sqrt((12*epsi ** 3*om2 ** 3-3*epsi ** 2*om2 ** 2*gam ** 2+36*epsi ** 2*om2 ** 2*ck2-60*epsi*om2*ck2*gam ** 2+12*ck2*gam ** 4+36*epsi*om2*ck2 ** 2+24*ck2 ** 2*gam ** 2+12*ck2 ** 3))) ** (-1./3.))
!- 
!-       reS= real(s)
!-       imS= aimag(s)
!- 
!-       ! check root:
!-       det = cabs((s**2 + ck2)*(s**2 + gam*s) + om2*s**2*epsi)
!-       write(*,'("*OLD* evalDisp: eps,omegap=",2e12.4," gam,k=",2e12.4," |det|=",e12.4)') eps,omegap,gam,k,det
!-       write(*,'("*OLD* evalDisp: reS,imS=",2(1P,e22.15))') reS,imS
!-       return
!-       end   





      !===========================================================================================
      ! Evaluate the complex index of refraction m=(mr,mi) for scattering off dispersive dielectrics
      !     m = mr + i*mi = sqrt{ mu_2 eps_2 (1+ Chi_2(s) ) / mu_1 eps_1 (1+Chi_1(s) ) }
      !===========================================================================================
      subroutine evalComplexIndexOfRefraction(mu1,eps1,chi1r,chi1i, 
     & mu2,eps2,chi2r,chi2i, mr,mi)

      implicit none
      real mu1,eps1,chi1r,chi1i, mu2,eps2,chi2r,chi2i, mr,mi
      ! local variables
      complex*16 chi1,chi2,mc

      chi1 = cmplx(chi1r,chi1i)
      chi2 = cmplx(chi2r,chi2i)

      mc = csqrt( ( mu2*eps2*(1.+chi2) )/( mu1*eps1*(1.+chi1) ) )

      mr=real(mc)
      mi=aimag(mc)

      return
      end



! =========================================================
!    Normalize a vector to have length 1
! =========================================================

! =========================================================
!    Normalize a complex vector to have length 1
! =========================================================

      ! ============================================================================================
      ! Evaluate the exact solution for a plane material interface in 3D
      ! Input:
      !   sr + I*si : incident complex frequency
      !   chiSum1 = chiSum1r + chiSum1i :  suseptibility on left
      !   chiSum2 = chiSum2r + chiSum2i : 
      !   (ax,ay,az) : incident amplitude vector 
      !   eps1,eps2    : magnetic permeabilities 
      !   mu1,mu2    : magnetic permeabilities 
      !   (nx,ny,nz) : normal to the interface
      !   (kx,ky,kz) : incident wave number 
      ! Output:
      !   krr(3) + I*kri(3) : refelected wave vector (complex: ktx=ktxr+ I*ktxi) etc.
      !   ktr(3) + I*kti(3) : transmitted wave vector (complex: ktx=ktxr+ I*ktxi) etc.
      !   arr(3) + I ari(3) : complex amplitide vector for reflected wave 
      !   atr(3) + I ati(3) : complex amplitide vector for transmitted wave 
      ! ============================================================================================
      subroutine evalMaterialInterfaceSolution3d( sr,si, chiSum1r, 
     & chiSum1i, chiSum2r, chiSum2i, av, nv, kv,eps1,eps2,mu1,mu2, 
     & krr,kri, ktr,kti, arr,ari,atr,ati )

      implicit none

      real av(3), nv(3), kv(3), krr(3), kri(3), ktr(3), kti(3)
      real eps1,eps2,mu1,mu2,sr,si,chiSum1r, chiSum1i, chiSum2r, 
     & chiSum2i
      real arr(3),ari(3), atr(3),ati(3)

      ! local
      real ax,ay,az, nx,ny,nz, kx,ky,kz
      real ktxr,ktxi,ktyr,ktyi,ktzr,ktzi
      real qx,qy,qz, gx,gy,gz, hx,hy,hz, px,py,pz
      real kDotN,aDotK,aNorm,qDotA,pDotA
      real krx,kry,krz
      complex*16 qDotQ,qDotH,qDotM,qDotG
      complex*16 pDotQ,pDotH,pDotM,pDotG
      ! real Eperp,Epar, Tperp,Tpar

      complex*16 qDotKrCrossQ,qDotKrCrossH,qDotKtCrossQ,qDotKtCrossM,
     & qDotKCrossA
      complex*16 pDotKrCrossQ,pDotKrCrossH,pDotKtCrossQ,pDotKtCrossM,
     & pDotKCrossA

      complex*16  ktx,kty,ktz,mx,my,mz,chiSum1,chiSum2,resid1,resid2
      complex*16 er(3),et(3)
      real ktSq,c1,c2
      complex*16 s, ktn

      complex*16 a4(4,4),b(4),x1,x2,x3,x4, a4s(4,4)
      integer n,nrhs,ipvt(4),info,m
      real eps
      integer debug

      eps = 1.e-12 ! fix me -- 100*REAL_EPS
      debug=0

      aDotk = ax*kx+ay*ky+az*kz
      if( debug.gt.0 )then
        write(*,'("ENTERING evalMaterialInterfaceSolution3d")')
        write(*,'(" s =",2(1pe10.2))') sr,si
        write(*,'(" av =",3(1pe10.2))') av(1),av(2),av(3)
        write(*,'(" kv =",3(1pe10.2))') kv(1),kv(2),kv(3)
        write(*,'(" nv =",3(1pe10.2))') nv(1),nv(2),nv(3)

        write(*,'(" aDotK =",(1pe10.2)," (Should be zero)")') aDotK
      end if
      if( abs(aDotk) .gt. eps )then
        write(*,'("evalMatInterfaceSolution3d: ERROR a.k ~= 0, a.k=",
     & e10.2)') aDotk
      end if


      ax=av(1)
      ay=av(2)
      az=av(3)

      kx=kv(1)
      ky=kv(2)
      kz=kv(3)

      nx=nv(1)
      ny=nv(2)
      nz=nv(3)

      ! complex transmission wave number kt 
      ktxr=ktr(1)
      ktyr=ktr(2)
      ktzr=ktr(3)

      ktxi=kti(1)
      ktyi=kti(2)
      ktzi=kti(3)

      ktx = cmplx(ktxr,ktxi)
      kty = cmplx(ktyr,ktyi)
      ktz = cmplx(ktzr,ktzi)


      ! q =(n X k)/ | n X k |
      qx = ny*kz - nz*ky
      qy = nz*kx - nx*kz
      qz = nx*ky - ny*kx
        aNorm = sqrt( qx**2+qy**2+qz**2)
        if( aNorm > 0. )then
          qx=qx/aNorm
          qy=qy/aNorm
          qz=qz/aNorm
        else
          write(*,'(" evalMatInterfaceSolution3d:ERROR: q is zero!")')
          stop 1111
        end if

      ! g = (q X k) / |q X k |
      gx = qy*kz - qz*ky
      gy = qz*kx - qx*kz
      gz = qx*ky - qy*kx
        aNorm = sqrt( gx**2+gy**2+gz**2)
        if( aNorm > 0. )then
          gx=gx/aNorm
          gy=gy/aNorm
          gz=gz/aNorm
        else
          write(*,'(" evalMatInterfaceSolution3d:ERROR: g is zero!")')
          stop 1111
        end if

      ! kr : reflected wave number 
      kDotN = kx*nx+ky*ny+kz*nz
      krx = kx - 2.*kDotN*nx
      kry = ky - 2.*kDotN*ny
      krz = kz - 2.*kDotN*nz
      krr(1)=krx
      kri(1)=0.
      krr(2)=kry
      kri(2)=0.
      krr(3)=krz
      kri(3)=0.

      c1 = sqrt(1./(eps1*mu1))
      c2 = sqrt(1./(eps2*mu2))

      s = cmplx(sr,si)
      chiSum1 = cmplx(chiSum1r,chiSum1i)
      chiSum2 = cmplx(chiSum2r,chiSum2i)

      !  s^2 + c1^2( kx^2 + ky^2 + kz^2 ) + s^2 chiSum1 = 0
      !  s^2 + c1^2( kn^2 + kt1^2 + kt2^2 ) + s^2 chiSum1 = 0
      ! Tranmitted:
      !  s^2 + c2^2( ktn^2 + kt1^2 + kt2^2 ) + s^2 chiSum2 = 0
      !
      ktSq = kx**2 + ky**2 + kz**2 - kDotN**2 ! kt1^2 + kt2^2
      ! ktn = normal component of the transmitted wave vector
      ktn = csqrt( -ktSq - (s**2)*(1.+chiSum2)/(c2**2) )

      if( .false. .and. aimag(ktn).lt. 0. )then
        ! If Im(ktn)<0 this mode will grow expoenntially in space
        ! In this case take the other 
        ktn=-ktn
      end if

      ! transmitted wave vector: 
      ktx = kx + (ktn -kDotN)*nx
      kty = ky + (ktn -kDotN)*ny
      ktz = kz + (ktn -kDotN)*nz
      ktr(1)= real(ktx)
      kti(1)=aimag(ktx)
      ktr(2)= real(kty)
      kti(2)=aimag(kty)
      ktr(3)= real(ktz)
      kti(3)=aimag(ktz)


      if( debug.gt.0 )then

        write(*,'(" transmitted wave: ktn=(",1pe10.2,",",1pe10.2,") (
     & could choose negative?)")') real(ktn),aimag(ktn)

        ! Check dispersion relations
        resid1 = (s**2)*(1.+chiSum1) + (c1**2)*( kx**2  + ky**2  + kz**
     & 2  )
        resid2 = (s**2)*(1.+chiSum2) + (c2**2)*( ktx**2 + kty**2 + ktz*
     & *2 )
        write(*,'(" residual in dispersion relation 1 =",2e10.2)') 
     & resid1
        write(*,'(" residual in dispersion relation 2 =",2e10.2)') 
     & resid2
      end if

      ! h = q X kr / | |
      hx = qy*krz - qz*kry
      hy = qz*krx - qx*krz
      hz = qx*kry - qy*krx
        aNorm = sqrt( hx**2+hy**2+hz**2)
        if( aNorm > 0. )then
          hx=hx/aNorm
          hy=hy/aNorm
          hz=hz/aNorm
        else
          write(*,'(" evalMatInterfaceSolution3d:ERROR: h is zero!")')
          stop 1111
        end if

      ! m = q X kt / | |   (may be complex)
      mx = qy*ktz - qz*kty
      my = qz*ktx - qx*ktz
      mz = qx*kty - qy*ktx
        aNorm = sqrt( cabs(mx)**2+cabs(my)**2+cabs(mz)**2)
        if( aNorm > 0. )then
          mx=mx/aNorm
          my=my/aNorm
          mz=mz/aNorm
        else
          write(*,'(" evalMatInterfaceSolution3d:ERROR: m is zero!")')
          stop 1111
        end if

      if( debug.gt.0 )then
        aNorm= sqrt( conjg(mx)*mx + conjg(my)*my + conjg(mz)*mz )
        write(*,'(" norm(m) = ",1pe10.2)') aNorm
      end if

      ! p = n X q : another vector in tangent plane of the interface
      px = ny*qz - nz*qy
      py = nz*qx - nx*qz
      pz = nx*qy - ny*qx
        aNorm = sqrt( px**2+py**2+pz**2)
        if( aNorm > 0. )then
          px=px/aNorm
          py=py/aNorm
          pz=pz/aNorm
        else
          write(*,'(" evalMatInterfaceSolution3d:ERROR: p is zero!")')
          stop 1111
        end if

      if( debug.gt.0 )then
        write(*,'(" kr =",3(1pe10.2))') krx,kry,krz

        write(*,'(" krr =",3(1pe10.2))') krr(1),krr(2),krr(3)
        write(*,'(" kri =",3(1pe10.2))') kri(1),kri(2),kri(3)

        write(*,'(" ktr =",3(1pe10.2))') ktr(1),ktr(2),ktr(3)
        write(*,'(" kti =",3(1pe10.2))') kti(1),kti(2),kti(3)

        write(*,'(" q = n X kh = ",3(1pe10.2))') qx,qy,qz
        write(*,'(" g = q X kh = ",3(1pe10.2))') gx,gy,gz

        write(*,'(" h = n X kh = ",3(1pe10.2))') hx,hy,hz
        write(*,'(" p = n X kh = ",3(1pe10.2))') px,py,pz
        write(*,'(" m = q X kth= ",6(1pe10.2))') mx,my,mz
      end if


      ! Eperp = ax*qx+ay*qy+az*qz
      ! Epar  = ax*gx+ay*gy+az*gz
      qDotA = qx*ax+qy*ay+qz*az
      pDotA = px*ax+py*ay+pz*az

      ! ---- Evaluate the four transmission and reflection coefficients -----
      !   See notes...

      qDotQ = qx*qx+qy*qy+qz*qz
      qDotH = qx*hx+qy*hy+qz*hz
      qDotM = qx*mx+qy*my+qz*mz
      qDotG = qx*gx+qy*gy+qz*gz

      pDotQ = px*qx+py*qy+pz*qz
      pDotH = px*hx+py*hy+pz*hz
      pDotM = px*mx+py*my+pz*mz
      pDotG = px*gx+py*gy+pz*gz

      a4(1,1) = qDotQ
      a4(1,2) = qDotH
      a4(1,3) =-qDotQ
      a4(1,4) =-qDotM
      b(1)    =-qDotA

      a4(2,1) = pDotQ
      a4(2,2) = pDotH
      a4(2,3) =-pDotQ
      a4(2,4) =-pDotM
      b(2)    =-pDotA

      ! t.[  x1*(kr X q)/mu1 + x2*( kr X h )/mu1 - ( x3*(kt X q)/mu2 + x4*( kt X m )/mu2 ] = - t.[ (k X a)/mu1 ]

      qDotKrCrossQ = qx*( kry*qz-krz*qy  ) + qy*( krz*qx-krx*qz ) + qz*
     & ( krx*qy-kry*qx )
      qDotKrCrossH = qx*( kry*hz-krz*hy  ) + qy*( krz*hx-krx*hz ) + qz*
     & ( krx*hy-kry*hx )
      qDotKtCrossQ = qx*( kty*qz-ktz*qy  ) + qy*( ktz*qx-ktx*qz ) + qz*
     & ( ktx*qy-kty*qx )
      qDotKtCrossM = qx*( kty*mz-ktz*my  ) + qy*( ktz*mx-ktx*mz ) + qz*
     & ( ktx*my-kty*mx )
      qDotKCrossA  = qx*( ky*az-kz*ay    ) + qy*( kz*ax-kx*az   ) + qz*
     & ( kx*ay -ky*ax  )

      a4(3,1) = qDotKrCrossQ/mu1
      a4(3,2) = qDotKrCrossH/mu1
      a4(3,3) =-qDotKtCrossQ/mu2
      a4(3,4) =-qDotKtCrossM/mu2
      b(3)    =-qDotKCrossA/mu1

      pDotKrCrossQ = px*( kry*qz-krz*qy  ) + py*( krz*qx-krx*qz ) + pz*
     & ( krx*qy-kry*qx )
      pDotKrCrossH = px*( kry*hz-krz*hy  ) + py*( krz*hx-krx*hz ) + pz*
     & ( krx*hy-kry*hx )
      pDotKtCrossQ = px*( kty*qz-ktz*qy  ) + py*( ktz*qx-ktx*qz ) + pz*
     & ( ktx*qy-kty*qx )
      pDotKtCrossM = px*( kty*mz-ktz*my  ) + py*( ktz*mx-ktx*mz ) + pz*
     & ( ktx*my-kty*mx )
      pDotKCrossA  = px*( ky*az-kz*ay    ) + py*( kz*ax-kx*az   ) + pz*
     & ( kx*ay -ky*ax  )

      a4(4,1) = pDotKrCrossQ/mu1
      a4(4,2) = pDotKrCrossH/mu1
      a4(4,3) =-pDotKtCrossQ/mu2
      a4(4,4) =-pDotKtCrossM/mu2
      b(4)    =-pDotKCrossA/mu1

      if( debug.gt.0 )then
       write(*,'(" a= (",2(1pe10.2),",",2(1pe10.2),",",2(1pe10.2),",",
     & 2(1pe10.2),")")') a4(1,1),a4(1,2),a4(1,3),a4(1,4)
       write(*,'(" a= (",2(1pe10.2),",",2(1pe10.2),",",2(1pe10.2),",",
     & 2(1pe10.2),")")') a4(2,1),a4(2,2),a4(2,3),a4(2,4)
       write(*,'(" a= (",2(1pe10.2),",",2(1pe10.2),",",2(1pe10.2),",",
     & 2(1pe10.2),")")') a4(3,1),a4(3,2),a4(3,3),a4(3,4)
       write(*,'(" a= (",2(1pe10.2),",",2(1pe10.2),",",2(1pe10.2),",",
     & 2(1pe10.2),")")') a4(4,1),a4(4,2),a4(4,3),a4(4,4)

       write(*,'(" b1=(",1pe10.2,",",1pe10.2,")")') real(b(1)),aimag(b(
     & 1))
       write(*,'(" b2=(",1pe10.2,",",1pe10.2,")")') real(b(2)),aimag(b(
     & 2))
       write(*,'(" b3=(",1pe10.2,",",1pe10.2,")")') real(b(3)),aimag(b(
     & 3))
       write(*,'(" b4=(",1pe10.2,",",1pe10.2,")")') real(b(4)),aimag(b(
     & 4))
      end if

      ! save matrix for checking
      do n=1,4
        do m=1,4
          a4s(m,n)=a4(m,n)
        end do
      end do
      n=4
      nrhs=1

      call zgesv( n, nrhs, a4, n, ipvt,b, n,info )
      if( info.ne.0 )then
        write(*,'("Error return from zgesv: info=",i6)') info
        stop 3333
      end if

      if( debug.gt.0 )then
       write(*,'(" x1=(",1pe10.2,",",1pe10.2,")")') real(b(1)),aimag(b(
     & 1))
       write(*,'(" x2=(",1pe10.2,",",1pe10.2,")")') real(b(2)),aimag(b(
     & 2))
       write(*,'(" x3=(",1pe10.2,",",1pe10.2,")")') real(b(3)),aimag(b(
     & 3))
       write(*,'(" x4=(",1pe10.2,",",1pe10.2,")")') real(b(4)),aimag(b(
     & 4))
      end if

      x1 = b(1)
      x2 = b(2)
      x3 = b(3)
      x4 = b(4)

      ! Complex reflected wave coefficient vector
      er(1) = x1*qx + x2*hx
      er(2) = x1*qy + x2*hy
      er(3) = x1*qz + x2*hz

      ! Complex transmitted wave coefficient vector
      et(1) = x3*qx + x4*mx
      et(2) = x3*qy + x4*my
      et(3) = x3*qz + x4*mz

      ! return real and imag. parts or er and et
      do m=1,3
        arr(m)= real(er(m))
        ari(m)=aimag(er(m))
        atr(m)= real(et(m))
        ati(m)=aimag(et(m))
      end do

      if( debug.gt.0 )then
        ! checks
        ! tm.( a + er - et )
        resid1 = a4s(1,1)*x1 + a4s(1,2)*x2 + a4s(1,3)*x3 + a4s(1,4)*x4 
     & + qDotA
        write(*,'(" check: 1st eqn: resid=",2(1pe10.2))') resid1

        resid1 = a4s(2,1)*x1 + a4s(2,2)*x2 + a4s(2,3)*x3 + a4s(2,4)*x4 
     & + pDotA
        write(*,'(" check: 2nd eqn: resid=",2(1pe10.2))') resid1

        resid1 = qx*( ax+er(1)-et(1) ) + qy*( ay+er(2)-et(2) ) + qz*( 
     & az+er(3)-et(3) )
        write(*,'(" check: q.( a + er - et )=",2(1pe10.2))') resid1

        resid1 = px*( ax+er(1)-et(1) ) + py*( ay+er(2)-et(2) ) + pz*( 
     & az+er(3)-et(3) )
        write(*,'(" check: p.( a + er - et )=",2(1pe10.2))') resid1

        resid2 = px*er(1) + py*er(2) + pz*er(3) - (pDotQ*x1 +pDotH*x2)
        write(*,'(" check: resid2=",2(1pe10.2))') resid2
        resid2 = px*et(1) + py*et(2) + pz*et(3) - (pDotQ*x3 +pDotM*x4)
        write(*,'(" check: resid2=",2(1pe10.2))') resid2
      end if


      return
      end
