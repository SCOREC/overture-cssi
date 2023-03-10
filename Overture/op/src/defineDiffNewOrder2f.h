c Define statement functions for difference approximations of order 2 
c Thsi file was generated by op/src/makeIncludeNew.p  
c Arguments: u,rx,dr,dx: names for the grid function, jacobian, unit square spacing and rectangular grid spacing
c To include derivatives of rx use OPTION=RX
#beginMacro defineDifferenceNewOrder2Components0(u,rx,dr,dx,OPTION)

#If #OPTION == "RX"
dr ## 12(kd) = 1./(2.*dr(kd))
dr ## 22(kd) = 1./(dr(kd)**2)
#End

u ## r2(i1,i2,i3)=(u(i1+1,i2,i3)-u(i1-1,i2,i3))*dr ## 12(0)
u ## s2(i1,i2,i3)=(u(i1,i2+1,i3)-u(i1,i2-1,i3))*dr ## 12(1)
u ## t2(i1,i2,i3)=(u(i1,i2,i3+1)-u(i1,i2,i3-1))*dr ## 12(2)

u ## rr2(i1,i2,i3)=(-2.*u(i1,i2,i3)+(u(i1+1,i2,i3)+u(i1-1,i2,i3)) )*dr ## 22(0)
u ## ss2(i1,i2,i3)=(-2.*u(i1,i2,i3)+(u(i1,i2+1,i3)+u(i1,i2-1,i3)) )*dr ## 22(1)
u ## rs2(i1,i2,i3)=(u ## r2(i1,i2+1,i3)-u ## r2(i1,i2-1,i3))*dr ## 12(1)
u ## tt2(i1,i2,i3)=(-2.*u(i1,i2,i3)+(u(i1,i2,i3+1)+u(i1,i2,i3-1)) )*dr ## 22(2)
u ## rt2(i1,i2,i3)=(u ## r2(i1,i2,i3+1)-u ## r2(i1,i2,i3-1))*dr ## 12(2)
u ## st2(i1,i2,i3)=(u ## s2(i1,i2,i3+1)-u ## s2(i1,i2,i3-1))*dr ## 12(2)
u ## rrr2(i1,i2,i3)=(-2.*(u(i1+1,i2,i3)-u(i1-1,i2,i3))+(u(i1+2,i2,i3)-u(i1-2,i2,i3)) )*dr ## 22(0)*dr ## 12(0)
u ## sss2(i1,i2,i3)=(-2.*(u(i1,i2+1,i3)-u(i1,i2-1,i3))+(u(i1,i2+2,i3)-u(i1,i2-2,i3)) )*dr ## 22(1)*dr ## 12(1)
u ## ttt2(i1,i2,i3)=(-2.*(u(i1,i2,i3+1)-u(i1,i2,i3-1))+(u(i1,i2,i3+2)-u(i1,i2,i3-2)) )*dr ## 22(2)*dr ## 12(2)

#If #OPTION == "RX"
rx ##r2(i1,i2,i3,m,n)=(rx ##(i1+1,i2,i3,m,n)-rx ##(i1-1,i2,i3,m,n))*dr ## 12(0)
rx ##s2(i1,i2,i3,m,n)=(rx ##(i1,i2+1,i3,m,n)-rx ##(i1,i2-1,i3,m,n))*dr ## 12(1)
rx ##t2(i1,i2,i3,m,n)=(rx ##(i1,i2,i3+1,m,n)-rx ##(i1,i2,i3-1,m,n))*dr ## 12(2)
rx ##rr2(i1,i2,i3,m,n)=(-2.*rx ##(i1,i2,i3,m,n)+(rx ##(i1+1,i2,i3,m,n)+rx ##(i1-1,i2,i3,m,n)) )*dr ## 22(0)
rx ##ss2(i1,i2,i3,m,n)=(-2.*rx ##(i1,i2,i3,m,n)+(rx ##(i1,i2+1,i3,m,n)+rx ##(i1,i2-1,i3,m,n)) )*dr ## 22(1)
rx ##rs2(i1,i2,i3,m,n)=(rx ## r2(i1,i2+1,i3,m,n)-rx ## r2(i1,i2-1,i3,m,n))*dr ## 12(1)
#End

u ## x21(i1,i2,i3)= rx(i1,i2,i3,0,0)*u ## r2(i1,i2,i3)
u ## y21(i1,i2,i3)=0
u ## z21(i1,i2,i3)=0

u ## x22(i1,i2,i3)= rx(i1,i2,i3,0,0)*u ## r2(i1,i2,i3)+rx(i1,i2,i3,1,0)*u ## s2(i1,i2,i3)
u ## y22(i1,i2,i3)= rx(i1,i2,i3,0,1)*u ## r2(i1,i2,i3)+rx(i1,i2,i3,1,1)*u ## s2(i1,i2,i3)
u ## z22(i1,i2,i3)=0
u ## x23(i1,i2,i3)=rx(i1,i2,i3,0,0)*u ## r2(i1,i2,i3)+rx(i1,i2,i3,1,0)*u ## s2(i1,i2,i3)+rx(i1,i2,i3,2,0)*u ## t2(i1,i2,i3)
u ## y23(i1,i2,i3)=rx(i1,i2,i3,0,1)*u ## r2(i1,i2,i3)+rx(i1,i2,i3,1,1)*u ## s2(i1,i2,i3)+rx(i1,i2,i3,2,1)*u ## t2(i1,i2,i3)
u ## z23(i1,i2,i3)=rx(i1,i2,i3,0,2)*u ## r2(i1,i2,i3)+rx(i1,i2,i3,1,2)*u ## s2(i1,i2,i3)+rx(i1,i2,i3,2,2)*u ## t2(i1,i2,i3)

#If #OPTION == "RX"
rx ## x21(i1,i2,i3,m,n)= rx(i1,i2,i3,0,0)*rx ##r2(i1,i2,i3,m,n)
rx ##x22(i1,i2,i3,m,n)= rx(i1,i2,i3,0,0)*rx ##r2(i1,i2,i3,m,n)+rx(i1,i2,i3,1,0)*rx ##s2(i1,i2,i3,m,n)
rx ##y22(i1,i2,i3,m,n)= rx(i1,i2,i3,0,1)*rx ##r2(i1,i2,i3,m,n)+rx(i1,i2,i3,1,1)*rx ##s2(i1,i2,i3,m,n)
rx ##x23(i1,i2,i3,m,n)=rx(i1,i2,i3,0,0)*rx ##r2(i1,i2,i3,m,n)+rx(i1,i2,i3,1,0)*rx ##s2(i1,i2,i3,m,n)+rx(i1,i2,i3,2,0)*rx ##t2(i1,i2,i3,m,n)
rx ##y23(i1,i2,i3,m,n)=rx(i1,i2,i3,0,1)*rx ##r2(i1,i2,i3,m,n)+rx(i1,i2,i3,1,1)*rx ##s2(i1,i2,i3,m,n)+rx(i1,i2,i3,2,1)*rx ##t2(i1,i2,i3,m,n)
rx ##z23(i1,i2,i3,m,n)=rx(i1,i2,i3,0,2)*rx ##r2(i1,i2,i3,m,n)+rx(i1,i2,i3,1,2)*rx ##s2(i1,i2,i3,m,n)+rx(i1,i2,i3,2,2)*rx ##t2(i1,i2,i3,m,n)
#End

u ## xx21(i1,i2,i3)=(rx(i1,i2,i3,0,0)**2)*u ## rr2(i1,i2,i3)+(rx ## x2 ## 2(i1,i2,i3,0,0))*u ## r2(i1,i2,i3)
u ## yy21(i1,i2,i3)=0
u ## xy21(i1,i2,i3)=0
u ## xz21(i1,i2,i3)=0
u ## yz21(i1,i2,i3)=0
u ## zz21(i1,i2,i3)=0
u ## laplacian21(i1,i2,i3)=u ## xx21(i1,i2,i3)
u ## xx22(i1,i2,i3)=(rx(i1,i2,i3,0,0)**2)*u ## rr2(i1,i2,i3)+2.*(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,0))*u ## rs2(i1,i2,i3)+(rx(i1,i2,i3,1,0)**2)*u ## ss2(i1,i2,i3)+(rx ## x2 ## 2(i1,i2,i3,0,0))*u ## r2(i1,i2,i3)+(rx ## x2 ## 2(i1,i2,i3,1,0))*u ## s2(i1,i2,i3)
u ## yy22(i1,i2,i3)=(rx(i1,i2,i3,0,1)**2)*u ## rr2(i1,i2,i3)+2.*(rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,1))*u ## rs2(i1,i2,i3)+(rx(i1,i2,i3,1,1)**2)*u ## ss2(i1,i2,i3)+(rx ## y2 ## 2(i1,i2,i3,0,1))*u ## r2(i1,i2,i3)+(rx ## y2 ## 2(i1,i2,i3,1,1))*u ## s2(i1,i2,i3)
u ## xy22(i1,i2,i3)=rx(i1,i2,i3,0,0)*rx(i1,i2,i3,0,1)*u ## rr2(i1,i2,i3)+(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,1)+rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,0))*u ## rs2(i1,i2,i3)+rx(i1,i2,i3,1,0)*rx(i1,i2,i3,1,1)*u ## ss2(i1,i2,i3)+rx ## x2 ## 2(i1,i2,i3,0,1)*u ## r2(i1,i2,i3)+rx ## x2 ## 2(i1,i2,i3,1,1)*u ## s2(i1,i2,i3)
u ## xz22(i1,i2,i3)=0
u ## yz22(i1,i2,i3)=0
u ## zz22(i1,i2,i3)=0
u ## laplacian22(i1,i2,i3)=(rx(i1,i2,i3,0,0)**2+rx(i1,i2,i3,0,1)**2)*u ## rr2(i1,i2,i3)+2.*(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,0)+ rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,1))*u ## rs2(i1,i2,i3)+(rx(i1,i2,i3,1,0)**2+rx(i1,i2,i3,1,1)**2)*u ## ss2(i1,i2,i3)+(rx ## x2 ## 2(i1,i2,i3,0,0)+rx ## y2 ## 2(i1,i2,i3,0,1))*u ## r2(i1,i2,i3)+(rx ## x2 ## 2(i1,i2,i3,1,0)+rx ## y2 ## 2(i1,i2,i3,1,1))*u ## s2(i1,i2,i3)
u ## xx23(i1,i2,i3)=rx(i1,i2,i3,0,0)**2*u ## rr2(i1,i2,i3)+rx(i1,i2,i3,1,0)**2*u ## ss2(i1,i2,i3)+rx(i1,i2,i3,2,0)**2*u ## tt2(i1,i2,i3)+2.*rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,0)*u ## rs2(i1,i2,i3)+2.*rx(i1,i2,i3,0,0)*rx(i1,i2,i3,2,0)*u ## rt2(i1,i2,i3)+2.*rx(i1,i2,i3,1,0)*rx(i1,i2,i3,2,0)*u ## st2(i1,i2,i3)+rx ## x2 ## 3(i1,i2,i3,0,0)*u ## r2(i1,i2,i3)+rx ## x2 ## 3(i1,i2,i3,1,0)*u ## s2(i1,i2,i3)+rx ## x2 ## 3(i1,i2,i3,2,0)*u ## t2(i1,i2,i3)
u ## yy23(i1,i2,i3)=rx(i1,i2,i3,0,1)**2*u ## rr2(i1,i2,i3)+rx(i1,i2,i3,1,1)**2*u ## ss2(i1,i2,i3)+rx(i1,i2,i3,2,1)**2*u ## tt2(i1,i2,i3)+2.*rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,1)*u ## rs2(i1,i2,i3)+2.*rx(i1,i2,i3,0,1)*rx(i1,i2,i3,2,1)*u ## rt2(i1,i2,i3)+2.*rx(i1,i2,i3,1,1)*rx(i1,i2,i3,2,1)*u ## st2(i1,i2,i3)+rx ## y2 ## 3(i1,i2,i3,0,1)*u ## r2(i1,i2,i3)+rx ## y2 ## 3(i1,i2,i3,1,1)*u ## s2(i1,i2,i3)+rx ## y2 ## 3(i1,i2,i3,2,1)*u ## t2(i1,i2,i3)
u ## zz23(i1,i2,i3)=rx(i1,i2,i3,0,2)**2*u ## rr2(i1,i2,i3)+rx(i1,i2,i3,1,2)**2*u ## ss2(i1,i2,i3)+rx(i1,i2,i3,2,2)**2*u ## tt2(i1,i2,i3)+2.*rx(i1,i2,i3,0,2)*rx(i1,i2,i3,1,2)*u ## rs2(i1,i2,i3)+2.*rx(i1,i2,i3,0,2)*rx(i1,i2,i3,2,2)*u ## rt2(i1,i2,i3)+2.*rx(i1,i2,i3,1,2)*rx(i1,i2,i3,2,2)*u ## st2(i1,i2,i3)+rx ## z2 ## 3(i1,i2,i3,0,2)*u ## r2(i1,i2,i3)+rx ## z2 ## 3(i1,i2,i3,1,2)*u ## s2(i1,i2,i3)+rx ## z2 ## 3(i1,i2,i3,2,2)*u ## t2(i1,i2,i3)
u ## xy23(i1,i2,i3)=rx(i1,i2,i3,0,0)*rx(i1,i2,i3,0,1)*u ## rr2(i1,i2,i3)+rx(i1,i2,i3,1,0)*rx(i1,i2,i3,1,1)*u ## ss2(i1,i2,i3)+rx(i1,i2,i3,2,0)*rx(i1,i2,i3,2,1)*u ## tt2(i1,i2,i3)+(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,1)+rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,0))*u ## rs2(i1,i2,i3)+(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,2,1)+rx(i1,i2,i3,0,1)*rx(i1,i2,i3,2,0))*u ## rt2(i1,i2,i3)+(rx(i1,i2,i3,1,0)*rx(i1,i2,i3,2,1)+rx(i1,i2,i3,1,1)*rx(i1,i2,i3,2,0))*u ## st2(i1,i2,i3)+rx ## x2 ## 3(i1,i2,i3,0,1)*u ## r2(i1,i2,i3)+rx ## x2 ## 3(i1,i2,i3,1,1)*u ## s2(i1,i2,i3)+rx ## x2 ## 3(i1,i2,i3,2,1)*u ## t2(i1,i2,i3)
u ## xz23(i1,i2,i3)=rx(i1,i2,i3,0,0)*rx(i1,i2,i3,0,2)*u ## rr2(i1,i2,i3)+rx(i1,i2,i3,1,0)*rx(i1,i2,i3,1,2)*u ## ss2(i1,i2,i3)+rx(i1,i2,i3,2,0)*rx(i1,i2,i3,2,2)*u ## tt2(i1,i2,i3)+(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,2)+rx(i1,i2,i3,0,2)*rx(i1,i2,i3,1,0))*u ## rs2(i1,i2,i3)+(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,2,2)+rx(i1,i2,i3,0,2)*rx(i1,i2,i3,2,0))*u ## rt2(i1,i2,i3)+(rx(i1,i2,i3,1,0)*rx(i1,i2,i3,2,2)+rx(i1,i2,i3,1,2)*rx(i1,i2,i3,2,0))*u ## st2(i1,i2,i3)+rx ## x2 ## 3(i1,i2,i3,0,2)*u ## r2(i1,i2,i3)+rx ## x2 ## 3(i1,i2,i3,1,2)*u ## s2(i1,i2,i3)+rx ## x2 ## 3(i1,i2,i3,2,2)*u ## t2(i1,i2,i3)
u ## yz23(i1,i2,i3)=rx(i1,i2,i3,0,1)*rx(i1,i2,i3,0,2)*u ## rr2(i1,i2,i3)+rx(i1,i2,i3,1,1)*rx(i1,i2,i3,1,2)*u ## ss2(i1,i2,i3)+rx(i1,i2,i3,2,1)*rx(i1,i2,i3,2,2)*u ## tt2(i1,i2,i3)+(rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,2)+rx(i1,i2,i3,0,2)*rx(i1,i2,i3,1,1))*u ## rs2(i1,i2,i3)+(rx(i1,i2,i3,0,1)*rx(i1,i2,i3,2,2)+rx(i1,i2,i3,0,2)*rx(i1,i2,i3,2,1))*u ## rt2(i1,i2,i3)+(rx(i1,i2,i3,1,1)*rx(i1,i2,i3,2,2)+rx(i1,i2,i3,1,2)*rx(i1,i2,i3,2,1))*u ## st2(i1,i2,i3)+rx ## y2 ## 3(i1,i2,i3,0,2)*u ## r2(i1,i2,i3)+rx ## y2 ## 3(i1,i2,i3,1,2)*u ## s2(i1,i2,i3)+rx ## y2 ## 3(i1,i2,i3,2,2)*u ## t2(i1,i2,i3)
u ## laplacian23(i1,i2,i3)=(rx(i1,i2,i3,0,0)**2+rx(i1,i2,i3,0,1)**2+rx(i1,i2,i3,0,2)**2)*u ## rr2(i1,i2,i3)+(rx(i1,i2,i3,1,0)**2+rx(i1,i2,i3,1,1)**2+rx(i1,i2,i3,1,2)**2)*u ## ss2(i1,i2,i3)+(rx(i1,i2,i3,2,0)**2+rx(i1,i2,i3,2,1)**2+rx(i1,i2,i3,2,2)**2)*u ## tt2(i1,i2,i3)+2.*(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,0)+ rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,1)+rx(i1,i2,i3,0,2)*rx(i1,i2,i3,1,2))*u ## rs2(i1,i2,i3)+2.*(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,2,0)+ rx(i1,i2,i3,0,1)*rx(i1,i2,i3,2,1)+rx(i1,i2,i3,0,2)*rx(i1,i2,i3,2,2))*u ## rt2(i1,i2,i3)+2.*(rx(i1,i2,i3,1,0)*rx(i1,i2,i3,2,0)+ rx(i1,i2,i3,1,1)*rx(i1,i2,i3,2,1)+rx(i1,i2,i3,1,2)*rx(i1,i2,i3,2,2))*u ## st2(i1,i2,i3)+(rx ## x2 ## 3(i1,i2,i3,0,0)+rx ## y2 ## 3(i1,i2,i3,0,1)+rx ## z2 ## 3(i1,i2,i3,0,2))*u ## r2(i1,i2,i3)+(rx ## x2 ## 3(i1,i2,i3,1,0)+rx ## y2 ## 3(i1,i2,i3,1,1)+rx ## z2 ## 3(i1,i2,i3,1,2))*u ## s2(i1,i2,i3)+(rx ## x2 ## 3(i1,i2,i3,2,0)+rx ## y2 ## 3(i1,i2,i3,2,1)+rx ## z2 ## 3(i1,i2,i3,2,2))*u ## t2(i1,i2,i3)
c============================================================================================
c Define derivatives for a rectangular grid
c
c============================================================================================
#If #OPTION == "RX"
dx ## 12(kd) = 1./(2.*dx(kd))
dx ## 22(kd) = 1./(dx(kd)**2)
#End

u ## x23r(i1,i2,i3)=(u(i1+1,i2,i3)-u(i1-1,i2,i3))*dx ## 12(0)
u ## y23r(i1,i2,i3)=(u(i1,i2+1,i3)-u(i1,i2-1,i3))*dx ## 12(1)
u ## z23r(i1,i2,i3)=(u(i1,i2,i3+1)-u(i1,i2,i3-1))*dx ## 12(2)

u ## xx23r(i1,i2,i3)=(-2.*u(i1,i2,i3)+(u(i1+1,i2,i3)+u(i1-1,i2,i3)) )*dx ## 22(0)
u ## yy23r(i1,i2,i3)=(-2.*u(i1,i2,i3)+(u(i1,i2+1,i3)+u(i1,i2-1,i3)) )*dx ## 22(1)
u ## xy23r(i1,i2,i3)=(u ## x23r(i1,i2+1,i3)-u ## x23r(i1,i2-1,i3))*dx ## 12(1)
u ## zz23r(i1,i2,i3)=(-2.*u(i1,i2,i3)+(u(i1,i2,i3+1)+u(i1,i2,i3-1)) )*dx ## 22(2)
u ## xz23r(i1,i2,i3)=(u ## x23r(i1,i2,i3+1)-u ## x23r(i1,i2,i3-1))*dx ## 12(2)
u ## yz23r(i1,i2,i3)=(u ## y23r(i1,i2,i3+1)-u ## y23r(i1,i2,i3-1))*dx ## 12(2)

u ## x21r(i1,i2,i3)= u ## x23r(i1,i2,i3)
u ## y21r(i1,i2,i3)= u ## y23r(i1,i2,i3)
u ## z21r(i1,i2,i3)= u ## z23r(i1,i2,i3)
u ## xx21r(i1,i2,i3)= u ## xx23r(i1,i2,i3)
u ## yy21r(i1,i2,i3)= u ## yy23r(i1,i2,i3)
u ## zz21r(i1,i2,i3)= u ## zz23r(i1,i2,i3)
u ## xy21r(i1,i2,i3)= u ## xy23r(i1,i2,i3)
u ## xz21r(i1,i2,i3)= u ## xz23r(i1,i2,i3)
u ## yz21r(i1,i2,i3)= u ## yz23r(i1,i2,i3)
u ## laplacian21r(i1,i2,i3)=u ## xx23r(i1,i2,i3)
u ## x22r(i1,i2,i3)= u ## x23r(i1,i2,i3)
u ## y22r(i1,i2,i3)= u ## y23r(i1,i2,i3)
u ## z22r(i1,i2,i3)= u ## z23r(i1,i2,i3)
u ## xx22r(i1,i2,i3)= u ## xx23r(i1,i2,i3)
u ## yy22r(i1,i2,i3)= u ## yy23r(i1,i2,i3)
u ## zz22r(i1,i2,i3)= u ## zz23r(i1,i2,i3)
u ## xy22r(i1,i2,i3)= u ## xy23r(i1,i2,i3)
u ## xz22r(i1,i2,i3)= u ## xz23r(i1,i2,i3)
u ## yz22r(i1,i2,i3)= u ## yz23r(i1,i2,i3)
u ## laplacian22r(i1,i2,i3)=u ## xx23r(i1,i2,i3)+u ## yy23r(i1,i2,i3)
u ## laplacian23r(i1,i2,i3)=u ## xx23r(i1,i2,i3)+u ## yy23r(i1,i2,i3)+u ## zz23r(i1,i2,i3)
u ## xxx22r(i1,i2,i3)=(-2.*(u ## (i1+1,i2,i3)-u ## (i1-1,i2,i3))+(u ## (i1+2,i2,i3)-u ## (i1-2,i2,i3)) )*dx ## 22(0)*dx ## 12(0)
u ## yyy22r(i1,i2,i3)=(-2.*(u ## (i1,i2+1,i3)-u ## (i1,i2-1,i3))+(u ## (i1,i2+2,i3)-u ## (i1,i2-2,i3)) )*dx ## 22(1)*dx ## 12(1)
u ## xxy22r(i1,i2,i3)=( u ## xx22r(i1,i2+1,i3)-u ## xx22r(i1,i2-1,i3))/(2.*dx(1))
u ## xyy22r(i1,i2,i3)=( u ## yy22r(i1+1,i2,i3)-u ## yy22r(i1-1,i2,i3))/(2.*dx(0))
u ## xxxx22r(i1,i2,i3)=(6.*u ## (i1,i2,i3)-4.*(u ## (i1+1,i2,i3)+u ## (i1-1,i2,i3)) +(u ## (i1+2,i2,i3)+u ## (i1-2,i2,i3)) )/(dx(0)**4)
u ## yyyy22r(i1,i2,i3)=(6.*u ## (i1,i2,i3)-4.*(u ## (i1,i2+1,i3)+u ## (i1,i2-1,i3)) +(u ## (i1,i2+2,i3)+u ## (i1,i2-2,i3)) )/(dx(1)**4)
u ## xxyy22r(i1,i2,i3)=( 4.*u ## (i1,i2,i3)-2.*(u ## (i1+1,i2,i3)+u ## (i1-1,i2,i3)+u ## (i1,i2+1,i3)+u ## (i1,i2-1,i3))+   (u ## (i1+1,i2+1,i3)+u ## (i1-1,i2+1,i3)+u ## (i1+1,i2-1,i3)+u ## (i1-1,i2-1,i3)) )/(dx(0)**2*dx(1)**2)
u ## LapSq22r(i1,i2,i3)= ( 6.*u ## (i1,i2,i3)- 4.*(u ## (i1+1,i2,i3)+u ## (i1-1,i2,i3))+(u ## (i1+2,i2,i3)+u ## (i1-2,i2,i3)) )/(dx(0)**4)+( 6.*u ## (i1,i2,i3)-4.*(u ## (i1,i2+1,i3)+u ## (i1,i2-1,i3)) +(u ## (i1,i2+2,i3)+u ## (i1,i2-2,i3)) )/(dx(1)**4)+( 8.*u ## (i1,i2,i3)-4.*(u ## (i1+1,i2,i3)+u ## (i1-1,i2,i3)+u ## (i1,i2+1,i3)+u ## (i1,i2-1,i3))+2.*(u ## (i1+1,i2+1,i3)+u ## (i1-1,i2+1,i3)+u ## (i1+1,i2-1,i3)+u ## (i1-1,i2-1,i3)) )/(dx(0)**2*dx(1)**2)
u ## xxx23r(i1,i2,i3)=(-2.*(u ## (i1+1,i2,i3)-u ## (i1-1,i2,i3))+(u ## (i1+2,i2,i3)-u ## (i1-2,i2,i3)) )*dx ## 22(0)*dx ## 12(0)
u ## yyy23r(i1,i2,i3)=(-2.*(u ## (i1,i2+1,i3)-u ## (i1,i2-1,i3))+(u ## (i1,i2+2,i3)-u ## (i1,i2-2,i3)) )*dx ## 22(1)*dx ## 12(1)
u ## zzz23r(i1,i2,i3)=(-2.*(u ## (i1,i2,i3+1)-u ## (i1,i2,i3-1))+(u ## (i1,i2,i3+2)-u ## (i1,i2,i3-2)) )*dx ## 22(1)*dx ## 12(2)
u ## xxy23r(i1,i2,i3)=( u ## xx22r(i1,i2+1,i3)-u ## xx22r(i1,i2-1,i3))/(2.*dx(1))
u ## xyy23r(i1,i2,i3)=( u ## yy22r(i1+1,i2,i3)-u ## yy22r(i1-1,i2,i3))/(2.*dx(0))
u ## xxz23r(i1,i2,i3)=( u ## xx22r(i1,i2,i3+1)-u ## xx22r(i1,i2,i3-1))/(2.*dx(2))
u ## yyz23r(i1,i2,i3)=( u ## yy22r(i1,i2,i3+1)-u ## yy22r(i1,i2,i3-1))/(2.*dx(2))
u ## xzz23r(i1,i2,i3)=( u ## zz22r(i1+1,i2,i3)-u ## zz22r(i1-1,i2,i3))/(2.*dx(0))
u ## yzz23r(i1,i2,i3)=( u ## zz22r(i1,i2+1,i3)-u ## zz22r(i1,i2-1,i3))/(2.*dx(1))
u ## xxxx23r(i1,i2,i3)=(6.*u ## (i1,i2,i3)-4.*(u ## (i1+1,i2,i3)+u ## (i1-1,i2,i3))+(u ## (i1+2,i2,i3)+u ## (i1-2,i2,i3)) )/(dx(0)**4)
u ## yyyy23r(i1,i2,i3)=(6.*u ## (i1,i2,i3)-4.*(u ## (i1,i2+1,i3)+u ## (i1,i2-1,i3))+(u ## (i1,i2+2,i3)+u ## (i1,i2-2,i3)) )/(dx(1)**4)
u ## zzzz23r(i1,i2,i3)=(6.*u ## (i1,i2,i3)-4.*(u ## (i1,i2,i3+1)+u ## (i1,i2,i3-1))+(u ## (i1,i2,i3+2)+u ## (i1,i2,i3-2)) )/(dx(2)**4)
u ## xxyy23r(i1,i2,i3)=( 4.*u ## (i1,i2,i3)-2.*(u ## (i1+1,i2,i3)+u ## (i1-1,i2,i3)+u ## (i1,i2+1,i3)+u ## (i1,i2-1,i3))+   (u ## (i1+1,i2+1,i3)+u ## (i1-1,i2+1,i3)+u ## (i1+1,i2-1,i3)+u ## (i1-1,i2-1,i3)) )/(dx(0)**2*dx(1)**2)
u ## xxzz23r(i1,i2,i3)=( 4.*u ## (i1,i2,i3)-2.*(u ## (i1+1,i2,i3)+u ## (i1-1,i2,i3)+u ## (i1,i2,i3+1)+u ## (i1,i2,i3-1))+   (u ## (i1+1,i2,i3+1)+u ## (i1-1,i2,i3+1)+u ## (i1+1,i2,i3-1)+u ## (i1-1,i2,i3-1)) )/(dx(0)**2*dx(2)**2)
u ## yyzz23r(i1,i2,i3)=( 4.*u ## (i1,i2,i3)-2.*(u ## (i1,i2+1,i3)  +u ## (i1,i2-1,i3)+  u ## (i1,i2  ,i3+1)+u ## (i1,i2  ,i3-1))+   (u ## (i1,i2+1,i3+1)+u ## (i1,i2-1,i3+1)+u ## (i1,i2+1,i3-1)+u ## (i1,i2-1,i3-1)) )/(dx(1)**2*dx(2)**2)

#endMacro
c Arguments: u,rx,dr,dx: names for the grid function, jacobian, unit square spacing and rectangular grid spacing
c To include derivatives of rx use OPTION=RX
#beginMacro defineDifferenceNewOrder2Components1(u,rx,dr,dx,OPTION)

#If #OPTION == "RX"
dr ## 12(kd) = 1./(2.*dr(kd))
dr ## 22(kd) = 1./(dr(kd)**2)
#End

u ## r2(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*dr ## 12(0)
u ## s2(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*dr ## 12(1)
u ## t2(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*dr ## 12(2)

u ## rr2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)) )*dr ## 22(0)
u ## ss2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd)) )*dr ## 22(1)
u ## rs2(i1,i2,i3,kd)=(u ## r2(i1,i2+1,i3,kd)-u ## r2(i1,i2-1,i3,kd))*dr ## 12(1)
u ## tt2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd)) )*dr ## 22(2)
u ## rt2(i1,i2,i3,kd)=(u ## r2(i1,i2,i3+1,kd)-u ## r2(i1,i2,i3-1,kd))*dr ## 12(2)
u ## st2(i1,i2,i3,kd)=(u ## s2(i1,i2,i3+1,kd)-u ## s2(i1,i2,i3-1,kd))*dr ## 12(2)
u ## rrr2(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*dr ## 22(0)*dr ## 12(0)
u ## sss2(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*dr ## 22(1)*dr ## 12(1)
u ## ttt2(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*dr ## 22(2)*dr ## 12(2)

#If #OPTION == "RX"
rx ##r2(i1,i2,i3,m,n)=(rx ##(i1+1,i2,i3,m,n)-rx ##(i1-1,i2,i3,m,n))*dr ## 12(0)
rx ##s2(i1,i2,i3,m,n)=(rx ##(i1,i2+1,i3,m,n)-rx ##(i1,i2-1,i3,m,n))*dr ## 12(1)
rx ##t2(i1,i2,i3,m,n)=(rx ##(i1,i2,i3+1,m,n)-rx ##(i1,i2,i3-1,m,n))*dr ## 12(2)
rx ##rr2(i1,i2,i3,m,n)=(-2.*rx ##(i1,i2,i3,m,n)+(rx ##(i1+1,i2,i3,m,n)+rx ##(i1-1,i2,i3,m,n)) )*dr ## 22(0)
rx ##ss2(i1,i2,i3,m,n)=(-2.*rx ##(i1,i2,i3,m,n)+(rx ##(i1,i2+1,i3,m,n)+rx ##(i1,i2-1,i3,m,n)) )*dr ## 22(1)
rx ##rs2(i1,i2,i3,m,n)=(rx ## r2(i1,i2+1,i3,m,n)-rx ## r2(i1,i2-1,i3,m,n))*dr ## 12(1)
#End

u ## x21(i1,i2,i3,kd)= rx(i1,i2,i3,0,0)*u ## r2(i1,i2,i3,kd)
u ## y21(i1,i2,i3,kd)=0
u ## z21(i1,i2,i3,kd)=0

u ## x22(i1,i2,i3,kd)= rx(i1,i2,i3,0,0)*u ## r2(i1,i2,i3,kd)+rx(i1,i2,i3,1,0)*u ## s2(i1,i2,i3,kd)
u ## y22(i1,i2,i3,kd)= rx(i1,i2,i3,0,1)*u ## r2(i1,i2,i3,kd)+rx(i1,i2,i3,1,1)*u ## s2(i1,i2,i3,kd)
u ## z22(i1,i2,i3,kd)=0
u ## x23(i1,i2,i3,kd)=rx(i1,i2,i3,0,0)*u ## r2(i1,i2,i3,kd)+rx(i1,i2,i3,1,0)*u ## s2(i1,i2,i3,kd)+rx(i1,i2,i3,2,0)*u ## t2(i1,i2,i3,kd)
u ## y23(i1,i2,i3,kd)=rx(i1,i2,i3,0,1)*u ## r2(i1,i2,i3,kd)+rx(i1,i2,i3,1,1)*u ## s2(i1,i2,i3,kd)+rx(i1,i2,i3,2,1)*u ## t2(i1,i2,i3,kd)
u ## z23(i1,i2,i3,kd)=rx(i1,i2,i3,0,2)*u ## r2(i1,i2,i3,kd)+rx(i1,i2,i3,1,2)*u ## s2(i1,i2,i3,kd)+rx(i1,i2,i3,2,2)*u ## t2(i1,i2,i3,kd)

#If #OPTION == "RX"
rx ## x21(i1,i2,i3,m,n)= rx(i1,i2,i3,0,0)*rx ##r2(i1,i2,i3,m,n)
rx ##x22(i1,i2,i3,m,n)= rx(i1,i2,i3,0,0)*rx ##r2(i1,i2,i3,m,n)+rx(i1,i2,i3,1,0)*rx ##s2(i1,i2,i3,m,n)
rx ##y22(i1,i2,i3,m,n)= rx(i1,i2,i3,0,1)*rx ##r2(i1,i2,i3,m,n)+rx(i1,i2,i3,1,1)*rx ##s2(i1,i2,i3,m,n)
rx ##x23(i1,i2,i3,m,n)=rx(i1,i2,i3,0,0)*rx ##r2(i1,i2,i3,m,n)+rx(i1,i2,i3,1,0)*rx ##s2(i1,i2,i3,m,n)+rx(i1,i2,i3,2,0)*rx ##t2(i1,i2,i3,m,n)
rx ##y23(i1,i2,i3,m,n)=rx(i1,i2,i3,0,1)*rx ##r2(i1,i2,i3,m,n)+rx(i1,i2,i3,1,1)*rx ##s2(i1,i2,i3,m,n)+rx(i1,i2,i3,2,1)*rx ##t2(i1,i2,i3,m,n)
rx ##z23(i1,i2,i3,m,n)=rx(i1,i2,i3,0,2)*rx ##r2(i1,i2,i3,m,n)+rx(i1,i2,i3,1,2)*rx ##s2(i1,i2,i3,m,n)+rx(i1,i2,i3,2,2)*rx ##t2(i1,i2,i3,m,n)
#End

u ## xx21(i1,i2,i3,kd)=(rx(i1,i2,i3,0,0)**2)*u ## rr2(i1,i2,i3,kd)+(rx ## x2 ## 2(i1,i2,i3,0,0))*u ## r2(i1,i2,i3,kd)
u ## yy21(i1,i2,i3,kd)=0
u ## xy21(i1,i2,i3,kd)=0
u ## xz21(i1,i2,i3,kd)=0
u ## yz21(i1,i2,i3,kd)=0
u ## zz21(i1,i2,i3,kd)=0
u ## laplacian21(i1,i2,i3,kd)=u ## xx21(i1,i2,i3,kd)
u ## xx22(i1,i2,i3,kd)=(rx(i1,i2,i3,0,0)**2)*u ## rr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,0))*u ## rs2(i1,i2,i3,kd)+(rx(i1,i2,i3,1,0)**2)*u ## ss2(i1,i2,i3,kd)+(rx ## x2 ## 2(i1,i2,i3,0,0))*u ## r2(i1,i2,i3,kd)+(rx ## x2 ## 2(i1,i2,i3,1,0))*u ## s2(i1,i2,i3,kd)
u ## yy22(i1,i2,i3,kd)=(rx(i1,i2,i3,0,1)**2)*u ## rr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,1))*u ## rs2(i1,i2,i3,kd)+(rx(i1,i2,i3,1,1)**2)*u ## ss2(i1,i2,i3,kd)+(rx ## y2 ## 2(i1,i2,i3,0,1))*u ## r2(i1,i2,i3,kd)+(rx ## y2 ## 2(i1,i2,i3,1,1))*u ## s2(i1,i2,i3,kd)
u ## xy22(i1,i2,i3,kd)=rx(i1,i2,i3,0,0)*rx(i1,i2,i3,0,1)*u ## rr2(i1,i2,i3,kd)+(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,1)+rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,0))*u ## rs2(i1,i2,i3,kd)+rx(i1,i2,i3,1,0)*rx(i1,i2,i3,1,1)*u ## ss2(i1,i2,i3,kd)+rx ## x2 ## 2(i1,i2,i3,0,1)*u ## r2(i1,i2,i3,kd)+rx ## x2 ## 2(i1,i2,i3,1,1)*u ## s2(i1,i2,i3,kd)
u ## xz22(i1,i2,i3,kd)=0
u ## yz22(i1,i2,i3,kd)=0
u ## zz22(i1,i2,i3,kd)=0
u ## laplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3,0,0)**2+rx(i1,i2,i3,0,1)**2)*u ## rr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,0)+ rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,1))*u ## rs2(i1,i2,i3,kd)+(rx(i1,i2,i3,1,0)**2+rx(i1,i2,i3,1,1)**2)*u ## ss2(i1,i2,i3,kd)+(rx ## x2 ## 2(i1,i2,i3,0,0)+rx ## y2 ## 2(i1,i2,i3,0,1))*u ## r2(i1,i2,i3,kd)+(rx ## x2 ## 2(i1,i2,i3,1,0)+rx ## y2 ## 2(i1,i2,i3,1,1))*u ## s2(i1,i2,i3,kd)
u ## xx23(i1,i2,i3,kd)=rx(i1,i2,i3,0,0)**2*u ## rr2(i1,i2,i3,kd)+rx(i1,i2,i3,1,0)**2*u ## ss2(i1,i2,i3,kd)+rx(i1,i2,i3,2,0)**2*u ## tt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,0)*u ## rs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3,0,0)*rx(i1,i2,i3,2,0)*u ## rt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3,1,0)*rx(i1,i2,i3,2,0)*u ## st2(i1,i2,i3,kd)+rx ## x2 ## 3(i1,i2,i3,0,0)*u ## r2(i1,i2,i3,kd)+rx ## x2 ## 3(i1,i2,i3,1,0)*u ## s2(i1,i2,i3,kd)+rx ## x2 ## 3(i1,i2,i3,2,0)*u ## t2(i1,i2,i3,kd)
u ## yy23(i1,i2,i3,kd)=rx(i1,i2,i3,0,1)**2*u ## rr2(i1,i2,i3,kd)+rx(i1,i2,i3,1,1)**2*u ## ss2(i1,i2,i3,kd)+rx(i1,i2,i3,2,1)**2*u ## tt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,1)*u ## rs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3,0,1)*rx(i1,i2,i3,2,1)*u ## rt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3,1,1)*rx(i1,i2,i3,2,1)*u ## st2(i1,i2,i3,kd)+rx ## y2 ## 3(i1,i2,i3,0,1)*u ## r2(i1,i2,i3,kd)+rx ## y2 ## 3(i1,i2,i3,1,1)*u ## s2(i1,i2,i3,kd)+rx ## y2 ## 3(i1,i2,i3,2,1)*u ## t2(i1,i2,i3,kd)
u ## zz23(i1,i2,i3,kd)=rx(i1,i2,i3,0,2)**2*u ## rr2(i1,i2,i3,kd)+rx(i1,i2,i3,1,2)**2*u ## ss2(i1,i2,i3,kd)+rx(i1,i2,i3,2,2)**2*u ## tt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3,0,2)*rx(i1,i2,i3,1,2)*u ## rs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3,0,2)*rx(i1,i2,i3,2,2)*u ## rt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3,1,2)*rx(i1,i2,i3,2,2)*u ## st2(i1,i2,i3,kd)+rx ## z2 ## 3(i1,i2,i3,0,2)*u ## r2(i1,i2,i3,kd)+rx ## z2 ## 3(i1,i2,i3,1,2)*u ## s2(i1,i2,i3,kd)+rx ## z2 ## 3(i1,i2,i3,2,2)*u ## t2(i1,i2,i3,kd)
u ## xy23(i1,i2,i3,kd)=rx(i1,i2,i3,0,0)*rx(i1,i2,i3,0,1)*u ## rr2(i1,i2,i3,kd)+rx(i1,i2,i3,1,0)*rx(i1,i2,i3,1,1)*u ## ss2(i1,i2,i3,kd)+rx(i1,i2,i3,2,0)*rx(i1,i2,i3,2,1)*u ## tt2(i1,i2,i3,kd)+(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,1)+rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,0))*u ## rs2(i1,i2,i3,kd)+(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,2,1)+rx(i1,i2,i3,0,1)*rx(i1,i2,i3,2,0))*u ## rt2(i1,i2,i3,kd)+(rx(i1,i2,i3,1,0)*rx(i1,i2,i3,2,1)+rx(i1,i2,i3,1,1)*rx(i1,i2,i3,2,0))*u ## st2(i1,i2,i3,kd)+rx ## x2 ## 3(i1,i2,i3,0,1)*u ## r2(i1,i2,i3,kd)+rx ## x2 ## 3(i1,i2,i3,1,1)*u ## s2(i1,i2,i3,kd)+rx ## x2 ## 3(i1,i2,i3,2,1)*u ## t2(i1,i2,i3,kd)
u ## xz23(i1,i2,i3,kd)=rx(i1,i2,i3,0,0)*rx(i1,i2,i3,0,2)*u ## rr2(i1,i2,i3,kd)+rx(i1,i2,i3,1,0)*rx(i1,i2,i3,1,2)*u ## ss2(i1,i2,i3,kd)+rx(i1,i2,i3,2,0)*rx(i1,i2,i3,2,2)*u ## tt2(i1,i2,i3,kd)+(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,2)+rx(i1,i2,i3,0,2)*rx(i1,i2,i3,1,0))*u ## rs2(i1,i2,i3,kd)+(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,2,2)+rx(i1,i2,i3,0,2)*rx(i1,i2,i3,2,0))*u ## rt2(i1,i2,i3,kd)+(rx(i1,i2,i3,1,0)*rx(i1,i2,i3,2,2)+rx(i1,i2,i3,1,2)*rx(i1,i2,i3,2,0))*u ## st2(i1,i2,i3,kd)+rx ## x2 ## 3(i1,i2,i3,0,2)*u ## r2(i1,i2,i3,kd)+rx ## x2 ## 3(i1,i2,i3,1,2)*u ## s2(i1,i2,i3,kd)+rx ## x2 ## 3(i1,i2,i3,2,2)*u ## t2(i1,i2,i3,kd)
u ## yz23(i1,i2,i3,kd)=rx(i1,i2,i3,0,1)*rx(i1,i2,i3,0,2)*u ## rr2(i1,i2,i3,kd)+rx(i1,i2,i3,1,1)*rx(i1,i2,i3,1,2)*u ## ss2(i1,i2,i3,kd)+rx(i1,i2,i3,2,1)*rx(i1,i2,i3,2,2)*u ## tt2(i1,i2,i3,kd)+(rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,2)+rx(i1,i2,i3,0,2)*rx(i1,i2,i3,1,1))*u ## rs2(i1,i2,i3,kd)+(rx(i1,i2,i3,0,1)*rx(i1,i2,i3,2,2)+rx(i1,i2,i3,0,2)*rx(i1,i2,i3,2,1))*u ## rt2(i1,i2,i3,kd)+(rx(i1,i2,i3,1,1)*rx(i1,i2,i3,2,2)+rx(i1,i2,i3,1,2)*rx(i1,i2,i3,2,1))*u ## st2(i1,i2,i3,kd)+rx ## y2 ## 3(i1,i2,i3,0,2)*u ## r2(i1,i2,i3,kd)+rx ## y2 ## 3(i1,i2,i3,1,2)*u ## s2(i1,i2,i3,kd)+rx ## y2 ## 3(i1,i2,i3,2,2)*u ## t2(i1,i2,i3,kd)
u ## laplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3,0,0)**2+rx(i1,i2,i3,0,1)**2+rx(i1,i2,i3,0,2)**2)*u ## rr2(i1,i2,i3,kd)+(rx(i1,i2,i3,1,0)**2+rx(i1,i2,i3,1,1)**2+rx(i1,i2,i3,1,2)**2)*u ## ss2(i1,i2,i3,kd)+(rx(i1,i2,i3,2,0)**2+rx(i1,i2,i3,2,1)**2+rx(i1,i2,i3,2,2)**2)*u ## tt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,1,0)+ rx(i1,i2,i3,0,1)*rx(i1,i2,i3,1,1)+rx(i1,i2,i3,0,2)*rx(i1,i2,i3,1,2))*u ## rs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3,0,0)*rx(i1,i2,i3,2,0)+ rx(i1,i2,i3,0,1)*rx(i1,i2,i3,2,1)+rx(i1,i2,i3,0,2)*rx(i1,i2,i3,2,2))*u ## rt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3,1,0)*rx(i1,i2,i3,2,0)+ rx(i1,i2,i3,1,1)*rx(i1,i2,i3,2,1)+rx(i1,i2,i3,1,2)*rx(i1,i2,i3,2,2))*u ## st2(i1,i2,i3,kd)+(rx ## x2 ## 3(i1,i2,i3,0,0)+rx ## y2 ## 3(i1,i2,i3,0,1)+rx ## z2 ## 3(i1,i2,i3,0,2))*u ## r2(i1,i2,i3,kd)+(rx ## x2 ## 3(i1,i2,i3,1,0)+rx ## y2 ## 3(i1,i2,i3,1,1)+rx ## z2 ## 3(i1,i2,i3,1,2))*u ## s2(i1,i2,i3,kd)+(rx ## x2 ## 3(i1,i2,i3,2,0)+rx ## y2 ## 3(i1,i2,i3,2,1)+rx ## z2 ## 3(i1,i2,i3,2,2))*u ## t2(i1,i2,i3,kd)
c============================================================================================
c Define derivatives for a rectangular grid
c
c============================================================================================
#If #OPTION == "RX"
dx ## 12(kd) = 1./(2.*dx(kd))
dx ## 22(kd) = 1./(dx(kd)**2)
#End

u ## x23r(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*dx ## 12(0)
u ## y23r(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*dx ## 12(1)
u ## z23r(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*dx ## 12(2)

u ## xx23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)) )*dx ## 22(0)
u ## yy23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd)) )*dx ## 22(1)
u ## xy23r(i1,i2,i3,kd)=(u ## x23r(i1,i2+1,i3,kd)-u ## x23r(i1,i2-1,i3,kd))*dx ## 12(1)
u ## zz23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd)) )*dx ## 22(2)
u ## xz23r(i1,i2,i3,kd)=(u ## x23r(i1,i2,i3+1,kd)-u ## x23r(i1,i2,i3-1,kd))*dx ## 12(2)
u ## yz23r(i1,i2,i3,kd)=(u ## y23r(i1,i2,i3+1,kd)-u ## y23r(i1,i2,i3-1,kd))*dx ## 12(2)

u ## x21r(i1,i2,i3,kd)= u ## x23r(i1,i2,i3,kd)
u ## y21r(i1,i2,i3,kd)= u ## y23r(i1,i2,i3,kd)
u ## z21r(i1,i2,i3,kd)= u ## z23r(i1,i2,i3,kd)
u ## xx21r(i1,i2,i3,kd)= u ## xx23r(i1,i2,i3,kd)
u ## yy21r(i1,i2,i3,kd)= u ## yy23r(i1,i2,i3,kd)
u ## zz21r(i1,i2,i3,kd)= u ## zz23r(i1,i2,i3,kd)
u ## xy21r(i1,i2,i3,kd)= u ## xy23r(i1,i2,i3,kd)
u ## xz21r(i1,i2,i3,kd)= u ## xz23r(i1,i2,i3,kd)
u ## yz21r(i1,i2,i3,kd)= u ## yz23r(i1,i2,i3,kd)
u ## laplacian21r(i1,i2,i3,kd)=u ## xx23r(i1,i2,i3,kd)
u ## x22r(i1,i2,i3,kd)= u ## x23r(i1,i2,i3,kd)
u ## y22r(i1,i2,i3,kd)= u ## y23r(i1,i2,i3,kd)
u ## z22r(i1,i2,i3,kd)= u ## z23r(i1,i2,i3,kd)
u ## xx22r(i1,i2,i3,kd)= u ## xx23r(i1,i2,i3,kd)
u ## yy22r(i1,i2,i3,kd)= u ## yy23r(i1,i2,i3,kd)
u ## zz22r(i1,i2,i3,kd)= u ## zz23r(i1,i2,i3,kd)
u ## xy22r(i1,i2,i3,kd)= u ## xy23r(i1,i2,i3,kd)
u ## xz22r(i1,i2,i3,kd)= u ## xz23r(i1,i2,i3,kd)
u ## yz22r(i1,i2,i3,kd)= u ## yz23r(i1,i2,i3,kd)
u ## laplacian22r(i1,i2,i3,kd)=u ## xx23r(i1,i2,i3,kd)+u ## yy23r(i1,i2,i3,kd)
u ## laplacian23r(i1,i2,i3,kd)=u ## xx23r(i1,i2,i3,kd)+u ## yy23r(i1,i2,i3,kd)+u ## zz23r(i1,i2,i3,kd)
u ## xxx22r(i1,i2,i3,kd)=(-2.*(u ## (i1+1,i2,i3,kd)-u ## (i1-1,i2,i3,kd))+(u ## (i1+2,i2,i3,kd)-u ## (i1-2,i2,i3,kd)) )*dx ## 22(0)*dx ## 12(0)
u ## yyy22r(i1,i2,i3,kd)=(-2.*(u ## (i1,i2+1,i3,kd)-u ## (i1,i2-1,i3,kd))+(u ## (i1,i2+2,i3,kd)-u ## (i1,i2-2,i3,kd)) )*dx ## 22(1)*dx ## 12(1)
u ## xxy22r(i1,i2,i3,kd)=( u ## xx22r(i1,i2+1,i3,kd)-u ## xx22r(i1,i2-1,i3,kd))/(2.*dx(1))
u ## xyy22r(i1,i2,i3,kd)=( u ## yy22r(i1+1,i2,i3,kd)-u ## yy22r(i1-1,i2,i3,kd))/(2.*dx(0))
u ## xxxx22r(i1,i2,i3,kd)=(6.*u ## (i1,i2,i3,kd)-4.*(u ## (i1+1,i2,i3,kd)+u ## (i1-1,i2,i3,kd)) +(u ## (i1+2,i2,i3,kd)+u ## (i1-2,i2,i3,kd)) )/(dx(0)**4)
u ## yyyy22r(i1,i2,i3,kd)=(6.*u ## (i1,i2,i3,kd)-4.*(u ## (i1,i2+1,i3,kd)+u ## (i1,i2-1,i3,kd)) +(u ## (i1,i2+2,i3,kd)+u ## (i1,i2-2,i3,kd)) )/(dx(1)**4)
u ## xxyy22r(i1,i2,i3,kd)=( 4.*u ## (i1,i2,i3,kd)-2.*(u ## (i1+1,i2,i3,kd)+u ## (i1-1,i2,i3,kd)+u ## (i1,i2+1,i3,kd)+u ## (i1,i2-1,i3,kd))+   (u ## (i1+1,i2+1,i3,kd)+u ## (i1-1,i2+1,i3,kd)+u ## (i1+1,i2-1,i3,kd)+u ## (i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
u ## LapSq22r(i1,i2,i3,kd)= ( 6.*u ## (i1,i2,i3,kd)- 4.*(u ## (i1+1,i2,i3,kd)+u ## (i1-1,i2,i3,kd))+(u ## (i1+2,i2,i3,kd)+u ## (i1-2,i2,i3,kd)) )/(dx(0)**4)+( 6.*u ## (i1,i2,i3,kd)-4.*(u ## (i1,i2+1,i3,kd)+u ## (i1,i2-1,i3,kd)) +(u ## (i1,i2+2,i3,kd)+u ## (i1,i2-2,i3,kd)) )/(dx(1)**4)+( 8.*u ## (i1,i2,i3,kd)-4.*(u ## (i1+1,i2,i3,kd)+u ## (i1-1,i2,i3,kd)+u ## (i1,i2+1,i3,kd)+u ## (i1,i2-1,i3,kd))+2.*(u ## (i1+1,i2+1,i3,kd)+u ## (i1-1,i2+1,i3,kd)+u ## (i1+1,i2-1,i3,kd)+u ## (i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
u ## xxx23r(i1,i2,i3,kd)=(-2.*(u ## (i1+1,i2,i3,kd)-u ## (i1-1,i2,i3,kd))+(u ## (i1+2,i2,i3,kd)-u ## (i1-2,i2,i3,kd)) )*dx ## 22(0)*dx ## 12(0)
u ## yyy23r(i1,i2,i3,kd)=(-2.*(u ## (i1,i2+1,i3,kd)-u ## (i1,i2-1,i3,kd))+(u ## (i1,i2+2,i3,kd)-u ## (i1,i2-2,i3,kd)) )*dx ## 22(1)*dx ## 12(1)
u ## zzz23r(i1,i2,i3,kd)=(-2.*(u ## (i1,i2,i3+1,kd)-u ## (i1,i2,i3-1,kd))+(u ## (i1,i2,i3+2,kd)-u ## (i1,i2,i3-2,kd)) )*dx ## 22(1)*dx ## 12(2)
u ## xxy23r(i1,i2,i3,kd)=( u ## xx22r(i1,i2+1,i3,kd)-u ## xx22r(i1,i2-1,i3,kd))/(2.*dx(1))
u ## xyy23r(i1,i2,i3,kd)=( u ## yy22r(i1+1,i2,i3,kd)-u ## yy22r(i1-1,i2,i3,kd))/(2.*dx(0))
u ## xxz23r(i1,i2,i3,kd)=( u ## xx22r(i1,i2,i3+1,kd)-u ## xx22r(i1,i2,i3-1,kd))/(2.*dx(2))
u ## yyz23r(i1,i2,i3,kd)=( u ## yy22r(i1,i2,i3+1,kd)-u ## yy22r(i1,i2,i3-1,kd))/(2.*dx(2))
u ## xzz23r(i1,i2,i3,kd)=( u ## zz22r(i1+1,i2,i3,kd)-u ## zz22r(i1-1,i2,i3,kd))/(2.*dx(0))
u ## yzz23r(i1,i2,i3,kd)=( u ## zz22r(i1,i2+1,i3,kd)-u ## zz22r(i1,i2-1,i3,kd))/(2.*dx(1))
u ## xxxx23r(i1,i2,i3,kd)=(6.*u ## (i1,i2,i3,kd)-4.*(u ## (i1+1,i2,i3,kd)+u ## (i1-1,i2,i3,kd))+(u ## (i1+2,i2,i3,kd)+u ## (i1-2,i2,i3,kd)) )/(dx(0)**4)
u ## yyyy23r(i1,i2,i3,kd)=(6.*u ## (i1,i2,i3,kd)-4.*(u ## (i1,i2+1,i3,kd)+u ## (i1,i2-1,i3,kd))+(u ## (i1,i2+2,i3,kd)+u ## (i1,i2-2,i3,kd)) )/(dx(1)**4)
u ## zzzz23r(i1,i2,i3,kd)=(6.*u ## (i1,i2,i3,kd)-4.*(u ## (i1,i2,i3+1,kd)+u ## (i1,i2,i3-1,kd))+(u ## (i1,i2,i3+2,kd)+u ## (i1,i2,i3-2,kd)) )/(dx(2)**4)
u ## xxyy23r(i1,i2,i3,kd)=( 4.*u ## (i1,i2,i3,kd)-2.*(u ## (i1+1,i2,i3,kd)+u ## (i1-1,i2,i3,kd)+u ## (i1,i2+1,i3,kd)+u ## (i1,i2-1,i3,kd))+   (u ## (i1+1,i2+1,i3,kd)+u ## (i1-1,i2+1,i3,kd)+u ## (i1+1,i2-1,i3,kd)+u ## (i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
u ## xxzz23r(i1,i2,i3,kd)=( 4.*u ## (i1,i2,i3,kd)-2.*(u ## (i1+1,i2,i3,kd)+u ## (i1-1,i2,i3,kd)+u ## (i1,i2,i3+1,kd)+u ## (i1,i2,i3-1,kd))+   (u ## (i1+1,i2,i3+1,kd)+u ## (i1-1,i2,i3+1,kd)+u ## (i1+1,i2,i3-1,kd)+u ## (i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
u ## yyzz23r(i1,i2,i3,kd)=( 4.*u ## (i1,i2,i3,kd)-2.*(u ## (i1,i2+1,i3,kd)  +u ## (i1,i2-1,i3,kd)+  u ## (i1,i2  ,i3+1,kd)+u ## (i1,i2  ,i3-1,kd))+   (u ## (i1,i2+1,i3+1,kd)+u ## (i1,i2-1,i3+1,kd)+u ## (i1,i2+1,i3-1,kd)+u ## (i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)

#endMacro


#beginMacro declareDifferenceNewOrder2(u,rx,dr,dx,OPTION)
real dr ##12
real dr ##22
real u ##r2
real u ##s2
real u ##t2
real u ##rr2
real u ##ss2
real u ##rs2
real u ##tt2
real u ##rt2
real u ##st2
real u ##rrr2
real u ##sss2
real u ##ttt2
real rx ##r2
real rx ##s2
real rx ##t2
real rx ##rr2
real rx ##ss2
real rx ##rs2
real u ##x21
real u ##y21
real u ##z21
real u ##x22
real u ##y22
real u ##z22
real u ##x23
real u ##y23
real u ##z23
real rx ##x21
real rx ##x22
real rx ##y22
real rx ##x23
real rx ##y23
real rx ##z23
real u ##xx21
real u ##yy21
real u ##xy21
real u ##xz21
real u ##yz21
real u ##zz21
real u ##laplacian21
real u ##xx22
real u ##yy22
real u ##xy22
real u ##xz22
real u ##yz22
real u ##zz22
real u ##laplacian22
real u ##xx23
real u ##yy23
real u ##zz23
real u ##xy23
real u ##xz23
real u ##yz23
real u ##laplacian23
real dx ##12
real dx ##22
real u ##x23r
real u ##y23r
real u ##z23r
real u ##xx23r
real u ##yy23r
real u ##xy23r
real u ##zz23r
real u ##xz23r
real u ##yz23r
real u ##x21r
real u ##y21r
real u ##z21r
real u ##xx21r
real u ##yy21r
real u ##zz21r
real u ##xy21r
real u ##xz21r
real u ##yz21r
real u ##laplacian21r
real u ##x22r
real u ##y22r
real u ##z22r
real u ##xx22r
real u ##yy22r
real u ##zz22r
real u ##xy22r
real u ##xz22r
real u ##yz22r
real u ##laplacian22r
real u ##laplacian23r
real u ##xxx22r
real u ##yyy22r
real u ##xxy22r
real u ##xyy22r
real u ##xxxx22r
real u ##yyyy22r
real u ##xxyy22r
real u ##LapSq22r
real u ##xxx23r
real u ##yyy23r
real u ##zzz23r
real u ##xxy23r
real u ##xyy23r
real u ##xxz23r
real u ##yyz23r
real u ##xzz23r
real u ##yzz23r
real u ##xxxx23r
real u ##yyyy23r
real u ##zzzz23r
real u ##xxyy23r
real u ##xxzz23r
real u ##yyzz23r
#endMacro
