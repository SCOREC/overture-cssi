// This file automatically generated from champBoundaryConditions.bC with bpp.
#include "Cgad.h"
#include "CompositeGridOperators.h"
#include "ParallelUtility.h"
//#include "AdamsPCData.h"
//#include "AdvanceOptions.h"
//#include "OGPulseFunction.h"
//#include "OGTrigFunction.h"
#include "Oges.h"

#include "SparseRep.h"

#include "gridFunctionNorms.h"
#include "Interface.h"  

// forward declarations:

int getChampParameters( int grid, int side, int axis, int grid2, int side2, int axis2, 
                                                real dt, real dx, Parameters & parameters, RealArray & champParameters );

// Macros to define 2nd-order difference approximations: 
#include "cgux2a.h"
//#include "champ4InterfaceStencilRectangularMacro.h" 

// // static real S = 1.*pow(10.,-2); // *note* also defined in AdParameters.C and userDefinedInitialConditions.C
// // static real G = 0.0;
// // #define PI 3.14159

// static real dCoupleBC=0.;
// static real cDecouple=1.;
// //static real cDecouple=0.;
// static real flux=1.;
// // static real h0= 1.5;
// // static real he = 0.;

// static real steady=0.;
// static real cpressure=-500.;


// static real tabserr = pow(10.,-8); //would be nice to pass this parameter
// static real trelerr = pow(10.,-10);

// static Ogshow *pshow=NULL;

    

#define ForBoundary(side,axis)   for( axis=0; axis<mg.numberOfDimensions(); axis++ ) for( side=0; side<=1; side++ )

#define FOR_3IJD(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); int J1Base =J1.getBase(),   J2Base =J2.getBase(),  J3Base =J3.getBase();  for(int i3=I3Base,j3=J3Base; i3<=I3Bound; i3++,j3++) for(int i2=I2Base,j2=J2Base; i2<=I2Bound; i2++,j2++) for(int i1=I1Base,j1=J1Base; i1<=I1Bound; i1++,j1++)

#define FOR_3IJ(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3) I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); J1Base =J1.getBase(),   J2Base =J2.getBase(),  J3Base =J3.getBase();  for(int i3=I3Base,j3=J3Base; i3<=I3Bound; i3++,j3++) for(int i2=I2Base,j2=J2Base; i2<=I2Bound; i2++,j2++) for(int i1=I1Base,j1=J1Base; i1<=I1Bound; i1++,j1++)
// Use this for indexing into coefficient matrices representing systems of equations
#define CE(c,e) (stencilSize*((c)+numberOfComponentsForCoefficients*(e)))
#define M123(m1,m2,m3) (m1+halfWidth1+width*(m2+halfWidth2+width*(m3+halfWidth3)))
#define M123CE(m1,m2,m3,c,e) (M123(m1,m2,m3)+CE(c,e))


#define ForStencil(m1,m2,m3)   for( m3=-halfWidth3; m3<=halfWidth3; m3++) for( m2=-halfWidth2; m2<=halfWidth2; m2++) for( m1=-halfWidth1; m1<=halfWidth1; m1++)

#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(int i3=I3Base; i3<=I3Bound; i3++) for(int i2=I2Base; i2<=I2Bound; i2++) for(int i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(int i3=I3Base; i3<=I3Bound; i3++) for(int i2=I2Base; i2<=I2Bound; i2++) for(int i1=I1Base; i1<=I1Bound; i1++)

// =======================================================================
// indexToEquation( n,i1,i2,i3 ) : defines the global index for each unknown in the system
//     n=component number (uc,vc,...)
//    (i1,i2,i3) = grid point
// =======================================================================
#define indexToEquation( n,i1,i2,i3 ) (n+1+ numberOfComponentsForCoefficients*(i1-equationNumberBase1+equationNumberLength1*(i2-equationNumberBase2+equationNumberLength2*(i3-equationNumberBase3))) + equationOffset)

// =======================================================================
// =======================================================================
#define setEquationNumber(m, ni,i1,i2,i3,  nj,j1,j2,j3 )equationNumber(m,i1,i2,i3)=indexToEquation( nj,j1,j2,j3)

// =======================================================================
// =======================================================================
#define setClassify(n,i1,i2,i3, type) classify(i1,i2,i3,n)=type

// =======================================================================
//  Macro to zero out the matrix coefficients for equations e1,e1+1,..,e2
// =======================================================================
#define zeroMatrixCoefficients( coeff,e1,e2, i1,i2,i3 )for( int m=CE(0,e1); m<=CE(0,e2+1)-1; m++ ) coeff(m,i1,i2,i3)=0.


// These macros define cI, cr, cs, crr, crs .. etc 
// Macros to evaluate matrix Coefficients for CHAMp conditions.
// File written by champGenCoeff.maple





// // =======================================================================================================================
// // Macro: Evaluate some variables that appear in the CHAMP conditions
// // =======================================================================================================================
// #beginMacro evalChampVariables()
//   const Real an1L = normal(i1,i2,i3,0), an2L =normal(i1,i2,i3,1);       // normal from left side, index with i1,i2,i3

//   // (b1L,b2L) are coefficients in the the normal derivative (left) = bL1*D_r + bL2*D_s 

//   const Real b1L = an1L*rxL(i1,i2,i3,0,0) + an2L*rxL(i1,i2,i3,0,1);     // note index left with i1,i2,i3
//   const Real b2L = an1L*rxL(i1,i2,i3,1,0) + an2L*rxL(i1,i2,i3,1,1);

//   // Macros needed when computing derivatives of b1 and b2 
//   #define b1Lm(i1,i2,i3,j1,j2,j3) ( normal(i1,i2,i3,0)*rxL(i1,i2,i3,0,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,0,1) )
//   #define b2Lm(i1,i2,i3,j1,j2,j3) ( normal(i1,i2,i3,0)*rxL(i1,i2,i3,1,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,1,1) )

//   // Macros needed when computing derivatives of c11,c12,c22,c1,c2
//   #define c11Lm(i1,i2,i3)  (SQR(rxL(i1,i2,i3,0,0)) + SQR(rxL(i1,i2,i3,0,1)))
//   #define c12Lm(i1,i2,i3)  (rxL(i1,i2,i3,0,0)*rxL(i1,i2,i3,1,0) + rxL(i1,i2,i3,0,1)*rxL(i1,i2,i3,1,1))
//   #define c21Lm(i1,i2,i3)  (rxL(i1,i2,i3,0,0)*rxL(i1,i2,i3,1,0) + rxL(i1,i2,i3,0,1)*rxL(i1,i2,i3,1,1))
//   #define c22Lm(i1,i2,i3)  (SQR(rxL(i1,i2,i3,1,0)) + SQR(rxL(i1,i2,i3,1,1)))
//   #define c1Lm(i1,i2,i3)  (RXX2(i1,i2,i3)+RYY2(i1,i2,i3))
//   #define c2Lm(i1,i2,i3)  (SXX2(i1,i2,i3)+SYY2(i1,i2,i3))
//   // (b1R,b2R) the coefficients in the normal derivative (right) = bR1*TR_r + bR2*TR_s 

//   const Real b1R = an1L*rxR(j1,j2,j3,0,0) + an2L*rxR(j1,j2,j3,0,1);     // note index right with j1,j2,j3
//   const Real b2R = an1L*rxR(j1,j2,j3,1,0) + an2L*rxR(j1,j2,j3,1,1);

//   // Macros needed when computing derivatives of b1 and b2 
//   #define b1Rm(i1,i2,i3,j1,j2,j3) ( normal(i1,i2,i3,0)*rxR(j1,j2,j3,0,0) + normal(i1,i2,i3,1)*rxR(j1,j2,j3,0,1) )
//   #define b2Rm(i1,i2,i3,j1,j2,j3) ( normal(i1,i2,i3,0)*rxR(j1,j2,j3,1,0) + normal(i1,i2,i3,1)*rxR(j1,j2,j3,1,1) )

//   #define c11Rm(j1,j2,j3)  (SQR(rxR(j1,j2,j3,0,0)) + SQR(rxR(j1,j2,j3,0,1)))
//   #define c12Rm(j1,j2,j3)  (rxR(j1,j2,j3,0,0)*rxR(j1,j2,j3,1,0) + rxR(j1,j2,j3,0,1)*rxR(j1,j2,j3,1,1))
//   #define c21Rm(j1,j2,j3)  (rxR(j1,j2,j3,0,0)*rxR(j1,j2,j3,1,0) + rxR(j1,j2,j3,0,1)*rxR(j1,j2,j3,1,1))
//   #define c22Rm(j1,j2,j3)  (SQR(rxR(j1,j2,j3,1,0)) + SQR(rxR(j1,j2,j3,1,1)))
//   #define c1Rm(j1,j2,j3)  (RXX2(j1,j2,j3)+RYY2(j1,j2,j3))
//   #define c2Rm(j1,j2,j3)  (SXX2(j1,j2,j3)+SYY2(j1,j2,j3))

//   // Use macros to compute derivatives of metrics 
//   #define Diffr2(name,i1,i2,i3,axis,dr)  ((name(i1+1,i2,i3)-name(i1-1,i2,i3))/(2.*dr[axis]))
//   #define Diffs2(name,i1,i2,i3,axis,dr)  ((name(i1,i2+1,i3)-name(i1,i2-1,i3))/(2.*dr[axis]))
//   #define Diffrr2(name,i1,i2,i3,axis,dr) ((name(i1+1,i2,i3)-2.*name(i1,i2,i3)+name(i1-1,i2,i3))/(SQR(dr[axis])))
//   #define Diffss2(name,i1,i2,i3,axis,dr) ((name(i1,i2+1,i3)-2.*name(i1,i2,i3)+name(i1,i2-1,i3))/(SQR(dr[axis])))
//   #define Diffr4(name,i1,i2,i3,axis,dr)  (((1./12.)*name(i1-2,i2,i3)-(2./3.)*name(i1-1,i2,i3)+(2./3.)*name(i1+1,i2,i3)-(1./12.)*name(i1+2,i2,i3))/(dr[axis]))
//   #define Diffs4(name,i1,i2,i3,axis,dr)  (((1./12.)*name(i1,i2-2,i3)-(2./3.)*name(i1,i2-1,i3)+(2./3.)*name(i1,i2+1,i3)-(1./12.)*name(i1,i2+2,i3))/(dr[axis]))
//   #define Diffrr4(name,i1,i2,i3,axis,dr)  ((-1.*name(i1-2,i2,i3)+16.*name(i1-1,i2,i3)-30.*name(i1,i2,i3)+16.*name(i1+1,i2,i3)-1.*name(i1+2,i2,i3))/(12.*SQR(dr[axis])))
//   #define Diffss4(name,i1,i2,i3,axis,dr)  ((-1.*name(i1,i2-2,i3)+16.*name(i1,i2-1,i3)-30.*name(i1,i2,i3)+16.*name(i1,i2+1,i3)-1.*name(i1,i2+2,i3))/(12.*SQR(dr[axis])))
//   #define Diffrrr2(name,i1,i2,i3,axis,dr)  (((-1./2.)*name(i1-2,i2,i3)+1.*name(i1-1,i2,i3)-1.*name(i1+1,i2,i3)+(1./2.)*name(i1+2,i2,i3))/(pow(dr[axis],3)))
//   #define Diffsss2(name,i1,i2,i3,axis,dr)  (((-1./2.)*name(i1,i2-2,i3)+1.*name(i1,i2-1,i3)-1.*name(i1,i2+1,i3)+(1./2.)*name(i1,i2+2,i3))/(pow(dr[axis],3)))
//   #define Diffrrrr2(name,i1,i2,i3,axis,dr)  ((name(i1-2,i2,i3)-4.*name(i1-1,i2,i3)+6.*name(i1,i2,i3)-4.*name(i1+1,i2,i3)+name(i1+2,i2,i3))/(pow(dr[axis],4)))
//   #define Diffssss2(name,i1,i2,i3,axis,dr)  ((name(i1,i2-2,i3)-4.*name(i1,i2-1,i3)+6.*name(i1,i2,i3)-4.*name(i1,i2+1,i3)+name(i1,i2+2,i3))/(pow(dr[axis],4)))
//   #define Diffrs4(name,i1,i2,i3,axis1,axis2,dr) ((64.*name(i1+1,i2+1,i3)-64.*name(i1+1,i2-1,i3)-8.*name(i1+1,i2+2,i3)+8.*name(i1+1,i2-2,i3)-64.*name(i1-1,i2+1,i3)+64.*name(i1-1,i2-1,i3)+8.*name(i1-1,i2+2,i3)-8.*name(i1-1,i2-2,i3)-8.*name(i1+2,i2+1,i3)+8.*name(i1+2,i2-1,i3)+name(i1+2,i2+2,i3)-name(i1+2,i2-2,i3)+8.*name(i1-2,i2+1,i3)-8.* name(i1-2,i2-1,i3)-name(i1-2,i2+2,i3)+name(i1-2,i2-2,i3))/(dr[axis2]*dr[axis1]*144.))
  
//   #define Diffbr2(name,i1,i2,i3,j1,j2,j3,axis,dr)  ((name(i1+1,i2,i3,j1+1,j2,j3)-name(i1-1,i2,i3,j1-1,j2,j3))/(2.*dr[axis]))
//   #define Diffbs2(name,i1,i2,i3,j1,j2,j3,axis,dr)  ((name(i1,i2+1,i3,j1,j2+1,j3)-name(i1,i2-1,i3,j1,j2-1,j3))/(2.*dr[axis]))
//   #define Diffbrr2(name,i1,i2,i3,j1,j2,j3,axis,dr) ((name(i1+1,i2,i3,j1+1,j2,j3)-2.*name(i1,i2,i3,j1,j2,j3)+name(i1-1,i2,i3,j1-1,j2,j3))/(SQR(dr[axis])))
//   #define Diffbss2(name,i1,i2,i3,j1,j2,j3,axis,dr) ((name(i1,i2+1,i3,j1,j2+1,j3)-2.*name(i1,i2,i3,j1,j2,j3)+name(i1,i2-1,i3,j1,j2-1,j3))/(SQR(dr[axis])))
//   #define Diffbr4(name,i1,i2,i3,j1,j2,j3,axis,dr)  (((1./12.)*name(i1-2,i2,i3,j1-2,j2,j3)-(2./3.)*name(i1-1,i2,i3,j1-1,j2,j3)+(2./3.)*name(i1+1,i2,i3,j1+1,j2,j3)-(1./12.)*name(i1+2,i2,i3,j1+2,j2,j3))/(dr[axis]))
//   #define Diffbs4(name,i1,i2,i3,j1,j2,j3,axis,dr)  (((1./12.)*name(i1,i2-2,i3,j1,j2-2,j3)-(2./3.)*name(i1,i2-1,i3,j1,j2-1,j3)+(2./3.)*name(i1,i2+1,i3,j1,j2+1,j3)-(1./12.)*name(i1,i2+2,i3,j1,j2+2,j3))/(dr[axis]))
//   #define Diffbrr4(name,i1,i2,i3,j1,j2,j3,axis,dr)  ((-name(i1-2,i2,i3,j1-2,j2,j3)+16.*name(i1-1,i2,i3,j1-1,j2,j3)-30.*name(i1,i2,i3,j1,j2,j3)+16.*name(i1+1,i2,i3,j1+1,j2,j3)-name(i1+2,i2,i3,j1+2,j2,j3))/(12.*SQR(dr[axis])))
//   #define Diffbss4(name,i1,i2,i3,j1,j2,j3,axis,dr)  ((-name(i1,i2-2,i3,j1,j2-2,j3)+16.*name(i1,i2-1,i3,j1,j2-1,j3)-30.*name(i1,i2,i3,j1,j2,j3)+16.*name(i1,i2+1,i3,j1,j2+1,j3)-name(i1,i2+2,i3,j1,j2+2,j3))/(12.*SQR(dr[axis])))
//   #define Diffbrrr2(name,i1,i2,i3,j1,j2,j3,axis,dr)  (((-1./2.)*name(i1-2,i2,i3,j1-2,j2,j3)+1.*name(i1-1,i2,i3,j1-1,j2,j3)-1.*name(i1+1,i2,i3,j1+1,j2,j3)+(1./2.)*name(i1+2,i2,i3,j1+2,j2,j3)/(pow(dr[axis],3)))
//   #define Diffbsss2(name,i1,i2,i3,j1,j2,j3,axis,dr)  (((-1./2.)*name(i1,i2-2,i3,j1,j2-2,j3)+1.*name(i1,i2-1,i3,j1,j2-1,j3)-1.*name(i1,i2+1,i3,j1,j2+1,j3)+(1./2.)*name(i1,i2+2,i3,j1,j2+2,j3))/(pow(dr[axis],3)))
//   #define Diffbrrrr2(name,i1,i2,i3,j1,j2,j3,axis,dr)  ((name(i1-2,i2,i3,j1-2,j2,j3)-4.*name(i1-1,i2,i3,j1-1,j2,j3)+6.*name(i1,i2,i3,j1,j2,j3)-4.*name(i1+1,i2,i3,j1+1,j2,j3)+name(i1+2,i2,i3,j1+2,j2,j3))/(pow(dr[axis],4)))
//   #define Diffbssss2(name,i1,i2,i3,j1,j2,j3,axis,dr)  ((name(i1,i2-2,i3,j1,j2-2,j3)-4.*name(i1,i2-1,i3,j1,j2-1,j3)+6.*name(i1,i2,i3,j1,j2,j3)-4.*name(i1,i2+1,i3,j1,j2+1,j3)+name(i1,i2+2,i3,j1,j2+2,j3))/(pow(dr[axis],4)))
//   #define Diffbrs4(name,i1,i2,i3,j1,j2,j3,axis1,axis2,dr) ((64.*name(i1+1,i2+1,i3,j1+1,j2+1,j3)-64.*name(i1+1,i2-1,i3,j1+1,j2-1,j3)-8.*name(i1+1,i2+2,i3,j1+1,j2+2,j3)+8.*name(i1+1,i2-2,i3,j1+1,j2-2,j3)-64.*name(i1-1,i2+1,i3,j1-1,j2+1,j3)+64.*name(i1-1,i2-1,i3,j1-1,j2-1,j3)+8.*name(i1-1,i2+2,i3,j1-1,j2+2,j3)-8.*name(i1-1,i2-2,i3,j1-1,j2-2,j3)-8.*name(i1+2,i2+1,i3,j1+2,j2+1,j3)+8.*name(i1+2,i2-1,i3,j1+2,j2-1,j3)+name(i1+2,i2+2,i3,j1+2,j2+2,j3)-name(i1+2,i2-2,i3,j1+2,j2-2,j3)+8.*name(i1-2,i2+1,i3,j1-2,j2+1,j3)-8.*name(i1-2,i2-1,i3,j1-2,j2-1,j3)-name(i1-2,i2+2,i3,j1-2,j2+2,j3)+name(i1-2,i2-2,i3,j1-2,j2-2,j3))/(dr[axis2]*dr[axis1]*144.))
  
      
//   // Use macros to compute derivatives of metrics 
//   // Macros in cgux2a need the following to be set:
//   //    dr12(axis) = 1./( 2*dr(axis) )  : for computing D0 
//   //    inverseVertexDerivative         : name of metric array

//   #define d12(axis) dr12L[axis]
//   #define inverseVertexDerivative rxLeft
//   Real r1xxL = RXX2(i1,i2,i3);  // r1 == R
//   Real r1yyL = RYY2(i1,i2,i3);
//   Real r2xxL = SXX2(i1,i2,i3);  // r2 == S 
//   Real r2yyL = SYY2(i1,i2,i3);
//   const Real c1L  = r1xxL + r1yyL ;                                  // coeff of D/D r_1 in Laplacian
//   const Real c2L  = r2xxL + r2yyL ;                                  // coeff of D/D r_2 in Laplacian

//   const Real c11L = SQR(rxL(i1,i2,i3,0,0)) + SQR(rxL(i1,i2,i3,0,1)); // coeff of D^2/D r_1^2 in Laplacian 
//   const Real c12L = rxL(i1,i2,i3,0,0)*rxL(i1,i2,i3,1,0) + rxL(i1,i2,i3,0,1)*rxL(i1,i2,i3,1,1); // coeff of D^2/(D r_1 D_r2) in Laplacian
//   const Real c21L = c12L;  
//   const Real c22L = SQR(rxL(i1,i2,i3,1,0)) + SQR(rxL(i1,i2,i3,1,1)); // coeff of D^2/D r_2^2 in Laplacian 
    
//   Real b1Lr2, b2Lr2;
//   if( axis2==0 )
//   {  // D0r2 
//     b1Lr2 = ( b1Lm(i1,i2+1,i3,j1,j2+1,j3) - b1Lm(i1,i2-1,i3,j1,j2-1,j3) )/(2.*drL[1]);
//     b2Lr2 = ( b2Lm(i1,i2+1,i3,j1,j2+1,j3) - b2Lm(i1,i2-1,i3,j1,j2-1,j3) )/(2.*drL[1]);
//   }
//   else
//   {
//     // D0r1 
//     b1Lr2 = ( b1Lm(i1+1,i2,i3,j1+1,j2,j3) - b1Lm(i1-1,i2,i3,j1-1,j2,j3) )/(2.*drL[0]);
//     b2Lr2 = ( b2Lm(i1+1,i2,i3,j1+1,j2,j3) - b2Lm(i1-1,i2,i3,j1-1,j2,j3) )/(2.*drL[0]);
//   }

//   //Real c11Lr2,c11Lr2Test;
//   //c11Lr2 = Diffr2(c11Lm,i1,i2,i3,0,drL);
//   //c11Lr2Test = (c11Lm(i1+1,i2,i3)-c11Lm(i1-1,i2,i3))/(drL[0]*2);

//   Real c11Lr4, c12Lr4, c22Lr4, c1Lr4, c2Lr4, b1Lr4, b2Lr4;
//   Real c11Ls4, c12Ls4, c22Ls4, c1Ls4, c2Ls4, b1Ls4, b2Ls4;
//   Real c11Lrr4, c12Lrr4, c22Lrr4, c1Lrr4, c2Lrr4, b1Lrr4, b2Lrr4;
//   Real c11Lss4, c12Lss4, c22Lss4, c1Lss4, c2Lss4, b1Lss4, b2Lss4;
//   Real c11Lrs4, c12Lrs4, c22Lrs4, c1Lrs4, c2Lrs4, b1Lrs4, b2Lrs4;
//   Real c11Lsss2, c12Lsss2, c22Lsss2, c1Lsss2, c2Lsss2, b1Lsss2, b2Lsss2;

//   if( orderOfAccuracy==4 ) // *wdh* 
//   {
//     if( axis2==0 )
//     {
//       c11Lr4=Diffr4(c11Lm,i1,i2,i3,0,drL); c12Lr4=Diffr4(c12Lm,i1,i2,i3,0,drL); c22Lr4=Diffr4(c22Lm,i1,i2,i3,0,drL);
//       c1Lr4 =Diffr4(c1Lm,i1,i2,i3,0,drL); c2Lr4 =Diffr4(c2Lm,i1,i2,i3,0,drL); 
//       //c22Lr4 = ( SQR(rxL(i1+1,i2,i3,1,0)) + SQR(rxL(i1+1,i2,i3,1,1)) - (SQR(rxL(i1-1,i2,i3,1,0)) + SQR(rxL(i1-1,i2,i3,1,1))));
//       //c22Lr4 = c22Lm(i1+1,i2,i3)-c22Lm(i1-1,i2,i3);

//       c11Ls4=Diffs4(c11Lm,i1,i2,i3,1,drL); c12Ls4=Diffs4(c12Lm,i1,i2,i3,1,drL); c22Ls4=Diffs4(c22Lm,i1,i2,i3,1,drL); 
//       c1Ls4 =Diffs4(c1Lm,i1,i2,i3,1,drL); c2Ls4 =Diffs4(c2Lm,i1,i2,i3,1,drL);

//       c11Lrr4=Diffrr4(c11Lm,i1,i2,i3,0,drL); c12Lrr4=Diffrr4(c12Lm,i1,i2,i3,0,drL); c22Lrr4=Diffrr4(c22Lm,i1,i2,i3,0,drL);
//       c1Lrr4 =Diffrr4(c1Lm,i1,i2,i3,0,drL); c2Lrr4 =Diffrr4(c2Lm,i1,i2,i3,0,drL);

//       c11Lss4=Diffss4(c11Lm,i1,i2,i3,1,drL); c12Lss4=Diffss4(c12Lm,i1,i2,i3,1,drL); c22Lss4=Diffss4(c22Lm,i1,i2,i3,1,drL); 
//       c1Lss4 =Diffss4(c1Lm,i1,i2,i3,1,drL); c2Lss4 =Diffss4(c2Lm,i1,i2,i3,1,drL);

//       c11Lsss2=Diffsss2(c11Lm,i1,i2,i3,1,drL); c12Lsss2=Diffsss2(c12Lm,i1,i2,i3,1,drL); c22Lsss2=Diffsss2(c22Lm,i1,i2,i3,1,drL); 
//       c1Lsss2 =Diffsss2(c1Lm,i1,i2,i3,1,drL); c2Lsss2 =Diffsss2(c2Lm,i1,i2,i3,1,drL);

//       c11Lrs4=Diffrs4(c11Lm,i1,i2,i3,0,1,drL); c12Lrs4=Diffrs4(c12Lm,i1,i2,i3,0,1,drL); c22Lrs4=Diffrs4(c22Lm,i1,i2,i3,0,1,drL);
//       c1Lrs4 =Diffrs4(c1Lm,i1,i2,i3,0,1,drL); c2Lrs4 =Diffrs4(c2Lm,i1,i2,i3,0,1,drL);

//       b1Lr4=Diffbr4(b1Lm,i1,i2,i3,j1,j2,j3,0,drL);
//       b1Ls4=Diffbs4(b1Lm,i1,i2,i3,j1,j2,j3,1,drL);
//       b1Lrr4=Diffbrr4(b1Lm,i1,i2,i3,j1,j2,j3,0,drL);
//       b1Lss4=Diffbss4(b1Lm,i1,i2,i3,j1,j2,j3,1,drL);
//       b1Lsss2=Diffbsss2(b1Lm,i1,i2,i3,j1,j2,j3,1,drL);
//       b1Lrs4=Diffbrs4(b1Lm,i1,i2,i3,j1,j2,j3,0,1,drL);
//       //b1Lrs4 = (64. * b1Lm(i1 + 1,i2 + 1,i3,j1 + 1,j2 + 1,j3) - 64. * b1Lm(i1 - 1,i2 + 1,i3,j1 - 1,j2 + 1,j3) - 8. * b1Lm(i1 + 2,i2 + 1,i3,j1 + 2,j2 + 1,j3) + 8. * b1Lm(i1 - 2,i2 + 1,i3,j1 - 2,j2 + 1,j3) - 64. * b1Lm(i1 + 1,i2 - 1,i3,j1 + 1,j2 - 1,j3) + 64. * b1Lm(i1 - 1,i2 - 1,i3,j1 - 1,j2 - 1,j3))//                + 8. * b1Lm(i1 + 2,i2 - 1,i3,j1 + 2,j2 - 1,j3) - 8. * b1Lm(i1 - 2,i2 - 1,i3,j1 - 2,j2 - 1,j3) - 8. * b1Lm(i1 + 1,i2 + 2,i3,j1 + 1,j2 + 2,j3) + 8. * b1Lm(i1 - 1,i2 + 2,i3,j1 - 1,j2 + 2,j3)//                + 8. * b1Lm(i1 + 1,i2 - 2,i3,j1 + 1,j2 - 2,j3) - 8. * b1Lm(i1 - 1,i2 - 2,i3,j1 - 1,j2 - 2,j3)+ 1.*b1Lm(i1 + 2,i2 + 2,i3,j1 + 2,j2 + 2,j3) - 1.*b1Lm(i1 - 2,i2 + 2,i3,j1 - 2,j2 + 2,j3);//  - b1Lm(i1 + 2,i2 - 2,i3,j1 + 2,j2 - 2,j3) + b1Lm(i1 - 2,i2 - 2,i3,j1 - 2,j2 - 2,j3));
//       //b1Lrs4 = ( normal(i1+ 2,i2+ 2,i3,0)*rxL(i1+ 2,i2+ 2,i3,0,0) );//+ normal(i1+ 2,i2+ 2,i3,1)*rxL(i1+ 2,i2+ 2,i3,0,1)) ;//b1Lm(i1 + 2,i2 + 2,i3,j1 + 2,j2 + 2,j3);// - b1Lm(i1 - 2,i2 + 2,i3,j1 - 2,j2 + 2,j3);

//       b2Lr4=Diffbr4(b2Lm,i1,i2,i3,j1,j2,j3,0,drL);
//       b2Ls4=Diffbs4(b2Lm,i1,i2,i3,j1,j2,j3,1,drL);
//       b2Lrr4=Diffbrr4(b2Lm,i1,i2,i3,j1,j2,j3,0,drL);
//       b2Lss4=Diffbss4(b2Lm,i1,i2,i3,j1,j2,j3,1,drL);
//       b2Lsss2=Diffbsss2(b2Lm,i1,i2,i3,j1,j2,j3,1,drL);
//       b2Lrs4=Diffbrs4(b2Lm,i1,i2,i3,j1,j2,j3,0,1,drL);
//     }
//     else
//     {
//       c11Lr4=Diffr4(c11Lm,i1,i2,i3,1,drL); c12Lr4=Diffr4(c12Lm,i1,i2,i3,1,drL); c22Lr4=Diffr4(c22Lm,i1,i2,i3,1,drL);
//       c1Lr4 =Diffr4(c1Lm,i1,i2,i3,1,drL); c2Lr4 =Diffr4(c2Lm,i1,i2,i3,1,drL);

//       c11Ls4=Diffs4(c11Lm,i1,i2,i3,0,drL); c12Ls4=Diffs4(c12Lm,i1,i2,i3,0,drL); c22Ls4=Diffs4(c22Lm,i1,i2,i3,0,drL); 
//       c1Ls4 =Diffs4(c1Lm,i1,i2,i3,0,drL); c2Ls4 =Diffs4(c2Lm,i1,i2,i3,0,drL);

//       c11Lrr4=Diffrr4(c11Lm,i1,i2,i3,1,drL); c12Lrr4=Diffrr4(c12Lm,i1,i2,i3,1,drL); c22Lrr4=Diffrr4(c22Lm,i1,i2,i3,1,drL);
//       c1Lrr4 =Diffrr4(c1Lm,i1,i2,i3,1,drL); c2Lrr4 =Diffrr4(c2Lm,i1,i2,i3,1,drL);

//       c11Lss4=Diffss4(c11Lm,i1,i2,i3,0,drL); c12Lss4=Diffss4(c12Lm,i1,i2,i3,0,drL); c22Lss4=Diffss4(c22Lm,i1,i2,i3,0,drL); 
//       c1Lss4 =Diffss4(c1Lm,i1,i2,i3,0,drL); c2Lss4 =Diffss4(c2Lm,i1,i2,i3,0,drL);

//       c11Lsss2=Diffsss2(c11Lm,i1,i2,i3,0,drL); c12Lsss2=Diffsss2(c12Lm,i1,i2,i3,0,drL); c22Lsss2=Diffsss2(c22Lm,i1,i2,i3,0,drL); 
//       c1Lsss2 =Diffsss2(c1Lm,i1,i2,i3,0,drL); c2Lsss2 =Diffsss2(c2Lm,i1,i2,i3,0,drL);

//       c11Lrs4=Diffrs4(c11Lm,i1,i2,i3,1,0,drL); c12Lrs4=Diffrs4(c12Lm,i1,i2,i3,1,0,drL); c22Lrs4=Diffrs4(c22Lm,i1,i2,i3,1,0,drL);
//       c1Lrs4 =Diffrs4(c1Lm,i1,i2,i3,1,0,drL); c2Lrs4 =Diffrs4(c2Lm,i1,i2,i3,1,0,drL);

//       b1Lr4=Diffbr4(b1Lm,i1,i2,i3,j1,j2,j3,1,drL);
//       b1Ls4=Diffbs4(b1Lm,i1,i2,i3,j1,j2,j3,0,drL);
//       b1Lrr4=Diffbrr4(b1Lm,i1,i2,i3,j1,j2,j3,1,drL);
//       b1Lss4=Diffbss4(b1Lm,i1,i2,i3,j1,j2,j3,0,drL);
//       b1Lsss2=Diffbsss2(b1Lm,i1,i2,i3,j1,j2,j3,0,drL);
//       b1Lrs4=Diffbrs4(b1Lm,i1,i2,i3,j1,j2,j3,1,0,drL);

//       b2Lr4=Diffbr4(b2Lm,i1,i2,i3,j1,j2,j3,1,drL);
//       b2Ls4=Diffbs4(b2Lm,i1,i2,i3,j1,j2,j3,0,drL);
//       b2Lrr4=Diffbrr4(b2Lm,i1,i2,i3,j1,j2,j3,1,drL);
//       b2Lss4=Diffbss4(b2Lm,i1,i2,i3,j1,j2,j3,0,drL);
//       b2Lsss2=Diffbsss2(b2Lm,i1,i2,i3,j1,j2,j3,0,drL);
//       b2Lrs4=Diffbrs4(b2Lm,i1,i2,i3,j1,j2,j3,1,0,drL);
//     }
//   } // end order=4

//   #define d12(axis) dr12R[axis]
//   #define inverseVertexDerivative rxRight
//   Real r1xxR = RXX2(j1,j2,j3);  // r1 == R
//   Real r1yyR = RYY2(j1,j2,j3);
//   Real r2xxR = SXX2(j1,j2,j3);  // r2 == S 
//   Real r2yyR = SYY2(j1,j2,j3);

//   const Real c1R  = r1xxR + r1yyR ;                                  // coeff of D/D r_1 in Laplacian
//   const Real c2R  = r2xxR + r2yyR ;                                  // coeff of D/D r_2 in Laplacian

//   const Real c11R = SQR(rxR(j1,j2,j3,0,0)) + SQR(rxR(j1,j2,j3,0,1)); // coeff of D^2/D r_1^2 in Laplacian 
//   const Real c12R = rxR(j1,j2,j3,0,0)*rxR(j1,j2,j3,1,0) + rxR(j1,j2,j3,0,1)*rxR(j1,j2,j3,1,1); // coeff of D^2/(D r_1 D_r2) in Laplacian
//   const Real c21R = c12R;  
//   const Real c22R = SQR(rxR(j1,j2,j3,1,0)) + SQR(rxR(j1,j2,j3,1,1)); // coeff of D^2/D r_2^2 in Laplacian 
    
//   Real c11Rr4, c12Rr4, c22Rr4, c1Rr4, c2Rr4, b1Rr4, b2Rr4;
//   Real c11Rs4, c12Rs4, c22Rs4, c1Rs4, c2Rs4, b1Rs4, b2Rs4;
//   Real c11Rrr4, c12Rrr4, c22Rrr4, c1Rrr4, c2Rrr4, b1Rrr4, b2Rrr4;
//   Real c11Rss4, c12Rss4, c22Rss4, c1Rss4, c2Rss4, b1Rss4, b2Rss4;
//   Real c11Rrs4, c12Rrs4, c22Rrs4, c1Rrs4, c2Rrs4, b1Rrs4, b2Rrs4;
//   Real c11Rsss2, c12Rsss2, c22Rsss2, c1Rsss2, c2Rsss2, b1Rsss2, b2Rsss2;

//   if( orderOfAccuracy==4 ) // *wdh* 
//   {
//     if( axis2==0 )
//     {
//       c11Rr4=Diffr4(c11Rm,j1,j2,j3,0,drR); c12Rr4=Diffr4(c12Rm,j1,j2,j3,0,drR); c22Rr4=Diffr4(c22Rm,j1,j2,j3,0,drR);
//       c1Rr4 =Diffr4(c1Rm,j1,j2,j3,0,drR); c2Rr4 =Diffr4(c2Rm,j1,j2,j3,0,drR);

//       c11Rs4=Diffs4(c11Rm,j1,j2,j3,1,drR); c12Rs4=Diffs4(c12Rm,j1,j2,j3,1,drR); c22Rs4=Diffs4(c22Rm,j1,j2,j3,1,drR); 
//       c1Rs4 =Diffs4(c1Rm,j1,j2,j3,1,drR); c2Rs4 =Diffs4(c2Rm,j1,j2,j3,1,drR);

//       c11Rrr4=Diffrr4(c11Rm,j1,j2,j3,0,drR); c12Rrr4=Diffrr4(c12Rm,j1,j2,j3,0,drR); c22Rrr4=Diffrr4(c22Rm,j1,j2,j3,0,drR);
//       c1Rrr4 =Diffrr4(c1Rm,j1,j2,j3,0,drR); c2Rrr4 =Diffrr4(c2Rm,j1,j2,j3,0,drR);

//       c11Rss4=Diffss4(c11Rm,j1,j2,j3,1,drR); c12Rss4=Diffss4(c12Rm,j1,j2,j3,1,drR); c22Rss4=Diffss4(c22Rm,j1,j2,j3,1,drR); 
//       c1Rss4 =Diffss4(c1Rm,j1,j2,j3,1,drR); c2Rss4 =Diffss4(c2Rm,j1,j2,j3,1,drR);

//       c11Rsss2=Diffsss2(c11Rm,j1,j2,j3,1,drR); c12Rsss2=Diffsss2(c12Rm,j1,j2,j3,1,drR); c22Rsss2=Diffsss2(c22Rm,j1,j2,j3,1,drR); 
//       c1Rsss2 =Diffsss2(c1Rm,j1,j2,j3,1,drR); c2Rsss2 =Diffsss2(c2Rm,j1,j2,j3,1,drR);

//       c11Rrs4=Diffrs4(c11Rm,j1,j2,j3,0,1,drR); c12Rrs4=Diffrs4(c12Rm,j1,j2,j3,0,1,drR); c22Rrs4=Diffrs4(c22Rm,j1,j2,j3,0,1,drR);
//       c1Rrs4 =Diffrs4(c1Rm,j1,j2,j3,0,1,drR); c2Rrs4 =Diffrs4(c2Rm,j1,j2,j3,0,1,drR);

//       b1Rr4=Diffbr4(b1Rm,i1,i2,i3,j1,j2,j3,0,drR);
//       b1Rs4=Diffbs4(b1Rm,i1,i2,i3,j1,j2,j3,1,drR);
//       b1Rrr4=Diffbrr4(b1Rm,i1,i2,i3,j1,j2,j3,0,drR);
//       b1Rss4=Diffbss4(b1Rm,i1,i2,i3,j1,j2,j3,1,drR);
//       b1Rsss2=Diffbsss2(b1Rm,i1,i2,i3,j1,j2,j3,1,drR);
//       b1Rrs4=Diffbrs4(b1Rm,i1,i2,i3,j1,j2,j3,0,1,drR);

//       b2Rr4=Diffbr4(b2Rm,i1,i2,i3,j1,j2,j3,0,drR);
//       b2Rs4=Diffbs4(b2Rm,i1,i2,i3,j1,j2,j3,1,drR);
//       b2Rrr4=Diffbrr4(b2Rm,i1,i2,i3,j1,j2,j3,0,drR);
//       b2Rss4=Diffbss4(b2Rm,i1,i2,i3,j1,j2,j3,1,drR);
//       b2Rsss2=Diffbsss2(b2Rm,i1,i2,i3,j1,j2,j3,1,drR);
//       b2Rrs4=Diffbrs4(b2Rm,i1,i2,i3,j1,j2,j3,0,1,drR);
//     }
//     else
//     {
//       c11Rr4=Diffr4(c11Rm,j1,j2,j3,1,drR); c12Rr4=Diffr4(c12Rm,j1,j2,j3,1,drR); c22Rr4=Diffr4(c22Rm,j1,j2,j3,1,drR);
//       c1Rr4 =Diffr4(c1Rm,j1,j2,j3,1,drR); c2Rr4 =Diffr4(c2Rm,j1,j2,j3,1,drR);

//       c11Rs4=Diffs4(c11Rm,j1,j2,j3,0,drR); c12Rs4=Diffs4(c12Rm,j1,j2,j3,0,drR); c22Rs4=Diffs4(c22Rm,j1,j2,j3,0,drR); 
//       c1Rs4 =Diffs4(c1Rm,j1,j2,j3,0,drR); c2Rs4 =Diffs4(c2Rm,j1,j2,j3,0,drR);

//       c11Rrr4=Diffrr4(c11Rm,j1,j2,j3,1,drR); c12Rrr4=Diffrr4(c12Rm,j1,j2,j3,1,drR); c22Rrr4=Diffrr4(c22Rm,j1,j2,j3,1,drR);
//       c1Rrr4 =Diffrr4(c1Rm,j1,j2,j3,1,drR); c2Rrr4 =Diffrr4(c2Rm,j1,j2,j3,1,drR);

//       c11Rss4=Diffss4(c11Rm,j1,j2,j3,0,drR); c12Rss4=Diffss4(c12Rm,j1,j2,j3,0,drR); c22Rss4=Diffss4(c22Rm,j1,j2,j3,0,drR); 
//       c1Rss4 =Diffss4(c1Rm,j1,j2,j3,0,drR); c2Rss4 =Diffss4(c2Rm,j1,j2,j3,0,drR);

//       c11Rsss2=Diffsss2(c11Rm,j1,j2,j3,0,drR); c12Rsss2=Diffsss2(c12Rm,j1,j2,j3,0,drR); c22Rsss2=Diffsss2(c22Rm,j1,j2,j3,0,drR); 
//       c1Rsss2 =Diffsss2(c1Rm,j1,j2,j3,0,drR); c2Rsss2 =Diffsss2(c2Rm,j1,j2,j3,0,drR);

//       c11Rrs4=Diffrs4(c11Rm,j1,j2,j3,1,0,drR); c12Rrs4=Diffrs4(c12Rm,j1,j2,j3,1,0,drR); c22Rrs4=Diffrs4(c22Rm,j1,j2,j3,1,0,drR);
//       c1Rrs4 =Diffrs4(c1Rm,j1,j2,j3,1,0,drR); c2Rrs4 =Diffrs4(c2Rm,j1,j2,j3,1,0,drR);

//       b1Rr4=Diffbr4(b1Rm,i1,i2,i3,j1,j2,j3,1,drR);
//       b1Rs4=Diffbs4(b1Rm,i1,i2,i3,j1,j2,j3,0,drR);
//       b1Rrr4=Diffbrr4(b1Rm,i1,i2,i3,j1,j2,j3,1,drR);
//       b1Rss4=Diffbss4(b1Rm,i1,i2,i3,j1,j2,j3,0,drR);
//       b1Rsss2=Diffbsss2(b1Rm,i1,i2,i3,j1,j2,j3,0,drR);
//       b1Rrs4=Diffbrs4(b1Rm,i1,i2,i3,j1,j2,j3,1,0,drR);

//       b2Rr4=Diffbr4(b2Rm,i1,i2,i3,j1,j2,j3,1,drR);
//       b2Rs4=Diffbs4(b2Rm,i1,i2,i3,j1,j2,j3,0,drR);
//       b2Rrr4=Diffbrr4(b2Rm,i1,i2,i3,j1,j2,j3,1,drR);
//       b2Rss4=Diffbss4(b2Rm,i1,i2,i3,j1,j2,j3,0,drR);
//       b2Rsss2=Diffbsss2(b2Rm,i1,i2,i3,j1,j2,j3,0,drR);
//       b2Rrs4=Diffbrs4(b2Rm,i1,i2,i3,j1,j2,j3,1,0,drR);
//     }
//   } // end order=4

//   Real b1Rr2, b2Rr2;
//   if( axis2==0 )
//   {  // D0r2 
//     b1Rr2 = ( b1Rm(i1,i2+1,i3,j1,j2+1,j3) - b1Rm(i1,i2-1,i3,j1,j2-1,j3) )/(2.*drR[1]);
//     b2Rr2 = ( b2Rm(i1,i2+1,i3,j1,j2+1,j3) - b2Rm(i1,i2-1,i3,j1,j2-1,j3) )/(2.*drR[1]);
//   }
//   else
//   {
//     // D0r1 
//     b1Rr2 = ( b1Rm(i1+1,i2,i3,j1+1,j2,j3) - b1Rm(i1-1,i2,i3,j1-1,j2,j3) )/(2.*drR[0]);
//     b2Rr2 = ( b2Rm(i1+1,i2,i3,j1+1,j2,j3) - b2Rm(i1-1,i2,i3,j1-1,j2,j3) )/(2.*drR[0]);
//   }

//   // --- we ned to evaluate tangential derivatives of
//   //          b2L/b1R
//   //          b2R/b1R 
//   # define b2Lb1R(i1,i2,i3,j1,j2,j3) (b2Lm(i1,i2,i3,j1,j2,j3)/b1Rm(i1,i2,i3,j1,j2,j3)) 
//   # define b2Rb1R(i1,i2,i3,j1,j2,j3) (b2Rm(i1,i2,i3,j1,j2,j3)/b1Rm(i1,i2,i3,j1,j2,j3)) 

//   Real b2Lb1Rr2, b2Rb1Rr2;
//   // ** CHECK ME **
//   if( axis2==0 )
//   {  // D0r2 
//     b2Lb1Rr2 = ( b2Lb1R(i1,i2+1,i3,j1,j2+1,j3) - b2Lb1R(i1,i2-1,i3,j1,j2-1,j3) )/(2.*drR[1]) ; //  D( (b2L/b1R) )/D( rR2)
//     b2Rb1Rr2 = ( b2Rb1R(i1,i2+1,i3,j1,j2+1,j3) - b2Rb1R(i1,i2-1,i3,j1,j2-1,j3) )/(2.*drR[1]) ; //  D( (b2L/b1R) )/D( rR2)
//   }
//   else
//   {
//     // D0r1 
//     b2Lb1Rr2 = ( b2Lb1R(i1+1,i2,i3,j1+1,j2,j3) - b2Lb1R(i1-1,i2,i3,j1-1,j2,j3) )/(2.*drR[0]) ;
//     b2Rb1Rr2 = ( b2Rb1R(i1+1,i2,i3,j1+1,j2,j3) - b2Rb1R(i1-1,i2,i3,j1-1,j2,j3) )/(2.*drR[0]) ;
//   }
// #endMacro

// ===============================================================================
// Macro: 
// Retrieve information about the domain solver on the opposite side:
//   int domain2, side2, axis2
//   DomainSolver solver2
//   GridFunction gf2
//   CompositeGrid cg2
//   MappedGrid mg2
//   Boundary index: J1,J2,J3
//   bool isRectaangular2
//   real dx2[]
//   CompositeGridOperators cgop2
//   MappedGridOperators mgop2
// ===============================================================================

// ==========================================================================
// Macro: fill the matrix with extrapolation for a given ghost=1,2,3,...
// ==========================================================================

//#Include "champ4InterfaceStencilCurvilinearMacro.h"

/*
*/



// ==========================================================================
// Macro: Evaluate the advection coefficients 
// ==========================================================================

// ====================================================================================================
//
//   Fill in the CHAMP implicit boundary conditions.
//   ------------------------------------------------
//
//  The coefficients in the matrix are defined by two arrays:
//
//    coeff(m,i1,i2,i3), m=0,1,2,...    : the coefficients associated with the equation(s) at the grid point (i1,i2,i3)
//                                        For a system of equations,
//                                                m=0,1,...,stencilDim-1             : holds the first eqn. in the system,
//                                                m=stencilDim,1,...,2*stencilDim-1  : holds the second eqn, etc.
//                                        (stencilDim is defined below)
//
//    equationNumber(m,i1,i2,i3)        : a global index that identifies the "unknown" in the system
//                                        (the global index of each unknown is defined by the indexToEquation macro)
//
// NOTES:
//   See cg/ins/src/insImp.h     : for fortran version macros
//       cg/ins/src/insImpINS.bf : for examples
// ====================================================================================================
int Cgad::champBoundaryConditions( realCompositeGridFunction & coeffcg, Parameters & parameters, Real dt )
{

    const int debug = parameters.dbase.get<int>("debug");

    CompositeGrid & cg = *coeffcg.getCompositeGrid();
    const int numberOfDimensions = cg.numberOfDimensions();
    int numberOfComponents = 1;
    
  // const int eq1=0, eq2=1;   // equation numbers
  // const int uc=0, vc=1;     // component numbers
    
    const std::vector<real> & kappa  = parameters.dbase.get<std::vector<real> >("kappa");
    const real & thermalConductivity = parameters.dbase.get<real>("thermalConductivity");

    const bool & variableAdvection = parameters.dbase.get<bool >("variableAdvection");

    const int & multiDomainProblem = parameters.dbase.get<int>("multiDomainProblem"); 
    const bool twilightZoneFlow = parameters.dbase.get<bool>("twilightZoneFlow");

    const int orderOfAccuracy = parameters.dbase.get<int>("orderOfAccuracy"); 

    const int extrapOrder = orderOfAccuracy+1;

    const Real extrapCoeff3[] = {1.,-3.,3.,-1.};
    const Real extrapCoeff4[] = {1.,-4.,6.,-4.,1.};
    const Real extrapCoeff5[] = {1.,-5.,10.,-10.,5.,-1.};
    const Real *extrapCoeff;
    if( extrapOrder==3 )
        extrapCoeff = extrapCoeff3;
    else if( extrapOrder==4 )
        extrapCoeff = extrapCoeff4;
    else if( extrapOrder==5 )
        extrapCoeff = extrapCoeff5;
    else
      {
        printF("Cgad::champBoundaryConditions: unexpected extrapOrder=%d\n",extrapOrder);
        OV_ABORT("ERROR");
      }


    int mapleOption = 0;  // switch between versions  0 = old,  1=Sijia's version
    if( orderOfAccuracy==4 )
        mapleOption=1;

  // Sl = optimized Schwartz parameter
  // theta = Ktc1/Ktc2
  // beta  = Dc1/Dc2
    if( !parameters.dbase.has_key("champParameters") )
    {
    // NOTE: This array may have already been created in CgAd: setupPdeParameters

        parameters.dbase.put<RealArray>("champParameters");
        RealArray & champParameters = parameters.dbase.get<RealArray>("champParameters");

    // ---- number of champ paramters we save-- 
        const int numChampPar=10; 

    //  pl    = champParameters(0,side,axis,grid);    // optimized Scwartz Parameter for side 1
    //  pr    = champParameters(1,side,axis,grid);    // optimized Scwartz Parameter for side 2
    //  theta = champParameters(2,side,axis,grid);    // K1/K2
    //  beta  = champParameters(3,side,axis,grid);    // D1/D2   
    //  Sl    = champParameters(4,side,axis,grid);  
    //  dxs   = champParameters(5,side,axis,grid);
    //  DL    = champParameters(6,side,axis,grid);
    //  KL    = champParameters(7,side,axis,grid);
    //  DR    = champParameters(8,side,axis,grid);
    //  KR    = champParameters(9,side,axis,grid);
    // 

        champParameters.redim(numChampPar,2,3,cg.numberOfComponentGrids());
        champParameters=-123456.; 



    }
    RealArray & champParameters = parameters.dbase.get<RealArray>("champParameters");

    const IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");

    printF("\n -- CGAD-- champBoundaryConditions: multiDomainProblem=%d,  kappa=%g, thermalConductivity=%g\n",multiDomainProblem,kappa[0],thermalConductivity);



    Range all;
    Range c0(0,0), c1(1,1);  // c0 = first component, c1 = second component
    
    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    Index Jv[3], &J1=Jv[0], &J2=Jv[1], &J3=Jv[2];  // opposite side 
    Index Igv[3], &Ig1=Igv[0], &Ig2=Igv[1], &Ig3=Igv[2];

    Index Ibv[3], &Ib1=Ibv[0], &Ib2=Ibv[1], &Ib3=Ibv[2];
    Index Jbv[3], &Jb1=Jbv[0], &Jb2=Jbv[1], &Jb3=Jbv[2];  // opposite side 
    int j1,j2,j3, i1m,i2m,i3m, m1,m2,m3;
    int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2];
    int isv[3], &is1=isv[0], &is2=isv[1], &is3=isv[2];
    int side,axis;
        
        
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        mg.update( MappedGrid::THEvertexBoundaryNormal );
        realArray & vertex = mg.vertex();
        intArray & mask = mg.mask();
        realMappedGridFunction & coeff = coeffcg[grid];
        MappedGridOperators & mgop = *coeff.getOperators();

        const int isRectangular=mg.isRectangular();

        real dx[3]={1.,1.,1.};
        if( isRectangular )
            mg.getDeltaX(dx);
        else
            mg.update(MappedGrid::THEinverseVertexDerivative);



        assert( coeff.sparse!=NULL );
        
        SparseRepForMGF & sparse = *coeff.sparse;
        int numberOfComponentsForCoefficients = sparse.numberOfComponents;  // size of the system of equations
        int numberOfGhostLines = sparse.numberOfGhostLines;
        int stencilSize = sparse.stencilSize;
        int stencilDim=stencilSize*numberOfComponentsForCoefficients; // number of coefficients per equation
        
        
        const int equationOffset=sparse.equationOffset;
        IntegerArray & equationNumber = sparse.equationNumber;
        IntegerArray & classify = sparse.classify;
        
        const int equationNumberBase1  =equationNumber.getBase(1);
        const int equationNumberLength1=equationNumber.getLength(1);
        const int equationNumberBase2  =equationNumber.getBase(2);
        const int equationNumberLength2=equationNumber.getLength(2);
        const int equationNumberBase3  =equationNumber.getBase(3);
        
        const int orderOfAccuracy2 = mgop.getOrderOfAccuracy(); 
        assert( orderOfAccuracy2==orderOfAccuracy );
        
    // stencil width's and half-width's :
        const int width = orderOfAccuracy+1;
        const int halfWidth1 = (width-1)/2;
        const int halfWidth2 = numberOfDimensions>1 ? halfWidth1 : 0;
        const int halfWidth3 = numberOfDimensions>2 ? halfWidth1 : 0;
        
        Range M0 = stencilSize;
        Range M = coeff.dimension(0);

        const int numGhost = halfWidth1; // number of ghost points 

    // For coefficients from Sijia: 
        Real coefB[5][5]={0.,0.,0.,0.,0.,
                                            0.,0.,0.,0.,0.,
                                            0.,0.,0.,0.,0.,
                                            0.,0.,0.,0.,0.,
                                            0.,0.,0.,0.,0.}; // initialize
        #define coefA(i,j) coefB[i+halfWidth1][j+halfWidth2]
        Real coefBg[5][5]={0.,0.,0.,0.,0.,
                                            0.,0.,0.,0.,0.,
                                            0.,0.,0.,0.,0.,
                                            0.,0.,0.,0.,0.,
                                            0.,0.,0.,0.,0.}; // initialize
        #define coefAg(i,j) coefBg[i+halfWidth1][j+halfWidth2]    
            
        const bool fillSecondGhostUsingExtrapolation=true; 

        ForBoundary(side,axis)
        {
            if( interfaceType(side,axis,grid) == Parameters::heatFluxInterface && mg.boundaryCondition(side,axis)>0 )
            {
                printF("\n +++++ CHAMP BC: FILL MATRIX BC FOR (grid,side,axis)=(%d,%d,%d) order=%d, numGhost=%d fillSecondGhostUsingExtrapolation=%d mapleOption=%d (heatFluxInterface) ++++++++\n",
                          grid,side,axis,orderOfAccuracy,numGhost,(int)fillSecondGhostUsingExtrapolation,mapleOption );

        // Set the index-shift for this side
                is1=is2=is3=0;
                isv[axis]=1-2*side;

        // *wdh* April 16, 2022 -- fix for periodic in tangnetial direction
        // getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                getBoundaryIndex(mg.indexRange(),side,axis,Ib1,Ib2,Ib3);

                OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal); 

              

        // Retrieve information about the domain solver on the opposite side:
        //   int domain2, side2, axis2
        //   DomainSolver solver2
        //   GridFunction gf2
        //   CompositeGrid cg2
        //   MappedGrid mg2
        //   Boundary index: J1,J2,J3
        //   CompositeGridOperators cgop2
        //   MappedGridOperators mgop2        
          // For testing the solver on the opposite side is just the the same
                    DomainSolver *pSolver2 = this;
                    int domain2 = 0, grid2 = grid, side2=1-side, axis2=axis; // right side of the left grid.
          //int domain2 = 0, grid2 = grid, side2=side, axis2=axis; 
                    if( multiDomainProblem )
                    {
           // --- GET THE SOLVER FOR THE OPPOSITE SIDE --- 
                        bool sameSide=false; 
                        GridFaceDescriptor & gfd = getInterfaceGridFaceDescriptor( grid, side, axis, parameters, sameSide );
                        domain2 = gfd.domain, grid2 = gfd.grid, side2=gfd.side, axis2=gfd.axis; 
            // Here is the Cgmp object:
                        DomainSolver *pCgmp = parameters.dbase.get<DomainSolver*>("multiDomainSolver");
                        assert( pCgmp->domainSolver[domain2]!=NULL );
            // Here is the solver in the other side: 
            // DomainSolver & solver2 = *(pCgmp->domainSolver[domain2]);
                        pSolver2 = pCgmp->domainSolver[domain2]; 
                    }
                    DomainSolver & solver2 = *pSolver2;
                    const int current2 = solver2.current;   // index of grid function that holds the current solution 
                    printF("FILL CHAMP MATRIX: OPPOSITE SIDE: domain2=%d, grid2=%i, side2=%i, axis2=%i className=%s, name=%s current2=%d\n",
                              domain2,grid2,side2,axis2,(const char*)solver2.getClassName(),(const char*)solver2.getName(),current2);
          // -- find the CompositeGrid and MappedGrid for the opposite side
                    GridFunction & gf2 = solver2.gf[current2];  
                    CompositeGrid & cg2 = gf2.cg; 
                    MappedGrid & mg2 = cg2[grid2];
          // *wdh* April 16, 2022 -- fix for periodic in tangnetial direction
          // getBoundaryIndex(mg2.gridIndexRange(),side2,axis2,Jb1,Jb2,Jb3);
                    getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3);
                    const bool isRectangular2 = mg2.isRectangular();
                    real dx2[3]={1.,1.,1.};
                    if( isRectangular2 )
                        mg2.getDeltaX(dx2);
                    else
                        mg2.update(MappedGrid::THEinverseVertexDerivative);  
          // --------- Operators on Opposite side -----------
                    CompositeGridOperators & cgop2 = *gf2.u.getOperators();
                    MappedGridOperators & mgop2 = cgop2[grid2]; 
          // ---- consistency checks -------
                    printF(" THIS SIDE    : Ib1=[%d,%d] Ib2=[%d,%d] Ib3=[%d,%d]\n",Ib1.getBase(),Ib1.getBound(), Ib2.getBase(),Ib2.getBound(), Ib3.getBase(),Ib3.getBound());
                    printF(" OPPOSITE SIDE: Jb1=[%d,%d] Jb2=[%d,%d] Jb3=[%d,%d]\n",Jb1.getBase(),Jb1.getBound(), Jb2.getBase(),Jb2.getBound(), Jb3.getBase(),Jb3.getBound());
                    int axisp1  = (axis  +1 ) % numberOfDimensions;
                    int axis2p1 = (axis2 +1 ) % numberOfDimensions;
          // For now we that both sides are indexed in a similiar way: Ib1 <-> Jb1 and Ib2 <-> Jb2
                    if(  !( Ibv[axis].getLength() == Jbv[axis2].getLength() && Ibv[axisp1] == Jbv[axis2p1] ) )
                    {
                        printF("CGAD:CHAMP:FILL MATRIX CURVILINEAR: ERROR: Indexes of interfaces do not match -- this is currently required\n");
                        printF(" THIS SIDE    : Ib1=[%d,%d] Ib2=[%d,%d] Ib3=[%d,%d]\n",Ib1.getBase(),Ib1.getBound(), Ib2.getBase(),Ib2.getBound(), Ib3.getBase(),Ib3.getBound());
                        printF(" OPPOSITE SIDE: Jb1=[%d,%d] Jb2=[%d,%d] Jb3=[%d,%d]\n",Jb1.getBase(),Jb1.getBound(), Jb2.getBase(),Jb2.getBound(), Jb3.getBase(),Jb3.getBound());
                        OV_ABORT("champBoundaryConditions::ERROR");
                    }
                    if( isRectangular != isRectangular2 )
                    {
                        printF("CGAD:CHAMP: isRectangular=%d is NOT EQUAL TO isRectagular2=%d -- this is currently required\n", 
                                      (int)isRectangular2,(int)isRectangular2);
                        OV_ABORT("champBoundaryConditions::ERROR");
                    }



        // NO: 
        // Champ h = (=/-) dx[axis] : sign depends on which side we are on 
        //  isc = 1 : left side:  h=-dx
        //      =-1 : right side: h=+dx
        // const Real isc = 1-2*side; 
        // const Real dxs = -isc*dx[axis]; 

        // const Real dxs = dx[axis]; 
                const Real dxs = dx2[axis2]; 

                if( isRectangular && isRectangular2 )
                    printF("+++++ CHAMP BC: dx[axis]=%9.3e. Opposite side: dx2[axis2]=%9.3e\n",dx[axis],dx2[axis2]);


        // --- Get grid spacing in normal direction for computing lambdaD = DL*dt/dxL*dxL  ---
                Real dnL = dx[axis];  

                if( !isRectangular )
                {
          // -- Estimate dn in normal direction ---
          //    dn  = dr /[ n.( dr/dx) ] 
          // Do this for now 
                    OV_GET_SERIAL_ARRAY(real, mg.inverseVertexDerivative(),rxLeft );  // metric on left 
          //  The rx array is only 4d, here is a macro to make it look 5d 
                    #define rxL(i1,i2,i3,m1,m2)  rxLeft(i1,i2,i3,(m1)+numberOfDimensions*(m2))          
                    const Real drL[3]   = { mg.gridSpacing(0),mg.gridSpacing(1),mg.gridSpacing(2) };     // grid spacing on left 

                    const int i1=Ib1.getBase(), i2=Ib2.getBase(), i3=Ib3.getBase();    // find dx at this point 
                    Real bL = -( normal(i1,i2,i3,0)*rxL(i1,i2,i3,axis,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,axis,1)) ;
          // Real b2L = -( normal(i1,i2,i3,0)*rxL(i1,i2,i3,1,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,1,1)) ;  
                                    
                    dnL = fabs( drL[axis]/bL );    // approx. normal grid spacing on left 
                    printF("+++++ CHAMP BC: curvilinear grid: estimate normal grid spacing: i1=%d i2=%d drL=%e, bL=%e, "
                                  "normal=[%e,%e] rx=[%e,%e] dnL=%e\n",
                                  i1,i2,drL[axis], bL,normal(i1,i2,i3,0),normal(i1,i2,i3,1),rxL(i1,i2,i3,axis,0),rxL(i1,i2,i3,axis,1),dnL);
          // OV_ABORT("stop here for now");
                }

        // --- get the champ parameters and store in the array champParameters ---
                if( multiDomainProblem )
                {
                    getChampParameters( grid, side, axis, grid2, side2, axis2, dt, dnL, parameters, champParameters );
                }
                else
                {
          // single domain problem : we must be tesing the CHAMP conditions

                    Real pl=1., pr=.1, theta=.5, beta=2.; 

                    champParameters(0,side,axis,grid)=pl;     // optimized Scwartz Parameter for side 1 : Sl = pl/dx
                    champParameters(1,side,axis,grid)=pr;     // optimized Scwartz Parameter for side 2 : Sr = pr/dx 
                    champParameters(2,side,axis,grid)=theta;  // K1/K2
                    champParameters(3,side,axis,grid)=beta;   // D1/D2                 
                }

                const Real pl    = champParameters(0,side,axis,grid);    // optimized Scwartz Parameter for side 1
        // pr is not used here and (currently) may not be correct: 
        // const Real pr    = champParameters(1,side,axis,grid);    // optimized Scwartz Parameter for side 2
                const Real theta = champParameters(2,side,axis,grid);    // K1/K2
                const Real beta  = champParameters(3,side,axis,grid);    // D1/D2  

                const Real DL    = champParameters(6,side,axis,grid);
                const Real KL    = champParameters(7,side,axis,grid);
                const Real DR    = champParameters(8,side,axis,grid);
                const Real KR    = champParameters(9,side,axis,grid);            


        // -- advection coefficients 
                std::vector<real> & aL = parameters.dbase.get<std::vector<real> >("a");
                std::vector<real> & bL = parameters.dbase.get<std::vector<real> >("b");
                std::vector<real> & cL = parameters.dbase.get<std::vector<real> >("c"); 

                Parameters & parametersR = solver2.parameters;
                std::vector<real> & aR = parametersR.dbase.get<std::vector<real> >("a");
                std::vector<real> & bR = parametersR.dbase.get<std::vector<real> >("b");
                std::vector<real> & cR = parametersR.dbase.get<std::vector<real> >("c"); 

                Real u1L=aL[0]; 
                Real u2L=bL[0]; 
                Real u1R=aR[0]; 
                Real u2R=bR[0];

        // advection coefficients scaled by 1/D 
                Real u1DL=aL[0]/DL; 
                Real u2DL=bL[0]/DL; 
                Real u1DR=aR[0]/DR; 
                Real u2DR=bR[0]/DR;         
                bool advectionIsOn=false;
                if( u1L!=0. || u2L!=0. || u1R!=0. || u2R!=0. || variableAdvection )
                {
                    advectionIsOn=true;

                    if( isRectangular )
                        OV_ABORT("champBC: finish me for rectangular + advection");

          // if( orderOfAccuracy==4 )
          //   OV_ABORT("champBC: finish me for order=4 + advection");
                }

                const int gridAdvect = variableAdvection ? grid : 0;
                realCompositeGridFunction * pAdvectVarL = variableAdvection ? parameters.dbase.get<realCompositeGridFunction*>("advectVar") : &coeffcg;
                OV_GET_SERIAL_ARRAY_CONDITIONAL(real,(*pAdvectVarL)[gridAdvect],advectVarL,variableAdvection); 

                const int grid2Advect = variableAdvection ? grid2 : 0;
                realCompositeGridFunction * pAdvectVarR= variableAdvection ? parametersR.dbase.get<realCompositeGridFunction*>("advectVar") : &coeffcg;
                OV_GET_SERIAL_ARRAY_CONDITIONAL(real,(*pAdvectVarR)[grid2Advect],advectVarR,variableAdvection); 

                printF("***** fillChampBC: (side,axis,grid)=(%d,%d,%d) : Kleft/Kright=%g, Dleft/Dright=%g, CHAMP parameter pl=%g ****\n",
                          side,axis,grid,theta,beta,pl); 
                printF("***** fillChampBC: DL=%g, KL=%g, u1DL=%g, u2DL=%g, advectionIsOn=%d, variableAdvection=%d\n",
                            DL,KL,u1DL,u2DL,(int)advectionIsOn,(int)variableAdvection );
                printF("***** fillChampBC: DR=%g, KR=%g, u1DR=%g, u2DR=%g\n",DR,KR,u1DR,u2DR);

                
        // --- Fill in the CHAMP interface conditions -----
        //   [ theta*D_x + h*L[beta] ] + Sl*[ I + theta*h*D_x + (h^2)/2 *L[beta] ] 
        //      L[beta] = beta*D_xx + (beta-1)*D_yy        
        //  LHS = a0*T + a1*D_n T + a2*T_xx + a3*T_yy 
        // NOTES:
        //     (1) The CHAMP condition is centered on the boundary point (i1,i2,i3)
        //         but is associated with the ghost-point (i1m,i2m,i3m)
                
                
                const int eqnStart=0, eqnEnd=0;       // equation numbers

                if( isRectangular )
                {
          // ------------------ RECTANGULAR GRID CASE ----------------------------
                    printF("\n &&&&&& CHAMP BC: FILL MATRIX BC FOR RECTANGULAR GRID &&&&&&\n");

          // store SL
                    const Real Sl = pl/dxs; 
                    champParameters(4,side,axis,grid)=Sl;  // save this           
                    champParameters(5,side,axis,grid)=dxs; // save this           
                
          // Evaluate the (single component) Laplace operator for points on the boundary
                    realSerialArray xxCoeff(M0,Ib1,Ib2,Ib3), yyCoeff(M0,Ib1,Ib2,Ib3), xCoeff(M0,Ib1,Ib2,Ib3), yCoeff(M0,Ib1,Ib2,Ib3), idCoeff(M0,Ib1,Ib2,Ib3);
                    mgop.assignCoefficients(MappedGridOperators::xxDerivative,xxCoeff,Ib1,Ib2,Ib3,0,0); //
                    mgop.assignCoefficients(MappedGridOperators::yyDerivative,yyCoeff,Ib1,Ib2,Ib3,0,0); //
                    mgop.assignCoefficients(MappedGridOperators::xDerivative ,xCoeff, Ib1,Ib2,Ib3,0,0);
                    mgop.assignCoefficients(MappedGridOperators::yDerivative ,yCoeff, Ib1,Ib2,Ib3,0,0);
                    mgop.assignCoefficients(MappedGridOperators::identityOperator,idCoeff,Ib1,Ib2,Ib3,0,0);
                        
                            
                    Real maxDiff=0; // for testing compute any difference between coeff from Sijia and WDH

                    FOR_3D(i1,i2,i3,Ib1,Ib2,Ib3) // loop over points on the boundary
                    {
            //printF("n1=%g,n2=%g\n",normal(i1,i2,i3,0),normal(i1,i2,i3,1));
                        const real a0 = Sl;
                        const real a1 = theta + Sl*dxs*theta;
            // // ONLY VALID FOR CARTESIAN : axis==0 
            // const real a2 = dxs*( beta      ) + Sl*( .5*SQR(dxs)*beta      );
            // const real a3 = dxs*( (beta-1.) ) + Sl*( .5*SQR(dxs)*(beta-1.) );
            // const real a4 = a3; // for 3D 
                        real axx, ayy, azz;
                        if( axis==0 )
                        {
                            axx = (beta   )*( dxs + Sl*( .5*SQR(dxs) ) );
                            ayy = (beta-1.)*( dxs + Sl*( .5*SQR(dxs) ) );
                            azz = ayy; // for 3D 
                        }
                        else if( axis==1 )
                        {
                            ayy = (beta   )*( dxs + Sl*( .5*SQR(dxs) ) );
                            axx = (beta-1.)*( dxs + Sl*( .5*SQR(dxs) ) );
                            azz = axx; // for 3D       
                        }
                        assert( numberOfDimensions==2 );            
                            
                        i1m=i1-is1, i2m=i2-is2, i3m=i3-is3; //  ghost point is (i1m,i2m,i3m)

                        if( false )
                            printF(" FILL MATRIX BC FOR GHOST PT (i1m,i2m,i3m)=(%d,%d,%d) Sl=%g, dxs=%g, (heatFluxInterface)\n",i1m,i2m,i3m,Sl,dxs);

            // --- Coefficients from Sijia: 
            //#include "champ2InterfaceStencilRectangular.h"

                        if( orderOfAccuracy==2 )
                        {
                            #include "champ2InterfaceStencilRectangular.h"
                        }
                        else if( orderOfAccuracy==4 )
                        {
              //#include "champ4InterfaceStencilRectangular.h" 
                            Real hI=dxs;
                                  	Real dxl,dyl;
                                  	Real a00,a01,a02,a03,a04,a10,a11,a12,a13,a14,a20,a21,a22,a23,a24,a30,a31,a32,a33,a34,a40,a41,a42,a43,a44;
              	//printF("========== Running Sijia's 4th-order Rectangular script ==========\n");
                                  	int is=1;
                                  	int shftI4=2;
                                  	dxl = dx[0];
                                  	dyl = dx[1];
                                  	if (side==0){
                                      	    is=-1;
                                  	}
                                  	else{
                                      	    is=1;
                                  	}
              	//printF("     side1=%i, side2=%i, axis1=%i, axis2=%i, is=%i\n", side,side2,axis,axis2,is );
              	//printF("     theta=%g, beta=%g, Sl=%g\n", theta,beta,Sl);
              	//printF("     dxl=%g, dyl=%g, dxs=%g\n", dxl, dyl,dxs);
                                  	Real Dcc[5][5] = {0.,0.,0.,0.,0.,
                                                    	                  0.,0.,0.,0.,0.,
                                                    	                  0.,0.,1.,0.,0.,
                                                    	                  0.,0.,0.,0.,0.,
                                                    	                  0.,0.,0.,0.,0.}; // initialize
                                  	Real Dx4cc[5][5]={0.,0., 1./(12.*dxl),0.,0.,
                                                    	                  0.,0.,-8./(12.*dxl),0.,0.,
                                                    	                  0.,0.,0.,0.,0.,
                                                    	                  0.,0., 8./(12.*dxl),0.,0.,
                                                    	                  0.,0.,-1./(12.*dxl),0.,0.}; 
                                  	Real Dy4cc[5][5]={0.,0.,0.,0.,0.,
                                                    	                  0.,0.,0.,0.,0.,
                                                    	                  1./(12.*dyl),-8./(12.*dyl),0.,8./(12.*dyl),-1./(12.*dyl),
                                                    	                  0.,0.,0.,0.,0.,
                                                    	                  0.,0.,0.,0.,0.};
                                  	Real Dxx4cc[5][5]={0.,0., -1./(12.*SQR(dxl)),0.,0.,
                                                     	                   0.,0., 16./(12.*SQR(dxl)),0.,0.,
                                                     	                   0.,0.,-30./(12.*SQR(dxl)),0.,0.,
                                                     	                   0.,0., 16./(12.*SQR(dxl)),0.,0.,
                                                     	                   0.,0., -1./(12.*SQR(dxl)),0.,0.};
                                  	Real Dyy4cc[5][5]={0.,0.,0.,0.,0.,
                                                    	                  0.,0.,0.,0.,0.,
                                                    	                  -1./(12.*SQR(dyl)),16./(12.*SQR(dyl)),-30./(12.*SQR(dyl)),16./(12.*SQR(dyl)),-1./(12.*SQR(dyl)),
                                                    	                  0.,0.,0.,0.,0.,
                                                    	                  0.,0.,0.,0.,0.};
                                  	Real Dxxx2cc[5][5]={0.,0.,(-1./2.)/(pow(dxl,3)),0.,0.,
                                                     	                   0.,0.,       1./(pow(dxl,3)),0.,0.,
                                                     	                   0.,0.,0.,0.,0.,
                                                     	                   0.,0.,      -1./(pow(dxl,3)),0.,0.,
                                                     	                   0.,0.,  (1./2.)/(pow(dxl,3)),0.,0.};
                                  	Real Dyyy2cc[5][5]={0.,0.,0.,0.,0.,
                                                    	                  0.,0.,0.,0.,0.,
                                                    	                  (-1./2.)/(pow(dyl,3)),1./(pow(dyl,3)),0.,-1./(pow(dyl,3)),(1./2.)/(pow(dyl,3)),
                                                    	                  0.,0.,0.,0.,0.,
                                                    	                  0.,0.,0.,0.,0.}; 
                                  	Real Dxxxx2cc[5][5]={0.,0., 1./(pow(dxl,4)),0.,0.,
                                                       	                     0.,0.,-4./(pow(dxl,4)),0.,0.,
                                                       	                     0.,0., 6./(pow(dxl,4)),0.,0.,
                                                       	                     0.,0.,-4./(pow(dxl,4)),0.,0.,
                                                       	                     0.,0., 1./(pow(dxl,4)),0.,0.};
                                  	Real Dyyyy2cc[5][5]= {0.,0.,0.,0.,0.,
                                                        	                      0.,0.,0.,0.,0.,
                                                        	                      1./(pow(dyl,4)),-4./(pow(dyl,4)),6./(pow(dyl,4)) ,-4./(pow(dyl,4)),1./(pow(dyl,4)),
                                                        	                      0.,0.,0.,0.,0.,
                                                        	                      0.,0.,0.,0.,0.};
                                  	Real Dxy4cc[5][5]={ 1./(144.*dxl*dyl), -8./(144.*dxl*dyl),0.,  8./(144.*dxl*dyl),-1./(144.*dxl*dyl),
                                                     	                   -8./(144.*dxl*dyl), 64./(144.*dxl*dyl),0.,-64./(144.*dxl*dyl), 8./(144.*dxl*dyl),
                                                      	                    0.,0.,0.,0.,0.,
                                                      	                    8./(144.*dxl*dyl),-64./(144.*dxl*dyl),0., 64./(144.*dxl*dyl),-8./(144.*dxl*dyl),
                                                     	                   -1./(144.*dxl*dyl),  8./(144.*dxl*dyl),0., -8./(144.*dxl*dyl), 1./(144.*dxl*dyl)};
                                  	Real Dxxy4cc[5][5]={ -1./(144.*SQR(dxl)*dyl),   8./(144.*SQR(dxl)*dyl),0.,  -8./(144.*SQR(dxl)*dyl),  1./(144.*SQR(dxl)*dyl),
                                                       	                     16./(144.*SQR(dxl)*dyl),-128./(144.*SQR(dxl)*dyl),0., 128./(144.*SQR(dxl)*dyl),-16./(144.*SQR(dxl)*dyl),
                                                      	                    -30./(144.*SQR(dxl)*dyl), 240./(144.*SQR(dxl)*dyl),0.,-240./(144.*SQR(dxl)*dyl), 30./(144.*SQR(dxl)*dyl),
                                                       	                     16./(144.*SQR(dxl)*dyl),-128./(144.*SQR(dxl)*dyl),0., 128./(144.*SQR(dxl)*dyl),-16./(144.*SQR(dxl)*dyl),
                                                       	                     -1./(144.*SQR(dxl)*dyl),   8./(144.*SQR(dxl)*dyl),0.,  -8./(144.*SQR(dxl)*dyl),  1./(144.*SQR(dxl)*dyl)};
                                  	Real Dxyy4cc[5][5]={ -1./(144.*dxl*SQR(dyl)),  16./(144.*dxl*SQR(dyl)), -30./(144.*dxl*SQR(dyl)),  16./(144.*dxl*SQR(dyl)),-1./(144.*dxl*SQR(dyl)),
                                                        	                      8./(144.*dxl*SQR(dyl)),-128./(144.*dxl*SQR(dyl)), 240./(144.*dxl*SQR(dyl)),-128./(144.*dxl*SQR(dyl)), 8./(144.*dxl*SQR(dyl)),
                                                        	                      0.,0.,0.,0.,0.,
                                                       	                     -8./(144.*dxl*SQR(dyl)), 128./(144.*dxl*SQR(dyl)),-240./(144.*dxl*SQR(dyl)), 128./(144.*dxl*SQR(dyl)),-8./(144.*dxl*SQR(dyl)),
                                                        	                      1./(144.*dxl*SQR(dyl)), -16./(144.*dxl*SQR(dyl)),  30./(144.*dxl*SQR(dyl)), -16./(144.*dxl*SQR(dyl)), 1./(144.*dxl*SQR(dyl))};
                                  	Real Dxxyy4cc[5][5]={  1./(144.*SQR(dxl)*SQR(dyl)), -16./(144.*SQR(dxl)*SQR(dyl)),  30./(144.*SQR(dxl)*SQR(dyl)), -16./(144.*SQR(dxl)*SQR(dyl)),  1./(144.*SQR(dxl)*SQR(dyl)),
                                                       	                     -16./(144.*SQR(dxl)*SQR(dyl)), 256./(144.*SQR(dxl)*SQR(dyl)),-480./(144.*SQR(dxl)*SQR(dyl)), 256./(144.*SQR(dxl)*SQR(dyl)),-16./(144.*SQR(dxl)*SQR(dyl)),
                                                        	                      30./(144.*SQR(dxl)*SQR(dyl)),-480./(144.*SQR(dxl)*SQR(dyl)), 900./(144.*SQR(dxl)*SQR(dyl)),-480./(144.*SQR(dxl)*SQR(dyl)), 30./(144.*SQR(dxl)*SQR(dyl)),
                                                       	                     -16./(144.*SQR(dxl)*SQR(dyl)), 256./(144.*SQR(dxl)*SQR(dyl)),-480./(144.*SQR(dxl)*SQR(dyl)), 256./(144.*SQR(dxl)*SQR(dyl)),-16./(144.*SQR(dxl)*SQR(dyl)),
                                                         	                       1./(144.*SQR(dxl)*SQR(dyl)), -16./(144.*SQR(dxl)*SQR(dyl)),  30./(144.*SQR(dxl)*SQR(dyl)), -16./(144.*SQR(dxl)*SQR(dyl)),  1./(144.*SQR(dxl)*SQR(dyl))};
                                  	Real Dxyyy2cc[5][5]={-1./(24.*dxl*pow(dyl,3)),  2./(24.*dxl*pow(dyl,3)),0., -2./(24.*dxl*pow(dyl,3)), 1./(24.*dxl*pow(dyl,3)),
                                                        	                      8./(24.*dxl*pow(dyl,3)),-16./(24.*dxl*pow(dyl,3)),0., 16./(24.*dxl*pow(dyl,3)),-8./(24.*dxl*pow(dyl,3)),
                                                        	                      0.,0.,0.,0.,0.,
                                                       	                     -8./(24.*dxl*pow(dyl,3)), 16./(24.*dxl*pow(dyl,3)),0.,-16./(24.*dxl*pow(dyl,3)), 8./(24.*dxl*pow(dyl,3)),
                                                        	                      1./(24.*dxl*pow(dyl,3)), -2./(24.*dxl*pow(dyl,3)),0.,  2./(24.*dxl*pow(dyl,3)),-1./(24.*dxl*pow(dyl,3))};
                                  	Real Dxxxy2cc[5][5]={-1./(24.*pow(dxl,3)*dyl),  8./(24.*pow(dxl,3)*dyl),0., -8./(24.*pow(dxl,3)*dyl), 1./(24.*pow(dxl,3)*dyl),
                                                        	                      2./(24.*pow(dxl,3)*dyl),-16./(24.*pow(dxl,3)*dyl),0., 16./(24.*pow(dxl,3)*dyl),-2./(24.*pow(dxl,3)*dyl),
                                                        	                      0.,0.,0.,0.,0.,
                                                       	                     -2./(24.*pow(dxl,3)*dyl), 16./(24.*pow(dxl,3)*dyl),0.,-16./(24.*pow(dxl,3)*dyl), 2./(24.*pow(dxl,3)*dyl),
                                                        	                      1./(24.*pow(dxl,3)*dyl), -8./(24.*pow(dxl,3)*dyl),0.,  8./(24.*pow(dxl,3)*dyl),-1./(24.*pow(dxl,3)*dyl)};
                                  	#define Dc(i,j) Dcc[i+shftI4][j+shftI4]
                                  	#define Dx4c(i,j) Dx4cc[i+shftI4][j+shftI4]
                                  	#define Dy4c(i,j) Dy4cc[i+shftI4][j+shftI4]
                                  	#define Dxx4c(i,j) Dxx4cc[i+shftI4][j+shftI4]
                                  	#define Dyy4c(i,j) Dyy4cc[i+shftI4][j+shftI4]
                                  	#define Dxxx2c(i,j) Dxxx2cc[i+shftI4][j+shftI4]
                                  	#define Dyyy2c(i,j) Dyyy2cc[i+shftI4][j+shftI4]
                                  	#define Dxxxx2c(i,j) Dxxxx2cc[i+shftI4][j+shftI4]
                                  	#define Dyyyy2c(i,j) Dyyyy2cc[i+shftI4][j+shftI4]
                                  	#define Dxy4c(i,j) Dxy4cc[i+shftI4][j+shftI4]
                                  	#define Dxxy4c(i,j) Dxxy4cc[i+shftI4][j+shftI4]    
                                  	#define Dxyy4c(i,j) Dxyy4cc[i+shftI4][j+shftI4]    
                                  	#define Dxxyy4c(i,j) Dxxyy4cc[i+shftI4][j+shftI4]    
                                  	#define Dxxxy2c(i,j) Dxxxy2cc[i+shftI4][j+shftI4]    
                                  	#define Dxyyy2c(i,j) Dxyyy2cc[i+shftI4][j+shftI4]        
              	// MAPLE input 
                                  	a00  = Sl;
                                  	a01  = 0;
                                  	a02  =  (hI * (beta - 1) * (Sl * hI + 2)) / 0.2e1;
                                  	a03  = 0;
                                  	a04  = pow(hI, 0.3e1) *   pow( (beta - 1),  2) * (Sl * hI + 0.4e1) / 0.24e2;
                                  	a10  = is * theta * (Sl * hI + 1);
                                  	a11  = 0;
                                  	a12  =  (hI * hI * is * theta * (beta - 1) * (Sl * hI + 3)) / 0.6e1;
                                  	a13  = 0;
                                  	a14  = 0;
                                  	a20  =  (beta * (Sl * hI + 2) * hI) / 0.2e1;
                                  	a21  = 0;
                                  	a22  =  beta * pow(hI, 0.3e1) *  (beta - 1) * (Sl * hI + 0.4e1) / 0.12e2;
                                  	a23  = 0;
                                  	a24  = 0;
                                  	a30  =  (beta * is * theta * (Sl * hI + 3) * hI * hI) / 0.6e1;
                                  	a31  = 0;
                                  	a32  = 0;
                                  	a33  = 0;
                                  	a34  = 0;
                                  	a40  = pow(hI, 0.3e1) * beta * beta * (Sl * hI + 0.4e1) / 0.24e2;
                                  	a41  = 0;
                                  	a42  = 0;
                                  	a43  = 0;
                                  	a44  = 0;
                                  	/* 
                                  	printF("    a00=%g, a01=%g, a02=%g, a03=%g, a04=%g \n", a00,a01,a02,a03,a04);
                                  	printF("    a10=%g, a11=%g, a12=%g, a13=%g, a14=%g \n", a10,a11,a12,a13,a14);
                                  	printF("    a20=%g, a21=%g, a22=%g, a23=%g, a24=%g \n", a20,a21,a22,a23,a24);
                                  	printF("    a30=%g, a31=%g, a32=%g, a33=%g, a34=%g \n", a30,a31,a32,a33,a34);
                                  	printF("    a40=%g, a41=%g, a42=%g, a43=%g, a44=%g \n", a40,a41,a42,a43,a44);
                                  	*/
                                  	if (axis==0){
                                      	    for (int i = -2; i < 3; i++) {
                                          	        for (int j = -2; j < 3; j++) {
                                              	            coefA(i,j) = Dc(i,j)*a00+Dx4c(i,j)*a10+Dxx4c(i,j)*a20+Dxxx2c(i,j)*a30+Dxxxx2c(i,j)*a40	                                    +Dy4c(i,j)*a01+Dyy4c(i,j)*a02+Dyyy2c(i,j)*a03+Dyyyy2c(i,j)*a04	                                    +Dxxy4c(i,j)*a21+Dxyy4c(i,j)*a12	                                    +Dxy4c(i,j)*a11+Dxxyy4c(i,j)*a22	                                    +Dxxxy2c(i,j)*a31+Dxyyy2c(i,j)*a13;
                                          	        }
                                      	    } 
                                  	}
                                   	 else{
                                      	    for (int i = -2; i < 3; i++) {
                                          	        for (int j = -2; j < 3; j++) {
                                              	            coefA(j,i) = Dc(j,i)*a00+Dy4c(j,i)*a10+Dyy4c(j,i)*a20+Dyyy2c(j,i)*a30+Dyyyy2c(j,i)*a40	                                    +Dx4c(j,i)*a01+Dxx4c(j,i)*a02+Dxxx2c(j,i)*a03+Dxxxx2c(j,i)*a04	                                    +Dxyy4c(j,i)*a21+Dxxy4c(j,i)*a12	                                    +Dxy4c(j,i)*a11+Dxxyy4c(j,i)*a22	                                    +Dxyyy2c(j,i)*a31+Dxxxy2c(j,i)*a13;
                                          	        }
                                      	    } 
                                   	 }
                        }
            //#include "champ4InterfaceStencilRectangular.h"
                        if( orderOfAccuracy==4 && debug & 4 ) 
                        {
                            printF("SH: order=4: coefA=(%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"
                                          "                    %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"
                                          "                    %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"
                                          "                    %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"
                                          "                    %10.3e,%10.3e,%10.3e,%10.3e,%10.3e )\n",coefA(-2,-2),coefA(-1,-2),coefA(0,-2),coefA(1,-2),coefA(2,-2),coefA(-2,-1),coefA(-1,-1),coefA(0,-1),coefA(1,-1),coefA(2,-1),coefA(-2, 0),coefA(-1, 0),coefA(0, 0),coefA(1, 0),coefA(2, 0),coefA(-2, 1),coefA(-1, 1),coefA(0, 1),coefA(1, 1),coefA(2, 1),coefA(-2, 2),coefA(-1, 2),coefA(0, 2),coefA(1, 2),coefA(2, 2));
                        }


            //#include "champ2InterfaceStencil.C"
            //printF(" My coefficients: coefA =(%g,%g,%g,%g,%g,%g,%g,%g,%g)\n", coefA(-1,-1),coefA(0,-1),coefA(1,-1),coefA(-1,0),coefA(0,0),coefA(1,0),coefA(-1,1),coefA(0,1),coefA(1,1));
                        for( int e=eqnStart; e<=eqnEnd; e++ ) // equation eq1, eq2, ...
                        {
                            const int c=0; // component number 
                            ForStencil(m1,m2,m3)
                            {
                                int m  = M123(m1,m2,m3);        // the single-component coeff-index
                                int mm = M123CE(m1,m2,m3,c,e);  // the system coeff-index
                                
                                if( mapleOption==0 )
                                {
                                    coeff(mm,i1m,i2m,i3m) = a0*idCoeff(m,i1,i2,i3) 
                                                                                + a1*( normal(i1,i2,i3,0)*xCoeff(m,i1,i2,i3) + normal(i1,i2,i3,1)*yCoeff(m,i1,i2,i3) )
                                                                                +axx*xxCoeff(m,i1,i2,i3) + ayy*yyCoeff(m,i1,i2,i3); 
                                }
                                else if( mapleOption==1 )
                                {
                  //printF("---->Use MAPLE generated coefficients.\n");
                                    coeff(mm,i1m,i2m,i3m) = coefA(m1,m2);
                                }
                //coeff(mm,i1m,i2m,i3m) =  a0*idCoeff(m,i1,i2,i3) 
                //                       + a1*( normal(i1,i2,i3,0)*xCoeff(m,i1,i2,i3) + normal(i1,i2,i3,1)*yCoeff(m,i1,i2,i3) )
                //                        +axx*xxCoeff(m,i1,i2,i3) + ayy*yyCoeff(m,i1,i2,i3); 
                                

                                if( multiDomainProblem==0 )
                                { // debug info 
                                    printF("champBC: (i1,i2,i3)=(%3d,%3d,%3d) (m1,m2,m3)=(%3d,%3d,%3d) coeff=%10.3e,  coeff=%10.3e (new) diff=%9.2e\n",
                                              i1,i2,i3,m1,m2,m3,coeff(mm,i1m,i2m,i3m), coefA(m1,m2), coeff(mm,i1m,i2m,i3m)-coefA(m1,m2));
                                } 
                                if( true )
                                {
                  // for testing compute any difference between Sijia and WDH
                                    maxDiff = max(maxDiff, fabs(coeff(mm,i1m,i2m,i3m)-coefA(m1,m2)));
                                }                                       

                // coeff(mm,i1m,i2m,i3m) = 100; // **TEST ****

                // Specify that the above coeff value is the coefficient of component c at the grid point (j1,j2,j3).
                                j1=i1+m1, j2=i2+m2, j3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
                                setEquationNumber(mm, e,i1m,i2m,i3m,  c,j1,j2,j3 );  // macro to set equationNumber

                // ---- Fill Ghost 2 ----

                                if( numGhost>1 )
                                {
                                    const int ghost=2;
                                    if( fillSecondGhostUsingExtrapolation )
                                    {
                                            const int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                      // --- fill in the coefficients of the extrapolation formula ---
                                            for( int me=0; me<=extrapOrder; me++ )
                                            {
                                                coeff(me,i1m,i2m,i3m) = extrapCoeff[me];
                                                const int j1=i1m + me*is1, j2=i2m + me*is2, j3=i3m + me*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                                setEquationNumber(me, e,i1m,i2m,i3m,  c,j1,j2,j3 );             // macro to set equationNumber
                                            }                
                                    }
                                    else
                                    {
                    //fillGhostChamp(ghost,dxs); 

                                        const int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
                                        real hI=ghost*dxs;
                                              	Real dxl,dyl;
                                              	Real a00,a01,a02,a03,a04,a10,a11,a12,a13,a14,a20,a21,a22,a23,a24,a30,a31,a32,a33,a34,a40,a41,a42,a43,a44;
                    	//printF("========== Running Sijia's 4th-order Rectangular script ==========\n");
                                              	int is=1;
                                              	int shftI4=2;
                                              	dxl = dx[0];
                                              	dyl = dx[1];
                                              	if (side==0){
                                                  	    is=-1;
                                              	}
                                              	else{
                                                  	    is=1;
                                              	}
                    	//printF("     side1=%i, side2=%i, axis1=%i, axis2=%i, is=%i\n", side,side2,axis,axis2,is );
                    	//printF("     theta=%g, beta=%g, Sl=%g\n", theta,beta,Sl);
                    	//printF("     dxl=%g, dyl=%g, dxs=%g\n", dxl, dyl,dxs);
                                              	Real Dcc[5][5] = {0.,0.,0.,0.,0.,
                                                                	                  0.,0.,0.,0.,0.,
                                                                	                  0.,0.,1.,0.,0.,
                                                                	                  0.,0.,0.,0.,0.,
                                                                	                  0.,0.,0.,0.,0.}; // initialize
                                              	Real Dx4cc[5][5]={0.,0., 1./(12.*dxl),0.,0.,
                                                                	                  0.,0.,-8./(12.*dxl),0.,0.,
                                                                	                  0.,0.,0.,0.,0.,
                                                                	                  0.,0., 8./(12.*dxl),0.,0.,
                                                                	                  0.,0.,-1./(12.*dxl),0.,0.}; 
                                              	Real Dy4cc[5][5]={0.,0.,0.,0.,0.,
                                                                	                  0.,0.,0.,0.,0.,
                                                                	                  1./(12.*dyl),-8./(12.*dyl),0.,8./(12.*dyl),-1./(12.*dyl),
                                                                	                  0.,0.,0.,0.,0.,
                                                                	                  0.,0.,0.,0.,0.};
                                              	Real Dxx4cc[5][5]={0.,0., -1./(12.*SQR(dxl)),0.,0.,
                                                                 	                   0.,0., 16./(12.*SQR(dxl)),0.,0.,
                                                                 	                   0.,0.,-30./(12.*SQR(dxl)),0.,0.,
                                                                 	                   0.,0., 16./(12.*SQR(dxl)),0.,0.,
                                                                 	                   0.,0., -1./(12.*SQR(dxl)),0.,0.};
                                              	Real Dyy4cc[5][5]={0.,0.,0.,0.,0.,
                                                                	                  0.,0.,0.,0.,0.,
                                                                	                  -1./(12.*SQR(dyl)),16./(12.*SQR(dyl)),-30./(12.*SQR(dyl)),16./(12.*SQR(dyl)),-1./(12.*SQR(dyl)),
                                                                	                  0.,0.,0.,0.,0.,
                                                                	                  0.,0.,0.,0.,0.};
                                              	Real Dxxx2cc[5][5]={0.,0.,(-1./2.)/(pow(dxl,3)),0.,0.,
                                                                 	                   0.,0.,       1./(pow(dxl,3)),0.,0.,
                                                                 	                   0.,0.,0.,0.,0.,
                                                                 	                   0.,0.,      -1./(pow(dxl,3)),0.,0.,
                                                                 	                   0.,0.,  (1./2.)/(pow(dxl,3)),0.,0.};
                                              	Real Dyyy2cc[5][5]={0.,0.,0.,0.,0.,
                                                                	                  0.,0.,0.,0.,0.,
                                                                	                  (-1./2.)/(pow(dyl,3)),1./(pow(dyl,3)),0.,-1./(pow(dyl,3)),(1./2.)/(pow(dyl,3)),
                                                                	                  0.,0.,0.,0.,0.,
                                                                	                  0.,0.,0.,0.,0.}; 
                                              	Real Dxxxx2cc[5][5]={0.,0., 1./(pow(dxl,4)),0.,0.,
                                                                   	                     0.,0.,-4./(pow(dxl,4)),0.,0.,
                                                                   	                     0.,0., 6./(pow(dxl,4)),0.,0.,
                                                                   	                     0.,0.,-4./(pow(dxl,4)),0.,0.,
                                                                   	                     0.,0., 1./(pow(dxl,4)),0.,0.};
                                              	Real Dyyyy2cc[5][5]= {0.,0.,0.,0.,0.,
                                                                    	                      0.,0.,0.,0.,0.,
                                                                    	                      1./(pow(dyl,4)),-4./(pow(dyl,4)),6./(pow(dyl,4)) ,-4./(pow(dyl,4)),1./(pow(dyl,4)),
                                                                    	                      0.,0.,0.,0.,0.,
                                                                    	                      0.,0.,0.,0.,0.};
                                              	Real Dxy4cc[5][5]={ 1./(144.*dxl*dyl), -8./(144.*dxl*dyl),0.,  8./(144.*dxl*dyl),-1./(144.*dxl*dyl),
                                                                 	                   -8./(144.*dxl*dyl), 64./(144.*dxl*dyl),0.,-64./(144.*dxl*dyl), 8./(144.*dxl*dyl),
                                                                  	                    0.,0.,0.,0.,0.,
                                                                  	                    8./(144.*dxl*dyl),-64./(144.*dxl*dyl),0., 64./(144.*dxl*dyl),-8./(144.*dxl*dyl),
                                                                 	                   -1./(144.*dxl*dyl),  8./(144.*dxl*dyl),0., -8./(144.*dxl*dyl), 1./(144.*dxl*dyl)};
                                              	Real Dxxy4cc[5][5]={ -1./(144.*SQR(dxl)*dyl),   8./(144.*SQR(dxl)*dyl),0.,  -8./(144.*SQR(dxl)*dyl),  1./(144.*SQR(dxl)*dyl),
                                                                   	                     16./(144.*SQR(dxl)*dyl),-128./(144.*SQR(dxl)*dyl),0., 128./(144.*SQR(dxl)*dyl),-16./(144.*SQR(dxl)*dyl),
                                                                  	                    -30./(144.*SQR(dxl)*dyl), 240./(144.*SQR(dxl)*dyl),0.,-240./(144.*SQR(dxl)*dyl), 30./(144.*SQR(dxl)*dyl),
                                                                   	                     16./(144.*SQR(dxl)*dyl),-128./(144.*SQR(dxl)*dyl),0., 128./(144.*SQR(dxl)*dyl),-16./(144.*SQR(dxl)*dyl),
                                                                   	                     -1./(144.*SQR(dxl)*dyl),   8./(144.*SQR(dxl)*dyl),0.,  -8./(144.*SQR(dxl)*dyl),  1./(144.*SQR(dxl)*dyl)};
                                              	Real Dxyy4cc[5][5]={ -1./(144.*dxl*SQR(dyl)),  16./(144.*dxl*SQR(dyl)), -30./(144.*dxl*SQR(dyl)),  16./(144.*dxl*SQR(dyl)),-1./(144.*dxl*SQR(dyl)),
                                                                    	                      8./(144.*dxl*SQR(dyl)),-128./(144.*dxl*SQR(dyl)), 240./(144.*dxl*SQR(dyl)),-128./(144.*dxl*SQR(dyl)), 8./(144.*dxl*SQR(dyl)),
                                                                    	                      0.,0.,0.,0.,0.,
                                                                   	                     -8./(144.*dxl*SQR(dyl)), 128./(144.*dxl*SQR(dyl)),-240./(144.*dxl*SQR(dyl)), 128./(144.*dxl*SQR(dyl)),-8./(144.*dxl*SQR(dyl)),
                                                                    	                      1./(144.*dxl*SQR(dyl)), -16./(144.*dxl*SQR(dyl)),  30./(144.*dxl*SQR(dyl)), -16./(144.*dxl*SQR(dyl)), 1./(144.*dxl*SQR(dyl))};
                                              	Real Dxxyy4cc[5][5]={  1./(144.*SQR(dxl)*SQR(dyl)), -16./(144.*SQR(dxl)*SQR(dyl)),  30./(144.*SQR(dxl)*SQR(dyl)), -16./(144.*SQR(dxl)*SQR(dyl)),  1./(144.*SQR(dxl)*SQR(dyl)),
                                                                   	                     -16./(144.*SQR(dxl)*SQR(dyl)), 256./(144.*SQR(dxl)*SQR(dyl)),-480./(144.*SQR(dxl)*SQR(dyl)), 256./(144.*SQR(dxl)*SQR(dyl)),-16./(144.*SQR(dxl)*SQR(dyl)),
                                                                    	                      30./(144.*SQR(dxl)*SQR(dyl)),-480./(144.*SQR(dxl)*SQR(dyl)), 900./(144.*SQR(dxl)*SQR(dyl)),-480./(144.*SQR(dxl)*SQR(dyl)), 30./(144.*SQR(dxl)*SQR(dyl)),
                                                                   	                     -16./(144.*SQR(dxl)*SQR(dyl)), 256./(144.*SQR(dxl)*SQR(dyl)),-480./(144.*SQR(dxl)*SQR(dyl)), 256./(144.*SQR(dxl)*SQR(dyl)),-16./(144.*SQR(dxl)*SQR(dyl)),
                                                                     	                       1./(144.*SQR(dxl)*SQR(dyl)), -16./(144.*SQR(dxl)*SQR(dyl)),  30./(144.*SQR(dxl)*SQR(dyl)), -16./(144.*SQR(dxl)*SQR(dyl)),  1./(144.*SQR(dxl)*SQR(dyl))};
                                              	Real Dxyyy2cc[5][5]={-1./(24.*dxl*pow(dyl,3)),  2./(24.*dxl*pow(dyl,3)),0., -2./(24.*dxl*pow(dyl,3)), 1./(24.*dxl*pow(dyl,3)),
                                                                    	                      8./(24.*dxl*pow(dyl,3)),-16./(24.*dxl*pow(dyl,3)),0., 16./(24.*dxl*pow(dyl,3)),-8./(24.*dxl*pow(dyl,3)),
                                                                    	                      0.,0.,0.,0.,0.,
                                                                   	                     -8./(24.*dxl*pow(dyl,3)), 16./(24.*dxl*pow(dyl,3)),0.,-16./(24.*dxl*pow(dyl,3)), 8./(24.*dxl*pow(dyl,3)),
                                                                    	                      1./(24.*dxl*pow(dyl,3)), -2./(24.*dxl*pow(dyl,3)),0.,  2./(24.*dxl*pow(dyl,3)),-1./(24.*dxl*pow(dyl,3))};
                                              	Real Dxxxy2cc[5][5]={-1./(24.*pow(dxl,3)*dyl),  8./(24.*pow(dxl,3)*dyl),0., -8./(24.*pow(dxl,3)*dyl), 1./(24.*pow(dxl,3)*dyl),
                                                                    	                      2./(24.*pow(dxl,3)*dyl),-16./(24.*pow(dxl,3)*dyl),0., 16./(24.*pow(dxl,3)*dyl),-2./(24.*pow(dxl,3)*dyl),
                                                                    	                      0.,0.,0.,0.,0.,
                                                                   	                     -2./(24.*pow(dxl,3)*dyl), 16./(24.*pow(dxl,3)*dyl),0.,-16./(24.*pow(dxl,3)*dyl), 2./(24.*pow(dxl,3)*dyl),
                                                                    	                      1./(24.*pow(dxl,3)*dyl), -8./(24.*pow(dxl,3)*dyl),0.,  8./(24.*pow(dxl,3)*dyl),-1./(24.*pow(dxl,3)*dyl)};
                                              	#define Dc(i,j) Dcc[i+shftI4][j+shftI4]
                                              	#define Dx4c(i,j) Dx4cc[i+shftI4][j+shftI4]
                                              	#define Dy4c(i,j) Dy4cc[i+shftI4][j+shftI4]
                                              	#define Dxx4c(i,j) Dxx4cc[i+shftI4][j+shftI4]
                                              	#define Dyy4c(i,j) Dyy4cc[i+shftI4][j+shftI4]
                                              	#define Dxxx2c(i,j) Dxxx2cc[i+shftI4][j+shftI4]
                                              	#define Dyyy2c(i,j) Dyyy2cc[i+shftI4][j+shftI4]
                                              	#define Dxxxx2c(i,j) Dxxxx2cc[i+shftI4][j+shftI4]
                                              	#define Dyyyy2c(i,j) Dyyyy2cc[i+shftI4][j+shftI4]
                                              	#define Dxy4c(i,j) Dxy4cc[i+shftI4][j+shftI4]
                                              	#define Dxxy4c(i,j) Dxxy4cc[i+shftI4][j+shftI4]    
                                              	#define Dxyy4c(i,j) Dxyy4cc[i+shftI4][j+shftI4]    
                                              	#define Dxxyy4c(i,j) Dxxyy4cc[i+shftI4][j+shftI4]    
                                              	#define Dxxxy2c(i,j) Dxxxy2cc[i+shftI4][j+shftI4]    
                                              	#define Dxyyy2c(i,j) Dxyyy2cc[i+shftI4][j+shftI4]        
                    	// MAPLE input 
                                              	a00  = Sl;
                                              	a01  = 0;
                                              	a02  =  (hI * (beta - 1) * (Sl * hI + 2)) / 0.2e1;
                                              	a03  = 0;
                                              	a04  = pow(hI, 0.3e1) *   pow( (beta - 1),  2) * (Sl * hI + 0.4e1) / 0.24e2;
                                              	a10  = is * theta * (Sl * hI + 1);
                                              	a11  = 0;
                                              	a12  =  (hI * hI * is * theta * (beta - 1) * (Sl * hI + 3)) / 0.6e1;
                                              	a13  = 0;
                                              	a14  = 0;
                                              	a20  =  (beta * (Sl * hI + 2) * hI) / 0.2e1;
                                              	a21  = 0;
                                              	a22  =  beta * pow(hI, 0.3e1) *  (beta - 1) * (Sl * hI + 0.4e1) / 0.12e2;
                                              	a23  = 0;
                                              	a24  = 0;
                                              	a30  =  (beta * is * theta * (Sl * hI + 3) * hI * hI) / 0.6e1;
                                              	a31  = 0;
                                              	a32  = 0;
                                              	a33  = 0;
                                              	a34  = 0;
                                              	a40  = pow(hI, 0.3e1) * beta * beta * (Sl * hI + 0.4e1) / 0.24e2;
                                              	a41  = 0;
                                              	a42  = 0;
                                              	a43  = 0;
                                              	a44  = 0;
                                              	/* 
                                              	printF("    a00=%g, a01=%g, a02=%g, a03=%g, a04=%g \n", a00,a01,a02,a03,a04);
                                              	printF("    a10=%g, a11=%g, a12=%g, a13=%g, a14=%g \n", a10,a11,a12,a13,a14);
                                              	printF("    a20=%g, a21=%g, a22=%g, a23=%g, a24=%g \n", a20,a21,a22,a23,a24);
                                              	printF("    a30=%g, a31=%g, a32=%g, a33=%g, a34=%g \n", a30,a31,a32,a33,a34);
                                              	printF("    a40=%g, a41=%g, a42=%g, a43=%g, a44=%g \n", a40,a41,a42,a43,a44);
                                              	*/
                                              	if (axis==0){
                                                  	    for (int i = -2; i < 3; i++) {
                                                      	        for (int j = -2; j < 3; j++) {
                                                          	            coefAg(i,j) = Dc(i,j)*a00+Dx4c(i,j)*a10+Dxx4c(i,j)*a20+Dxxx2c(i,j)*a30+Dxxxx2c(i,j)*a40	                                    +Dy4c(i,j)*a01+Dyy4c(i,j)*a02+Dyyy2c(i,j)*a03+Dyyyy2c(i,j)*a04	                                    +Dxxy4c(i,j)*a21+Dxyy4c(i,j)*a12	                                    +Dxy4c(i,j)*a11+Dxxyy4c(i,j)*a22	                                    +Dxxxy2c(i,j)*a31+Dxyyy2c(i,j)*a13;
                                                      	        }
                                                  	    } 
                                              	}
                                               	 else{
                                                  	    for (int i = -2; i < 3; i++) {
                                                      	        for (int j = -2; j < 3; j++) {
                                                          	            coefAg(j,i) = Dc(j,i)*a00+Dy4c(j,i)*a10+Dyy4c(j,i)*a20+Dyyy2c(j,i)*a30+Dyyyy2c(j,i)*a40	                                    +Dx4c(j,i)*a01+Dxx4c(j,i)*a02+Dxxx2c(j,i)*a03+Dxxxx2c(j,i)*a04	                                    +Dxyy4c(j,i)*a21+Dxxy4c(j,i)*a12	                                    +Dxy4c(j,i)*a11+Dxxyy4c(j,i)*a22	                                    +Dxyyy2c(j,i)*a31+Dxxxy2c(j,i)*a13;
                                                      	        }
                                                  	    } 
                                               	 }
                                        coeff(mm,i1m,i2m,i3m) = coefAg(m1,m2);
                                        const int j1=i1m + mm*is1, j2=i2m + mm*is2, j3=i3m + mm*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                                        setEquationNumber(mm, e,i1m,i2m,i3m,  c,j1,j2,j3 ); 
                                    }
                                }
                                
                          }
                          if( orderOfAccuracy==2 && (maxDiff > 10.*REAL_EPSILON/(SQR(dxs)) ) )
                          {
                              printF("\n **** champBC:ERROR: rectangular: There is a difference between coeff from Sijia and WDH, maxDiff=%9.2e *** \n",maxDiff);
                              ForStencil(m1,m2,m3)  
                              {
                                  int mm = M123CE(m1,m2,m3,c,e);  // the system coeff-index
                                
                                  printF("champBC: (i1,i2,i3)=(%3d,%3d,%3d) (m1,m2,m3)=(%3d,%3d,%3d) coeff=%10.3e (WDH),  coeff=%10.3e (Sijia) diff=%9.2e\n",
                                              i1,i2,i3,m1,m2,m3,coeff(mm,i1m,i2m,i3m), coefA(m1,m2), coeff(mm,i1m,i2m,i3m)-coefA(m1,m2));
                              }             
               // OV_ABORT("Stop here for now");
                          }
                      }  // end for e

                    }   // end FOR_3D
                    if( debug & 8 )
                    {
                        printF("My ghost coefficients: coefAg=(%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n                               %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n                               %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n                               %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n                               %10.3e,%10.3e,%10.3e,%10.3e,%10.3e )\n",coefAg(-2,-2),coefAg(-1,-2),coefAg(0,-2),coefAg(1,-2),coefAg(2,-2),coefAg(-2,-1),coefAg(-1,-1),coefAg(0,-1),coefAg(1,-1),coefAg(2,-1),coefAg(-2, 0),coefAg(-1, 0),coefAg(0, 0),coefAg(1, 0),coefAg(2, 0),coefAg(-2, 1),coefAg(-1, 1),coefAg(0, 1),coefAg(1, 1),coefAg(2, 1),coefAg(-2, 2),coefAg(-1, 2),coefAg(0, 2),coefAg(1, 2),coefAg(2, 2));
                    }
                    if( true )
                        printF(">>>>>>> champBC::INFO max-diff in coeff between Sijia and WDH = %9.2e <<<<<<<<<\n",maxDiff);

                }
                else
                {

          // ---------------------- CURVILINEAR GRID --------------------------
                    printF("\n &&&&&& CHAMP BC: FILL MATRIX BC FOR CURVILINEAR GRID &&&&&&\n");

          // ******* WE NEED THE METRICS FROM THE OPPOSITE SIDE *******
          // getOppositeSideDomainInfo(); 

          // // Evaluate the (single component) Laplace operator for points on the boundary
          // realSerialArray lapCoeff(M0,Ib1,Ib2,Ib3), xCoeff(M0,Ib1,Ib2,Ib3), yCoeff(M0,Ib1,Ib2,Ib3), idCoeff(M0,Ib1,Ib2,Ib3);
          // // realSerialArray xxCoeff(M0,Ib1,Ib2,Ib3), yyCoeff(M0,Ib1,Ib2,Ib3);
          // // mgop.assignCoefficients(MappedGridOperators::xxDerivative,xxCoeff,Ib1,Ib2,Ib3,0,0); //
          // // mgop.assignCoefficients(MappedGridOperators::yyDerivative,yyCoeff,Ib1,Ib2,Ib3,0,0); //
          // mgop.assignCoefficients(MappedGridOperators::laplacianOperator,lapCoeff,Ib1,Ib2,Ib3,0,0); //

          // mgop.assignCoefficients(MappedGridOperators::xDerivative ,xCoeff, Ib1,Ib2,Ib3,0,0);
          // mgop.assignCoefficients(MappedGridOperators::yDerivative ,yCoeff, Ib1,Ib2,Ib3,0,0);
          // mgop.assignCoefficients(MappedGridOperators::identityOperator,idCoeff,Ib1,Ib2,Ib3,0,0);

          // // r1r2Coeff : left side
          // RealArray r1r2Coeff(M0,Ib1,Ib2,Ib3);
          // mgop.assignCoefficients(MappedGridOperators::r1r2Derivative,r1r2Coeff,Ib1,Ib2,Ib3,0,0); //


          // // --- Matrix elements from the opposite side: 
          // RealArray r2Coeff(M0,Jb1,Jb2,Jb3), r2r2Coeff(M0,Jb1,Jb2,Jb3);

          // // These are from the right-side: 
          // if( axis2==0 )
          // {
          //   mgop2.assignCoefficients(MappedGridOperators::r2Derivative  ,r2Coeff  ,Jb1,Jb2,Jb3,0,0); // tangential direction is r2 
          //   mgop2.assignCoefficients(MappedGridOperators::r2r2Derivative,r2r2Coeff,Jb1,Jb2,Jb3,0,0); //
          // }
          // else
          // {
          //   mgop2.assignCoefficients(MappedGridOperators::r1Derivative  ,r2Coeff  ,Jb1,Jb2,Jb3,0,0); // tangential direction is r1 
          //   mgop2.assignCoefficients(MappedGridOperators::r1r1Derivative,r2r2Coeff,Jb1,Jb2,Jb3,0,0); //            
          // }


                    OV_GET_SERIAL_ARRAY(real, mg.inverseVertexDerivative(),rxLeft );  // metric on left 
                    OV_GET_SERIAL_ARRAY(real,mg2.inverseVertexDerivative(),rxRight);  // metric on right
          //  The rx array is only 4d, here is a macro to make it look 5d 
                    #define rxL(i1,i2,i3,m1,m2)  rxLeft(i1,i2,i3,(m1)+numberOfDimensions*(m2))          
                    #define rxR(i1,i2,i3,m1,m2) rxRight(i1,i2,i3,(m1)+numberOfDimensions*(m2))          

                    RealArray NCoeff(M0), LCoeff(M0); // Matrix representations of operators Nlr and Llr (see champ4/notes)

          // put grid spacings (dr) in c-arrays for convenience
                    const Real drL[3]   = { mg.gridSpacing(0),mg.gridSpacing(1),mg.gridSpacing(2) };     // grid spacing on left 
                    const Real dr12L[3] = { 1./(2.*drL[0]), 1./(2.*drL[1]), 1./(2.*drL[2]) };            // for D0 

                    const Real drR[3]   = { mg2.gridSpacing(0),mg2.gridSpacing(1),mg2.gridSpacing(2) };  // grid spacing on right 
                    const Real dr12R[3] = { 1./(2.*drR[0]), 1./(2.*drR[1]), 1./(2.*drR[2]) };            // for D0 


          // Real uDotGradLr=0., uDotGradLs=0.;
          // Real uDotGradRr=0., uDotGradRs=0.;

          // assert( axis1==0 && axis2==0 ); // for now we assume this

                    Real Sl = -1.; // Set below, first time through the loops 

          // --- For twilightzone we save some coefficients ----
          // printF("Save some coefficients for TZ\n");

                    realArray *ccPointer= & NCoeff; // not used by default
                    if( multiDomainProblem )
                    {
                        if( twilightZoneFlow )
                        {
                            bool sameSide=true; 
                            GridFaceDescriptor & myGFD = getInterfaceGridFaceDescriptor( grid, side, axis, parameters, sameSide );
                            if( !myGFD.dbase.has_key("ccChamp") )
                                myGFD.dbase.put<RealArray>("ccChamp");

              // -- for TZ we save some coefficients in the champ formula ---
                            RealArray & cc = myGFD.dbase.get<RealArray>("ccChamp");
                            const int numCC = orderOfAccuracy==2 ? 6 : 15;
                            Index Id1, Id2, Id3;
                            getBoundaryIndex(mg.gridIndexRange(),side,axis,Id1,Id2,Id3); // dimension this size for applyBC
                            cc.redim(Id1,Id2,Id3,numCC); 

                            cc=0.;           
                            ccPointer = &cc; 
                        }
                    }
                    else
                    { // single domain testing case 
                        if( twilightZoneFlow )
                        {
                            RealArray & ccChamp = parameters.dbase.put<RealArray>("ccChamp");
                            ccChamp.redim(Ib1,Ib2,Ib3,6); 
                            ccPointer = &ccChamp;
                        }
                    }
          // RealArray & cc = twilightZoneFlow ? myGFD.dbase.get<RealArray>("ccChamp") : NCoeff;
                    RealArray & cc = *ccPointer;

                    Real maxDiff=0., maxCoeff=0.; // for testing compute any difference between coeff from Sijia and WDH

                    bool fillMatrixWDH=false; 
                    if( orderOfAccuracy==2 || 
                            orderOfAccuracy==4 )
                    {
            // *new* way April 10, 2022
                        fillMatrixWDH=true; 

            // assert( axis==0 && axis2==0 ); // ** FIX ME***

            // Real Sl=-1.;
            // const Real hR = (1-2*side2)*drR[axis2]; // scaled grid spacing on right
            // if( Sl<0 )
            // {
            //   // We should probably always take Sl to be positive and consitently take the normal from the current side.


            //   // Define Sl from the first grid point encountered   *** Do this for now ***  **** FIX ME ***
            //   Sl = pl/fabs(hR);  // note hR may be negative 
            //   // Sl = pl/hR;  // note hR may be negative 
            //   // champParameters(4,side,axis,grid)=Sl; // save Sl
            // }

                        if( axis==0 && axis2==0 )
                        {
                            if( orderOfAccuracy==2 )
                            {
// Coefficients for CHAMP curvilinear order 2, NORMAL DIRECTION=0.
// File written by champGenCode.maple
const Real KLR = theta;
const Real DLR = beta;
printF("WDH: KLR=%9.3e, DLR=%9.3e\n", KLR,DLR);
// Derivatives of the metrics: order 4
// Derivatives of the metrics: order 2
// Fourth-order one-sided approximations
// Here are fully one-sided: *check me* *wdh* April 15, 2022
// Macro : fill end ghost points on b1, b2, by periodicity or extrapolation
                  
// Macro : fill end ghost points on crr, crs, ... by periodicity or extrapolation
                  
const int axis1L=axis,  axis2L=(axis1L+1) % numberOfDimensions;
const int axis1R=axis2, axis2R=(axis1R+1) % numberOfDimensions;
const IntegerArray & gir1 = mg.indexRange();
const IntegerArray & gir2 = mg2.indexRange();
int ksv[3], &ks1=ksv[0], &ks2=ksv[1], &ks3=ksv[2];
int axist;
//  ----- We need b1,b2 at 4 extra tangential ghosts (order=4)---- 
//  since L1 is evaluated at 4 extra ghost 
int extra=2*numGhost;
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray n1L(Ib1,Ib2,Ib3), n2L(Ib1,Ib2,Ib3)  ; // save the normal, it needs to be extended too;
RealArray b1L(Ib1,Ib2,Ib3), b2L(Ib1,Ib2,Ib3);
RealArray b1R(Jb1,Jb2,Jb3), b2R(Jb1,Jb2,Jb3);
// Evaluate the normal derivative coefficients at points on the boundary. Use right-normal = -left-normal.
extra=numGhost;  // we can only directly evaluate 2 ghost using metrics 
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
FOR_3IJD(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    n1L(i1,i2,i3) = normal(i1,i2,i3,0);
    n2L(i1,i2,i3) = normal(i1,i2,i3,1);
    b1L(i1,i2,i3) = -( normal(i1,i2,i3,0)*rxL(i1,i2,i3,0,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,0,1)) ;
    b2L(i1,i2,i3) = -( normal(i1,i2,i3,0)*rxL(i1,i2,i3,1,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,1,1)) ;
    b1R(j1,j2,j3) = -( normal(i1,i2,i3,0)*rxR(j1,j2,j3,0,0) + normal(i1,i2,i3,1)*rxR(j1,j2,j3,0,1)) ;
    b2R(j1,j2,j3) = -( normal(i1,i2,i3,0)*rxR(j1,j2,j3,1,0) + normal(i1,i2,i3,1)*rxR(j1,j2,j3,1,1)) ;
} // end FOR_3IJD

// Fill 3rd and 4th ghost in tangential directions for b1, b2 
const int ghostToFillStart = numGhost+1;
const int ghostToFillEnd   = 2*numGhost;
    axist = (axis + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir1(0,0),gir1(0,0)); 
        I2=Range(gir1(0,1),gir1(0,1)); 
        I3=Range(gir1(0,2),gir1(0,2)); 
        Iv[axis]=Range(gir1(side,axis),gir1(side,axis));
        Iv[axist]=gir1(sidet,axist);
        for( int ghost=ghostToFillStart; ghost<=ghostToFillEnd; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir1(1,0)-gir1(0,0))*ks1;
                Index Ip2 = Ig2 + (gir1(1,1)-gir1(0,1))*ks2;
                Index Ip3 = Ig3 + (gir1(1,2)-gir1(0,2))*ks3;
                n1L(Ig1,Ig2,Ig3) = n1L(Ip1,Ip2,Ip3); 
                n2L(Ig1,Ig2,Ig3) = n2L(Ip1,Ip2,Ip3); 
            }
            else;
            {
                n1L(Ig1,Ig2,Ig3) = (5.*n1L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*n1L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*n1L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*n1L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+n1L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                n2L(Ig1,Ig2,Ig3) = (5.*n2L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*n2L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*n2L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*n2L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+n2L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t
    axist = (axis + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir1(0,0),gir1(0,0)); 
        I2=Range(gir1(0,1),gir1(0,1)); 
        I3=Range(gir1(0,2),gir1(0,2)); 
        Iv[axis]=Range(gir1(side,axis),gir1(side,axis));
        Iv[axist]=gir1(sidet,axist);
        for( int ghost=ghostToFillStart; ghost<=ghostToFillEnd; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir1(1,0)-gir1(0,0))*ks1;
                Index Ip2 = Ig2 + (gir1(1,1)-gir1(0,1))*ks2;
                Index Ip3 = Ig3 + (gir1(1,2)-gir1(0,2))*ks3;
                b1L(Ig1,Ig2,Ig3) = b1L(Ip1,Ip2,Ip3); 
                b2L(Ig1,Ig2,Ig3) = b2L(Ip1,Ip2,Ip3); 
            }
            else;
            {
                b1L(Ig1,Ig2,Ig3) = (5.*b1L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b1L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b1L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b1L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b1L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                b2L(Ig1,Ig2,Ig3) = (5.*b2L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b2L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b2L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b2L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b2L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t
    axist = (axis2 + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir2(0,0),gir2(0,0)); 
        I2=Range(gir2(0,1),gir2(0,1)); 
        I3=Range(gir2(0,2),gir2(0,2)); 
        Iv[axis2]=Range(gir2(side2,axis2),gir2(side2,axis2));
        Iv[axist]=gir2(sidet,axist);
        for( int ghost=ghostToFillStart; ghost<=ghostToFillEnd; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg2.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir2(1,0)-gir2(0,0))*ks1;
                Index Ip2 = Ig2 + (gir2(1,1)-gir2(0,1))*ks2;
                Index Ip3 = Ig3 + (gir2(1,2)-gir2(0,2))*ks3;
                b1R(Ig1,Ig2,Ig3) = b1R(Ip1,Ip2,Ip3); 
                b2R(Ig1,Ig2,Ig3) = b2R(Ip1,Ip2,Ip3); 
            }
            else;
            {
                b1R(Ig1,Ig2,Ig3) = (5.*b1R(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b1R(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b1R(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b1R(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b1R(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                b2R(Ig1,Ig2,Ig3) = (5.*b2R(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b2R(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b2R(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b2R(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b2R(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t

// --- Set the optimized Schwarz parameter: Sl = pl/dx (now that we know b1,b2)  ---- 
//  For now set a constant value on the face 
if( Sl<0 )
{
    j1=Jb1.getBase(); j2=Jb2.getBase(); j3=Jb3.getBase(); // take b1R from this point 
    const Real dnR = fabs( drR[axis2]/b1R(j1,j2,j3) );    // approx. normal grid spacing on right
    Sl = pl/dnR; 
    champParameters(4,side,axis,grid)=Sl; // save Sl
}
extra=0;
getBoundaryIndex(mg.indexRange(),  side,axis, I1,I2,I3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3,extra);
int extraNormal=numGhost-1;
Iv[axis1L]=Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]=Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
RealArray crrL(I1,I2,I3), crsL(I1,I2,I3), cssL(I1,I2,I3), crL(I1,I2,I3), csL(I1,I2,I3);
RealArray crrR(J1,J2,J3), crsR(J1,J2,J3), cssR(J1,J2,J3), crR(J1,J2,J3), csR(J1,J2,J3);
RealArray rxxL(3,3), rxxR(3,3);
RealArray rxyL(3,3), rxyR(3,3);
// We can fill in crr, crs, and ss to 2 ghost since they do not involve terms using rxx, sxx, etc. 
extra=0;
getBoundaryIndex(mg.indexRange(),  side,axis, I1,I2,I3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3,extra);
Iv[axis1L]=Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]=Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
FOR_3IJ(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3)
{
  // Evaluate coefficients of the Laplacian at points near the boundary
  // uxx = (rx*Dr + sx*Ds)*[ (rx)*ur + sx*us ]
  // uxx = (rx)^2 urr + 2*(rx*sx)*urs + (sx)^2 *uss + (rxx)*ur + (sxx)*us
  // uyy = (ry)^2 urr + 2*(ry*sy)*urs + (sy)^2 *uss + (ryy)*ur + (syy)*us
    crrL(i1,i2,i3) = SQR(rxL(i1,i2,i3,0,0)) + SQR(rxL(i1,i2,i3,0,1));                                // (rx)^2 + (ry)^2 
    crsL(i1,i2,i3) = 2.*(rxL(i1,i2,i3,0,0)*rxL(i1,i2,i3,1,0) + rxL(i1,i2,i3,0,1)*rxL(i1,i2,i3,1,1)); // 2( rx*sx + ry*sy ) 
    cssL(i1,i2,i3) = SQR(rxL(i1,i2,i3,1,0)) + SQR(rxL(i1,i2,i3,1,1));                                // (sx)^2 + (sy)^2
  // These are done below
  //crL(i1,i2,i3)  = rxxL(0,0) + rxyL(0,1);  // rxx + ryy
  //csL(i1,i2,i3)  = rxxL(1,0) + rxyL(1,1);  // sxx + syy
  // Evaluate coefficients of the Laplacian at points near the boundary
    crrR(j1,j2,j3) = SQR(rxR(j1,j2,j3,0,0)) + SQR(rxR(j1,j2,j3,0,1));
    crsR(j1,j2,j3) = 2.*( rxR(j1,j2,j3,0,0)*rxR(j1,j2,j3,1,0) + rxR(j1,j2,j3,0,1)*rxR(j1,j2,j3,1,1) );
    cssR(j1,j2,j3) = SQR(rxR(j1,j2,j3,1,0)) + SQR(rxR(j1,j2,j3,1,1));
  // These are done below
  //crR(j1,j2,j3)  = rxxR(0,0) + rxyR(0,1);  // rxx + ryy
  //csR(j1,j2,j3)  = rxxR(1,0) + rxyR(1,1);  // sxx + syy
  // printF("WDH: crrL=%9.3e, crsL=%9.3e, cssL=%9.3e\n",crrL(i1,i2,i3),crsL(i1,i2,i3),cssL(i1,i2,i3) );
  // printF("WDH: crrR=%9.3e, crsR=%9.3e, cssR=%9.3e\n",crrR(j1,j2,j3),crsR(j1,j2,j3),cssR(j1,j2,j3) );
} // end FOR_3IJ

// We can only compute cr and cs up to the boundary since we need two ghost points to evaluate rxx, rxy, sxx, sxy
getBoundaryIndex(mg.indexRange(),side,axis,I1,I2,I3);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3);
Iv[axis1L]= Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]= Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
Real rxr,rxs;
FOR_3IJ(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3)
{
    for( int m2=0; m2<numberOfDimensions; m2++ )
    {
        for( int m1=0; m1<numberOfDimensions; m1++ )
        {
        // Evaluate derivatives of metrics to order 2 .. we are computing more than we need here *fix me*
                rxxL(m1,m2) = rxL(i1,i2,i3,0,0)*(-rxL(i1-1,i2,i3,m1,m2)+rxL(i1+1,i2,i3,m1,m2))/(2.*drL[0])+rxL(i1,i2,i3,1,0)*(-rxL(i1,i2-1,i3,m1,m2)+rxL(i1,i2+1,i3,m1,m2))/(2.*drL[1]);
                rxyL(m1,m2) = rxL(i1,i2,i3,0,1)*(-rxL(i1-1,i2,i3,m1,m2)+rxL(i1+1,i2,i3,m1,m2))/(2.*drL[0])+rxL(i1,i2,i3,1,1)*(-rxL(i1,i2-1,i3,m1,m2)+rxL(i1,i2+1,i3,m1,m2))/(2.*drL[1]);
                                                                                                            
                rxxR(m1,m2) = rxR(j1,j2,j3,0,0)*(-rxR(j1-1,j2,j3,m1,m2)+rxR(j1+1,j2,j3,m1,m2))/(2.*drR[0])+rxR(j1,j2,j3,1,0)*(-rxR(j1,j2-1,j3,m1,m2)+rxR(j1,j2+1,j3,m1,m2))/(2.*drR[1]);
                rxyR(m1,m2) = rxR(j1,j2,j3,0,1)*(-rxR(j1-1,j2,j3,m1,m2)+rxR(j1+1,j2,j3,m1,m2))/(2.*drR[0])+rxR(j1,j2,j3,1,1)*(-rxR(j1,j2-1,j3,m1,m2)+rxR(j1,j2+1,j3,m1,m2))/(2.*drR[1]);
        }
    }
  // Evaluate coefficients of the Laplacian at points near the boundary
  // uxx = (rx*Dr + sx*Ds)*[ (rx)*ur + sx*us ]
  // uxx = (rx)^2 urr + 2*(rx*sx)*urs + (sx)^2 *uss + (rxx)*ur + (sxx)*us
  // uyy = (ry)^2 urr + 2*(ry*sy)*urs + (sy)^2 *uss + (ryy)*ur + (syy)*us
    crL(i1,i2,i3)  = rxxL(0,0) + rxyL(0,1);  // rxx + ryy
    csL(i1,i2,i3)  = rxxL(1,0) + rxyL(1,1);  // sxx + syy
  // Evaluate coefficients of the Laplacian at points near the boundary
    crR(j1,j2,j3)  = rxxR(0,0) + rxyR(0,1);  // rxx + ryy
    csR(j1,j2,j3)  = rxxR(1,0) + rxyR(1,1);  // sxx + syy
    if( advectionIsOn )  // include advection terms scaled by 1/D : (u/D).grad 
    {
            if( variableAdvection )
            {
                u1DL = advectVarL(i1,i2,i3,0)/DL;
                u2DL = advectVarL(i1,i2,i3,1)/DL;
                u1DR = advectVarR(j1,j2,j3,0)/DR;
                u2DR = advectVarR(j1,j2,j3,1)/DR;
            }
        crL(i1,i2,i3) = crL(i1,i2,i3) - ( u1DL*rxL(i1,i2,i3,0,0) + u2DL*rxL(i1,i2,i3,0,1) ); 
        csL(i1,i2,i3) = csL(i1,i2,i3) - ( u1DL*rxL(i1,i2,i3,1,0) + u2DL*rxL(i1,i2,i3,1,1) ); 
        crR(j1,j2,j3) = crR(j1,j2,j3) - ( u1DR*rxR(j1,j2,j3,0,0) + u2DR*rxR(j1,j2,j3,0,1) ); 
        csR(j1,j2,j3) = csR(j1,j2,j3) - ( u1DR*rxR(j1,j2,j3,1,0) + u2DR*rxR(j1,j2,j3,1,1) ); 
    }
  // printF("WDH: (i1,i2)=(%4d,%4d) crL=%9.3e csL=%9.3e\n",i1,i2,crL(i1,i2,i3),csL(i1,i2,i3) );
  // printF("WDH: (j1,j2)=(%4d,%4d) crR=%9.3e csR=%9.3e\n",j1,j2,crR(j1,j2,j3),csR(j1,j2,j3) );
} // end FOR_3IJ


// ----- Evaluate coefficients of L1 on the boundary -----

// 2*numGhost = 2 extra ghost needed for 2nd order, true ? 
extra=2*numGhost;
getBoundaryIndex(mg.indexRange(), side ,axis ,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L1ICoeff(Ib1,Ib2,Ib3), L1ICoeffs(Ib1,Ib2,Ib3);
RealArray L1rCoeff(Ib1,Ib2,Ib3), L1rCoeffs(Ib1,Ib2,Ib3);
RealArray L1sCoeff(Ib1,Ib2,Ib3), L1sCoeffs(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    const Real n1R = -n1L(i1,i2,i3), n2R = -n2L(i1,i2,i3);    // for advection, use right normal 
        if( variableAdvection )
        {
            u1DL = advectVarL(i1,i2,i3,0)/DL;
            u2DL = advectVarL(i1,i2,i3,1)/DL;
            u1DR = advectVarR(j1,j2,j3,0)/DR;
            u2DR = advectVarR(j1,j2,j3,1)/DR;
        }
    L1ICoeff(i1,i2,i3) = ((-1.*KLR*u2DL+1.*u2DR)*n2R+(-1.*KLR*u1DL+1.*u1DR)*n1R)/b1R(j1,j2,j3);
    L1rCoeff(i1,i2,i3) = 1.*KLR*b1L(i1,i2,i3)/b1R(j1,j2,j3);
    L1sCoeff(i1,i2,i3) = (1.*KLR*b2L(i1,i2,i3)-1.*b2R(j1,j2,j3))/b1R(j1,j2,j3);
} // end FOR_3D

// ----- Evaluate coefficients of L2 on the boundary -----
extra=0;
getBoundaryIndex(mg.indexRange(), side,  axis,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L2rrCoeff(Ib1,Ib2,Ib3);
RealArray L2rsCoeff(Ib1,Ib2,Ib3);
RealArray L2ssCoeff(Ib1,Ib2,Ib3);
RealArray L2rCoeff(Ib1,Ib2,Ib3);
RealArray L2sCoeff(Ib1,Ib2,Ib3);
RealArray L2ICoeff(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  // Compute tangential derivatives of L1 - label tangential directions a and b (for 3D)
    L1rCoeffs(i1,i2,i3) = ((L1rCoeff(i1,i2+1,i3)-L1rCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L1sCoeffs(i1,i2,i3) = ((L1sCoeff(i1,i2+1,i3)-L1sCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L1ICoeffs(i1,i2,i3) = ((L1ICoeff(i1,i2+1,i3)-L1ICoeff(i1,i2-1,i3))/(2.*drL[1]));
  // Compute L2 coefficients
    L2rrCoeff(i1,i2,i3) = 1.*DLR*crrL(i1,i2,i3)/crrR(j1,j2,j3);
    L2rsCoeff(i1,i2,i3) = (1.*DLR*crsL(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1rCoeff(i1,i2,i3))/crrR(j1,j2,j3);
    L2ssCoeff(i1,i2,i3) = 1.*(DLR*cssL(i1,i2,i3)-crsR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-cssR(j1,j2,j3))/crrR(j1,j2,j3);
    L2rCoeff(i1,i2,i3) = 1.*(DLR*crL(i1,i2,i3)-crR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-crsR(j1,j2,j3)*L1rCoeffs(i1,i2,i3))/crrR(j1,j2,j3);
    L2sCoeff(i1,i2,i3) = ((-1.*L1ICoeff(i1,i2,i3)-1.*L1sCoeffs(i1,i2,i3))*crsR(j1,j2,j3)+1.*DLR*csL(i1,i2,i3)-1.*crR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*csR(j1,j2,j3))/crrR(j1,j2,j3);
    L2ICoeff(i1,i2,i3) = (-1.*crR(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1ICoeffs(i1,i2,i3))/crrR(j1,j2,j3);
} // end FOR_3D

// ----- Define coefficients in difference operators -----
Range R4(-1,1);
RealArray iCoeff(R4);
RealArray rCoeff(R4), rrCoeff(R4), rrrCoeff(R4), rrrrCoeff(R4);
RealArray sCoeff(R4), ssCoeff(R4), sssCoeff(R4), ssssCoeff(R4);
      iCoeff(-1)=  0.;    iCoeff(0)=  1.;    iCoeff(1)= 0.;   
      rCoeff(-1)= -1.;    rCoeff(0)=  0.;    rCoeff(1)= 1.;    rCoeff    /=(2.*drL[0])  ;
    rrCoeff(-1)=  1.;   rrCoeff(0)= -2.;   rrCoeff(1)= 1.;    rrCoeff   /=(SQR(drL[0]));
      sCoeff(-1)= -1.;    sCoeff(0)=  0.;    sCoeff(1)= 1.;    sCoeff    /=(2.*drL[1])  ;
    ssCoeff(-1)=  1.;   ssCoeff(0)= -2.;   ssCoeff(1)= 1.;    ssCoeff   /=(SQR(drL[1]));

// ----- Fill in the Matrix Coefficients for CHAMP -----
const Real dxs = (1-2*side2)*drR[axis2];
const Real h = dxs;
printF("WDH: grid=%d, (side,axis)=(%d,%d) (side2,axis2)=(%d,%d) dxs=%9.3e, Sl=%9.3e\n",grid, side,axis, side2,axis2, dxs,Sl);
RealArray coeff4(R4,R4);
const int e=0, c=0; // eqn number and component number
// ------- FILL CHAMP conditions into the matrix -----
// NOTE: skip top boundary if periodic (use indexRange) 
// NOTE: skip adjacent boundaries if Dirichlet BC *finish me**
extra=0; 
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
int axisp = (axis+1) % numberOfDimensions;
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    i1m=i1-is1, i2m=i2-is2, i3m=i3-is3; //  ghost point is (i1m,i2m,i3m)
    const Real bn = b1R(j1,j2,j3), bt = b2R(j1,j2,j3);  // b in normal and tangential directions
  // The next macro defines cI, cr, cs, crr, ...
  // ------  ORDER = 2 , NormalDirection = 0,------
    const Real h2By2=h*h/2., h3By6=h2By2*h/3., h4By24=h3By6*h/4.;
    Real L1IC = L1ICoeff(i1,i2,i3);
    Real L1ICs = L1ICoeffs(i1,i2,i3);
    Real L1rC = L1rCoeff(i1,i2,i3);
    Real L1rCs = L1rCoeffs(i1,i2,i3);
    Real L1sC = L1sCoeff(i1,i2,i3);
    Real L1sCs = L1sCoeffs(i1,i2,i3);
    Real L2IC = L2ICoeff(i1,i2,i3);
    Real L2rC = L2rCoeff(i1,i2,i3);
    Real L2sC = L2sCoeff(i1,i2,i3);
    Real L2rrC = L2rrCoeff(i1,i2,i3);
    Real L2rsC = L2rsCoeff(i1,i2,i3);
    Real L2ssC = L2ssCoeff(i1,i2,i3);
    Real cI  = (Sl*L1IC-bn*L2IC-bt*L1ICs)*h+Sl*h2By2*L2IC-bn*L1IC+Sl;
    Real cr  = (Sl*L1rC-bn*L2rC-bt*L1rCs)*h+Sl*h2By2*L2rC-bn*L1rC;
    Real cs  = ((-L1sCs-L1IC)*h-L2IC*h2By2-2)*bt+(Sl*L1sC-bn*L2sC)*h+Sl*h2By2*L2sC-bn*L1sC;
    Real crr  = L2rrC*(Sl*h2By2-bn*h);
    Real crs  = (-bn*L2rsC-bt*L1rC)*h+(Sl*L2rsC-bt*L2rC)*h2By2;
    Real css  = (-bn*L2ssC-bt*L1sC)*h+(Sl*L2ssC-bt*L2sC)*h2By2;
//  if( orderOfAccuracy==4 )
//  {
//    printF("domain2=%d: cI   =%12.4e, cr   =%12.4e, cs   =%12.4e, crr  =%12.4e, crs=%12.4e, css=%12.4e\n",domain2,cI,cr,cs,crr,crs,css);
//    printF("            crrr =%12.4e, crrs =%12.4e, crss =%12.4e, csss =%12.4e\n",crrr,crrs,crss,csss);
//    printF("            crrrr=%12.4e, crrrs=%12.4e, crrss=%12.4e, crsss=%12.4e, cssss=%12.4e\n",crrrr,crrrs,crrss,crsss,cssss);
//    printF("            L4ssssC=%12.4e, (Sl*h4By24-bn*h3By6)=%12.4e h4By24=%12.4e h3By6=%12.4e, Sl=%12.4e, bn=%12.4e, h=%12.4e, bt=%12.4e\n",L4ssssC,(Sl*h4By24-bn*h3By6),h4By24,h3By6,Sl,bn,h,bt);
//    if( domain2==0 )
//    {
//      OV_ABORT("stop here for now");
//    }
//  }
    ForStencil(m1,m2,m3)
    {
        int m  = M123(m1,m2,m3);        // the single-component coeff-index
        int mm = M123CE(m1,m2,m3,c,e);  // the system coeff-index
        coeff4(m1,m2) = 
                          + cI   *   iCoeff(m1) * iCoeff(m2) 
                          + cr   *   rCoeff(m1) * iCoeff(m2) 
                          + cs   *   iCoeff(m1) * sCoeff(m2) 
                          + crr  *  rrCoeff(m1) * iCoeff(m2) 
                          + crs  *   rCoeff(m1) * sCoeff(m2) 
                          + css  *   iCoeff(m1) *ssCoeff(m2) 
                          ;
        if( fillMatrixWDH )
        {
            coeff(mm,i1m,i2m,i3m) = coeff4(m1,m2);
      // Specify that the above coeff value is the coefficient of component c at the grid point (j1,j2,j3).
            const int k1=i1+m1, k2=i2+m2, k3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)    
            setEquationNumber(mm, e,i1m,i2m,i3m,  c,k1,k2,k3 );  // macro to set equationNumber                  
                                                                                                                                                                                                                      
       // Fill Ghost 2 -- extrapolation for now                                                            
              const int ghost=2;                                                                                  
                  const int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
         // --- fill in the coefficients of the extrapolation formula ---
                  for( int me=0; me<=extrapOrder; me++ )
                  {
                      coeff(me,i1m,i2m,i3m) = extrapCoeff[me];
                      const int j1=i1m + me*is1, j2=i2m + me*is2, j3=i3m + me*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                      setEquationNumber(me, e,i1m,i2m,i3m,  c,j1,j2,j3 );             // macro to set equationNumber
                  }                
        }
        if( twilightZoneFlow )
        {
      // --- For twilightZone we save some coefficients that go into the CHAMP matrix ---
            cc(i1,i2,i3, 0) = cI;   // coeff of I
            cc(i1,i2,i3, 1) = cr;   // coeff of ur
            cc(i1,i2,i3, 2) = cs;   // coeff of us
            cc(i1,i2,i3, 3) = crr;  // coeff of urr
            cc(i1,i2,i3, 4) = crs;  // coeff of urs
            cc(i1,i2,i3, 5) = css;  // coeff of uss
        } 
    } // end ForStencil
  if( debug & 4 ) 
  {
      printF("WDH:order=2: (i1,i2,i3)=(%3d,%3d,%3d)\n",i1,i2,i3);     
      printF("    coeff =(%10.3e,%10.3e,%10.3e,\n"    
                    "            %10.3e,%10.3e,%10.3e,\n"    
                    "            %10.3e,%10.3e,%10.3e )\n",  
                      coeff4(-1,-1),coeff4(0,-1),coeff4(1,-1),  
                      coeff4(-1, 0),coeff4(0, 0),coeff4(1, 0),  
                      coeff4(-1, 1),coeff4(0, 1),coeff4(1, 1)); 
  }
} // end FOR_3D
                            }
                            else
                            {
// Coefficients for CHAMP curvilinear order 4, NORMAL DIRECTION=0.
// File written by champGenCode.maple
const Real KLR = theta;
const Real DLR = beta;
printF("WDH: KLR=%9.3e, DLR=%9.3e\n", KLR,DLR);
// Derivatives of the metrics: order 4
// Derivatives of the metrics: order 2
// Fourth-order one-sided approximations
// Here are fully one-sided: *check me* *wdh* April 15, 2022
// Macro : fill end ghost points on b1, b2, by periodicity or extrapolation
                  
// Macro : fill end ghost points on crr, crs, ... by periodicity or extrapolation
                  
const int axis1L=axis,  axis2L=(axis1L+1) % numberOfDimensions;
const int axis1R=axis2, axis2R=(axis1R+1) % numberOfDimensions;
const IntegerArray & gir1 = mg.indexRange();
const IntegerArray & gir2 = mg2.indexRange();
int ksv[3], &ks1=ksv[0], &ks2=ksv[1], &ks3=ksv[2];
int axist;
//  ----- We need b1,b2 at 4 extra tangential ghosts (order=4)---- 
//  since L1 is evaluated at 4 extra ghost 
int extra=2*numGhost;
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray n1L(Ib1,Ib2,Ib3), n2L(Ib1,Ib2,Ib3)  ; // save the normal, it needs to be extended too;
RealArray b1L(Ib1,Ib2,Ib3), b2L(Ib1,Ib2,Ib3);
RealArray b1R(Jb1,Jb2,Jb3), b2R(Jb1,Jb2,Jb3);
// Evaluate the normal derivative coefficients at points on the boundary. Use right-normal = -left-normal.
extra=numGhost;  // we can only directly evaluate 2 ghost using metrics 
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
FOR_3IJD(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    n1L(i1,i2,i3) = normal(i1,i2,i3,0);
    n2L(i1,i2,i3) = normal(i1,i2,i3,1);
    b1L(i1,i2,i3) = -( normal(i1,i2,i3,0)*rxL(i1,i2,i3,0,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,0,1)) ;
    b2L(i1,i2,i3) = -( normal(i1,i2,i3,0)*rxL(i1,i2,i3,1,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,1,1)) ;
    b1R(j1,j2,j3) = -( normal(i1,i2,i3,0)*rxR(j1,j2,j3,0,0) + normal(i1,i2,i3,1)*rxR(j1,j2,j3,0,1)) ;
    b2R(j1,j2,j3) = -( normal(i1,i2,i3,0)*rxR(j1,j2,j3,1,0) + normal(i1,i2,i3,1)*rxR(j1,j2,j3,1,1)) ;
} // end FOR_3IJD

// Fill 2nd ghost in tangential directions for b1, b2, and nL 
const int ghostToFillStart = numGhost+1;
const int ghostToFillEnd   = 2*numGhost;
    axist = (axis + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir1(0,0),gir1(0,0)); 
        I2=Range(gir1(0,1),gir1(0,1)); 
        I3=Range(gir1(0,2),gir1(0,2)); 
        Iv[axis]=Range(gir1(side,axis),gir1(side,axis));
        Iv[axist]=gir1(sidet,axist);
        for( int ghost=ghostToFillStart; ghost<=ghostToFillEnd; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir1(1,0)-gir1(0,0))*ks1;
                Index Ip2 = Ig2 + (gir1(1,1)-gir1(0,1))*ks2;
                Index Ip3 = Ig3 + (gir1(1,2)-gir1(0,2))*ks3;
                n1L(Ig1,Ig2,Ig3) = n1L(Ip1,Ip2,Ip3); 
                n2L(Ig1,Ig2,Ig3) = n2L(Ip1,Ip2,Ip3); 
            }
            else;
            {
                n1L(Ig1,Ig2,Ig3) = (5.*n1L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*n1L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*n1L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*n1L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+n1L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                n2L(Ig1,Ig2,Ig3) = (5.*n2L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*n2L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*n2L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*n2L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+n2L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t
    axist = (axis + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir1(0,0),gir1(0,0)); 
        I2=Range(gir1(0,1),gir1(0,1)); 
        I3=Range(gir1(0,2),gir1(0,2)); 
        Iv[axis]=Range(gir1(side,axis),gir1(side,axis));
        Iv[axist]=gir1(sidet,axist);
        for( int ghost=ghostToFillStart; ghost<=ghostToFillEnd; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir1(1,0)-gir1(0,0))*ks1;
                Index Ip2 = Ig2 + (gir1(1,1)-gir1(0,1))*ks2;
                Index Ip3 = Ig3 + (gir1(1,2)-gir1(0,2))*ks3;
                b1L(Ig1,Ig2,Ig3) = b1L(Ip1,Ip2,Ip3); 
                b2L(Ig1,Ig2,Ig3) = b2L(Ip1,Ip2,Ip3); 
            }
            else;
            {
                b1L(Ig1,Ig2,Ig3) = (5.*b1L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b1L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b1L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b1L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b1L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                b2L(Ig1,Ig2,Ig3) = (5.*b2L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b2L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b2L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b2L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b2L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t
    axist = (axis2 + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir2(0,0),gir2(0,0)); 
        I2=Range(gir2(0,1),gir2(0,1)); 
        I3=Range(gir2(0,2),gir2(0,2)); 
        Iv[axis2]=Range(gir2(side2,axis2),gir2(side2,axis2));
        Iv[axist]=gir2(sidet,axist);
        for( int ghost=ghostToFillStart; ghost<=ghostToFillEnd; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg2.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir2(1,0)-gir2(0,0))*ks1;
                Index Ip2 = Ig2 + (gir2(1,1)-gir2(0,1))*ks2;
                Index Ip3 = Ig3 + (gir2(1,2)-gir2(0,2))*ks3;
                b1R(Ig1,Ig2,Ig3) = b1R(Ip1,Ip2,Ip3); 
                b2R(Ig1,Ig2,Ig3) = b2R(Ip1,Ip2,Ip3); 
            }
            else;
            {
                b1R(Ig1,Ig2,Ig3) = (5.*b1R(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b1R(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b1R(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b1R(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b1R(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                b2R(Ig1,Ig2,Ig3) = (5.*b2R(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b2R(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b2R(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b2R(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b2R(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t

// --- Set the optimized Schwarz parameter: Sl = pl/dx (now that we know b1,b2)  ---- 
//  For now set a constant value on the face 
if( Sl<0 )
{
    j1=Jb1.getBase(); j2=Jb2.getBase(); j3=Jb3.getBase(); // take b1R from this point 
    const Real dnR = fabs( drR[axis2]/b1R(j1,j2,j3) );    // approx. normal grid spacing on right
    Sl = pl/dnR; 
    champParameters(4,side,axis,grid)=Sl; // save Sl
}
// *NOTE* WE NEED AN EXTRA GHOST IN TANGENTIAL DIRECTION FOR CURVILINEAR 
// Index bounds for points near the boundary: add 3 extra points in tangential directions and 1 in normal direction.
extra=numGhost+1;
getBoundaryIndex(mg.indexRange(),  side,axis, I1,I2,I3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3,extra);
int extraNormal=numGhost-1;
Iv[axis1L]=Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]=Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
RealArray crrL(I1,I2,I3), crsL(I1,I2,I3), cssL(I1,I2,I3), crL(I1,I2,I3), csL(I1,I2,I3);
RealArray crrR(J1,J2,J3), crsR(J1,J2,J3), cssR(J1,J2,J3), crR(J1,J2,J3), csR(J1,J2,J3);
RealArray rxxL(3,3), rxxR(3,3);
RealArray rxyL(3,3), rxyR(3,3);
// We can fill in crr, crs, and ss to 2 ghost since they do not involve terms using rxx, sxx, etc. 
extra=numGhost;
getBoundaryIndex(mg.indexRange(),  side,axis, I1,I2,I3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3,extra);
Iv[axis1L]=Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]=Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
FOR_3IJ(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3)
{
  // Evaluate coefficients of the Laplacian at points near the boundary
  // uxx = (rx*Dr + sx*Ds)*[ (rx)*ur + sx*us ]
  // uxx = (rx)^2 urr + 2*(rx*sx)*urs + (sx)^2 *uss + (rxx)*ur + (sxx)*us
  // uyy = (ry)^2 urr + 2*(ry*sy)*urs + (sy)^2 *uss + (ryy)*ur + (syy)*us
    crrL(i1,i2,i3) = SQR(rxL(i1,i2,i3,0,0)) + SQR(rxL(i1,i2,i3,0,1));                                // (rx)^2 + (ry)^2 
    crsL(i1,i2,i3) = 2.*(rxL(i1,i2,i3,0,0)*rxL(i1,i2,i3,1,0) + rxL(i1,i2,i3,0,1)*rxL(i1,i2,i3,1,1)); // 2( rx*sx + ry*sy ) 
    cssL(i1,i2,i3) = SQR(rxL(i1,i2,i3,1,0)) + SQR(rxL(i1,i2,i3,1,1));                                // (sx)^2 + (sy)^2
  // These are done below
  //crL(i1,i2,i3)  = rxxL(0,0) + rxyL(0,1);  // rxx + ryy
  //csL(i1,i2,i3)  = rxxL(1,0) + rxyL(1,1);  // sxx + syy
  // Evaluate coefficients of the Laplacian at points near the boundary
    crrR(j1,j2,j3) = SQR(rxR(j1,j2,j3,0,0)) + SQR(rxR(j1,j2,j3,0,1));
    crsR(j1,j2,j3) = 2.*( rxR(j1,j2,j3,0,0)*rxR(j1,j2,j3,1,0) + rxR(j1,j2,j3,0,1)*rxR(j1,j2,j3,1,1) );
    cssR(j1,j2,j3) = SQR(rxR(j1,j2,j3,1,0)) + SQR(rxR(j1,j2,j3,1,1));
  // These are done below
  //crR(j1,j2,j3)  = rxxR(0,0) + rxyR(0,1);  // rxx + ryy
  //csR(j1,j2,j3)  = rxxR(1,0) + rxyR(1,1);  // sxx + syy
  // printF("WDH: crrL=%9.3e, crsL=%9.3e, cssL=%9.3e\n",crrL(i1,i2,i3),crsL(i1,i2,i3),cssL(i1,i2,i3) );
  // printF("WDH: crrR=%9.3e, crsR=%9.3e, cssR=%9.3e\n",crrR(j1,j2,j3),crsR(j1,j2,j3),cssR(j1,j2,j3) );
} // end FOR_3IJ

// We can only compute cr and cs up to the boundary since we need two ghost points to evaluate rxx, rxy, sxx, sxy
getBoundaryIndex(mg.indexRange(),side,axis,I1,I2,I3);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3);
Iv[axis1L]= Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]= Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
Real rxr,rxs;
FOR_3IJ(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3)
{
    for( int m2=0; m2<numberOfDimensions; m2++ )
    {
        for( int m1=0; m1<numberOfDimensions; m1++ )
        {
        // Evaluate derivatives of metrics to order 4 .. we are computing more than we need here *fix me*
                if( i1-2 >= mg.dimension(0,0) && i1+2 <= mg.dimension(1,0) )
                    rxr = (((1./12.)*rxL(i1-2,i2,i3,m1,m2)-(2./3.)*rxL(i1-1,i2,i3,m1,m2)+(2./3.)*rxL(i1+1,i2,i3,m1,m2)-(1./12.)*rxL(i1+2,i2,i3,m1,m2))/(drL[0]));              // centered 
                else if( i1-2 < mg.dimension(0,0) )
                    rxr = (((-25./12.)*rxL(i1,i2,i3,m1,m2)+(4.)*rxL(i1+1,i2,i3,m1,m2)-(3.)*rxL(i1+2,i2,i3,m1,m2)+(4./3.)*rxL(i1+3,i2,i3,m1,m2)-(1./4.)*rxL(i1+4,i2,i3,m1,m2))/(drL[0]));   // one-sided to right
                else 
                    rxr = (((+25./12.)*rxL(i1,i2,i3,m1,m2)-(4.)*rxL(i1-1,i2,i3,m1,m2)+(3.)*rxL(i1-2,i2,i3,m1,m2)-(4./3.)*rxL(i1-3,i2,i3,m1,m2)+(1./4.)*rxL(i1-4,i2,i3,m1,m2))/(drL[0]));  // one-sided to left
                if( i2-2 >= mg.dimension(0,1) && i2+2 <= mg.dimension(1,1) )
                    rxs = (((1./12.)*rxL(i1,i2-2,i3,m1,m2)-(2./3.)*rxL(i1,i2-1,i3,m1,m2)+(2./3.)*rxL(i1,i2+1,i3,m1,m2)-(1./12.)*rxL(i1,i2+2,i3,m1,m2))/(drL[1]));              // centered 
                else if( i2-2 < mg.dimension(0,1) )
                    rxs = (((-25./12.)*rxL(i1,i2,i3,m1,m2)+(4.)*rxL(i1,i2+1,i3,m1,m2)-(3.)*rxL(i1,i2+2,i3,m1,m2)+(4./3.)*rxL(i1,i2+3,i3,m1,m2)-(1./4.)*rxL(i1,i2+4,i3,m1,m2))/(drL[1]));   // one-sided to right
                else 
                    rxs = (((+25./12.)*rxL(i1,i2,i3,m1,m2)-(4.)*rxL(i1,i2-1,i3,m1,m2)+(3.)*rxL(i1,i2-2,i3,m1,m2)-(4./3.)*rxL(i1,i2-3,i3,m1,m2)+(1./4.)*rxL(i1,i2-4,i3,m1,m2))/(drL[1]));  // one-sided to left
                rxxL(m1,m2) = rxL(i1,i2,i3,0,0)*rxr + rxL(i1,i2,i3,1,0)*rxs;
                rxyL(m1,m2) = rxL(i1,i2,i3,0,1)*rxr + rxL(i1,i2,i3,1,1)*rxs;
                                                                                                            
                if( j1-2 >= mg2.dimension(0,0) && j1+2 <= mg2.dimension(1,0) )
                    rxr = (((1./12.)*rxR(j1-2,j2,j3,m1,m2)-(2./3.)*rxR(j1-1,j2,j3,m1,m2)+(2./3.)*rxR(j1+1,j2,j3,m1,m2)-(1./12.)*rxR(j1+2,j2,j3,m1,m2))/(drR[0]));  // centered 
                else if( j1-2 < mg2.dimension(0,0) )
                    rxr = (((-25./12.)*rxR(j1,j2,j3,m1,m2)+(4.)*rxR(j1+1,j2,j3,m1,m2)-(3.)*rxR(j1+2,j2,j3,m1,m2)+(4./3.)*rxR(j1+3,j2,j3,m1,m2)-(1./4.)*rxR(j1+4,j2,j3,m1,m2))/(drR[0]));  // one-sided to right
                else 
                    rxr = (((+25./12.)*rxR(j1,j2,j3,m1,m2)-(4.)*rxR(j1-1,j2,j3,m1,m2)+(3.)*rxR(j1-2,j2,j3,m1,m2)-(4./3.)*rxR(j1-3,j2,j3,m1,m2)+(1./4.)*rxR(j1-4,j2,j3,m1,m2))/(drR[0]));  // one-sided to left
                if( j2-2 >= mg2.dimension(0,1) && j2+2 <= mg2.dimension(1,1) )
                    rxs = (((1./12.)*rxR(j1,j2-2,j3,m1,m2)-(2./3.)*rxR(j1,j2-1,j3,m1,m2)+(2./3.)*rxR(j1,j2+1,j3,m1,m2)-(1./12.)*rxR(j1,j2+2,j3,m1,m2))/(drR[1]));  // centered 
                else if( j2-2 < mg2.dimension(0,1) )
                    rxs = (((-25./12.)*rxR(j1,j2,j3,m1,m2)+(4.)*rxR(j1,j2+1,j3,m1,m2)-(3.)*rxR(j1,j2+2,j3,m1,m2)+(4./3.)*rxR(j1,j2+3,j3,m1,m2)-(1./4.)*rxR(j1,j2+4,j3,m1,m2))/(drR[1]));  // one-sided to right
                else 
                    rxs = (((+25./12.)*rxR(j1,j2,j3,m1,m2)-(4.)*rxR(j1,j2-1,j3,m1,m2)+(3.)*rxR(j1,j2-2,j3,m1,m2)-(4./3.)*rxR(j1,j2-3,j3,m1,m2)+(1./4.)*rxR(j1,j2-4,j3,m1,m2))/(drR[1]));  // one-sided to left
                rxxR(m1,m2) = rxR(j1,j2,j3,0,0)*rxr + rxR(j1,j2,j3,1,0)*rxs;
                rxyR(m1,m2) = rxR(j1,j2,j3,0,1)*rxr + rxR(j1,j2,j3,1,1)*rxs;
                                                                                                                                        
        }
    }
  // Evaluate coefficients of the Laplacian at points near the boundary
  // uxx = (rx*Dr + sx*Ds)*[ (rx)*ur + sx*us ]
  // uxx = (rx)^2 urr + 2*(rx*sx)*urs + (sx)^2 *uss + (rxx)*ur + (sxx)*us
  // uyy = (ry)^2 urr + 2*(ry*sy)*urs + (sy)^2 *uss + (ryy)*ur + (syy)*us
    crL(i1,i2,i3)  = rxxL(0,0) + rxyL(0,1);  // rxx + ryy
    csL(i1,i2,i3)  = rxxL(1,0) + rxyL(1,1);  // sxx + syy
  // Evaluate coefficients of the Laplacian at points near the boundary
    crR(j1,j2,j3)  = rxxR(0,0) + rxyR(0,1);  // rxx + ryy
    csR(j1,j2,j3)  = rxxR(1,0) + rxyR(1,1);  // sxx + syy
    if( advectionIsOn )  // include advection terms scaled by 1/D : (u/D).grad 
    {
            if( variableAdvection )
            {
                u1DL = advectVarL(i1,i2,i3,0)/DL;
                u2DL = advectVarL(i1,i2,i3,1)/DL;
                u1DR = advectVarR(j1,j2,j3,0)/DR;
                u2DR = advectVarR(j1,j2,j3,1)/DR;
            }
        crL(i1,i2,i3) = crL(i1,i2,i3) - ( u1DL*rxL(i1,i2,i3,0,0) + u2DL*rxL(i1,i2,i3,0,1) ); 
        csL(i1,i2,i3) = csL(i1,i2,i3) - ( u1DL*rxL(i1,i2,i3,1,0) + u2DL*rxL(i1,i2,i3,1,1) ); 
        crR(j1,j2,j3) = crR(j1,j2,j3) - ( u1DR*rxR(j1,j2,j3,0,0) + u2DR*rxR(j1,j2,j3,0,1) ); 
        csR(j1,j2,j3) = csR(j1,j2,j3) - ( u1DR*rxR(j1,j2,j3,1,0) + u2DR*rxR(j1,j2,j3,1,1) ); 
    }
  // printF("WDH: (i1,i2)=(%4d,%4d) crL=%9.3e csL=%9.3e\n",i1,i2,crL(i1,i2,i3),csL(i1,i2,i3) );
  // printF("WDH: (j1,j2)=(%4d,%4d) crR=%9.3e csR=%9.3e\n",j1,j2,crR(j1,j2,j3),csR(j1,j2,j3) );
} // end FOR_3IJ

// --- Now fill in one extra ghost for crr,crs,css and two extra ghost for cr,cs --- 
// --- Fill ghost points on ends through extrapolation or periodicity --- 
// Fill ghost in tangential directions for crr, crs, etc.
    axist = (axis + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir1(0,0),gir1(0,0)); 
        I2=Range(gir1(0,1),gir1(0,1)); 
        I3=Range(gir1(0,2),gir1(0,2)); 
        Iv[axis]=Range(gir1(side,axis)-1,gir1(side,axis)+1);
        Iv[axist]=gir1(sidet,axist);
        for( int ghost=1; ghost<=numGhost+1; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir1(1,0)-gir1(0,0))*ks1;
                Index Ip2 = Ig2 + (gir1(1,1)-gir1(0,1))*ks2;
                Index Ip3 = Ig3 + (gir1(1,2)-gir1(0,2))*ks3;
                crrL(Ig1,Ig2,Ig3) = crrL(Ip1,Ip2,Ip3); 
                crsL(Ig1,Ig2,Ig3) = crsL(Ip1,Ip2,Ip3); 
                cssL(Ig1,Ig2,Ig3) = cssL(Ip1,Ip2,Ip3); 
                crL(Ig1,Ig2,Ig3)  =  crL(Ip1,Ip2,Ip3); 
                csL(Ig1,Ig2,Ig3)  =  csL(Ip1,Ip2,Ip3); 
            }
            else;
            {
                if( ghost==3 )
                { 
                    crrL(Ig1,Ig2,Ig3) = (5.*crrL(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*crrL(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*crrL(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*crrL(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+crrL(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                    crsL(Ig1,Ig2,Ig3) = (5.*crsL(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*crsL(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*crsL(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*crsL(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+crsL(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                    cssL(Ig1,Ig2,Ig3) = (5.*cssL(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*cssL(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*cssL(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*cssL(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+cssL(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                } 
                crL(Ig1,Ig2,Ig3) = (5.*crL(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*crL(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*crL(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*crL(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+crL(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                csL(Ig1,Ig2,Ig3) = (5.*csL(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*csL(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*csL(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*csL(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+csL(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t
    axist = (axis2 + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir2(0,0),gir2(0,0)); 
        I2=Range(gir2(0,1),gir2(0,1)); 
        I3=Range(gir2(0,2),gir2(0,2)); 
        Iv[axis2]=Range(gir2(side2,axis2)-1,gir2(side2,axis2)+1);
        Iv[axist]=gir2(sidet,axist);
        for( int ghost=1; ghost<=numGhost+1; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg2.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir2(1,0)-gir2(0,0))*ks1;
                Index Ip2 = Ig2 + (gir2(1,1)-gir2(0,1))*ks2;
                Index Ip3 = Ig3 + (gir2(1,2)-gir2(0,2))*ks3;
                crrR(Ig1,Ig2,Ig3) = crrR(Ip1,Ip2,Ip3); 
                crsR(Ig1,Ig2,Ig3) = crsR(Ip1,Ip2,Ip3); 
                cssR(Ig1,Ig2,Ig3) = cssR(Ip1,Ip2,Ip3); 
                crR(Ig1,Ig2,Ig3)  =  crR(Ip1,Ip2,Ip3); 
                csR(Ig1,Ig2,Ig3)  =  csR(Ip1,Ip2,Ip3); 
            }
            else;
            {
                if( ghost==3 )
                { 
                    crrR(Ig1,Ig2,Ig3) = (5.*crrR(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*crrR(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*crrR(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*crrR(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+crrR(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                    crsR(Ig1,Ig2,Ig3) = (5.*crsR(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*crsR(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*crsR(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*crsR(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+crsR(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                    cssR(Ig1,Ig2,Ig3) = (5.*cssR(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*cssR(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*cssR(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*cssR(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+cssR(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                } 
                crR(Ig1,Ig2,Ig3) = (5.*crR(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*crR(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*crR(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*crR(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+crR(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                csR(Ig1,Ig2,Ig3) = (5.*csR(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*csR(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*csR(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*csR(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+csR(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t
                  
// ---- Evaluate derivatives of coefficients : (crr).rr, (csr).rs etc. 
//  L3 : needs (crr).r (crr).s etc at one tangential ghost. 
extra=numGhost-1;
getBoundaryIndex(mg.indexRange(),  side,axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray b1Lrr(Ib1,Ib2,Ib3), b1Lrs(Ib1,Ib2,Ib3), b1Lss(Ib1,Ib2,Ib3), b1Lr(Ib1,Ib2,Ib3), b1Ls(Ib1,Ib2,Ib3);
RealArray b2Lrr(Ib1,Ib2,Ib3), b2Lrs(Ib1,Ib2,Ib3), b2Lss(Ib1,Ib2,Ib3), b2Lr(Ib1,Ib2,Ib3), b2Ls(Ib1,Ib2,Ib3);
RealArray crrLrr(Ib1,Ib2,Ib3), crrLrs(Ib1,Ib2,Ib3), crrLss(Ib1,Ib2,Ib3), crrLr(Ib1,Ib2,Ib3), crrLs(Ib1,Ib2,Ib3);
RealArray crsLrr(Ib1,Ib2,Ib3), crsLrs(Ib1,Ib2,Ib3), crsLss(Ib1,Ib2,Ib3), crsLr(Ib1,Ib2,Ib3), crsLs(Ib1,Ib2,Ib3);
RealArray cssLrr(Ib1,Ib2,Ib3), cssLrs(Ib1,Ib2,Ib3), cssLss(Ib1,Ib2,Ib3), cssLr(Ib1,Ib2,Ib3), cssLs(Ib1,Ib2,Ib3);
RealArray crLrr(Ib1,Ib2,Ib3), crLrs(Ib1,Ib2,Ib3), crLss(Ib1,Ib2,Ib3), crLr(Ib1,Ib2,Ib3), crLs(Ib1,Ib2,Ib3);
RealArray csLrr(Ib1,Ib2,Ib3), csLrs(Ib1,Ib2,Ib3), csLss(Ib1,Ib2,Ib3), csLr(Ib1,Ib2,Ib3), csLs(Ib1,Ib2,Ib3);
// Evaluate derivatives of coefficients at points on the boundary (left)
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  // tangential directives of b1,b2 only 
    b1Lss(i1,i2,i3) = ((b1L(i1,i2+1,i3)-2.*b1L(i1,i2,i3)+b1L(i1,i2-1,i3))/(SQR(drL[1])));
      b1Ls(i1,i2,i3) =  ((b1L(i1,i2+1,i3)-b1L(i1,i2-1,i3))/(2.*drL[1]));
    b2Lss(i1,i2,i3) = ((b2L(i1,i2+1,i3)-2.*b2L(i1,i2,i3)+b2L(i1,i2-1,i3))/(SQR(drL[1])));
      b2Ls(i1,i2,i3) =  ((b2L(i1,i2+1,i3)-b2L(i1,i2-1,i3))/(2.*drL[1]));
  // normal and tangential directives of c 
    crrLrr(i1,i2,i3) = ((crrL(i1+1,i2,i3)-2.*crrL(i1,i2,i3)+crrL(i1-1,i2,i3))/(SQR(drL[0])));
    crrLss(i1,i2,i3) = ((crrL(i1,i2+1,i3)-2.*crrL(i1,i2,i3)+crrL(i1,i2-1,i3))/(SQR(drL[1])));
    crrLrs(i1,i2,i3) = ((((crrL(i1+1,i2+1,i3)-crrL(i1+1,i2-1,i3))/(2.*drL[1]))-((crrL(i1-1,i2+1,i3)-crrL(i1-1,i2-1,i3))/(2.*drL[1])))/(2.*drL[0]));
      crrLr(i1,i2,i3) =  ((crrL(i1+1,i2,i3)-crrL(i1-1,i2,i3))/(2.*drL[0]));
      crrLs(i1,i2,i3) =  ((crrL(i1,i2+1,i3)-crrL(i1,i2-1,i3))/(2.*drL[1]));
    crsLrr(i1,i2,i3) = ((crsL(i1+1,i2,i3)-2.*crsL(i1,i2,i3)+crsL(i1-1,i2,i3))/(SQR(drL[0])));
    crsLss(i1,i2,i3) = ((crsL(i1,i2+1,i3)-2.*crsL(i1,i2,i3)+crsL(i1,i2-1,i3))/(SQR(drL[1])));
    crsLrs(i1,i2,i3) = ((((crsL(i1+1,i2+1,i3)-crsL(i1+1,i2-1,i3))/(2.*drL[1]))-((crsL(i1-1,i2+1,i3)-crsL(i1-1,i2-1,i3))/(2.*drL[1])))/(2.*drL[0]));
      crsLr(i1,i2,i3) =  ((crsL(i1+1,i2,i3)-crsL(i1-1,i2,i3))/(2.*drL[0]));
      crsLs(i1,i2,i3) =  ((crsL(i1,i2+1,i3)-crsL(i1,i2-1,i3))/(2.*drL[1]));
    cssLrr(i1,i2,i3) = ((cssL(i1+1,i2,i3)-2.*cssL(i1,i2,i3)+cssL(i1-1,i2,i3))/(SQR(drL[0])));
    cssLss(i1,i2,i3) = ((cssL(i1,i2+1,i3)-2.*cssL(i1,i2,i3)+cssL(i1,i2-1,i3))/(SQR(drL[1])));
    cssLrs(i1,i2,i3) = ((((cssL(i1+1,i2+1,i3)-cssL(i1+1,i2-1,i3))/(2.*drL[1]))-((cssL(i1-1,i2+1,i3)-cssL(i1-1,i2-1,i3))/(2.*drL[1])))/(2.*drL[0]));
      cssLr(i1,i2,i3) =  ((cssL(i1+1,i2,i3)-cssL(i1-1,i2,i3))/(2.*drL[0]));
      cssLs(i1,i2,i3) =  ((cssL(i1,i2+1,i3)-cssL(i1,i2-1,i3))/(2.*drL[1]));
    crLrr(i1,i2,i3) = ((crL(i1+1,i2,i3)-2.*crL(i1,i2,i3)+crL(i1-1,i2,i3))/(SQR(drL[0])));
    crLss(i1,i2,i3) = ((crL(i1,i2+1,i3)-2.*crL(i1,i2,i3)+crL(i1,i2-1,i3))/(SQR(drL[1])));
    crLrs(i1,i2,i3) = ((((crL(i1+1,i2+1,i3)-crL(i1+1,i2-1,i3))/(2.*drL[1]))-((crL(i1-1,i2+1,i3)-crL(i1-1,i2-1,i3))/(2.*drL[1])))/(2.*drL[0]));
      crLr(i1,i2,i3) =  ((crL(i1+1,i2,i3)-crL(i1-1,i2,i3))/(2.*drL[0]));
      crLs(i1,i2,i3) =  ((crL(i1,i2+1,i3)-crL(i1,i2-1,i3))/(2.*drL[1]));
    csLrr(i1,i2,i3) = ((csL(i1+1,i2,i3)-2.*csL(i1,i2,i3)+csL(i1-1,i2,i3))/(SQR(drL[0])));
    csLss(i1,i2,i3) = ((csL(i1,i2+1,i3)-2.*csL(i1,i2,i3)+csL(i1,i2-1,i3))/(SQR(drL[1])));
    csLrs(i1,i2,i3) = ((((csL(i1+1,i2+1,i3)-csL(i1+1,i2-1,i3))/(2.*drL[1]))-((csL(i1-1,i2+1,i3)-csL(i1-1,i2-1,i3))/(2.*drL[1])))/(2.*drL[0]));
      csLr(i1,i2,i3) =  ((csL(i1+1,i2,i3)-csL(i1-1,i2,i3))/(2.*drL[0]));
      csLs(i1,i2,i3) =  ((csL(i1,i2+1,i3)-csL(i1,i2-1,i3))/(2.*drL[1]));
} // end FOR_3D
RealArray b1Rrr(Jb1,Jb2,Jb3), b1Rrs(Jb1,Jb2,Jb3), b1Rss(Jb1,Jb2,Jb3), b1Rr(Jb1,Jb2,Jb3), b1Rs(Jb1,Jb2,Jb3);
RealArray b2Rrr(Jb1,Jb2,Jb3), b2Rrs(Jb1,Jb2,Jb3), b2Rss(Jb1,Jb2,Jb3), b2Rr(Jb1,Jb2,Jb3), b2Rs(Jb1,Jb2,Jb3);
RealArray crrRrr(Jb1,Jb2,Jb3), crrRrs(Jb1,Jb2,Jb3), crrRss(Jb1,Jb2,Jb3), crrRr(Jb1,Jb2,Jb3), crrRs(Jb1,Jb2,Jb3);
RealArray crsRrr(Jb1,Jb2,Jb3), crsRrs(Jb1,Jb2,Jb3), crsRss(Jb1,Jb2,Jb3), crsRr(Jb1,Jb2,Jb3), crsRs(Jb1,Jb2,Jb3);
RealArray cssRrr(Jb1,Jb2,Jb3), cssRrs(Jb1,Jb2,Jb3), cssRss(Jb1,Jb2,Jb3), cssRr(Jb1,Jb2,Jb3), cssRs(Jb1,Jb2,Jb3);
RealArray crRrr(Jb1,Jb2,Jb3), crRrs(Jb1,Jb2,Jb3), crRss(Jb1,Jb2,Jb3), crRr(Jb1,Jb2,Jb3), crRs(Jb1,Jb2,Jb3);
RealArray csRrr(Jb1,Jb2,Jb3), csRrs(Jb1,Jb2,Jb3), csRss(Jb1,Jb2,Jb3), csRr(Jb1,Jb2,Jb3), csRs(Jb1,Jb2,Jb3);
// Evaluate derivatives of coefficients at points on the boundary (right)
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
 // tangential directives of b1,b2 only 
    b1Rss(j1,j2,j3) = ((b1R(j1,j2+1,j3)-2.*b1R(j1,j2,j3)+b1R(j1,j2-1,j3))/(SQR(drR[1])));
      b1Rs(j1,j2,j3) =  ((b1R(j1,j2+1,j3)-b1R(j1,j2-1,j3))/(2.*drR[1]));
    b2Rss(j1,j2,j3) = ((b2R(j1,j2+1,j3)-2.*b2R(j1,j2,j3)+b2R(j1,j2-1,j3))/(SQR(drR[1])));
      b2Rs(j1,j2,j3) =  ((b2R(j1,j2+1,j3)-b2R(j1,j2-1,j3))/(2.*drR[1]));
  // normal and tangential directives of c 
    crrRrr(j1,j2,j3) = ((crrR(j1+1,j2,j3)-2.*crrR(j1,j2,j3)+crrR(j1-1,j2,j3))/(SQR(drR[0])));
    crrRss(j1,j2,j3) = ((crrR(j1,j2+1,j3)-2.*crrR(j1,j2,j3)+crrR(j1,j2-1,j3))/(SQR(drR[1])));
    crrRrs(j1,j2,j3) = ((((crrR(j1+1,j2+1,j3)-crrR(j1+1,j2-1,j3))/(2.*drR[1]))-((crrR(j1-1,j2+1,j3)-crrR(j1-1,j2-1,j3))/(2.*drR[1])))/(2.*drR[0]));
      crrRr(j1,j2,j3) =  ((crrR(j1+1,j2,j3)-crrR(j1-1,j2,j3))/(2.*drR[0]));
      crrRs(j1,j2,j3) =  ((crrR(j1,j2+1,j3)-crrR(j1,j2-1,j3))/(2.*drR[1]));
    crsRrr(j1,j2,j3) = ((crsR(j1+1,j2,j3)-2.*crsR(j1,j2,j3)+crsR(j1-1,j2,j3))/(SQR(drR[0])));
    crsRss(j1,j2,j3) = ((crsR(j1,j2+1,j3)-2.*crsR(j1,j2,j3)+crsR(j1,j2-1,j3))/(SQR(drR[1])));
    crsRrs(j1,j2,j3) = ((((crsR(j1+1,j2+1,j3)-crsR(j1+1,j2-1,j3))/(2.*drR[1]))-((crsR(j1-1,j2+1,j3)-crsR(j1-1,j2-1,j3))/(2.*drR[1])))/(2.*drR[0]));
      crsRr(j1,j2,j3) =  ((crsR(j1+1,j2,j3)-crsR(j1-1,j2,j3))/(2.*drR[0]));
      crsRs(j1,j2,j3) =  ((crsR(j1,j2+1,j3)-crsR(j1,j2-1,j3))/(2.*drR[1]));
    cssRrr(j1,j2,j3) = ((cssR(j1+1,j2,j3)-2.*cssR(j1,j2,j3)+cssR(j1-1,j2,j3))/(SQR(drR[0])));
    cssRss(j1,j2,j3) = ((cssR(j1,j2+1,j3)-2.*cssR(j1,j2,j3)+cssR(j1,j2-1,j3))/(SQR(drR[1])));
    cssRrs(j1,j2,j3) = ((((cssR(j1+1,j2+1,j3)-cssR(j1+1,j2-1,j3))/(2.*drR[1]))-((cssR(j1-1,j2+1,j3)-cssR(j1-1,j2-1,j3))/(2.*drR[1])))/(2.*drR[0]));
      cssRr(j1,j2,j3) =  ((cssR(j1+1,j2,j3)-cssR(j1-1,j2,j3))/(2.*drR[0]));
      cssRs(j1,j2,j3) =  ((cssR(j1,j2+1,j3)-cssR(j1,j2-1,j3))/(2.*drR[1]));
    crRrr(j1,j2,j3) = ((crR(j1+1,j2,j3)-2.*crR(j1,j2,j3)+crR(j1-1,j2,j3))/(SQR(drR[0])));
    crRss(j1,j2,j3) = ((crR(j1,j2+1,j3)-2.*crR(j1,j2,j3)+crR(j1,j2-1,j3))/(SQR(drR[1])));
    crRrs(j1,j2,j3) = ((((crR(j1+1,j2+1,j3)-crR(j1+1,j2-1,j3))/(2.*drR[1]))-((crR(j1-1,j2+1,j3)-crR(j1-1,j2-1,j3))/(2.*drR[1])))/(2.*drR[0]));
      crRr(j1,j2,j3) =  ((crR(j1+1,j2,j3)-crR(j1-1,j2,j3))/(2.*drR[0]));
      crRs(j1,j2,j3) =  ((crR(j1,j2+1,j3)-crR(j1,j2-1,j3))/(2.*drR[1]));
    csRrr(j1,j2,j3) = ((csR(j1+1,j2,j3)-2.*csR(j1,j2,j3)+csR(j1-1,j2,j3))/(SQR(drR[0])));
    csRss(j1,j2,j3) = ((csR(j1,j2+1,j3)-2.*csR(j1,j2,j3)+csR(j1,j2-1,j3))/(SQR(drR[1])));
    csRrs(j1,j2,j3) = ((((csR(j1+1,j2+1,j3)-csR(j1+1,j2-1,j3))/(2.*drR[1]))-((csR(j1-1,j2+1,j3)-csR(j1-1,j2-1,j3))/(2.*drR[1])))/(2.*drR[0]));
      csRr(j1,j2,j3) =  ((csR(j1+1,j2,j3)-csR(j1-1,j2,j3))/(2.*drR[0]));
      csRs(j1,j2,j3) =  ((csR(j1,j2+1,j3)-csR(j1,j2-1,j3))/(2.*drR[1]));
} // end FOR_3IJ

// ----- Evaluate coefficients of L1 on the boundary -----

// 2*numGhost = 4 extra ghost needed for fourth order 
extra=2*numGhost;
getBoundaryIndex(mg.indexRange(), side ,axis ,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L1ICoeff(Ib1,Ib2,Ib3), L1ICoeffs(Ib1,Ib2,Ib3), L1ICoeffss(Ib1,Ib2,Ib3), L1ICoeffsss(Ib1,Ib2,Ib3);
RealArray L1rCoeff(Ib1,Ib2,Ib3), L1rCoeffs(Ib1,Ib2,Ib3), L1rCoeffss(Ib1,Ib2,Ib3), L1rCoeffsss(Ib1,Ib2,Ib3);
RealArray L1sCoeff(Ib1,Ib2,Ib3), L1sCoeffs(Ib1,Ib2,Ib3), L1sCoeffss(Ib1,Ib2,Ib3), L1sCoeffsss(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    const Real n1R = -n1L(i1,i2,i3), n2R = -n2L(i1,i2,i3);    // for advection, use right normal 
        if( variableAdvection )
        {
            u1DL = advectVarL(i1,i2,i3,0)/DL;
            u2DL = advectVarL(i1,i2,i3,1)/DL;
            u1DR = advectVarR(j1,j2,j3,0)/DR;
            u2DR = advectVarR(j1,j2,j3,1)/DR;
        }
    L1ICoeff(i1,i2,i3) = ((-1.*KLR*u2DL+1.*u2DR)*n2R+(-1.*KLR*u1DL+1.*u1DR)*n1R)/b1R(j1,j2,j3);
    L1rCoeff(i1,i2,i3) = 1.*KLR*b1L(i1,i2,i3)/b1R(j1,j2,j3);
    L1sCoeff(i1,i2,i3) = (1.*KLR*b2L(i1,i2,i3)-1.*b2R(j1,j2,j3))/b1R(j1,j2,j3);
} // end FOR_3D

// ----- Evaluate coefficients of L2 on the boundary -----
extra=numGhost;
getBoundaryIndex(mg.indexRange(), side,  axis,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L2rrCoeff(Ib1,Ib2,Ib3), L2rrCoeffs(Ib1,Ib2,Ib3), L2rrCoeffss(Ib1,Ib2,Ib3);
RealArray L2rsCoeff(Ib1,Ib2,Ib3), L2rsCoeffs(Ib1,Ib2,Ib3), L2rsCoeffss(Ib1,Ib2,Ib3);
RealArray L2ssCoeff(Ib1,Ib2,Ib3), L2ssCoeffs(Ib1,Ib2,Ib3), L2ssCoeffss(Ib1,Ib2,Ib3);
RealArray L2rCoeff(Ib1,Ib2,Ib3), L2rCoeffs(Ib1,Ib2,Ib3), L2rCoeffss(Ib1,Ib2,Ib3);
RealArray L2sCoeff(Ib1,Ib2,Ib3), L2sCoeffs(Ib1,Ib2,Ib3), L2sCoeffss(Ib1,Ib2,Ib3);
RealArray L2ICoeff(Ib1,Ib2,Ib3), L2ICoeffs(Ib1,Ib2,Ib3), L2ICoeffss(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  // Compute tangential derivatives of L1 - label tangential directions a and b (for 3D)
    L1rCoeffs(i1,i2,i3) = (((1./12.)*L1rCoeff(i1,i2-2,i3)-(2./3.)*L1rCoeff(i1,i2-1,i3)+(2./3.)*L1rCoeff(i1,i2+1,i3)-(1./12.)*L1rCoeff(i1,i2+2,i3))/(drL[1]));
    L1sCoeffs(i1,i2,i3) = (((1./12.)*L1sCoeff(i1,i2-2,i3)-(2./3.)*L1sCoeff(i1,i2-1,i3)+(2./3.)*L1sCoeff(i1,i2+1,i3)-(1./12.)*L1sCoeff(i1,i2+2,i3))/(drL[1]));
    L1ICoeffs(i1,i2,i3) = (((1./12.)*L1ICoeff(i1,i2-2,i3)-(2./3.)*L1ICoeff(i1,i2-1,i3)+(2./3.)*L1ICoeff(i1,i2+1,i3)-(1./12.)*L1ICoeff(i1,i2+2,i3))/(drL[1]));
  // Compute L2 coefficients
    L2rrCoeff(i1,i2,i3) = 1.*DLR*crrL(i1,i2,i3)/crrR(j1,j2,j3);
    L2rsCoeff(i1,i2,i3) = (1.*DLR*crsL(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1rCoeff(i1,i2,i3))/crrR(j1,j2,j3);
    L2ssCoeff(i1,i2,i3) = 1.*(DLR*cssL(i1,i2,i3)-crsR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-cssR(j1,j2,j3))/crrR(j1,j2,j3);
    L2rCoeff(i1,i2,i3) = 1.*(DLR*crL(i1,i2,i3)-crR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-crsR(j1,j2,j3)*L1rCoeffs(i1,i2,i3))/crrR(j1,j2,j3);
    L2sCoeff(i1,i2,i3) = ((-1.*L1ICoeff(i1,i2,i3)-1.*L1sCoeffs(i1,i2,i3))*crsR(j1,j2,j3)+1.*DLR*csL(i1,i2,i3)-1.*crR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*csR(j1,j2,j3))/crrR(j1,j2,j3);
    L2ICoeff(i1,i2,i3) = (-1.*crR(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1ICoeffs(i1,i2,i3))/crrR(j1,j2,j3);
} // end FOR_3D

// ----- Evaluate coefficients of L3 on the boundary -----
extra=numGhost-1;
getBoundaryIndex(mg.indexRange(),side,axis,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L3rrrCoeff(Ib1,Ib2,Ib3), L3rrrCoeffs(Ib1,Ib2,Ib3);
RealArray L3rrsCoeff(Ib1,Ib2,Ib3), L3rrsCoeffs(Ib1,Ib2,Ib3);
RealArray L3rssCoeff(Ib1,Ib2,Ib3), L3rssCoeffs(Ib1,Ib2,Ib3);
RealArray L3sssCoeff(Ib1,Ib2,Ib3), L3sssCoeffs(Ib1,Ib2,Ib3);
RealArray L3rrCoeff(Ib1,Ib2,Ib3), L3rrCoeffs(Ib1,Ib2,Ib3);
RealArray L3rsCoeff(Ib1,Ib2,Ib3), L3rsCoeffs(Ib1,Ib2,Ib3);
RealArray L3ssCoeff(Ib1,Ib2,Ib3), L3ssCoeffs(Ib1,Ib2,Ib3);
RealArray L3rCoeff(Ib1,Ib2,Ib3), L3rCoeffs(Ib1,Ib2,Ib3);
RealArray L3sCoeff(Ib1,Ib2,Ib3), L3sCoeffs(Ib1,Ib2,Ib3);
RealArray L3ICoeff(Ib1,Ib2,Ib3), L3ICoeffs(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    const Real n1R = -n1L(i1,i2,i3), n2R = -n2L(i1,i2,i3);    // for advection, use right normal 
        if( variableAdvection )
        {
            u1DL = advectVarL(i1,i2,i3,0)/DL;
            u2DL = advectVarL(i1,i2,i3,1)/DL;
            u1DR = advectVarR(j1,j2,j3,0)/DR;
            u2DR = advectVarR(j1,j2,j3,1)/DR;
        }
  // Compute tangential derivatives of L1 
    L1rCoeffss(i1,i2,i3) = ((L1rCoeff(i1,i2+1,i3)-2.*L1rCoeff(i1,i2,i3)+L1rCoeff(i1,i2-1,i3))/(SQR(drL[1])));
    L1sCoeffss(i1,i2,i3) = ((L1sCoeff(i1,i2+1,i3)-2.*L1sCoeff(i1,i2,i3)+L1sCoeff(i1,i2-1,i3))/(SQR(drL[1])));
    L1ICoeffss(i1,i2,i3) = ((L1ICoeff(i1,i2+1,i3)-2.*L1ICoeff(i1,i2,i3)+L1ICoeff(i1,i2-1,i3))/(SQR(drL[1])));
  // Compute tangential derivatives of L2
    L2rrCoeffs(i1,i2,i3) = ((L2rrCoeff(i1,i2+1,i3)-L2rrCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L2rsCoeffs(i1,i2,i3) = ((L2rsCoeff(i1,i2+1,i3)-L2rsCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L2ssCoeffs(i1,i2,i3) = ((L2ssCoeff(i1,i2+1,i3)-L2ssCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L2rCoeffs(i1,i2,i3)  = ((L2rCoeff(i1,i2+1,i3)-L2rCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L2sCoeffs(i1,i2,i3)  = ((L2sCoeff(i1,i2+1,i3)-L2sCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L2ICoeffs(i1,i2,i3)  = ((L2ICoeff(i1,i2+1,i3)-L2ICoeff(i1,i2-1,i3))/(2.*drL[1]));
  // Compute L3 coefficients
    L3rrrCoeff(i1,i2,i3) = 1.*DLR*KLR*b1L(i1,i2,i3)*crrL(i1,i2,i3)/b1R(j1,j2,j3)/crrR(j1,j2,j3);
    L3rrsCoeff(i1,i2,i3) = ((-1.*crsR(j1,j2,j3)*b1R(j1,j2,j3)-1.*crrR(j1,j2,j3)*b2R(j1,j2,j3))*L2rrCoeff(i1,i2,i3)+(1.*crsL(i1,i2,i3)*b1L(i1,i2,i3)+1.*b2L(i1,i2,i3)*crrL(i1,i2,i3))*KLR*DLR)/b1R(j1,j2,j3)/crrR(j1,j2,j3);
    L3rssCoeff(i1,i2,i3) = 1.*(DLR*KLR*b1L(i1,i2,i3)*cssL(i1,i2,i3)+DLR*KLR*b2L(i1,i2,i3)*crsL(i1,i2,i3)-b1R(j1,j2,j3)*crsR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)-b1R(j1,j2,j3)*cssR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-b2R(j1,j2,j3)*crrR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)-b2R(j1,j2,j3)*crsR(j1,j2,j3)*L1rCoeff(i1,i2,i3))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
    L3sssCoeff(i1,i2,i3) = ((-1.*L2ssCoeff(i1,i2,i3)*crsR(j1,j2,j3)-1.*L1sCoeff(i1,i2,i3)*cssR(j1,j2,j3))*b1R(j1,j2,j3)+(-1.*crrR(j1,j2,j3)*L2ssCoeff(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*cssR(j1,j2,j3))*b2R(j1,j2,j3)+1.*DLR*KLR*b2L(i1,i2,i3)*cssL(i1,i2,i3))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
    L3rrCoeff(i1,i2,i3) = (((1.*n1R*u1DR+1.*n2R*u2DR)*crrR(j1,j2,j3)+(-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*b1R(j1,j2,j3)-1.*crrRs(j1,j2,j3)*b2R(j1,j2,j3))*L2rrCoeff(i1,i2,i3)+((-1.*n1R*u1DL-1.*n2R*u2DL)*crrL(i1,i2,i3)+(1.*crL(i1,i2,i3)+1.*crrLr(i1,i2,i3))*b1L(i1,i2,i3)+1.*b2L(i1,i2,i3)*crrLs(i1,i2,i3))*KLR*DLR-1.*b1R(j1,j2,j3)*crsR(j1,j2,j3)*L2rrCoeffs(i1,i2,i3)-1.*b2R(j1,j2,j3)*crrR(j1,j2,j3)*L2rrCoeffs(i1,i2,i3))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
    L3rsCoeff(i1,i2,i3) = (((-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1rCoeff(i1,i2,i3)+(-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*L2rsCoeff(i1,i2,i3)+(-1.*L2rsCoeffs(i1,i2,i3)-1.*L2rCoeff(i1,i2,i3))*crsR(j1,j2,j3)-2.*cssR(j1,j2,j3)*L1rCoeffs(i1,i2,i3))*b1R(j1,j2,j3)+((-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*b2R(j1,j2,j3)+crsR(j1,j2,j3)*(n1R*u1DR+n2R*u2DR))*L1rCoeff(i1,i2,i3)+(-1.*crsL(i1,i2,i3)*n1R*u1DL-1.*crsL(i1,i2,i3)*n2R*u2DL+(crsLr(i1,i2,i3)+csL(i1,i2,i3))*b1L(i1,i2,i3)+b2L(i1,i2,i3)*(crsLs(i1,i2,i3)+crL(i1,i2,i3)))*KLR*DLR+(-1.*crrRs(j1,j2,j3)*L2rsCoeff(i1,i2,i3)+(-1.*L2rsCoeffs(i1,i2,i3)-1.*L2rCoeff(i1,i2,i3))*crrR(j1,j2,j3)-2.*crsR(j1,j2,j3)*L1rCoeffs(i1,i2,i3))*b2R(j1,j2,j3)+crrR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)*(n1R*u1DR+n2R*u2DR))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
    L3ssCoeff(i1,i2,i3) = (((-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1sCoeff(i1,i2,i3)+(-1.*L2ssCoeffs(i1,i2,i3)-1.*L2sCoeff(i1,i2,i3))*crsR(j1,j2,j3)+(-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*L2ssCoeff(i1,i2,i3)+(-2.*L1sCoeffs(i1,i2,i3)-1.*L1ICoeff(i1,i2,i3))*cssR(j1,j2,j3)-1.*cssRr(j1,j2,j3))*b1R(j1,j2,j3)+((-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1sCoeff(i1,i2,i3)+(-2.*L1sCoeffs(i1,i2,i3)-1.*L1ICoeff(i1,i2,i3))*crsR(j1,j2,j3)-1.*crrRs(j1,j2,j3)*L2ssCoeff(i1,i2,i3)+(-1.*L2ssCoeffs(i1,i2,i3)-1.*L2sCoeff(i1,i2,i3))*crrR(j1,j2,j3)-1.*cssRs(j1,j2,j3)-1.*csR(j1,j2,j3))*b2R(j1,j2,j3)+crsR(j1,j2,j3)*(n1R*u1DR+n2R*u2DR)*L1sCoeff(i1,i2,i3)+crrR(j1,j2,j3)*(n1R*u1DR+n2R*u2DR)*L2ssCoeff(i1,i2,i3)+(-1.*cssL(i1,i2,i3)*n1R*u1DL-1.*cssL(i1,i2,i3)*n2R*u2DL+(cssLs(i1,i2,i3)+csL(i1,i2,i3))*b2L(i1,i2,i3)+cssLr(i1,i2,i3)*b1L(i1,i2,i3))*KLR*DLR+cssR(j1,j2,j3)*(n1R*u1DR+n2R*u2DR))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
    L3rCoeff(i1,i2,i3) = (((-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1rCoeffs(i1,i2,i3)+(-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*L2rCoeff(i1,i2,i3)-1.*crsR(j1,j2,j3)*L2rCoeffs(i1,i2,i3)-1.*cssR(j1,j2,j3)*L1rCoeffss(i1,i2,i3)-1.*crRr(j1,j2,j3)*L1rCoeff(i1,i2,i3))*b1R(j1,j2,j3)+((-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*b2R(j1,j2,j3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crsR(j1,j2,j3))*L1rCoeffs(i1,i2,i3)+(-1.*crsR(j1,j2,j3)*L1rCoeffss(i1,i2,i3)-1.*crRs(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L2rCoeff(i1,i2,i3)-1.*L2rCoeffs(i1,i2,i3)*crrR(j1,j2,j3))*b2R(j1,j2,j3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crrR(j1,j2,j3)*L2rCoeff(i1,i2,i3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crR(j1,j2,j3)*L1rCoeff(i1,i2,i3)+(-1.*crL(i1,i2,i3)*n1R*u1DL-1.*crL(i1,i2,i3)*n2R*u2DL+1.*crLr(i1,i2,i3)*b1L(i1,i2,i3)+1.*crLs(i1,i2,i3)*b2L(i1,i2,i3))*KLR*DLR)/b1R(j1,j2,j3)/crrR(j1,j2,j3);
    L3sCoeff(i1,i2,i3) = (((-1.*L2ICoeff(i1,i2,i3)-1.*L2sCoeffs(i1,i2,i3))*crsR(j1,j2,j3)+(-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1sCoeffs(i1,i2,i3)+(-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1ICoeff(i1,i2,i3)+(-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*L2sCoeff(i1,i2,i3)-1.*csRr(j1,j2,j3)-1.*cssR(j1,j2,j3)*L1sCoeffss(i1,i2,i3)-2.*cssR(j1,j2,j3)*L1ICoeffs(i1,i2,i3)-1.*crRr(j1,j2,j3)*L1sCoeff(i1,i2,i3))*b1R(j1,j2,j3)+((-1.*L1sCoeffss(i1,i2,i3)-2.*L1ICoeffs(i1,i2,i3))*crsR(j1,j2,j3)+(-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1sCoeffs(i1,i2,i3)+(-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1ICoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L2sCoeff(i1,i2,i3)+(-1.*L2ICoeff(i1,i2,i3)-1.*L2sCoeffs(i1,i2,i3))*crrR(j1,j2,j3)-1.*csRs(j1,j2,j3)-1.*crRs(j1,j2,j3)*L1sCoeff(i1,i2,i3))*b2R(j1,j2,j3)+(L1ICoeff(i1,i2,i3)+L1sCoeffs(i1,i2,i3))*(n1R*u1DR+n2R*u2DR)*crsR(j1,j2,j3)+(u1DR*L2sCoeff(i1,i2,i3)*crrR(j1,j2,j3)+(crR(j1,j2,j3)*L1sCoeff(i1,i2,i3)+csR(j1,j2,j3))*u1DR-1.*csL(i1,i2,i3)*u1DL*KLR*DLR)*n1R+(u2DR*L2sCoeff(i1,i2,i3)*crrR(j1,j2,j3)+(crR(j1,j2,j3)*L1sCoeff(i1,i2,i3)+csR(j1,j2,j3))*u2DR-1.*csL(i1,i2,i3)*u2DL*KLR*DLR)*n2R+DLR*KLR*(b1L(i1,i2,i3)*csLr(i1,i2,i3)+b2L(i1,i2,i3)*csLs(i1,i2,i3)))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
    L3ICoeff(i1,i2,i3) = (((-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1ICoeffs(i1,i2,i3)+(-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*L2ICoeff(i1,i2,i3)-1.*crRr(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*crsR(j1,j2,j3)*L2ICoeffs(i1,i2,i3)-1.*cssR(j1,j2,j3)*L1ICoeffss(i1,i2,i3))*b1R(j1,j2,j3)+((-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*b2R(j1,j2,j3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crsR(j1,j2,j3))*L1ICoeffs(i1,i2,i3)+(-1.*L1ICoeff(i1,i2,i3)*crRs(j1,j2,j3)-1.*L2ICoeffs(i1,i2,i3)*crrR(j1,j2,j3)-1.*crrRs(j1,j2,j3)*L2ICoeff(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1ICoeffss(i1,i2,i3))*b2R(j1,j2,j3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crrR(j1,j2,j3)*L2ICoeff(i1,i2,i3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crR(j1,j2,j3)*L1ICoeff(i1,i2,i3))/b1R(j1,j2,j3)/crrR(j1,j2,j3);
} // end FOR_3D

// ----- Evaluate coefficients of L4 on the boundary -----
extra=0;
getBoundaryIndex(mg.indexRange(), side,axis  ,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L4rrrrCoeff(Ib1,Ib2,Ib3);
RealArray L4rrrsCoeff(Ib1,Ib2,Ib3);
RealArray L4rrssCoeff(Ib1,Ib2,Ib3);
RealArray L4rsssCoeff(Ib1,Ib2,Ib3);
RealArray L4ssssCoeff(Ib1,Ib2,Ib3);
RealArray L4rrrCoeff(Ib1,Ib2,Ib3);
RealArray L4rrsCoeff(Ib1,Ib2,Ib3);
RealArray L4rssCoeff(Ib1,Ib2,Ib3);
RealArray L4sssCoeff(Ib1,Ib2,Ib3);
RealArray L4rrCoeff(Ib1,Ib2,Ib3);
RealArray L4rsCoeff(Ib1,Ib2,Ib3);
RealArray L4ssCoeff(Ib1,Ib2,Ib3);
RealArray L4rCoeff(Ib1,Ib2,Ib3);
RealArray L4sCoeff(Ib1,Ib2,Ib3);
RealArray L4ICoeff(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  // Compute tangential derivatives of L1
    L1rCoeffsss(i1,i2,i3) = (((-1./2.)*L1rCoeff(i1,i2-2,i3)+1.*L1rCoeff(i1,i2-1,i3)-1.*L1rCoeff(i1,i2+1,i3)+(1./2.)*L1rCoeff(i1,i2+2,i3))/(pow(drL[1],3)));
    L1sCoeffsss(i1,i2,i3) = (((-1./2.)*L1sCoeff(i1,i2-2,i3)+1.*L1sCoeff(i1,i2-1,i3)-1.*L1sCoeff(i1,i2+1,i3)+(1./2.)*L1sCoeff(i1,i2+2,i3))/(pow(drL[1],3)));
    L1ICoeffsss(i1,i2,i3) = (((-1./2.)*L1ICoeff(i1,i2-2,i3)+1.*L1ICoeff(i1,i2-1,i3)-1.*L1ICoeff(i1,i2+1,i3)+(1./2.)*L1ICoeff(i1,i2+2,i3))/(pow(drL[1],3)));
  // Compute tangential derivatives of L2
    L2rrCoeffss(i1,i2,i3) = ((L2rrCoeff(i1,i2+1,i3)-2.*L2rrCoeff(i1,i2,i3)+L2rrCoeff(i1,i2-1,i3))/(SQR(drL[1])));
    L2rsCoeffss(i1,i2,i3) = ((L2rsCoeff(i1,i2+1,i3)-2.*L2rsCoeff(i1,i2,i3)+L2rsCoeff(i1,i2-1,i3))/(SQR(drL[1])));
    L2ssCoeffss(i1,i2,i3) = ((L2ssCoeff(i1,i2+1,i3)-2.*L2ssCoeff(i1,i2,i3)+L2ssCoeff(i1,i2-1,i3))/(SQR(drL[1])));
    L2rCoeffss(i1,i2,i3)  = ((L2rCoeff(i1,i2+1,i3)-2.*L2rCoeff(i1,i2,i3)+L2rCoeff(i1,i2-1,i3))/(SQR(drL[1])));
    L2sCoeffss(i1,i2,i3)  = ((L2sCoeff(i1,i2+1,i3)-2.*L2sCoeff(i1,i2,i3)+L2sCoeff(i1,i2-1,i3))/(SQR(drL[1])));
    L2ICoeffss(i1,i2,i3)  = ((L2ICoeff(i1,i2+1,i3)-2.*L2ICoeff(i1,i2,i3)+L2ICoeff(i1,i2-1,i3))/(SQR(drL[1])));
  // Compute tangential derivatives of L3
    L3rrrCoeffs(i1,i2,i3) = ((L3rrrCoeff(i1,i2+1,i3)-L3rrrCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L3rrsCoeffs(i1,i2,i3) = ((L3rrsCoeff(i1,i2+1,i3)-L3rrsCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L3rssCoeffs(i1,i2,i3) = ((L3rssCoeff(i1,i2+1,i3)-L3rssCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L3sssCoeffs(i1,i2,i3) = ((L3sssCoeff(i1,i2+1,i3)-L3sssCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L3rrCoeffs(i1,i2,i3)  = ((L3rrCoeff(i1,i2+1,i3)-L3rrCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L3rsCoeffs(i1,i2,i3)  = ((L3rsCoeff(i1,i2+1,i3)-L3rsCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L3ssCoeffs(i1,i2,i3)  = ((L3ssCoeff(i1,i2+1,i3)-L3ssCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L3rCoeffs(i1,i2,i3)   = ((L3rCoeff(i1,i2+1,i3)-L3rCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L3sCoeffs(i1,i2,i3)   = ((L3sCoeff(i1,i2+1,i3)-L3sCoeff(i1,i2-1,i3))/(2.*drL[1]));
    L3ICoeffs(i1,i2,i3)   = ((L3ICoeff(i1,i2+1,i3)-L3ICoeff(i1,i2-1,i3))/(2.*drL[1]));
  // Compute L4 coefficients
    L4rrrrCoeff(i1,i2,i3) = 1.*pow(DLR,2.)*pow(crrL(i1,i2,i3),2.)/pow(crrR(j1,j2,j3),2.);
    L4rrrsCoeff(i1,i2,i3) = (2.*pow(DLR,2.)*crrL(i1,i2,i3)*crsL(i1,i2,i3)-2.*crrR(j1,j2,j3)*crsR(j1,j2,j3)*L3rrrCoeff(i1,i2,i3))/pow(crrR(j1,j2,j3),2.);
    L4rrssCoeff(i1,i2,i3) = (2.*pow(DLR,2.)*crrL(i1,i2,i3)*cssL(i1,i2,i3)+pow(DLR,2.)*pow(crsL(i1,i2,i3),2.)-2.*crrR(j1,j2,j3)*crsR(j1,j2,j3)*L3rrsCoeff(i1,i2,i3)-2.*crrR(j1,j2,j3)*cssR(j1,j2,j3)*L2rrCoeff(i1,i2,i3)-pow(crsR(j1,j2,j3),2.)*L2rrCoeff(i1,i2,i3))/pow(crrR(j1,j2,j3),2.);
    L4rsssCoeff(i1,i2,i3) = (2.*pow(DLR,2.)*crsL(i1,i2,i3)*cssL(i1,i2,i3)-2.*crrR(j1,j2,j3)*crsR(j1,j2,j3)*L3rssCoeff(i1,i2,i3)-2.*crrR(j1,j2,j3)*cssR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)-pow(crsR(j1,j2,j3),2.)*L2rsCoeff(i1,i2,i3)-2.*crsR(j1,j2,j3)*cssR(j1,j2,j3)*L1rCoeff(i1,i2,i3))/pow(crrR(j1,j2,j3),2.);
    L4ssssCoeff(i1,i2,i3) = (pow(DLR,2.)*pow(cssL(i1,i2,i3),2.)-2.*crrR(j1,j2,j3)*crsR(j1,j2,j3)*L3sssCoeff(i1,i2,i3)-2.*crrR(j1,j2,j3)*cssR(j1,j2,j3)*L2ssCoeff(i1,i2,i3)-pow(crsR(j1,j2,j3),2.)*L2ssCoeff(i1,i2,i3)-2.*crsR(j1,j2,j3)*cssR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-pow(cssR(j1,j2,j3),2.))/pow(crrR(j1,j2,j3),2.);
    L4rrrCoeff(i1,i2,i3) = (((-2.*crR(j1,j2,j3)-2.*crrRr(j1,j2,j3))*L3rrrCoeff(i1,i2,i3)-2.*crsR(j1,j2,j3)*L3rrrCoeffs(i1,i2,i3))*crrR(j1,j2,j3)-1.*crrRs(j1,j2,j3)*crsR(j1,j2,j3)*L3rrrCoeff(i1,i2,i3)+((2.*crL(i1,i2,i3)+2.*crrLr(i1,i2,i3))*crrL(i1,i2,i3)+crrLs(i1,i2,i3)*crsL(i1,i2,i3))*pow(DLR,2.))/pow(crrR(j1,j2,j3),2.);
    L4rrsCoeff(i1,i2,i3) = (((-2.*L3rrCoeff(i1,i2,i3)-2.*L3rrsCoeffs(i1,i2,i3))*crsR(j1,j2,j3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2rrCoeff(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crrRr(j1,j2,j3))*L3rrsCoeff(i1,i2,i3)-4.*cssR(j1,j2,j3)*L2rrCoeffs(i1,i2,i3))*crrR(j1,j2,j3)-2.*pow(crsR(j1,j2,j3),2.)*L2rrCoeffs(i1,i2,i3)+((-2.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2rrCoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3rrsCoeff(i1,i2,i3))*crsR(j1,j2,j3)-2.*crrRs(j1,j2,j3)*cssR(j1,j2,j3)*L2rrCoeff(i1,i2,i3)+((2.*crL(i1,i2,i3)+crrLr(i1,i2,i3)+crsLs(i1,i2,i3))*crsL(i1,i2,i3)+(2.*crsLr(i1,i2,i3)+2.*csL(i1,i2,i3))*crrL(i1,i2,i3)+2.*crrLs(i1,i2,i3)*cssL(i1,i2,i3))*pow(DLR,2.))/pow(crrR(j1,j2,j3),2.);
    L4rssCoeff(i1,i2,i3) = ((-1.*L2rCoeff(i1,i2,i3)-2.*L2rsCoeffs(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3rsCoeff(i1,i2,i3)-2.*L3rssCoeffs(i1,i2,i3))*crrR(j1,j2,j3)+(-2.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2rsCoeff(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*L1rCoeff(i1,i2,i3)-1.*L3rssCoeff(i1,i2,i3)*crrRs(j1,j2,j3)-6.*cssR(j1,j2,j3)*L1rCoeffs(i1,i2,i3))*crsR(j1,j2,j3)+((-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2rsCoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L1rCoeff(i1,i2,i3)+(-2.*L2rCoeff(i1,i2,i3)-4.*L2rsCoeffs(i1,i2,i3))*cssR(j1,j2,j3)+(-2.*crR(j1,j2,j3)-2.*crrRr(j1,j2,j3))*L3rssCoeff(i1,i2,i3))*crrR(j1,j2,j3)-2.*crrRs(j1,j2,j3)*cssR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3))*cssR(j1,j2,j3)*L1rCoeff(i1,i2,i3)+((crsLr(i1,i2,i3)+cssLs(i1,i2,i3)+2.*csL(i1,i2,i3))*crsL(i1,i2,i3)+(2.*crL(i1,i2,i3)+2.*crsLs(i1,i2,i3))*cssL(i1,i2,i3)+2.*crrL(i1,i2,i3)*cssLr(i1,i2,i3))*pow(DLR,2.))/pow(crrR(j1,j2,j3),2.);
    L4sssCoeff(i1,i2,i3) = ((-1.*L2sCoeff(i1,i2,i3)-2.*L2ssCoeffs(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3ssCoeff(i1,i2,i3)-2.*L3sssCoeffs(i1,i2,i3))*crrR(j1,j2,j3)+(-6.*L1sCoeffs(i1,i2,i3)-2.*L1ICoeff(i1,i2,i3))*cssR(j1,j2,j3)+(-2.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2ssCoeff(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*L1sCoeff(i1,i2,i3)-1.*cssRr(j1,j2,j3)-1.*crrRs(j1,j2,j3)*L3sssCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-2.*L2sCoeff(i1,i2,i3)-4.*L2ssCoeffs(i1,i2,i3))*cssR(j1,j2,j3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2ssCoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L1sCoeff(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crrRr(j1,j2,j3))*L3sssCoeff(i1,i2,i3))*crrR(j1,j2,j3)+(-2.*crrRs(j1,j2,j3)*L2ssCoeff(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3))*L1sCoeff(i1,i2,i3)-2.*cssRs(j1,j2,j3)-2.*csR(j1,j2,j3))*cssR(j1,j2,j3)+((2.*csL(i1,i2,i3)+2.*cssLs(i1,i2,i3))*cssL(i1,i2,i3)+crsL(i1,i2,i3)*cssLr(i1,i2,i3))*pow(DLR,2.))/pow(crrR(j1,j2,j3),2.);
    L4rrCoeff(i1,i2,i3) = (((-2.*crRr(j1,j2,j3)-1.*crrRrr(j1,j2,j3))*L2rrCoeff(i1,i2,i3)-2.*L3rrCoeffs(i1,i2,i3)*crsR(j1,j2,j3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2rrCoeffs(i1,i2,i3)-2.*crR(j1,j2,j3)*L3rrCoeff(i1,i2,i3)-2.*crrRr(j1,j2,j3)*L3rrCoeff(i1,i2,i3)-2.*cssR(j1,j2,j3)*L2rrCoeffss(i1,i2,i3))*crrR(j1,j2,j3)+((-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*crsR(j1,j2,j3)-1.*crrRr(j1,j2,j3)*crR(j1,j2,j3)-1.*cssR(j1,j2,j3)*crrRss(j1,j2,j3)-1.*csR(j1,j2,j3)*crrRs(j1,j2,j3)-1.*pow(crR(j1,j2,j3),2.))*L2rrCoeff(i1,i2,i3)-1.*pow(crsR(j1,j2,j3),2.)*L2rrCoeffss(i1,i2,i3)+((-2.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2rrCoeffs(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3rrCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((crrLrr(i1,i2,i3)+2.*crLr(i1,i2,i3))*crrL(i1,i2,i3)+(crrLrs(i1,i2,i3)+crLs(i1,i2,i3))*crsL(i1,i2,i3)+cssL(i1,i2,i3)*crrLss(i1,i2,i3)+crrLr(i1,i2,i3)*crL(i1,i2,i3)+csL(i1,i2,i3)*crrLs(i1,i2,i3)+pow(crL(i1,i2,i3),2.))*pow(DLR,2.)-2.*cssR(j1,j2,j3)*crrRs(j1,j2,j3)*L2rrCoeffs(i1,i2,i3))/pow(crrR(j1,j2,j3),2.);
    L4rsCoeff(i1,i2,i3) = ((-2.*L2rCoeffs(i1,i2,i3)-1.*L2rsCoeffss(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3rsCoeffs(i1,i2,i3)-2.*L3rCoeff(i1,i2,i3))*crrR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1rCoeff(i1,i2,i3)-6.*cssR(j1,j2,j3)*L1rCoeffss(i1,i2,i3)+(-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*L2rsCoeff(i1,i2,i3)+(-2.*L2rsCoeffs(i1,i2,i3)-2.*L2rCoeff(i1,i2,i3))*crR(j1,j2,j3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2rsCoeffs(i1,i2,i3)+(-2.*crsRr(j1,j2,j3)-2.*cssRs(j1,j2,j3)-4.*csR(j1,j2,j3))*L1rCoeffs(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2rCoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3rsCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1rCoeff(i1,i2,i3)+(-2.*L2rsCoeffss(i1,i2,i3)-4.*L2rCoeffs(i1,i2,i3))*cssR(j1,j2,j3)+(-2.*crRr(j1,j2,j3)-1.*crrRrr(j1,j2,j3))*L2rsCoeff(i1,i2,i3)-2.*crR(j1,j2,j3)*L3rsCoeff(i1,i2,i3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2rsCoeffs(i1,i2,i3)-4.*cssRr(j1,j2,j3)*L1rCoeffs(i1,i2,i3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2rCoeff(i1,i2,i3)-2.*crrRr(j1,j2,j3)*L3rsCoeff(i1,i2,i3))*crrR(j1,j2,j3)+((-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*cssR(j1,j2,j3)+(-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3)*csR(j1,j2,j3))*L1rCoeff(i1,i2,i3)+(-4.*crR(j1,j2,j3)*L1rCoeffs(i1,i2,i3)-1.*crrRss(j1,j2,j3)*L2rsCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2rCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2rsCoeffs(i1,i2,i3)-4.*crsRs(j1,j2,j3)*L1rCoeffs(i1,i2,i3))*cssR(j1,j2,j3)+((crsLrs(i1,i2,i3)+crLr(i1,i2,i3)+csLs(i1,i2,i3))*crsL(i1,i2,i3)+(crsLrr(i1,i2,i3)+2.*csLr(i1,i2,i3))*crrL(i1,i2,i3)+(crsLss(i1,i2,i3)+2.*crLs(i1,i2,i3))*cssL(i1,i2,i3)+(crsLr(i1,i2,i3)+2.*csL(i1,i2,i3))*crL(i1,i2,i3)+csL(i1,i2,i3)*crsLs(i1,i2,i3))*pow(DLR,2.)+(-1.*crrRr(j1,j2,j3)*crR(j1,j2,j3)-1.*csR(j1,j2,j3)*crrRs(j1,j2,j3)-1.*pow(crR(j1,j2,j3),2.))*L2rsCoeff(i1,i2,i3))/pow(crrR(j1,j2,j3),2.);
    L4ssCoeff(i1,i2,i3) = ((-1.*L2ssCoeffss(i1,i2,i3)-2.*L2sCoeffs(i1,i2,i3)-1.*L2ICoeff(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3ssCoeffs(i1,i2,i3)-2.*L3sCoeff(i1,i2,i3))*crrR(j1,j2,j3)+(-6.*L1ICoeffs(i1,i2,i3)-6.*L1sCoeffss(i1,i2,i3))*cssR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1sCoeff(i1,i2,i3)+(-2.*L2sCoeff(i1,i2,i3)-2.*L2ssCoeffs(i1,i2,i3))*crR(j1,j2,j3)+(-4.*L1sCoeffs(i1,i2,i3)-2.*L1ICoeff(i1,i2,i3))*csR(j1,j2,j3)+(-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*L2ssCoeff(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2ssCoeffs(i1,i2,i3)+(-2.*cssRs(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L1sCoeffs(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3))*L1ICoeff(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2sCoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3ssCoeff(i1,i2,i3)-1.*cssRrs(j1,j2,j3)-1.*csRr(j1,j2,j3))*crsR(j1,j2,j3)+((-2.*L2ICoeff(i1,i2,i3)-2.*L2ssCoeffss(i1,i2,i3)-4.*L2sCoeffs(i1,i2,i3))*cssR(j1,j2,j3)+(-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1sCoeff(i1,i2,i3)-2.*crR(j1,j2,j3)*L3ssCoeff(i1,i2,i3)+(-2.*L2sCoeff(i1,i2,i3)-2.*L2ssCoeffs(i1,i2,i3))*csR(j1,j2,j3)+(-2.*crRr(j1,j2,j3)-1.*crrRrr(j1,j2,j3))*L2ssCoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L1ICoeff(i1,i2,i3)-4.*cssRr(j1,j2,j3)*L1sCoeffs(i1,i2,i3)-1.*cssRrr(j1,j2,j3)-2.*crrRr(j1,j2,j3)*L3ssCoeff(i1,i2,i3)-2.*crsRr(j1,j2,j3)*L2ssCoeffs(i1,i2,i3)-2.*crsRr(j1,j2,j3)*L2sCoeff(i1,i2,i3))*crrR(j1,j2,j3)+((-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*L1sCoeff(i1,i2,i3)+(-4.*L1sCoeffs(i1,i2,i3)-2.*L1ICoeff(i1,i2,i3))*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3)*L1ICoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2ssCoeffs(i1,i2,i3)-2.*csRs(j1,j2,j3)-4.*crsRs(j1,j2,j3)*L1sCoeffs(i1,i2,i3)-1.*crrRss(j1,j2,j3)*L2ssCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2sCoeff(i1,i2,i3)-1.*cssRss(j1,j2,j3))*cssR(j1,j2,j3)+((-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3)*csR(j1,j2,j3))*L1sCoeff(i1,i2,i3)-1.*pow(crR(j1,j2,j3),2.)*L2ssCoeff(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)*L2ssCoeff(i1,i2,i3)-1.*cssRr(j1,j2,j3))*crR(j1,j2,j3)-1.*pow(csR(j1,j2,j3),2.)+(-1.*crrRs(j1,j2,j3)*L2ssCoeff(i1,i2,i3)-1.*cssRs(j1,j2,j3))*csR(j1,j2,j3)+((cssLrs(i1,i2,i3)+csLr(i1,i2,i3))*crsL(i1,i2,i3)+(cssLss(i1,i2,i3)+2.*csLs(i1,i2,i3))*cssL(i1,i2,i3)+cssLr(i1,i2,i3)*crL(i1,i2,i3)+csL(i1,i2,i3)*cssLs(i1,i2,i3)+crrL(i1,i2,i3)*cssLrr(i1,i2,i3)+pow(csL(i1,i2,i3),2.))*pow(DLR,2.))/pow(crrR(j1,j2,j3),2.);
    L4rCoeff(i1,i2,i3) = (-1.*pow(crsR(j1,j2,j3),2.)*L2rCoeffss(i1,i2,i3)+(-2.*L3rCoeffs(i1,i2,i3)*crrR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1rCoeffs(i1,i2,i3)-2.*cssR(j1,j2,j3)*L1rCoeffsss(i1,i2,i3)+(-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*L2rCoeff(i1,i2,i3)-2.*crR(j1,j2,j3)*L2rCoeffs(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*L1rCoeffss(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2rCoeffs(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3rCoeff(i1,i2,i3)-1.*crRrs(j1,j2,j3)*L1rCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1rCoeffs(i1,i2,i3)-2.*cssR(j1,j2,j3)*L2rCoeffss(i1,i2,i3)+(-2.*crRr(j1,j2,j3)-1.*crrRrr(j1,j2,j3))*L2rCoeff(i1,i2,i3)-2.*crR(j1,j2,j3)*L3rCoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L1rCoeffss(i1,i2,i3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2rCoeffs(i1,i2,i3)-2.*crrRr(j1,j2,j3)*L3rCoeff(i1,i2,i3)-1.*crRrr(j1,j2,j3)*L1rCoeff(i1,i2,i3))*crrR(j1,j2,j3)+((-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*cssR(j1,j2,j3)+(-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3)*csR(j1,j2,j3))*L1rCoeffs(i1,i2,i3)+(-2.*crR(j1,j2,j3)*L1rCoeffss(i1,i2,i3)-1.*crrRss(j1,j2,j3)*L2rCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2rCoeffs(i1,i2,i3)-2.*crsRs(j1,j2,j3)*L1rCoeffss(i1,i2,i3)-1.*crRss(j1,j2,j3)*L1rCoeff(i1,i2,i3))*cssR(j1,j2,j3)+(-1.*crrRr(j1,j2,j3)*crR(j1,j2,j3)-1.*csR(j1,j2,j3)*crrRs(j1,j2,j3)-1.*pow(crR(j1,j2,j3),2.))*L2rCoeff(i1,i2,i3)-1.*crRr(j1,j2,j3)*crR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*csR(j1,j2,j3)*crRs(j1,j2,j3)*L1rCoeff(i1,i2,i3)+pow(DLR,2.)*(crL(i1,i2,i3)*crLr(i1,i2,i3)+crLrr(i1,i2,i3)*crrL(i1,i2,i3)+crLrs(i1,i2,i3)*crsL(i1,i2,i3)+crLs(i1,i2,i3)*csL(i1,i2,i3)+crLss(i1,i2,i3)*cssL(i1,i2,i3)))/pow(crrR(j1,j2,j3),2.);
    L4sCoeff(i1,i2,i3) = ((-1.*L2sCoeffss(i1,i2,i3)-2.*L2ICoeffs(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3sCoeffs(i1,i2,i3)-2.*L3ICoeff(i1,i2,i3))*crrR(j1,j2,j3)+(-6.*L1ICoeffss(i1,i2,i3)-2.*L1sCoeffsss(i1,i2,i3))*cssR(j1,j2,j3)+(-2.*L2ICoeff(i1,i2,i3)-2.*L2sCoeffs(i1,i2,i3))*crR(j1,j2,j3)+(-2.*L1sCoeffss(i1,i2,i3)-4.*L1ICoeffs(i1,i2,i3))*csR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1sCoeffs(i1,i2,i3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1ICoeff(i1,i2,i3)+(-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*L2sCoeff(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3))*L1sCoeffss(i1,i2,i3)+(-2.*cssRs(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L1ICoeffs(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2sCoeffs(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2ICoeff(i1,i2,i3)-1.*csRrs(j1,j2,j3)-1.*crRrs(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3sCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-2.*L2sCoeffss(i1,i2,i3)-4.*L2ICoeffs(i1,i2,i3))*cssR(j1,j2,j3)-2.*crR(j1,j2,j3)*L3sCoeff(i1,i2,i3)+(-2.*L2ICoeff(i1,i2,i3)-2.*L2sCoeffs(i1,i2,i3))*csR(j1,j2,j3)+(-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1sCoeffs(i1,i2,i3)+(-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1ICoeff(i1,i2,i3)+(-2.*crRr(j1,j2,j3)-1.*crrRrr(j1,j2,j3))*L2sCoeff(i1,i2,i3)-1.*csRrr(j1,j2,j3)-2.*cssRr(j1,j2,j3)*L1sCoeffss(i1,i2,i3)-4.*cssRr(j1,j2,j3)*L1ICoeffs(i1,i2,i3)-1.*crRrr(j1,j2,j3)*L1sCoeff(i1,i2,i3)-2.*crsRr(j1,j2,j3)*L2ICoeff(i1,i2,i3)-2.*crrRr(j1,j2,j3)*L3sCoeff(i1,i2,i3)-2.*crsRr(j1,j2,j3)*L2sCoeffs(i1,i2,i3))*crrR(j1,j2,j3)+((-2.*L1sCoeffss(i1,i2,i3)-4.*L1ICoeffs(i1,i2,i3))*crR(j1,j2,j3)+(-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*L1sCoeffs(i1,i2,i3)+(-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*L1ICoeff(i1,i2,i3)-1.*csRss(j1,j2,j3)-1.*crrRss(j1,j2,j3)*L2sCoeff(i1,i2,i3)-1.*crRss(j1,j2,j3)*L1sCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2sCoeffs(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2ICoeff(i1,i2,i3)-2.*crsRs(j1,j2,j3)*L1sCoeffss(i1,i2,i3)-4.*crsRs(j1,j2,j3)*L1ICoeffs(i1,i2,i3))*cssR(j1,j2,j3)-1.*pow(crR(j1,j2,j3),2.)*L2sCoeff(i1,i2,i3)+((-2.*L1sCoeffs(i1,i2,i3)-2.*L1ICoeff(i1,i2,i3))*csR(j1,j2,j3)-1.*csRr(j1,j2,j3)-1.*crsRr(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*crsRr(j1,j2,j3)*L1sCoeffs(i1,i2,i3)-1.*crRr(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*crrRr(j1,j2,j3)*L2sCoeff(i1,i2,i3))*crR(j1,j2,j3)+(-1.*csRs(j1,j2,j3)-1.*crrRs(j1,j2,j3)*L2sCoeff(i1,i2,i3)-1.*crRs(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*crsRs(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*crsRs(j1,j2,j3)*L1sCoeffs(i1,i2,i3))*csR(j1,j2,j3)+pow(DLR,2.)*(crL(i1,i2,i3)*csLr(i1,i2,i3)+crrL(i1,i2,i3)*csLrr(i1,i2,i3)+crsL(i1,i2,i3)*csLrs(i1,i2,i3)+csL(i1,i2,i3)*csLs(i1,i2,i3)+csLss(i1,i2,i3)*cssL(i1,i2,i3)))/pow(crrR(j1,j2,j3),2.);
    L4ICoeff(i1,i2,i3) = (-1.*pow(crsR(j1,j2,j3),2.)*L2ICoeffss(i1,i2,i3)+(-2.*L3ICoeffs(i1,i2,i3)*crrR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1ICoeffs(i1,i2,i3)-2.*cssR(j1,j2,j3)*L1ICoeffsss(i1,i2,i3)+(-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*L2ICoeff(i1,i2,i3)-2.*crR(j1,j2,j3)*L2ICoeffs(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*L1ICoeffss(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L2ICoeffs(i1,i2,i3)-1.*crrRs(j1,j2,j3)*L3ICoeff(i1,i2,i3)-1.*crRrs(j1,j2,j3)*L1ICoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1ICoeffs(i1,i2,i3)-2.*cssR(j1,j2,j3)*L2ICoeffss(i1,i2,i3)+(-2.*crRr(j1,j2,j3)-1.*crrRrr(j1,j2,j3))*L2ICoeff(i1,i2,i3)-2.*crR(j1,j2,j3)*L3ICoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L1ICoeffss(i1,i2,i3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L2ICoeffs(i1,i2,i3)-2.*crrRr(j1,j2,j3)*L3ICoeff(i1,i2,i3)-1.*crRrr(j1,j2,j3)*L1ICoeff(i1,i2,i3))*crrR(j1,j2,j3)+((-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*cssR(j1,j2,j3)+(-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3)*csR(j1,j2,j3))*L1ICoeffs(i1,i2,i3)+(-1.*crrRss(j1,j2,j3)*L2ICoeff(i1,i2,i3)-2.*crsRs(j1,j2,j3)*L1ICoeffss(i1,i2,i3)-1.*crRss(j1,j2,j3)*L1ICoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L2ICoeffs(i1,i2,i3)-2.*crR(j1,j2,j3)*L1ICoeffss(i1,i2,i3))*cssR(j1,j2,j3)+(-1.*crrRr(j1,j2,j3)*crR(j1,j2,j3)-1.*csR(j1,j2,j3)*crrRs(j1,j2,j3)-1.*pow(crR(j1,j2,j3),2.))*L2ICoeff(i1,i2,i3)-1.*crRr(j1,j2,j3)*crR(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*csR(j1,j2,j3)*crRs(j1,j2,j3)*L1ICoeff(i1,i2,i3))/pow(crrR(j1,j2,j3),2.);
} // end FOR_3D

// ----- Define coefficients in difference operators -----
Range R4(-2,2);
RealArray iCoeff(R4);
RealArray rCoeff(R4), rrCoeff(R4), rrrCoeff(R4), rrrrCoeff(R4);
RealArray sCoeff(R4), ssCoeff(R4), sssCoeff(R4), ssssCoeff(R4);
      iCoeff(-2)= 0.;    iCoeff(-1)=  0.;    iCoeff(0)=  1.;    iCoeff(1)= 0.;    iCoeff(2)= 0.;
      rCoeff(-2)= 1.;    rCoeff(-1)= -8.;    rCoeff(0)=  0.;    rCoeff(1)= 8.;    rCoeff(2)=-1.; rCoeff    /=(    12.*drL[0]);
    rrCoeff(-2)=-1.;   rrCoeff(-1)= 16.;   rrCoeff(0)=-30.;   rrCoeff(1)=16.;   rrCoeff(2)=-1.; rrCoeff   /=(12.*SQR(drL[0]));
  rrrCoeff(-2)=-1.;  rrrCoeff(-1)=  2.;  rrrCoeff(0)=  0.;  rrrCoeff(1)=-2.;  rrrCoeff(2)= 1.; rrrCoeff  /=( 2.*pow(drL[0],3));
rrrrCoeff(-2)= 1.; rrrrCoeff(-1)= -4.; rrrrCoeff(0)=  6.; rrrrCoeff(1)=-4.; rrrrCoeff(2)= 1.; rrrrCoeff /=(    pow(drL[0],4));
      sCoeff(-2)= 1.;    sCoeff(-1)= -8.;    sCoeff(0)=  0.;    sCoeff(1)= 8.;    sCoeff(2)=-1.; sCoeff    /=(    12.*drL[1]);
    ssCoeff(-2)=-1.;   ssCoeff(-1)= 16.;   ssCoeff(0)=-30.;   ssCoeff(1)=16.;   ssCoeff(2)=-1.; ssCoeff   /=(12.*SQR(drL[1]));
  sssCoeff(-2)=-1.;  sssCoeff(-1)=  2.;  sssCoeff(0)=  0.;  sssCoeff(1)=-2.;  sssCoeff(2)= 1.; sssCoeff  /=( 2.*pow(drL[1],3));
ssssCoeff(-2)= 1.; ssssCoeff(-1)= -4.; ssssCoeff(0)=  6.; ssssCoeff(1)=-4.; ssssCoeff(2)= 1.; ssssCoeff /=(    pow(drL[1],4));

// ----- Fill in the Matrix Coefficients for CHAMP -----
const Real dxs = (1-2*side2)*drR[axis2];
const Real h = dxs;
printF("WDH: grid=%d, (side,axis)=(%d,%d) (side2,axis2)=(%d,%d) dxs=%9.3e, Sl=%9.3e\n",grid, side,axis, side2,axis2, dxs,Sl);
RealArray coeff4(R4,R4);
const int e=0, c=0; // eqn number and component number
// ------- FILL CHAMP conditions into the matrix -----
// NOTE: skip top boundary if periodic (use indexRange) 
// NOTE: skip adjacent boundaries if Dirichlet BC *finish me**
extra=0; 
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
int axisp = (axis+1) % numberOfDimensions;
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    i1m=i1-is1, i2m=i2-is2, i3m=i3-is3; //  ghost point is (i1m,i2m,i3m)
    const Real bn = b1R(j1,j2,j3), bt = b2R(j1,j2,j3);  // b in normal and tangential directions
  // The next macro defines cI, cr, cs, crr, ...
  // ------  ORDER = 4 , NormalDirection = 0,------
    const Real h2By2=h*h/2., h3By6=h2By2*h/3., h4By24=h3By6*h/4.;
    Real L1IC = L1ICoeff(i1,i2,i3);
    Real L1ICs = L1ICoeffs(i1,i2,i3);
    Real L1ICss = L1ICoeffss(i1,i2,i3);
    Real L1ICsss = L1ICoeffsss(i1,i2,i3);
    Real L1rC = L1rCoeff(i1,i2,i3);
    Real L1rCs = L1rCoeffs(i1,i2,i3);
    Real L1rCss = L1rCoeffss(i1,i2,i3);
    Real L1rCsss = L1rCoeffsss(i1,i2,i3);
    Real L1sC = L1sCoeff(i1,i2,i3);
    Real L1sCs = L1sCoeffs(i1,i2,i3);
    Real L1sCss = L1sCoeffss(i1,i2,i3);
    Real L1sCsss = L1sCoeffsss(i1,i2,i3);
    Real L2IC = L2ICoeff(i1,i2,i3);
    Real L2ICs = L2ICoeffs(i1,i2,i3);
    Real L2ICss = L2ICoeffss(i1,i2,i3);
    Real L2rC = L2rCoeff(i1,i2,i3);
    Real L2rCs = L2rCoeffs(i1,i2,i3);
    Real L2rCss = L2rCoeffss(i1,i2,i3);
    Real L2sC = L2sCoeff(i1,i2,i3);
    Real L2sCs = L2sCoeffs(i1,i2,i3);
    Real L2sCss = L2sCoeffss(i1,i2,i3);
    Real L2rrC = L2rrCoeff(i1,i2,i3);
    Real L2rrCs = L2rrCoeffs(i1,i2,i3);
    Real L2rrCss = L2rrCoeffss(i1,i2,i3);
    Real L2rsC = L2rsCoeff(i1,i2,i3);
    Real L2rsCs = L2rsCoeffs(i1,i2,i3);
    Real L2rsCss = L2rsCoeffss(i1,i2,i3);
    Real L2ssC = L2ssCoeff(i1,i2,i3);
    Real L2ssCs = L2ssCoeffs(i1,i2,i3);
    Real L2ssCss = L2ssCoeffss(i1,i2,i3);
    Real L3IC = L3ICoeff(i1,i2,i3);
    Real L3ICs = L3ICoeffs(i1,i2,i3);
    Real L3rC = L3rCoeff(i1,i2,i3);
    Real L3rCs = L3rCoeffs(i1,i2,i3);
    Real L3sC = L3sCoeff(i1,i2,i3);
    Real L3sCs = L3sCoeffs(i1,i2,i3);
    Real L3rrC = L3rrCoeff(i1,i2,i3);
    Real L3rrCs = L3rrCoeffs(i1,i2,i3);
    Real L3rsC = L3rsCoeff(i1,i2,i3);
    Real L3rsCs = L3rsCoeffs(i1,i2,i3);
    Real L3ssC = L3ssCoeff(i1,i2,i3);
    Real L3ssCs = L3ssCoeffs(i1,i2,i3);
    Real L3rrrC = L3rrrCoeff(i1,i2,i3);
    Real L3rrrCs = L3rrrCoeffs(i1,i2,i3);
    Real L3rrsC = L3rrsCoeff(i1,i2,i3);
    Real L3rrsCs = L3rrsCoeffs(i1,i2,i3);
    Real L3rssC = L3rssCoeff(i1,i2,i3);
    Real L3rssCs = L3rssCoeffs(i1,i2,i3);
    Real L3sssC = L3sssCoeff(i1,i2,i3);
    Real L3sssCs = L3sssCoeffs(i1,i2,i3);
    Real L4IC = L4ICoeff(i1,i2,i3);
    Real L4rC = L4rCoeff(i1,i2,i3);
    Real L4sC = L4sCoeff(i1,i2,i3);
    Real L4rrC = L4rrCoeff(i1,i2,i3);
    Real L4rsC = L4rsCoeff(i1,i2,i3);
    Real L4ssC = L4ssCoeff(i1,i2,i3);
    Real L4rrrC = L4rrrCoeff(i1,i2,i3);
    Real L4rrsC = L4rrsCoeff(i1,i2,i3);
    Real L4rssC = L4rssCoeff(i1,i2,i3);
    Real L4sssC = L4sssCoeff(i1,i2,i3);
    Real L4rrrrC = L4rrrrCoeff(i1,i2,i3);
    Real L4rrrsC = L4rrrsCoeff(i1,i2,i3);
    Real L4rrssC = L4rrssCoeff(i1,i2,i3);
    Real L4rsssC = L4rsssCoeff(i1,i2,i3);
    Real L4ssssC = L4ssssCoeff(i1,i2,i3);
    Real cI  = Sl*(h*L1IC+h2By2*L2IC+h3By6*L3IC+h4By24*L4IC+1)+(-h*L2IC-h2By2*L3IC-h3By6*L4IC-L1IC)*bn-bt*(h*L1ICs+h2By2*L2ICs+h3By6*L3ICs);
    Real cr  = Sl*(h*L1rC+h2By2*L2rC+h3By6*L3rC+h4By24*L4rC)+(-h*L2rC-h2By2*L3rC-h3By6*L4rC-L1rC)*bn-bt*(h*L1rCs+h2By2*L2rCs+h3By6*L3rCs);
    Real cs  = ((-L2IC-L2sCs)*h2By2+(-L3IC-L3sCs)*h3By6+(-L1sCs-L1IC)*h-L4IC*h4By24-2)*bt+(Sl*L2sC-bn*L3sC)*h2By2+(Sl*L3sC-bn*L4sC)*h3By6+(h*L1sC+h4By24*L4sC)*Sl-bn*(h*L2sC+L1sC);
    Real crr  = (Sl*L2rrC-bn*L3rrC-bt*L2rrCs)*h2By2+(Sl*L3rrC-bn*L4rrC-bt*L3rrCs)*h3By6-bn*L2rrC*h+Sl*h4By24*L4rrC;
    Real crs  = ((-L2rsCs-L2rC)*h2By2+(-L3rC-L3rsCs)*h3By6-L1rC*h-L4rC*h4By24)*bt+(Sl*L2rsC-bn*L3rsC)*h2By2+(Sl*L3rsC-bn*L4rsC)*h3By6-bn*L2rsC*h+Sl*h4By24*L4rsC;
    Real css  = ((-L2ssCs-L2sC)*h2By2+(-L3sC-L3ssCs)*h3By6-L1sC*h-L4sC*h4By24)*bt+(Sl*L2ssC-bn*L3ssC)*h2By2+(Sl*L3ssC-bn*L4ssC)*h3By6-bn*L2ssC*h+Sl*h4By24*L4ssC;
    Real crrr  = (Sl*L3rrrC-bn*L4rrrC-bt*L3rrrCs)*h3By6+Sl*h4By24*L4rrrC-bn*h2By2*L3rrrC;
    Real crrs  = ((-L3rrsCs-L3rrC)*bt+Sl*L3rrsC-L4rrsC*bn)*h3By6+(-h2By2*L2rrC-h4By24*L4rrC)*bt+Sl*h4By24*L4rrsC-bn*h2By2*L3rrsC;
    Real crss  = ((-L3rssCs-L3rsC)*bt+Sl*L3rssC-L4rssC*bn)*h3By6+(-h2By2*L2rsC-h4By24*L4rsC)*bt+Sl*h4By24*L4rssC-bn*h2By2*L3rssC;
    Real csss  = ((-L3sssCs-L3ssC)*bt+Sl*L3sssC-L4sssC*bn)*h3By6+(-h2By2*L2ssC-h4By24*L4ssC)*bt+Sl*h4By24*L4sssC-bn*h2By2*L3sssC;
    Real crrrr  = L4rrrrC*(Sl*h4By24-bn*h3By6);
    Real crrrs  = (-h3By6*L3rrrC-h4By24*L4rrrC)*bt+(Sl*h4By24-bn*h3By6)*L4rrrsC;
    Real crrss  = (-h3By6*L3rrsC-h4By24*L4rrsC)*bt+(Sl*h4By24-bn*h3By6)*L4rrssC;
    Real crsss  = (-h3By6*L3rssC-h4By24*L4rssC)*bt+(Sl*h4By24-bn*h3By6)*L4rsssC;
    Real cssss  = (-h3By6*L3sssC-h4By24*L4sssC)*bt+(Sl*h4By24-bn*h3By6)*L4ssssC;
//  if( orderOfAccuracy==4 )
//  {
//    printF("domain2=%d: cI   =%12.4e, cr   =%12.4e, cs   =%12.4e, crr  =%12.4e, crs=%12.4e, css=%12.4e\n",domain2,cI,cr,cs,crr,crs,css);
//    printF("            crrr =%12.4e, crrs =%12.4e, crss =%12.4e, csss =%12.4e\n",crrr,crrs,crss,csss);
//    printF("            crrrr=%12.4e, crrrs=%12.4e, crrss=%12.4e, crsss=%12.4e, cssss=%12.4e\n",crrrr,crrrs,crrss,crsss,cssss);
//    printF("            L4ssssC=%12.4e, (Sl*h4By24-bn*h3By6)=%12.4e h4By24=%12.4e h3By6=%12.4e, Sl=%12.4e, bn=%12.4e, h=%12.4e, bt=%12.4e\n",L4ssssC,(Sl*h4By24-bn*h3By6),h4By24,h3By6,Sl,bn,h,bt);
//    if( domain2==0 )
//    {
//      OV_ABORT("stop here for now");
//    }
//  }
    ForStencil(m1,m2,m3)
    {
        int m  = M123(m1,m2,m3);        // the single-component coeff-index
        int mm = M123CE(m1,m2,m3,c,e);  // the system coeff-index
        coeff4(m1,m2) = 
                          + cI   *   iCoeff(m1) * iCoeff(m2) 
                          + cr   *   rCoeff(m1) * iCoeff(m2) 
                          + cs   *   iCoeff(m1) * sCoeff(m2) 
                          + crr  *  rrCoeff(m1) * iCoeff(m2) 
                          + crs  *   rCoeff(m1) * sCoeff(m2) 
                          + css  *   iCoeff(m1) *ssCoeff(m2) 
                          + crrr   * rrrCoeff(m1)  *   iCoeff(m2) 
                          + crrs   *  rrCoeff(m1)  *   sCoeff(m2) 
                          + crss   *   rCoeff(m1)  *  ssCoeff(m2) 
                          + csss   *   iCoeff(m1)  * sssCoeff(m2) 
                          + crrrr  *rrrrCoeff(m1)  *   iCoeff(m2) 
                          + crrrs  * rrrCoeff(m1)  *   sCoeff(m2) 
                          + crrss  *  rrCoeff(m1)  *  ssCoeff(m2) 
                          + crsss  *   rCoeff(m1)  * sssCoeff(m2) 
                          + cssss  *   iCoeff(m1)  *ssssCoeff(m2) 
                          ;
        if( fillMatrixWDH )
        {
            coeff(mm,i1m,i2m,i3m) = coeff4(m1,m2);
      // Specify that the above coeff value is the coefficient of component c at the grid point (j1,j2,j3).
            const int k1=i1+m1, k2=i2+m2, k3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)    
            setEquationNumber(mm, e,i1m,i2m,i3m,  c,k1,k2,k3 );  // macro to set equationNumber                  
                                                                                                                                                                                                                      
       // Fill Ghost 2 -- extrapolation for now                                                            
              const int ghost=2;                                                                                  
                  const int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
         // --- fill in the coefficients of the extrapolation formula ---
                  for( int me=0; me<=extrapOrder; me++ )
                  {
                      coeff(me,i1m,i2m,i3m) = extrapCoeff[me];
                      const int j1=i1m + me*is1, j2=i2m + me*is2, j3=i3m + me*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                      setEquationNumber(me, e,i1m,i2m,i3m,  c,j1,j2,j3 );             // macro to set equationNumber
                  }                
        }
        if( twilightZoneFlow )
        {
      // --- For twilightZone we save some coefficients that go into the CHAMP matrix ---
            cc(i1,i2,i3, 0) = cI;   // coeff of I
            cc(i1,i2,i3, 1) = cr;   // coeff of ur
            cc(i1,i2,i3, 2) = cs;   // coeff of us
            cc(i1,i2,i3, 3) = crr;  // coeff of urr
            cc(i1,i2,i3, 4) = crs;  // coeff of urs
            cc(i1,i2,i3, 5) = css;  // coeff of uss
            cc(i1,i2,i3, 6) = crrr;    // coeff of urrr
            cc(i1,i2,i3, 7) = crrs;    // coeff of urrs
            cc(i1,i2,i3, 8) = crss;    // coeff of urss
            cc(i1,i2,i3, 9) = csss;    // coeff of usss
            cc(i1,i2,i3,10) = crrrr;   // coeff of urrrr
            cc(i1,i2,i3,11) = crrrs;   // coeff of urrrs
            cc(i1,i2,i3,12) = crrss;   // coeff of urrss
            cc(i1,i2,i3,13) = crsss;   // coeff of ursss
            cc(i1,i2,i3,14) = cssss;   // coeff of ussss
        } 
    } // end ForStencil
  if( debug & 4 ) 
  {
      printF("WDH:order=4: (i1,i2,i3)=(%3d,%3d,%3d)\n",i1,i2,i3);     
      printF("    coeff4=(%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"    
                    "            %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"    
                    "            %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"    
                    "            %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"    
                    "            %10.3e,%10.3e,%10.3e,%10.3e,%10.3e )\n",  
                      coeff4(-2,-2),coeff4(-1,-2),coeff4(0,-2),coeff4(1,-2),coeff4(2,-2),  
                      coeff4(-2,-1),coeff4(-1,-1),coeff4(0,-1),coeff4(1,-1),coeff4(2,-1),  
                      coeff4(-2, 0),coeff4(-1, 0),coeff4(0, 0),coeff4(1, 0),coeff4(2, 0),  
                      coeff4(-2, 1),coeff4(-1, 1),coeff4(0, 1),coeff4(1, 1),coeff4(2, 1),  
                      coeff4(-2, 2),coeff4(-1, 2),coeff4(0, 2),coeff4(1, 2),coeff4(2, 2)); 
  }
} // end FOR_3D
                            }
                        }
                        else if( axis==1 && axis2==1 )
                        {
                          if( orderOfAccuracy==2 )
                            {
// Coefficients for CHAMP curvilinear order 2, NORMAL DIRECTION=1.
// File written by champGenCode.maple
const Real KLR = theta;
const Real DLR = beta;
printF("WDH: KLR=%9.3e, DLR=%9.3e\n", KLR,DLR);
// Derivatives of the metrics: order 4
// Derivatives of the metrics: order 2
// Fourth-order one-sided approximations
// Here are fully one-sided: *check me* *wdh* April 15, 2022
// Macro : fill end ghost points on b1, b2, by periodicity or extrapolation
                  
// Macro : fill end ghost points on crr, crs, ... by periodicity or extrapolation
                  
const int axis1L=axis,  axis2L=(axis1L+1) % numberOfDimensions;
const int axis1R=axis2, axis2R=(axis1R+1) % numberOfDimensions;
const IntegerArray & gir1 = mg.indexRange();
const IntegerArray & gir2 = mg2.indexRange();
int ksv[3], &ks1=ksv[0], &ks2=ksv[1], &ks3=ksv[2];
int axist;
//  ----- We need b1,b2 at 4 extra tangential ghosts (order=4)---- 
//  since L1 is evaluated at 4 extra ghost 
int extra=2*numGhost;
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray n1L(Ib1,Ib2,Ib3), n2L(Ib1,Ib2,Ib3)  ; // save the normal, it needs to be extended too;
RealArray b1L(Ib1,Ib2,Ib3), b2L(Ib1,Ib2,Ib3);
RealArray b1R(Jb1,Jb2,Jb3), b2R(Jb1,Jb2,Jb3);
// Evaluate the normal derivative coefficients at points on the boundary. Use right-normal = -left-normal.
extra=numGhost;  // we can only directly evaluate 2 ghost using metrics 
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
FOR_3IJD(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    n1L(i1,i2,i3) = normal(i1,i2,i3,0);
    n2L(i1,i2,i3) = normal(i1,i2,i3,1);
    b1L(i1,i2,i3) = -( normal(i1,i2,i3,0)*rxL(i1,i2,i3,0,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,0,1)) ;
    b2L(i1,i2,i3) = -( normal(i1,i2,i3,0)*rxL(i1,i2,i3,1,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,1,1)) ;
    b1R(j1,j2,j3) = -( normal(i1,i2,i3,0)*rxR(j1,j2,j3,0,0) + normal(i1,i2,i3,1)*rxR(j1,j2,j3,0,1)) ;
    b2R(j1,j2,j3) = -( normal(i1,i2,i3,0)*rxR(j1,j2,j3,1,0) + normal(i1,i2,i3,1)*rxR(j1,j2,j3,1,1)) ;
} // end FOR_3IJD

// Fill 3rd and 4th ghost in tangential directions for b1, b2 
const int ghostToFillStart = numGhost+1;
const int ghostToFillEnd   = 2*numGhost;
    axist = (axis + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir1(0,0),gir1(0,0)); 
        I2=Range(gir1(0,1),gir1(0,1)); 
        I3=Range(gir1(0,2),gir1(0,2)); 
        Iv[axis]=Range(gir1(side,axis),gir1(side,axis));
        Iv[axist]=gir1(sidet,axist);
        for( int ghost=ghostToFillStart; ghost<=ghostToFillEnd; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir1(1,0)-gir1(0,0))*ks1;
                Index Ip2 = Ig2 + (gir1(1,1)-gir1(0,1))*ks2;
                Index Ip3 = Ig3 + (gir1(1,2)-gir1(0,2))*ks3;
                n1L(Ig1,Ig2,Ig3) = n1L(Ip1,Ip2,Ip3); 
                n2L(Ig1,Ig2,Ig3) = n2L(Ip1,Ip2,Ip3); 
            }
            else;
            {
                n1L(Ig1,Ig2,Ig3) = (5.*n1L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*n1L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*n1L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*n1L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+n1L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                n2L(Ig1,Ig2,Ig3) = (5.*n2L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*n2L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*n2L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*n2L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+n2L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t
    axist = (axis + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir1(0,0),gir1(0,0)); 
        I2=Range(gir1(0,1),gir1(0,1)); 
        I3=Range(gir1(0,2),gir1(0,2)); 
        Iv[axis]=Range(gir1(side,axis),gir1(side,axis));
        Iv[axist]=gir1(sidet,axist);
        for( int ghost=ghostToFillStart; ghost<=ghostToFillEnd; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir1(1,0)-gir1(0,0))*ks1;
                Index Ip2 = Ig2 + (gir1(1,1)-gir1(0,1))*ks2;
                Index Ip3 = Ig3 + (gir1(1,2)-gir1(0,2))*ks3;
                b1L(Ig1,Ig2,Ig3) = b1L(Ip1,Ip2,Ip3); 
                b2L(Ig1,Ig2,Ig3) = b2L(Ip1,Ip2,Ip3); 
            }
            else;
            {
                b1L(Ig1,Ig2,Ig3) = (5.*b1L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b1L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b1L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b1L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b1L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                b2L(Ig1,Ig2,Ig3) = (5.*b2L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b2L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b2L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b2L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b2L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t
    axist = (axis2 + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir2(0,0),gir2(0,0)); 
        I2=Range(gir2(0,1),gir2(0,1)); 
        I3=Range(gir2(0,2),gir2(0,2)); 
        Iv[axis2]=Range(gir2(side2,axis2),gir2(side2,axis2));
        Iv[axist]=gir2(sidet,axist);
        for( int ghost=ghostToFillStart; ghost<=ghostToFillEnd; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg2.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir2(1,0)-gir2(0,0))*ks1;
                Index Ip2 = Ig2 + (gir2(1,1)-gir2(0,1))*ks2;
                Index Ip3 = Ig3 + (gir2(1,2)-gir2(0,2))*ks3;
                b1R(Ig1,Ig2,Ig3) = b1R(Ip1,Ip2,Ip3); 
                b2R(Ig1,Ig2,Ig3) = b2R(Ip1,Ip2,Ip3); 
            }
            else;
            {
                b1R(Ig1,Ig2,Ig3) = (5.*b1R(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b1R(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b1R(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b1R(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b1R(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                b2R(Ig1,Ig2,Ig3) = (5.*b2R(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b2R(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b2R(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b2R(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b2R(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t

// --- Set the optimized Schwarz parameter: Sl = pl/dx (now that we know b1,b2)  ---- 
//  For now set a constant value on the face 
if( Sl<0 )
{
    j1=Jb1.getBase(); j2=Jb2.getBase(); j3=Jb3.getBase(); // take b1R from this point 
    const Real dnR = fabs( drR[axis2]/b2R(j1,j2,j3) );    // approx. normal grid spacing on right
    Sl = pl/dnR; 
    champParameters(4,side,axis,grid)=Sl; // save Sl
}
extra=0;
getBoundaryIndex(mg.indexRange(),  side,axis, I1,I2,I3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3,extra);
int extraNormal=numGhost-1;
Iv[axis1L]=Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]=Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
RealArray crrL(I1,I2,I3), crsL(I1,I2,I3), cssL(I1,I2,I3), crL(I1,I2,I3), csL(I1,I2,I3);
RealArray crrR(J1,J2,J3), crsR(J1,J2,J3), cssR(J1,J2,J3), crR(J1,J2,J3), csR(J1,J2,J3);
RealArray rxxL(3,3), rxxR(3,3);
RealArray rxyL(3,3), rxyR(3,3);
// We can fill in crr, crs, and ss to 2 ghost since they do not involve terms using rxx, sxx, etc. 
extra=0;
getBoundaryIndex(mg.indexRange(),  side,axis, I1,I2,I3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3,extra);
Iv[axis1L]=Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]=Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
FOR_3IJ(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3)
{
  // Evaluate coefficients of the Laplacian at points near the boundary
  // uxx = (rx*Dr + sx*Ds)*[ (rx)*ur + sx*us ]
  // uxx = (rx)^2 urr + 2*(rx*sx)*urs + (sx)^2 *uss + (rxx)*ur + (sxx)*us
  // uyy = (ry)^2 urr + 2*(ry*sy)*urs + (sy)^2 *uss + (ryy)*ur + (syy)*us
    crrL(i1,i2,i3) = SQR(rxL(i1,i2,i3,0,0)) + SQR(rxL(i1,i2,i3,0,1));                                // (rx)^2 + (ry)^2 
    crsL(i1,i2,i3) = 2.*(rxL(i1,i2,i3,0,0)*rxL(i1,i2,i3,1,0) + rxL(i1,i2,i3,0,1)*rxL(i1,i2,i3,1,1)); // 2( rx*sx + ry*sy ) 
    cssL(i1,i2,i3) = SQR(rxL(i1,i2,i3,1,0)) + SQR(rxL(i1,i2,i3,1,1));                                // (sx)^2 + (sy)^2
  // These are done below
  //crL(i1,i2,i3)  = rxxL(0,0) + rxyL(0,1);  // rxx + ryy
  //csL(i1,i2,i3)  = rxxL(1,0) + rxyL(1,1);  // sxx + syy
  // Evaluate coefficients of the Laplacian at points near the boundary
    crrR(j1,j2,j3) = SQR(rxR(j1,j2,j3,0,0)) + SQR(rxR(j1,j2,j3,0,1));
    crsR(j1,j2,j3) = 2.*( rxR(j1,j2,j3,0,0)*rxR(j1,j2,j3,1,0) + rxR(j1,j2,j3,0,1)*rxR(j1,j2,j3,1,1) );
    cssR(j1,j2,j3) = SQR(rxR(j1,j2,j3,1,0)) + SQR(rxR(j1,j2,j3,1,1));
  // These are done below
  //crR(j1,j2,j3)  = rxxR(0,0) + rxyR(0,1);  // rxx + ryy
  //csR(j1,j2,j3)  = rxxR(1,0) + rxyR(1,1);  // sxx + syy
  // printF("WDH: crrL=%9.3e, crsL=%9.3e, cssL=%9.3e\n",crrL(i1,i2,i3),crsL(i1,i2,i3),cssL(i1,i2,i3) );
  // printF("WDH: crrR=%9.3e, crsR=%9.3e, cssR=%9.3e\n",crrR(j1,j2,j3),crsR(j1,j2,j3),cssR(j1,j2,j3) );
} // end FOR_3IJ

// We can only compute cr and cs up to the boundary since we need two ghost points to evaluate rxx, rxy, sxx, sxy
getBoundaryIndex(mg.indexRange(),side,axis,I1,I2,I3);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3);
Iv[axis1L]= Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]= Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
Real rxr,rxs;
FOR_3IJ(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3)
{
    for( int m2=0; m2<numberOfDimensions; m2++ )
    {
        for( int m1=0; m1<numberOfDimensions; m1++ )
        {
        // Evaluate derivatives of metrics to order 2 .. we are computing more than we need here *fix me*
                rxxL(m1,m2) = rxL(i1,i2,i3,0,0)*(-rxL(i1-1,i2,i3,m1,m2)+rxL(i1+1,i2,i3,m1,m2))/(2.*drL[0])+rxL(i1,i2,i3,1,0)*(-rxL(i1,i2-1,i3,m1,m2)+rxL(i1,i2+1,i3,m1,m2))/(2.*drL[1]);
                rxyL(m1,m2) = rxL(i1,i2,i3,0,1)*(-rxL(i1-1,i2,i3,m1,m2)+rxL(i1+1,i2,i3,m1,m2))/(2.*drL[0])+rxL(i1,i2,i3,1,1)*(-rxL(i1,i2-1,i3,m1,m2)+rxL(i1,i2+1,i3,m1,m2))/(2.*drL[1]);
                                                                                                            
                rxxR(m1,m2) = rxR(j1,j2,j3,0,0)*(-rxR(j1-1,j2,j3,m1,m2)+rxR(j1+1,j2,j3,m1,m2))/(2.*drR[0])+rxR(j1,j2,j3,1,0)*(-rxR(j1,j2-1,j3,m1,m2)+rxR(j1,j2+1,j3,m1,m2))/(2.*drR[1]);
                rxyR(m1,m2) = rxR(j1,j2,j3,0,1)*(-rxR(j1-1,j2,j3,m1,m2)+rxR(j1+1,j2,j3,m1,m2))/(2.*drR[0])+rxR(j1,j2,j3,1,1)*(-rxR(j1,j2-1,j3,m1,m2)+rxR(j1,j2+1,j3,m1,m2))/(2.*drR[1]);
        }
    }
  // Evaluate coefficients of the Laplacian at points near the boundary
  // uxx = (rx*Dr + sx*Ds)*[ (rx)*ur + sx*us ]
  // uxx = (rx)^2 urr + 2*(rx*sx)*urs + (sx)^2 *uss + (rxx)*ur + (sxx)*us
  // uyy = (ry)^2 urr + 2*(ry*sy)*urs + (sy)^2 *uss + (ryy)*ur + (syy)*us
    crL(i1,i2,i3)  = rxxL(0,0) + rxyL(0,1);  // rxx + ryy
    csL(i1,i2,i3)  = rxxL(1,0) + rxyL(1,1);  // sxx + syy
  // Evaluate coefficients of the Laplacian at points near the boundary
    crR(j1,j2,j3)  = rxxR(0,0) + rxyR(0,1);  // rxx + ryy
    csR(j1,j2,j3)  = rxxR(1,0) + rxyR(1,1);  // sxx + syy
    if( advectionIsOn )  // include advection terms scaled by 1/D : (u/D).grad 
    {
            if( variableAdvection )
            {
                u1DL = advectVarL(i1,i2,i3,0)/DL;
                u2DL = advectVarL(i1,i2,i3,1)/DL;
                u1DR = advectVarR(j1,j2,j3,0)/DR;
                u2DR = advectVarR(j1,j2,j3,1)/DR;
            }
        crL(i1,i2,i3) = crL(i1,i2,i3) - ( u1DL*rxL(i1,i2,i3,0,0) + u2DL*rxL(i1,i2,i3,0,1) ); 
        csL(i1,i2,i3) = csL(i1,i2,i3) - ( u1DL*rxL(i1,i2,i3,1,0) + u2DL*rxL(i1,i2,i3,1,1) ); 
        crR(j1,j2,j3) = crR(j1,j2,j3) - ( u1DR*rxR(j1,j2,j3,0,0) + u2DR*rxR(j1,j2,j3,0,1) ); 
        csR(j1,j2,j3) = csR(j1,j2,j3) - ( u1DR*rxR(j1,j2,j3,1,0) + u2DR*rxR(j1,j2,j3,1,1) ); 
    }
  // printF("WDH: (i1,i2)=(%4d,%4d) crL=%9.3e csL=%9.3e\n",i1,i2,crL(i1,i2,i3),csL(i1,i2,i3) );
  // printF("WDH: (j1,j2)=(%4d,%4d) crR=%9.3e csR=%9.3e\n",j1,j2,crR(j1,j2,j3),csR(j1,j2,j3) );
} // end FOR_3IJ


// ----- Evaluate coefficients of L1 on the boundary -----

// 2*numGhost = 2 extra ghost needed for 2nd order, true ? 
extra=2*numGhost;
getBoundaryIndex(mg.indexRange(), side ,axis ,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L1ICoeff(Ib1,Ib2,Ib3), L1ICoeffr(Ib1,Ib2,Ib3);
RealArray L1rCoeff(Ib1,Ib2,Ib3), L1rCoeffr(Ib1,Ib2,Ib3);
RealArray L1sCoeff(Ib1,Ib2,Ib3), L1sCoeffr(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    const Real n1R = -n1L(i1,i2,i3), n2R = -n2L(i1,i2,i3);    // for advection, use right normal 
        if( variableAdvection )
        {
            u1DL = advectVarL(i1,i2,i3,0)/DL;
            u2DL = advectVarL(i1,i2,i3,1)/DL;
            u1DR = advectVarR(j1,j2,j3,0)/DR;
            u2DR = advectVarR(j1,j2,j3,1)/DR;
        }
    L1ICoeff(i1,i2,i3) = ((-1.*KLR*u2DL+1.*u2DR)*n2R+(-1.*KLR*u1DL+1.*u1DR)*n1R)/b2R(j1,j2,j3);
    L1rCoeff(i1,i2,i3) = (1.*KLR*b1L(i1,i2,i3)-1.*b1R(j1,j2,j3))/b2R(j1,j2,j3);
    L1sCoeff(i1,i2,i3) = 1.*b2L(i1,i2,i3)*KLR/b2R(j1,j2,j3);
} // end FOR_3D

// ----- Evaluate coefficients of L2 on the boundary -----
extra=0;
getBoundaryIndex(mg.indexRange(), side,  axis,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L2rrCoeff(Ib1,Ib2,Ib3);
RealArray L2rsCoeff(Ib1,Ib2,Ib3);
RealArray L2ssCoeff(Ib1,Ib2,Ib3);
RealArray L2rCoeff(Ib1,Ib2,Ib3);
RealArray L2sCoeff(Ib1,Ib2,Ib3);
RealArray L2ICoeff(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  // Compute tangential derivatives of L1 - label tangential directions a and b (for 3D)
    L1rCoeffr(i1,i2,i3) = ((L1rCoeff(i1+1,i2,i3)-L1rCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L1sCoeffr(i1,i2,i3) = ((L1sCoeff(i1+1,i2,i3)-L1sCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L1ICoeffr(i1,i2,i3) = ((L1ICoeff(i1+1,i2,i3)-L1ICoeff(i1-1,i2,i3))/(2.*drL[0]));
  // Compute L2 coefficients
    L2rrCoeff(i1,i2,i3) = 1.*(DLR*crrL(i1,i2,i3)-crsR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-crrR(j1,j2,j3))/cssR(j1,j2,j3);
    L2rsCoeff(i1,i2,i3) = (1.*DLR*crsL(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1sCoeff(i1,i2,i3))/cssR(j1,j2,j3);
    L2ssCoeff(i1,i2,i3) = 1.*cssL(i1,i2,i3)*DLR/cssR(j1,j2,j3);
    L2rCoeff(i1,i2,i3) = ((-1.*L1ICoeff(i1,i2,i3)-1.*L1rCoeffr(i1,i2,i3))*crsR(j1,j2,j3)+1.*DLR*crL(i1,i2,i3)-1.*csR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*crR(j1,j2,j3))/cssR(j1,j2,j3);
    L2sCoeff(i1,i2,i3) = 1.*(DLR*csL(i1,i2,i3)-crsR(j1,j2,j3)*L1sCoeffr(i1,i2,i3)-csR(j1,j2,j3)*L1sCoeff(i1,i2,i3))/cssR(j1,j2,j3);
    L2ICoeff(i1,i2,i3) = (-1.*crsR(j1,j2,j3)*L1ICoeffr(i1,i2,i3)-1.*csR(j1,j2,j3)*L1ICoeff(i1,i2,i3))/cssR(j1,j2,j3);
} // end FOR_3D

// ----- Define coefficients in difference operators -----
Range R4(-1,1);
RealArray iCoeff(R4);
RealArray rCoeff(R4), rrCoeff(R4), rrrCoeff(R4), rrrrCoeff(R4);
RealArray sCoeff(R4), ssCoeff(R4), sssCoeff(R4), ssssCoeff(R4);
      iCoeff(-1)=  0.;    iCoeff(0)=  1.;    iCoeff(1)= 0.;   
      rCoeff(-1)= -1.;    rCoeff(0)=  0.;    rCoeff(1)= 1.;    rCoeff    /=(2.*drL[0])  ;
    rrCoeff(-1)=  1.;   rrCoeff(0)= -2.;   rrCoeff(1)= 1.;    rrCoeff   /=(SQR(drL[0]));
      sCoeff(-1)= -1.;    sCoeff(0)=  0.;    sCoeff(1)= 1.;    sCoeff    /=(2.*drL[1])  ;
    ssCoeff(-1)=  1.;   ssCoeff(0)= -2.;   ssCoeff(1)= 1.;    ssCoeff   /=(SQR(drL[1]));

// ----- Fill in the Matrix Coefficients for CHAMP -----
const Real dxs = (1-2*side2)*drR[axis2];
const Real h = dxs;
printF("WDH: grid=%d, (side,axis)=(%d,%d) (side2,axis2)=(%d,%d) dxs=%9.3e, Sl=%9.3e\n",grid, side,axis, side2,axis2, dxs,Sl);
RealArray coeff4(R4,R4);
const int e=0, c=0; // eqn number and component number
// ------- FILL CHAMP conditions into the matrix -----
// NOTE: skip top boundary if periodic (use indexRange) 
// NOTE: skip adjacent boundaries if Dirichlet BC *finish me**
extra=0; 
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
int axisp = (axis+1) % numberOfDimensions;
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    i1m=i1-is1, i2m=i2-is2, i3m=i3-is3; //  ghost point is (i1m,i2m,i3m)
    const Real bn = b2R(j1,j2,j3), bt = b1R(j1,j2,j3);  // b in normal and tangential directions
  // The next macro defines cI, cr, cs, crr, ...
  // ------  ORDER = 2 , NormalDirection = 1,------
    const Real h2By2=h*h/2., h3By6=h2By2*h/3., h4By24=h3By6*h/4.;
    Real L1IC = L1ICoeff(i1,i2,i3);
    Real L1ICr = L1ICoeffr(i1,i2,i3);
    Real L1rC = L1rCoeff(i1,i2,i3);
    Real L1rCr = L1rCoeffr(i1,i2,i3);
    Real L1sC = L1sCoeff(i1,i2,i3);
    Real L1sCr = L1sCoeffr(i1,i2,i3);
    Real L2IC = L2ICoeff(i1,i2,i3);
    Real L2rC = L2rCoeff(i1,i2,i3);
    Real L2sC = L2sCoeff(i1,i2,i3);
    Real L2rrC = L2rrCoeff(i1,i2,i3);
    Real L2rsC = L2rsCoeff(i1,i2,i3);
    Real L2ssC = L2ssCoeff(i1,i2,i3);
    Real cI  = (Sl*L1IC-bn*L2IC-bt*L1ICr)*h+Sl*h2By2*L2IC-bn*L1IC+Sl;
    Real cr  = ((-L1rCr-L1IC)*h-L2IC*h2By2-2)*bt+(Sl*L1rC-bn*L2rC)*h+Sl*h2By2*L2rC-bn*L1rC;
    Real cs  = (Sl*L1sC-bn*L2sC-bt*L1sCr)*h+Sl*h2By2*L2sC-bn*L1sC;
    Real crr  = (-bn*L2rrC-bt*L1rC)*h+(Sl*L2rrC-bt*L2rC)*h2By2;
    Real crs  = (-bn*L2rsC-bt*L1sC)*h+(Sl*L2rsC-bt*L2sC)*h2By2;
    Real css  = L2ssC*(Sl*h2By2-bn*h);
//  if( orderOfAccuracy==4 )
//  {
//    printF("domain2=%d: cI   =%12.4e, cr   =%12.4e, cs   =%12.4e, crr  =%12.4e, crs=%12.4e, css=%12.4e\n",domain2,cI,cr,cs,crr,crs,css);
//    printF("            crrr =%12.4e, crrs =%12.4e, crss =%12.4e, csss =%12.4e\n",crrr,crrs,crss,csss);
//    printF("            crrrr=%12.4e, crrrs=%12.4e, crrss=%12.4e, crsss=%12.4e, cssss=%12.4e\n",crrrr,crrrs,crrss,crsss,cssss);
//    printF("            L4ssssC=%12.4e, (Sl*h4By24-bn*h3By6)=%12.4e h4By24=%12.4e h3By6=%12.4e, Sl=%12.4e, bn=%12.4e, h=%12.4e, bt=%12.4e\n",L4ssssC,(Sl*h4By24-bn*h3By6),h4By24,h3By6,Sl,bn,h,bt);
//    if( domain2==0 )
//    {
//      OV_ABORT("stop here for now");
//    }
//  }
    ForStencil(m1,m2,m3)
    {
        int m  = M123(m1,m2,m3);        // the single-component coeff-index
        int mm = M123CE(m1,m2,m3,c,e);  // the system coeff-index
        coeff4(m1,m2) = 
                          + cI   *   iCoeff(m1) * iCoeff(m2) 
                          + cr   *   rCoeff(m1) * iCoeff(m2) 
                          + cs   *   iCoeff(m1) * sCoeff(m2) 
                          + crr  *  rrCoeff(m1) * iCoeff(m2) 
                          + crs  *   rCoeff(m1) * sCoeff(m2) 
                          + css  *   iCoeff(m1) *ssCoeff(m2) 
                          ;
        if( fillMatrixWDH )
        {
            coeff(mm,i1m,i2m,i3m) = coeff4(m1,m2);
      // Specify that the above coeff value is the coefficient of component c at the grid point (j1,j2,j3).
            const int k1=i1+m1, k2=i2+m2, k3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)    
            setEquationNumber(mm, e,i1m,i2m,i3m,  c,k1,k2,k3 );  // macro to set equationNumber                  
                                                                                                                                                                                                                      
       // Fill Ghost 2 -- extrapolation for now                                                            
              const int ghost=2;                                                                                  
                  const int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
         // --- fill in the coefficients of the extrapolation formula ---
                  for( int me=0; me<=extrapOrder; me++ )
                  {
                      coeff(me,i1m,i2m,i3m) = extrapCoeff[me];
                      const int j1=i1m + me*is1, j2=i2m + me*is2, j3=i3m + me*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                      setEquationNumber(me, e,i1m,i2m,i3m,  c,j1,j2,j3 );             // macro to set equationNumber
                  }                
        }
        if( twilightZoneFlow )
        {
      // --- For twilightZone we save some coefficients that go into the CHAMP matrix ---
            cc(i1,i2,i3, 0) = cI;   // coeff of I
            cc(i1,i2,i3, 1) = cr;   // coeff of ur
            cc(i1,i2,i3, 2) = cs;   // coeff of us
            cc(i1,i2,i3, 3) = crr;  // coeff of urr
            cc(i1,i2,i3, 4) = crs;  // coeff of urs
            cc(i1,i2,i3, 5) = css;  // coeff of uss
        } 
    } // end ForStencil
  if( debug & 4 ) 
  {
      printF("WDH:order=2: (i1,i2,i3)=(%3d,%3d,%3d)\n",i1,i2,i3);     
      printF("    coeff =(%10.3e,%10.3e,%10.3e,\n"    
                    "            %10.3e,%10.3e,%10.3e,\n"    
                    "            %10.3e,%10.3e,%10.3e )\n",  
                      coeff4(-1,-1),coeff4(0,-1),coeff4(1,-1),  
                      coeff4(-1, 0),coeff4(0, 0),coeff4(1, 0),  
                      coeff4(-1, 1),coeff4(0, 1),coeff4(1, 1)); 
  }
} // end FOR_3D
                            }
                            else
                            {
// Coefficients for CHAMP curvilinear order 4, NORMAL DIRECTION=1.
// File written by champGenCode.maple
const Real KLR = theta;
const Real DLR = beta;
printF("WDH: KLR=%9.3e, DLR=%9.3e\n", KLR,DLR);
// Derivatives of the metrics: order 4
// Derivatives of the metrics: order 2
// Fourth-order one-sided approximations
// Here are fully one-sided: *check me* *wdh* April 15, 2022
// Macro : fill end ghost points on b1, b2, by periodicity or extrapolation
                  
// Macro : fill end ghost points on crr, crs, ... by periodicity or extrapolation
                  
const int axis1L=axis,  axis2L=(axis1L+1) % numberOfDimensions;
const int axis1R=axis2, axis2R=(axis1R+1) % numberOfDimensions;
const IntegerArray & gir1 = mg.indexRange();
const IntegerArray & gir2 = mg2.indexRange();
int ksv[3], &ks1=ksv[0], &ks2=ksv[1], &ks3=ksv[2];
int axist;
//  ----- We need b1,b2 at 4 extra tangential ghosts (order=4)---- 
//  since L1 is evaluated at 4 extra ghost 
int extra=2*numGhost;
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray n1L(Ib1,Ib2,Ib3), n2L(Ib1,Ib2,Ib3)  ; // save the normal, it needs to be extended too;
RealArray b1L(Ib1,Ib2,Ib3), b2L(Ib1,Ib2,Ib3);
RealArray b1R(Jb1,Jb2,Jb3), b2R(Jb1,Jb2,Jb3);
// Evaluate the normal derivative coefficients at points on the boundary. Use right-normal = -left-normal.
extra=numGhost;  // we can only directly evaluate 2 ghost using metrics 
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
FOR_3IJD(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    n1L(i1,i2,i3) = normal(i1,i2,i3,0);
    n2L(i1,i2,i3) = normal(i1,i2,i3,1);
    b1L(i1,i2,i3) = -( normal(i1,i2,i3,0)*rxL(i1,i2,i3,0,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,0,1)) ;
    b2L(i1,i2,i3) = -( normal(i1,i2,i3,0)*rxL(i1,i2,i3,1,0) + normal(i1,i2,i3,1)*rxL(i1,i2,i3,1,1)) ;
    b1R(j1,j2,j3) = -( normal(i1,i2,i3,0)*rxR(j1,j2,j3,0,0) + normal(i1,i2,i3,1)*rxR(j1,j2,j3,0,1)) ;
    b2R(j1,j2,j3) = -( normal(i1,i2,i3,0)*rxR(j1,j2,j3,1,0) + normal(i1,i2,i3,1)*rxR(j1,j2,j3,1,1)) ;
} // end FOR_3IJD

// Fill 2nd ghost in tangential directions for b1, b2, and nL 
const int ghostToFillStart = numGhost+1;
const int ghostToFillEnd   = 2*numGhost;
    axist = (axis + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir1(0,0),gir1(0,0)); 
        I2=Range(gir1(0,1),gir1(0,1)); 
        I3=Range(gir1(0,2),gir1(0,2)); 
        Iv[axis]=Range(gir1(side,axis),gir1(side,axis));
        Iv[axist]=gir1(sidet,axist);
        for( int ghost=ghostToFillStart; ghost<=ghostToFillEnd; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir1(1,0)-gir1(0,0))*ks1;
                Index Ip2 = Ig2 + (gir1(1,1)-gir1(0,1))*ks2;
                Index Ip3 = Ig3 + (gir1(1,2)-gir1(0,2))*ks3;
                n1L(Ig1,Ig2,Ig3) = n1L(Ip1,Ip2,Ip3); 
                n2L(Ig1,Ig2,Ig3) = n2L(Ip1,Ip2,Ip3); 
            }
            else;
            {
                n1L(Ig1,Ig2,Ig3) = (5.*n1L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*n1L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*n1L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*n1L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+n1L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                n2L(Ig1,Ig2,Ig3) = (5.*n2L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*n2L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*n2L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*n2L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+n2L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t
    axist = (axis + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir1(0,0),gir1(0,0)); 
        I2=Range(gir1(0,1),gir1(0,1)); 
        I3=Range(gir1(0,2),gir1(0,2)); 
        Iv[axis]=Range(gir1(side,axis),gir1(side,axis));
        Iv[axist]=gir1(sidet,axist);
        for( int ghost=ghostToFillStart; ghost<=ghostToFillEnd; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir1(1,0)-gir1(0,0))*ks1;
                Index Ip2 = Ig2 + (gir1(1,1)-gir1(0,1))*ks2;
                Index Ip3 = Ig3 + (gir1(1,2)-gir1(0,2))*ks3;
                b1L(Ig1,Ig2,Ig3) = b1L(Ip1,Ip2,Ip3); 
                b2L(Ig1,Ig2,Ig3) = b2L(Ip1,Ip2,Ip3); 
            }
            else;
            {
                b1L(Ig1,Ig2,Ig3) = (5.*b1L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b1L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b1L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b1L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b1L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                b2L(Ig1,Ig2,Ig3) = (5.*b2L(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b2L(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b2L(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b2L(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b2L(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t
    axist = (axis2 + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir2(0,0),gir2(0,0)); 
        I2=Range(gir2(0,1),gir2(0,1)); 
        I3=Range(gir2(0,2),gir2(0,2)); 
        Iv[axis2]=Range(gir2(side2,axis2),gir2(side2,axis2));
        Iv[axist]=gir2(sidet,axist);
        for( int ghost=ghostToFillStart; ghost<=ghostToFillEnd; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg2.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir2(1,0)-gir2(0,0))*ks1;
                Index Ip2 = Ig2 + (gir2(1,1)-gir2(0,1))*ks2;
                Index Ip3 = Ig3 + (gir2(1,2)-gir2(0,2))*ks3;
                b1R(Ig1,Ig2,Ig3) = b1R(Ip1,Ip2,Ip3); 
                b2R(Ig1,Ig2,Ig3) = b2R(Ip1,Ip2,Ip3); 
            }
            else;
            {
                b1R(Ig1,Ig2,Ig3) = (5.*b1R(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b1R(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b1R(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b1R(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b1R(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                b2R(Ig1,Ig2,Ig3) = (5.*b2R(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*b2R(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*b2R(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*b2R(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+b2R(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t

// --- Set the optimized Schwarz parameter: Sl = pl/dx (now that we know b1,b2)  ---- 
//  For now set a constant value on the face 
if( Sl<0 )
{
    j1=Jb1.getBase(); j2=Jb2.getBase(); j3=Jb3.getBase(); // take b1R from this point 
    const Real dnR = fabs( drR[axis2]/b2R(j1,j2,j3) );    // approx. normal grid spacing on right
    Sl = pl/dnR; 
    champParameters(4,side,axis,grid)=Sl; // save Sl
}
// *NOTE* WE NEED AN EXTRA GHOST IN TANGENTIAL DIRECTION FOR CURVILINEAR 
// Index bounds for points near the boundary: add 3 extra points in tangential directions and 1 in normal direction.
extra=numGhost+1;
getBoundaryIndex(mg.indexRange(),  side,axis, I1,I2,I3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3,extra);
int extraNormal=numGhost-1;
Iv[axis1L]=Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]=Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
RealArray crrL(I1,I2,I3), crsL(I1,I2,I3), cssL(I1,I2,I3), crL(I1,I2,I3), csL(I1,I2,I3);
RealArray crrR(J1,J2,J3), crsR(J1,J2,J3), cssR(J1,J2,J3), crR(J1,J2,J3), csR(J1,J2,J3);
RealArray rxxL(3,3), rxxR(3,3);
RealArray rxyL(3,3), rxyR(3,3);
// We can fill in crr, crs, and ss to 2 ghost since they do not involve terms using rxx, sxx, etc. 
extra=numGhost;
getBoundaryIndex(mg.indexRange(),  side,axis, I1,I2,I3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3,extra);
Iv[axis1L]=Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]=Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
FOR_3IJ(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3)
{
  // Evaluate coefficients of the Laplacian at points near the boundary
  // uxx = (rx*Dr + sx*Ds)*[ (rx)*ur + sx*us ]
  // uxx = (rx)^2 urr + 2*(rx*sx)*urs + (sx)^2 *uss + (rxx)*ur + (sxx)*us
  // uyy = (ry)^2 urr + 2*(ry*sy)*urs + (sy)^2 *uss + (ryy)*ur + (syy)*us
    crrL(i1,i2,i3) = SQR(rxL(i1,i2,i3,0,0)) + SQR(rxL(i1,i2,i3,0,1));                                // (rx)^2 + (ry)^2 
    crsL(i1,i2,i3) = 2.*(rxL(i1,i2,i3,0,0)*rxL(i1,i2,i3,1,0) + rxL(i1,i2,i3,0,1)*rxL(i1,i2,i3,1,1)); // 2( rx*sx + ry*sy ) 
    cssL(i1,i2,i3) = SQR(rxL(i1,i2,i3,1,0)) + SQR(rxL(i1,i2,i3,1,1));                                // (sx)^2 + (sy)^2
  // These are done below
  //crL(i1,i2,i3)  = rxxL(0,0) + rxyL(0,1);  // rxx + ryy
  //csL(i1,i2,i3)  = rxxL(1,0) + rxyL(1,1);  // sxx + syy
  // Evaluate coefficients of the Laplacian at points near the boundary
    crrR(j1,j2,j3) = SQR(rxR(j1,j2,j3,0,0)) + SQR(rxR(j1,j2,j3,0,1));
    crsR(j1,j2,j3) = 2.*( rxR(j1,j2,j3,0,0)*rxR(j1,j2,j3,1,0) + rxR(j1,j2,j3,0,1)*rxR(j1,j2,j3,1,1) );
    cssR(j1,j2,j3) = SQR(rxR(j1,j2,j3,1,0)) + SQR(rxR(j1,j2,j3,1,1));
  // These are done below
  //crR(j1,j2,j3)  = rxxR(0,0) + rxyR(0,1);  // rxx + ryy
  //csR(j1,j2,j3)  = rxxR(1,0) + rxyR(1,1);  // sxx + syy
  // printF("WDH: crrL=%9.3e, crsL=%9.3e, cssL=%9.3e\n",crrL(i1,i2,i3),crsL(i1,i2,i3),cssL(i1,i2,i3) );
  // printF("WDH: crrR=%9.3e, crsR=%9.3e, cssR=%9.3e\n",crrR(j1,j2,j3),crsR(j1,j2,j3),cssR(j1,j2,j3) );
} // end FOR_3IJ

// We can only compute cr and cs up to the boundary since we need two ghost points to evaluate rxx, rxy, sxx, sxy
getBoundaryIndex(mg.indexRange(),side,axis,I1,I2,I3);
getBoundaryIndex(mg2.indexRange(),side2,axis2,J1,J2,J3);
Iv[axis1L]= Range( Iv[axis1L].getBase()-extraNormal,Iv[axis1L].getBound()+extraNormal );
Jv[axis1R]= Range( Jv[axis1R].getBase()-extraNormal,Jv[axis1R].getBound()+extraNormal );
Real rxr,rxs;
FOR_3IJ(i1,i2,i3,I1,I2,I3,j1,j2,j3,J1,J2,J3)
{
    for( int m2=0; m2<numberOfDimensions; m2++ )
    {
        for( int m1=0; m1<numberOfDimensions; m1++ )
        {
        // Evaluate derivatives of metrics to order 4 .. we are computing more than we need here *fix me*
                if( i1-2 >= mg.dimension(0,0) && i1+2 <= mg.dimension(1,0) )
                    rxr = (((1./12.)*rxL(i1-2,i2,i3,m1,m2)-(2./3.)*rxL(i1-1,i2,i3,m1,m2)+(2./3.)*rxL(i1+1,i2,i3,m1,m2)-(1./12.)*rxL(i1+2,i2,i3,m1,m2))/(drL[0]));              // centered 
                else if( i1-2 < mg.dimension(0,0) )
                    rxr = (((-25./12.)*rxL(i1,i2,i3,m1,m2)+(4.)*rxL(i1+1,i2,i3,m1,m2)-(3.)*rxL(i1+2,i2,i3,m1,m2)+(4./3.)*rxL(i1+3,i2,i3,m1,m2)-(1./4.)*rxL(i1+4,i2,i3,m1,m2))/(drL[0]));   // one-sided to right
                else 
                    rxr = (((+25./12.)*rxL(i1,i2,i3,m1,m2)-(4.)*rxL(i1-1,i2,i3,m1,m2)+(3.)*rxL(i1-2,i2,i3,m1,m2)-(4./3.)*rxL(i1-3,i2,i3,m1,m2)+(1./4.)*rxL(i1-4,i2,i3,m1,m2))/(drL[0]));  // one-sided to left
                if( i2-2 >= mg.dimension(0,1) && i2+2 <= mg.dimension(1,1) )
                    rxs = (((1./12.)*rxL(i1,i2-2,i3,m1,m2)-(2./3.)*rxL(i1,i2-1,i3,m1,m2)+(2./3.)*rxL(i1,i2+1,i3,m1,m2)-(1./12.)*rxL(i1,i2+2,i3,m1,m2))/(drL[1]));              // centered 
                else if( i2-2 < mg.dimension(0,1) )
                    rxs = (((-25./12.)*rxL(i1,i2,i3,m1,m2)+(4.)*rxL(i1,i2+1,i3,m1,m2)-(3.)*rxL(i1,i2+2,i3,m1,m2)+(4./3.)*rxL(i1,i2+3,i3,m1,m2)-(1./4.)*rxL(i1,i2+4,i3,m1,m2))/(drL[1]));   // one-sided to right
                else 
                    rxs = (((+25./12.)*rxL(i1,i2,i3,m1,m2)-(4.)*rxL(i1,i2-1,i3,m1,m2)+(3.)*rxL(i1,i2-2,i3,m1,m2)-(4./3.)*rxL(i1,i2-3,i3,m1,m2)+(1./4.)*rxL(i1,i2-4,i3,m1,m2))/(drL[1]));  // one-sided to left
                rxxL(m1,m2) = rxL(i1,i2,i3,0,0)*rxr + rxL(i1,i2,i3,1,0)*rxs;
                rxyL(m1,m2) = rxL(i1,i2,i3,0,1)*rxr + rxL(i1,i2,i3,1,1)*rxs;
                                                                                                            
                if( j1-2 >= mg2.dimension(0,0) && j1+2 <= mg2.dimension(1,0) )
                    rxr = (((1./12.)*rxR(j1-2,j2,j3,m1,m2)-(2./3.)*rxR(j1-1,j2,j3,m1,m2)+(2./3.)*rxR(j1+1,j2,j3,m1,m2)-(1./12.)*rxR(j1+2,j2,j3,m1,m2))/(drR[0]));  // centered 
                else if( j1-2 < mg2.dimension(0,0) )
                    rxr = (((-25./12.)*rxR(j1,j2,j3,m1,m2)+(4.)*rxR(j1+1,j2,j3,m1,m2)-(3.)*rxR(j1+2,j2,j3,m1,m2)+(4./3.)*rxR(j1+3,j2,j3,m1,m2)-(1./4.)*rxR(j1+4,j2,j3,m1,m2))/(drR[0]));  // one-sided to right
                else 
                    rxr = (((+25./12.)*rxR(j1,j2,j3,m1,m2)-(4.)*rxR(j1-1,j2,j3,m1,m2)+(3.)*rxR(j1-2,j2,j3,m1,m2)-(4./3.)*rxR(j1-3,j2,j3,m1,m2)+(1./4.)*rxR(j1-4,j2,j3,m1,m2))/(drR[0]));  // one-sided to left
                if( j2-2 >= mg2.dimension(0,1) && j2+2 <= mg2.dimension(1,1) )
                    rxs = (((1./12.)*rxR(j1,j2-2,j3,m1,m2)-(2./3.)*rxR(j1,j2-1,j3,m1,m2)+(2./3.)*rxR(j1,j2+1,j3,m1,m2)-(1./12.)*rxR(j1,j2+2,j3,m1,m2))/(drR[1]));  // centered 
                else if( j2-2 < mg2.dimension(0,1) )
                    rxs = (((-25./12.)*rxR(j1,j2,j3,m1,m2)+(4.)*rxR(j1,j2+1,j3,m1,m2)-(3.)*rxR(j1,j2+2,j3,m1,m2)+(4./3.)*rxR(j1,j2+3,j3,m1,m2)-(1./4.)*rxR(j1,j2+4,j3,m1,m2))/(drR[1]));  // one-sided to right
                else 
                    rxs = (((+25./12.)*rxR(j1,j2,j3,m1,m2)-(4.)*rxR(j1,j2-1,j3,m1,m2)+(3.)*rxR(j1,j2-2,j3,m1,m2)-(4./3.)*rxR(j1,j2-3,j3,m1,m2)+(1./4.)*rxR(j1,j2-4,j3,m1,m2))/(drR[1]));  // one-sided to left
                rxxR(m1,m2) = rxR(j1,j2,j3,0,0)*rxr + rxR(j1,j2,j3,1,0)*rxs;
                rxyR(m1,m2) = rxR(j1,j2,j3,0,1)*rxr + rxR(j1,j2,j3,1,1)*rxs;
                                                                                                                                        
        }
    }
  // Evaluate coefficients of the Laplacian at points near the boundary
  // uxx = (rx*Dr + sx*Ds)*[ (rx)*ur + sx*us ]
  // uxx = (rx)^2 urr + 2*(rx*sx)*urs + (sx)^2 *uss + (rxx)*ur + (sxx)*us
  // uyy = (ry)^2 urr + 2*(ry*sy)*urs + (sy)^2 *uss + (ryy)*ur + (syy)*us
    crL(i1,i2,i3)  = rxxL(0,0) + rxyL(0,1);  // rxx + ryy
    csL(i1,i2,i3)  = rxxL(1,0) + rxyL(1,1);  // sxx + syy
  // Evaluate coefficients of the Laplacian at points near the boundary
    crR(j1,j2,j3)  = rxxR(0,0) + rxyR(0,1);  // rxx + ryy
    csR(j1,j2,j3)  = rxxR(1,0) + rxyR(1,1);  // sxx + syy
    if( advectionIsOn )  // include advection terms scaled by 1/D : (u/D).grad 
    {
            if( variableAdvection )
            {
                u1DL = advectVarL(i1,i2,i3,0)/DL;
                u2DL = advectVarL(i1,i2,i3,1)/DL;
                u1DR = advectVarR(j1,j2,j3,0)/DR;
                u2DR = advectVarR(j1,j2,j3,1)/DR;
            }
        crL(i1,i2,i3) = crL(i1,i2,i3) - ( u1DL*rxL(i1,i2,i3,0,0) + u2DL*rxL(i1,i2,i3,0,1) ); 
        csL(i1,i2,i3) = csL(i1,i2,i3) - ( u1DL*rxL(i1,i2,i3,1,0) + u2DL*rxL(i1,i2,i3,1,1) ); 
        crR(j1,j2,j3) = crR(j1,j2,j3) - ( u1DR*rxR(j1,j2,j3,0,0) + u2DR*rxR(j1,j2,j3,0,1) ); 
        csR(j1,j2,j3) = csR(j1,j2,j3) - ( u1DR*rxR(j1,j2,j3,1,0) + u2DR*rxR(j1,j2,j3,1,1) ); 
    }
  // printF("WDH: (i1,i2)=(%4d,%4d) crL=%9.3e csL=%9.3e\n",i1,i2,crL(i1,i2,i3),csL(i1,i2,i3) );
  // printF("WDH: (j1,j2)=(%4d,%4d) crR=%9.3e csR=%9.3e\n",j1,j2,crR(j1,j2,j3),csR(j1,j2,j3) );
} // end FOR_3IJ

// --- Now fill in one extra ghost for crr,crs,css and two extra ghost for cr,cs --- 
// --- Fill ghost points on ends through extrapolation or periodicity --- 
// Fill ghost in tangential directions for crr, crs, etc.
    axist = (axis + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir1(0,0),gir1(0,0)); 
        I2=Range(gir1(0,1),gir1(0,1)); 
        I3=Range(gir1(0,2),gir1(0,2)); 
        Iv[axis]=Range(gir1(side,axis)-1,gir1(side,axis)+1);
        Iv[axist]=gir1(sidet,axist);
        for( int ghost=1; ghost<=numGhost+1; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir1(1,0)-gir1(0,0))*ks1;
                Index Ip2 = Ig2 + (gir1(1,1)-gir1(0,1))*ks2;
                Index Ip3 = Ig3 + (gir1(1,2)-gir1(0,2))*ks3;
                crrL(Ig1,Ig2,Ig3) = crrL(Ip1,Ip2,Ip3); 
                crsL(Ig1,Ig2,Ig3) = crsL(Ip1,Ip2,Ip3); 
                cssL(Ig1,Ig2,Ig3) = cssL(Ip1,Ip2,Ip3); 
                crL(Ig1,Ig2,Ig3)  =  crL(Ip1,Ip2,Ip3); 
                csL(Ig1,Ig2,Ig3)  =  csL(Ip1,Ip2,Ip3); 
            }
            else;
            {
                if( ghost==3 )
                { 
                    crrL(Ig1,Ig2,Ig3) = (5.*crrL(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*crrL(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*crrL(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*crrL(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+crrL(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                    crsL(Ig1,Ig2,Ig3) = (5.*crsL(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*crsL(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*crsL(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*crsL(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+crsL(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                    cssL(Ig1,Ig2,Ig3) = (5.*cssL(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*cssL(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*cssL(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*cssL(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+cssL(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                } 
                crL(Ig1,Ig2,Ig3) = (5.*crL(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*crL(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*crL(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*crL(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+crL(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                csL(Ig1,Ig2,Ig3) = (5.*csL(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*csL(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*csL(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*csL(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+csL(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t
    axist = (axis2 + 1) % numberOfDimensions; // tangential direction
    for( int sidet=0; sidet<=1; sidet++ ) // top and bottom of face
    { 
        ks1=ks2=ks3=0;  ksv[axist]=1-2*sidet;
        I1=Range(gir2(0,0),gir2(0,0)); 
        I2=Range(gir2(0,1),gir2(0,1)); 
        I3=Range(gir2(0,2),gir2(0,2)); 
        Iv[axis2]=Range(gir2(side2,axis2)-1,gir2(side2,axis2)+1);
        Iv[axist]=gir2(sidet,axist);
        for( int ghost=1; ghost<=numGhost+1; ghost++ )
        { 
            Ig1=I1-ks1*ghost; Ig2=I2-ks2*ghost; Ig3=I3-ks3*ghost;
            if( mg2.boundaryCondition(sidet,axist)<0 )
            { // periodic
                Index Ip1 = Ig1 + (gir2(1,0)-gir2(0,0))*ks1;
                Index Ip2 = Ig2 + (gir2(1,1)-gir2(0,1))*ks2;
                Index Ip3 = Ig3 + (gir2(1,2)-gir2(0,2))*ks3;
                crrR(Ig1,Ig2,Ig3) = crrR(Ip1,Ip2,Ip3); 
                crsR(Ig1,Ig2,Ig3) = crsR(Ip1,Ip2,Ip3); 
                cssR(Ig1,Ig2,Ig3) = cssR(Ip1,Ip2,Ip3); 
                crR(Ig1,Ig2,Ig3)  =  crR(Ip1,Ip2,Ip3); 
                csR(Ig1,Ig2,Ig3)  =  csR(Ip1,Ip2,Ip3); 
            }
            else;
            {
                if( ghost==3 )
                { 
                    crrR(Ig1,Ig2,Ig3) = (5.*crrR(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*crrR(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*crrR(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*crrR(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+crrR(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                    crsR(Ig1,Ig2,Ig3) = (5.*crsR(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*crsR(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*crsR(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*crsR(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+crsR(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                    cssR(Ig1,Ig2,Ig3) = (5.*cssR(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*cssR(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*cssR(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*cssR(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+cssR(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                } 
                crR(Ig1,Ig2,Ig3) = (5.*crR(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*crR(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*crR(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*crR(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+crR(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
                csR(Ig1,Ig2,Ig3) = (5.*csR(Ig1+ks1,Ig2+ks2,Ig3+ks3)-10.*csR(Ig1+ks1+ks1,Ig2+ks2+ks2,Ig3+ks3+ks3)+10.*csR(Ig1+ks1+2*ks1,Ig2+ks2+2*ks2,Ig3+ks3+2*ks3)-5.*csR(Ig1+ks1+3*ks1,Ig2+ks2+3*ks2,Ig3+ks3+3*ks3)+csR(Ig1+ks1+4*ks1,Ig2+ks2+4*ks2,Ig3+ks3+4*ks3)); 
            }                  
        } // end for ghost 
    } // end for side t
                  
// ---- Evaluate derivatives of coefficients : (crr).rr, (csr).rs etc. 
//  L3 : needs (crr).r (crr).s etc at one tangential ghost. 
extra=numGhost-1;
getBoundaryIndex(mg.indexRange(),  side,axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray b1Lrr(Ib1,Ib2,Ib3), b1Lrs(Ib1,Ib2,Ib3), b1Lss(Ib1,Ib2,Ib3), b1Lr(Ib1,Ib2,Ib3), b1Ls(Ib1,Ib2,Ib3);
RealArray b2Lrr(Ib1,Ib2,Ib3), b2Lrs(Ib1,Ib2,Ib3), b2Lss(Ib1,Ib2,Ib3), b2Lr(Ib1,Ib2,Ib3), b2Ls(Ib1,Ib2,Ib3);
RealArray crrLrr(Ib1,Ib2,Ib3), crrLrs(Ib1,Ib2,Ib3), crrLss(Ib1,Ib2,Ib3), crrLr(Ib1,Ib2,Ib3), crrLs(Ib1,Ib2,Ib3);
RealArray crsLrr(Ib1,Ib2,Ib3), crsLrs(Ib1,Ib2,Ib3), crsLss(Ib1,Ib2,Ib3), crsLr(Ib1,Ib2,Ib3), crsLs(Ib1,Ib2,Ib3);
RealArray cssLrr(Ib1,Ib2,Ib3), cssLrs(Ib1,Ib2,Ib3), cssLss(Ib1,Ib2,Ib3), cssLr(Ib1,Ib2,Ib3), cssLs(Ib1,Ib2,Ib3);
RealArray crLrr(Ib1,Ib2,Ib3), crLrs(Ib1,Ib2,Ib3), crLss(Ib1,Ib2,Ib3), crLr(Ib1,Ib2,Ib3), crLs(Ib1,Ib2,Ib3);
RealArray csLrr(Ib1,Ib2,Ib3), csLrs(Ib1,Ib2,Ib3), csLss(Ib1,Ib2,Ib3), csLr(Ib1,Ib2,Ib3), csLs(Ib1,Ib2,Ib3);
// Evaluate derivatives of coefficients at points on the boundary (left)
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  // tangential directives of b1,b2 only 
    b1Lrr(i1,i2,i3) = ((b1L(i1+1,i2,i3)-2.*b1L(i1,i2,i3)+b1L(i1-1,i2,i3))/(SQR(drL[0])));
      b1Lr(i1,i2,i3) =  ((b1L(i1+1,i2,i3)-b1L(i1-1,i2,i3))/(2.*drL[0]));
    b2Lrr(i1,i2,i3) = ((b2L(i1+1,i2,i3)-2.*b2L(i1,i2,i3)+b2L(i1-1,i2,i3))/(SQR(drL[0])));
      b2Lr(i1,i2,i3) =  ((b2L(i1+1,i2,i3)-b2L(i1-1,i2,i3))/(2.*drL[0]));
  // normal and tangential directives of c 
    crrLrr(i1,i2,i3) = ((crrL(i1+1,i2,i3)-2.*crrL(i1,i2,i3)+crrL(i1-1,i2,i3))/(SQR(drL[0])));
    crrLss(i1,i2,i3) = ((crrL(i1,i2+1,i3)-2.*crrL(i1,i2,i3)+crrL(i1,i2-1,i3))/(SQR(drL[1])));
    crrLrs(i1,i2,i3) = ((((crrL(i1+1,i2+1,i3)-crrL(i1+1,i2-1,i3))/(2.*drL[1]))-((crrL(i1-1,i2+1,i3)-crrL(i1-1,i2-1,i3))/(2.*drL[1])))/(2.*drL[0]));
      crrLr(i1,i2,i3) =  ((crrL(i1+1,i2,i3)-crrL(i1-1,i2,i3))/(2.*drL[0]));
      crrLs(i1,i2,i3) =  ((crrL(i1,i2+1,i3)-crrL(i1,i2-1,i3))/(2.*drL[1]));
    crsLrr(i1,i2,i3) = ((crsL(i1+1,i2,i3)-2.*crsL(i1,i2,i3)+crsL(i1-1,i2,i3))/(SQR(drL[0])));
    crsLss(i1,i2,i3) = ((crsL(i1,i2+1,i3)-2.*crsL(i1,i2,i3)+crsL(i1,i2-1,i3))/(SQR(drL[1])));
    crsLrs(i1,i2,i3) = ((((crsL(i1+1,i2+1,i3)-crsL(i1+1,i2-1,i3))/(2.*drL[1]))-((crsL(i1-1,i2+1,i3)-crsL(i1-1,i2-1,i3))/(2.*drL[1])))/(2.*drL[0]));
      crsLr(i1,i2,i3) =  ((crsL(i1+1,i2,i3)-crsL(i1-1,i2,i3))/(2.*drL[0]));
      crsLs(i1,i2,i3) =  ((crsL(i1,i2+1,i3)-crsL(i1,i2-1,i3))/(2.*drL[1]));
    cssLrr(i1,i2,i3) = ((cssL(i1+1,i2,i3)-2.*cssL(i1,i2,i3)+cssL(i1-1,i2,i3))/(SQR(drL[0])));
    cssLss(i1,i2,i3) = ((cssL(i1,i2+1,i3)-2.*cssL(i1,i2,i3)+cssL(i1,i2-1,i3))/(SQR(drL[1])));
    cssLrs(i1,i2,i3) = ((((cssL(i1+1,i2+1,i3)-cssL(i1+1,i2-1,i3))/(2.*drL[1]))-((cssL(i1-1,i2+1,i3)-cssL(i1-1,i2-1,i3))/(2.*drL[1])))/(2.*drL[0]));
      cssLr(i1,i2,i3) =  ((cssL(i1+1,i2,i3)-cssL(i1-1,i2,i3))/(2.*drL[0]));
      cssLs(i1,i2,i3) =  ((cssL(i1,i2+1,i3)-cssL(i1,i2-1,i3))/(2.*drL[1]));
    crLrr(i1,i2,i3) = ((crL(i1+1,i2,i3)-2.*crL(i1,i2,i3)+crL(i1-1,i2,i3))/(SQR(drL[0])));
    crLss(i1,i2,i3) = ((crL(i1,i2+1,i3)-2.*crL(i1,i2,i3)+crL(i1,i2-1,i3))/(SQR(drL[1])));
    crLrs(i1,i2,i3) = ((((crL(i1+1,i2+1,i3)-crL(i1+1,i2-1,i3))/(2.*drL[1]))-((crL(i1-1,i2+1,i3)-crL(i1-1,i2-1,i3))/(2.*drL[1])))/(2.*drL[0]));
      crLr(i1,i2,i3) =  ((crL(i1+1,i2,i3)-crL(i1-1,i2,i3))/(2.*drL[0]));
      crLs(i1,i2,i3) =  ((crL(i1,i2+1,i3)-crL(i1,i2-1,i3))/(2.*drL[1]));
    csLrr(i1,i2,i3) = ((csL(i1+1,i2,i3)-2.*csL(i1,i2,i3)+csL(i1-1,i2,i3))/(SQR(drL[0])));
    csLss(i1,i2,i3) = ((csL(i1,i2+1,i3)-2.*csL(i1,i2,i3)+csL(i1,i2-1,i3))/(SQR(drL[1])));
    csLrs(i1,i2,i3) = ((((csL(i1+1,i2+1,i3)-csL(i1+1,i2-1,i3))/(2.*drL[1]))-((csL(i1-1,i2+1,i3)-csL(i1-1,i2-1,i3))/(2.*drL[1])))/(2.*drL[0]));
      csLr(i1,i2,i3) =  ((csL(i1+1,i2,i3)-csL(i1-1,i2,i3))/(2.*drL[0]));
      csLs(i1,i2,i3) =  ((csL(i1,i2+1,i3)-csL(i1,i2-1,i3))/(2.*drL[1]));
} // end FOR_3D
RealArray b1Rrr(Jb1,Jb2,Jb3), b1Rrs(Jb1,Jb2,Jb3), b1Rss(Jb1,Jb2,Jb3), b1Rr(Jb1,Jb2,Jb3), b1Rs(Jb1,Jb2,Jb3);
RealArray b2Rrr(Jb1,Jb2,Jb3), b2Rrs(Jb1,Jb2,Jb3), b2Rss(Jb1,Jb2,Jb3), b2Rr(Jb1,Jb2,Jb3), b2Rs(Jb1,Jb2,Jb3);
RealArray crrRrr(Jb1,Jb2,Jb3), crrRrs(Jb1,Jb2,Jb3), crrRss(Jb1,Jb2,Jb3), crrRr(Jb1,Jb2,Jb3), crrRs(Jb1,Jb2,Jb3);
RealArray crsRrr(Jb1,Jb2,Jb3), crsRrs(Jb1,Jb2,Jb3), crsRss(Jb1,Jb2,Jb3), crsRr(Jb1,Jb2,Jb3), crsRs(Jb1,Jb2,Jb3);
RealArray cssRrr(Jb1,Jb2,Jb3), cssRrs(Jb1,Jb2,Jb3), cssRss(Jb1,Jb2,Jb3), cssRr(Jb1,Jb2,Jb3), cssRs(Jb1,Jb2,Jb3);
RealArray crRrr(Jb1,Jb2,Jb3), crRrs(Jb1,Jb2,Jb3), crRss(Jb1,Jb2,Jb3), crRr(Jb1,Jb2,Jb3), crRs(Jb1,Jb2,Jb3);
RealArray csRrr(Jb1,Jb2,Jb3), csRrs(Jb1,Jb2,Jb3), csRss(Jb1,Jb2,Jb3), csRr(Jb1,Jb2,Jb3), csRs(Jb1,Jb2,Jb3);
// Evaluate derivatives of coefficients at points on the boundary (right)
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
 // tangential directives of b1,b2 only 
    b1Rrr(j1,j2,j3) = ((b1R(j1+1,j2,j3)-2.*b1R(j1,j2,j3)+b1R(j1-1,j2,j3))/(SQR(drR[0])));
      b1Rr(j1,j2,j3) =  ((b1R(j1+1,j2,j3)-b1R(j1-1,j2,j3))/(2.*drR[0]));
    b2Rrr(j1,j2,j3) = ((b2R(j1+1,j2,j3)-2.*b2R(j1,j2,j3)+b2R(j1-1,j2,j3))/(SQR(drR[0])));
      b2Rr(j1,j2,j3) =  ((b2R(j1+1,j2,j3)-b2R(j1-1,j2,j3))/(2.*drR[0]));
  // normal and tangential directives of c 
    crrRrr(j1,j2,j3) = ((crrR(j1+1,j2,j3)-2.*crrR(j1,j2,j3)+crrR(j1-1,j2,j3))/(SQR(drR[0])));
    crrRss(j1,j2,j3) = ((crrR(j1,j2+1,j3)-2.*crrR(j1,j2,j3)+crrR(j1,j2-1,j3))/(SQR(drR[1])));
    crrRrs(j1,j2,j3) = ((((crrR(j1+1,j2+1,j3)-crrR(j1+1,j2-1,j3))/(2.*drR[1]))-((crrR(j1-1,j2+1,j3)-crrR(j1-1,j2-1,j3))/(2.*drR[1])))/(2.*drR[0]));
      crrRr(j1,j2,j3) =  ((crrR(j1+1,j2,j3)-crrR(j1-1,j2,j3))/(2.*drR[0]));
      crrRs(j1,j2,j3) =  ((crrR(j1,j2+1,j3)-crrR(j1,j2-1,j3))/(2.*drR[1]));
    crsRrr(j1,j2,j3) = ((crsR(j1+1,j2,j3)-2.*crsR(j1,j2,j3)+crsR(j1-1,j2,j3))/(SQR(drR[0])));
    crsRss(j1,j2,j3) = ((crsR(j1,j2+1,j3)-2.*crsR(j1,j2,j3)+crsR(j1,j2-1,j3))/(SQR(drR[1])));
    crsRrs(j1,j2,j3) = ((((crsR(j1+1,j2+1,j3)-crsR(j1+1,j2-1,j3))/(2.*drR[1]))-((crsR(j1-1,j2+1,j3)-crsR(j1-1,j2-1,j3))/(2.*drR[1])))/(2.*drR[0]));
      crsRr(j1,j2,j3) =  ((crsR(j1+1,j2,j3)-crsR(j1-1,j2,j3))/(2.*drR[0]));
      crsRs(j1,j2,j3) =  ((crsR(j1,j2+1,j3)-crsR(j1,j2-1,j3))/(2.*drR[1]));
    cssRrr(j1,j2,j3) = ((cssR(j1+1,j2,j3)-2.*cssR(j1,j2,j3)+cssR(j1-1,j2,j3))/(SQR(drR[0])));
    cssRss(j1,j2,j3) = ((cssR(j1,j2+1,j3)-2.*cssR(j1,j2,j3)+cssR(j1,j2-1,j3))/(SQR(drR[1])));
    cssRrs(j1,j2,j3) = ((((cssR(j1+1,j2+1,j3)-cssR(j1+1,j2-1,j3))/(2.*drR[1]))-((cssR(j1-1,j2+1,j3)-cssR(j1-1,j2-1,j3))/(2.*drR[1])))/(2.*drR[0]));
      cssRr(j1,j2,j3) =  ((cssR(j1+1,j2,j3)-cssR(j1-1,j2,j3))/(2.*drR[0]));
      cssRs(j1,j2,j3) =  ((cssR(j1,j2+1,j3)-cssR(j1,j2-1,j3))/(2.*drR[1]));
    crRrr(j1,j2,j3) = ((crR(j1+1,j2,j3)-2.*crR(j1,j2,j3)+crR(j1-1,j2,j3))/(SQR(drR[0])));
    crRss(j1,j2,j3) = ((crR(j1,j2+1,j3)-2.*crR(j1,j2,j3)+crR(j1,j2-1,j3))/(SQR(drR[1])));
    crRrs(j1,j2,j3) = ((((crR(j1+1,j2+1,j3)-crR(j1+1,j2-1,j3))/(2.*drR[1]))-((crR(j1-1,j2+1,j3)-crR(j1-1,j2-1,j3))/(2.*drR[1])))/(2.*drR[0]));
      crRr(j1,j2,j3) =  ((crR(j1+1,j2,j3)-crR(j1-1,j2,j3))/(2.*drR[0]));
      crRs(j1,j2,j3) =  ((crR(j1,j2+1,j3)-crR(j1,j2-1,j3))/(2.*drR[1]));
    csRrr(j1,j2,j3) = ((csR(j1+1,j2,j3)-2.*csR(j1,j2,j3)+csR(j1-1,j2,j3))/(SQR(drR[0])));
    csRss(j1,j2,j3) = ((csR(j1,j2+1,j3)-2.*csR(j1,j2,j3)+csR(j1,j2-1,j3))/(SQR(drR[1])));
    csRrs(j1,j2,j3) = ((((csR(j1+1,j2+1,j3)-csR(j1+1,j2-1,j3))/(2.*drR[1]))-((csR(j1-1,j2+1,j3)-csR(j1-1,j2-1,j3))/(2.*drR[1])))/(2.*drR[0]));
      csRr(j1,j2,j3) =  ((csR(j1+1,j2,j3)-csR(j1-1,j2,j3))/(2.*drR[0]));
      csRs(j1,j2,j3) =  ((csR(j1,j2+1,j3)-csR(j1,j2-1,j3))/(2.*drR[1]));
} // end FOR_3IJ

// ----- Evaluate coefficients of L1 on the boundary -----

// 2*numGhost = 4 extra ghost needed for fourth order 
extra=2*numGhost;
getBoundaryIndex(mg.indexRange(), side ,axis ,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L1ICoeff(Ib1,Ib2,Ib3), L1ICoeffr(Ib1,Ib2,Ib3), L1ICoeffrr(Ib1,Ib2,Ib3), L1ICoeffrrr(Ib1,Ib2,Ib3);
RealArray L1rCoeff(Ib1,Ib2,Ib3), L1rCoeffr(Ib1,Ib2,Ib3), L1rCoeffrr(Ib1,Ib2,Ib3), L1rCoeffrrr(Ib1,Ib2,Ib3);
RealArray L1sCoeff(Ib1,Ib2,Ib3), L1sCoeffr(Ib1,Ib2,Ib3), L1sCoeffrr(Ib1,Ib2,Ib3), L1sCoeffrrr(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    const Real n1R = -n1L(i1,i2,i3), n2R = -n2L(i1,i2,i3);    // for advection, use right normal 
        if( variableAdvection )
        {
            u1DL = advectVarL(i1,i2,i3,0)/DL;
            u2DL = advectVarL(i1,i2,i3,1)/DL;
            u1DR = advectVarR(j1,j2,j3,0)/DR;
            u2DR = advectVarR(j1,j2,j3,1)/DR;
        }
    L1ICoeff(i1,i2,i3) = ((-1.*KLR*u2DL+1.*u2DR)*n2R+(-1.*KLR*u1DL+1.*u1DR)*n1R)/b2R(j1,j2,j3);
    L1rCoeff(i1,i2,i3) = (1.*KLR*b1L(i1,i2,i3)-1.*b1R(j1,j2,j3))/b2R(j1,j2,j3);
    L1sCoeff(i1,i2,i3) = 1.*b2L(i1,i2,i3)*KLR/b2R(j1,j2,j3);
} // end FOR_3D

// ----- Evaluate coefficients of L2 on the boundary -----
extra=numGhost;
getBoundaryIndex(mg.indexRange(), side,  axis,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L2rrCoeff(Ib1,Ib2,Ib3), L2rrCoeffr(Ib1,Ib2,Ib3), L2rrCoeffrr(Ib1,Ib2,Ib3);
RealArray L2rsCoeff(Ib1,Ib2,Ib3), L2rsCoeffr(Ib1,Ib2,Ib3), L2rsCoeffrr(Ib1,Ib2,Ib3);
RealArray L2ssCoeff(Ib1,Ib2,Ib3), L2ssCoeffr(Ib1,Ib2,Ib3), L2ssCoeffrr(Ib1,Ib2,Ib3);
RealArray L2rCoeff(Ib1,Ib2,Ib3), L2rCoeffr(Ib1,Ib2,Ib3), L2rCoeffrr(Ib1,Ib2,Ib3);
RealArray L2sCoeff(Ib1,Ib2,Ib3), L2sCoeffr(Ib1,Ib2,Ib3), L2sCoeffrr(Ib1,Ib2,Ib3);
RealArray L2ICoeff(Ib1,Ib2,Ib3), L2ICoeffr(Ib1,Ib2,Ib3), L2ICoeffrr(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  // Compute tangential derivatives of L1 - label tangential directions a and b (for 3D)
    L1rCoeffr(i1,i2,i3) = (((1./12.)*L1rCoeff(i1-2,i2,i3)-(2./3.)*L1rCoeff(i1-1,i2,i3)+(2./3.)*L1rCoeff(i1+1,i2,i3)-(1./12.)*L1rCoeff(i1+2,i2,i3))/(drL[0]));
    L1sCoeffr(i1,i2,i3) = (((1./12.)*L1sCoeff(i1-2,i2,i3)-(2./3.)*L1sCoeff(i1-1,i2,i3)+(2./3.)*L1sCoeff(i1+1,i2,i3)-(1./12.)*L1sCoeff(i1+2,i2,i3))/(drL[0]));
    L1ICoeffr(i1,i2,i3) = (((1./12.)*L1ICoeff(i1-2,i2,i3)-(2./3.)*L1ICoeff(i1-1,i2,i3)+(2./3.)*L1ICoeff(i1+1,i2,i3)-(1./12.)*L1ICoeff(i1+2,i2,i3))/(drL[0]));
  // Compute L2 coefficients
    L2rrCoeff(i1,i2,i3) = 1.*(DLR*crrL(i1,i2,i3)-crsR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-crrR(j1,j2,j3))/cssR(j1,j2,j3);
    L2rsCoeff(i1,i2,i3) = (1.*DLR*crsL(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1sCoeff(i1,i2,i3))/cssR(j1,j2,j3);
    L2ssCoeff(i1,i2,i3) = 1.*cssL(i1,i2,i3)*DLR/cssR(j1,j2,j3);
    L2rCoeff(i1,i2,i3) = ((-1.*L1ICoeff(i1,i2,i3)-1.*L1rCoeffr(i1,i2,i3))*crsR(j1,j2,j3)+1.*DLR*crL(i1,i2,i3)-1.*csR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*crR(j1,j2,j3))/cssR(j1,j2,j3);
    L2sCoeff(i1,i2,i3) = 1.*(DLR*csL(i1,i2,i3)-crsR(j1,j2,j3)*L1sCoeffr(i1,i2,i3)-csR(j1,j2,j3)*L1sCoeff(i1,i2,i3))/cssR(j1,j2,j3);
    L2ICoeff(i1,i2,i3) = (-1.*crsR(j1,j2,j3)*L1ICoeffr(i1,i2,i3)-1.*csR(j1,j2,j3)*L1ICoeff(i1,i2,i3))/cssR(j1,j2,j3);
} // end FOR_3D

// ----- Evaluate coefficients of L3 on the boundary -----
extra=numGhost-1;
getBoundaryIndex(mg.indexRange(),side,axis,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L3rrrCoeff(Ib1,Ib2,Ib3), L3rrrCoeffr(Ib1,Ib2,Ib3);
RealArray L3rrsCoeff(Ib1,Ib2,Ib3), L3rrsCoeffr(Ib1,Ib2,Ib3);
RealArray L3rssCoeff(Ib1,Ib2,Ib3), L3rssCoeffr(Ib1,Ib2,Ib3);
RealArray L3sssCoeff(Ib1,Ib2,Ib3), L3sssCoeffr(Ib1,Ib2,Ib3);
RealArray L3rrCoeff(Ib1,Ib2,Ib3), L3rrCoeffr(Ib1,Ib2,Ib3);
RealArray L3rsCoeff(Ib1,Ib2,Ib3), L3rsCoeffr(Ib1,Ib2,Ib3);
RealArray L3ssCoeff(Ib1,Ib2,Ib3), L3ssCoeffr(Ib1,Ib2,Ib3);
RealArray L3rCoeff(Ib1,Ib2,Ib3), L3rCoeffr(Ib1,Ib2,Ib3);
RealArray L3sCoeff(Ib1,Ib2,Ib3), L3sCoeffr(Ib1,Ib2,Ib3);
RealArray L3ICoeff(Ib1,Ib2,Ib3), L3ICoeffr(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    const Real n1R = -n1L(i1,i2,i3), n2R = -n2L(i1,i2,i3);    // for advection, use right normal 
        if( variableAdvection )
        {
            u1DL = advectVarL(i1,i2,i3,0)/DL;
            u2DL = advectVarL(i1,i2,i3,1)/DL;
            u1DR = advectVarR(j1,j2,j3,0)/DR;
            u2DR = advectVarR(j1,j2,j3,1)/DR;
        }
  // Compute tangential derivatives of L1 
    L1rCoeffrr(i1,i2,i3) = ((L1rCoeff(i1+1,i2,i3)-2.*L1rCoeff(i1,i2,i3)+L1rCoeff(i1-1,i2,i3))/(SQR(drL[0])));
    L1sCoeffrr(i1,i2,i3) = ((L1sCoeff(i1+1,i2,i3)-2.*L1sCoeff(i1,i2,i3)+L1sCoeff(i1-1,i2,i3))/(SQR(drL[0])));
    L1ICoeffrr(i1,i2,i3) = ((L1ICoeff(i1+1,i2,i3)-2.*L1ICoeff(i1,i2,i3)+L1ICoeff(i1-1,i2,i3))/(SQR(drL[0])));
  // Compute tangential derivatives of L2
    L2rrCoeffr(i1,i2,i3) = ((L2rrCoeff(i1+1,i2,i3)-L2rrCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L2rsCoeffr(i1,i2,i3) = ((L2rsCoeff(i1+1,i2,i3)-L2rsCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L2ssCoeffr(i1,i2,i3) = ((L2ssCoeff(i1+1,i2,i3)-L2ssCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L2rCoeffr(i1,i2,i3)  = ((L2rCoeff(i1+1,i2,i3)-L2rCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L2sCoeffr(i1,i2,i3)  = ((L2sCoeff(i1+1,i2,i3)-L2sCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L2ICoeffr(i1,i2,i3)  = ((L2ICoeff(i1+1,i2,i3)-L2ICoeff(i1-1,i2,i3))/(2.*drL[0]));
  // Compute L3 coefficients
    L3rrrCoeff(i1,i2,i3) = ((-1.*crsR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*cssR(j1,j2,j3)*L2rrCoeff(i1,i2,i3)-1.*crrR(j1,j2,j3))*b1R(j1,j2,j3)+(-1.*L1rCoeff(i1,i2,i3)*crrR(j1,j2,j3)-1.*L2rrCoeff(i1,i2,i3)*crsR(j1,j2,j3))*b2R(j1,j2,j3)+1.*DLR*KLR*b1L(i1,i2,i3)*crrL(i1,i2,i3))/b2R(j1,j2,j3)/cssR(j1,j2,j3);
    L3rrsCoeff(i1,i2,i3) = 1.*(DLR*KLR*b1L(i1,i2,i3)*crsL(i1,i2,i3)+DLR*KLR*b2L(i1,i2,i3)*crrL(i1,i2,i3)-b1R(j1,j2,j3)*crsR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-b1R(j1,j2,j3)*cssR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)-b2R(j1,j2,j3)*crrR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-b2R(j1,j2,j3)*crsR(j1,j2,j3)*L2rsCoeff(i1,i2,i3))/b2R(j1,j2,j3)/cssR(j1,j2,j3);
    L3rssCoeff(i1,i2,i3) = ((-1.*cssR(j1,j2,j3)*b1R(j1,j2,j3)-1.*b2R(j1,j2,j3)*crsR(j1,j2,j3))*L2ssCoeff(i1,i2,i3)+(1.*cssL(i1,i2,i3)*b1L(i1,i2,i3)+1.*b2L(i1,i2,i3)*crsL(i1,i2,i3))*KLR*DLR)/b2R(j1,j2,j3)/cssR(j1,j2,j3);
    L3sssCoeff(i1,i2,i3) = 1.*cssL(i1,i2,i3)*b2L(i1,i2,i3)*KLR*DLR/b2R(j1,j2,j3)/cssR(j1,j2,j3);
    L3rrCoeff(i1,i2,i3) = (((-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1rCoeff(i1,i2,i3)+(-1.*L2rCoeff(i1,i2,i3)-1.*L2rrCoeffr(i1,i2,i3))*crsR(j1,j2,j3)+(-1.*csR(j1,j2,j3)-1.*cssRs(j1,j2,j3))*L2rrCoeff(i1,i2,i3)+(-1.*L1ICoeff(i1,i2,i3)-2.*L1rCoeffr(i1,i2,i3))*crrR(j1,j2,j3)-1.*crrRs(j1,j2,j3))*b2R(j1,j2,j3)+((-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1rCoeff(i1,i2,i3)+(-1.*L1ICoeff(i1,i2,i3)-2.*L1rCoeffr(i1,i2,i3))*crsR(j1,j2,j3)-1.*cssRr(j1,j2,j3)*L2rrCoeff(i1,i2,i3)+(-1.*L2rCoeff(i1,i2,i3)-1.*L2rrCoeffr(i1,i2,i3))*cssR(j1,j2,j3)-1.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3))*b1R(j1,j2,j3)+crsR(j1,j2,j3)*(n1R*u1DR+n2R*u2DR)*L1rCoeff(i1,i2,i3)+cssR(j1,j2,j3)*(n1R*u1DR+n2R*u2DR)*L2rrCoeff(i1,i2,i3)+(-1.*crrL(i1,i2,i3)*n1R*u1DL-1.*crrL(i1,i2,i3)*n2R*u2DL+(crL(i1,i2,i3)+crrLr(i1,i2,i3))*b1L(i1,i2,i3)+b2L(i1,i2,i3)*crrLs(i1,i2,i3))*KLR*DLR+crrR(j1,j2,j3)*(n1R*u1DR+n2R*u2DR))/cssR(j1,j2,j3)/b2R(j1,j2,j3);
    L3rsCoeff(i1,i2,i3) = (((-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1sCoeff(i1,i2,i3)+(-1.*csR(j1,j2,j3)-1.*cssRs(j1,j2,j3))*L2rsCoeff(i1,i2,i3)+(-1.*L2sCoeff(i1,i2,i3)-1.*L2rsCoeffr(i1,i2,i3))*crsR(j1,j2,j3)-2.*crrR(j1,j2,j3)*L1sCoeffr(i1,i2,i3))*b2R(j1,j2,j3)+((-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*b1R(j1,j2,j3)+crsR(j1,j2,j3)*(n1R*u1DR+n2R*u2DR))*L1sCoeff(i1,i2,i3)+(-1.*crsL(i1,i2,i3)*n1R*u1DL-1.*crsL(i1,i2,i3)*n2R*u2DL+(crsLr(i1,i2,i3)+csL(i1,i2,i3))*b1L(i1,i2,i3)+b2L(i1,i2,i3)*(crsLs(i1,i2,i3)+crL(i1,i2,i3)))*KLR*DLR+(-1.*cssRr(j1,j2,j3)*L2rsCoeff(i1,i2,i3)-2.*crsR(j1,j2,j3)*L1sCoeffr(i1,i2,i3)+(-1.*L2sCoeff(i1,i2,i3)-1.*L2rsCoeffr(i1,i2,i3))*cssR(j1,j2,j3))*b1R(j1,j2,j3)+cssR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)*(n1R*u1DR+n2R*u2DR))/b2R(j1,j2,j3)/cssR(j1,j2,j3);
    L3ssCoeff(i1,i2,i3) = (((1.*n1R*u1DR+1.*n2R*u2DR)*cssR(j1,j2,j3)+(-1.*csR(j1,j2,j3)-1.*cssRs(j1,j2,j3))*b2R(j1,j2,j3)-1.*cssRr(j1,j2,j3)*b1R(j1,j2,j3))*L2ssCoeff(i1,i2,i3)+((-1.*n1R*u1DL-1.*n2R*u2DL)*cssL(i1,i2,i3)+(1.*csL(i1,i2,i3)+1.*cssLs(i1,i2,i3))*b2L(i1,i2,i3)+1.*cssLr(i1,i2,i3)*b1L(i1,i2,i3))*KLR*DLR-1.*b1R(j1,j2,j3)*cssR(j1,j2,j3)*L2ssCoeffr(i1,i2,i3)-1.*b2R(j1,j2,j3)*crsR(j1,j2,j3)*L2ssCoeffr(i1,i2,i3))/cssR(j1,j2,j3)/b2R(j1,j2,j3);
    L3rCoeff(i1,i2,i3) = (((-1.*L2rCoeffr(i1,i2,i3)-1.*L2ICoeff(i1,i2,i3))*crsR(j1,j2,j3)+(-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1rCoeffr(i1,i2,i3)+(-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1ICoeff(i1,i2,i3)+(-1.*csR(j1,j2,j3)-1.*cssRs(j1,j2,j3))*L2rCoeff(i1,i2,i3)-1.*crrR(j1,j2,j3)*L1rCoeffrr(i1,i2,i3)-1.*csRs(j1,j2,j3)*L1rCoeff(i1,i2,i3)-2.*crrR(j1,j2,j3)*L1ICoeffr(i1,i2,i3)-1.*crRs(j1,j2,j3))*b2R(j1,j2,j3)+((-2.*L1ICoeffr(i1,i2,i3)-1.*L1rCoeffrr(i1,i2,i3))*crsR(j1,j2,j3)+(-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1rCoeffr(i1,i2,i3)+(-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*L1ICoeff(i1,i2,i3)-1.*cssRr(j1,j2,j3)*L2rCoeff(i1,i2,i3)+(-1.*L2rCoeffr(i1,i2,i3)-1.*L2ICoeff(i1,i2,i3))*cssR(j1,j2,j3)-1.*csRr(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*crRr(j1,j2,j3))*b1R(j1,j2,j3)+(L1ICoeff(i1,i2,i3)+L1rCoeffr(i1,i2,i3))*(n1R*u1DR+n2R*u2DR)*crsR(j1,j2,j3)+(u1DR*L2rCoeff(i1,i2,i3)*cssR(j1,j2,j3)+(csR(j1,j2,j3)*L1rCoeff(i1,i2,i3)+crR(j1,j2,j3))*u1DR-1.*crL(i1,i2,i3)*u1DL*KLR*DLR)*n1R+(u2DR*L2rCoeff(i1,i2,i3)*cssR(j1,j2,j3)+(csR(j1,j2,j3)*L1rCoeff(i1,i2,i3)+crR(j1,j2,j3))*u2DR-1.*crL(i1,i2,i3)*u2DL*KLR*DLR)*n2R+DLR*KLR*(b1L(i1,i2,i3)*crLr(i1,i2,i3)+b2L(i1,i2,i3)*crLs(i1,i2,i3)))/b2R(j1,j2,j3)/cssR(j1,j2,j3);
    L3sCoeff(i1,i2,i3) = (((-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1sCoeffr(i1,i2,i3)+(-1.*csR(j1,j2,j3)-1.*cssRs(j1,j2,j3))*L2sCoeff(i1,i2,i3)-1.*csRs(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*crrR(j1,j2,j3)*L1sCoeffrr(i1,i2,i3)-1.*crsR(j1,j2,j3)*L2sCoeffr(i1,i2,i3))*b2R(j1,j2,j3)+((-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*b1R(j1,j2,j3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crsR(j1,j2,j3))*L1sCoeffr(i1,i2,i3)+(-1.*cssRr(j1,j2,j3)*L2sCoeff(i1,i2,i3)-1.*csRr(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*crsR(j1,j2,j3)*L1sCoeffrr(i1,i2,i3)-1.*cssR(j1,j2,j3)*L2sCoeffr(i1,i2,i3))*b1R(j1,j2,j3)+(1.*n1R*u1DR+1.*n2R*u2DR)*cssR(j1,j2,j3)*L2sCoeff(i1,i2,i3)+(1.*n1R*u1DR+1.*n2R*u2DR)*csR(j1,j2,j3)*L1sCoeff(i1,i2,i3)+(-1.*csL(i1,i2,i3)*n1R*u1DL-1.*csL(i1,i2,i3)*n2R*u2DL+1.*b1L(i1,i2,i3)*csLr(i1,i2,i3)+1.*b2L(i1,i2,i3)*csLs(i1,i2,i3))*KLR*DLR)/b2R(j1,j2,j3)/cssR(j1,j2,j3);
    L3ICoeff(i1,i2,i3) = (((-1.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1ICoeffr(i1,i2,i3)+(-1.*csR(j1,j2,j3)-1.*cssRs(j1,j2,j3))*L2ICoeff(i1,i2,i3)-1.*L1ICoeffrr(i1,i2,i3)*crrR(j1,j2,j3)-1.*L2ICoeffr(i1,i2,i3)*crsR(j1,j2,j3)-1.*csRs(j1,j2,j3)*L1ICoeff(i1,i2,i3))*b2R(j1,j2,j3)+((-1.*crsRr(j1,j2,j3)-1.*csR(j1,j2,j3))*b1R(j1,j2,j3)+(1.*n1R*u1DR+1.*n2R*u2DR)*crsR(j1,j2,j3))*L1ICoeffr(i1,i2,i3)+(-1.*crsR(j1,j2,j3)*L1ICoeffrr(i1,i2,i3)-1.*csRr(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*cssR(j1,j2,j3)*L2ICoeffr(i1,i2,i3)-1.*cssRr(j1,j2,j3)*L2ICoeff(i1,i2,i3))*b1R(j1,j2,j3)+(1.*n1R*u1DR+1.*n2R*u2DR)*cssR(j1,j2,j3)*L2ICoeff(i1,i2,i3)+(1.*n1R*u1DR+1.*n2R*u2DR)*csR(j1,j2,j3)*L1ICoeff(i1,i2,i3))/b2R(j1,j2,j3)/cssR(j1,j2,j3);
    printF("L3s=%9.2e L1sCoeff=%9.2e L1sCoeffr=%9.2e L1sCoeffrr=%9.2e\n",L3sCoeff(i1,i2,i3),L1sCoeff(i1,i2,i3),L1sCoeffr(i1,i2,i3),L1sCoeffrr(i1,i2,i3));
    printF("L2sCoeff=%9.2e L2sCoeffr=%9.2e \n",L2sCoeff(i1,i2,i3),L2sCoeffr(i1,i2,i3));
    printF("b1R=%9.2e, b2r=%9.2e, crR=%9.2e csR=%9.2e cssR=%9.2e csRs=%9.2e, csR=%9.2e\n",b1R(j1,j2,j3),b2R(j1,j2,j3),crR(j1,j2,j3),csR(j1,j2,j3),cssRs(j1,j2,j3),csRs(j1,j2,j3),csR(j1,j2,j3));
} // end FOR_3D

// ----- Evaluate coefficients of L4 on the boundary -----
extra=0;
getBoundaryIndex(mg.indexRange(), side,axis  ,Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
RealArray L4rrrrCoeff(Ib1,Ib2,Ib3);
RealArray L4rrrsCoeff(Ib1,Ib2,Ib3);
RealArray L4rrssCoeff(Ib1,Ib2,Ib3);
RealArray L4rsssCoeff(Ib1,Ib2,Ib3);
RealArray L4ssssCoeff(Ib1,Ib2,Ib3);
RealArray L4rrrCoeff(Ib1,Ib2,Ib3);
RealArray L4rrsCoeff(Ib1,Ib2,Ib3);
RealArray L4rssCoeff(Ib1,Ib2,Ib3);
RealArray L4sssCoeff(Ib1,Ib2,Ib3);
RealArray L4rrCoeff(Ib1,Ib2,Ib3);
RealArray L4rsCoeff(Ib1,Ib2,Ib3);
RealArray L4ssCoeff(Ib1,Ib2,Ib3);
RealArray L4rCoeff(Ib1,Ib2,Ib3);
RealArray L4sCoeff(Ib1,Ib2,Ib3);
RealArray L4ICoeff(Ib1,Ib2,Ib3);
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
  // Compute tangential derivatives of L1
    L1rCoeffrrr(i1,i2,i3) = (((-1./2.)*L1rCoeff(i1-2,i2,i3)+1.*L1rCoeff(i1-1,i2,i3)-1.*L1rCoeff(i1+1,i2,i3)+(1./2.)*L1rCoeff(i1+2,i2,i3))/(pow(drL[0],3)));
    L1sCoeffrrr(i1,i2,i3) = (((-1./2.)*L1sCoeff(i1-2,i2,i3)+1.*L1sCoeff(i1-1,i2,i3)-1.*L1sCoeff(i1+1,i2,i3)+(1./2.)*L1sCoeff(i1+2,i2,i3))/(pow(drL[0],3)));
    L1ICoeffrrr(i1,i2,i3) = (((-1./2.)*L1ICoeff(i1-2,i2,i3)+1.*L1ICoeff(i1-1,i2,i3)-1.*L1ICoeff(i1+1,i2,i3)+(1./2.)*L1ICoeff(i1+2,i2,i3))/(pow(drL[0],3)));
  // Compute tangential derivatives of L2
    L2rrCoeffrr(i1,i2,i3) = ((L2rrCoeff(i1+1,i2,i3)-2.*L2rrCoeff(i1,i2,i3)+L2rrCoeff(i1-1,i2,i3))/(SQR(drL[0])));
    L2rsCoeffrr(i1,i2,i3) = ((L2rsCoeff(i1+1,i2,i3)-2.*L2rsCoeff(i1,i2,i3)+L2rsCoeff(i1-1,i2,i3))/(SQR(drL[0])));
    L2ssCoeffrr(i1,i2,i3) = ((L2ssCoeff(i1+1,i2,i3)-2.*L2ssCoeff(i1,i2,i3)+L2ssCoeff(i1-1,i2,i3))/(SQR(drL[0])));
    L2rCoeffrr(i1,i2,i3)  = ((L2rCoeff(i1+1,i2,i3)-2.*L2rCoeff(i1,i2,i3)+L2rCoeff(i1-1,i2,i3))/(SQR(drL[0])));
    L2sCoeffrr(i1,i2,i3)  = ((L2sCoeff(i1+1,i2,i3)-2.*L2sCoeff(i1,i2,i3)+L2sCoeff(i1-1,i2,i3))/(SQR(drL[0])));
    L2ICoeffrr(i1,i2,i3)  = ((L2ICoeff(i1+1,i2,i3)-2.*L2ICoeff(i1,i2,i3)+L2ICoeff(i1-1,i2,i3))/(SQR(drL[0])));
  // Compute tangential derivatives of L3
    L3rrrCoeffr(i1,i2,i3) = ((L3rrrCoeff(i1+1,i2,i3)-L3rrrCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L3rrsCoeffr(i1,i2,i3) = ((L3rrsCoeff(i1+1,i2,i3)-L3rrsCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L3rssCoeffr(i1,i2,i3) = ((L3rssCoeff(i1+1,i2,i3)-L3rssCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L3sssCoeffr(i1,i2,i3) = ((L3sssCoeff(i1+1,i2,i3)-L3sssCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L3rrCoeffr(i1,i2,i3)  = ((L3rrCoeff(i1+1,i2,i3)-L3rrCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L3rsCoeffr(i1,i2,i3)  = ((L3rsCoeff(i1+1,i2,i3)-L3rsCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L3ssCoeffr(i1,i2,i3)  = ((L3ssCoeff(i1+1,i2,i3)-L3ssCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L3rCoeffr(i1,i2,i3)   = ((L3rCoeff(i1+1,i2,i3)-L3rCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L3sCoeffr(i1,i2,i3)   = ((L3sCoeff(i1+1,i2,i3)-L3sCoeff(i1-1,i2,i3))/(2.*drL[0]));
    L3ICoeffr(i1,i2,i3)   = ((L3ICoeff(i1+1,i2,i3)-L3ICoeff(i1-1,i2,i3))/(2.*drL[0]));
  // Compute L4 coefficients
    L4rrrrCoeff(i1,i2,i3) = (pow(DLR,2.)*pow(crrL(i1,i2,i3),2.)-2.*crrR(j1,j2,j3)*crsR(j1,j2,j3)*L1rCoeff(i1,i2,i3)-2.*crrR(j1,j2,j3)*cssR(j1,j2,j3)*L2rrCoeff(i1,i2,i3)-pow(crsR(j1,j2,j3),2.)*L2rrCoeff(i1,i2,i3)-2.*crsR(j1,j2,j3)*cssR(j1,j2,j3)*L3rrrCoeff(i1,i2,i3)-pow(crrR(j1,j2,j3),2.))/pow(cssR(j1,j2,j3),2.);
    L4rrrsCoeff(i1,i2,i3) = (2.*pow(DLR,2.)*crrL(i1,i2,i3)*crsL(i1,i2,i3)-2.*crrR(j1,j2,j3)*crsR(j1,j2,j3)*L1sCoeff(i1,i2,i3)-2.*crrR(j1,j2,j3)*cssR(j1,j2,j3)*L2rsCoeff(i1,i2,i3)-pow(crsR(j1,j2,j3),2.)*L2rsCoeff(i1,i2,i3)-2.*crsR(j1,j2,j3)*cssR(j1,j2,j3)*L3rrsCoeff(i1,i2,i3))/pow(cssR(j1,j2,j3),2.);
    L4rrssCoeff(i1,i2,i3) = (2.*pow(DLR,2.)*crrL(i1,i2,i3)*cssL(i1,i2,i3)+pow(DLR,2.)*pow(crsL(i1,i2,i3),2.)-2.*crrR(j1,j2,j3)*cssR(j1,j2,j3)*L2ssCoeff(i1,i2,i3)-pow(crsR(j1,j2,j3),2.)*L2ssCoeff(i1,i2,i3)-2.*crsR(j1,j2,j3)*cssR(j1,j2,j3)*L3rssCoeff(i1,i2,i3))/pow(cssR(j1,j2,j3),2.);
    L4rsssCoeff(i1,i2,i3) = (2.*pow(DLR,2.)*crsL(i1,i2,i3)*cssL(i1,i2,i3)-2.*crsR(j1,j2,j3)*cssR(j1,j2,j3)*L3sssCoeff(i1,i2,i3))/pow(cssR(j1,j2,j3),2.);
    L4ssssCoeff(i1,i2,i3) = 1.*pow(DLR,2.)*pow(cssL(i1,i2,i3),2.)/pow(cssR(j1,j2,j3),2.);
    L4rrrCoeff(i1,i2,i3) = ((-1.*L2rCoeff(i1,i2,i3)-2.*L2rrCoeffr(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3rrrCoeffr(i1,i2,i3)-2.*L3rrCoeff(i1,i2,i3))*cssR(j1,j2,j3)+(-2.*L1ICoeff(i1,i2,i3)-6.*L1rCoeffr(i1,i2,i3))*crrR(j1,j2,j3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*L2rrCoeff(i1,i2,i3)+(-2.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1rCoeff(i1,i2,i3)-1.*crrRs(j1,j2,j3)-1.*cssRr(j1,j2,j3)*L3rrrCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-2.*L2rCoeff(i1,i2,i3)-4.*L2rrCoeffr(i1,i2,i3))*crrR(j1,j2,j3)+(-2.*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3))*L2rrCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L1rCoeff(i1,i2,i3)+(-2.*csR(j1,j2,j3)-2.*cssRs(j1,j2,j3))*L3rrrCoeff(i1,i2,i3))*cssR(j1,j2,j3)+(-2.*cssRr(j1,j2,j3)*L2rrCoeff(i1,i2,i3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*L1rCoeff(i1,i2,i3)-2.*crR(j1,j2,j3)-2.*crrRr(j1,j2,j3))*crrR(j1,j2,j3)+((2.*crL(i1,i2,i3)+2.*crrLr(i1,i2,i3))*crrL(i1,i2,i3)+crrLs(i1,i2,i3)*crsL(i1,i2,i3))*pow(DLR,2.))/pow(cssR(j1,j2,j3),2.);
    L4rrsCoeff(i1,i2,i3) = ((-2.*L2rsCoeffr(i1,i2,i3)-1.*L2sCoeff(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3rrsCoeffr(i1,i2,i3)-2.*L3rsCoeff(i1,i2,i3))*cssR(j1,j2,j3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*L2rsCoeff(i1,i2,i3)+(-2.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1sCoeff(i1,i2,i3)-6.*crrR(j1,j2,j3)*L1sCoeffr(i1,i2,i3)-1.*cssRr(j1,j2,j3)*L3rrsCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-2.*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3))*L2rsCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L1sCoeff(i1,i2,i3)+(-4.*L2rsCoeffr(i1,i2,i3)-2.*L2sCoeff(i1,i2,i3))*crrR(j1,j2,j3)+(-2.*csR(j1,j2,j3)-2.*cssRs(j1,j2,j3))*L3rrsCoeff(i1,i2,i3))*cssR(j1,j2,j3)-2.*crrR(j1,j2,j3)*cssRr(j1,j2,j3)*L2rsCoeff(i1,i2,i3)+(-2.*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3))*crrR(j1,j2,j3)*L1sCoeff(i1,i2,i3)+((2.*crL(i1,i2,i3)+crrLr(i1,i2,i3)+crsLs(i1,i2,i3))*crsL(i1,i2,i3)+(2.*crsLr(i1,i2,i3)+2.*csL(i1,i2,i3))*crrL(i1,i2,i3)+2.*crrLs(i1,i2,i3)*cssL(i1,i2,i3))*pow(DLR,2.))/pow(cssR(j1,j2,j3),2.);
    L4rssCoeff(i1,i2,i3) = (((-2.*L3rssCoeffr(i1,i2,i3)-2.*L3ssCoeff(i1,i2,i3))*crsR(j1,j2,j3)+(-2.*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3))*L2ssCoeff(i1,i2,i3)+(-2.*csR(j1,j2,j3)-2.*cssRs(j1,j2,j3))*L3rssCoeff(i1,i2,i3)-4.*L2ssCoeffr(i1,i2,i3)*crrR(j1,j2,j3))*cssR(j1,j2,j3)-2.*pow(crsR(j1,j2,j3),2.)*L2ssCoeffr(i1,i2,i3)+((-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*L2ssCoeff(i1,i2,i3)-1.*cssRr(j1,j2,j3)*L3rssCoeff(i1,i2,i3))*crsR(j1,j2,j3)-2.*crrR(j1,j2,j3)*cssRr(j1,j2,j3)*L2ssCoeff(i1,i2,i3)+((crsLr(i1,i2,i3)+cssLs(i1,i2,i3)+2.*csL(i1,i2,i3))*crsL(i1,i2,i3)+(2.*crL(i1,i2,i3)+2.*crsLs(i1,i2,i3))*cssL(i1,i2,i3)+2.*crrL(i1,i2,i3)*cssLr(i1,i2,i3))*pow(DLR,2.))/pow(cssR(j1,j2,j3),2.);
    L4sssCoeff(i1,i2,i3) = (((-2.*csR(j1,j2,j3)-2.*cssRs(j1,j2,j3))*L3sssCoeff(i1,i2,i3)-2.*crsR(j1,j2,j3)*L3sssCoeffr(i1,i2,i3))*cssR(j1,j2,j3)-1.*crsR(j1,j2,j3)*cssRr(j1,j2,j3)*L3sssCoeff(i1,i2,i3)+((2.*csL(i1,i2,i3)+2.*cssLs(i1,i2,i3))*cssL(i1,i2,i3)+crsL(i1,i2,i3)*cssLr(i1,i2,i3))*pow(DLR,2.))/pow(cssR(j1,j2,j3),2.);
    L4rrCoeff(i1,i2,i3) = ((-1.*L2ICoeff(i1,i2,i3)-1.*L2rrCoeffrr(i1,i2,i3)-2.*L2rCoeffr(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3rCoeff(i1,i2,i3)-2.*L3rrCoeffr(i1,i2,i3))*cssR(j1,j2,j3)+(-6.*L1rCoeffrr(i1,i2,i3)-6.*L1ICoeffr(i1,i2,i3))*crrR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1rCoeff(i1,i2,i3)+(-2.*L2rrCoeffr(i1,i2,i3)-2.*L2rCoeff(i1,i2,i3))*csR(j1,j2,j3)+(-2.*L1ICoeff(i1,i2,i3)-4.*L1rCoeffr(i1,i2,i3))*crR(j1,j2,j3)+(-1.*cssRrs(j1,j2,j3)-1.*csRr(j1,j2,j3))*L2rrCoeff(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3))*L2rrCoeffr(i1,i2,i3)+(-2.*crrRr(j1,j2,j3)-2.*crsRs(j1,j2,j3))*L1rCoeffr(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1ICoeff(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3))*L2rCoeff(i1,i2,i3)-1.*cssRr(j1,j2,j3)*L3rrCoeff(i1,i2,i3)-1.*crrRrs(j1,j2,j3)-1.*crRs(j1,j2,j3))*crsR(j1,j2,j3)+((-2.*L2ICoeff(i1,i2,i3)-2.*L2rrCoeffrr(i1,i2,i3)-4.*L2rCoeffr(i1,i2,i3))*crrR(j1,j2,j3)+(-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*L1rCoeff(i1,i2,i3)-2.*csR(j1,j2,j3)*L3rrCoeff(i1,i2,i3)+(-2.*L2rrCoeffr(i1,i2,i3)-2.*L2rCoeff(i1,i2,i3))*crR(j1,j2,j3)+(-1.*cssRss(j1,j2,j3)-2.*csRs(j1,j2,j3))*L2rrCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L1ICoeff(i1,i2,i3)-2.*cssRs(j1,j2,j3)*L3rrCoeff(i1,i2,i3)-4.*crrRs(j1,j2,j3)*L1rCoeffr(i1,i2,i3)-1.*crrRss(j1,j2,j3)-2.*crsRs(j1,j2,j3)*L2rrCoeffr(i1,i2,i3)-2.*crsRs(j1,j2,j3)*L2rCoeff(i1,i2,i3))*cssR(j1,j2,j3)+((-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1rCoeff(i1,i2,i3)+(-2.*L1ICoeff(i1,i2,i3)-4.*L1rCoeffr(i1,i2,i3))*csR(j1,j2,j3)-2.*crsRr(j1,j2,j3)*L1ICoeff(i1,i2,i3)-4.*crsRr(j1,j2,j3)*L1rCoeffr(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L2rrCoeffr(i1,i2,i3)-1.*crrRrr(j1,j2,j3)-2.*crRr(j1,j2,j3)-1.*cssRrr(j1,j2,j3)*L2rrCoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L2rCoeff(i1,i2,i3))*crrR(j1,j2,j3)+((-2.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*csR(j1,j2,j3)-1.*crsRr(j1,j2,j3)*crR(j1,j2,j3))*L1rCoeff(i1,i2,i3)-1.*pow(csR(j1,j2,j3),2.)*L2rrCoeff(i1,i2,i3)+(-1.*crrRs(j1,j2,j3)-1.*cssRs(j1,j2,j3)*L2rrCoeff(i1,i2,i3))*csR(j1,j2,j3)-1.*pow(crR(j1,j2,j3),2.)+(-1.*crrRr(j1,j2,j3)-1.*cssRr(j1,j2,j3)*L2rrCoeff(i1,i2,i3))*crR(j1,j2,j3)+((crrLrr(i1,i2,i3)+2.*crLr(i1,i2,i3))*crrL(i1,i2,i3)+(crrLrs(i1,i2,i3)+crLs(i1,i2,i3))*crsL(i1,i2,i3)+cssL(i1,i2,i3)*crrLss(i1,i2,i3)+crrLr(i1,i2,i3)*crL(i1,i2,i3)+csL(i1,i2,i3)*crrLs(i1,i2,i3)+pow(crL(i1,i2,i3),2.))*pow(DLR,2.))/pow(cssR(j1,j2,j3),2.);
    L4rsCoeff(i1,i2,i3) = ((-1.*L2rsCoeffrr(i1,i2,i3)-2.*L2sCoeffr(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3rsCoeffr(i1,i2,i3)-2.*L3sCoeff(i1,i2,i3))*cssR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1sCoeff(i1,i2,i3)-6.*crrR(j1,j2,j3)*L1sCoeffrr(i1,i2,i3)+(-1.*cssRrs(j1,j2,j3)-1.*csRr(j1,j2,j3))*L2rsCoeff(i1,i2,i3)+(-2.*L2sCoeff(i1,i2,i3)-2.*L2rsCoeffr(i1,i2,i3))*csR(j1,j2,j3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3))*L2rsCoeffr(i1,i2,i3)+(-2.*crsRs(j1,j2,j3)-2.*crrRr(j1,j2,j3)-4.*crR(j1,j2,j3))*L1sCoeffr(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3))*L2sCoeff(i1,i2,i3)-1.*cssRr(j1,j2,j3)*L3rsCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*L1sCoeff(i1,i2,i3)+(-2.*L2rsCoeffrr(i1,i2,i3)-4.*L2sCoeffr(i1,i2,i3))*crrR(j1,j2,j3)+(-1.*cssRss(j1,j2,j3)-2.*csRs(j1,j2,j3))*L2rsCoeff(i1,i2,i3)-2.*csR(j1,j2,j3)*L3rsCoeff(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3))*L2rsCoeffr(i1,i2,i3)-4.*crrRs(j1,j2,j3)*L1sCoeffr(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3))*L2sCoeff(i1,i2,i3)-2.*cssRs(j1,j2,j3)*L3rsCoeff(i1,i2,i3))*cssR(j1,j2,j3)+((-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*crrR(j1,j2,j3)+(-2.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*csR(j1,j2,j3)-1.*crsRr(j1,j2,j3)*crR(j1,j2,j3))*L1sCoeff(i1,i2,i3)+(-4.*csR(j1,j2,j3)*L1sCoeffr(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L2sCoeff(i1,i2,i3)-4.*crsRr(j1,j2,j3)*L1sCoeffr(i1,i2,i3)-1.*cssRrr(j1,j2,j3)*L2rsCoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L2rsCoeffr(i1,i2,i3))*crrR(j1,j2,j3)+((crsLrs(i1,i2,i3)+crLr(i1,i2,i3)+csLs(i1,i2,i3))*crsL(i1,i2,i3)+(crsLrr(i1,i2,i3)+2.*csLr(i1,i2,i3))*crrL(i1,i2,i3)+(crsLss(i1,i2,i3)+2.*crLs(i1,i2,i3))*cssL(i1,i2,i3)+(crsLr(i1,i2,i3)+2.*csL(i1,i2,i3))*crL(i1,i2,i3)+csL(i1,i2,i3)*crsLs(i1,i2,i3))*pow(DLR,2.)+(-1.*pow(csR(j1,j2,j3),2.)-1.*cssRr(j1,j2,j3)*crR(j1,j2,j3)-1.*csR(j1,j2,j3)*cssRs(j1,j2,j3))*L2rsCoeff(i1,i2,i3))/pow(cssR(j1,j2,j3),2.);
    L4ssCoeff(i1,i2,i3) = (((-1.*cssRss(j1,j2,j3)-2.*csRs(j1,j2,j3))*L2ssCoeff(i1,i2,i3)-2.*crsR(j1,j2,j3)*L3ssCoeffr(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3))*L2ssCoeffr(i1,i2,i3)-2.*cssRs(j1,j2,j3)*L3ssCoeff(i1,i2,i3)-2.*csR(j1,j2,j3)*L3ssCoeff(i1,i2,i3)-2.*crrR(j1,j2,j3)*L2ssCoeffrr(i1,i2,i3))*cssR(j1,j2,j3)+((-1.*cssRrs(j1,j2,j3)-1.*csRr(j1,j2,j3))*crsR(j1,j2,j3)-1.*pow(csR(j1,j2,j3),2.)-1.*crrR(j1,j2,j3)*cssRrr(j1,j2,j3)-1.*csR(j1,j2,j3)*cssRs(j1,j2,j3)-1.*cssRr(j1,j2,j3)*crR(j1,j2,j3))*L2ssCoeff(i1,i2,i3)-1.*pow(crsR(j1,j2,j3),2.)*L2ssCoeffrr(i1,i2,i3)+((-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3)-2.*csR(j1,j2,j3))*L2ssCoeffr(i1,i2,i3)-1.*cssRr(j1,j2,j3)*L3ssCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((cssLrs(i1,i2,i3)+csLr(i1,i2,i3))*crsL(i1,i2,i3)+(cssLss(i1,i2,i3)+2.*csLs(i1,i2,i3))*cssL(i1,i2,i3)+cssLr(i1,i2,i3)*crL(i1,i2,i3)+csL(i1,i2,i3)*cssLs(i1,i2,i3)+crrL(i1,i2,i3)*cssLrr(i1,i2,i3)+pow(csL(i1,i2,i3),2.))*pow(DLR,2.)-2.*crrR(j1,j2,j3)*cssRr(j1,j2,j3)*L2ssCoeffr(i1,i2,i3))/pow(cssR(j1,j2,j3),2.);
    L4rCoeff(i1,i2,i3) = ((-1.*L2rCoeffrr(i1,i2,i3)-2.*L2ICoeffr(i1,i2,i3))*pow(crsR(j1,j2,j3),2.)+((-2.*L3rCoeffr(i1,i2,i3)-2.*L3ICoeff(i1,i2,i3))*cssR(j1,j2,j3)+(-6.*L1ICoeffrr(i1,i2,i3)-2.*L1rCoeffrrr(i1,i2,i3))*crrR(j1,j2,j3)+(-2.*L2ICoeff(i1,i2,i3)-2.*L2rCoeffr(i1,i2,i3))*csR(j1,j2,j3)+(-2.*L1rCoeffrr(i1,i2,i3)-4.*L1ICoeffr(i1,i2,i3))*crR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1rCoeffr(i1,i2,i3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1ICoeff(i1,i2,i3)+(-1.*cssRrs(j1,j2,j3)-1.*csRr(j1,j2,j3))*L2rCoeff(i1,i2,i3)+(-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1rCoeffrr(i1,i2,i3)+(-2.*crrRr(j1,j2,j3)-2.*crsRs(j1,j2,j3))*L1ICoeffr(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3))*L2rCoeffr(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3))*L2ICoeff(i1,i2,i3)-1.*crRrs(j1,j2,j3)-1.*csRrs(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*cssRr(j1,j2,j3)*L3rCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-4.*L2ICoeffr(i1,i2,i3)-2.*L2rCoeffrr(i1,i2,i3))*crrR(j1,j2,j3)-2.*csR(j1,j2,j3)*L3rCoeff(i1,i2,i3)+(-2.*L2ICoeff(i1,i2,i3)-2.*L2rCoeffr(i1,i2,i3))*crR(j1,j2,j3)+(-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*L1rCoeffr(i1,i2,i3)+(-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*L1ICoeff(i1,i2,i3)+(-1.*cssRss(j1,j2,j3)-2.*csRs(j1,j2,j3))*L2rCoeff(i1,i2,i3)-1.*crRss(j1,j2,j3)-2.*crsRs(j1,j2,j3)*L2ICoeff(i1,i2,i3)-2.*crsRs(j1,j2,j3)*L2rCoeffr(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L1rCoeffrr(i1,i2,i3)-4.*crrRs(j1,j2,j3)*L1ICoeffr(i1,i2,i3)-2.*cssRs(j1,j2,j3)*L3rCoeff(i1,i2,i3)-1.*csRss(j1,j2,j3)*L1rCoeff(i1,i2,i3))*cssR(j1,j2,j3)+((-2.*L1rCoeffrr(i1,i2,i3)-4.*L1ICoeffr(i1,i2,i3))*csR(j1,j2,j3)+(-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1rCoeffr(i1,i2,i3)+(-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*L1ICoeff(i1,i2,i3)-1.*crRrr(j1,j2,j3)-2.*cssRr(j1,j2,j3)*L2ICoeff(i1,i2,i3)-2.*crsRr(j1,j2,j3)*L1rCoeffrr(i1,i2,i3)-4.*crsRr(j1,j2,j3)*L1ICoeffr(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L2rCoeffr(i1,i2,i3)-1.*csRrr(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*cssRrr(j1,j2,j3)*L2rCoeff(i1,i2,i3))*crrR(j1,j2,j3)-1.*pow(csR(j1,j2,j3),2.)*L2rCoeff(i1,i2,i3)+((-2.*L1rCoeffr(i1,i2,i3)-2.*L1ICoeff(i1,i2,i3))*crR(j1,j2,j3)-1.*crRs(j1,j2,j3)-1.*crsRs(j1,j2,j3)*L1rCoeffr(i1,i2,i3)-1.*cssRs(j1,j2,j3)*L2rCoeff(i1,i2,i3)-1.*csRs(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*crsRs(j1,j2,j3)*L1ICoeff(i1,i2,i3))*csR(j1,j2,j3)+(-1.*crRr(j1,j2,j3)-1.*crsRr(j1,j2,j3)*L1rCoeffr(i1,i2,i3)-1.*cssRr(j1,j2,j3)*L2rCoeff(i1,i2,i3)-1.*csRr(j1,j2,j3)*L1rCoeff(i1,i2,i3)-1.*crsRr(j1,j2,j3)*L1ICoeff(i1,i2,i3))*crR(j1,j2,j3)+pow(DLR,2.)*(crL(i1,i2,i3)*crLr(i1,i2,i3)+crLrr(i1,i2,i3)*crrL(i1,i2,i3)+crLrs(i1,i2,i3)*crsL(i1,i2,i3)+crLs(i1,i2,i3)*csL(i1,i2,i3)+crLss(i1,i2,i3)*cssL(i1,i2,i3)))/pow(cssR(j1,j2,j3),2.);
    L4sCoeff(i1,i2,i3) = (-1.*pow(crsR(j1,j2,j3),2.)*L2sCoeffrr(i1,i2,i3)+(-2.*L3sCoeffr(i1,i2,i3)*cssR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1sCoeffr(i1,i2,i3)-2.*crrR(j1,j2,j3)*L1sCoeffrrr(i1,i2,i3)+(-1.*cssRrs(j1,j2,j3)-1.*csRr(j1,j2,j3))*L2sCoeff(i1,i2,i3)-2.*csR(j1,j2,j3)*L2sCoeffr(i1,i2,i3)+(-2.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1sCoeffrr(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3))*L2sCoeffr(i1,i2,i3)-1.*cssRr(j1,j2,j3)*L3sCoeff(i1,i2,i3)-1.*csRrs(j1,j2,j3)*L1sCoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*L1sCoeffr(i1,i2,i3)-2.*crrR(j1,j2,j3)*L2sCoeffrr(i1,i2,i3)+(-1.*cssRss(j1,j2,j3)-2.*csRs(j1,j2,j3))*L2sCoeff(i1,i2,i3)-2.*csR(j1,j2,j3)*L3sCoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L1sCoeffrr(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3))*L2sCoeffr(i1,i2,i3)-2.*cssRs(j1,j2,j3)*L3sCoeff(i1,i2,i3)-1.*csRss(j1,j2,j3)*L1sCoeff(i1,i2,i3))*cssR(j1,j2,j3)+((-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*crrR(j1,j2,j3)+(-2.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*csR(j1,j2,j3)-1.*crsRr(j1,j2,j3)*crR(j1,j2,j3))*L1sCoeffr(i1,i2,i3)+(-1.*cssRrr(j1,j2,j3)*L2sCoeff(i1,i2,i3)-2.*crsRr(j1,j2,j3)*L1sCoeffrr(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L2sCoeffr(i1,i2,i3)-2.*csR(j1,j2,j3)*L1sCoeffrr(i1,i2,i3)-1.*csRrr(j1,j2,j3)*L1sCoeff(i1,i2,i3))*crrR(j1,j2,j3)+(-1.*pow(csR(j1,j2,j3),2.)-1.*cssRr(j1,j2,j3)*crR(j1,j2,j3)-1.*csR(j1,j2,j3)*cssRs(j1,j2,j3))*L2sCoeff(i1,i2,i3)-1.*csR(j1,j2,j3)*csRs(j1,j2,j3)*L1sCoeff(i1,i2,i3)-1.*crR(j1,j2,j3)*csRr(j1,j2,j3)*L1sCoeff(i1,i2,i3)+pow(DLR,2.)*(crL(i1,i2,i3)*csLr(i1,i2,i3)+crrL(i1,i2,i3)*csLrr(i1,i2,i3)+crsL(i1,i2,i3)*csLrs(i1,i2,i3)+csL(i1,i2,i3)*csLs(i1,i2,i3)+csLss(i1,i2,i3)*cssL(i1,i2,i3)))/pow(cssR(j1,j2,j3),2.);
    L4ICoeff(i1,i2,i3) = (-1.*pow(crsR(j1,j2,j3),2.)*L2ICoeffrr(i1,i2,i3)+(-2.*L3ICoeffr(i1,i2,i3)*cssR(j1,j2,j3)+(-1.*crsRrs(j1,j2,j3)-1.*crRr(j1,j2,j3)-1.*csRs(j1,j2,j3))*L1ICoeffr(i1,i2,i3)-2.*crrR(j1,j2,j3)*L1ICoeffrrr(i1,i2,i3)+(-1.*cssRrs(j1,j2,j3)-1.*csRr(j1,j2,j3))*L2ICoeff(i1,i2,i3)-2.*csR(j1,j2,j3)*L2ICoeffr(i1,i2,i3)+(-2.*crR(j1,j2,j3)-1.*crrRr(j1,j2,j3)-1.*crsRs(j1,j2,j3))*L1ICoeffrr(i1,i2,i3)+(-1.*cssRs(j1,j2,j3)-1.*crsRr(j1,j2,j3))*L2ICoeffr(i1,i2,i3)-1.*csRrs(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*cssRr(j1,j2,j3)*L3ICoeff(i1,i2,i3))*crsR(j1,j2,j3)+((-1.*crsRss(j1,j2,j3)-2.*crRs(j1,j2,j3))*L1ICoeffr(i1,i2,i3)-2.*crrR(j1,j2,j3)*L2ICoeffrr(i1,i2,i3)+(-1.*cssRss(j1,j2,j3)-2.*csRs(j1,j2,j3))*L2ICoeff(i1,i2,i3)-2.*csR(j1,j2,j3)*L3ICoeff(i1,i2,i3)-2.*crrRs(j1,j2,j3)*L1ICoeffrr(i1,i2,i3)+(-2.*crR(j1,j2,j3)-2.*crsRs(j1,j2,j3))*L2ICoeffr(i1,i2,i3)-2.*cssRs(j1,j2,j3)*L3ICoeff(i1,i2,i3)-1.*csRss(j1,j2,j3)*L1ICoeff(i1,i2,i3))*cssR(j1,j2,j3)+((-1.*crsRrr(j1,j2,j3)-2.*csRr(j1,j2,j3))*crrR(j1,j2,j3)+(-2.*crR(j1,j2,j3)-1.*crsRs(j1,j2,j3))*csR(j1,j2,j3)-1.*crsRr(j1,j2,j3)*crR(j1,j2,j3))*L1ICoeffr(i1,i2,i3)+(-1.*csRrr(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*cssRrr(j1,j2,j3)*L2ICoeff(i1,i2,i3)-2.*cssRr(j1,j2,j3)*L2ICoeffr(i1,i2,i3)-2.*crsRr(j1,j2,j3)*L1ICoeffrr(i1,i2,i3)-2.*csR(j1,j2,j3)*L1ICoeffrr(i1,i2,i3))*crrR(j1,j2,j3)+(-1.*pow(csR(j1,j2,j3),2.)-1.*cssRr(j1,j2,j3)*crR(j1,j2,j3)-1.*csR(j1,j2,j3)*cssRs(j1,j2,j3))*L2ICoeff(i1,i2,i3)-1.*crR(j1,j2,j3)*csRr(j1,j2,j3)*L1ICoeff(i1,i2,i3)-1.*csR(j1,j2,j3)*csRs(j1,j2,j3)*L1ICoeff(i1,i2,i3))/pow(cssR(j1,j2,j3),2.);
} // end FOR_3D

// ----- Define coefficients in difference operators -----
Range R4(-2,2);
RealArray iCoeff(R4);
RealArray rCoeff(R4), rrCoeff(R4), rrrCoeff(R4), rrrrCoeff(R4);
RealArray sCoeff(R4), ssCoeff(R4), sssCoeff(R4), ssssCoeff(R4);
      iCoeff(-2)= 0.;    iCoeff(-1)=  0.;    iCoeff(0)=  1.;    iCoeff(1)= 0.;    iCoeff(2)= 0.;
      rCoeff(-2)= 1.;    rCoeff(-1)= -8.;    rCoeff(0)=  0.;    rCoeff(1)= 8.;    rCoeff(2)=-1.; rCoeff    /=(    12.*drL[0]);
    rrCoeff(-2)=-1.;   rrCoeff(-1)= 16.;   rrCoeff(0)=-30.;   rrCoeff(1)=16.;   rrCoeff(2)=-1.; rrCoeff   /=(12.*SQR(drL[0]));
  rrrCoeff(-2)=-1.;  rrrCoeff(-1)=  2.;  rrrCoeff(0)=  0.;  rrrCoeff(1)=-2.;  rrrCoeff(2)= 1.; rrrCoeff  /=( 2.*pow(drL[0],3));
rrrrCoeff(-2)= 1.; rrrrCoeff(-1)= -4.; rrrrCoeff(0)=  6.; rrrrCoeff(1)=-4.; rrrrCoeff(2)= 1.; rrrrCoeff /=(    pow(drL[0],4));
      sCoeff(-2)= 1.;    sCoeff(-1)= -8.;    sCoeff(0)=  0.;    sCoeff(1)= 8.;    sCoeff(2)=-1.; sCoeff    /=(    12.*drL[1]);
    ssCoeff(-2)=-1.;   ssCoeff(-1)= 16.;   ssCoeff(0)=-30.;   ssCoeff(1)=16.;   ssCoeff(2)=-1.; ssCoeff   /=(12.*SQR(drL[1]));
  sssCoeff(-2)=-1.;  sssCoeff(-1)=  2.;  sssCoeff(0)=  0.;  sssCoeff(1)=-2.;  sssCoeff(2)= 1.; sssCoeff  /=( 2.*pow(drL[1],3));
ssssCoeff(-2)= 1.; ssssCoeff(-1)= -4.; ssssCoeff(0)=  6.; ssssCoeff(1)=-4.; ssssCoeff(2)= 1.; ssssCoeff /=(    pow(drL[1],4));

// ----- Fill in the Matrix Coefficients for CHAMP -----
const Real dxs = (1-2*side2)*drR[axis2];
const Real h = dxs;
printF("WDH: grid=%d, (side,axis)=(%d,%d) (side2,axis2)=(%d,%d) dxs=%9.3e, Sl=%9.3e\n",grid, side,axis, side2,axis2, dxs,Sl);
RealArray coeff4(R4,R4);
const int e=0, c=0; // eqn number and component number
// ------- FILL CHAMP conditions into the matrix -----
// NOTE: skip top boundary if periodic (use indexRange) 
// NOTE: skip adjacent boundaries if Dirichlet BC *finish me**
extra=0; 
getBoundaryIndex(mg.indexRange(), side, axis, Ib1,Ib2,Ib3,extra);
getBoundaryIndex(mg2.indexRange(),side2,axis2,Jb1,Jb2,Jb3,extra);
int axisp = (axis+1) % numberOfDimensions;
FOR_3IJ(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3)
{
    i1m=i1-is1, i2m=i2-is2, i3m=i3-is3; //  ghost point is (i1m,i2m,i3m)
    const Real bn = b2R(j1,j2,j3), bt = b1R(j1,j2,j3);  // b in normal and tangential directions
  // The next macro defines cI, cr, cs, crr, ...
  // ------  ORDER = 4 , NormalDirection = 1,------
    const Real h2By2=h*h/2., h3By6=h2By2*h/3., h4By24=h3By6*h/4.;
    Real L1IC = L1ICoeff(i1,i2,i3);
    Real L1ICr = L1ICoeffr(i1,i2,i3);
    Real L1ICrr = L1ICoeffrr(i1,i2,i3);
    Real L1ICrrr = L1ICoeffrrr(i1,i2,i3);
    Real L1rC = L1rCoeff(i1,i2,i3);
    Real L1rCr = L1rCoeffr(i1,i2,i3);
    Real L1rCrr = L1rCoeffrr(i1,i2,i3);
    Real L1rCrrr = L1rCoeffrrr(i1,i2,i3);
    Real L1sC = L1sCoeff(i1,i2,i3);
    Real L1sCr = L1sCoeffr(i1,i2,i3);
    Real L1sCrr = L1sCoeffrr(i1,i2,i3);
    Real L1sCrrr = L1sCoeffrrr(i1,i2,i3);
    Real L2IC = L2ICoeff(i1,i2,i3);
    Real L2ICr = L2ICoeffr(i1,i2,i3);
    Real L2ICrr = L2ICoeffrr(i1,i2,i3);
    Real L2rC = L2rCoeff(i1,i2,i3);
    Real L2rCr = L2rCoeffr(i1,i2,i3);
    Real L2rCrr = L2rCoeffrr(i1,i2,i3);
    Real L2sC = L2sCoeff(i1,i2,i3);
    Real L2sCr = L2sCoeffr(i1,i2,i3);
    Real L2sCrr = L2sCoeffrr(i1,i2,i3);
    Real L2rrC = L2rrCoeff(i1,i2,i3);
    Real L2rrCr = L2rrCoeffr(i1,i2,i3);
    Real L2rrCrr = L2rrCoeffrr(i1,i2,i3);
    Real L2rsC = L2rsCoeff(i1,i2,i3);
    Real L2rsCr = L2rsCoeffr(i1,i2,i3);
    Real L2rsCrr = L2rsCoeffrr(i1,i2,i3);
    Real L2ssC = L2ssCoeff(i1,i2,i3);
    Real L2ssCr = L2ssCoeffr(i1,i2,i3);
    Real L2ssCrr = L2ssCoeffrr(i1,i2,i3);
    Real L3IC = L3ICoeff(i1,i2,i3);
    Real L3ICr = L3ICoeffr(i1,i2,i3);
    Real L3rC = L3rCoeff(i1,i2,i3);
    Real L3rCr = L3rCoeffr(i1,i2,i3);
    Real L3sC = L3sCoeff(i1,i2,i3);
    Real L3sCr = L3sCoeffr(i1,i2,i3);
    Real L3rrC = L3rrCoeff(i1,i2,i3);
    Real L3rrCr = L3rrCoeffr(i1,i2,i3);
    Real L3rsC = L3rsCoeff(i1,i2,i3);
    Real L3rsCr = L3rsCoeffr(i1,i2,i3);
    Real L3ssC = L3ssCoeff(i1,i2,i3);
    Real L3ssCr = L3ssCoeffr(i1,i2,i3);
    Real L3rrrC = L3rrrCoeff(i1,i2,i3);
    Real L3rrrCr = L3rrrCoeffr(i1,i2,i3);
    Real L3rrsC = L3rrsCoeff(i1,i2,i3);
    Real L3rrsCr = L3rrsCoeffr(i1,i2,i3);
    Real L3rssC = L3rssCoeff(i1,i2,i3);
    Real L3rssCr = L3rssCoeffr(i1,i2,i3);
    Real L3sssC = L3sssCoeff(i1,i2,i3);
    Real L3sssCr = L3sssCoeffr(i1,i2,i3);
    Real L4IC = L4ICoeff(i1,i2,i3);
    Real L4rC = L4rCoeff(i1,i2,i3);
    Real L4sC = L4sCoeff(i1,i2,i3);
    Real L4rrC = L4rrCoeff(i1,i2,i3);
    Real L4rsC = L4rsCoeff(i1,i2,i3);
    Real L4ssC = L4ssCoeff(i1,i2,i3);
    Real L4rrrC = L4rrrCoeff(i1,i2,i3);
    Real L4rrsC = L4rrsCoeff(i1,i2,i3);
    Real L4rssC = L4rssCoeff(i1,i2,i3);
    Real L4sssC = L4sssCoeff(i1,i2,i3);
    Real L4rrrrC = L4rrrrCoeff(i1,i2,i3);
    Real L4rrrsC = L4rrrsCoeff(i1,i2,i3);
    Real L4rrssC = L4rrssCoeff(i1,i2,i3);
    Real L4rsssC = L4rsssCoeff(i1,i2,i3);
    Real L4ssssC = L4ssssCoeff(i1,i2,i3);
    Real cI  = Sl*(h*L1IC+h2By2*L2IC+h3By6*L3IC+h4By24*L4IC+1)+(-h*L2IC-h2By2*L3IC-h3By6*L4IC-L1IC)*bn-bt*(h*L1ICr+h2By2*L2ICr+h3By6*L3ICr);
    Real cr  = ((-L2IC-L2rCr)*h2By2+(-L3IC-L3rCr)*h3By6+(-L1rCr-L1IC)*h-L4IC*h4By24-2)*bt+(Sl*L2rC-bn*L3rC)*h2By2+(Sl*L3rC-bn*L4rC)*h3By6+(h*L1rC+h4By24*L4rC)*Sl-bn*(h*L2rC+L1rC);
    Real cs  = Sl*(h*L1sC+h2By2*L2sC+h3By6*L3sC+h4By24*L4sC)+(-h*L2sC-h2By2*L3sC-h3By6*L4sC-L1sC)*bn-bt*(h*L1sCr+h2By2*L2sCr+h3By6*L3sCr);
    Real crr  = ((-L2rrCr-L2rC)*h2By2+(-L3rC-L3rrCr)*h3By6-L1rC*h-L4rC*h4By24)*bt+(Sl*L2rrC-bn*L3rrC)*h2By2+(Sl*L3rrC-bn*L4rrC)*h3By6-bn*L2rrC*h+Sl*h4By24*L4rrC;
    Real crs  = ((-L2rsCr-L2sC)*h2By2+(-L3sC-L3rsCr)*h3By6-L1sC*h-L4sC*h4By24)*bt+(Sl*L2rsC-bn*L3rsC)*h2By2+(Sl*L3rsC-bn*L4rsC)*h3By6-bn*L2rsC*h+Sl*h4By24*L4rsC;
    Real css  = (Sl*L2ssC-bn*L3ssC-bt*L2ssCr)*h2By2+(Sl*L3ssC-bn*L4ssC-bt*L3ssCr)*h3By6-bn*L2ssC*h+Sl*h4By24*L4ssC;
    Real crrr  = ((-L3rrrCr-L3rrC)*bt+Sl*L3rrrC-bn*L4rrrC)*h3By6+(-h2By2*L2rrC-h4By24*L4rrC)*bt+Sl*h4By24*L4rrrC-bn*h2By2*L3rrrC;
    Real crrs  = ((-L3rrsCr-L3rsC)*bt+Sl*L3rrsC-L4rrsC*bn)*h3By6+(-h2By2*L2rsC-h4By24*L4rsC)*bt+Sl*h4By24*L4rrsC-bn*h2By2*L3rrsC;
    Real crss  = ((-L3rssCr-L3ssC)*bt+Sl*L3rssC-L4rssC*bn)*h3By6+(-h2By2*L2ssC-h4By24*L4ssC)*bt+Sl*h4By24*L4rssC-bn*h2By2*L3rssC;
    Real csss  = (Sl*L3sssC-bn*L4sssC-bt*L3sssCr)*h3By6+Sl*h4By24*L4sssC-bn*h2By2*L3sssC;
    Real crrrr  = (-h3By6*L3rrrC-h4By24*L4rrrC)*bt+L4rrrrC*(Sl*h4By24-bn*h3By6);
    Real crrrs  = (-h3By6*L3rrsC-h4By24*L4rrsC)*bt+(Sl*h4By24-bn*h3By6)*L4rrrsC;
    Real crrss  = (-h3By6*L3rssC-h4By24*L4rssC)*bt+(Sl*h4By24-bn*h3By6)*L4rrssC;
    Real crsss  = (-h3By6*L3sssC-h4By24*L4sssC)*bt+(Sl*h4By24-bn*h3By6)*L4rsssC;
    Real cssss  = (Sl*h4By24-bn*h3By6)*L4ssssC;
//  if( orderOfAccuracy==4 )
//  {
//    printF("domain2=%d: cI   =%12.4e, cr   =%12.4e, cs   =%12.4e, crr  =%12.4e, crs=%12.4e, css=%12.4e\n",domain2,cI,cr,cs,crr,crs,css);
//    printF("            crrr =%12.4e, crrs =%12.4e, crss =%12.4e, csss =%12.4e\n",crrr,crrs,crss,csss);
//    printF("            crrrr=%12.4e, crrrs=%12.4e, crrss=%12.4e, crsss=%12.4e, cssss=%12.4e\n",crrrr,crrrs,crrss,crsss,cssss);
//    printF("            L4ssssC=%12.4e, (Sl*h4By24-bn*h3By6)=%12.4e h4By24=%12.4e h3By6=%12.4e, Sl=%12.4e, bn=%12.4e, h=%12.4e, bt=%12.4e\n",L4ssssC,(Sl*h4By24-bn*h3By6),h4By24,h3By6,Sl,bn,h,bt);
//    if( domain2==0 )
//    {
//      OV_ABORT("stop here for now");
//    }
//  }
    ForStencil(m1,m2,m3)
    {
        int m  = M123(m1,m2,m3);        // the single-component coeff-index
        int mm = M123CE(m1,m2,m3,c,e);  // the system coeff-index
        coeff4(m1,m2) = 
                          + cI   *   iCoeff(m1) * iCoeff(m2) 
                          + cr   *   rCoeff(m1) * iCoeff(m2) 
                          + cs   *   iCoeff(m1) * sCoeff(m2) 
                          + crr  *  rrCoeff(m1) * iCoeff(m2) 
                          + crs  *   rCoeff(m1) * sCoeff(m2) 
                          + css  *   iCoeff(m1) *ssCoeff(m2) 
                          + crrr   * rrrCoeff(m1)  *   iCoeff(m2) 
                          + crrs   *  rrCoeff(m1)  *   sCoeff(m2) 
                          + crss   *   rCoeff(m1)  *  ssCoeff(m2) 
                          + csss   *   iCoeff(m1)  * sssCoeff(m2) 
                          + crrrr  *rrrrCoeff(m1)  *   iCoeff(m2) 
                          + crrrs  * rrrCoeff(m1)  *   sCoeff(m2) 
                          + crrss  *  rrCoeff(m1)  *  ssCoeff(m2) 
                          + crsss  *   rCoeff(m1)  * sssCoeff(m2) 
                          + cssss  *   iCoeff(m1)  *ssssCoeff(m2) 
                          ;
        if( fillMatrixWDH )
        {
            coeff(mm,i1m,i2m,i3m) = coeff4(m1,m2);
      // Specify that the above coeff value is the coefficient of component c at the grid point (j1,j2,j3).
            const int k1=i1+m1, k2=i2+m2, k3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)    
            setEquationNumber(mm, e,i1m,i2m,i3m,  c,k1,k2,k3 );  // macro to set equationNumber                  
                                                                                                                                                                                                                      
       // Fill Ghost 2 -- extrapolation for now                                                            
              const int ghost=2;                                                                                  
                  const int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
         // --- fill in the coefficients of the extrapolation formula ---
                  for( int me=0; me<=extrapOrder; me++ )
                  {
                      coeff(me,i1m,i2m,i3m) = extrapCoeff[me];
                      const int j1=i1m + me*is1, j2=i2m + me*is2, j3=i3m + me*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
                      setEquationNumber(me, e,i1m,i2m,i3m,  c,j1,j2,j3 );             // macro to set equationNumber
                  }                
        }
        if( twilightZoneFlow )
        {
      // --- For twilightZone we save some coefficients that go into the CHAMP matrix ---
            cc(i1,i2,i3, 0) = cI;   // coeff of I
            cc(i1,i2,i3, 1) = cr;   // coeff of ur
            cc(i1,i2,i3, 2) = cs;   // coeff of us
            cc(i1,i2,i3, 3) = crr;  // coeff of urr
            cc(i1,i2,i3, 4) = crs;  // coeff of urs
            cc(i1,i2,i3, 5) = css;  // coeff of uss
            cc(i1,i2,i3, 6) = crrr;    // coeff of urrr
            cc(i1,i2,i3, 7) = crrs;    // coeff of urrs
            cc(i1,i2,i3, 8) = crss;    // coeff of urss
            cc(i1,i2,i3, 9) = csss;    // coeff of usss
            cc(i1,i2,i3,10) = crrrr;   // coeff of urrrr
            cc(i1,i2,i3,11) = crrrs;   // coeff of urrrs
            cc(i1,i2,i3,12) = crrss;   // coeff of urrss
            cc(i1,i2,i3,13) = crsss;   // coeff of ursss
            cc(i1,i2,i3,14) = cssss;   // coeff of ussss
        } 
    } // end ForStencil
  if( debug & 4 ) 
  {
      printF("WDH:order=4: (i1,i2,i3)=(%3d,%3d,%3d)\n",i1,i2,i3);     
      printF("    coeff4=(%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"    
                    "            %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"    
                    "            %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"    
                    "            %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"    
                    "            %10.3e,%10.3e,%10.3e,%10.3e,%10.3e )\n",  
                      coeff4(-2,-2),coeff4(-1,-2),coeff4(0,-2),coeff4(1,-2),coeff4(2,-2),  
                      coeff4(-2,-1),coeff4(-1,-1),coeff4(0,-1),coeff4(1,-1),coeff4(2,-1),  
                      coeff4(-2, 0),coeff4(-1, 0),coeff4(0, 0),coeff4(1, 0),coeff4(2, 0),  
                      coeff4(-2, 1),coeff4(-1, 1),coeff4(0, 1),coeff4(1, 1),coeff4(2, 1),  
                      coeff4(-2, 2),coeff4(-1, 2),coeff4(0, 2),coeff4(1, 2),coeff4(2, 2)); 
  }
} // end FOR_3D
                            }

                        }
                        else
                        {
                            printF("champBC: ERROR: axis=%d and axis2=%d not currently supported\n",axis,axis2);
                            OV_ABORT("ERROR: FINISH ME BILL!");
                        }
            // fillMatrixWDH=false; 
                    }
                    else
                    {
                        OV_ABORT("champBC: ERROR");
                    }

          // FOR_3IJD(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3) // loop over points on BOTH SIDES of the interface 
          // {
                            
          //   i1m=i1-is1, i2m=i2-is2, i3m=i3-is3; //  ghost point is (i1m,i2m,i3m)

          //   // ----- Evaluate variables in the champ conditions -----
          //   // *WDH This next macro is only for order=4 *** FIX ME **************************
          //   evalChampVariables(); 


          //   //printF("SL=%g,\n", Sl);

          //   const Real hR = (1-2*side2)*drR[axis2]/b1R; // scaled grid spacing on right
                    
          //   if( Sl<0 )
          //   {
          //     // Define Sl from the first grid point encountered   *** Do this for now ***
          //     // Sl = pl/fabs(hR);  // note hR may be negative 
          //     Sl = pl/hR;  // note hR may be negative 
          //     champParameters(4,side,axis,grid)=Sl; // save Sl
          //   }

          //   const real a0 = Sl;
          //   const real a1 = 1. + hR*Sl;
          //   const real a2 = hR*(1.0 + .5*hR*Sl);

          //   if( false ) 
          //     printF(" FILL MATRIX BC FOR GHOST PT (i1m,i2m,i3m)=(%d,%d,%d) CURVILINEAR Sl=%g, hR=%g, (heatFluxInterface)\n",i1m,i2m,i3m,Sl,hR);

          //   //#include "testCHAMP4.h"
                        
          //   if( orderOfAccuracy==2 )
          //   {
          //     #include "champ2InterfaceStencilCurvilinear.h"
          //   }
          //   else if( orderOfAccuracy==4 )
          //   {
          //     //#include "champ4InterfaceStencilCurvilinear.h"

          //     Real hI=drR[0];
          //     #Include "champ4InterfaceStencilCurvilinearMacro.h"
          //     evalChamp4StencilCurvilinear(hI,coefA);
          //   }


          //   //Real hI=drr;
          //   //#Include "champ4InterfaceStencilCurvilinearMacro.h"
          //   //evalChamp4StencilCurvilinear(hI,coefA);

          //   //#include "champ2InterfaceStencilCurvilinear.h"

          //   if( debug & 4 && orderOfAccuracy==4 )
          //   {
          //     printF("SH: order=4 coefA=(%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"
          //            "                   %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"
          //            "                   %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"
          //            "                   %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n"
          //            "                   %10.3e,%10.3e,%10.3e,%10.3e,%10.3e )\n",
          //             coefA(-2,-2),coefA(-1,-2),coefA(0,-2),coefA(1,-2),coefA(2,-2),//             coefA(-2,-1),coefA(-1,-1),coefA(0,-1),coefA(1,-1),coefA(2,-1),//             coefA(-2, 0),coefA(-1, 0),coefA(0, 0),coefA(1, 0),coefA(2, 0),//             coefA(-2, 1),coefA(-1, 1),coefA(0, 1),coefA(1, 1),coefA(2, 1),//             coefA(-2, 2),coefA(-1, 2),coefA(0, 2),coefA(1, 2),coefA(2, 2));

          //     printF("===========================================\n");
          //     //#Include "champ2InterfaceStencilCurvilinear.h"
          //     //
          //     //printF("    My coefficients: coefA =(%g,%g,%g,%g,%g,%g,%g,%g,%g)\n", coefA(-1,-1),coefA(0,-1),coefA(1,-1),coefA(-1,0),coefA(0,0),coefA(1,0),coefA(-1,1),coefA(0,1),coefA(1,1));

                            
          //     printF("Domain2= %i:  i1=%i, i2=%i, i3=%i, j1=%i, j2=%i, j3=%i \n",domain2,i1,i2,i3,j1,j2,j3);
          //     printF("-------------------------------------------\n");
          //     printF("    an1L=%g, an2L=%g \n",an1L,an2L);
          //     printF("    b1L=%g, b2L=%g, c11L=%g, c12L=%g, c21L=%g, c22L=%g, c1L=%g, c2L=%g\n",
          //             b1L,b2L,c11L,c12L,c21L,c22L,c1L,c2L); 
          //     printF("    b1R=%g, b2R=%g, c11R=%g, c12R=%g, c21R=%g, c22R=%g, c1R=%g, c2R=%g\n",
          //             b1R,b2R,c11R,c12R,c21R,c22R,c1R,c2R); 
          //     printF("    b1Lr2=%g, b2Lr2=%g, b1Rr2=%g, b2Rr2=%g \n",b1Lr2,b2Lr2,b1Rr2,b2Rr2);
                            
          //     printF("     c2L =[%g,%g,%g,%g,%g;\n           //                        %g,%g,%g,%g,%g;\n           //                       %g,%g,%g,%g,%g;\n           //                       %g,%g,%g,%g,%g;\n           //                       %g,%g,%g,%g,%g];\n",//             c2Lm(i1-2,i2-2,i3),c2Lm(i1-2,i2-1,i3),c2Lm(i1-2,i2,i3),c2Lm(i1-2,i2+1,i3),c2Lm(i1-2,i2+2,i3),//             c2Lm(i1-1,i2-2,i3),c2Lm(i1-1,i2-1,i3),c2Lm(i1-1,i2,i3),c2Lm(i1-1,i2+1,i3),c2Lm(i1-1,i2+2,i3), //             c2Lm(i1  ,i2-2,i3),c2Lm(i1  ,i2-1,i3),c2Lm(i1  ,i2,i3),c2Lm(i1,  i2+1,i3),c2Lm(i1  ,i2+2,i3), //             c2Lm(i1+1,i2-2,i3),c2Lm(i1+1,i2-1,i3),c2Lm(i1+1,i2,i3),c2Lm(i1+1,i2+1,i3),c2Lm(i1+1,i2+2,i3), \         
          //             c2Lm(i1+2,i2-2,i3),c2Lm(i1+2,i2-1,i3),c2Lm(i1+2,i2,i3),c2Lm(i1+2,i2+1,i3),c2Lm(i1+2,i2+2,i3));
                            
          //     printF("     b2L =[%g,%g,%g,%g,%g;\n           //                        %g,%g,%g,%g,%g;\n           //                       %g,%g,%g,%g,%g;\n           //                       %g,%g,%g,%g,%g;\n           //                       %g,%g,%g,%g,%g];\n",//           b2Lm(i1-2,i2-2,i3,j1-2,j2-2,j3),b2Lm(i1-2,i2-1,i3,j1-2,j2-1,j3),b2Lm(i1-2,i2,i3,j1-2,j2,j3),b2Lm(i1-2,i2+1,i3,j1-2,j2+1,j3),b2Lm(i1-2,i2+2,i3,j1-2,j2+2,j3),//           b2Lm(i1-1,i2-2,i3,j1-1,j2-2,j3),b2Lm(i1-1,i2-1,i3,j1-1,j2-1,j3),b2Lm(i1-1,i2,i3,j1-1,j2,j3),b2Lm(i1-1,i2+1,i3,j1-1,j2+1,j3),b2Lm(i1-1,i2+2,i3,j1-1,j2+2,j3), //           b2Lm(i1  ,i2-2,i3,j1  ,j2-2,j3),b2Lm(i1  ,i2-1,i3,j1  ,j2-1,j3),b2Lm(i1  ,i2,i3,j1  ,j2,j3),b2Lm(i1,  i2+1,i3,j1  ,j2+1,j3),b2Lm(i1  ,i2+2,i3,j1  ,j2+2,j3), //           b2Lm(i1+1,i2-2,i3,j1+1,j2-2,j3),b2Lm(i1+1,i2-1,i3,j1+1,j2-1,j3),b2Lm(i1+1,i2,i3,j1+1,j2,j3),b2Lm(i1+1,i2+1,i3,j1+1,j2+1,j3),b2Lm(i1+1,i2+2,i3,j1+1,j2+2,j3), \         
          //           b2Lm(i1+2,i2-2,i3,j1+2,j2-2,j3),b2Lm(i1+2,i2-1,i3,j1+2,j2-1,j3),b2Lm(i1+2,i2,i3,j1+2,j2,j3),b2Lm(i1+2,i2+1,i3,j1+2,j2+1,j3),b2Lm(i1+2,i2+2,i3,j1+2,j2+2,j3));
                                                          

          //     printF("-------------------------------------------\n");
          //     printF("    c11Lr4=%10.3e; c11Lrr4=%10.3e; c11Ls4=%10.3e; c11Lss4=%10.3e; c11Lsss2=%10.3e; c11Lrs4=%10.3e;\n",c11Lr4,c11Lrr4,c11Ls4,c11Lss4,c11Lsss2,c11Lrs4);
          //     printF("    c12Lr4=%10.3e; c12Lrr4=%10.3e; c12Ls4=%10.3e; c12Lss4=%10.3e; c12Lsss2=%10.3e; c12Lrs4=%10.3e;\n",c12Lr4,c12Lrr4,c12Ls4,c12Lss4,c12Lsss2,c12Lrs4);
          //     printF("    c22Lr4=%10.3e; c22Lrr4=%10.3e; c22Ls4=%10.3e; c22Lss4=%10.3e; c22Lsss2=%10.3e; c22Lrs4=%10.3e;\n",c22Lr4,c22Lrr4,c22Ls4,c22Lss4,c22Lsss2,c22Lrs4);
          //     printF("    c1Lr4=%10.3e; c1Lrr4=%10.3e; c1Ls4=%10.3e; c1Lss4=%10.3e; c1Lsss2=%10.3e; c1Lrs4=%10.3e;\n",c1Lr4,c1Lrr4,c1Ls4,c1Lss4,c1Lsss2,c1Lrs4);
          //     printF("    c2Lr4=%10.3e; c2Lrr4=%10.3e; c2Ls4=%10.3e; c2Lss4=%10.3e; c2Lsss2=%10.3e; c2Lrs4=%10.3e;\n",c2Lr4,c2Lrr4,c2Ls4,c2Lss4,c2Lsss2,c2Lrs4);
          //     printF("    b1Lr4=%10.3e; b1Lrr4=%10.3e; b1Ls4=%10.3e; b1Lss4=%10.3e; b1Lsss2=%10.3e; b1Lrs4=%10.3e;\n",b1Lr4,b1Lrr4,b1Ls4,b1Lss4,b1Lsss2,b1Lrs4);
          //     printF("    b2Lr4=%10.3e; b2Lrr4=%10.3e; b2Ls4=%10.3e; b2Lss4=%10.3e; b2Lsss2=%10.3e; b2Lrs4=%10.3e;\n",b2Lr4,b2Lrr4,b2Ls4,b2Lss4,b2Lsss2,b2Lrs4);
          //     printF("-------------------------------------------\n");
          //     printF("    c11Rr4=%10.3e; c11Rrr4=%10.3e; c11Rs4=%10.3e; c11Rss4=%10.3e; c11Rsss2=%10.3e; c11Rrs4=%10.3e;\n",c11Rr4,c11Rrr4,c11Rs4,c11Rss4,c11Rsss2,c11Rrs4);
          //     printF("    c12Rr4=%10.3e; c12Rrr4=%10.3e; c12Rs4=%10.3e; c12Rss4=%10.3e; c12Rsss2=%10.3e; c12Rrs4=%10.3e;\n",c12Rr4,c12Rrr4,c12Rs4,c12Rss4,c12Rsss2,c12Rrs4);
          //     printF("    c22Rr4=%10.3e; c22Rrr4=%10.3e; c22Rs4=%10.3e; c22Rss4=%10.3e; c22Rsss2=%10.3e; c22Rrs4=%10.3e;\n",c22Rr4,c22Rrr4,c22Rs4,c22Rss4,c22Rsss2,c22Rrs4);
          //     printF("    c1Rr4=%10.3e; c1Rrr4=%10.3e; c1Rs4=%10.3e; c1Rss4=%10.3e; c1Rsss2=%10.3e; c1Rrs4=%10.3e;\n",c1Rr4,c1Rrr4,c1Rs4,c1Rss4,c1Rsss2,c1Rrs4);
          //     printF("    c2Rr4=%10.3e; c2Rrr4=%10.3e; c2Rs4=%10.3e; c2Rss4=%10.3e; c2Rsss2=%10.3e; c2Rrs4=%10.3e;\n",c2Rr4,c2Rrr4,c2Rs4,c2Rss4,c2Rsss2,c2Rrs4);
          //     printF("    b1Rr4=%10.3e; b1Rrr4=%10.3e; b1Rs4=%10.3e; b1Rss4=%10.3e; b1Rsss2=%10.3e; b1Rrs4=%10.3e;\n",b1Rr4,b1Rrr4,b1Rs4,b1Rss4,b1Rsss2,b1Rrs4);
          //     printF("    b2Rr4=%10.3e; b2Rrr4=%10.3e; b2Rs4=%10.3e; b2Rss4=%10.3e; b2Rsss2=%10.3e; b2Rrs4=%10.3e;\n",b2Rr4,b2Rrr4,b2Rs4,b2Rss4,b2Rsss2,b2Rrs4);
          //     printF("-------------------------------------------\n");
                            
          //   }

          //   ForStencil(m1,m2,m3)
          //   {
          //     const int m  = M123(m1,m2,m3);        // the single-component coeff-index
          //     // const int mm = M123CE(m1,m2,m3,c,e);  // the system coeff-index

          //     //printF(" (i1,i2)=(%3d,%3d) (j1,j2)=(%3d,%3d) anl=(%g,%g) m=%d: xCoeff=%g yCoeff=%g r2Coeff=%g r2r2Coeff=%g \n",
          //     //        i1,i2,j1,j2,an1L,an2L,m,xCoeff(m,i1,i2,i3),yCoeff(m,i1,i2,i3),r2Coeff(m,j1,j2,j3),r2r2Coeff(m,j1,j2,j3));
                            
          //     //printF("    b1L=%g, b2L=%g, c11L=%g, c12L=%g, c21L=%g, c22L=%g c1L=%g c2L=%g\n",
          //     //        b1L,b2L,c11L,c12L,c21L,c22L,c1L,c2L); 
          //     //printF("    b1R=%g, b2R=%g, c11R=%g, c12R=%g, c21R=%g, c22R=%g c1R=%g c2R=%g\n",
          //     //        b1R,b2R,c12R,c21R,c22R,c1R,c2R); 
          //     //printF(" b1Lr2=%g, b2Lr2=%g, b1Rr2=%g, b2Rr2=%g \n",b1Lr2,b2Lr2,b1Rr2,b2Rr2);
                            
          //     //printF("    b2Lb1Rr2=%g, b2Rb1Rr2=%g \n",b2Lb1Rr2,b2Rb1Rr2); 
          //                    // Nlr = theta*nL*grad_L - b2*D_s
          //     //printF("lapCoeff=%g \n",lapCoeff(m,i1,i2,i3));
                            
          //     NCoeff(m) = theta*( an1L*xCoeff(m,i1,i2,i3) + an2L*yCoeff(m,i1,i2,i3) ) - b2R*r2Coeff(m,j1,j2,j3);
                            
          //     // NOTE: Here we assume the tangential r derivatives are the same on both sides! (see champ4/notes) 
                            
          //     LCoeff(m) = (b1R*b1R/c11R)*( beta*lapCoeff(m,i1,i2,i3) 
          //                                 - (c12R+c21R)*( theta*( (b1L/b1R)*r1r2Coeff(m,i1,i2,i3) + (b2L/b1R)*r2r2Coeff(m,j1,j2,j3) + b2Lb1Rr2*r2Coeff(m,j1,j2,j3) )
          //                                                 - b2Rb1Rr2*r2Coeff(m,j1,j2,j3) - (b2R/b1R)*r2r2Coeff(m,j1,j2,j3) )
          //                                 - c22R*r2r2Coeff(m,j1,j2,j3) 
          //                                 - (c1R/b1R)*NCoeff(m) 
          //                                 - c2R*r2Coeff(m,j1,j2,j3) );

          //     // LCoeff(m) = (b1R*b1R/c11R)*( beta*( xxCoeff(m,i1,i2,i3)  + yyCoeff(m,i1,i2,i3) ) ); // TEST 
          //   }
          //   //printF(" NCoeff=(%g,%g,%g,%g,%g,%g,%g,%g,%g)\n", NCoeff(0),NCoeff(1),NCoeff(2),NCoeff(3),NCoeff(4),NCoeff(5),NCoeff(6),NCoeff(7),NCoeff(8));
          //   //printF(" LCoeff=(%g,%g,%g,%g,%g,%g,%g,%g,%g)\n", LCoeff(0),LCoeff(1),LCoeff(2),LCoeff(3),LCoeff(4),LCoeff(5),LCoeff(6),LCoeff(7),LCoeff(8));
                        
          //   //printF(" a0=%g, a1=%g, a2=%g\n",a0,a1,a2);

          //   // if( orderOfAccuracy==2 && twilightZoneFlow && !fillMatrixWDH )
          //   if( orderOfAccuracy==2 && twilightZoneFlow  )
          //   {
          //     // --- For twilightZone we save some coefficients that go into the CHAMP matrix ---

          //     // The next expressions are from champMacros.h:  addTwilightZoneCorrectionForChamp()
          //     // 
          //     // Nc = theta*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + cc(Ib1,Ib2,Ib3,0)*uer2; 
          //     // 
          //     // interfaceData.u(Ib1,Ib2,Ib3,n) += Sl*ue   
          //     //                                   + cc(Ib1,Ib2,Ib3,1)*( uexx + ueyy ) 
          //     //                                   + cc(Ib1,Ib2,Ib3,2)*uer1r2 
          //     //                                   + cc(Ib1,Ib2,Ib3,3)*uer2r2
          //     //                                   + cc(Ib1,Ib2,Ib3,4)*uer2
          //     //                                   + cc(Ib1,Ib2,Ib3,5)*Nc;

          //     cc(i1,i2,i3,0) =    ( - b2R);                                                         // coeff of u2r in NCoeff
          //     cc(i1,i2,i3,1) = a2*(b1R*b1R/c11R)*( beta                                        );   // coeff of Lap(u) in a2*LCoeff
          //     cc(i1,i2,i3,2) = a2*(b1R*b1R/c11R)*( -(c12R+c21R)*theta*(b1L/b1R)                );   // coeff of Dr1r2 in a2*LCoeff
          //     cc(i1,i2,i3,3) = a2*(b1R*b1R/c11R)*( -(c12R+c21R)*( (b2L/b1R) + (b2R/b1R)) -c22R );   // coeff of Dr2r2 in a2*LCoeff
          //     cc(i1,i2,i3,4) = a2*(b1R*b1R/c11R)*( -(c12R+c21R)*(  b2Lb1Rr2 + b2Rb1Rr2 ) + c2R );   // coeff of Dr2 in a2*LCoeff
          //     cc(i1,i2,i3,5) = a2*(b1R*b1R/c11R)*( (c1R/b1R) ) + a1;                                // coeff of NCoeff in  a1*NCoeff + a2*LCoeff

          //   }

          //   for( int e=eqnStart; e<=eqnEnd; e++ ) // equation eq1, eq2, ...
          //   {
          //     int c=0; // component number 
          //     ForStencil(m1,m2,m3)
          //     {
          //       int m  = M123(m1,m2,m3);        // the single-component coeff-index
          //       int mm = M123CE(m1,m2,m3,c,e);  // the system coeff-index
          //       // n.grad = n1*D_x + n2*D_y =
          //       //        = ( n1*r.x + n2*r.y ) D_r + ( n1*s.x + n2*s.y ) D_s
          //       //        = b1*D_r + b2*D_s 

          //       // Matrix for the Champ condition: 
          //       if( !fillMatrixWDH )
          //       {
          //         if( mapleOption==0 )
          //         {
          //           coeff(mm,i1m,i2m,i3m) =  a0*idCoeff(m,i1,i2,i3) + a1*NCoeff(mm) + a2*LCoeff(mm); 
          //         }
          //         else if( mapleOption==1 )
          //         {
          //           //printF("---->Use MAPLE generated coefficients.\n");
          //           coeff(mm,i1m,i2m,i3m) = coefA(m1,m2);
          //         }
          //       }
                                
          //       //coeff(mm,i1m,i2m,i3m) =  a0*idCoeff(m,i1,i2,i3) + a1*NCoeff(mm) + a2*LCoeff(mm); 
          //       //printF("idCoef: (i1,i2,i3)=(%d,%d,%d) (m1,m2,m3)=(%d,%d,%d) coeff = %e\n",i1,i2,i3,m1,m2,m3,idCoeff(m,i1,i2,i3));
                                
          //       //printF("LapCoef: (i1,i2,i3)=(%d,%d,%d) (m1,m2,m3)=(%d,%d,%d) coeff = %e\n",i1,i2,i3,m1,m2,m3,lapCoeff(m,i1,i2,i3));
                                
          //       if( multiDomainProblem==0 )
          //       { // debug info 
          //         printF("champBC: (i1,i2,i3)=(%3d,%3d,%d) (m1,m2,m3)=(%3d,%3d,%3d) coeff=%10.3e, coeff=%10.3e (new), diff=%9.2e\n",
          //              i1,i2,i3,m1,m2,m3,coeff(mm,i1m,i2m,i3m),coefA(m1,m2),coeff(mm,i1m,i2m,i3m)-coefA(m1,m2));
          //       }
          //       if( true )
          //       {
          //         // for testing compute any difference between Sijia and WDH
          //         maxDiff = max(maxDiff, fabs(coeff(mm,i1m,i2m,i3m)-coefA(m1,m2)));
          //         maxCoeff = max( maxCoeff,abs(coeff(mm,i1m,i2m,i3m)) );
          //       }                     

          //       // Specify that the above coeff value is the coefficient of component c at the grid point (j1,j2,j3).
          //       const int k1=i1+m1, k2=i2+m2, k3=i3+m3;   // the stencil is centred on the boundary pt (i1,i2,i3)
          //       setEquationNumber(mm, e,i1m,i2m,i3m,  c,k1,k2,k3 );  // macro to set equationNumber

          //       // Fill Ghost 2 -- extrapolation for now 
          //       if( numGhost>1 )
          //       {
          //         const int ghost=2;
          //         if( fillSecondGhostUsingExtrapolation )
          //         {
          //           fillGhostExtrapolation(ghost); 
          //         }
          //         else
          //         {
          //           //fillGhostChamp(ghost,drR[0]); 
          //           const int i1m=i1-is1*ghost, i2m=i2-is2*ghost, i3m=i3-is3*ghost; //  ghost point is (i1m,i2m,i3m)
          //           real hI=ghost*drR[0];
          //           evalChamp4StencilCurvilinear(hI,coefAg);
          //           coeff(mm,i1m,i2m,i3m) = coefAg(m1,m2);
          //           const int j1=i1m + mm*is1, j2=i2m + mm*is2, j3=i3m + mm*is3;     // index of point "m" in extrap formula is shifted in the direction (is1,is2,is3)
          //           setEquationNumber(mm, e,i1m,i2m,i3m,  c,j1,j2,j3 ); 
          //         }
          //       }                
          //    }


          //    if( orderOfAccuracy==2 && (maxDiff > 100.*REAL_EPSILON*maxCoeff) )
          //    {
          //      printF("\n **** champBC:ERROR: CURVILINEAR: There is a difference between coeff from Sijia and WDH, maxDiff=%9.2e *** \n",maxDiff);
          //      ForStencil(m1,m2,m3)  
          //      {
          //        int mm = M123CE(m1,m2,m3,c,e);  // the system coeff-index
                                
          //        printF("champBC: (i1,i2,i3)=(%3d,%3d,%3d) (m1,m2,m3)=(%3d,%3d,%3d) coeff=%10.3e (WDH),  coeff=%10.3e (Sijia) diff=%9.2e\n",
          //              i1,i2,i3,m1,m2,m3,coeff(mm,i1m,i2m,i3m), coefA(m1,m2), coeff(mm,i1m,i2m,i3m)-coefA(m1,m2));
          //      }             
          //      // OV_ABORT("Stop here for now");
          //    }

          //  }  // end for e

          // } // end FOR_3IJD

          // if( debug& 4 && orderOfAccuracy==4 )
          // {
          //   printF("My ghost coefficients: coefAg=(%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n                               //                                          %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n                               //                                          %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n                               //                                          %10.3e,%10.3e,%10.3e,%10.3e,%10.3e,\n                               //                                          %10.3e,%10.3e,%10.3e,%10.3e,%10.3e )\n",//             coefAg(-2,-2),coefAg(-1,-2),coefAg(0,-2),coefAg(1,-2),coefAg(2,-2),//             coefAg(-2,-1),coefAg(-1,-1),coefAg(0,-1),coefAg(1,-1),coefAg(2,-1),//             coefAg(-2, 0),coefAg(-1, 0),coefAg(0, 0),coefAg(1, 0),coefAg(2, 0),//             coefAg(-2, 1),coefAg(-1, 1),coefAg(0, 1),coefAg(1, 1),coefAg(2, 1),//             coefAg(-2, 2),coefAg(-1, 2),coefAg(0, 2),coefAg(1, 2),coefAg(2, 2));
          // }

          // if( true )
          //   printF(">>>>>>> champBC::INFO max-diff in coeff between Sijia and WDH = %9.2e <<<<<<<<<\n",maxDiff);

          // // -- Save a copy fo the CHAMP coefficient matrix for computing the residual in champResidualOpt.bf90 ---
          // sameSide=true; 
          // GridFaceDescriptor & myGFD = getInterfaceGridFaceDescriptor( grid, side, axis, parameters, sameSide );
          // myGFD.dbase.put<RealArray>("coeffChamp");
          // RealArray & coeffChamp = myGFD.dbase.get<RealArray>("coeffChamp");
          // coeffChamp.redim(M,Ib1,Ib2,Ib3); 

          // coeffChamp = coeff(M,Ib1-is1,Ib2-is2,Ib3-is3);


          // if( false ) // --> residual is now computed else-where
          // {
          //   // ---- DOUBLE CHECK THE CHAMP COEFFICIENTS ------
          //   // Compute the residual in using the matrix coefficients:
          //   //      A*u - ( difference stencil )
          //   // 
          //   // OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));

          //   Index I1,I2,I3;
          //   // Index Ig1,Ig2,Ig3;
          //   getIndex(mg.dimension(),I1,I2,I3);  // no need to do so many points 

          //   RealArray u(I1,I2,I3,1);
          //   OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
          //   // int rectangular=0;
          //   // // Set u to exact solution (not really necessary)
          //   // Real t=.1; /// jusr choose a value
          //   // Range N=numberOfComponents;
          //   // e.gd( u ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,I1,I2,I3,N,t);

          //   // Choose some u: 
          //   u = 1. + xLocal(I1,I2,I3,0)*( .5 + .25*xLocal(I1,I2,I3,0) )
          //          + xLocal(I1,I2,I3,1)*( .5 + .25*xLocal(I1,I2,I3,1) );  

          //   // *wdh* April 16, 2022 -- fix for periodic in tangnetial direction
          //   getBoundaryIndex(mg.indexRange(),side,axis,Ib1,Ib2,Ib3);
          //   // getGhostIndex(   mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);

          //   RealArray res(Ib1,Ib2,Ib3);
          //   res =0.;  
          //   FOR_3(i1,i2,i3,Ib1,Ib2,Ib3)
          //   {
          //     int i1m=i1-is1, i2m=i2-is2, i3m=i3-is3; //  ghost point is (i1m,i2m,i3m)
          //     ForStencil(m1,m2,m3)
          //     {
          //        const int m  = M123(m1,m2,m3);        // the single-component coeff-index
          //        // printF("(i1,i2)=(%3d,%3d) m=%4d (m1,m2)=(%2d,%2d) coeff=%10.2e\n ",i1,i2,m,m1,m2,coeff(m,i1m,i2m,i3m));

          //        res(i1,i2,i3) += coeff(m,i1m,i2m,i3m)*u(i1+m1,i2+m2,i3+m3);
          //     }
          //   }
          //   // RealArray uexx(Ib1,Ib2,Ib3), ueyy(Ib1,Ib2,Ib3);
          //   // e.gd( uex ,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,N,t);
          //   // e.gd( uey ,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,N,t);
          //   // e.gd( uexx,xLocal,mg.numberOfDimensions(),rectangular,0,2,0,0,Ib1,Ib2,Ib3,N,t);
          //   // e.gd( ueyy,xLocal,mg.numberOfDimensions(),rectangular,0,0,2,0,Ib1,Ib2,Ib3,N,t); 


          //   // Define d12 and d22 to use dr from the left: 
          //   const Real dr22L[3] = { 1./(SQR(drL[0])), 1./(SQR(drL[1])), 1./(SQR(drL[2])) };        // for D+D- 
          //   #define d12(dir) dr12L[dir]
          //   #define d22(dir) dr22L[dir]

          //   // Make some definitions for the macros in cgux2a.h 
          //   #define U u 
          //   #define UR2(I1,I2,I3,KD) ( (U(I1+1,I2,I3,KD)-U(I1-1,I2,I3,KD))*d12(0) )
          //   #define US2(I1,I2,I3,KD) ( (U(I1,I2+1,I3,KD)-U(I1,I2-1,I3,KD))*d12(1) )
          //   #define UT2(I1,I2,I3,KD) ( (U(I1,I2,I3+1,KD)-U(I1,I2,I3-1,KD))*d12(2) )


          //   // const Real dr22R[3] = { 1./(SQR(drR[0])), 1./(SQR(drR[1])), 1./(SQR(drR[2])) };        // for D+D- 

          //   FOR_3IJD(i1,i2,i3,Ib1,Ib2,Ib3,j1,j2,j3,Jb1,Jb2,Jb3) // loop over points on BOTH SIDES of the interface 
          //   {

          //     // Evaluate variables in the champ conditions:
          //     evalChampVariables();

          //     const Real hR = (1-2*side2)*drR[axis2]/b1R; // scaled grid spacing on right
          //     const real a0 = Sl;
          //     const real a1 = 1. + hR*Sl;
          //     const real a2 = hR*(1.0 + .5*hR*Sl);

          //     // 
          //     // NCoeff(m) = theta*( an1L*xCoeff(m,i1,i2,i3) + an2L*yCoeff(m,i1,i2,i3) ) - b2R*r2Coeff(m,j1,j2,j3);
          //     const Real ur2   =  US2(i1,i2,i3,0); // d(u)/(d r_2)
          //     const Real ur1r2 = URS2(i1,i2,i3,0); // d^2(u)/(d r_1 d r_2 )
          //     const Real ur2r2 = USS2(i1,i2,i3,0); // d^2(u)/(d (r_2)^2 )

          //     // NCoeff(m) = theta*( an1L*xCoeff(m,i1,i2,i3) + an2L*yCoeff(m,i1,i2,i3) ) - b2R*r2Coeff(m,j1,j2,j3); 
          //     // LCoeff(m) = (b1R*b1R/c11R)*( beta*lapCoeff(m,i1,i2,i3) 
          //     //                             - (c12R+c21R)*( theta*( (b1L/b1R)*r1r2Coeff(m,i1,i2,i3) + (b2L/b1R)*r2r2Coeff(m,j1,j2,j3) + b2Lb1Rr2*r2Coeff(m,j1,j2,j3) )
          //     //                                             - b2Rb1Rr2*r2Coeff(m,j1,j2,j3) - (b2R/b1R)*r2r2Coeff(m,j1,j2,j3) )
          //     //                             - c22R*r2r2Coeff(m,j1,j2,j3) 
          //     //                             - (c1R/b1R)*NCoeff(m) 
          //     //                             - c2R*r2Coeff(m,j1,j2,j3) ); 

          //     // -- switch to rxLeft for computing ux and uy etc. on left              
          //     #define inverseVertexDerivative rxLeft 
          //     Real Nlr =  theta*( an1L*UX22(i1,i2,i3,0) + an2L*UY22(i1,i2,i3,0) ) - b2R*ur2 ;

          //     Real Llr  = (b1R*b1R/c11R)*( beta*( UXX22(i1,i2,i3,0) + UYY22(i1,i2,i3,0) )
          //                                 - (c12R+c21R)*( theta*( (b1L/b1R)*ur1r2 + (b2L/b1R)*ur2r2 + b2Lb1Rr2*ur2  ) 
          //                                                - b2Rb1Rr2*ur2 - (b2R/b1R)*ur2r2 )
          //                                 - c22R*ur2r2
          //                                 - (c1R/b1R)*Nlr 
          //                                 - c2R*ur2
          //                                );
          //     // Llr= (b1R*b1R/c11R)*( beta*( UXX22(i1,i2,i3,0) + UYY22(i1,i2,i3,0) ) ); // TEST 

          //     res(i1,i2,i3) -= a0*u(i1,i2,i3) + a1*Nlr + a2*Llr; 
          //   }

          //   Real maxRes = max(fabs(res));
          //   printF("CHAMP-MATRIX: max-residual=%8.2e\n",maxRes);
          //   if( maxRes >1.e-10 )
          //   {
          //     OV_ABORT("residual is too large ??")
          //   }

          // }

                    

                } // end curvilinear grid

                if( !multiDomainProblem )
                {
                    OV_ABORT("champBoundaryConditions: TESTING ON SINGLE DOMAIN -- STOP HERE FOR NOW");

                }

        // -- Save a copy of the CHAMP coefficient matrix for computing the residual in champResidualOpt.bf90 ---
                const bool sameSide=true; 
                GridFaceDescriptor & myGFD = getInterfaceGridFaceDescriptor( grid, side, axis, parameters, sameSide );
                if( !myGFD.dbase.has_key("coeffChamp") )
                    myGFD.dbase.put<RealArray>("coeffChamp");

                RealArray & coeffChamp = myGFD.dbase.get<RealArray>("coeffChamp");
                coeffChamp.redim(M,Ib1,Ib2,Ib3); 

                coeffChamp = coeff(M,Ib1-is1,Ib2-is2,Ib3-is3);          
                    
            } // end if( mg.boundaryCondition(side,axis)>0 )
        } // end ForBoundary 
        
        if( debug & 8 )
            ::displayCoeff(coeff,"coeff");
            
    } // end for grid 

    if( orderOfAccuracy!=2 && orderOfAccuracy!=4 )
    {
        printf("champBoundaryConditions: orderOfAccuracy=%d. FINISH ME\n",orderOfAccuracy);
        OV_ABORT("ERROR: FINISH ME");

    }


  // OV_ABORT("champBoundaryConditions : END: STOP HERE FOR NOW");
    
    return 0;

}

