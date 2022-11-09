// ---- Two concentric annulii : Eigenvalues for eigenfunctions of the heat equation CHT problem ----
// See documenation for CgAd in cgDoc/mp/tex/doubleAnnulus.tex for derivation and formulae. 
// File written by cg/ad/codes/concentricAnnulusEigenvalues.maple 
// numBesselOrderDBC  : Bessel orders are m=0,1,2,...,numBesselOrderDBC-1 
// numRootDBC         : number of roots
// ----------- DIRCHLET BCs at r=ra and r=rc, INTERFACE at r=rb ----- 
// Solution is u = [ c1J*Jn(lam1*r) + c1Y*Yn(lam1*r) ] cos(n*theta) exp(- s^2*t ),   a < r < b
//             u = [ c2J*Jn(lam2*r) + c2Y*Yn(lam2*r) ] cos(n*theta) exp(- s^2*t),   b < r < c
// where s_{mn} is the eigenvalue for n=0,1,2,3... and m=,1,2,3,... 
//  c1J = Jn(lam2*b) - (Jn(lam2*c)/(Yn(lam2*c)*Yn(lam2*b) 
//  c2J = Jn(lam1*b) - (Jn(lam1*a)/(Yn(lam1*a)*Yn(lam1*b) 
//  c1Y = -Jn(lam1*a)/Yn(lam1*a) * c1J 
//  c2Y = -Jn(lam2*c)/Yn(lam2*c) * c2J 
//  lam1=alpha1*s, lam2=alpha2*s, alpha1=1/sqrt{D1} alpha2=1/sqrt{D2} 
 // Note: Kr = K2/K1 : ratio of thermal conductivitites (K1=inner, K2=outer).
 const Real ra=5.00000000000000e-01, rb=1.00000000000000e+00, rc=1.50000000000000e+00, Kr=5.00000000000000e-01; 
 const Real D1=9.00000000000000e-01, D2=8.00000000000000e-01; 
 const int numBesselOrder=3, numRoot=3;
 Real concentricAnnulusEigs[numBesselOrder][numRoot]={
 { // Roots of determinant condition for n=0
    2.77023150510382e+00,
    5.71364452131061e+00,
    8.70518643026619e+00}, // end n=0
 { // Roots of determinant condition for n=1
    2.95217791037101e+00,
    5.79224396583258e+00,
    8.77892976493773e+00}, // end n=1
 { // Roots of determinant condition for n=2
    3.42942342921976e+00,
    6.02341926474434e+00,
    8.99649908763766e+00}, // end n=2
                               };
