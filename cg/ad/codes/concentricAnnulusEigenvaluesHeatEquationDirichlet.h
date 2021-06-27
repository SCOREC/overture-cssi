// ---- Two concentric annulii : Eigenvalues for eigenfunctions of the heat equation ----
// See documenation for CgAd in cgDoc/ad for derivation and formulae. 
// File written by cg/ad/codes/concentricAannulusEigenvalues.maple 
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
 const Real ra=5.00000000000000e-01, rb=1.00000000000000e+00, rc=1.50000000000000e+00, Kr=2.00000000000000e+00; 
 const Real D1=1.00000000000000e-01, D2=5.00000000000000e-02; 
 const int numBesselOrder=3, numRoot=4;
 Real concentricAnnulusEigs[numBesselOrder][numRoot]={
 { // Roots of determinant condition for n=0
    1.75976787679971e+00,
    2.33858965449530e+00,
    3.39615006714118e+00,
    4.05226914831230e+00}, // end n=0
 { // Roots of determinant condition for n=1
    1.79216409479020e+00,
    2.35841530708184e+00,
    3.40764898946815e+00,
    4.06912736649571e+00}, // end n=1
 { // Roots of determinant condition for n=2
    1.88161256990436e+00,
    2.42053214361087e+00,
    3.44109725016957e+00,
    4.11998596651301e+00}, // end n=2
                               };
