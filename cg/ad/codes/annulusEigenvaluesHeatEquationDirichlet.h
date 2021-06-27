// ---- Eigenvalues for eigenfunctions of the heat equation on an annulus ----
// File written by cg/ad/codes/annulusEigenvalues.maple 
// numBesselOrderDBC  : Bessel orders are m=0,1,2,...,numBesselOrderDBC-1 
// numRootDBC         : number of roots
// ----------- DIRCHLET BCs at r=ra and r=rb ----- 
// Solution is u = cJ*Jn(lambda*r) + cY*Yn(lambda*r)
//  cJ = amp*Yn(lambda*a), cY = -amp*Jn(lambda*a) 
 const Real ra=5.00000000000000e-01, rb=1.00000000000000e+00; 
 const int numBesselOrder=3, numRoot=4;
 Real annulusEigs[numBesselOrder][numRoot]={
 { // Roots for f(z,n) = Jn(a*z)*Yn(a*z) - Jn(b*z)*Yn(a*z), n=0
    6.24606183919138e+00,
    1.25468714279844e+01,
    1.88364150845032e+01,
    2.51228463710507e+01}, // end n=0
 { // Roots for f(z,n) = Jn(a*z)*Yn(a*z) - Jn(b*z)*Yn(a*z), n=1
    6.39315676162127e+00,
    1.26246990207465e+01,
    1.88889298509645e+01,
    2.51624056202082e+01}, // end n=1
 { // Roots for f(z,n) = Jn(a*z)*Yn(a*z) - Jn(b*z)*Yn(a*z), n=2
    6.81384285313505e+00,
    1.28555318451921e+01,
    1.90457045399067e+01,
    2.52807623389876e+01}, // end n=2
                               };
