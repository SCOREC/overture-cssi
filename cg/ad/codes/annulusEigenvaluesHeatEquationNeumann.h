// ---- Eigenvalues for eigenfunctions of the heat equation on an annulus ----
// File written by cg/ad/codes/annulusEigenvalues.maple 
// numBesselOrderDBC  : Bessel orders are m=0,1,2,...,numBesselOrderDBC-1 
// numRootDBC         : number of roots
// ----------- NEUMANN BCs at r=ra and r=rb ----- 
// Solution is u = cJ*Jn(lambda*r) + cY*Yn(lambda*r)
//  cJ = amp*Yn(lambda*a), cY = -amp*Jn(lambda*a) 
 const Real ra=5.00000000000000e-01, rb=1.00000000000000e+00; 
 const int numBesselOrder=5, numRoot=5;
 Real annulusEigs[numBesselOrder][numRoot]={
 { // Roots for f(z,n) = Jn(a*z)*Yn(a*z) - Jn(b*z)*Yn(a*z), n=0
    6.39315676162127e+00,
    1.26246990207465e+01,
    1.88889298509645e+01,
    2.51624056202082e+01,
    3.14397085388595e+01}, // end n=0
 { // Roots for f(z,n) = Jn(a*z)*Yn(a*z) - Jn(b*z)*Yn(a*z), n=1
    1.35467201027317e+00,
    6.56494238232276e+00,
    1.27064223370974e+01,
    1.89426593061038e+01,
    2.52024870402830e+01}, // end n=1
 { // Roots for f(z,n) = Jn(a*z)*Yn(a*z) - Jn(b*z)*Yn(a*z), n=2
    7.06258161604745e+00,
    1.29494113826463e+01,
    1.91031584965659e+01,
    2.53224312291551e+01,
    3.15674823440921e+01}, // end n=2
 { // Roots for f(z,n) = Jn(a*z)*Yn(a*z) - Jn(b*z)*Yn(a*z), n=3
    7.84010909785816e+00,
    1.33475998904036e+01,
    1.93684288354592e+01,
    2.55213543943649e+01,
    3.17266192967601e+01}, // end n=3
 { // Roots for f(z,n) = Jn(a*z)*Yn(a*z) - Jn(b*z)*Yn(a*z), n=4
    8.83644071386484e+00,
    1.38922797245307e+01,
    1.97353588717768e+01,
    2.57978368957929e+01,
    3.19483474572928e+01}, // end n=4
                               };
