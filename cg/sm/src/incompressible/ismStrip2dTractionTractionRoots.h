// --- Eigenvalues of a 2D periodic strip for incompressible elasticity -----
//   BC = TT
//   File written by ismStrip2d.mw on 2022-06-21
//   Solution is of the form  : u(x,y,t) = uHat(x) *exp( I*k*y )*exp( I*omega t) , k=2*Pi*m 
//   uHat1 =  A1*cos(beta*x) + B1*sin(beta*x) + pHat'(x)/(mu*(beta^2+k^2) 
//   uHat2 =  from divergence
//   pHat  =  A3*exp(k*x) + B3*exp(-k*x) 
//   beta^2 = (omega/c)^2 - k^2 
//   Note that beta is generally real except for surface waves where it is pure imaginary 
const int mxd=4, myd=3; // max number of roots in the table
const Real omegaArray[]={
    // m=1, n=1,2,..,4
    5.7766060906888207e+00, 6.6005868537057281e+00, 8.8857658763167325e+00, 1.1003069256654355e+01,
    // m=2, n=1,2,..,4
    1.2106218072213127e+01, 1.3753223189861312e+01, 1.5624812844669151e+01, 1.7771531752633465e+01,
    // m=3, n=1,2,..,4
    1.8027501109917280e+01, 1.9472778749087922e+01, 2.0798488827365205e+01, 2.2516638126074782e+01
};
#define omegaRoot(mx,my) omegaArray[(mx-1)+mxd*(my-1)]
omega = omegaRoot(mx,my);
ky = twoPi*my; kx=ky; k=kx;
B3 = scaleFactor;  // scaleFactor should be set before this point 
Real betaSq = SQR(omega/c) - SQR(kx);
if( betaSq >0 )
{
  isSurfaceWave=false;
  beta = sqrt( betaSq );
  A1 = -0.4e1 / (0.4e1 * exp(k) * pow(beta, 0.5e1) * pow(k, 0.3e1) - 0.4e1 * exp(k) * pow(k, 0.7e1) * beta - sin(beta) * pow(beta, 0.8e1) + 0.2e1 * sin(beta) * pow(beta, 0.6e1) * k * k - 0.2e1 * sin(beta) * beta * beta * pow(k, 0.6e1) + sin(beta) * pow(k, 0.8e1) - 0.4e1 * cos(beta) * pow(beta, 0.5e1) * pow(k, 0.3e1) + 0.4e1 * cos(beta) * beta * pow(k, 0.7e1)) * pow(k, 0.3e1) * (0.2e1 * pow(k, 0.3e1) * exp(k) * beta - sin(beta) * pow(beta, 0.4e1) + 0.2e1 * sin(beta) * beta * beta * k * k - sin(beta) * pow(k, 0.4e1) - 0.2e1 * pow(k, 0.3e1) * exp(-k) * beta) * B3;
  A3 = (sin(beta) * pow(beta, 0.4e1) - 0.2e1 * sin(beta) * beta * beta * k * k + sin(beta) * pow(k, 0.4e1) - 0.4e1 * cos(beta) * beta * pow(k, 0.3e1) + 0.4e1 * pow(k, 0.3e1) * exp(-k) * beta) * B3 / (0.4e1 * pow(k, 0.3e1) * exp(k) * beta - sin(beta) * pow(beta, 0.4e1) + 0.2e1 * sin(beta) * beta * beta * k * k - sin(beta) * pow(k, 0.4e1) - 0.4e1 * cos(beta) * beta * pow(k, 0.3e1));
  B1 = 0.2e1 * (exp(k) * beta * beta - exp(k) * k * k - 0.2e1 * cos(beta) * beta * beta + 0.2e1 * cos(beta) * k * k + exp(-k) * beta * beta - exp(-k) * k * k) * pow(k, 0.3e1) * B3 / (beta * beta + k * k) / (0.4e1 * pow(k, 0.3e1) * exp(k) * beta - sin(beta) * pow(beta, 0.4e1) + 0.2e1 * sin(beta) * beta * beta * k * k - sin(beta) * pow(k, 0.4e1) - 0.4e1 * cos(beta) * beta * pow(k, 0.3e1));
  // uHat2 = -uHat1.x/ky;
  A2 = -beta*B1/ky;
  B2 =  beta*A1/ky;
}
else
{
  isSurfaceWave=true;
  alpha = sqrt( -betaSq );
  A1 = 0.4e1 * B3 * (0.2e1 * pow(k, 0.3e1) * exp(-k) * alpha + pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.2e1 * pow(k, 0.3e1) * exp(k) * alpha) * pow(k, 0.3e1) / (alpha * alpha + k * k) / (alpha + k) / (pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.4e1 * pow(k, 0.3e1) * alpha * (exp(k) - cosh(alpha))) / (k - alpha);
  B1 = 0.2e1 * (alpha * alpha + k * k) * (exp(k) - 0.2e1 * cosh(alpha) + exp(-k)) * pow(k, 0.3e1) * B3 / (alpha + k) / (pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.4e1 * pow(k, 0.3e1) * alpha * (exp(k) - cosh(alpha))) / (k - alpha);
  A3 = -B3 * (0.4e1 * pow(k, 0.3e1) * exp(-k) * alpha + pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.4e1 * cosh(alpha) * alpha * pow(k, 0.3e1)) / (pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.4e1 * pow(k, 0.3e1) * alpha * (exp(k) - cosh(alpha)));
  // uHat2 = -uHat1.x/ky;
  A2 = -alpha*B1/ky;
  B2 = -alpha*A1/ky;
}
P1p=  A3*kx/(rho*SQR(omega));
P1m= -B3*kx/(rho*SQR(omega));
P2p= -kx*P1p/ky;
P2m= +kx*P1m/ky;
