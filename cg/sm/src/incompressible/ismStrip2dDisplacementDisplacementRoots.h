// --- Eigenvalues of a 2D periodic strip for incompressible elasticity -----
//   BC = DD
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
    7.6465512957455341e+00, 1.0252755871737867e+01, 1.3170257593343554e+01, 1.6177532569133238e+01,
    // m=2, n=1,2,..,4
    1.3104495532075780e+01, 1.4553250335532520e+01, 1.6601752731685955e+01, 1.9009502048418872e+01,
    // m=3, n=1,2,..,4
    1.9173530688079666e+01, 2.0105183319631333e+01, 2.1545698437568630e+01, 2.3380429851173121e+01
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
  A1 = 0.1e1 / (cos(beta) * pow(beta, 0.3e1) + cos(beta) * k * k * beta + beta * beta * k * sin(beta) + sin(beta) * pow(k, 0.3e1) - exp(k) * pow(beta, 0.3e1) - k * k * exp(k) * beta) * (0.2e1 * sin(beta) * k - beta * exp(k) + beta * exp(-k)) * k * B3;
  A3 = (beta * cos(beta) - sin(beta) * k - beta * exp(-k)) * B3 / (beta * cos(beta) + sin(beta) * k - beta * exp(k));
  B1 = -(-exp(k) - exp(-k) + 0.2e1 * cos(beta)) * k * k * B3 / (beta * beta + k * k) / (beta * cos(beta) + sin(beta) * k - beta * exp(k));
  // uHat2 = -uHat1.x/ky;
  A2 = -beta*B1/ky;
  B2 =  beta*A1/ky;
}
else
{
  isSurfaceWave=true;
  alpha = sqrt( -betaSq );
  A1 = (-alpha * exp(k) + alpha * exp(-k) + 0.2e1 * sinh(alpha) * k) * k * B3 / (-alpha * exp(k) + alpha * cosh(alpha) + sinh(alpha) * k) / (-alpha + k) / (alpha + k);
  B1 = -(exp(k) + exp(-k) - 0.2e1 * cosh(alpha)) * k * k * B3 / (-alpha * exp(k) + alpha * cosh(alpha) + sinh(alpha) * k) / (-alpha * alpha + k * k);
  A3 = (alpha * exp(-k) - alpha * cosh(alpha) + sinh(alpha) * k) / (alpha * exp(k) - alpha * cosh(alpha) - sinh(alpha) * k) * B3;
  // uHat2 = -uHat1.x/ky;
  A2 = -alpha*B1/ky;
  B2 = -alpha*A1/ky;
}
P1p=  A3*kx/(rho*SQR(omega));
P1m= -B3*kx/(rho*SQR(omega));
P2p= -kx*P1p/ky;
P2m= +kx*P1m/ky;
