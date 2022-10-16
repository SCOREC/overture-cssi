// --- Eigenvalues of a 2D periodic strip for incompressible elasticity -----
//   BC = TD
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
    6.1046733125517700e+00, 8.2872475749655039e+00, 1.0694170271326492e+01, 1.3041586135065603e+01,
    // m=2, n=1,2,..,4
    1.2008573987474925e+01, 1.3377455285429044e+01, 1.5071190101197233e+01, 1.7199046325461177e+01,
    // m=3, n=1,2,..,4
    1.8007351701186442e+01, 1.9295616953990335e+01, 2.0422956693880925e+01, 2.2014459535174920e+01
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
  A1 = -0.4e1 / (sin(beta) * pow(beta, 0.6e1) - sin(beta) * pow(beta, 0.4e1) * k * k - sin(beta) * beta * beta * pow(k, 0.4e1) + sin(beta) * pow(k, 0.6e1) + 0.4e1 * cos(beta) * pow(beta, 0.3e1) * pow(k, 0.3e1) + 0.4e1 * cos(beta) * beta * pow(k, 0.5e1) + 0.2e1 * exp(k) * pow(beta, 0.5e1) * k - 0.2e1 * exp(k) * beta * pow(k, 0.5e1)) * pow(k, 0.3e1) * (sin(beta) * beta * beta - sin(beta) * k * k + k * exp(k) * beta - k * exp(-k) * beta) * B3;
  A3 = -(sin(beta) * pow(beta, 0.4e1) - 0.2e1 * sin(beta) * beta * beta * k * k + sin(beta) * pow(k, 0.4e1) - 0.4e1 * cos(beta) * beta * pow(k, 0.3e1) - 0.2e1 * exp(-k) * pow(beta, 0.3e1) * k + 0.2e1 * exp(-k) * beta * pow(k, 0.3e1)) * B3 / (sin(beta) * pow(beta, 0.4e1) - 0.2e1 * sin(beta) * beta * beta * k * k + sin(beta) * pow(k, 0.4e1) + 0.4e1 * cos(beta) * beta * pow(k, 0.3e1) + 0.2e1 * exp(k) * pow(beta, 0.3e1) * k - 0.2e1 * exp(k) * beta * pow(k, 0.3e1));
  B1 = (0.4e1 * cos(beta) * k * k + exp(k) * beta * beta - k * k * exp(k) + exp(-k) * beta * beta - k * k * exp(-k)) * k * B3 * (beta * beta - k * k) / (beta * beta + k * k) / (sin(beta) * pow(beta, 0.4e1) - 0.2e1 * sin(beta) * beta * beta * k * k + sin(beta) * pow(k, 0.4e1) + 0.4e1 * cos(beta) * beta * pow(k, 0.3e1) + 0.2e1 * exp(k) * pow(beta, 0.3e1) * k - 0.2e1 * exp(k) * beta * pow(k, 0.3e1));
  // uHat2 = -uHat1.x/ky;
  A2 = -beta*B1/ky;
  B2 =  beta*A1/ky;
}
else
{
  isSurfaceWave=true;
  alpha = sqrt( -betaSq );
  A1 = 0.4e1 * pow(k, 0.3e1) * B3 * (k * exp(-k) * alpha + (alpha * alpha + k * k) * sinh(alpha) - k * exp(k) * alpha) / (pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.2e1 * alpha * k * ((alpha * alpha + k * k) * exp(k) - 0.2e1 * cosh(alpha) * k * k)) / (k - alpha) / (alpha + k);
  B1 = k * ((alpha * alpha + k * k) * exp(-k) + (alpha * alpha + k * k) * exp(k) - 0.4e1 * cosh(alpha) * k * k) * B3 * (alpha * alpha + k * k) / ((-0.2e1 * pow(alpha, 0.3e1) * k - 0.2e1 * alpha * pow(k, 0.3e1)) * exp(k) + pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) + 0.4e1 * cosh(alpha) * alpha * pow(k, 0.3e1)) / (k - alpha) / (alpha + k);
  A3 = -((0.2e1 * pow(alpha, 0.3e1) * k + 0.2e1 * alpha * pow(k, 0.3e1)) * exp(-k) + pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.4e1 * cosh(alpha) * alpha * pow(k, 0.3e1)) * B3 / (pow(alpha * alpha + k * k, 0.2e1) * sinh(alpha) - 0.2e1 * alpha * k * ((alpha * alpha + k * k) * exp(k) - 0.2e1 * cosh(alpha) * k * k));
  // uHat2 = -uHat1.x/ky;
  A2 = -alpha*B1/ky;
  B2 = -alpha*A1/ky;
}
P1p=  A3*kx/(rho*SQR(omega));
P1m= -B3*kx/(rho*SQR(omega));
P2p= -kx*P1p/ky;
P2m= +kx*P1m/ky;
