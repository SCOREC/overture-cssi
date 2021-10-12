// ------------------------------------------------------------------------------------------------ 
//    DiskDispBC: Exact solution for incompressible linear elasticity 
// Solid disk of radius R=1; r=0 bounded solution; r=1 Displacement BC. 
// File Written by Disk_Annulus_ExactSoln.mw (maple) 
// Set C4 to scale the solution 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omegaRoot,betaRoot,Pi,k,C1,C2,C3,ra,rb; 

rho=1;
mu=1;
betaRoot=6.38016189592398350623661;
k=2;
C1 = 0.0e0;
C2 = 0.607082650035880738991908e1 * C4;
omegaRoot = betaRoot * sqrt(mu / rho);
Pi=3.14159265358979323846264;
#beginMacro diskMacro(r,theta,t,ur,vr,pr)
  ur = 0.1e1 / r * jn(k, 0.638016189592398350623661e1 * r) * C4 * cos(k * theta) * cos(t * omegaRoot) + 0.1e1 / r * jn(k, 0.638016189592398350623661e1 * r) * C4 * sin(k * theta) * sin(t * omegaRoot) - pow(r, k) / r * jn(k, 0.638016189592398350623661e1) * C4 * cos(k * theta) * cos(t * omegaRoot) - pow(r, k) / r * jn(k, 0.638016189592398350623661e1) * C4 * sin(k * theta) * sin(t * omegaRoot);
  vr = 0.638016189592398350623661e1 / (double) k * jn(k + 1, 0.638016189592398350623661e1 * r) * C4 * sin((double) (k * theta)) * cos(t * omegaRoot) - 0.638016189592398350623661e1 / (double) k * jn(k + 1, 0.638016189592398350623661e1 * r) * C4 * cos((double) (k * theta)) * sin(t * omegaRoot) - 0.100000000000000000000000e1 / r * jn(k, 0.638016189592398350623661e1 * r) * C4 * sin((double) (k * theta)) * cos(t * omegaRoot) + 0.100000000000000000000000e1 / r * jn(k, 0.638016189592398350623661e1 * r) * C4 * cos((double) (k * theta)) * sin(t * omegaRoot) + 0.100000000000000000000000e1 * pow(r, (double) k) / r * jn(k, 0.638016189592398350623661e1) * C4 * sin((double) (k * theta)) * cos(t * omegaRoot) - 0.100000000000000000000000e1 * pow(r, (double) k) / r * jn(k, 0.638016189592398350623661e1) * C4 * cos((double) (k * theta)) * sin(t * omegaRoot);
  pr = -0.407064658182003197420524e2 * jn(k, 0.638016189592398350623661e1) * C4 / k * mu * pow(r, k) * cos(k * theta) * cos(t * omegaRoot) - 0.407064658182003197420524e2 * jn(k, 0.638016189592398350623661e1) * C4 / k * mu * pow(r, k) * sin(k * theta) * sin(t * omegaRoot);
#endMacro
