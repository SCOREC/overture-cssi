// ------------------------------------------------------------------------------------------------ 
//    AnnulusDispDispBCs: Exact solution for incompressible linear elasticity 
// Annulus of outer radius Rb=2, inner radius Ra = 1; r=1 Displacement BC; r=2 Displacement BC. 
// File Written by Disk_Annulus_ExactSoln.mw (maple) 
// Set C4 to scale the solution 
// ------------------------------------------------------------------------------------------------ 
Real rho,c,mu,omegaRoot,betaRoot,Pi,k,C1,C2,C3,ra,rb; 
ra=.5;
rb=1;
rho=1;
mu=1;
betaRoot=12.5105206473408546981549;
k=1;
C1 = 0.192992035693309510535441e2 * C4;
C2 = -0.304935551313522923733292e2 * C4;
C3 = 0.992579994490862297289017e0 * C4;
omegaRoot = betaRoot * sqrt(mu / rho);
Pi=3.14159265358979323846264;
#beginMacro annulusMacro(r,theta,t,ur,vr,pr)
  ur = 0.1e1 / r * jn(1, 0.125105206473408546981549e2 * r) * C4 * cos(t * omegaRoot) * cos(theta) + 0.1e1 / r * jn(1, 0.125105206473408546981549e2 * r) * C4 * sin(t * omegaRoot) * sin(theta) + 0.992579994490862297289017e0 / r * yn(1, 0.125105206473408546981549e2 * r) * C4 * cos(t * omegaRoot) * cos(theta) + 0.992579994490862297289017e0 / r * yn(1, 0.125105206473408546981549e2 * r) * C4 * sin(t * omegaRoot) * sin(theta) + 0.123307251957620161779820e0 * C4 * cos(t * omegaRoot) * cos(theta) + 0.123307251957620161779820e0 * C4 * sin(t * omegaRoot) * sin(theta) + 0.194830655687808236749749e0 * pow(r, -0.2e1) * C4 * cos(t * omegaRoot) * cos(theta) + 0.194830655687808236749749e0 * pow(r, -0.2e1) * C4 * sin(t * omegaRoot) * sin(theta);
  vr = -0.100000000000000000000000e1 / r * jn(1, 0.125105206473408546981549e2 * r) * C4 * sin(t * omegaRoot) * cos(theta) + 0.100000000000000000000000e1 / r * jn(1, 0.125105206473408546981549e2 * r) * C4 * cos(t * omegaRoot) * sin(theta) + 0.124176925152154045774849e2 * yn(0, 0.125105206473408546981549e2 * r) * C4 * sin(t * omegaRoot) * cos(theta) - 0.124176925152154045774849e2 * yn(0, 0.125105206473408546981549e2 * r) * C4 * cos(t * omegaRoot) * sin(theta) + 0.125105206473408546981549e2 * jn(0, 0.125105206473408546981549e2 * r) * C4 * sin(t * omegaRoot) * cos(theta) - 0.125105206473408546981549e2 * jn(0, 0.125105206473408546981549e2 * r) * C4 * cos(t * omegaRoot) * sin(theta) - 0.992579994490862297289018e0 / r * yn(1, 0.125105206473408546981549e2 * r) * C4 * sin(t * omegaRoot) * cos(theta) + 0.992579994490862297289018e0 / r * yn(1, 0.125105206473408546981549e2 * r) * C4 * cos(t * omegaRoot) * sin(theta) + 0.123307251957620161779820e0 * C4 * sin(t * omegaRoot) * cos(theta) - 0.123307251957620161779820e0 * C4 * cos(t * omegaRoot) * sin(theta) - 0.194830655687808236749749e0 * pow(r, -0.2e1) * C4 * sin(t * omegaRoot) * cos(theta) + 0.194830655687808236749749e0 * pow(r, -0.2e1) * C4 * cos(t * omegaRoot) * sin(theta);
  pr = 0.192992035693309510535441e2 * C4 * mu * pow(r, k) * cos(k * theta) * cos(t * omegaRoot) + 0.192992035693309510535441e2 * C4 * mu * pow(r, k) * sin(k * theta) * sin(t * omegaRoot) - 0.304935551313522923733292e2 * C4 * mu / pow(r, k) * cos(k * theta) * cos(t * omegaRoot) - 0.304935551313522923733292e2 * C4 * mu / pow(r, k) * sin(k * theta) * sin(t * omegaRoot);
#endMacro
