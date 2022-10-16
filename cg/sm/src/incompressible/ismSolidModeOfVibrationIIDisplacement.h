// File written by AMP/ims/maple-DWS/solidVibrationModes.mpl
// =============================================================================================
// Macro: evaluate the solid mode of vibration II for an incompressible solid, Displacement BC
// =============================================================================================
#beginMacro evalSolidMode2Traction(r,theta,phi, u1,u2,u3,p)
  // Lp = the associated Legendre function
    Lp = 0.3e1 * (cos(phi) - 0.1e1) * (cos(phi) + 0.1e1);
  // phip = Pi - phi
  if( r>eps && phi>eps && phip>eps )
  { // r, phi and phip=Pi=phi are big
    u1 = -0.10e2 / 0.3e1 * (cos(phi) + 0.1e1) * ((((double) (-6 * r * r * k * k + 15) * pow(cos(phi), 0.2e1) + (double) (3 * r * r * k * k) - 0.9e1) * sin((double) (k * r)) + (double) k * ((0.9e1 + (double) (r * r * k * k - 15) * pow(cos(phi), 0.2e1)) * cos((double) (k * r)) - 0.3e1 / 0.5e1 * (double) (k * k) * (double)   pow((double) r, (double) 4) * (double) OmegaFact) * (double) r) * cos(theta) * cos(0.2e1 * theta) + ((double) (-3 * r * r * k * k + 6) * sin((double) (k * r)) + (double) k * (double) r * ((double) (r * r * k * k - 6) * cos((double) (k * r)) - 0.3e1 / 0.5e1 * (double) (k * k) * (double)   pow((double) r, (double) 4) * (double) OmegaFact)) * sin(theta) * sin(0.2e1 * theta)) * Lp * (cos(phi) - 0.1e1) * pow(sin(phi), -0.3e1) * (double)   pow((double) k, (double) (-5)) * (double)   pow((double) r, (double) (-4));
    u2 = -0.10e2 / 0.3e1 * (cos(phi) + 0.1e1) * ((((double) (-6 * r * r * k * k + 15) * pow(cos(phi), 0.2e1) + (double) (3 * r * r * k * k) - 0.9e1) * sin((double) (k * r)) + (double) k * ((0.9e1 + (double) (r * r * k * k - 15) * pow(cos(phi), 0.2e1)) * cos((double) (k * r)) - 0.3e1 / 0.5e1 * (double) (k * k) * (double)   pow((double) r, (double) 4) * (double) OmegaFact) * (double) r) * sin(theta) * cos(0.2e1 * theta) - ((double) (-3 * r * r * k * k + 6) * sin((double) (k * r)) + (double) k * (double) r * ((double) (r * r * k * k - 6) * cos((double) (k * r)) - 0.3e1 / 0.5e1 * (double) (k * k) * (double)   pow((double) r, (double) 4) * (double) OmegaFact)) * sin(0.2e1 * theta) * cos(theta)) * Lp * (cos(phi) - 0.1e1) * pow(sin(phi), -0.3e1) * (double)   pow((double) k, (double) (-5)) * (double)   pow((double) r, (double) (-4));
    u3  = 0.10e2 / 0.3e1 * cos(phi) * cos((double) (2 * theta)) * (cos(k * r) * pow(k, 0.3e1) * pow(r, 0.3e1) - 0.6e1 * sin(k * r) * k * k * r * r - 0.15e2 * cos(k * r) * k * r + 0.15e2 * sin(k * r)) * Lp * pow(sin(phi), -0.2e1) * pow(k, -0.5e1) * pow(r, -0.4e1) * (pow(cos(phi), 0.2e1) - 0.1e1);
    p = -OmegaFact * r * r * Lp * cos((double) (2 * theta));
  }
  else if( r>eps )
  { // small phi or phip but not small r -- series error is O(phi^3)
    if( phi<eps ) // small phi but not small r
    { 
    u1 = -0.10e2 * (((double) (-3 * r * r * k * k + 6) * sin((double) (k * r)) + (double) k * (double) r * ((double) (r * r * k * k - 6) * cos((double) (k * r)) - 0.3e1 / 0.5e1 * (double) (k * k) * (double)   pow((double) r, (double) 4) * (double) OmegaFact)) * cos(theta) * cos(0.2e1 * theta) + ((double) (-3 * r * r * k * k + 6) * sin((double) (k * r)) + (double) k * (double) r * ((double) (r * r * k * k - 6) * cos((double) (k * r)) - 0.3e1 / 0.5e1 * (double) (k * k) * (double)   pow((double) r, (double) 4) * (double) OmegaFact)) * sin(theta) * sin(0.2e1 * theta)) * (double)   pow((double) k, (double) (-5)) * (double)   pow((double) r, (double) (-4)) * phi;
    u2 = -0.10e2 * (((double) (-3 * r * r * k * k + 6) * sin((double) (k * r)) + (double) k * (double) r * ((double) (r * r * k * k - 6) * cos((double) (k * r)) - 0.3e1 / 0.5e1 * (double) (k * k) * (double)   pow((double) r, (double) 4) * (double) OmegaFact)) * sin(theta) * cos(0.2e1 * theta) - ((double) (-3 * r * r * k * k + 6) * sin((double) (k * r)) + (double) k * (double) r * ((double) (r * r * k * k - 6) * cos((double) (k * r)) - 0.3e1 / 0.5e1 * (double) (k * k) * (double)   pow((double) r, (double) 4) * (double) OmegaFact)) * sin(0.2e1 * theta) * cos(theta)) * (double)   pow((double) k, (double) (-5)) * (double)   pow((double) r, (double) (-4)) * phi;
    u3  = 0.10e2 * cos((double) (2 * theta)) * (cos(k * r) * pow(k, 0.3e1) * pow(r, 0.3e1) - 0.6e1 * sin(k * r) * k * k * r * r - 0.15e2 * cos(k * r) * k * r + 0.15e2 * sin(k * r)) * pow(k, -0.5e1) * pow(r, -0.4e1) * phi * phi - 0.25e2 / 0.3e1 * cos((double) (2 * theta)) * (cos(k * r) * pow(k, 0.3e1) * pow(r, 0.3e1) - 0.6e1 * sin(k * r) * k * k * r * r - 0.15e2 * cos(k * r) * k * r + 0.15e2 * sin(k * r)) * pow(k, -0.5e1) * pow(r, -0.4e1) * pow(phi, 0.4e1);
    p = 0.3e1 * OmegaFact * r * r * cos((double) (2 * theta)) * phi * phi - OmegaFact * r * r * cos((double) (2 * theta)) * pow(phi, 0.4e1);
    } 
    else if( phip<eps )// small phip, big r
    { 
    u1 = 0.2e1 * pow(k, -0.5e1) * pow(r, -0.4e1) * (-0.3e1 * pow(k, 0.3e1) * pow(r, 0.5e1) * OmegaFact + 0.5e1 * cos(k * r) * pow(k, 0.3e1) * pow(r, 0.3e1) - 0.15e2 * sin(k * r) * k * k * r * r - 0.30e2 * cos(k * r) * k * r + 0.30e2 * sin(k * r)) * (cos(theta) * cos(0.2e1 * theta) + sin(theta) * sin(0.2e1 * theta)) * phip;
    u2 = -0.2e1 * pow(k, -0.5e1) * pow(r, -0.4e1) * (-0.3e1 * pow(k, 0.3e1) * pow(r, 0.5e1) * OmegaFact + 0.5e1 * cos(k * r) * pow(k, 0.3e1) * pow(r, 0.3e1) - 0.15e2 * sin(k * r) * k * k * r * r - 0.30e2 * cos(k * r) * k * r + 0.30e2 * sin(k * r)) * (sin((double) (2 * theta)) * cos((double) theta) - sin((double) theta) * cos((double) (2 * theta))) * phip;
    u3  = -0.10e2 * cos((double) (2 * theta)) * ((double) (-6 * r * r * k * k + 15) * sin((double) (k * r)) + cos((double) (k * r)) * (double) k * (double) r * (double) (r * r * k * k - 15)) * (double)   pow((double) k, (double) (-5)) * (double)   pow((double) r, (double) (-4)) * phip * phip + 0.25e2 / 0.3e1 * cos((double) (2 * theta)) * ((double) (-6 * r * r * k * k + 15) * sin((double) (k * r)) + cos((double) (k * r)) * (double) k * (double) r * (double) (r * r * k * k - 15)) * (double)   pow((double) k, (double) (-5)) * (double)   pow((double) r, (double) (-4)) * pow(phip, 0.4e1);
    p = 0.3e1 * OmegaFact * r * r * cos((double) (2 * theta)) * phip * phip - OmegaFact * r * r * cos((double) (2 * theta)) * pow(phip, 0.4e1);
    } 
    else 
    { // this case should not happen 
      OV_ABORT("ERROR"); 
    } 
  }
  else if( phi>eps && phip>eps )
  { // small r, big phi, big phip -- series error is O(r^3) 
    u1 = 0.2e1 * pow(k, -0.2e1) * pow(cos(phi) + 0.1e1, 0.2e1) * cos(theta) * (k * k + (double) (3 * OmegaFact)) * pow(cos(phi) - 0.1e1, 0.2e1) * pow(sin(phi), -0.3e1) * r;
    u2 = -0.2e1 * pow(k, -0.2e1) * pow(cos(phi) + 0.1e1, 0.2e1) * sin(theta) * (k * k + (double) (3 * OmegaFact)) * pow(cos(phi) - 0.1e1, 0.2e1) * pow(sin(phi), -0.3e1) * r;
    u3  = 0;
    p = (-0.3e1 * OmegaFact * pow(cos(phi), 0.2e1) + 0.3e1 * OmegaFact) * cos((double) (2 * theta)) * r * r;
  }
  else if( phi<eps )
  { // small r,  small phi 
    u1 = 0.2e1 * pow(k, -0.2e1) * cos(theta) * (k * k + (double) (3 * OmegaFact)) * r * phi;
    u2 = -0.2e1 * pow(k, -0.2e1) * sin(theta) * (k * k + (double) (3 * OmegaFact)) * r * phi;
    u3  = 0;
    p = 0.3e1 * OmegaFact * r * r * cos((double) (2 * theta)) * phi * phi - OmegaFact * r * r * cos((double) (2 * theta)) * pow(phi, 0.4e1);
  }
  else if( phip<eps )
  { // small r,  small phip 
    u1 = -0.2e1 * pow(k, -0.2e1) * cos(theta) * (k * k + (double) (3 * OmegaFact)) * r * phip;
    u2 = 0.2e1 * pow(k, -0.2e1) * sin(theta) * (k * k + (double) (3 * OmegaFact)) * r * phip;
    u3  = 0;
    p = 0.3e1 * OmegaFact * r * r * cos((double) (2 * theta)) * phip * phip - OmegaFact * r * r * cos((double) (2 * theta)) * pow(phip, 0.4e1);
  }
  else 
  { // this case should not happen 
    OV_ABORT("ERROR"); 
  } 
#endMacro

// ---- solid vibrational mode parameters ------
m=2, n=2;
cs=1.0000000000000000e+00; // Shear wave speed 
a=1.0000000000000000e+00;  // Sphere radius 
k=6.9879320005005200e+00;
omega=6.9879320005005200e+00;
OmegaFact=6.6906715092055407e-01;
