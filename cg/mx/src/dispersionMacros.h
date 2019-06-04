// 
//   Macros to adjust non-dispersive known solutions (e.g. scattering from a cylinder)
//   to the dispersive case 
//
//#define exDpw(x,y,t,dpwExp) sin(twoPi*(kx*(x)+ky*(y))-omegaDpwRe*(t))*(pwc[0]*(dpwExp))
//#define eyDpw(x,y,t,dpwExp) sin(twoPi*(kx*(x)+ky*(y))-omegaDpwRe*(t))*(pwc[1]*(dpwExp))
//#define hzDpw(x,y,t,dpwExp) sin(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc[5]

// ====================================================================================
/// Macro: Return the time-dependent coefficients for a known solution
// 
// NOTE: This next section is used in getInitialConditions.bC,
//        getErrors.bC and assignBoundaryConditions.bC 
// ====================================================================================
#beginMacro getKnownSolutionTimeCoefficients();
  real cost,sint,costm,sintm,dcost,dsint;
  real phiPc,phiPs, phiPcm,phiPsm;
  	  
  // define coefficients of the real(E), Im(E), Re(H) and Im(H)
  real cEr, cEi, cHr, cHi;  // coefficients for time t
  real cErm, cEim, cHrm, cHim;  // coefficients for time t-dt

  // coefficients of Polarization vector
  real cPr[10],cPi[10], cPrm[10], cPim[10];


  if( localDispersionModel==noDispersion )
  {
    //     --- NON-DISPERSIVE ---
    const real tm=t-dt;

    cost = cos(-twoPi*cc0*t); // *wdh* 040626 add "-"
    sint = sin(-twoPi*cc0*t); // *wdh* 040626 add "-"
  
    costm= cos(-twoPi*cc0*tm); // *wdh* 040626 add "-"
    sintm= sin(-twoPi*cc0*tm); // *wdh* 040626 add "-"
  
    // SOSUP needs the time derivative
    dcost =  twoPi*cc0*sint;  // d(cos(..))/dt 
    dsint = -twoPi*cc0*cost;  // d(sin(..))/dt 

    // new way: 
    //   E = Im( (Er + i*Ei)*exp(-i*omega*t) )
    //     = Im( (Er + i*Ei)*( cos(-i*omega*t) + i*sin(-i*omega*t) ) )
    //     = Er*sin(-i*omega*t) + Ei*
    cEr = sin(-twoPi*cc0*t);  cErm = sin(-twoPi*cc0*tm);  
    cEi = cos(-twoPi*cc0*t);  cEim = cos(-twoPi*cc0*tm);

    cHr = sin(-twoPi*cc0*t);  cHrm = sin(-twoPi*cc0*tm);
    cHi = cos(-twoPi*cc0*t);  cHim = cos(-twoPi*cc0*tm);
  }
  else
  {
    // -- DISPERSIVE MODEL -- *CHECK ME*
    const real tm=t-dt;
  
    // ---Evaluate the dispersion relation for "s" using outside grid  -------------------
    int outerGrid=0;  // use this grid *wdh* June 11, 2018
    DispersiveMaterialParameters & dmp = getDispersiveMaterialParameters(outerGrid);
    const real kk = twoPi*cc0;  // Parameter c*k in dispersion relation **check me**
    // *new way*
    real sr,si,chir[10],chii[10],chiSumr,chiSumi;
    dmp.evaluateDispersionRelation( cGrid(outerGrid),kk, sr, si, chir,chii,chiSumr,chiSumi ); 

    //  chi comes from the actual grid 
    DispersiveMaterialParameters & dmp2 = getDispersiveMaterialParameters(grid);
    real kr,ki;
    dmp2.evaluateComplexWaveNumber( cGrid(grid),sr,si, kr,ki,chir,chii,chiSumr,chiSumi ); // eval (kr,ki), chi given (sr,si)
    const real alphaP = dmp2.alphaP;

    // real reS, imS;
    // dmp.computeDispersionRelation( c,eps,mu,kk, reS, imS );
    // real expS = exp(reS*t), expSm=exp(reS*(t-dt));

    // si=-si;  // flip sign    **** FIX ME ****

    if( t<=3.*dt ) 
    {
      printF("--MX--GIC dispersion: s=(%12.4e,%12.4e) chi=(%12.4e,%12.4e)\n",sr,si,chir[0],chii[0]);
      printF("--MX--GIC grid=%i: complex k=(kr,ki)=(%e,%e)\n",grid,kr,ki);
      //printF("--MX--GIC scatCyl si/(twoPi*cc0)=%g\n",si/twoPi*cc0);
    }
    
    // Im( (Er+i*Ei)*( ct + i*st ) )
    //   st*Er + ct*Ei 
    real expt=exp(sr*t), exptm=exp(sr*tm);
    real ct =cos( si*t  )*expt,  st =sin( si*t )*expt;
    real ctm=cos( si*tm )*exptm, stm=sin( si*tm )*exptm;
    
    sint =  st;      // Coeff of Er
    cost =  ct;      // Coeff of Ei
  
    sintm =  stm;
    costm =  ctm;
  
    // SOSUP needs the time derivative **FIX ME**
    if( method==sosup )
    {
      OV_ABORT("finish me");
    }
    
    dcost =  -si*st + sr*ct;  //  d/dt of "cost" 
    dsint =  (si*ct + sr*st);  //  d/dt of "cost" 

  
    // P/eps = Im{ chi(s)*E } = Im{ (chir+i*chi)*( Er + i*Ei)( ct+i*st ) }
    //       =  Im{ (chir+i*chi)*( Er*ct- Ei*st + i( Er*st + Ei*ct ) )
    //       =  chir*(  Er*st + Ei*ct ) + chii*( Er*ct- Ei*st )
    //       =  (chir*st + chii*ct)*Er + (chir*ct-chii*st )*Ei 
  
    phiPs = chir[0]*sint+chii[0]*cost;  // Coeff of Er 
    phiPc = chir[0]*cost-chii[0]*sint;  // coeff of Ei
  	    
    phiPsm = chir[0]*sintm+chii[0]*costm;  // Coeff of Er(t-dt)
    phiPcm = chir[0]*costm-chii[0]*sintm;  // coeff of Ei(t-dt)

    // *new* 
    cEr = st;  cErm = stm;
    cEi = ct;  cEim = ctm;

    // *wdh* May 26, 2019  -- do not do this for dieletric disk/sphere solutions 
    if( !(knownSolutionOption==scatteringFromADielectricDiskKnownSolution ||
          knownSolutionOption==scatteringFromADielectricSphereKnownSolution ) ) 
    {
      // -- This computes the coefficient of H given E from
      //        H = Hat*e^{s t}
      //        Ehat = eta(s)/(s*mu0) * curl( Hhat)
      // 
      //     --> assumes E amplitude of incident wave is =1 I think
      // 

      // ----- Compute the coefficients of H for the scattering of a GDM wave from a cylinder ----
      // Coefficients generated by cg/mx/codes/gdm.maple

      real coeffH=1./kk;  // is this right? -- probably should be +/- (eps*cc*kk)
      
      // Coefficients generated by cg/mx/codes/gdm.maple
      // real chirSum=chir[0], chiiSum=chii[0];
      real chirSum=chiSumr, chiiSum=chiSumi;  // *wdh* May 7, 2019
      cHr=-coeffH*(((-1.*chirSum*sr+si*chiiSum)*ct+st*(si*chirSum+sr*chiiSum))*alphaP-1.*ct*sr+si*st);
      cHi=-coeffH*(((si*chirSum+sr*chiiSum)*ct+chirSum*sr*st-1.*si*st*chiiSum)*alphaP+ct*si+sr*st);
      cHrm=-coeffH*(((-1.*chirSum*sr+si*chiiSum)*ctm+stm*(si*chirSum+sr*chiiSum))*alphaP-1.*ctm*sr+si*stm);
      cHim=-coeffH*(((si*chirSum+sr*chiiSum)*ctm+chirSum*sr*stm-1.*si*stm*chiiSum)*alphaP+ctm*si+sr*stm);

    }
    else
    {
      cHr = st;  cHrm = stm;
      cHi = ct;  cHim = ctm;
    }
    if( numberOfPolarizationVectors>10 )
    {
      OV_ABORT("ERROR: FIX ME numberOfPolarizationVectors>10");
    }
    

    for( int jv=0; jv<numberOfPolarizationVectors; jv++ )
    {
      cPr[jv] = eps*( chir[jv]*sint+chii[jv]*cost);  // Coeff of Er for P 
      cPi[jv] = eps*( chir[jv]*cost-chii[jv]*sint);  // coeff of Eifor P 
  	    
      cPrm[jv] = eps*(chir[jv]*sintm+chii[jv]*costm);  // Coeff of Er(t-dt)
      cPim[jv] = eps*(chir[jv]*costm-chii[jv]*sintm);  // coeff of Ei(t-dt)
    }
    
    
  }
#endMacro
