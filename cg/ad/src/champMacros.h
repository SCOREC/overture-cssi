// ====================================================================================
// Macro: Add twilght-zone corrections for CHAMP
// ====================================================================================
#beginMacro addTwilightZoneCorrectionForChamp()              

  OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));

  const RealArray & champParameters = parameters.dbase.get<RealArray>("champParameters");
  // const real Sl    = champParameters(0,side,axis,grid);
  // const real theta = champParameters(2,side,axis,grid);
  // const real beta  = champParameters(3,side,axis,grid);  
  const Real pl    = champParameters(0,side,axis,grid);    // optimized Scwartz Parameter for side 1
  const Real pr    = champParameters(1,side,axis,grid);    // optimized Scwartz Parameter for side 2
  const Real theta = champParameters(2,side,axis,grid);    // K1/K2
  const Real beta  = champParameters(3,side,axis,grid);    // D1/D2    
  const Real Sl    = champParameters(4,side,axis,grid);    // pl/dxs; 

  // printF("**** ADD TZ Correction to results from getData for CHAMP:  Sl=%g, theta=%g, beta=%g *********************************\n",Sl,theta,beta);                  

  OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);

  RealArray ue(Ib1,Ib2,Ib3,N), uex(Ib1,Ib2,Ib3,N), uey(Ib1,Ib2,Ib3,N), uexx(Ib1,Ib2,Ib3,N), ueyy(Ib1,Ib2,Ib3,N);
  int rectangular=0;
  e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,N,t);
  e.gd( uex ,xLocal,mg.numberOfDimensions(),rectangular,0,1,0,0,Ib1,Ib2,Ib3,N,t);
  e.gd( uey ,xLocal,mg.numberOfDimensions(),rectangular,0,0,1,0,Ib1,Ib2,Ib3,N,t);
  e.gd( uexx,xLocal,mg.numberOfDimensions(),rectangular,0,2,0,0,Ib1,Ib2,Ib3,N,t);
  e.gd( ueyy,xLocal,mg.numberOfDimensions(),rectangular,0,0,2,0,Ib1,Ib2,Ib3,N,t);


  if( isRectangular )
  {

    // const Real dxs = dx[axis];   // *** FIX ME need dx[opposite-side]
    // We use dx fro the opposite side
    const Real dxs = champParameters(5,side,axis,grid); // value saved in champBoundaryConditions.bC 

    const real a0 = Sl;
    const real a1 = theta + Sl*dxs*theta;

    // const real a2 = dxs*( beta      ) + Sl*( .5*SQR(dxs)*beta      );
    // const real a3 = dxs*( (beta-1.) ) + Sl*( .5*SQR(dxs)*(beta-1.) );
    // const real a4 = a3; // for 3D 

    real axx, ayy, azz;
    if( axis==0 )
    {
      axx = (beta   )*( dxs + Sl*( .5*SQR(dxs) ) );
      ayy = (beta-1.)*( dxs + Sl*( .5*SQR(dxs) ) );
      azz = ayy; // for 3D 
    }
    else if( axis==1 )
    {
      ayy = (beta   )*( dxs + Sl*( .5*SQR(dxs) ) );
      axx = (beta-1.)*( dxs + Sl*( .5*SQR(dxs) ) );
      azz = axx; // for 3D       
    }
    assert( numberOfDimensions==2 );

    int n=0;
    assert( N.getLength()==1 );
    interfaceData.u(Ib1,Ib2,Ib3,n) += a0*ue + a1*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + axx*uexx + ayy*ueyy;
  }
  else
  {
    // -------- CURVILINEAR GRID ---------


    const Real dxs = mg.gridSpacing(axis);

    const real a0 = Sl;
    const real a1 = theta + Sl*dxs*theta;


    const real a2 = dxs*( beta      ) + Sl*( .5*SQR(dxs)*beta      );
    const real a3 = dxs*( (beta-1.) ) + Sl*( .5*SQR(dxs)*(beta-1.) );
    const real a4 = a3; // for 3D 

    //  printF("**** ADD TZ Correction to results from getData for CHAMP:  Sl=%g, theta=%g, beta=%g **************** FINISH ME *****************\n",Sl,theta,beta); 

    // NCoeff(m) = theta*( an1L*xCoeff(m,i1,i2,i3) + an2L*yCoeff(m,i1,i2,i3) ) - b2R*r2Coeff(m,j1,j2,j3); 

    // // NOTE: Here we assume the tangential r derivatives are the same on both sides! (see champ4/notes) 
    // LCoeff(m) = (b1R*b1R/c11R)*( beta*lapCoeff(m,i1,i2,i3) 
    //                             - (c12R+c21R)*( theta*( (b1L/b1R)*r1r2Coeff(m,i1,i2,i3) + (b2L/b1R)*r2r2Coeff(m,j1,j2,j3) + b2Lb1Rr2*r2Coeff(m,j1,j2,j3) )
    //                                             + b2Rb1Rr2*r2Coeff(m,j1,j2,j3) + (b2R/b1R)*r2r2Coeff(m,j1,j2,j3) )
    //                             - c22R*r2r2Coeff(m,j1,j2,j3) 
    //                             + (c1R/b1R)*NCoeff(m) 
    //                             + c2R*r2Coeff(m,j1,j2,j3) );

    const Real dr1 = mg.gridSpacing(axis); 

    // -- compute tangential derivatives w.r.t r2  ---
    const int axisp1 = (axis+1) % mg.numberOfDimensions();
    const Real dr2 = mg.gridSpacing(axisp1); 

    int jsv[3]={0,0,0}, &js1=jsv[0], &js2=jsv[1], &js3=jsv[2];
    jsv[axisp1]=1; 

    RealArray uep(Ib1+js1,Ib2+js2,Ib3+js3,N), uem(Ib1-js1,Ib2-js2,Ib3-js3,N);
    e.gd( uep,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1+js1,Ib2+js2,Ib3+js3,N,t);
    e.gd( uem,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1-js1,Ib2-js2,Ib3-js3,N,t);

    RealArray uer2(Ib1,Ib2,Ib3,N), uer2r2(Ib1,Ib2,Ib3,N);
    uer2   = ( uep - uem )*(1./(2.*dr2));
    uer2r2 = ( uep -2.*ue + uem )*(1./(dr2*dr2));

    // -- compute the mixed derivative w.r.t r1 and r2 ----
    RealArray uepp(Ib1+1,Ib2+1,Ib3,N), uemm(Ib1-1,Ib2-1,Ib3,N), uepm(Ib1+1,Ib2-1,Ib3,N), uemp(Ib1-1,Ib2+1,Ib3,N);
    e.gd( uemm,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1-1,Ib2-1,Ib3,N,t);
    e.gd( uepm,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1+1,Ib2-1,Ib3,N,t);
    e.gd( uemp,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1-1,Ib2+1,Ib3,N,t);
    e.gd( uepp,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1+1,Ib2+1,Ib3,N,t);
    RealArray uer1r2(Ib1,Ib2,Ib3,N);
    uer1r2 = ( uepp - uemp -  uepm + uemm )*(1./(4.*dr1*dr2));   // D0r1 * D0r2 

    int n=0;
    assert( N.getLength()==1 );
    // interfaceData.u(Ib1,Ib2,Ib3,n) += a0*ue + a1*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + a2*uexx + a3*ueyy;

    // Some coefficients in the CHAMP conditions were saved when the matrix was created (champBoundaryConditions.bC)
    RealArray & cc = myGFD.dbase.get<RealArray>("ccChamp");

    // RealArray cc(Ib1,Ib2,Ib3,6); // need to save these when matrix is formed.  **FIX ME***
    // cc=0.;
    // cc(Ib1,Ib2,Ib3,1) = beta; 
    // cc(Ib1,Ib2,Ib3,3) = -1.; 
    // OV_ABORT("finish me");

    RealArray Nc(Ib1,Ib2,Ib3,N);
    Nc = theta*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + cc(Ib1,Ib2,Ib3,0)*uer2; 

    interfaceData.u(Ib1,Ib2,Ib3,n) += Sl*ue   
                                      + cc(Ib1,Ib2,Ib3,1)*( uexx + ueyy ) 
                                      + cc(Ib1,Ib2,Ib3,2)*uer1r2 
                                      + cc(Ib1,Ib2,Ib3,3)*uer2r2
                                      + cc(Ib1,Ib2,Ib3,4)*uer2
                                      + cc(Ib1,Ib2,Ib3,5)*Nc;


  }
#endMacro
