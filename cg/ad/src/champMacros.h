// ====================================================================================
// Macro: Add twilght-zone corrections for CHAMP
// ====================================================================================
#beginMacro addTwilightZoneCorrectionForChamp()              

  Index Ib1,Ib2,Ib3; 
  getBoundaryIndex(mg.indexRange(),side,axis,Ib1,Ib2,Ib3);

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

  const Real KLR = theta;
  const Real DLR = beta;

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

    if( orderOfAccuracy==2 )
    {

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
    else if( orderOfAccuracy==4 )
    {
      if( true )
        printF("Add TZ correction for CHAMP4 at t=%9.3e, KLR=%g\n",t,KLR);

      int n=0;
      assert( N.getLength()==1 );
      if( axis==0 )
      {
        RealArray uexxx(Ib1,Ib2,Ib3,N), uexyy(Ib1,Ib2,Ib3,N), uexxxx(Ib1,Ib2,Ib3,N), uexxyy(Ib1,Ib2,Ib3,N), ueyyyy(Ib1,Ib2,Ib3,N);
        e.gd( uexxx ,xLocal,mg.numberOfDimensions(),rectangular,0,3,0,0,Ib1,Ib2,Ib3,N,t);
        e.gd( uexyy ,xLocal,mg.numberOfDimensions(),rectangular,0,1,2,0,Ib1,Ib2,Ib3,N,t);
        e.gd( uexxxx,xLocal,mg.numberOfDimensions(),rectangular,0,4,0,0,Ib1,Ib2,Ib3,N,t);
        e.gd( uexxyy,xLocal,mg.numberOfDimensions(),rectangular,0,2,2,0,Ib1,Ib2,Ib3,N,t);
        e.gd( ueyyyy,xLocal,mg.numberOfDimensions(),rectangular,0,0,4,0,Ib1,Ib2,Ib3,N,t);

        RealArray L1(Ib1,Ib2,Ib3,N),L2(Ib1,Ib2,Ib3,N),L3(Ib1,Ib2,Ib3,N),L4(Ib1,Ib2,Ib3,N);

        const Real nSign = 2*side-1;
        
        L1 = (KLR*nSign)*uex;
        L2 = DLR*uexx + (DLR-1.)*ueyy;
        L3 = (KLR*DLR*nSign)*( uexxx ) + ((KLR*DLR-KLR)*nSign)*( uexyy );
        L4 = (DLR*DLR)*( uexxxx + 2.*uexxyy + ueyyyy ) - ueyyyy - (2.*DLR)*( uexxyy + ueyyyy) + 2.*ueyyyy;

        interfaceData.u(Ib1,Ib2,Ib3,n) += Sl*( ue + dxs*L1 + (dxs*dxs*.5)*L2 + (dxs*dxs*dxs/6.)*L3 + (dxs*dxs*dxs*dxs/24.)*L4 )
                                             + L1 + dxs*L2 + (dxs*dxs*.5)*L3 + (dxs*dxs*dxs/6.)*L4;

      }
      else if( axis==1 )
      {
        OV_ABORT("finish me");
      }
    }
    else
    {
      OV_ABORT("error -- orderOfAccuracy");
    }
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

    int n=0;
    assert( N.getLength()==1 );
    // interfaceData.u(Ib1,Ib2,Ib3,n) += a0*ue + a1*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + a2*uexx + a3*ueyy;

    // Some coefficients in the CHAMP conditions were saved when the matrix was created (champBoundaryConditions.bC)
    RealArray & cc = myGFD.dbase.get<RealArray>("ccChamp");

#defineMacro diffr2(u,i1,i2,i3,axis,dr)  ((u(i1+1,i2,i3)-u(i1-1,i2,i3))/(2.*dr[axis]))
#defineMacro diffs2(u,i1,i2,i3,axis,dr)  ((u(i1,i2+1,i3)-u(i1,i2-1,i3))/(2.*dr[axis]))
#defineMacro diffrr2(u,i1,i2,i3,axis,dr) ((u(i1+1,i2,i3)-2.*u(i1,i2,i3)+u(i1-1,i2,i3))/(SQR(dr[axis])))
#defineMacro diffss2(u,i1,i2,i3,axis,dr) ((u(i1,i2+1,i3)-2.*u(i1,i2,i3)+u(i1,i2-1,i3))/(SQR(dr[axis])))
#defineMacro diffr4(u,i1,i2,i3,axis,dr)  (((1./12.)*u(i1-2,i2,i3)-(2./3.)*u(i1-1,i2,i3)+(2./3.)*u(i1+1,i2,i3)-(1./12.)*u(i1+2,i2,i3))/(dr[axis]))
#defineMacro diffs4(u,i1,i2,i3,axis,dr)  (((1./12.)*u(i1,i2-2,i3)-(2./3.)*u(i1,i2-1,i3)+(2./3.)*u(i1,i2+1,i3)-(1./12.)*u(i1,i2+2,i3))/(dr[axis]))
#defineMacro diffrr4(u,i1,i2,i3,axis,dr)  ((-1.*u(i1-2,i2,i3)+16.*u(i1-1,i2,i3)-30.*u(i1,i2,i3)+16.*u(i1+1,i2,i3)-1.*u(i1+2,i2,i3))/(12.*SQR(dr[axis])))
#defineMacro diffss4(u,i1,i2,i3,axis,dr)  ((-1.*u(i1,i2-2,i3)+16.*u(i1,i2-1,i3)-30.*u(i1,i2,i3)+16.*u(i1,i2+1,i3)-1.*u(i1,i2+2,i3))/(12.*SQR(dr[axis])))
#defineMacro diffrrr2(u,i1,i2,i3,axis,dr)  (((-1./2.)*u(i1-2,i2,i3)+1.*u(i1-1,i2,i3)-1.*u(i1+1,i2,i3)+(1./2.)*u(i1+2,i2,i3))/(pow(dr[axis],3)))
#defineMacro diffsss2(u,i1,i2,i3,axis,dr)  (((-1./2.)*u(i1,i2-2,i3)+1.*u(i1,i2-1,i3)-1.*u(i1,i2+1,i3)+(1./2.)*u(i1,i2+2,i3))/(pow(dr[axis],3)))
#defineMacro diffrrrr2(u,i1,i2,i3,axis,dr)  ((u(i1-2,i2,i3)-4.*u(i1-1,i2,i3)+6.*u(i1,i2,i3)-4.*u(i1+1,i2,i3)+u(i1+2,i2,i3))/(pow(dr[axis],4)))
#defineMacro diffssss2(u,i1,i2,i3,axis,dr)  ((u(i1,i2-2,i3)-4.*u(i1,i2-1,i3)+6.*u(i1,i2,i3)-4.*u(i1,i2+1,i3)+u(i1,i2+2,i3))/(pow(dr[axis],4)))
#defineMacro diffrs4(u,i1,i2,i3,axis1,axis2,dr) ((64.*u(i1+1,i2+1,i3)-64.*u(i1+1,i2-1,i3)-8.*u(i1+1,i2+2,i3)+8.*u(i1+1,i2-2,i3)-64.*u(i1-1,i2+1,i3)+64.*u(i1-1,i2-1,i3)+8.*u(i1-1,i2+2,i3)-8.*u(i1-1,i2-2,i3)-8.*u(i1+2,i2+1,i3)+8.*u(i1+2,i2-1,i3)+u(i1+2,i2+2,i3)-u(i1+2,i2-2,i3)+8.*u(i1-2,i2+1,i3)-8.* u(i1-2,i2-1,i3)-u(i1-2,i2+2,i3)+u(i1-2,i2-2,i3))/(dr[axis2]*dr[axis1]*144.))

#defineMacro diffrs2(u,i1,i2,i3,axis1,axis2,dr) ((diffs2(u,i1+1,i2,i3,axis2,dr)-diffs2(u,i1-1,i2,i3,axis2,dr))/(2.*dr[axis1]))

#defineMacro diffrs4(u,i1,i2,i3,axis1,axis2,dr)   (((1./12.)*  diffs4(u,i1-2,i2,i3,axis2,dr)-(2./3.)*  diffs4(u,i1-1,i2,i3,axis2,dr)+(2./3.)*  diffs4(u,i1+1,i2,i3,axis2,dr)-(1./12.)*  diffs4(u,i1+2,i2,i3,axis2,dr))/(dr[axis1]))
#defineMacro diffrss4(u,i1,i2,i3,axis1,axis2,dr)  (((1./12.)* diffss4(u,i1-2,i2,i3,axis2,dr)-(2./3.)* diffss4(u,i1-1,i2,i3,axis2,dr)+(2./3.)* diffss4(u,i1+1,i2,i3,axis2,dr)-(1./12.)* diffss4(u,i1+2,i2,i3,axis2,dr))/(dr[axis1]))
#defineMacro diffrsss2(u,i1,i2,i3,axis1,axis2,dr) (((1./12.)*diffsss2(u,i1-2,i2,i3,axis2,dr)-(2./3.)*diffsss2(u,i1-1,i2,i3,axis2,dr)+(2./3.)*diffsss2(u,i1+1,i2,i3,axis2,dr)-(1./12.)*diffsss2(u,i1+2,i2,i3,axis2,dr))/(dr[axis1]))

#defineMacro diffrrs4(u,i1,i2,i3,axis1,axis2,dr)  ((-1.* diffs4(u,i1-2,i2,i3,axis2,dr)+16.* diffs4(u,i1-1,i2,i3,axis2,dr)-30.* diffs4(u,i1,i2,i3,axis2,dr)+16.* diffs4(u,i1+1,i2,i3,axis2,dr)-1.* diffs4(u,i1+2,i2,i3,axis2,dr))/(12.*SQR(dr[axis1])))
#defineMacro diffrrss4(u,i1,i2,i3,axis1,axis2,dr) ((-1.*diffss4(u,i1-2,i2,i3,axis2,dr)+16.*diffss4(u,i1-1,i2,i3,axis2,dr)-30.*diffss4(u,i1,i2,i3,axis2,dr)+16.*diffss4(u,i1+1,i2,i3,axis2,dr)-1.*diffss4(u,i1+2,i2,i3,axis2,dr))/(12.*SQR(dr[axis1])))

#defineMacro diffrrrs2(u,i1,i2,i3,axis1,axis2,dr) (((-1./2.)*diffs4(u,i1-2,i2,i3,axis2,dr)+1.*diffs4(u,i1-1,i2,i3,axis2,dr)-1.*diffs4(u,i1+1,i2,i3,axis2,dr)+(1./2.)*diffs4(u,i1+2,i2,i3,axis2,dr)/(pow(dr[axis1],3))))


    if( orderOfAccuracy==2 )
    {
      const RealArray & dra = mg.gridSpacing(); 
      real dr[3]={dra(0),dra(1),dra(2)};

      // evaluate the exact solution at points on and near the boundary
      Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
      I1=Ib1; I2=Ib2; I3=Ib3;
      for( int dir=0; dir<numberOfDimensions; dir++ )
      {
        Iv[dir] = Range( Iv[dir].getBase()-2, Iv[dir].getBound()+2 );  // we have a 5-pt stencil
      }

      RealArray ue(I1,I2,I3,N);
      e.gd( ue,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,I1,I2,I3,N,t);

      RealArray    uer(Ib1,Ib2,Ib3),    ues(Ib1,Ib2,Ib3);
      RealArray   uerr(Ib1,Ib2,Ib3),   uers(Ib1,Ib2,Ib3),   uess(Ib1,Ib2,Ib3);

      uer  = diffr2(ue,Ib1,Ib2,Ib3,0,dr);
      ues  = diffs2(ue,Ib1,Ib2,Ib3,1,dr);

      uerr = diffrr2(ue,Ib1,Ib2,Ib3,0,dr);
      uers = diffrs2(ue,Ib1,Ib2,Ib3,0,1,dr);
      uess = diffss2(ue,Ib1,Ib2,Ib3,1,dr);

      // interfaceData.u(Ib1,Ib2,Ib3,n) += Sl*ue(Ib1,Ib2,Ib3); 
      // interfaceData.u(Ib1,Ib2,Ib3,n) +=  cc(Ib1,Ib2,Ib3, 0)*ue(Ib1,Ib2,Ib3);
      interfaceData.u(Ib1,Ib2,Ib3,n) +=  cc(Ib1,Ib2,Ib3, 0)*ue(Ib1,Ib2,Ib3)
                                        +cc(Ib1,Ib2,Ib3, 1)*uer
                                        +cc(Ib1,Ib2,Ib3, 2)*ues
                                        +cc(Ib1,Ib2,Ib3, 3)*uerr
                                        +cc(Ib1,Ib2,Ib3, 4)*uers
                                        +cc(Ib1,Ib2,Ib3, 5)*uess;


      // RealArray uep(Ib1+js1,Ib2+js2,Ib3+js3,N), uem(Ib1-js1,Ib2-js2,Ib3-js3,N);
      // e.gd( uep,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1+js1,Ib2+js2,Ib3+js3,N,t);
      // e.gd( uem,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1-js1,Ib2-js2,Ib3-js3,N,t);

      // RealArray uer2(Ib1,Ib2,Ib3,N), uer2r2(Ib1,Ib2,Ib3,N);
      // uer2   = ( uep - uem )*(1./(2.*dr2));
      // uer2r2 = ( uep -2.*ue + uem )*(1./(dr2*dr2));

      // // -- compute the mixed derivative w.r.t r1 and r2 ----
      // RealArray uepp(Ib1+1,Ib2+1,Ib3,N), uemm(Ib1-1,Ib2-1,Ib3,N), uepm(Ib1+1,Ib2-1,Ib3,N), uemp(Ib1-1,Ib2+1,Ib3,N);
      // e.gd( uemm,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1-1,Ib2-1,Ib3,N,t);
      // e.gd( uepm,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1+1,Ib2-1,Ib3,N,t);
      // e.gd( uemp,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1-1,Ib2+1,Ib3,N,t);
      // e.gd( uepp,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1+1,Ib2+1,Ib3,N,t);
      // RealArray uer1r2(Ib1,Ib2,Ib3,N);
      // uer1r2 = ( uepp - uemp -  uepm + uemm )*(1./(4.*dr1*dr2));   // D0r1 * D0r2       
      // RealArray Nc(Ib1,Ib2,Ib3,N);
      // Nc = theta*( normal(Ib1,Ib2,Ib3,0)*uex + normal(Ib1,Ib2,Ib3,1)*uey ) + cc(Ib1,Ib2,Ib3,0)*uer2; 

      // interfaceData.u(Ib1,Ib2,Ib3,n) += Sl*ue   
      //                                   + cc(Ib1,Ib2,Ib3,1)*( uexx + ueyy ) 
      //                                   + cc(Ib1,Ib2,Ib3,2)*uer1r2 
      //                                   + cc(Ib1,Ib2,Ib3,3)*uer2r2
      //                                   + cc(Ib1,Ib2,Ib3,4)*uer2
      //                                   + cc(Ib1,Ib2,Ib3,5)*Nc;
    }
    else
    {


      const RealArray & dra = mg.gridSpacing(); 
      real dr[3]={dra(0),dra(1),dra(2)};

      // evaluate the exact solution at points on and near the boundary
      Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
      I1=Ib1; I2=Ib2; I3=Ib3;
      for( int dir=0; dir<numberOfDimensions; dir++ )
      {
        Iv[dir] = Range( Iv[dir].getBase()-2, Iv[dir].getBound()+2 );  // we have a 5-pt stencil
      }

      RealArray ue(I1,I2,I3,N);
      e.gd( ue,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,I1,I2,I3,N,t);

      RealArray    uer(Ib1,Ib2,Ib3),    ues(Ib1,Ib2,Ib3);
      RealArray   uerr(Ib1,Ib2,Ib3),   uers(Ib1,Ib2,Ib3),   uess(Ib1,Ib2,Ib3);
      RealArray  uerrr(Ib1,Ib2,Ib3),  uerrs(Ib1,Ib2,Ib3),  uerss(Ib1,Ib2,Ib3),  uesss(Ib1,Ib2,Ib3);
      RealArray uerrrr(Ib1,Ib2,Ib3), uerrrs(Ib1,Ib2,Ib3), uerrss(Ib1,Ib2,Ib3), uersss(Ib1,Ib2,Ib3), uessss(Ib1,Ib2,Ib3);

      uer  = diffr4(ue,Ib1,Ib2,Ib3,0,dr);
      ues  = diffs4(ue,Ib1,Ib2,Ib3,1,dr);

      uerr = diffrr4(ue,Ib1,Ib2,Ib3,0,dr);
      uers = diffrs4(ue,Ib1,Ib2,Ib3,0,1,dr);
      uess = diffss4(ue,Ib1,Ib2,Ib3,1,dr);

      uerrr = diffrrr2(ue,Ib1,Ib2,Ib3,0,dr);
      uerrs = diffrrs4(ue,Ib1,Ib2,Ib3,0,1,dr);
      uerss = diffrss4(ue,Ib1,Ib2,Ib3,0,1,dr);
      uesss = diffsss2(ue,Ib1,Ib2,Ib3,1,dr);

      uerrrr = diffrrrr2(ue,Ib1,Ib2,Ib3,0,dr);
      uerrrs = diffrrrs2(ue,Ib1,Ib2,Ib3,0,1,dr);
      uerrss = diffrrss4(ue,Ib1,Ib2,Ib3,0,1,dr);
      uersss = diffrsss2(ue,Ib1,Ib2,Ib3,0,1,dr);
      uessss = diffssss2(ue,Ib1,Ib2,Ib3,1,dr);

      // printF("TZ correction: champ4: Sl=%9.3e \n",Sl);
      // ::display(cc(Ib1,Ib2,Ib3, 0),"cc(Ib1,Ib2,Ib3, 0)","%9.3e");

      // printF("TZ correction: ||uer||=%8.2e\n",max(fabs(uer)));
      // printF("TZ correction: ||ues||=%8.2e\n",max(fabs(ues)));

      // printF("TZ correction: ||uerr||=%8.2e\n",max(fabs(uerr)));
      // printF("TZ correction: ||uers||=%8.2e\n",max(fabs(uers)));
      // printF("TZ correction: ||uess||=%8.2e\n",max(fabs(uess)));

      // printF("TZ correction: ||uerrr||=%8.2e\n",max(fabs(uerrr)));
      // printF("TZ correction: ||uerrs||=%8.2e\n",max(fabs(uerrs)));
      // printF("TZ correction: ||uerss||=%8.2e\n",max(fabs(uerss)));
      // printF("TZ correction: ||uesss||=%8.2e\n",max(fabs(uesss)));

      // printF("TZ correction: ||uerrrr||=%8.2e\n",max(fabs(uerrrr)));
      // printF("TZ correction: ||uerrrs||=%8.2e\n",max(fabs(uerrrs)));
      // printF("TZ correction: ||uerrss||=%8.2e\n",max(fabs(uerrss)));
      // printF("TZ correction: ||uersss||=%8.2e\n",max(fabs(uersss)));
      // printF("TZ correction: ||uessss||=%8.2e\n",max(fabs(uessss)));

      // interfaceData.u(Ib1,Ib2,Ib3,n) += Sl*ue(Ib1,Ib2,Ib3); 
      // interfaceData.u(Ib1,Ib2,Ib3,n) +=  cc(Ib1,Ib2,Ib3, 0)*ue(Ib1,Ib2,Ib3);
      interfaceData.u(Ib1,Ib2,Ib3,n) +=  cc(Ib1,Ib2,Ib3, 0)*ue(Ib1,Ib2,Ib3)
                                        +cc(Ib1,Ib2,Ib3, 1)*uer
                                        +cc(Ib1,Ib2,Ib3, 2)*ues
                                        +cc(Ib1,Ib2,Ib3, 3)*uerr
                                        +cc(Ib1,Ib2,Ib3, 4)*uers
                                        +cc(Ib1,Ib2,Ib3, 5)*uess
                                        +cc(Ib1,Ib2,Ib3, 6)*uerrr
                                        +cc(Ib1,Ib2,Ib3, 7)*uerrs
                                        +cc(Ib1,Ib2,Ib3, 8)*uerss
                                        +cc(Ib1,Ib2,Ib3, 9)*uesss
                                        +cc(Ib1,Ib2,Ib3,10)*uerrrr
                                        +cc(Ib1,Ib2,Ib3,11)*uerrrs
                                        +cc(Ib1,Ib2,Ib3,12)*uerrss
                                        +cc(Ib1,Ib2,Ib3,13)*uersss
                                        +cc(Ib1,Ib2,Ib3,14)*uessss;


    }

  }
#endMacro
