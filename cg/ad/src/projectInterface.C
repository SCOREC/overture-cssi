#include "Cgad.h"
#include "AdParameters.h"
#include "ParallelUtility.h"
// #include "ParallelGridUtility.h"
#include "Interface.h"  

#define ForBoundary(side,axis)   for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) \
                                 for( int side=0; side<=1; side++ )

#define FOR_3D(i1,i2,i3,I1,I2,I3) \
int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  \
int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
for(int i3=I3Base; i3<=I3Bound; i3++) \
for(int i2=I2Base; i2<=I2Bound; i2++) \
for(int i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) \
I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  \
I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
for(int i3=I3Base; i3<=I3Bound; i3++) \
for(int i2=I2Base; i2<=I2Bound; i2++) \
for(int i1=I1Base; i1<=I1Bound; i1++)


/// --------------------------------------------------------------------------------------
/// \brief Project interface values such as Temperature for conjugate heat transfer.
/// --------------------------------------------------------------------------------------
int Cgad::
projectInterface(real t, real dt, GridFunction & cgf )
{


  const int & multiDomainProblem = parameters.dbase.get<int>("multiDomainProblem"); 

  const int projectInterfaceTemperature = parameters.dbase.get<bool>("projectInterfaceTemperature");

  printP("\n ##### Entering CgAd::projectInterface t=%9.3e projectInterfaceTemperature=%d ######\n",t,projectInterfaceTemperature);

  if( ! (multiDomainProblem && projectInterfaceTemperature ) )
  {
     return 0;
  }

  realCompositeGridFunction & u = cgf.u;

  CompositeGridOperators & cgop =  *cgf.u.getOperators();


  CompositeGrid & cg = *u.getCompositeGrid();
  const int numberOfDimensions = cg.numberOfDimensions();
  // int numberOfComponents = 1;
  Range N = parameters.dbase.get<Range >("Rt");   // time dependent variables
  
  const bool twilightZoneFlow    = parameters.dbase.get<bool>("twilightZoneFlow");
  const int & debug              = parameters.dbase.get<int>("debug");
  const int & myid               = parameters.dbase.get<int>("myid");
  // const int & orderOfAccuracy    = parameters.dbase.get<int>("orderOfAccuracy");
  // assert( orderOfAccuracy==2 || orderOfAccuracy==4 );

  if( debug & 1 )
    printP("+++++ Entering CgAd::projectInterface t=%9.3e +++++\n",t);


  // RealArray & champParameters = parameters.dbase.get<RealArray>("champParameters");

  const IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");

  // printF("\n -- CGAD-- champBoundaryConditions: multiDomainProblem=%d,  kappa=%g, thermalConductivity=%g\n",multiDomainProblem,kappa[0],thermalConductivity);

  
  Index Ibv[3], &Ib1=Ibv[0], &Ib2=Ibv[1], &Ib3=Ibv[2];
  // Index Igv[3], &Ig1=Igv[0], &Ig2=Igv[1], &Ig3=Igv[2];

  // int i1,i2,i3, j1,j2,j3, i1m,i2m,i3m, m1,m2,m3;
  // int isv[3], &is1=isv[0], &is2=isv[1], &is3=isv[2];

  // real maxInterfaceResidual=0.; // max residual over all interfaces 
  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {

    MappedGrid & mg = cg[grid];
    OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);

    // const int isRectangular=mg.isRectangular();

    // real dx[3]={1.,1.,1.};
    // if( isRectangular )
    //   mg.getDeltaX(dx);
    // else
    //   mg.update(MappedGrid::THEinverseVertexDerivative);


    // --- loop over faces ----      
    ForBoundary(side,axis)
    {
      if( interfaceType(side,axis,grid) == Parameters::heatFluxInterface && mg.boundaryCondition(side,axis)>0 )
      {
        // printF("+++++ projectInterface for (grid,side,axis)=(%d,%d,%d) (heatFluxInterface)\n",grid,side,axis);

        // // Set the index-shift for this side
        // is1=is2=is3=0;
        // isv[axis]=1-2*side;

        // const Real KLR = champParameters(2,side,axis,grid);    // K1/K2
        // const Real beta  = champParameters(3,side,axis,grid);    // D1/D2    

        // const Real & kThermalLeft  = parameters.dbase.get<Real>("thermalConductivity");
        // const Real & kThermalRight = kThermalLeft/KLR; 

       // Retrieve the parameters from the opposite side of the interface:
       Parameters & paramRight = getInterfaceParameters( grid, side, axis, parameters);
     
       // This only works if the opposite side is a Cgad object *** FIX ME ***
       std::vector<real> & kappaLeft  = parameters.dbase.get<std::vector<real> >("kappa");
       std::vector<real> & kappaRight = paramRight.dbase.get<std::vector<real> >("kappa");
       const Real & kThermalLeft  = parameters.dbase.get<Real>("thermalConductivity");
       const Real & kThermalRight = paramRight.dbase.get<Real>("thermalConductivity");


        if( debug & 1 )
          printP(" CgAd::projectInterface (grid,side,axis)=(%d,%d,%d) : kThermalLeft=%g, kThermalRight=%g\n",grid,side,axis,kThermalLeft,kThermalRight);

        OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal); 

        InterfaceData interfaceData;
        // Range Rx=numberOfDimensions;
        getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
        // getGhostIndex(   mg.gridIndexRange(),side,axis,Ig1,Ig2,Ig3);

        interfaceData.u.redim(Ib1,Ib2,Ib3,N); // TRight is returned here 
        interfaceData.t=t;
        interfaceData.u=0;
        // ::display(interfaceData.u(Ib1,Ib2,Ib3,N),"BEFORE REQUEST DATA FOR INTERFACE","%7.4f ");

        // GridFaceDescriptor & myGFD = getInterfaceGridFaceDescriptor( grid, side, axis, parameters );

        GridFaceDescriptor targetInfo(-1,grid,side,axis);
        targetInfo.a[0]=1.; targetInfo.a[1]=0; targetInfo.a[2]=0.; // get T = 1.*T + 0*T.n 

        int interfaceDataOptions = Parameters::heatFluxInterfaceData;
        bool saveTimeHistory=false;

        if( debug & 8 )
        {
          printF("\n");
          printF(" >>>> projectInterface:: REQUEST INTERFACE DATA for (grid,side,axis)=(%d,%d,%d) t=%9.3e.\n",
                 grid,side,axis,t);
        }

        getInterfaceData( t, grid, side, axis, 
                          interfaceDataOptions,
                          interfaceData.u,
                          parameters,saveTimeHistory,
                          &targetInfo );

        // ----- PROJECT THE TEMPERATURE -----
        // 
        //           Kleft*uLeft + Kright*uRight 
        //   uAve = -----------------------------
        //               Kleft + Kright 

        const Real alpha1 =  kThermalLeft/(kThermalLeft+kThermalRight);
        const Real alpha2 = kThermalRight/(kThermalLeft+kThermalRight);
        // RealArray uAve(Ib1,Ib2,Ib3,N);
        RealArray & uRight = interfaceData.u; 

        if( twilightZoneFlow  )
        {
          OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
          OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);    
          RealArray ue(Ib1,Ib2,Ib3,N);
          int rectangular=0;
          e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,N,t);        
          Real err = max(fabs(uLocal(Ib1,Ib2,Ib3,N)-ue));
          printP("projectInterface: t=%9.3e, error in interface T BEFORE projection=%9.3e\n",t,err);
        }
        else
        {
          
          // -- advection coefficients 
          std::vector<real> & aL = parameters.dbase.get<std::vector<real> >("a");
          std::vector<real> & bL = parameters.dbase.get<std::vector<real> >("b");
          std::vector<real> & cL = parameters.dbase.get<std::vector<real> >("c"); 
  
          std::vector<real> & aR = paramRight.dbase.get<std::vector<real> >("a");
          std::vector<real> & bR = paramRight.dbase.get<std::vector<real> >("b");
          std::vector<real> & cR = paramRight.dbase.get<std::vector<real> >("c");

          const Real DL = kappaLeft[0];
          const Real DR = kappaRight[0];
          const Real KL = kThermalLeft;
          const Real KR = kThermalRight;

          // Real u1L=aL[0]; 
          // Real u2L=bL[0]; 
          // Real u1R=aR[0]; 
          // Real u2R=bR[0];  
          // advection coefficinets scaled by 1/D 
          const Real u1DL=aL[0]/DL; 
          const Real u2DL=bL[0]/DL; 
          const Real u1DR=aR[0]/DR; 
          const Real u2DR=bR[0]/DR;                             

          RealArray uR, unR;
          uR = interfaceData.u; // save a copy

          // -- get K T.n ----
          // Note: flip sign of normal 
          targetInfo.a[0]=0.; targetInfo.a[1]=-KR; targetInfo.a[2]=0.; // get T = 0*T + kThermalRight*T.n 
          getInterfaceData( t, grid, side, axis, 
                           interfaceDataOptions,
                           interfaceData.u,
                           parameters,saveTimeHistory,
                           &targetInfo );
          unR = interfaceData.u;           

          Real jumpT =max(fabs(uLocal(Ib1,Ib2,Ib3,N)-uR));

          // Jump in flux:
          RealArray ux(Ib1,Ib2,Ib3,N), uy(Ib1,Ib2,Ib3,N);
          MappedGridOperators & op = cgop[grid];

          op.derivative(MappedGridOperators::xDerivative,uLocal,ux,Ib1,Ib2,Ib3,N);
          op.derivative(MappedGridOperators::yDerivative,uLocal,uy,Ib1,Ib2,Ib3,N);

          // ----- Flux = K*( T.n - (n.u/D) T ) ------

          OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);
          ux = KL*( normal(Ib1,Ib2,Ib3,0)*ux + normal(Ib1,Ib2,Ib3,1)*uy  - (normal(Ib1,Ib2,Ib3,0)*u1DL + normal(Ib1,Ib2,Ib3,1)*u2DL )*uLocal(Ib1,Ib2,Ib3,N) );

          unR -= KR*( (normal(Ib1,Ib2,Ib3,0)*u1DR + normal(Ib1,Ib2,Ib3,1)*u2DR )*uR );

          if( mg.numberOfDimensions()==3 )
          {
            OV_ABORT("finishe me -- 3D");
          }         

          Real jumpFlux = max(fabs(ux-unR)); // jump in flux 

          printP("projectInterface: t=%9.3e, BEFORE projection, [T]=%9.3e, [flux]=[K*( T.n - (n.u/D) T )]=%9.3e\n",t,jumpT,jumpFlux);

          // Real jumpT =max(fabs(uLocal(Ib1,Ib2,Ib3,N)-uRight));
          // printP("projectInterface: t=%9.3e, BEFORE projection, [T]=%9.3e\n",t,jumpT);

          uRight = uR; // reset for use below 



        }
        // Here is the projection: thermal conductivity weighted average
        uLocal(Ib1,Ib2,Ib3,N) = alpha1*uLocal(Ib1,Ib2,Ib3,N) + alpha2*uRight;
         
        if( twilightZoneFlow  )
        {
          OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
          OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);    
          RealArray ue(Ib1,Ib2,Ib3,N);
          int rectangular=0;
          e.gd( ue  ,xLocal,mg.numberOfDimensions(),rectangular,0,0,0,0,Ib1,Ib2,Ib3,N,t);
          
          if( true )
          {
            uLocal(Ib1,Ib2,Ib3,N) += (1.-alpha1)*ue(Ib1,Ib2,Ib3,N);
          }
          Real err = max(fabs(uLocal(Ib1,Ib2,Ib3,N)-ue));
          printP("projectInterface: t=%9.3e, error in interface T AFTER projection=%9.3e\n",t,err);
        }


      }
    }
  }
        







  return 1; // Temperature was projected
}