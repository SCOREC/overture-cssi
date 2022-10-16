// This file automatically generated from interfaceBoundaryConditions.bC with bpp.
#include "DomainSolver.h"
#include "Interface.h"
#include "ParallelUtility.h"
#include "ArrayEvolution.h"

// ****** NOTE: this base class implementation assumes that interface BC is for the Temperature only *****

// ===========================================================================
/// \brief Setup an interface boundary condition.
/// \details This function is used when solving the interface equations 
///           by iteration. It will setup the interface conditions that should be
///           used. For example, on a heatFlux interface the interface BC may 
///           be Dirichlet, Neumann or mixed. This choice is determined by cgmp.
/// \param info (input) : contains the info on which interface to set. 
// ===========================================================================
int
DomainSolver::
setInterfaceBoundaryCondition( GridFaceDescriptor & info )
{

    CompositeGrid & cg = gf[0].cg;
    
    const IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");

  // *** As a start, here is how we save the interface BC info ***
  //   interfaceCondition(side,axis,grid) = [dirichletInterface/neumannInterface]
  //
    if (!parameters.dbase.has_key("interfaceCondition")) 
    {
        parameters.dbase.put<IntegerArray>("interfaceCondition");
        IntegerArray & interfaceCondition = parameters.dbase.get<IntegerArray>("interfaceCondition");

        interfaceCondition.redim(2,3,cg.numberOfComponentGrids());
        interfaceCondition=Parameters::dirichletInterface;
    }
    
    IntegerArray & interfaceCondition = parameters.dbase.get<IntegerArray>("interfaceCondition");

    const int grid=info.grid, side=info.side, axis=info.axis;

    if( grid<0 || grid>=cg.numberOfComponentGrids() ||
            side<0 || side>1 || axis<0 || axis>=cg.numberOfDimensions() )
    {
        printF("DomainSolver::setInterfaceBoundaryCondition:ERROR: invalid values: (side,axis,grid)=(%i,%i,%i)\n",
                      side,axis,grid);
        OV_ABORT("DomainSolver::setInterfaceBoundaryCondition:ERROR");
    }

    if( interfaceType(side,axis,grid)==Parameters::heatFluxInterface )
    {
    // ****************************************
    // ********* Heat flux interface **********
    // ****************************************

        printP("DomainSolver::setInterfaceBC: set heat-flux interface BC=%i for (side,axis,grid)=(%i,%i,%i)"
                      " a0=%6.2f a1=%6.2f\n",
                      info.interfaceBC,side,axis,grid,info.a[0],info.a[1]);

        interfaceCondition(side,axis,grid)=info.interfaceBC;

    // from insp.C
#define mixedRHS(component,side,axis,grid)         bcData(component+numberOfComponents*(0),side,axis,grid)
#define mixedCoeff(component,side,axis,grid)       bcData(component+numberOfComponents*(1),side,axis,grid)
#define mixedNormalCoeff(component,side,axis,grid) bcData(component+numberOfComponents*(2),side,axis,grid)

        RealArray & bcData = parameters.dbase.get<RealArray>("bcData");
        const int & numberOfComponents = parameters.dbase.get<int >("numberOfComponents");
        if( bcData.getLength(0)<3*numberOfComponents || bcData.getLength(1)!=2 )
        {
            bcData.display("error");
            OV_ABORT("error");
        }

        const int tc = parameters.dbase.get<int >("tc");   
        assert( tc>=0 );  

        mixedCoeff(tc,side,axis,grid)=info.a[0];
        mixedNormalCoeff(tc,side,axis,grid)=info.a[1];
        
    }
    else if( interfaceType(side,axis,grid)==Parameters::tractionInterface )
    {
    // *******************************************
    // ********** Traction Interface *************
    // *******************************************

        printP("DomainSolver::setInterfaceBC: do nothing for interfaceType(%i,%i,%i)=%i \n",side,axis,grid,
                      interfaceCondition(side,axis,grid));
    }
    else
    {
        printP("DomainSolver::setInterfaceBC: do nothing for interfaceType(%i,%i,%i)=%i \n",side,axis,grid,
                      interfaceCondition(side,axis,grid));
    }

}

// ===================================================================================
/// \brief Return the interface data required for a given type of interface.
/// \param info (input) : the descriptor for the interface.
/// \param interfaceDataOptions (output) : a list of items from Parameters::InterfaceDataEnum that define
///                    which data to get (or which data were set).  Multiple items are
///                     chosen by bit-wise or of the different options
/// \note: this function should be over-loaded.
// ===================================================================================
int DomainSolver::
getInterfaceDataOptions( GridFaceDescriptor & info, int & interfaceDataOptions ) const
{

    const int grid=info.grid, side=info.side, axis=info.axis;

    IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");
    if( grid<0 || grid>interfaceType.getBound(2) ||
            side<0 || side>1 || axis<0 || axis>interfaceType.getBound(1) )
    {
        printP("DomainSolver::getInterfaceDataOptions:ERROR: invalid values: (grid,side,axis)=(%i,%i,%i)\n",
                      grid,side,axis);
        OV_ABORT("DomainSolver::getInterfaceDataOptions:ERROR");
    }
    if( interfaceType(side,axis,grid)==Parameters::heatFluxInterface )
    {
        interfaceDataOptions=Parameters::heatFluxInterfaceData;
    }
    else if( interfaceType(side,axis,grid)==Parameters::tractionInterface ) 
    {
    // *** do this for now ** this function should be over-loaded by Cgad, Cgins, Cgcns, and Cgsm
        if( getClassName()=="Cgins" || getClassName()=="Cgcns" )
        {
            interfaceDataOptions=Parameters::positionInterfaceData;
        }
        else if( getClassName()=="Cgsm" )
        {
            interfaceDataOptions=Parameters::tractionInterfaceData;
        }
        else
        {
            printP("DomainSolver::getInterfaceDataOptions:ERROR: unknown class using tractionInterface\n");
            OV_ABORT("DomainSolver::getInterfaceDataOptions:ERROR");
        }
    }
    else
    {
        printP("DomainSolver::getInterfaceDataOptions:ERROR: interfaceType(grid=%i,side=%i,axis=%i)=%i\n",
                      grid,side,axis,interfaceType(side,axis,grid));
        OV_ABORT("DomainSolver::getInterfaceDataOptions:ERROR");
    }
    
    

    return 0;
}

// =====================================================================================================
// Macro: Evaluate the heatflux condition in the twilight-zone
//       ue <- a[0]*ue + a1*n.grad( ue )
// K1,K2,K3 : evaluate the exact solution at x(K1,K2,K3,*)
// N1,N2,N3 : evaluate the normal at normal(N1,N2,N3,*)
// =====================================================================================================


// =====================================================================================================
// Macro: Evaluate the heatflux condition
//       f = a[0]*u + a1*n.grad( u )
// =====================================================================================================




// =====================================================================================================
// Macro: Evaluate the CHAMP right-hand-side
//       f = a[0]*u + a[1]*n.grad( u )  at x=+dx on right side, or -dx, left side
// NOTE:
//   Evaluate the solution and it's derivatives at one grid line inside
//   Evaluate the normal on the boundary  
// =====================================================================================================

// ========================================================================================================
/// \brief Set or get the right-hand-side for an interface boundary condition.
/// \details This function is used when solving the interface equations 
///           by iteration.
/// \param option (input) : option=getInterfaceRightHandSide : get the RHS, 
///                         option=setInterfaceRightHandSide : set the RHS
/// \param interfaceDataOptions (input) : a list of items from Parameters::InterfaceDataEnum that define
////                    which data to get (or which data were set).  Multiple items are
///                     chosen by bit-wise or of the different options   
/// \param info (input) : contains the GridFaceDescriptor info used to set the right-hand-side.
/// \param gfd (input) : the master GridFaceDescriptor. 
/// \param gfIndex (input) : use the solution from gf[gfIndex]
/// \param t (input) : current time.
/// \param saveTimeHistory (input) : if true, save a time-history of the requested data. This is the
///    new way to save a time-history when interfaceCommunicationMode==requestInterfaceDataWhenNeeded
// ==========================================================================================================
int
DomainSolver::
interfaceRightHandSide( InterfaceOptionsEnum option, 
                                                int interfaceDataOptions,
                                                GridFaceDescriptor & info, 
                                                GridFaceDescriptor & gfd,
                                                int gfIndex, real t,
                                                bool saveTimeHistory /* = false */ )
{
  // gfIndex <0 means choose the grid function corresponding to time t

  // *** THIS NEXT CODE IS DUPLICATED IN interfaceRightHandSide in CGINS AND CGSM **FIX ME**
    bool getValuesByTimeExtrapolation = false;
    if( gfIndex==-1 )
    {
    // Find the solution that matches time=t
        const int & currentGF = parameters.dbase.get<int>("currentGF");
        const int & nextGF    = parameters.dbase.get<int>("nextGF");

        assert( current>=0 );
        if( gf[current].t == t )
        {
            gfIndex=current;
        }
        else if( currentGF<0 )   // do this for now
        {
            gfIndex=current;
            printP("interfaceRightHandSide:WARNING cannot find gfIndex to match t=%9.3e, using current...\n",t);
        }
        else
        {
        
            if( !(currentGF>=0 && nextGF>=0) )
            {
                printP("interfaceRightHandSide:ERROR: t=%9.2e, current=%i gf[current].t=%9.2e, currentGF=%i, nextGF=%i\n",
                              t,current,gf[current].t,currentGF,nextGF);
                OV_ABORT("FIX ME");
            }

            if( gf[currentGF].t == t )
            {
                gfIndex=currentGF;
            }
            else if( gf[nextGF].t == t )
            {
                gfIndex=nextGF;
            }
            else 
            {
        // ************** FIX ME ************
                printP("interfaceRightHandSide:WARNING cannot find gfIndex to match t=%9.3e\n"
                              "      currentGF=%i, gf[currentGF].t=%9.3e, nextGF=%i, gf[nextGF].t=%9.3e, saveTimeHistory=%d\n",
                              t,currentGF,gf[currentGF].t,nextGF,gf[nextGF].t,(int)saveTimeHistory);
                if( fabs(gf[currentGF].t-t) <  fabs(gf[nextGF].t-t) )
                    gfIndex=currentGF; 
                else
                    gfIndex=nextGF; 

                getValuesByTimeExtrapolation=true;
                if( !gfd.dbase.has_key("heatFluxHistory") )
                    getValuesByTimeExtrapolation=false;       // there is no time history yet

        // OV_ABORT("fix me");
            }
        }
        if( debug() & 4 )
        {
            printF("\n");
            printP(" >>>>> INFO: DomainSolver::interfaceRHS: Setting gfIndex=%d, gf[gfIndex].t=%9.3e, t=%9.3e\n\n",gfIndex,gf[gfIndex].t,t);
        }
    }


    CompositeGrid & cg = gf[gfIndex].cg;
    const int numberOfDimensions = cg.numberOfDimensions();
    
    const IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");

    const int grid=info.grid, side=info.side, axis=info.axis;

    if( grid<0 || grid>=cg.numberOfComponentGrids() ||
            side<0 || side>1 || axis<0 || axis>=cg.numberOfDimensions() )
    {
        printP("DomainSolver::interfaceRightHandSide:ERROR: invalid values: (grid,side,axis)=(%i,%i,%i)\n",
                      grid,side,axis);
        printP("gfIndex=%d, numberOfComponentGrids=%d, numberOfDimensions=%d\n",gfIndex,cg.numberOfComponentGrids(),cg.numberOfDimensions());
        OV_ABORT("DomainSolver::interfaceRightHandSide:ERROR");
    }

    MappedGrid & mg = cg[grid];
    RealArray & bd = parameters.getBoundaryData(side,axis,grid,mg);

    const int numberOfComponents = parameters.dbase.get<int >("numberOfComponents") - parameters.dbase.get<int >("numberOfExtraVariables");
    const bool & twilightZoneFlow = parameters.dbase.get<bool >("twilightZoneFlow");

    assert( info.u != NULL );
    RealArray & f = *info.u;

  // *wdh* April 8, 2021 -- This was re-worked to support new communicationMode
  // Index I1=f.dimension(0),I2=f.dimension(1),I3=f.dimension(2);
    Index J1=f.dimension(0), J2=f.dimension(1), J3=f.dimension(2); // target(get) or source (set) index bounds

    Index I1,I2,I3;                                                // source (get) or source(set) index bounds
    getBoundaryIndex(mg.gridIndexRange(),side,axis,I1,I2,I3);

    if( I1.getLength()!=J1.getLength() || I2.getLength()!=J2.getLength() || I3.getLength()!=J3.getLength() )
    {
        printP("interfaceRightHandSide:ERROR: (grid,side,axis) = (%d,%d,%d) \n"
                      " source and target Index ranges are not conformable:\n"
                      "   I1=[%d,%d] I2=[%d,%d] I3=[%d,%d] : source (get) or source(set) \n" 
                      "   J1=[%d,%d] J2=[%d,%d] J3=[%d,%d] : target (get) or source (set) \n",
                      grid,side,axis,
                      I1.getBase(),I1.getBound(), I2.getBase(),I2.getBound(), I3.getBase(),I3.getBound(),
                      J1.getBase(),J1.getBound(), J2.getBase(),J2.getBound(), J3.getBase(),J3.getBound()
                      );
        OV_ABORT("ERROR FIX ME BILL!");
    }

    if( interfaceType(side,axis,grid)==Parameters::heatFluxInterface )
    {
    // ****************************************
    // ********* Heat flux interface **********
    // ****************************************

        real *a = info.a;

        if( debug() & 4 )
        {
            printP("DomainSolver::interfaceRHS:heatFlux %s RHS for (side,axis,grid)=(%i,%i,%i) a=[%5.2f,%5.2f]"
                          " t=%9.3e gfIndex=%i (current=%i)\n",
                          (option==0 ? "get" : "set"),side,axis,grid,a[0],a[1],t,gfIndex,current);
        }

        const int tc = parameters.dbase.get<int >("tc");   
        assert( tc>=0 );
        Range N(tc,tc);

    // We could optimize this for rectangular grids 
        mg.update(MappedGrid::THEvertexBoundaryNormal);
#ifdef USE_PPP
        const realSerialArray & normal = mg.vertexBoundaryNormalArray(side,axis);
#else
        const realSerialArray & normal = mg.vertexBoundaryNormal(side,axis);
#endif


        if( option==setInterfaceRightHandSide )
        {
      // ***************************************************************************
      // ********************** SET HEAT-FLUX DATA *********************************
      // ***************************************************************************

      // **** set the RHS *****

            bd(I1,I2,I3,tc)=f(J1,J2,J3,tc);

            if( false )
                ::display(bd(I1,I2,I3,tc)," RHS values","%4.2f ");

            if(  twilightZoneFlow ) // turn this off for testing the case where the same TZ holds across all domains
            {
        // ---add forcing for twilight-zone flow---
        //   ue <- a[0]*ue + a1*n.grad( ue )
                    OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                    const bool isRectangular = false; // ** do this for now ** mg.isRectangular();
                    if( !isRectangular )
                        mg.update(MappedGrid::THEcenter);
                    realArray & x= mg.center();
                    #ifdef USE_PPP
                        realSerialArray xLocal; 
                        if( !isRectangular ) 
                            getLocalArrayWithGhostBoundaries(x,xLocal);
                    #else
                        const realSerialArray & xLocal = x;
                    #endif
                    realSerialArray ue(I1,I2,I3,N);
                    if( a[0]!=0. )
                    {
                        e.gd( ue ,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,N,t);  // exact solution 
                        ue(I1,I2,I3,N) = a[0]*ue(I1,I2,I3,N);
                    }
                    else
                    {
                        ue(I1,I2,I3,N) =0.;
                    }
                    if( a[1]!=0. )
                    {
                        realSerialArray uex(I1,I2,I3,N), uey(I1,I2,I3,N);
                        e.gd( uex ,xLocal,numberOfDimensions,isRectangular,0,1,0,0,I1,I2,I3,N,t);
                        e.gd( uey ,xLocal,numberOfDimensions,isRectangular,0,0,1,0,I1,I2,I3,N,t);
                        if( numberOfDimensions==2 )
                        {
                            ue(I1,I2,I3,N) += a[1]*( normal(I1,I2,I3,0)*uex + normal(I1,I2,I3,1)*uey );
                        }
                        else
                        {
                            realSerialArray uez(I1,I2,I3,N);
                            e.gd( uez ,xLocal,numberOfDimensions,isRectangular,0,0,0,1,I1,I2,I3,N,t);
                            ue(I1,I2,I3,N) += a[1]*( normal(I1,I2,I3,0)*uex + normal(I1,I2,I3,1)*uey + normal(I1,I2,I3,2)*uez ); 
                        }
                    }

        //   add on TZ flow:
        //   bd <- bd + a[0]*ue + a[1]*( nu*ue.n )
                bd(I1,I2,I3,tc) += ue(I1,I2,I3,N);

                if( false )
                    bd(I1,I2,I3,tc) = ue(I1,I2,I3,N);
            }      
            
        }
        else if( option==getInterfaceRightHandSide )
        {

      // ***************************************************************************
      // ********************** GET HEAT-FLUX DATA *********************************
      // ***************************************************************************

      // **** get the RHS ****

            realMappedGridFunction & u = gf[gfIndex].u[grid];
            OV_GET_SERIAL_ARRAY(real,gf[gfIndex].u[grid],uLocal); 
      // #ifdef USE_PPP
      //       realSerialArray uLocal; getLocalArrayWithGhostBoundaries(u,uLocal);
      // #else
      //       realSerialArray & uLocal = gf[gfIndex].u[grid];
      // #endif

            bool getChampData=false;        
            if( option==getInterfaceRightHandSide && info.dbase.has_key("getChampData") )  
            {
                    getChampData = info.dbase.get<int>("getChampData");
            }    
            if( getChampData )
            {
        // --- get interface RHS for the CHAMP CHT conditions  -----
                if( getValuesByTimeExtrapolation )
                {
          // const int orderOfAccuracyInSpace = parameters.dbase.get<int>("orderOfAccuracyInSpace");
                    ArrayEvolution & heatFluxHistory = gfd.dbase.get<ArrayEvolution>("heatFluxHistory");

          // const int orderOfAccuracy = parameters.dbase.get<int>("orderOfAccuracy");
                    const int orderOfTimeAccuracy = parameters.dbase.get<int >("orderOfTimeAccuracy");

                    const int orderOfTimeExtrapolation = orderOfTimeAccuracy+1; // Use this many time-levels to extrapolate the time data
                    const int numTimeLevelsInHistory = heatFluxHistory.getNumberOfTimeLevels();
                    printP(" >>>> interfaceRightHandSide: EVAL CHAMP RHS USING TIME EXTRAPOLATION TO t=%9.3e, orderOfTimeExtrapolation=%d, orderOfTimeAccuracy=%d, numTimeLevelsInHistory=%d <<<<<\n\n",
                                t,orderOfTimeExtrapolation,orderOfTimeAccuracy,numTimeLevelsInHistory);

                    const int numberOfTimeDerivatives=0;
                    heatFluxHistory.eval( t, f, numberOfTimeDerivatives, orderOfTimeExtrapolation );

                    saveTimeHistory=false; // do NOT save this value in the time history  
                }
                else
                {
                        if( debug() & 4 )
                            printF("DomainSolver::interfaceRightHandSide:CHAMP RHS: S = a[0]=%g,   a[1]=%g (= (+/-) 1 ?) \n",a[0],a[1]);
                        const bool isRectangular=mg.isRectangular();
            // --- CHAMP RHS is evaluated at a distance h off the boundary 
            // First line inside: 
                        Index K1,K2,K3;
                        getGhostIndex(mg.gridIndexRange(),side,axis,K1,K2,K3,-1);
            // printF(" axis=%d, side=%d : boundary (I1,I2)=[%d,%d][%d,%d], First-line inside (K1,K2)=[%d,%d][%d,%d]\n",
            //           axis,side, I1.getBase(),I1.getBound(), I2.getBase(),I2.getBound(), K1.getBase(),K1.getBound(), K2.getBase(),K2.getBound() );
                        if( true || isRectangular )  // use this always *wdh* April 12, 2022
                        {
              // -------------------------- RECTANGULAR --------------
                            real dx[3]={1.,1.,1.};
                            mg.getDeltaX(dx);
                            f(J1,J2,J3,tc) = a[0]*uLocal(K1,K2,K3,tc);
                            if( a[1]!=0. )
                            {
                // add on a[1]*( nu*u.n ) on the boundary 
                // **be careful** -- the normal changes sign on the two sides of the interface ---
                                MappedGridOperators & op = *(u.getOperators());
                                realSerialArray ux(K1,K2,K3,N), uy(K1,K2,K3,N);
                                op.derivative(MappedGridOperators::xDerivative,uLocal,ux,K1,K2,K3,N);
                                op.derivative(MappedGridOperators::yDerivative,uLocal,uy,K1,K2,K3,N);
                                if( cg.numberOfDimensions()==2 )
                                {
                                    f(J1,J2,J3,tc) += a[1]*( normal(I1,I2,I3,0)*ux + normal(I1,I2,I3,1)*uy );
                                }
                                else
                                {
                                    realSerialArray uz(K1,K2,K3);
                                    op.derivative(MappedGridOperators::zDerivative,uLocal,uz,K1,K2,K3,N);
                                    f(J1,J2,J3,tc) += a[1]*( normal(I1,I2,I3,0)*ux + normal(I1,I2,I3,1)*uy + normal(I1,I2,I3,2)*uz );
                                }
                            }
                            if( debug() & 8 )
                            {
                                printP("interfaceRightHandSide:CHAMP-RHS (grid,side,axis) = (%d,%d,%d) \n"
                                  "   I1=[%d,%d] I2=[%d,%d] I3=[%d,%d] : source (get) or source(set) \n" 
                                  "   J1=[%d,%d] J2=[%d,%d] J3=[%d,%d] : target (get) or source (set) \n",
                                  grid,side,axis,
                                  I1.getBase(),I1.getBound(), I2.getBase(),I2.getBound(), I3.getBase(),I3.getBound(),
                                  J1.getBase(),J1.getBound(), J2.getBase(),J2.getBound(), J3.getBase(),J3.getBound()
                                  );
                                ::display(f(J1,J2,J3,tc) ,sPrintF("interfaceRHS: CHAMP RHS(h):  %f*u + %f*( u.n ) (before TZ corrections)",a[0],a[1]));
                            }
                            if( twilightZoneFlow ) // turn this off for testing the case where the same TZ holds across all domains
                            {
                // ---add forcing for twilight-zone flow---
                //   ue <- a[0]*ue + a1*n.grad( ue )
                                    OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                                    const bool isRectangular = false; // ** do this for now ** mg.isRectangular();
                                    if( !isRectangular )
                                        mg.update(MappedGrid::THEcenter);
                                    realArray & x= mg.center();
                                    #ifdef USE_PPP
                                        realSerialArray xLocal; 
                                        if( !isRectangular ) 
                                            getLocalArrayWithGhostBoundaries(x,xLocal);
                                    #else
                                        const realSerialArray & xLocal = x;
                                    #endif
                                    realSerialArray ue(K1,K2,K3,N);
                                    if( a[0]!=0. )
                                    {
                                        e.gd( ue ,xLocal,numberOfDimensions,isRectangular,0,0,0,0,K1,K2,K3,N,t);  // exact solution 
                                        ue(K1,K2,K3,N) = a[0]*ue(K1,K2,K3,N);
                                    }
                                    else
                                    {
                                        ue(K1,K2,K3,N) =0.;
                                    }
                                    if( a[1]!=0. )
                                    {
                                        realSerialArray uex(K1,K2,K3,N), uey(K1,K2,K3,N);
                                        e.gd( uex ,xLocal,numberOfDimensions,isRectangular,0,1,0,0,K1,K2,K3,N,t);
                                        e.gd( uey ,xLocal,numberOfDimensions,isRectangular,0,0,1,0,K1,K2,K3,N,t);
                                        if( numberOfDimensions==2 )
                                        {
                                            ue(K1,K2,K3,N) += a[1]*( normal(I1,I2,I3,0)*uex + normal(I1,I2,I3,1)*uey );
                                        }
                                        else
                                        {
                                            realSerialArray uez(K1,K2,K3,N);
                                            e.gd( uez ,xLocal,numberOfDimensions,isRectangular,0,0,0,1,K1,K2,K3,N,t);
                                            ue(K1,K2,K3,N) += a[1]*( normal(I1,I2,I3,0)*uex + normal(I1,I2,I3,1)*uey + normal(I1,I2,I3,2)*uez ); 
                                        }
                                    }
                //   subtract off TZ flow:
                //   f <- f - ( a[0]*ue + a[1]*( nu*ue.n ) )
                                if( false )
                                {
                                    ::display(f(J1,J2,J3,tc) ," a[0]*u + a[1]*( k u.n )");
                                    ::display(ue(K1,K2,K3,tc)," a[0]*ue + a[1]*( k ue.n )");
                                }
                                f(J1,J2,J3,tc) -= ue(K1,K2,K3,N);
                                if( debug() & 8  )
                                {
                                    printP("RHS FOR CHAMP AFTER ADDING TZ CORRECTION, t=%9.3e\n",t);
                                    ::display(f(J1,J2,J3,tc) ," a[0]*u + a[1]*( k u.n ) - [a[0]*ue + a[1]*( k ue.n )]");
                                    e.gd( ue ,xLocal,numberOfDimensions,isRectangular,0,0,0,0,K1,K2,K3,N,t);  // exact solution 
                                    ::display(ue(K1,K2,K3,tc) ," ue on the first line inside");
                                    ::display(uLocal(K1,K2,K3,tc) ," uLocal on the first line inside");
                                }
                            }
                        }
                        else
                        {
               // ---------------- CURVILINEAR CASE --------------------
               // CHAMP RHS: (see champ4/notes)
               //   ( SL + b1^R * D_{ r_1 } ) Tr( (+/-)dr_1 )
               // b1 = an1*rx + an2*ry + an3*rz 
                            OV_GET_SERIAL_ARRAY( real,mg.inverseVertexDerivative(),rx );
              // printF(" rx : [%d,%d,%d,%d]\n",rx.getLength(0),rx.getLength(1),rx.getLength(2),rx.getLength(3));
              // The rx array is only 4d, here is a macro to make it look 5d 
                            #define RX(i1,i2,i3,m1,m2) rx(i1,i2,i3,(m1)+numberOfDimensions*(m2))
              // --- b1 is evaluated on the boundary ---
                            RealArray b1(I1,I2,I3); 
                            if( numberOfDimensions==2 )
                                b1 = normal(I1,I2,I3,0)*RX(I1,I2,I3,axis,0) + normal(I1,I2,I3,1)*RX(I1,I2,I3,axis,1);
                            else
                                b1 = normal(I1,I2,I3,0)*RX(I1,I2,I3,axis,0) + normal(I1,I2,I3,1)*RX(I1,I2,I3,axis,1) + normal(I1,I2,I3,2)*RX(I1,I2,I3,axis,2);
              // Delta index in the axis direction: 
                            int idv[3]={0,0,0}, &id1=idv[0], &id2=idv[1], &id3=idv[2];
                            idv[axis]=1;  // NOTE: We always compute D0r 
                            const Real dr = mg.gridSpacing(axis); 
                            const Real factor = (a[1]/(2.*dr));   // DO include a[1] since we really want the normal from the other side
              // const Real factor = (1./(2.*dr));   // Do not include the change of sign in a[1] for the curvilinear case
              // Champ RHS is centered at the first grid line: 
                            f(J1,J2,J3,tc) = a[0]*uLocal(K1,K2,K3,tc) + b1*( uLocal(K1+id1,K2+id2,K3+id3,tc) - uLocal(K1-id1,K2-id2,K3-id3,tc) )*factor;
                            if( twilightZoneFlow ) 
                            {
                // ---add forcing for twilight-zone flow---
                                OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                                OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
                                realSerialArray ue(K1,K2,K3,N), uep(K1+id1,K2+id2,K3+id3,N), uem(K1-id1,K2-id2,K3-id3,N);
                                e.gd( uem ,xLocal,numberOfDimensions,isRectangular,0,0,0,0,K1-id1,K2-id2,K3-id3,N,t);  //  T(0)   
                                e.gd( ue  ,xLocal,numberOfDimensions,isRectangular,0,0,0,0,K1,    K2,    K3,    N,t);  //  T(dr)
                                e.gd( uep ,xLocal,numberOfDimensions,isRectangular,0,0,0,0,K1+id1,K2+id2,K3+id3,N,t);  //  T(2*dr)
                                f(J1,J2,J3,tc) -= a[0]*ue + b1*( uep - uem )*factor;
                            }
                        }
            // OV_ABORT("DomainSolver::interfaceRightHandSide: get RHS for CHAMP -- FINISH ME");  
                }

            }
            else 
            {
        //  f = a[0]*u + a1*n.grad( u )
                    f(J1,J2,J3,tc) = a[0]*uLocal(I1,I2,I3,tc);
                    if( a[1]!=0. )
                    {
            // add on a[1]*( nu*u.n ) on the boundary 
            // **be careful** -- the normal changes sign on the two sides of the interface ---
                        MappedGridOperators & op = *(u.getOperators());
                        realSerialArray ux(I1,I2,I3,N), uy(I1,I2,I3,N);
                        op.derivative(MappedGridOperators::xDerivative,uLocal,ux,I1,I2,I3,N);
                        op.derivative(MappedGridOperators::yDerivative,uLocal,uy,I1,I2,I3,N);
                        if( cg.numberOfDimensions()==2 )
                        {
                            f(J1,J2,J3,tc) += a[1]*( normal(I1,I2,I3,0)*ux + normal(I1,I2,I3,1)*uy );
                        }
                        else
                        {
                            realSerialArray uz(I1,I2,I3);
                            op.derivative(MappedGridOperators::zDerivative,uLocal,uz,I1,I2,I3,N);
                            f(J1,J2,J3,tc) += a[1]*( normal(I1,I2,I3,0)*ux + normal(I1,I2,I3,1)*uy + normal(I1,I2,I3,2)*uz );
                        }
                    }
                    if( debug() & 8 )
                    {
                        printP("interfaceRightHandSide: (grid,side,axis) = (%d,%d,%d) \n"
                          "   I1=[%d,%d] I2=[%d,%d] I3=[%d,%d] : source (get) or source(set) \n" 
                          "   J1=[%d,%d] J2=[%d,%d] J3=[%d,%d] : target (get) or source (set) \n",
                          grid,side,axis,
                          I1.getBase(),I1.getBound(), I2.getBase(),I2.getBound(), I3.getBase(),I3.getBound(),
                          J1.getBase(),J1.getBound(), J2.getBase(),J2.getBound(), J3.getBase(),J3.getBound()
                          );
                        ::display(f(J1,J2,J3,tc) ,sPrintF("interfaceRHS: get heat-flux data:  %f*u + %f*( u.n ) (before TZ corrections)",a[0],a[1]));
                    }
                    if( twilightZoneFlow ) // turn this off for testing the case where the same TZ holds across all domains
                    {
            // ---add forcing for twilight-zone flow---
            //   ue <- a[0]*ue + a1*n.grad( ue )
                            OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                            const bool isRectangular = false; // ** do this for now ** mg.isRectangular();
                            if( !isRectangular )
                                mg.update(MappedGrid::THEcenter);
                            realArray & x= mg.center();
                            #ifdef USE_PPP
                                realSerialArray xLocal; 
                                if( !isRectangular ) 
                                    getLocalArrayWithGhostBoundaries(x,xLocal);
                            #else
                                const realSerialArray & xLocal = x;
                            #endif
                            realSerialArray ue(I1,I2,I3,N);
                            if( a[0]!=0. )
                            {
                                e.gd( ue ,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,N,t);  // exact solution 
                                ue(I1,I2,I3,N) = a[0]*ue(I1,I2,I3,N);
                            }
                            else
                            {
                                ue(I1,I2,I3,N) =0.;
                            }
                            if( a[1]!=0. )
                            {
                                realSerialArray uex(I1,I2,I3,N), uey(I1,I2,I3,N);
                                e.gd( uex ,xLocal,numberOfDimensions,isRectangular,0,1,0,0,I1,I2,I3,N,t);
                                e.gd( uey ,xLocal,numberOfDimensions,isRectangular,0,0,1,0,I1,I2,I3,N,t);
                                if( numberOfDimensions==2 )
                                {
                                    ue(I1,I2,I3,N) += a[1]*( normal(I1,I2,I3,0)*uex + normal(I1,I2,I3,1)*uey );
                                }
                                else
                                {
                                    realSerialArray uez(I1,I2,I3,N);
                                    e.gd( uez ,xLocal,numberOfDimensions,isRectangular,0,0,0,1,I1,I2,I3,N,t);
                                    ue(I1,I2,I3,N) += a[1]*( normal(I1,I2,I3,0)*uex + normal(I1,I2,I3,1)*uey + normal(I1,I2,I3,2)*uez ); 
                                }
                            }
            //   subtract off TZ flow:
            //   f <- f - ( a[0]*ue + a[1]*( nu*ue.n ) )
                        if( false )
                        {
                            ::display(f(J1,J2,J3,tc) ," a[0]*u + a[1]*( k u.n )");
                            ::display(ue(I1,I2,I3,tc)," a[0]*ue + a[1]*( k ue.n )");
                        }
                        f(J1,J2,J3,tc) -= ue(I1,I2,I3,N);
                        if( false )
                            ::display(f(J1,J2,J3,tc) ," a[0]*u + a[1]*( k u.n ) - [a[0]*ue + a[1]*( k ue.n )]");
                    }
            }

        }
        else
        {
            printF("DomainSolver::interfaceRightHandSide:ERROR: unknown option=%i\n",option);
            OV_ABORT("error");
        }

    // *wdh* Dec 11, 2021
        if( saveTimeHistory )
        {
      // -- save a time history of the heat flux
            if( !gfd.dbase.has_key("heatFluxHistory") )
            {
                gfd.dbase.put<ArrayEvolution>("heatFluxHistory");
                printP("\n @@@@@@ --interfaceRightHandSide-- Create ArrayEvolution heatFluxTimeHistory at t=%9.3e (for time extrap of champ RHS data) @@@@@@@ \n",t);
            }
            ArrayEvolution & heatFluxHistory = gfd.dbase.get<ArrayEvolution>("heatFluxHistory");

            heatFluxHistory.add( t, f );
            printP("\n @@@@@@ --interfaceRightHandSide-- Save heat flux time history at t=%9.3e @@@@@@@ \n",t);
        }



        
    // This next block was split and moved up above *wdh* April 22, 2021
    // if( // false &&  // turn this off for testing the case where the same TZ holds across all domains
    //     twilightZoneFlow )
    // {
    //   // ---add forcing for twilight-zone flow---
    //   //   ue <- a[0]*ue + a1*n.grad( ue )
    //   getTwilightZoneHeatFluxMacro()


    //   if( option==getInterfaceRightHandSide )
    //   { // get 
    //     //   subtract off TZ flow:
    //     //   f <- f - ( a[0]*ue + a[1]*( nu*ue.n ) )
    //     if( false )
    //     {
    //       ::display(f(J1,J2,J3,tc) ," a[0]*u + a[1]*( k u.n )");
    //       ::display(ue(I1,I2,I3,tc)," a[0]*ue + a[1]*( k ue.n )");
    //     }
    //     f(J1,J2,J3,tc) -= ue(I1,I2,I3,N);
    //     if( false )
    //     {
    //       ::display(f(J1,J2,J3,tc) ," a[0]*u + a[1]*( k u.n ) - [a[0]*ue + a[1]*( k ue.n )]");
    //     }
    //   }
    //   else if( option==setInterfaceRightHandSide )
    //   { // set 
    //     //   add on TZ flow:
    //     //   bd <- bd + a[0]*ue + a[1]*( nu*ue.n )
    //     bd(I1,I2,I3,tc) += ue(I1,I2,I3,N);

    //     if( false )
    //     {
    //       bd(I1,I2,I3,tc) = ue(I1,I2,I3,N);
    //     }
                

    //   }
    //   else
    //   {
    //     OV_ABORT("error");
    //   }
        
    // } // end if TZ 
        
    }
    else if( interfaceType(side,axis,grid)==Parameters::tractionInterface )
    {
    // *******************************************
    // ********** Traction Interface *************
    // *******************************************

        if( debug() & 2 )
            printP("DomainSolver::iterativeInterfaceRHS:traction: %s RHS for (grid,side,axis)=(%i,%i,%i) "
                          " t=%9.3e gfIndex=%i (current=%i)\n",
                          (option==0 ? "get" : "set"),grid,side,axis,t,gfIndex,current);


        const int uc = parameters.dbase.get<int >("uc");
        Range V(uc,uc+numberOfDimensions-1);

        if( option==setInterfaceRightHandSide )
        {
      // ***************************************************************************
      // ********************** SET TRACTION DATA **********************************
      // ***************************************************************************

      // **** set the RHS *****
      //   (TZ is done below) <- todo 

            bd=0.;
            bd(I1,I2,I3,V)=f(J1,J2,J3,V);

            if( debug() & 8 )
            {
                bd(I1,I2,I3,V).display("setInterfaceRightHandSide: Here is the RHS");
            }
        
        }


        else if( option==getInterfaceRightHandSide )
        {
      // ***************************************************************************
      // ********************** GET TRACTION DATA **********************************
      // ***************************************************************************

            realMappedGridFunction & u = gf[gfIndex].u[grid];
#ifdef USE_PPP
            realSerialArray uLocal; getLocalArrayWithGhostBoundaries(u,uLocal);
#else
            realSerialArray & uLocal = gf[gfIndex].u[grid];
#endif


      // **** we need to call a virtual function to evaluate the traction or the position of the interface 


      // ************ do this for now as a test ***************
            if( getClassName()=="Cgins" || getClassName()=="Cgcns" )
            {
        // We could optimize this for rectangular grids 
                mg.update(MappedGrid::THEvertexBoundaryNormal);
#ifdef USE_PPP
                const realSerialArray & normal = mg.vertexBoundaryNormalArray(side,axis);
#else
                const realSerialArray & normal = mg.vertexBoundaryNormal(side,axis);
#endif

        // *new* way : Use this: 

                Index Ib1,Ib2,Ib3;
        // getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
                Ib1=I1, Ib2=I2, Ib3=I3;
                realSerialArray traction(Ib1,Ib2,Ib3,numberOfDimensions);  // this should be a serial array

                gf[gfIndex].conservativeToPrimitive();

                int ipar[] = {grid,side,axis,gf[gfIndex].form}; // 
                real rpar[] = { gf[gfIndex].t }; // 

                parameters.getNormalForce( gf[gfIndex].u,traction,ipar,rpar );

                Range D=numberOfDimensions;

                #ifndef USE_PPP
                    f(J1,J2,J3,V)=traction(I1,I2,I3,D);
                #else
                    OV_ABORT("ERROR: finish me for parallel");
                #endif

        // f(J1,J2,J3,V)=-traction(I1,I2,I3,D); // &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



                if( debug() & 8 )
                {
                    f(J1,J2,J3,V).display("getInterfaceRightHandSide: Here is the RHS (traction=normalForce)");
                }

                if( interfaceDataOptions & Parameters::tractionRateInterfaceData )
                {
                
                    OV_ABORT("DomainSolver::getInterfaceRightHandSide: ERROR: unable to compute the time derivative of the traction!");
                }

//      const int pc = parameters.dbase.get<int >("pc");
//      realSerialArray p(I1,I2,I3);

//         p=1.;  // do this for now ******************************

        // ******** into what components should we fill the force? +++++++++++++++++++++++++++++++++++++=

//      f(J1,J2,J3,uc  )=normal(I1,I2,I3,0)*p(I1,I2,I3);
//      f(J1,J2,J3,uc+1)=normal(I1,I2,I3,1)*p(I1,I2,I3);
//         if( numberOfDimensions>2 )
//        f(J1,J2,J3,uc+2)=normal(I1,I2,I3,2)*p(I1,I2,I3);
                
//      f(J1,J2,J3,V)=0.;
//      f(J1,J2,J3,uc)=1.;

//      Range R4=f.dimension(3);
//      f(J1,J2,J3,R4)=0.;
//      f(J1,J2,J3,0)=1.;
                
            }
            else if( getClassName()=="Cgsm" )
            {
        // We could optimize this for rectangular grids 
                mg.update(MappedGrid::THEvertex);
#ifdef USE_PPP
                realSerialArray vertex; getLocalArrayWithGhostBoundaries(mg.vertex(),vertex);
#else
                const realSerialArray & vertex = mg.vertex();
#endif
        // return the position of the boundary 
                Range D(0,numberOfDimensions-1);

        // ******** into what components should we fill the positions? +++++++++++++++++++++++++++++++++++++=
                bool methodComputesDisplacements=true;
                if( parameters.dbase.has_key("methodComputesDisplacements") )
                {
                    methodComputesDisplacements=parameters.dbase.get<bool>("methodComputesDisplacements");
                }
        // Some solid-mechanics solvers compute the displacements, others the full deformation. 
                if( methodComputesDisplacements )
                    f(J1,J2,J3,V)=vertex(I1,I2,I3,D)+uLocal(I1,I2,I3,V);
                else
                {  // for now this is for Jeff's hemp code -- fix me --
                      if( true )
                      {
             // for linear elasticity, the displacement is saved in components u1c,u2c,u3c
                          if( !parameters.dbase.has_key("u1c") )
                          {
                              printP("interfaceRightHandSide:ERROR: unable to find variable u1c in the data-base\n");
                              OV_ABORT("error");
                          }
                          const int u1c = parameters.dbase.get<int >("u1c");
                          assert( u1c>=0 );
                          Range V1(u1c,u1c+numberOfDimensions-1);
                          f(J1,J2,J3,V)=vertex(I1,I2,I3,D)+uLocal(I1,I2,I3,V1);
                      }
                      else
                      {
                          f(J1,J2,J3,V)=uLocal(I1,I2,I3,V);
                      }
                      
                }
                
                if( debug() & 8 )
                {
                    f(J1,J2,J3,V).display("getInterfaceRightHandSide: Here is the RHS (vertex+displacement)");
                }
            }
            else
            {
                printP("ERROR: unknown className=[%s]\n",(const char*)getClassName());
                OV_ABORT("error");
            }
            
                
//       f(J1,J2,J3,V)=0.;
//       f(J1,J2,J3,uc)=1.;
            
        }
        



    }
    else
    {
        printP("DomainSolver::interfaceRightHandSide:unexpected interfaceType=%i for (grid,side,axis)=(%d,%d,%d)\n",
                      interfaceType(side,axis,grid),grid,side,axis);
        OV_ABORT("error");
    }
    

    return 0;
}
