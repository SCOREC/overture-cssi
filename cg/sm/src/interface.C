// This file automatically generated from interface.bC with bpp.
#include "Cgsm.h"
#include "SmParameters.h"
#include "ParallelUtility.h"
#include "Interface.h"  
#include "ArrayEvolution.h"    


// include some interface bpp macros
//------------------------------------------------------------------------------------
// This file contains macros used by the interface routines.
//------------------------------------------------------------------------------------



// ===========================================================================
// Get/set the interface RHS for a heat flux interface
// ===========================================================================

// ===================================================================================
/// \brief Return the interface data required for a given type of interface.
/// \param info (input) : the descriptor for the interface.
/// \param interfaceDataOptions (output) : a list of items from Parameters::InterfaceDataEnum that define
///                    which data to get (or which data were set).  Multiple items are
///                     chosen by bit-wise or of the different options
/// \return : the number of interface data items required. (use this to dimension arrays).
/// \note: this function should be over-loaded.
// ===================================================================================
int Cgsm::
getInterfaceDataOptions( GridFaceDescriptor & info, int & interfaceDataOptions ) const
{
    int numberOfItems=0;
    const int numberOfDimensions = parameters.dbase.get<int>("numberOfDimensions");

    SmParameters::PDEVariation & pdeVariation = parameters.dbase.get<SmParameters::PDEVariation>("pdeVariation");
    
    const int grid=info.grid, side=info.side, axis=info.axis;

    IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");
    if( grid<0 || grid>interfaceType.getBound(2) ||
            side<0 || side>1 || axis<0 || axis>interfaceType.getBound(1) )
    {
        printP("Cgsm::getInterfaceDataOptions:ERROR: invalid values: (grid,side,axis)=(%i,%i,%i)\n",
         	   grid,side,axis);
        OV_ABORT("Cgsm::getInterfaceDataOptions:ERROR");
    }
    if( interfaceType(side,axis,grid)==Parameters::heatFluxInterface )
    {
        interfaceDataOptions=Parameters::heatFluxInterfaceData;
        numberOfItems=1;
        
        printP("Cgsm::getInterfaceDataOptions:ERROR: not expecting a heatFluxInterface!\n");
        OV_ABORT("Cgsm::getInterfaceDataOptions:ERROR");
    }
    else if( interfaceType(side,axis,grid)==Parameters::tractionInterface ) 
    {
        const bool projectInterface = parameters.dbase.get<bool>("projectInterface");
        const bool useAddedMassAlgorithm = parameters.dbase.get<bool>("useAddedMassAlgorithm");
        
        if( true )
        {
            printP("--SM--getInterfaceDataOptions: useAddedMassAlgorithm = %i ***\n",(int)useAddedMassAlgorithm);

        }
        
        if( debug() & 2 )
            printP("*** Cgsm:getInterfaceDataOptions: projectInterface = %i ***\n",projectInterface);

        if( !projectInterface && !useAddedMassAlgorithm )
        {
      // -- the standard interface approximation requires the traction from the fluid

      // We need the traction (boundary force) at a tractionInterface:
            interfaceDataOptions=Parameters::tractionInterfaceData;
            numberOfItems+=numberOfDimensions;
        }
        else if( useAddedMassAlgorithm )
        {
      // *wdh* May 26 2017
      // We want the interface traction and velocity for the INS-SM amp-scheme
            interfaceDataOptions = ( Parameters::velocityInterfaceData     |
                         			       Parameters::tractionInterfaceData    
                                                          );

            numberOfItems+=2*numberOfDimensions;    
        }
        else
        {
      // *** HAS THIS OPTION EVER BEEN USED?? *wdh* May 26, 2017 ??

      // if we project the interface values then we need the following from the fluid:

            interfaceDataOptions = ( Parameters::positionInterfaceData     |
                         			       Parameters::velocityInterfaceData     |
                         			       Parameters::tractionInterfaceData    
                                                          );

            numberOfItems+=3*numberOfDimensions;
            
      // Parameters::accelerationInterfaceData 
      // Parameters::tractionRateInterfaceData 

        }
        
        if( pdeVariation==SmParameters::godunov )
        {
      // The godunov solver also needs the time derivative of the traction: 
            interfaceDataOptions = interfaceDataOptions | Parameters::tractionRateInterfaceData;
            numberOfItems+=numberOfDimensions;
        }


    }
    else
    {
        printP("Cgsm::getInterfaceDataOptions:ERROR: interfaceType(grid=%i,side=%i,axis=%i)=%i\n",
         	   grid,side,axis,interfaceType(side,axis,grid));
        OV_ABORT("Cgsm::getInterfaceDataOptions:ERROR");
    }

    return numberOfItems;
}

// =======================================================================================================
// Macro: Adjust the traction for TZ
//
//   param: OPTION (input): OPTION=[GET|SET]
//
// =======================================================================================================




// =======================================================================================================
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
Cgsm::
interfaceRightHandSide( InterfaceOptionsEnum option, 
                                                int interfaceDataOptions,
                                                GridFaceDescriptor & info, 
                                                GridFaceDescriptor & gfd, 
                  			int gfIndex, real t,
                                                bool saveTimeHistory /* = false */ )
{

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
            printF("Cgsm: interfaceRightHandSide:WARNING cannot find gfIndex to match t=%9.3e, using current with t=%.3e"
                          " (currentGF=%i nextGF=%i)\n",t,gf[current].t,currentGF,nextGF);
        }
        else
        {
        
            if( !(currentGF>=0 && nextGF>=0) )
            {
      	printF("Cgsm: interfaceRightHandSide:ERROR: t=%9.2e, current=%i gf[current].t=%9.2e, currentGF=%i, nextGF=%i\n",
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
      	printF("Cgsm: interfaceRightHandSide:WARNING cannot find gfIndex to match t=%9.3e\n"
             	       "      currentGF=%i, gf[currentGF].t=%9.3e, nextGF=%i, gf[nextGF].t=%9.3e\n",
             	       t,currentGF,gf[currentGF].t,nextGF,gf[nextGF].t);
      	if( fabs(gf[currentGF].t-t) <  fabs(gf[nextGF].t-t) )
        	  gfIndex=currentGF; 
      	else
        	  gfIndex=nextGF; 
	// OV_ABORT("fix me");
            }
        }
        
    }
    

  // *wdh* 081212 CompositeGrid & cg = gf[0].cg;
    CompositeGrid & cg = gf[gfIndex].cg;
    const int numberOfDimensions = cg.numberOfDimensions();
    const real epsT = REAL_EPSILON*100.*(1+fabs(t));
    const int & correctionStage = parameters.dbase.get<int>("correctionStage");
    
    const IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");

    const int grid=info.grid, side=info.side, axis=info.axis;

    if( grid<0 || grid>=cg.numberOfComponentGrids() ||
            side<0 || side>1 || axis<0 || axis>=cg.numberOfDimensions() )
    {
        printP("Cgsm::interfaceRightHandSide:ERROR: invalid values: (grid,side,axis)=(%i,%i,%i)\n",
         	   grid,side,axis);
        OV_ABORT("Cgsm::interfaceRightHandSide:ERROR");
    }

    MappedGrid & mg = cg[grid];
    RealArray & bd = parameters.getBoundaryData(side,axis,grid,mg);

    const int numberOfComponents = parameters.dbase.get<int >("numberOfComponents")- parameters.dbase.get<int >("numberOfExtraVariables");
    SmParameters::PDEVariation & pdeVariation = parameters.dbase.get<SmParameters::PDEVariation>("pdeVariation");
    const bool & twilightZoneFlow = parameters.dbase.get<bool >("twilightZoneFlow");

    assert( info.u != NULL );
    RealArray & f = *info.u;
    Index I1=f.dimension(0),I2=f.dimension(1),I3=f.dimension(2);
    Range Rx=numberOfDimensions;

  // For TZ we keep traction of where quanities are saved in the output array "f"
    int baseVelocity=-1, baseAcceleration=-1, baseTraction=-1, baseTractionRate=-1;
    
    FILE *&debugFile = parameters.dbase.get<FILE* >("debugFile");
    FILE *&pDebugFile = parameters.dbase.get<FILE* >("pDebugFile");

    if( interfaceType(side,axis,grid)==Parameters::heatFluxInterface )
    {
    // ****************************************
    // ********* Heat flux interface **********
    // ****************************************

        real *a = info.a;
        if( debug() & 4 )
        {
            printP("Cgsm::interfaceRHS:heatFlux %s RHS for (side,axis,grid)=(%i,%i,%i) a=[%5.2f,%5.2f]"
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
      // **** set the RHS *****
      //   (TZ is done below)
            bd(I1,I2,I3,tc)=f(I1,I2,I3,tc);
            if( false )
            {
                ::display(bd(I1,I2,I3,tc)," RHS values","%4.2f ");
            }
        }
        else if( option==getInterfaceRightHandSide )
        {
      // **** get the RHS ****
            realMappedGridFunction & u = gf[gfIndex].u[grid];
        #ifdef USE_PPP
            realSerialArray uLocal; getLocalArrayWithGhostBoundaries(u,uLocal);
        #else
            realSerialArray & uLocal = gf[gfIndex].u[grid];
        #endif
            f(I1,I2,I3,tc) = a[0]*uLocal(I1,I2,I3,tc);
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
                    f(I1,I2,I3,tc) += a[1]*( normal(I1,I2,I3,0)*ux + normal(I1,I2,I3,1)*uy );
                }
                else
                {
                    realSerialArray uz(I1,I2,I3);
                    op.derivative(MappedGridOperators::zDerivative,uLocal,uz,I1,I2,I3,N);
                    f(I1,I2,I3,tc) += a[1]*( normal(I1,I2,I3,0)*ux + normal(I1,I2,I3,1)*uy + normal(I1,I2,I3,2)*uz );
                }
            }
            if( debug() & 4 )
            {
                ::display(f(I1,I2,I3,tc) ,sPrintF("getRHS:  %f*u + %f*( u.n ) ",a[0],a[1]));
            }
        }
        else
        {
            printF("Cgsm::interfaceRightHandSide:ERROR: unknown option=%i\n",option);
            Overture::abort("error");
        }
        if( // false &&  // turn this off for testing the case where the same TZ holds across all domains
            parameters.dbase.get<bool >("twilightZoneFlow") )
        {
      // ---add forcing for twlight-zone flow---
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
            if( option==getInterfaceRightHandSide )
            { // get 
        //   subtract off TZ flow:
        //   f <- f - ( a[0]*ue + a[1]*( nu*ue.n ) )
                if( false )
                {
                    ::display(f(I1,I2,I3,tc) ," a[0]*u + a[1]*( k u.n )");
                    ::display(ue(I1,I2,I3,tc)," a[0]*ue + a[1]*( k ue.n )");
                }
                f(I1,I2,I3,tc) -= ue(I1,I2,I3,N);
                if( false )
                {
                    ::display(f(I1,I2,I3,tc) ," a[0]*u + a[1]*( k u.n ) - [a[0]*ue + a[1]*( k ue.n )]");
                }
            }
            else if( option==setInterfaceRightHandSide )
            { // set 
        //   add on TZ flow:
        //   bd <- bd + a[0]*ue + a[1]*( nu*ue.n )
                bd(I1,I2,I3,tc) += ue(I1,I2,I3,N);
                if( false )
                {
                    bd(I1,I2,I3,tc) = ue(I1,I2,I3,N);
                }
            }
            else
            {
                Overture::abort("error");
            }
        } // end if TZ 
        
    }
    else if( interfaceType(side,axis,grid)==Parameters::tractionInterface )
    {
    // *******************************************
    // ********** Traction Interface *************
    // *******************************************

        if( debug() & 2 )
            printP("Cgsm::interfaceRightHandSide:traction: %s RHS for (grid,side,axis)=(%i,%i,%i) "
           	     " t=%9.3e gfIndex=%i (current=%i)\n",
           	     (option==0 ? "get" : "set"),grid,side,axis,t,gfIndex,current);


        const int uc = parameters.dbase.get<int >("uc");
        const int v1c = parameters.dbase.get<int >("v1c");

        const bool methodComputesVelocity = v1c>=0;  // not all methods compute velocity directly

        Range Dc(uc,uc+numberOfDimensions-1);  // displacement components
        Range C, Cd;
        if( option==setInterfaceRightHandSide )
        {
      // -----------------------------------------------------
      // ------------ Interface data has been provided -------
      // ------------ Fill in the bd array             -------
      // -----------------------------------------------------
      // **** set the RHS *****

            int numSaved=0; // keeps track of how many things we have saved in f 

            bd=0.;

            if( interfaceDataOptions & Parameters::positionInterfaceData )
            {
        // --- interface position is given ---
                if( debug() & 4 )
        	  printP("Cgsm:interfaceRightHandSide: interface position provided t=%9.3e\n",t);

      	C=Range(numSaved,numSaved+numberOfDimensions-1);

        // the interface position is currently not used.
	// bd(I1,I2,I3,Dc)=f(I1,I2,I3,C);        // where should we put this?        

                numSaved+=numberOfDimensions;
            }
            
            if( interfaceDataOptions & Parameters::velocityInterfaceData )
            {
        // --- interface velocity is given ---
                if( debug() & 4 )
        	  printP("Cgsm:interfaceRightHandSide: interface velocity provided t=%9.3e\n",t);

      	C=Range(numSaved,numSaved+numberOfDimensions-1);

                if( true || debug() & 4 )
      	{
        	  printP(">>> --SM-- interfaceRHS: interface velocity provided t=%9.3e  in C=[%i,%i]<<<\n",
                                  t,C.getBase(),C.getBound());
      	}
	// -- for now save the velocity data in the dbase *wdh* May 26, 2017
                aString velocityDataName;
      	sPrintF(velocityDataName,"velocityG%iS%iA%i",grid,side,axis);
      	if( !parameters.dbase.has_key(velocityDataName) )
      	{
        	  InterfaceData & interfaceData = parameters.dbase.put<InterfaceData>(velocityDataName);
        	  interfaceData.u.redim(bd.dimension(0),bd.dimension(1),bd.dimension(2),numberOfDimensions);
        	  interfaceData.u=0;
      	}
      	InterfaceData & interfaceData = parameters.dbase.get<InterfaceData>(velocityDataName);
      	interfaceData.t=t;
      	interfaceData.u(I1,I2,I3,Rx)=f(I1,I2,I3,C);  


        // the interface velocity is currently not used.
	// bd(I1,I2,I3,Dc)=f(I1,I2,I3,C);     // where should we put this?          

                numSaved+=numberOfDimensions;
            }

            if( interfaceDataOptions & Parameters::accelerationInterfaceData )
            {
        // --- interface acceleration is given ---
                if( debug() & 4 )
        	  printP("Cgsm:interfaceRightHandSide: interface acceleration provided t=%9.3e\n",t);

      	C=Range(numSaved,numSaved+numberOfDimensions-1);

        // the interface acceleration is currently not used.
	// bd(I1,I2,I3,Dc)=f(I1,I2,I3,C);    // where should we put this?

                numSaved+=numberOfDimensions;
            }

            if( interfaceDataOptions & Parameters::tractionInterfaceData )
            {
        // --- interface traction is given ---
                if( debug() & 4 )
        	  printP("Cgsm:interfaceRightHandSide: interface traction provided t=%9.3e\n",t);

      	C=Range(numSaved,numSaved+numberOfDimensions-1);

	// bd(I1,I2,I3,Dc)=f(I1,I2,I3,Dc);   // old
      	bd(I1,I2,I3,Dc)=f(I1,I2,I3,C);    // ***** save traction data here for now -- FIX ME 

	// ::display(bd(I1,I2,I3,Dc),"Cgsm::interface traction from fluid");
      	

                numSaved+=numberOfDimensions;
            }

            if( interfaceDataOptions & Parameters::tractionRateInterfaceData )
            {
	// save the time derivative of the traction:
                if( debug() & 4 )
        	  printP("Cgsm:interfaceRightHandSide: the tractionRate is provided at t=%9.3e\n",t);

                C=Range(numSaved,numSaved+numberOfDimensions-1);

        // traction rate is filled into the velocity components (v1c,...)
                Range Vc(v1c,v1c+numberOfDimensions-1);
        // bd(I1,I2,I3,Vc)=f(I1,I2,I3,Dc+numberOfDimensions);      // ************** fix me bases
                bd(I1,I2,I3,Vc)=f(I1,I2,I3,C);

      	if( debug() & 8 )
      	{
        	  ::display(bd(I1,I2,I3,Vc),sPrintF("Cgsm::interfaceRightHandSide: time derivative of the traction, t=%9.3e",t),"%8.2e ");
      	}

      	numSaved+=numberOfDimensions;
            }
            

            if( debug() & 8 )
            {
                bd(I1,I2,I3,Dc).display("Cgsm::interfaceRightHandSide:set: Here is the RHS (traction)");
            }
        
        }

        else if( option==getInterfaceRightHandSide )
        {

      // ****************************************************************************
      // *********************Return the interface data******************************
      // ****************************************************************************

      // For TESTING return exact values (if known)
            bool useExactInterfacePosition=false;
            bool useExactInterfaceVelocity=false; 
            bool useExactInterfaceAcceleration=false;
            bool useExactInterfaceTraction=false;
            
            bool evalExactTraction=true;  // evaluate (but do not use) exact traction (for comparing to the computed)
            
      // Should we save the predicted values in the time history
            bool savePredictedVelocityInHistory=true;  // *default* true
            bool savePredictedTractionInHistory=true;  // *default* true
            
      // always eval traction from the time-history:
            bool useTractionTimeHistory=false; // *default* false

            if( useExactInterfacePosition )
                printF("--SM--IRHS: USING EXACT INTERFACE POSITIONS t=%9.3e ***TEMP***\n",t);
            
            if( useExactInterfaceVelocity )
                printF("--SM--IRHS: USING EXACT INTERFACE VELOCITIES t=%9.3e ***TEMP***\n",t);
            
            if( useExactInterfaceAcceleration )
                printF("--SM--IRHS: USING EXACT INTERFACE ACCELERATION t=%9.3e ***TEMP***\n",t);
            if( useExactInterfaceTraction )
                printF("--SM--IRHS: USING EXACT INTERFACE TRACTION t=%9.3e ***TEMP***\n",t);
            
            
            realMappedGridFunction & u = gf[gfIndex].u[grid];
            OV_GET_SERIAL_ARRAY(real,u,uLocal);

      // We could optimize this for rectangular grids 
            mg.update(MappedGrid::THEvertex);
            OV_GET_SERIAL_ARRAY(real,mg.vertex(),vertex);

            int numSaved=0; // keeps track of how many things we have saved in f 

            if( interfaceDataOptions & Parameters::positionInterfaceData )
            {
	// -- return the position of the boundary --

      	C=Range(numSaved,numSaved+numberOfDimensions-1); // save displacement in these components of f
                Cd=C;  // Save me for TZ below

      	if( !parameters.dbase.has_key("u1c") )
      	{
        	  printP("interfaceRightHandSide:ERROR: unable to find displacement component u1c in the data-base\n");
        	  OV_ABORT("error");
      	}
      	const int u1c = parameters.dbase.get<int >("u1c");
      	assert( u1c>=0 );
      	Range Dc(u1c,u1c+numberOfDimensions-1);  // displacement components
                RealArray displacement(I1,I2,I3,Rx);
                
                bool initialStateAssigned=false;
                if( t<=epsT || useExactInterfacePosition )
                {
          // -- get initial or past-time state of the interface  --
                    if( !twilightZoneFlow )
                    {
                        initialStateAssigned=getKnownInterfaceState(Parameters::boundaryPosition,t, side,axis,grid,
                                                                                                                mg,I1,I2,I3,Rx,displacement );
            // if( true ||  t<=2.*dt & debug() & 2 )
            // {
            //   display(displacement,sPrintF("--SM-- displacement from known solution at t=%9.3e",t),"%.2e ");
            // }

                    }
                    else
                    {
            // --- set past time position from TZ solution ---- *wdh* April 17, 2018
                        initialStateAssigned=true;
                        
                        printP("interfaceRightHandSide:set past-time interface position from TZ solution, t=%9.3e\n",t);
                        
                        OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                        bool isRectangular=false;
                        if( !isRectangular )
                            mg.update(MappedGrid::THEcenter | MappedGrid::THEvertexBoundaryNormal );
                        OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);

                        Range Dc(uc,uc+numberOfDimensions-1); // displacements
                        realSerialArray ue(I1,I2,I3,Dc);
                        e.gd( ue,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,Dc,t);  // ue : exact solution        

                        displacement=xLocal(I1,I2,I3,Rx)+ue(I1,I2,I3,Dc); 
                    }
                    
                }
                if( !initialStateAssigned )
                {

          // Some solid-mechanics solvers compute the displacements, others the full deformation.	
                    bool methodComputesDisplacements=true;
                    if( parameters.dbase.has_key("methodComputesDisplacements") )
                    {
                        methodComputesDisplacements=parameters.dbase.get<bool>("methodComputesDisplacements");
                    }
                    if( methodComputesDisplacements )
                    {
                        
                        displacement=vertex(I1,I2,I3,Rx)+uLocal(I1,I2,I3,Dc);  

                        if( debug() & 4 && t <= 3.*dt )
                        {
                            printP("interfaceRightHandSide:get interface position and save in components [%i,%i] t=%9.3e\n",
                                          C.getBase(),C.getBound(),t);
                                  
                            ::display(uLocal(I1,I2,I3,Dc),sPrintF("Cgsm:interfaceRightHandSide: uLocal: displacement t=%9.3e (before any TZ fix)",t),debugFile);
                            if( twilightZoneFlow )
                            {
                                OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
                                bool isRectangular=false;
                                if( !isRectangular )
                                    mg.update(MappedGrid::THEcenter | MappedGrid::THEvertexBoundaryNormal );
                                OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);

                                Range Dc(uc,uc+numberOfDimensions-1); // displacements
                                realSerialArray ue(I1,I2,I3,Dc);
                                e.gd( ue,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,Dc,t);  // ue : exact solution        

                                ::display(ue(I1,I2,I3,Dc),sPrintF("Cgsm:interfaceRightHandSide: exactDisplacement, t=%9.3e",t),debugFile);
                                ::display(fabs(ue(I1,I2,I3,Dc)-uLocal(I1,I2,I3,Dc)),sPrintF("Cgsm:interfaceRightHandSide: error in displacement, t=%9.3e",t),debugFile);
                            }
                            
                        }

                    }
                    else
                    {  // for now this is for Jeff's hemp code -- fix me --
                        if( true )
                        {
              // for linear elasticity, the displacement is saved in components u1c,u2c,u3c
                            displacement=vertex(I1,I2,I3,Rx)+uLocal(I1,I2,I3,Dc);
                        }
                        else
                        {
                            displacement=uLocal(I1,I2,I3,Dc);
                        }
         	   
                    }
                }
                
      	f(I1,I2,I3,C)=displacement;
                if( t<=2.*dt & debug() & 4 )
                {
                    display(f(I1,I2,I3,C),sPrintF("--SM-- interface positions at t=%9.3e",t),debugFile,"%.2e ");
                }
                
                
                if( !methodComputesVelocity && saveTimeHistory )
                {
          // If the method does not compute the velocity then we save the time history of the displacement
          // so that we can compute the velocity and acceleration 
          // -- save a time history of the displacement
                    if( !gfd.dbase.has_key("displacementHistory") )
                    {
                        gfd.dbase.put<ArrayEvolution>("displacementHistory");

            // *wdh* May 1, 2018 -- increase order of accuracy to 3 (3-pt)
                        ArrayEvolution & displacementHistory = gfd.dbase.get<ArrayEvolution>("displacementHistory");
                        const int timeOrderOfAccuracy=parameters.dbase.get<int>("interfaceArrayEvolutionTimeAccuracy"); 
                        displacementHistory.setOrderOfAccuracy(timeOrderOfAccuracy);
                    }
                    ArrayEvolution & displacementHistory = gfd.dbase.get<ArrayEvolution>("displacementHistory");

                    printF("--SM----IRHS-- Save displacement time history at t=%9.3e\n",t);
                    displacementHistory.add( t, displacement);  

                }

      	if( debug() & 8 )
      	{
        	  f(I1,I2,I3,C).display("interfaceRightHandSide:get: Here is the RHS (vertex+displacement)");
      	}

      	numSaved+=numberOfDimensions;
            }
            
            
            if( interfaceDataOptions & Parameters::velocityInterfaceData )
            {
        // ------------------------------------
	// --- save the interface velocity ----
        // ------------------------------------

                bool useTimeHistoryForPredictorStage=false;

      	if( debug() & 2 )
        	  printP("Cgsm:interfaceRightHandSide: Save the interface velocity.\n");
      	
                C=Range(numSaved,numSaved+numberOfDimensions-1);   // save velocity in these components of f

                baseVelocity=C.getBase();  // save for TZ correction

      	const int v1c = parameters.dbase.get<int >("v1c");
                Range Rx=numberOfDimensions;
                RealArray velocity(I1,I2,I3,Rx);

                bool initialStateAssigned=false;
                if( t<=epsT || useExactInterfaceVelocity )
                {
          // -- get initial or past-time state of the interface  --
                    initialStateAssigned=getKnownInterfaceState(Parameters::boundaryVelocity,t, side,axis,grid,
                                                                                                            mg,I1,I2,I3,Rx,velocity );

                }
                if( !initialStateAssigned )
                {
                    if( v1c>=0 )
                    {
                        Range V(v1c,v1c+numberOfDimensions-1);           
                        if( !useTimeHistoryForPredictorStage || correctionStage!=0 )
                        {
                            velocity=uLocal(I1,I2,I3,V);
                        }
                        else
                        { // ** TEST** Use time extrap for predictor stage
                            assert( correctionStage==0 );
                            ArrayEvolution & velocityHistory = gfd.dbase.get<ArrayEvolution>("velocityHistory");
                            int numberOfTimeDerivatives=0;
                            const int timeOrderOfAccuracy=parameters.dbase.get<int>("interfaceArrayEvolutionTimeAccuracy"); 
                            int order=timeOrderOfAccuracy;
                            printF("--SM--interface: Use time-history for velocity on predictor stage, correctionStage=%i, order=%i,  t=%.3e\n",
                                          correctionStage,order,t);
                            velocityHistory.eval( t, velocity, numberOfTimeDerivatives, order );  
                        }
                        
                    }
                    else
                    {
            // Scheme does not compute velocities --- approximate by differences:
            // const int u1c = parameters.dbase.get<int >("u1c");
            // Range U(u1c,u1c+numberOfDimensions-1);         

                        if( gfd.dbase.has_key("displacementHistory") )
                        {
                            ArrayEvolution & displacementHistory = gfd.dbase.get<ArrayEvolution>("displacementHistory");

              // Do we need this to higher order so we can get the acceleration to 2nd order?
              // const int timeOrderOfAccuracy=parameters.dbase.get<int>("interfaceArrayEvolutionTimeAccuracy"); 
              // int orderOfAccuracy=3;
                            printF("--SM----IRHS-- EVAL the velocity as d(u)/dt from displacement time history, t=%9.3e\n",t);
                            const int numberOfTimeDerivatives=1;
                            displacementHistory.eval( t,velocity,numberOfTimeDerivatives ); // should be 3rd order accuracy now 
                        }
                        else
                        {
                            OV_ABORT("ERROR: displacementHistory not found. finish me");
                        }

                    }
                }
                
        // **** Comment: If we want to under-relax the solid velocity it would probably be done here *wdh* March 27, 2019
                const bool & relaxCorrectionSteps = parameters.dbase.get<int>("relaxCorrectionSteps");
                if( relaxCorrectionSteps && t>0. )
                {
                    const real & omega = parameters.dbase.get<real>("addedMassRelaxationFactor");
        	  const real & tol = parameters.dbase.get<real>("subIterationConvergenceTolerance");
                    const int correctionStage = parameters.dbase.get<int>("correctionStage");

                    ArrayEvolution & velocityHistory = gfd.dbase.get<ArrayEvolution>("velocityHistory");
                    RealArray velocityCurrent(I1,I2,I3,Rx);
                    velocityHistory.eval(t,velocityCurrent); // eval current guess for the velocity
                    real maxDiff = max(fabs(velocity-velocityCurrent));
                    bool isConverged = maxDiff < tol;
                    parameters.dbase.get<int>("correctionIterationsHaveConverged")=isConverged;  // Will be used by advanceSteps
                    
                    printF("--SM----IRHS-- relaxCorrectionSteps: t=%9.3e, omega=%g, tol=%g, correction=%i, maxDiff=%8.2e, isConverged=%i\n",
                                  t,omega,tol,correctionStage,maxDiff,(int)isConverged);

                    
                    
                    f(I1,I2,I3,C)=  omega*velocity + (1.-omega)*velocityCurrent;
                }
                else
                {
                      f(I1,I2,I3,C)=velocity;
                }
                
      	numSaved+=numberOfDimensions;


        // Save a time history of the velocity on the interface that can be used to
        // predict the acceleration
                if( saveTimeHistory )
                {
          // -- save a time history of the velocity
                    if( !gfd.dbase.has_key("velocityHistory") )
                    {
                        gfd.dbase.put<ArrayEvolution>("velocityHistory");
                        ArrayEvolution & velocityHistory = gfd.dbase.get<ArrayEvolution>("velocityHistory");
            // *wdh* May 1, 2018 -- increase order of accuracy to 3 (3-pt)
                        const int timeOrderOfAccuracy=parameters.dbase.get<int>("interfaceArrayEvolutionTimeAccuracy"); 

                        int order=timeOrderOfAccuracy;
            // int order=2;  // *wdh* May 2, 2018
                        velocityHistory.setOrderOfAccuracy(order);
                    }
                    ArrayEvolution & velocityHistory = gfd.dbase.get<ArrayEvolution>("velocityHistory");

          // *wdh* May 1, 2018-- do not save velocity from the predictor step -- this is bad for some reason
                    if( t<=0. || correctionStage!=0 || savePredictedVelocityInHistory )
                    {
                        if( debug() & 4 )
                            printF("--SM----IRHS-- Save interface velocity time history at t=%9.3e in ArrayEvolution - correctionStage=%i\n",
                                          t,correctionStage);
                        velocityHistory.add( t, velocity);  
                            
                        if( max(fabs(velocity))>1e5 ) /// **TEMP**
                        {
                            ::display(velocity,"--SM-- velocity","%.2e ");
                            ::display(uLocal,"--SM-- uLocal","%.2e ");
                            
                            OV_ABORT("ERROR");
                        }

                    }
                    
                }

                if( FALSE ) // **************** OLD *************
                {
                    const int & numberOfInterfaceTimeLevels=parameters.dbase.get<int>("numberOfInterfaceTimeLevels");
                    aString velocityTimeHistoryName;
                    sPrintF(velocityTimeHistoryName,"velocityTimeHistoryG%iS%iA%i",grid,side,axis);
                    if( !parameters.dbase.has_key(velocityTimeHistoryName) )
                    {
                        std::vector<InterfaceData> & velocityTimeHistory =  
                            parameters.dbase.put<std::vector<InterfaceData> >(velocityTimeHistoryName);
                        velocityTimeHistory.resize(numberOfInterfaceTimeLevels); // allocate space 
                        for( int i=0; i<numberOfInterfaceTimeLevels; i++ )
                        {
                            velocityTimeHistory[i].t=-REAL_MAX*.1;  // set negative bogus value for time 
                        }
                  
                    }
                    std::vector<InterfaceData> & velocityTimeHistory =
                        parameters.dbase.get<std::vector<InterfaceData> >(velocityTimeHistoryName);
                
                    int & currentTimeLevel = parameters.dbase.get<int>("currentInterfaceTimeLevel");
                    if( currentTimeLevel==-1 )
                    {
                        currentTimeLevel=0;
                        velocityTimeHistory[currentTimeLevel].t =t;
                    }
                    else if( t > velocityTimeHistory[currentTimeLevel].t )
                    {
            // add a new time level in the history
                        currentTimeLevel = ovmod(currentTimeLevel + 1, numberOfInterfaceTimeLevels);
                    }
              
                    assert( currentTimeLevel>=0 && currentTimeLevel<velocityTimeHistory.size() );
                    velocityTimeHistory[currentTimeLevel].t=t;
          // ** June 28, 2017 velocityTimeHistory[currentTimeLevel].u=uLocal(I1,I2,I3,V);
                    velocityTimeHistory[currentTimeLevel].u=velocity;
                    printF("--SM-- *OLD* SET INTERFACE VELOCITY TIME HISTORY current=%i, t=%9.3e\n",currentTimeLevel,t);
                }
                
            }
        
            if( interfaceDataOptions & Parameters::accelerationInterfaceData )
            {
	// -- save the interface acceleration --
      	C=Range(numSaved,numSaved+numberOfDimensions-1);   // save acceleration in these components of f
                baseAcceleration=C.getBase();  // save for TZ correction

      	if( debug() & 2 )
                {
        	  printP("Cgsm:interfaceRightHandSide: return the interface acceleration, t=%9.3e.\n",t);
                    fprintf(debugFile,"Cgsm:interfaceRightHandSide: return the interface acceleration t=%9.3e.\n",t);
                }
                
                bool initialStateAssigned=false;
                if( t<=epsT || useExactInterfaceAcceleration )
                {
          // -- get initial or past-time state of the interface  --
                    initialStateAssigned=getKnownInterfaceState(Parameters::boundaryAcceleration,t, side,axis,grid,
                                                                                                            mg,I1,I2,I3,C,f );

                }
                if( !initialStateAssigned )
                {
                    getInterfaceAcceleration( gfd, t, side,axis,grid, mg, I1,I2,I3,C, f,useExactInterfaceAcceleration );
                }
                

	// f(I1,I2,I3,C)=uLocal(I1,I2,I3,V);

      	numSaved+=numberOfDimensions;
            }
        
            if( interfaceDataOptions & Parameters::tractionInterfaceData )
            {
	// -- save the interface traction --
      	C=Range(numSaved,numSaved+numberOfDimensions-1);   // save traction in these components of f
                baseTraction=C.getBase();  // save for TZ correction

      	if( debug() & 2 )
        	  printP(">>>Cgsm:interfaceRHS: Save interface traction.. t=%9.3e in C=[%i,%i] <<<\n",
                                    t,C.getBase(),C.getBound());

                bool initialStateAssigned=false;
                if( t<=epsT || useExactInterfaceTraction )
                {
          // -- get initial or past-time state of the interface  --
                    initialStateAssigned=getKnownInterfaceState(Parameters::boundaryTraction,t, side,axis,grid,
                                                                                                            mg,I1,I2,I3,C,f );

                }
                if( !initialStateAssigned )
                {

                    const int pdeTypeForGodunovMethod = parameters.dbase.get<int >("pdeTypeForGodunovMethod");
                    if(  pdeTypeForGodunovMethod!=0 ) // nonlinear solid
                    {
                        printF("--SM-- get interface traction: WARNING: Not implemented for a nonlinear solid\n");
                    }
                    else if(  pdeTypeForGodunovMethod==0 ) // linear-elasticity
                    {
            // ------ get the traction  --------
                        const int s11c = parameters.dbase.get<int >("s11c"); assert( s11c>=0 );
                        const int s12c = parameters.dbase.get<int >("s12c"); assert( s12c>=0 );
                        const int s22c = parameters.dbase.get<int >("s22c"); assert( s22c>=0 );

            // -- here is the normal to the un-deformed surface -- do this for now
                        if( t <= 2.*dt )
                            printF("\n --SM-- IRHS WARNING Using REFERENCE NORMAL FOR TRACTION ***FIX ME*** ***TEMP*** \n\n");
                        
                        mg.update(MappedGrid::THEvertexBoundaryNormal);
                        OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);
                    
            // 
            // traction:  nv^T sigma
            // We need the normal and Cauchy stress tensor
                        if( numberOfDimensions==2 )
                        {
          	    
                            int c0=numSaved, c1=c0+1;
                            f(I1,I2,I3,c0)=( normal(I1,I2,I3,0)*uLocal(I1,I2,I3,s11c) +
                                                              normal(I1,I2,I3,1)*uLocal(I1,I2,I3,s12c) );
      	
                            f(I1,I2,I3,c1)=( normal(I1,I2,I3,0)*uLocal(I1,I2,I3,s12c) +
                                                              normal(I1,I2,I3,1)*uLocal(I1,I2,I3,s22c) );
                        }
                        else
                        {
                            OV_ABORT("finish me for 3D");
                        }
                    }
                    
                    
                    RealArray traction2(I1,I2,I3,Rx);
                    if( useTractionTimeHistory && correctionStage==0 )
                    {
                        if( gfd.dbase.has_key("tractionHistory") )
                        {
                            ArrayEvolution & tractionHistory = gfd.dbase.get<ArrayEvolution>("tractionHistory");

              // Do we need this to higher order so we can get the acceleration to 2nd order?
              // const int timeOrderOfAccuracy=parameters.dbase.get<int>("interfaceArrayEvolutionTimeAccuracy"); 
              // int orderOfAccuracy=3;
                            printF("--SM----IRHS-- EVAL the traction from the traction time history, t=%9.3e\n",t);
                            tractionHistory.eval( t,traction2 ); 

                        }
                        else
                        {
                            OV_ABORT("ERROR: tractionHistory not found. finish me");
                        }
                    }
                    
                    if( t<=3.*dt && debug() & 4 )
                    {
                        Range Ct=Range(baseTraction,baseTraction+numberOfDimensions-1);
            // ::display(f(I1,I2,I3,Ct),sPrintF("--SM--IRHS: computed traction t=%9.3e",t), debugFile,"%6.2f ");

                        if( evalExactTraction )
                        {
                            RealArray traction(I1,I2,I3,Rx);
                            getKnownInterfaceState(Parameters::boundaryTraction,t, side,axis,grid,
                                                                          mg,I1,I2,I3,Rx,traction );
                            real diff=max(fabs(f(I1,I2,I3,C)-traction));
                            printF("-SM--> Max error in traction = %9.3e, |traction|=%8.2e (correctionStage=%i)\n",
                                          diff,max(fabs(traction)),correctionStage);

                            if( useTractionTimeHistory && correctionStage==0 )
                            {
                                real diff=max(fabs(traction2-traction));
                                printF("-SM--> Max error in traction (time history) = %9.3e, |traction|=%8.2e (correctionStage=%i)\n",
                                              diff,max(fabs(traction)),correctionStage);
                            }
              // if( false )
              // {
              //   printF("--SM-- interface ****TEMP*** use exact traction\n");
              //   f(I1,I2,I3,C)=traction; // ************* TEMP *************
              // }
                            
                        }
                        

                    }
                    if( useTractionTimeHistory && correctionStage==0 )
                    {
                        printF("--SM-- interface: use traction from time-history for predictor stage \n");
                        f(I1,I2,I3,C)=traction2;
                    }

                    
                } // end if !initialStateAssigned
                

        // Save a time history of the traction on the interface *wdh* May 1, 2018
                if( !saveTimeHistory )
                {
                    OV_ABORT("ERROR : traction: saveTimeHistory=false");
                }
                
                if( saveTimeHistory )
                {
          // -- save a time history of the traction
                    if( !gfd.dbase.has_key("tractionHistory") )
                    {
                        gfd.dbase.put<ArrayEvolution>("tractionHistory");
                        ArrayEvolution & tractionHistory = gfd.dbase.get<ArrayEvolution>("tractionHistory");
            // *wdh* May 1, 2018 -- increase order of accuracy to 3 (3-pt)
                        const int timeOrderOfAccuracy=parameters.dbase.get<int>("interfaceArrayEvolutionTimeAccuracy"); 
                        int order=timeOrderOfAccuracy;
            // int order=2;  // *wdh* May 2, 2018
                        tractionHistory.setOrderOfAccuracy(order);
                    }
                    ArrayEvolution & tractionHistory = gfd.dbase.get<ArrayEvolution>("tractionHistory");

          // *wdh* May 1, 2018-- do not save traction from the predictor step -- this is bad for some reason
                    if( t<=0. || correctionStage!=0 || savePredictedTractionInHistory )
                    {
                        if( debug() & 4 )
                            printF("--SM----IRHS-- Save interface traction time history at t=%9.3e in ArrayEvolution - correctionStage=%i\n",
                                          t,correctionStage);
                        Range Rx=numberOfDimensions;
                        RealArray traction(I1,I2,I3,Rx);
                        traction=f(I1,I2,I3,C);
                        tractionHistory.add( t, traction );
                    }
                } // end if saveTimeHistory 
                
      	

//- // compute the solid normal (n1s,n2s)
//- #beginMacro getSolidNormal()
//-  rx2=rsxy2(j1,j2,j3,axis2,0)
//-  ry2=rsxy2(j1,j2,j3,axis2,1)
//-  r2Norm=normalSign2*max(epsx,sqrt(rx2**2+ry2**2))
//-  n1s=rx2/r2Norm
//-  n2s=ry2/r2Norm
//- #endMacro
//- 
//-         // -- linear elasticity: 
//- 	s11s =u2(j1,j2,j3,s11c);
//- 	s12s =u2(j1,j2,j3,s12c);
//- 	s22s =u2(j1,j2,j3,s22c);
//-         // compute the solid normal (n1s,n2s)
//- 	getSolidNormal();
//-         //  solid traction is ns.sigmas:
//-         traction1 = n1s*s11s + n2s*s12s  
//-         traction2 = n1s*s12s + n2s*s22s


        // OV_ABORT("Finish me: save the traction");
      	
	// f(I1,I2,I3,C)=uLocal(I1,I2,I3,V);

      	numSaved+=numberOfDimensions;
            }
        
            if( interfaceDataOptions & Parameters::tractionRateInterfaceData )
            {
	// -- save the interface traction rate --
      	if( debug() & 2 )
        	  printP("Cgsm:interfaceRightHandSide: Save the interface traction rate. FINISH ME \n");

      	C=Range(numSaved,numSaved+numberOfDimensions-1);   // save traction rate in these components of f
                baseTractionRate=C.getBase();  // save for TZ correction

                OV_ABORT("ERROR");
                
	// f(I1,I2,I3,C)=uLocal(I1,I2,I3,V);

      	numSaved+=numberOfDimensions;
            }
        


        } // end option==getInterfaceRightHandSide
        


    // *****************************************************************************
    // ************************ Twilight Zone Forcing ******************************
    // *****************************************************************************
    // if( t< 2.*dt && debug() & 4 )
    // {
    //   // ::display(ue(I1,I2,I3,Dc),"Exact displacements\n");
    //   Range Ct=Range(baseTraction,baseTraction+numberOfDimensions-1);
    //   ::display(f(I1,I2,I3,Ct),sPrintF("interface: traction *I* (t=%9.3e)\n",t),debugFile);
    // }

        if( parameters.dbase.get<bool >("twilightZoneFlow") )
        {
      // ---add forcing for twilight-zone flow---

            OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));

            const bool isRectangular = false; // ** do this for now ** mg.isRectangular();

            if( !isRectangular )
      	mg.update(MappedGrid::THEcenter | MappedGrid::THEvertexBoundaryNormal );

            realArray & x= mg.center();
#ifdef USE_PPP
            realSerialArray xLocal; 
            if( !isRectangular ) 
      	getLocalArrayWithGhostBoundaries(x,xLocal);
#else
            const realSerialArray & xLocal = x;
#endif

            const int rc=parameters.dbase.get<int >("rc");
            const int uc=parameters.dbase.get<int >("uc");
            const int vc=parameters.dbase.get<int >("vc");
            const int wc=parameters.dbase.get<int >("wc");


            const int s11c=parameters.dbase.get<int >("s11c");
            const int s12c=parameters.dbase.get<int >("s12c");
            const int s13c=parameters.dbase.get<int >("s13c");
            const int s21c=parameters.dbase.get<int >("s21c");
            const int s22c=parameters.dbase.get<int >("s22c");
            const int s23c=parameters.dbase.get<int >("s23c");
            const int s31c=parameters.dbase.get<int >("s31c");
            const int s32c=parameters.dbase.get<int >("s32c");
            const int s33c=parameters.dbase.get<int >("s33c");

            if( option==setInterfaceRightHandSide )
            { // set 
	//   add on TZ flow:
	//   bd <- bd + (true boundary traction)

                
        // -- adjust the traction and/or tractionRate by adding on the TZ traction
                    #ifdef USE_PPP
                    const realSerialArray & normal = mg.vertexBoundaryNormalArray(side,axis);
                    #else
                    const realSerialArray & normal = mg.vertexBoundaryNormal(side,axis);
                    #endif
                    Range Dc(uc,uc+numberOfDimensions-1); // displacements
                    RealArray & muGrid = parameters.dbase.get<RealArray>("muGrid");
                    RealArray & lambdaGrid = parameters.dbase.get<RealArray>("lambdaGrid");
                    real lambda = lambdaGrid(grid);
                    real mu = muGrid(grid);
                    real alpha = lambda+2.*mu;
                    if( pdeVariation==SmParameters::nonConservative || pdeVariation==SmParameters::conservative )
                    {
            // TZ forcing for traction for SOS codes: compute traction from displacements 
                        realSerialArray uex(I1,I2,I3,Dc),uey(I1,I2,I3,Dc);
                        e.gd( uex,xLocal,numberOfDimensions,isRectangular,0,1,0,0,I1,I2,I3,Dc,t);  // u.x : exact solution 
                        e.gd( uey,xLocal,numberOfDimensions,isRectangular,0,0,1,0,I1,I2,I3,Dc,t);  // u.y : exact solution 
                        if( interfaceDataOptions & Parameters::tractionInterfaceData  )
                        {
                                bd(I1,I2,I3,0) += (normal(I1,I2,I3,0)*( alpha*uex(I1,I2,I3,uc)+lambda*uey(I1,I2,I3,vc)) +
                                                                      normal(I1,I2,I3,1)*( mu*( uey(I1,I2,I3,uc)+uex(I1,I2,I3,vc) ) ) );
                                bd(I1,I2,I3,1) += (normal(I1,I2,I3,0)*( mu*( uey(I1,I2,I3,uc)+uex(I1,I2,I3,vc) )) +
                                                                      normal(I1,I2,I3,1)*( alpha*uey(I1,I2,I3,vc)+lambda*uex(I1,I2,I3,uc) ) );
                        }
                        if( interfaceDataOptions & Parameters::tractionRateInterfaceData )
                        {
              // tractionRate: time derivative of the traction
              // compute time derivatives of velocity spatial derivatives
               // traction rate is filled into the velocity components (v1c,...)
                            e.gd( uex,xLocal,numberOfDimensions,isRectangular,1,1,0,0,I1,I2,I3,Dc,t);  // v.x : exact solution 
                            e.gd( uey,xLocal,numberOfDimensions,isRectangular,1,0,1,0,I1,I2,I3,Dc,t);  // v.y : exact solution 
              // *note* uex is really vex
                                bd(I1,I2,I3,v1c+0) += (normal(I1,I2,I3,0)*( alpha*uex(I1,I2,I3,uc)+lambda*uey(I1,I2,I3,vc)) +
                                                       			       normal(I1,I2,I3,1)*( mu*( uey(I1,I2,I3,uc)+uex(I1,I2,I3,vc) ) ) );
                                bd(I1,I2,I3,v1c+1) += (normal(I1,I2,I3,0)*( mu*( uey(I1,I2,I3,uc)+uex(I1,I2,I3,vc) )) +
                                             			       normal(I1,I2,I3,1)*( alpha*uey(I1,I2,I3,vc)+lambda*uex(I1,I2,I3,uc) ) );
                        }
                    }
                    else if( pdeVariation==SmParameters::godunov )
                    {
            // TZ forcing for traction for FOS codes: compute traction from stress	  
                        Range Sc(s11c,s22c);
                        realSerialArray ue(I1,I2,I3,Sc);
                        e.gd( ue,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,Sc,t);   // exact values for the stress
                        if( interfaceDataOptions & Parameters::tractionInterfaceData  )
                        {
              // if( t<3.*dt && debug() & 4 )
              // {
              //   Range Ct=Range(baseTraction,baseTraction+1);
              //   ::display(f(I1,I2,I3,Ct),sPrintF("--SM--IRHS: before adjust traction for TZ t=%9.3e",t),
              //             debugFile,"%6.2f ");
              // }
                                bd(I1,I2,I3,uc) += normal(I1,I2,I3,0)*ue(I1,I2,I3,s11c) + normal(I1,I2,I3,1)*ue(I1,I2,I3,s21c);
                                bd(I1,I2,I3,vc) += normal(I1,I2,I3,0)*ue(I1,I2,I3,s12c) + normal(I1,I2,I3,1)*ue(I1,I2,I3,s22c);
                        }
            // if( t<3.*dt && debug() & 4 )
            // {
            //   RealArray te(I1,I2,I3,2);
            //   te(I1,I2,I3,0) = normal(I1,I2,I3,0)*ue(I1,I2,I3,s11c) + normal(I1,I2,I3,1)*ue(I1,I2,I3,s21c);
            //   te(I1,I2,I3,1) = normal(I1,I2,I3,0)*ue(I1,I2,I3,s12c) + normal(I1,I2,I3,1)*ue(I1,I2,I3,s22c);
            //   ::display(te,sPrintF("--SM--IRHS: exact traction t=%9.3e",t),
            //                   debugFile,"%6.2f ");
            //   Range Ct=Range(baseTraction,baseTraction+1);
            //   ::display(f(I1,I2,I3,Ct),sPrintF("--SM--IRHS: after adjust traction for TZ t=%9.3e",t),
            //                   debugFile,"%6.2f ");
            // }
            // ::display(ue(I1,I2,I3,s11c),"TZ ue(I1,I2,I3,s11c)");
            // ::display(bd(I1,I2,I3,uc),"TZ traction on solid boundary");
                        if( interfaceDataOptions & Parameters::tractionRateInterfaceData )
                        {
              // tractionRate: time derivative of the traction
              // compute time derivatives of velocity spatial derivatives
               // traction rate is filled into the velocity components (v1c,...)
                            e.gd( ue,xLocal,numberOfDimensions,isRectangular,1,0,0,0,I1,I2,I3,Sc,t);   // exact values for the d(stress)/dt 
                                bd(I1,I2,I3,v1c+0) += normal(I1,I2,I3,0)*ue(I1,I2,I3,s11c) + normal(I1,I2,I3,1)*ue(I1,I2,I3,s21c);
                                bd(I1,I2,I3,v1c+1) += normal(I1,I2,I3,0)*ue(I1,I2,I3,s12c) + normal(I1,I2,I3,1)*ue(I1,I2,I3,s22c);
                        }
                    }
                    else
                    {
                        printF("interface:ERROR: pdeVariation=%i not implemented yet. \n",(int)pdeVariation);
                    }
                    if( numberOfDimensions==3 )
                    {
                        printP("interface:ERROR: adjust traction forcing in 3D for TZ. Finish me!\n");
                        OV_ABORT("error");
                    }

      	
        // bd(I1,I2,I3,V) += xe(I1,I2,I3,V);
      	
            }
            else if( option==getInterfaceRightHandSide )
            { // get 
	//   subtract off TZ flow:
        //   f <- f - ( sm-boundary-position )  +  ( TZ-interface-position)

                if( interfaceDataOptions & Parameters::positionInterfaceData ) // *wdh* April 20, 2018
                {
                    Range D=numberOfDimensions;
                    Range Dc(uc,uc+numberOfDimensions-1); // displacements
                    realSerialArray ue(I1,I2,I3,Dc);
                    e.gd( ue,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,Dc,t);  // u : exact solution        


                    bool methodComputesDisplacements=true;
                    if( parameters.dbase.has_key("methodComputesDisplacements") )
                        methodComputesDisplacements=parameters.dbase.get<bool>("methodComputesDisplacements");

          // Some solid-mechanics solvers compute the displacements, others the full deformation.	
                    if( methodComputesDisplacements )
                    {
                        f(I1,I2,I3,Cd) -= xLocal(I1,I2,I3,D) + ue(I1,I2,I3,Dc); 

                        const real *v0 = parameters.dbase.get<real [3]>("tzInterfaceVelocity");  // this should match vg0 in the Deforming body motion 
                        const real *a0 = parameters.dbase.get<real [3]>("tzInterfaceAcceleration");         // this should match ag0 in the Deforming body motion 
                        if( t<3.*dt )
                        {
                            printP("interface:INFO: set TZ interface at t=%9.3e for (grid,side,axis)=(%i,%i,%i) to be x0 + v0*t + a0*t^2\n"
                                          " v0=(%g,%g,%g), a0=(%g,%g,%g) (NOTE: these should match the values for the DeformingBodyMotion)\n",t,grid,side,axis,
                                          v0[0],v0[1],v0[2], a0[0],a0[1],a0[2]);
                            if( debug() & 4 )
                            {
                // ::display(ue(I1,I2,I3,Dc),"Exact displacements\n");
                                ::display(f(I1,I2,I3,Cd),sPrintF(" Displacements before adding x0 + v0*t + a0*t^2 (t=%9.3e)",t),
                                                  debugFile);
                            }
                        
                        }
        	  
                        for( int axis=0; axis<numberOfDimensions; axis++ )
                        {
                            f(I1,I2,I3,Cd.getBase()+axis) += xLocal(I1,I2,I3,axis) + t*( v0[axis] + t*( a0[axis] ) );
                        }

                        if( debug() & 4 &&  t<3.*dt )
                        {
                            ::display(f(I1,I2,I3,Cd)," positions xLocal+Displacements AFTER adding x0 + v0*t + a0*t^2\n",debugFile);
                        }

                    }
                    else
                    {
                        printP("interface:ERROR: subtract boundary position for TZ with !methodComputesDisplacements. Finish me!\n");
                        OV_ABORT("error");
                    }
                }
                
                if( interfaceDataOptions & Parameters::tractionInterfaceData ||
                        interfaceDataOptions & Parameters::tractionRateInterfaceData )
                {
          // ADJUST REQUESTED TRACTION FOR TZ 

                    if( t<3.*dt )
                        printP("interface:INFO: adjust requested traction for TZ, t=%8.2e.\n",t);

          // OV_ABORT("Cgsm:interface: interface traction requested -- fix for TZ");

          // -- adjust the traction by subtractingthe TZ traction
          // FIX ME -- NEED TO ADJUST CERTAIN COMPONENTS OF f, not bd
                        #ifdef USE_PPP
                        const realSerialArray & normal = mg.vertexBoundaryNormalArray(side,axis);
                        #else
                        const realSerialArray & normal = mg.vertexBoundaryNormal(side,axis);
                        #endif
                        Range Dc(uc,uc+numberOfDimensions-1); // displacements
                        RealArray & muGrid = parameters.dbase.get<RealArray>("muGrid");
                        RealArray & lambdaGrid = parameters.dbase.get<RealArray>("lambdaGrid");
                        real lambda = lambdaGrid(grid);
                        real mu = muGrid(grid);
                        real alpha = lambda+2.*mu;
                        if( pdeVariation==SmParameters::nonConservative || pdeVariation==SmParameters::conservative )
                        {
              // TZ forcing for traction for SOS codes: compute traction from displacements 
                            realSerialArray uex(I1,I2,I3,Dc),uey(I1,I2,I3,Dc);
                            e.gd( uex,xLocal,numberOfDimensions,isRectangular,0,1,0,0,I1,I2,I3,Dc,t);  // u.x : exact solution 
                            e.gd( uey,xLocal,numberOfDimensions,isRectangular,0,0,1,0,I1,I2,I3,Dc,t);  // u.y : exact solution 
                            if( interfaceDataOptions & Parameters::tractionInterfaceData  )
                            {
                                    assert( baseTraction>=0 );
                                    f(I1,I2,I3,baseTraction+0) -= (normal(I1,I2,I3,0)*( alpha*uex(I1,I2,I3,uc)+lambda*uey(I1,I2,I3,vc)) +
                                                                                                  normal(I1,I2,I3,1)*( mu*( uey(I1,I2,I3,uc)+uex(I1,I2,I3,vc) ) ) );
                                    f(I1,I2,I3,baseTraction+1) -= (normal(I1,I2,I3,0)*( mu*( uey(I1,I2,I3,uc)+uex(I1,I2,I3,vc) )) +
                                                                                                  normal(I1,I2,I3,1)*( alpha*uey(I1,I2,I3,vc)+lambda*uex(I1,I2,I3,uc) ) );
                            }
                            if( interfaceDataOptions & Parameters::tractionRateInterfaceData )
                            {
                // tractionRate: time derivative of the traction
                // compute time derivatives of velocity spatial derivatives
                 // traction rate is filled into the velocity components (v1c,...)
                                e.gd( uex,xLocal,numberOfDimensions,isRectangular,1,1,0,0,I1,I2,I3,Dc,t);  // v.x : exact solution 
                                e.gd( uey,xLocal,numberOfDimensions,isRectangular,1,0,1,0,I1,I2,I3,Dc,t);  // v.y : exact solution 
                // *note* uex is really vex
                                    assert( baseTractionRate>=0 );
                                    f(I1,I2,I3,baseTractionRate+0) -= (normal(I1,I2,I3,0)*( alpha*uex(I1,I2,I3,uc)+lambda*uey(I1,I2,I3,vc)) +
                                                                       			                   normal(I1,I2,I3,1)*( mu*( uey(I1,I2,I3,uc)+uex(I1,I2,I3,vc) ) ) );
                                    f(I1,I2,I3,baseTractionRate+1) -= (normal(I1,I2,I3,0)*( mu*( uey(I1,I2,I3,uc)+uex(I1,I2,I3,vc) )) +
                                                             			                   normal(I1,I2,I3,1)*( alpha*uey(I1,I2,I3,vc)+lambda*uex(I1,I2,I3,uc) ) );
                            }
                        }
                        else if( pdeVariation==SmParameters::godunov )
                        {
              // TZ forcing for traction for FOS codes: compute traction from stress	  
                            Range Sc(s11c,s22c);
                            realSerialArray ue(I1,I2,I3,Sc);
                            e.gd( ue,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,Sc,t);   // exact values for the stress
                            if( interfaceDataOptions & Parameters::tractionInterfaceData  )
                            {
                // if( t<3.*dt && debug() & 4 )
                // {
                //   Range Ct=Range(baseTraction,baseTraction+1);
                //   ::display(f(I1,I2,I3,Ct),sPrintF("--SM--IRHS: before adjust traction for TZ t=%9.3e",t),
                //             debugFile,"%6.2f ");
                // }
                                    assert( baseTraction>=0 );
                                    f(I1,I2,I3,baseTraction+0) -= normal(I1,I2,I3,0)*ue(I1,I2,I3,s11c) + normal(I1,I2,I3,1)*ue(I1,I2,I3,s21c);
                                    f(I1,I2,I3,baseTraction+1) -= normal(I1,I2,I3,0)*ue(I1,I2,I3,s12c) + normal(I1,I2,I3,1)*ue(I1,I2,I3,s22c);
                            }
              // if( t<3.*dt && debug() & 4 )
              // {
              //   RealArray te(I1,I2,I3,2);
              //   te(I1,I2,I3,0) = normal(I1,I2,I3,0)*ue(I1,I2,I3,s11c) + normal(I1,I2,I3,1)*ue(I1,I2,I3,s21c);
              //   te(I1,I2,I3,1) = normal(I1,I2,I3,0)*ue(I1,I2,I3,s12c) + normal(I1,I2,I3,1)*ue(I1,I2,I3,s22c);
              //   ::display(te,sPrintF("--SM--IRHS: exact traction t=%9.3e",t),
              //                   debugFile,"%6.2f ");
              //   Range Ct=Range(baseTraction,baseTraction+1);
              //   ::display(f(I1,I2,I3,Ct),sPrintF("--SM--IRHS: after adjust traction for TZ t=%9.3e",t),
              //                   debugFile,"%6.2f ");
              // }
              // ::display(ue(I1,I2,I3,s11c),"TZ ue(I1,I2,I3,s11c)");
              // ::display(bd(I1,I2,I3,uc),"TZ traction on solid boundary");
                            if( interfaceDataOptions & Parameters::tractionRateInterfaceData )
                            {
                // tractionRate: time derivative of the traction
                // compute time derivatives of velocity spatial derivatives
                 // traction rate is filled into the velocity components (v1c,...)
                                e.gd( ue,xLocal,numberOfDimensions,isRectangular,1,0,0,0,I1,I2,I3,Sc,t);   // exact values for the d(stress)/dt 
                                    assert( baseTractionRate>=0 );
                                    bd(I1,I2,I3,baseTractionRate+0) -= normal(I1,I2,I3,0)*ue(I1,I2,I3,s11c) + normal(I1,I2,I3,1)*ue(I1,I2,I3,s21c);
                                    bd(I1,I2,I3,baseTractionRate+1) -= normal(I1,I2,I3,0)*ue(I1,I2,I3,s12c) + normal(I1,I2,I3,1)*ue(I1,I2,I3,s22c);
                            }
                        }
                        else
                        {
                            printF("interface:ERROR: pdeVariation=%i not implemented yet. \n",(int)pdeVariation);
                        }
                        if( numberOfDimensions==3 )
                        {
                            printP("interface:ERROR: adjust traction forcing in 3D for TZ. Finish me!\n");
                            OV_ABORT("error");
                        }

                    if( t< 2.*dt && debug() & 4 )
                    {
            // ::display(ue(I1,I2,I3,Dc),"Exact displacements\n");
                        Range Ct=Range(baseTraction,baseTraction+numberOfDimensions-1);
                        ::display(f(I1,I2,I3,Ct),sPrintF("interface:get traction (after adjust for TZ) (t=%9.3e)",t),debugFile);
                    }

                }

                if( interfaceDataOptions & Parameters::velocityInterfaceData )
                {
                    if( v1c>=0 )
                    {
                        Range V(v1c,v1c+numberOfDimensions-1);           
                        RealArray vex(I1,I2,I3,V);
                        e.gd( vex,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,V,t); // exact velocity

            // OV_ABORT("Cgsm:interface: interface velocity requested -- fix for TZ");

                        assert( baseVelocity>=0 );
                        Range Rx=numberOfDimensions;
                        f(I1,I2,I3,Rx+baseVelocity) -= vex(I1,I2,I3,Rx+v1c);
                        
                    }
                    else
                    {
                          OV_ABORT("Cgsm:interface: interface velocity requested -- fix for TZ and SOS scheme");
                    }
                    
                }

                if( interfaceDataOptions & Parameters::accelerationInterfaceData )
                {
          // OV_ABORT("Cgsm:interface: interface acceleration requested -- fix for TZ");
                    if( v1c>=0 )
                    {
                        Range V(v1c,v1c+numberOfDimensions-1);           
                        RealArray aex(I1,I2,I3,V);
                        e.gd( aex,xLocal,numberOfDimensions,isRectangular,1,0,0,0,I1,I2,I3,V,t); // exact acceleration = v.t 

            // OV_ABORT("Cgsm:interface: interface velocity requested -- fix for TZ");

                        assert( baseAcceleration>=0 );
                        Range Rx=numberOfDimensions;
                        f(I1,I2,I3,Rx+baseAcceleration) -= aex(I1,I2,I3,Rx+v1c);

                        if( t< 2.*dt && debug() & 4 )
                        {
              // ::display(ue(I1,I2,I3,Dc),"Exact displacements\n");
                            Range Ct=Range(baseTraction,baseTraction+numberOfDimensions-1);
                            ::display(f(I1,I2,I3,Rx+baseAcceleration),
                                      sPrintF("interface:get acceleration (after adjust for TZ) (t=%9.3e)",t),debugFile);
                        }
                    }
                    else
                    {
                          OV_ABORT("Cgsm:interface: interface velocity requested -- fix for TZ and SOS scheme");
                    }
                }
                
                if( interfaceDataOptions & Parameters::tractionRateInterfaceData )
                {
                    OV_ABORT("Cgsm:interface: interface tractionRate requested -- fix for TZ");
                }
                

            }
            else
            {
      	OV_ABORT("error");
            }
        
        } // end if TZ 




    }
    else
    {
        printF("Cgsm::interfaceRightHandSide:unexpected interfaceType=%i\n",interfaceType(side,axis,grid));
        OV_ABORT("error");
    }
    

    return 0;
}




// ================================================================================================
/// \brief project interface values for the AMP scheme.
/// \notes: Currently this routine projects the velocity on the interface for the INS-SM AMP scheme.

//  -- This routine assigns the velocity on the interface using the value from the
//     the fluid which should be equal to the projected interface velocity, the interface
//     traction from the fluid will be assigned following tis routine in the assignBC routine.
// ================================================================================================
int Cgsm::
projectInterface( int grid, real dt, int current )
{
    const bool useAddedMassAlgorithm = parameters.dbase.get<bool>("useAddedMassAlgorithm");

  // if( FALSE && useAddedMassAlgorithm )
  // {
  //   printF("\n--SM-- projectInterfaceVelocity DO NOT PROJECT VELOCITY **TEMP***\n");
  //   return 0;
  // }
    
    if( !useAddedMassAlgorithm )
        return 0;
    
    FILE *&debugFile = parameters.dbase.get<FILE* >("debugFile");
    FILE *&pDebugFile = parameters.dbase.get<FILE* >("pDebugFile");

    const bool projectInterface = parameters.dbase.get<bool>("projectInterface");
    GridFunction & cgf = gf[current];
    const real t= cgf.t; 
    if( t<2.*dt )
    {
        printF("--SM-- projectInterface: t=%9.2e, grid=%i useAddedMassAlgorithm=%i projectInterface=%i\n",
         	   t,grid,useAddedMassAlgorithm,projectInterface);
        if( debug() & 4 )
        {
            fPrintF(debugFile,"\n--SM-- projectInterface: t=%9.2e, grid=%i useAddedMassAlgorithm=%i projectInterface=%i\n",
         	   t,grid,useAddedMassAlgorithm,projectInterface);
        }
        
    }
    
    realMappedGridFunction & u = gf[current].u[grid];
    OV_GET_SERIAL_ARRAY(real,u,uLocal);
    MappedGrid & mg = cgf.cg[grid];
    const int numberOfDimensions=cg.numberOfDimensions();
    const Parameters::InterfaceCommunicationModeEnum & interfaceCommunicationMode= 
        parameters.dbase.get<Parameters::InterfaceCommunicationModeEnum>("interfaceCommunicationMode");

    const Parameters::KnownSolutionsEnum & knownSolution = 
        parameters.dbase.get<Parameters::KnownSolutionsEnum >("knownSolution"); 

    IntegerArray & interfaceType = parameters.dbase.get<IntegerArray >("interfaceType");
    const int & correctionStage = parameters.dbase.get<int>("correctionStage");

    Index Ibv[3], &Ib1=Ibv[0], &Ib2=Ibv[1], &Ib3=Ibv[2];
    const int uc   = parameters.dbase.get<int >("uc");
    const int v1c  = parameters.dbase.get<int >("v1c");
    const int s11c = parameters.dbase.get<int >("s11c");
    const int s12c = parameters.dbase.get<int >("s12c");
    const int s21c = parameters.dbase.get<int >("s21c");
    const int s22c = parameters.dbase.get<int >("s22c");
    Range Dc(uc,uc+numberOfDimensions-1);    // displacement components
    Range Vc(v1c,v1c+numberOfDimensions-1);
    Range Rx = numberOfDimensions;
      	
    const real rho   =parameters.dbase.get<real>("rho");
    const real lambda=parameters.dbase.get<real>("lambda");
    const real mu    =parameters.dbase.get<real>("mu");
    const real cp    =sqrt((lambda+2.*mu)/rho);
    const real cs    =sqrt(mu/rho);
    const real zp    =rho*cp;
    const real zs    =rho*cs;

    InterfaceData interfaceData;

    for( int axis=axis1; axis<mg.numberOfDimensions(); axis++ )
    {
        for( int side=Start; side<=End; side++ )
        {
            if( interfaceType(side,axis,grid)==Parameters::tractionInterface )
            {
      	if( t <= 2.*dt && debug() & 4 )
                {
        	  fPrintF(debugFile,"--SM-- projectInterface for (side,axis,grid)=(%i,%i,%i) bc=%i\n",
             		 side,axis,grid,mg.boundaryCondition(side,axis));
                }
                
                int extra=1;
                getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3,extra);  
                bool ok=ParallelUtility::getLocalArrayBounds(u,uLocal,Ib1,Ib2,Ib3,1);

                if( interfaceCommunicationMode==Parameters::requestInterfaceDataWhenNeeded )
                {
          // *new* way June 30, 2017 -- explicitly request interface data from Cgins, ...
                    if( t <= 2.*dt )
                    {
                        printF("--SM-- projectInterface: REQUEST interface velocity from fluid at t=%9.3e\n",t);
                        fPrintF(debugFile,"--SM-- projectInterface: REQUEST interface velocity from fluid at t=%9.3e\n",t);
                    }
                    
                    if( ok )
                        interfaceData.u.redim(Ib1,Ib2,Ib3,Rx);
                    interfaceData.t=t;
                    interfaceData.u=0;

                    int interfaceDataOptions = Parameters::velocityInterfaceData;
                    bool saveTimeHistory=false; // what should this be ?
                    getInterfaceData( t, grid, side, axis, 
                                                        interfaceDataOptions,
                                                        interfaceData.u,
                                                        parameters,saveTimeHistory );

                    RealArray & interfaceVelocity =interfaceData.u;
                    mg.update(MappedGrid::THEvertexBoundaryNormal);
                    OV_GET_VERTEX_BOUNDARY_NORMAL(mg,side,axis,normal);         


                    if( ok )
                    {

            // Project the velocity and traction on the interface for the INS-SM AMP scheme.

                        bool useImpedanceProjection=true; // ***********************************************


                        real alphaV=1., alphaS=0.;  // set below for impedance averages
                        
                        RealArray & bd = parameters.getBoundaryData(side,axis,grid,mg);

                        if( !useImpedanceProjection || correctionStage <0 )
                        {
              // old 
                            uLocal(Ib1,Ib2,Ib3,Vc)=interfaceVelocity(Ib1,Ib2,Ib3,Rx);
                        }
                        else
                        {
              // -------------------------------------------------
              // -- Project the interface velocity and traction --
              // -------------------------------------------------

              // Retrieve the parameters from the fluid
              // FIX ME -- lookup first time and then save locally 
                            Parameters & fluidParams = getInterfaceParameters( grid,side,axis,parameters );
                            const real nuFluid=fluidParams.dbase.get<real>("nu");
                            const real muFluid=fluidParams.dbase.get<real>("mu");
                            const real fluidDensity=fluidParams.dbase.get<real>("fluidDensity");
                            assert( fabs(fluidDensity*nuFluid/muFluid-1.) < REAL_EPSILON*100. );
                        

                            Index Ib1,Ib2,Ib3, Ip1,Ip2,Ip3;
                            getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3); // boundary line 

              //  --------------- Compute the normal distance of the first grid line  ---------------------
              //   *** SAVE THIS VALUE TO AVOID RECOMPUTING ****
                            real dn;
                            parameters.getNormalGridSpacing( mg,side,axis, dn );
    

              // mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
              // OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
              // Index Ib1,Ib2,Ib3, Ip1,Ip2,Ip3;
              // getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3); // boundary line 
              // getGhostIndex(mg.gridIndexRange(),side,axis,Ip1,Ip2,Ip3,-1); // first line in 

              // int includeGhost=0;
              // bool ok;
              // ok=ParallelUtility::getLocalArrayBounds(mg.vertex(),xLocal,Ib1,Ib2,Ib3,includeGhost); // restrict bounds to this processor
              // ok=ParallelUtility::getLocalArrayBounds(mg.vertex(),xLocal,Ip1,Ip2,Ip3,includeGhost);

              // // Assume grid is nearly orthogonal -- we could check using the mask 
              // real ndist=1.;
              // if( ok )
              // {
              //   if( numberOfDimensions==2 )
              //     ndist = sqrt( min(SQR(xLocal(Ip1,Ip2,Ip3,0)-xLocal(Ib1,Ib2,Ib3,0))+
              //                      SQR(xLocal(Ip1,Ip2,Ip3,1)-xLocal(Ib1,Ib2,Ib3,1))) );
              //   else
              //     ndist = sqrt( min(SQR(xLocal(Ip1,Ip2,Ip3,0)-xLocal(Ib1,Ib2,Ib3,0))+
              //                      SQR(xLocal(Ip1,Ip2,Ip3,1)-xLocal(Ib1,Ib2,Ib3,1))+
              //                      SQR(xLocal(Ip1,Ip2,Ip3,2)-xLocal(Ib1,Ib2,Ib3,2))) );
              // }
              // const real dn =ParallelUtility::getMinValue(ndist);


              // traction is saved here: bd(I1,I2,I3,Dc)=ui(I1,I2,I3,Rx);
                            RealArray fluidTraction(Ib1,Ib2,Ib3,Dc);
                            fluidTraction=bd(Ib1,Ib2,Ib3,Dc); // fluid traction has aleady been requested 


              // For added-mass (FIBR)
              //   zf = zfMuByH*(mu/h) + zfRhoHByDt*(rho*h/dt)
                            const real & zfMuByH = parameters.dbase.get<real>("zfMuByH");  
                            const real & zfRhoHByDt = parameters.dbase.get<real>("zfRhoHByDt");
                            const real & zfMono = parameters.dbase.get<real>("zfMono");  
                            const real & cfl = parameters.dbase.get<real>("cfl");  

              // zf = rho*nu/h = mu/h 
              //   zp*dt/rho = cp*dt approx h (approx)
              // => zf approx  (mu)/(cp*dt) 

                            real zf1 = muFluid/dn;
                            real zf2 = fluidDensity*dn/dt; 

              // real zfNew = 2.*(fluidDensity*muFluid)/(zp*dt);

              // const real fluidAddedMassLengthScale =  fluidParams.dbase.get<real>("fluidAddedMassLengthScale");
              // const real zfFluid = fluidDensity*fluidAddedMassLengthScale/dt; 

                            real zfWeighted= zfMuByH*zf1 + zfRhoHByDt*zf2; 
                            real zf=zfWeighted;
                            
                            bool useMonolithicImpedance=zfMono>0. ; // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                            real zfm=1., k=1.;
                            if( true )
                            {
                // --- Impedance from the monolithic algorithm ---
                                real rhof = fluidDensity;
                                real muf=muFluid;

                                const int axisp1 = (axis + 1) % mg.numberOfDimensions();
                                const real ds = mg.gridSpacing(axisp1);  
                                
                // k = 1./(2.*ds);
                                k = 1./dn;

                                real beta=sqrt( rhof/(muf*dt) + k*k );
                                zfm = (rhof*rhof + ( rhof*beta*zp + 4.*rhof*muf*k*k )*dt - ( 4.*(beta-k)*muf*muf*k*k*k )*dt*dt )/( (rhof+ zp*dt*(beta-k) )*dt*k );
                                
                                printF("SSSSSSSS sm: ");
                                printF("old zfm=%.2e, ",zfm);

                // this is asymptotically equivalent to old zfm
                                zfm = 2.*muf*k+1./(k*dt);
                                printF("new zfm=%.2e\n",zfm);

                            }
                            if( useMonolithicImpedance )
                            {
                                zf=zfMono*zfm;
                            }
                            else
                            {
                                printF("***** WARNING: NOT USING MONOLITHIC ZF, (zfMono=%.2e) *****\n",zfMono);
                            }
                            
                            /* ---- OLD
                            if( correctionStage>0 )
                            {
                                zf=1e8;
                                if( t<= 2.*dt )
                                    printF("SMSMSMSMS  CGSM correctionStage=%i: Set zf=%.2e\n",correctionStage,zf);
                            }
                            else
                            {
                                if( t<= 2.*dt )
                                    printF("SMSMSMSMS  CGSM correctionStage=%i, Use zf=%.2e\n",correctionStage,zf);

                            }
                              ---- */

              // coeff of [traction] in AMP interface condition:
                            const real & ampSigmaJumpCoeff = parameters.dbase.get<real>("ampSigmaJumpCoeff"); 
              // coeff of [velocity] in AMP interface condition
                            const real & ampVelocityJumpCoeff = parameters.dbase.get<real>("ampVelocityJumpCoeff"); 

                            real Cmu   = useMonolithicImpedance? 0. :  zfMuByH;
                            real Cmono = useMonolithicImpedance? zfMono :  0.;
                            if( t<= 2.*dt )
                                printF("SSSSSSSS sm: project zf =(%.2f)*mu/h + (%.2f)*rho*h/dt =%.2e, mu/h=%.2e, rho*h/dt=%.2e zMonolithic=%.2e zf=%.2e "
                                              "zfMono=%.2g, Cmu=%.2g Cmono=%.2g sigmaJumpCoeff=%.1f velocityJumpCoeff=%.1f \n",
                                              zfMuByH,zfRhoHByDt,zfWeighted,zf1,zf2,zfm, zf,zfMono,Cmu,Cmono,ampSigmaJumpCoeff,ampVelocityJumpCoeff);

                            real zfHat=muFluid*fluidDensity/(zp*dt);
                            real viscousCFL=(muFluid/fluidDensity)*dt/(dn*dn);
                        
                            
              // printF("SSSSSSSS sm: project zfOld=muFluid/dn = %.2e  zfNew=2*(rhoF*muF)/(zp*dt) = %.2e, "
              //     " zfFluid=rho*L/dt= %.2e, use: zf=%.2e\n",
              //        zfOld,zfNew,zfFluid,zf);

              //       ", rhoF=%.2e muFluid=%.2e nuFluid=%.2e zp=%.2e dt=%.2e\n",zfOld,zfNew,fluidDensity,muFluid,nuFluid,zp,dt);

                            alphaV = zf/( zf + zp );
                            alphaS = 1.-alphaV; // (1./zf)/( 1./zf + 1./zp );  // = zp/(zf+zp) = 1-alphaV
              // this is a multiplier for the jump terms
              // during the correction, this will be zero
                            real jumpMult = 1.; 
                            if( correctionStage>0 )
                            {
                                alphaV=1.; alphaS=0.; jumpMult = 0.;
                                printF("SMSMSMSMS  CGSM correctionStage=%i: Set alphaV=%.2e, alphaS=%.2e, jumpMult=%.2e\n",correctionStage,alphaV,alphaS,jumpMult);
                            }
                            else
                            {
                                printF("SMSMSMSMS  CGSM correctionStage=%i, Use zf=%.2e\n",correctionStage,zf);
                            }

              // printF("\\delta=%.2e \\zp=%.2e \\mu=%.2e \\zf=\\f{\\mu}{h}=%.2e \\zfh=\\f{\\rho\\dt}{h}=%.2e "
              //        "\\zfHat=%.2e \\zfm=%.2e (k=%.2e) \\f{\\nu\\dt}{h^2}=%.2e (h=%.2e) \\zfMuByH=%.2e \\alphaV=%.2e\n",
              //        rho,zp,muFluid,zf1,zf2,zfHat,zfm,k,viscousCFL,dn,zfMuByH,alphaV);
                            const real & tFinal = parameters.dbase.get<real>("tFinal"); 
                            if( t<=2.*dt )
                                printF(" %.2g & %.2g & %.2g & %.2g & %.2g & %.2g & %.2g & %.2g & %.2g & %.2g & %.2g & %.2g & %.2g & %.2g \\\\ \n",
                                              rho,muFluid,dn, viscousCFL,cfl,zp,zf1,zf2,zfHat,zfm,Cmu,Cmono,alphaV,tFinal);
                            

                            
                            if( t<=2.*dt )
                            {
                                printF("--SM-- projectInterface: muFluid=%.2e, dn=%.3e, zf=%.2e, zp=%.2e, alphaV=zf/(zf+zp)=%.2e alphaS=%.2e"
                                              " correctionStage=%i\n",
                                              muFluid,dn,zf,zp,alphaV,alphaS,correctionStage);
                            }
                            if( debug() & 4 )
                            {
                                fPrintF(debugFile,"--SM-- projectInterface: muFluid=%.2e, dn=%.3e, zf=%.2e, zp=%.2e, alphaV=zf/(zf+zp)=%.2e alphaS=%.2e"
                                              " correctionStage=%i\n",
                                              muFluid,dn,zf,zp,alphaV,alphaS,correctionStage);
                                
                            }
                        
                            RealArray solidTraction(Ib1,Ib2,Ib3,Rx);
                            assert( numberOfDimensions==2 );
                            solidTraction(Ib1,Ib2,Ib3,0) = (uLocal(Ib1,Ib2,Ib3,s11c)*normal(Ib1,Ib2,Ib3,0) + 
                                                                                            uLocal(Ib1,Ib2,Ib3,s12c)*normal(Ib1,Ib2,Ib3,1));
                            solidTraction(Ib1,Ib2,Ib3,1) = (uLocal(Ib1,Ib2,Ib3,s21c)*normal(Ib1,Ib2,Ib3,0) + 
                                                                                            uLocal(Ib1,Ib2,Ib3,s22c)*normal(Ib1,Ib2,Ib3,1));
                        


                            if( false )
                            {
                                ::display(interfaceVelocity(Ib1,Ib2,Ib3,Rx)," Fluid velocity","%.3e ");
                                ::display(uLocal(Ib1,Ib2,Ib3,Vc)," Solid velocity","%.3e ");
                                ::display(solidTraction," Solid traction","%.3e ");
                                ::display(fluidTraction," Fluid traction","%.3e ");
                            }

              // // Do this for now: *fix me* for multiple interfaces
              // if( !parameters.dbase.has_key("solidVelocityPredicted") )
              // {
              //   parameters.dbase.put<RealArray>("solidVelocityPredicted");
              //   parameters.dbase.put<RealArray>("solidTractionPredicted");
              // }
              // RealArray & solidVelocityPredicted = parameters.dbase.get<RealArray>("solidVelocityPredicted");
              // RealArray & solidTractionPredicted = parameters.dbase.get<RealArray>("solidTractionPredicted");
              // solidVelocityPredicted=uLocal(Ib1,Ib2,Ib3,Vc);
              // solidTractionPredicted=solidTraction;


              // project velocity - impedance average: 
              // ** check sign of the jump term:
              // ** fix me for normal and tangential components -- use zs for tangential
                            uLocal(Ib1,Ib2,Ib3,Vc)=( alphaV*interfaceVelocity(Ib1,Ib2,Ib3,Rx) + (1.-alphaV)*uLocal(Ib1,Ib2,Ib3,Vc)
                                                                              + jumpMult*ampSigmaJumpCoeff*(1./(zf+zp))*( fluidTraction - solidTraction)
                                );
              // project traction - inverse impedance average:
                            bd(Ib1,Ib2,Ib3,Dc) = ( alphaS*fluidTraction + (1.-alphaS)*solidTraction
                                                                          + jumpMult*ampVelocityJumpCoeff*(1./(1./zf+1./zp))*( interfaceVelocity(Ib1,Ib2,Ib3,Rx) -uLocal(Ib1,Ib2,Ib3,Vc) )
                                );

                            const int & fluidSolidCornerFix = parameters.dbase.get<int>("fluidSolidCornerFix");
                            if( fluidSolidCornerFix != 0 )
                            {
                                RealArray cornerFix;
                                printF("FFFFFFFF CGSM: projIntr: cornerFixup for sigmasI: fluidSolidCornerFix=%i\n",fluidSolidCornerFix);
                        
                                parameters.getFluidSolidCornerFixFunction( mg,grid,side,axis,cornerFix,SmParameters::displacementBC );
                                if( false )
                                    ::display(cornerFix,"CGSM: corner fix function","%9.2e ");

                                if( true )
                                {
                                    for( int d=Dc.getBase(); d<=Dc.getBound(); d++ )
                                    {
                    // ramp displacement to zero near adjacent displacement corners
                                        uLocal(Ib1,Ib2,Ib3,d) *=cornerFix(Ib1,Ib2,Ib3); 
                                    }
                                    for( int d=Vc.getBase(); d<=Vc.getBound(); d++ )
                                    {
                    // ramp velocity to zero near adjacent displacement corners
                                        uLocal(Ib1,Ib2,Ib3,d) *=cornerFix(Ib1,Ib2,Ib3); 
                                    }
                                }
                                
                                for( int d=Dc.getBase(); d<=Dc.getBound(); d++ )
                                {
                  // ramp traction to zero near adjacent displacement corners
                                    bd(Ib1,Ib2,Ib3,d) *=cornerFix(Ib1,Ib2,Ib3); 
                                }

                            }

              // else
              // {
              //   // Save the solid velocity and traction for the CORRECTOR stage

              //   RealArray & solidVelocityPredicted = parameters.dbase.get<RealArray>("solidVelocityPredicted");
              //   RealArray & solidTractionPredicted = parameters.dbase.get<RealArray>("solidTractionPredicted");

              //   // project velocity - impedance average: 
              //   uLocal(Ib1,Ib2,Ib3,Vc)=( alphaV*interfaceVelocity(Ib1,Ib2,Ib3,Rx) + (1.-alphaV)*solidVelocityPredicted
              //                                       + (1./(zf+zp))*( fluidTraction - solidTractionPredicted )
              //     );
              //   // project traction - inverse impedance average:
              //   // ** DO NOT** use solidVelocityPredicted in jump ??
              //   bd(Ib1,Ib2,Ib3,Dc) = ( alphaS*fluidTraction + (1.-alphaS)*solidTractionPredicted
              //                          //       + (1./(1./zf+1./zp))*( interfaceVelocity(Ib1,Ib2,Ib3,Rx) - uLocal(Ib1,Ib2,Ib3,Vc) )
              //     );
              // }

              // -----------------------------------------------------
              // --- save a time history of the interface traction ---
              // -----------------------------------------------------
              //    *wdh* July 19, 2018 
                            int orderOfAccuracyForTractionRate=2; // This requires 3 time levels of the traction.
              // **** FIX ME FOR MULTIPLE INTERFACES ****
                            if( !parameters.dbase.has_key("interfaceTractionHistory") )
                            {
                                parameters.dbase.put<ArrayEvolution>("interfaceTractionHistory");
                                ArrayEvolution & interfaceTractionHistory = parameters.dbase.get<ArrayEvolution>("interfaceTractionHistory");
                // const int timeOrderOfAccuracy=parameters.dbase.get<int>("interfaceArrayEvolutionTimeAccuracy"); 
                                interfaceTractionHistory.setOrderOfAccuracy(orderOfAccuracyForTractionRate);
                            }
                            ArrayEvolution & interfaceTractionHistory = parameters.dbase.get<ArrayEvolution>("interfaceTractionHistory");

                            RealArray interfaceTraction;
                            interfaceTraction=bd(Ib1,Ib2,Ib3,Dc);
                            interfaceTractionHistory.add( t, interfaceTraction);    // save the interface traction

              // --- evaluate the interface traction-rate ----- 
                            RealArray tractionRate(Ib1,Ib2,Ib3,Rx);
                            if( interfaceTractionHistory.getNumberOfTimeLevels() > orderOfAccuracyForTractionRate ) 
                            {
                // There are enough time history values save to compute the traction rate
                                const int numberOfTimeDerivatives=1;
                                interfaceTractionHistory.eval( t,tractionRate,numberOfTimeDerivatives);

                                if( t<= 4.*dt || debug() & 4 )
                                {
                                    printF(" --SM-- TTTTTTTTTTT projectInterface: eval traction-rate from traction time history order=%i t=%.3e\n",
                                                    orderOfAccuracyForTractionRate,t);
                                    fPrintF(debugFile," projectInterface: eval traction-rate from traction time history order=%i t=%.3e\n",
                                                    orderOfAccuracyForTractionRate,t);
                                }
                                
                            }
                            else
                            {
                                if( knownSolution==Parameters::userDefinedKnownSolution )
                                {
                                    int body=0;
                                    parameters.getUserDefinedDeformingBodyKnownSolution( body,Parameters::boundaryTractionRate,
                                                                                                                                              t, grid, mg, Ib1,Ib2,Ib3,Rx,tractionRate );
                                    printF("--SM-- TTTTTTTTTTTTTTTTTTTT Use EXACT traction-rate at t=%.3e (not enough time-history)\n",t);
                                    fPrintF(debugFile,"--SM-- TTTTTTTTTTTTTTTTTTTT Use EXACT traction-rate at t=%.3e (not enough time-history)\n",t);
                                }

                            }
                            
                            bd(Ib1,Ib2,Ib3,Vc) = tractionRate;  // set projected traction-rate


                        } // end impedance projection 
                        
                        

                        if( parameters.dbase.get<bool >("twilightZoneFlow") )
                        {
              // ---add forcing for twilight-zone flow--- 

                            OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));

                            const bool isRectangular = false; // ** do this for now ** mg.isRectangular();
                            mg.update(MappedGrid::THEcenter | MappedGrid::THEvertex);
                            OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
                            realSerialArray ve(Ib1,Ib2,Ib3,Vc);

                    	      e.gd( ve,xLocal,numberOfDimensions,isRectangular,0,0,0,0,Ib1,Ib2,Ib3,Vc,t);  

                            if( !useImpedanceProjection )
                            {
                                uLocal(Ib1,Ib2,Ib3,Vc) += ve;   // fluid velocity is corrected by TZ to be zero so set solid velocity to exact 
                            }
                            else
                            {
                // TZ fix for impedance average 
                                RealArray solidTractione(Ib1,Ib2,Ib3,Rx);
                                RealArray s11e(Ib1,Ib2,Ib3),s12e(Ib1,Ib2,Ib3),s21e(Ib1,Ib2,Ib3),s22e(Ib1,Ib2,Ib3);
                                
                                assert( numberOfDimensions==2 );

                                e.gd( s11e,xLocal,numberOfDimensions,isRectangular,0,0,0,0,Ib1,Ib2,Ib3,s11c,t); 
                                e.gd( s12e,xLocal,numberOfDimensions,isRectangular,0,0,0,0,Ib1,Ib2,Ib3,s21c,t); 
                                e.gd( s21e,xLocal,numberOfDimensions,isRectangular,0,0,0,0,Ib1,Ib2,Ib3,s12c,t); 
                                e.gd( s22e,xLocal,numberOfDimensions,isRectangular,0,0,0,0,Ib1,Ib2,Ib3,s22c,t); 

                                solidTractione(Ib1,Ib2,Ib3,0) = (s11e*normal(Ib1,Ib2,Ib3,0) + 
                                                                                                  s12e*normal(Ib1,Ib2,Ib3,1));
                                solidTractione(Ib1,Ib2,Ib3,1) = (s21e*normal(Ib1,Ib2,Ib3,0) + 
                                                                                                  s22e*normal(Ib1,Ib2,Ib3,1));
                                
                                
                                uLocal(Ib1,Ib2,Ib3,Vc) += alphaV*ve;
                                
                                bd(Ib1,Ib2,Ib3,Dc) += alphaS*solidTractione;
                                
                            }
                            

                        }



                        if( t <= 2.*dt && debug() & 4 )
                        {
                            ::display(uLocal(Ib1,Ib2,Ib3,Vc),sPrintF("--SM--PIV-- Projected interface velocity t=%.3e",t),debugFile,"%6.3f ");
                            ::display(bd(Ib1,Ib2,Ib3,Dc)    ,sPrintF("--SM--PIV-- Projected interface traction t=%.3e",t),debugFile,"%6.3f ");
                        }
                        
                    }
                    
                }
                else  // *OLD WAY
                {
                    aString velocityDataName;
                    sPrintF(velocityDataName,"velocityG%iS%iA%i",grid,side,axis);
                    if( !parameters.dbase.has_key(velocityDataName) )
                    {
            // -- At t=0 the interfaceData may not have been set yet, in this case call getInterfaceData

                        printF("--SM-- projectInterface::WARNING: interface data: [%s] not found, t=%9.3e!\n",
                                      (const char*)velocityDataName,t);
                        OV_ABORT("finish me");


                    }

                    InterfaceData & interfaceData = parameters.dbase.get<InterfaceData>(velocityDataName);
                    aString buff;
                    if( t <= 2.*dt )
                        ::display(interfaceData.u,sPrintF(buff,"--SM-- projectInterface: interface velocity [%s] at ti=%9.3e (t=%9.3e)",
                                                                                            (const char*)velocityDataName,interfaceData.t,t),"%6.3f ");
                    RealArray & interfaceVelocity =interfaceData.u;

                    if( ok )
                    {
            // Project the velocity and traction on the interface for the INS-SM AMP scheme.
                        uLocal(Ib1,Ib2,Ib3,Vc)=interfaceVelocity(Ib1,Ib2,Ib3,Rx);


                    }
                }
                
            }
        }
    }
    
  // *wdh* May 8, 2018 -- do a periodic update -- what about non-periodic ?? **FIX ME**
  // really only need to do Vc
    u.periodicUpdate();
    

    return 0;
    
}
    
//============================================================================================================
/// \brief Return the acceleration on the interface
/// \param t (input) : time to eval the acceleration
/// \param side,axis,grid : face
/// \param mg (input) : MappedGrid
/// \param I1,I2,I3,C (input) : assign points f(I1,I2,I3,C) 
/// \param f (output) : put acceleration here 
//============================================================================================================
int Cgsm::
getInterfaceAcceleration( GridFaceDescriptor & gfd, const real t, const int side, const int axis, const int grid, 
                                                    MappedGrid & mg, const Index & I1, const Index & I2, const Index & I3, const Range & C, 
                                                    RealArray & f, bool useExactInterfaceAcceleration )
{
    const bool & twilightZoneFlow = parameters.dbase.get<bool >("twilightZoneFlow");

    const bool evalExactAcceleration=false; // eval exact accel. for testing 
  // const bool useExactAcceleration=false && !twilightZoneFlow;
    bool useTimeHistory=!useExactInterfaceAcceleration;
    const int & correctionStage = parameters.dbase.get<int>("correctionStage");
    
    RealArray accel;
    ArrayEvolution & velocityHistory = gfd.dbase.get<ArrayEvolution>("velocityHistory");
    
    if( evalExactAcceleration || (velocityHistory.getNumberOfTimeLevels()<= 1 && twilightZoneFlow) )
    {
    // --- evaluate the exact acceleration for testing ---    
        if( twilightZoneFlow )
        {
      // get exact acceleration from TZ
            OGFunction & e = *(parameters.dbase.get<OGFunction* >("exactSolution"));
            mg.update(MappedGrid::THEcenter | MappedGrid::THEvertex );
            OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
            const int numberOfDimensions=mg.numberOfDimensions();
            const int v1c = parameters.dbase.get<int >("v1c");
            Range Vc(v1c,v1c+numberOfDimensions-1);
            accel.redim(I1,I2,I3,Vc);
            const bool isRectangular = false; // ** do this for now ** mg.isRectangular();
            e.gd( accel,xLocal,numberOfDimensions,isRectangular,1,0,0,0,I1,I2,I3,Vc,t);  // (v)_t : exact solution 
            
        }
        else
        {
            Range Rx=mg.numberOfDimensions();
            accel.redim(I1,I2,I3,Rx);

      // get exact accel from a known solution 
            const Parameters::KnownSolutionsEnum & knownSolution = 
                parameters.dbase.get<Parameters::KnownSolutionsEnum >("knownSolution");  
      // printF("--SM-- getInterfaceAcceleration: knownSolution=%i\n",(int)knownSolution);

            int bodyNumber=0; // fix me for multiple deforming bodies
          	    
      // RealArray state(I1,I2,I3,Rx);   // Fix me -- could avoid a copy by passing C to the next function:
            parameters.getUserDefinedDeformingBodyKnownSolution( 
                bodyNumber,
                Parameters::boundaryAcceleration, 
                t, grid, mg,I1,I2,I3,Rx, accel );
        }

    }
    if( useExactInterfaceAcceleration )
    {
        f(I1,I2,I3,C)=accel;
        
        if( true || debug() & 2 )
        {
            printF("--SM-- getInterfaceAcceleration: Setting acceleration from EXACT known solution ***TEMP**"
                          " at t=%9.3e.\n",t);
            ::display(f(I1,I2,I3,C),sPrintF("--SM-- GIA EXACT acceleration, t=%9.3e",t),"%7.5f ");
        }
    }
    else if( useTimeHistory )
    {
    
    // --- Compute the acceleration from the velocity time history ----
        if( true ) // *** NEW* June 29, 2017
        {
            Range Rx=mg.numberOfDimensions();
            
            if( velocityHistory.getNumberOfTimeLevels()<= 1 && twilightZoneFlow )
            {
        // ***** TEMP FIX FOR TZ -- use exact acceleration if we don't have enough time levels *wdh* March 14, 2018
                f(I1,I2,I3,C)=accel;
        
                if( true || debug() & 2 )
                {
                    printF("--SM-- getInterfaceAcceleration: Setting acceleration from EXACT known solution ***TEMP**"
                                  " at t=%9.3e.\n",t);
                }
            }
            else
            {
                RealArray acceleration(I1,I2,I3,Rx);
                const int numberOfTimeDerivatives=1;
                const int timeOrderOfAccuracy=parameters.dbase.get<int>("interfaceArrayEvolutionTimeAccuracy"); 
                int order= timeOrderOfAccuracy;
        // order=2;  // *wdh* May 8, 2018
                if( velocityHistory.getNumberOfTimeLevels()<= 1 )
                {
                    printF("--SM-- getInterfaceAcceleration: WARNING: setting acceleration ZERO at t=%.2e since there is "
                                  "only 1 time level in the velocity history\n",t);
                    f(I1,I2,I3,C)=0.;
                }
                else
                {
                    velocityHistory.eval( t, acceleration, numberOfTimeDerivatives,order );  
                    f(I1,I2,I3,C)=acceleration;
                }
                
                if( debug() & 2 )
                {
                    printF("--SM-- getInterfaceAccel: ACCELERATION computed from VELOCITY TIME HISTORY t=%9.3e order=%i correctionStage=%i "
                                  "*NEW*\n",order,correctionStage,t);
                    if( evalExactAcceleration )
                    {
                        real diff=max(fabs(acceleration-accel));
                        printF("--> Max error in accel = %9.3e, |accel|=%8.2e\n",diff,max(fabs(accel)));
                        if( diff>1e5 )
                        {
                            ::display(acceleration,"acceleration","%.2e ");
                        }
                        
                    }
                
                }
            }
            
        }
        else  // *OLD TIME HISTORY
        {
            
      // Look for the interface velocity 
            aString velocityTimeHistoryName;
            sPrintF(velocityTimeHistoryName,"velocityTimeHistoryG%iS%iA%i",grid,side,axis);
            std::vector<InterfaceData> & velocityTimeHistory =
                parameters.dbase.get<std::vector<InterfaceData> >(velocityTimeHistoryName);
                
            const int & numberOfInterfaceTimeLevels=parameters.dbase.get<int>("numberOfInterfaceTimeLevels");
            const int & cur = parameters.dbase.get<int>("currentInterfaceTimeLevel");
            assert (cur>=0 && cur<numberOfInterfaceTimeLevels );

            int prev  = ovmod(cur-1,numberOfInterfaceTimeLevels);
            int prev2 = ovmod(cur-2,numberOfInterfaceTimeLevels);
            real t2 = velocityTimeHistory[cur  ].t;
            real t1 = velocityTimeHistory[prev ].t;
            real t0 = velocityTimeHistory[prev2].t;
        
      // Compute the time derivative of the Lagrange polynomial: 
      //    v(t) = l0(t)*v0 + l1(t)*v1 + l2(t)*v2 
      // where 
      //   l0 = (t-t1)*(t-t2)/( (t0-t1)*(t0-t2) );
      //   l1 = (t-t2)*(t-t0)/( (t1-t2)*(t1-t0) );
      //   l2 = (t-t0)*(t-t1)/( (t2-t0)*(t2-t1) );

      // const real t0=time(prev2), t1=time(prev), t2=time(current);
            if( t0>=0.0 )
            {
        // enough previous levels exist for a 2nd-order approximation
                real dt0= t1-t0;
                real dt1= t2-t1;
                assert( dt0>0. && dt1>0. );
                assert( fabs(dt0-dt1)<dt0*1.e-6 );

                assert( fabs(t2-t) < REAL_EPSILON*10.*(1.+t) );


                RealArray & v2 = velocityTimeHistory[cur  ].u;
                RealArray & v1 = velocityTimeHistory[prev ].u;
                RealArray & v0 = velocityTimeHistory[prev2].u;
                Range V = v0.dimension(3);

                real l0t = (2.*t-(t1+t2))/( (t0-t1)*(t0-t2) );
                real l1t = (2.*t-(t2+t0))/( (t1-t2)*(t1-t0) );
                real l2t = (2.*t-(t0+t1))/( (t2-t0)*(t2-t1) );

                RealArray aI(I1,I2,I3,V);  //**TEMP**
                int orderOfApproximation=2;
                if( orderOfApproximation==2 )
                {
          // second-order approx.
                    aI(I1,I2,I3,V) = l0t*v0(I1,I2,I3,V) + l1t*v1(I1,I2,I3,V) + l2t*v2(I1,I2,I3,V);
                }
                else
                { // first order accurate approx. 
                    aI(I1,I2,I3,V) = (v2(I1,I2,I3,V) - v1(I1,I2,I3,V))/dt1;
                }
            
                f(I1,I2,I3,C)= aI(I1,I2,I3,V);

                if( true || debug() & 2 )
                    printF("--SM-- getInterfaceAccel: acceleration computed from velocity time history t=%9.3e (order=%i)",
                                  t,orderOfApproximation);
            
                if( debug() & 4 )
                {
                    ::display(aI,sPrintF("--SM-- GIA acceleration from velocity time history t=%9.3e (order=%i)",t,orderOfApproximation),"%7.5f ");
                    ::display(accel,sPrintF("--SM-- GIA EXACT acceleration                      t=%9.3e ",t,orderOfApproximation),"%7.5f ");
                }

            }
            else
            {

                if( evalExactAcceleration )
                {
                    f(I1,I2,I3,C)=accel;
                }
                else
                {
                    OV_ABORT("getInterfaceAccel - finish me");
                }
            
                if( true || debug() & 2 )
                {
                    printF("--SM-- getInterfaceAcceleration: Setting acceleration from EXACT known solution ***TEMP**"
                                  " at t=%9.3e.\n",t);
                    ::display(f(I1,I2,I3,C),sPrintF("GIA EXACT acceleration, t=%9.3e",t),"%7.5f ");
                }
            }
        
        
        
            int numLevelsSet=0;
            for( int i=0; i<numberOfInterfaceTimeLevels; i++ )
            {
                if( velocityTimeHistory[i].t >=0. ) numLevelsSet++;
            }
        
            printF("--SM-- GIA: velocity time history: t=%9.3e, cur=%i velocityTimeHistory[cur].t=%9.3e numLevelsSet=%i\n",
                          t,cur,velocityTimeHistory[cur].t,numLevelsSet);

      // assert( current>=0 && current<velocityTimeHistory.size() );
      // velocityTimeHistory[current].t=t;
      // velocityTimeHistory[current].u=uLocal(I1,I2,I3,V);

        }
    
    }
    else
    {
        OV_ABORT("--SM-- getInterfaceAcceleration -- FINISH ME ");
    }
    

    return 0;
}

//============================================================================================================
/// \brief Return the state (position, velocity,acceleration) of the interface for known solutions
///
/// \param stateOption (input) : specify which information to return
/// \param t (input) : time to eval the state
/// \param side,axis,grid : face
/// \param mg (input) : MappedGrid
/// \param I1,I2,I3,C (input) : assign points f(I1,I2,I3,C) 
/// \param f (output) : put acceleration here 
///
/// \Return value: 1= known solution was assigned, 0=no values were assigned
//============================================================================================================
int Cgsm::
getKnownInterfaceState( Parameters::DeformingBodyStateOptionEnum stateOption,
                                                const real t, const int side, const int axis, const int grid, 
                                                MappedGrid & mg, const Index & I1, const Index & I2, const Index & I3, const Range & C, 
                                                RealArray & state )
{
    int stateAssigned=0;
    
    const bool & twilightZoneFlow = parameters.dbase.get<bool >("twilightZoneFlow");

    const Parameters::KnownSolutionsEnum & knownSolution = 
        parameters.dbase.get<Parameters::KnownSolutionsEnum >("knownSolution"); 

    if( twilightZoneFlow )
    {
        printF("--SM-- getKnownInterfaceState: finish me\n");
    }
    else if( knownSolution==Parameters::userDefinedKnownSolution )
    {
        printF("--SM-- getKnownInterfaceState: stateOption=%i setting knownSolution=%i at t=%9.3e\n",
                      (int)stateOption,(int)knownSolution,t);

        int bodyNumber=0; // fix me for multiple deforming bodies
          	    
        parameters.getUserDefinedDeformingBodyKnownSolution( 
            bodyNumber,
            stateOption,
            t, grid, mg,I1,I2,I3,C, state );

        stateAssigned=1;
    }
    else
    {
    // Return default state
        printF("--SM-- getKnownInterfaceState: finish me\n");

    }
    

    return stateAssigned;
}
