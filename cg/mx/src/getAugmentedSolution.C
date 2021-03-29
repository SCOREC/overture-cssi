// This file automatically generated from getAugmentedSolution.bC with bpp.
//#define BOUNDS_CHECK
//#define OV_DEBUG

#include "Maxwell.h"
#include "PlotStuff.h"
#include "GL_GraphicsInterface.h"
#include "DialogData.h"
#include "ParallelUtility.h"
#include "display.h"
#include "DispersiveMaterialParameters.h"
// #include "BodyForce.h"

#define FOR_3D(i1,i2,i3,I1,I2,I3) int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)

#define FOR_3(i1,i2,i3,I1,I2,I3) I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase();  I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); for(i3=I3Base; i3<=I3Bound; i3++) for(i2=I2Base; i2<=I2Bound; i2++) for(i1=I1Base; i1<=I1Bound; i1++)

// Define macros for forcing functions:
//   (Ex).t = (1/eps)*[  (Hz).y ]
//   (Ey).t = (1/eps)*[ -(Hz).x ]
//   (Hz).t = (1/mu) *[ (Ex).y - (Ey).x ]

#define exTrue(x,y,t) sin(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc[0]
#define eyTrue(x,y,t) sin(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc[1]
#define hzTrue(x,y,t) sin(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc[5]

#define extTrue(x,y,t) (-twoPi*cc)*cos(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc[0]
#define eytTrue(x,y,t) (-twoPi*cc)*cos(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc[1]
#define hztTrue(x,y,t) (-twoPi*cc)*cos(twoPi*(kx*(x)+ky*(y)-cc*(t)))*pwc[5]

#define exLaplacianTrue(x,y,t) sin(twoPi*(kx*(x)+ky*(y)-cc*(t)))*(-(twoPi*twoPi*(kx*kx+ky*ky))*pwc[0])
#define eyLaplacianTrue(x,y,t) sin(twoPi*(kx*(x)+ky*(y)-cc*(t)))*(-(twoPi*twoPi*(kx*kx+ky*ky))*pwc[1])
#define hzLaplacianTrue(x,y,t) sin(twoPi*(kx*(x)+ky*(y)-cc*(t)))*(-(twoPi*twoPi*(kx*kx+ky*ky))*pwc[5])

// Here is a plane wave with the shape of a Gaussian
// xi = kx*(x)+ky*(y)-cc*(t)
// cc=  c*sqrt( kx*kx+ky*ky );
#define hzGaussianPulse(xi)  exp(-betaGaussianPlaneWave*((xi)*(xi)))
#define exGaussianPulse(xi)  hzGaussianPulse(xi)*(-ky/(eps*cc))
#define eyGaussianPulse(xi)  hzGaussianPulse(xi)*( kx/(eps*cc))

#define hzLaplacianGaussianPulse(xi)  ((4.*betaGaussianPlaneWave*betaGaussianPlaneWave*(kx*kx+ky*ky))*xi*xi-(2.*betaGaussianPlaneWave*(kx*kx+ky*ky)))*exp(-betaGaussianPlaneWave*((xi)*(xi)))
#define exLaplacianGaussianPulse(xi)  hzLaplacianGaussianPulse(xi,t)*(-ky/(eps*cc))
#define eyLaplacianGaussianPulse(xi)  hzLaplacianGaussianPulse(xi,t)*( kx/(eps*cc))

// 3D
//
//   (Ex).t = (1/eps)*[ (Hz).y - (Hy).z ]
//   (Ey).t = (1/eps)*[ (Hx).z - (Hz).x ]
//   (Ez).t = (1/eps)*[ (Hy).x - (Hx).y ]
//   (Hx).t = (1/mu) *[ (Ey).z - (Ez).y ]
//   (Hy).t = (1/mu) *[ (Ez).x - (Ex).z ]
//   (Hz).t = (1/mu) *[ (Ex).y - (Ey).x ]

// ****************** finish this -> should `rotate' the 2d solution ****************

#define exTrue3d(x,y,z,t) sin(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*pwc[0]
#define eyTrue3d(x,y,z,t) sin(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*pwc[1]
#define ezTrue3d(x,y,z,t) sin(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*pwc[2]

#define extTrue3d(x,y,z,t) (-twoPi*cc)*cos(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*pwc[0]
#define eytTrue3d(x,y,z,t) (-twoPi*cc)*cos(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*pwc[1]
#define eztTrue3d(x,y,z,t) (-twoPi*cc)*cos(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*pwc[2]



#define hxTrue3d(x,y,z,t) sin(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*pwc[3]
#define hyTrue3d(x,y,z,t) sin(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*pwc[4]
#define hzTrue3d(x,y,z,t) sin(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*pwc[5]

#define exLaplacianTrue3d(x,y,z,t) sin(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*(-(twoPi*twoPi*(kx*kx+ky*ky+kz*kz))*pwc[0])
#define eyLaplacianTrue3d(x,y,z,t) sin(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*(-(twoPi*twoPi*(kx*kx+ky*ky+kz*kz))*pwc[1])
#define ezLaplacianTrue3d(x,y,z,t) sin(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*(-(twoPi*twoPi*(kx*kx+ky*ky+kz*kz))*pwc[2])

#define hxLaplacianTrue3d(x,y,z,t) sin(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*(-(twoPi*twoPi*(kx*kx+ky*ky+kz*kz))*pwc[3])
#define hyLaplacianTrue3d(x,y,z,t) sin(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*(-(twoPi*twoPi*(kx*kx+ky*ky+kz*kz))*pwc[4])
#define hzLaplacianTrue3d(x,y,z,t) sin(twoPi*(kx*(x)+ky*(y)+kz*(z)-cc*(t)))*(-(twoPi*twoPi*(kx*kx+ky*ky+kz*kz))*pwc[5])



//==================================================================================================
// Evaluate Tom Hagstom's exact solution defined as an integral of Guassian sources
// 
// OPTION: OPTION=solution or OPTION=error OPTION=bounary to compute the solution or the error or
//     the boundary condition
//
//==================================================================================================


//==================================================================================================
// The DEFINE_GF_MACRO is a helper for EXTRACT_GFP and sets up the cpp macros for a given field
//
//==================================================================================================
//==================================================================================================

//==================================================================================================
// The EXTRACT_GFP macro extracts gridfunction pointers and array bounds.
//                 We use these to write code that works in 2/3D for both the nfdtd and
//                 dsi schemes.
//
// The macro expects the user to have the following variables defined in the enclosing scope:
//    Index I1,I2,I3 - used to get the index ranges for each grid
//    CompositeGrid &cg - the composite grid used by the fields
//    int grid      - the current grid to setup the pointers for
//    int i1,i2,i3  - grid indices into the arrays at the appropriate centering
//    MappedGrid *Maxwell:mgp - or -
//          CompositeGrid *Maxwell cgfields -or- CompositeGrid *Maxwell dsi_cgfields
//
// The macro defines:
//    cpp macros:
//               UH{XYZ}(i0,i1,i2) - access the current h field at i1,i2,i3 with the appropriate centering
//               UE{XYZ}(i0,i1,i2) - access the current e field at i1,i2,i3 with the appropriate centering        
//               UMH{XYZ}(i0,i1,i2) - the h field at the previous timestep
//               UME{XYZ}(i0,i1,i2) - the e field at the previous timestep
//               UNH{XYZ}(i0,i1,i2) - the h field at the next timestep
//               UNE{XYZ}(i0,i1,i2) - the h field at the next timestep
//
//               ERRH{XYZ}(i0,i1,i2) - acces the h field error gridfunction
//               ERRE{XYZ}(i0,i1,i2) - acces the e field error gridfunction
//
//               XEP(i0,i1,i2,i3) - coordinates of e centering
//               XHP(i0,i1,i2,i3) - coordinates of h centering
// 
//    variables:
//               MappedGrid &mg - the current mapped grid (cg[grid])
//
//               const bool isStructured - true for structured grids
//               const bool isRectangular - true for rectangular grids
//
//               realMappedGridFunction uh - view of the current h or h.n field
//               realMappedGridFunction ue - view of the current e or e.n field
//               realMappedGridFunction umh - view of the previous h or h.n field
//               realMappedGridFunction ume - view of the previous e or e.n field
//               realMappedGridFunction unh - view of the next h or h.n field
//               realMappedGridFunction une - view of the next e or e.n field
//               realMappedGridFunction errh - view of the h field error
//               realMappedGridFunction erre - view of the e field error
//               realArray xe - view of the x coordinates at the e centering
//               realArray xh - view of the x coordinates at the h centering
//               realArray ye - view of the y coordinates at the e centering
//               realArray yh - view of the y coordinates at the h centering
//               realArray ze - view of the z coordinates at the e centering
//               realArray zh - view of the z coordinates at the h centering
//               realArray xce - coordinates of the e centering
//               realArray xch - coordinates of the h centering
//               realArray emptyArray - used for setting references to things we don't need
//               real *uhp - data pointer for the current h or h.n field
//               real *uep - data pointer for the current h or e.n field
//               real *umhp - data pointer for the previous h or h.n field
//               real *umep - data pointer for the previous h or e.n field
//               real *unhp - data pointer for the next h or h.n field
//               real *unep - data pointer for the next h or e.n field
//               real *xep - data pointer for the coordinates at the e centering
//               real *xhp - data pointer for the coordinates at the h centering
//
//               real dx[3] - dx in each direction for rectangular grids  (={0,0,0} if !isRectangular)
//               real xab[2][3] - coordinate bounds for rectangular grids (={ {0,0},.. } if !isRectangular)
//
//               int uhDim0,uhDim1,uhDim2 - array dimensions for the e gridfunctions
//               int ueDim0,ueDim1,ueDim2 - array dimensions for the h gridfunctions
//               int xeDim0,xeDim1,xeDim2 - array dimensions for the e centering coordinates
//               int xhDim0,xhDim1,xhDim2 - array dimensions for the h centering coordinates
//
// KNOWN ASSUMPTIONS:  * gridFunctions for the same variable at different time levels have the same
//                                   raw data sizes
//                     * there are unrecognized and perhaps subtle assumptions being made
//         
// OPTION: 
//==================================================================================================
//==================================================================================================
//==================================================================================================



//==================================================================================================
// This bpp macro undefs the cpp macros defined by EXTRACT_GFP
//               UH(i0,i1,i2) - access the current h field at i1,i2,i3 with the appropriate centering
//               UE(i0,i1,i2) - access the current e field at i1,i2,i3 with the appropriate centering        
//               UMH(i0,i1,i2) - the h field at the previous timestep
//               UME(i0,i1,i2) - the e field at the previous timestep
//               UNH(i0,i1,i2) - the h field at the next timestep
//               UNE(i0,i1,i2) - the h field at the next timestep
//               XEP(i0,i1,i2,i3) - coordinates of e centering
//               XHP(i0,i1,i2,i3) - coordinates of h centering
// OPTION: 
//==================================================================================================


// =========================================================================================
// Macro: assign the field names in the grid function
// =========================================================================================

// =========================================================================================
// Macro: initialize the variables for plotting the polarization
// =========================================================================================


// =========================================================================================
// Macro: assign the names for the polarization variables 
// =========================================================================================


// =========================================================================================
// Macro: save the polarization in the augmented solution
// =========================================================================================


// =========================================================================================
// Macro: Initialize nonlinear variables for plotting. 
// =========================================================================================


// =========================================================================================
// Macro: assign the names for the nonlinear variables 
// =========================================================================================


// =========================================================================================
// Macro: save the polarization in the augmented solution
// =========================================================================================



// =============================================================================================
/// \brief Create a grid function that holds all the things we can plot.
/// \param t (input) : if t<0 then only fill the component names into the grid function v.
// =============================================================================================
realCompositeGridFunction& Maxwell::
getAugmentedSolution(int current, realCompositeGridFunction & v, const real t)
{

    assert( cgp!=NULL );
    CompositeGrid & cg = *cgp;
    const int numberOfDimensions = cg.numberOfDimensions();
    
  // const int numberOfComponents=mgp==NULL ? cgfields[current][0].getLength(3) : fields[current].getLength(3);

    int numberOfComponents;
        numberOfComponents= dbase.get<int>("numberOfComponents"); // *new way* Sept 01, 2017 *wdh*

    const bool saveErrors = plotErrors && !(errp==NULL && cgerrp==NULL);
    const bool saveDissipation =  plotDissipation && (( (artificialDissipation>0. || artificialDissipationCurvilinear>0.) 
                                                                && (method==nfdtd || method==bamx)) || dissipation || cgdissipation);

    const int & solveForAllFields = dbase.get<int>("solveForAllFields");

    Range all;
    const bool saveDsiDiss= saveDissipation && method==dsiMatVec && cg.numberOfDimensions()==3;
    
  // Determine the number of components to plot and the component numbers for the errors, etc.
  //    nErr : component where the error is stored
  //    ndd  : component where the dissipation is stored
    int numberToPlot=numberOfComponents;                  // save fields
    int nErr=numberToPlot;    numberToPlot += numberOfComponents*int(saveErrors);
    int ndd=numberToPlot;     numberToPlot += numberOfComponents*int(saveDissipation);
                                                        numberToPlot += cg.numberOfDimensions()*int( saveDsiDiss ); 
    int nVarDis=numberToPlot; numberToPlot += int(useVariableDissipation);
    int nDivE=numberToPlot;   numberToPlot += int(plotDivergence); 

    int nDivH=-1;
    int plotDivH = plotDivergence && (method==yee || solveForMagneticField ) && numberOfDimensions==3;
    plotDivH = plotDivH || ( method==bamx && (solveForAllFields==1 || numberOfDimensions==3) );
    if( plotDivH ) 
    { // plot div(H) or div(B) too
        nDivH=numberToPlot;  numberToPlot +=1;
    }
    if( method!=nfdtd && method!=yee && method!=sosup && method!=bamx )
    {
        numberToPlot += 2;  // something for Kyle
    }
    int nEdiss=numberToPlot;  numberToPlot += (cg.numberOfDimensions()+1)*int(e_cgdissipation ? 1 : 0);
    int nRho=numberToPlot;    numberToPlot += int(plotRho); 
    int nEnergyDensity=numberToPlot; numberToPlot += int(plotEnergyDensity);
    int nIntensity=numberToPlot; numberToPlot += int(plotIntensity);

  // There are 2 components of the harmonic field, Er and Ei for each component of E. 
    int nHarmonicE=numberToPlot;  numberToPlot += 2*(cg.numberOfDimensions())*int(plotHarmonicElectricFieldComponents);

    bool plotCurlE=false;     // for testing plot curl( E_known )
    if( method==yee && false )
        plotCurlE=true;
    int nCurlE = numberToPlot; numberToPlot += 2*(1 + 2*(numberOfDimensions-2))*int(plotCurlE);
    
    const int & maxNumberOfPolarizationComponents = parameters.dbase.get<int>("maxNumberOfPolarizationComponents");
  // total number of polarization components per grid 
    const IntegerArray & totalNumberOfPolarizationComponents =
        parameters.dbase.get<IntegerArray>("totalNumberOfPolarizationComponents");



  // Initialize the variables for plotting the polarization
    // bool plotPolarization = false && method==bamx && solveForAllFields && dispersionModel!=noDispersion; // *** FIX ME ***
    // const int numPolarizationVectors=6;  // [Px,Py,Pz, Mx,My,Mz]
        bool plotPolarization = dispersionModel!=noDispersion; 
        bool plotMagnetization = false;
        int numPolarizationVectors=-1;
        if( method==bamx )
        {
            plotMagnetization=true;
            if( solveForAllFields )
                numPolarizationVectors=6;   // [Px,Py,Pz, Mx,My,Mz]
            else
                numPolarizationVectors=3;
        }
        else
        {
            numPolarizationVectors=numberOfDimensions;  // [Px,Py,Pz]
        }
        bool plotPolarizationErrors = false && plotPolarization && saveErrors;  // finish me 
        const int nPolarization=numberToPlot;
        const int nPolarizationErr=nPolarization+numPolarizationVectors;
        if(  plotPolarization )
            numberToPlot += numPolarizationVectors; // Plot polarization and magnetization vectors
        if( plotPolarizationErrors )
            numberToPlot += numPolarizationVectors;
    // Plot each individual compnent of the polarization  
        bool plotPolarizationComponents = dbase.get<bool>("plotPolarizationComponents") && dispersionModel!=noDispersion; 
    // total number of polarization components per grid 
    // ::display(totalNumberOfPolarizationComponents,"totalNumberOfPolarizationComponents");
    // const int maxNumberOfPolarizationVectors = max(totalNumberOfPolarizationComponents);
        int maxNumberOfPolarizationVectors=0;
        const int nPolarComponent= numberToPlot; // start saving polarization components in this index 
        if( plotPolarizationComponents )
        {
            for( int domain=0; domain<cg.numberOfDomains(); domain++ )
            {
                const DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
                const int numberOfPolarizationVectors=dmp.getNumberOfPolarizationVectors();
                maxNumberOfPolarizationVectors=max(maxNumberOfPolarizationVectors,numberOfPolarizationVectors);
            }
            numberToPlot += maxNumberOfPolarizationVectors*numberOfDimensions;
        }
        bool plotPolarizationComponentErrors = plotPolarizationComponents && saveErrors;
        const int nPolarizationComponentErrors=numberToPlot;
        if( plotPolarizationComponentErrors )
        {
      // plot errors in the individual polarization components
            numberToPlot += maxNumberOfPolarizationVectors*numberOfDimensions;
        }
        

  // Initialize nonlinear variables for plotting
        bool plotNonlinearComponents = dbase.get<bool>("plotNonlinearComponents");
        const int nNonlinear = numberToPlot; // start saving nonlinear components in this index 
        bool isNonlinear=false;
        int maxNumberOfAtomicLevels=0;
        if( plotNonlinearComponents )
        {
      // ---- plot components of the nonlinear model ----
            for( int domain=0; domain<cg.numberOfDomains(); domain++ )
            {
                const DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
                isNonlinear = isNonlinear || dmp.isNonlinearMaterial();
                const int numberOfAtomicLevels = dmp.getNumberOfAtomicLevels();
                maxNumberOfAtomicLevels=max(maxNumberOfAtomicLevels,numberOfAtomicLevels);
            }
            if( isNonlinear )
                numberToPlot += maxNumberOfAtomicLevels;
            else
                plotNonlinearComponents=false;   // there are non nonlinear material domains 
        }
        bool plotNonlinearComponentErrors = plotNonlinearComponents && saveErrors;
        const int nNonlinearComponentErrors=numberToPlot;
        if( plotNonlinearComponentErrors )
        {
            numberToPlot += maxNumberOfAtomicLevels;
        }



  // we build a grid function with more components (errors, dissipation) for plotting
    v.updateToMatchGrid(cg,all,all,all,numberToPlot);

    v=0;

  // Assign the field names in the augemented solution grid function
        if( method==nfdtd || method==yee || method==sosup || method==bamx ) 
        {
            for( int n=0; n<numberOfComponents; n++ )
            {
                if( mgp!=NULL )
                {
                    MappedGrid & mg = *mgp;
                    if( mg.getGridType()==MappedGrid::structuredGrid )
                    {
                        if ( method==nfdtd || method==sosup || method==bamx || ( n<fields[current].getLength(3) ) )
                            v.setName(fields[current].getName(n),n);
                        else if ( n<fields[current].getLength(3) ) 
                            v.setName(fields[current+numberOfTimeLevels].getName(n-fields[current].getLength(3)),n);
                        if( saveErrors )
                            v.setName(errp->getName(n),n+numberOfComponents);
                    }
                }
                else
                {
          // *wdh* v.setName(cgfields[current].getName(n),n);
                    v.setName(getCGField(HField,current).getName(n),n);
                    if( saveErrors )
                        v.setName(cgerrp->getName(n),n+numberOfComponents);
                    if( saveDissipation )
                        v.setName(cgdissipation->getName(n),n+ndd);
                }
            }
        }
        else
        {
            if( cg.numberOfDimensions()==2 )
            {
                int i=3;
                v.setName("Hz",0);
                v.setName("Ex",1);
                v.setName("Ey",2);
                if( method==dsiMatVec )
                {
                    v.setName("E.n",3);
                    i=4;
                }
                if ( (dissipation ||cgdissipation )&& plotDissipation)
                {
                    v.setName("Hz dissp",i++);
                    v.setName("E.n dissp",i++);
                    v.setName("Ex dissp",i++);
                    v.setName("Ey dissp",i++);
                }
                if( saveErrors )
                {
                    v.setName("Hz-err",i++);
                    v.setName("Ex-err",i++);
                    v.setName("Ey-err",i);
                }
            }
            else
            {
                v.setName("Hx",hx);
                v.setName("Hy",hy);
                v.setName("Hz",hz);
                v.setName("Ex",ex+3);
                v.setName("Ey",ey+3);
                v.setName("Ez",ez+3);
                v.setName("H.n",ez+4);
                v.setName("E.n",ez+5);
                int i=8;
                if ( (dissipation ||cgdissipation) && plotDissipation)
                {
                    v.setName("H.n dissp",i++);
                    v.setName("Hx dissp",i++);
                    v.setName("Hy dissp",i++);
                    v.setName("Hz dissp",i++);
                    v.setName("E.n dissp",i++);
                    v.setName("Ex dissp",i++);
                    v.setName("Ey dissp",i++);
                    v.setName("Ez dissp",i++);
                }
                if( saveErrors )
                {
                    v.setName("Hx-err",i++);
                    v.setName("Hy-err",i++);
                    v.setName("Hz-err",i++);
                    v.setName("Ex-err",i++);
                    v.setName("Ey-err",i++);
                    v.setName("Ez-err",i);
                }
            }
        }
        if( plotDivergence && (method==nfdtd || method==yee || method==sosup || method==bamx) )
        {
            if( method==bamx )
            {
                v.setName("div(D)",nDivE);
                if( nDivH>=0 )
                    v.setName("div(B)",nDivH);
            }
            else
            {
                v.setName("div(E)",nDivE);
                if( nDivH>=0 )
                    v.setName("div(H)",nDivH);
            }
        }
        if( plotCurlE && (method==nfdtd || method==yee || method==sosup || method==bamx) )
        {
            if( numberOfDimensions==3 )
            {
                v.setName("curlExr",nCurlE  );
                v.setName("curlEyr",nCurlE+1);
                v.setName("curlEzr",nCurlE+2);
                v.setName("curlExi",nCurlE+3);
                v.setName("curlEyi",nCurlE+4);
                v.setName("curlEzi",nCurlE+5);
            }
        }
      

  // Assign the names for the polarization variables 
        if( plotPolarization )
        {
            v.setName("Px",nPolarization+0);
            v.setName("Py",nPolarization+1);
            if( numberOfDimensions==3 || method==bamx )   // ** fix me ***
                v.setName("Pz",nPolarization+2);              
            if( plotMagnetization )
            {
                v.setName("Mx",nPolarization+3);
                v.setName("My",nPolarization+4);
                v.setName("Mz",nPolarization+5);
            }
        }
        if( plotPolarizationErrors )
        {
            v.setName("Px error",nPolarizationErr+0);
            v.setName("Py error",nPolarizationErr+1);
            if( numberOfDimensions==3 || method==bamx )   // ** fix me ***
                v.setName("Pz error",nPolarizationErr+2);
            if( plotMagnetization )
            {
                v.setName("Mx error",nPolarizationErr+3);
                v.setName("My error",nPolarizationErr+4);
                v.setName("Mz error",nPolarizationErr+5);
            }
        }
        if( plotPolarizationComponents )
        {
      // --- name individual polarization components ---
            for( int iv=0; iv<maxNumberOfPolarizationVectors; iv++ )
            {
                for( int dir=0; dir<numberOfDimensions; dir++ )
                {
          // Names: P1x, P1y, P1z, P2x, P2y, P2z, ....
                    v.setName(sPrintF("P%d%s",iv+1,(dir==0 ? "x" : dir==1 ? "y" : "z")), nPolarComponent+dir+numberOfDimensions*(iv));
                }
            }
        }
        if( plotPolarizationComponentErrors )
        {
      // --- name individual polarization components ---
            for( int iv=0; iv<maxNumberOfPolarizationVectors; iv++ )
            {
                for( int dir=0; dir<numberOfDimensions; dir++ )
                {
          // Names: P1xErr, P1yErr, P1zErr, P2xErr, P2yErr, P2zErr, ....
                    v.setName(sPrintF("P%d%sErr",iv+1,(dir==0 ? "x" : dir==1 ? "y" : "z")),
                                        nPolarizationComponentErrors +dir+numberOfDimensions*(iv));
                }
            }
        }
      
  // Assign names for the nonlinear variables 
        if( plotNonlinearComponents )
        {
      // --- name components of the nonlinear model ---
            for( int n=0; n<maxNumberOfAtomicLevels; n++ )
            {
        // Names: N0,N1,N2,...
                v.setName(sPrintF("N%d",n), nNonlinear+n);
            }
        }
        if( plotNonlinearComponentErrors )
        {
            for( int n=0; n<maxNumberOfAtomicLevels; n++ )
            {
        // Names: N0Err, N1Err, N2Err,...
                v.setName(sPrintF("N%dErr",n), nNonlinearComponentErrors+n);
            }
        }
    
    


    if( useVariableDissipation )
        v.setName("varDis",nVarDis);
    if( plotRho )
        v.setName("rho",nRho);
    if( plotEnergyDensity )
        v.setName("energyDensity",nEnergyDensity);

    if( plotIntensity )
    {
        v.setName("intensity",nIntensity);
    }
        
    if( plotHarmonicElectricFieldComponents )
    {
        v.setName("Exr",nHarmonicE+0);
        v.setName("Exi",nHarmonicE+1);
        v.setName("Eyr",nHarmonicE+2);
        v.setName("Eyi",nHarmonicE+3);
        if( numberOfDimensions==3 )
        {
            v.setName("Ezr",nHarmonicE+4);
            v.setName("Ezi",nHarmonicE+5);
        }
    }
        
    if( t<0. )
    {
    // in this case we only assign the component names and return 
        return v;
    }


//   if( plotIntensity || plotHarmonicElectricFieldComponents )
//   {
//     if( false && intensityOption==1 )
//     {
//       // compute the intensity using current and prev values
//       int stepNumber=0;
//       real nextTimeToPlot=0.;
//       real dt=deltaT; // check this 
//       computeIntensity(current,t,dt,stepNumber,nextTimeToPlot);
//     }
//   }
  // printF(" plot: cg.numberOfComponentGrids() = %i \n",cg.numberOfComponentGrids());
    

    divEMax=0.;

    if( method==yee )
    {
    // compute node centered fields for plotting -- this will fill in v ---
        int option=3;
        int iparam[5] = { nDivE,nDivH,0,0,0 }; // 
        getValuesFDTD( option, iparam, current, t, deltaT, &v );
    // ::display(v[0],"v after getValuesFDTD","%5.2f");
        if( plotDivergence )
        {
            option=2; // compute div(E) ( and div(H) in 3D)
            getValuesFDTD( option, iparam, current, t, deltaT, &v );
        }
        if( plotCurlE )
        {
            option=4; // compute curl(E)
            iparam[0]=nCurlE;
            getValuesFDTD( option, iparam, current, t, deltaT, &v );
        }
        
    }


    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = mgp!=NULL ? *mgp : cg[grid];

        realMappedGridFunction & u = mgp!=NULL ? fields[current] : getCGField(HField,current)[grid];
        realMappedGridFunction & vg = v[grid];
    
        OV_GET_SERIAL_ARRAY(real,u,uLocal);
        OV_GET_SERIAL_ARRAY(real,vg,vLocal);

        Range N=numberOfComponents;
        if( method==yee )
        {
      // this is done above
        }
        else if( method==nfdtd  || method==sosup || method==bamx )
        {
            vLocal(all,all,all,N)=uLocal(all,all,all,N); // for now make a copy *** fix this **
        }
        else
        {
            Range N1 = fields[current].getLength(3);
            Range N2 =fields[current+numberOfTimeLevels].getLength(3);
#ifdef USE_PPP
            realSerialArray f1Local; getLocalArrayWithGhostBoundaries(fields[current],f1Local);
            realSerialArray f2Local; getLocalArrayWithGhostBoundaries(fields[current+numberOfTimeLevels],f2Local);
#else
            const realSerialArray & f1Local = fields[current];
            const realSerialArray & f2Local = fields[current+numberOfTimeLevels];
#endif
            vLocal(all,all,all,N1) = f1Local(all,all,all,N1);
            vLocal(all,all,all,N2) = f2Local(all,all,all,N2);
        }


        const bool & solveForScatteredField = dbase.get<bool>("solveForScatteredField");
        
        if( ( solveForScatteredField && plotTotalField ) &&
                cg.domainNumber(grid)==0 )  //   assumes domain 0 is the exterior domain
        {
      // *** NOTE: only add plane wave to the outer domain

            const real cc= c*sqrt( kx*kx+ky*ky+kz*kz );
                
            mg.update(MappedGrid::THEcenter | MappedGrid::THEvertex);  // *** fix for rectangular ***
                
      // // subtract off or add on the the incident field
      // const real pm = plotScatteredField ? 1. : -1.;

            const real pm = 1.;  // add on the incident field
            
            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
            getIndex(mg.dimension(),I1,I2,I3);
            const int includeGhost=1;
            bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3,includeGhost); 


#ifdef USE_PPP
            realSerialArray xLocal; getLocalArrayWithGhostBoundaries(mg.center(),xLocal);
#else
            const realSerialArray & xLocal = mg.center();
#endif

            if( ok )
            {
                if( mg.numberOfDimensions()==2 )
                {
                    const realSerialArray & x = xLocal(I1,I2,I3,0);
                    const realSerialArray & y = xLocal(I1,I2,I3,1);

                    vLocal(I1,I2,I3,ex)-=pm*exTrue(x,y,t);
                    vLocal(I1,I2,I3,ey)-=pm*eyTrue(x,y,t);
                    vLocal(I1,I2,I3,hz)-=pm*hzTrue(x,y,t);
                    if( method==sosup )
                    {
                        vLocal(I1,I2,I3,ext)-=pm*extTrue(x,y,t);
                        vLocal(I1,I2,I3,eyt)-=pm*eytTrue(x,y,t);
                        vLocal(I1,I2,I3,hzt)-=pm*hztTrue(x,y,t);
                    }
                    
                }
                else
                {
                    const realSerialArray & x = xLocal(I1,I2,I3,0);
                    const realSerialArray & y = xLocal(I1,I2,I3,1);
                    const realSerialArray & z = xLocal(I1,I2,I3,2);

                    if( solveForElectricField )
                    {
                        vLocal(I1,I2,I3,ex)-=pm*exTrue3d(x,y,z,t);
                        vLocal(I1,I2,I3,ey)-=pm*eyTrue3d(x,y,z,t);
                        vLocal(I1,I2,I3,ez)-=pm*ezTrue3d(x,y,z,t);
                        if( method==sosup )
                        {
                            vLocal(I1,I2,I3,ext)-=pm*extTrue3d(x,y,z,t);
                            vLocal(I1,I2,I3,eyt)-=pm*eytTrue3d(x,y,z,t);
                            vLocal(I1,I2,I3,ezt)-=pm*eztTrue3d(x,y,z,t);
                        }
                        
                    }

                }
            } // end if ok 
                
        } // end if( plotScatteredField || plotTotalField )
            
        if( plotEnergyDensity )
        {
            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
            getIndex(mg.dimension(),I1,I2,I3);
            const int includeGhost=1;
            bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3,includeGhost); 

            if( ok )
            {
                c = cGrid(grid);
                eps = epsGrid(grid);
                mu = muGrid(grid);

                if( mg.numberOfDimensions()==2 )
                {
          // vLocal(I1,I2,I3,nEnergyDensity)= eps*( SQR(uLocal(I1,I2,I3,ex))+SQR(uLocal(I1,I2,I3,ey)) );
                    vLocal(I1,I2,I3,nEnergyDensity)= ( (.5*eps)*( SQR(uLocal(I1,I2,I3,ex))+SQR(uLocal(I1,I2,I3,ey)) )+
                                                                                          (.5*mu )*( SQR(uLocal(I1,I2,I3,hz)) ) );
                }
                else
                {
                    Overture::abort("finish me -- we need H here");
                    vLocal(I1,I2,I3,nEnergyDensity)= eps*( SQR(uLocal(I1,I2,I3,ex))+SQR(uLocal(I1,I2,I3,ey))+
                                                                                                  SQR(uLocal(I1,I2,I3,ez)) );
                }
            } // end if ok 
                
        } // end if( plotEnergyDensity )
            
        
    // --- save the polarization in the augmented solution ---
        {
            if( plotPolarization || plotPolarizationComponents )
            {
        // --- plot P and M (sums) ----
                Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                getIndex(mg.dimension(),I1,I2,I3);
                const int includeGhost=1;
                if( method!=bamx )
                {
          // --- Isotropic Maxwell : plot total polarization vector and individual components  ----
                    const int domain = cg.domainNumber(grid);
                    const DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
                    if( dmp.isDispersiveMaterial() )
                    {
                        realMappedGridFunction & p = getDispersionModelMappedGridFunction( grid,current );
                        OV_GET_SERIAL_ARRAY(real,p,pLocal);
                        bool ok = ParallelUtility::getLocalArrayBounds(p,pLocal,I1,I2,I3,includeGhost);
                        if( ok )
                        {
                            const int numberOfPolarizationVectors = dmp.getNumberOfPolarizationVectors();
                            printF("getAugmentedSolution: grid=%d, save polarization, numberOfPolarizationVectors=%d\n",grid,numberOfPolarizationVectors);
                            for( int dir=0; dir<numberOfDimensions; dir++ )
                            {
                                if( plotPolarization )
                                    vLocal(I1,I2,I3,nPolarization+dir) = 0.;  // initialize sum 
                                for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
                                {
                                    const int pc = dir+ iv*numberOfDimensions; // component 
                                    if( plotPolarizationComponents )
                                        vLocal(I1,I2,I3,nPolarComponent+pc) = pLocal(I1,I2,I3,pc);  // individual component 
                  // --- sum the component polarization vectors ---
                                    if( plotPolarization )
                                        vLocal(I1,I2,I3,nPolarization+dir) += pLocal(I1,I2,I3,pc);   // sum 
                                }
                            }
                            if( plotPolarizationComponentErrors )
                            {
                // --- plot errors in individual components of P ---
                                bool getErrorGridFunction=true;
                                realMappedGridFunction & pErr = getDispersionModelMappedGridFunction( grid,current,getErrorGridFunction  );
                                OV_GET_SERIAL_ARRAY(real,pErr,pErrLocal);
                                for( int dir=0; dir<numberOfDimensions; dir++ )
                                {
                                    for( int iv=0; iv<numberOfPolarizationVectors; iv++ )
                                    {
                                        const int pc = dir+ iv*numberOfDimensions; // component 
                                        vLocal(I1,I2,I3,nPolarizationComponentErrors+pc) = pErrLocal(I1,I2,I3,pc);
                                    }
                                }
                            }
                            if( plotPolarizationComponents )
                            {
                // zero out any components that do not live on this grid 
                                for( int pc=numberOfPolarizationVectors; pc<maxNumberOfPolarizationComponents*numberOfDimensions; pc++ )
                                    vLocal(I1,I2,I3,nPolarComponent+pc)=0.;
                                if( plotPolarizationComponentErrors )
                                {
                                    for( int pc=numberOfPolarizationVectors; pc<maxNumberOfPolarizationComponents*numberOfDimensions; pc++ )
                                        vLocal(I1,I2,I3,nPolarizationComponentErrors+pc)=0.;
                                }
                            }
                        }
                    }
                    else
                    {
            // Set P=0 in a non-dispersive region
                        bool ok = ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3,includeGhost);
                        if( ok )
                        {
                            if( plotPolarization )
                            {
                                for( int ec=0; ec<numPolarizationVectors; ec++ )
                                    vLocal(I1,I2,I3,nPolarization+ec) = 0.;
                            }
                            if( plotPolarizationComponents )
                            {
                                for( int pc=0; pc<maxNumberOfPolarizationComponents*numberOfDimensions; pc++ )
                                    vLocal(I1,I2,I3,nPolarComponent+pc)=0.;
                            }
                        }
                    }
                }
                else if( method==bamx ) 
                {
          // ------ BA MAXWELL -----
          // printF("BAMX: Evaluate [Px,Py,Pz] [Mx,My,Mz] for plotting\n");
                    realMappedGridFunction & p = getDispersionModelMappedGridFunction( grid,current );
                    OV_GET_SERIAL_ARRAY(real,p,pLocal);
                    bool ok = ParallelUtility::getLocalArrayBounds(p,pLocal,I1,I2,I3,includeGhost); 
          // --- Errors in P ----
                    bool getErrorGridFunction=plotPolarizationErrors;
                    realMappedGridFunction & pErr = getDispersionModelMappedGridFunction( grid,current,getErrorGridFunction );
                    OV_GET_SERIAL_ARRAY(real,pErr,pErrLocal);
          // printF("pLocal: [%d,%d][%d,%d][%d,%d][%d,%d]\n",
          //        pLocal.getBase(0),pLocal.getBound(0),
          //        pLocal.getBase(1),pLocal.getBound(1),
          //        pLocal.getBase(2),pLocal.getBound(2),
          //        pLocal.getBase(3),pLocal.getBound(3));
          // printF("vLocal: [%d,%d][%d,%d][%d,%d][%d,%d]\n",
          //        vLocal.getBase(0),vLocal.getBound(0),
          //        vLocal.getBase(1),vLocal.getBound(1),
          //        vLocal.getBase(2),vLocal.getBound(2),
          //        vLocal.getBase(3),vLocal.getBound(3));
                    if( ok )
                    {
                        std::vector<DispersiveMaterialParameters> & dmpVector = 
                            dbase.get<std::vector<DispersiveMaterialParameters> >("materialRegionParameters");
            // printF("numberOfMaterialRegions=%d, dmpVector.size()=%d\n",numberOfMaterialRegions,dmpVector.size());
                        assert( numberOfMaterialRegions==dmpVector.size() );
                        if( numberOfMaterialRegions>1 )
                        {
              // assert( pBodyMask!=NULL );
              // const IntegerArray & matMask = *pBodyMask;  // material index 
                            intCompositeGridFunction & materialMask = parameters.dbase.get<intCompositeGridFunction>("materialMask");
                            OV_GET_SERIAL_ARRAY(int,materialMask[grid],matMask);
              // ------ make an array of Np(6,6) for the different materials ----
                            IntegerArray Npa(6,6,numberOfMaterialRegions);
                            Range M6=6;
                            for( int mr=0; mr<numberOfMaterialRegions; mr++ )
                            {
                                DispersiveMaterialParameters & dmp = dmpVector[mr];             
                                const IntegerArray & Np = dmp.getBianisotropicNp();
                // ::display(Np,"Np","%3i ");
                                Npa(M6,M6,mr) = Np(M6,M6);
                            }
                            int i1,i2,i3;
                            FOR_3D(i1,i2,i3,I1,I2,I3)
                            {
                                const int mr = matMask(i1,i2,i3);
                                assert( mr>=0 && mr<numberOfMaterialRegions );
                                DispersiveMaterialParameters & dmp = dmpVector[mr]; 
                                int pc=0;
                                for( int k1=0; k1<6; k1++ )
                                {
                                    int ec=k1;
                                    vLocal(i1,i2,i3,nPolarization+ec)=0.;
                                    if( plotPolarizationErrors )
                                        vLocal(i1,i2,i3,nPolarizationErr+ec)=0.;
                                    for( int k2=0; k2<6; k2++ )
                                    {
                                        for( int n=0; n<Npa(k1,k2,mr); n++ )
                                        {
                                            vLocal(i1,i2,i3,nPolarization+ec) += pLocal(i1,i2,i3,pc);
                                            if( plotPolarizationErrors )
                                                vLocal(i1,i2,i3,nPolarizationErr+ec) += pErrLocal(i1,i2,i3,pc);
                                            pc+=2;   // we store p and pt so increment by 2
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
              // single material region
                            for( int mr=0; mr<numberOfMaterialRegions; mr++ )
                            {
                                DispersiveMaterialParameters & dmp = dmpVector[mr]; 
                                const IntegerArray & Np = dmp.getBianisotropicNp();
                                int pc=0;
                                for( int k1=0; k1<6; k1++ )
                                {
                                    int ec=k1;
                                    vLocal(I1,I2,I3,nPolarization+ec)=0.;
                                    if( plotPolarizationErrors )
                                        vLocal(I1,I2,I3,nPolarizationErr+ec)=0.;
                                    for( int k2=0; k2<6; k2++ )
                                    {
                                        for( int n=0; n<Np(k1,k2); n++ )
                                        {
                                            vLocal(I1,I2,I3,nPolarization+ec) += pLocal(I1,I2,I3,pc);
                                            if( plotPolarizationErrors )
                                                vLocal(I1,I2,I3,nPolarizationErr+ec) += pErrLocal(I1,I2,I3,pc);
                                            pc+=2;   // we store p and pt so increment by 2
                                        }
                                    }
                                }
                            }
                        }
                    }
                } // end BA Maxwell 
            } // end if plotPolarization
        }

    // --- save the nonlinear variables in the augmented solution ---
        {
            if( plotNonlinearComponents )
            {
        // --- plot nonlinear components ----
                Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
                getIndex(mg.dimension(),I1,I2,I3);
                const int includeGhost=1;
                const int domain = cg.domainNumber(grid);
                const DispersiveMaterialParameters & dmp = getDomainDispersiveMaterialParameters(domain);
                if( dmp.isNonlinearMaterial() )
                {
                    realMappedGridFunction & q = getNonlinearModelMappedGridFunction( grid,current );
                    OV_GET_SERIAL_ARRAY(real,q,qLocal);
                    bool ok = ParallelUtility::getLocalArrayBounds(q,qLocal,I1,I2,I3,includeGhost);
                    if( ok )
                    {
                        const int numberOfAtomicLevels = dmp.getNumberOfAtomicLevels();
                        for( int n=0; n<numberOfAtomicLevels; n++ )
                            vLocal(I1,I2,I3,nNonlinear+n) = qLocal(I1,I2,I3,n); 
            // zero out any components that do not live on this grid 
                        for( int n=numberOfAtomicLevels; n<maxNumberOfAtomicLevels; n++ )
                            vLocal(I1,I2,I3,nNonlinear+n) = 0.;
                        if( plotNonlinearComponentErrors )
                        {
              // --- plot errors in the nonlinear variables
                            bool getErrorGridFunction=true;
                            realMappedGridFunction & qErr = getNonlinearModelMappedGridFunction( grid,current,getErrorGridFunction  );
                            OV_GET_SERIAL_ARRAY(real,qErr,qErrLocal);
                            for( int n=0; n<numberOfAtomicLevels; n++ )
                                vLocal(I1,I2,I3,nNonlinearComponentErrors+n) = qErrLocal(I1,I2,I3,n);
                            for( int n=numberOfAtomicLevels; n<maxNumberOfAtomicLevels; n++ )
                                vLocal(I1,I2,I3,nNonlinearComponentErrors+n) = 0.;
                        }
                    } // end if ok 
                }
                else
                {
          // Set Q=0 in a non-dispersive region
                    bool ok = ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3,includeGhost);
                    if( ok )
                    {
                        for( int n=0; n<maxNumberOfAtomicLevels; n++ )
                            vLocal(I1,I2,I3,nNonlinear+n) = 0.;
                        if( plotNonlinearComponentErrors )
                        {
                            for( int n=0; n<maxNumberOfAtomicLevels; n++ )
                                vLocal(I1,I2,I3,nNonlinearComponentErrors+n) = 0.;
                        }
                    }
                }
            } // end if plotNonlinearComponents
        }






        if( plotIntensity )
        {
            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
            getIndex(mg.dimension(),I1,I2,I3);
            const int includeGhost=1;
            bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3,includeGhost); 

            if( pIntensity!=NULL )
            {
                realCompositeGridFunction & intensity = *pIntensity;
                #ifdef USE_PPP
                    realSerialArray intensityLocal; getLocalArrayWithGhostBoundaries(intensity[grid],intensityLocal);
                #else
                    const realSerialArray & intensityLocal = intensity[grid];
                #endif
                if( ok )
                {
                    vLocal(I1,I2,I3,nIntensity)= intensityLocal(I1,I2,I3);
                } // end if ok 
            }
            else
            {
                vLocal(I1,I2,I3,nIntensity)=0.;
            }
        } // end if( plotIntensity )

        if( plotHarmonicElectricFieldComponents )
        {
      // plot Er and Ei assuming : E(x,t) = Er(x)*cos(w*t) + Ei(x)*sin(w*t)

            Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
            getIndex(mg.dimension(),I1,I2,I3);
            const int includeGhost=1;
            bool ok = ParallelUtility::getLocalArrayBounds(u,uLocal,I1,I2,I3,includeGhost); 

            Range Rx2=2*numberOfDimensions;
            if( pHarmonicElectricField!=NULL )
            {
                realCompositeGridFunction & hef = *pHarmonicElectricField;
                #ifdef USE_PPP
                    realSerialArray hefLocal; getLocalArrayWithGhostBoundaries(hef[grid],hefLocal);
                #else
                    const realSerialArray & hefLocal = hef[grid];
                #endif
                if( ok )
                {
                    vLocal(I1,I2,I3,Rx2+nHarmonicE)= hefLocal(I1,I2,I3,Rx2);
                } // end if ok 
            }
            else
            {
                vLocal(I1,I2,I3,Rx2+nHarmonicE)=0.;
            }
        }


        if( saveErrors )
        {
            realMappedGridFunction & err = errp!=NULL ? *errp : cgerrp!=NULL ? (*cgerrp)[grid] : u;    
#ifdef USE_PPP
            realSerialArray errLocal; getLocalArrayWithGhostBoundaries(err,errLocal);
#else
            const realSerialArray & errLocal = err;
#endif
            vLocal(all,all,all,N+numberOfComponents)=errLocal(all,all,all,N);
        }
            
        if( useVariableDissipation )
        {
#ifdef USE_PPP
            realSerialArray varDissLocal; getLocalArrayWithGhostBoundaries((*variableDissipation)[grid],varDissLocal);
#else
            const realSerialArray & varDissLocal = (*variableDissipation)[grid];
#endif
            vLocal(all,all,all,nVarDis)=varDissLocal;
        }
            
        if( saveDissipation )
        {
#ifdef USE_PPP
            realSerialArray dissLocal; getLocalArrayWithGhostBoundaries((*cgdissipation)[grid],dissLocal);
#else
            const realSerialArray & dissLocal = (*cgdissipation)[grid];
#endif
            vLocal(all,all,all,N+ndd)=dissLocal(all,all,all,N);
        }

        if( plotRho )
        {
            getChargeDensity( current,t,v,nRho );
        }


    }  // end for grid 
      
  // Moved outside grid loop ! *wdh* Dec 13, 2020
    if( method==nfdtd  || method==sosup || method==bamx )
    {
    // printF(" $$$$ plot: call getMaxDivergence $$$$\n");

        if( plotDivergence )
        {
            getMaxDivergence( current,t, &v,nDivE, &v,nRho);
        }
        else
        {
            getMaxDivergence( current,t );
        }   
    }

    if( plotDivergence && mgp==NULL ) 
    {
    // we need to interpolate the divergence to give values at the interp. pts. for plotting
        v.interpolate(Range(nDivE,nDivE));
    }
    
    return v;
    
}

