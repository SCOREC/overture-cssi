// ========================================================================================================
/// \class Cgad
/// \brief Solver for advection-diffusion (AD) equations.
// ========================================================================================================

#include "Cgad.h"
#include "AdParameters.h"
#include "Ogshow.h"

Cgad::
Cgad(CompositeGrid & cg_, 
      GenericGraphicsInterface *ps /* =NULL */, 
      Ogshow *show /* =NULL */ , 
      const int & plotOption_ /* =1 */) 
   : DomainSolver(*(new AdParameters),cg_,ps,show,plotOption_)
// ===================================================================================================
// Notes:
//   AdParameters (passed to the DomainSolver constructor above) replaces the base class Parameters
// ===================================================================================================
{
  className="Cgad";
  name="ad";

  // should this be somewhere else? setup?
  if( realPartOfEigenvalue.size() != cg.numberOfComponentGrids() )
    realPartOfEigenvalue.resize(cg.numberOfComponentGrids(),-1.);
  if( imaginaryPartOfEigenvalue.size() != cg.numberOfComponentGrids() )
    imaginaryPartOfEigenvalue.resize(cg.numberOfComponentGrids(),-1.);


  // ------- Set the default order of accuracy from the grid parameters  ------
  //  Added April 16, 2021 *wdh*
  
  int minDiscretizationWidth=INT_MAX;
  int minInterpolationWidth=INT_MAX;
  Range R=cg.numberOfDimensions();
  const IntegerArray & iw = cg.interpolationWidth;
  // iw.display("iw");
  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    MappedGrid & mg = cg[grid];
    const IntegerArray & dw = mg.discretizationWidth();
      
    // dw.display("dw");
      
    minDiscretizationWidth=min(minDiscretizationWidth,min(dw(R)));
      
    for( int grid2=0; grid2<cg.numberOfComponentGrids(); grid2++ )
    {
      if( grid!=grid2 )
        minInterpolationWidth=min( minInterpolationWidth,min(iw(R,grid,grid2)));
    }
  }
  if( minInterpolationWidth==INT_MAX ) minInterpolationWidth=minDiscretizationWidth;
  printF("Cgad::constructor: minDiscretizationWidth=%i, minInterpolationWidth=%i.\n",minDiscretizationWidth,
         minInterpolationWidth);

  const int maxOrderOfAccuracy=8;  // *************
    
  int & orderOfAccuracyInSpace = parameters.dbase.get<int>("orderOfAccuracy");
  int & orderOfAccuracyInTime  = parameters.dbase.get<int>("orderOfTimeAccuracy");
  

  orderOfAccuracyInSpace=min(maxOrderOfAccuracy,minDiscretizationWidth-1,minInterpolationWidth-1);
  if( orderOfAccuracyInSpace%2 ==1 )
    orderOfAccuracyInSpace--;   // must be even
 

  printP("Cgad::constructor: INFO: setting default order of accuracy=%i based on the input grid parameters\n",
         orderOfAccuracyInSpace);
  
}



Cgad::
~Cgad()
{
  delete & parameters;
}




int Cgad::
updateToMatchGrid(CompositeGrid & cg)
{
  printF("\n $$$$$$$$$$$$$$$ Cgad: updateToMatchGrid(CompositeGrid & cg) $$$$$$$$$$$$\n\n");
  
  int returnValue =DomainSolver::updateToMatchGrid(cg);

  if( realPartOfEigenvalue.size() != cg.numberOfComponentGrids() )
    realPartOfEigenvalue.resize(cg.numberOfComponentGrids(),-1.);
  if( imaginaryPartOfEigenvalue.size() != cg.numberOfComponentGrids() )
    imaginaryPartOfEigenvalue.resize(cg.numberOfComponentGrids(),-1.);

  return returnValue;
  
}

// int Cgad::
// formImplicitTimeSteppingMatrix(realMappedGridFunction & coeff,
//                             const real & dt0, 
//                             int scalarSystem, 
//                             realMappedGridFunction & uL,
//                             const int & grid )
// {
//   Overture::abort("Cgad::formImplicitTimeSteppingMatrix:ERROR: not implemented");
//   return 0;
// }


int Cgad::
updateGeometryArrays(GridFunction & cgf)
{
  if( debug() & 4 ) printF(" --- Cgad::updateGeometryArrays ---\n");

  real cpu0=getCPU();

  int grid;
  for( grid=0; grid<cgf.cg.numberOfComponentGrids(); grid++ )
  {
    if( !cgf.cg[grid].isRectangular() || twilightZoneFlow() ||  parameters.gridIsMoving(grid) )
      cgf.cg[grid].update(MappedGrid::THEcenter | MappedGrid::THEvertex );  
  }
  parameters.dbase.get<RealArray>("timing")(parameters.dbase.get<int>("timeForUpdatePressureEquation"))+=getCPU()-cpu0;

  if( realPartOfEigenvalue.size() != cgf.cg.numberOfComponentGrids() )
    realPartOfEigenvalue.resize(cgf.cg.numberOfComponentGrids(),-1.);
  
  if( imaginaryPartOfEigenvalue.size() != cgf.cg.numberOfComponentGrids() )
    imaginaryPartOfEigenvalue.resize(cgf.cg.numberOfComponentGrids(),-1.);

  return DomainSolver::updateGeometryArrays(cgf);
}
// ===================================================================================================================
/// \brief provide titles for show file output
///
// ===================================================================================================================
void Cgad::
saveShowFileComments( Ogshow &show )
{
    // save comments that go at the top of each plot
  char buffer[80]; 
  // save comments that go at the top of each plot
  aString timeLine="";
  if(  parameters.dbase.has_key("timeLine") )
    timeLine=parameters.dbase.get<aString>("timeLine");

  std::vector<real> & kappa = parameters.dbase.get<std::vector<real> >("kappa");
  std::vector<real> & a = parameters.dbase.get<std::vector<real> >("a");
  std::vector<real> & b = parameters.dbase.get<std::vector<real> >("b");
  std::vector<real> & c = parameters.dbase.get<std::vector<real> >("c");   
  
  aString showFileTitle[5];
  if( pdeName =="thinFilmEquations" )
  {
     showFileTitle[0]=sPrintF(buffer,"Thin-Film");
  }
  else
  {
    if( parameters.dbase.get<int>("numberOfDimensions")==2 )
      showFileTitle[0]=sPrintF(buffer,"Convection Diffusion, a=%g, b=%g, kappa=%g",a[0],b[0],kappa[0]);
    else
      showFileTitle[0]=sPrintF(buffer,"Convection Diffusion, a=%g, b=%g, c=%g, kappa=%g",
                               a[0],b[0],c[0],kappa[0]);
  }
  
  showFileTitle[1]=timeLine;
  showFileTitle[2]="";  // marks end of titles
  
  for( int i=0; showFileTitle[i]!=""; i++ )
    show.saveComment(i,showFileTitle[i]);
}


// ===================================================================================================================
/// \brief Output run-time parameters for the header.
/// \param file (input) : write values to this file.
///
// ===================================================================================================================
void 
Cgad::
writeParameterSummary( FILE * file )
{
  DomainSolver::writeParameterSummary( file );

  const int & numberOfComponents = parameters.dbase.get<int>("numberOfComponents");
  std::vector<real> & kappa = parameters.dbase.get<std::vector<real> >("kappa");
  std::vector<real> & a = parameters.dbase.get<std::vector<real> >("a");
  std::vector<real> & b = parameters.dbase.get<std::vector<real> >("b");
  std::vector<real> & c = parameters.dbase.get<std::vector<real> >("c");

  if ( file==parameters.dbase.get<FILE* >("checkFile") )
  {
    // -- check file header ---
    if( pdeName=="advection diffusion" || pdeName=="convection diffusion" )
    {
      fPrintF(file,"\\caption{advection-diffusion, gridName, $\\kappa=%3.2g$, $t=%2.1f$, ",
              kappa[0],parameters.dbase.get<real >("tFinal"));
    }
    else
    {
      fPrintF(file,"\\caption{%s, gridName, $t=%2.1f$, ",(const char*)pdeName,parameters.dbase.get<real >("tFinal"));
    }
    
     return;
  }


  real & thermalConductivity = parameters.dbase.get<real>("thermalConductivity");

  fPrintF(file," numberOfComponents=%i:\n",numberOfComponents);
  aString *componentName = parameters.dbase.get<aString* >("componentName");
  for( int m=0; m<numberOfComponents; m++ )
    fPrintF(file,"   component %i: %s\n",m,(const char*)componentName[m]);

  if( pdeName=="advection diffusion" || pdeName=="convection diffusion" )
  {
    if( parameters.dbase.get<bool >("variableDiffusivity") )
    {
      fPrintF(file," The coefficients of diffusivity are variable.\n");
    }
    else
    {
      fPrintF(file," The coefficients of diffusivity are constant:\n  ");
      aString name = "kappa";
      for( int m=0; m<numberOfComponents; m++ )
      {
        if( numberOfComponents==1 )
          fPrintF(file," %s=%g",(const char*)name,kappa[m]);
        else
          fPrintF(file," %s[%i]=%g,",(const char*)name,m,kappa[m]);
      }
      fPrintF(file,"\n");

    
    }
    if( parameters.dbase.get<bool >("variableAdvection") )
    {
      fPrintF(file," The advection coefficients are variable.\n");
    }
    else
    {
      fPrintF(file," The advection coefficients are constant:\n");
      for( int n=1; n<4; n++ )
      {
        std::vector<real> & par = n==1 ? a : n==2 ? b : c;
        aString name = n==1 ? "a" : n==2 ? "b" : "c";
        for( int m=0; m<numberOfComponents; m++ )
        {
          if( numberOfComponents==1 )
            fPrintF(file,"   %s=%g",(const char*)name,par[m]);
          else
            fPrintF(file,"   %s[%i]=%g,",(const char*)name,m,par[m]);
        }
        fPrintF(file,"\n");
      }

      fPrintF(file," applyChampInterfaceConditions=%d (for multi-domain problems).\n",parameters.dbase.get<int>("applyChampInterfaceConditions"));
      fPrintF(file," champOption = %d (0=use old, 1=use new implementation).\n",parameters.dbase.get<int>("champOption"));
    }

  }
  else if( pdeName=="thinFilmEquations" )
  {
    const real & S  = parameters.dbase.get<real>("inverseCapillaryNumber");
    const real & G  = parameters.dbase.get<real>("scaledStokesNumber");
    const real & h0 = parameters.dbase.get<real>("thinFilmBoundaryThickness");
    const real & he = parameters.dbase.get<real>("thinFilmLidThickness");

    fPrintF(file," thin film parameters: S=%.2g G=%.2g, h0=%.2g, he=%.2g .\n",S,G,h0,he);
  }
  
  const bool & implicitAdvection = parameters.dbase.get<bool >("implicitAdvection");
  if( implicitAdvection )
    fPrintF(file," Treat advection terms implicitly (when using implicit time-stepping).\n");
  else
    fPrintF(file," Treat advection terms explicitly (when using implicit time-stepping).\n");
  
  fPrintF(file," thermalConductivity=%g\n",thermalConductivity);

}

