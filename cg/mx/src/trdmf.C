// ====================================================================================
///  \file trdmf.C
///  \brief test program to read a Dispersive Material File
// ===================================================================================


#include "Overture.h"
#include "PlotStuff.h"
#include "DispersiveMaterialParameters.h"
#include "display.h"



int 
main(int argc, char *argv[]) 
{
  Overture::start(argc,argv);  // initialize Overture and A++/P++

  printF("Usage: trdmf -file=<data file> -matType=[isotropic|bianisotropic] -noplot -debug=<> ... \n" );

  aString fileName="numMaterial.txt";
  aString materialType="isotropic";
  
  int plotOption=true;
  
  const real nm = 1e-9;         // nanometers  (meter-per-nm)
  const real um = 1e-6;         // micrometers
  const real c0 = 299792458;    // the speed of light, [m/c]
  const real L0 = 100*nm;       // length scale -- appropriate for visible light 

  const real velocityScale = c0;  // velocity scale 
  const real lengthScale   = L0;  // length scale 
  real V0=velocityScale;

  
  char buff[180];
  int len=0;
  if( argc > 1 )
  { 
    aString line;
    for( int i=1; i<argc; i++ )
    {
      line=argv[i];
      if( line=="-noplot" || line=="noplot" )
        plotOption=false;
      else if( line=="-plot" )
        plotOption=true;
      else if( len=line.matches("-file=") )
      {
        fileName=line(len,line.length()-1);
        printF("trdmf: reading dispersive material data from file [%s]\n",(const char*)fileName);
      }
      else if( len=line.matches("-matType=") )
      {
        materialType=line(len,line.length()-1);
        printF("trdmf: materialType=[%s]\n",(const char*)materialType);
      }
      else
      {
        printF("Unknown command=[%s]\n",(const char*)line);
	OV_ABORT("error");
      }
      
    }
  }

  DispersiveMaterialParameters dmp;


  if( true )
  {
    if( materialType == "isotropic" )
      dmp.setMaterialType( DispersiveMaterialParameters::isotropic );
    else
      dmp.setMaterialType( DispersiveMaterialParameters::bianisotropic );
    dmp.setScales( velocityScale,lengthScale );
  }
  else
  {
    real V0=1., L0=1.;

    if( materialType == "isotropic" )
    {
      dmp.setMaterialType( DispersiveMaterialParameters::isotropic );

      V0=c0; // velocity scale
      // L0=200 * 1e-9;  // length scale 
      L0=100 * 1e-9;  // length scale 
      dmp.setScales( V0,L0 );
    }
    else
    {
      dmp.setMaterialType( DispersiveMaterialParameters::bianisotropic );
      V0=1.;  // velocity scale
      L0=1.;  // length scale 
      dmp.setScales( V0,L0 );
    }
  }
  
  int npRequested=-1; // -1 means use largest Np avaiable
  dmp.readFromFile(fileName, npRequested);


  if( materialType == "bianisotropic" )
  {
    RealArray K0, bianisotropicParameters;
    IntegerArray Np;
    dmp.getBianisotropicParameters( K0, bianisotropicParameters, Np );
    
    printF("TEST READ MATERIAL PARAMETERS FILE\n");
    printF("BIANISOTROPIC PARAMETERS:\n");
    ::display(K0,"K0");
    ::display(Np,"Np");
    if( false )
    {
      ::display(bianisotropicParameters,"bianisotropicParameters");
    }
    
    

  }


  GenericGraphicsInterface & gi = *Overture::getGraphicsInterface("trdmf",plotOption,argc,argv);


  // ------- Plot some quantities --------

  dmp.update( gi);


  if( true )
  {
    
    const real PHz = 1e15;
    const real nm = 1.e-9;

    int Np;
    real epsInf;
    RealArray modelParams;
  
    dmp.getIsotropicParameters( Np,  epsInf, modelParams );
    printF("Np=%d, epsInf=%g\n",Np,epsInf);
    for( int j=0; j<Np; j++ )
      printF("j=%d: [a0,a1,b0,b1]=[%g,%g,%g,%g]\n",j,modelParams(0,j),modelParams(1,j),modelParams(2,j),modelParams(3,j));



    int nw=101;
    real wMin=.15, wMax=1.1;
    RealArray w(nw), epsHatr(nw), epsHati(nw), nHatr(nw), nHati(nw), om(nw), lam(nw);
  
    real OmegaScale = 9.4182578365442600e+15; // From Gold approx 10 PHz  *** FIX ME ***
    real Omega0 = 2.*Pi*V0/L0;
    printF(" L0=%9.2e (nm), Omega0 = 2.*Pi*V0/L = %9.3e, OmegaScale=%9.3e (from file)\n",L0/nm, Omega0,OmegaScale);

    for( int i=0; i<nw; i++ )
    {
      w(i) = wMin + (wMax-wMin)*(i)/real(nw-1);

      real omega = w(i)*OmegaScale;          // omega in MKS units (Hz)
      real lambda = (2.*Pi*V0/(omega)) / nm; // lambda in nm 

      w(i) *= OmegaScale/Omega0;  // non-dimensional omega

      dmp.evalEpsAndN( w(i), epsHatr(i), epsHati(i), nHatr(i), nHati(i) );

    

      om(i)=omega/PHz;
      lam(i)=lambda;
      // real lambda = (c0/w(i)) * L0;
    
      printF(" w=%7.3g (= %7.3g PHz) lambda=%9.2e (nm) epsHat=(%9.2e,%9.2e) n=(%9.2e,%9.2e)\n",
             w(i),omega/PHz,lambda,epsHatr(i), epsHati(i), nHatr(i), nHati(i) );
    }

    RealArray fields(nw,4);
    Range all;
    fields(all,0) = epsHatr;
    fields(all,1) = epsHati;
    fields(all,2) = nHatr;
    fields(all,3) = nHati;


    // x-axis: 0=omega, 1=lambda
    int xAxisType=0;

    PlotStuffParameters psp;

    GUIState dialog;

    dialog.setWindowTitle("Read Material File tester");
    dialog.setExitCommand("exit", "exit");

    aString cmds[] = {"plot",
                      ""};

    int numberOfPushButtons=0;  // number of entries in cmds
    while( cmds[numberOfPushButtons]!="" ){numberOfPushButtons++;}; // 
    int numRows=(numberOfPushButtons+1)/2;
    dialog.setPushButtons( cmds, cmds, numRows ); 

    aString opCommand1[] = {"omega",
                            "lambda",
                            ""};

    dialog.setOptionMenuColumns(1);
    dialog.addOptionMenu( "x-axis:", opCommand1, opCommand1, xAxisType );

    gi.pushGUI(dialog);

  
    aString answer;
  
    for(;;)
    {
    
      gi.getAnswer(answer,"");  

      if( answer=="exit" )
      {
        break;
      }
      else if( answer=="omega" )
      {
        xAxisType=0;
      }
      else if( answer=="lambda" )
      {
        xAxisType=1;
      }

      else if( answer=="plot" )
      {


        psp.set(GI_TOP_LABEL,sPrintF("Test Read Material File"));
        const aString names[]={"epsr","epsi","nr","ni"};

        if( xAxisType==0 )
          PlotIt::plot(gi, om, fields, "Dispersive Material","omega (PHz)", names,psp );
        else
          PlotIt::plot(gi, lam, fields, "Dispersive Material","lambda (nm)", names,psp );

      }
    
    }
  

    gi.popGUI(); // restore the previous GUI  
  }
  
  Overture::finish(); 
  return 0;
}
