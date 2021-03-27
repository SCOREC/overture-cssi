#
#  ---- Include file for cgmx cmd scripts: get common options  ------
# 
#     -method=[fd|Yee|sosup|bamx] 
if( $method          eq "" ){ $method="fd"; }
if( $grid            eq "" ){ $grid="square20.order2.hdf"; }
#     -go=[run/halt/og]
if( $go              eq "" ){ $go="go"; }
# User supplied extra commands in an include file to execute at the end of the script (before -go) (e.g. for plotting)
if( $extraCommandsEnd eq "" ){ $extraCommandsEnd=""; }
#    Dispersive model: -dm=[none|gdm]
if( $dm              eq "" ){ $dm="none"; }
#    Nonlinear model:  -nm=[none|mla]
#       mla = multi-level atomic model 
if( $nm              eq "" ){ $nm="none"; }
if( $plotNonlinearComponents eq "" ){ $plotNonlinearComponents=1; }
# 
if( $tFinal          eq "" ){ $tFinal=5; }
if( $tPlot           eq "" ){ $tPlot=.1; }
if( $cfl             eq "" ){ $cfl=.9; }
#    [kx,ky,kz] : incident wave vector (scaled later by 2*pi)
if( $kx              eq "" ){ $kx=1; }
if( $ky              eq "" ){ $ky=0; }
if( $kz              eq "" ){ $kz=0; }
if( $diss            eq "" ){ $diss=1.; }
if( $dissOrder       eq "" ){ $dissOrder=-1; }
if( $filter          eq "" ){ $filter=0; }
if( $varDiss         eq "" ){ $varDiss=0; }
if( $varDissSmooths  eq "" ){$varDissSmooths=20; }
#   cons=1 : use conservative differences
if( $cons            eq "" ){ $cons=0; }
if( $debug           eq "" ){ $debug=0; }
if( $checkErrors     eq "" ){ $checkErrors=0; }
if( $maxIterationsForImplicitInterpolation eq "" ){ $maxIterationsForImplicitInterpolation=10; }
#
#  ------ initial conditions ----
#   -ic=[pw|gp|zeroInitialCondition]
#      pw : plane wave
#      gp : Gaussian plane wave  
if( $ic              eq "" ){ $ic="pw";  }
#  ------ known solution ----
#   -known=[planeWave|...]
if( $known           eq "" ){ $known="#";  }
# ----- initial bounding box -----
if( $xabb         eq "" ){ $xabb=-100.; }
if( $xbbb         eq "" ){ $xbbb= -1.5; }
if( $yabb         eq "" ){ $yabb=-100.; }
if( $ybbb         eq "" ){ $ybbb= 100.; }
if( $zabb         eq "" ){ $zabb=-100.; }
if( $zbbb         eq "" ){ $zbbb=-100.; }
#
#  ------- Plane wave parameters --------
#  [ax,ay,az] : plane wave coeffs. all zero -> use default
if( $ax              eq "" ){ $ax=0.; }
if( $ay              eq "" ){ $ay=0.;  }
if( $az              eq "" ){ $az=0.; }
#
#  ----- Gaussian plane wave parameters --------
#    [x0,y0,z0], beta, k0 : for Gaussian plane wave
if( $x0              eq "" ){ $x0=.5; }
if( $y0              eq "" ){ $y0=0; }
if( $z0              eq "" ){ $z0=0; }
if( $beta            eq "" ){ $beta=50; }
if( $k0              eq "" ){ $k0=0; }
#
# ----- show file parameters ----
if( $show              eq "" ){ $show=" "; }
if( $compareToShowFile eq "" ){ $compareToShowFile=""; }
if( $flushFrequency    eq "" ){ $flushFrequency=10; }
#
# ---- interface parameters -------
if( $interfaceEquationOption eq "" ){ $interfaceEquationOption=1; }
if( $interfaceIterations     eq "" ){ $interfaceIterations=10; }
if( $interfaceOmega          eq "" ){ $interfaceOmega=.5; }
if( $useNewInterface         eq "" ){ $useNewInterface=1; }
if( $tallCellRatioBound      eq "" ){ $tallCellRatioBound=1.25; }
#
# ---- plotting options ------
if( $plotPolarizationComponents eq "" ){ $plotPolarizationComponents=1; }
if( $plotIntensity              eq "" ){ $plotIntensity=0; }
if( $intensityOption            eq "" ){ $intensityOption=1;  }
#
# ---- Upwind (sosup) dissipation parameters ----
if( $useSosupDissipation eq ""       ){ $useSosupDissipation=0; }
if( $sosupParameter eq ""            ){ $sosupParameter=1.; }
if( $sosupDissipationOption eq ""    ){ $sosupDissipationOption=1; }
if( $sosupDissipationFrequency eq "" ){ $sosupDissipationFrequency=1; }
if( $selectiveDissipation eq ""      ){ $selectiveDissipation=0; }
#
#  ----- boundary conditions ------
#
# bcCmds = array of user supplied BC commands
@bcCmds = ();
#
# --- radiation boundary conditions ---
#  -rbc=[abcEM2|rbcNonLocal|abcPML] 
if( $rbc          eq "" ){ $rbc="abcEM2"; }
if( $pmlLines     eq "" ){ $pmlLines=11; }
if( $pmlPower     eq "" ){ $pmlPower=6; }
if( $pmlStrength  eq "" ){ $pmlStrength=50.; }
#
#  ----- probes ---------
if( $probeFrequency    eq "" ){ $probeFrequency=1; }
if( $probeFileName     eq "" ){ $probeFileName=""; }
if( $xLeftProbe        eq "" ){ $xLeftProbe=-1.5; }
if( $xRightProbe       eq "" ){ $xRightProbe=1.5; }
if( $yLeftProbe        eq "" ){ $yLeftProbe=0; }
if( $yRightProbe       eq "" ){ $yRightProbe=0; }
#
# ------ material properities -------
#  File conatining info on dispersive and non-dispersive materials: 
if( $materialFile  eq "" ){ $materialFile="diskMaterials.h"; }
# optionally supply material files on the command line to over-ride those in $materialFile : 
@matFileArray = ();         
# number of bodies to set materials for 
if( $numBodies         eq "" ){ $numBodies=1; }
# 
# Arrays to hold GDM parameters:  (old way)
#
$alphaP = ();
@npv=();
if( $modeGDM eq "" ){ $modeGDM=-1; }
@a01 = (); @a11=(); @b01=(); @b11=(); # these must be null for GetOptions to work, defaults are given below
@a02 = (); @a12=(); @b02=(); @b12=(); # for a second GDM domain 
# 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"diss=f"=>\$diss,"tp=f"=>\$tPlot,"show=s"=>\$show,"debug=i"=>\$debug, \
 "cfl=f"=>\$cfl, "bg=s"=>\$backGround,"bcn=s"=>\$bcn,"go=s"=>\$go,"noplot=s"=>\$noplot,\
  "plotIntensity=i"=>\$plotIntensity,"ax=f"=>\$ax,"ay=f"=>\$ay,"az=f"=>\$az,"intensityOption=i"=>\$intensityOption,\
  "dtMax=f"=>\$dtMax,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz, "numBlocks=i"=>\$numBlocks,\
  "ii=i"=>\$interfaceIterations,"varDiss=i"=>\$varDiss ,"varDissSmooths=i"=>\$varDissSmooths,\
  "xb=f"=>\$xb,"yb=f"=>\$yb,"cons=i"=>\$cons,"compareToShowFile=s"=>\$compareToShowFile,\
  "probeFileName=s"=>\$probeFileName,"xLeftProbe=f"=>\$xLeftProbe,"xRightProbe=f"=>\$xRightProbe,\
  "yLeftProbe=f"=>\$yLeftProbe,"yRightProbe=f"=>\$yRightProbe,"flushFrequency=f"=>\$flushFrequency,\
  "checkErrors=i"=>\$checkErrors,"sidebc=s"=>\$sidebc,"dissOrder=i"=>\$dissOrder,"method=s"=>\$method,\
  "filter=i"=>\$filter, "backGround=s"=>\$backGround,"rbc=s"=>\$rbc,"pmlLines=i"=>\$pmlLines,\
  "pmlPower=i"=>\$pmlPower,"pmlStrength=f"=>\$pmlStrength,"leftBC=s"=>\$leftBC,"bcBody=s"=>\$bcBody,\
  "eps0=f"=>\$eps0,"eps1=f"=>\$eps1,"eps2=f"=>\$eps2,"eps3=f"=>\$eps3,"eps4=f"=>\$eps4,\
  "xabb=f"=>\$xabb,"xbbb=f"=>\$xbbb,"yabb=f"=>\$yabb,"ybbb=f"=>\$ybbb,"zabb=f"=>\$zabb,"zbbb=f"=>\$zbbb, \
  "xar=f"=>\$xar,"xbr=f"=>\$xbr,"xat=f"=>\$xat,"xbt=f"=>\$xbt,"dm=s"=>\$dm,"ic=s"=>\$ic,\
  "npv=i{1,}"=>\@npv,"alphaP=f{1,}"=>\@alphaP,"a01=f{1,}"=>\@a01,"a11=f{1,}"=>\@a11,"b01=f{1,}"=>\@b01,"b11=f{1,}"=>\@b11,\
  "a02=f{1,}"=>\@a02,"a12=f{1,}"=>\@a12,"b02=f{1,}"=>\@b02,"b12=f{1,}"=>\@b12,\
  "useSosupDissipation=i"=>\$useSosupDissipation,"sosupParameter=f"=>\$sosupParameter,\
  "sosupDissipationOption=i"=>\$sosupDissipationOption,"sosupDissipationFrequency=i"=>\$sosupDissipationFrequency,\
  "selectiveDissipation=i"=>\$selectiveDissipation,"x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"k0=f"=>\$k0,\
  "probeFrequency=i"=>\$probeFrequency,\
  "materialFile=s"=>\$materialFile,"numBodies=i"=>\$numBodies,"matFileArray=s{1,}"=>\@matFileArray,\
  "bcCmds=s{1,}"=>\@bcCmds,"nm=s"=>\$nm,"plotNonlinearComponents=f"=>\$plotNonlinearComponents,\
  "tallCellRatioBound=f"=>\$tallCellRatioBound,"extraCommandsEnd=s"=>\$extraCommandsEnd, \
  "maxIterationsForImplicitInterpolation=i"=>\$maxIterationsForImplicitInterpolation );
# -------------------------------------------------------------------------------------------------
#
if( $dm eq "none" ){ $dm="no dispersion"; }
if( $dm eq "gdm" ){ $dm="GDM"; $cons=0; } # Turn off conservative for GDM
if( $nm eq "none" ){ $nm="#"; }
if( $nm eq "mla" ){ $nm="multilevelAtomic"; }
#
if( $method eq "sosup" ){ $diss=0.; }
if( $method eq "fd" ){ $method="nfdtd"; }
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
#
# Give defaults here for array arguments: 
if( $alphaP[0] eq "" ){ @alphaP=(-1,-1); } # default -1 means use 1/eps
if( $interfaceNormal[0] eq "" ){ @interfaceNormal=(1,0,0); }
if( $interfacePoint[0] eq "" ){ @interfacePoint=(0,0,0); }
if( $npv[0] eq "" ){ @npv=(0,0); }
if( $a01[0] eq "" ){ @a01=(1,0,0,0); }
if( $a11[0] eq "" ){ @a11=(0,0,0,0); }
if( $b01[0] eq "" ){ @b01=(0,0,0,0); }
if( $b11[0] eq "" ){ @b11=(0,0,0,0); }
#
if( $a02[0] eq "" ){ @a02=(1,0,0,0); }
if( $a12[0] eq "" ){ @a12=(0,0,0,0); }
if( $b02[0] eq "" ){ @b02=(0,0,0,0); }
if( $b12[0] eq "" ){ @b12=(0,0,0,0); }
