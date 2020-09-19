#================================================================================================
#
#  cgmx example:  MULTI-LEVEL ATOMIC NONLINEAR MODEL 
echo to terminal 1
#
# Usage:
#   
#  cgmx [-noplot] mla -g=<name> -tf=<tFinal> -tp=<tPlot> -kx=<num> -ky=<num> -kz=<num> ...
#                 -plotIntensity=[0|1] -diss=<>  -filter=[0|1] -debug=<num> -cons=[0/1] -varDiss=<0|1> ...
#                 -rbc=[abcEM2|rbcNonLocal|abcPML] -leftBC=[rbc|planeWave] -method=[fd|Yee|sosup] 
#                 -probeFileName=<s> -xLeftProbe=<f> -xRightProbe=<f> ...
#                 -useSosupDissipation=[0|1] -dm=[none|gdm] -ic=[slabs|pw|gp|pulse] -go=[run/halt/og]
#
# Arguments:
#  -kx= -ky= -kz= : wave number of the incident wave
#  -dm : dispersion model
#  -ic : initial condition, pw=plane-wave, gp=Gaussian-pulse
# 
#================================================================================================
# 
$tFinal=5; $tPlot=.1; $diss=1.; $filter=0; $dissOrder=-1; $cfl=.9; $varDiss=0; $varDissSmooths=20; $sidebc="symmetry"; 
$kx=1; $ky=0; $kz=0; $plotIntensity=0; $intensityOption=1; $checkErrors=0; $method="NFDTD"; $dm="none"; $ic="slabs"; 
$ax=0.; $ay=0.; $az=0.; # plane wave coeffs. all zero -> use default
$numSlabs=0; # 0 = default case of scattering from a "innerDomain" 
$x0=.5; $y0=0; $z0=0; $beta=50; # for Gaussian plane wave IC
# --- old: 
$eps0=1.; $mu0=1.; # outer domain 
$eps1=1.; $mu1=1.; # block 1 
$eps2=1.; $mu2=1.; # block 2 
$eps3=1.; $mu3=1.; # block 3 
$eps4=1.; $mu4=1.; # block 4 
# -- new: 
$matFileLeft ="baMatIsoEps1.txt";
$matFileRight="baMatIsoEps1.txt";
$matFile1="baMatIsoEps2.txt";
$matFile2="baMatIsoEps3.txt";
$matFile3="baMatIsoEps8.txt";
$matFile4="baMatIsoEps4.txt";
#
$show=" "; $backGround="leftBackGround"; $compareToShowFile=""; $flushFrequency=100; 
$interfaceEquationOption=1; $interfaceIterations=10;  $interfaceOmega=.5; $useNewInterface=1; 
$grid="slabsGrid2d2Slabse4.order4.hdf";
$cons=1; $go="halt"; 
$xa=-100.; $xb=-1.5; $ya=-100.; $yb=100.; $za=-100.; $zb=100.;  # initial condition bounding box
$leftBC="rbc"; $bcBody=""; 
$probeFileName="probeFile"; $xLeftProbe=-1.5; $xRightProbe=1.5; $yLeftProbe=0; $yRightProbe=0; $probeFrequency=1; 
$xar=-2.; $xbr=-1.; # reflection probe x-bounds
$xat= 1.; $xbt= 2.; # transmission probe x-bounds 
$rbc="abcEM2"; $pmlLines=11; $pmlPower=6; $pmlStrength=50.; 
# -- sosup dissipation --
$useSosupDissipation=0; $sosupParameter=1.;  $sosupDissipationOption=1; $sosupDissipationFrequency=1;
$selectiveDissipation=0;
# 
# $alphaP=-1.; $a0=1.; $a1=0.; $b0=0.; $b1=1.;  # GDM parameters
$dm="none"; $alphaP = (); @npv=();  $modeGDM=-1; 
@a01 = (); @a11=(); @b01=(); @b11=(); # these must be null for GetOptions to work, defaults are given below
@a02 = (); @a12=(); @b02=(); @b12=(); # for a second GDM domain 
$dmFile=""; # old way
$matFile=""; # "SilverJCDispersionFits.txt"; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"diss=f"=>\$diss,"tp=f"=>\$tPlot,"show=s"=>\$show,"debug=i"=>\$debug, \
 "cfl=f"=>\$cfl, "bg=s"=>\$backGround,"bcn=s"=>\$bcn,"go=s"=>\$go,"noplot=s"=>\$noplot,\
   "plotIntensity=i"=>\$plotIntensity,"ax=f"=>\$ax,"ay=f"=>\$ay,"az=f"=>\$az,"intensityOption=i"=>\$intensityOption,\
  "dtMax=f"=>\$dtMax,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz, "numSlabs=i"=>\$numSlabs,\
   "ii=i"=>\$interfaceIterations,"varDiss=i"=>\$varDiss ,"varDissSmooths=i"=>\$varDissSmooths,\
   "xb=f"=>\$xb,"yb=f"=>\$yb,"cons=i"=>\$cons,"compareToShowFile=s"=>\$compareToShowFile,\
   "probeFileName=s"=>\$probeFileName,"xLeftProbe=f"=>\$xLeftProbe,"xRightProbe=f"=>\$xRightProbe,\
   "yLeftProbe=f"=>\$yLeftProbe,"yRightProbe=f"=>\$yRightProbe,"flushFrequency=f"=>\$flushFrequency,\
   "checkErrors=i"=>\$checkErrors,"sidebc=s"=>\$sidebc,"dissOrder=i"=>\$dissOrder,"method=s"=>\$method,\
   "filter=i"=>\$filter, "backGround=s"=>\$backGround,"rbc=s"=>\$rbc,"pmlLines=i"=>\$pmlLines,\
   "pmlPower=i"=>\$pmlPower,"pmlStrength=f"=>\$pmlStrength,"leftBC=s"=>\$leftBC,"bcBody=s"=>\$bcBody,\
   "eps0=f"=>\$eps0,"eps1=f"=>\$eps1,"eps2=f"=>\$eps2,"eps3=f"=>\$eps3,"eps4=f"=>\$eps4,\
   "xar=f"=>\$xar,"xbr=f"=>\$xbr,"xat=f"=>\$xat,"xbt=f"=>\$xbt,"dm=s"=>\$dm,"ic=s"=>\$ic,\
   "npv=i{1,}"=>\@npv,"alphaP=f{1,}"=>\@alphaP,"a01=f{1,}"=>\@a01,"a11=f{1,}"=>\@a11,"b01=f{1,}"=>\@b01,"b11=f{1,}"=>\@b11,\
   "a02=f{1,}"=>\@a02,"a12=f{1,}"=>\@a12,"b02=f{1,}"=>\@b02,"b12=f{1,}"=>\@b12,\
  "useSosupDissipation=i"=>\$useSosupDissipation,"sosupParameter=f"=>\$sosupParameter,\
  "sosupDissipationOption=i"=>\$sosupDissipationOption,"sosupDissipationFrequency=i"=>\$sosupDissipationFrequency,\
  "selectiveDissipation=i"=>\$selectiveDissipation,"x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,\
  "dmFile=s"=>\$dmFile,"matFile1=s"=>\$matFile1,"matFile2=s"=>\$matFile2,"matFile3=s"=>\$matFile3,"matFile4=s"=>\$matFile4,\
  "matFileLeft=s"=>\$matFileLeft,"matFileRight=s"=>\$matFileRight,"probeFrequency=i"=>\$probeFrequency );
# -------------------------------------------------------------------------------------------------
if( $dm eq "none" ){ $dm="no dispersion"; }
if( $dm eq"drude" || $dm eq "Drude" ){ $dm="Drude"; }
if( $dm eq "gdm" ){ $dm="GDM"; $cons=0; } # Turn off conservative for GDM
#
if( $dmFile ne "" ){ $matFile=$dmFile; } # for backward compat
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
# 
$grid
#
$method
# ---- dispersion model ---
$dm
# 
#  --- nonlinear model ----
$nonLinearModel="multilevelAtomic";
$nonLinearModel
# pause
#
$kxa= abs($kx);
if( $kxa > 1. ){ $xb = int( $xb*$kxa +.5 )/$kxa; }  # we need to clip the plane wave on a period
if( $kx < 0 ){ $xa=$xb; $xb=100.; }
if( $leftBC eq "rbc" && $ic ne "gp" ){ $cmd="initial condition bounding box $xa $xb $ya $yb $za $zb"; }else{ $cmd="#"; }
$cmd
# initial condition bounding box $xa $xb $ya $yb $za $zb
#  damp initial conditions at face (side,axis)=(0,1) of the box
bounding box decay face 1 0 
$betaBB=5.; # exponent in tanh function for smooth transition to zero outside the bounding box
bounding box decay exponent $betaBB
# 
$epsPW=$eps0; $muPW=$mu0; # parameters for the incident plane wave
plane wave coefficients $ax $ay $az $epsPW $muPW
#
# if( $leftBC eq "rbc" ){ $cmd = "planeWaveInitialCondition"; }else{ $cmd="zeroInitialCondition"; }
# if( $ic eq "gp" ){ $cmd="Gaussian plane wave: $beta $x0 $y0 0 (beta,x0,y0,z0)\n gaussianPlaneWave"; }
if( $ic eq "gpw" ){ $cmd="gaussianPlaneWave\n Gaussian plane wave: $beta $x0 0 0 (beta,x0,y0,z0)"; }
if( $ic eq "slabs" ){ $cmd="user defined known solution\n  slabs\n  oneSlabIsoEps2\n done\n userDefinedKnownSolutionInitialCondition"; }
if( $ic eq "pulse" ){ $cmd="user defined known solution\n nonlinear pulse\n soliton\n done\n userDefinedKnownSolutionInitialCondition"; }
# Execute the initial condition commands: 
$cmd
# 
# if( $checkErrors ){ $known="planeWaveKnownSolution"; }else{ $known="#"; }
# $known
#
use new interface routines $useNewInterface
# zeroInitialCondition
# ====
# planeWaveScatteredFieldInitialCondition
# ====
# twilightZone
#  degreeSpace, degreeTime  1 1
#
kx,ky,kz $kx $ky $kz
#
# bc: all=dirichlet
# bc: all=perfectElectricalConductor
bc: all=$rbc
#
#- if( $leftBC eq "planeWave" ){ $cmd="bc: $backGround(0,0)=planeWaveBoundaryCondition"; }else{ $cmd="#"; }
#- $cmd 
#- if( $bcBody eq "pec" ){ $cmd="bc: annulus=perfectElectricalConductor"; }else{ $cmd="#"; }
#- $cmd
#
pml width,strength,power $pmlLines $pmlStrength $pmlPower
# 
# -- we need to subtract out the incident field on the "inflow" boundary before
#    applying the radiation boundary condition: 
if( $ic eq "gpw" ){ $adjustFields=1; }else{ $adjustFields=0; }
#- adjust boundaries for incident field $adjustFields all
if( $kx > 0 ){ $backGround="leftBackGround"; }else{ $backGround="rightBackGround"; }
adjust boundaries for incident field $adjustFields $backGround
#
# 
# NOTE: material interfaces have share>=100
#
#  -- Assign material parameters ---
#  For multi-block grids we assume domain names of blockDomain0, blockDomain1, etc. 
#
#  for single domains:
number of polarization vectors: $npv[1]
if( $numSlabs eq 0 ){ $cmd="GDM domain name: all\n material file: $matFile1"; }else{ $cmd="#"; }
$cmd
# 
GDM domain name: leftDomain
material file: $matFileLeft
GDM domain name: rightDomain
material file: $matFileRight
#
# GDM domain name: blockDomain0
# number of polarization vectors: $npv[1]
# material file: $matFile1
if( $numSlabs ge 1 ){ $cmd="GDM domain name: blockDomain0\n material file: $matFile1"; }else{ $cmd="#"; }
$cmd 
#
if( $numSlabs ge 2 ){ $cmd="GDM domain name: blockDomain1\n material file: $matFile2"; }else{ $cmd="#"; }
$cmd
if( $numSlabs ge 3 ){ $cmd="GDM domain name: blockDomain2\n material file: $matFile3"; }else{ $cmd="#"; }
$cmd 
if( $numSlabs ge 4 ){ $cmd="GDM domain name: blockDomain2\n material file: $matFile4"; }else{ $cmd="#"; }
$cmd 
#
interface BC iterations $interfaceIterations
# interfaceEquationsOption=0 : use extrap for 2nd ghost, 1=use eqns
interface equations option $interfaceEquationOption
omega for interface iterations $interfaceOmega
#
# bc: Annulus=perfectElectricalConductor
tFinal $tFinal
tPlot  $tPlot
#
apply filter $filter
order of dissipation $dissOrder
dissipation $diss
#
# -- sosup dissipation
#
if( $selectiveDissipation eq "1" ){ $cmd="selective dissipation...\n  turn off rectangular\n continue"; }else{ $cmd="#"; }
$cmd 
#
use sosup dissipation $useSosupDissipation
sosup parameter $sosupParameter
sosup dissipation option $sosupDissipationOption
sosup dissipation frequency $sosupDissipationFrequency
#
#*********************************
show file options...
  MXSF:compressed
  MXSF:open
    $show
  # MXSF:frequency to save 
  MXSF:frequency to flush $flushFrequency 
exit
#**********************************
#
use variable dissipation $varDiss
number of variable dissipation smooths $varDissSmooths
use conservative difference $cons
# order of dissipation 4
debug 0
#
cfl $cfl 
plot divergence 0
printf("checkErrors=$checkErrors\n");
plot errors $checkErrors
check errors $checkErrors
plot intensity $plotIntensity
intensity option $intensityOption
$intensityAveragingInterval=1.;
intensity averaging interval $intensityAveragingInterval
# $c1 = 1./sqrt($eps1*$mu1); $omega= $c1*sqrt( $kx*$kx + $ky*$ky + $kz*$kz );
# time harmonic omega $omega (omega/(2pi), normally c*|k|
#
# We can compute errors compared to another solution (e.g. to test far-field BC's)
if( $compareToShowFile ne "" ){ $cmd ="compare to show file 1\n reference show file: $compareToShowFile\n plot errors 1\n check errors 1"; }else{ $cmd="#"; }
$cmd 
#
# output probes every time step: 
probe frequency $probeFrequency
#
# Point probe: 
create a probe...
  $leftProbeFileName="left$probeFileName.dat"; 
  file name $leftProbeFileName
  probe name leftProbe
  nearest grid point to $xLeftProbe $yLeftProbe 0
exit
create a probe...
  $rightProbeFileName="right$probeFileName.dat"; 
  file name $rightProbeFileName
  probe name rightProbe
  nearest grid point to $xRightProbe $yRightProbe 0
exit
# -- surface integral probe, "average" option means scale by surface area  --
if( $probeFileName ne "" ){ $intProbeName = "integralLeft$probeFileName"; }
if( $intProbeName ne "" ){ $cmd="create a probe...\n   probe name $intProbeName\n  file name  $intProbeName.dat\n coordinate plane probe \n integral\n average\n  grid coordinate plane 0 $xLeftProbe $yLeftProbe 0 (axis, x,y,z)\n all components\n exit"; }else{ $cmd="#"; }
$cmd 
# 
if( $probeFileName ne "" ){ $intProbeName = "integralRight$probeFileName"; }
if( $intProbeName ne "" ){ $cmd="create a probe...\n   probe name $intProbeName\n  file name  $intProbeName.dat\n coordinate plane probe \n integral\n average\n grid coordinate plane 0 $xRightProbe $yRightProbe 0 (axis, x,y,z)\n all components\n exit"; }else{ $cmd="#"; }
$cmd 
# OLD: 
# probe file: probeFile.dat
# specify probes
#   -1.5 0 0 
#    1.5 .0 0
# done
#
continue
#
#
if( $az==0 ){ $cmd="plot:Ey"; }else{ $cmd="plot:Ez"; }
$cmd 
contour
  plot contour lines (toggle)
  # vertical scale factor 0.2
  # min max -1.1 1.1
  # plot a contour plane in 3d 
  if( $grid =~ /3d/ ){ $cmd="delete contour plane 2\n delete contour plane 1\n delete contour plane 0\n add contour plane  0.00000e+00  0.00000e+00  1.00000e+00 0 0 0"; }else{ $cmd="#"; }
  $cmd
exit
$go
