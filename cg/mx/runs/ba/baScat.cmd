#================================================================================================
#
#  cgmx example:  test the BA Maxwell solver
echo to terminal 1
#
# Usage:
#   
#  cgmx [-noplot] baScat -g=<name> -tf=<tFinal> -tp=<tPlot> -kx=<num> -ky=<num> -kz=<num> ...
#                 -plotIntensity=[0|1] -diss=<>  -filter=[0|1] -debug=<num> -cons=[0/1] -varDiss=<0|1> ...
#                 -rbc=[abcEM2|rbcNonLocal|abcPML|absorbing] -leftBC=[rbc|planeWave] -method=[fd|Yee|sosup] 
#                 -probeFileName=<s> -xLeftProbe=<f> -xRightProbe=<f> ...
#                 -useSosupDissipation=[0|1] -dm=[none|gdm] -ic=[pw|gp] -go=[run/halt/og]
#
# Arguments:
#  -kx= -ky= -kz= : integer wave number of the incident wave
#  -varDiss :  if 1, use variable dissipation (only add dissipation near interpolation points)
#  -dm : dispersion model
#  -ic : initial condition, pw=plane-wave, gp=Gaussian-pulse
# 
# Notes: 
# refractive index from vacuum to glass is 1.5
#       $epsGlass = 1.5^2 = 2.25
# 
#   pts/wave-length = (1/ds)*(1/ky) = factor*100/ky 
# 
# NOTE: diss=.5 was too small for fine grid runs;  diss=5. works but a smaller value may also be ok. 
# 
#     glass  n=1.5 -> eps1= (1.5)^2 
#     silicon  n=3.45 -> eps1= (3.45)^2 =11.9025
# Examples:
#
#================================================================================================
# 
$tFinal=5; $tPlot=.1; $diss=1.; $filter=0; $dissOrder=-1; $cfl=.9; $varDiss=0; $varDissSmooths=20; $sidebc="symmetry"; 
$kx=1; $ky=0; $kz=0; $plotIntensity=0; $intensityOption=1; $checkErrors=0; $method="NFDTD"; $dm="none"; $ic="pw"; 
$ax=0.; $ay=0.; $az=0.; # plane wave coeffs. all zero -> use default
$numBlocks=0; # 0 = default case of scattering from a "innerDomain" 
$x0=.5; $y0=0; $z0=0; $beta=50; # for Gaussian plane wave IC
$omega=5.;  $amp=1.; $rampTime=-1.;  # Gaussian souce
$xc=.5; $yc=.5; $zc=0.; $radius=.25; # for cylinder material region
$ae=.25; $be=.25; $ce=.25;  # for sphere or ellipsoid material region
$eps0=1.; $mu0=1.; # outer domain 
$eps1=1.; $mu1=1.; # block 1 
$eps2=1.; $mu2=1.; # block 2 
$eps3=1.; $mu3=1.; # block 3 
$eps4=1.; $mu4=1.; # block 4 
$ts="bamx"; $matFile="";  $numMatRegions=1; $matFile2="";  $matFile3=""; $matFile4="";
$solveForAllFields=0; $regionFile="boxRegion.h"; 
$show=" "; $backGround="backGround"; $compareToShowFile=""; $flushFrequency=10; 
$interfaceEquationOption=1; $interfaceIterations=10;  $interfaceOmega=.5; $useNewInterface=1; 
$grid="square32.order2";
$gridCmd="rectangeArg.cmd"; # comannd file for building a grid inline with ogen
$cons=1; $go="halt"; 
$bbxa=-100.; $bbxb=-1.5; $bbya=-100.; $bbyb=100.; $zabb=-100.; $zbbb=100.;  # initial condition bounding box
$leftBC="rbc"; $bcBody=""; $bc1=""; $bc2=""; $bc3=""; $bc4=""; $bc5=""; $bc6=""; 
$probeFileName="";
$xLeftProbe=-1.5; $xRightProbe=1.5; $yLeftProbe=0; $yRightProbe=0; $probeFrequency=1; 
$intProbeName="";  # integral probe
$xar=-2.; $xbr=-1.; # reflection probe x-bounds
$xat= 1.; $xbt= 2.; # transmission probe x-bounds 
$rbc="abcEM2"; $pmlLines=11; $pmlPower=6; $pmlStrength=50.; 
# -- sosup dissipation --
$useSosupDissipation=0; $sosupParameter=1.;  $sosupDissipationOption=1; $sosupDissipationFrequency=1;
$selectiveDissipation=0;
$useSuperGrid=0; $superGridWidth=.2; 
# 
# $alphaP=-1.; $a0=1.; $a1=0.; $b0=0.; $b1=1.;  # GDM parameters
$dm="none"; $alphaP = (); @npv=();  $modeGDM=-1; 
@a01 = (); @a11=(); @b01=(); @b11=(); # these must be null for GetOptions to work, defaults are given below
@a02 = (); @a12=(); @b02=(); @b12=(); # for a second GDM domain 
$dmFile=""; # "SilverJCDispersionFits.txt"; 
# -- new: boxRegionArray.h defines regions
$nxBox=2; $nyBox=2; $nzBox=1;
$xWidth=.25; $yWidth=.25; $zWidth=.25; # width of boxes 
$xSep=.25; $ySep=.25; $zSep=.25;       # separation between the boxes
$xc=.125; $yc=.125; $zc=.125;  # lower left corner of the boxes
@matFileArray = (); # holds list of materials for boxes
$materialFile = ""; # defines materials in a file *new way* March 1, 2021
# ----------------------------- get command line arguments ---------------------------------------
#
#  NOTE: To pass -ya and -yb to ogen cmd files do NOT include -yabb and -ybbb in this list (?!)
GetOptions( "g=s"=>\$grid,"gridCmd=s"=>\$gridCmd,"tf=f"=>\$tFinal,"diss=f"=>\$diss,"tp=f"=>\$tPlot, \
 "show=s"=>\$show,"debug=i"=>\$debug,"matFileArray=s{1,}"=>\@matFileArray, \
 "cfl=f"=>\$cfl, "bg=s"=>\$backGround,"bcn=s"=>\$bcn,"go=s"=>\$go,"noplot=s"=>\$noplot,\
   "plotIntensity=i"=>\$plotIntensity,"ax=f"=>\$ax,"ay=f"=>\$ay,"az=f"=>\$az,"intensityOption=i"=>\$intensityOption,\
  "dtMax=f"=>\$dtMax,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz, "numBlocks=i"=>\$numBlocks,\
   "ii=i"=>\$interfaceIterations,"varDiss=i"=>\$varDiss ,"varDissSmooths=i"=>\$varDissSmooths,\
   "cons=i"=>\$cons,"compareToShowFile=s"=>\$compareToShowFile,\
   "bbxa=f"=>\$bbxa,"bbxb=f"=>\$bbxb,"bbya=f"=>\$bbya,"bbyb=f"=>\$bbyb,"zabb=f"=>\$zabb,"zbbb=f"=>\$zbbb, \
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
  "dmFile=s"=>\$dmFile,"probeFrequency=i"=>\$probeFrequency,"ts=s"=>\$ts,\
  "matFile=s"=>\$matFile,"matFile2=s"=>\$matFile2,"matFile3=s"=>\$matFile3,"matFile4=s"=>\$matFile4,\
  "numMatRegions=i"=>\$numMatRegions,"bc1=s"=>\$bc1,"bc2=s"=>\$bc2,"bc3=s"=>\$bc3,"bc4=s"=>\$bc4,"bc5=s"=>\$bc5,"bc6=s"=>\$bc6,\
  "solveForAllFields=i"=>\$solveForAllFields,"regionFile=s"=>\$regionFile,"xc=f"=>\$xc,"yc=f"=>\$yc,"zc=f"=>\$zc,\
  "radius=f"=>\$radius,"ae=f"=>\$ae,"be=f"=>\$be,"ce=f"=>\$ce,"useSuperGrid=i"=>\$useSuperGrid,\
  "superGridWidth=f"=>\$superGridWidth,"intProbeName=s"=>\$intProbeName,"amp=f"=>\$amp,"omega=f"=>\$omega,\
  "rampTime=f"=>\$rampTime,"nxBox=i"=>\$nxBox,"nyBox=i"=>\$nyBox,"nzBox=i"=>\$nzBox,\
  "xWidth=f"=>\$xWidth,"yWidth=f"=>\$yWidth,"zWidth=f"=>\$zWidth,"xc=f"=>\$xc,"yc=f"=>\$yc,"zc=f"=>\$zc,\
  "xSep=f"=>\$xSep,"ySep=f"=>\$ySep,"zSep=f"=>\$zSep,"materialFile=s"=>\$materialFile );
# -------------------------------------------------------------------------------------------------
# GetOptions( "yb=f"=>\$yb ); 
# printf(" baScat: ya=$ya yb=$yb\n");
# pause
#printf(" num mats = $#matFileArray+1\n");
#printf(" matFileArray[0]=$matFileArray[0]\n");
#printf(" matFileArray[1]=$matFileArray[1]\n");
#pause
#
if( $ts eq "me" ){ $ts="modifiedEquationTimeStepping"; }
if( $rbc eq "c" || $rbc eq "char" ){ $rbc="characteristic"; }
if( $rbc eq "a" ){ $rbc="absorbing"; }
$orderOfRungeKutta=4;
if( $ts eq "rk1" ){ $ts="rungeKutta"; $orderOfRungeKutta=1; }
if( $ts eq "rk2" ){ $ts="rungeKutta"; $orderOfRungeKutta=2; }
if( $ts eq "rk3" ){ $ts="rungeKutta"; $orderOfRungeKutta=3; }
if( $ts eq "rk4" ){ $ts="rungeKutta"; $orderOfRungeKutta=4; }
# printf(" TS = [$ts], order=[$orderOfRungeKutta] \n");
#
if( $dm eq "none" ){ $dm="no dispersion"; }
if( $dm eq"drude" || $dm eq "Drude" ){ $dm="Drude"; }
if( $dm eq "gdm" ){ $dm="GDM"; $cons=0; } # Turn off conservative for GDM
#
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
if( $npv[0] eq "" ){ @npv=(1,0); }
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
# If $grid="ogen" then we build a grid on the fly
$grid = ($grid eq "ogen") ? "ogen\n read command file\n $gridCmd" : $grid;
# 
$grid
# pause
#
$method
# time-stepping method
$ts
order of Runge Kutta $orderOfRungeKutta
#
solve for all fields $solveForAllFields
# dispersion model:
$dm
#
use super-grid absorbing layers $useSuperGrid 
super-grid width $superGridWidth
# printf("useSuperGrid = $useSuperGrid\n");
# pause
# 
#-# Drude params 1 1 all (gamma,omegap,domain-name)
#-if( $numBlocks eq 0 ){ $cmd="GDM params $a0 $a1 $b0 $b1 innerDomain (a0,a1,b0,b1,domain-name)\n"; }\
#-else{\
#-  $cmd="#"; \
#-  for( $i=0; $i<$numBlocks; $i++ ){ $cmd .= "\n GDM params $a0 $a1 $b0 $b1 blockDomain$i (a0,a1,b0,b1,domain-name)"; }\
#-}
#-$cmd 
#
# planeWaveInitialCondition
if( $leftBC eq "rbc" ){ $cmd = "planeWaveInitialCondition"; }else{ $cmd="zeroInitialCondition"; }
if( $ic eq "gp" ){ $cmd="Gaussian plane wave: $beta $x0 $y0 0 (beta,x0,y0,z0)\n gaussianPlaneWave"; }
if( $ic eq "gpw" ){ $cmd="gaussianPlaneWave\n Gaussian plane wave: $beta $x0 0 0 (beta,x0,y0,z0)"; }
if( $ic eq "gs"  ){ $cmd="gaussianSource\n Gaussian source: $beta $omega $x0 $y0 $z0 $amp $rampTime"; }
if( $ic eq "mgs"  ){ $a=1.; $t0=0; $p=1; $cmd="userDefinedForcing\n my source\n $a $beta $omega $p $x0 $&y0 $z0 $t0\n exit"; }
$cmd 
if( $checkErrors && $ic eq "pw" ){ $known="planeWaveKnownSolution"; }else{ $known="#"; }
$known
$kxa= abs($kx);
if( $kxa > 1. ){ $bbxb = int( $bbxb*$kxa +.5 )/$kxa; }  # we need to clip the plane wave on a period
if( $kx < 0 ){ $bbxa=$bbxb; $bbxb=100.; }
if( $leftBC eq "rbc" && $ic ne "gp" && $ic ne "gpw" ){ $cmd="initial condition bounding box $bbxa $bbxb $bbya $bbyb $zabb $zbbb"; }else{ $cmd="#"; }
$cmd
# initial condition bounding box $bbxa $bbxb $bbya $bbyb $zabb $zbbb
#  damp initial conditions at face (side,axis)=(0,1) of the box
bounding box decay face 1 0 
$betaBB=5.; # exponent in tanh function for smooth transition to zero outside the bounding box
bounding box decay exponent $betaBB
# 
$epsPW=$eps0; $muPW=$mu0; # parameters for the incident plane wave
plane wave coefficients $ax $ay $az $epsPW $muPW
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
# TEST: 
$cmd="#"; 
if( $bc1 ne "" ){  $cmd .= "\nbc: square(0,0)=$bc1"; }
if( $bc2 ne "" ){  $cmd .= "\nbc: square(1,0)=$bc2"; }
if( $bc3 ne "" ){  $cmd .= "\nbc: square(0,1)=$bc3"; }
if( $bc4 ne "" ){  $cmd .= "\nbc: square(0,1)=$bc4"; }
if( $bc5 ne "" ){  $cmd .= "\nbc: square(0,2)=$bc5"; }
if( $bc6 ne "" ){  $cmd .= "\nbc: square(0,2)=$bc6"; }
$cmd 
# bc: square(0,0)=symmetry
# bc: square(1,0)=absorbing
# bc: square(0,1)=symmetry
# bc: square(1,1)=symmetry
#
if( $leftBC eq "planeWave" ){ $cmd="bc: $backGround(0,0)=planeWaveBoundaryCondition"; }else{ $cmd="#"; }
$cmd 
if( $bcBody eq "pec" ){ $cmd="bc: annulus=perfectElectricalConductor"; }else{ $cmd="#"; }
$cmd
#
pml width,strength,power $pmlLines $pmlStrength $pmlPower
# 
# -- we need to subtract out the incident field on the "inflow" boundary before
#    applying the radiation boundary condition: 
if( $leftBC eq "planeWave" ){ $adjustFields=0; }else{ $adjustFields=1; }
adjust boundaries for incident field $adjustFields all
adjust boundaries for incident field $adjustFields $backGround
# 
# NOTE: material interfaces have share>=100
#
#  -- Assign material parameters ---
#  For multi-block grids we assume domain names of blockDomain0, blockDomain1, etc. 
#
@epsv = ( $eps1, $eps2, $eps3, $eps4 );
@muv = ( $mu1, $mu2, $mu3, $mu4 );
if( $numBlocks eq 0 ){ $cmd="coefficients $eps1 $mu1 innerDomain   (eps,mu,grid/domain-name)\n"; }\
else{\
  $cmd="#"; \
  for( $i=0; $i<$numBlocks; $i++ ){ $cmd .= "\n coefficients $epsv[$i] $muv[$i] blockDomain$i   (eps,mu,grid/domain-name)"; }\
}
$cmd 
coefficients $eps0 $mu0 outerDomain   (eps,mu,grid/domain-name)
GDM domain name: all
number of polarization vectors: $npv[0]
# -- specify material regions or material file
$cmd="#"; 
if( $materialFile ne "" ){ $cmd="include $materialFile"; }   # use materialFile values if provided 
$cmd
# printf(" dmfile[0]=$dmFile[0]\n");
# pause
#
if( $matFile ne "" ){ $cmd="material file: $matFile"; }else{ $cmd="#"; }
if( $materialFile ne "" ){ $cmd="material file: $dmFile[0]"; }   # use materialFile values if provided 
$cmd
#
#
$cmd="#"; 
if( $numMatRegions>1 ){\
$cmd = "include $regionFile";\
}
# 
# box: .75 1. .25 .75 0 1 (xa,xb,ya,yb,za,zb)
# material file: baAir.txt
# box: .25 .5 .5 .75 0 1 (xa,xb,ya,yb,za,zb)
# material file: baAir.txt
# cylinder: 0 .5 0 .25 0 1 (x0,y0,z0, radius, za,zb)
$cmd 
# 
# *****************
# for Yee we define the material regions
# if( $method eq "Yee" ){ $cmds = "define embedded bodies\n PEC cylinder\n $rad $x0 $y0 $z0\n exit"; }else{ $cmds="#"; }
if( $method eq "Yee" ){ $cmds = "define embedded bodies\n plane material interface\n 1 0 0 0 0 0\n $eps1 $mu1 0 0\n exit"; }else{ $cmds="#"; }
$cmds 
# ****************
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
debug $debug
#
cfl $cfl 
plot divergence 1
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
# --- *NEW* USER DEFINED PROBE
#- create user defined probe...
#-   $probeFile ="reflection$probeFileName.dat"; 
#-   file name $probeFile
#-   probe box $xar $xbr -2 2 -1 1 (xa,xb, ya,yb, za,zb)
#-   probe name reflectionProbe
#-   $L=-.25; # reflection coefficient is centered here 
#-   R/T offset $L
#-   incident amplitude 1.
#-   incident phase 0.
#-   reflection probe
#- exit
#- # 
#- create user defined probe...
#-   $probeFile ="transmission$probeFileName.dat"; 
#-   file name $probeFile
#-   probe box $xat $xbt -2 2 -1 1 (xa,xb, ya,yb, za,zb)
#-   probe name reflectionProbe
#-   $L=.25; # reflection coefficient is centered here 
#-   R/T offset $L
#-   transmission probe
#- exit
# 
#- # box probe
#- create a probe... 
#- open graphics
#- 
#- region probe
#- 
#- bounding box...
#-   bounding box grid 0
#-   bounding box 10 15 0 40 (i1a,i1b, i2a,i2b, i3a,i3b)
#-   number of layers 2
#- exit
#- probe name boxProbe
#- file name boxProbeFile.dat
#- sum
#- exit
# 
# -- surface integral probe --
#-- if( $intProbeName ne "" ){ $cmd="create a probe...\n   probe name $intProbeName\n  file name  $intProbeName.dat\n coordinate plane probe \n grid coordinate plane 0 .2512 0 0 (axis, x,y,z)\n   integral\n  all components\n exit"; }else{ $cmd="#"; }
#-- $cmd 
#-- #
#-- # Point probe: 
#-- create a probe...
#--   $leftProbeFileName="left$probeFileName.dat"; 
#--   file name $leftProbeFileName
#--   probe name leftProbe
#--   nearest grid point to $xLeftProbe $yLeftProbe 0
#-- exit
#-- create a probe...
#--   $rightProbeFileName="right$probeFileName.dat"; 
#--   file name $rightProbeFileName
#--   probe name rightProbe
#--   nearest grid point to $xRightProbe $yRightProbe 0
#-- exit
# **NEW** May 15, 2020
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
plot:Ey
contour
  plot contour lines (toggle)
  # vertical scale factor 0.2
  # min max -1.1 1.1
  if( $grid =~ /box/ ){ $cmd="contour shift 0.5\n -shift contour planes"; }else{ $cmd="#"; }
  $cmd 
exit
$go


plot:intensity
contour
  plot contour lines (toggle)
 # vertical scale factor 0.
 # min max .9 1.1
exit
echo to terminal 1
$go


#
plot:Ex
contour
plot contour lines (toggle)
vertical scale factor 0.
# min max -1.1 1.1
exit
$go
