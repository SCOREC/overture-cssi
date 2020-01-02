#================================================================================================
#  cgmx example:  scattering from a dielectric cylinder or sphere (and compare to the exact solution)
#
# Usage:
#   
#  cgmx [-noplot] dielectricCyl  -g=<name> -tf=<tFinal> -tp=<tPlot> -kx=<num> -ky=<num> -kz=<num> -show=<name> ...
#                                -eps1=<> -eps2=<> -interfaceIts=<> -diss=<> -filter=[0|1] -dissc=<> -debug=<num> ...
#                                -cons=[0/1] -plotIntensity=[0|1] ...
#        -useSosupDissipation=[0|1] -sosupDissipationOption=[0|1] -sosupDissipationFrequency=<i> 
#                                -method=[nfdtd|Yee|sosup|bamx] -errorNorm=[0|1|2] -go=[run/halt/og]
# Arguments:
#  -kx= -ky= -kz= : integer wave numbers of the incident wave
#  -interit : number of iterations to solve the interface equations 
#  -errorNorm:  set to 1 or 2 to show L1 and L2 norm errors
#  -diss, -dissc : coefficients of art. dissipation. If dissc>=0 then use this for curvilinear grids. 
#
echo to terminal 0
# Examples: (see Makefile to build grids)
#  -- NEW 8th order filter: smaller errors and more stable for fewer -interit 
#   cgmx dielectricCyl -g=innerOutere8.order4 -kx=2 -eps1=.25 -eps2=1. -go=halt -filter=1 -tp=.5 -tf=10 -interit=1
#   cgmx dielectricCyl -g=innerOutere4.order2 -kx=2 -eps1=.25 -eps2=1. -go=halt -dissOrder=4 -tp=.5 -tf=10 
#   cgmx dielectricCyl -g=innerOutere4.order4 -kx=1.25 -eps1=2.25 -eps2=1. -go=halt -diss=0. -filter=1 -tp=.1 -tf=10 -interit=1
#
#  -- broken:
#   cgmx dielectricCyl -g=innerOutere2.order4.hdf -kx=2 -eps1=.25 -eps2=1. -go=halt -tp=.01
# 
#  --NEW: build grids with domains
#   cgmx dielectricCyl -g=innerOutere2.order4.hdf -kx=2 -eps1=.25 -eps2=1. -go=halt
#   cgmx dielectricCyl -g=innerOutere4.order4.hdf -kx=2 -eps1=.25 -eps2=1. -go=halt
#   cgmx dielectricCyl -g=innerOutere8.order4.hdf -kx=2 -eps1=.25 -eps2=1. -go=halt
# 
#   cgmx dielectricCyl -g=innerOutere4.order4.hdf -kx=2 -eps1=2.25 -eps2=1. -go=halt
# 
#   mpirun -np 2 $cgmxp dielectricCyl -g=innerOutere2.order2.hdf -kx=2 -eps1=.25 -eps2=1. -useNewInterface=1 -go=halt
#   mpirun -np 2 $cgmxp dielectricCyl -g=$ov/ogen.p/innerOuter4.order2parallel.hdf -kx=2 -eps1=.25 -eps2=1. -useNewInterface=1 -go=halt
# 
#   srun -N1 -n1 -ppdebug $cgmxp noplot dielectricCyl -g=innerOutere2.order2.hdf -kx=2 -eps1=.25 -eps2=1. -tf=.1 -useNewInterface=1 -go=go
# 
# -- Yee scheme : dielectric cylinder:
#   cgmx dielectricCyl -g=bigSquareSize1f4.hdf -kx=2 -eps1=.25 -eps2=1. -method=Yee -errorNorm=2 -tp=.01 -go=halt
#   cgmx dielectricCyl -g=bigSquareSize1f8.hdf -kx=2 -eps1=.25 -eps2=1. -method=Yee -errorNorm=2 -go=halt
#   cgmx dielectricCyl -g=bigSquareSize1f16.hdf -kx=2 -eps1=.25 -eps2=1. -method=Yee -errorNorm=2 -go=halt
#   cgmx dielectricCyl -g=bigSquareSize1f32.hdf -kx=2 -eps1=.25 -eps2=1. -method=Yee -errorNorm=2 -go=halt
# -- Yee scheme : dielectric sphere
#   cgmx dielectricCyl -cyl=0 -g=bigBox1.order2 -kx=1 -eps1=.25 -eps2=1. -method=Yee -errorNorm=2 -tp=.01 -go=halt
#   cgmx dielectricCyl -cyl=0 -g=bigBox2.order2 -kx=1 -eps1=.25 -eps2=1. -method=Yee -errorNorm=2 -tp=.01 -go=halt
#   cgmx dielectricCyl -cyl=0 -g=bigBox8.order2 -kx=1 -eps1=.25 -eps2=1. -method=Yee -errorNorm=2 -cfl=.5 -tp=.01 -go=halt
#
# -- sosup
#   cgmx dielectricCyl -g=innerOutere4.order4.hdf -kx=2 -eps1=.25 -eps2=1. -method=sosup -go=halt
#
# -- OLD: 
#   cgmx dielectricCyl -g=innerOuter4.order2.hdf -kx=2 -eps1=.25 -eps2=1. -go=halt
#   cgmx dielectricCyl -g=innerOuter4.order4.hdf -kx=2 -eps1=.25 -eps2=1. -go=halt
#   cgmx dielectricCyl -g=innerOuter8.order4.hdf -kx=2 -eps1=2.  -eps2=1. -go=halt
#   cgmx dielectricCyl -g=innerOuter16.order4.hdf -kx=2 -eps1=2. -eps2=1. -go=halt
#
#  -- dielectric sphere: eps1=inside eps2=outside (should be 1)
#   cgmx dielectricCyl -cyl=0 -g=solidSphereInABoxi1.order2 -kx=1 -eps1=.25 -eps2=1. -go=halt
#   cgmx dielectricCyl -cyl=0 -g=solidSphereInABoxi2.order2 -kx=1 -eps1=.25 -eps2=1. -go=halt
#   cgmx dielectricCyl -cyl=0 -g=solidSphereInABoxi1.order2 -kx=1 -eps1=1. -eps2=1. -go=halt -tp=.01
# 
#   cgmx dielectricCyl -cyl=0 -g=solidSphereInABoxe2.order2 -kx=1 -eps1=.25 -eps2=1. -go=halt 
#   cgmx dielectricCyl -cyl=0 -g=solidSphereInABoxe4.order2 -kx=1 -eps1=.25 -eps2=1. -go=halt
# 
# -- fourth order:
#   cgmx dielectricCyl -cyl=0 -g=solidSphereInABoxi2.order4 -kx=1 -eps1=.25 -eps2=1. -go=halt  [ok]
#   cgmx noplot dielectricCyl -cyl=0 -g=solidSphereInABoxi4.order4 -kx=1 -eps1=.25 -eps2=1. -tf=.2 -go=go 
#  -- fix grid dimensions for convergence tests:
#   cgmx noplot dielectricCyl -cyl=0 -g=solidSphereInABoxFixedi1.order2 -kx=1 -eps1=.25 -eps2=1. -go=go -tf=.5
#   cgmx noplot dielectricCyl -cyl=0 -g=solidSphereInABoxFixedi2.order2 -kx=1 -eps1=.25 -eps2=1. -go=go -tf=.5
# --- new sphere grids
#   cgmx noplot dielectricCyl -cyl=0 -g=solidSphereInABoxNewi1.order4 -kx=1 -eps1=.25 -eps2=1. -go=og
#   cgmx dielectricCyl -cyl=0 -g=solidSphereInABoxe2.order4 -kx=1 -eps1=.25 -eps2=1. -go=halt -tp=.01
#   cgmx dielectricCyl -cyl=0 -g=solidSphereInABoxe4.order4 -kx=1 -eps1=.25 -eps2=1. -go=halt -tp=.01 
# 
# parallel: 
#  srun -N1 -n1 -ppdebug $cgmxp noplot dielectricCyl -cyl=0 -g=solidSphereInABoxe2.order2 -kx=1 -eps1=.25 -eps2=1. -tp=.1 -tf=.2 -go=go
#  srun -N1 -n1 -ppdebug $cgmxp noplot dielectricCyl -cyl=0 -g=solidSphereInABoxe2.order4 -kx=1 -eps1=.25 -eps2=1. -tp=.1 -tf=.2 -go=go
#  totalview srun -a -N1 -n1 -ppdebug $cgmxp noplot ...
#  srun -ppdebug -N2 -n2 memcheck_all $cgmxp noplot ...
# diss=.5 not enough for:
#  srun -N1 -n1 -ppdebug $cgmxp noplot dielectricCyl -cyl=0 -g=solidSphereInABoxe8.order4 -kx=1 -eps1=.25 -eps2=1. -tp=.1 -tf=.5 -diss=2. -go=go
#================================================================================================
# 
$tFinal=2.; $tPlot=.1;  $show=" "; $method="NFDTD"; $ts="me"; 
$cfl = .8; $diss=.5; $dissOrder=-1; $filter=0; $dissc=-1.; $debug=0;  $plotIntensity=0;
$cyl=1;   # set to 0 for a sphere 
$kx=2; $ky=0; $kz=0;
$eps1=.25; $mu1=1.; # inner
$eps2=1.;  $mu2=1.; # outer 
$show=" "; $backGround="backGround"; $useNewInterface=1; 
$interfaceEquationOption=1; $interfaceIterations=5;  $interfaceOmega=.5;
$grid="innerOuter4.order4.hdf";
$cons=0; $go="halt";  $errorNorm=0;
$flushFrequency=8; 
$dm="none"; @npv=();  $modeGDM=-1; 
$alphaP = (); 
@a01 = (); @a11=(); @b01=(); @b11=(); # these must be null for GetOptions to work, defaults are given below
@a02 = (); @a12=(); @b02=(); @b12=(); 
$dmFile=""; # "SilverJCDispersionFits.txt"; 
$lengthScale=1.e-7; # length-scale = 100 nm 
$sphereRadius=1.; 
#
$stageOption ="default";
$useSosupDissipation=0; $sosupParameter=1.;  $sosupDissipationOption=1; $sosupDissipationFrequency=1;
$selectiveDissipation=0;
#
$drOption="computeFrequency"; # dispersion relation option 
$matFile=""; $solveForAllFields=0; $numMatRegions=1; $matFile2="";
$regionFile="cylRegion.h"; $xc=0.; $yc=0.; $zc=0; $radius=.4; # case of sphere is handled below
# 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"diss=f"=>\$diss,"dissc=f"=>\$dissc,"tp=f"=>\$tPlot,"show=s"=>\$show,"debug=i"=>\$debug, \
 "cfl=f"=>\$cfl, "bg=s"=>\$backGround,"bcn=s"=>\$bcn,"go=s"=>\$go,"noplot=s"=>\$noplot,"plotIntensity=i"=>\$plotIntensity,\
 "interfaceIts=i"=>\$interfaceIterations,"cyl=i"=>\$cyl,"useNewInterface=i"=>\$useNewInterface,"errorNorm=i"=>\$errorNorm,\
 "dtMax=f"=>\$dtMax,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"eps1=f"=>\$eps1,"eps2=f"=>\$eps2, "cons=i"=>\$cons,\
 "method=s"=>\$method,"dissOrder=i"=>\$dissOrder,"filter=i"=>\$filter,"flushFrequency=i"=>\$flushFrequency,\
 "interfaceEquationOption=i"=>\$interfaceEquationOption,"interfaceOmega=f"=>\$interfaceOmega,"lengthScale=f"=>\$lengthScale,\
  "useSosupDissipation=i"=>\$useSosupDissipation,"sosupParameter=f"=>\$sosupParameter,\
  "sosupDissipationOption=i"=>\$sosupDissipationOption,"sosupDissipationFrequency=i"=>\$sosupDissipationFrequency,\
  "selectiveDissipation=i"=>\$selectiveDissipation,"modeGDM=i"=>\$modeGDM,"stageOption=s"=>\$stageOption,\
  "dm=s"=>\$dm,"npv=i{1,}"=>\@npv,"alphaP=f{1,}"=>\@alphaP,\
  "a01=f{1,}"=>\@a01,"a11=f{1,}"=>\@a11,"b01=f{1,}"=>\@b01,"b11=f{1,}"=>\@b11,\
  "a02=f{1,}"=>\@a02,"a12=f{1,}"=>\@a12,"b02=f{1,}"=>\@b02,"b12=f{1,}"=>\@b12,\
   "dmFile=s"=>\$dmFile,"ii=i"=>\$interfaceIterations,"sphereRadius=f"=>\$sphereRadius,"ts=s"=>\$ts,\
   "matFile=s"=>\$matFile,"matFile2=s"=>\$matFile2,"solveForAllFields=i"=>\$solveForAllFields, "drOption=s"=>\$drOption,\
   "numMatRegions=i"=>\$numMatRegions,"regionFile=s"=>\$regionFile );
# -------------------------------------------------------------------------------------------------
if( $method eq "bamx" ){ $numMatRegions=2; }
if( $method eq "bamx" && $cyl eq 0  ){ $solveForAllFields=1; $regionFile="sphereRegion.h"; $radius=1.; $ae=$radius; $be=$radius; $ce=$radius; }
if( $method eq "sosup" ){ $diss=0.; }
if( $method eq "fd" ){ $method="nfdtd"; }
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
# 
if( $ts eq "me" ){ $ts="modifiedEquationTimeStepping"; }
$orderOfRungeKutta=4;
if( $ts eq "rk1" ){ $ts="rungeKutta"; $orderOfRungeKutta=1; }
if( $ts eq "rk2" ){ $ts="rungeKutta"; $orderOfRungeKutta=2; }
if( $ts eq "rk3" ){ $ts="rungeKutta"; $orderOfRungeKutta=3; }
if( $ts eq "rk4" ){ $ts="rungeKutta"; $orderOfRungeKutta=4; }
# 
#
if( $dm eq "none" ){ $dm="no dispersion"; }
if( $dm eq"gdm" ){ $dm="GDM"; }
# Give defaults here for array arguments: 
if( $alphaP[0] eq "" ){ @alphaP=(-1,-1); } # default -1 means use 1/eps
# 
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
# 
echo to terminal 1
#
$grid
#
$method
# time-stepping method
$ts
order of Runge Kutta $orderOfRungeKutta
#
solve for all fields $solveForAllFields
# Set length scale (used by dispersion models)
length scale: $lengthScale
# dispersion model:
$dm
$drOption
# printf(" dm=$dm\n");
#
GDM mode: $modeGDM
#
# ----- SET eps and mu BEFORE GDM parameters so alphaP=1/eps by default ------
#      innerAnnulus
#      innerSquare
#      outerAnnulus
#      outerSquare
# NOTE: material interfaces have share>=100
$cmd="#";
if( $cyl eq 1 && $method ne "Yee" ){ $cmd = \
  "coefficients $eps1 1. innerAnnulus* (eps,mu,grid-name)\n" .\
  "coefficients $eps1 1. innerSquare (eps,mu,grid-name)\n" .\
  "coefficients $eps2 1. outerAnnulus* (eps,mu,grid-name)\n" .\
  "coefficients $eps2 1. outerSquare (eps,mu,grid-name)\n"; }
if( $cyl eq 0 && $method ne "Yee" ){ $cmd = \
  "coefficients $eps1 1. inner* (eps,mu,grid-name)\n" .\
  "coefficients $eps2 1. outer* (eps,mu,grid-name)\n";}
$cmd
# 
$domain="all"; 
# ------------ Set GDM parameters in the outer domain -----------
GDM domain name: outerDomain
  number of polarization vectors: $npv[0]
  GDM alphaP: $alphaP[0]
$cmd="#"; 
if( $npv[0] == 1 ){ \
   $cmd = " GDM coeff: 0 $a01[0] $a11[0] $b01[0] $b11[0] (eqn, a0,a1,b0,b1)\n"; \
 }
if( $npv[0] == 2 ){ \
   $cmd  = " GDM coeff: 0 $a01[0] $a11[0] $b01[0] $b11[0] (eqn, a0,a1,b0,b1)\n"; \
   $cmd .= " GDM coeff: 1 $a01[1] $a11[1] $b01[1] $b11[1] (eqn, a0,a1,b0,b1)"; \
      }
$cmd
# ------------ Set GDM parameters in the inner domain -----------
GDM domain name: innerDomain
  number of polarization vectors: $npv[1]
  GDM alphaP: $alphaP[1]
$cmd="#"; 
if( $npv[1] == 1 ){ \
   $cmd = " GDM coeff: 0 $a02[0] $a12[0] $b02[0] $b12[0] (eqn, a0,a1,b0,b1)\n"; \
 }
if( $npv[1] == 2 ){ \
   $cmd  = " GDM coeff: 0 $a02[0] $a12[0] $b02[0] $b12[0] (eqn, a0,a1,b0,b1)\n"; \
   $cmd .= " GDM coeff: 1 $a02[1] $a12[1] $b02[1] $b12[1] (eqn, a0,a1,b0,b1)"; \
      }
$cmd
#
# -- specify material regions or material file
$cmd="#"; 
if( $matFile ne "" ){ $cmd="GDM domain name: all\nmaterial file: $matFile"; }else{ $cmd="#"; }
$cmd
#
$cmd="#"; 
if( $numMatRegions>1 ){\
$cmd = "include $regionFile";\
}
$cmd
# 
# 
#* planeWaveInitialCondition
# ++ zeroInitialCondition
# ==== in 3d solve for the full field ===
$ic = $cyl ? "planeWaveScatteredFieldInitialCondition" : "planeWaveInitialCondition";
# $ic
planeWaveScatteredFieldInitialCondition
$known = $cyl ? "scatteringFromADielectricDiskKnownSolution" : "scatteringFromADielectricSphereKnownSolution";
$known
# *****************
# for Yee we define the cylinder as a masked stair step region
if( $cyl eq 1 ){ $rad=.4; }else{ $rad=$sphereRadius; }
$x0=0.; $y0=0; $z0=0;  
$cmds="#";
if( $cyl eq 1 && $method eq "Yee" ){ $cmds = "define embedded bodies\n dielectric cylinder\n $rad $x0 $y0 $z0\n $eps1 $mu1 0. 0. \nexit"; }
if( $cyl eq 0 && $method eq "Yee" ){ $cmds = "define embedded bodies\n dielectric sphere\n $rad $x0 $y0 $z0\n $eps1 $mu1 0. 0. \nexit"; }
$cmds 
scattering radius $rad
# ****************
# ====
# twilightZone
#  degreeSpace, degreeTime  1 1
#
kx,ky,kz $kx $ky $kz
# 
use new interface routines $useNewInterface
#
bc: all=dirichlet
# --
#bc: all=abcEM2
#bc: outerBox(0,0)=planeWaveBoundaryCondition
# --
# ++ bc: all=perfectElectricalConductor
# ++ bc: outerSquare(0,0)=planeWaveBoundaryCondition
# 
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
# this is broken: order of dissipation 6
order of dissipation $dissOrder
apply filter $filter
dissipation $diss
dissipation (curvilinear) $dissc
#
if( $selectiveDissipation eq "1" ){ $cmd="selective dissipation...\n  turn off rectangular\n continue"; }else{ $cmd="#"; }
$cmd 
#
use sosup dissipation $useSosupDissipation
sosup parameter $sosupParameter
sosup dissipation option $sosupDissipationOption
sosup dissipation frequency $sosupDissipationFrequency
# -- set default stage options
if( $stageOption eq "default" && $useSosupDissipation eq 0 ){ $stageOption="IDB"; }
if( $stageOption eq "default" && $useSosupDissipation eq 1 ){ $stageOption="D-IB"; }
#
# --- Define multi-stage time-step: 
if( $stageOption eq "IDB" ){ $stages="updateInterior,addDissipation,applyBC"; }
if( $stageOption eq "D-IB" ){ $stages="addDissipation\n updateInterior,applyBC"; }
if( $stageOption eq "DB-IB" ){ $stages="addDissipation,applyBC\n updateInterior,applyBC"; }
if( $stageOption eq "IB-DB" ){ $stages="updateInterior,applyBC\n addDissipation,applyBC"; }
if( $stageOption eq "IB-D" ){ $stages="updateInterior,applyBC\n addDissipation"; }
# -- options to precompute V=uDot used in the dissipation
if( $stageOption eq "IVDB" ){ $stages="updateInterior,computeUt,addDissipation,applyBC"; }
if( $stageOption eq "VD-IB" ){ $stages="computeUt,addDissipation\n updateInterior,applyBC"; }
if( $stageOption eq "VDB-IB" ){ $stages="computeUt,addDissipation,applyBC\n updateInterior,applyBC"; }
if( $stageOption eq "IB-VDB" ){ $stages="updateInterior,applyBC\n computeUt,addDissipation,applyBC"; }
if( $stageOption eq "IB-VD" ){ $stages="updateInterior,applyBC\n computeUt,addDissipation"; }
if( $method ne "bamx" ){ $cmd="set stages...\n  $stages\n done"; }else{ $cmd="#"; }
$cmd
#
#
use conservative difference $cons 
debug $debug
#
cfl $cfl 
plot errors 1
check errors 1
error norm $errorNorm
plot intensity $plotIntensity
# 
#*********************************
show file options...
  MXSF:compressed
  MXSF:open
    $show
  MXSF:frequency to flush $flushFrequency
exit
#**********************************
continue
#
plot:Ey
# Move contour planes in 3D 
if( $cyl eq 0 ){ $cmd="contour\n delete contour plane 0\n contour shift 0.5\n -shift contour planes\n exit"; }else{ $cmd="#"; }
$cmd
# 
$go
