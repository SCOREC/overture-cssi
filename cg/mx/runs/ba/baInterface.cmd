# 
# cgmx: test for BA material interfaces
#
# Usage:
#    cgmx [-noplot] baInterface -g=<name> -tf=<tFinal> -tp=<tPlot> -cfl=<> -diss=<> -eps1=<> -eps2=<> -ic=<> ...
#                   -bc=[abcEM2|rbcNonLocal|abcPML|perfectElectricalConductor|symmetry] ...
#                   -useNewInterface=[0|1] -method=[nfdtd|Yee|bmax] -errorNorm=[0|1|2] -interfaceIts=<i> ...
#                   -dm=[none|gdm] -npv=i i  -dispersionOption=[s|k] ...
#                   -caseName=[surfaceWaveCase1SpacePeriodic|surfaceWaveCase1TimePeriodic]
# 
# Examples:
# 
#  
# -- set default values for parameters ---
$kx=2; $ky=0; $kz=0; $left="leftSquare*"; $right="rightSquare*"; $degreex=2; $degreet=2; $method="NFDTD";
$tFinal=5.; $tPlot=.2; $cfl=.9; $show=" "; $interfaceIts=3; $debug=0; $diss=.1; $dissOrder=-1;
$useNewInterface=1; $errorNorm=1; $interfaceEquationOption=1; $interfaceOmega=.5; $setDivergenceAtInterfaces=0; 
$useImpedanceInterfaceProjection=1; $cons=1; $frequencyToFlush=10; 
$useSosupDissipation=0; $sosupParameter=1.; 
$dispersionOption="s"; 
$ts="me"; $matFile=""; $solveForAllFields=0;  $matFile2="baAir.txt";  $numMatRegions=1; $regionFile="";
$caseName="surfaceWaveCase1SpacePeriodic"; # case name for surface wave 
#
$eps1=1.; $mu1=1.;
$eps2=1.; $mu2=1.;
$ax=0.; $ay=0.; $az=0.; # plane wave coeffs. all zero -> use default
$ic = "zeroInitialCondition"; $bc1=""; $bc2=""; $bc3=""; $bc4=""; $bc5=""; $bc6=""; $bc7=""; $bc8=""; 
$pmic = "planeMaterialInterfaceInitialCondition";
$bc = "perfectElectricalConductor";
$x0=.5; $y0=0; $z0=0; $beta=50; # for Gaussian plane wave IC
@interfaceNormal = (); @interfacePoint = ();
$tz = "#"; $go="halt";
$dm="none"; $alphaP = (); @npv=();  $modeGDM=-1; 
@a01 = (); @a11=(); @b01=(); @b11=(); # these must be null for GetOptions to work, defaults are given below
@a02 = (); @a12=(); @b02=(); @b12=(); 
# ----------------------------- get command line arguments ---------------------------------------
#  -- first get any commonly used options: (requires the env variable CG to be set)
# $getCommonOptions = "$ENV{'CG'}/mp/cmd/getCommonOptions.h";
# include $getCommonOptions
#  -- now get additional options: 
GetOptions("bc=s"=>\$bc,"cfl=f"=>\$cfl,"debug=i"=>\$debug,"diss=f"=>\$diss,"eps1=f"=>\$eps1,"eps2=f"=>\$eps2,\
           "go=s"=>\$go,"g=s"=>\$grid,"ic=s"=>\$ic,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"degreex=f"=>\$degreex,\
           "degreet=f"=>\$degreet,"method=s"=>\$method,"errorNorm=i"=>\$errorNorm,"dissOrder=i"=>\$dissOrder,\
           "interfacePoint=f{1,}"=>\@interfacePoint,"interfaceNormal=f{1,}"=>\@interfaceNormal,"ax=f"=>\$ax,"ay=f"=>\$ay,"az=f"=>\$az,\
           "left=s"=>\$left,"right=s"=>\$right,"restart=s"=>\$restart,"interfaceIts=i"=>\$interfaceIts,\
           "show=s"=>\$show,"tf=f"=>\$tFinal,"tp=f"=>\$tPlot,"tz=s"=>\$tz,"useNewInterface=i"=>\$useNewInterface,\
           "interfaceEquationOption=i"=>\$interfaceEquationOption,"interfaceOmega=f"=>\$interfaceOmega,"tz=s"=>\$tz,\
           "bc1=s"=>\$bc1,"bc2=s"=>\$bc2,"bc3=s"=>\$bc3,"bc4=s"=>\$bc4,"bc5=s"=>\$bc5,"bc6=s"=>\$bc6,\
           "bc7=s"=>\$bc7,"bc8=s"=>\$bc8,"setDivergenceAtInterfaces=s"=>\$setDivergenceAtInterfaces,"cons=i"=>\$cons,\
           "useImpedanceInterfaceProjection=s"=>\$useImpedanceInterfaceProjection,"modeGDM=i"=>\$modeGDM,"alphaP=f{1,}"=>\@alphaP,\
           "dm=s"=>\$dm,"npv=i{1,}"=>\@npv,"useSosupDissipation=i"=>\$useSosupDissipation,"sosupParameter=f"=>\$sosupParameter,\
           "a01=f{1,}"=>\@a01,"a11=f{1,}"=>\@a11,"b01=f{1,}"=>\@b01,"b11=f{1,}"=>\@b11,\
           "a02=f{1,}"=>\@a02,"a12=f{1,}"=>\@a12,"b02=f{1,}"=>\@b02,"b12=f{1,}"=>\@b12,\
	   "x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"frequencyToFlush=i"=>\$frequencyToFlush,\
	   "matFile=s"=>\$matFile,"matFile2=s"=>\$matFile2,"numMatRegions=i"=>\$numMatRegions,"regionFile=s"=>\$regionFile,\
	   "ts=s"=>\$ts,"solveForAllFields=i"=>\$solveForAllFields, "dispersionOption=s"=>\$dispersionOption,\
	   "caseName=s"=>\$caseName );
# -------------------------------------------------------------------------------------------------
#
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
if( $ic eq "pmic" ){ $ic="planeMaterialInterfaceInitialCondition"; }
if( $ic eq "tz" ){ $ic="twilightZone"; }
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $tz eq "" ){ $tz="#"; }
if( $tz eq "poly" ){ $tz="polynomial"; }
if( $tz eq "trig" ){ $tz="trigonometric"; }
if( $method eq "sosup" ){ $diss=0.; }
# 
if( $ts eq "me" ){ $ts="modifiedEquationTimeStepping"; }
$orderOfRungeKutta=4;
if( $ts eq "rk1" ){ $ts="rungeKutta"; $orderOfRungeKutta=1; }
if( $ts eq "rk2" ){ $ts="rungeKutta"; $orderOfRungeKutta=2; }
if( $ts eq "rk3" ){ $ts="rungeKutta"; $orderOfRungeKutta=3; }
if( $ts eq "rk4" ){ $ts="rungeKutta"; $orderOfRungeKutta=4; }
# printf(" TS = [$ts], order=[$orderOfRungeKutta] \n");
#
if( $dm eq "none" ){ $dm="no dispersion"; }
if( $dm eq"gdm" ){ $dm="GDM"; }
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
# time-stepping method
$ts
order of Runge Kutta $orderOfRungeKutta
#
solve for all fields $solveForAllFields
#
# dispersion model:
$dm
# printf(" dm=$dm\n");
#
kx,ky,kz $kx $ky $kz
plane wave coefficients $ax $ay $az $eps1 $mu1
#
use new interface routines $useNewInterface
set divergence at interfaces $setDivergenceAtInterfaces
use impedance interface projection $useImpedanceInterfaceProjection
# These next parameters define the exact solution for a material interface: 
material interface point $interfacePoint[0] $interfacePoint[1] $interfacePoint[2] 
material interface normal $interfaceNormal[0] $interfaceNormal[1] $interfaceNormal[2] 
#
# ****************
bc: all=$bc
#** bc: all=dirichlet
#* bc: leftSquare(0,0)=planeWaveBoundaryCondition
# ****************
if( $bc1 ne "" ){ $cmd="bc: bcNumber1=$bc1"; }else{ $cmd="#"; }
$cmd
if( $bc2 ne "" ){ $cmd="bc: bcNumber2=$bc2"; }else{ $cmd="#"; }
$cmd
if( $bc3 ne "" ){ $cmd="bc: bcNumber3=$bc3"; }else{ $cmd="#"; }
$cmd
if( $bc4 ne "" ){ $cmd="bc: bcNumber4=$bc4"; }else{ $cmd="#"; }
$cmd
if( $bc5 ne "" ){ $cmd="bc: bcNumber5=$bc5"; }else{ $cmd="#"; }
$cmd
if( $bc6 ne "" ){ $cmd="bc: bcNumber6=$bc6"; }else{ $cmd="#"; }
$cmd
if( $bc7 ne "" ){ $cmd="bc: bcNumber7=$bc7"; }else{ $cmd="#"; }
$cmd
if( $bc8 ne "" ){ $cmd="bc: bcNumber8=$bc8"; }else{ $cmd="#"; }
$cmd
# 
#
# option: 1=extrapolate as initial guess for material interface ghost values
interface option 1
# interfaceEquationsOption=0 : use extrap for 2nd ghost, 1=use eqns
interface equations option $interfaceEquationOption
omega for interface iterations $interfaceOmega
#
interface BC iterations $interfaceIts
#
#
# *****************
#   --- Define eps and mu in the difference domains -----
$cmds="#";
if( $method eq "Yee" ){ $cmds = "define embedded bodies\n plane material interface\n $interfaceNormal $interfacePoint\n $eps2 $mu2 0. 0. \nexit"; }
$cmds 
$cmd="#";
if( $method ne "Yee" ){ $cmd = \
  "coefficients $eps1 1. leftDomain (eps,mu,grid-name)\n" . \
  "coefficients $eps2 1. rightDomain (eps,mu,grid-name)\n"; }
$cmd 
#
GDM mode: $modeGDM
GDM domain name: all
# -- specify material regions or material file
material file: $matFile
edit material
 # Specify option when solving the dispersion relation
 if( $dispersionOption eq "s" ){ $cmd="compute s given k"; }else{ $cmd="compute k given s"; }
 $cmd
exit
# 
forcing options...
define material region...
material file: $matFile2
# NOTE: material box must be big enough to include ghost points
if( $ic ne "isw" ){ $cmd="box: .00000001 1.1 -.55 .55 -.55 .55 (xa,xb,ya,yb,za,zb)"; }else{ $cmd="box: -2.1 2.1 -1.05 -.0000001 -.55 .55 (xa,xb,ya,yb,za,zb)"; }
$cmd 
continue
# 
#- $domain="all"; 
#- # ------------ Set GDM parameters on the left domain -----------
#- GDM domain name: leftDomain
#-   number of polarization vectors: $npv[0]
#-   GDM alphaP: $alphaP[0]
#- $cmd="#"; 
#- if( $npv[0] == 1 ){ \
#-    $cmd = " GDM coeff: 0 $a01[0] $a11[0] $b01[0] $b11[0] (eqn, a0,a1,b0,b1)\n"; \
#-  }
#- if( $npv[0] == 2 ){ \
#-    $cmd  = " GDM coeff: 0 $a01[0] $a11[0] $b01[0] $b11[0] (eqn, a0,a1,b0,b1)\n"; \
#-    $cmd .= " GDM coeff: 1 $a01[1] $a11[1] $b01[1] $b11[1] (eqn, a0,a1,b0,b1)"; \
#-       }
#- if( $npv[0] == 3 ){ \
#-    $cmd  = " GDM coeff: 0 $a01[0] $a11[0] $b01[0] $b11[0] (eqn, a0,a1,b0,b1)\n"; \
#-    $cmd .= " GDM coeff: 1 $a01[1] $a11[1] $b01[1] $b11[1] (eqn, a0,a1,b0,b1)\n"; \
#-    $cmd .= " GDM coeff: 2 $a01[2] $a11[2] $b01[2] $b11[2] (eqn, a0,a1,b0,b1)"; \
#-       }
#- $cmd
#- # ------------ Set GDM parameters on the right domain -----------
#- GDM domain name: rightDomain
#-   number of polarization vectors: $npv[1]
#-   GDM alphaP: $alphaP[1]
#- $cmd="#"; 
#- if( $npv[1] == 1 ){ \
#-    $cmd = " GDM coeff: 0 $a02[0] $a12[0] $b02[0] $b12[0] (eqn, a0,a1,b0,b1)\n"; \
#-  }
#- if( $npv[1] == 2 ){ \
#-    $cmd  = " GDM coeff: 0 $a02[0] $a12[0] $b02[0] $b12[0] (eqn, a0,a1,b0,b1)\n"; \
#-    $cmd .= " GDM coeff: 1 $a02[1] $a12[1] $b02[1] $b12[1] (eqn, a0,a1,b0,b1)"; \
#-       }
#- if( $npv[1] == 3 ){ \
#-    $cmd  = " GDM coeff: 0 $a02[0] $a12[0] $b02[0] $b12[0] (eqn, a0,a1,b0,b1)\n"; \
#-    $cmd .= " GDM coeff: 1 $a02[1] $a12[1] $b02[1] $b12[1] (eqn, a0,a1,b0,b1)\n"; \
#-    $cmd .= " GDM coeff: 2 $a02[2] $a12[2] $b02[2] $b12[2] (eqn, a0,a1,b0,b1)"; \
#-       }
#- $cmd
#
# 
#  The dispersive case is handled by a user defined known solution
# if( $ic ne "gp" && $dm ne "no dispersion" && $tz eq "#" ){
if( $ic ne "gp" && $tz eq "#" ){\
 $ic = "user defined known solution\n  dispersive plane wave interface\n done\n userDefinedKnownSolutionInitialCondition"; }
# Surface wave: 
if( $ic ne "isw" && $tz eq "#" ){\
 $ic = "user defined known solution\n interface surface wave\n $caseName\n done\n userDefinedKnownSolutionInitialCondition"; }
# 
bc: all=$bc
if( $ic eq "gp" ){ $ic="Gaussian plane wave: $beta $x0 $y0 0 (beta,x0,y0,z0)\n gaussianPlaneWave"; }
$ic
# printf(" ic=$ic\n");
#
$tz 
** trigonometric
 degreeSpace, degreeTime  $degreex $degreet
#
# bc: Annulus=perfectElectricalConductor
tFinal $tFinal
tPlot $tPlot
#
order of dissipation $dissOrder
dissipation $diss
#
use sosup dissipation $useSosupDissipation
sosup parameter $sosupParameter
#
use conservative difference $cons
debug $debug
#
cfl $cfl 
plot errors 1
check errors 1
error norm $errorNorm
#*********************************
show file options...
# MXSF:compressed
MXSF:open
  $show
# MXSF:frequency to save 1
MXSF:frequency to flush $frequencyToFlush
exit
#**********************************
continue
#
plot:Ey
# Move contour planes in 3D 
if( $grid =~ /^box3d/ ){ $cmd="contour\n contour shift 0.5\n -shift contour planes\n exit"; }else{ $cmd="#"; }
$cmd
# 
$go


movie mode
finish
