#================================================================================================
#
#  Master cgmx command file script    (Created March 2021)
#
echo to terminal 1
#
# Usage:
#   
#  cgmx [-noplot] mxScript OPTIONS 
#
# OPTIONS:
#   -g=<name>
#   -method=[fd|Yee|sosup] 
#   -ic : initial condition, pw=plane-wave, gp=Gaussian-pulse
#   -dm=[none|gdm] : dispersion model
#   -tf=<tFinal> : final time 
#   -tp=<tPlot>  : plot times 
#   -kx=<num> -ky=<num> -kz=<num> : wave number used in incident field and for other purposes
#   -bcCmds="b1" "b2" "bc2" ...   : array of bc commands
#   -plotIntensity=[0|1]
#   -diss=<>  -filter=[0|1]
#   -debug=<num> -cons=[0/1]
#   -varDiss=<0|1> ...
#   -rbc=[abcEM2|rbcNonLocal|abcPML]
#   -leftBC=[rbc|planeWave]
#   -probeFileName=<s> -xLeftProbe=<f> -xRightProbe=<f> ...
#   -useSosupDissipation=[0|1]
#   -extraCommandsEnd=<file.h>  : extra commands to include at the end
#
#================================================================================================
# 
$eps0=1.; $mu0=1.; # outer domain 
$backGround="backGround";
$leftBC="rbc"; $bcBody=""; 
#
# ---- specify materials in a file: -----
# 
$materialFile="diskMaterials.h";
$materialFile="diskMaterials.h";
# number of bodies to set materials for:
$numBodies=1; 
#
# --------------------------------------------------
# ------- get common command line arguments --------
include getCommonOptions.h
# --------------------------------------------------
include $extraCommandsStart
# 
$grid
#
$method
# dispersion model:
$dm
# nonlinear model:
$nm
# 
# planeWaveInitialCondition
if( $leftBC eq "rbc" ){ $cmd = "planeWaveInitialCondition"; }else{ $cmd="zeroInitialCondition"; }
if( $ic eq "gp" ){ $cmd="Gaussian plane wave: $beta $x0 $y0 0 $k0 (beta,x0,y0,z0,k0)\n gaussianPlaneWave"; }
if( $ic eq "gpw" ){ $cmd="gaussianPlaneWave\n Gaussian plane wave: $beta $x0 0 0 $k0 (beta,x0,y0,z0,k0)"; }
$cmd 
# 
if( $checkErrors ){ $known="planeWaveKnownSolution"; }else{ $known="#"; }
$known
# 
$kxa= abs($kx);
if( $kxa > 1. ){ $xbbb = int( $xbbb*$kxa +.5 )/$kxa; }  # we need to clip the plane wave on a period
if( $kx < 0 ){ $xabb=$xbbb; $xbbb=100.; }
if( $leftBC eq "rbc" && $ic ne "gp"  && $ic ne "gpw" ){ $cmd="initial condition bounding box $xa $xb $ya $yb $za $zb"; }else{ $cmd="#"; }
$cmd
#- # initial condition bounding box $xabb $xbbb $yabb $ybbb $zabb $zbbb
#- #  damp initial conditions at face (side,axis)=(0,1) of the box
#- bounding box decay face 1 0 
#- $betaBB=5.; # exponent in tanh function for smooth transition to zero outside the bounding box
#- bounding box decay exponent $betaBB
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
# ------- boundary conditions ----
#
# bc: all=dirichlet
# bc: all=perfectElectricalConductor
if( $rbc eq "pec" ){ $rbc="perfectElectricalConductor"; }
bc: all=$rbc
#
if( $leftBC eq "planeWave" ){ $cmd="bc: $backGround(0,0)=planeWaveBoundaryCondition"; }else{ $cmd="#"; }
$cmd 
if( $bcBody eq "pec" ){ $cmd="bc: annulus=perfectElectricalConductor"; }else{ $cmd="#"; }
$cmd
#
# User supplied BCs (array of string commands)
$cmd="#"; for( $i=0; $i<= $#bcCmds; $i++ ){ $cmd .= "\n" . $bcCmds[$i]; } 
# printf("cmd=$cmd\n");
# pause
$cmd
# 
# -- we need to subtract out the incident field on the "inflow" boundary before
#    applying the radiation boundary condition: 
if( $leftBC eq "planeWave" ){ $adjustFields=0; }else{ $adjustFields=1; }
adjust boundaries for incident field $adjustFields all
adjust boundaries for incident field $adjustFields $backGround
# 
#
# specify material properties for domains a file  (wdh: Feb 9, 2021)
# 
if( $materialFile ne "none" ){ $cmd="include $materialFile"; }else{ $cmd="#"; }
$cmd
#
# printf("materialFile=$materialFile\n");
# pause
# over-ride default materials with command line arguments:
# Periodic wrap the materials to set all bodies
if( $#matFileArray >=0 ){ @dmFile = (); $numMaterials =  $#matFileArray+1;  \
for( $i=0; $i<= $numBodies; $i++ ){ $mat = $i % $numMaterials; $dmFile[$i] = $matFileArray[$mat]; } \
}
# printf(" dmFile=@dmFile\n");
# pause
# 
# --- set material properties -----
#    $i=0 : set properties for $outerDomain (default - "outerDomain")
#    NOTE: material interfaces have share>=100
#
$cmd=""; 
for( $i=0; $i<=$numBodies; $i++ ){ \
  $cmd .= "GDM domain name: $domainName[$i]\n";  \
  $cmd .= "number of polarization vectors: $npv[$i]\n";  \
  $cmd .= "material file: $dmFile[$i]\n"; \
  }
  $cmd .= "#";
  # printf("cmd=\n$cmd\n");
$cmd
#
#
#  ------------- Set common options ----------
#
include setCommonOptions.h
#
# --- Probes ---
#
# output probes every this many time steps: 
probe frequency $probeFrequency
#
# ---- Optionally create left and right point probes ----
#
$leftProbeFileName="left$probeFileName.dat"; 
# 
if( $probeFileName ne "" && $probeFileName ne "none" ){ \
$cmd = "create a probe...\n" . \
       "  file name $leftProbeFileName\n" .  \
       "  probe name leftProbe\n" .  \
       "  nearest grid point to $xLeftProbe $yLeftProbe 0\n" . \
       "exit"; \
   }else{ $cmd="#"; }
$cmd
#
$rightProbeFileName="right$probeFileName.dat"; 
if( $probeFileName ne "" && $probeFileName ne "none" ){ \
$cmd = "create a probe...\n" . \
       "  file name $rightProbeFileName\n" . \
       "  probe name rightProbe\n" . \
       "  nearest grid point to $xRightProbe $yRightProbe 0\n" . \
       "exit"; \
   }else{ $cmd="#"; }
$cmd
#
# -- surface integral probe, "average" option means scale by surface area  --
if( $probeFileName ne ""  && $probeFileName ne "none" ){ $intProbeName = "integralLeft$probeFileName"; }
if( $intProbeName ne "" ){ $cmd="create a probe...\n   probe name $intProbeName\n  file name  $intProbeName.dat\n coordinate plane probe \n integral\n average\n  grid coordinate plane 0 $xLeftProbe $yLeftProbe 0 (axis, x,y,z)\n all components\n exit"; }else{ $cmd="#"; }
$cmd 
# 
if( $probeFileName ne ""  && $probeFileName ne "none" ){ $intProbeName = "integralRight$probeFileName"; }
if( $intProbeName ne "" ){ $cmd="create a probe...\n   probe name $intProbeName\n  file name  $intProbeName.dat\n coordinate plane probe \n integral\n average\n grid coordinate plane 0 $xRightProbe $yRightProbe 0 (axis, x,y,z)\n all components\n exit"; }else{ $cmd="#"; }
$cmd 
#
continue
#
#
# if( $az==0 ){ $cmd="plot:Ey"; }else{ $cmd="plot:Ez"; }
# $cmd 
# contour
#   plot contour lines (toggle)
#   # vertical scale factor 0.2
#   # min max -1.1 1.1
#   # plot a contour plane in 3d 
#   if( $grid =~ /3d/ || $grid =~ /Sphere/ || $grid =~ /Ellipsoid/ ){ $cmd="delete contour plane 2\n delete contour plane 1\n delete contour plane 0\n add contour plane  0.00000e+00  0.00000e+00  1.00000e+00 0 0 0"; }else{ $cmd="#"; }
#   $cmd
# exit
#
if( $extraCommandsEnd ne "" ){ $cmd="include $extraCommandsEnd"; }else{ $cmd="#"; }
$cmd
#
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
