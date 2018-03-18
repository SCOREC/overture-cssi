#================================================================================================
#  cgmx example:  plane wave moving through a domain
#
# Usage:
#   
#  cgmx [-noplot] gdm  -g=<name> -tf=<tFinal> -tp=<tPlot> -kx=<num> -ky=<num> -kz=<num> -show=<name> ...
#        -plotIntensity=[0|1] -eps1=<> -eps2=<> -interit=<> -diss=<> -filter=[0|1] -debug=<num> -cons=[0/1] ...
#        -method=[nfdtd|Yee|sosup] -bcn=[default|d|abc] -plotHarmonicComponents=[0|1] ] -dm=[none|gdm]
#        -useSosupDissipation=[0|1] -sosupParameter=[0-1] -sosupDissipationOption=[0|1] ...
#        -stageOption=[IDB|IBDB|D-IB|...] -a0= f f .. -a1=f f .. -b0= f f.. -b1=f f ..    -go=[run/halt/og]
#
# Arguments:
#  -kx= -ky= -kz= : integer wave numbers of the incident wave
#  -interit : number of iterations to solve the interface equations 
#  -filter=1 : add the high order filter
#  -plotHarmonicComponents : plot the harmonic components of the E field
#  -dm : dispersion model
#  -a0 : array values can be specified as -a0= 1 1.5 3.   (assign 3 array values)
#
echo to terminal 0
# Examples: 
#
#   cgmx gdm -g=square10 -kx=1 -go=halt
#
#================================================================================================
# 
$tFinal=2.; $tPlot=.1;  $show=" "; $method="NFDTD"; $bcn="default"; $matRegion=0; $bg="square";  $dm="none"; 
$cfl = .9; $diss=.5; $dissOrder=-1; $filter=0; $divClean=0; $divCleanCoeff=1; $solveForH=0; $projectInterp=0;
$debug=0; $divDamping=.0; $plotIntensity=0; $intensityOption=0; $abcDir=0; $abcSide=1; $plotHarmonicComponents=0; 
$cyl=1;   # set to 0 for a sphere 
$kx=2; $ky=0; $kz=0; 
$ax=0.; $ay=0.; $az=0.; # plane wave coeffs. all zero -> use default
$x0=.5; $y0=.5; $z0=0.; # Gaussian pulse center
$x1=-.5; $y1=-.5; $z1=0.; # Gaussian pulse center
$omega=-1.; # time harmonic omega for intensity computations 
$eps=1.; $mu=1.;
$show=" "; $backGround="backGround"; $useNewInterface=0; $checkErrors=0; 
$interfaceIterations=3;
$grid="innerOuter4.order4.hdf";
$cons=1; $go="halt";  $useSosupDissipation=0; $sosupParameter=1.;  $sosupDissipationOption=0;
$stageOption ="default";
$probeFileName="ProbeFile"; $xLeftProbe=-.5; $xRightProbe=-.25; $yLeftProbe=.25; $yRightProbe=.5;
# GDM parameters
$npv=1; $alphaP=1.; $modeGDM=-1; 
@a0 = (); @a1=(); @b0=(); @b1=(); # these must be null for GetOptions to work, defaults are given below 
# ----------------------------- get command line arguments ---------------------------------------
# GetOptions('a=s{2}' => \@opt_a, 'b=s{2}' => \@opt_b  );
# printf(" opt_a[0]=[%s] opt_a[1]=[%s]\n",$opt_a[0],$opt_a[1]);
# GetOptions("a0=f{2}"=>\@a0 );
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"diss=f"=>\$diss,"tp=f"=>\$tPlot,"show=s"=>\$show,"debug=i"=>\$debug, \
    "cfl=f"=>\$cfl, "bg=s"=>\$backGround,"bcn=s"=>\$bcn,"go=s"=>\$go,"noplot=s"=>\$noplot,"bg=s"=>\$bg,\
    "divDamping=f"=>\$divDamping,"plotIntensity=i"=>\$plotIntensity,"intensityOption=i"=>\$intensityOption,\
    "interit=i"=>\$interfaceIterations,"cyl=i"=>\$cyl,"useNewInterface=i"=>\$useNewInterface,\
    "dtMax=f"=>\$dtMax,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"ax=f"=>\$ax,"ay=f"=>\$ay,"az=f"=>\$az,\
    "eps=f"=>\$eps,"mu=f"=>\$mu, "cons=i"=>\$cons,"method=s"=>\$method,"bcn=s"=>\$bcn,"matRegion=i"=>\$matRegion,\
    "abcDir=i"=>\$abcDir,"abcSide=i"=>\$abcSide,"dissOrder=i"=>\$dissOrder,"filter=i"=>\$filter,\
    "divClean=i"=>\$divClean,"divCleanCoeff=f"=>\$divCleanCoeff,"solveForH=i"=>\$solveForH,\
    "probeFileName=s"=>\$probeFileName,"xLeftProbe=f"=>\$xLeftProbe,"xRightProbe=f"=>\$xRightProbe,\
    "yLeftProbe=f"=>\$yLeftProbe,"yRightProbe=f"=>\$yRightProbe,\
    "projectInterp=i"=>\$projectInterp,"plotHarmonicComponents=i"=>\$plotHarmonicComponents,\
    "useSosupDissipation=i"=>\$useSosupDissipation,"sosupParameter=f"=>\$sosupParameter,\
    "sosupDissipationOption=i"=>\$sosupDissipationOption,"modeGDM=i"=>\$modeGDM,\
    "checkErrors=i"=>\$checkErrors,"dm=s"=>\$dm,"stageOption=s"=>\$stageOption,\
    "alphaP=f"=>\$alphaP,"a0=f{1,}"=>\@a0,"a1=f{1,}"=>\@a1,"b0=f{1,}"=>\@b0,"b1=f{1,}"=>\@b1,"npv=i"=>\$npv,\
    "x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"x1=f"=>\$x1,"y1=f"=>\$y1,"z1=f"=>\$z1  );
# -------------------------------------------------------------------------------------------------
# printf(" opt_a[0]=[%s] opt_a[1]=[%s]\n",$opt_a[0],$opt_a[1]);
# printf(" opt_b[0]=[%s] opt_b[1]=[%s]\n",$opt_b[0],$opt_b[1]);
# Give defaults here for array arguments: 
if( $a0[0] eq "" ){ @a0=(1,0,0,0); }
if( $a1[0] eq "" ){ @a1=(0,0,0,0); }
if( $b0[0] eq "" ){ @b0=(0,0,0,0); }
if( $b1[0] eq "" ){ @b1=(0,0,0,0); }
printf(" a0[0]=%f, a0[1]=%f\n",$a0[0],$a0[1]);
printf(" b1[0]=%f, b1[1]=%f\n",$b1[0],$b1[1]);
#
if( $dm eq "none" ){ $dm="no dispersion"; }
if( $dm eq"drude" || $dm eq "Drude" ){ $dm="Drude"; }
if( $dm eq"gdm" ){ $dm="GDM"; }
#
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
# 
echo to terminal 1
#
$grid
#
$method
# dispersion model:
$dm
# 
# Drude params 1 1 all (gamma,omegap,domain-name)
#GDM params $a0 $a1 $b0 $b1 all (a0,a1,b0,b1,domain-name)
GDM mode: $modeGDM
$domain="all"; 
$cmd="#"; 
if( $npv == 1 ){ $cmd = "GDM params $a0[0] $a1[0] $b0[0] $b1[0] all (a0,a1,b0,b1,domain-name)"; }
if( $npv == 2 ){ \
   $cmd  = "GDM domain name: $domain\n"; \
   $cmd .= " number of polarization vectors: $npv\n"; \
   $cmd .= " GDM coeff: 0 $a0[0] $a1[0] $b0[0] $b1[0] (eqn, a0,a1,b0,b1)\n"; \
   $cmd .= " GDM coeff: 1 $a0[1] $a1[1] $b0[1] $b1[1] (eqn, a0,a1,b0,b1)"; \
      }
$cmd
## planeWaveInitialCondition
## planeWaveKnownSolution
#
#- gaussianPlaneWaveKnownSolution
#- gaussianPlaneWave
#- Gaussian plane wave: 50 0.5 0 0 (beta,x0,y0,z0)
# 
#               " g(x,y,t) = a*cos(2*pi*omega*(t-t0) )*exp( -beta*[ (x-x0)^2 + (y-y0)^2 ]^p )\n"
#	       " F(Ex) = -(y-y0)*g(x,y,t) \n"
#	       " F(Ey) =  (x-x0)*g(x,y,t) \n"
#	       " F(Hz) =         g(x,y,t) \n"
userDefinedInitialConditions
gaussian pulses
  2
  $amp0=50.; $beta0=40; $p0=1.; 
  $amp1=50.; $beta1=40; $p1=1.; 
  #  a,beta,omega,p,x0,y0,z0,t0
   $amp0 $beta0 0. $p0 $x0 $y0 $z0 0.
   $amp1 $beta1 0. $p1 $x1 $y1 $z1 0.
exit
#-- 
#-- gaussianPulseInitialCondition
#--    # H ??  scale*exp( -pow(beta*(xe*xe+ye*ye),exponent) );
#--   $scale=15; $beta=50.; $exponent=2; 
#--   Gaussian pulse: $beta $scale $exponent $x0 $y0 $z0 (beta,scale,exponent,x0,y0,z0)
#
# ++ zeroInitialCondition
# ====
# twilightZone
#  degreeSpace, degreeTime  1 1
#
kx,ky,kz $kx $ky $kz
plane wave coefficients $ax $ay $az $eps $mu
# 
# *****************
# *****************
# for Yee we define the cylinder as a masked stair step region
$rad=.5; $x0=0.; $y0=0; $z0=0; $eps1=2.; $mu1=1.; 
if( $method eq "Yee" && $matRegion==1 ){ $cmds = "define embedded bodies\n dielectric cylinder\n $rad $x0 $y0 $z0\n $eps1 $mu1 0. 0. \nexit"; }else{ $cmds="#"; }
$cmds 
# ****************
#
bc: all=perfectElectricalConductor
# --------------- fix me --------------
if( $bcn eq "d" ){ $bcn = "bc: all=dirichlet"; }
if( $bcn eq "s" ){ $bcn = "bc: all=symmetry"; }else{ $bcn = "#"; }
$bcn
# 
#
# bc: Annulus=perfectElectricalConductor
tFinal $tFinal
tPlot  $tPlot
#
solve for magnetic field $solveForH
#
coefficients $eps $mu all (eps,mu,grid-name)
#
dissipation $diss
apply filter $filter
order of dissipation $dissOrder
divergence damping $divDamping
project interpolation points $projectInterp
#
use sosup dissipation $useSosupDissipation
# sosup parameter $sosupParameter
# sosup dissipation option $sosupDissipationOption
#
#
use divergence cleaning $divClean
div cleaning coefficient $divCleanCoeff
#
use conservative difference $cons 
debug $debug
#
cfl $cfl 
plot errors $checkErrors
check errors $checkErrors
plot intensity $plotIntensity
intensity option $intensityOption
plot harmonic E field $plotHarmonicComponents
# $c=1./sqrt($eps*$mu);
# $omega= $c*sqrt( $kx*$kx + $ky*$ky + $kz*$kz );
# time harmonic omega $omega (omega/(2pi), normally c*|k|
# 
#*********************************
show file options...
  MXSF:compressed
  MXSF:open
  $show
  MXSF:frequency to flush 5
exit
#**********************************
create user defined probe...
 $leftProbeName="leftProbe";
 $leftProbeFileName="left$probeFileName.dat";
 probe name $leftProbeName
 file name $leftProbeFileName
 grid point probe
  $xLeftProbe $yLeftProbe 0
exit
# 
create user defined probe...
 $rightProbeName="rightProbe";
 $rightProbeFileName="right$probeFileName.dat";
 probe name $rightProbeName
 file name $rightProbeFileName
 grid point probe
  $xRightProbe $yRightProbe 0
exit
# # Point probe: 
# create a probe...
#   $leftProbeFileName="left$probeFileName.dat"; 
#   file name $leftProbeFileName
#   probe name leftProbe
#   nearest grid point to $xLeftProbe $yLeftProbe 0
# exit
# create a probe...
#   $rightProbeFileName="right$probeFileName.dat"; 
#   file name $rightProbeFileName
#   probe name rightProbe
#   nearest grid point to $xRightProbe $yRightProbe 0
# exit
# 
continue
#
plot:Ey
# plot:intensity
# 
$go
