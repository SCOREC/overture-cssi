#================================================================================================
#  cgmx example:  plane wave moving through a domain
#
# Usage:
#   
#  cgmx [-noplot] planeWave  -g=<name> -tf=<tFinal> -tp=<tPlot> -kx=<num> -ky=<num> -kz=<num> -show=<name> ...
#        -plotIntensity=[0|1] -eps1=<> -eps2=<> -interit=<> -diss=<> -filter=[0|1] -debug=<num> -cons=[0/1] ...
#        -method=[nfdtd|Yee|sosup|bamx] -bcn=[default|d|abc] -plotHarmonicComponents=[0|1] ] -dm=[none|gdm]
#        -useSosupDissipation=[0|1] -sosupParameter=[0-1] -sosupDissipationOption=[0|1] ...
#        -stageOption=[IDB|IBDB|D-IB|...] -a0= f f .. -a1=f f .. -b0= f f.. -b1=f f ..  -dmFile=<s> -dmSaveFile=<s> ...
#        -ts=[me|rk4] -solveForAllFields=[0|1] -drOption=[computeFrequency|computeWaveNumber]
#        -go=[run/halt/og]
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
#   cgmx planeWave -g=square10 -kx=1 -go=halt
#   cgmx planeWave -g=box10 -kx=1 -go=halt
#   cgmx planeWave -g=square10 -kx=2 -tp=.5 -tf=5. -plotIntensity=1 -intensityOption=1 -go=halt 
#   cgmx planeWave -g=square128.order4 -kx=4 -tp=.5 -tf=5. -plotIntensity=1 -intensityOption=1 -go=halt 
#
#  -- sosup:
#    cgmx planeWave -g=square32.order2 -kx=2 -ky=2 -tp=.1 -tf=5. -method=sosup -diss=0. -checkErrors=1 -go=halt  
#    cgmx planeWave -g=square40.order4.ng3 -kx=2 -ky=2 -tp=.1 -tf=5. -method=sosup -diss=0. -checkErrors=1 -go=halt [OK
#    cgmx planeWave -g=square40.order6.ng4 -kx=2 -ky=2 -tp=.1 -tf=5. -method=sosup -diss=0. -checkErrors=1 -go=halt
#
#  -- div clean:
#   cgmx planeWave -g=square40p -kx=2 -ky=1 -divClean=1 -divCleanCoeff=10. -tf=50 -tp=.5 -go=halt 
#   cgmx planeWave -g=square40p -kx=2 -ky=1 -divClean=1 -divCleanCoeff=20. -tf=50 -tp=.5 -go=halt 
#   cgmx planeWave -g=square80p -kx=2 -ky=1 -divClean=1 -divCleanCoeff=20. -tf=50 -tp=.5 -go=halt 
#   cgmx planeWave -g=nonSquare128p -kx=2 -ky=1 -divClean=1 -divCleanCoeff=10. -tf=50 -tp=.5 -go=halt 
# 
#   cgmx planeWave -g=square32p.order4 -kx=2 -ky=1 -divClean=1 -divCleanCoeff=10. -diss=0. -tf=50 -tp=.5 -go=halt
#   cgmx planeWave -g=square64p.order4 -kx=2 -ky=1 -divClean=1 -divCleanCoeff=20. -diss=0. -tf=50 -tp=.5 -go=halt 
# 
#    -- div clean and project interp sort of works: -- div does not decay very fast --
#   cgmx planeWave -g=sis8p -kx=2 -ky=1 -divClean=1 -divCleanCoeff=40. -projectInterp=1 -tf=50 -tp=.1 -go=halt
#     ... trouble with projectInterp=1 .. large errors 
#   cgmx planeWave -g=rsise4.order2p -kx=2 -ky=1 -divClean=1 -divCleanCoeff=5. -projectInterp=1 -filter=1 -tf=50 -tp=.1 -go=halt 
#     ... filter=1 looks good
#   cgmx planeWave -g=rsise4.order2p -kx=2 -ky=1 -divClean=1 -divCleanCoeff=5. -projectInterp=0 -diss=1. -tf=50 -tp=.1 -go=halt 
#   cgmx planeWave -g=rsise8.order2p -kx=2 -ky=1 -divClean=1 -divCleanCoeff=10. -projectInterp=0 -filter=1 -diss=0. -tf=50 -tp=.1 -go=halt 
#   cgmx planeWave -g=rsise4.order4p -kx=2 -ky=1 -divClean=1 -divCleanCoeff=10. -filter=1 -diss=0. -tf=50 -tp=.1 -go=halt 
# 
#      ... 3D compute E and H
#   cgmx planeWave -g=box40.order2p -ax=.5 -ay=-.5 -az=1. -kx=2 -ky=1 -kz=0 -divClean=1 -divCleanCoeff=10. -solveForH=1 -diss=0. -tf=50 -tp=.01 -go=halt
#   cgmx planeWave -g=nonBox40.order2p -ax=.5 -ay=-.5 -az=-.25 -kx=2 -ky=1 -kz=2 -divClean=1 -divCleanCoeff=10. -solveForH=1 -diss=0. -tf=50 -tp=.01 -go=halt
# 
# -- test Yee
#   cgmx planeWave -g=square10 -kx=1 -method=Yee -bcn=d -go=halt
#   cgmx noplot planeWave -g=square128 -kx=1 -ky=1 -method=Yee -bcn=d -go=og
#   cgmx planeWave -g=box64 -kx=1 -ky=1 -kz=1 -method=Yee -bcn=d -go=halt
#     -- abc: 
#   cgmx planeWave -g=square128 -kx=5 -method=Yee -bcn=abc -go=halt
#   cgmx noplot planeWave -g=box64 -kx=2 -ay=1. -method=Yee -bcn=abc -abcDir=0 -bg=box -tp=.01 -go=og
#   cgmx noplot planeWave -g=box64 -kx=-2 -ay=1. -method=Yee -bcn=abc -abcDir=0 -abcSide=0 -bg=box -tp=.01 -go=og
#   cgmx noplot planeWave -g=box64 -kx=0 -ky=2 -az=1. -method=Yee -bcn=abc -abcDir=1 -bg=box -tp=.01 -go=og
#   cgmx noplot planeWave -g=box64 -kx=0 -ky=-2 -az=1. -method=Yee -bcn=abc -abcDir=1 -abcSide=0 -bg=box -tp=.01 -go=og
#   cgmx noplot planeWave -g=box64 -kx=0 -kz=2 -ax=1. -method=Yee -bcn=abc -abcDir=2 -bg=box -tp=.01 -go=og
#   cgmx noplot planeWave -g=box64 -kx=0 -kz=-2 -ax=1. -method=Yee -bcn=abc -abcDir=2 -abcSide=0 -bg=box -tp=.01 -go=og
# 
#    -- Yee with an embedded material region: 
#   cgmx planeWave -g=bigSquare4.order2 -kx=2 -method=Yee -matRegion=1 -bcn=d -go=halt
# 
# -- test abc: 
#   cgmx planeWave -g=square64 -kx=2 -tp=.1 -tf=5. -bcn=abc -plotIntensity=1 -intensityOption=1 -go=halt 
#   cgmx noplot planeWave -g=box64 -kx=0. -kz=2. -ax=1. -ay=1. -tp=.1 -tf=5. -bcn=abc -plotIntensity=1 -intensityOption=1 -go=og
# 
#   cgmx planeWave -g=sise4.order2 -kx=1 -ky=1 -plotIntensity=1 -intensityOption=1 -go=halt 
#   cgmx planeWave -g=sise4.order4 -kx=1 -ky=1 -plotIntensity=1 -intensityOption=1 -go=halt 
#   cgmx planeWave -g=rsise4.order4 -kx=1 -ky=1 -divDamping=1. -go=halt
#
#   cgmx planeWave -g=cice3.order2 -kx=1 -ky=1 -plotIntensity=1 -intensityOption=1 -go=halt 
#   cgmx planeWave -g=cice3.order4 -kx=1 -ky=1 -plotIntensity=1 -intensityOption=1 -go=halt 
# 
#   cgmx planeWave -g=box16.order4 -kx=1 -ky=1 -kz=0 -plotIntensity=1 -intensityOption=1 -go=halt 
#   cgmx planeWave -g=box64.order4 -kx=1 -ky=1 -kz=0 -plotIntensity=1 -intensityOption=1 -go=halt 
#   cgmx planeWave -g=box16.order4 -kx=0 -ky=0 -kz=1 -ax=1. -plotIntensity=1 -intensityOption=1 -go=halt 
#   cgmx planeWave -g=bib2.order4 -kx=1 -ky=1 -kz=0 -plotIntensity=1 -intensityOption=1 -go=halt 
# 
# -- new 8th order dissipation
#    cgmx planeWave -g=sise4.order4 -kx=1 -ky=1 -diss=.0 -filter=1 -go=halt 
#      diss=0 has a weak instablity at about t=10: 
#    cgmx planeWave -g=rsise4.order4 -kx=1 -ky=1 -diss=.0 -filter=1 -tp=.5 -tf=50 -go=halt 
#      ok with    ad8 = -adc*dt/(256.*numberOfDimensions);  ad8*h^8 (D+D-)^4 
#    cgmx planeWave -g=rsise4.order4 -kx=1 -ky=1 -diss=0. -filter=1 -tp=.1 -tf=50 -go=halt -cfl=.9
#      -- ok with new filter
#      cgmx planeWave -g=rsise4.order4 -kx=1 -ky=1 -diss=0. -filter=1 -tp=.1 -tf=50 -go=halt -cfl=.9
#      cgmx planeWave -g=rsise8.order4 -kx=1 -ky=1 -diss=0. -filter=1 -tp=.1 -tf=50 -go=halt -cfl=.9
#      cgmx planeWave -g=rsise16.order4 -kx=2 -ky=2 -diss=0. -filter=1 -tp=.1 -tf=50 -go=halt -cfl=.9
# 
#      cgmx planeWave -g=cice3.order4 -kx=1 -ky=1 -diss=0. -filter=1 -tp=.1 -tf=50 -go=halt -cfl=.9
#      cgmx planeWave -g=cice6.order4 -kx=1 -ky=1 -diss=0. -filter=1 -tp=.1 -tf=50 -go=halt -cfl=.9
#   
#      cgmx planeWave -g=rbibe2.order4 -kx=2 -ky=2 -kz=2. -diss=0. -filter=1 -tp=.1 -tf=50 -go=halt
#        -- this starts to go unstable around t=15-20 with diss=0 : 
#      cgmx planeWave -g=rbibe4.order4 -kx=2 -ky=2 -kz=2. -diss=0. -filter=1 -tp=1. -tf=50. -go=halt
# 
#    srun -N1 -n2 -ppdebug memcheck_all $cgmxp noplot planeWave -g=square16.order4 -kx=1 -ky=1 -diss=.0 -filter=1 -go=go
#    totalview srun -N1 -n2 -ppdebug $cgmxp planeWave -g=square16.order4 -kx=1 -ky=1 -diss=.0 -filter=1 -go=halt 
#    mpirun -np 2 $cgmxp planeWave -g=sise4.order4 -kx=1 -ky=1 -diss=.0 -filter=1 -go=halt 
# ----
#      Intensity = .25*c*( eps a^2 + mu b^2 ) = .5*c*eps*a^2
#        k=(0,0,1), a=(0,0,1) -> b=(0,1,0)  I =.5
#        k=(0,0,1), a=(0,1,1) -> b=(-1,1,0) I =1.
#   cgmx planeWave -g=rbibe2.order2 -kx=0 -ky=0 -kz=1 -ax=1. -plotIntensity=1 -intensityOption=1 -go=halt  
#   cgmx planeWave -g=rbibe2.order4 -kx=0 -ky=0 -kz=1 -ax=1. -plotIntensity=1 -intensityOption=1 -go=halt 
# 
#   cgmx planeWave -g=rbibe2.order2 -eps=2.25 -kx=0 -ky=0 -kz=1 -ax=.81649658 -ay=.81649658 -plotIntensity=1 -intensityOption=1 -go=halt  
#
#  mpirun -np 2 $cgmxp noplot planeWave -g=square10 -kx=1 -show="planeWave.show" -tf=.1 -go=go
#================================================================================================
# 
$tFinal=2.; $tPlot=.1;  $show=" "; $method="NFDTD"; $bcn="default"; $matRegion=0; $bg="square";  $dm="none"; 
$cfl = .9; $diss=.5; $dissOrder=-1; $filter=0; $divClean=0; $divCleanCoeff=1; $solveForH=0; $projectInterp=0;
$ts="me"; $matFile=""; $solveForAllFields=0;
$dtMax=1e10; 
$debug=0; $divDamping=.0; $plotIntensity=0; $intensityOption=0; $abcDir=0; $abcSide=1; $plotHarmonicComponents=0; 
$cyl=1;   # set to 0 for a sphere 
$kx=2; $ky=0; $kz=0;
$drOption="computeFrequency"; # dispersion relation option 
$ax=0.; $ay=0.; $az=0.; # plane wave coeffs. all zero -> use default
$omega=-1.; # time harmonic omega for intensity computations 
$eps=1.; $mu=1.;
$show=" "; $backGround="backGround"; $useNewInterface=0; $checkErrors=1; 
$interfaceIterations=3;
$grid="innerOuter4.order4.hdf";
$cons=1; $go="halt";  $useSosupDissipation=0; $sosupParameter=1.;  $sosupDissipationOption=0;
$stageOption ="default";
# GDM parameters
$npv=1; $alphaP=-1.; $modeGDM=-1; 
@a0 = (); @a1=(); @b0=(); @b1=(); # these must be null for GetOptions to work, defaults are given below 
$dmFile=""; # "SilverJCDispersionFits.txt";
$dmSaveFile=""; # for saving GDM parameters to a file
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
    "projectInterp=i"=>\$projectInterp,"plotHarmonicComponents=i"=>\$plotHarmonicComponents,\
    "useSosupDissipation=i"=>\$useSosupDissipation,"sosupParameter=f"=>\$sosupParameter,\
    "sosupDissipationOption=i"=>\$sosupDissipationOption,"modeGDM=i"=>\$modeGDM,\
    "checkErrors=i"=>\$checkErrors,"dm=s"=>\$dm,"stageOption=s"=>\$stageOption,\
    "alphaP=f"=>\$alphaP,"a0=f{1,}"=>\@a0,"a1=f{1,}"=>\@a1,"b0=f{1,}"=>\@b0,"b1=f{1,}"=>\@b1,"npv=i"=>\$npv,\
    "dmFile=s"=>\$dmFile, "dmSaveFile=s"=>\$dmSaveFile,"matFile=s"=>\$matFile,"ts=s"=>\$ts,\
    "solveForAllFields=i"=>\$solveForAllFields, "drOption=s"=>\$drOption );
# -------------------------------------------------------------------------------------------------
# printf(" opt_a[0]=[%s] opt_a[1]=[%s]\n",$opt_a[0],$opt_a[1]);
# printf(" opt_b[0]=[%s] opt_b[1]=[%s]\n",$opt_b[0],$opt_b[1]);
if( $ts eq "me" ){ $ts="modifiedEquationTimeStepping"; }
$orderOfRungeKutta=4;
if( $ts eq "rk1" ){ $ts="rungeKutta"; $orderOfRungeKutta=1; }
if( $ts eq "rk2" ){ $ts="rungeKutta"; $orderOfRungeKutta=2; }
if( $ts eq "rk3" ){ $ts="rungeKutta"; $orderOfRungeKutta=3; }
if( $ts eq "rk4" ){ $ts="rungeKutta"; $orderOfRungeKutta=4; }
# printf(" TS = [$ts], order=[$orderOfRungeKutta] \n");
#
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
# time-stepping method
$ts
order of Runge Kutta $orderOfRungeKutta
#
solve for all fields $solveForAllFields
#
# dispersion model:
$dm
$drOption
# 
coefficients $eps $mu all (eps,mu,grid-name)
#
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
if( $npv > 2 ){ \
   $cmd  = "GDM domain name: $domain\n"; \
   $cmd .= " number of polarization vectors: $npv\n"; \
      }
$cmd
#
# -- read material parameters from a file 
if( $dmFile ne "" ) { $cmd="GDM domain name: $domain\n number of polarization vectors: $npv\n material file: $dmFile" }else{ $cmd="#"; }
$cmd 
# -- test: save material parameters to a file
if( $dmSaveFile ne "" ) { $cmd="GDM domain name: $domain\nsave material file: $dmSaveFile" }else{ $cmd="#"; }
$cmd 
# -- specify material regions or material file
$cmd="#"; 
if( $matFile ne "" ){ $cmd="material file: $matFile"; }else{ $cmd="#"; }
$cmd
$cmd="#"; 
# *********************** TEST ********************
# $numMatRegions=2; 
# if( $numMatRegions>1 ){\
# $matFile2="baAir.txt"; \
# $regionFile="sphereRegion.h"; $radius=1.; $xc=0; $yc=0; $zc=0; $ae=$radius; $be=$radius; $ce=$radius;\
# $cmd = "include $regionFile";\
# }
# $cmd
# 
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
planeWaveInitialCondition
planeWaveKnownSolution
#
#- gaussianPlaneWaveKnownSolution
#- gaussianPlaneWave
#- Gaussian plane wave: 50 0.5 0 0 (beta,x0,y0,z0)
# 
#- gaussianPulseInitialCondition
#- Gaussian pulse: 10 2 10 .5 .5 0 (beta,scale,exponent,x0,y0,z0)
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
bc: all=dirichlet
# ------------ ABC: ------------------------
# " bc: square(0,0)=planeWaveBoundaryCondition\n"
$abcSide2=1-$abcSide;
if( $bcn eq "abc" ){ \
$cmds = "bc: all=symmetry\n" .\
" bc: $bg($abcSide2,$abcDir)=dirichlet\n" .\
" bc: $bg($abcSide,$abcDir)=abcEM2\n" .\
" adjust boundaries for incident field 0 all\n" .\
" adjust boundaries for incident field 1 $bg\n"; }else{ $cmds ="#"; }
$cmds
# bc: $bg(1,0)=dirichlet
# ------------ ABC: ------------------------
#bc: all=abcEM2
#bc: box(0,2)=planeWaveBoundaryCondition
#bc: box(0,0)=symmetry
#bc: box(1,0)=symmetry
#bc: box(0,1)=symmetry
#bc: box(1,1)=symmetry
#adjust boundaries for incident field 0 all
#adjust boundaries for incident field 1 box
# -------------------------------------------
# 
# ++ bc: all=perfectElectricalConductor
# ++ bc: outerSquare(0,0)=planeWaveBoundaryCondition
#
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
dissipation $diss
apply filter $filter
order of dissipation $dissOrder
divergence damping $divDamping
project interpolation points $projectInterp
#
use sosup dissipation $useSosupDissipation
sosup parameter $sosupParameter
sosup dissipation option $sosupDissipationOption
#
#
use divergence cleaning $divClean
div cleaning coefficient $divCleanCoeff
#
use conservative difference $cons 
debug $debug
#
cfl $cfl 
dtMax $dtMax 
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
continue
#
plot:Ey
# plot:intensity
# 
$go
