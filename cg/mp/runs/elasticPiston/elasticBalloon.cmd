#
# cgmp:   INS + Elasticity: elastic "balloon" problem
# 
# Usage:
#    cgmp [-noplot] elasticBalloon -g=<name> -method=[ins|cns] -nu=<> -mu=<> -kappa=<num> -tf=<tFinal> -tp=<tPlot> ...
#           -solver=[yale|best] -psolver=[yale|best] -ktcFluid=<> -ktcFluid=<> -tz=[poly/trig/none] -bg=<backGroundGrid> ...
#           -degreeSpace=<num> -degreeTime=<num> -ts=[fe|be|im|pc] -nc=[] -d1=<> -d2=<> -smVariation=[nc|c|g|h] ...
#           -known=[piston|shear] -sideBC=[noSlipWall|slipWall|dirichlet] -useImplicitAmpBCs=[0|1]
# 
#  -ktcFluid -ktcSolid : thermal conductivities 
#  -ts = time-stepping-method, be=backward-Euler, fe=forward-Euler, im=implicit-multistep
#  -d1, -d2 : names for domains
# 
# Examples:
# 
# --- set default values for parameters ---
# 
$grid="deformingChannelGrid4.order2"; $domain1="fluidDomain"; $domain2="solidDomain";
$method="ins"; $probeFile="probeFile"; $multiDomainAlgorithm=1;  $pi=0; $pOffset=0.; 
$tFinal=20.; $tPlot=.1;  $cfl=.9; $show="";  $pdebug=0; $debug=0; $go="halt"; $cdv=""; $cDt=""; 
$muFluid=0.; $rhoFluid=1.4; $pFluid=1.; $TFluid=$pFluid/$rhoFluid; 
$nu=.1; $rhoSolid=1.; $prandtl=.72; $cnsVariation="jameson"; $ktcFluid=-1.; $u0=0.; $xShock=-1.5; $uShock=1.25; 
$p0=1.; 
$cnsEOS="ideal"; 
$cnsGammaStiff=1.4; $cnsPStiff=0.;   # for stiffened EOS -- by default make it look like an ideal gas
$lambdaSolid=1.; $muSolid=1.;
$thetad=0.; # rotation of domain (degrees)
## $stressRelaxation=1; $relaxAlpha=0.1; $relaxDelta=0.1; 
$stressRelaxation=4; $relaxAlpha=.5; $relaxDelta=.5; 
# $stressRelaxation=0;  # *wdh* turn off stress-relaxation for testing TZ
## $displacementDissipation=0.; # for CgSm -- turn off for TZ
$scf=1.; # solidScaleFactor : scale rho,mu and lambda by this amount 
$thermalExpansivity=1.; $T0=1.; $Twall=1.;  $kappa=.01; $ktcSolid=-1.; 
$diss=.5;   # 2nd-order linear dissipation for cgsm 
## $diss=.0;   # TURN OFF FOR TZ *wdh* April 19, 2018
$smVariation = "g"; 
$setGhostByExtrapolation=0; 
$tsSM="modifiedEquationTimeStepping";
$tz="none"; $degreeSpace=1; $degreeTime=1;
$degreeSpaceSM=""; $degreeTimeSM="";  # if set, use this as the degree for cgsm
$gravity = "0 0. 0."; $boundaryPressureOffset=0.; $cnsGodunovOrder=2; 
$fic = "uniform";  # fluid initial condition
$backGround="outerSquare"; $deformingGrid="interface"; 
#
$ts="pc";   # MP solver
$tsINS="pc"; # INS time-stepping method 
$numberOfCorrections=1;  # cgmp and cgins 
$coupled=0; $iTol=1.e-3; $iOmega=1.; $flushFrequency=10; $useNewInterfaceTransfer=0; 
$useTP=0; # 1=use traditional partitioned scheme with sub-iterations
$tol=1.e-3; $omega=.5; # for sub-iterations
$projectMultiDomainInitialConditions=0; 
$useNewTimeSteppingStartup=1;  # *NEW* July 1, 2017
$freqFullUpdate=1; # frequency for using full ogen update in moving grids 
$useExactPressureBC=0; # use exact RHS for pressure wall BC (variable used in insDomain.h)
$useCurlFormOfTraction=0; # use div(v)=0 to alter viscous traction
#
$smoothInterface=0;  # smooth the interface (in DeformingBodyMotion.C )
$numberOfInterfaceSmooths=4; 
$startCurve="nurbs"; # start curve for the deforming body 
#
# $option="beamUnderPressure"; # this currently means ramp the inflow
$option="bulkSolidPiston"; # define pressure BC from known solution
$known="piston"; #  or "shear"
#
$sideBC="slipWall"; 
#
$bcOption=4;   # does this do anything ? I thibnk this is for cgcns
$orderOfExtrapForOutflow=3; $orderOfExtrapForGhost2=2; $orderOfExtrapForInterpNeighbours=2; 
$projectInitialConditions=0; # for INS
# 
$psolver="yale"; 
$solver="yale"; 
$ksp="bcgs"; $pc="bjacobi"; $subksp="preonly"; $subpc="ilu"; $iluLevels=3;
# -- p-wave strength: don't make too big or else solid may become inverted in the deformed space
$append=0; 
#
#  Normal force = sin(Pi*t/timeInterval)^2 * pAmp * sin(freq*theta) 
$pAmp=.05; $freq=3.; $timeInterval=.25; 
# ------------------------- turn on added mass here ----------------
$addedMass=0;  $addedMassVelocityBC=0; $zfMuByH=5.; $zfRhoHByDt=0.; 
$addedMassLengthScale=1.; # default is 1 
$useImplicitAmpBCs=0; # set to 1 to use new implicit AMP BC's -- do this for now, make default later
# $predictedBoundaryPressureNeeded=1; # predict pressure for velocity BC *wdh* Dec 25, 2017
$predictedBoundaryPressureNeeded=0; # WDH: WHY IS THIS NEEDED ?? TURN OFF FOR NOW - April 19, 2018
# ---- piston parameters:  choose t0=1/(4*k) to make yI(0)=0 
$Pi=4.*atan2(1.,1.);
$amp=.1; $k=.5; $t0=1./(4*$k);  $H=1.; $Hbar=.5; $rho=1.; 
$rampOrder=2;  # number of zero derivatives at start and end of the ramp
$ra=-10.; $rb=-9.; # ramp interval -- actual interval shifted by Hbar/cp 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "g=s"=>\$grid,"tf=f"=>\$tFinal,"nu=f"=>\$nu,"muFluid=f"=>\$muFluid,"kappa=f"=>\$kappa, "bg=s"=>\$backGround,\
 "tp=f"=>\$tPlot, "solver=s"=>\$solver, "psolver=s"=>\$psolver,"useTP=i"=> \$useTP,\
 "tz=s"=>\$tz,"degreeSpace=i"=>\$degreeSpace, "degreeTime=i"=>\$degreeTime,\
 "degreeSpaceSM=i"=>\$degreeSpaceSM,"degreeTimeSM=i"=>\$degreeTimeSM,\
 "show=s"=>\$show,"method=s"=>\$method,"ts=s"=>\$ts,"tsSM=s"=>\$tsSM,"noplot=s"=>\$noplot,"ktcFluid=f"=>\$ktcFluid,\
  "ktcSolid=f"=>\$ktcSolid,"muSolid=f"=>\$muSolid,"lambdaSolid=f"=>\$lambdaSolid, "T0=f"=>\$T0,"Twall=f"=>\$Twall,\
  "nc=i"=> \$numberOfCorrections, "numberOfCorrections=i"=> \$numberOfCorrections,"coupled=i"=>\$coupled,\
  "d1=s"=>\$domain1,"d2=s"=>\$domain2,"dg=s"=>\$deformingGrid,"debug=i"=>\$debug,"kThermalFluid=f"=>\$kThermalFluid,\
  "cfl=f"=>\$cfl,"rhoSolid=f"=>\$rhoSolid,"cnsVariation=s"=>\$cnsVariation,"diss=f"=>\$diss,"fic=s"=>\$fic,"go=s"=>\$go,\
   "smVariation=s"=>\$smVariation,"scf=f"=>\$scf,"probeFile=s"=>\$probeFile,"pOffset=f"=>\$boundaryPressureOffset,\
   "cnsGodunovOrder=f"=>\$cnsGodunovOrder,"flushFrequency=i"=>\$flushFrequency,\
   "cnsEOS=s"=>\$cnsEOS,"cnsGammaStiff=f"=>\$cnsGammaStiff,"cnsPStiff=f"=>\$cnsPStiff,"u0=f"=>\$u0,\
   "useNewInterfaceTransfer=i"=>\$useNewInterfaceTransfer,"multiDomainAlgorithm=i"=>\$multiDomainAlgorithm,\
   "pi=i"=>\$pi,"xShock=f"=>\$xShock,"uShock=f"=>\$uShock,"bcOption=i"=>\$bcOption,"option=s"=>\$option,\
   "stressRelaxation=f"=>\$stressRelaxation,"relaxAlpha=f"=>\$relaxAlpha,"relaxDelta=f"=>\$relaxDelta,\
   "displacementDissipation=f"=>\$displacementDissipation,\
   "p0=f"=>\$p0,"sideBC=s"=>\$sideBC,"iOmega=f"=>\$iOmega,"iTol=f"=>\$iTol,"addedMass=f"=>\$addedMass,\
   "projectInitialConditions=f"=>\$projectInitialConditions,"restart=s"=>\$restart,"append=i"=>\$append,\
   "projectMultiDomainInitialConditions=f"=>\$projectMultiDomainInitialConditions,"known=s"=>\$known,\
   "amp=f"=>\$amp,"rampOrder=i"=>\$rampOrder,"ra=f"=>\$ra,"rb=f"=>\$rb,"cdv=f"=>\$cdv,"cDt=f"=>\$cDt,\
   "omega=f"=>\$omega,"tol=f"=>\$tol,"zfMuByH=f"=>\$zfMuByH,"zfRhoHByDt=f"=>\$zfRhoHByDt,\
   "useNewTimeSteppingStartup=i"=> \$useNewTimeSteppingStartup,"tsINS=s"=>\$tsINS,"addedMassVelocityBC=i"=>\$addedMassVelocityBC,\
   "freqFullUpdate=i"=>\$freqFullUpdate,"smoothInterface=i"=>\$smoothInterface,"addedMassLengthScale=f"=>\$addedMassLengthScale,\
   "numberOfInterfaceSmooths=i"=>\$numberOfInterfaceSmooths,"useImplicitAmpBCs=i"=>\$useImplicitAmpBCs,\
   "useCurlFormOfTraction=i"=>\$useCurlFormOfTraction,"dtMax=f"=>\$dtMax,"thetad=f"=>\$thetad,\
   "useExactPressureBC=i"=>\$useExactPressureBC,"startCurve=s"=>\$startCurve,\
   "setGhostByExtrapolation=i"=>\$setGhostByExtrapolation,\
   "pAmp=f"=>\$pAmp,"freq=f"=>\$freq,"timeInterval=f"=>\$timeInterval );
# -------------------------------------------------------------------------------------------------
if( $solver eq "best" ){ $solver="choose best iterative solver"; }
if( $psolver eq "best" ){ $psolver="choose best iterative solver"; }
if( $ts eq "fe" ){ $ts="forward Euler"; }
if( $ts eq "be" ){ $ts="backward Euler"; }
if( $ts eq "im" ){ $ts="implicit"; }
if( $ts eq "pc" ){ $ts="adams PC"; }
#
if( $tsINS eq "fe" ){ $tsINS="forward Euler";}
if( $tsINS eq "be" ){ $tsINS="backward Euler"; }
if( $tsINS eq "im" ){ $tsINS="implicit"; }
if( $tsINS eq "bdf" ){ $tsINS="implicit BDF"; }
if( $tsINS eq "imex" ){ $tsINS="implicit explicit multistep"; }
if( $tsINS eq "pc" ){ $tsINS="adams PC"; }
if( $tsINS eq "pc4" ){ $tsINS="adams PC order 4"; $useNewImp=0; } # NOTE: turn off new implicit for fourth order
if( $tsINS eq "mid"){ $tsINS="midpoint"; }  
if( $tsINS eq "afs"){ $tsINS="approximate factorization"; $newts=1;  $implicitVariation="full"; }
#
if( $tz eq "none" ){ $tz="turn off twilight zone"; }
if( $tz eq "poly" ){ $tz="turn on twilight zone\n turn on polynomial"; $cdv=0.; }
if( $tz eq "trig" ){ $tz="turn on twilight zone\n turn on trigonometric"; $cdv=0.; }
if( $go eq "halt" ){ $go = "break"; }
if( $go eq "og" ){ $go = "open graphics"; }
if( $go eq "run" || $go eq "go" ){ $go = "movie mode\n finish"; }
#
if( $projectInitialConditions eq "1" ){ $projectInitialConditions = "project initial conditions"; }else{ $projectInitialConditions = "do not project initial conditions"; }
# 
if( $smVariation eq "nc" ){ $smVariation = "non-conservative"; }
if( $smVariation eq "c" ){ $smVariation = "conservative"; $cons=1; }
if( $smVariation eq "g" ){ $smVariation = "godunov"; }
if( $smVariation eq "h" ){ $smVariation = "hemp"; }
#
if( $method eq "ins" && $kThermalFluid eq "" ){ $kThermalFluid=$nu/$prandtl; }
if( $method eq "cns" && $kThermalFluid eq "" ){ $kThermalFluid=$muFluid/$prandtl; }
if( $ktcFluid < 0. ){ $ktcFluid=$kThermalFluid;} if( $ktcSolid < 0. ){ $ktcSolid=$kappa; }
# 
$grid
# ----------  define deforming bodies by a share flag of 100 ----
# ----------  NOTE: we parameterize the boundary by index so grid points match! ---
# **** NEW WAY TO SPECIFY DEFORMING BODY FOR A BULK SOLID 
# $vInitial=-.54414;
$numberOfPastTimeLevels=3; 
$gridEvolutionVelocityAccuracy=3; 
$gridEvolutionAccelerationAccuracy=2; 
# choose the start curve type:
$startCurveCmd="#";
if( $startCurve eq "nurbs"  ){ $startCurveCmd="nurbs start curve"; }
if( $startCurve eq "spline" ){ $startCurveCmd="spline start curve"; }
if( $tz eq "turn off twilight zone" ){ $useKnown=1; }else{ $useKnown=0; }
$moveCmds = \
  "turn on moving grids\n" . \
  "specify grids to move\n" . \
  "    deforming body\n" . \
  "      $startCurveCmd\n" . \
  "      bulk solid\n" . \
  "        debug\n $debug \n" . \
  "      velocity order of accuracy\n $gridEvolutionVelocityAccuracy\n" . \
  "      acceleration order of accuracy\n $gridEvolutionAccelerationAccuracy\n" . \
  "      number of past time levels: $numberOfPastTimeLevels\n" . \
  "      smooth surface $smoothInterface \n" . \
  "      number of surface smooths: $numberOfInterfaceSmooths \n" . \
  "      relax correction steps $useTP \n" . \
  "      sub iteration convergence tolerance\n $tol \n" . \
  "      added mass relaxation factor\n $omega\n" . \
  "     done\n" . \
  "     choose grids by share flag\n" . \
  "        100 \n" . \
  "   done\n" . \
  "done";
# 
#$probeFileName = $probeFile . "Fluid.dat";
#$extraCmds = \
#    "frequency to save probes 1\n" . \
#    "create a probe\n" . \
#    "  file name $probeFileName\n" . \
#    "  nearest grid point to 0. .5 0.\n" . \
#    "  exit";
# ------- specify fluid domain ----------
#  Cgins:
$domainName=$domain1; $solverName="fluid"; 
$modelNameINS="none"; 
#
$T0=0.; 
$bc = "all=$sideBC\n bcNumber100=noSlipWall uniform(u=.0,T=$T0)\n bcNumber100=tractionInterface";
#
$ic="uniform flow\n" . "p=0., u=$u0, T=$T0";
$rhoBar=$rhoSolid*$scf; $lambdaBar=$lambdaSolid*$scf; $muBar=$muSolid*$scf;
$thetaR=$thetad*$Pi/180.;
$knownCmds ="#";
if( $tz ne "turn off twilight zone" ){ $ic="#"; }
# For sub-time-step iterations: 
$numberOfTimeStepCorrections=$numberOfCorrections; 
#
echo to terminal 0
# $extraCmds="check error on ghost\n 1";
include $ENV{CG}/mp/cmd/insDomain.h
$extraCmds="#"; 
echo to terminal 1
# 
# ------- specify elastic solid domain ----------
$domainName=$domain2; $solverName="solid"; 
#  Normal force = sin(Pi*t/timeInterval)^2 * pAmp * sin(freq*theta) 
$bcCommands="bcNumber4=tractionBC, userDefinedBoundaryData\n sinusoidal pressure force\n $pAmp $freq $timeInterval\n done\n bcNumber100=tractionBC\n bcNumber100=tractionInterface"; 
#
$initialConditionCommands="zeroInitialCondition";
#
$smCheckErrors=0;
# 
echo to terminal 0
if( $degreeSpaceSM ne "" ){ $degreeSpace=$degreeSpaceSM; }
if( $degreeTimeSM ne "" ){ $degreeTime=$degreeTimeSM; }
include $ENV{CG}/mp/cmd/smDomain.h
echo to terminal 1
# 
continue
#
# -- set parameters for cgmp ---
# 
  final time $tFinal
  times to plot $tPlot
  cfl $cfl
  dtMax $dtMax
  $ts
  number of PC corrections $numberOfCorrections
  OBPDE:interface tolerance $iTol
  OBPDE:interface omega $iOmega
  OBPDE:solve coupled interface equations $coupled
  # project multi-domain initial condtions
  OBPDE:project initial conditions $projectMultiDomainInitialConditions
  OBPDE:use new interface transfer $useNewInterfaceTransfer
  # relax correction steps for TP scheme: 
  OBPDE:relax correction steps $useTP
 # -- for testing solve the domains in reverse order: 
 # OBPDE:domain order 1 0
  OBPDE:project interface $pi
  if( $multiDomainAlgorithm eq 1 ){ $cmd="OBPDE:step all then match advance"; }else{ $cmd="#"; }
  $cmd 
  #
  $tz
  # DEFINE THE MULTI_STAGE ALGORITHM --
  OBPDE:multi-stage
  actions=takeStep classNames=Cgsm
  actions=takeStep,applyBC classNames=Cgins
  actions=applyBC classNames=Cgsm
  #actions=takeStep classNames=Cgsm,Cgins
  #actions=takeStep,applyBC classNames=Cgins
  #actions=takeStep,applyBC domainNames=fluidDomain
  #actions=takeStep,applyBC domainNames=fluidDomain,solidDomain
  done
#
  debug $debug
  show file options
    if( $append eq 0 ){ $cmd="OBPSF:create new show file"; }else{ $cmd="OBPSF:append to old show file"; }
    $cmd
    compressed
      open
       $show
    frequency to flush
      $flushFrequency
    exit
  continue
#
continue
# --
        erase
        plot domain: fluid
        contour
          vertical scale factor 0.
           # ghost lines 1
           if( $known eq "shear" ){ $cmd="plot:u"; }else{ $cmd="plot:p"; }
           $cmd
           ##  wire frame
          exit
        plot domain: solid
        contour
          vertical scale factor 0.
          adjust grid for displacement 1
        exit
        if( $known eq "shear" ){ $cmd="plot:solid : v1"; }else{ $cmd="plot:solid : v2"; }
        $cmd
        plot all
$go


   erase
   plot domain: fluid
   grid
     exit this menu
   plot domain: solid
   displacement
     exit this menu
   plot all




          OBIC:user defined...
            bubble with shock
            r=1 T=1 u=0 v=0
            .2 .5 0.
            r=2. T=2. 
            -5.
            r=2 T=2
            exit

