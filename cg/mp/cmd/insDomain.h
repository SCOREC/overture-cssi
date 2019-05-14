# -*- mode: text; -*-
# file is included by other command files to define a Cgins domain
#
# Here are optional parameters that can be set before including this file: 
#   $domainName [required] : name of the domain to assign
#   $modelNameINS : "none", "Boussinesq"
#   $solverName : name given to the domain (e.g. "fluid")
#   $tz, $degreeSpace, $degreeTime : 
#   $tsINS -- time-stepping method 
#   $dtMax
#   $nu, $kThermal, $gravity, $thermalExpansivity, $ktc
#   $ic : specify initial condition commands
#   $bc : specify boundary condition commands
#   $gravity : 
#   $solver, $atoli, $rtoli : implicit solver
#   $psolver, $atolp, $rtolp : pressure solver
#   $pdebug 
#   $axisymmetric : if non-null turn on the axisymmetric flow option
#   $implicitFactor : for implicit time-stepping: .5=CN, 1.=BE, 0.=FE
#   $predictedBoundaryPressureNeeded : some BC's need a predicted pressure on the boundary
#   $projectInitialConditions : 
#   $numberOfTimeStepCorrections : numnber of correction steps for predictor-corrector time-stepping
#   $extraCmds : extra commands
#   $checkForInflowAtOutflow = [0|1]   
#   $useNeumannAtOutflow = [0|1]
#   $useNewTimeSteppingStartup = [0|1]
#   $orderOfExtrapForOutflow
#   $addedMassVelocityBC=[0|1]
#
# ------- start new domain ----------
if( $tsINS eq "" ){ $tsINS=$ts; }
if( $thermalExpansivity eq "" ){ $thermalExpansivity=1.; }
if( $modelNameINS eq "" ){ $modelNameINS="Boussinesq model"; }elsif( $modelNameINS eq "none" ){ $modelNameINS="#"; } #
if( $adcBoussinesq eq "" ){ $adcBoussinesq=0.; }
if( $dtMax eq "" ){ $dtMax=1.e20; }
if( $solver eq "" ){ $solver="choose best iterative solver"; }
if( $psolver eq "" ){ $psolver=$solver; }
if( $rtolp eq "" ){ $rtolp=1.e-5; }
if( $atolp eq "" ){ $atolp=1.e-7; }
if( $rtoli eq "" ){ $rtoli=1.e-5; }
if( $atoli eq "" ){ $atoli=1.e-7; }
if( $pdebug eq "" ){ $pdebug=0; }
if( $idebug eq "" ){ $idebug=0; }
if( $tz eq "" ){ $tz="turn off twilight zone"; }
if( $fx eq "" ){ $fx=1.; }
if( $fy eq "" ){ $fy=1.; }
if( $fz eq "" ){ $fz=1.; }
if( $ft eq "" ){ $ft=1.; }
if( $u0 eq "" ){ $u0=0.; }
if( $T0 eq "" ){ $T0=0.; }
if( $cdv eq "" ){ $cdv=1.; }
if( $cDt eq "" ){ $cDt=.25; }  # the div damping is at most cDt/dt
if( $ad2 eq "" ){ $ad2=0; }
if( $ad21 eq "" ){ $ad21=1.; }
if( $ad22 eq "" ){ $ad22=1.; }
if( $implicitFactor eq "" ){ $implicitFactor=.5; }
if( $implicitVariation eq "" ){ $implicitVariation = "implicitViscous"; }
if( $predictedPressureNeeded eq "" ){ $predictedPressureNeeded=0; } 
if( $predictedBoundaryPressureNeeded eq "" ){ $predictedBoundaryPressureNeeded=0; } 
if( $projectInitialConditions eq "" ){ $projectInitialConditions="#";} 
if( $extraCmds eq "" ){ $extraCmds="#"; }
if( $moveCmds eq "" ){ $moveCmds="#"; }
if( $useNewTimeSteppingStartup eq "" ){ $useNewTimeSteppingStartup=0; }
if( $checkForInflowAtOutflow eq "" ){ $checkForInflowAtOutflow=0; }
if( $useNeumannAtOutflow eq "" ){ $useNeumannAtOutflow=0; }
if( $fluidDensity eq "" ){ $fluidDensity=1.; }
if( $useImplicitAmpBCs eq "" ){ $useImplicitAmpBCs=0; }
if( $ogesDtol eq "" ){ $ogesDtol=1e20; }
if( $orderOfExtrapForOutflow eq "" ){ $orderOfExtrapForOutflow=2; }
if( $useExactPressureBC eq "" ){ $useExactPressureBC=0; }
if( $useCurlFormOfTraction eq "" ){ $useCurlFormOfTraction=0; }
if( $addedMassLengthScale eq "" ){ $addedMassLengthScale=1.; }
if( $addedMassVelocityBC eq "" ){ $addedMassVelocityBC=0; }
if( $useTP eq "" ){ $useTP=0; }
if( $fluidSolidCornerFix eq "" ){ $fluidSolidCornerFix=0; }
if( $zfMuByH eq "" ){ $zfMuByH=2.;}
if( $zfRhoHByDt eq "" ){ $zfRhoHByDt=2.;}
if( $zfMono eq "" ){ $zfMono=0.;}
if( $zfMuRhoByZpDt eq "" ){ $zfMuRhoByZpDt=2.;}
if( $surfaceTension eq "" ){ $surfaceTension=0.;}
if( $pAtmosphere eq "" ){ $pAtmosphere=0.;}
# 
setup $domainName
 set solver Cgins
 solver name $solverName
 solver parameters
# 
  incompressible Navier Stokes
  # #wdh# 2012/05/30 Boussinesq model
  $modelNameINS
#   define real parameter kappa $kThermal
  define real parameter thermalExpansivity $thermalExpansivity
  define real parameter adcBoussinesq $adcBoussinesq 
  define real parameter surfaceTension $surfaceTension
  define real parameter pAtmosphere $pAtmosphere
  continue
# 
  $tz
  degree in space $degreeSpace
  degree in time $degreeTime
  OBTZ:frequencies (x,y,z,t) $fx, $fy, $fz, $ft
# 
  $tsINS
#
  if( $numberOfTimeStepCorrections ne "" ){ $cmd="number of PC corrections $numberOfTimeStepCorrections"; }else{ $cmd="#"; }
  $cmd
#
##  useNewImplicitMethod
  $implicitVariation
  implicit factor $implicitFactor 
  refactor frequency $refactorFrequency
# 
  dtMax $dtMax
  debug $debug 
#
  use added mass algorithm $addedMass
  # for now we let the solver know that the added mass algorithm needed predicted values for the pressure:
  # *wdh* April 22, 2018 predicted pressure needed $addedMass
  predicted pressure needed $predictedPressureNeeded
  use new time-stepping startup $useNewTimeSteppingStartup
  useImplicitAmpBCs $useImplicitAmpBCs
  # This next option may be temporary, until we figure out the right thing to do: June, 25, 2018
  added mass velocity BC: $addedMassVelocityBC
  use moving grid sub-iterations $useTP
#
  zfMuByH $zfMuByH 
  zfRhoHByDt $zfRhoHByDt
  zfMono $zfMono 
  zfMuRhoByZpDt $zfMuRhoByZpDt
  fluid solid corner fix: $fluidSolidCornerFix
#
  predicted boundary pressure needed $predictedBoundaryPressureNeeded
# 
  frequency for full grid gen update $freqFullUpdate
#
  pde parameters
    nu  $nu
    kThermal $kThermal
    thermal conductivity $ktc
    fluid density
      $fluidDensity
    gravity
      $gravity
    if( $checkForInflowAtOutflow eq "1" ){ $cmd="OBPDE:check for inflow at outflow"; }else{ $cmd="#"; }
    $cmd 
    if( $useNeumannAtOutflow eq "1" ){ $cmd = "use Neumann BC at outflow"; }else{ $cmd="#"; }
    $cmd 
   done
# 
  $moveCmds
# 
    OBPDE:second-order artificial diffusion $ad2
    OBPDE:ad21,ad22  $ad21, $ad22
    OBPDE:divergence damping  $cdv 
    OBPDE:cDt div damping $cDt
    # Use exact RHS for pressure BC on a wall for testing known solutions and TZ: 
    OBPDE:use exact pressure BC $useExactPressureBC 
    OBPDE:use curl form of the traction $useCurlFormOfTraction
    OBPDE:added mass length scale $addedMassLengthScale
## finish me    OBPDE:added mass velocity BC: $addedMassVelocityBC
# 
$setAxi = $axisymmetric ? "turn on axisymmetric flow" : "#";
$setAxi
$setAxi="";
  if( $commands eq "" ){ $commands_ins="debug $debug";}else{$commands_ins=$commands; } #
  $commands_ins
  $commands_ins=""; 
# 
    maximum number of iterations for implicit interpolation
      10
  pressure solver options
     $psolver
     # yale
     # these tolerances are chosen for PETSc
     relative tolerance
      $rtolp
     absolute tolerance
      $atolp
    maximum allowable increase in the residual
       $ogesDtol
     debug 
       $pdebug
    exit
# 
  implicit time step solver options
     $solver
     # these tolerances are chosen for PETSc
     relative tolerance
      $rtoli
     absolute tolerance
      $atoli
     debug 
       $idebug
    exit
#
  boundary conditions...
   order of extrap for outflow $orderOfExtrapForOutflow (-1=default)
  done
#
  boundary conditions
   $bc 
  done
# 
  initial conditions
  if( $tz eq "turn off twilight zone" && $ic eq "" ){ $ic="uniform flow\n" . "p=1., u=$u0, T=$T0"; }elsif( $ic eq "" ){ $ic="#";}
if( $restart eq "" ){ $icCmds = $ic; }\
  else{ $icCmds = "use grid from show file 1\n always interpolate from show file 1\n read from a show file\n $restart\n -1"; }
    $icCmds
  continue
  $projectInitialConditions
  $extraCmds
# 
  continue
done
# kkc 080328 reset commands so that it can be used in subsequent includes (commands was created in this header
##$commands="";
