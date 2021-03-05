#
#  ---- Include file for cgmx cmd scripts: set common options  ------
#
pml width,strength,power $pmlLines $pmlStrength $pmlPower
interface BC iterations $interfaceIterations
# interfaceEquationsOption=0 : use extrap for 2nd ghost, 1=use eqns
interface equations option $interfaceEquationOption
omega for interface iterations $interfaceOmega
#
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
plot divergence 0
plot errors $checkErrors
check errors $checkErrors
plot polarization components $plotPolarizationComponents
plot intensity $plotIntensity
intensity option $intensityOption
$intensityAveragingInterval=1.;
intensity averaging interval $intensityAveragingInterval
#
# We can compute errors compared to another solution (e.g. to test far-field BC's)
if( $compareToShowFile ne "" ){ $cmd ="compare to show file 1\n reference show file: $compareToShowFile\n plot errors 1\n check errors 1"; }else{ $cmd="#"; }
$cmd 
#
