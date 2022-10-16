# 
# Plot results from two domain CHT runs
#
# doubleAnnulus: 
#   plotStuff plotSolution -show=doubleAnnlusG4.show -name=doubleAnnulusCHT -solution=3 -Tmin=-.2 -Tmax=.2 -vsf=.6
# Order=4: 
#   plotStuff plotSolution -show=doubleAnnulusO4G4.show -name=doubleAnnulusO4G4 -solution=3 -Tmin=-.08 -Tmax=.08 -errMin=-4.4e-6 -errMax=4.4e-6
#   plotStuff plotSolution -show=doubleAnnulusO4G8.show -name=doubleAnnulusO4G8 -solution=3 -Tmin=-.08 -Tmax=.08 -vsf=.6
#
#   plotStuff plotSolution -show=diskArray2x2yG4.show -name=diskArray2x2yG4T -solution=6 -Tmin=-1.0 -Tmax=1.0 -errMin=-3.0e-6 -errMax=2.6e-6
# 
#   plotStuff plotSolution -show=fourDisksDomainG2.show -name=fourDisksDomainT -solution=6 -numDomains=5 -domainName=outerDomain diskDomain1 diskDomain2 diskDomain3 diskDomain4 -Tmin=-1.0 -Tmax=1.0 
#   plotStuff plotSolution -show=fourDisksDomainG4.show -name=fourDisksDomainT -solution=6 -numDomains=5 -domainName=outerDomain diskDomain1 diskDomain2 diskDomain3 diskDomain4 -Tmin=-1.0 -Tmax=1.0 -errMin=-3.0e-7 -errMax=4e-7
#   plotStuff plotSolution -show=fourDisksDomainG4ui.show -name=fourDisksDomainT -solution=6 -numDomains=5 -domainName=outerDomain diskDomain1 diskDomain2 diskDomain3 diskDomain4 -Tmin=-1.0 -Tmax=1.0 -errMin=-6.0e-7 -errMax=8e-7
# 
$show="fourHeatedDisks.show"; $name="fourHeatedDisks"; $Tmin=0; $Tmax=-1; $numToPlot=1; $vsf=0.; 
$errMin=0; $errMax=-1; 
$numDomains=2; 
# $domain1="outerDomain"; $domain2="innerDomain"; 
@domainName = ();  # these must be null for GetOptions to work, defaults are given below 
* ----------------------------- get command line arguments --------------------------------------- 
GetOptions( "show=s"=>\$show, "name=s"=>\$name,"solution=i"=>\$solution,"Tmin=f"=>\$Tmin,"Tmax=f"=>\$Tmax ,\
            "numToPlot=i"=>\$numToPlot,"domain1=s"=>\$domain1,"domain2=s"=>\$domain2,"vsf=f"=>\$vsf,\
            "errMin=f"=>\$errMin,"errMax=f"=>\$errMax,"numDomains=i"=>\$numDomains,"domainName=s{1,}"=>\@domainName ); 
# 
if( $domainName[0] eq "" ){ @domainName=("outerDomain","innerDomain"); }# defaults if not set 
#
$show 
# 
# 
solution: $solution
#
plot:T
#
if( $Tmax>$Tmin ){ $minMaxCmd="min max $Tmin $Tmax"; }else{ $minMaxCmd="#"; }
$cmd="#"; 
for( $d=0; $d<$numDomains; $d++ ){\
  $cmd .= "\nframe series:$domainName[$d]\n"  \
          . "contour \n"                        \
          . "  # plot contour lines (toggle)\n" \
          . "  vertical scale factor $vsf\n"    \
          . "  $minMaxCmd\n"                    \
          . "exit";                             \
}
$cmd
# 
DISPLAY AXES:0 0
# side view
# set view:0 -0.0758308 -0.0151057 0 1.05751 1 0 0 0 0.5 -0.866025 0 0.866025 0.5
pause
plot
$plotName = $name . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
# plot Errors 
#
#
if( $errMax>$errMin ){ $minMaxCmd="min max $errMin $errMax"; }else{ $minMaxCmd="#"; }
$cmd="#"; 
for( $d=0; $d<$numDomains; $d++ ){\
  $cmd .= "\nframe series:$domainName[$d]\n"    \
          . "contour \n"                        \
          . "plot:T (error)\n"                  \
          . "  # plot contour lines (toggle)\n" \
          . "  vertical scale factor $vsf\n"    \
          . "  $minMaxCmd\n"                    \
          . "exit";                             \
}
$cmd
$plotName = $name . "Error.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0  


frame series:$domainName[0]
contour 
  # plot contour lines (toggle)
  vertical scale factor $vsf
  $cmd
exit
frame series:$domainName[1]
contour
  # plot contour lines (toggle)
  vertical scale factor $vsf
  if( $Tmax>$Tmin ){ $cmd="min max $Tmin $Tmax"; }else{ $cmd="#"; }
  $cmd
exit
# 
DISPLAY AXES:0 0
# side view
# set view:0 -0.0758308 -0.0151057 0 1.05751 1 0 0 0 0.5 -0.866025 0 0.866025 0.5
pause
plot
$plotName = $name . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
# plot Errors 
#
#
frame series:innerDomain
plot:T (error)
contour
  min max $errMin $errMax
exit
frame series:outerDomain
plot:T (error)
contour
  min max $errMin $errMax
exit
$plotName = $name . "Error.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0  


# --- movie:
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
bigger:0
DISPLAY COLOUR BAR:0 0
movie file name: $name
save movie files 1
show movie

# -- save one hardcopy
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0
# 
plot
$plotName = $name . $label . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
