# 
# Plot results from two domain CHT runs
#
# doubleAnnulus: 
#   plotStuff plotSolution -show=doubleAnnlusG4.show -name=doubleAnnulusCHT -solution=3 -Tmin=-.2 -Tmax=.2 -vsf=.6
# 
$show="fourHeatedDisks.show"; $name="fourHeatedDisks"; $Tmin=0; $Tmax=-1; $numToPlot=1; $vsf=0.; 
$domain1="outerDomain"; $domain2="innerDomain"; 
* ----------------------------- get command line arguments --------------------------------------- 
GetOptions( "show=s"=>\$show, "name=s"=>\$name,"solution=i"=>\$solution,"Tmin=f"=>\$Tmin,"Tmax=f"=>\$Tmax ,\
            "numToPlot=i"=>\$numToPlot,"domain1=s"=>\$domain1,"domain2=s"=>\$domain2,"vsf=f"=>\$vsf ); 
# 
$show 
# 
# -- set displacement components for Godunov 
frame series:$domain1
# derived types 
# speed 
# exit 
# 
# previous 
# plot:stressNorm 
solution: $solution
#
plot:T
contour 
  # plot contour lines (toggle)
  vertical scale factor $vsf
  if( $Tmax>$Tmin ){ $cmd="min max $Tmin $Tmax"; }else{ $cmd="#"; }
  $cmd
exit
frame series:$domain2
contour
  # plot contour lines (toggle)
  vertical scale factor $vsf
  if( $Tmax>$Tmin ){ $cmd="min max $Tmin $Tmax"; }else{ $cmd="#"; }
  $cmd
exit
DISPLAY AXES:0 0
set view:0 -0.0758308 -0.0151057 0 1.05751 1 0 0 0 0.5 -0.866025 0 0.866025 0.5
pause
plot
$plotName = $name . ".ps"; 
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
