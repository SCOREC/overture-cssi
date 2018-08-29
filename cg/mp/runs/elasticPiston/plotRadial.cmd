#
#  plotStuff plotRadial.cmd -show=radialG2.show -name=radialG2 -solution=6 -vMax=.015
# 
#  plotStuff plotRadial.cmd -show=repG4scf1.show -name=radialG4 -solution=10 -vMax=.25
#
#  plotStuff plotRadial.cmd -show=repG8Scf1.show -name=radialG8 -solution=7 -tp=0p3 -cs=v1 -cf=u
#  plotStuff plotRadial.cmd -show=repG8Scf1.show -name=radialG8 -solution=19 -tp=0p9 -cs=v2 -cf=v
#  plotStuff plotRadial.cmd -show=repG8Scf1.show -name=radialG8 -solution=19 -tp=0p9 -cs=stressNorm -cf=p 
#
#  plotStuff plotRadial.cmd -show=repG8Scf1.show -name=radialG8 -solution=15 -tp=0p7 -cs=vr -cf=vr -cMin=-.21 -cMax=-.1
#
# Radial shear:
#  plotStuff plotRadial.cmd -show=radialShearG4Scf10.show -name=radialShearG4 -solution=7 -tp=0p3 -cs=vTheta -cf=vTheta
#
#  plotStuff plotRadial.cmd -show=radialShearG8Scf10.show -name=radialShearG8 -solution=3 -tp=0p1 -cs=vTheta -cf=vTheta
#  plotStuff plotRadial.cmd -show=radialShearG8Scf10.show -name=radialShearG8 -solution=7 -tp=0p3 -cs=vTheta -cf=vTheta
#  plotStuff plotRadial.cmd -show=radialShearG8Scf10.show -name=radialShearG8 -solution=11 -tp=0p5 -cs=vTheta -cf=vTheta
# 
#  plotStuff plotRadial.cmd -show=radialShearG8Scf10.show -name=radialShearG8 -solution=11 -tp=0p5 -cs=stressNorm -cf=p 
#
$show="radialG2.hdf"; $cMin=""; $cMax=""; $solution=1; $cs="v1"; $cf="u"; $tp="0p0"; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"cMin=f"=>\$cMin,"cMax=f"=>\$cMax,"solution=i"=>\$solution,\
      "cf=s"=>\$cf,"cs=s"=>\$cs ,"tp=s"=>\$tp );
#
$show
#
solution: $solution
frame series:fluidDomain
* 
derived types
  speed
  user defined
    velocity (cylindrical coordinates)
    0 0 0 0 0 1
  exit
exit
contour
  plot:$cf
#  plot:speed
  vertical scale factor 0.
  plot contour lines (toggle)
  if( $cMax ne "" ){ $cmd="min max $cMin $cMax"; }else{ $cmd="#"; }
  $cmd
  exit
frame series:solidDomain
* 
derived types
  speed
 specify velocity components
  0 1 2
 stressNorm
  user defined
    velocity (cylindrical coordinates)
    0 0 0 0 0 1
  exit
exit
contour
  adjust grid for displacement 1
  * plot:vorz
  # plot:speed
  plot:$cs 
  vertical scale factor 0.
  if( $cMax ne "" ){ $cmd="min max $cMin $cMax"; }else{ $cmd="#"; }
  $cmd
  plot contour lines (toggle)
exit
pause
#
DISPLAY COLOUR BAR:0 0
DISPLAY AXES:0 0
$plotName = $name . "$cs$cf" ."t$tp.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0


#
$plotName = $name . "SLt0p5.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0



  # bigger
  DISPLAY SQUARES:0 0
#
 colour boundaries by grid number
if( $name eq "elasticPistonGrid" ){ \
  $cmd="colour grid lines from chosen name\n grid colour 0 RED\n grid colour 1 BLUE\n grid colour 2 GREEN"; }else{ $cmd="#"; }
  $cmd
#
  line width scale factor:0 3
  plot interpolation points 1
  # colour interpolation points 1
#
  hardcopy file name:0 $name.ps
  hardcopy save:0
  hardcopy save:0
