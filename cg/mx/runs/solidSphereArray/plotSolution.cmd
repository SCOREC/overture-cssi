#
# 
#   plotStuff plotSolution.cmd -show=fourSolidEllipsoids.show -name=fourSolidEllipsoids -solution=3
#   plotStuff plotSolution.cmd -show=eightSolidEllipsoids.show -name=eightSolidEllipsoids -solution=3
#   plotStuff plotSolution.cmd -show=36SolidEllipsoids.show -name=36SolidEllipsoids -solution=7 -y1=-1.3 -y2=1.3 -z1=-1.95 
#   plotStuff plotSolution.cmd -show=125SolidEllipsoids.show -name=125SolidEllipsoids -solution=9 -y1=-1.3 -y2=1.3 -z1=-1.95 -eMax=.5 -view=2
#
# MLA: 
#   plotStuff plotSolution.cmd -show=fourMlaEllipsoids.show -name=fourSolidMlaEllipsoids -solution=5
#   plotStuff plotSolution.cmd -show=eightMlaEllipsoids.show -name=eightSolidMlaEllipsoids -solution=6
#   plotStuff plotSolution.cmd -show=36MlaEllipsoids.show -name=36SolidMlaEllipsoids -solution=6 -y1=-1.3 -y2=1.3 -z1=-1.95 
#
$show="fourSolidEllipsoids.hdf"; $solution="-1"; $name="fourSolidEllipsoids"; $field="Ey"; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
$y1=-.65; $y2=.65; $z1=-.65; $eMax=""; $view=1;
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field,\
      "y1=f"=>\$y1,"y2=f"=>\$y2,"z1=f"=>\$z1,"eMax=f"=>\$eMax,"view=i"=>\$view );
#
$show
#
derived types
  E field norm
exit
solution: $solution
contour
  plot:eFieldNorm
  plot contour lines (toggle)
  # for 8 ellipsoid:
  delete contour plane 2
  delete contour plane 1
  delete contour plane 0
  add contour plane  0 0 1  0   0 $z1
  add contour plane  0 1 0  0 $y1   0 
  add contour plane  0 1 0  0 $y2   0 
  # 125:
  if( $eMax ne "" ){ $cmd="min max 0.0 $eMax"; }else{ $cmd="#"; }
  $cmd 
exit
grid
  toggle grid 0 0
  plot block boundaries 0
  plot shaded surfaces (3D) 0
  coarsening factor 3
exit
# 36 : 
$cmd="set view:0 -0.00241692 -0.00483384 0 1.27701 0.866025 0.17101 -0.469846 0 0.939693 0.34202 0.5 -0.296198 0.813798";
# 125 : 
if( $view eq 2 ){ $cmd="set view:0 0.0465255 -0.0997615 0 1.06004 0.866025 0.17101 -0.469846 0 0.939693 0.34202 0.5 -0.296198 0.813798"; }
$cmd
# pause
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
# line width scale factor:0 2
# make lines faint 
line width scale factor:0 3
Foreground colour:0 white
plot 
#
$plotName = $name . "EfieldNorm.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
pause
#
plot:Py
$plotName = $name . "Py.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
pause
#
plot:N3
$plotName = $name . "N3.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0



  plot contour lines (toggle)
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  # min max 0 0.8
exit
# --- movie ---
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0
bigger:0
movie file name: $name
save movie files 1
show movie


pause
#
Foreground colour:0 white
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 4
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0
# 
plot
$plotName = $name . "EfieldNorm.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0

