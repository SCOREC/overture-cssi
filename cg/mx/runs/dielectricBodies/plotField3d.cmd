#
#   plotStuff plotField3d.cmd -show=sphereG4.show -name=sphereG4
#   plotStuff plotField3d.cmd -show=solidPillBox.show -name=solidPillBox
#
#   plotStuff plotField3d.cmd -show=solidSphereEps4G6.show -name=solidSphere
#
#   plotStuff plotField3d.cmd -show=solidEllipsoidG4.show -name=solidEllipsoid
# 
#   plotStuff plotField3d.cmd -show=IBeam3dG8a.show -name=IBeam3da 
# 
#   plotStuff plotField3d.cmd -show=IBeam3dG8.show -name=IBeam3d
#
# plotStuff plotField3d.cmd -show=ellipsoid543G4 -solution=5 -name=ellipsoid543G4
#
# plotStuff plotField3d.cmd -show=ellipsoid543Eps4G4 -solution=5 -name=ellipsoid543Eps4G4
# plotStuff plotField3d.cmd -show=ellipsoid1aEpsG4 -solution=5 -name=ellipsoid1aEpsG4
# plotStuff plotField3d.cmd -show=ellipsoid1aEpsG8 -solution=5 -name=ellipsoid1aEpsG8
# plotStuff plotField3d.cmd -show=solidSphereRad0p5G4i -name=solidSphereRad0p5G4i
# 
# plotStuff plotField3d.cmd -show=solidEllipsoid546 -name=solidEllipsoid546 -sol=3
#
$show="sphereG4.hdf"; $sol=21; 
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "sol=s"=>\$sol );
#
$show
#
# previous
solution: $sol
derived types
  E field norm
exit
plot:eFieldNorm
# plot:Ey
Foreground colour:0 white 
contour
  contour lines 0
  # min max -1 1
  #
  plot the grid
    toggle grid 0 0
    plot block boundaries 0
    plot shaded surfaces (3D) 0
    coarsening factor 3
  exit this menu
 exit
  DISPLAY COLOUR BAR:0 0
  DISPLAY AXES:0 0
  DISPLAY LABELS:0 0
# ellispoid: 
set view:0 -0.0060423 0.00302115 0 1.05751 0.939693 0.116978 -0.321394 0 0.939693 0.34202 0.34202 -0.321394 0.883022  
pause
  $plotName = $name . "Enorm.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0





# --------------  I-BEAM 3D  ------------
$sol=19;
set view:0 -0.00302115 -0.00302115 0 1.04747 0.939693 -0.0593912 0.336824 0 0.984808 0.173648 -0.34202 -0.163176 0.925417
solution: $sol
derived types
  E field norm
exit
plot:eFieldNorm
contour
  contour lines 0
  # delete contour plane 0
  # min max -1 1
  #
  plot the grid
    toggle grid 0 0
    plot block boundaries 0
    # plot shaded surfaces (3D) 0
    grid colour 1 BRASS
    grid colour 3 BRASS
    plot grid lines 0
    # coarsening factor 3
  exit this menu
  DISPLAY AXES:0 0
pause
  $plotName = $name . "Enorm.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
  

# ---solidSphereEps4G6 : 
set view:0 -0.0469789 -0.0060423 0 1.0216 0.939693 0.0593912 -0.336824 0 0.984808 0.173648 0.34202 -0.163176 0.925417
# ellipsoid: 
set view:0 -0.0969789 0 0 1.0216 0.939693 -0.296198 0.17101 0 0.5 0.866025 -0.34202 -0.813798 0.469846
solution: 21
plot:Ey
contour
  contour lines 0
  # min max -1 1
  #
  plot the grid
    toggle grid 0 0
    plot block boundaries 0
    plot shaded surfaces (3D) 0
    coarsening factor 2
  exit this menu
  DISPLAY AXES:0 0
pause
#
  line width scale factor:0 3
  plot
  $plotName = $name . "Ey.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
  


solution: 31
plot:Ey
set view:0 0.00302115 -0.0120846 0 1.02795 0.866025 0.25 -0.433013 0 0.866025 0.5 0.5 -0.433013 0.75
contour
  contour lines 0
  min max -1 1
  #
  plot the grid
    toggle grid 0 0
    plot block boundaries 0
    plot shaded surfaces (3D) 0
    coarsening factor 4
  exit this menu
  DISPLAY AXES:0 0
pause
#
  line width scale factor:0 3
  plot
  $plotName = $name . "Ey.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
# 
  set view:0 -0.00712762 -0.0201969 0 3.20091 1 0 0 0 1 0 0 0 1
  DISPLAY COLOUR BAR:0 0
  $plotName = $name . "EyZoom.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0

contour
  plot contour lines (toggle)
  exit
solution: 31
plot:Ey
#
  line width scale factor:0 3
  plot
  $plotName = $name . "Ey.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
# 
erase
stream lines
  streamline density 100
  arrow size 0.05
exit
  $plotName = $name . "SL.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
