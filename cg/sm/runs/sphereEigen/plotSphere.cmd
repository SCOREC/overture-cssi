#
# plotStuff plotSphere -show=sphere2ClassIIn2m1.show
# 
# Plot errors over time:
# plotStuff plotSphere -show=../sphere2/sphere2c.show -matlab=sphere2cErrors.m
# plotStuff plotSphere -show=sphere2g.show -matlab=sphere2gErrors.m
#
# plotStuff plotSphere -show=../sphere2/sphere4ca.show -matlab=sphere4cErrors.m
# plotStuff plotSphere -show=sphere4ga.show -matlab=sphere4gErrors.m
#
# cp {sphere2cErrors.m,sphere2gErrors.m} $smog/sphereEigen   
# cp {sphere4cErrors.m,sphere4gErrors.m} $smog/sphereEigen
#
# Movie:
#  plotStuff plotSphere -show=sphere2
#
$matlab="sphere2cErrors.m";
GetOptions( "show=s"=> \$show, "name=s"=> \$name, "matlab=s"=> \$matlab );
#
$show
# 
displacement
  y-r:0
  y-r:0
  y-r:0
  x+r:0
  x+r:0
  plot block boundaries 0
  DISPLAY SQUARES:0 0
  displacement scale factor .05
  exit this menu
# 
plot sequence:errors
  reset:0
  uError
  add vError
  add wError
  save results to a matlab file
  $matlab

# ++++++++++++ MOVIE ++++++++++++++
set view:0 0.0120846 0.00302115 0 1.03115 0.866025 -0.17101 0.469846 0 0.939693 0.34202 -0.5 -0.296198 0.813798
bigger 1.1
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
movie file name: deformingSphere
solution: 1
save movie files 1
show movie



set view:0 -0.221792 -0.113275 0 1.23036 0.866025 -0.17101 0.469846 0 0.939693 0.34202 -0.5 -0.296198 0.813798
DISPLAY AXES:0 0


# 
line width scale factor:0 4
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
# 
$cmds = ""; 
for( $i=1; $i < 22; $i=$i+2 )\
{ \
$cmds .= "solution $i\n"; \
$cmds .= "hardcopy file name:0 sphereClassIIn2m1t$i.ps\n hardcopy save:0\n"; \
}; $cmds .= "#"; 
$cmds


displacement
  set plot bounds
  -1.5 1.5 -1.5 1.5 1.5 -.15
  displacement scale factor .08
  exit this menu