# 
# plotStuff plotSphere -show=sphere2ClassIIn2m1.show 
# 
GetOptions( "show=s"=> \$show, "name=s"=> \$name ); 
# 
$show 
# 
displacement 
  set plot bounds 
  -1.5 1.5 -1.5 1.5 1.5 -.15 
  displacement scale factor .08 
  exit this menu 
contour
  adjust grid for displacement 1
  plot the grid
    plot shaded surfaces (3D) 0
    y-r:0
    y-r:0
    x+r:0
    x+r:0
    x+r:0
    exit this menu
  pick to delete contour planes
  delete contour plane 2
  delete contour plane 1
  delete contour plane -1
  delete contour plane -1
  bigger:0
  bigger:0
  x-r:0
  x-r:0
  y-r:0
  y-r:0
  delete contour plane -1
  delete contour plane -1
  delete contour plane 0
  pick to add boundary surface
  add boundary surface 1 1 2 
  add boundary surface 3 1 2 
  reset:0
  y-r:0
  y-r:0
  y-r:0
  y-r:0
  x+r:0
  x+r:0
  x+r:0
  plot:s11
  exit
erase
contour
  exit
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
derived types
stressNorm
exit
plot:stressNorm
y-r:0
y-r:0
y-r:0
y-r:0
y-r:0
y-r:0
y-r:0
y-r:0
y+r:0
y+r:0
y+r:0
y+r:0
y+r:0
y+r:0
y+r:0
y+r:0
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
previous
plot:v1
derived types
specify velocity components
6 7 8
speed
exit
plot:speed
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
plot:s12
plot:s22
plot:s11
next
next
next
next
next
plot:u
reset:0
y-r:0
y-r:0
y-r:0
x+r:0
x+r:0
plot:w
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
next
plot:v
plot:s33
plot:s11
