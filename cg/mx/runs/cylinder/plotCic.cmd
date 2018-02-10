#
# plotStuff plotCic.cmd
#
$show="cicSosup.show";
# 
$show
#
DISPLAY AXES:0 0
plot:Ey
contour
  # min max -0.8 0.8
  min max -0.85 1.
  exit
save movie files 1
plot:Ex error
movie file name: cicSosupExError
show movie


movie file name: cicSosupEy
show movie

movie file name: cicSosupEx
show movie
