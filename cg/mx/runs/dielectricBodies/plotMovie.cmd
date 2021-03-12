#
#   plotStuff plotMovie.cmd -show=ellipseL10A30G8movie.show -name=ellipseL10A30G8EyMovie
#
$show="ellipseG8.hdf"; $solution="-1"; $name="plot"; $field="Ey"; 
$eMin=-.25; $eMax=1.; 
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution );
#
$show
#
plot:$field
contour
  plot contour lines (toggle)
  coarsening factor 1 (<0 : adaptive)
  # min max $eMin $eMax
 vertical scale factor 0.
exit
DISPLAY AXES:0 0
DISPLAY COLOUR BAR:0 0
DISPLAY LABELS:0 0
bigger:0
movie file name: $name
save movie files 1


