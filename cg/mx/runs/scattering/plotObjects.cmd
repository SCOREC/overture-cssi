#
#   plotStuff plotObjects.cmd -show=objectsG4.show -name=objectsG4t0p3 -solution=6 
#   plotStuff plotObjects.cmd -show=objectsG8.show -name=objectsG8t0p3 -solution=6 
#
#   plotStuff plotObjects.cmd -show=objectsG4.show -name=objectsG4t1p0 -solution=21
#   plotStuff plotObjects.cmd -show=objectsG8.show -name=objectsG8t1p0 -solution=21
#
#   plotStuff plotObjects.cmd -show=objectsG4.show -name=objectsG4t5p0 -solution=101
#   plotStuff plotObjects.cmd -show=objectsG8.show -name=objectsG8t5p0 -solution=101
#
# new grids: 
#   plotStuff plotObjects.cmd -show=objectsG4a.show -name=objectsG4at1p0 -solution=21
#   plotStuff plotObjects.cmd -show=objectsG8a.show -name=objectsG8at1p0 -solution=21
#
#   plotStuff plotObjects.cmd -show=objectsG4a.show -name=objectsG4at5p0 -solution=101
#   plotStuff plotObjects.cmd -show=objectsG8a.show -name=objectsG8at5p0 -solution=101
#
# Resonant mode: 
#   plotStuff plotObjects.cmd -show=objectsG4M0.show -name=objectsG4M0t1p0 -solution=21
#   plotStuff plotObjects.cmd -show=objectsG8M0.show -name=objectsG8M0t1p0 -solution=21
#
#   plotStuff plotObjects.cmd -show=objectsG4M0.show -name=objectsG4M0t5p0 -solution=101
#   plotStuff plotObjects.cmd -show=objectsG8M0.show -name=objectsG8M0t5p0 -solution=101
#
# -- new grids 
#   plotStuff plotObjects.cmd -show=objectsG4M0a.show -name=objectsG4M0at1p0 -solution=21
#   plotStuff plotObjects.cmd -show=objectsG8M0a.show -name=objectsG8M0at1p0 -solution=21
#
#   plotStuff plotObjects.cmd -show=objectsG4M0a.show -name=objectsG4M0at5p0 -solution=101
#   plotStuff plotObjects.cmd -show=objectsG8M0a.show -name=objectsG8M0at5p0 -solution=101
#
# -- longer time: 
#   plotStuff plotObjects.cmd -show=objectsG4b.show 
#   plotStuff plotObjects.cmd -show=objectsG4M1.show
#
# Movie:
#   plotStuff plotObjects.cmd -show=objectsG4m.show -EyMin=-1.7 -EyMax=1.1 -name=cgmxScat12ObjectsGDM
#
$show="ellipseG8.hdf"; $solution="-1"; $EyMin=-1.5; $EyMax=1.5; $ExMin=-1.5; $ExMax=1.5; $HzMin=-1.; $HzMax=1.;
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,\
   "ExMin=f"=>\$ExMin,"ExMax=f"=>\$ExMax, "EyMin=f"=>\$EyMin,"EyMax=f"=>\$EyMax  );
#
$show
# 
contour
  plot contour lines (toggle)
  coarsening factor 1 (<0 : adaptive)
  # Set max and min for each component 
  plot:Hz
    min max $HzMin $HzMax
  plot:Ex
    min max $ExMin $ExMax
  plot:Ey
    min max $EyMin $EyMax
  vertical scale factor 0.
exit
solution: $solution
# movie: 
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0
set view:0 0.0155536 -0.021148 0 1.83688 1 0 0 0 1 0 0 0 1
movie file name: $name
save movie files 1


# zoom: 
  set view:0 -0.0307779 -0.0121165 0 2.08435 1 0 0 0 1 0 0 0 1

pause
# hardcopy vertical resolution:0 2048
# hardcopy horizontal resolution:0 2048
# line width scale factor:0 3
DISPLAY AXES:0 0
#
  plot
  $plotName = $name . "Ey.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
# 
  plot:Ex
  $plotName = $name . "Ex.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
# 
exit
exit


erase
stream lines
  streamline density 100
  arrow size 0.05
exit
  $plotName = $name . "SL.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0