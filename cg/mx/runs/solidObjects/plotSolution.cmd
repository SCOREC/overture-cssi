#
#   plotStuff plotSolution.cmd -show=solidsO2G4.show -name=solidsO2G4 -solution=1 
# 
#   solutions 21 or 17 are nice:
#   plotStuff plotSolution.cmd -show=solidsO4G6.show -name=sixSolidObjects4G6Enorm -solution=17
#
# RPI: 
#   plotStuff plotSolution.cmd -show=rpiG8O4.show -name=rpiG8O4Enormt2p5 -solution=6
#
# Movies:
#   plotStuff plotSolution.cmd -show=rpiG8O4a.show -name=rpiScat
#
#   plotStuff plotSolution.cmd -show=solidsO4G6a.show -name=solidScat
#
$show="ellipseG8.hdf"; $solution="-1"; $EyMin=-1.5; $EyMax=1.5; $ExMin=-1.5; $ExMax=1.5; $HzMin=-1.; $HzMax=1.;
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,\
   "ExMin=f"=>\$ExMin,"ExMax=f"=>\$ExMax, "EyMin=f"=>\$EyMin,"EyMax=f"=>\$EyMax  );
#
$show
# 
solution: $solution
derived types
  E field norm
exit
plot:eFieldNorm
# 
contour
  plot contour lines (toggle)
  coarsening factor 1 (<0 : adaptive)
  min max 0 1
  vertical scale factor 0.
exit


#
Foreground colour:0 white
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0
Foreground colour:0 wheat
bigger
solution: 1
# --------------- movie  - -lower resolution 
# hardcopy vertical resolution:0 512
# hardcopy horizontal resolution:0 512
line width scale factor:0 2
movie file name: $name
save movie files 1


solution: $solution
pause
$plotName = $name . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0


#
$numPerTime=2; $numToSave=7; $tSave=1.;
$field="eFieldNorm"; 
$fieldNoBlanks="Enorm"; 
$cmd="#";
for( $i=0; $i<$numToSave; $i++ ){\
  $solution = $numPerTime*$i +1; $t = $i*$tSave; $t =~ s/\./p/g;  \
  $cmd .= "solution: $solution\n"; \
  $plotName = $name . $fieldNoBlanks . "t$t.ps"; \
  $cmd .= "hardcopy file name:0 $plotName\n"; \
  $cmd .= "hardcopy save:0\n"; \
}
$cmd



# -- movie: 
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
