#
#   plotStuff plotSolution.cmd -show=MD4.show -name=multiDisk4 -solution=1 
#   plotStuff plotSolution.cmd -show=MD4a.show -name=multiDisk4a -solution=1
#   plotStuff plotSolution.cmd -show=MD8a.show -name=multiDisk8a -solution=1
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
# Foreground colour:0 white
solution: 1
derived types
  E field norm
exit
plot:eFieldNorm
# 
contour
  plot contour lines (toggle)
  coarsening factor 1 (<0 : adaptive)
exit
DISPLAY AXES:0 0
pause
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

