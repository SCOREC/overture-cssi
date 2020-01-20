#
# 
#   plotStuff plotSurfaceWave.cmd -show=baSurfaceWaveG16.show -name=baSurfaceWaveG32 -solution=11
#   plotStuff plotSurfaceWave.cmd -show=baSurfaceWaveG32.show -name=baSurfaceWaveG32 -solution=11
#
#   plotStuff plotSurfaceWave.cmd -show=baSurfaceWaveTimePeriodicG8.show -name=baSurfaceWaveTimePeriodicG8 -solution=11
#
$show="baMatIntGDM.show.hdf"; $solution="-1"; $name="plot"; $field="Ey"; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field );
#
$show
show forcing regions 1
#
derived types
  E field norm
exit
contour
  plot contour lines (toggle)
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  # min max 0 0.8
exit
x-
solution: $solution
pause
DISPLAY AXES:0 0
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 4
# 
plot:Ex
$plotName = $name . "Ex.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
plot:Ey
$plotName = $name . "Ey.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
