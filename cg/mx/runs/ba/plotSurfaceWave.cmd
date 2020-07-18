#
# 
#   plotStuff plotSurfaceWave.cmd -show=baSurfaceWaveG16.show -name=baSurfaceWaveG16 -solution=11
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
#
@comp = ( "Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "Ex error", "Ey error", "Ez error", "Hx error", "Hy error", "Hz error");
@compName = ( "Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "ExErr", "EyErr", "EzErr", "HxErr", "HyErr", "HzErr");
$cmd="#"; 
for( $i=0; $i<@comp; $i++ ){\
   $plotName = $name . $compName[$i] . ".ps"; \
   $cmd .= "\n plot:$comp[$i]\n hardcopy file name:0 $plotName\n hardcopy save:0";\
 }
# printf("cmd=$cmd\n");
$cmd
# -------- plot rotated surfaces
plot:Ey error
contour
  vertical scale factor 0.2
  plot contour lines (toggle)
#
exit
set view:0 -0.0361329 0 0 1.59135 0.173648 -0.17101 0.969846 0.984808 0.0301537 -0.17101 -1.09336e-16 0.984808 0.173648
@comp = ( "Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "Ex error", "Ey error", "Ez error", "Hx error", "Hy error", "Hz error");
@compName = ( "Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "ExErr", "EyErr", "EzErr", "HxErr", "HyErr", "HzErr");
$cmd="#"; 
for( $i=0; $i<@comp; $i++ ){\
   $plotName = $name . $compName[$i] . "Rotated.ps"; \
   $cmd .= "\n plot:$comp[$i]\n hardcopy file name:0 $plotName\n hardcopy save:0";\
 }
$cmd

 hardcopy file name:0 baSurfaceWaveG32EyErrRotatedB.ps
 hardcopy save:0


plot:Ex
$plotName = $name . "Ex.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
plot:Ey
$plotName = $name . "Ey.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
plot:Ez
$plotName = $name . "Ez.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
plot:Hx
$plotName = $name . "Hx.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
plot:Hy
$plotName = $name . "Hy.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
plot:Hz
$plotName = $name . "Hz.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
