#   ---- Plot errors over time in the super-grid ----
# 
# BA-GDM -- different super-grid widths 
#   plotStuff plotSuperGridErrors.cmd -show=baMatGDMGenNDSGWidth0p05.show -name=baMatGDMGenNDSGWidth0p05
#   plotStuff plotSuperGridErrors.cmd -show=baMatGDMGenNDSGWidth0p1.show -name=baMatGDMGenNDSGWidth0p1
#   plotStuff plotSuperGridErrors.cmd -show=baMatGDMGenNDSGWidth0p2.show -name=baMatGDMGenNDSGWidth0p2
#
# THREE-DIMENSIONS
#   plotStuff plotSuperGridErrors.cmd -show=cboxSmall8Air.show -name=cboxSmall8AirWidth0p1
#
$show="baMatIntGDM.show.hdf"; $solution="-1"; $name="plot"; $field="Ey"; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field );
#
$show
# show forcing regions 1
#
solution: $solution
plot sequence:errors
  add Ey_error
  add Ez_error
  add Hx_error
  add Hy_error
  add Hz_error
  $fileName = $name . ".m"; 
  save results to a matlab file
  $fileName 


contour
  plot:Hz
  # plot contour lines (toggle)
exit
DISPLAY AXES:0 0
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 4
#
plot:Hz
$plotName = $name . "Hz.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
pause
# 
contour
  plot:Hz error
  ghost lines -20
#
$plotName = $name . "HzErr.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0



derived types
  E field norm
exit
contour
  plot:eFieldNorm

  plot contour lines (toggle)
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  # min max 0 0.8
exit
solution: $solution
pause
DISPLAY AXES:0 0
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 4
plot 
$plotName = $name . "EfieldNorm.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0


pause
#
Foreground colour:0 white
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 4
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0
# 
plot
$plotName = $name . "EfieldNorm.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0


plot:$field
contour
  plot contour lines (toggle)
  coarsening factor 1 (<0 : adaptive)
  # min max -1 1
 vertical scale factor 0.
exit
pause
#
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 3
DISPLAY AXES:0 0
#
#
solution: $solution
pause
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
DISPLAY AXES:0 0
#
  line width scale factor:0 3
  plot
  $plotName = $name . "Ey.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
# 
plot:Ex
  $plotName = $name . "Ex.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0





$fieldNoBlanks=$field; $fieldNoBlanks =~ s/ //g; 
$cmd="#";
for( $i=0; $i<$numToSave; $i++ ){\
  $solution = $numPerTime*$i +1; $t = $i*$tSave; $t =~ s/\./p/g;  \
  $cmd .= "solution: $solution\n"; \
  $plotName = $name . $fieldNoBlanks . "t$t.ps"; \
  $cmd .= "plot:$field\n"; \
  $cmd .= "hardcopy file name:0 $plotName\n"; \
  $cmd .= "hardcopy save:0\n"; \
}
$cmd 
# save a matlab file with the errors
plot sequence:errors
  Ey_error
  add Ex_error
  save results to a matlab file
  $matlabFile = $name . "Err.m"; 
  $matlabFile
# 



plot:Hz
  $plotName = $name . "Hz.ps"; 
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
