#
# 
#   plotStuff plotSolution.cmd -show=waveyDiskb.show -name=waveyDiskb -solution=21
#   plotStuff plotSolution.cmd -show=waveyDiskc.show -name=waveyDiskc -solution=21
# 
#   plotStuff plotSolution.cmd -show=CBodyEps4G8.show -name=splitRing -solution=1
# 
#   plotStuff plotSolution.cmd -show=64DiskG4.show -name=64Disk -solution=16 -emin=0 -emax=.35 
#
#   plotStuff plotSolution.cmd -show=oneDiskTwoLayersMlaGdm.show -name=oneDiskTwoLayersMlaGdm -solution=13 -emin=0 -emax=.85 
#
$show="ellipseG8.hdf"; $solution="-1"; $name="plot"; $field="Ey"; $emin=0; $emax=-1; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax );
#
$show
#
derived types
  E field norm
exit
contour
  plot:eFieldNorm
  plot contour lines (toggle)
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  if( $emax > $emin ){ $cmd="min max $emin $emax"; }else{ $cmd="#"; }
  $cmd
exit
solution: $solution
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

# --- movie ---
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0


bigger:0
movie file name: $name
save movie files 1
show movie


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
