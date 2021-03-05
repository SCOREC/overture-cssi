#
# 
#   plotStuff plotSolution.cmd -show=baMatIntGDM.show -name=baMatIntGDM -solution=21
#   plotStuff plotSolution.cmd -show=baSixBox.show -name=baSixBox -solution=11
#   plotStuff plotSolution.cmd -show=cylBox.show -name=baCylBox -solution=7
#   plotStuff plotSolution.cmd -show=baGI.show -name=baGI -solution=6
# 
#   plotStuff plotSolution.cmd -show=baDieCyl.show -name=baDieCyl -solution=6
#
# 2D BA interface
#   plotStuff plotSolution.cmd -show=baMatIntG8.show -name=baMatInt2d -solution=11
#   plotStuff plotSolution.cmd -show=baMatIntG32.show -name=baMatInt2d -solution=11
#   plotStuff plotSolution.cmd -show=baMatIntG32NoDiss.show -name=baMatInt2dNoDiss -solution=11
#
# BA co-centric squares and cyls
#   plotStuff plotSolution.cmd -show=baBoxCylArrayG32.show -name=baBoxCylArrayG32t1p5 -eMax=.9 -solution=16
#   plotStuff plotSolution.cmd -show=baBoxCylArrayG32.show -name=baBoxCylArrayG32t1p0 -eMax=.9 -solution=11
#   plotStuff plotSolution.cmd -show=baBoxCylArrayG64.show -name=baBoxCylArrayG64t1p0 -eMax=.9 -solution=11
#   plotStuff plotSolution.cmd -show=baBoxCylArrayG64.show -name=baBoxCylArrayG64t1p5 -eMax=.9 -solution=16
#
#   plotStuff plotSolution.cmd -show=baBoxCylArrayG128.show -name=baBoxCylArrayG128t0p8 -eMax=.9 -solution=9
#   plotStuff plotSolution.cmd -show=baBoxCylArrayG128.show -name=baBoxCylArrayG128t1p0 -eMax=.9 -solution=11
#   plotStuff plotSolution.cmd -show=baBoxCylArrayG128.show -name=baBoxCylArrayG128t1p2 -eMax=.9 -solution=13
#   plotStuff plotSolution.cmd -show=baBoxCylArrayG128.show -name=baBoxCylArrayG128t1p5 -eMax=.9 -solution=16
#
#   plotStuff plotSolution.cmd -show=fourDiskEps4BAH -name=fourDiskEps4BAHSolution -field=Ey -solution=4
#   plotStuff plotSolution.cmd -show=16DiskSlabsEps4Eps8Eps6Eps3 -name=16DiskSlabsEps4Eps8Eps6Eps3 -field=Ey -solution=9
# 
$show="baMatIntGDM.show.hdf"; $solution="-1"; $name="plot"; $field="eFieldNorm"; $eMin=0; $eMax=0; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field,"eMin=f"=>\$eMin,"eMax=f"=>\$eMax );
#
$show
show forcing regions 1
#
derived types
  E field norm
exit
contour
  plot:$field
  plot contour lines (toggle)
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  if( $eMax>$eMin ){ $cmd="min max $eMin $eMax"; }else{ $cmd="#"; }
  $cmd
  # min max 0 0.8
exit
solution: $solution
pause
DISPLAY AXES:0 0
# ---
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0
# --
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 4
plot 
$plotName = $name . $field . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0

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
  vertical scale factor 0.3
  plot contour lines (toggle)
#
exit
set view:0 -0.0969789 -0.00906344 0 1.02795 0.984808 0.0593912 -0.163176 -0.173648 0.336824 -0.925417 2.40985e-17 0.939693 0.34202
DISPLAY COLOUR BAR:0 0
x+:0
# x+:0
smaller 1.075
@comp = ( "Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "Ex error", "Ey error", "Ez error", "Hx error", "Hy error", "Hz error");
@compName = ( "Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "ExErr", "EyErr", "EzErr", "HxErr", "HyErr", "HzErr");
$cmd="#"; 
for( $i=0; $i<@comp; $i++ ){\
   $plotName = $name . $compName[$i] . "Rotated.ps"; \
   $cmd .= "\n plot:$comp[$i]\n hardcopy file name:0 $plotName\n hardcopy save:0";\
 }
$cmd


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
