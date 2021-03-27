#
# Plot solutions from the multi-level atomic model 
# 
#   plotStuff plotMLA.cmd -show=fourDiskMLA.show -name=fourDiskMLA -solution=6 -timeLabel="t1p0"  [tp=.2]
#
#  plotStuff plotMLA.cmd -show=64DiskMLA.show -name=64DiskMLA -emax=.5 -solution=21 -timeLabel="t10"    [ tp=.5 ]
#
#  plotStuff plotMLA.cmd -show=rpiMLA.show -name=rpiMLA -emax=-1 -solution=7 -timeLabel="t3"  [tp=.5]
#  plotStuff plotMLA.cmd -show=rpiMLAO4.show -name=rpiMLAO4 -emax=-1 -solution=7 -timeLabel="t3"  [tp=.5]
#
#  plotStuff plotMLA.cmd -show=oneDiskMLAO4G8.show -name=oneDiskMLA -emax=-1 -solution=5 -timeLabel="t2"  [tp=.5]
#
$show="fourDiskMLA.hdf"; $solution="-1"; $name="plot"; $field="Ey"; $emin=0; $emax=-1; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
$timeLabel=""; 
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,"timeLabel=s"=>\$timeLabel,\
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
$plotName = $name . "EfieldNorm$timeLabel.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
pause
plot:Py
$plotName = $name . "Py$timeLabel.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
plot:N3
$plotName = $name . "N3$timeLabel.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
plot:N0
$plotName = $name . "N0$timeLabel.ps"; 
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