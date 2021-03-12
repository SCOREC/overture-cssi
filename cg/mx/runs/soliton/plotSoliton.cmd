#
#  plot the soliton solution
#
# plotStuff plotSoliton.cmd -show=solitonO2 -name=solitonO2 -solution=21
#
# plotStuff plotSoliton.cmd -show=solitonO4 -name=solitonO4 -solution=21
#
# fine grid:
# plotStuff plotSoliton.cmd -show=solitonO4G4 -name=solitonO4G4 -solution=21
#
# Curved channel:
#  plotStuff plotSoliton.cmd -show=curvedChannelSoliton -name=curvedChannelSoliton -solution=61
#
$show="solitonO4.hdf"; $solution="-1"; $name="plot"; $field="Ey"; $emin=0; $emax=-1; $tPlot=10;
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,"tPlot=f"=>\$tPlot,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax );
#
$show
#
derived types
  E field norm
exit
solution: $solution
#
contour
  plot:eFieldNorm
  plot contour lines (toggle)
#
  line plots 
    specify a boundary
      4001 0 0 1
#    specify lines
#      1 4001
#      0. .1 400 .1
    Ey
    add Py
    add N0
    add eFieldNorm
pause    
    # x0 holds the x values 
    add x0
#       
    save results to a matlab file
      $time = ($solution-1)*$tPlot;; 
      $matlab = $name . "t$time" . ".m"; 
      $matlab
    exit this menu
  exit this menu
exit
exit


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

