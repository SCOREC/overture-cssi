#
#  plotStuff plotDeformingChannel -show=<name> -option=[contour|grid]
#
#  -- slipWalls : elastic piston
#  plotStuff plotDeformingChannel -show=dcG4.show -option=grid -name=deformingChannelG4
#
#  -- noSlipWalls on ends 
#  plotStuff plotDeformingChannel -show=dcPinnedG4.show
#  plotStuff plotDeformingChannel -show=dcPinnedG4.show -option=grid -name=deformingChannelPinnedG4
#
$show="dc.show"; $option="contour"; $name="dc"; 
* ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show, "name=s"=>\$name,"solution=i"=>\$solution,"vorMin=f"=>\$vorMin,"vorMax=f"=>\$vorMax,\
            "option=s"=>\$option, "name=s"=>\$name, "matlab=s"=>\$matlab );
#
#
$show
#
# -- plot grids 
if( $option eq "grid" ){ $cmd="frame series:fluidDomain\n grid\n exit this menu\n frame series:solidDomain\n displacement\n exit this menu";}else{ $cmd="#"; }
$cmd
#
pause
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
DISPLAY AXES:0 0
DISPLAY SQUARES:0 0
line width scale factor:0 3
#  plotGrids(num,timeLabel)
sub plotGrids\
{ local($num,$label)=@_; \
  $plotName = $name . "t$label" . "Grids.ps"; \
  $cmds = "solution: $num \n" . \
   "hardcopy file name:0 $plotName\n" . \
   "hardcopy save:0\n"; \
}
plotGrids(1,"0p0");
$cmds
plotGrids(6,"1p0");
$cmds
plotGrids(11,"2p0");
$cmds
plotGrids(16,"3p0");
$cmds
plotGrids(21,"4p0");
$cmds
plotGrids(26,"5p0");
$cmds
plotGrids(31,"6p0");
$cmds
plotGrids(36,"7p0");
$cmds
plotGrids(41,"8p0");
$cmds
plotGrids(46,"9p0");
$cmds
plotGrids(51,"10p0");
$cmds


#
contour
  vertical scale factor 0.
  plot:v2
  adjust grid for displacement 1
  exit
frame series:fluidDomain
contour
  vertical scale factor 0.
  plot:v
exit

next
next
next
next
next
