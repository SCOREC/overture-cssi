#
# 
#   plotStuff plotSolution3d.cmd -show=baMI3dG8.show -name=baMI3dG8 -solution=5
# 
#   plotStuff plotSolution3d.cmd -show=baMatInt3DG4.show -name=baMatInt3D -solution=5
#   plotStuff plotSolution3d.cmd -show=baMatInt3DG4Diss.show -name=baMatInt3D -solution=5
#   plotStuff plotSolution3d.cmd -show=baMatInt3DG8Diss.show -name=baMatInt3D -solution=11
#
#   plotStuff plotSolution3d.cmd -show=baMatInt3DG8O4.show -name=baMatInt3DO4 -solution=11
# 
#   plotStuff plotSolution3d.cmd -show=baDieSphere.show -name=baDieSphere -solution=5 -emax=1.5
#
#  BA 3D Many boxes:
#   plotStuff plotSolution3d.cmd -show=baBox27G8.show -name=baBox27G8 -solution=5 
#   plotStuff plotSolution3d.cmd -show=baBox27G16.show -name=baBox27G16 -solution=7 -emax=.9 
#
$show="baMI3d.show.hdf"; $solution="-1"; $name="plot"; $field="Ey"; $emax=1.2; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field,"emax=f"=>\$emax );
#
$show
show forcing regions 1
#
derived types
  E field norm
exit
contour
  plot:eFieldNorm
  delete contour plane 0
  contour lines 0
  min max 0 $emax
exit
DISPLAY AXES:0 0
set view:0 -0.0879154 0.00302115 0 1.03115 0.766044 0.219846 -0.604023 0 0.939693 0.34202 0.642788 -0.262003 0.719846
x- 
solution: $solution
pause
# 
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 4
plot 
$plotName = $name . "EfieldNorm.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0



# Plot lots of things: 
@comp = ( "Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "Ex error", "Ey error", "Ez error", "Hx error", "Hy error", "Hz error");
@compName = ( "Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "ExErr", "EyErr", "EzErr", "HxErr", "HyErr", "HzErr");
$cmd="#"; 
for( $i=0; $i<@comp; $i++ ){\
   $plotName = $name . $compName[$i] . ".ps"; \
   $cmd .= "\n plot:$comp[$i]\n hardcopy file name:0 $plotName\n hardcopy save:0";\
 }
# printf("cmd=$cmd\n");
$cmd
