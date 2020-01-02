#
# 
#   plotStuff plotSolution3d.cmd -show=baMI3dG8.show -name=baMI3dG8 -solution=5
#   plotStuff plotSolution3d.cmd -show=baMI3dG4.show -name=baMI3dG4 -solution=5
# 
#   plotStuff plotSolution3d.cmd -show=baDieSphere.show -name=baDieSphere -solution=5 -emax=1.5
#
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
pause 
solution: $solution
pause
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 4
plot 
$plotName = $name . "EfieldNorm.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
