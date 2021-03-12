#
#  comp compDiskMLA.cmd -name=compDiskMLAO2 -show=oneDiskMLA -solution=5
#  comp compDiskMLA.cmd -name=compDiskMLAO4 -show=oneDiskMLAO4 -solution=5
#
$solution=2; $name="compDiskMLA"; $show="oneDiskMLA"; 
# 
GetOptions( "solution=i"=>\$solution, "name=s"=>\$name, "show=s"=>\$show );
#      
output file name: $name.out
matlab file name: $name.m
#
specify files
  $show1 = $show . "G2.show"; 
  $show2 = $show . "G4.show"; 
  $show3 = $show . "G8.show"; 
  $show1
  $show2
  $show3
exit
choose a solution
  $solution
#
# Ex Ey Hz Px Py P1x P1y P2x P2y N0 N1 N2 N3
#  0  1  2  3  4  5  6   7   8   9  10 11 12
define a vector component
Ev
  0 1 
done
define a vector component
Pv
  3 4 
done
define a vector component
Nv
  9 10 11 12
done
#
compute errors



define a vector component
Hv
  3 4 5
done
# -- zero out errors near physical boundaries (for SUPER-GRID)
boundary error offset: 40 
compute errors  
#
plot differences
  plot contour lines (toggle)
  vertical scale factor 0.
  exit
  min max -0.005 0.005
  DISPLAY AXES:0 0
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  hardcopy file name:0 baBoxCylArrayDiffG128G64t1p5.ps
  hardcopy save:0
