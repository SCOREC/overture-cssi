#
# Configuration file for solidEllipsoidArrayGrid.cmd
#
#    3 ellipsoids in a line
#     
#  Ellipsoid volume: (4*pi/3) * a * b * c 
#
# Box bounds: 
$xa=-2.5; $xb=2.5; $ya=-.75; $yb=.75; $za=-.75; $zb=.75;
#
# ------ list of ellipsoid parameters: --------
#   Note: probably best to have c=longest semi-axis so end patches are smaller
#
$numEllipsoids=3;
$i=0;
# thin 
$ae[$i]= .375;   $be[$i]= .375; $ce[$i]=.5;   $xe[$i]=-1.50; $ye[$i]=0.; $ze[$i]=0.; $rot1[$i]=  0; $rot2[$i]= 90; $i=$i+1; 
# flat
$ae[$i]= .5;     $be[$i]= .5;   $ce[$i]=.35;  $xe[$i]=  0.0; $ye[$i]=0.; $ze[$i]=0.; $rot1[$i]=  0; $rot2[$i]= 90; $i=$i+1; 
$ae[$i]= .55;    $be[$i]= .4;   $ce[$i]=.4;   $xe[$i]= 1.50; $ye[$i]=0.; $ze[$i]=0.; $rot1[$i]=  0; $rot2[$i]=  0; $i=$i+1;  
#
