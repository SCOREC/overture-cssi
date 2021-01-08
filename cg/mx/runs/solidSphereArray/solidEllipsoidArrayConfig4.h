#
# Configuration file for solidEllipsoidArrayGrid.cmd
#
#    4 ellipsoids 
#     
# Box bounds: 
$xa=-1.5; $xb=1.5; $ya=-1.5; $yb=1.5; $za=-.75; $zb=.75;
#
# ------ list of ellipsoid parameters: --------
#   Note: probably best to have c=longest semi-axis so end patches are smaller
#
$numEllipsoids=4;
$i=0;
$ae[$i]= .5;   $be[$i]= .5; $ce[$i]=.5; $xe[$i]=-0.65; $ye[$i]=-.65; $ze[$i]=0.; $rot1[$i]=  0; $rot2[$i]=  0; $i=$i+1; 
$ae[$i]= .375; $be[$i]= .5; $ce[$i]=.5; $xe[$i]= 0.65; $ye[$i]=-.65; $ze[$i]=0.; $rot1[$i]= 20; $rot2[$i]=  0; $i=$i+1; 
$ae[$i]= .55;  $be[$i]= .4; $ce[$i]=.4; $xe[$i]=-0.65; $ye[$i]=0.65; $ze[$i]=0.; $rot1[$i]= 45; $rot2[$i]= 45; $i=$i+1;  
$ae[$i]= .55;  $be[$i]= .4; $ce[$i]=.4; $xe[$i]= 0.65; $ye[$i]=0.65; $ze[$i]=0.; $rot1[$i]=  0; $rot2[$i]=  0; $i=$i+1;  
#
