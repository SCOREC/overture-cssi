#
# Configuration file for solidEllipsoidArrayGrid.cmd
#
#    1 rotated ellipsoid
#     
# Box bounds: 
$xa=-1.; $xb=1.; $ya=-.75; $yb=.75; $za=-.75; $zb=.75;
#
# ------ list of ellipsoid parameters: --------
#   Note: probably best to have c=longest semi-axis so end patches are smaller
#
$numEllipsoids=1;
$i=0;
$ae[$i]= .5; $be[$i]= .4; $ce[$i]=.45; $xe[$i]=0; $ye[$i]=0.; $ze[$i]=0.; $rot1[$i]= 0; $rot2[$i]= 0; $i=$i+1; 
# $ae[$i]= .5; $be[$i]= .375; $ce[$i]=.65; $xe[$i]=0; $ye[$i]=0.; $ze[$i]=0.; $rot1[$i]= 20; $rot2[$i]= 20; $i=$i+1; 
#
