#
# Configuration file for solidEllipsoidArrayGrid.cmd
#
#    8 ellipsoids 
#     
# Box bounds: 
$xa=-1.5; $xb=1.5; $ya=-1.5; $yb=1.5; $za=-1.5; $zb=1.5;
#
# ------ list of ellipsoid parameters: --------
#   Note: probably best to have c=longest semi-axis so end patches are smaller
#
$numEllipsoids=8;
$i=0;
$de=.65; # delta between ellipsoids 
# flat ellipsoid 
$ae[$i]= .55;  $be[$i]= .5; $ce[$i]=.40; $xe[$i]=-$de; $ye[$i]=-$de; $ze[$i]=-$de; $rot1[$i]=  0; $rot2[$i]= 90; $i=$i+1; 
$ae[$i]= .375; $be[$i]= .5; $ce[$i]=.5;  $xe[$i]= $de; $ye[$i]=-$de; $ze[$i]=-$de; $rot1[$i]= 20; $rot2[$i]=  0; $i=$i+1; 
$ae[$i]= .55;  $be[$i]= .4; $ce[$i]=.4;  $xe[$i]=-$de; $ye[$i]= $de; $ze[$i]=-$de; $rot1[$i]= 45; $rot2[$i]= 45; $i=$i+1;  
$ae[$i]= .55;  $be[$i]= .4; $ce[$i]=.4;  $xe[$i]= $de; $ye[$i]= $de; $ze[$i]=-$de; $rot1[$i]=  0; $rot2[$i]=  0; $i=$i+1;  
#
# thin ellipsoid
$ae[$i]= .375; $be[$i]= .375; $ce[$i]=.5;  $xe[$i]=-$de; $ye[$i]=-$de; $ze[$i]= $de; $rot1[$i]=  0; $rot2[$i]=  0; $i=$i+1; 
$ae[$i]= .375; $be[$i]= .5;   $ce[$i]=.5;  $xe[$i]= $de; $ye[$i]=-$de; $ze[$i]= $de; $rot1[$i]= 20; $rot2[$i]=  0; $i=$i+1; 
# flat ellipsoid 
$ae[$i]= .50;  $be[$i]= .5;   $ce[$i]=.35; $xe[$i]=-$de; $ye[$i]= $de; $ze[$i]= $de; $rot1[$i]=  0; $rot2[$i]=  0; $i=$i+1;  
$ae[$i]= .55;  $be[$i]= .4;   $ce[$i]=.4;  $xe[$i]= $de; $ye[$i]= $de; $ze[$i]= $de; $rot1[$i]=  0; $rot2[$i]=  0; $i=$i+1;  
