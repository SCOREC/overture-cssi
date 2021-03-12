#  -*-perl-*-
#
# Define an array of box shaped regions for the BA MX solver
#   Input:
#      $nxBox, $nyBox, $nzBox    : number of boxes in each direction
#      $xWidth, $yWidth, $zWidth : width of boxes 
#      $xSep, $ySep, $zSep       : separation between the boxes
#      $xc, $yc, $zc             : lower left corner of the boxes
#      $matFileArray[$mat]       : list of material .txt files (periodic wrap if there are more boxes than materials)
#
# NOTE: materials are re-used for different boxes, if possible, to reduce the number of GDM vectors
#
$numRegions=$nxBox*$nyBox*$nzBox;
if( $xSep eq "" ){ $xSep=0; }
if( $ySep eq "" ){ $ySep=0; }
if( $zSep eq "" ){ $zSep=0; }
#  use materialFile values if provided, not $dmFile[0] = background material 
if( $materialFile ne "" ){ @matFileArray = (); for($i=0; $i<$#dmFile; $i++ ){ $matFileArray[$i]=$dmFile[$i+1]; } }
#
$numMaterials= $#matFileArray+1; # number materials supplied
#
# NOTE: For efficiency keep track on different boxes that belong to the same material
#   We thus first make a list of all boxes belonging to a given material
#   We also search previous materials to see if a material can be re-used
#
# $numBox[$mat] : counts number of boxes for each material
for( $mat=0; $mat<$numMaterials; $mat++ ){ $numBox[$mat]=0; }  
#
# Assign boxes to materials 
#   $matBox[$mat][$j] = list of boxes associated with material $mat
#
$i=0;  # counts boxes
for( $iz=0; $iz<$nzBox; $iz++ ){\
for( $iy=0; $iy<$nyBox; $iy++ ){\
for( $ix=0; $ix<$nxBox; $ix++ ){\
  $xa[$i]=$xc+$ix*($xWidth+$xSep); $ya[$i]=$yc+$iy*($yWidth+$ySep); $za[$i]=$zc+$iz*($zWidth+$zSep); \
  $xb[$i]=$xa[$i]+$xWidth;         $yb[$i]=$ya[$i]+$yWidth;         $zb[$i]=$za[$i]+$zWidth; \
  $mat = $i % $numMaterials; \
  for( $j=0; $j<$i; $j++ ){ if( $matFileArray[$j] eq $matFileArray[$mat] ){ $mat=$j; last; }} \
  $matBox[$mat][$numBox[$mat]]=$i; $numBox[$mat]++;    \
  $i++;  \
}}}
#
$cmd="";
for( $mat=0; $mat<$numMaterials; $mat++ ){\
  if( $numBox[$mat]>0 ){\
    $cmd .= "material file: $matFileArray[$mat]\n";  \
    for( $j=0; $j<$numBox[$mat]; $j++ ){ \
       $i = $matBox[$mat][$j]; \
       $cmd .=  "box: $xa[$i] $xb[$i] $ya[$i] $yb[$i] $za[$i] $zb[$i] (xa,xb,ya,yb,za,zb)\n"; \
    } \
  } \
}
$cmd .="continue";
# printf("cmd=$cmd\n");
#  
forcing options...
define material region...
# Now execute the commands: 
$cmd 
#

