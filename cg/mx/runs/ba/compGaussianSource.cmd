#
#  comp compGaussianSource.cmd -option=squareIso -label="BAMX: Gaussian source 2D square, super-grid"
#  comp compGaussianSource.cmd -option=squareGdm -label="BAMX: Gaussian source 2D square, baMatGDMgenND, super-grid"
#
#  comp compGaussianSource.cmd -option=boxIso -label="BAMX: Gaussian source 3D box, baMatIsoEps1, super-grid"
#  comp compGaussianSource.cmd -option=boxGdm -label="BAMX: Gaussian source 3D box, baMatGDMgenND, super-grid"
#
$option="square"; $sgWidth=.25; $nd=2;
$label="BAMX: Gaussian source"; 
GetOptions( "option=s"=>\$option, "label=s"=>\$label, "sgWidth=f"=>\$sgWidth  );
# 
output file name: comp.out
#
caption label: $label
#
specify files
 $files="pause";
 if( $option eq "squareIso" ){ $files="baMatIsoEps1G8.show\n baMatIsoEps1G16.show\n baMatIsoEps1G32.show\n #  baMatIsoEps1G64.show "; }
 if( $option eq "squareGdm" ){ $files="baMatGDMgenNDG8.show\n baMatGDMgenNDG16.show\n baMatGDMgenNDG32.show\n baMatGDMgenNDG64.show "; }
 if( $option eq "boxIso" ){ $nd=3; $files="baMatIsoEps1BoxG4.show\n baMatIsoEps1BoxG8.show\n baMatIsoEps1BoxG16.show"; }
 if( $option eq "boxGdm" ){ $nd=3; $files="baMatGDMgenNDBoxG4.show\n baMatGDMgenNDBoxG8.show\n baMatGDMgenNDBoxG16.show "; }
#  if( $option eq "boxGdm" ){ $nd=3; $files="baMatGDMgenNDBoxG4.show\n baMatGDMgenNDBoxG8.show\n baMatGDMgenNDBoxG16.show\n baMatGDMgenNDBoxG32.show "; }
 $files 
exit
choose a solution
  21
#
define a vector component
Ev
  if( $nd eq "2" ){ $components="0 1"; }else{ $components="0 1 2"; }
  $components
  # 0 1 2
done
#
# -- zero out errors near physical boundaries (for SUPER-GRID)
# boundary error offset: 40
absorbing layer width: .25 
#
compute errors  
