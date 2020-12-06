#
# Convert $mapName to a nurbs and rotate and/or shift  
#   Input:
#     $angle
#     $xShift, $yShift 
#     $ng = number of ghost points 
#
$numGhost=$ng+1; 
sub convertToNurbs\
{ local($old,$new,$angle,$xShift,$yShift)=@_; \
  $commands = "nurbs (surface)\n" . \
              "interpolate from mapping with options\n" . "$old\n" . "parameterize by index (uniform)\n" . \
              " number of ghost points to include\n $numGhost\n" . \
              "done\n" . \
              "rotate\n" . "$angle 1\n" . "0 0 0\n" . \
              "shift\n" . "$xShift $yShift\n" . \
              "mappingName\n" . "$new\n" . "exit\n"; \
}
convertToNurbs($mapName,$gridName,$angle,$xShift,$yShift);
$commands
