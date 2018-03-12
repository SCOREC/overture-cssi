#
# plotStuff plotGrid.cmd
#
$grid="tubeGride1.order2.hdf";
#
$grid
#
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
#
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  line width scale factor:0 3
#
  PIC:brass
#  grid colour 1 BRASS
  grid colour 1 TURQUOISE
  grid colour 0 MEDIUMBLUE
  set view:0 0.0259138 0.00320958 0 1.19434 0.642788 0.262003 -0.719846 0 0.939693 0.34202 0.766044 -0.219846 0.604023
#
  hardcopy file name:0 tubeGrid1.ps
  hardcopy save:0
