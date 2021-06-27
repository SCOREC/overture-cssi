#
#  Plot the grid for the heat exchanger - pipe in a box
#
#    plotStuff plotGrid3d.cmd -show=heatExchangerGride2.order2.hdf -name=pipeInABox
#
$show=""; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show, "name=s"=>\$name );
# -------------------------------------------------------------------------------------------------
$show
#
DISPLAY AXES:0 0
set view:0 0.000216774 -0.00711932 0 0.900403 0.866025 0.321394 -0.383022 -0.5 0.55667 -0.663414 6.93889e-17 0.766044 0.642788
#
  hardcopy file name:0 heatExchangerGrid.ps
  hardcopy save:0


# 
  toggle grid 0 0
  pick closest 1
  grid colour 1 BRASS
  toggle grid 4 0
  toggle grid 5 0
  toggle grid 6 0
  toggle grid 7 0
  grid colour 3 COPPER
  grid colour 2 CORNFLOWERBLUE
  plot block boundaries 0
  colour block boundaries black
  set view:0 0.00154683 -0.00309366 0 2.00772 -0.766044 0.321394 -0.55667 0 0.866025 0.5 0.642788 0.383022 -0.663414
  add coordinate plane 2 1 7 (grid dir index)
  add coordinate plane 1 1 28 (grid dir index)
# pause
  add coordinate plane 1 1 15 (grid dir index)
# pause
#  add coordinate plane 1 1 37 (grid dir index)
#pause
  add coordinate plane 1 0 12 (grid dir index)
# pause
  add coordinate plane 2 0 7 (grid dir index)
#
$plotName = "$name.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
