#
# plotStuff plotSolidFuelAssemblyGrids.cmd -show=solidFuelAssemblyGrid3pinse2.order2.hdf -name=solidFuelAssembly7RodGrid
# plotStuff plotSolidFuelAssemblyGrids.cmd -show=sfaOneRodHeight1Gride2.order2.hdf -name=solidFuelAssembly1RodH1Grid
# plotStuff plotSolidFuelAssemblyGrids.cmd -show=sfa19RodsGride2.order2.hdf -name=solidFuelAssembly19RodsH2Grid
#
$show = "/home/henshaw.0/Overture/ogen/solidFuelAssembly3pinsl4e2.hdf"; 
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "opt=s"=>\$opt );
#
# 
$show
# 
  colour boundaries by domain
  plot block boundaries 0
  coarsening factor 2
  set view:0 0.00225445 -0.0407627 0 1.17399 0.939693 -0.17101 0.296198 0.34202 0.469846 -0.813798 -2.09427e-17 0.866025 0.5
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  #
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  line width scale factor:0 3
  plot
  $plotName = $name . ".ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0


  x-r 90
  y-r:0
  y-r:0
  x+r:0
  x+r:0
  x+r:0

  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  #
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  hardcopy file name:0 solidFuelAssembly7RodGrid.ps
  line width scale factor:0 3
  plot
  hardcopy file name:0 solidFuelAssembly7RodGrid.ps
  hardcopy save:0



  colour boundaries by domain
  plot block boundaries 0
  x-r 90
  y-r:0
  y-r:0
  x+r:0
  x+r:0
  x+r:0
  plot grid lines 0
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  bigger 1.05
  hardcopy file name:0 solidFuelAssemblyDomains.ps
  hardcopy save:0
  plot grid lines 1
  coarsening factor 1
  coarsening factor 4
  hardcopy file name:0 solidFuelAssemblyGrids.ps
  colour boundaries by grid number
  colour boundaries by domain
  hardcopy file name:0 solidFuelAssemblyGrids.ps
  hardcopy save:0
  colour boundaries by grid number
  hardcopy file name:0 solidFuelAssemblyGrids2.ps
  hardcopy save:0
  plot grid lines 0
  plot boundary condition (toggle) 5 0
  plot boundary condition (toggle) 5 1
  plot boundary condition (toggle) 3 0
  plot boundary condition (toggle) 2 0
  plot boundary condition (toggle) 4 0
  plot boundary condition (toggle) 3 1
  plot boundary condition (toggle) 2 1
  pick colour...
  PIC:brass
  grid colour 0 BRASS
  grid colour 1 BRASS
  grid colour 2 BRASS
  grid colour 3 BRASS
  grid colour 4 BRASS
  grid colour 5 BRASS
  grid colour 6 BRASS
  grid colour 7 BRASS
  grid colour 8 BRASS
  grid colour 9 BRASS
  grid colour 10 BRASS
  grid colour 11 BRASS
  grid colour 12 BRASS
  grid colour 13 BRASS
  grid colour 14 BRASS
  grid colour 15 BRASS
  grid colour 16 BRASS
  grid colour 17 BRASS
  grid colour 18 BRASS
  grid colour 19 BRASS
  grid colour 20 BRASS
  grid colour 21 BRASS
  grid colour 22 BRASS
  grid colour 23 BRASS
  PIC:blue
  grid colour 1 BLUE
  grid colour 23 BLUE
  PIC:turquoise
  grid colour 23 TURQUOISE
  hardcopy file name:0 solidFuelAssemblyDomains2.ps
  hardcopy save:0
