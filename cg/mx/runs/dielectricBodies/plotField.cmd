#
#   plotStuff plotField.cmd -show=ellipseG8.show -name=ellipseG8
#   plotStuff plotField.cmd -show=rodG8.show -name=rodG8
#   plotStuff plotField.cmd -show=diskInBoxG8.show -name=diskInBoxG8
#   plotStuff plotField.cmd -show=crossG16.show -name=cross
#
#   plotStuff plotField.cmd -show=starFishG128.show -name=starFish
#   plotStuff plotField.cmd -show=starFishMovie.show -name=starFish
# 
#   plotStuff plotField.cmd -show=starFisha.show -name=starFish -solution=21
# 
#   plotStuff plotField.cmd -show=waveyDiskb.show -name=waveyDiskb -solution=21
#
#    plotStuff plotField.cmd -show=dieBlockG8.show -solution=101 -name=dieBlock
#    plotStuff plotField.cmd -show=dieBlockG8.show -solution=1  -name=dieBlockt0p0
#    plotStuff plotField.cmd -show=dieBlockG8.show -solution=11 -name=dieBlockt1p0
#    plotStuff plotField.cmd -show=dieBlockG8.show -solution=31 -name=dieBlockt3p0
#    plotStuff plotField.cmd -show=dieBlockG8.show -solution=51 -name=dieBlockt5p0
# 
#    plotStuff plotField.cmd -show=dieBlockG8Eps4.show -solution=101 -name=dieBlockG8Eps4
#
#    plotStuff plotField.cmd -show=pecDiskG4Eps4.show -solution=201 -name=pecDiskG4Eps4
#
#    plotStuff plotField.cmd -show=rod16kx0p5.show -solution=201 -name=rod16kx0p5
#    plotStuff plotField.cmd -show=rod16kx0p25.show -solution=201 -name=rod16kx0p25
# 
#    plotStuff plotField.cmd -show=cross16kx0p5.show -solution=201 -name=cross16kx0p5
#    plotStuff plotField.cmd -show=cross16kx0p25.show -solution=201 -name=cross16kx0p25
#
#    plotStuff plotField.cmd -show=ellipse16kx0p5.show -solution=201 -name=ellipse16kx0p5
#    plotStuff plotField.cmd -show=ellipse8kx0p25.show -solution=201 -name=ellipse8kx0p25
#
# compare solutions for L=20, 20 and 4 (kx=.5)
#    plotStuff plotField.cmd -show=ellipseL20A30G8.show -solution=41 -name=ellipseL20A30
#    plotStuff plotField.cmd -show=ellipseL10A30G8.show -solution=21 -name=ellipseL10A30
#    plotStuff plotField.cmd -show=ellipseL4A30G8.show -solution=21 -name=ellipseL4A30
#
# compare solutions for L=20, 20 and 4 (kx=.25)
#    plotStuff plotField.cmd -show=ellipseL20A30G8kx0p25.show -solution=41 -name=ellipseL20A30kx0p25
#    plotStuff plotField.cmd -show=ellipseL10A30G8kx0p25.show -solution=41 -name=ellipseL10A30kx0p25
#    plotStuff plotField.cmd -show=ellipseL4A30G8kx0p25.show -solution=41 -name=ellipseL4A30kx0p25
#
#   plotStuff plotField.cmd -show=ellGL4gp.show -solution=21 -name=ellGL4gp
#   plotStuff plotField.cmd -show=ellGL10gp.show -solution=21 -name=ellGL10gp
#
# Long domain G8 reference solution: ( ellipse-scat)
#   plotStuff plotField.cmd -show=ellipseL20A30G8.show -name=ellL20G8 -tSave=1 -numPerTime=10 -numToSave=5
#
# L2 NL-RBC short-domain  compared to long G8
#  plotStuff plotField.cmd -show=ellL2G4NLvL20G8.show -name=ellL2nl -tSave=10 -numPerTime=20 -numToSave=6 -field="Ey error"
# L4 NL-RBC short-domain  compared to long G8
#  plotStuff plotField.cmd -show=ellL4G4NLvL20G8.show -name=ellL4nl -tSave=10 -numPerTime=20 -numToSave=6 -field="Ey error"
# L2  PML  compared to long G8
#  plotStuff plotField.cmd -show=ellL2G4PMLvL20G8.show -name=ellL2pml -tSave=10 -numPerTime=20 -numToSave=6 -field="Ey error"
# L4 PML  compared to long G8
#  plotStuff plotField.cmd -show=ellL4G4PMLvL20G8.show -name=ellL4pml -tSave=10 -numPerTime=20 -numToSave=6 -field="Ey error" 
# L10 PML  compared to long G8
#  plotStuff plotField.cmd -show=ellL10G4PMLvL20G8.show -name=ellL10pml -tSave=10 -numPerTime=20 -numToSave=6 -field="Ey error" 
# L4 PML41  compared to long G8
#  plotStuff plotField.cmd -show=ellL4G4PML41vL20G8.show -name=ellL4pml41 -tSave=10 -numPerTime=20 -numToSave=6 -field="Ey error" 
# L10 PML41  compared to long G8
#  plotStuff plotField.cmd -show=ellL10G4PML41vL20G8.show -name=ellL10pml41 -tSave=10 -numPerTime=20 -numToSave=6 -field="Ey error"
# L10 PML81  compared to long G8
#  plotStuff plotField.cmd -show=ellL10G4PM81LvL20G8.show -name=ellL10pml81 -tSave=10 -numPerTime=20 -numToSave=6 -field="Ey error" 
#   +++ Order 2:
# L2 NL-RBC short-domain  compared to long G8
#  plotStuff plotField.cmd -show=ellL2G4NLOrder2vL20G8.show -name=ellL2nlOrder2 -tSave=10 -numPerTime=20 -numToSave=6 -field="Ey error"
# L4 NL-RBC short-domain  compared to long G8
#  plotStuff plotField.cmd -show=ellL4G4NLOrder2vL20G8.show -name=ellL4nlOrder2 -tSave=10 -numPerTime=20 -numToSave=6 -field="Ey error"
# L4 PML  compared to long G8
#  plotStuff plotField.cmd -show=ellL4G4PMLOrder2vL20G8.show -name=ellL4pmlOrder2 -tSave=10 -numPerTime=20 -numToSave=6 -field="Ey error"
# L10 PML41  compared to long G8
#  plotStuff plotField.cmd -show=ellL10G4PML41Order2vL20G8.show -name=ellL10pml41Order2 -tSave=10 -numPerTime=20 -numToSave=6 -field="Ey error"
#
$show="ellipseG8.hdf"; $solution="-1"; $name="plot"; $field="Ey"; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field );
#
$show
#
plot:$field
contour
  plot contour lines (toggle)
  coarsening factor 1 (<0 : adaptive)
  # min max -1 1
 vertical scale factor 0.
exit
pause
#
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 3
DISPLAY AXES:0 0
#
#
solution: $solution
pause
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
DISPLAY AXES:0 0
#
  line width scale factor:0 3
  plot
  $plotName = $name . "Ey.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
# 
plot:Ex
  $plotName = $name . "Ex.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0





$fieldNoBlanks=$field; $fieldNoBlanks =~ s/ //g; 
$cmd="#";
for( $i=0; $i<$numToSave; $i++ ){\
  $solution = $numPerTime*$i +1; $t = $i*$tSave; $t =~ s/\./p/g;  \
  $cmd .= "solution: $solution\n"; \
  $plotName = $name . $fieldNoBlanks . "t$t.ps"; \
  $cmd .= "plot:$field\n"; \
  $cmd .= "hardcopy file name:0 $plotName\n"; \
  $cmd .= "hardcopy save:0\n"; \
}
$cmd 
# save a matlab file with the errors
plot sequence:errors
  Ey_error
  add Ex_error
  save results to a matlab file
  $matlabFile = $name . "Err.m"; 
  $matlabFile
# 



plot:Hz
  $plotName = $name . "Hz.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
# 
erase
stream lines
  streamline density 100
  arrow size 0.05
exit
  $plotName = $name . "SL.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
