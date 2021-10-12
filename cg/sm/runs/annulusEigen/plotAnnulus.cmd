#
# plotStuff plotAnnulus
#
$show="annulusMode1.show";
$show="annulus4Mode1g.show"; 
#
$show
#
contour
  plot:s11
  adjust grid for displacement 1
  set plot bounds
    # -.75 .75 -.75 .75
    # -.85 .85 -.85 .85
    -1.1 1.1 -1.1 1.1 
  displacement scale factor .25
  exit
x-:0
x-:0
bigger 1.1
DISPLAY AXES:0 0
# 
solution 1
plot:s11
hardcopy file name:0 annulusS11t0p0.ps
hardcopy save:0
pause
# 
solution 7
plot:s12
hardcopy file name:0 annulusS12t0.6p0.ps
hardcopy save:0
pause
# 
solution 21
plot:s22
hardcopy file name:0 annulusS22t2p0.ps
hardcopy save:0



displacement
  displacement scale factor .5
  exit this menu
displacement
  exit this menu
