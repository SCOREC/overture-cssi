#
# Channel grid for the periodic strip surface wave 
#
#
# usage: ogen [noplot] channelGrid -factor=<num> -order=[2/4/6/8] -xa= -xb= -ya= -yb= -prefix=<> -name=<> -periodic=[p|np|pn]
#
#  ogen -noplot channelGrid -prefix=case1ChannelGrid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=2
# 
# examples:
#
$prefix="channelGrid"; $xa=-1.; $xb=1.; $ya=-1.; $yb=1.;
$interp="e"; $interpType = "explicit for all grids";
$Nx=-1; # if set, use this many grid points in x 
$adjustCC=0; # apply cell-=centered adjustment for painted interface 
$order=2; $factor=1; $ds0=.1; # default values
$orderOfAccuracy = "second order"; $ng=2;  $periodic=""; $name=""; 
$numGhost=-1;  # if this value is set, then use this number of ghost points
$extraLines=0; 
# 
# get command line arguments
# Getopt::Long::Configure("prefix_pattern=(--channelGrid|--|-)");
GetOptions( "order=i"=>\$order,"factor=f"=>\$factor,"xa=f"=>\$xa,"xb=f"=>\$xb,"ds0=f"=>\$ds0,\
            "ya=f"=>\$ya,"yb=f"=>\$yb,"ybx=f"=>\$ybx,\
            "periodic=s"=>\$periodic,"name=s"=>\$name,"prefix=s"=>\$prefix,"numGhost=i"=> \$numGhost, \
	    "extraLines=i"=> \$extraLines,"Nx=i"=>\$Nx,"adjustCC=i"=>\$adjustCC );
# printf("channelGrid: factor=$factor xa=$xa xb=$xb ya=$ya yb=$yb ybx=$ybx\n");
# pause
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix="";
if( $periodic eq "p" ){ $suffix = "p"; }
if( $periodic eq "np" ){ $suffix = "np"; }
if( $periodic eq "pn" ){ $suffix = "pn"; }
$suffix .= ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; }
if( $name eq "" ){$name = $prefix . $interp . $factor . $suffix . ".hdf";}
# 
$ds=$ds0/$factor;
# 
create mappings
#
# left boundary grid
# 
rectangle
  set corners
    $xaLeft=$xa; $xbLeft=$xa + .25*($xb-$xa); 
    $xaLeft $xbLeft $ya $yb
  lines
    $nx = int( ($xbLeft-$xaLeft)/$ds +1.5 + $extraLines );
    $ny = int( ($yb-$ya)/$ds +1.5 + $extraLines );
    if( $Nx > 0 ){ $nx=$Nx; }
    $nx $ny
  boundary conditions
     if( $periodic eq "np" ){ $bc ="1 0 -1 -1"; }else{ $bc="1 0 3 4"; }
    $bc
  share   
     0 0 3 4
  mappingName
    leftGrid
exit
#
# middle grid
# 
$dsm=$ds/1.3; # change middle grid spacing 
rectangle
  set corners
    $xaMid=$xbLeft-$ds*.5; $xbMid=$xa + .75*($xb-$xa); 
    $xaMid $xbMid $ya $yb
  lines
    $nx = int( ($xbMid-$xaMid)/$dsm +1.5 + $extraLines );
    $ny = int( ($yb-$ya)/$dsm +1.5 + $extraLines );
    if( $Nx > 0 ){ $nx=$Nx; }
    $nx $ny
  boundary conditions
     if( $periodic eq "np" ){ $bc ="0 0 -1 -1"; }else{ $bc="0 0 3 4"; }
    $bc
  share   
     0 0 3 4
  mappingName
    middleGrid
exit
#
# right boundary grid
# 
rectangle
  set corners
    $xaRight=$xbMid-$ds*.5; $xbRight=$xb;
    $xaRight $xbRight $ya $yb
  lines
    $nx = int( ($xbRight-$xaRight)/$ds +1.5 + $extraLines );
    $ny = int( ($yb-$ya)/$ds +1.5 + $extraLines );
    if( $Nx > 0 ){ $nx=$Nx; }
    $nx $ny
  boundary conditions
     if( $periodic eq "np" ){ $bc ="0 2 -1 -1"; }else{ $bc="0 2 3 4"; }
    $bc
  share   
     0 0 3 4    
  mappingName
    rightGrid
exit
#
exit
generate an overlapping grid
    leftGrid
    middleGrid
    rightGrid
  done
  change parameters
    order of accuracy 
      $orderOfAccuracy
    interpolation type
      $interpType
    ghost points
      all
      $ngp=$ng+1; 
      $ng $ng $ng $ngp $ng $ng 
  exit
  #  display intermediate results
  # open graphics
  compute overlap
# 
  # display computed geometry
  exit
#
save an overlapping grid
$name
rectangle
exit

