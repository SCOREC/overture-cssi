#
# Two co-centric annulii for a two-domain example (e.g. for conjugate heat transfer)
#
# NOTE: THIS VERSION HAS THE ANNLUUS WITH r1=r and r2=theta (reversed from normal)
#
# usage: ogen [noplot] doubleAnnulusGrid-factor=<num> -order=[2/4/6/8] -interp=[e/i]
# 
# examples:
#     ogen -noplot doubleAnnulusGrid -order=2 -interp=e -factor=0.5
#     ogen -noplot doubleAnnulusGrid -order=2 -interp=e -factor=1
#     ogen -noplot doubleAnnulusGrid -order=2 -interp=e -factor=2
#     ogen -noplot doubleAnnulusGrid -order=2 -interp=e -factor=4
#     ogen -noplot doubleAnnulusGrid -order=2 -interp=e -factor=8
#
# Double quarter annulus
#    ogen -noplot doubleAnnulusGrid -prefix=doubleQuarterAnnulusGrid -theta2=.25 -order=2 -interp=e -factor=1
#
$prefix="doubleAnnulusGrid"; 
$theta1=0; $theta2=1.; # range of theta (normalized to [0,1])
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=>\$factor,"interp=s"=>\$interp,"prefix=s"=> \$prefix,"theta1=f"=>\$theta1,"theta2=f"=>\$theta2 );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
$name = $prefix . "$interp$factor" . $suffix . ".hdf";
# 
$ds=.1/$factor;
$pi=atan2(1.,1.)*4; 
# 
# Here is the radius of the circular boundary:
$innerRad=.5;
$middleRad=1.;
$outerRad=1.5;
# 
create mappings
#
# --- INNER ANNULUS ---
#
Annulus
 # keep the number of radial points on the annulus fixed:
  inner and outer radii
    $innerRad $middleRad
  angles: $theta1, $theta2
  lines
    $nTheta = int( 2.*$pi*($theta2-$theta1)*($middleRad)/$ds + 1.5 );
    $nr = int( ($middleRad-$innerRad)/$ds + 2.5 );
    $nTheta $nr
  mappingName
   innerAnnulus0
exit
# FLIP r1 <-> r2 
reparameterize
  transform which mapping?
    innerAnnulus0
  reorient domain coordinates
    1 0
  boundary conditions
    if( ($theta2-$theta1) == 1 ){ $cmd="1 100 -1 -1 "; }else{ $cmd="1 100 3 4"; }
    $cmd
  share
    0 100 0 0 
  mappingName
   innerAnnulus    
exit
#
# ---- OUTER ANNULUS ----
Annulus
 # keep the number of radial points on the annulus fixed:
  inner and outer radii
    $middleRad $outerRad
  angles: $theta1, $theta2    
  lines
    $nr = int( ($outerRad-$middleRad)/$ds + 2.5 );
    $nTheta $nr
  mappingName
   outerAnnulus0
exit
# FLIP r1 <-> r2 
reparameterize
  transform which mapping?
    outerAnnulus0
  reorient domain coordinates
    1 0
  boundary conditions
    if( ($theta2-$theta1) == 1 ){ $cmd="100 2 -1 -1 "; }else{ $cmd="100 2 3 4"; }
    $cmd
  share
    100 2  0 0    
  mappingName
   outerAnnulus 
exit
#
exit
generate an overlapping grid
    innerAnnulus
    outerAnnulus
  done
  change parameters
 # choose implicit or explicit interpolation
    specify a domain
 # domain name:
      outerDomain 
 # grids in the domain:
      outerAnnulus
      done
    specify a domain
 # domain name:
      innerDomain 
 # grids in the domain:
      innerAnnulus
      done
# 
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
  exit
#  display intermediate results
 # open graphics
  compute overlap
# 
  exit
#
save an overlapping grid
$name
doubleAnnulus
exit
