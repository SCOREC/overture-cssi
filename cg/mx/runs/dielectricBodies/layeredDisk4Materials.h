#
#  Specify material files for the array of 4 solid disks with one layer
#
#          eps=8 |  eps=4
#         -------------------
#          eps=4 |  eps=8
#
if( $outerDomain eq "" ){ $outerDomain="outerDomain"; } # default name of outer-domain
$baseDomainName = "innerDomain";
#  
$d=0; $disk=1; # $d=domain counter, $disk=disk counter
# Outer-domain:
$domainName[$d]="$outerDomain"; $npv[$d]=0; $dmFile[$d]="baIsoEps1.txt"; $d=$d+1;
# 
# $numLayers=1; # fix me 
for( $i2=0; $i2<2; $i2++ ){\
for( $i1=0; $i1<2; $i1++ ){\
  $domainName[$d]="$baseDomainName$disk";  $npv[$d]=0; $dmFile[$d]="baIsoEps2.txt"; \
  if( $nm eq "multilevelAtomic" ){ $npv[$d]=2; $dmFile[$d]="mlaMat4levels.txt"; } \
  $d=$d+1;   \
  for( $layer=1; $layer<=$numLayers; $layer++ )\
  {\
    $layerDomainName = "disk$disk" . "Layer$layer" . "Domain"; \
    $domainName[$d]="$layerDomainName";   $npv[$d]=0; \
    if( $layer == 1 ){ $dmFile[$d]="baIsoEps3.txt"; }else{ $dmFile[$d]="baIsoEps4.txt"; }    \
    if( $layer == 1  && $dm eq "GDM" ){ $npv[$d]=1; $dmFile[$d]="gdmMaterial1.txt"; } \
    $d=$d+1;                     \
  }\
 $disk=$disk+1; \
}}



