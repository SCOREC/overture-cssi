#
#  Specify material files for the array of 64 solid MLA disks
#  Disks are numbered left to right, bottom to top
#
#       57  58                      64
#       49                          56
#       41                          48
#       33                          40
#       25                          32
#       17                          24
#        9  10  11  12  13  14  15  16
#        1   2   3   4   5   6   7   8
#
if( $outerDomain eq "" ){ $outerDomain="outerDomain"; } # default name of outer-domain
$baseDomainName = "innerDomain";
#  
$i=0; # (outer domain is domain 0)
# Outer-domain:
$domainName[$i]="$outerDomain"; $npv[$i]=0; $dmFile[$i]="baIsoEps1.txt"; $i=$i+1;
#
for( $i2=0; $i2<8; $i2++ ){  \
for( $i1=0; $i1<8; $i1++ )\
{\
  $i1m = $i1 % 4; $i2m = $i2 % 4; \
  $dmFile[$i]="baIsoEps1.txt"; \
  if( ($i1m+$i2m) % 4  == 0 ) { $npv[$i]=2; $dmFile[$i]="mlaMat4levels.txt"; }  \
  if( ($i1m+$i2m) % 4  == 1 ) { $npv[$i]=2; $dmFile[$i]="mlaMat4levels.txt"; }   \
  if( ($i1m+$i2m) % 4  == 2 ) { $npv[$i]=2; $dmFile[$i]="mlaMat4levels.txt"; } \
  if( ($i1m+$i2m) % 4  == 3 ) { $npv[$i]=2; $dmFile[$i]="mlaMat4levels.txt"; } \
  $domainName[$i]="$baseDomainName$i";  $i=$i+1; \
}}
