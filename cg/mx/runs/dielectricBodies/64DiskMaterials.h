#
#  Specify material files for the array of 64 solid disks
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
  if( ($i1m+$i2m) % 4  == 0 ) { $dmFile[$i]="baIsoEps4.txt"; }  \
  if( ($i1m+$i2m) % 4  == 1 ) { $dmFile[$i]="baIsoEps8.txt"; }   \
  if( ($i1m+$i2m) % 4  == 2 ) { $dmFile[$i]="baIsoEps2.txt"; } \
  if( ($i1m+$i2m) % 4  == 3 ) { $dmFile[$i]="baIsoEps5.txt"; } \
  $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $i=$i+1; \
}}
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps4.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps8.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps2.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps5.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps4.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps8.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps2.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps5.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps5.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps2.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps8.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps4.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps5.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps2.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps8.txt";  $i=$i+1; \
#-   $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps4.txt";  $i=$i+1; \
#- 
#
#- # 
#- # non-dispersive materials: 
#- # $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps2.txt"; $i=$i+1;
#- # $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps3.txt"; $i=$i+1;
#- #
#- # dispersive materials
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1;
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1;
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1;
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1;
#- #
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1;
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1;
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1;
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1;
#- #
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1;
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1;
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1;
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1;
#- #
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1;
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1;
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1;
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1;
#- #
#- for( $j=0; $j<200; $j++ )\
#- {\
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=1; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1; \
#- $domainName[$i]="$baseDomainName$i"; $npv[$i]=2; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1; \
#- }
