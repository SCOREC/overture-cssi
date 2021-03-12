#
#  Specify material files for the array of 16 solid disks
#
if( $outerDomain eq "" ){ $outerDomain="outerDomain"; } # default name of outer-domain
$baseDomainName = "innerDomain";
#  
$i=0; # (outer domain is domain 0)
# Outer-domain:
$domainName[$i]="$outerDomain"; $npv[$i]=0; $dmFile[$i]="baIsoEps1.txt"; $i=$i+1;
# 
for( $i2=0; $i2<4; $i2++ ){\
for( $i1=0; $i1<4; $i1++ ){\
  if( ($i2 % 4) == 0 ){ \
    if( ($i1 % 4) == 0 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1; } \
    if( ($i1 % 4) == 1 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1; } \
    if( ($i1 % 4) == 2 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps6.txt"; $i=$i+1; } \
    if( ($i1 % 4) == 3 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps3.txt"; $i=$i+1; } \
  }\
  elsif( ($i2 % 4) == 1 ) { \
    if( ($i1 % 4) == 0 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1; } \
    if( ($i1 % 4) == 1 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1; } \
    if( ($i1 % 4) == 2 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps3.txt"; $i=$i+1; } \
    if( ($i1 % 4) == 3 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps6.txt"; $i=$i+1; } \
  }\
  elsif( ($i2 % 4) == 2 ) { \
    if( ($i1 % 4) == 0 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps6.txt"; $i=$i+1; } \
    if( ($i1 % 4) == 1 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps3.txt"; $i=$i+1; } \
    if( ($i1 % 4) == 2 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1; } \
    if( ($i1 % 4) == 3 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1; } \
  }\
  else { \
    if( ($i1 % 4) == 0 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps3.txt"; $i=$i+1; } \
    if( ($i1 % 4) == 1 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps6.txt"; $i=$i+1; } \
    if( ($i1 % 4) == 2 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1; } \
    if( ($i1 % 4) == 3 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1; } \
  }\
} \
}


