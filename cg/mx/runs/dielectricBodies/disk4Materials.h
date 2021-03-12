#
#  Specify material files for the array of 4 solid disks
#
#          eps=8 |  eps=4
#         -------------------
#          eps=4 |  eps=8
#
if( $outerDomain eq "" ){ $outerDomain="outerDomain"; } # default name of outer-domain
$baseDomainName = "innerDomain";
#  
$i=0; # (outer domain is domain 0)
# Outer-domain:
$domainName[$i]="$outerDomain"; $npv[$i]=0; $dmFile[$i]="baIsoEps1.txt"; $i=$i+1;
# 
for( $i2=0; $i2<2; $i2++ ){\
for( $i1=0; $i1<2; $i1++ ){\
  if( ($i2 % 2) == 0 ){ \
    if( ($i1 % 2) == 0 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1; } \
    if( ($i1 % 2) == 1 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1; } \
  }\
  else { \
    if( ($i1 % 2) == 0 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps8.txt"; $i=$i+1; } \
    if( ($i1 % 2) == 1 ){ $domainName[$i]="$baseDomainName$i"; $npv[$i]=0; $dmFile[$i]="baIsoEps4.txt"; $i=$i+1; } \
  }\
} \
}


