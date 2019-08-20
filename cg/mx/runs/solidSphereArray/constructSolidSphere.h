#  -------------- Construct a SOLID SPHERE ---------------
# 
$xSphere=$xv[$count];
$ySphere=$yv[$count];
$zSphere=$zv[$count];
$count = $count+1; 
# 
include buildSolidSphere.h
# 
$domainName[$count] = "sphereDomain$count"; 
$domainGridNames[$count] = $innerNames; 
$domainGridNames[0] .= "\n" . $outerNames; 
