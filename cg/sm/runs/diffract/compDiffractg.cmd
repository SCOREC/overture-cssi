#
#   (comp for FOS-G results)
#
# $ovg/bin/comp compDiffractg.cmd -option=noAmr -outputFile="compDiffractFOSnoAmr.log"
# $ovg/bin/comp compDiffractg.cmd -option=amr -outputFile="compDiffractFOSAmr.log"
#
# Use limiter:
# $ovg/bin/comp compDiffractg.cmd -lim=1 -option=noAmr -outputFile="compDiffractFOSnoAmrLim.log"
# $ovg/bin/comp compDiffractg.cmd -lim=1 -option=amr -outputFile="compDiffractFOSAmrLim.log"
#
# -- 110416 : redo for smog revision
#  $gf/comp -noplot compDiffractg.cmd -option=noAmr -outputFile="compDiffractFOSnoAmrNew.log"
#  $gf/comp -noplot compDiffractg.cmd -option=amr -outputFile="compDiffractFOSAmrNew.log"
#  $gf/comp -noplot compDiffractg.cmd -option=amr -outputFile="compDiffractFOSAmrNew.log" -sol=5
#
$option="noAmr"; $outputFile="compDiffractFOSnoAmr.log"; $lim=0; $sol=9; 
# $option="amr"; $outputFile="compDiffractFOSAmr.log"; 
GetOptions( "option=s"=>\$option,"outputFile=s"=>\$outputFile,"lim=i"=>\$lim,"sol=i"=>\$sol );
#
output file name
  $outputFile
specify files (coarse to fine)
# if( $option eq "noAmr" ){ $files = "diffract16g.show\n diffract32g.show\n diffract64g.show"; }\
#                     else{ $files = "diffract8l2r2g.show\n diffract8l2r4g.show\n diffract64g.show"; }
# New: use finer fine grid
if( $option eq "noAmr" ){ $files = "diffract16g.show\n diffract32g.show\n diffract128g.show"; }\
                    else{ $files = "diffract8l2r2g.show\n diffract8l2r4g.show\n diffract128g.show"; }
# TEMP: use coarse grids: 
# if( $option eq "noAmr" ){ $files = "diffract8g.show\n diffract16g.show\n diffract64g.show"; }
#
if( $lim ne 0 ){\
 if( $option eq "noAmr" ){ $files = "diffract16glim.show\n diffract32glim.show\n diffract64glim.show"; }\
                     else{ $files = "diffract8l2r2glim.show\n diffract8l2r4glim.show\n diffract64glim.show"; } }
$files
exit
# 
choose a solution
# 9 : t=1.6 
$sol
# 3
define a vector component
uNorm
  6 7 
 done
define a vector component
vNorm
  0 1 
 done
define a vector component
sNorm
  2 3 4 5 
 done
compute errors

exit
