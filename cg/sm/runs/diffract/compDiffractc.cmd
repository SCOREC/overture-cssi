#
#   (comp for SOS-C results)
# $ovg/bin/comp compDiffractc.cmd
# $ovg/bin/comp -noplot compDiffractc.cmd
#
# $ovg/bin/comp compDiffractc.cmd -option=noAmr -outputFile="compDiffractSOSnoAmr.log"
# $ovg/bin/comp compDiffractc.cmd -option=amr -outputFile="compDiffractSOSAmr.log"
#
# -- 110416 : redo for smog revision
# $gf/comp -noplot compDiffractc.cmd -option=noAmr -outputFile="compDiffractSOSnoAmrNew.log"
# $gf/comp -noplot compDiffractc.cmd -option=amr   -outputFile="compDiffractSOSAmrNew.log"
#
$option="noAmr"; $outputFile="compDiffractSOSnoAmr.log"; 
# $option="amr"; $outputFile="compDiffractSOSAmr.log"; 
GetOptions( "option=s"=>\$option,"outputFile=s"=>\$outputFile );
#
output file name
  $outputFile
specify files (coarse to fine)
# if( $option eq "noAmr" ){ $files = "diffract16c.show\n diffract32c.show\n diffract64c.show"; }\
#                    else{ $files = "diffract8l2r2.show\n diffract8l2r4.show\n diffract64c.show"; }
# NEW
if( $option eq "noAmr" ){ $files = "diffract16cf.show\n diffract32cf.show\n diffract128cf.show"; }\
                    else{ $files = "diffract8l2r2f4.show\n diffract8l2r4f4.show\n diffract128cf.show"; }
#if( $option eq "noAmr" ){ $files = "diffract16cf.show\n diffract32cf.show\n diffract64cf.show"; }\
#                    else{ $files = "diffract8l2r2f4.show\n diffract8l2r4f4.show\n diffract64cf.show"; }
# if( $option eq "noAmr" ){ $files = "diffract16cf.show\n diffract32cf.show\n diffract128cf.show"; }\
#                     else{ $files = "diffract8l2r2f4.show\n diffract32cf.show\n diffract64cf.show"; }
# if( $option eq "noAmr" ){ $files = "diffract16cf.show\n diffract32cf.show\n diffract128cf.show"; }
$files
exit
# 
choose a solution
# 9 : t=1.6 
9
compute errors
# 3
define a vector component
uNorm
  0 1 
 done
define a vector component
vNorm
  4 5 
 done
define a vector component
sNorm
  6 7 8
 done
compute errors

exit
