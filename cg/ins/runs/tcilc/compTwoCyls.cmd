#
#  Self convergence for flow past two cyls
#
#  IMEX22, nu=0.1, : ramped p 
#  comp compTwoCyls.cmd -icase=IMEX22nu0p1 -sol=3     
#  IMEX22, nu=0.1, : ramped v 
#  comp compTwoCyls.cmd -icase=IMEX22nu0p1v -sol=3    
# 
#  IMEX44, nu=0.1, : ramped p 
#  comp compTwoCyls.cmd -icase=IMEX44nu0p1 -sol=2
#  comp compTwoCyls.cmd -icase=IMEX44nu0p1a -sol=3     (reduced dtMax)
# 
#  IMEX44, nu=0.1, : No-BEWNO ramped p 
#  comp compTwoCyls.cmd -icase=IMEX44nu0p1nb -sol=3
#
#  IMEX24, nu=0.1, : BEWNO ramped p 
#  comp compTwoCyls.cmd -icase=IMEX24nu0p1 -sol=3
#
#  IMEX44, nu=0.1, : ramped v 
#  comp compTwoCyls.cmd -icase=IMEX44nu0p1v -sol=2
# 
# -- ONE CYL
#  comp compTwoCyls.cmd -icase=oneCylIMEX22nu0p1 -sol=3
#  comp compTwoCyls.cmd -icase=oneCylIMEX44nu0p1 -sol=3
#  comp compTwoCyls.cmd -icase=oneCylPC44nu0p1 -sol=3
#
$files="pause"; $sol=-1; $icase="none"; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "icase=s"=>\$icase,"sol=s"=>\$sol );
#
if( $icase eq "IMEX22nu0p1" ){ $files ="twoCyl1IMEX2nu0p1.show\n twoCyl2IMEX2nu0p1.show\n twoCyl4IMEX2nu0p1.show\n twoCyl8IMEX2nu0p1.show"; }
if( $icase eq "IMEX22nu0p1v" ){ $files ="twoCyl1IMEX2nu0p1v.show\n twoCyl2IMEX2nu0p1v.show\n twoCyl4IMEX2nu0p1v.show\n twoCyl8IMEX2nu0p1v.show"; }
#
## if( $icase eq "IMEX44nu0p1" ){ $files ="twoCyl1IMEX44nu0p1.show\n twoCyl2IMEX44nu0p1.show\n twoCyl4IMEX44nu0p1.show\n twoCyl8IMEX44nu0p1.show\n twoCyl16IMEX44nu0p1.show"; }
if( $icase eq "IMEX44nu0p1" ){ $files ="twoCyl4IMEX44nu0p1.show\n twoCyl8IMEX44nu0p1.show\n twoCyl16IMEX44nu0p1.show"; }
#
if( $icase eq "IMEX44nu0p1a" ){ $files ="twoCyl1IMEX44nu0p1a.show\n twoCyl2IMEX44nu0p1a.show\n twoCyl4IMEX44nu0p1a.show\n twoCyl8IMEX44nu0p1a.show"; }
#
if( $icase eq "IMEX24nu0p1" ){ $files ="twoCyl1IMEX24nu0p1.show\n twoCyl2IMEX24nu0p1.show\n twoCyl4IMEX24nu0p1.show\n twoCyl8IMEX24nu0p1.show"; }
if( $icase eq "IMEX44nu0p1nb" ){ $files ="twoCyl1IMEX44nu0p1nb.show\n twoCyl2IMEX44nu0p1nb.show\n twoCyl4IMEX44nu0p1nb.show\n twoCyl8IMEX44nu0p1nb.show"; }
if( $icase eq "IMEX44nu0p1v" ){ $files ="twoCyl1IMEX44nu0p1v.show\n twoCyl2IMEX44nu0p1v.show\n twoCyl4IMEX44nu0p1v.show\n twoCyl8IMEX44nu0p1v.show"; }
#
if( $icase eq "oneCylIMEX22nu0p1" ){ $files ="oneCyl2IMEX2nu0p1.show\n oneCyl4IMEX2nu0p1.show\n oneCyl8IMEX2nu0p1.show\n oneCyl16IMEX2nu0p1.show"; }
#
#3 if( $icase eq "oneCylIMEX44nu0p1" ){ $files ="oneCyl2IMEX44nu0p1.show\n oneCyl4IMEX44nu0p1.show\n oneCyl8IMEX44nu0p1.show"; }
if( $icase eq "oneCylIMEX44nu0p1" ){ $files ="oneCyl2IMEX44nu0p1.show\n oneCyl4IMEX44nu0p1.show\n oneCyl8IMEX44nu0p1.show\n oneCyl16IMEX44nu0p1.show"; }
#
if( $icase eq "oneCylPC44nu0p1" ){ $files ="oneCyl2PC44nu0p1.show\n oneCyl4PC44nu0p1.show\n oneCyl8PC44nu0p1.show"; }
# 
specify files
 $files
exit
choose a solution
  $sol
#
$output = $icase . ".log"; 
output file name: $outputFile
#
$matlabFile = $icase . ".m"; 
matlab file name: $matlabFile
compute errors

