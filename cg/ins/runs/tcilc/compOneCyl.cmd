#
#  Self convergence for one cylinder rotating in  square
#
#  comp compOneCyl.cmd -icase=oneCyl -sol=6
#  comp compOneCyl.cmd -icase=oneCylIMEX22 -sol=6
#  comp compOneCyl.cmd -icase=oneCylIMEX24 -sol=6 
#  comp compOneCyl.cmd -icase=oneCylIMEX24a -sol=6   [ larger dt
#  comp compOneCyl.cmd -icase=oneCylIMEX24b -sol=6   [ dt=.05, ...
#
#  comp compOneCyl.cmd -icase=oneCylIMEX44 -sol=6   (rampOrder=6)
#  comp compOneCyl.cmd -icase=oneCylIMEX44a -sol=6  (rampOrder=5)
#  comp compOneCyl.cmd -icase=oneCylIMEX44aa -sol=6  (bigger dt=.02,...)
#  comp compOneCyl.cmd -icase=oneCylIMEX44ab -sol=11  (bigger dt=.05,...)
#
#  comp compOneCyl.cmd -icase=oneCylIMEX44b -sol=6  (rampOrder=4)
#  comp compOneCyl.cmd -icase=oneCylIMEX44c -sol=6  (rampOrder=5, nu=.05)
#  comp compOneCyl.cmd -icase=oneCylIMEX44d -sol=6  (rampOrder=3)
# 
#  comp compOneCyl.cmd -icase=oneCylIMEX44s -sol=6  (stretched grid)
#
#  comp compOneCyl.cmd -icase=oneCylPC44 -sol=2
#
# PAPER: self-convergence for rotating cyl   -nu=.05, dt=[.05,.025,.0125]
#  comp compOneCyl.cmd -icase=oneCylIMEX44 -sol=11
#  comp compOneCyl.cmd -icase=oneCylIMEX44f -sol=11
#  comp compOneCyl.cmd -icase=oneCylIMEX24 -sol=11
#  comp compOneCyl.cmd -icase=oneCylIMEX24f -sol=11
#  comp compOneCyl.cmd -icase=oneCylIMEX22 -sol=11
#  comp compOneCyl.cmd -icase=oneCylIMEX22f -sol=11
#  comp compOneCyl.cmd -icase=oneCylPC44 -sol=11
#  comp compOneCyl.cmd -icase=oneCylPC44f -sol=11
#  comp compOneCyl.cmd -icase=oneCylPC22 -sol=11
#  comp compOneCyl.cmd -icase=oneCylPC22f -sol=11
#  
$files="pause"; $sol=-1; $icase="none"; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "icase=s"=>\$icase,"sol=s"=>\$sol );
#
if( $icase eq "oneCyl" ){ $files ="oneCyl2.show\n oneCyl4.show\n oneCyl8.show"; }
#
if( $icase eq "oneCylIMEX44" ){ $files ="oneCyl2IMEX44.show\n oneCyl4IMEX44.show\n oneCyl8IMEX44.show"; }
if( $icase eq "oneCylIMEX44f" ){ $files ="oneCyl4IMEX44.show\n oneCyl8IMEX44.show\n oneCyl16IMEX44.show"; }
if( $icase eq "oneCylIMEX22" ){ $files ="oneCyl2IMEX22.show\n oneCyl4IMEX22.show\n oneCyl8IMEX22.show"; }
if( $icase eq "oneCylIMEX22f" ){ $files ="oneCyl4IMEX22.show\n oneCyl8IMEX22.show\n oneCyl16IMEX22.show"; }
if( $icase eq "oneCylIMEX24" ){ $files ="oneCyl2IMEX24.show\n oneCyl4IMEX24.show\n oneCyl8IMEX24.show"; }
if( $icase eq "oneCylIMEX24f" ){ $files ="oneCyl4IMEX24.show\n oneCyl8IMEX24.show\n oneCyl16IMEX24.show"; }
if( $icase eq "oneCylPC44" ){ $files ="oneCyl2PC44.show\n oneCyl4PC44.show\n oneCyl8PC44.show"; }
if( $icase eq "oneCylPC44f" ){ $files ="oneCyl4PC44.show\n oneCyl8PC44.show\n oneCyl16PC44.show"; }
if( $icase eq "oneCylPC22" ){ $files ="oneCyl2PC22.show\n oneCyl4PC22.show\n oneCyl8PC22.show"; }
if( $icase eq "oneCylPC22f" ){ $files ="oneCyl4PC22.show\n oneCyl8PC22.show\n oneCyl16PC22.show"; }
#
if( $icase eq "oneCylIMEX24a" ){ $files ="oneCyl2IMEX24a.show\n oneCyl4IMEX24a.show\n oneCyl8IMEX24a.show"; }
if( $icase eq "oneCylIMEX24b" ){ $files ="oneCyl2IMEX24b.show\n oneCyl4IMEX24b.show\n oneCyl8IMEX24b.show"; }
if( $icase eq "oneCylIMEX44a" ){ $files ="oneCyl2IMEX44a.show\n oneCyl4IMEX44a.show\n oneCyl8IMEX44a.show"; }
if( $icase eq "oneCylIMEX44aa" ){ $files ="oneCyl2IMEX44aa.show\n oneCyl4IMEX44aa.show\n oneCyl8IMEX44aa.show"; }
if( $icase eq "oneCylIMEX44ab" ){ $files ="oneCyl2IMEX44ab.show\n oneCyl4IMEX44ab.show\n oneCyl8IMEX44ab.show"; }
if( $icase eq "oneCylIMEX44b" ){ $files ="oneCyl2IMEX44b.show\n oneCyl4IMEX44b.show\n oneCyl8IMEX44b.show"; }
if( $icase eq "oneCylIMEX44c" ){ $files ="oneCyl2IMEX44c.show\n oneCyl4IMEX44c.show\n oneCyl8IMEX44c.show"; }
if( $icase eq "oneCylIMEX44d" ){ $files ="oneCyl2IMEX44d.show\n oneCyl4IMEX44d.show\n oneCyl8IMEX44d.show"; }
if( $icase eq "oneCylIMEX44s" ){ $files ="oneCylStretched2IMEX44s.show\n oneCylStretched4IMEX44s.show\n oneCylStretched8IMEX44s.show"; }
if( $icase eq "oneCylPC44" ){ $files ="oneCyl2PC44.show\n oneCyl4PC44.show\n oneCyl8PC44.show"; }
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
$matlabFile = $icase . "sol$sol.m"; 
matlab file name: $matlabFile
compute errors

