#
#   Scattering of a dipsersive plane wave from many disks
#
#   comp compMultiDisk.cmd -order=2 -solution=7 -logFile=compMultiDiskOrder2.log
# 
#   comp compMultiDisk.cmd -order=4 -solution=3 -logFile=compMultiDiskOrder4.log
#
$solution=1; $name="compMultiDisk"; $logFile="compMultiDiskOrder2.log"; $order=2; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"solution=i"=>\$solution,"logFile=s"=>\$logFile,"order=i"=>\$order  );
# 
specify files
  if( $order eq 2 ){ $files="MD2.show\n MD4.show\n  MD8.show"; }else{  $files="MD2a.show\n  MD4a.show\n  MD8a.show"; }
  $files 
exit
#
output file name: $logFile
#
choose a solution
  $solution
#
compute errors
