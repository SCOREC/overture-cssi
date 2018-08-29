#
#   DISK IN A CHANNEL -- SELF CONVERGENCE
#
#   comp compDisk.cmd -solution=16 -logFile=compDiskInAChannelScf10t1p5.log
#
$solution=1; $name="diskInAChannel"; $logFile="compDiskInAChannelScf10.log"; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"solution=i"=>\$solution,"logFile=s"=>\$logFile  );
# 
specify files
  diskInAChannelScf10G2.show
  diskInAChannelScf10G4.show
  diskInAChannelScf10G8.show
exit
#
output file name: $logFile
#
# --- solid domain ------
choose a frame series (domain) per file
 0
 0
 0
define a vector component
  usv
  6 7 
done
define a vector component
  vsv
  0 1
done
define a vector component
  sigmav
  2 3 4 5
done
choose a solution
  $solution
compute errors
#
# --- fluid domain ------
choose a frame series (domain) per file
 1
 1
 1
delete all vector components
# -- currently only vector components are saved if any are specified
define a vector component
  pp
  0
done
define a vector component
  vv
  1 2
done
choose a solution
  $solution
compute errors