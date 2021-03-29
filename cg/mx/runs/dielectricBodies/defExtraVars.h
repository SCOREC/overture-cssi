# This is included in mxScript.cmd using -extraCommandsStart=defExtraVars.h 
if( $numLayers eq "" ){ $numLayers=1; }
#
GetOptions( "numLayers=i"=>\$numLayers );
# printf("numLayers=$numLayers\n");
# pause
