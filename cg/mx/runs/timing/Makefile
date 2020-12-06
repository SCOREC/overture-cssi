
boxPerf:
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=box -tf=1. -scheme=FD    -order=2 -numRes=5
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=box -tf=1. -scheme=FDA   -order=2 -numRes=5
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=box -tf=1. -scheme=UPC   -order=2 -numRes=5
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=box -tf=1. -scheme=SOSUP -order=2 -numRes=5
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=box -tf=1. -scheme=FD    -order=4 -numRes=4
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=box -tf=1. -scheme=FDA   -order=4 -numRes=4
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=box -tf=1. -scheme=UPC   -order=4 -numRes=4
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=box -tf=1. -scheme=SOSUP -order=4 -numRes=4

nonBoxPerf:
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonBox -tf=1. -scheme=FD    -order=2 -numRes=5
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonBox -tf=1. -scheme=FDA   -order=2 -numRes=5
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonBox -tf=1. -scheme=UPC   -order=2 -numRes=5
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonBox -tf=1. -scheme=SOSUP -order=2 -numRes=5
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonBox -tf=1. -scheme=FD    -order=4 -numRes=4
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonBox -tf=1. -scheme=FDA   -order=4 -numRes=4
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonBox -tf=1. -scheme=UPC   -order=4 -numRes=4
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonBox -tf=1. -scheme=SOSUP -order=4 -numRes=4

squarePerf:
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=square -tf=1. -scheme=FD    -order=2 -numRes=8
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=square -tf=1. -scheme=FDA   -order=2 -numRes=8
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=square -tf=1. -scheme=UPC   -order=2 -numRes=8
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=square -tf=1. -scheme=SOSUP -order=2 -numRes=8
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=square -tf=1. -scheme=FD    -order=4 -numRes=4
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=square -tf=1. -scheme=FDA   -order=4 -numRes=4
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=square -tf=1. -scheme=UPC   -order=4 -numRes=4
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=square -tf=1. -scheme=SOSUP -order=4 -numRes=4

nonSquarePerf:
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonSquare -tf=1. -scheme=FD    -order=2 -numRes=6
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonSquare -tf=1. -scheme=FDA   -order=2 -numRes=6
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonSquare -tf=1. -scheme=UPC   -order=2 -numRes=6
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonSquare -tf=1. -scheme=SOSUP -order=2 -numRes=6
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonSquare -tf=1. -scheme=FD    -order=4 -numRes=4
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonSquare -tf=1. -scheme=FDA   -order=4 -numRes=4
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonSquare -tf=1. -scheme=UPC   -order=4 -numRes=4
	runCgmx.p -runCase=performance -cmd=boxEigen.cmd -g=nonSquare -tf=1. -scheme=SOSUP -order=4 -numRes=4


square256:
	runCgmx.p -cmd=boxEigen.cmd -g=square256.order2 -tf=.2 -matlabFile="cpuSquare256Order2" -title="Square 256, Order 2"
	runCgmx.p -cmd=boxEigen.cmd -g=square256.order4.ng3 -tf=.2 -matlabFile="cpuSquare256Order4" -title="Square 256, Order 4"
	runCgmx.p -cmd=boxEigen.cmd -g=nonSquare256.order2 -tf=.2 -matlabFile="cpuNonSquare256Order2" -title="NonSquare 256, Order 2"
	runCgmx.p -cmd=boxEigen.cmd -g=nonSquare256.order4.ng3 -tf=.2 -matlabFile="cpuNonSquare256Order4" -title="NonSquare 256, Order 4"

square1024:
	runCgmx.p -cmd=boxEigen.cmd -g=square1024.order2 -tf=.05 -matlabFile="cpuSquare1024Order2" -title="Square 1024, Order 2"
	runCgmx.p -cmd=boxEigen.cmd -g=square1024.order4.ng3 -tf=.05 -matlabFile="cpuSquare1024Order4" -title="Square 1024, Order 4"
	runCgmx.p -cmd=boxEigen.cmd -g=nonSquare1024.order2 -tf=.05 -matlabFile="cpuNonSquare1024Order2" -title="NonSquare 1024, Order 2"
	runCgmx.p -cmd=boxEigen.cmd -g=nonSquare1024.order4.ng3 -tf=.05 -matlabFile="cpuNonSquare1024Order4" -title="NonSquare 1024, Order 4"



box4:
	runCgmx.p -cmd=boxEigen.cmd -g=box4.order2        -tf=.2 -matlabFile="cpuBox40Order2"    -title="Box 40, Order 2"
	runCgmx.p -cmd=boxEigen.cmd -g=nonBox4.order2     -tf=.2 -matlabFile="cpuNonBox40Order2" -title="NonBox 40, Order 2"
	runCgmx.p -cmd=boxEigen.cmd -g=box4.order4.ng3    -tf=.2 -matlabFile="cpuBox40Order4"    -title="Box 40, Order 4"
	runCgmx.p -cmd=boxEigen.cmd -g=nonBox4.order4.ng3 -tf=.2 -matlabFile="cpuNonBox40Order4" -title="NonBox 40, Order 4"


box16:
	runCgmx.p -cmd=boxEigen.cmd -g=box16.order2 -tf=.2 -matlabFile="cpuBox160Order2" -title="Box 160, Order 2"
	runCgmx.p -cmd=boxEigen.cmd -g=box16.order4.ng3 -tf=.2 -matlabFile="cpuBox160Order4" -title="Box 160, Order 4"




nonBox16:
	runCgmx.p -cmd=boxEigen.cmd -g=nonBox16.order2 -tf=.2 -matlabFile="cpuNonBox160Order2" -title="NonBox 160, Order 2"
	runCgmx.p -cmd=boxEigen.cmd -g=nonBox16.order4.ng3 -tf=.2 -matlabFile="cpuNonBox160Order4" -title="NonBox 160, Order 4"
