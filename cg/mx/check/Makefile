#
# Type "make check" to run the regression tests
# Note: the main set of regression tests are found in file "tests"
# 
# Example: To run one check:
#    check.p twoSquaresInterfaceGDM22r
#    check.p twoBoxesInterface32r
#
#
common = ../../common

check:  buildGrids
	ln -sf $(common)/src/check.p .
	./check.p -solver=cgmx

check-quiet:  buildGrids
	ln -sf $(common)/src/check.p .
	./check.p -solver=cgmx -verbose=0

# Interface checks alone
checkInterface: buildGrids
	ln -sf $(common)/src/check.p .
	./check.p -solver=cgmx -testsFile=interfaceTests

# ------ build grids used by check -----------------
tube1.order4.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/tube1.order4.cmd

square16.order4.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/square16.order4.cmd

cici1.order4.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/cicArg -factor=1 interp=i -order=4

box16.order4.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/box16.order4.cmd

box1.order2.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/boxArg.cmd -order=2 -factor=1 

nonBox1.order2.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/boxArg.cmd -prefix=nonBox -order=2 -factor=1

box1.order4.ng3.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/boxArg.cmd  -order=4 -factor=1 -numGhost=3

nonBox1.order4.ng3.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/boxArg.cmd -prefix=nonBox -order=4 -factor=1 -numGhost=3

innerOuter4.order4.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/innerOuter -factor=4 -order=4 -deltaRad=.5 -interp=e -name="innerOuter4.order4.hdf"

# ---- twoSquaresInterface grids -----
twoSquaresInterfacee1.order2.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoSquaresInterface -order=2 -interp=e -factor=1

twoSquaresInterfacee2.order2.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoSquaresInterface -order=2 -interp=e -factor=2

twoNonSquaresInterfacee1.order2.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoSquaresInterface -prefix=twoNonSquaresInterface -angle=0.0 -order=2 -interp=e -factor=1

twoNonSquaresInterfacee2.order2.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoSquaresInterface -prefix=twoNonSquaresInterface -angle=0.0 -order=2 -interp=e -factor=2

twoSquaresInterfacee1.order4.ng3.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoSquaresInterface -order=4 -interp=e -factor=1 -numGhost=3

twoSquaresInterfacee2.order4.ng3.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoSquaresInterface -order=4 -interp=e -factor=2 -numGhost=3

twoNonSquaresInterfacee1.order4.ng3.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoSquaresInterface -prefix=twoNonSquaresInterface -angle=0.0 -order=4 -interp=e -factor=1 -numGhost=3

twoNonSquaresInterfacee2.order4.ng3.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoSquaresInterface -prefix=twoNonSquaresInterface -angle=0.0 -order=4 -interp=e -factor=2 -numGhost=3

# ---- twoBoxesInterface grids -----
twoBoxesInterfacee1.order2.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoBoxesInterface -order=2 -interp=e -factor=1

twoBoxesInterfacee2.order2.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoBoxesInterface -order=2 -interp=e -factor=2

twoNonBoxesInterfacee1.order2.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoBoxesInterface -prefix=twoNonBoxesInterface -angle=0.0 -order=2 -interp=e -factor=1

twoNonBoxesInterfacee2.order2.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoBoxesInterface -prefix=twoNonBoxesInterface -angle=0.0 -order=2 -interp=e -factor=2

twoBoxesInterfacee1.order4.ng3.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoBoxesInterface -order=4 -interp=e -factor=1 -numGhost=3

twoBoxesInterfacee2.order4.ng3.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoBoxesInterface -order=4 -interp=e -factor=2 -numGhost=3

twoNonBoxesInterfacee1.order4.ng3.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoBoxesInterface -prefix=twoNonBoxesInterface -angle=0.0 -order=4 -interp=e -factor=1 -numGhost=3

twoNonBoxesInterfacee2.order4.ng3.hdf:
	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoBoxesInterface -prefix=twoNonBoxesInterface -angle=0.0 -order=4 -interp=e -factor=2 -numGhost=3



square20.order2.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/squareArg -order=2 -nx=20 

nonSquare20.order2.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/nonSquare -order=2 -n=20 

# nore ghost points for sosup: 
square20.order4.ng3.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/squareArg -order=4 -nx=20 -numGhost=3

nonSquare20.order4.ng3.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/nonSquare -order=4 -n=20 -numGhost=3

square20.order6.ng4.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/squareArg -order=6 -nx=20 -numGhost=4

solidSphereInABoxe2.order2.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/solidSphereInABox -order=2 -interp=e -factor=2

solidSphereInABoxe2.order4.ng3.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/solidSphereInABox -order=4 -interp=e -numGhost=3 -factor=2


# build grids that are used by the regression tests (and not already in Overture/sampleGrids)
buildGrids: tube1.order4.hdf square16.order4.hdf cici1.order4.hdf box16.order4.hdf nonBox1.order4.ng3.hdf box1.order4.ng3.hdf \
            innerOuter4.order4.hdf \
            square20.order4.ng3.hdf square20.order6.ng4.hdf \
            solidSphereInABoxe2.order2.hdf solidSphereInABoxe2.order4.ng3.hdf square20.order2.hdf \
            nonSquare20.order2.hdf nonSquare20.order4.ng3.hdf box1.order2.hdf nonBox1.order2.hdf \
            twoSquaresInterfacee1.order4.ng3.hdf twoSquaresInterfacee2.order4.ng3.hdf \
            twoNonSquaresInterfacee1.order4.ng3.hdf twoNonSquaresInterfacee2.order4.ng3.hdf \
            twoBoxesInterfacee1.order4.ng3.hdf twoBoxesInterfacee2.order4.ng3.hdf \
            twoNonBoxesInterfacee1.order4.ng3.hdf twoNonBoxesInterfacee2.order4.ng3.hdf 
	touch buildGrids


checkCheckFiles:	
	@cd $(common); $(MAKE) -$(MAKEFLAGS) $(BINDIR)checkCheckFiles

cgmx:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) cgmx

clean:
	rm -f *.hdf



