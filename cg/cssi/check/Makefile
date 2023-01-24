# Type "make check" to run the regression tests
# NOTE: The file tests contains the list of test cases

common = ../../common
check: buildGrids
	ln -sf $(common)/src/check.p .
	./check.p -solver=cgcssi
check-quiet: buildGrids
	ln -sf $(common)/src/check.p .
	./check.p -solver=cgcssi -verbose=0

cgcssi:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) cgcssi

# build grids that are used by the regression tests (and not already in Overture/sampleGrids)
buildGrids: sicFixede2.order2.hdf multifluidChannelGrid.hdf


sicFixede2.order2.hdf:
	${Overture}/bin/ogen -noplot ${Overture}/sampleGrids/sicArg -order=2 -fixedRadius=.25 -interp=e -factor=2

multifluidChannelGrid.hdf:
	${Overture}/bin/ogen -noplot ../runs/cavityCollapse/multifluidChannelGrid.cmd

