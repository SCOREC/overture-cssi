
ogen = ${Overture}/bin/ogen

grids:
	$(ogen) CreateGrid.cmd -factor=1
	$(ogen) CreateGrid.cmd -factor=2
	$(ogen) CreateGrid.cmd -factor=4
	$(ogen) CreateGrid.cmd -factor=8
	$(ogen) radialElasticPistonGrid -interp=e -factor=1
	$(ogen) radialElasticPistonGrid -interp=e -factor=2
	$(ogen) radialElasticPistonGrid -interp=e -factor=4
	$(ogen) radialElasticPistonGrid -interp=e -factor=8
	$(ogen) radialTravelingWaveGrid -interp=e -factor=1
	$(ogen) radialTravelingWaveGrid -interp=e -factor=2
	$(ogen) radialTravelingWaveGrid -interp=e -factor=4
	$(ogen) radialTravelingWaveGrid -interp=e -factor=8
	$(ogen) radialTravelingWaveGrid -interp=e -factor=16


rotated:
	$(ogen) -noplot elasticPistonGridRotated -factor=1 -per=1 -thetad=30
	$(ogen) -noplot elasticPistonGridRotated -factor=2 -per=1 -thetad=30
	$(ogen) -noplot elasticPistonGridRotated -factor=4 -per=1 -thetad=30
	$(ogen) -noplot elasticPistonGridRotated -factor=8 -per=1 -thetad=30
	$(ogen) -noplot elasticPistonGridRotated -factor=16 -per=1 -thetad=30


ins:
	cd ../../; make
