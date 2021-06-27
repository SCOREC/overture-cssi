#
#   Makefile for the elastic piston
#


ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp
conv = $(CG)/common/src/conv.p

# radial elastic piston example (light solid):
radialElasticPiston: radialElasticPistonGride4.order2.hdf
	$(cgmp) -noplot radialElasticPiston -g=radialElasticPistonGride4.order2.hdf -addedMass=1 -useImplicitAmpBCs=1 -scf=1.0 -nu=.05 -tp=.05 -tf=1.0 -smVariation=g -debug=1 -tsINS=im -stressRelaxation=4 -diss=.5  -setGhostByExtrapolation=1 -show=repG4Scf1.show -go=go >! repG4Scf1.out

# annular SHEAR example (medium solid):
annularShearMedium: radialElasticPistonGride4.order2.hdf
	$(cgmp) -noplot radialElasticPiston -g=radialElasticPistonGride4.order2.hdf -known=shear -addedMass=1 -numberOfCorrections=1 -useImplicitAmpBCs=1 -scf=1 -nu=0.1 -tp=.05 -tf=.5 -smVariation=g -debug=1 -tsINS=im -stressRelaxation=4 -diss=.5  -setGhostByExtrapolation=1 -show=radialShearG4Scf1.show -go=go >! radialShearG4Scf1.out

# grid
radialElasticPistonGride4.order2.hdf:
	$(ogen) -noplot radialElasticPistonGrid -interp=e -factor=4



grids: elasticPiston radial rotated

elasticPiston:
	$(ogen) -noplot elasticPistonGrid -factor=1
	$(ogen) -noplot elasticPistonGrid -factor=2
	$(ogen) -noplot elasticPistonGrid -factor=4
	$(ogen) -noplot elasticPistonGrid -factor=8
	$(ogen) -noplot elasticPistonGrid -factor=16
	$(ogen) -noplot elasticPistonGrid -factor=1 -per=1
	$(ogen) -noplot elasticPistonGrid -factor=2 -per=1
	$(ogen) -noplot elasticPistonGrid -factor=4 -per=1
	$(ogen) -noplot elasticPistonGrid -factor=8 -per=1
	$(ogen) -noplot elasticPistonGrid -factor=16 -per=1
	$(ogen) -noplot elasticPistonGrid -factor=32 -per=1
	$(ogen) -noplot elasticPistonGrid -factor=64 -per=1


# Narrow grid in x-direction
elasticPistonN:
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -factor=1 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -factor=2 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -factor=4 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -factor=8 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -per=1 -factor=1 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -per=1 -factor=2 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -per=1 -factor=4 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -per=1 -factor=8 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -per=1 -factor=16

# coarse grid in x: 
elasticPistonX:
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -factor=1 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -factor=2 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -factor=4 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -factor=8 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -per=1 -factor=1 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -per=1 -factor=2 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -per=1 -factor=4 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -per=1 -factor=8 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -per=1 -factor=16


radial:
	$(ogen) -noplot radialElasticPistonGrid -interp=e -factor=1
	$(ogen) -noplot radialElasticPistonGrid -interp=e -factor=2
	$(ogen) -noplot radialElasticPistonGrid -interp=e -factor=4
	$(ogen) -noplot radialElasticPistonGrid -interp=e -factor=8
	$(ogen) -noplot radialElasticPistonGrid -interp=e -factor=16
	$(ogen) -noplot radialTravelingWaveGrid -interp=e -factor=1
	$(ogen) -noplot radialTravelingWaveGrid -interp=e -factor=2
	$(ogen) -noplot radialTravelingWaveGrid -interp=e -factor=4
	$(ogen) -noplot radialTravelingWaveGrid -interp=e -factor=8
	$(ogen) -noplot radialTravelingWaveGrid -interp=e -factor=16

rotated:
	$(ogen) -noplot elasticPistonGridRotated -factor=1 -per=1 -thetad=0
	$(ogen) -noplot elasticPistonGridRotated -factor=2 -per=1 -thetad=0
	$(ogen) -noplot elasticPistonGridRotated -factor=4 -per=1 -thetad=0
	$(ogen) -noplot elasticPistonGridRotated -factor=8 -per=1 -thetad=0
	$(ogen) -noplot elasticPistonGridRotated -factor=1 -per=1 -thetad=30
	$(ogen) -noplot elasticPistonGridRotated -factor=2 -per=1 -thetad=30
	$(ogen) -noplot elasticPistonGridRotated -factor=4 -per=1 -thetad=30
	$(ogen) -noplot elasticPistonGridRotated -factor=8 -per=1 -thetad=30
	# $(ogen) -noplot elasticPistonGridRotated -factor=1 -per=1 -thetad=0
	# $(ogen) -noplot elasticPistonGridRotated -factor=2 -per=1 -thetad=30
	# $(ogen) -noplot elasticPistonGridRotated -factor=4 -per=1 -thetad=30
	# $(ogen) -noplot elasticPistonGridRotated -factor=8 -per=1 -thetad=30
	# $(ogen) -noplot elasticPistonGridRotated -factor=1 -per=1 -thetad=45
	# $(ogen) -noplot elasticPistonGridRotated -factor=2 -per=1 -thetad=45
	# $(ogen) -noplot elasticPistonGridRotated -factor=4 -per=1 -thetad=45
	# $(ogen) -noplot elasticPistonGridRotated -factor=8 -per=1 -thetad=45


# 
# NOTE: To run the following convergence tests:
#      (1) Make a symbolic link to cgmp/conv/insElasticPiston.conv
#      (2) Make a sub-directory results to hold results
radialWaveConv:
	# radial traveling wave
	$(conv) insElasticPiston.conv -known=radialTravelingWave -cg=rtw -tp=0.3 -tf=0.3 -showAllComponents=0 -tsINS=im -smVariation=g -caseid=7 -scf=0.1 -nu=0.1  -zfMuByH=1 -amp=1e-7 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -ng=4 -zfOption=m -zfMono=1
	$(conv) insElasticPiston.conv -known=radialTravelingWave -cg=rtw -tp=0.3 -tf=0.3 -showAllComponents=0 -tsINS=im -smVariation=g -caseid=0 -scf=1 -nu=0.1  -zfMuByH=1 -amp=1e-7 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -ng=4 -zfOption=m -zfMono=1
	$(conv) insElasticPiston.conv -known=radialTravelingWave -cg=rtw -tp=0.3 -tf=0.3 -showAllComponents=0 -tsINS=im -smVariation=g -caseid=8 -scf=10 -nu=0.1  -zfMuByH=1 -amp=1e-7 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -ng=4 -zfOption=m -zfMono=1

radialShearConv:
	# radial shear
	$(conv) insElasticPiston.conv -known=shear -nu=0.1 -tp=0.6 -tf=0.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfOption=m -zfMono=1 -zfMuByH=1 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -ng=4
	$(conv) insElasticPiston.conv -known=shear -nu=0.1 -tp=0.6 -tf=0.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfOption=m -zfMono=1 -zfMuByH=1 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -ng=4
	$(conv) insElasticPiston.conv -known=shear -nu=0.1 -tp=0.6 -tf=0.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfOption=m -zfMono=1 -zfMuByH=1 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -ng=4

radialPistonConv:
	# radial piston
	$(conv) insElasticPiston.conv -nu=0.05 -tp=0.6 -tf=0.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfOption=m -zfMono=1 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -ng=4
	$(conv) insElasticPiston.conv -nu=0.05 -tp=0.6 -tf=0.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfOption=m -zfMono=1 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -ng=4
	$(conv) insElasticPiston.conv -nu=0.05 -tp=0.6 -tf=0.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfOption=m -zfMono=1 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -ng=4

blockShearConv:
	# block shear
	$(conv) insElasticPiston.conv -known=shear -nu=0.1 -tp=0.6 -tf=0.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=1 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=4 -setGhostByExtrapolation=1 -ng=4
	$(conv) insElasticPiston.conv -known=shear -nu=0.1 -tp=0.6 -tf=0.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfMuByH=1 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=4 -setGhostByExtrapolation=1 -ng=4
	$(conv) insElasticPiston.conv -known=shear -nu=0.1 -tp=0.6 -tf=0.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=1 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=4 -setGhostByExtrapolation=1 -ng=4

blockPistonConv:
	# block piston
	$(conv) insElasticPiston.conv -nu=0.05 -tp=0.6 -tf=0.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=1 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=0 -setGhostByExtrapolation=1 -ng=4
	$(conv) insElasticPiston.conv -nu=0.05 -tp=0.6 -tf=0.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfMuByH=1 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=0 -setGhostByExtrapolation=1 -ng=4
	$(conv) insElasticPiston.conv -nu=0.05 -tp=0.6 -tf=0.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=1 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=0 -setGhostByExtrapolation=1 -ng=4

blockPistonConvTP:
	#$(conv) insElasticPiston.conv -nu=0.05 -tp=0.6 -tf=0.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=0.001 -addedMass=0 -tsSM=improvedEuler -useTP=1 -cfl=0.9 -diss=0.5 -stressRelaxation=0 -setGhostByExtrapolation=1 -ng=4
	$(conv) insElasticPiston.conv -nu=0.05 -tp=0.6 -tf=0.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -addedMass=0 -tsSM=improvedEuler -useTP=1 -cfl=0.9 -diss=0.5 -stressRelaxation=0 -setGhostByExtrapolation=1 -ng=4
	#$(conv) insElasticPiston.conv -nu=0.05 -tp=0.6 -tf=0.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1000 -addedMass=0 -tsSM=improvedEuler -useTP=1 -cfl=0.9 -diss=0.5 -stressRelaxation=0 -setGhostByExtrapolation=1 -ng=4

blockPistonConvAMP:
	$(conv) insElasticPiston.conv -nu=0.05 -tp=0.6 -tf=0.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfMuByH=0 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=0 -setGhostByExtrapolation=1 -ng=4
	$(conv) insElasticPiston.conv -nu=0.05 -tp=0.6 -tf=0.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1000 -zfMuByH=0 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=0 -setGhostByExtrapolation=1 -ng=4
	$(conv) insElasticPiston.conv -nu=0.05 -tp=0.6 -tf=0.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=0.001 -zfMuByH=0 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=0.9 -diss=0.5 -stressRelaxation=0 -setGhostByExtrapolation=1 -ng=4


mp:
	cd ../../; \
	make

