common = ../../common

# The convergence rate script is called conv.p
all:  
	ln -sf $(common)/src/conv.p .
	ln -sf ../cmd/cnsDomain.h .
	ln -sf ../cmd/insDomain.h .
	ln -sf ../cmd/adDomain.h .


cgmp:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) cgmp

DropBox = $(HOME)/Dropbox

ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp


# double annulus grids
doubleAnnulusGrids:
	$(ogen) -noplot doubleAnnulusGrid -order=2 -interp=e -factor=1
	$(ogen) -noplot doubleAnnulusGrid -order=2 -interp=e -factor=2
	$(ogen) -noplot doubleAnnulusGrid -order=2 -interp=e -factor=4
	$(ogen) -noplot doubleAnnulusGrid -order=2 -interp=e -factor=8
	$(ogen) -noplot doubleAnnulusGrid -order=2 -interp=e -factor=16
	$(ogen) -noplot doubleAnnulusGrid -order=2 -interp=e -factor=32
	$(ogen) -noplot doubleAnnulusGrid -order=2 -interp=e -factor=64
	$(ogen) -noplot doubleAnnulusGrid -order=4 -interp=e -factor=1
	$(ogen) -noplot doubleAnnulusGrid -order=4 -interp=e -factor=2
	$(ogen) -noplot doubleAnnulusGrid -order=4 -interp=e -factor=4
	$(ogen) -noplot doubleAnnulusGrid -order=4 -interp=e -factor=8
	$(ogen) -noplot doubleAnnulusGrid -order=4 -interp=e -factor=16
	$(ogen) -noplot doubleAnnulusGrid -order=4 -interp=e -factor=32
	$(ogen) -noplot doubleAnnulusGrid -order=4 -interp=e -factor=64	

# convergence for concentric annulii, CHAMP and DN 
doubleAnnulus:
	conv.p twoDomainsCHT.conv -problem=doubleAnnulusChamp -cg=doubleAnnulusGrid -known=concentricCylinders -useChamp=1 -tp=.3 -tf=.3 -nMode=1 -mRoot=0 -ng=6
	conv.p twoDomainsCHT.conv -problem=doubleAnnulus -cg=doubleAnnulusGrid -known=concentricCylinders -useChamp=0 -tp=.3 -tf=.3 -nMode=1 -mRoot=0 -ng=6



# Radial traveling wave
# caseid=0 scf=1 and nu=.1
# caseid=2 : scf=.001 nu=.1
# caseid=4 : scf=10 nu=.1
#
# RESULTS FOR PAPER
radialTravelingWave:
	conv.p insElasticPiston.conv -known=radialTravelingWave -cg=rtw -tp=.05 -tf=.2 -showAllComponents=0 -tsINS=im -smVariation=g -caseid=0 -scf=1 -nu=.1  -zfMuByH=1 -amp=1e-5 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=5.e-2 -ng=4
	cp insElasticPiston.rtw.radialTravelingWaveAmpCaseid0SCF1zf1+0Imtp2.ConvTable.tex $(DropBox)/CG/fibr/tables
	cp insElasticPiston_rtw_radialTravelingWaveAmpCaseid0SCF1zf1p0Imtp2.m $(DropBox)/CG/fibr/matlab
	conv.p insElasticPiston.conv -known=radialTravelingWave -cg=rtw -tp=.05 -tf=.3 -showAllComponents=0 -tsINS=im -smVariation=g -caseid=2 -scf=1e-3 -nu=.1  -zfMuByH=1 -amp=1e-5 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=5.e-2 -ng=4
	cp insElasticPiston.rtw.radialTravelingWaveAmpCaseid2SCF1e-3zf1+0Imtp3.ConvTable.tex $(DropBox)/CG/fibr/tables
	cp insElasticPiston_rtw_radialTravelingWaveAmpCaseid2SCF1em3zf1p0Imtp3.m $(DropBox)/CG/fibr/matlab
	conv.p insElasticPiston.conv -known=radialTravelingWave -cg=rtw -tp=.05 -tf=.3 -showAllComponents=0 -tsINS=im -smVariation=g -caseid=4 -scf=10 -nu=.1  -zfMuByH=1 -amp=1e-5 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=5.e-2 -ng=4
	cp insElasticPiston.rtw.radialTravelingWaveAmpCaseid4SCF10zf1+0Imtp3.ConvTable.tex $(DropBox)/CG/fibr/tables 
	cp insElasticPiston_rtw_radialTravelingWaveAmpCaseid4SCF10zf1p0Imtp3.m $(DropBox)/CG/fibr/matlab

# case 2 run longer
radialTravelingWave2:
	conv.p insElasticPiston.conv -known=radialTravelingWave -cg=rtw -tp=.05 -tf=.5 -showAllComponents=0 -tsINS=im -smVariation=g -caseid=2 -scf=1e-3 -nu=.1  -zfMuByH=1 -amp=1e-5 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=5.e-2 -ng=4
	cp insElasticPiston.rtw.radialTravelingWaveAmpCaseid2SCF1e-3zf1+0Imtp5.ConvTable.tex $(DropBox)/CG/fibr/tables
	cp insElasticPiston_rtw_radialTravelingWaveAmpCaseid2SCF1em3zf1p0Imtp5.m $(DropBox)/CG/fibr/matlab

# case 2 run even longer
radialTravelingWave2a:
	conv.p insElasticPiston.conv -known=radialTravelingWave -cg=rtw -tp=.05 -tf=.9 -showAllComponents=0 -tsINS=im -smVariation=g -caseid=2 -scf=1e-3 -nu=.1  -zfMuByH=1 -amp=1e-5 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=5.e-2 -ng=4
	cp insElasticPiston.rtw.radialTravelingWaveAmpCaseid2SCF1e-3zf1+0Imtp9.ConvTable.tex $(DropBox)/CG/fibr/tables
	cp insElasticPiston_rtw_radialTravelingWaveAmpCaseid2SCF1em3zf1p0Imtp9.m $(DropBox)/CG/fibr/matlab



# Elastic piston 
elasticPiston: blockElasticPiston shearElasticPiston


#  Block elastic piston -- monolithic zf ---
blockElasticPistonm:
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfMuByH=1 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1zfmImtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=1 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1e3zfmImtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=1 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1e-3zfmImtp6.ConvTable.tex $(DropBox)/CG/fibr/tables

shearElasticPistonm:
	conv.p insElasticPiston.conv -known=shear -nu=10 -tp=.1 -tf=.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=1 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgp.shearAmpSCF1e3zfmImtp6.ConvTable.tex  $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -known=shear -nu=.1 -tp=.1 -tf=.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=1 -zfRhoHByDt=0 -zfOption=m -zfMono=1 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgp.shearAmpSCF1e-3zfmImtp6.ConvTable.tex  $(DropBox)/CG/fibr/tables

radialElasticPistonm:
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfOption=m -zfMono=1 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1zfmImtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfOption=m -zfMono=1 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1e3zfmImtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfOption=m -zfMono=1 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1e-3zfmImtp6.ConvTable.tex $(DropBox)/CG/fibr/tables

radialShearElasticPistonm:
	conv.p insElasticPiston.conv -known=shear -nu=10. -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=10 -zfOption=m -zfMono=1 -zfMuByH=1 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.shearAmpSCF10zfmImtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -known=shear -nu=10. -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfOption=m -zfMono=1 -zfMuByH=1 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.shearAmpSCF1e3zfmImtp6.ConvTable.tex $(DropBox)/CG/fibr/tables




#  Block elastic piston
blockElasticPiston100:
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfMuByH=100 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1zf100+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=100 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1e3zf100+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=100 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1e-3zf100+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables


blockElasticPiston5:
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfMuByH=5 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1zf5+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=5 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1e3zf5+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=5 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1e-3zf5+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables


blockElasticPiston2:
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfMuByH=2 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1zf2+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=2 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1e3zf2+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=2 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1e-3zf2+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables


blockElasticPiston1:
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1zf1+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1e3zf1+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=epgpN -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=0 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgpN.AmpSCF1e-3zf1+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables


# block shear elastic piston
shearElasticPiston100:
	conv.p insElasticPiston.conv -known=shear -nu=10 -tp=.1 -tf=.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=100 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgp.shearAmpSCF1e3zf100+0Imtp6.ConvTable.tex  $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -known=shear -nu=.1 -tp=.1 -tf=.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=100 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgp.shearAmpSCF1e-3zf100+0Imtp6.ConvTable.tex  $(DropBox)/CG/fibr/tables


shearElasticPiston5:
	conv.p insElasticPiston.conv -known=shear -nu=10 -tp=.1 -tf=.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=5 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgp.shearAmpSCF1e3zf5+0Imtp6.ConvTable.tex  $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -known=shear -nu=.1 -tp=.1 -tf=.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=5 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgp.shearAmpSCF1e-3zf5+0Imtp6.ConvTable.tex  $(DropBox)/CG/fibr/tables


shearElasticPiston2:
	conv.p insElasticPiston.conv -known=shear -nu=10 -tp=.1 -tf=.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=2 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgp.shearAmpSCF1e3zf2+0Imtp6.ConvTable.tex  $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -known=shear -nu=.1 -tp=.1 -tf=.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=2 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgp.shearAmpSCF1e-3zf2+0Imtp6.ConvTable.tex  $(DropBox)/CG/fibr/tables

shearElasticPiston1:
	conv.p insElasticPiston.conv -known=shear -nu=10 -tp=.1 -tf=.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgp.shearAmpSCF1e3zf1+0Imtp6.ConvTable.tex  $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -known=shear -nu=.1 -tp=.1 -tf=.6 -cg=epgp -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -dtMax=2.e-2 -setGhostByExtrapolation=1 -ng=4
	cp insElasticPiston.epgp.shearAmpSCF1e-3zf1+0Imtp6.ConvTable.tex  $(DropBox)/CG/fibr/tables


# radial elastic piston
radialElasticPiston100:
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1  -zfMuByH=100 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1zf100+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3  -zfMuByH=100 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1e3zf100+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3  -zfMuByH=100 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1e-3zf100+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables

radialElasticPiston5:
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1  -zfMuByH=5 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1zf5+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3  -zfMuByH=5 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1e3zf5+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3  -zfMuByH=5 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1e-3zf5+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables


radialElasticPiston2:
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfMuByH=2 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1zf2+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=2 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1e3zf2+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=2 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1e-3zf2+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables

# RESULTS FOR PAPER
radialElasticPiston1:
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1zf1+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	cp insElasticPiston_rep_AmpSCF1zf1p0Imtp6.m $(DropBox)/CG/fibr/matlab
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1e3zf1+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	cp insElasticPiston_rep_AmpSCF1e3zf1p0Imtp6.m $(DropBox)/CG/fibr/matlab
	conv.p insElasticPiston.conv -nu=.05 -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e-3 -zfMuByH=1 -zfRhoHByDt=0 -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.AmpSCF1e-3zf1+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	cp insElasticPiston_rep_AmpSCF1em3zf1p0Imtp6.m $(DropBox)/CG/fibr/matlab

# radial shear solution
radialShearElasticPiston100:
	conv.p insElasticPiston.conv -known=shear -nu=10. -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=10 -zfMuByH=100 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.shearAmpSCF10zf100+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -known=shear -nu=10. -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=100 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.shearAmpSCF1e3zf100+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables

radialShearElasticPiston5:
	conv.p insElasticPiston.conv -known=shear -nu=10. -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=10 -zfMuByH=5 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.shearAmpSCF10zf5+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -known=shear -nu=10. -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=5 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.shearAmpSCF1e3zf5+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables


radialShearElasticPiston2:
	conv.p insElasticPiston.conv -known=shear -nu=10. -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=10 -zfMuByH=2 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.shearAmpSCF10zf2+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	conv.p insElasticPiston.conv -known=shear -nu=10. -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=2 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.shearAmpSCF1e3zf2+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables

# RESULTS FOR PAPER
radialShearElasticPiston1:
	conv.p insElasticPiston.conv -known=shear -nu=10. -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=10 -zfMuByH=1 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.shearAmpSCF10zf1+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	cp insElasticPiston_rep_shearAmpSCF10zf1p0Imtp6.m $(DropBox)/CG/fibr/matlab
	conv.p insElasticPiston.conv -known=shear -nu=10. -tp=.1 -tf=.6 -cg=rep -showAllComponents=0 -tsINS=im -smVariation=g -scf=1e3 -zfMuByH=1 -zfRhoHByDt=0  -addedMass=1 -useImplicitAmpBCs=1 -numberOfCorrections=1 -cfl=.9 -diss=0.5 -stressRelaxation=4 -startCurve=spline  -setGhostByExtrapolation=1 -dtMax=2.e-2 -ng=4
	cp insElasticPiston.rep.shearAmpSCF1e3zf1+0Imtp6.ConvTable.tex $(DropBox)/CG/fibr/tables
	cp insElasticPiston_rep_shearAmpSCF1e3zf1p0Imtp6.m $(DropBox)/CG/fibr/matlab
