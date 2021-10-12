# CgMp examples -- elastic disk in a fluid channel FSI-INS-AD

ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp

# ------- Elastic disk in a fluid channel ------
#  Self-convergence runs (from fibr2019 JCP 2019)

diskInAChannelG2: diskInAChannelGride2.order2.hdf
	$(cgmp) -noplot solidDiskInAChannel -method=ins -smVariation=g  -g=diskInAChannelGride2.order2  -addedMass=1 -scf=10.0 -zfMuByH=1 -p0=1 -tp=.1 -tf=1.5 -tsINS=im -useImplicitAmpBCs=1 -nc=1 -nu=.05 -pOffset=.0  -diss=.5 -debug=1 -rampedInflow=1 -flushFrequency=4 -multiDomainAlgorithm=1 -cfl=.9 -setGhostByExtrapolation=1 -show=diskInAChannelScf10G2.show -dtMax=5e-3 -go=go >! diskInAChannelScf10G2.out

diskInAChannelG4: diskInAChannelGride4.order2.hdf
	$(cgmp) -noplot solidDiskInAChannel -method=ins -smVariation=g  -g=diskInAChannelGride4.order2  -addedMass=1 -scf=10.0 -zfMuByH=1 -p0=1 -tp=.1 -tf=1.5 -tsINS=im -useImplicitAmpBCs=1 -nc=1 -nu=.05 -pOffset=.0  -diss=.5 -debug=1 -rampedInflow=1 -flushFrequency=4 -multiDomainAlgorithm=1 -cfl=.9 -setGhostByExtrapolation=1 -show=diskInAChannelScf10G4.show -dtMax=2.5e-3 -go=go >! diskInAChannelScf10G4.out

diskInAChannelG8: diskInAChannelGride8.order2.hdf
	$(cgmp) -noplot -noplot solidDiskInAChannel -method=ins -smVariation=g  -g=diskInAChannelGride8.order2  -addedMass=1 -scf=10.0 -zfMuByH=1 -p0=1 -tp=.1 -tf=1.5 -tsINS=im -useImplicitAmpBCs=1 -nc=1 -nu=.05 -pOffset=.0  -diss=.5 -debug=1 -rampedInflow=1 -flushFrequency=4 -multiDomainAlgorithm=1 -cfl=.9 -setGhostByExtrapolation=1 -show=diskInAChannelScf10G8.show -dtMax=1.25e-3 -go=go >! diskInAChannelScf10G8.out

# grids 
diskInAChannelGride2.order2.hdf:
	$(ogen) -noplot diskInAChannelGrid -interp=e -factor=2

diskInAChannelGride4.order2.hdf:
	$(ogen) -noplot diskInAChannelGrid -interp=e -factor=4

diskInAChannelGride8.order2.hdf:
	$(ogen) -noplot diskInAChannelGrid -interp=e -factor=8
