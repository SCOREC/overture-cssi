
cgmx = $(CGBUILDPREFIX)/mx/bin/cgmx

plotStuff = $(Overture)/bin/plotStuff

# GDM scattering from multiple bodies
runmb:
	$(cgmx) -noplot scattering -g=objectsGride4.order4.hdf  -bg=backGround -rbc=abcPML -boundaryForcing=1 -plotTotalField=0 -kx=2 -tf=5 -tp=.05 -dm=gdm -a0=.9 -a1=0. -b0=1 -b1=0  -modeGDM=1 -useSosupDissipation=1 -show=objectsG4M0a.show -go=go > objectsG4M0a.out
	$(cgmx) -noplot scattering -g=objectsGride4.order4.hdf  -bg=backGround -rbc=abcPML -boundaryForcing=1 -plotTotalField=0 -kx=2 -tf=5 -tp=.05 -dm=gdm -a0=.9 -a1=0. -b0=1 -b1=0  -useSosupDissipation=1 -show=objectsG4a.show -go=go > objectsG4a.out
	$(cgmx) -noplot scattering -g=objectsGride8.order4.hdf  -bg=backGround -rbc=abcPML -boundaryForcing=1 -plotTotalField=0 -kx=2 -tf=5 -tp=.05 -dm=gdm -a0=.9 -a1=0. -b0=1 -b1=0  -useSosupDissipation=1 -show=objectsG8a.show -go=go > objectsG8a.out
	$(cgmx) -noplot scattering -g=objectsGride8.order4.hdf  -bg=backGround -rbc=abcPML -boundaryForcing=1 -plotTotalField=0 -kx=2 -tf=5 -tp=.05 -dm=gdm -a0=.9 -a1=0. -b0=1 -b1=0  -modeGDM=1 -useSosupDissipation=1 -show=objectsG8M0a.show -go=go > objectsG8M0a.out


# plot results for scattering from multiple bodies
plotmb:
	$(plotStuff) -nopause plotObjects.cmd -show=objectsG4a.show -name=objectsG4at1p0 -solution=21
	$(plotStuff) -nopause plotObjects.cmd -show=objectsG8a.show -name=objectsG8at1p0 -solution=21
	$(plotStuff) -nopause plotObjects.cmd -show=objectsG4a.show -name=objectsG4at5p0 -solution=101
	$(plotStuff) -nopause plotObjects.cmd -show=objectsG8a.show -name=objectsG8at5p0 -solution=101
	$(plotStuff) -nopause plotObjects.cmd -show=objectsG4M0a.show -name=objectsG4M0at1p0 -solution=21
	$(plotStuff) -nopause plotObjects.cmd -show=objectsG8M0a.show -name=objectsG8M0at1p0 -solution=21
	$(plotStuff) -nopause plotObjects.cmd -show=objectsG4M0a.show -name=objectsG4M0at5p0 -solution=101
	$(plotStuff) -nopause plotObjects.cmd -show=objectsG8M0a.show -name=objectsG8M0at5p0 -solution=101
