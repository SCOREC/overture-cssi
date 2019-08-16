#
#  Scattering of a dipsersive plane wave from many disks
#

cgmx = $(CGBUILDPREFIX)/mx/bin/cgmx
ogen = $(Overture)/bin/ogen

# Make grids: 
grids:
	$(ogen) -noplot multiDiskGrid -interp=e -order=2 -factor=2
	$(ogen) -noplot multiDiskGrid -interp=e -order=2 -factor=4
	$(ogen) -noplot multiDiskGrid -interp=e -order=2 -factor=8
	$(ogen) -noplot multiDiskGrid -interp=e -order=2 -factor=16
	$(ogen) -noplot multiDiskGrid -interp=e -order=4 -numGhost=3 -factor=2 
	$(ogen) -noplot multiDiskGrid -interp=e -order=4 -numGhost=3 -factor=4 
	$(ogen) -noplot multiDiskGrid -interp=e -order=4 -numGhost=3 -factor=8
	$(ogen) -noplot multiDiskGrid -interp=e -order=4 -numGhost=3 -factor=16


# runs for self convergence
# SECOND-ORDER
runsOrder2:
	$(cgmx) -noplot multiDomain -g=multiDiskGride2.order2 -dm=gdm -dmFile=gdmMaterial1.txt gdmMaterial2.txt gdmMaterial3.txt gdmMaterial4.txt -backGround=backGround -rbc=rbcNonLocal -kx=.5 -eps0=1 -useSosupDissipation=1 -tf=5. -tp=.5 -ic=gp -x0=-3. -beta=120. -show=MD2.show -go=go > MD2.out
	$(cgmx) -noplot multiDomain -g=multiDiskGride4.order2 -dm=gdm -dmFile=gdmMaterial1.txt gdmMaterial2.txt gdmMaterial3.txt gdmMaterial4.txt -backGround=backGround -rbc=rbcNonLocal -kx=.5 -eps0=1 -useSosupDissipation=1 -tf=5. -tp=.5 -ic=gp -x0=-3. -beta=120. -show=MD4.show -go=go > MD4.out
	$(cgmx) -noplot multiDomain -g=multiDiskGride8.order2 -dm=gdm -dmFile=gdmMaterial1.txt gdmMaterial2.txt gdmMaterial3.txt gdmMaterial4.txt -backGround=backGround -rbc=rbcNonLocal -kx=.5 -eps0=1 -useSosupDissipation=1 -tf=5. -tp=.5 -ic=gp -x0=-3. -beta=120. -show=MD8.show -go=go > MD8.out


# FOURTH-ORDER
runsOrder4:
	$(cgmx) -noplot multiDomain -g=multiDiskGride2.order4.ng3 -dm=gdm -dmFile=gdmMaterial1.txt gdmMaterial2.txt gdmMaterial3.txt gdmMaterial4.txt -backGround=backGround -rbc=rbcNonLocal -kx=.5 -eps0=1 -useSosupDissipation=1 -tf=5. -tp=.5 -ic=gp -x0=-3. -beta=120. -show=MD2a.show -go=go > MD2a.out
	$(cgmx) -noplot multiDomain -g=multiDiskGride4.order4.ng3 -dm=gdm -dmFile=gdmMaterial1.txt gdmMaterial2.txt gdmMaterial3.txt gdmMaterial4.txt -backGround=backGround -rbc=rbcNonLocal -kx=.5 -eps0=1 -useSosupDissipation=1 -tf=5. -tp=.5 -ic=gp -x0=-3. -beta=120. -show=MD4a.show -go=go > MD4a.out
	$(cgmx) -noplot multiDomain -g=multiDiskGride8.order4.ng3 -dm=gdm -dmFile=gdmMaterial1.txt gdmMaterial2.txt gdmMaterial3.txt gdmMaterial4.txt -backGround=backGround -rbc=rbcNonLocal -kx=.5 -eps0=1 -useSosupDissipation=1 -tf=5. -tp=.5 -ic=gp -x0=-3. -beta=120. -show=MD8a.show -go=go > MD8a.out


#	$(cgmx) noplot multiDomain -g=multiDiskGride16.order4.ng3 -dm=gdm -dmFile=gdmMaterial1.txt gdmMaterial2.txt gdmMaterial3.txt gdmMaterial4.txt -backGround=backGround -rbc=rbcNonLocal -kx=.5 -eps0=1 -useSosupDissipation=1 -tf=1. -tp=.5 -ic=gp -x0=-2.5 -beta=120. -show=MD16a.show -go=go >! MD16a.out
