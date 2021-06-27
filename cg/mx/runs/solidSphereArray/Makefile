#
#   Makefile for solidSphereArray
#


ogen = $(Overture)/bin/ogen

cgmx = $(CGBUILDPREFIX)/mx/bin/cgmx

# 64 solid spheres
64spheres:
	 $(ogen) -noplot solidSphereArrayGrid -prefix=64SolidSphereGrid -xa=-5.5 -xb=5.5 -ya=-5.5 -yb=5.5 -za=-5.5 -zb=5.5 -nsx=4 -nsy=4 -nsz=4 -order=2 -factor=2
	 $(ogen) -noplot solidSphereArrayGrid -prefix=64SolidSphereGrid -xa=-5.5 -xb=5.5 -ya=-5.5 -yb=5.5 -za=-5.5 -zb=5.5 -nsx=4 -nsy=4 -nsz=4 -order=4 -numGhost=3 -factor=2


# 27 solid spheres
27spheres:
	$(ogen) -noplot solidSphereArrayGrid -prefix=27SolidSphereGrid -xa=-4.5 -xb=4.5 -ya=-5 -yb=5 -za=-5 -zb=5 -nsx=3 -nsy=3 -nsz=3 -interp=e -order=2 -factor=2
	$(ogen) -noplot solidSphereArrayGrid -prefix=27SolidSphereGrid -xa=-4.5 -xb=4.5 -ya=-5 -yb=5 -za=-5 -zb=5 -nsx=3 -nsy=3 -nsz=3 -interp=e -order=4 -numGhost=3 -factor=2


# 8 solid spheres
8spheres:
	$(ogen) -noplot solidSphereArrayGrid -prefix=eightSolidSphereGrid -xa=-3 -xb=3 -ya=-3 -yb=3 -za=-3 -zb=3 - -nsx=2 -nsy=2 -nsz=2 -interp=e -order=2 -factor=2
	$(ogen) -noplot solidSphereArrayGrid -prefix=eightSolidSphereGrid -xa=-3 -xb=3 -ya=-3 -yb=3 -za=-3 -zb=3 - -nsx=2 -nsy=2 -nsz=2 -interp=e -order=4 -numGhost=3 -factor=2



# =================== Ellipsoid grids ==========================

#  ---- One ellipsoid ----
ellipsoidGrid1Order2G2:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=1SolidEllipsoidsGrid -nsx=1 -nsy=1 -nsz=1 -configFile=solidEllipsoidArrayConfig3d.h -order=2 -per=npp -interp=i -factor=2  -rFactor=.65 -nrExtra=3 -orthographicPatchParameter=.7
ellipsoidGrid1Order2G4:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=1SolidEllipsoidsGrid -nsx=1 -nsy=1 -nsz=1 -configFile=solidEllipsoidArrayConfig3d.h -order=2 -per=npp -interp=e -rFactor=.65 -nrExtra=3 -orthographicPatchParameter=.7 -factor=4
ellipsoidGrid1Order2G8:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=1SolidEllipsoidsGrid -nsx=1 -nsy=1 -nsz=1 -configFile=solidEllipsoidArrayConfig3d.h -order=2 -per=npp -interp=e -rFactor=.65 -nrExtra=3 -orthographicPatchParameter=.7 -factor=8
ellipsoidGrid1Order2G16:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=1SolidEllipsoidsGrid -nsx=1 -nsy=1 -nsz=1 -configFile=solidEllipsoidArrayConfig3d.h -order=2 -per=npp -interp=e -rFactor=.65 -nrExtra=3 -orthographicPatchParameter=.7 -factor=16

ellipsoidGrid1Order4G4:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=1SolidEllipsoidsGrid -nsx=1 -nsy=1 -nsz=1 -configFile=solidEllipsoidArrayConfig3d.h -order=4 -per=npp -rFactor=.7 -nrExtra=4 -orthographicPatchParameter=.7 -interp=e -numGhost=3 -factor=4 
ellipsoidGrid1Order4G8:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=1SolidEllipsoidsGrid -nsx=1 -nsy=1 -nsz=1 -configFile=solidEllipsoidArrayConfig3d.h -order=4 -per=npp -rFactor=.7 -nrExtra=5 -orthographicPatchParameter=.7 -interp=e -numGhost=3 -factor=8
ellipsoidGrid1Order4G16:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=1SolidEllipsoidsGrid -nsx=1 -nsy=1 -nsz=1 -configFile=solidEllipsoidArrayConfig3d.h -order=4 -per=npp -rFactor=.7 -nrExtra=6 -orthographicPatchParameter=.7 -interp=e -numGhost=3 -factor=16

# --- four ellipsoids ----
ellipsoidGrid4:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=fourSolidEllipsoidsGrid -configFile=solidEllipsoidArrayConfig4.h -order=2 -per=npp -rFactor=.65 -nrExtra=3 -orthographicPatchParameter=.7 -factor=4

# --- eight ellipsoids ----
ellipsoidGrid8:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=eightSolidEllipsoidsGrid -configFile=solidEllipsoidArrayConfig8.h -order=4 -per=npp -rFactor=.65 -nrExtra=4 -orthographicPatchParameter=.6 -factor=4 -matchGridsAtInterfaces=1 -interp=i -numGhost=3
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=eightSolidEllipsoidsGrid -configFile=solidEllipsoidArrayConfig8.h -order=2 -per=npp -rFactor=.65 -nrExtra=3 -orthographicPatchParameter=.7 -factor=4

# ---- 32 ellipsoids ----
ellipsoidGrid32:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=36SolidEllipsoidsGrid -nsx=3 -nsy=3 -nsz=4 -configFile=solidEllipsoidArrayConfig3d.h -order=2 -per=npp -rFactor=.65 -nrExtra=3 -orthographicPatchParameter=.7 -factor=4

# ================== ONE ELLIPSOID MLA ================================
# Note: coarse grid needs CFL=.5 for tf>2 (maybe nonlinear terms affect time-step ?)
runEllipsoid1Order2G2:
	$(cgmx) -noplot mxScript -g=1SolidEllipsoidsGridi2.order2 -domainBaseName=ellipsoid -numBodies=1 -dm=gdm -nm=mla -materialFile=ellipsoidMlaMaterials.h -backGround=backGround -cfl=.5 -rbc=abcEM2 -kx=.5 -eps0=1 -eps1=2 -useSosupDissipation=1 -tf=2 -tp=.1 -ic=gp -x0=-1.0 -beta=120 -flushFrequency=4 -yContour=0 -zContour=0 -extraCommandsEnd=plotSpheres.h -show=oneSolidEllipsoidMLAOrder2G2.show -go=go > oneSolidEllipsoidMLAOrder2G2.out

runEllipsoid1Order2G4:
	$(cgmx) -noplot mxScript -g=1SolidEllipsoidsGride4.order2 -domainBaseName=ellipsoid -numBodies=1 -dm=gdm -nm=mla -materialFile=ellipsoidMlaMaterials.h -backGround=backGround -rbc=abcEM2 -kx=.5 -eps0=1 -eps1=2 -useSosupDissipation=1 -tf=2 -tp=.1 -ic=gp -x0=-1.0 -beta=120 -flushFrequency=4  -yContour=0 -zContour=0 -extraCommandsEnd=plotSpheres.h -show=oneSolidEllipsoidMLAOrder2G4.show -go=go > oneSolidEllipsoidMLAOrder2G4.out 

runEllipsoid1Order2G8: 
	$(cgmx) -noplot mxScript -g=1SolidEllipsoidsGride8.order2 -domainBaseName=ellipsoid -numBodies=1 -dm=gdm -nm=mla -materialFile=ellipsoidMlaMaterials.h -backGround=backGround -rbc=abcEM2 -kx=.5 -eps0=1 -eps1=2 -useSosupDissipation=1 -tf=2 -tp=.1 -ic=gp -x0=-1.0 -beta=120 -flushFrequency=4  -yContour=0 -zContour=0 -extraCommandsEnd=plotSpheres.h -show=oneSolidEllipsoidMLAOrder2G8.show -go=go > oneSolidEllipsoidMLAOrder2G8.out	

runEllipsoid1Order2G16:
	$(cgmx) -noplot mxScript -g=1SolidEllipsoidsGride16.order2 -domainBaseName=ellipsoid -numBodies=1 -dm=gdm -nm=mla -materialFile=ellipsoidMlaMaterials.h -backGround=backGround -rbc=abcEM2 -kx=.5 -eps0=1 -eps1=2 -useSosupDissipation=1 -tf=2 -tp=.1 -ic=gp -x0=-1.0 -beta=120 -flushFrequency=4 -yContour=0 -zContour=0 -extraCommandsEnd=plotSpheres.h -show=oneSolidEllipsoidMLAOrder2G16.show -go=go > oneSolidEllipsoidMLAOrder2G16.out


runEllipsoid1Order4G4:
	$(cgmx) -noplot mxScript -g=1SolidEllipsoidsGride4.order4.ng3.hdf -domainBaseName=ellipsoid -numBodies=1 -dm=gdm -nm=mla -materialFile=ellipsoidMlaMaterials.h -backGround=backGround -rbc=rbcNonLocal -kx=.5 -eps0=1 -eps1=2 -useSosupDissipation=1 -tf=2 -tp=.1 -ic=gp -x0=-1.0 -beta=120 -flushFrequency=4 -yContour=0 -zContour=0 -extraCommandsEnd=plotSpheres.h -show=oneSolidEllipsoidMLAOrder4G4.show -go=go > oneSolidEllipsoidMLAOrder4G4.out

runEllipsoid1Order4G8:
	setsid $(cgmx) -noplot mxScript -g=1SolidEllipsoidsGride8.order4.ng3.hdf -domainBaseName=ellipsoid -numBodies=1 -dm=gdm -nm=mla -materialFile=ellipsoidMlaMaterials.h -backGround=backGround -rbc=rbcNonLocal -kx=.5 -eps0=1 -eps1=2 -useSosupDissipation=1 -tf=2 -tp=.1 -ic=gp -x0=-1.0 -beta=120 -flushFrequency=4 -yContour=0 -zContour=0 -extraCommandsEnd=plotSpheres.h -show=oneSolidEllipsoidMLAOrder4G8.show -go=go > oneSolidEllipsoidMLAOrder4G8.out

runEllipsoid1Order4G16:
	setsid $(cgmx) -noplot mxScript -g=1SolidEllipsoidsGride16.order4.ng3.hdf -domainBaseName=ellipsoid -numBodies=1 -dm=gdm -nm=mla -materialFile=ellipsoidMlaMaterials.h -backGround=backGround -rbc=rbcNonLocal -kx=.5 -eps0=1 -eps1=2 -useSosupDissipation=1 -tf=2 -tp=.1 -ic=gp -x0=-1.0 -beta=120 -flushFrequency=4 -yContour=0 -zContour=0 -extraCommandsEnd=plotSpheres.h -show=oneSolidEllpsoidMLAOrder4G16.show -go=go > oneSolidEllpsoidMLAOrder44G16.out 






