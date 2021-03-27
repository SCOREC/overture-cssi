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



# Ellipsoid grids:
ellipsoidGrid1:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=1SolidEllipsoidsGrid -nsx=1 -nsy=1 -nsz=1 -configFile=solidEllipsoidArrayConfig3d.h -order=2 -per=npp -factor=4 -rFactor=.65 -nrExtra=3 -orthographicPatchParameter=.7
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=1SolidEllipsoidsGrid -nsx=1 -nsy=1 -nsz=1 -configFile=solidEllipsoidArrayConfig3d.h -order=4 -per=npp -factor=4 -rFactor=.7 -nrExtra=4 -orthographicPatchParameter=.7 -interp=e

ellipsoidGrid4:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=fourSolidEllipsoidsGrid -configFile=solidEllipsoidArrayConfig4.h -order=2 -per=npp -rFactor=.65 -nrExtra=3 -orthographicPatchParameter=.7 -factor=4

ellipsoidGrid8:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=eightSolidEllipsoidsGrid -configFile=solidEllipsoidArrayConfig8.h -order=4 -per=npp -rFactor=.65 -nrExtra=4 -orthographicPatchParameter=.6 -factor=4 -matchGridsAtInterfaces=1 -interp=i -numGhost=3
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=eightSolidEllipsoidsGrid -configFile=solidEllipsoidArrayConfig8.h -order=2 -per=npp -rFactor=.65 -nrExtra=3 -orthographicPatchParameter=.7 -factor=4

ellipsoidGrid32:
	$(ogen) -noplot solidEllipsoidArrayGrid -prefix=36SolidEllipsoidsGrid -nsx=3 -nsy=3 -nsz=4 -configFile=solidEllipsoidArrayConfig3d.h -order=2 -per=npp -rFactor=.65 -nrExtra=3 -orthographicPatchParameter=.7 -factor=4


