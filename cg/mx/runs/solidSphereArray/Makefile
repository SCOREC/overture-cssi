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





