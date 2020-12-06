#
#   Makefile for cgmx interface 
#


ogen = $(Overture)/bin/ogen


twoBoxesInterface:
	$(ogen) -noplot twoBoxesInterface -factor=1 -order=2 -interp=e
	$(ogen) -noplot twoBoxesInterface -factor=2 -order=2 -interp=e
	$(ogen) -noplot twoBoxesInterface -factor=4 -order=2 -interp=e
	$(ogen) -noplot twoBoxesInterface -factor=8 -order=2 -interp=e
	$(ogen) -noplot twoBoxesInterface -factor=1 -order=4 -interp=e
	$(ogen) -noplot twoBoxesInterface -factor=2 -order=4 -interp=e
	$(ogen) -noplot twoBoxesInterface -factor=4 -order=4 -interp=e
	$(ogen) -noplot twoBoxesInterface -factor=8 -order=4 -interp=e

# Shorter domain 
tbiGrid:
	$(ogen) -noplot twoBoxesInterface -prefix=tbiGrid -xa=-.25 -xb=.25 -factor=1 -order=2 -interp=e
	$(ogen) -noplot twoBoxesInterface -prefix=tbiGrid -xa=-.25 -xb=.25 -factor=2 -order=2 -interp=e
	$(ogen) -noplot twoBoxesInterface -prefix=tbiGrid -xa=-.25 -xb=.25 -factor=4 -order=2 -interp=e
	$(ogen) -noplot twoBoxesInterface -prefix=tbiGrid -xa=-.25 -xb=.25 -factor=8 -order=2 -interp=e
	$(ogen) -noplot twoBoxesInterface -prefix=tbiGrid -xa=-.25 -xb=.25 -factor=1 -order=4 -interp=e -numGhost=3
	$(ogen) -noplot twoBoxesInterface -prefix=tbiGrid -xa=-.25 -xb=.25 -factor=2 -order=4 -interp=e -numGhost=3
	$(ogen) -noplot twoBoxesInterface -prefix=tbiGrid -xa=-.25 -xb=.25 -factor=4 -order=4 -interp=e -numGhost=3
	$(ogen) -noplot twoBoxesInterface -prefix=tbiGrid -xa=-.25 -xb=.25 -factor=8 -order=4 -interp=e -numGhost=3
