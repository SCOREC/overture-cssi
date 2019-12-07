#
#  Makefile for solidObjects
#

# Grid generator:
ogen = $(Overture)/bin/ogen

# Note: grids need to be initially fine enough to handle corners
grids:
	$(ogen) -noplot solidObjectsGrid -interp=e -order=2 -factor=4 
	$(ogen) -noplot solidObjectsGrid -interp=e -numGhost=3 -order=4 -factor=6

# to avoid negative volumes at order=4, start at factor=8 *fix me*
rpiGrids:
	$(ogen) -noplot rpiGrid -direction=vertical -interp=e -order=2 -factor=4
	$(ogen) -noplot rpiGrid -direction=vertical -interp=e -numGhost=3 -order=4 -factor=8
