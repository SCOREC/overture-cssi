# CgCns example -- shock hitting a deforming beam

ogen        = $(Overture)/bin/ogen
cgcns       = $(CGBUILDPREFIX)/cns/bin/cgcns
beamGridDir = $(CG)/ins/runs/beamInAChannel


# beam in a channel grid:
grids:
	$(ogen) -noplot $(beamGridDir)/beamInAChannelGrid -interp=e -xb=2 -yb=2 -factor=2
	$(ogen) -noplot $(beamGridDir)/beamInAChannelGrid -interp=e -xb=2 -yb=2 -factor=4
	$(ogen) -noplot $(beamGridDir)/beamInAChannelGrid -interp=e -xb=2 -yb=2 -factor=8
	$(ogen) -noplot $(beamGridDir)/beamInAChannelGrid -interp=e -xb=2 -yb=2 -factor=16
	$(ogen) -noplot $(beamGridDir)/beamInAChannelGrid -interp=e -xb=2 -yb=2 -factor=32

# medium grid run
shockBeam8:
	$(cgcns) -noplot beamInAChannel -g=beamInAChannelGride8.order2.hdf -delta=10.0 -tf=2.0 -tp=.1 -show=beam8.show -go=go > beam8.out 

# Fine grid run
shockBeam32:
	$(cgcns) -noplot beamInAChannel -g=beamInAChannelGride32.order2.hdf -delta=10.0 -tf=2.0 -tp=.1 -show=beam32.show -go=go > beam32.out 

# # ------- two flat plates in 2D ------
# #  Implicit time-stepping with sub-iterations for interface
# flatPlates: 
# 	$(cgmp) flatPlates -H=1 -L=1 -factor=1 -T_top=2 -u_top=1 -ts=imp -imp=1 -nc=10 -tf=50 -tp=1 -dtMax=.1 -fluid_kappa=.02 -fluid_k=.02 -coupled=0 -iTol=1.e-3 -iOmega=.75  -show=flatPlates2d.show

