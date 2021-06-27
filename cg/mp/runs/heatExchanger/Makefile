# CgMp examples -- heated exchanger

ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp

# ------- 3D pipe in a box ------
#  Explicit time-stepping (predictor corrector):
#  Solve CHT interface equations using "coupled" compatibility conditions
hexExplicit: heatExchangerGride2.order2.hdf
	$(cgmp) -noplot heatExchanger.cmd -g=heatExchangerGride2.order2.hdf -solver=best -show=hexExplicit.show -nu=.025 -tp=.1 -tf=1.0 -go=go > hexExplicit.out

#  Partitioned implicit time-stepping:
#  Solve partitioned CHT interface equations using an iteration (DN ?)
hexImplicit: heatExchangerGride2.order2.hdf
	$(cgmp) heatExchanger.cmd -g=heatExchangere2.order2.hdf -solver=best -ts=im -nc=20 -dtMax=.05 -tp=.05 -coupled=0 -iOmega=.5  -debug=3

# grid for 3d pipe in a box
heatExchangerGride2.order2.hdf:
	$(ogen) -noplot heatExchangerGrid -factor=2 -interp=e

