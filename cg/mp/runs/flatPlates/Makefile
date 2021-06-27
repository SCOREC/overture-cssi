# CgMp examples -- flat plates CHT examples

ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp

# ------- two flat plates in 2D ------
#  Implicit time-stepping with sub-iterations for interface
flatPlates: 
	$(cgmp) flatPlates -H=1 -L=1 -factor=1 -T_top=2 -u_top=1 -ts=imp -imp=1 -nc=10 -tf=50 -tp=1 -dtMax=.1 -fluid_kappa=.02 -fluid_k=.02 -coupled=0 -iTol=1.e-3 -iOmega=.75  -show=flatPlates2d.show


# ------- two flat plates in 2D =, add GRAVITY ------
#  --> u develops a parabolic profile
#  Implicit time-stepping with sub-iterations for interface
flatPlatesWithGravity: 
	$(cgmp) flatPlates -H=1 -L=1 -factor=1 -T_top=5 -u_top=0 -ts=imp -imp=1 -nc=10 -tf=50 -tp=1 -dtMax=.1 -solid_kappa=1. -solid_k=1. -fluid_kappa=.2 -fluid_k=.2 -coupled=0 -iTol=1.e-3 -iOmega=.75 -gravity="-1. 0 0." -debug=1 -show=flatPlates1.show

# ------- two flat plates in 3D ------
#  Implicit time-stepping with sub-iterations for interface
flatPlates3d:
	$(cgmp) flatPlates -nd=3 -H=1 -L=1 -factor=1 -T_top=5 -u_top=0 -ts=imp -imp=1 -nc=10 -tf=10 -tp=1 -dtMax=.1 -solid_kappa=1. -solid_k=1. -fluid_kappa=.2 -fluid_k=.2 -coupled=0 -iTol=1.e-3 -iOmega=.75 -gravity="-1. 0 0." -debug=1 -show="flatPlates3d1.show"