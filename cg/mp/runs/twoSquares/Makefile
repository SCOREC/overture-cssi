# CgMp examples -- heat transfer in multiple solid domains

ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp


# double annulus grids
doubleAnnulusGrids:
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleAnnulusGrid -order=2 -interp=e -factor=1
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleAnnulusGrid -order=2 -interp=e -factor=2
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleAnnulusGrid -order=2 -interp=e -factor=4
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleAnnulusGrid -order=2 -interp=e -factor=8
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleAnnulusGrid -order=2 -interp=e -factor=16
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleAnnulusGrid -order=2 -interp=e -factor=32
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleAnnulusGrid -order=4 -interp=e -factor=1
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleAnnulusGrid -order=4 -interp=e -factor=2
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleAnnulusGrid -order=4 -interp=e -factor=4
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleAnnulusGrid -order=4 -interp=e -factor=8
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleAnnulusGrid -order=4 -interp=e -factor=16
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleAnnulusGrid -order=4 -interp=e -factor=32	


# double annulus, just a sector
doubleEighthAnnulusGrids:
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleEighthAnnulusGrid -theta2=.125 -order=2 -interp=e -factor=1
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleEighthAnnulusGrid -theta2=.125 -order=2 -interp=e -factor=2
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleEighthAnnulusGrid -theta2=.125 -order=2 -interp=e -factor=4
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleEighthAnnulusGrid -theta2=.125 -order=2 -interp=e -factor=8
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleEighthAnnulusGrid -theta2=.125 -order=2 -interp=e -factor=16
	$(ogen) -noplot doubleAnnulusGrid -prefix=doubleEighthAnnulusGrid -theta2=.125 -order=2 -interp=e -factor=32

# two squares or non-squares in a horizontal orientation
twoSquaresInterface:
	$(ogen) -noplot twoSquaresInterface  -interp=e -order=2 -factor=1 
	$(ogen) -noplot twoSquaresInterface  -interp=e -order=2 -factor=2 
	$(ogen) -noplot twoSquaresInterface  -interp=e -order=2 -factor=4
	$(ogen) -noplot twoSquaresInterface  -interp=e -angle=0 -order=2 -factor=1 
	$(ogen) -noplot twoSquaresInterface  -interp=e -angle=0 -order=2 -factor=2 
	$(ogen) -noplot twoSquaresInterface  -interp=e -angle=0 -order=2 -factor=4
	$(ogen) -noplot twoSquaresInterface  -interp=e -order=4 -factor=1 
	$(ogen) -noplot twoSquaresInterface  -interp=e -order=4 -factor=2 
	$(ogen) -noplot twoSquaresInterface  -interp=e -order=4 -factor=4
	$(ogen) -noplot twoSquaresInterface  -interp=e -angle=0 -order=4 -factor=1 
	$(ogen) -noplot twoSquaresInterface  -interp=e -angle=0 -order=4 -factor=2 
	$(ogen) -noplot twoSquaresInterface  -interp=e -angle=0 -order=4 -factor=4	

# two horizontal squares, right-side finer in dx
twoHorizontalSquares:
	$(ogen) -noplot twoSquaresInterface -prefix=twoHorizontalSquares -xFactorRight=2 -interp=e -order=2 -factor=1 
	$(ogen) -noplot twoSquaresInterface -prefix=twoHorizontalSquares -xFactorRight=2 -interp=e -order=2 -factor=2 
	$(ogen) -noplot twoSquaresInterface -prefix=twoHorizontalSquares -xFactorRight=2 -interp=e -order=2 -factor=4
	$(ogen) -noplot twoSquaresInterface -prefix=twoHorizontalNonSquares -xFactorRight=2 -interp=e -angle=0.0 -order=2 -factor=1 
	$(ogen) -noplot twoSquaresInterface -prefix=twoHorizontalNonSquares -xFactorRight=2 -interp=e -angle=0.0 -order=2 -factor=2 
	$(ogen) -noplot twoSquaresInterface -prefix=twoHorizontalNonSquares -xFactorRight=2 -interp=e -angle=0.0 -order=2 -factor=4	

# two squares or non-squares in a vertical orientation
twoSquaresVertical:
	$(ogen) -noplot twoSquaresInterface -prefix=twoVerticalSquares -interp=e -orientation=vertical -order=2 -factor=1 
	$(ogen) -noplot twoSquaresInterface -prefix=twoVerticalSquares -interp=e -orientation=vertical -order=2 -factor=2 
	$(ogen) -noplot twoSquaresInterface -prefix=twoVerticalSquares -interp=e -orientation=vertical -order=2 -factor=4
	$(ogen) -noplot twoSquaresInterface -prefix=twoVerticalNonSquares -interp=e -orientation=vertical -angle=0.0 -order=2 -factor=1 
	$(ogen) -noplot twoSquaresInterface -prefix=twoVerticalNonSquares -interp=e -orientation=vertical -angle=0.0 -order=2 -factor=2 
	$(ogen) -noplot twoSquaresInterface -prefix=twoVerticalNonSquares -interp=e -orientation=vertical -angle=0.0 -order=2 -factor=4 





# # ------- two flat plates in 2D ------
# #  Implicit time-stepping with sub-iterations for interface
# flatPlates: 
# 	$(cgmp) flatPlates -H=1 -L=1 -factor=1 -T_top=2 -u_top=1 -ts=imp -imp=1 -nc=10 -tf=50 -tp=1 -dtMax=.1 -fluid_kappa=.02 -fluid_k=.02 -coupled=0 -iTol=1.e-3 -iOmega=.75  -show=flatPlates2d.show


# # ------- two flat plates in 2D =, add GRAVITY ------
# #  --> u develops a parabolic profile
# #  Implicit time-stepping with sub-iterations for interface
# flatPlatesWithGravity: 
# 	$(cgmp) flatPlates -H=1 -L=1 -factor=1 -T_top=5 -u_top=0 -ts=imp -imp=1 -nc=10 -tf=50 -tp=1 -dtMax=.1 -solid_kappa=1. -solid_k=1. -fluid_kappa=.2 -fluid_k=.2 -coupled=0 -iTol=1.e-3 -iOmega=.75 -gravity="-1. 0 0." -debug=1 -show=flatPlates1.show

# # ------- two flat plates in 3D ------
# #  Implicit time-stepping with sub-iterations for interface
# flatPlates3d:
# 	$(cgmp) flatPlates -nd=3 -H=1 -L=1 -factor=1 -T_top=5 -u_top=0 -ts=imp -imp=1 -nc=10 -tf=10 -tp=1 -dtMax=.1 -solid_kappa=1. -solid_k=1. -fluid_kappa=.2 -fluid_k=.2 -coupled=0 -iTol=1.e-3 -iOmega=.75 -gravity="-1. 0 0." -debug=1 -show="flatPlates3d1.show"