# CgMp examples -- heated bodies

ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp

# --- two solid squares (coupled heat equations), manufactured solutions ----
twoSquaresHeat: twoSquaresInterfacee1.order2.hdf
	$(cgmp) aa -g=twoSquaresInterfacee1.order2.hdf -kappa1=1. -ktc1=1. -kappa2=.5 -ktc2=.5 -tf=.05 -tp=.01

twoSquaresInterfacee1.order2.hdf:
	$(ogen) -noplot twoSquaresInterface -interp=e -order=2 -factor=2 

# ---- One hot disk in an INS fluid -----
oneDiskINS: innerOutere2.order2.hdf
	$(cgmp) io.cmd -g=innerOutere2.order2.hdf -nu=.1 -kappa=.05 -tf=1.0 -solver=yale
# grid for one disk 
innerOutere2.order2.hdf:
	$(ogen) -noplot io -order=2 -interp=e -factor=2


# ----- Four hot disks, INS fluid -----
fourDiskINS: diskArray2x2ye1.order2.hdf
	$(cgmp) io.cmd -g=diskArray2x2ye1.order2.hdf -nu=.05 -kappa=.01 -tp=.05 -solver=yale
# grid for 4 disks
diskArray2x2ye1.order2.hdf:
	$(ogen) -noplot diskArray -factor=1 -interp=e -nCylx=2 -nCyly=2

# ----- 16 hot disks, INS fluid -----
16DiskINS: diskArray4x4ye2.order2.hdf
	$(cgmp) -noplot io.cmd -g=diskArray4x4ye2.order2.hdf -nu=.01 -kappa=.005 -tf=3 -tp=.05 -show=16HeatedDisks.show -go=go > 16HeatedDisks.out
# grid for 16 disks
diskArray4x4ye2.order2.hdf.hdf:
	$(ogen) -noplot diskArray -interp=e -nCylx=4 -nCyly=4 -factor=2

