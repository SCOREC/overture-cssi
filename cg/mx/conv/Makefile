common = ../../common



# The convergence rate script is called conv.p
all:  
	ln -sf $(common)/src/conv.p .

cgmx:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) cgmx

ogen = $(Overture)/bin/ogen

# directory to store interface tables:
interfaceTables = ./interfaceTables

#  -- GDM interface convergence runs
#   Order = 2, 4
#   modeGDM = 
#   Np = 


# ---- GDM DISK in a SQUARE ---
#   Orders = 2,4
#   Modes = -1 (default NR-mode) and 1 (R-mode)
gdmDisk: gdmDiskOrder2  gdmDiskOrder4

gdmDiskOrder2:
	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.5 -tp=.1 -order=2 -name="innerOuterFixede" -dm=gdm -modeGDM=-1 -npv=1 2  -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09  -diss=0 -useSosupDissipation=1 -ng=4
	mv dielectricCylGDMNp1Np2Order2maxNormConvTable.tex $(interfaceTables)
	mv dielectricCylNFDTDOrder2GDMNp1Np2maxNorm.m  $(interfaceTables)
	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.5 -tp=.1 -order=2 -name="innerOuterFixede" -dm=gdm -modeGDM=1 -npv=1 2  -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09  -diss=0 -useSosupDissipation=1 -ng=4
	mv dielectricCylGDMNp1Np2Mode1Order2maxNormConvTable.tex $(interfaceTables)
	mv dielectricCylNFDTDOrder2GDMNp1Np2GdmMode1maxNorm.m $(interfaceTables)

gdmDiskOrder4:
	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.5 -tp=.1 -order=4 -name="innerOuterFixede" -dm=gdm -modeGDM=-1 -npv=1 2  -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09  -diss=0 -useSosupDissipation=1 -ng=4
	mv dielectricCylGDMNp1Np2Order4maxNormConvTable.tex $(interfaceTables)
	mv dielectricCylNFDTDOrder4GDMNp1Np2maxNorm.m  $(interfaceTables)
	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.5 -tp=.1 -order=4 -name="innerOuterFixede" -dm=gdm -modeGDM=1 -npv=1 2  -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09  -diss=0 -useSosupDissipation=1 -ng=4
	mv dielectricCylGDMNp1Np2Mode1Order4maxNormConvTable.tex $(interfaceTables)
	mv dielectricCylNFDTDOrder4GDMNp1Np2GdmMode1maxNorm.m $(interfaceTables)


# -------------------- 2D PLANAR INTERFACE ---------------------------------------
#   Orders = 2,4
#   Modes = -1 (default NR-mode) and 1 (R-mode)
interfaceGDM2d :interfaceGDM2dOrder2 interfaceGDM2dOrder4

interfaceGDM2dOrder2:
	conv.p interface.conv -dm=gdm -eps1=2. -eps2=3. -npv=3 2 -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09 -kx=2 -ky=1 -ic=pmic -bc=d -cg=twoSquaresInterface -order=2 -useSosupDissipation=1  -tf=.5 -tp=.1 -modeGDM=-1 -ng=4
	mv interfaceNFDTDtwoSquaresInterfaceGDMNp3Np2Order2maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoSquaresInterfaceOrder2GDMNp3Np2maxNorm.m $(interfaceTables)
	conv.p interface.conv -dm=gdm -eps1=2. -eps2=3. -npv=3 2 -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09 -kx=2 -ky=1 -ic=pmic -bc=d -cg=twoSquaresInterface -order=2 -useSosupDissipation=1  -tf=.5 -tp=.1 -modeGDM=1 -ng=4
	mv interfaceNFDTDtwoSquaresInterfaceGDMNp3Np2Mode1Order2maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoSquaresInterfaceOrder2GDMNp3Np2Mode1maxNorm.m $(interfaceTables)

interfaceGDM2dOrder4:
	conv.p interface.conv -dm=gdm -eps1=2. -eps2=3. -npv=3 2 -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09 -kx=2 -ky=1 -ic=pmic -bc=d -cg=twoSquaresInterface -order=4 -useSosupDissipation=1  -tf=.5 -tp=.1 -modeGDM=-1 -ng=4
	mv interfaceNFDTDtwoSquaresInterfaceGDMNp3Np2Order4maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoSquaresInterfaceOrder4GDMNp3Np2maxNorm.m $(interfaceTables)
	conv.p interface.conv -dm=gdm -eps1=2. -eps2=3. -npv=3 2 -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09 -kx=2 -ky=1 -ic=pmic -bc=d -cg=twoSquaresInterface -order=4 -useSosupDissipation=1  -tf=.5 -tp=.1 -modeGDM=1 -ng=4
	mv interfaceNFDTDtwoSquaresInterfaceGDMNp3Np2Mode1Order4maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoSquaresInterfaceOrder4GDMNp3Np2Mode1maxNorm.m $(interfaceTables)

# -------------------- 3D PLANAR INTERFACE ---------------------------------------
interfaceGDM3d: interfaceGDM3dOrder2 interfaceGDM3dOrder4


interfaceGDM3dOrder2:
	conv.p interface.conv -kx=2 -ky=1 -kz=.5 -ic=pmic -eps1=1 -eps2=2. -dm=gdm -npv=2 1 -a01=1 1.5 -a11=.1 .125 -b01=1 1.5 -b11=.2 .1 -a02=2 -a12=.05 -b02=2 -b12=.1 -ic=pmic -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=twoBoxesInterface -tf=.5 -tp=.1 -order=2 -useSosupDissipation=1 -modeGDM=-1 -ng=4
	mv interfaceNFDTDtwoBoxesInterfaceGDMNp2Np1Order2maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoBoxesInterfaceOrder2GDMNp2Np1maxNorm.m  $(interfaceTables)
	conv.p interface.conv -kx=2 -ky=1 -kz=.5 -ic=pmic -eps1=1 -eps2=2. -dm=gdm -npv=2 1 -a01=1 1.5 -a11=.1 .125 -b01=1 1.5 -b11=.2 .1 -a02=2 -a12=.05 -b02=2 -b12=.1 -ic=pmic -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=twoBoxesInterface -tf=.5 -tp=.1 -order=2 -useSosupDissipation=1 -modeGDM=1 -ng=4
	mv interfaceNFDTDtwoBoxesInterfaceGDMNp2Np1Mode1Order2maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoBoxesInterfaceOrder2GDMNp2Np1Mode1maxNorm.m  $(interfaceTables)

interfaceGDM3dOrder4:
	conv.p interface.conv -kx=2 -ky=1 -kz=.5 -ic=pmic -eps1=1 -eps2=2. -dm=gdm -npv=2 1 -a01=1 1.5 -a11=.1 .125 -b01=1 1.5 -b11=.2 .1 -a02=2 -a12=.05 -b02=2 -b12=.1 -ic=pmic -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=twoBoxesInterface -tf=.5 -tp=.1 -order=4 -useSosupDissipation=1 -modeGDM=-1 -ng=4
	mv interfaceNFDTDtwoBoxesInterfaceGDMNp2Np1Order4maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoBoxesInterfaceOrder4GDMNp2Np1maxNorm.m  $(interfaceTables)
	conv.p interface.conv -kx=2 -ky=1 -kz=.5 -ic=pmic -eps1=1 -eps2=2. -dm=gdm -npv=2 1 -a01=1 1.5 -a11=.1 .125 -b01=1 1.5 -b11=.2 .1 -a02=2 -a12=.05 -b02=2 -b12=.1 -ic=pmic -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=twoBoxesInterface -tf=.5 -tp=.1 -order=4 -useSosupDissipation=1 -modeGDM=1 -ng=4
	mv interfaceNFDTDtwoBoxesInterfaceGDMNp2Np1Mode1Order4maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoBoxesInterfaceOrder4GDMNp2Np1Mode1maxNorm.m  $(interfaceTables)

# OLD: 
# --- DOUBLE-BOX -- grid is slightly finer in x-direction to over-come problem with 3d order 4 interface
# wdh: This issue has now been fixed
interfaceGDM3dOrder2old:
	conv.p interface.conv -dm=gdm -npv=1 1 -kx=2 -ky=1 -kz=0 -ic=pmic -eps1=1 -eps2=2.  -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=doubleBoxInterface -tf=.5 -tp=.1 -order=2 -useSosupDissipation=1 -modeGDM=-1 -ng=4
	mv interfaceNFDTDdoubleBoxInterfaceGDMNp1Np1Order2maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDdoubleBoxInterfaceOrder2GDMNp1Np1maxNorm.m  $(interfaceTables)
	conv.p interface.conv -dm=gdm -npv=1 1 -kx=2 -ky=1 -kz=0 -ic=pmic -eps1=1 -eps2=2.  -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=doubleBoxInterface -tf=.5 -tp=.1 -order=2 -useSosupDissipation=1 -modeGDM=0 -ng=4
	mv interfaceNFDTDdoubleBoxInterfaceGDMNp1Np1Mode0Order2maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDdoubleBoxInterfaceOrder2GDMNp1Np1Mode0maxNorm.m  $(interfaceTables)

interfaceGDM3dOrder4old:
	conv.p interface.conv -dm=gdm -npv=1 1 -kx=2 -ky=1 -kz=0 -ic=pmic -eps1=1 -eps2=2.  -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=doubleBoxInterface -tf=.5 -tp=.1 -order=4 -useSosupDissipation=1 -modeGDM=-1 -ng=4
	mv interfaceNFDTDdoubleBoxInterfaceGDMNp1Np1Order4maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDdoubleBoxInterfaceOrder4GDMNp1Np1maxNorm.m  $(interfaceTables)
	conv.p interface.conv -dm=gdm -npv=1 1 -kx=2 -ky=1 -kz=0 -ic=pmic -eps1=1 -eps2=2.  -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=doubleBoxInterface -tf=.5 -tp=.1 -order=4 -useSosupDissipation=1 -modeGDM=0 -ng=4
	mv interfaceNFDTDdoubleBoxInterfaceGDMNp1Np1Mode0Order4maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDdoubleBoxInterfaceOrder4GDMNp1Np1Mode0maxNorm.m  $(interfaceTables)


# -------------- SCATTERING FROM A DIELETRIC SPHERE ------------
# TEST ME  
dielectricSphere2:
	conv.p dielectric.conv -tp=.1 -tf=.5 -order=2 -name="solidSphereInABoxe" -kx=1 -ky=0 -eps1=2 -eps2=1 -dm=gdm -modeGDM=-1 -npv=1 2  -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09  -cfl=.9  -useSosupDissipation=1 -ng=3
	mv dielectricSphereGDMNp1Np2Order2maxNormConvTable.tex $(interfaceTables)
	mv dielectricSphereNFDTDOrder2GDMNp1Np2maxNorm.m $(interfaceTables)
	conv.p dielectric.conv -tp=.1 -tf=.5 -order=2 -name="solidSphereInABoxe" -kx=1 -ky=0 -eps1=2 -eps2=1 -dm=gdm -modeGDM=1 -npv=1 2  -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09  -cfl=.9  -useSosupDissipation=1 -ng=3
	mv dielectricSphereGDMNp1Np2Mode1Order2maxNormConvTable.tex $(interfaceTables)
	mv dielectricSphereNFDTDOrder2GDMNp1Np2GdmMode1maxNorm.m $(interfaceTables) 


dielectricSphere4:
	conv.p dielectric.conv -tp=.1 -tf=.5 -order=4 -name="solidSphereInABoxe" -kx=1 -ky=0 -eps1=2 -eps2=1 -dm=gdm -modeGDM=-1 -npv=1 2  -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09  -cfl=.9  -useSosupDissipation=1 -ng=3
	mv dielectricSphereGDMNp1Np2Order4maxNormConvTable.tex $(interfaceTables)
	mv dielectricSphereNFDTDOrder4GDMNp1Np2maxNorm.m $(interfaceTables)
	conv.p dielectric.conv -tp=.1 -tf=.5 -order=4 -name="solidSphereInABoxe" -kx=1 -ky=0 -eps1=2 -eps2=1 -dm=gdm -modeGDM=1 -npv=1 2  -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09  -cfl=.9  -useSosupDissipation=1 -ng=3
	mv dielectricSphereGDMNp1Np2Mode1Order4maxNormConvTable.tex $(interfaceTables)
	mv dielectricSphereNFDTDOrder4GDMNp1Np2GdmMode1maxNorm.m $(interfaceTables) 



# -------- BA MAXWELL INTERFACE *WARNING* COPIES FILES TO bamx paper directory
bamxp = ~/Dropbox/DARPA/BAMX/papers/bamx
matFile = baMatGDMCase1
matFile2 = baMatGDM4
baInterface: baInterface2d baInterface3d

baInterface2d:
	conv.p interface.conv -cg=rect -method=bamx -kx=1 -ky=1 -dm=gdm -dispersionOption=k -matFile=$(matFile).txt -matFile2=$(matFile2).txt  -ic=pmic -tf=1 -tp=.5 -order=2 -saveVectorErrors=1  -errorNorm=1 -option=interfaceWithDiss -diss=1 -alwaysCompute=0 -ng=7
	cp interfaceWithDissbamxrect$(matFile)$(matFile2)rk4Order2L1NormConvTable.tex $(bamxp)/tables
	cp interfaceWithDissbamxrectrk4Order2$(matFile)$(matFile2)L1Norm.m $(bamxp)/matlab 
	conv.p interface.conv -cg=rect -method=bamx -kx=1 -ky=1 -dm=gdm -dispersionOption=k -matFile=$(matFile).txt -matFile2=$(matFile2).txt  -ic=pmic -tf=1 -tp=.5 -order=4 -saveVectorErrors=1  -errorNorm=1 -option=interfaceWithDiss -diss=1 -alwaysCompute=0 -ng=7
	cp interfaceWithDissbamxrect$(matFile)$(matFile2)rk4Order4L1NormConvTable.tex $(bamxp)/tables
	cp interfaceWithDissbamxrectrk4Order4$(matFile)$(matFile2)L1Norm.m $(bamxp)/matlab 

baInterface2dNoDiss:
	conv.p interface.conv -cg=rect -method=bamx -kx=1 -ky=1 -dm=gdm -dispersionOption=k -matFile=$(matFile).txt -matFile2=$(matFile2).txt  -ic=pmic -tf=1 -tp=.5 -order=2 -saveVectorErrors=1  -errorNorm=1  -alwaysCompute=0 -ng=7
	cp interfacebamxrect$(matFile)$(matFile2)rk4Order2L1NormConvTable.tex $(bamxp)/tables
	cp interfacebamxrectrk4Order2$(matFile)$(matFile2)L1Norm.m $(bamxp)/matlab 
	conv.p interface.conv -cg=rect -method=bamx -kx=1 -ky=1 -dm=gdm -dispersionOption=k -matFile=$(matFile).txt -matFile2=$(matFile2).txt  -ic=pmic -tf=1 -tp=.5 -order=4 -saveVectorErrors=1  -errorNorm=1 -alwaysCompute=0 -ng=7
	cp interfacebamxrect$(matFile)$(matFile2)rk4Order4L1NormConvTable.tex $(bamxp)/tables
	cp interfacebamxrectrk4Order4$(matFile)$(matFile2)L1Norm.m $(bamxp)/matlab 

baInterface3d:
	conv.p interface.conv -cg=box3d -method=bamx -kx=1 -ky=.5 -kz=.5 -dm=gdm -dispersionOption=k -matFile=$(matFile).txt -matFile2=$(matFile2).txt -ic=pmic -tf=.5 -tp=.5 -order=2 -saveVectorErrors=1  -errorNorm=1 -alwaysCompute=0 -ng=4
	cp interfacebamxbox3d$(matFile)$(matFile2)rk4Order2L1NormConvTable.tex $(bamxp)/tables
	cp interfacebamxbox3drk4Order2$(matFile)$(matFile2)L1Norm.m $(bamxp)/matlab 
	conv.p interface.conv -cg=box3d -method=bamx -kx=1 -ky=.5 -kz=.5 -dm=gdm -dispersionOption=k -matFile=$(matFile).txt -matFile2=$(matFile2).txt -ic=pmic -tf=.5 -tp=.5 -order=2 -saveVectorErrors=1  -errorNorm=1 -option=interfaceWithDiss -diss=1 -alwaysCompute=0 -ng=4
	cp interfaceWithDissbamxbox3d$(matFile)$(matFile2)rk4Order2L1NormConvTable.tex $(bamxp)/tables
	cp interfaceWithDissbamxbox3drk4Order2$(matFile)$(matFile2)L1Norm.m $(bamxp)/matlab 
	conv.p interface.conv -cg=box3d -method=bamx -kx=1 -ky=.5 -kz=.5 -dm=gdm -dispersionOption=k -matFile=$(matFile).txt -matFile2=$(matFile2).txt -ic=pmic -tf=.5 -tp=.5 -order=4 -saveVectorErrors=1  -errorNorm=1 -alwaysCompute=0 -ng=4
	cp interfacebamxbox3d$(matFile)$(matFile2)rk4Order4L1NormConvTable.tex $(bamxp)/tables
	cp interfacebamxbox3drk4Order4$(matFile)$(matFile2)L1Norm.m $(bamxp)/matlab 
	conv.p interface.conv -cg=box3d -method=bamx -kx=1 -ky=.5 -kz=.5 -dm=gdm -dispersionOption=k -matFile=$(matFile).txt -matFile2=$(matFile2).txt -ic=pmic -tf=.5 -tp=.5 -order=4 -saveVectorErrors=1  -errorNorm=1 -option=interfaceWithDiss -diss=1 -alwaysCompute=0 -ng=4
	cp interfaceWithDissbamxbox3d$(matFile)$(matFile2)rk4Order4L1NormConvTable.tex $(bamxp)/tables
	cp interfaceWithDissbamxbox3drk4Order4$(matFile)$(matFile2)L1Norm.m $(bamxp)/matlab 


# ---------------------------- BA SURFACE WAVE -----------------------
#    **************WARNING THIS WILL OVERWITE FILES IN THE BAMX PAPER DIRECTORY*****************
baSurfaceWave: baSurfaceWaveOrder2 baSurfaceWaveOrder4

baSurfaceWaveOrder2:
	conv.p interface.conv -cg=rect -gridSuffix=pn -method=bamx -dm=gdm -xa=-1.05 -xb=1.05 -ya=-.55 -yb=1e-10 -matFile=baMatIsotropic1.txt -matFile2=baMatBiIsotropic1.txt  -ic=isw -tf=.5 -tp=.5 -order=2 -dispersionOption=k  -saveVectorErrors=1 -errorNorm=1 -ng=5
	cp interfaceSurfaceWavebamxrectbaMatIsotropic1baMatBiIsotropic1rk4Order2L1NormConvTable.tex $(bamxp)/tables
	cp interfaceSurfaceWavebamxrectrk4Order2baMatIsotropic1baMatBiIsotropic1L1Norm.m  $(bamxp)/matlab 
	conv.p interface.conv -cg=rect -gridSuffix=pn -method=bamx -dm=gdm -xa=-1.05 -xb=1.05 -ya=-.55 -yb=1e-10 -matFile=baMatIsotropic1.txt -matFile2=baMatBiIsotropic1.txt  -ic=isw -tf=.5 -tp=.5 -order=2 -dispersionOption=k  -saveVectorErrors=1 -errorNorm=1 -option=interfaceWithDiss -diss=1 -ng=5
	cp interfaceWithDissSurfaceWavebamxrectbaMatIsotropic1baMatBiIsotropic1rk4Order2L1NormConvTable.tex  $(bamxp)/tables
	cp interfaceWithDissSurfaceWavebamxrectrk4Order2baMatIsotropic1baMatBiIsotropic1L1Norm.m $(bamxp)/matlab


baSurfaceWaveOrder4:
	conv.p interface.conv -cg=rect -gridSuffix=pn -method=bamx -dm=gdm -xa=-1.05 -xb=1.05 -ya=-.55 -yb=1e-10 -matFile=baMatIsotropic1.txt -matFile2=baMatBiIsotropic1.txt  -ic=isw -tf=.5 -tp=.5 -order=4 -dispersionOption=k  -saveVectorErrors=1 -errorNorm=1 -ng=5
	cp interfaceSurfaceWavebamxrectbaMatIsotropic1baMatBiIsotropic1rk4Order4L1NormConvTable.tex $(bamxp)/tables
	cp interfaceSurfaceWavebamxrectrk4Order4baMatIsotropic1baMatBiIsotropic1L1Norm.m  $(bamxp)/matlab 
	conv.p interface.conv -cg=rect -gridSuffix=pn -method=bamx -dm=gdm -xa=-1.05 -xb=1.05 -ya=-.55 -yb=1e-10 -matFile=baMatIsotropic1.txt -matFile2=baMatBiIsotropic1.txt  -ic=isw -tf=.5 -tp=.5 -order=4 -dispersionOption=k  -saveVectorErrors=1 -errorNorm=1 -option=interfaceWithDiss -diss=1 -ng=5
	cp interfaceWithDissSurfaceWavebamxrectbaMatIsotropic1baMatBiIsotropic1rk4Order4L1NormConvTable.tex  $(bamxp)/tables
	cp interfaceWithDissSurfaceWavebamxrectrk4Order4baMatIsotropic1baMatBiIsotropic1L1Norm.m $(bamxp)/matlab


# --------------------------- SCATTERING FROM A CYLINDER WITH BAMX --------------------
baCyl: baCylOrder2 baCylOrder4

baCylOrder2:
	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.6 -tp=.6 -order=2 -name="cSquaree" -method=bamx -ts=rk4 -matFile=baAir.txt -matFile2=baIsoEps0p25.txt  -errorNorm=1 -diss=0 -alwaysCompute=0 -ng=7
	cp dielectricCylbaAirbaIsoEps0p25Order2NormConvTable.tex $(bamxp)/tables
	cp dielectricCylbamxOrder2L1Norm.m $(bamxp)/matlab 
	conv.p dielectric.conv -prefix="WithDiss" -cyl=1 -kx=3 -tf=.6 -tp=.6 -order=2 -name="cSquaree" -method=bamx -ts=rk4 -matFile=baAir.txt -matFile2=baIsoEps0p25.txt  -errorNorm=1 -diss=1 -alwaysCompute=0 -ng=7
	cp dielectricCylWithDissbaAirbaIsoEps0p25Order2NormConvTable.tex $(bamxp)/tables
	cp dielectricCylWithDissbamxOrder2L1Norm.m $(bamxp)/matlab 

baCylOrder4:
	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.6 -tp=.6 -order=4 -name="cSquaree" -method=bamx -ts=rk4 -matFile=baAir.txt -matFile2=baIsoEps0p25.txt  -errorNorm=1 -diss=0 -alwaysCompute=0 -ng=7
	cp dielectricCylbaAirbaIsoEps0p25Order4NormConvTable.tex $(bamxp)/tables
	cp dielectricCylbamxOrder4L1Norm.m $(bamxp)/matlab 
	conv.p dielectric.conv -prefix="WithDiss" -cyl=1 -kx=3 -tf=.6 -tp=.6 -order=4 -name="cSquaree" -method=bamx -ts=rk4 -matFile=baAir.txt -matFile2=baIsoEps0p25.txt  -errorNorm=1 -diss=1 -alwaysCompute=0 -ng=7
	cp dielectricCylWithDissbaAirbaIsoEps0p25Order4NormConvTable.tex $(bamxp)/tables
	cp dielectricCylWithDissbamxOrder4L1Norm.m $(bamxp)/matlab 



# --------------------------- SCATTERING FROM A SPHERE WITH BAMX --------------------
baSphere: baSphereOrder2 baSphereOrder4

baSphereOrder2:
	conv.p dielectric.conv -cyl=0 -kx=1 -tf=.6 -tp=.6 -order=2 -name="cboxe" -method=bamx -ts=rk4 -matFile=baAir.txt -matFile2=baIsoEps0p25.txt -errorNorm=1 -diss=0 -saveVectorErrors=1 -solveForAllFields=1 -alwaysCompute=0 -ng=4
	cp dielectricSpherebaAirbaIsoEps0p25Order2NormConvTable.tex $(bamxp)/tables
	cp dielectricSpherebamxOrder2L1Norm.m $(bamxp)/matlab
	conv.p dielectric.conv -cyl=0 -kx=1 -tf=.6 -tp=.6 -order=2 -name="cboxe" -method=bamx -ts=rk4 -matFile=baAir.txt -matFile2=baIsoEps0p25.txt -errorNorm=1 -prefix="WithDiss" -diss=1 -saveVectorErrors=1 -solveForAllFields=1 -alwaysCompute=0 -ng=4
	cp dielectricSphereWithDissbaAirbaIsoEps0p25Order2NormConvTable.tex $(bamxp)/tables
	cp dielectricSphereWithDissbamxOrder2L1Norm.m $(bamxp)/matlab


baSphereOrder4:
	conv.p dielectric.conv -cyl=0 -kx=1 -tf=.6 -tp=.6 -order=4 -name="cboxe" -method=bamx -ts=rk4 -matFile=baAir.txt -matFile2=baIsoEps0p25.txt -errorNorm=1 -diss=0 -saveVectorErrors=1 -solveForAllFields=1 -alwaysCompute=0 -ng=4
	cp dielectricSpherebaAirbaIsoEps0p25Order4NormConvTable.tex $(bamxp)/tables
	cp dielectricSpherebamxOrder4L1Norm.m $(bamxp)/matlab
	conv.p dielectric.conv -cyl=0 -kx=1 -tf=.6 -tp=.6 -order=4 -name="cboxe" -method=bamx -ts=rk4 -matFile=baAir.txt -matFile2=baIsoEps0p25.txt -errorNorm=1 -prefix="WithDiss" -diss=1 -saveVectorErrors=1 -solveForAllFields=1 -alwaysCompute=0 -ng=4
	cp dielectricSphereWithDissbaAirbaIsoEps0p25Order4NormConvTable.tex $(bamxp)/tables
	cp dielectricSphereWithDissbamxOrder4L1Norm.m $(bamxp)/matlab


# ------------------------ MLA INTERFACE -----------------

interfaceMLA: interfaceMLA2d interfaceMLA3d copyInterfaceMLA

interfaceMLA2d:
	conv.p interface.conv -dm=gdm -nm=mla -matFile1=mlaMat2.txt -matFile2=mlaMat3.txt -ic=tz -tz=trig -bc=d -cg=twoSquaresInterface -order=2 -useSosupDissipation=1 -tf=.3 -tp=.3 -ng=4
	processCheckFiles.p -file=./interfaceNFDTDtwoSquaresInterfaceOrder2MLAmaxNorm.twoSquaresInterface.master.check -useGridNames=0 -c=3,4,5 -startingGridResolution=20
	conv.p interface.conv -dm=gdm -nm=mla -matFile1=mlaMat2.txt -matFile2=mlaMat3.txt -ic=tz -tz=trig -bc=d -cg=twoSquaresInterface -order=4 -useSosupDissipation=1 -tf=.3 -tp=.3 -ng=4
	processCheckFiles.p -file=./interfaceNFDTDtwoSquaresInterfaceOrder4MLAmaxNorm.twoSquaresInterface.master.check -useGridNames=0 -c=3,4,5 -startingGridResolution=20

interfaceMLA3d:
	conv.p interface.conv -dm=gdm -nm=mla -matFile1=mlaMat2.txt -matFile2=mlaMat3.txt -ic=tz -tz=trig -bc=d -cg=twoNonBoxesInterface -order=2 -useSosupDissipation=1 -tf=.3 -tp=.3 -ng=4
	processCheckFiles.p -file=./interfaceNFDTDtwoNonBoxesInterfaceOrder2MLAmaxNorm.twoNonBoxesInterface.master.check -useGridNames=0 -c=3,4,5 -startingGridResolution=20
	conv.p interface.conv -dm=gdm -nm=mla -matFile1=mlaMat2.txt -matFile2=mlaMat3.txt -ic=tz -tz=trig -bc=d -cg=twoNonBoxesInterface -order=4 -useSosupDissipation=1 -tf=.3 -tp=.3 -ng=4
	processCheckFiles.p -file=./interfaceNFDTDtwoNonBoxesInterfaceOrder4MLAmaxNorm.twoNonBoxesInterface.master.check -useGridNames=0 -c=3,4,5 -startingGridResolution=20

# Copy files to Maxwell-Block paper
mbe = ~/Dropbox/DARPA/RPI/Maxwell-Bloch/papers/mbe
copyInterfaceMLA:
	cp interfaceNFDTDtwoSquaresInterfaceOrder2MLAmaxNorm.twoSquaresInterface.master.ConvTable.tex $(mbe)/tables
	cp interfaceNFDTDtwoSquaresInterfaceOrder2MLAmaxNorm.m $(mbe)/codes
	cp interfaceNFDTDtwoSquaresInterfaceOrder4MLAmaxNorm.twoSquaresInterface.master.ConvTable.tex $(mbe)/tables
	cp interfaceNFDTDtwoSquaresInterfaceOrder4MLAmaxNorm.m $(mbe)/codes
	cp interfaceNFDTDtwoNonBoxesInterfaceOrder2MLAmaxNorm.twoNonBoxesInterface.master.ConvTable.tex $(mbe)/tables
	cp interfaceNFDTDtwoNonBoxesInterfaceOrder2MLAmaxNorm.m $(mbe)/codes
	cp interfaceNFDTDtwoNonBoxesInterfaceOrder4MLAmaxNorm.twoNonBoxesInterface.master.ConvTable.tex $(mbe)/tables
	cp interfaceNFDTDtwoNonBoxesInterfaceOrder4MLAmaxNorm.m $(mbe)/codes

# ------------------------ CURVED MLA INTERFACE -----------------

# one disk for MLA self convergence
oneDiskGrids:
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=2 -factor=1
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=2 -factor=2
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=2 -factor=4
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=2 -factor=8
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=2 -factor=16
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=4 -numGhost=3 -factor=1
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=4 -numGhost=3 -factor=2
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=4 -numGhost=3 -factor=4
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=4 -numGhost=3 -factor=8
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=4 -numGhost=3 -factor=16

# one-disk MLA (trig)
oneDiskMLA:
	conv.p interface.conv -dm=gdm -nm=mla -matFile1=mlaMat2.txt -matFile2=mlaMat3.txt -ic=tz -tz=trig -bc=d -cg=oneDiskMLAGridFixed -leftDomain=outerDomain -rightDomain=innerDomain1 -order=2 -useSosupDissipation=1 -tf=.3 -tp=.3 -ng=3
	processCheckFiles.p -file=./interfaceNFDTDoneDiskMLAGridFixedOrder2MLAmaxNorm.oneDiskMLAGridFixed.master.check -useGridNames=0 -c=3,4,5 -startingGridResolution=20
	conv.p interface.conv -dm=gdm -nm=mla -matFile1=mlaMat2.txt -matFile2=mlaMat3.txt -ic=tz -tz=trig -bc=d -cg=oneDiskMLAGridFixed -leftDomain=outerDomain -rightDomain=innerDomain1 -order=4 -useSosupDissipation=1 -tf=.3 -tp=.3 -ng=3
	processCheckFiles.p -file=./interfaceNFDTDoneDiskMLAGridFixedOrder4MLAmaxNorm.oneDiskMLAGridFixed.master.check -useGridNames=0 -c=3,4,5 -startingGridResolution=20


copyOneDiskMLA:
	cp interfaceNFDTDoneDiskMLAGridFixedOrder2MLAmaxNorm.oneDiskMLAGridFixed.master.ConvTable.tex $(mbe)/tables
	cp interfaceNFDTDoneDiskMLAGridFixedOrder2MLAmaxNorm.m $(mbe)/codes
	cp interfaceNFDTDoneDiskMLAGridFixedOrder4MLAmaxNorm.oneDiskMLAGridFixed.master.ConvTable.tex $(mbe)/tables
	cp interfaceNFDTDoneDiskMLAGridFixedOrder4MLAmaxNorm.m $(mbe)/codes
