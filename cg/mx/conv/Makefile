common = ../../common

# The convergence rate script is called conv.p
all:  
	ln -sf $(common)/src/conv.p .

cgmx:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) cgmx


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
