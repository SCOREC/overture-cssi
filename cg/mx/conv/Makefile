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
gdmDisk: gdmDiskOrder2  gdmDiskOrder4

gdmDiskOrder2:
	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.2 -tp=.1 -order=2 -name="innerOuterFixede" -dm=gdm -npv=1 1 -diss=0 -useSosupDissipation=1 -ng=4
	mv dielectricCylGDMNp1Np1Order2maxNormConvTable.tex $(interfaceTables)
	mv dielectricCylNFDTDOrder2GDMNp1Np1maxNorm.m  $(interfaceTables)
	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.2 -tp=.1 -order=2 -name="innerOuterFixede" -dm=gdm -npv=1 1 -diss=0 -useSosupDissipation=1 -modeGDM=0 -ng=4
	mv dielectricCylGDMNp1Np1Mode0Order2maxNormConvTable.tex $(interfaceTables)
	mv dielectricCylNFDTDOrder2GDMNp1Np1GdmMode0maxNorm.m $(interfaceTables)

gdmDiskOrder4:
	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.2 -tp=.1 -order=4 -name="innerOuterFixede" -dm=gdm -npv=1 1 -diss=0 -useSosupDissipation=1 -ng=4
	mv dielectricCylGDMNp1Np1Order4maxNormConvTable.tex $(interfaceTables)
	mv dielectricCylNFDTDOrder4GDMNp1Np1maxNorm.m  $(interfaceTables)
	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.2 -tp=.1 -order=4 -name="innerOuterFixede" -dm=gdm -npv=1 1 -diss=0 -useSosupDissipation=1 -modeGDM=0 -ng=4
	mv dielectricCylGDMNp1Np1Mode0Order4maxNormConvTable.tex $(interfaceTables)
	mv dielectricCylNFDTDOrder4GDMNp1Np1GdmMode0maxNorm.m $(interfaceTables)



interfaceGDM2d :interfaceGDM2dOrder2 interfaceGDM2dOrder4

interfaceGDM2dOrder2:
	conv.p interface.conv -dm=gdm -eps1=1 -eps2=2. -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -kx=2 -ky=1 -ic=pmic -eps1=1 -eps2=2. -bc=d -cg=twoSquaresInterface -order=2 -useSosupDissipation=1  -tf=.5 -tp=.1 -modeGDM=-1 -ng=4
	mv interfaceNFDTDtwoSquaresInterfaceGDMNp1Np1Order2maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoSquaresInterfaceOrder2GDMNp1Np1maxNorm.m $(interfaceTables)
	conv.p interface.conv -dm=gdm -eps1=1 -eps2=2. -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -kx=2 -ky=1 -ic=pmic -eps1=1 -eps2=2. -bc=d -cg=twoSquaresInterface -order=2 -useSosupDissipation=1  -tf=.5 -tp=.1 -modeGDM=0 -ng=4
	mv interfaceNFDTDtwoSquaresInterfaceGDMNp1Np1Mode0Order2maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoSquaresInterfaceOrder2GDMNp1Np1Mode0maxNorm.m $(interfaceTables)

interfaceGDM2dOrder4:
	conv.p interface.conv -dm=gdm -eps1=1 -eps2=2. -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -kx=2 -ky=1 -ic=pmic -eps1=1 -eps2=2. -bc=d -cg=twoSquaresInterface -order=4 -useSosupDissipation=1  -tf=.5 -tp=.1 -modeGDM=-1 -ng=4
	mv interfaceNFDTDtwoSquaresInterfaceGDMNp1Np1Order4maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoSquaresInterfaceOrder4GDMNp1Np1maxNorm.m $(interfaceTables)
	conv.p interface.conv -dm=gdm -eps1=1 -eps2=2. -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -kx=2 -ky=1 -ic=pmic -eps1=1 -eps2=2. -bc=d -cg=twoSquaresInterface -order=4 -useSosupDissipation=1  -tf=.5 -tp=.1 -modeGDM=0 -ng=4
	mv interfaceNFDTDtwoSquaresInterfaceGDMNp1Np1Mode0Order4maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoSquaresInterfaceOrder4GDMNp1Np1Mode0maxNorm.m $(interfaceTables)


interfaceGDM3d: interfaceGDM3dOrder2 interfaceGDM3dOrder4

# --- DOUBLE-BOX -- grid is slightly finer in x-direction to over-come problem with 3d order 4 interface
interfaceGDM3dOrder2:
	conv.p interface.conv -dm=gdm -npv=1 1 -kx=2 -ky=1 -kz=0 -ic=pmic -eps1=1 -eps2=2.  -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=doubleBoxInterface -tf=.5 -tp=.1 -order=2 -useSosupDissipation=1 -modeGDM=-1 -ng=4
	mv interfaceNFDTDdoubleBoxInterfaceGDMNp1Np1Order2maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDdoubleBoxInterfaceOrder2GDMNp1Np1maxNorm.m  $(interfaceTables)
	conv.p interface.conv -dm=gdm -npv=1 1 -kx=2 -ky=1 -kz=0 -ic=pmic -eps1=1 -eps2=2.  -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=doubleBoxInterface -tf=.5 -tp=.1 -order=2 -useSosupDissipation=1 -modeGDM=0 -ng=4
	mv interfaceNFDTDdoubleBoxInterfaceGDMNp1Np1Mode0Order2maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDdoubleBoxInterfaceOrder2GDMNp1Np1Mode0maxNorm.m  $(interfaceTables)

interfaceGDM3dOrder4:
	conv.p interface.conv -dm=gdm -npv=1 1 -kx=2 -ky=1 -kz=0 -ic=pmic -eps1=1 -eps2=2.  -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=doubleBoxInterface -tf=.5 -tp=.1 -order=4 -useSosupDissipation=1 -modeGDM=-1 -ng=4
	mv interfaceNFDTDdoubleBoxInterfaceGDMNp1Np1Order4maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDdoubleBoxInterfaceOrder4GDMNp1Np1maxNorm.m  $(interfaceTables)
	conv.p interface.conv -dm=gdm -npv=1 1 -kx=2 -ky=1 -kz=0 -ic=pmic -eps1=1 -eps2=2.  -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=doubleBoxInterface -tf=.5 -tp=.1 -order=4 -useSosupDissipation=1 -modeGDM=0 -ng=4
	mv interfaceNFDTDdoubleBoxInterfaceGDMNp1Np1Mode0Order4maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDdoubleBoxInterfaceOrder4GDMNp1Np1Mode0maxNorm.m  $(interfaceTables)


interfaceGDM3dOrder2a:
	conv.p interface.conv -dm=gdm -npv=1 1 -kx=2 -ky=1 -kz=0 -ic=pmic -eps1=1 -eps2=2.  -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=twoBoxesInterface -tf=.5 -tp=.1 -order=2 -useSosupDissipation=1 -modeGDM=-1 -ng=3
	mv interfaceNFDTDtwoBoxesInterfaceGDMNp1Np1Order2maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoBoxesInterfaceOrder2GDMNp1Np1maxNorm.m  $(interfaceTables)
	conv.p interface.conv -dm=gdm -npv=1 1 -kx=2 -ky=1 -kz=0 -ic=pmic -eps1=1 -eps2=2.  -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=twoBoxesInterface -tf=.5 -tp=.1 -order=2 -useSosupDissipation=1 -modeGDM=0 -ng=3
	mv interfaceNFDTDtwoBoxesInterfaceGDMNp1Np1Mode0Order2maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoBoxesInterfaceOrder2GDMNp1Np1Mode0maxNorm.m  $(interfaceTables)

interfaceGDM3dOrder4a:
	conv.p interface.conv -dm=gdm -npv=1 1 -kx=2 -ky=1 -kz=0 -ic=pmic -eps1=1 -eps2=2.  -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=twoBoxesInterface -tf=.5 -tp=.1 -order=4 -useSosupDissipation=1 -modeGDM=-1 -ng=3
	mv interfaceNFDTDtwoBoxesInterfaceGDMNp1Np1Order4maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoBoxesInterfaceOrder4GDMNp1Np1maxNorm.m  $(interfaceTables)
	conv.p interface.conv -dm=gdm -npv=1 1 -kx=2 -ky=1 -kz=0 -ic=pmic -eps1=1 -eps2=2.  -npv=1 1 -a01=1 -a11=.1 -b01=1 -b11=.2 -a02=1 -a12=.2 -b02=1 -b12=.1 -interfaceNormal=1 0 0 -interfacePoint=0 .0 0 -bc=d -cg=twoBoxesInterface -tf=.5 -tp=.1 -order=4 -useSosupDissipation=1 -modeGDM=0 -ng=3
	mv interfaceNFDTDtwoBoxesInterfaceGDMNp1Np1Mode0Order4maxNormConvTable.tex $(interfaceTables)
	mv interfaceNFDTDtwoBoxesInterfaceOrder4GDMNp1Np1Mode0maxNorm.m  $(interfaceTables)
