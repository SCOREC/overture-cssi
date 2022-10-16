common = ../../common

# The convergence rate script is called conv.p
all:  
	ln -sf $(common)/src/conv.p .

cgsm:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) cgsm


ogen = $(Overture)/bin/ogen

tablesISM = ~/Dropbox/AMP/ism/tables
matlabISM = ~/Dropbox/AMP/ism/matlab

annulusILE: annulusILEpc annulusILEpcCopy annulusILEme annulusILEmeCopy

tf = 1.0
numRes = 4

# Annulus eigenmode, ILE, ts=pc2
annulusILEpcOld:
	conv.p incompressible.conv -known=annulus -icase=1 -n=1 -m=1 -tf=$(tf) -tp=$(tf) -cdv=.25 -ts=pc2 -useCurlCurl=1 -ng=$(numRes)
	conv.p incompressible.conv -known=annulus -icase=2 -n=1 -m=1 -tf=$(tf) -tp=$(tf) -cdv=.25 -ts=pc2 -useCurlCurl=1 -ng=$(numRes)
	conv.p incompressible.conv -known=annulus -icase=3 -n=1 -m=1 -tf=$(tf) -tp=$(tf) -cdv=.25 -ts=pc2 -useCurlCurl=1 -ng=$(numRes)

# copy files to ism paper folder
annulusILEpcCopyOld:
	cp annulusIncompressibleElasticityTSpc2Case1BCddModen1m1.ConvTable.tex	$(tablesISM)
	cp annulusIncompressibleElasticityTSpc2Case1BCddModen1m1.m $(matlabISM)
	cp annulusIncompressibleElasticityTSpc2Case2BCtdModen1m1.ConvTable.tex	$(tablesISM)
	cp annulusIncompressibleElasticityTSpc2Case2BCtdModen1m1.m $(matlabISM)	
	cp annulusIncompressibleElasticityTSpc2Case3BCttModen1m1.ConvTable.tex	$(tablesISM)
	cp annulusIncompressibleElasticityTSpc2Case3BCttModen1m1.m $(matlabISM)		

# New modes: 
annulusILEpc:
	conv.p incompressible.conv -known=annulus -icase=1 -n=1 -m=1 -tf=$(tf) -tp=$(tf) -cdv=.25 -ts=pc2 -useCurlCurl=1 -ng=$(numRes)
	conv.p incompressible.conv -known=annulus -icase=2 -n=2 -m=4 -tf=$(tf) -tp=$(tf) -cdv=.25 -ts=pc2 -useCurlCurl=1 -ng=$(numRes)
	conv.p incompressible.conv -known=annulus -icase=3 -n=3 -m=2 -tf=$(tf) -tp=$(tf) -cdv=.25 -ts=pc2 -useCurlCurl=1 -ng=$(numRes)

# copy files to ism paper folder
annulusILEpcCopy:
	cp annulusIncompressibleElasticityTSpc2Case1BCddModen1m1.ConvTable.tex	$(tablesISM)
	cp annulusIncompressibleElasticityTSpc2Case1BCddModen1m1.m $(matlabISM)
	cp annulusIncompressibleElasticityTSpc2Case2BCtdModen2m4.ConvTable.tex	$(tablesISM)
	cp annulusIncompressibleElasticityTSpc2Case2BCtdModen2m4.m $(matlabISM)	
	cp annulusIncompressibleElasticityTSpc2Case3BCttModen3m2.ConvTable.tex	$(tablesISM)
	cp annulusIncompressibleElasticityTSpc2Case3BCttModen3m2.m $(matlabISM)			

# Annulus eigenmode, ILE, ts=me
annulusILEmeOld:
	conv.p incompressible.conv -known=annulus -icase=1 -n=1 -m=1 -tf=$(tf) -tp=$(tf) -cdv=.25 -ts=me -useCurlCurl=1 -ng=$(numRes)
	conv.p incompressible.conv -known=annulus -icase=2 -n=1 -m=1 -tf=$(tf) -tp=$(tf) -cdv=.25 -ts=me -useCurlCurl=1 -ng=$(numRes)
	conv.p incompressible.conv -known=annulus -icase=3 -n=1 -m=1 -tf=$(tf) -tp=$(tf) -cdv=.25 -ts=me -useCurlCurl=1 -ng=$(numRes)	

# copy files to ism paper folder
annulusILEmeCopyOld:
	cp annulusIncompressibleElasticityTSmeCase1BCddModen1m1.ConvTable.tex	$(tablesISM)
	cp annulusIncompressibleElasticityTSmeCase1BCddModen1m1.m $(matlabISM)
	cp annulusIncompressibleElasticityTSmeCase2BCtdModen1m1.ConvTable.tex	$(tablesISM)
	cp annulusIncompressibleElasticityTSmeCase2BCtdModen1m1.m $(matlabISM)	
	cp annulusIncompressibleElasticityTSmeCase3BCttModen1m1.ConvTable.tex	$(tablesISM)
	cp annulusIncompressibleElasticityTSmeCase3BCttModen1m1.m $(matlabISM)		


# new modes:
annulusILEme:
	conv.p incompressible.conv -known=annulus -icase=1 -n=1 -m=1 -tf=$(tf) -tp=$(tf) -cdv=.25 -ts=me -useCurlCurl=1 -ng=$(numRes)
	conv.p incompressible.conv -known=annulus -icase=2 -n=2 -m=4 -tf=$(tf) -tp=$(tf) -cdv=.25 -ts=me -useCurlCurl=1 -ng=$(numRes)
	conv.p incompressible.conv -known=annulus -icase=3 -n=3 -m=2 -tf=$(tf) -tp=$(tf) -cdv=.25 -ts=me -useCurlCurl=1 -ng=$(numRes)	
# copy files to ism paper folder
annulusILEmeCopy:
	cp annulusIncompressibleElasticityTSmeCase1BCddModen1m1.ConvTable.tex	$(tablesISM)
	cp annulusIncompressibleElasticityTSmeCase1BCddModen1m1.m $(matlabISM)
	cp annulusIncompressibleElasticityTSmeCase2BCtdModen2m4.ConvTable.tex	$(tablesISM)
	cp annulusIncompressibleElasticityTSmeCase2BCtdModen2m4.m $(matlabISM)	
	cp annulusIncompressibleElasticityTSmeCase3BCttModen3m2.ConvTable.tex	$(tablesISM)
	cp annulusIncompressibleElasticityTSmeCase3BCttModen3m2.m $(matlabISM)		



tfDisk = 2
numResDisk = 4
diskILE: diskILEpc diskILEpcCopy diskILEme diskILEmeCopy

# Disk eigenmode, ILE, ts=me
diskILEmeOld:
	conv.p incompressible.conv -known=disk -icase=1 -n=1 -m=1 -tf=$(tfDisk) -tp=$(tfDisk) -cdv=.25 -ts=me -useCurlCurl=1 -ng=$(numResDisk)
	conv.p incompressible.conv -known=disk -icase=2 -n=1 -m=1 -tf=$(tfDisk) -tp=$(tfDisk) -cdv=.25 -ts=me -useCurlCurl=1 -ng=$(numResDisk)

# copy files to ism paper folder
diskILEmeCopyOld:
	cp diskIncompressibleElasticityTSmeCase1BCdModen1m1.ConvTable.tex	$(tablesISM)
	cp diskIncompressibleElasticityTSmeCase1BCdModen1m1.m $(matlabISM)
	cp diskIncompressibleElasticityTSmeCase2BCtModen1m1.ConvTable.tex	$(tablesISM)
	cp diskIncompressibleElasticityTSmeCase2BCtModen1m1.m $(matlabISM)	

# Disk eigenmode, ILE, ts=me *NEW MODES*
diskILEme:
	conv.p incompressible.conv -known=disk -icase=1 -n=1 -m=2 -tf=$(tfDisk) -tp=$(tfDisk) -cdv=.25 -ts=me -useCurlCurl=1 -ng=$(numResDisk)
	conv.p incompressible.conv -known=disk -icase=2 -n=2 -m=3 -tf=$(tfDisk) -tp=$(tfDisk) -cdv=.25 -ts=me -useCurlCurl=1 -ng=$(numResDisk)

# copy files to ism paper folder
diskILEmeCopy:
	cp diskIncompressibleElasticityTSmeCase1BCdModen1m2.ConvTable.tex	$(tablesISM)
	cp diskIncompressibleElasticityTSmeCase1BCdModen1m2.m $(matlabISM)
	cp diskIncompressibleElasticityTSmeCase2BCtModen2m3.ConvTable.tex	$(tablesISM)
	cp diskIncompressibleElasticityTSmeCase2BCtModen2m3.m $(matlabISM)	


# Disk eigenmode, ILE, ts=pc2
diskILEpcOld:
	conv.p incompressible.conv -known=disk -icase=1 -n=1 -m=1 -tf=$(tfDisk) -tp=$(tfDisk) -cdv=.25 -ts=pc2 -useCurlCurl=1 -ng=$(numResDisk)
	conv.p incompressible.conv -known=disk -icase=2 -n=1 -m=1 -tf=$(tfDisk) -tp=$(tfDisk) -cdv=.25 -ts=pc2 -useCurlCurl=1 -ng=$(numResDisk)

# copy files to ism paper folder 
diskILEpcCopyOld:
	cp diskIncompressibleElasticityTSpc2Case1BCdModen1m1.ConvTable.tex	$(tablesISM)
	cp diskIncompressibleElasticityTSpc2Case1BCdModen1m1.m $(matlabISM)
	cp diskIncompressibleElasticityTSpc2Case2BCtModen1m1.ConvTable.tex	$(tablesISM)
	cp diskIncompressibleElasticityTSpc2Case2BCtModen1m1.m $(matlabISM)	


# Disk eigenmode, ILE, ts=pc2 *NEW MODES*
diskILEpc:
	conv.p incompressible.conv -known=disk -icase=1 -n=1 -m=2 -tf=$(tfDisk) -tp=$(tfDisk) -cdv=.25 -ts=pc2 -useCurlCurl=1 -ng=$(numResDisk)
	conv.p incompressible.conv -known=disk -icase=2 -n=2 -m=3 -tf=$(tfDisk) -tp=$(tfDisk) -cdv=.25 -ts=pc2 -useCurlCurl=1 -ng=$(numResDisk)

# copy files to ism paper folder
diskILEpcCopy:
	cp diskIncompressibleElasticityTSpc2Case1BCdModen1m2.ConvTable.tex	$(tablesISM)
	cp diskIncompressibleElasticityTSpc2Case1BCdModen1m2.m $(matlabISM)
	cp diskIncompressibleElasticityTSpc2Case2BCtModen2m3.ConvTable.tex	$(tablesISM)
	cp diskIncompressibleElasticityTSpc2Case2BCtModen2m3.m $(matlabISM)	


# ----------------------- PERIODIC STRIP ILE --------------------------
tfStrip=1.5
numResStrip = 5
stripILE: stripILEme stripILEmeCopy stripILEpc stripILEpcCopy

# Strip eigenmode, ILE, ts=me
stripILEme:
	conv.p incompressible.conv -known=strip -icase=1 -tf=1.5 -tp=$(tfStrip) -cdv=.25 -ts=me -ng=$(numResStrip)
	conv.p incompressible.conv -known=strip -icase=2 -tf=1.5 -tp=$(tfStrip) -cdv=.25 -ts=me -ng=$(numResStrip)

# copy files to ism paper folder
stripILEmeCopy:
	cp stripIncompressibleElasticityTSmeCase1BCtdModen1m1.ConvTable.tex	$(tablesISM)
	cp stripIncompressibleElasticityTSmeCase1BCtdModen1m1.m $(matlabISM)
	cp stripIncompressibleElasticityTSmeCase2BCttModen1m1.ConvTable.tex	$(tablesISM)
	cp stripIncompressibleElasticityTSmeCase2BCttModen1m1.m $(matlabISM)	

# Strip eigenmode, ILE, ts=pc2
stripILEpc:
	conv.p incompressible.conv -known=strip -icase=1 -tf=1.5 -tp=$(tfStrip) -cdv=.25 -ts=pc2 -ng=$(numResStrip)
	conv.p incompressible.conv -known=strip -icase=2 -tf=1.5 -tp=$(tfStrip) -cdv=.25 -ts=pc2 -ng=$(numResStrip)

# copy files to ism paper folder
stripILEpcCopy:
	cp stripIncompressibleElasticityTSpc2Case1BCtdModen1m1.ConvTable.tex	$(tablesISM)
	cp stripIncompressibleElasticityTSpc2Case1BCtdModen1m1.m $(matlabISM)
	cp stripIncompressibleElasticityTSpc2Case2BCttModen1m1.ConvTable.tex	$(tablesISM)
	cp stripIncompressibleElasticityTSpc2Case2BCttModen1m1.m $(matlabISM)	


# ----------------------- HOLLOW CYLINDER --------------------------
cylGrids:
	$(ogen) -noplot -noplot cylinder -periodic=p -za=0 -zb=1 -factor=1
	$(ogen) -noplot -noplot cylinder -periodic=p -za=0 -zb=1 -factor=2
	$(ogen) -noplot -noplot cylinder -periodic=p -za=0 -zb=1 -factor=4
	$(ogen) -noplot -noplot cylinder -periodic=p -za=0 -zb=1 -factor=8

hollowCyl: 
	conv.p incompressible.conv -known=hollowCyl -caseName=hollowCylinderDD -tf=.5 -tp=.5 -cdv=.25 -ts=me -ng=3






