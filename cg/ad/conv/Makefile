common = ../../common

ogen = $(Overture)/bin/ogen
cgad = $(CGBUILDPREFIX)/ad/bin/cgad
convDir = $(CG)/ad/conv

# The convergence rate script is called conv.p
all:  
	ln -sf $(common)/src/conv.p .

cgad:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) cgad


# --- square eigenfunction ---
squareEigPC2FD2:
	conv.p exact.conv -known=squareEigen -cg=square -ts=pc2 -order=2 -tf=.5 -tp=.5 -kappa=.1 -kx=2 -ky=2 -ng=4

squareEigPC2FD4:	
	conv.p exact.conv -known=squareEigen -cg=square -ts=pc2 -order=4 -tf=.5 -tp=.5 -kappa=.1 -kx=2 -ky=2 -ng=4

squareEigBDF4FD4:	
	conv.p exact.conv -known=squareEigen -cg=square -ts=bdf -bdfOrder=4 -order=4 -tf=.5 -tp=.5 -kappa=.1 -kx=2 -ky=2 -dtMax=.05 -ng=5	

# --- disk eigen function ----
diskEigCNFD2:
		conv.p exact.conv -known=diskEigen -ts=im  -bdfOrder=2 -order=2 -cg=sic -n=1 -m=1 -tf=.5 -tp=.5 -kappa=.1 -bc=d -dtMax=.05 -ng=4

diskEigBDF2FD2:
		conv.p exact.conv -known=diskEigen -ts=bdf -bdfOrder=2 -order=2 -cg=sic -n=1 -m=1 -tf=.5 -tp=.5 -kappa=.1 -bc=d -dtMax=.05 -ng=4

diskEigBDF4FD4:
		conv.p exact.conv -known=diskEigen -ts=bdf -bdfOrder=4 -order=4 -cg=sic -n=1 -m=1 -tf=.5 -tp=.5 -kappa=.1 -bc=d -dtMax=.05 -ng=4

# --- annulus eigenfunction ---
annulusEigCNFD2:
	conv.p  exact.conv -known=annulusEigen -cg=annulus -ts=im -order=2 -n=0 -m=1 -tf=.5 -tp=.5 -kappa=.1 -bc=d -dtMax=.05 -ng=4

annulusEigBDF4FD4:
	conv.p  exact.conv -known=annulusEigen -cg=annulus -bdfOrder=4 -order=4 -n=0 -m=1 -tf=.5 -tp=.5 -kappa=.1 -bc=d -dtMax=.05 -ng=4	
