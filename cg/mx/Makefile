#
# Main Makefile for cgmx -- Maxwell's Equations
#
THISPKG := mx

# check for targets that do not need source file dependencies:
NODEP := $(strip $(filter usage doc clean cleano cleanc cleandep,$(MAKECMDGOALS)))

# ==== optional capabilities: use on/off ===

# set to "on" to build 6th and 8th order (but general BC's not implemented yet)
buildHigherOrder := off
# buildHigherOrder := on

# ---- sparse solver packages ----
# These are sometimes needed by CgMx: e.g. for projections, for implicit interpolation in parallel
# To compile with PETSc you should set: 
# setenv usePETSc  
ifeq ($(usePETSc),)
  usePETSc   := off
else
  usePETSc   := on
endif
# To compile with SuperLU you should set: 
# setenv useSuperLU 
ifeq ($(useSuperLU),)
  useSuperLU := off
else
  useSuperLU := on
endif


include configFile
include $(CONFIGDIR)/make.dirs
include $(CONFIGDIR)/make.$(OS)
include $(CONFIGDIR)/make.rules
# Over-ride default compile flags here: (or setenv COMPILE [dbg/opt/optdbg] from the shell)
#    compile:=[dbg/opt/optdbg]
#    profile:=[on/off]   
# compile:=dbg
include $(CONFIGDIR)/make.options
include make

# gfortran has trouble with interface3dOrder4.bf and -O3 ?
FC_FLAGS_OPT = -O

# memory checking for gcc or clang: 
# CXX_FLAGS += -fsanitize=address

# List of packages to build
PKGS := src

# include Overture/Oges to get the PETSc solver and Ogmg files
vpath % $(PKGS) $(GUIDIR)/src $(Overture)/Oges $(OvertureCheckout)/oges
VPATH = $(PKGS) $(GUIDIR)/src $(Overture)/Oges $(OvertureCheckout)/oges

SRC_FILES := 
TESTS := 
PKG_INCLUDE := $(foreach p, $(PKGS), $(p)/files)

all: $(BINDIR)cgmx

cgmx: $(BINDIR)cgmx

# Test: read dispersion material file NOTE TYPE "make tr" to build this 
trdmf: $(BINDIR)trdmf

usage: 
	-@echo "usage : "
	-@echo "        % make " $(COPTS) "<all|cgmx>"

config_os:
	-@echo "using configuration for $(OS)"

include $(PKG_INCLUDE)

# This next include defines how to compile the lists: 
#  SRC_FILES : list of C++, C, and fortran files that are to be compiled in the default way.
#  CSRCO : list of C++ files that are to be compiled optimized by default
#  FSCRO : list of Fortran files that are to be compiled optimized by default
#  FSRC   : list of Fortran files that are to be compiled without auto-dbl
#  FSCRSO : list of Fortran files that are to be compiled optimized and without auto-dbl
include $(CONFIGDIR)/make.src

# add dependencies for main programs that are not in the previous source file lists
DEP += cgmx.d

# objco:; -@echo "OBJCO=$(OBJCO)"

.PHONY: doc cgmx

# Cgmx needs the data base files (but does not link to libCommon) so make a link to them here:
# src/dbAccess.C : ../../common/dataBase/dbAccess.C
# 	ln -sf ../../common/dataBase/dbAccess.C src/dbAccess.C
# 	ln -sf ../../common/dataBase/dbFunctions.f src/dbFunctions.f
# 	ln -sf ../../common/dataBase/kk_ptr.cc src/kk_ptr.cc

libCommon:=$(CGBUILDPREFIX)common/libCgCommon

# make cgmx in the bin dir:
PETSC_FILES_OBJ := $(addprefix $(OBJPREFIX),$(PETSC_FILES))
CGMX_EXTRA_FILES = $(PETSC_FILES_OBJ)

# compile PETSc files optimized
$(PETSC_FILES_OBJ) : $(OBJPREFIX)%.o : %.C
	$(CXX) $(CXX_FLAGSO) -o $@ -c $<
	

fftmpi:
	-@echo "USE_FFTMPI = $(USE_FFTMPI)"

# For parallel FFTs you should "setenv USE_FFTMPI" 
ifneq ($(USE_FFTMPI),)
  # Note: if you change the version of fftmpi, recompile FourierTransform.o 

  # For now this is where fftmpi is found:
  ifeq ( $(USE_FFTMPI),KISS )
    # Version of fftMPI using KISS FFT: 
    FFTDIR = /home/henshw/software/fftmpi-1Oct18
  else 
    # Version of fftMPI using FFTW:
    FFTDIR = /home/henshw/software/fftmpiFFTW
    FFT1D = -L/home/henshw/software/fftw-3.3.8/lib -lfftw3
  endif 

  FFT_INC =       -DUSE_FFTMPI -I$(FFTDIR)/src
  FFT_PATH =      -L$(FFTDIR)/src
  FFT2D_LIB =   -Wl,-rpath,$(FFTDIR)/src $(FFT_PATH) -lfft2dmpi $(FFT1D)
  FFT3D_LIB =   -Wl,-rpath,$(FFTDIR)/src $(FFT_PATH) -lfft3dmpi $(FFT1D)
endif 

$(OBJDIR)FourierTransform.o: ./src/FourierTransform.C; $(CXX) $(CXX_FLAGSO) $(FFT_INC) -o $(OBJDIR)$*.o -c $<


# $(BINDIR)cgmx: $(OBJDIR)cgmx.o $(OBJALL) $(CGMX_EXTRA_FILES) src/dbAccess.C
$(BINDIR)cgmx: $(OBJDIR)cgmx.o $(OBJALL) $(CGMX_EXTRA_FILES) $(libCommon)
	-@mkdir -p $(BINDIR)
	$(CXX) $(CXX_FLAGS) -o $@ $(OBJDIR)cgmx.o $(OBJALL) $(CGMX_EXTRA_FILES)  $(FFT2D_LIB) $(FFT3D_LIB) $(LIBS)

#	$(CXX) $(CXX_FLAGS) -o $@ $(OBJDIR)cgmx.o $(CGMX) $(OBJALL) $(CGMX_EXTRA_FILES) $(LIBS)


# test program to read a Dispersive Material File
TESTS += $(BINDIR)trdmf
tr: $(OBJDIR)trdmf.o $(OBJDIR)DispersiveMaterialParameters.o $(OBJDIR)BAMaterialFunctions.o $(OBJDIR)BAPlaneInterface.o \
        $(OBJDIR)PlaneWaveExactSolution.o
	$(CXX) $(CXX_FLAGS) -o $(BINDIR)trdmf $(OBJDIR)trdmf.o $(OBJALL) $(LIBS)

# Test program for exact solutions
TESTS += $(BINDIR)texact
texact: $(OBJDIR)texact.o $(OBJDIR)SphereExactSolutions.o $(OBJDIR)scatSphere.o $(OBJDIR)PlaneInterfaceExactSolution.o $(OBJDIR)evalDispersionRelation.o $(OBJDIR)DispersiveMaterialParameters.o $(OBJDIR)BAMaterialFunctions.o $(OBJDIR)BAPlaneInterface.o $(OBJDIR)SlabsExactSolution.o $(OBJDIR)ComplexArray.o
	$(CXX) $(CXX_FLAGS) -o $(BINDIR)texact $(OBJDIR)texact.o $(OBJALL) $(LIBS)

# Test program for non-local radiation BC
TESTS += $(BINDIR)trad
# trad: $(OBJDIR)trad.o $(OBJDIR)RadiationBoundaryCondition.o $(OBJDIR)RadiationKernel.o $(OBJDIR)gaussianIntegral.o $(OBJDIR)FourierTransform1.o $(OBJDIR)FourierTransform.o 
trad: $(OBJDIR)trad.o $(OBJDIR)RadiationBoundaryCondition.o $(OBJDIR)RadiationKernel.o $(OBJDIR)gaussianIntegral.o $(OBJDIR)FourierTransform.o 
	$(CXX) $(CXX_FLAGS) -o $(BINDIR)trad $(OBJDIR)trad.o $(OBJALL) $(LIBS) $(FFT2D_LIB) $(FFT3D_LIB)



# tb: $(OBJDIR)tbm.o $(OBJDIR)DispersiveMaterialParameters.o
#	$(CXX) $(CXX_FLAGS) -o $(BINDIR)tbm $(OBJDIR)tbm.o $(OBJALL) $(LIBS)

# build the library in cg/common
$(libCommon):; @cd ../common; $(MAKE) lib

$(DEPDIR)depdep: $(DEPDIR) $(BINDIR) $(LIBDIR) $(OBJDIR)
	touch $(DEPDIR)depdep

$(BINDIR) $(LIBDIR) $(DEPDIR) $(OBJDIR):
	@-mkdir -p $@

doc:
	-@mkdir doc

check: 

# remove all c and c++ .o files, dependencies and executables (but leave fortran .o files)
cleanc: 
	-rm -f $(DEPPREFIX)*.d $(COBJ) $(OBJCO) $(BINDIR)cgmx $(LIBDIR)*.a

cleano:
	-rm -f $(OBJPREFIX)*.o
cleandep:
	-rm -f $(OBJPREFIX)*.d
clean:
	-rm -f *~ */*~ $(OBJPREFIX)*.o $(DEPPREFIX)*.d $(BINDIR)cgmx  \
         cmd/*.show* cmd/*.debug cmd/*.check cmd/*.log cmd/*.move cmd/cgmx.cmd 

# include dependencies except for a list of commands 
ifeq ($(NODEP),)
-include $(DEPDIR)depdep
-include $(DEP)
endif

