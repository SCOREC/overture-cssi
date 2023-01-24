#
# Top level Makefile for cg solvers
#
# NOTE: 
# setenv CGBUILDPREFIX /some/other/dir will put the obj/lib/bin in /some/other/dir/{common,cgad...}/{lib,obj,dep,bin}
#

all : libCommon cgad cgcns cgins cgasf cgsm cgmp unitTests cgmx message

everything: all 

message : cgad cgcns cgins cgasf cgsm cgmp unitTests cgmx
	@echo "CG:INFO: type \"make check-quiet\" or \"make -j8 check-quiet\" or \"make check\" to run regression tests"	

libCommon:
	@cd common; $(MAKE) 

# Astronomical plasma calculations
cgcssi: libCommon
	@cd cssi;    $(MAKE) 

# Unit tests
unitTests: libCommon
	@cd unitTests;     $(MAKE) 

# run regression tests:
check:
	-@cd unitTests/check;  $(MAKE)  check
	-@cd cssi/check; $(MAKE)  check

check-unitTests-quiet:
	-@cd unitTests/check;  $(MAKE) -s  check-quiet
check-cssi-quiet:
	-@cd cssi/check;  $(MAKE) -s check-quiet

# run regression tests with minimal output:  [ make -j8 check-quiet ]
check-quiet: check-unitTests-quiet \
             check-cssi-quiet

# clean out obj and bin directories and temporary files.
clean:
	rm -f *~ */*~ */*/*~ */*/*/*~
	-@cd cssi;            $(MAKE) clean
	-@cd common;        $(MAKE) clean
	-@cd common/moving; $(MAKE) clean
	-@cd user;          $(MAKE) clean
	-@cd unitTests;     $(MAKE) clean

