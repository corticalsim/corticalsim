# Standard makefile v1.2 by Simon Tindemans
# 02-05-2008
#
# This makefile is built to allow for easy extension and addition of make
# targets and keeps track of the last used configuration. Calling make with
# default parameters will use the previous configuration.
#
# There are three hard-coded build targets:
# - all	[default]	: compile the program with the previously used configuration
# - clean			: removes all known object files
# - cleanmore		: removes all configuration subdirectories and build flags
#
# Additional targets can be defined in CONFIG_LIST below. A combination of targets
# can be used as well, such as 'make clean release' to force a full rebuild.
#
# NOTES:
# - The makefile is designed to be used with GNU make. It is *not* compatible
# with the nmake program that comes with Microsoft Visual C++. VC++ users should 
# install GNU make (freeware).
# - Under windows, subdirectories for the make configurations (./release/, etc.) 
# are not created by the makefile. They may be made by the compiler. If not, this
# should be done manually.
# 


################## User: Program info ########################

# name of the executable (on Windows, .exe will be added by the script)
PROG_NAME = corticalSim.$(CONFIG)

# Compiler type. Current options are: gnu / pg / ms / intel
COMPILER=gnu
# Platform; should be 'win' for Windows or anything else for unix/linux, etc.
PLATFORM=anything_else

# list all necessary components. (header files, c++ files (.cpp compulsory) and object files
# NO directory names allowed
H_DEPS=corticalSim.h DLList.h eig3.h
C_DEPS=corticalSim.cpp geometry-generic.cpp geometry-special.cpp microtubule.cpp system.cpp queue.cpp parameters.cpp IO-analysis.cpp eig3.cpp
L_DEPS=

# header directory and library directory
HDIR=.
LDIR=.

# list of include directories, separated by white space
IDIRS=. $(BOOST_ROOT) ./libs


################# User: build configuration ##############

# list of possible configurations and default target for 'all'
CONFIG_LIST = release debug warning
DEFAULT_CONFIG = release

# name to use for flag files (FLAG_NAME.CONFIG)
FLAG_NAME = makeConfigurationFlag

############## User: compiler names and options for each build configuration ##########
# NOTE: the CC field can be left empty and absorbed in the _CFLAGS fields
gnu_CC=g++
#gnu_release_CFLAGS=-O3 -std=gnu++0x
#gnu_debug_CFLAGS=-g -std=gnu++0x
#gnu_warning_CFLAGS=-g -Wall -std=gnu++0x
gnu_release_CFLAGS=-O3 -std=gnu++11
gnu_debug_CFLAGS=-g -std=gnu++11
gnu_warning_CFLAGS=-g -Wall -std=gnu++11
ms_CC=cl
ms_release_CFLAGS=/Ox /EHsc /arch:SSE2 /nologo
ms_debug_CFLAGS=/EHsc /nologo /Zi
ms_warning_CFLAGS=/Wall /nologo /Zi /EHsc
intel_CC=icpc
intel_release_CFLAGS=-xW -O3 -ipo -no-prec-div -Wall -static-libcxa
# NOTE: for release, -static should be faster than -static-libcxa, but this causes issues on fedora linux...
intel_debug_CFLAGS=-g -Wall -w2
intel_warning_CFLAGS=-g -Wall -w2
pg_CC=pgCC
pg_release_CFLAGS=-tp k8-64 -fastsse -O3
pg_debug_CFLAGS=-g
pg_warning_CFLAGS=-g -Wall






################### no user changes below here #######################

#default build target is 'all'
all:

# based on the configurations, construct all possible flags and objects
ifeq ($(COMPILER),ms)
O_DEPS = $(patsubst %.cpp,%.obj,$(C_DEPS))
else
O_DEPS = $(patsubst %.cpp,%.o,$(C_DEPS))
endif
POS_FLAGS = $(addprefix $(FLAG_NAME).,$(CONFIG_LIST))
POS_OBJ = $(foreach cfg,$(CONFIG_LIST),$(patsubst %,$(cfg)/%,$(O_DEPS)))


##### construct build configuration overrides #####

# define generic rule for build configurations
define buildconfig
.PHONY: $(1)
$(1):
# set proper flags and recursively call make again
ifeq ($(PLATFORM), win)
	@copy /y makefile $(FLAG_NAME).$(1)
	@del /F $(filter-out $(FLAG_NAME).$(1),$(POS_FLAGS))
else
	@touch $(FLAG_NAME).$(1)
	@rm -f $(filter-out $(FLAG_NAME).$(1),$(POS_FLAGS))
endif
	@$(MAKE)
endef
# construct rules for all build configurations
$(foreach cfg,$(CONFIG_LIST),$(eval $(call buildconfig,$(cfg))))


##### Look for set build configuration flags ######

# look for configuration flag files of type 'FLAG_NAME.*' and extract type
FLAG_FILES := $(wildcard $(FLAG_NAME).*)
FLAG_TYPE := $(subst $(FLAG_NAME).,,$(FLAG_FILES))

ifeq ($(words $(FLAG_TYPE)),0)
 # no flags are found: default config
 $(info no configuration flag: using default [$(DEFAULT_CONFIG)])
 CONFIG := $(DEFAULT_CONFIG)
else 
ifeq ($(words $(FLAG_TYPE)),1)
 $(info found flag: [$(FLAG_TYPE)])
 # one flag found; check if configuration is valid
 #ifeq ($(findstring $(FLAG_TYPE),$(CONFIG_LIST)),)
 ifeq ($(words $(filter $(FLAG_TYPE),$(CONFIG_LIST))),0)
  $(error Invalid flag found: '$(FLAG_NAME).$(FLAG_TYPE)'. Please remove manually)
 endif
 #set corect default ; not 100% safe when flag_type is a subset of a proper config
 CONFIG := $(FLAG_TYPE)
else
 $(error Multiple build configurations flagged. Please remove '$(FLAG_NAME).*' manually.)
endif
endif


##### expand compiler flags and file paths
# NOTE delayed, by using '=' instead of ':='

# load compiler flags based on configuration
CFLAGS=$($(COMPILER)_$(CONFIG)_CFLAGS)
# load object directory based on configuration
ODIR=$(CONFIG)

ifeq ($(COMPILER),ms)
INCLUDELIST := $(addprefix /I,$(IDIRS))
else
INCLUDELIST := $(addprefix -I,$(IDIRS))
endif

# construct file paths
HDRLIST = $(patsubst %,$(HDIR)/%,$(H_DEPS))
OBJLIST = $(patsubst %,$(ODIR)/%,$(O_DEPS))
LIBLIST = $(patsubst %,$(LDIR)/%,$(L_DEPS))


###### create compilation and linking rules

# define generic rule for object compilation
define create_object
ifeq ($(PLATFORM),win)
$(1): $(2) $(HDRLIST)
	$($(COMPILER)_CC) $(CFLAGS) -c /Fo$(1) $(2) $(INCLUDELIST)
else
$(1): $(2) $(HDRLIST) | createdir
	$($(COMPILER)_CC) $(CFLAGS) -c -o $(1) $(2) $(INCLUDELIST)
endif
endef

# safety rule to construct configuration directories when necessary
# only called on non-windows systems
.PHONY: createdir
createdir:
	@mkdir -p $(CONFIG)

# create rules for all possible configurations
ifeq ($(COMPILER),ms)
$(foreach ofile,$(POS_OBJ), $(eval $(call create_object,$(ofile),$(patsubst %.obj,%.cpp,$(notdir $(ofile))))))
else
$(foreach ofile,$(POS_OBJ), $(eval $(call create_object,$(ofile),$(patsubst %.o,%.cpp,$(notdir $(ofile))))))
endif

# add the '.exe' extension to the executable on Windows systems
ifeq ($(PLATFORM), win)
PROG_NAME := $(strip $(PROG_NAME).exe)
endif

# rule for the main program compilation [default]
.PHONY: all
all: $(OBJLIST)
ifeq ($(COMPILER),ms)
	$($(COMPILER)_CC) $(CFLAGS) /Fe$(PROG_NAME) $^ $(LIBLIST) $(INCLUDELIST)
else
	$($(COMPILER)_CC) $(CFLAGS) -o $(PROG_NAME) $^ $(LIBLIST) -L$(LDIR) $(INCLUDELIST)
endif

# standard clean; only removes known files; leave flags alone for full rebuild
.PHONY: clean
clean:
ifeq ($(PLATFORM),win)
ifeq ($(COMPILER),ms)
	del *.obj /s /q 
else
	del *.o /s /q 
endif
else
	rm -f $(POS_OBJ)
endif

# extensive clean; removes configuration directories and flag files
.PHONY: cleanmore	
cleanmore:
ifeq ($(PLATFORM),win)
	del /Q /S $(CONFIG_LIST) $(FLAG_NAME).*
else
	rm -fR $(CONFIG_LIST) $(FLAG_NAME).*
endif
