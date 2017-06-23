# Makefile for Stan.
##


# The default target of this Makefile is...
help:

## Disable implicit rules.
SUFIXES:

##
# Users should only need to set these three variables for use.
# - CC: The compiler to use. Expecting g++ or clang++.
# - O: Optimization level. Valid values are {0, 1, 2, 3}.
# - AR: archiver (must specify for cross-compiling)
# - OS_TYPE: {mac, win, linux}
# - C++11: Compile with C++11 extensions, Valid values: {true, false}. 
##
CC = clang++
O = 0
O_STANC = 0
AR = ar
C++11 = true

##
# Set default compiler options.
## 
CFLAGS = -I src -isystem $(EIGEN) -isystem $(BOOST) -isystem $(CVODES)/include -isystem $(MATH) -Wall -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -DFUSION_MAX_VECTOR_SIZE=12 -DNO_FPRINTF_OUTPUT -pipe -Wno-unused-local-typedefs -Wno-c++11-extensions
CFLAGS_GTEST = -DGTEST_USE_OWN_TR1_TUPLE
LDLIBS = 
LDLIBS_STANC = -Lbin -lstanc
EXE = 
WINE =

-include $(HOME)/.config/stan/make.local  # define local variables
-include make/local                       # overwrite local variables


##
# Library locations
##
STAN ?= stan
MATH ?= $(STAN)/lib/stan_math/
-include $(MATH)make/libraries

##
# Get information about the compiler used.
# - CC_TYPE: {g++, clang++, mingw32-g++, other}
# - CC_MAJOR: major version of CC
# - CC_MINOR: minor version of CC
##
-include make/detect_cc

# OS_TYPE is set automatically by this script
##
# These includes should update the following variables
# based on the OS:
#   - CFLAGS
#   - CFLAGS_GTEST
#   - EXE
##
-include make/detect_os

model1: src/model1.cpp src/util.hpp src/gmo.hpp
	$(CC) $(CFLAGS) -O$(O) -o model1 src/model1.cpp -I src/

model2: src/model2.cpp src/util.hpp src/gmo.hpp
	$(CC) $(CFLAGS) -O$(O) -o model2 src/model2.cpp -I src/

test: src/test.cpp src/util.hpp src/gmo.hpp
	$(CC) $(CFLAGS) -O$(O) -o test src/test.cpp -I src/
	./test

%.hpp: %.stan
	~/scm/cmdstan/bin/stanc $< --o=$@
