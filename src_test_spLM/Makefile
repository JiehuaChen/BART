#
# compiler
#
CC := gcc 
CXX := c++
AR := ar

### Customising
#
# Adjust the following if necessary: 
# EXECUTABLE is the target executable's filename; 
# LIBDIR is a list of library path of project; 
# INCDIR is a list of include path of project; 
# LIBS is a list of libraries to link in (e.g. alleg, stdcx, iostr, stlport_gcc etc);
# COMOPT is a list of compile options.
# 

EXECUTABLE := testsplm
LIBDIR := /usr/lib/R/lib
INCDIR := /usr/share/R/include
#INCDIR := nc
LIBS := blas lapack R
COMOPT := -D_LINUX -D_DEBUG -MD -g

# Now alter any implicit rules' variables if you like, e.g.:

# compiler flags -- command-line options are: 
# make mode=release (default)
# make mode=development
# make mode=debug
ifeq ($(mode),release)
	CFLAGS := -O3 -Wall
else ifeq ($(mode),development)
	CFLAGS := -O2 -g -p -Wall
else 
	CFLAGS := -O0 -g -p -Wall
endif

CXXFLAGS := $(CFLAGS)
POSTBUILD := 

# You shouldn't need to change anything below this point.

RM-F := rm -f

SOURCE := $(wildcard *.c) $(wildcard *.cc) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cc,%.o,$(patsubst %.cpp,%.o,$(SOURCE))))
DEPS := $(patsubst %.o,%.d,$(OBJS))
MISSING_DEPS := $(filter-out $(wildcard $(DEPS)),$(DEPS))
MISSING_DEPS_SOURCES := $(wildcard $(patsubst %.d,%.c,$(MISSING_DEPS)) \
$(patsubst %.d,%.cc,$(MISSING_DEPS)) \
$(patsubst %.d,%.cpp,$(MISSING_DEPS)))
CPPFLAGS += $(COMOPT) $(addprefix -I, $(INCDIR))

.PHONY : everything deps objs clean veryclean rebuild

everything : $(EXECUTABLE)

deps : $(DEPS)

objs : $(OBJS)

clean :
	@$(RM-F)	*.o
	@$(RM-F)	*.d

veryclean: clean
	@$(RM-F) $(EXECUTABLE)

rebuild: veryclean everything

ifneq ($(MISSING_DEPS),)
$(MISSING_DEPS) :
	@$(RM-F) $(patsubst %.d,%.o,$@)
endif

-include $(DEPS)


$(EXECUTABLE) : $(OBJS)
	$(CXX) -o $(EXECUTABLE) $(CXXFLAGS) $(OBJS) $(addprefix -L, $(LIBDIR)) $(addprefix -l, $(LIBS))
	$(POSTBUILD)
