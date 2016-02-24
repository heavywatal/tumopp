## Directories and Files
PACKAGE := $(notdir ${CURDIR})
SRCDIR := .
OBJDIR := build
INCLUDEDIR := -isystem ${HOME}/local/include
PROGRAM := a.out


## Directories and Files
SRCS := $(wildcard ${SRCDIR}/*.cpp)
OBJS := $(notdir $(SRCS:.cpp=.o))
OBJS := $(addprefix ${OBJDIR}/,${OBJS})
vpath %.cpp ${SRCDIR}


## Options
GXX := $(notdir $(firstword $(foreach x,g++-5 g++,$(shell which $x))))
CXX_ARRAY := clang++ ${GXX}
CXX := $(firstword $(foreach x,${CXX_ARRAY},$(shell which $x)))
CC := $(CXX)
CPPFLAGS := -Wall -Wextra -Wno-unused-parameter -fno-strict-aliasing ${INCLUDEDIR} ${CPPDBG} -ftemplate-depth=512
CXXFLAGS := -std=c++14 -O3 ${CXXDBG}
LDFLAGS := -L${HOME}/local/lib
LDLIBS := -lsfmt -lboost_program_options -lboost_filesystem -lboost_system -lboost_iostreams -lz
TARGET_ARCH := -m64 -msse -msse2 -msse3

ifneq (,$(filter $(CXX), ${GXX}))
  CXXFLAGS += -mfpmath=sse
  CPPFLAGS += -pthread -D_GLIBCXX_USE_NANOSLEEP
  LDFLAGS += -L${HOME}/local/boost-gcc/lib -static-libstdc++
else
  CXXFLAGS += -stdlib=libc++
  LDFLAGS += -L${HOME}/local/boost-clang/lib
  ifeq ($(shell uname), Linux)
    LDLIBS += -lpthread -lsupc++
  endif
endif

export CXX CC TARGET_ARCH

## Targets
.PHONY: all clean run test help
.DEFAULT_GOAL := all

all:
	${MAKE} -j3 ${PROGRAM}

${PROGRAM}: ${OBJS}
	${LINK.cpp} ${OUTPUT_OPTION} $^ ${LOADLIBES} ${LDLIBS}

clean:
	${RM} ${OBJS} ${PROGRAM}

run: ${PROGRAM}
	@./$<

test: ${PROGRAM}
	./$< --test

help: ${PROGRAM}
	./$< --help


.PHONY: debug release instruments
debug:
	${MAKE} CXXDBG="-g" all

release:
	${MAKE} CPPDBG="-DNDEBUG" all

instruments: ${PROGRAM}
	instruments -t "Time Profiler" -D profile$$(date +%Y%m%d) ${PROGRAM}

${OBJDIR}/%.o: | ${OBJDIR}
	$(COMPILE.cpp) ${OUTPUT_OPTION} $<

${OBJDIR}:
	mkdir $@

## Dependencies
-include Dependfile

.PHONY: Depend
Depend:
	${CXX} -MM ${CPPFLAGS} ${CXXFLAGS} ${TARGET_ARCH} ${SRCS} | sed 's|\(\w*\.o:\)|${OBJDIR}/\1|' > Dependfile
