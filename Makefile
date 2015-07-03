## Directories and Files
PACKAGE := $(notdir ${CURDIR})
SRCDIR := .
OBJDIR := build
INCLUDEDIR := -isystem /usr/local/include -iquote ${HOME}/local/include
PROGRAM := a.out


## Directories and Files
SRCS := $(wildcard ${SRCDIR}/*.cpp)
OBJS := $(notdir $(SRCS:.cpp=.o))
OBJS := $(addprefix ${OBJDIR}/,${OBJS})
vpath %.cpp ${SRCDIR}


## Options
GXX := $(firstword $(foreach x,g++-4.9 g++-4.8 g++,$(shell which $x)))
CXX_ARRAY := clang++ ${GXX}
CXX := $(firstword $(foreach x,${CXX_ARRAY},$(shell which $x)))
CC := $(CXX)
CPPFLAGS := -Wall -Wextra -Wno-unused-parameter -fno-strict-aliasing ${INCLUDEDIR} ${CPPDBG} -ftemplate-depth=512
CXXFLAGS := -std=c++11 -O3 ${CXXDBG}
LDFLAGS := -L${HOME}/local/lib -L/usr/local/lib
LDLIBS := -lsfmt -lboost_program_options -lboost_filesystem -lboost_system -lboost_iostreams -lz
TARGET_ARCH := -march=core2 -m64 -msse -msse2 -msse3

ifneq (,$(filter $(CXX), ${GXX}))
  CXXFLAGS += -mfpmath=sse
  CPPFLAGS += -pthread -isystem /usr/local/boost-gcc/include -D_GLIBCXX_USE_NANOSLEEP
  LDFLAGS += -L/usr/local/boost-gcc/lib -static-libstdc++
else
  CPPFLAGS += -isystem /usr/local/boost-clang/include
  CXXFLAGS += -stdlib=libc++
  LDFLAGS += -L/usr/local/boost-clang/lib
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

run:
	@./${PROGRAM}

test:
	./${PROGRAM} --test

help:
	./${PROGRAM} --help


.PHONY: debug release instruments
debug:
	${MAKE} CXXDBG="-g" all

release:
	${MAKE} CPPDBG="-DNDEBUG" all

instruments: release
	@echo $(shell gdate +%F_%T)
	instruments -t "/Applications/Xcode.app/Contents/Applications/Instruments.app/Contents/Resources/templates/Time Profiler.tracetemplate" -D ~/tmp/profile$(shell gdate +%F_%T) ${PROGRAM} -N 50000

${OBJDIR}/%.o: | ${OBJDIR}
	$(COMPILE.cpp) ${OUTPUT_OPTION} $<

${OBJDIR}:
	mkdir $@

## Dependencies
-include Dependfile

.PHONY: Depend
Depend:
	${CXX} -MM ${CPPFLAGS} ${CXXFLAGS} ${TARGET_ARCH} ${SRCS} | sed 's|\(\w*\.o:\)|${OBJDIR}/\1|' > Dependfile

## misc.
.PHONY: open

mainsrcs := $(addprefix ${SRCDIR}/,\
gland.h \
gland.cpp \
tissue.h \
tissue.cpp \
simulation.h \
simulation.cpp \
main.cpp \
)

open:
	open -a mi ${mainsrcs}
