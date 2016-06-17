## Directories and Files
PACKAGE := $(notdir ${CURDIR})
SRCDIR := src
OBJDIR := build
DESTDIR := ${HOME}/local
INCLUDEDIR := -isystem ${HOME}/local/include
ARCHIVE := lib${PACKAGE}.a
SHARED_OBJ := lib${PACKAGE}.so
PROGRAM := a.out
INSTALL := install -p -D

## Directories and Files
SRCS := $(wildcard ${SRCDIR}/*.cpp)
HEADERS := $(wildcard ${SRCDIR}/*.hpp)
OBJS := $(notdir $(SRCS:.cpp=.o))
OBJS := $(addprefix ${OBJDIR}/,${OBJS})
vpath %.cpp ${SRCDIR}


## Options
GXX := $(notdir $(firstword $(foreach x,g++-6 g++-5 g++,$(shell which $x))))
CXX := $(notdir $(firstword $(foreach x,clang++ ${GXX},$(shell which $x 2>/dev/null))))
CC := $(CXX)
CPPFLAGS := -Wall -Wextra -Wno-unused-parameter -fno-strict-aliasing ${INCLUDEDIR} ${CPPDBG} -ftemplate-depth=512
CXXFLAGS := -std=c++14 -O2 -fPIC ${CXXDBG}
LDFLAGS := -L${HOME}/local/lib
LDLIBS := -lsfmt -lboost_program_options -lboost_filesystem -lboost_system -lboost_iostreams -lboost_zlib #-lz
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

${PROGRAM}: main.cpp
	${MAKE} install
	${LINK.cpp} ${OUTPUT_OPTION} $^ -l${PACKAGE}

${ARCHIVE}: ${OBJS}
	$(AR) -rcs $@ $^
	ranlib $@

${SHARED_OBJ}: ${OBJS}
	${LINK.cpp} ${OUTPUT_OPTION} -dynamiclib -install_name ${DESTDIR}/lib/$@ $^ ${LDLIBS}

install: ${SHARED_OBJ}
	$(INSTALL) -m 644 -t ${DESTDIR}/include/${PACKAGE} ${HEADERS}
	$(INSTALL) -t ${DESTDIR}/lib $^

clean:
	${RM} ${OBJS} ${SHARED_OBJ} ${ARCHIVE} ${PROGRAM}

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
