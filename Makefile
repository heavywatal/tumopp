## Directories and Files
PACKAGE := $(notdir ${CURDIR})
SRCDIR := src
OBJDIR := build
DESTDIR := ${HOME}/local
INCLUDEDIR := -I${HOME}/local/include
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
GXX := $(notdir $(firstword $(foreach x,g++-6 g++-5 g++,$(shell which $x 2>/dev/null))))
ifeq ($(shell uname -s),Darwin)
		CXX := clang++
  SHLIBFLAGS = -dynamiclib -install_name ${DESTDIR}/lib/$@ -Wl,-rpath,${BOOST}/lib
else
		CXX := ${GXX}
  SHLIBFLAGS = -shared -Wl,-soname,${DESTDIR}/lib/$@ -Wl,-rpath,${DESTDIR}/lib -Wl,-rpath,${BOOST}/lib
endif
export CXX
CC := $(CXX)
CPPFLAGS := -Wall -Wextra -Wno-unused-parameter -fno-strict-aliasing ${INCLUDEDIR} ${CPPDBG} -ftemplate-depth=512
CXXFLAGS := -std=c++14 -O2 -fPIC ${CXXDBG}
LDFLAGS = -L${DESTDIR}/lib -L${BOOST}/lib
LDLIBS := -lsfmt -lboost_program_options-mt -lboost_filesystem-mt -lboost_system-mt -lboost_iostreams-mt #-lboost_zlib
TARGET_ARCH := -m64 -msse -msse2 -msse3
ARFLAGS := -rcs

ifneq (,$(filter $(CXX), ${GXX}))
  CXXFLAGS += -mfpmath=sse
		BOOST := ${HOME}/local/boost-gcc
else
  CXXFLAGS += -stdlib=libc++
		BOOST := ${HOME}/local/boost-clang
  ifeq ($(shell uname -s), Linux)
    LDLIBS += -lsupc++
  endif
endif

## Targets
.PHONY: all archive shared clean run test help
.DEFAULT_GOAL := all

all:
	$(MAKE) -j3 ${PROGRAM}

${PROGRAM}: main.cpp ${SHARED_OBJ}
	$(MAKE) install
	$(LINK.cpp) ${OUTPUT_OPTION} $< -l${PACKAGE}

archive: ${ARCHIVE}
	@:

shared: ${SHARED_OBJ}
	@:

${ARCHIVE}: ${OBJS}
	$(AR) ${ARFLAGS} $@ $^
	ranlib $@

${SHARED_OBJ}: ${OBJS}
	$(LINK.cpp) ${OUTPUT_OPTION} ${SHLIBFLAGS} $^ ${LDLIBS}

install: ${SHARED_OBJ}
	$(INSTALL) -m 644 -t ${DESTDIR}/include/${PACKAGE} ${HEADERS}
	$(INSTALL) -t ${DESTDIR}/lib $^

clean:
	$(RM) ${OBJS} ${SHARED_OBJ} ${ARCHIVE} ${PROGRAM}

run: ${PROGRAM}
	@./$<

test: ${PROGRAM}
	./$< --test

help: ${PROGRAM}
	./$< --help


.PHONY: debug release profile
debug:
	$(MAKE) CXXDBG="-g" all

release:
	$(MAKE) CPPDBG="-DNDEBUG" all

profile: ${PROGRAM}
	instruments -t "Time Profiler" -D profile$$(date +%Y%m%d) ${PROGRAM}

${OBJDIR}/%.o: ${SRCDIR}/%.cpp | ${OBJDIR}
	$(COMPILE.cpp) ${OUTPUT_OPTION} $<

${OBJDIR}:
	mkdir $@

## Dependencies
-include Dependfile

.PHONY: Depend
Depend:
	$(CXX) -MM ${CPPFLAGS} ${CXXFLAGS} ${TARGET_ARCH} ${SRCS} | sed 's|\(\w*\.o:\)|${OBJDIR}/\1|' > Dependfile
