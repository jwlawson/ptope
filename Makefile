NAME = ptope
MAJOR = 0
MINOR = 6
VERSION = $(MAJOR).$(MINOR)

ifeq ($(CXX),g++)
#One of these causes a problem with overflowing hashes
CXXFLAGS += -Wall -Wextra -march=native\
	-fno-signed-zeros\
	-fno-math-errno\
	-fno-rounding-math\
	-fno-signaling-nans -fno-trapping-math -ffinite-math-only\
	-Wno-misleading-indentation
OPT += -O3 -g
AR = gcc-ar
else
CXXFLAGS += -Wall -xHOST
OPT += -O3 -ipo
AR = xiar
endif
CXXFLAGS += -DARMA_DONT_USE_WRAPPER -DARMA_NO_DEBUG -DNDEBUG
B_OPT += $(OPT)

uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
uname_O := $(shell sh -c 'uname -o 2>/dev/null || echo not')

# Using cygwin -std=gnu++11 should be used rather than -std=c++11
ifeq ($(uname_O),Cygwin)
	CXXFLAGS += -std=gnu++11 -DCYGWIN_STOI
endif
ifeq ($(uname_S),Linux)
	CXXFLAGS += -std=c++11
endif

LIB = lib$(NAME).so.$(VERSION)
STATIC = lib$(NAME).a
TEST = test$(NAME)

LDFLAGS = -shared -Wl,-soname,$(LIB)

# Specify base directory
BASE_DIR = .

# Specify source directory
SRC_DIR = $(BASE_DIR)/src

# Specify test directory
TEST_DIR = $(BASE_DIR)/test

# define the output directory for .o
OBJ_DIR = $(BASE_DIR)/build

INC_DIR = $(BASE_DIR)/include

# Install directory
PREFIX = $(HOME)

# define any directories containing header files other than /usr/include
INCLUDES = -I$(HOME)/include -I$(BASE_DIR)/include \
           -I$(BASE_DIR)/lib/include

# define library paths in addition to /usr/lib
LFLAGS = -L$(HOME)/lib -L$(BASE_DIR)/lib

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) use -lx
LIBS =
TEST_LIBS = -lgtest -lgtest_main -lopenblas -llapack -pthread \
						-lboost_system

# define the C source files
SRCS = $(filter-out $(SRC_DIR)/bench.cc,$(wildcard $(SRC_DIR)/*.cc))
TEST_SRCS = $(wildcard $(TEST_DIR)/*.cc)

# define the C object files
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .cc of all words in the macro SRCS
# with the .o suffix
#
_OBJS = $(SRCS:.cc=.o)
_TEST_OBJS = $(TEST_SRCS:.cc=.o)

# Puts objs in obj_dir
OBJS = $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(_OBJS))
TEST_OBJS = $(patsubst $(TEST_DIR)/%,$(OBJ_DIR)/%,$(_TEST_OBJS))

PROF = #-fprofile-use

.PHONY: clean

all:	$(LIB)

$(TEST): CXXFLAGS += -flto -fuse-linker-plugin
$(TEST): OPT = -O3
$(TEST): $(OBJS) $(TEST_OBJS)
	$(CXX) $(PROF) $(CXXFLAGS) $(B_OPT) $(INCLUDES) -o $(TEST) $(TEST_OBJS) $(OBJS) $(LFLAGS) $(TEST_LIBS)

test: $(TEST)
	@echo Running tests
	@./$(TEST)

bench: CXXFLAGS += -flto -fuse-linker-plugin
bench: OPT = -O3
bench: $(STATIC) src/bench.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OPT) -c src/bench.cc -o build/bench.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OPT) build/bench.o -L. -Wl,-Bstatic -lptope -lbenchmark -Wl,-Bdynamic $(LFLAGS) -lopenblas -llapack -lboost_system -pthread -o bench

lib:	$(LIB)
static:	$(STATIC)

$(LIB): CXXFLAGS += -fPIC
$(LIB):	$(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

$(STATIC): CXXFLAGS += -flto -fuse-linker-plugin -ffat-lto-objects
$(STATIC): OPT = -O3
$(STATIC): $(OBJS)
	$(AR) rcs $(STATIC) $(OBJS)

install:	$(LIB)
	cp $(LIB) $(PREFIX)/lib
	ldconfig -v -n $(PREFIX)/lib
	ln -fs $(PREFIX)/lib/$(LIB) $(PREFIX)/lib/lib$(NAME).so
	mkdir -p $(PREFIX)/include/$(NAME)
	cp -ru include/* $(PREFIX)/include/$(NAME)/

uninstall:
	rm $(PREFIX)/lib/lib$(NAME).*
	rm -r $(PREFIX)/include/$(NAME)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file)
# (see the gnu make manual section about automatic variables)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(INC_DIR)/%.h
	$(CXX) $(PROF) $(CXXFLAGS) $(OPT) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cc
	$(CXX) $(PROF) $(CXXFLAGS) $(OPT) $(INCLUDES) -c $< -o $@

$(OBJS): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	$(RM) *.o *~ $(MAIN) $(OBJ_DIR)/*.o $(LIB) $(STATIC) $(TEST)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
