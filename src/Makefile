# -----------------------------------------------------------------
#   Makefile for tifosi
# --------------------------------------------------------------------

FORCE_DYNAMIC = 

# Set this variable to either UNIX, MAC or WIN
SYS = UNIX

# Put C++ compiler here. It is platform specific.
#CXX_UNIX = g++
# Older versions of the compiler:
CXX_UNIX = g++
CXX_WIN = c:\bin\mingw\bin\mingw32-g++.exe

# Any other compiler flags here
CXXFLAGS = 
LIB = 

# --------------------------------------------------------------------
# Do not edit below this line!!!
# --------------------------------------------------------------------

CXXFLAGS += -O3
# DEBUG:
#CXXFLAGS += -g
OUTPUT = tifosi

# Define some specific flags

ifeq ($(SYS),WIN)
 CXXFLAGS += -DWIN
 CRM = del
 CXX = $(CXX_WIN)
 ifndef FORCE_DYNAMIC
  CXXFLAGS += -static
 endif
 CXXFLAGS += -I $(HITPP) -L $(LITPP) -I $(HARMA)
endif

ifeq ($(SYS),UNIX)
 CXXFLAGS += -DUNIX
 CRM = rm -f
 CXX = $(CXX_UNIX)
 ifndef FORCE_DYNAMIC
  CXXFLAGS += -static
 endif
endif

ifeq ($(SYS),MAC)
 CXXFLAGS += -DUNIX
 CRM = rm -f
 CXX = $(CXX_UNIX)
endif

SRC = aresta.cpp celula.cpp especies.cpp main.cpp poblacio.cpp vertex.cpp aleatoris.cpp
HDR = aleatoris.h aresta.h auxiliar.h celula.h especies.h main.h poblacio.h vertex.h
OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT)

$(OUTPUT) :
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean:
	$(CRM) *.o *~ tifosi

