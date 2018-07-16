ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs) -lm -ldl
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -lm -ldl

#ARCH          = -m32
OPT           = -O2 -march=amdfam10 -mfpmath=sse -lgomp
CXX           = /usr/bin/g++
CXXFLAGS      = -Wall $(OPT) $(ARCH) -fPIC
INCS          =  -I$(ROOTSYS)/include  -I../ds
LD            = /usr/bin/g++
LDFLAGS       = $(OPT)
SOFLAGS       = -Wl,-soname -shared 
LIB_EXTRA     = 
LIBS          = $(ROOTGLIBS) $(LIB_EXTRA) -L/usr/X11R6/lib -lXpm -lX11
LOCALLIB = ../ds/libEvent.so

exe:   ProcessRun EventDisplay EventDisplay_Cyl

%.o:	%.cxx
	$(CXX) -c $(INCS) $(CXXFLAGS) $<

%.o:	%.C
	$(CXX) -c $(INCS) $(CXXFLAGS) $<

EventDisplay_Cyl: EventDisplay_Cyl.o
	$(LD) $(LDFLAGS) -o EventDisplay_Cyl EventDisplay_Cyl.o $(LIBS) $(LOCALLIB)

EventDisplay: EventDisplay.o
	$(LD) $(LDFLAGS) -o EventDisplay EventDisplay.o $(LIBS) $(LOCALLIB)

ProcessRun: ProcessRun.o
	$(LD) $(LDFLAGS) -o ProcessRun ProcessRun.o $(LIBS) $(LOCALLIB)
