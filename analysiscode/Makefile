ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs) -lm -ldl
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -lm -ldl

#ARCH          = -m32
OPT           = -O2 -march=amdfam10 -mfpmath=sse
CXXFLAGS      = -Wall $(OPT) $(ARCH) -fPIC -std=c++11
INCS          =  -I$(ROOTSYS)/include  -I../ds -I../include
LD            =/usr/bin/g++
LDFLAGS       = $(OPT)
SOFLAGS       = -Wl,-soname -shared 
LIB_EXTRA     = 
LIBS          = $(ROOTGLIBS) $(LIB_EXTRA) -L/usr/X11R6/lib -lXpm -lX11
LOCALLIB = /home/tmcelroy/lolxsim/ds/libEvent.so
#LOCALLIB = -lEvent

exe:   ProcessRun EventDisplay EventDisplay_Cyl LXeElectronCherenkov

%.o:	%.cxx
	$(LD) -c $(INCS) $(CXXFLAGS) $<

%.o:	%.C
	$(LD) -c $(INCS) $(CXXFLAGS) $<

LOLXReadData.o: ../src/LOLXReadData.cc
	$(LD) -c $(INCS) $(CXXFLAGS) $<

EventDisplay_Cyl: EventDisplay_Cyl.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) $(LOCALLIB)

EventDisplay: EventDisplay.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) $(LOCALLIB)

ProcessRun: ProcessRun.o LOLXReadData.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) $(LOCALLIB) $(INCS)

LXeElectronCherenkov: LXeElectronCherenkov.o LOLXReadData.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) $(LOCALLIB) $(INCS)
