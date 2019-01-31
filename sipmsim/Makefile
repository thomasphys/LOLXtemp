ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs) -lm -ldl
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -lm -ldl

#ARCH          = -m32
OPT           = -O2 -march=amdfam10 -mfpmath=sse -lgomp
CXX           = /usr/bin/g++
CXXFLAGS      = -Wall $(OPT) $(ARCH) -fPIC
INCS          =  -I$(ROOTSYS)/include  
LD            = /usr/bin/g++
LDFLAGS       = $(OPT)
SOFLAGS       = -Wl,-soname -shared 
LIB_EXTRA     = 
LIBS          = $(ROOTGLIBS) $(LIB_EXTRA) -L/usr/X11R6/lib -lXpm -lX11

exe: crosstalksim

%.o:	%.cxx
	$(CXX) -c $(INCS) $(CXXFLAGS) $<

%.o:	%.C
	$(CXX) -c $(INCS) $(CXXFLAGS) $<

crosstalksim:  crosstalksim.o
	$(LD) $(LDFLAGS) -o crosstalksim crosstalksim.o $(LIBS)
    
clean:
	rm -f *.o
	rm -f *~

distclean: clean
	rm -f crosstalksim

