ROOTLIBS = $(shell $(ROOTSYS)/bin/root-config --libs) -lm -ldl
ROOTGLIBS = $(shell $(ROOTSYS)/bin/root-config --glibs) -lm -ldl

INC =  -I$(ROOTSYS)/include
LIBS = $(ROOTGLIBS)
GCC = g++ -O2
exe = libEvent.so
HEADERS = Event.h  MPPC.h

libEvent.so: Dict.cxx
	g++ -shared -o$@ `root-config --ldflags` -fPIC -std=c++11 $(INC) $(LIBS) $^

Dict.cxx: $(HEADERS) LinkDef.h 
	rootcint -f $@ -c  -p $(INC) $^
