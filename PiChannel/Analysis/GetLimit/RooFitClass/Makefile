GCC = g++
FLAGS = -O3 -D DEFECT_NO_EXPLICIT_QUALIFICATION
INCLUDES =  $(ROOTSYS)/include/
SPECIALFLAGS=-O -fPIC
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLDFLAGS=$(shell root-config --ldflags)
ROOTLIBS=$(shell root-config --libs)
CFLAGS = $(SPECIALFLAGS) -I.
LFLAGS = -O2
RCXX=$(CFLAGS) $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS) -l pthread

MYLIBNAME=RooTauLeptonInvisible
SRC1=TPseudoRestFrame
SRC2=$(MYLIBNAME)
OBJECTS=$(SRC1).o $(SRC2).o
DICNAME=MyDict

all: $(MYLIBNAME).so

$(MYLIBNAME).so: $(DICNAME).cxx $(OBJECTS)
	g++ -shared  $(CFLAGS) -o  $(MYLIBNAME).so $(ROOTLDFLAGS) $(CXXFLAGS) -I $(RCXX) $(DICNAME).cxx $(OBJECTS)
$(SRC2).o: $(SRC2).cxx  
	$(GCC) $(FLAGS) $(CFLAGS) -o $(SRC2).o -c $(SRC2).cxx -I $(RCXX)
$(SRC1).o:
	$(GCC) $(FLAGS) $(CFLAGS) -o $(SRC1).o -c $(SRC1).cxx -I $(RCXX)
$(DICNAME).cxx: $(SRC1).h $(SRC2).h $(MYLIBNAME)_Linkdef.h
	rootcling -f $@ -p $^ #-c $(ROOTCFLAGS) -p $^
clean:
	rm -f $(OBJECTS) $(MYLIBNAME).so $(DICNAME).cxx 




