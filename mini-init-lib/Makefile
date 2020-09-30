CFLAGS=-I. -I../maplecomsps_lrb_vsids_18/MapleCOMSPS_LRB_VSIDS_2/

CSRCS=init.cpp
COBJS=$(CSRCS:.cpp=.o)


%.o: %.cpp
	g++ $(CFLAGS) -c -o $@ $<

libmininit.a: $(COBJS)
	ar -rcsv $@ $^


all: libmininit.a

clean:
	rm -f libmininit.a *.o
