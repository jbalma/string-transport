CXX=CC
CFLAGS=-O3 -g -hgnu
LDFLAGS= -lstdc++

INCLUDE= -D__USE_OPENMP -D__USE_MPI

SRCS_CXX=cs.cpp Distribution.cpp Vex.cpp Points.cpp
OBJS_CXX=$(subst .cpp,.o,$(SRCS_CXX))


all: MAIN

MAIN: $(SRCS_CXX)
	$(CXX) $(CFLAGS) $(INCLUDE) $(SRCS_CXX)  -o CS.x

LINK: $(OBJS_CXX)
	$(LD) $(LD_FLAGS) $(INCLUDE) $(OBJS_CXX) -o CS.x

depend: .depend

.depend: $(OBJS_CXX)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;


MISC_REM = CS.x *.ptx *.cub

clean:
	$(RM) $(OBJS_CXX) $(MISC_REM)

dist-clean: clean
	$(RM) *~ .dependtool

