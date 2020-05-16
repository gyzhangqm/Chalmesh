include ../makefile.defs

SRCS= chalmesh.c welcome.c run_tutorial.c save_overlapping_grid.c volume_face.c

OBJS= $(SRCS:.c=.o)

chalmesh = $(OBJS)

../bin/chalmesh.$(SYSTEM): $(chalmesh) $(CLIBS)
	if [ ! -d ../bin ]; then mkdir ../bin; fi
	$(CC) -o $@ $(CFLAGS) $(chalmesh) $(CLIBS) $(GLLIBS) $(XLIBS) $(HDFLIBS) \
        $(EXTRALIBS) -lm

../bin/chalmesh_mesa.$(SYSTEM): $(chalmesh) $(CLIBS)
	if [ ! -d ../bin ]; then mkdir ../bin; fi
	$(CC) -o $@ $(CFLAGS) $(chalmesh) $(CLIBS) -L$(MESA_HOME)/lib $(GLLIBS_MESA) \
              $(XLIBS) $(HDFLIBS) $(EXTRALIBS) -lm

clean:
	rm -f *.o *.a
