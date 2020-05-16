include ../makefile.defs

SRCS= surface.c plot.c sphere_patch.c disc_point.c save_surface.c

OBJS= $(SRCS:.c=.o)

surface.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@	
	echo "Don't forget to remake the executable!"

clean:
	rm -f *.o *.a
