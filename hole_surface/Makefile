include ../makefile.defs

SRCS= hole_surface.c newton_search.c plot_sgrid.c administer_hs.c bi_linear.c

OBJS= $(SRCS:.c=.o) 

hole_surface.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@	
	echo "Don't forget to remake the executable!"

clean:
	rm -f *.o *.a
