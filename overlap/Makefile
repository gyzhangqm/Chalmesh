include ../makefile.defs

SRCS= administer_overlap.c over3d.c oct_tree.c \
      overlap_3d_parameters.c tri_linear.c plot_overlap.c plot_mode.c 

OBJS= $(SRCS:.c=.o)

overlap.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@	
	echo "Don't forget to remake the executable!"

clean:
	rm -f *.o *.a
