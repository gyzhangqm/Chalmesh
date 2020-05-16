include ../makefile.defs

SRCS= administer_vol.c cartesian.c normal_surface.c plot_vol.c \
      set_surf_label.c volume_bc.c volume_point.c cylindrical.c \
      volume_plot_mode.c hyperbolic.c normal_spheroid.c

OBJS= $(SRCS:.c=.o)

volume.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@	
	echo "Don't forget to remake the executable!"

clean:
	rm -f *.o *.a

