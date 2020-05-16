include ../makefile.defs

SRCS = ogl_plot.c
OBJS = $(SRCS:.c=.o)

ogl_plot.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@	

clean:
	rm -f *.o *.a
