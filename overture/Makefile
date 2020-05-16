include ../makefile.defs

SRCS= ogen_interface.C overmap.C mogl.C

OBJS= $(SRCS:.C=.o) 

ogen_interface.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@	
	echo "Don't forget to remake the executable!"

ogen: calling_ogen.o $(OBJS) 
	$(CCxx) $(CCFLAGS) -o $@ $(ogen) $(CLIBS2) $(GLIBS) 
  
