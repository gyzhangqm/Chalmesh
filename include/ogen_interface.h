#ifndef ogen_interface_h
#define ogen_interface_h

void *
initializeOvertureMappings(void);
void
cleanupOvertureMappings(void *overture_data);
void
makeMappingsAndPlotThem(void *overture_data);

#endif
