#include "ogen_interface.h"

int 
main(int argc, char *argv[])
{
  void *mappingInfo_ptr;
 
  mappingInfo_ptr = initializeOvertureMappings();

  makeMappingsAndPlotThem(mappingInfo_ptr);

  cleanupOvertureMappings(mappingInfo_ptr);

  return 0;
}
