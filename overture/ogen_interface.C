//-----------------------------------------------------------------------------------------------
//  ogen: Overlapping Grid Generation
//
// This program can be used to create an overlapping grid. 
// The first step in making an overlapping grid is to define some
// Mappings that define the component grids. The second step is 
// to call the grid generator which will determine how the component
// grids overlap. The third step is to save the overlapping grid 
// in a file.
//
//-----------------------------------------------------------------------------------------------

//=========================================================================================
// Here is the driver program for `ogen' - the overlapping grid generator
//
//   Usage: type 
//       ogen
//  to run with graphics, or type
//       ogen noplot
//   to run without graphics, or
//       ogen file.cmd
//   to run ogen with graphics and read in a command file, or
//       ogen noplot file.cmd
//   to run ogen without graphics and read in a command file.
//
//  By default user commands will be saved in the file "ogen.cmd"
//
//  You can add to the driver any nonstandard Mapping's that you want to use.
//  See the example below where (if the macro ADD_USERMAPPINGS is defined) an AirfoilMapping
//  is created and added to a list. The list is then passed to ogen. The Mapping
//  can be subsequently changed within ogen, if required.
//
//  Thus, for example, your compile line should look something like:
//      CC -DADD_USERMAPPINGS .... ogenDriver.C 
//
//===========================================================================================

#include "GL_GraphicsInterface.h"
#include "MappingInformation.h"
#include "HDF_DataBase.h"
#include "Overture.h"
#include "PlotStuff.h"

int 
createMappings( MappingInformation & mapInfo );

extern "C" {
  static PlotStuff *ps_ptr=NULL;

void *
initializeOvertureMappings(void){
  MappingInformation *mappingInfo_ptr;

  ios::sync_with_stdio();     // Synchronize C++ and C I/O subsystems
  Index::setBoundsCheck(on);  // Turn on A++ array bounds checking
    
  // allocate storage for user defined mappings
  mappingInfo_ptr = new MappingInformation;

  return mappingInfo_ptr;
}

void
cleanupOvertureMappings(void *overture_data){
  MappingInformation *mappingInfo_ptr = (MappingInformation *) overture_data;
// delete the mappings
  delete mappingInfo_ptr;
  // delete the plotting context
  delete ps_ptr;
}

void
makeMappingsAndPlotThem(void *overture_data){
  MappingInformation *mappingInfo_ptr = (MappingInformation *) overture_data;
  int plotOption=TRUE;

  // create a graphics interface:
  if (ps_ptr == NULL){
    ps_ptr = new PlotStuff(plotOption,"ogen: Overlapping Grid Generator");
    mappingInfo_ptr->graphXInterface=ps_ptr;
    ps_ptr->defaultPrompt="chalmesh-overture-interface>";
  }

  // create mappings
  createMappings(*mappingInfo_ptr);

}
}
