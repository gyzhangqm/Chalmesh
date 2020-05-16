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

#include "Overture.h"
#include "MappingInformation.h"
#include "PlotStuff.h"

// Here are some user defined mappings
#ifdef ADD_USER_MAPPINGS
#include "AirfoilMapping.h"
int addToMappingList(Mapping & map);
#endif

int 
ogen(MappingInformation & mappingInfo, PlotStuff & ps, const aString & commandFileName );

int 
main(int argc, char *argv[])
{
  int plotOption=TRUE;
  String commandFileName="";
  if( argc > 1 )
  { // look at arguments for "noplot" or some other name
    String line;
    for( int i=1; i<argc; i++ )
    {
      line=argv[i];
      if( line=="noplot" )
        plotOption=FALSE;
      else if( commandFileName=="" )
        commandFileName=line;    
    }
  }
  else
    cout << "Usage: `ogen [noplot][file.cmd]' \n"
            "          noplot:   run without graphics \n" 
            "          file.cmd: read this command file \n";

  // --- create user defined mappings ----
  MappingInformation mappingInfo;
#ifdef ADD_USER_MAPPINGS
  AirfoilMapping airfoil;
  mappingInfo.mappingList.addElement(airfoil);
  // Do this so we can read the airfoil mapping from a data-base file
  addToMappingList(airfoil);
#endif
  

  // Graphics interface:
  //  PlotStuff ps(plotOption,"ogen: Overlapping Grid Generator"); 

  PlotStuff *ps_ptr;
  ps_ptr = new PlotStuff(plotOption,"ogen: Overlapping Grid Generator"); 

  // By default start saving the command file called "ogen.cmd"
  String logFile="ogen.cmd";
  //  ps.saveCommandFile(logFile);
  ps_ptr->saveCommandFile(logFile);
  cout << "User commands are being saved in the file `" << (const char *)logFile << "'\n";

  // create more mappings and/or make an overlapping grid
  //  ogen( mappingInfo,ps,commandFileName);
  ogen( mappingInfo,*ps_ptr,commandFileName);

  delete ps_ptr;
  return 0;
}
