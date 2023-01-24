#include "OvertureTypes.h"
#include <list>
#include "CssiParameters.h"

#define KK_DEBUG
#include "DBase.hh"

using namespace std;
using namespace DBase;

extern "C"
{

#define getIntFromDataBase EXTERN_C_NAME(getintfromdatabase)
int getIntFromDataBase( DataBase *pdb, char *name_, int & value, int & nameLength);

#define getIntFromDataBaseCgcssi EXTERN_C_NAME(getintfromdatabasecgcssi)
int getIntFromDataBaseCgcssi( DataBase *pdb, char *name_, int & value, int & nameLength)
// =======================================================================
/// \brief Use this routine from fortran to look-up a int or enum variables in the Cgcssi data base
///
/// \return values:
///                1=found, 
///                0=not found,
///               -1=name found but not the correct type
///
/// \note: value is left unchanged if it was not found. 
// =======================================================================
{

  DataBase & dbase = *pdb;

  string name(name_,0,nameLength);
  // remove trailing blanks
  int i= name.find_last_not_of(" "); // position of last non-blank character
  name.erase(i+1,name.size()-i);


  // -- look for Cgcssi enums --
  if( name=="pde" && dbase.has_key("pde") )
  {
    value = (int) dbase.get<CssiParameters::PDE>("pde");
    // printF("getIntFromDataBaseCgcssi: pde found!, value=%i\n",value);
    return 1;
  }
  else if( name=="pdeVariation" && dbase.has_key("pdeVariation") )
  {
    value = (int) dbase.get<CssiParameters::PDEVariation>("pdeVariation");
    // printF("getIntFromDataBaseCgcssi: pdeVariation found!, value=%i\n",value);
    return 1;
  }
  else if( name=="conservativeGodunovMethod" && dbase.has_key("conservativeGodunovMethod") )
  {
    value = (int)dbase.get<CssiParameters::GodunovVariation>("conservativeGodunovMethod");
    return 1;
  }
  else if( name=="equationOfState" && dbase.has_key("equationOfState") )
  {
    value = (int)dbase.get<CssiParameters::EquationOfStateEnum>("equationOfState");
    return 1;
  }
  else 
  {
    return getIntFromDataBase( pdb, name_, value,nameLength);
  }
  
  return 0; // not found 
}
}
