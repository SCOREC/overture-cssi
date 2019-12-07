#ifndef _ReferenceCounting
#define _ReferenceCounting

//
// Who to blame:  Geoff Chesshire
//

#include "OvertureDefine.h"
#ifndef OV_USE_OLD_STL_HEADERS
#include <iostream>
#else
#include <iostream.h>
#endif
#include "aString.H"
#include "OvertureTypes.h"

class ReferenceCounting {
  public:
    inline ReferenceCounting() {
        className = "ReferenceCounting";
        globalID = incrementGlobalIDCounter();
        referenceCount = uncountedReferences = 0;
    }
    inline ReferenceCounting
      (const ReferenceCounting& x, const CopyType ct = DEEP) {
        className = "ReferenceCounting";
        globalID = incrementGlobalIDCounter();
        referenceCount = uncountedReferences = 0;
    }
    inline virtual ~ReferenceCounting() {
        if (referenceCount) {
            cerr << "ReferenceCounting::~ReferenceCounting():  "
                 << "referenceCount != 0 for globalID = "
                 << getGlobalID() << "." << endl;
            abort();
        } // end if
    }
    inline virtual ReferenceCounting& operator=(const ReferenceCounting& x)
     { return *this; }
    inline virtual void reference(const ReferenceCounting& x) {}
    inline virtual void breakReference() { }
    inline virtual ReferenceCounting* virtualConstructor
      (const CopyType ct = DEEP) const
      { return new ReferenceCounting(*this, ct); }
    inline Integer incrementReferenceCount()
      { return ++referenceCount + uncountedReferences; }
    inline Integer decrementReferenceCount()
      { return --referenceCount + uncountedReferences; }
    inline Integer getReferenceCount() const
      { return referenceCount + uncountedReferences; }
    inline Logical uncountedReferencesMayExist() {
        if (!referenceCount) uncountedReferences = 1;
        return uncountedReferences;
    }
    inline virtual aString getClassName() const { return className; }
    Integer getGlobalID() const { return globalID; }
    // change the global ID -- not a wise thing to normally do -- used by plotting routines.
    void setGlobalID(int id) { globalID=id; } 
    virtual void consistencyCheck() const;
  private:
    Integer referenceCount;
    Integer uncountedReferences;
    Integer globalID;
    aString className;
    static Integer globalIDCounter;
    static Integer incrementGlobalIDCounter();
};
ostream& operator<<(ostream& s, const ReferenceCounting& x);

#endif // _ReferenceCounting
