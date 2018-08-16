/*  -*-Mode: c++; -*-  */
#ifndef GRID_COLLECTION_FUNCTION
#define GRID_COLLECTION_FUNCTION "GridCollectionFunction.h"

#include "MappedGridFunction.h" 
#include "GenericGridCollectionFunction.h"
#include "ListOfMappedGridFunction.h"
#include "ListOfGridCollectionFunction.h"


// Note that FABS is converted to abs in the int case

class GridCollection;           // forward declaration
class GridCollectionData;       // forward declaration
class Interpolant;              // forward declaration
class GenericGridCollectionOperators;  // forward declaration
class SparseRepForGCF;          // forward declaration
class GenericCollectionOperators;
class GridCollectionOperators;
class GridFunctionParameters;
class intGridCollectionFunction;
class DisplayParameters;

//===================================================================
//  GridCollectionFunction
//
//  Define a grid function to be used with a composite grid.
//  This class contains a list of MappedGridFunction's
//  and a gridCollection
//
//  This is a reference counted class so that there is no need
//  to keep a pointer to a grid function. Use the reference
//  member function to make one grid function reference another.
//
//  Usage:
//     GridCollection cg(...);  // here is a composite grid
//     doubleGridCollectionFunction u(cg),v;
//     u[0]=5.;                // mapped grid function for grid 0
//     Index I(0,10);
//     u[1](I,I)=3.;           // mapped grid function for grid 1
//     v.reference(u);         // v is referenced to u
//     v[1]=7.;                // changes u as well
//     v.breakReference();     // v is no longer referenced to u
//     ...
//
//==================================================================

class GridCollectionFunction : public GenericGridCollectionFunction
{
 public:

  enum
  {
    undefinedValue=MappedGridFunction::undefinedValue,
    defaultValue=MappedGridFunction::defaultValue,
    forAll=MappedGridFunction::forAll,
    forAllGrids=forAll+1,
    maximumNumberOfIndicies=MappedGridFunction::maximumNumberOfIndicies,
    maximumNumberOfComponents=MappedGridFunction::maximumNumberOfComponents
  };
  
  enum updateReturnValue  // the return value from updateToMatchGrid is a mask of the following values
  {
    updateNoChange          = MappedGridFunction::updateNoChange, // no changes made
    updateReshaped          = MappedGridFunction::updateReshaped, // grid function was reshaped
    updateResized           = MappedGridFunction::updateResized,  // grid function was resized
    updateComponentsChanged = MappedGridFunction::updateComponentsChanged  // component dimensions may 
                                                   //   have changed (but grid was not resized or reshaped)
  };

//  IntegerArray positionOfComponent;
//  IntegerArray positionOfCoordinate;
  // positions of the component Index's
  inline int positionOfComponent(int i) const {return rcData->positionOfComponent[i];} 
  // positions of the coordinate Index's
  inline int positionOfCoordinate(int i) const {return rcData->positionOfCoordinate[i];} 
                                 //   takes on values found in the faceCenteringType enumerator
  // numberOfComponentGrids : number of grids including refinements but excluding multigrid levels
  int numberOfComponentGrids() const;
  // numberOfGrids : total number of grids including refinements and multigrid levels
  int numberOfGrids() const;

  int numberOfRefinementLevels() const;
  int numberOfMultigridLevels() const;
  inline int positionOfFaceCentering() const { return rcData->positionOfFaceCentering; } 
  inline int faceCentering() const { return rcData->faceCentering; } 

  ListOfMappedGridFunction mappedGridFunctionList;

  //  Collections of grid functions having the same base grid.
  ListOfGridCollectionFunction baseGrid;       

  //  Collections of grid functions having the same multigrid level.
  ListOfGridCollectionFunction multigridLevel; 

  //  Collections of grid functions having the same refinement level.
  ListOfGridCollectionFunction refinementLevel;

  // the next two pointers could be in rcData!
  GridCollectionData *gridCollectionData;  // for normal access to the gridCollection
  GridCollection     *gridCollection;      // pointer to the GridCollection
  GenericCollectionOperators *operators;      // object to define derivatives and BC's
  SparseRepForGCF *sparse;                     // pointer to info on sparse representation for coefficients

  //-----------------------------------------------------------------------------------------
  //----------------------------Constructors-------------------------------------------------
  //-----------------------------------------------------------------------------------------
  GridCollectionFunction();

  // This constructor takes ranges, the first 3 "nullRange" values are taken to be the
  // coordinate directions in the grid function.
  GridCollectionFunction(GridCollection & grid, 
			 const Range & R0, 
			 const Range & R1=nullRange,
			 const Range & R2=nullRange,
			 const Range & R3=nullRange,
			 const Range & R4=nullRange,
			 const Range & R5=nullRange,
			 const Range & R6=nullRange,
			 const Range & R7=nullRange );

  GridCollectionFunction(GridCollectionData & gcData, 
			 const Range & R0, 
			 const Range & R1=nullRange,
			 const Range & R2=nullRange,
			 const Range & R3=nullRange,
			 const Range & R4=nullRange,
			 const Range & R5=nullRange,
			 const Range & R6=nullRange,
			 const Range & R7=nullRange );

  GridCollectionFunction(GridCollection & grid, 
			 const int   & i0, 
			 const Range & R1=nullRange,
			 const Range & R2=nullRange,
			 const Range & R3=nullRange,
			 const Range & R4=nullRange,
			 const Range & R5=nullRange,
			 const Range & R6=nullRange,
			 const Range & R7=nullRange );

  GridCollectionFunction(GridCollectionData & gcData, 
			 const int   & i0, 
			 const Range & R1=nullRange,
			 const Range & R2=nullRange,
			 const Range & R3=nullRange,
			 const Range & R4=nullRange,
			 const Range & R5=nullRange,
			 const Range & R6=nullRange,
			 const Range & R7=nullRange );

  //
  // This constructor takes a GridFunctionType
  // 
  GridCollectionFunction(GridCollection & gc,
			 const GridFunctionParameters::GridFunctionType & type, 
			 const Range & component0=nullRange,       // defaults to Range(0,0)
			 const Range & component1=nullRange,
			 const Range & component2=nullRange,
			 const Range & component3=nullRange,
			 const Range & component4=nullRange );

  GridCollectionFunction(GridCollection & gridCollection);
  GridCollectionFunction(GridCollectionData & gridCollection);

  GridCollectionFunction(const GridCollectionFunction & cgf, const CopyType copyType=DEEP);

  virtual ~GridCollectionFunction();

  
  virtual GridCollectionFunction & operator= ( const GridCollectionFunction & cgf );

  virtual MappedGridFunction & operator[]( const int grid ) const
  {
    return mappedGridFunctionList[grid];
  }

  // Return a link to some specfied components ** cannot be virtual**
  GridCollectionFunction operator()(const Range & component0,
				    const Range & component1=nullRange,
				    const Range & component2=nullRange,
				    const Range & component3=nullRange,
				    const Range & component4=nullRange );

  virtual void consistencyCheck() const;

  virtual aString getClassName() const;

  int getComponentBound( int component ) const;       // get the bound of the given component
  int getComponentBase( int component ) const;        // get the base of the given component
  int getComponentDimension( int component ) const;   // get the dimension of the given component

  int getCoordinateBound( int coordinate ) const;       // get the bound of the given coordinate
  int getCoordinateBase( int coordinate ) const;        // get the base of the given coordinate
  int getCoordinateDimension( int coordinate ) const;   // get the dimension of the given coordinate

  GridFunctionParameters::GridFunctionType 
    getGridFunctionType(const Index & component0=nullIndex,  // return the type of the grid function
			const Index & component1=nullIndex,
			const Index & component2=nullIndex,
			const Index & component3=nullIndex,
			const Index & component4=nullIndex ) const;     

  GridFunctionParameters::GridFunctionTypeWithComponents 
    getGridFunctionTypeWithComponents(
				      const Index & component0=nullIndex,  // return the type of the grid function
				      const Index & component1=nullIndex,
				      const Index & component2=nullIndex,
				      const Index & component3=nullIndex,
				      const Index & component4=nullIndex ) const;     

  int getNumberOfComponents() const;                  // number of components

  bool isNull();  // TRUE is this is a null grid function (no grid)

  virtual void setIsACoefficientMatrix(const bool trueOrFalse=true, 
				       const int stencilSize=defaultValue, 
				       const int numberOfGhostLines=1,
				       const int numberOfComponentsForCoefficients=1,
				       const int offset=0 );

  virtual void setIsACoefficientMatrix(const bool trueOrFalse, 
				       const IntegerArray & stencilSize, 
				       const IntegerArray & numberOfGhostLines,
				       const int numberOfComponentsForCoefficients=1,
				       const int offset=0 );

  bool getIsACoefficientMatrix() const;

  virtual int dataCopy( const GridCollectionFunction & gcf );  // copy the array data only

  virtual void display( const aString & label=nullString, FILE *file = NULL, const aString & format=nullString ) const;
  virtual void display( const aString & label, const aString & format  ) const;
  virtual void display( const aString & label, const DisplayParameters & displayParameters) const; 

  virtual int interpolate(Interpolant & interpolant, 
			  const Range & C0 = nullRange,
			  const Range & C1 = nullRange,
			  const Range & C2 = nullRange );
  virtual int interpolate(const Range & C0 = nullRange,
			  const Range & C1 = nullRange,
			  const Range & C2 = nullRange  );

  // Clean up a grid function, release the memory
  virtual int destroy();

  // link to a component of another grid Function
  virtual void link(const GridCollectionFunction & gcf,
	    const Range & R0,
	    const Range & R1=nullRange,
	    const Range & R2=nullRange,
	    const Range & R3=nullRange,
	    const Range & R4=nullRange );

  void getRanges(Range & R0,     // return the current values for the Ranges
		 Range & R1,
		 Range & R2,
		 Range & R3,
		 Range & R4,
		 Range & R5,
		 Range & R6,
		 Range & R7 ) const;

  virtual void link(const GridCollectionFunction & gcf, 
            const int componentToLinkTo=0,
            const int numberOfComponents=1 );
  

  virtual void reference( const GridCollectionFunction & cgf );
  virtual void breakReference();

  // set the name of the grid function or of a component
  void setName(const aString & name, 
               const int & component0=defaultValue,
               const int & component1=defaultValue,
               const int & component2=defaultValue,
               const int & component3=defaultValue,
               const int & component4=defaultValue );

  // get the name of the grid function or a component
  aString getName(const int & component0=defaultValue,
                 const int & component1=defaultValue,
                 const int & component2=defaultValue,
                 const int & component3=defaultValue,
                 const int & component4=defaultValue ) const ;

  virtual GridCollection* getGridCollection(const bool abortIfNull=TRUE) const;   // return a pointer to the grid collection

  virtual int get( const GenericDataBase & dir, const aString & name);    // get from a database file
  virtual int put( GenericDataBase & dir, const aString & name) const;    // put to a database file

  // Update edges of periodic grids
  virtual void periodicUpdate(const Range & C0=nullRange,
			      const Range & C1=nullRange,
			      const Range & C2=nullRange,
			      const Range & C3=nullRange,
			      const Range & C4=nullRange,
			      const bool & derivativePeriodic=FALSE);

  virtual updateReturnValue updateToMatchGrid();
  virtual updateReturnValue updateToMatchGrid(GridCollectionData & gridData, 
			const Range & R0, 
			const Range & R1=nullRange,
			const Range & R2=nullRange,
			const Range & R3=nullRange,
			const Range & R4=nullRange,
			const Range & R5=nullRange,
			const Range & R6=nullRange,
			const Range & R7=nullRange );
  
  virtual updateReturnValue updateToMatchGrid(GridCollection & grid, 
			const Range & R0, 
			const Range & R1=nullRange,
			const Range & R2=nullRange,
			const Range & R3=nullRange,
			const Range & R4=nullRange,
			const Range & R5=nullRange,
			const Range & R6=nullRange,
			const Range & R7=nullRange );

  // define this version to avoid overloading ambiguities
  virtual updateReturnValue updateToMatchGrid(GridCollection & grid, 
			const int  & i0, 
			const Range & R1=nullRange,
			const Range & R2=nullRange,
			const Range & R3=nullRange,
			const Range & R4=nullRange,
			const Range & R5=nullRange,
			const Range & R6=nullRange,
			const Range & R7=nullRange );

  virtual updateReturnValue updateToMatchGrid(const Range & R0, 
			const Range & R1=nullRange,
			const Range & R2=nullRange,
			const Range & R3=nullRange,
			const Range & R4=nullRange,
			const Range & R5=nullRange,
			const Range & R6=nullRange,
			const Range & R7=nullRange );
  
  virtual updateReturnValue updateToMatchGrid(GridCollection & grid, 
				      const GridFunctionParameters::GridFunctionType & type, 
				      const Range & component0,
				      const Range & component1=nullRange,
				      const Range & component2=nullRange,
				      const Range & component3=nullRange,
				      const Range & component4=nullRange );

  virtual updateReturnValue updateToMatchGrid(const GridFunctionParameters::GridFunctionType & type, 
				      const Range & component0,
				      const Range & component1=nullRange,
				      const Range & component2=nullRange,
				      const Range & component3=nullRange,
				      const Range & component4=nullRange );

  virtual updateReturnValue updateToMatchGrid(GridCollection & grid, 
				      const GridFunctionParameters::GridFunctionType & type);
  virtual updateReturnValue updateToMatchGrid(const GridFunctionParameters::GridFunctionType & type);


  virtual updateReturnValue updateToMatchGrid(GridCollection & gridCollection);
  virtual updateReturnValue updateToMatchGrid(GridCollectionData & gridCollection);


  // update this grid function to match another grid function
  virtual updateReturnValue updateToMatchGridFunction(const GridCollectionFunction & cgf, 
				const Range & R0, 
				const Range & R1=nullRange,
				const Range & R2=nullRange,
				const Range & R3=nullRange,
				const Range & R4=nullRange,
				const Range & R5=nullRange,
				const Range & R6=nullRange,
				const Range & R7=nullRange );

  virtual updateReturnValue updateToMatchGridFunction(const GridCollectionFunction & gcf);

  // make sure the number of mappedGridFunction's is correct
  virtual updateReturnValue updateToMatchNumberOfGrids(GridCollection& gc);
  // Update the grid collection to match the component grids
  virtual updateReturnValue updateToMatchComponentGrids();
  

  void setInterpolant(Interpolant *interpolant );
  // inquire cell centredness:
  bool getIsCellCentered(const Index & axis=nullIndex, 
                         const Index & component0=nullIndex, 
                         const Index & component1=nullIndex, 
                         const Index & component2=nullIndex, 
                         const Index & component3=nullIndex, 
                         const Index & component4=nullIndex, 
                         const Index & grid=nullIndex ) const;
  // change cell centredness:
  void setIsCellCentered(const bool trueOrFalse, 
                         const Index & axis=nullIndex, 
                         const Index & component0=nullIndex, 
                         const Index & component1=nullIndex, 
                         const Index & component2=nullIndex, 
                         const Index & component3=nullIndex, 
                         const Index & component4=nullIndex, 
                         const Index & grid=nullIndex ); 
  // set a component to be face centred along a given axis
  void setIsFaceCentered(const int   & axis=forAll,
                         const Index & component0=nullIndex, 
                         const Index & component1=nullIndex, 
                         const Index & component2=nullIndex, 
                         const Index & component3=nullIndex, 
                         const Index & component4=nullIndex, 
			 const Index & grid=nullIndex );
  
  // inquire whether the grid function is face-centred
  bool getIsFaceCentered(const int   & axis=forAll,
                         const Index & component0=nullIndex, 
                         const Index & component1=nullIndex, 
                         const Index & component2=nullIndex, 
                         const Index & component3=nullIndex, 
                         const Index & component4=nullIndex, 
                         const Index & grid=nullIndex ) const;

  // This next function tells you whether the grid function is face centred in a standard way
  GridFunctionParameters::faceCenteringType getFaceCentering() const;

  // Set the type of face centering, the behaviour of this function depends on whether the
  //  argument "axis" has been specified ore else the current value for getFaceCentering():
  //  (1) if "axis" is given then make all components face centred in direction=axis
  //  (2) if getFaceCentering()==all : make components face centered in all directions, the
  //      grid function should have been contructed or updated using the faceRange to specify
  //      which Index is to be used for the "directions"
  void setFaceCentering( const int & axis=defaultValue );

  // Make sure the the faceCentering of the grid collection is consistent with the
  // face centering of the MappedGridFunctions
  int updateFaceCentering();

  // specify options on hold the array data is allocated.
  void setDataAllocationOption(int option );

  // return size of this object  
  virtual real sizeOf(FILE *file = NULL ) const;

  // fixup unused points 
  virtual int fixupUnusedPoints(const RealArray & value =Overture::nullRealArray(), 
			int numberOfGhostlines=1 );

//   // zero out unused points based on the classify array in a coefficient matrix
//   int zeroUnusedPoints(GridCollectionFunction & coeff, 
//                        double value = 0,
// 		       const Index & component0=nullIndex,
// 		       const Index & component1=nullIndex,
// 		       const Index & component2=nullIndex,
// 		       const Index & component3=nullIndex,
// 		       const Index & component4=nullIndex );


  // ************************************************
  // ***** DIFFERENTIATION CLASS FUNCTIONS **********
  // ************************************************


  // *wdh* 020423 Index I1,I2,I3;

  GridCollectionOperators* getOperators() const;    // return a pointer to the operators **not virtual**
  virtual void setOperators(GenericCollectionOperators & operators );  // supply a derivative object to use

  friend GridCollectionFunction & FABS( const GridCollectionFunction & cgf );

#define COLLECTION_FUNCTION GridCollectionFunction
#define INT_COLLECTION_FUNCTION intGridCollectionFunction
#define COLLECTION GridCollection
#define INTEGRAL_TYPE double
#include "derivativeDefinitions.h"
#undef COLLECTION_FUNCTION
#undef INT_COLLECTION_FUNCTION
#undef COLLECTION 
#undef INTEGRAL_TYPE

  // ***non-member friend functions ***
  friend GridCollectionFunction & multiply( const GridCollectionFunction & a, const GridCollectionFunction & coeff );

  // use this function to  "evaluate and expression" when calling a function with an expression
  // involving grid collection functions.
  friend GridCollectionFunction evaluate( GridCollectionFunction & u );

 private:
  
  // These are used by list's of ReferenceCounting objects
  virtual void reference( const ReferenceCounting & mgf )
  { GridCollectionFunction::reference( (GridCollectionFunction&) mgf ); }
  virtual ReferenceCounting & operator=( const ReferenceCounting & mgf )
  { return GridCollectionFunction::operator=( (GridCollectionFunction&) mgf ); }
  virtual ReferenceCounting* virtualConstructor( const CopyType ct = DEEP ) const
  { return ::new GridCollectionFunction(*this,ct); }  

 protected:
  void constructor(GridCollection *gridCollection0,
                   GridCollectionData *gridCollectionData0,
		   const Range & R0, 
		   const Range & R1,
		   const Range & R2,
		   const Range & R3,
		   const Range & R4,
		   const Range & R5,
		   const Range & R6,
		   const Range & R7,
                   const bool createMappedGridFunctionList=TRUE );
  void initialize();
  void dimensionName(); // make sure name is long enough
  virtual int updateCollections(); // update refinementLevel etc.
  void setPositionOfFaceCentering(const int number);

 protected:
  aString className;
  bool temporary;                  // TRUE if this is a temporary used internally
  int dataAllocationOption;        // ==1 means do not allocate data for rectangular grids.
  
 // holds data to be reference counted
 public:
    class RCData : public ReferenceCounting
   {
    public:
     friend class GridCollectionFunction;
     RCData(); 
     virtual ~RCData();
     RCData & operator=( const RCData & rcdata );
//     IntegerArray positionOfComponent;
//     IntegerArray positionOfCoordinate;

     int positionOfComponent[maximumNumberOfIndicies];  // positions of the components
     int positionOfCoordinate[maximumNumberOfIndicies]; // positions of the coordinate directions 

     int positionOfFaceCentering;
     int faceCentering;
     Range R[MappedGridFunction::maximumNumberOfIndicies+1];
     aString *name;              // pointer to an array of names for the grid function and components
     int numberOfNames;         // length of the name array
     Interpolant *interpolant;  // pointer to an interpolant
     ListOfGridCollectionFunction multigridLevel; //  Collections of grid functions having the same multigrid level.
     ListOfGridCollectionFunction refinementLevel; //  Collections of grid functions having the same refinement level.
    private:
  
      // These are used by list's of ReferenceCounting objects
      virtual void reference( const ReferenceCounting & mgf )
      { RCData::reference( (RCData&) mgf ); }
      virtual ReferenceCounting & operator=( const ReferenceCounting & mgf )
      { return RCData::operator=( (RCData&) mgf ); }
      virtual ReferenceCounting* virtualConstructor( const CopyType ) const
      { return ::new RCData(); }  
   };
   
 protected:
  RCData *rcData;

};  



#endif 
