// This file automatically generated from CheckParallel.bC with bpp.
#include "CheckParallel.h"
#include "HDF_DataBase.h"
#include "display.h"

CheckParallel::CheckParallel()
{
    counter=0;
    dbase.put<FILE*>("file")=stdout;  // default file for output 
}

CheckParallel::~CheckParallel()
{
    if( dbase.has_key("dataBase") )
    {
        HDF_DataBase & dataBase = dbase.get<HDF_DataBase>("dataBase");

        dataBase.unmount();
    }
    

}

/// ===============================================================================================
/// \brief set the default file to write to
/// ===============================================================================================
int CheckParallel::setDefaultFile( FILE *file )
{
    dbase.get<FILE*>("file")=file;

    return 0;
}



/// ===============================================================================================
/// \brief check the difference in an array between two parallel runs
///
/// \param x (input) : array to save or check
/// \param label (input) : label
/// \param file (input) : write info to this file (normal there should be a separate file for each processor)
///
/// Usage:
///
///   (1) Run your code first with np=1 processor, calling checkDiff with the array you  want to check.
///   Check diff will save a database file with the array information.
///
///   (2) Run your code again in parallel. The call to check diff will compute the difference
///      in the array values (over those points that are common)
/// ===============================================================================================

real CheckParallel::
checkDiff( const RealArray & x, const aString & label, FILE *file_ /* = NULL */ )
{
    real maxDiff=0.;
    counter++;        // keeps track of different calls
    const int myid = max(0,Communication_Manager::My_Process_Number);
    const int np   = max(1,Communication_Manager::Number_Of_Processors);
    FILE *file=file_;
    if( file_==NULL )
        file = dbase.get<FILE*>("file");
    if( !dbase.has_key("dataBase") )
    {
        dbase.put<HDF_DataBase>("dataBase");
        HDF_DataBase & dataBase = dbase.get<HDF_DataBase>("dataBase");
        aString fileName = "checkParallel.hdf";
        if( np==1 )
            dataBase.mount(fileName,"I");  // open for writing 
        else
            dataBase.mount(fileName,"R");  // open for reading 
        dataBase.setMode(HDF_DataBase::noStreamMode);
    }
    HDF_DataBase & dataBase = dbase.get<HDF_DataBase>("dataBase");
    aString name;
    sPrintF(name,"array%d",counter);
    if( np==1 )
    {
        printF("checkParallel: np=1, save [%s] label=[%s]\n",(const char*)name,(const char*)label);
        dataBase.put(x,name);
    }
    else
    {
        RealArray x1;
        dataBase.get(x1,name);  // here is the array saved for np=1
        Range Iv[4], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2], &I4=Iv[3];
    // Could leave off parallel ghost on far left and right
        for( int dir=0; dir<4; dir++ )
            Iv[dir]=x.dimension(dir);
        RealArray xdiff(I1,I2,I3,I4);
            xdiff = fabs(x(I1,I2,I3,I4)-x1(I1,I2,I3,I4));
        maxDiff = max(xdiff);
        fprintf(file,"checkParallel: np=%d, myid=%d, compare [%s] label=[%s] max-diff=%9.3e\n",
                      np,myid,(const char*)name,(const char*)label,maxDiff);
        fprintf(file,"               check points [%d,%d][%d,%d][%d,%d][%d,%d] from source [%d,%d][%d,%d][%d,%d][%d,%d] (np=1)\n",
                        I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),I4.getBase(),I4.getBound(),
                        x1.getBase(0),x1.getBound(0),x1.getBase(1),x1.getBound(1),x1.getBase(2),x1.getBound(2),x1.getBase(3),x1.getBound(3));
        if( maxDiff > 0. )
        {
            ::display(xdiff,"|x-x(np=1)|",file,"%8.1e");
        }
    }
    return maxDiff;
}

real CheckParallel::
checkDiff( const IntegerArray & x, const aString & label, FILE *file_ /* = NULL */ )
{
    real maxDiff=0.;
    counter++;        // keeps track of different calls
    const int myid = max(0,Communication_Manager::My_Process_Number);
    const int np   = max(1,Communication_Manager::Number_Of_Processors);
    FILE *file=file_;
    if( file_==NULL )
        file = dbase.get<FILE*>("file");
    if( !dbase.has_key("dataBase") )
    {
        dbase.put<HDF_DataBase>("dataBase");
        HDF_DataBase & dataBase = dbase.get<HDF_DataBase>("dataBase");
        aString fileName = "checkParallel.hdf";
        if( np==1 )
            dataBase.mount(fileName,"I");  // open for writing 
        else
            dataBase.mount(fileName,"R");  // open for reading 
        dataBase.setMode(HDF_DataBase::noStreamMode);
    }
    HDF_DataBase & dataBase = dbase.get<HDF_DataBase>("dataBase");
    aString name;
    sPrintF(name,"array%d",counter);
    if( np==1 )
    {
        printF("checkParallel: np=1, save [%s] label=[%s]\n",(const char*)name,(const char*)label);
        dataBase.put(x,name);
    }
    else
    {
        IntegerArray x1;
        dataBase.get(x1,name);  // here is the array saved for np=1
        Range Iv[4], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2], &I4=Iv[3];
    // Could leave off parallel ghost on far left and right
        for( int dir=0; dir<4; dir++ )
            Iv[dir]=x.dimension(dir);
        IntegerArray xdiff(I1,I2,I3,I4);
            xdiff = abs(x(I1,I2,I3,I4)-x1(I1,I2,I3,I4));
        maxDiff = max(xdiff);
        fprintf(file,"checkParallel: np=%d, myid=%d, compare [%s] label=[%s] max-diff=%9.3e\n",
                      np,myid,(const char*)name,(const char*)label,maxDiff);
        fprintf(file,"               check points [%d,%d][%d,%d][%d,%d][%d,%d] from source [%d,%d][%d,%d][%d,%d][%d,%d] (np=1)\n",
                        I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),I4.getBase(),I4.getBound(),
                        x1.getBase(0),x1.getBound(0),x1.getBase(1),x1.getBound(1),x1.getBase(2),x1.getBound(2),x1.getBase(3),x1.getBound(3));
        if( maxDiff > 0. )
        {
            ::display(xdiff,"|x-x(np=1)|",file,"%8.1e");
        }
    }
    return maxDiff;
}



// =======================================================================
// Use this routine from fortran to check an array 
//
// =======================================================================

real checkParallelArrayReal( CheckParallel *pcp, char *label_,
                                                  int & nd1a, int & nd1b, int & nd2a, int & nd2b, int & nd3a, int & nd3b, int & nd4a, int & nd4b, 
                                                  int & n1a, int & n1b, int & n2a, int & n2b, int & n3a, int & n3b, int & n4a, int & n4b,
                                                  real & x_, 
                                                  int & labelLength )
{
    if( n1a>n1b || n2a>n2b || n3a>n3b || n4a>n4b )
    {
        return 0.;
    }
    Range Iv[4], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2], &I4=Iv[3];
    I1 = Range(n1a,n1b);
    I2 = Range(n2a,n2b);
    I3 = Range(n3a,n3b);
    I4 = Range(n4a,n4b);
  // In order to use CheckParallel we create an A++ array to hold the Fortran data
    RealArray xd(I1,I2,I3,I4);
    real *px=&x_;
    const int nd1=nd1b-nd1a+1;
    const int nd2=nd2b-nd2a+1;
    const int nd3=nd3b-nd3a+1;
    #define x(i1,i2,i3,i4) px[(i1-nd1a) + nd1*( i2-nd2a + nd2*( i3-nd3a + nd3*( i4-nd4a) ) ) ]
    for( int i4=n4a; i4<=n4b; i4++ )
    for( int i3=n3a; i3<=n3b; i3++ )
    for( int i2=n2a; i2<=n2b; i2++ )
    for( int i1=n1a; i1<=n1b; i1++ )
    {
        xd(i1,i2,i3,i4) = x(i1,i2,i3,i4);  // copy fortran data
    }
    string label(label_,0,labelLength);
  // remove trailing blanks
    int i= label.find_last_not_of(" "); // position of last non-blank character
    label.erase(i+1,label.size()-i);
    assert( pcp !=NULL );
    return pcp->checkDiff( xd,label );
}

real checkParallelArrayInt( CheckParallel *pcp, char *label_,
                                                  int & nd1a, int & nd1b, int & nd2a, int & nd2b, int & nd3a, int & nd3b, int & nd4a, int & nd4b, 
                                                  int & n1a, int & n1b, int & n2a, int & n2b, int & n3a, int & n3b, int & n4a, int & n4b,
                                                  int & x_, 
                                                  int & labelLength )
{
    if( n1a>n1b || n2a>n2b || n3a>n3b || n4a>n4b )
    {
        return 0.;
    }
    Range Iv[4], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2], &I4=Iv[3];
    I1 = Range(n1a,n1b);
    I2 = Range(n2a,n2b);
    I3 = Range(n3a,n3b);
    I4 = Range(n4a,n4b);
  // In order to use CheckParallel we create an A++ array to hold the Fortran data
    IntegerArray xd(I1,I2,I3,I4);
    int *px=&x_;
    const int nd1=nd1b-nd1a+1;
    const int nd2=nd2b-nd2a+1;
    const int nd3=nd3b-nd3a+1;
    #define x(i1,i2,i3,i4) px[(i1-nd1a) + nd1*( i2-nd2a + nd2*( i3-nd3a + nd3*( i4-nd4a) ) ) ]
    for( int i4=n4a; i4<=n4b; i4++ )
    for( int i3=n3a; i3<=n3b; i3++ )
    for( int i2=n2a; i2<=n2b; i2++ )
    for( int i1=n1a; i1<=n1b; i1++ )
    {
        xd(i1,i2,i3,i4) = x(i1,i2,i3,i4);  // copy fortran data
    }
    string label(label_,0,labelLength);
  // remove trailing blanks
    int i= label.find_last_not_of(" "); // position of last non-blank character
    label.erase(i+1,label.size()-i);
    assert( pcp !=NULL );
    return pcp->checkDiff( xd,label );
}

