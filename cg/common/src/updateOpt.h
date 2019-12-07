#ifndef UPDATE_OPT_H
#define UPDATE_OPT_H

#define updateOpt EXTERN_C_NAME(updateopt)
#define updateOptNew EXTERN_C_NAME(updateoptnew)
extern "C"
{
   void updateOpt(const int &nd1a,const int &nd1b,const int &nd2a,const int &nd2b,
                  const int &nd3a,const int &nd3b,const int &nd4a,const int &nd4b, \
		  const int &mask,real &u1, const real&u2,  
                  const real&ut1, const real&ut2, const real&ut3, const real&ut4, 
                  const int &ipar, const real& rpar, int & ierr );

   void updateOptNew(const int &nd1a,const int &nd1b,const int &nd2a,const int &nd2b,
                  const int &nd3a,const int &nd3b,const int &nd4a,const int &nd4b, \
		  const int &mask, real &uNew,
                  const real&u1, const real&u2, const real&u3, const real&u4, const real&u5, 
                  const real&u6, const real&u7, const real&u8, const real&u9, const real&u10, 
                  const int &ipar, const real& rpar, int & ierr );
}


#endif
