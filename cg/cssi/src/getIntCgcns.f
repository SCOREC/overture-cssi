      ! This function can be called from fortran to lookup int's or enums in 
      ! a Cgcns data base
      integer function getIntCgcns(pdb,name,x)
      ! -- add the length of the string as an extra arg --
      implicit none
      double precision pdb
      character *(*) name
      integer x,getIntFromDataBaseCgcns

      getIntCgcns = getIntFromDataBaseCgcns(pdb,name,x,len(name))
      return 
      end
