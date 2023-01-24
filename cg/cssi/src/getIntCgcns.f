      ! This function can be called from fortran to lookup int's or enums in 
      ! a Cgcssi data base
      integer function getIntCgcssi(pdb,name,x)
      ! -- add the length of the string as an extra arg --
      implicit none
      double precision pdb
      character *(*) name
      integer x,getIntFromDataBaseCgcssi

      getIntCgcssi = getIntFromDataBaseCgcssi(pdb,name,x,len(name))
      return 
      end
