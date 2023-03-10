      subroutine cgespc2( nd,ng,nv,neq,nze,nqs,ia,ja,a,peqn,bc,pkr,
     & ndrs3,mrsab,lratio,epsz,nfict, neqb,nzeb,iab,jab,ab,
     & nze1,ia1,ja1,a1,
     & neqp,nzep,iep,iap,jap,ap,iepi,debug, id,rd,ierr )
c===========================================================
c    Boundary Pre-conditionner
c
c Purpose
c   Locally invert the equations for a boundary point and fictitous
c   points. This should remove the problem of zero pivots so that
c   the Yale solver can work on a wider class of matrices. This
c   routine should also help the convergence of interative solvers.
c
c
c Input
c  nd,ng,nv,neq
c  ia,ja,a - matrix generated by cgesg
c  peqn(ng) - for the mapping from (n,i1,i2,i3,k) to equation number
c             (defined in cgesg)
c  lratio : =1 single precision, 2= double precision
c  iab(neq+1),jab(nze),ab(nze) : work space
c  nfict : number of fictitous points (1 for 2nd order, 2 for 4th order)
c Output
c
c===========================================================
      integer ia(neq+1),ja(*),bc(nd,2,ng),
     & pkr(ng),ndrs3(3,2,ng),mrsab(nd,2,ng),lratio,peqn(ng),
     & iab(neqb+1),jab(nzeb),ia1(*),ja1(*),debug,id(*)
      real epsz,a(nze),ab(nzeb),a1(*),rd(*)
      integer iep(*),iap(*),jap(*), iepi(*)
      real ap(*)

c........local
      parameter( nvd=10,nbd=7*nvd )
      real b(nbd,nbd),work(nbd),rcond
      integer eqn,ipvt(nbd),ie(nbd),kds(3,2,3,2),ifp(3,3,2),iv(3)
      logical d
      data kds/1,1,1, 1,2,2,  1,1,1, 2,1,2,  1,1,1, 2,2,1,
     &         2,1,1, 2,2,2,  1,2,1, 2,2,2,  1,1,2, 2,2,2/
      data ifp/-1,0,0, 0,-1,0, 0,0,-1, 1,0,0, 0,1,0, 0,0,1/
c......start statement functions
      ndr(kd,k)=ndrs3(kd,2,k)-ndrs3(kd,1,k)+1
      eqn(n,i1,i2,i3,k)=n+ nv*(i1-ndrs3(1,1,k)+
     &               ndr(1,k)*(i2-ndrs3(2,1,k)+
     &               ndr(2,k)*(i3-ndrs3(3,1,k))))  + peqn(k)
      kr(i1,i2,i3,k)=id( pkr(k)+i1-ndrs3(1,1,k)+ndr(1,k)*(
     &                          i2-ndrs3(2,1,k)+ndr(2,k)*(
     &                          i3-ndrs3(3,1,k))) )
      d(i)=mod(debug/2**i,2).eq.1
c......end statement functions

      if( d(3) )then
        write(1,'('' *****Entering CGESBPC2 *****'')')
      end if

c     ---sort entries in a() so that column indices are increasing
      do i=1,neq
        call bsort( ia(i+1)-ia(i),ja(ia(i)),a(ia(i)) )
      end do

      if( nv.gt.nvd )then
        write(*,'('' CGES:CGESPER: dimension error nv > nvd'')')
        write(*,'('' CGES:CGESPER: ...no pivoting performed'')')
        return
      end if

c  initialize iab(1)=1
      ie0=0      ! current equation in (iab,jab,ab)
      iab(1)=1   ! holds temporary copy of transformed boundary equations
      iap(1)=1   ! holds sparse matrix used to transform the rhs

c     =====For each side of each grid Do ======
      do k=1,ng
        do kd=1,nd
          do ks=1,2
            if( bc(kd,ks,k).gt.0 )then
c             ===boundary condition side
c             ...loop over points on this side
              if( nd.eq.2 )then
                i3a=ndrs3(3,1,k)
                i3b=ndrs3(3,2,k)
              else
                i3a=mrsab(3,kds(3,1,kd,ks),k)
                i3b=mrsab(3,kds(3,2,kd,ks),k)
              end if
              do i3=i3a,i3b
                do i2=mrsab(2,kds(2,1,kd,ks),k),
     &                mrsab(2,kds(2,2,kd,ks),k)
                  do i1=mrsab(1,kds(1,1,kd,ks),k),
     &                  mrsab(1,kds(1,2,kd,ks),k)

                    if( kr(i1,i2,i3,k).gt.0 )then
c                     ...get eqn numbers for this boundary point
c                     ie(m) m=1,...,nb : equation numbers for the
c                        boundary point and fictitious points
                      nb=0
                      do n=1,nv
                        nb=nb+1
                        ie(nb)=eqn(n,i1,i2,i3,k)
                      end do
                      do i=1,ie0
                        if( ie(1).eq.iep(i) )then
c                         skip this point - we have already done it
                          goto 100
                        end if
                      end do
c                     ...here are the neighbouring fictitous point(s):
                      do ifict=1,nfict
                        i1n=i1+ifp(1,kd,ks)*ifict
                        i2n=i2+ifp(2,kd,ks)*ifict
                        i3n=i3+ifp(3,kd,ks)*ifict
                        do n=1,nv
                          nb=nb+1
                          ie(nb)=eqn(n,i1n,i2n,i3n,k)
                        end do
                      end do
c                     ...check for corners, add points in other directions
                      iv(1)=i1
                      iv(2)=i2
                      iv(3)=i3
                      do kdd=0,nd-2
                        kd2=mod(kd+kdd,nd)+1     ! kd2 = other directions
                        do ks2=1,2
c****
*                           if( nd.eq.2 .and.
*      &                        iv(kd2).eq.mrsab(kd2,ks2,k) )then
                          if( iv(kd2).eq.mrsab(kd2,ks2,k) )then
c                           ...here are the neighbouring fictitous point(s):
*       write(*,'('' CGESBPC2: PC corner, kd,ks,kd2,ks2='',4i4)')
*      & kd,ks,kd2,ks2
                            do ifict=1,nfict
                              i1n=i1+ifp(1,kd2,ks2)*ifict
                              i2n=i2+ifp(2,kd2,ks2)*ifict
                              i3n=i3+ifp(3,kd2,ks2)*ifict
                              do n=1,nv
                                nb=nb+1
                                ie(nb)=eqn(n,i1n,i2n,i3n,k)
                              end do
                            end do
                          end if
                        end do
                      end do
                      if( nb.gt.nbd )then
                        stop 'cgespc2: ERROR nb > nbd'
                      end if
c                     ---load the matrix with the coefficients
                      do ib=1,nb
                        do jb=1,nb
                          b(ib,jb)=0.
                        end do
                        i=ie(ib)
                        do ii=ia(i),ia(i+1)-1
                          jj=ja(ii)
                          do jb=1,nb
                            if( jj.eq.ie(jb) )then
                              b(ib,jb)=a(ii)
                            end if
                          end do
                        end do
                      end do

      if( d(3) )then
      write(1,9000) kd,ks,i1,i2,i3
 9000 format(1x,' CGESPER: kd,ks,i1,i2,i3 =',2i2,3i6)
      do n=1,nb
        write(1,9100) n,ie(n),(ja(j),a(j),j=ia(ie(n)),ia(ie(n)+1)-1)
      end do
 9100 format(1x,' n,ie(n),(j,a) =',i2,i8,(1x,10(i6,e9.1)))
      write(1,'('' b:'')')
      do ib=1,nb
        write(1,9200) (b(ib,jb),jb=1,nb)
      end do
 9200 format((1x,10e9.1))
      end if

c                     --Factor the matrix : b(nb,nb)
                      if( lratio.eq.1 )then
                        call sgeco( b,nbd,nb,ipvt,rcond,work )
                      else
                        call dgeco( b,nbd,nb,ipvt,rcond,work )
                      end if
                      if( rcond.lt.epsz )then
                        write(*,'('' CGESBPC2 Warning rcond='',e12.4'//
     &                   ','' i1,i2,i3,k='',3i6,i3)') rcond,i1,i2,i3,k
                        if( rcond.eq.0. )then
                          stop 'CGESBPC2'
                        end if
                      end if
c                     --Invert the matrix
                      if( lratio.eq.1 )then
                        call sgedi( b,nbd,nb,ipvt,det,work,1 )
                      else
                        call dgedi( b,nbd,nb,ipvt,det,work,1 )
                      end if
                      if( d(3) )then
                        write(1,'('' rcond,ipvt ='',e12.4,99i4)')
     &                   rcond,(ipvt(i),i=1,nb)

                        write(1,'('' b (inverse):'')')
                        do ib=1,nb
                          write(1,9200) (b(ib,jb),jb=1,nb)
                        end do
                      end if
c                     ...Save the inverse in the sparse array (iap,jap,ap)
c                     ...Change the equations in (ia,ja,a) by multiplying
c                        through by the inverse -> (iab,jab,ab)
                      do ib=1,nb
                        ie0=ie0+1
                        if( ie0.gt.neqp )then
                          stop 'CGESPC2: dimension error neqp too small'
                        end if
                        iep(ie0)=ie(ib)  ! equation number
                        iap(ie0+1)=iap(ie0)+nb
                        if( iap(ie0+1)-1.gt.nzep )then
                          stop 'CGESPC2: dimension error nzep too small'
                        end if
                        do j=iap(ie0),iap(ie0+1)-1
                          jb=j-iap(ie0)+1
                          jap(j)=ie(jb)
                          ap(j)=b(ib,jb)
c**                       c(jb)=b(ib,jb)
                        end do
c                       ...ab(i0,j) <- sum c(n)*a(ie(n),j)  n=1,...,ne
c                          c(n)=ap(iap(ie0)+n-1)
                        call cgessra( ie0,nb,ie,ap(iap(ie0)),
     &                   neq,ia,ja,a,iab,jab,ab,epsz )
                        if( iab(ie0+1)-1.gt.nzeb )then
                          write(*,'('' CGESPC2: nzeb too small '')')
                          write(*,'('' nzeb= '',i8)') nzeb
                          write(*,'('' CGESPC2: increase zratio '')')
                          stop 'CGESPC2'
                        end if
      if( d(4) )then
        write(1,9300) ie0,(jab(j),ab(j),j=iab(ie0),iab(ie0+1)-1)
      end if
                      end do
 9300 format(1x,' ie0,(jab,ab) =',i8,(1x,10(i6,e9.1)))

                    end if
 100                continue
                  end do
                end do
              end do
            end if
          end do
        end do
      end do

      if( d(3) )then
        write(*,'('' CGESPC2: neqb, neqb(true) ='',2i8)') neqb,ie0
        write(*,'('' CGESPC2: nzeb, nzeb(true) ='',2i8)') nzeb,
     &   iab(ie0+1)-1
        write(*,'('' CGESPC2: neqp, neqp(true) ='',2i8)') neqp,ie0
        write(*,'('' CGESPC2: nzep, nzep(true) ='',2i8)') nzep,
     &   iap(ie0+1)-1
      end if
      neqp=ie0 ! number of equations in sparse matrix for rhs
      nzep=iap(neqp+1)-1

c     ---sort the iep array into increasing order
      do i=1,neqp
        iepi(i)=i
      end do
*       if( d(1) ) write(*,'('' CGESBPC2: begin sort(2)...'')')
*       call bsorti( neqp,iep,iepi )
      call sortii( neqp,iep,iepi )
*       if( d(1) ) write(*,'('' CGESBPC2: end sort(2)...'')')
*       write(1,'('' neqp,nzep = '',2i6)') neqp,nzep
*       write(1,'('' Sorted iep:'',/,(1x,20i6))') (iep(i),i=1,neqp)
*       write(1,'('' Sorted iepi:'',/,(1x,20i6))') (iepi(i),i=1,neqp)


c     === form new sparse matrix with new equations at boundary replacing
c         the old equations
      m=1
      ia1(1)=1
      iep(neqp+1)=neq+1 ! special value at end
      do n=1,neq
        if( n.lt.iep(m) )then
          ia1(n+1)=ia1(n)+ia(n+1)-ia(n)
          if( ia1(n+1)-1.gt.nze1 )then
            write(*,'('' CGES:CGESPC: error nze1 too small '')')
            write(*,'('' nze1= '',i8)') nze1
            stop 'CGESPC2'
          end if
          do j=ia1(n),ia1(n+1)-1
            ja1(j)=ja(j-ia1(n)+ia(n))
            a1(j)=  a(j-ia1(n)+ia(n) )
          end do
        else
c         ---here is a new boundary equation
          m0=iepi(m)
          ia1(n+1)=ia1(n)+iab(m0+1)-iab(m0)
          if( ia1(n+1)-1.gt.nze1 )then
            write(*,'('' CGES:CGESPC: error nze1 too small '')')
            write(*,'('' nze1= '',i8)') nze1
            stop 'CGESPC2'
          end if
          do j=ia1(n),ia1(n+1)-1
            ja1(j)=jab(j-ia1(n)+iab(m0))
            a1(j)=  ab(j-ia1(n)+iab(m0) )
          end do
c         ---increment m
          m=m+1
        end if
      end do
      if( d(2) )then
        write(1,'('' CGESBPC: Old nze, new nze ='',2i8)') nze,
     &   ia1(neq+1)-1
      end if
      nzenew=ia1(neq+1)-1
      if( nzenew.gt.nqs )then
        write(*,'('' CGESPC2: old nze, new nze ='',2i9)') nze,nzenew
        write(*,'('' CGESPC2: ERROR nze > nqs, nze,nqs='',2i9)')
     &    nzenew,nqs
        write(*,'('' CGESPC2: Increase zratio to >='',f12.4)')
     &   (nzenew+1.)/neq
        stop
      end if

      nze=ia1(neq+1)-1  ! new number of non-zeroes

c     ---copy back into ia,ja,a
      do i=1,neq+1
        ia(i)=ia1(i)
      end do
      do i=1,nze
        ja(i)=ja1(i)
        a(i)=a1(i)
      end do
*       write(1,*) ' Matrix at end of CGESPC2...'
*       do n=1,neq
*         write(1,9700) n,(ja(j),a(j),j=ia(n),ia(n+1)-1)
*       end do
*  9700 format(1x,' n,(j,a) =',i8,(1x,10(i6,e9.1)))

c****** unsort the iep array - should really sort (iap,jap,ap)
      do i=1,neqp
        ia1(i)=iep(i)
      end do
      do i=1,neqp
        iep(iepi(i))=ia1(i)
      end do
      end

      subroutine cgespc3( neqp,nzep,iep0,iap0,jap0,ap0,
     & iep,iap,jap,ap )
c====================================================================
c Copy the (iep,iap,jap,ap) arrays from their temporary spot
c====================================================================
      integer iep0(*),iap0(*),jap0(*),iep(*),iap(*),jap(*)
      real ap0(*),ap(*)

      do i=1,neqp
        iep(i)=iep0(i)
      end do
      do i=1,neqp+1
        iap(i)=iap0(i)
      end do
      do i=1,nzep
        jap(i)=jap0(i)
        ap(i)=ap0(i)
      end do
*       write(1,*) ' Matrix at end of CGESPC3...'
*       do n=1,neqp
*         write(1,9700) n,iep(n),(jap(j),ap(j),j=iap(n),iap(n+1)-1)
*       end do
*  9700 format(1x,' n,iep(n),(jap,ap) =',2i8,(1x,10(i6,e9.1)))
      end

      subroutine cgespc4( neq,nze,ia,iabpc )
c====================================================================
c Compress the ia() array
c    *** Assumes that the rows in ia() are in order ***
c Input
c  neq,nze,ia
c Output
c  iabpc()
c====================================================================
      integer ia(nze),iabpc(neq+1)

      iabpc(1)=1
      m=1
      do i=1,nze
        if( ia(i).gt.m )then
          iabpc(m+1)=i
          m=m+1
        end if
      end do
      iabpc(neq+1)=nze+1

      end

      subroutine cgespc5( neq,nze,ia,iabpc )
c====================================================================
c Un-Compress the ia() array
c Input
c  neq,nze,iabpc
c Output
c  ia()
c====================================================================
      integer ia(nze),iabpc(neq+1)

      do m=1,neq
        do j=iabpc(m),iabpc(m+1)-1
          ia(j)=m
        end do
      end do

      end

      subroutine sort(rw, n)
c====================================================================
c     --- supplied by WDH
c====================================================================
      integer rw(n),temp
c
      do i=1,n-1
        ichange=0
        do j=1,n-i
          if( rw(j).gt.rw(j+1) ) then
            ichange=ichange+1
            temp=rw(j)
            rw(j)=rw(j+1)
            rw(j+1)=temp
          end if
        end do
        if( ichange.eq.0 ) goto 100
      end do
 100  continue

      end
      subroutine bsort( n,iperm,aperm )
c****
c**       BUBBLE SORT FOR INTEGER ARRAYS
c** Sort the elements of iperm(i),i=1,n into ascending order.
c** Permute the elements of aperm(i),i=1,n as well.
c****
c
      integer n,iperm(n)
      real aperm(n)
c
      do i=1,n-1
        ichange=0
        do j=1,n-i
          if( iperm(j).gt.iperm(j+1) ) then
            ichange=ichange+1
            itemp=iperm(j)
            iperm(j)=iperm(j+1)
            iperm(j+1)=itemp
            atemp=aperm(j)
            aperm(j)=aperm(j+1)
            aperm(j+1)=atemp
          end if
        end do
        if( ichange.eq.0 ) goto 100
      end do
 100  continue

      end

      subroutine bsorti( n,iperm,aperm )
c****
c**       BUBBLE SORT FOR INTEGER ARRAYS
c** Sort the elements of iperm(i),i=1,n into ascending order.
c** Permute the elements of aperm(i),i=1,n as well.
c****
c
      integer n,iperm(n)
      integer aperm(n),atemp
c
      do i=1,n-1
        ichange=0
        do j=1,n-i
          if( iperm(j).gt.iperm(j+1) ) then
            ichange=ichange+1
            itemp=iperm(j)
            iperm(j)=iperm(j+1)
            iperm(j+1)=itemp
            atemp=aperm(j)
            aperm(j)=aperm(j+1)
            aperm(j+1)=atemp
          end if
        end do
        if( ichange.eq.0 ) goto 100
      end do
 100  continue

      end

      subroutine sortir( la,ia,ra )
c===============================================================
c Quicksort
c
c   Sort the integer array ia into increasing order, permute
c   the real array ra as well
c
c Input -
c    la,ia,ra
c
c
c===============================================================
      integer            la
      real               ra(la),rat,ratt
      integer            ia(la)
      integer            iu(21),il(21),i,m,j,k,ij,l
      integer            t,tt,r
      if (la.le.0) return
      m = 1
      i = 1
      j = la
      r = .375
    5 if (i.eq.j) go to 45
      if (r.gt..5898437) go to 10
      r = r+3.90625e-2
      go to 15
   10 r = r-.21875
   15 k = i
c                                  select a central element of the
c                                  array and save it in location t
      ij = i+(j-i)*r
      t = ia(ij)
      rat = ra(ij)
c                                  if first element of array is greater
c                                  than t, interchange with t
      if (ia(i).le.t) go to 20
      ia(ij) = ia(i)
      ia(i) = t
      t = ia(ij)
      ra(ij) = ra(i)
      ra(i) = rat
      rat = ra(ij)
   20 l = j
c                                  if last element of array is less than
c                                  t, interchange with t
      if (ia(j).ge.t) go to 30
      ia(ij) = ia(j)
      ia(j) = t
      t = ia(ij)
      ra(ij) = ra(j)
      ra(j) = rat
      rat = ra(ij)
c                                  if first element of array is greater
c                                  than t, interchange with t
      if (ia(i).le.t) go to 30
      ia(ij) = ia(i)
      ia(i) = t
      t = ia(ij)
      ra(ij) = ra(i)
      ra(i) = rat
      rat = ra(ij)
      go to 30
   25 if (ia(l).eq.ia(k)) go to 30
      tt = ia(l)
      ia(l) = ia(k)
      ia(k) = tt
      ratt = ra(l)
      ra(l) = ra(k)
      ra(k) = ratt
c                                  find an element in the second half of
c                                  the array which is smaller than t
   30 l = l-1
      if (ia(l).gt.t) go to 30
c                                  find an element in the first half of
c                                  the array which is greater than t
   35 k = k+1
      if (ia(k).lt.t) go to 35
c                                  interchange these elements
      if (k.le.l) go to 25
c                                  save upper and lower subscripts of
c                                  the array yet to be sorted
      if (l-i.le.j-k) go to 40
      il(m) = i
      iu(m) = l
      i = k
      m = m+1
      go to 50
   40 il(m) = k
      iu(m) = j
      j = l
      m = m+1
      go to 50
c                                  begin again on another portion of
c                                  the unsorted array
   45 m = m-1
      if (m.eq.0) return
      i = il(m)
      j = iu(m)
   50 if (j-i.ge.11) go to 15
      if (i.eq.1) go to 5
      i = i-1
   55 i = i+1
      if (i.eq.j) go to 45
      t = ia(i+1)
      rat = ra(i+1)
      if (ia(i).le.t) go to 55
      k = i
   60 ia(k+1) = ia(k)
      ra(k+1) = ra(k)
      k = k-1
      if (t.lt.ia(k)) go to 60
      ia(k+1) = t
      ra(k+1) = rat
      go to 55
      end
c
c   purpose             - sorting of arrays by algebraic value
c                       - permutations returned
c
c   routine base on VSRTR from IMSL
c
c   usage               - call sorti (ia,la,ir)
c
c   arguments    ia     - on input, a contains the array to be sorted.
c                         on output, a contains the sorted array.
c                la     - input variable containing the number of
c                           elements in the array to be sorted.
c                ir     - vector of length la.
c                         on input, ir contains the integer values
c                           1,2,...,la. see remarks.
c                         on output, ir contains a record of the
c                           permutations made on the vector a.
c
c   remarks      the vector ir must be initialized before entering
c                sorti.  ordinarily, ir(1)=1, ir(2)=2, ...,
c                ir(la)=la.  for wider applicability, any integer
c                that is to be associated with ia(i) for i=1,2,...,la
c                may be entered into ir(i).
c
      subroutine sortii(la,ia,ir)
c                                  specifications for arguments
      integer            la,ir(la)
      integer            ia(la)
c                                  specifications for local variables
      integer            iu(21),il(21),i,m,j,k,ij,it,l,itt
      integer            t,tt,r
c                                  first executable statement
      if (la.le.0) return
      m = 1
      i = 1
      j = la
      r = .375
    5 if (i.eq.j) go to 45
      if (r.gt..5898437) go to 10
      r = r+3.90625e-2
      go to 15
   10 r = r-.21875
   15 k = i
c                                  select a central element of the
c                                  array and save it in location t
      ij = i+(j-i)*r
      t = ia(ij)
      it = ir(ij)
c                                  if first element of array is greater
c                                  than t, interchange with t
      if (ia(i).le.t) go to 20
      ia(ij) = ia(i)
      ia(i) = t
      t = ia(ij)
      ir(ij) = ir(i)
      ir(i) = it
      it = ir(ij)
   20 l = j
c                                  if last element of array is less than
c                                  t, interchange with t
      if (ia(j).ge.t) go to 30
      ia(ij) = ia(j)
      ia(j) = t
      t = ia(ij)
      ir(ij) = ir(j)
      ir(j) = it
      it = ir(ij)
c                                  if first element of array is greater
c                                  than t, interchange with t
      if (ia(i).le.t) go to 30
      ia(ij) = ia(i)
      ia(i) = t
      t = ia(ij)
      ir(ij) = ir(i)
      ir(i) = it
      it = ir(ij)
      go to 30
   25 if (ia(l).eq.ia(k)) go to 30
      tt = ia(l)
      ia(l) = ia(k)
      ia(k) = tt
      itt = ir(l)
      ir(l) = ir(k)
      ir(k) = itt
c                                  find an element in the second half of
c                                  the array which is smaller than t
   30 l = l-1
      if (ia(l).gt.t) go to 30
c                                  find an element in the first half of
c                                  the array which is greater than t
   35 k = k+1
      if (ia(k).lt.t) go to 35
c                                  interchange these elements
      if (k.le.l) go to 25
c                                  save upper and lower subscripts of
c                                  the array yet to be sorted
      if (l-i.le.j-k) go to 40
      il(m) = i
      iu(m) = l
      i = k
      m = m+1
      go to 50
   40 il(m) = k
      iu(m) = j
      j = l
      m = m+1
      go to 50
c                                  begin again on another portion of
c                                  the unsorted array
   45 m = m-1
      if (m.eq.0) return
      i = il(m)
      j = iu(m)
   50 if (j-i.ge.11) go to 15
      if (i.eq.1) go to 5
      i = i-1
   55 i = i+1
      if (i.eq.j) go to 45
      t = ia(i+1)
      it = ir(i+1)
      if (ia(i).le.t) go to 55
      k = i
   60 ia(k+1) = ia(k)
      ir(k+1) = ir(k)
      k = k-1
      if (t.lt.ia(k)) go to 60
      ia(k+1) = t
      ir(k+1) = it
      go to 55
      end
