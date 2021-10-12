subroutine mbe1d(rpar,ipar,n,uout)
! subroutine mbe1d(epsilon,dt,xa,xb,n,tfinal,uout)
    ! 1D single domain second order (full model) in domain [-l,l] with hardcoded 1 polarization and 2 atomic levels (hardcoded coefficients and thus 1 difference variable)
    ! E_tt - E_xx = epsilon*alphaP*P_tt
    ! P_tt + P = epsilon*(N0-N1)*E
    ! D_t = -epsilon*E*P_t (D=N0-N1)
    implicit none

    integer n,i,iOrder,ng,tz
    character(len=32) arg,filename1,filename2,filename3
    real*8 tfinal,epsilon,l
    real*8 dx,t,dt,cfl,c,alphaP,errE,errP,errN,eE,eP,eN
    integer it,Nt
    integer nP,nN
    integer writeFile
    real*8 xa,xb
    real*8 um(-2:n+2,3),uout(-2:n+2,3)
    real*8, save, allocatable :: un(:,:)
    real*8, save, allocatable :: u(:,:)
    real*8, save :: tcurrent
    real*8 rpar(0:*)
    integer ipar(0:*)
    
    real*8, allocatable :: xvec(:)
    ! real*8, allocatable :: um(:,:),u(:,:)!,un(:,:)
    real*8 b0,b1,pnec,peptc,etmp,ptmp,dtmp,et,ett,ettt,etttt,pt,ptt,pttt,ptttt

    common tz

    if (.not.allocated(u)) allocate(u(-2:n+2,3))
    if (.not.allocated(un)) allocate(un(-2:n+2,3))

    ! output file or not
    writeFile = 0

    ! call getarg(1,arg)
    ! read(arg,*) n
    ! n = 2**n-1

    tfinal = rpar(0)
    dt = rpar(1)
    xa = rpar(2)
    xb = rpar(3)
    epsilon = rpar(4)

    ! n = ipar(0)
    iOrder = ipar(1)

    ! iOrder=4 ! currently only support second order version
    ng = 2!iOrder/2

    ! print *, 'mbe1d--------------'
    ! print *, t,dt,xa,xb,epsilon,iOrder,n

    ! ! tz
    ! call getarg(2,arg)
    ! read(arg,*) tz ! tz=[0|1] with 0 for "exact" solution and 1 for manufactured solution
    tz = 0

    ! call getarg(3,arg)
    ! read(arg,*) tfinal

    ! file output
    if (writeFile.eq.1) then
        write (filename1,  "('Embe1dN',I4.4,'.txt')") n+1
        write (filename2,  "('Pmbe1dN',I4.4,'.txt')") n+1
        write (filename3,  "('Dmbe1dN',I4.4,'.txt')") n+1
        open(unit=111,file=filename1)
        open(unit=222,file=filename2)
        open(unit=333,file=filename3)
    endif

    ! hard coded for 1 polarization and 2 atomic levels
    nP=1
    nN=1 ! 1 for N_0-N_1

    ! 
    ! l=1000.
    ! print *, "time step 1", dt
    dx = (xb-xa)/n ! domain [-l,l]
    allocate(xvec(-ng:n+ng))
    do i=-ng,n+ng
        xvec(i)=xa+dx*i
    enddo
    c = 1.
    ! cfl=0.9
    ! dt = cfl*dx/c
    Nt = ceiling(abs(tfinal)/dt)
    if (tfinal<0.) then
        dt = -dt
    endif
    ! if (Nt.ne.0) then
    !     dt = tfinal/Nt
    ! endif
    ! print *, "Fortran 1D code, tfinal diff 1", Nt*dt-tfinal
    if (abs(Nt*dt-tfinal)>1.d-10) Nt = floor(abs(tfinal)/dt)
    ! print *, "Fortran 1D code, tfinal diff 2", Nt*dt-tfinal
    ! print *, "time step 2", dt

    ! print *, dx, dt, Nt, tfinal
    ! print *, tfinal, dt,'Number of time steps is:', Nt
    
    ! hardcoded coefficients
    ! E
    ! epsilon = 0.01d0
    alphaP = epsilon
    ! P
    pnec= epsilon ! \Delta N computed as N0-N1
    b0 = 1.d0
    b1 = 0.d0 ! no damping
    ! N
    peptc = -epsilon
    
    ! call setupInitials(n,ng,nP,nN,xvec,u,um,dt,dx,epsilon)

    eE = 0.d0
    eP = 0.d0
    eN = 0.d0

    if (Nt.eq.0) then
        call setupInitials(n,ng,nP,nN,xvec,u,um,dt,dx,epsilon,c,alphaP,b0,b1,pnec,peptc)
        t = 0.
        un = u ! value at t=0.
        u = um

        uout = un
        do i=-ng,n+ng
            uout(i,2) = epsilon*un(i,2)
        enddo

        !
        tcurrent = t
        ! print *, 'Case 1',t

    elseif (Nt.eq.1 .and. dt<0.) then
        t = dt
        call setupInitials(n,ng,nP,nN,xvec,u,um,-dt,dx,epsilon,c,alphaP,b0,b1,pnec,peptc)
        un = um ! value at t=0.

        uout = un
        do i=-ng,n+ng
            uout(i,2) = epsilon*un(i,2)
        enddo

        !
        tcurrent = t
        ! print *, 'Case 2',t

    elseif (Nt.gt.0 .and. dt>0. .and. abs(tfinal-tcurrent)<1.d-10) then ! still at the same time level
        t = Nt*dt

        ! maybe do nothing?
        ! u = u
        ! un = un
        uout = un
        do i=-ng,n+ng
            uout(i,2) = epsilon*un(i,2)
        enddo

        ! store current time
        tcurrent = t

        ! print *, 'Case 3',tcurrent, tfinal

    else ! march 1 time step to the next time level
        ! do it=1,Nt
            t = Nt*dt

            ! previous values
            ! do i=-ng,n+ng
            !     do j=1,3
            !         um(i,j) = u(i,j)
            !     enddo
            ! enddo

            ! do i=-ng,n+ng
            !     do j=1,3
            !         u(i,j) = un(i,j)
            !     enddo
            ! enddo
            um = u
            u = un

            ! march one time step
            call timeStep1D(pnec,peptc,ng,n,alphaP,dx,dt,xvec,t,b0,b1,c,um,u,un,iOrder)

            ! do i = 0,n
            !     print *, un(i,1),un(i,2),un(i,3)
            ! enddo

            call setupBC(dx,dt,t,ng,n,nP,nN,xvec,u,un,epsilon,iOrder,c,alphaP,b0,b1,pnec,peptc)

            ! compute error
            ! if (tz.eq.0) then
            !     errE = 0.
            !     errP = 0.
            !     errN = 0.
            !     do i=0,n
            !         call sechWave(etmp,et,ett,ptmp,pt,ptt,dtmp,xvec(i),t,epsilon)
            !         ! approximate errors (fast vs slow)
            !         errE = max(errE,abs(un(i,1)-etmp))
            !         errP = max(errP,abs(un(i,2)-ptmp))
            !         errN = max(errN,abs(un(i,3)-dtmp))
            !     enddo
            !     eE = max(eE,errE)
            !     eP = max(eP,errP)
            !     eN = max(eN,errN)
            !     ! print *, 'Time Step:',it, eE,eP,eN

            ! elseif (tz.eq.1) then
            !     call getError(n,errE,errP,errN,t,xvec,ng,nP,nN,un)
            !     ! print *, it,errE,errP,errN
            !     eE = max(eE,errE)
            !     eP = max(eP,errP)
            !     eN = max(eN,errN)
            ! endif   

            !
            ! um = u
            ! u = un
            ! output results for signaling problem
            if (tz.eq.0 .and. writeFile.eq.1) then
                do i=0,n
                    write(111,*) un(i,1)
                    write(222,*) un(i,2)
                    write(333,*) un(i,3)
                enddo
            endif
            tcurrent = tfinal
        ! enddo

        ! print *, 'Case 4',tcurrent,tfinal

        ! adjust P
        uout = un
        do i=-ng,n+ng
            uout(i,2) = epsilon*un(i,2)
        enddo


        ! print *, 'Max errors of (E,P,D) in space and time at time:', t,Nt
        ! print *, eE,eP,eN

    endif

    if (tz.eq.0) then
        errE = 0.
        errP = 0.
        errN = 0.
        do i=0,n
            call sechWave(etmp,et,ett,ettt,etttt,ptmp,pt,ptt,pttt,ptttt,dtmp,xvec(i),t,epsilon,c,alphaP,b0,b1,pnec,peptc)
            ! approximate errors (fast vs slow)
            errE = max(errE,abs(un(i,1)-etmp))
            errP = max(errP,abs(un(i,2)-ptmp))
            errN = max(errN,abs(un(i,3)-dtmp))
        enddo
        eE = max(eE,errE)
        eP = max(eP,errP)
        eN = max(eN,errN)
        ! print *, 'Time Step:',it, eE,eP,eN

    elseif (tz.eq.1) then
        call getError(n,errE,errP,errN,t,xvec,ng,nP,nN,un)
        ! print *, it,errE,errP,errN
        eE = max(eE,errE)
        eP = max(eP,errP)
        eN = max(eN,errN)
    endif   

    ! print *, 'Max errors of (E,P,D) in space and time at time:', t,Nt
    ! print *, eE,eP,eN

    ! output results for signaling problem
    ! if (tz.eq.0 .and. writeFile.eq.1) then
    !     do i=0,n
    !         write(111,*) u(i,1)
    !         write(222,*) u(i,2)
    !         write(333,*) u(i,3)
    !     enddo
    ! endif

    ! clean
    deallocate(xvec)
    if (writeFile.eq.1) then
        close(111)
        close(222)
        close(333)
    endif

! end subroutine mbe1d
contains

    subroutine setupInitials(n,ng,nP,nN,xvec,u,um,dt,dx,epsilon,c,alphaP,b0,b1,pnec,peptc)
        implicit none
        integer i,j,ng,nP,nN,n,tz
        real*8, allocatable :: xvec(:)
        ! real*8, allocatable :: um(:,:),u(:,:)
        real*8 um(-2:n+2,3),u(-2:n+2,3)
        real*8 dt,tmp1,tmp2,dx,epsilon
        real*8 etmp1,ptmp1,dtmp1,etmp2,ptmp2,dtmp2,et,ett,ettt,etttt,pt,ptt,pttt,ptttt
        real*8 c,alphaP,b0,b1,pnec,peptc
        common tz
           
        do i=-ng,n+ng
            
            if (tz.eq.0) then 
                ! call sechWave(etmp1,et,ett,ptmp1,pt,ptt,dtmp1,xvec(i),-dt,epsilon) ! filling initial conditions
                call sechWave(etmp2,et,ett,ettt,etttt,ptmp2,pt,ptt,pttt,ptttt,dtmp2,xvec(i),0.d0,epsilon,c,alphaP,b0,b1,pnec,peptc)
                ! um(i,1) = etmp1
                um(i,1) = etmp2-et*dt+dt**2/2.*ett-dt**3/6.*ettt+dt**4/24.*etttt
                u(i,1)  = etmp2
                ! um(i,2) = ptmp1
                um(i,2) = ptmp2-pt*dt+dt**2/2.*ptt-dt**3/6.*pttt+dt**4/24.*ptttt
                u(i,2)  = ptmp2
                ! um(i,3) = dtmp1
                u(i,3)  = dtmp2
                ! print *, etmp1-(etmp2-et*dt-dt**2/2.*ett),ptmp1-(ptmp2-pt*dt-dt**2/2.*ptt)
            elseif (tz.eq.1) then ! twilight zone, polynomial solution
                do j=1,1+nP+nN
                    call getExactSolution(tmp1,xvec(i),-dt)
                    call getExactSolution(tmp2,xvec(i),0.d0)
                    um(i,j) = tmp1
                    u(i,j)  = tmp2
                enddo
            endif 
            
        enddo

    end subroutine setupInitials

    subroutine sechWave(e,et,ett,ettt,etttt,p,pt,ptt,pttt,ptttt,d,x,t,epsilon,c,alphaP,b0,b1,pnec,peptc) ! incoming wave
        implicit none
        real*8 x,t,eta,eps,U,x0,x2,U2,epsilon,epsilon2,eps2,x1,eps1,U1
        real*8 e,p,d,et,ett,ettt,etttt,pt,ptt,pttt,ptttt
        real*8 e2,p2,d2,et2,ett2,ettt2,etttt2,pt2,ptt2,pttt2,ptttt2
        real*8 c,alphaP,b0,b1,pnec,peptc

        ! first soliton
        eta = 1.
        U = 3./4.
        x0 = 300.
        eps = epsilon/2.*sqrt(eta/(U-U*U))
        U1 = U
        x1 = x0
        eps1 = eps
        if (abs(x-x0)>500.) then
            e=2.*sqrt(eta*U/(1.-U))/cosh(eps*(x-x0-U*t))*sin(x-t)
            et=0.
            ett=0.
            ettt=0.
            etttt=0.
            p=2.*tanh(eps*(x-x0-U*t))/cosh(eps*(x-x0-U*t))*cos(x-t)
            pt=0.
            ptt=0.
            pttt=0.
            ptttt=0.
            ! d=1-2./(cosh(eps*(x-x0-U*t)))**2
        else
            call getDerivs(e,et,ett,ettt,etttt,p,pt,ptt,pttt,ptttt,d,eta,U,x,x0,t,eps,c,alphaP,b0,b1,pnec,peptc)
        endif
        ! print *, 'soliton1', eta,U,x,x0,t,eps,c,alphaP,b0,b1,pnec,peptc
        ! print *, 'soliton1 at t: ',t,'with values: ',e,et,ett,ettt,etttt,p,pt,ptt,pttt,ptttt,d

        ! second soliton
        U2 = 1./4.
        ! epsilon2 = epsilon*sqrt((U2-U2*U2)/(U-U*U))
        x0 = 900.
        eps = epsilon/2.*sqrt(eta/(U2-U2*U2))
        x2 = x0
        eps2 = eps
        if (abs(x-x0)>500.) then
            e2=2.*sqrt(eta*U2/(1.-U2))/cosh(eps2*(x-x0-U2*t))*sin(x-t)
            et2=0.
            ett2=0.
            ettt2=0.
            etttt2=0.
            p2=2.*tanh(eps2*(x-x0-U2*t))/cosh(eps2*(x-x0-U2*t))*cos(x-t)
            pt2=0.
            ptt2=0.
            pttt2=0.
            ptttt2=0.
            ! d2=1-2./(cosh(eps*(x-x0-U2*t)))**2
        else
            call getDerivs(e2,et2,ett2,ettt2,etttt2,p2,pt2,ptt2,pttt2,ptttt2,d2,eta,U2,x,x0,t,eps,c,alphaP,b0,b1,pnec,peptc)
        endif
        ! print *, 'soliton2', eta,U2,x,x0,t,eps,c,alphaP,b0,b1,pnec,peptc
        ! print *, 'soliton2 at t: ',t,'with values: ',e2,et2,ett2,ettt2,etttt2,p2,pt2,ptt2,pttt2,ptttt2,d2

        ! sum
        e = e+e2
        et = et + et2
        ett = ett + ett2
        ettt = ettt + ettt2
        etttt = etttt + etttt2

        p = p + p2
        pt = pt + pt2
        ptt = ptt + ptt2
        pttt = pttt + pttt2
        ptttt = ptttt + ptttt2

        ! d = d + d2
        d=1-2./(cosh(eps1*(x-x1-U1*t)))**2-2./(cosh(eps2*(x-x2-U2*t)))**2


    end subroutine sechWave

    subroutine getDerivs(e,et,ett,ettt,etttt,p,pt,ptt,pttt,ptttt,d,eta,U,x,x0,t,eps,c,alphaP,b0,b1,pnec,peptc)
        implicit none
        real*8 x,t,eta,eps,U,x0,e,p,d
        real*8 et,ett,ettt,etttt,pt,ptt,pttt,ptttt
        real*8 c,alphaP,b0,b1,pnec,peptc
        real*8 ex,exx,exxt,exxtt,exxxx
        real*8 pxx,ptxx,pttxx
        real*8 dt,dtt,dx,dxx
        ! E
        e = 2.*sqrt(eta*U/(1.-U))/cosh(eps*(x-x0-U*t))*sin(x-t)
        ! et = -2*sqrt(eta*U/(1 - U))*eps*U/cosh(eps*(-U*t + x - x0))*tanh(eps*(-U*t + x - x0))*sin(-x + t) & 
        !     - 2*sqrt(eta*U/(1 - U))/cosh(eps*(-U*t + x - x0))*cos(-x + t)
        ! ett = -2*sqrt(eta*U/(1 - U))*eps**2*U**2/cosh(eps*(-U*t + x - x0))*(tanh(eps*(-U*t + x - x0)))**2*sin(-x + t) &
        !      + 2*sqrt(eta*U/(1 - U))*eps**2*U**2/cosh(eps*(-U*t + x - x0))*(1 - (tanh(eps*(-U*t + x - x0)))**2)*sin(-x + t) &
        !      - 4*sqrt(eta*U/(1 - U))*eps*U/cosh(eps*(-U*t + x - x0))*tanh(eps*(-U*t + x - x0))*cos(-x + t) &
        !      + 2*sqrt(eta*U/(1 - U))/cosh(eps*(-U*t + x - x0))*sin(-x + t)
        et = - 2*sqrt(eta*U/(1. - U))*sin(-x + t)*eps*U*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
             - 2*sqrt(eta*U/(1. - U))*cos(-x + t)/cosh(eps*(-U*t + x - x0))
        exxt = -12.*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**3*sinh(eps*(-U*t + x - x0))**3*U/cosh(eps*(-U*t + x - x0))**4 &
               - 4.*sqrt(eta*U/(1. - U))*cos(-x + t)*eps**2*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
               + 10.*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**3*sinh(eps*(-U*t + x - x0))*U/cosh(eps*(-U*t + x - x0))**2 &
               - 8.*sqrt(eta*U/(1. - U))*cos(-x + t)*eps**2*sinh(eps*(-U*t + x - x0))**2*U/cosh(eps*(-U*t + x - x0))**3 &
               + 4.*sqrt(eta*U/(1. - U))*sin(-x + t)*eps*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
               + 4.*sqrt(eta*U/(1. - U))*cos(-x + t)*eps**2*U/cosh(eps*(-U*t + x - x0)) &
               + 2.*sqrt(eta*U/(1. - U))*cos(-x + t)*eps**2/cosh(eps*(-U*t + x - x0)) &
               + 2.*sqrt(eta*U/(1. - U))*sin(-x + t)*eps*U*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
               + 2.*sqrt(eta*U/(1. - U))*cos(-x + t)/cosh(eps*(-U*t + x - x0))
        ex = 2.*sqrt(eta*U/(1. - U))*sin(-x + t)*eps*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
             + 2.*sqrt(eta*U/(1. - U))*cos(-x + t)/cosh(eps*(-U*t + x - x0))
        exx = - 4.*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**2*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
              - 4.*sqrt(eta*U/(1. - U))*cos(-x + t)*eps*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
              + 2.*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**2/cosh(eps*(-U*t + x - x0)) &
              + 2.*sqrt(eta*U/(1. - U))*sin(-x + t)/cosh(eps*(-U*t + x - x0))
        exxxx = - 48.*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**4*sinh(eps*(-U*t + x - x0))**4/cosh(eps*(-U*t + x - x0))**5 &
                - 48.*sqrt(eta*U/(1. - U))*cos(-x + t)*eps**3*sinh(eps*(-U*t + x - x0))**3/cosh(eps*(-U*t + x - x0))**4 &
                + 56.*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**4*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
                + 24.*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**2*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
                + 40.*sqrt(eta*U/(1. - U))*cos(-x + t)*eps**3*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
                - 10.*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**4/cosh(eps*(-U*t + x - x0)) &
                + 8.*sqrt(eta*U/(1. - U))*cos(-x + t)*eps*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
                - 12.*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**2/cosh(eps*(-U*t + x - x0)) &
                - 2.*sqrt(eta*U/(1. - U))*sin(-x + t)/cosh(eps*(-U*t + x - x0))
        ! print *, 'e,ex,exx,exxxx',e,ex,exx,exxxx
        ! print *, 'eps*(-U*t + x - x0)',eps*(-U*t + x - x0)
        ! print *, sqrt(eta*U/(1. - U)),sinh(eps*(-U*t + x - x0)),cosh(eps*(-U*t + x - x0)),sin(-x + t),cos(-x + t)


        ! ett = - 4*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**2*U**2*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
        !       - 4*sqrt(eta*U/(1. - U))*cos(-x + t)*eps*U*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !       + 2*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**2*U**2/cosh(eps*(-U*t + x - x0)) &
        !       + 2*sqrt(eta*U/(1. - U))*sin(-x + t)/cosh(eps*(-U*t + x - x0))

        ! ettt = - 12*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**3*U**3*sinh(eps*(-U*t + x - x0))**3/cosh(eps*(-U*t + x - x0))**4 &
        !        - 12*sqrt(eta*U/(1. - U))*cos(-x + t)*eps**2*U**2*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
        !        + 10*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**3*U**3*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !        + 6*sqrt(eta*U/(1. - U))*sin(-x + t)*eps*U*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !        + 6*sqrt(eta*U/(1. - U))*cos(-x + t)*eps**2*U**2/cosh(eps*(-U*t + x - x0)) &
        !        + 2*sqrt(eta*U/(1. - U))*cos(-x + t)/cosh(eps*(-U*t + x - x0))

        ! etttt = - 48*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**4*U**4*sinh(eps*(-U*t + x - x0))**4/cosh(eps*(-U*t + x - x0))**5 &
        !         - 48*sqrt(eta*U/(1. - U))*cos(-x + t)*eps**3*U**3*sinh(eps*(-U*t + x - x0))**3/cosh(eps*(-U*t + x - x0))**4 &
        !         + 56*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**4*U**4*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
        !         + 24*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**2*U**2*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
        !         + 40*sqrt(eta*U/(1. - U))*cos(-x + t)*eps**3*U**3*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !         - 10*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**4*U**4/cosh(eps*(-U*t + x - x0)) &
        !         +  8*sqrt(eta*U/(1. - U))*cos(-x + t)*eps*U*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !         - 12*sqrt(eta*U/(1. - U))*sin(-x + t)*eps**2*U**2/cosh(eps*(-U*t + x - x0)) &
        !         -  2*sqrt(eta*U/(1. - U))*sin(-x + t)/cosh(eps*(-U*t + x - x0))

        ! P
        p = 2.*tanh(eps*(x-x0-U*t))/cosh(eps*(x-x0-U*t))*cos(x-t)
        pt = - 2*eps*U*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
             + 2*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps*U*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
             - 2*tanh(eps*(-U*t + x - x0))*sin(-x + t)/cosh(eps*(-U*t + x - x0))
        pxx = -4.*eps**2*tanh(eps*(-U*t + x - x0))*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
              - 4.*eps**2*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
              + 4.*eps*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)/cosh(eps*(-U*t + x - x0)) &
              + 4.*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**2*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
              - 4.*tanh(eps*(-U*t + x - x0))*sin(-x + t)*eps*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
              - 2.*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**2/cosh(eps*(-U*t + x - x0)) &
              - 2.*tanh(eps*(-U*t + x - x0))*cos(-x + t)/cosh(eps*(-U*t + x - x0))
        ptxx = 2.*tanh(eps*(-U*t + x - x0))*sin(-x + t)*eps**2/cosh(eps*(-U*t + x - x0)) &
               - 12.*eps**3*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)*sinh(eps*(-U*t + x - x0))**2*U/cosh(eps*(-U*t + x - x0))**3 &
               + 8.*eps**2*tanh(eps*(-U*t + x - x0))*U*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)/cosh(eps*(-U*t + x - x0)) &
               + 8.*eps**2*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)*U*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
               + 12.*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**3*sinh(eps*(-U*t + x - x0))**3*U/cosh(eps*(-U*t + x - x0))**4 &
               - 10.*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**3*sinh(eps*(-U*t + x - x0))*U/cosh(eps*(-U*t + x - x0))**2 &
               - 8.*tanh(eps*(-U*t + x - x0))*sin(-x + t)*eps**2*sinh(eps*(-U*t + x - x0))**2*U/cosh(eps*(-U*t + x - x0))**3 &
               - 2.*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps*U*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
               - 8.*eps**3*tanh(eps*(-U*t + x - x0))**2*U*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
               + 2.*tanh(eps*(-U*t + x - x0))*sin(-x + t)/cosh(eps*(-U*t + x - x0)) &
               + 2.*eps*U*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
               + 4.*eps**2*tanh(eps*(-U*t + x - x0))*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)/cosh(eps*(-U*t + x - x0)) &
               + 4.*eps**2*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
               + 6.*eps**3*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)*U/cosh(eps*(-U*t + x - x0)) &
               - 4.*tanh(eps*(-U*t + x - x0))*sin(-x + t)*eps**2*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
               + 4.*tanh(eps*(-U*t + x - x0))*sin(-x + t)*eps**2*U/cosh(eps*(-U*t + x - x0)) &
               - 4.*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
               + 4.*eps**3*U*(1 - tanh(eps*(-U*t + x - x0))**2)**2*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
               - 12.*eps**3*tanh(eps*(-U*t + x - x0))*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)*U*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
               + 4.*eps*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)/cosh(eps*(-U*t + x - x0))

        ! ptt = - 4*eps**2*U**2*tanh(eps*(-U*t + x - x0))*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !       - 4*eps**2*U**2*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !       + 4*eps*U*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !       + 4*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**2*U**2*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
        !       - 4*tanh(eps*(-U*t + x - x0))*sin(-x + t)*eps*U*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !       - 2*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**2*U**2/cosh(eps*(-U*t + x - x0)) &
        !       - 2*tanh(eps*(-U*t + x - x0))*cos(-x + t)/cosh(eps*(-U*t + x - x0))

        ! pttt =  4*eps**3*U**3*(1 - tanh(eps*(-U*t + x - x0))**2)**2*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !       + 6*eps**3*U**3*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !       + 6*tanh(eps*(-U*t + x - x0))*sin(-x + t)*eps**2*U**2/cosh(eps*(-U*t + x - x0)) &
        !       + 2*tanh(eps*(-U*t + x - x0))*sin(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !       - 12*eps**3*U**3*tanh(eps*(-U*t + x - x0))*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !       + 6*eps*U*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !       - 6*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps*U*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !       - 8*eps**3*U**3*tanh(eps*(-U*t + x - x0))**2*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !       + 12*eps**2*U**2*tanh(eps*(-U*t + x - x0))*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !       - 12*eps**3*U**3*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
        !       + 12*eps**2*U**2*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !       + 12*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**3*U**3*sinh(eps*(-U*t + x - x0))**3/cosh(eps*(-U*t + x - x0))**4 &
        !       - 12*tanh(eps*(-U*t + x - x0))*sin(-x + t)*eps**2*U**2*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
        !       - 10*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**3*U**3*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2

        ! ptttt =   2*tanh(eps*(-U*t + x - x0))*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !         - 16*eps**3*U**3*(1 - tanh(eps*(-U*t + x - x0))**2)**2*sin(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !         - 24*eps**3*U**3*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !         + 24*eps**2*U**2*tanh(eps*(-U*t + x - x0))*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !         + 24*eps**2*U**2*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !         - 8*eps*U*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !         + 12*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**2*U**2/cosh(eps*(-U*t + x - x0)) &
        !         - 24*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**2*U**2*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
        !         + 8*tanh(eps*(-U*t + x - x0))*sin(-x + t)*eps*U*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !         + 10*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**4*U**4/cosh(eps*(-U*t + x - x0)) &
        !         - 32*eps**4*U**4*tanh(eps*(-U*t + x - x0))**2*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !         - 48*eps**4*U**4*tanh(eps*(-U*t + x - x0))*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
        !         + 48*eps**3*U**3*tanh(eps*(-U*t + x - x0))*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !         - 16*eps**4*U**4*tanh(eps*(-U*t + x - x0))**3*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !         + 32*eps**3*U**3*tanh(eps*(-U*t + x - x0))**2*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !         - 48*eps**4*U**4*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)*sinh(eps*(-U*t + x - x0))**3/cosh(eps*(-U*t + x - x0))**4 &
        !         + 48*eps**3*U**3*(1 - tanh(eps*(-U*t + x - x0))**2)*sin(-x + t)*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
        !         + 48*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**4*U**4*sinh(eps*(-U*t + x - x0))**4/cosh(eps*(-U*t + x - x0))**5 &
        !         - 48*tanh(eps*(-U*t + x - x0))*sin(-x + t)*eps**3*U**3*sinh(eps*(-U*t + x - x0))**3/cosh(eps*(-U*t + x - x0))**4 &
        !         - 56*tanh(eps*(-U*t + x - x0))*cos(-x + t)*eps**4*U**4*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**3 &
        !         + 32*eps**4*U**4*(1 - tanh(eps*(-U*t + x - x0))**2)**2*cos(-x + t)*tanh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0)) &
        !         + 16*eps**4*U**4*(1 - tanh(eps*(-U*t + x - x0))**2)**2*cos(-x + t)*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !         + 24*eps**4*U**4*tanh(eps*(-U*t + x - x0))*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)/cosh(eps*(-U*t + x - x0)) &
        !         + 40*eps**4*U**4*(1 - tanh(eps*(-U*t + x - x0))**2)*cos(-x + t)*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2 &
        !         + 40*tanh(eps*(-U*t + x - x0))*sin(-x + t)*eps**3*U**3*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**2
        ! pt = -2*eps*U*(1 - (tanh(eps*(-U*t + x - x0)))**2)/cosh(eps*(-U*t + x - x0))*cos(-x + t) &
        !     + 2*(tanh(eps*(-U*t + x - x0)))**2*eps*U/cosh(eps*(-U*t + x - x0))*cos(-x + t) &
        !     - 2*tanh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))*sin(-x + t)
        ! ptt = -10*eps**2*U**2*tanh(eps*(-U*t + x - x0))*(1 - (tanh(eps*(-U*t + x - x0)))**2)/cosh(eps*(-U*t + x - x0))*cos(-x + t) &
        !      + 4*eps*U*(1 - (tanh(eps*(-U*t + x - x0)))**2)/cosh(eps*(-U*t + x - x0))*sin(-x + t) &
        !      + 2*(tanh(eps*(-U*t + x - x0)))**3*eps**2*U**2/cosh(eps*(-U*t + x - x0))*cos(-x + t) &
        !      - 4*(tanh(eps*(-U*t + x - x0)))**2*eps*U/cosh(eps*(-U*t + x - x0))*sin(-x + t) &
        !      - 2*tanh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))*cos(-x + t)
        ! D
        d = 1-2./(cosh(eps*(x-x0-U*t)))**2
        dx = 4.*eps*sinh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))**3
        dxx = -12.*eps**2*sinh(eps*(-U*t + x - x0))**2/cosh(eps*(-U*t + x - x0))**4 + 4.*eps**2/cosh(eps*(-U*t + x - x0))**2

        ! derived from pde
        ptt = -b0*p-b1*pt+pnec*d*e
        ett = c**2*exx-alphaP*ptt
        dt = peptc*e*pt
        pttxx = -b0*pxx-b1*ptxx+pnec*(dxx*e+2.*dx*ex+d*exx)

        pttt = -b0*pt-b1*ptt+pnec*dt*e+pnec*d*et
        ettt = c**2*exxt-alphaP*pttt
        dtt = peptc*et*pt+peptc*e*ptt

        ptttt = -b0*ptt-b1*pttt+pnec*dtt*e+2.*pnec*dt*et+pnec*d*ett
        exxtt = c**2*exxxx-alphaP*pttxx
        ! print *, 'here1',exxtt, c, exxxx, alphaP,pttxx
        etttt = c**2*exxtt-alphaP*ptttt
        ! print *, 'here2',etttt, c, exxtt, alphaP,ptttt
    end subroutine getDerivs

    ! subroutine setupBCOld(dx,dt,t,ng,n,nP,nN,xvec,u,un,epsilon) ! set up bc for E
    !     implicit none
    !     integer ng,n,nP,nN,ig,tz,c
    !     real*8, allocatable :: xvec(:)
    !     ! real*8, allocatable :: u(:,:)
    !     real*8, allocatable :: En(:)
    !     real*8 t,tmp1,tmp2,dx,dt,xp,yp
    !     real*8 etmp1,etmp2
    !     real*8 ptmp1,ptmp2
    !     real*8 dtmp1,dtmp2,epsilon
    !     real*8 et,ett,ptt
    !     common tz
    !     real*8, allocatable :: u(:,:)
    !     real*8, allocatable :: un(:,:)

    !     ! speed of light
    !     c = 1.d0

    !     allocate(En(-ng:n+ng))

    !     En = un(-ng:n+ng,1)
    !     if (tz.eq.0) then ! initial boundary value problem
    !         ! left boundary: incomingWave
    !         call sechWave(etmp1,et,ett,ptmp1,pt,ptt,dtmp1,xvec(0),t,epsilon)
    !         En(0) = etmp1

    !         ! right boundary
    !         call sechWave(etmp2,et,ett,ptmp2,pt,ptt,dtmp2,xvec(n),t,epsilon)
    !         En(n) = etmp2

    !         ! extrapolation
    !         ! third order [3 -3 1]
    !         do ig=1,ng
    !             En(-ig) = 3.*En(-ig+1)-3.*En(-ig+2)+1.*En(-ig+3)
    !             En(n+ig) = 3.*En(n+ig-1)-3.*En(n+ig-2)+1.*En(n+ig-3)
    !         enddo
    !         ! fourth order accurate extrap coef = [4 -6 4 -1]
    !         ! do ig=1,ng
    !         !     En(-ig) = 4.*En(-ig+1)-6.*En(-ig+2)+4.*En(-ig+3)-1.*En(-ig+4)
    !         !     En(n+ig) = 4.*En(n+ig-1)-6.*En(n+ig-2)+4.*En(n+ig-3)-1.*En(n+ig-4)
    !         ! enddo

    !     elseif (tz.eq.1) then! twilight zone
    !         call getExactSolution(tmp1,xvec(0),t)
    !         call getExactSolution(tmp2,xvec(n),t)
    !         En(0) = tmp1
    !         En(n) = tmp2

    !         ! cheat ghost points, need to be replaced by extrapolation
    !         ! do ig=1,ng
    !         !     call getExactSolution(tmp1,xvec(-ig),t)
    !         !     En(-ig) = tmp1
    !         !     call getExactSolution(tmp2,xvec(n+ig),t)
    !         !     En(n+ig) = tmp2
    !         ! enddo
            
    !         ! extrapolation
    !         ! third order [3 -3 1]
    !         do ig=1,ng
    !             En(-ig) = 3.*En(-ig+1)-3.*En(-ig+2)+1.*En(-ig+3)
    !             En(n+ig) = 3.*En(n+ig-1)-3.*En(n+ig-2)+1.*En(n+ig-3)
    !         enddo
    !         ! fourth order accurate extrap coef = [4 -6 4 -1]
    !         ! do ig=1,ng
    !         !     En(-ig) = 4.*En(-ig+1)-6.*En(-ig+2)+4.*En(-ig+3)-1.*En(-ig+4)
    !         !     En(n+ig) = 4.*En(n+ig-1)-6.*En(n+ig-2)+4.*En(n+ig-3)-1.*En(n+ig-4)
    !         ! enddo
    !     endif

    !     un(-ng:n+ng,1) = En

    !     deallocate(En)
        
    ! end subroutine setupBCOld

    subroutine setupBC(dx,dt,t,ng,n,nP,nN,xvec,u,un,epsilon,iOrder,c,alphaP,b0,b1,pnec,peptc) ! set up bc for E
        implicit none
        integer ng,n,nP,nN,ig,tz,iOrder
        real*8, allocatable :: xvec(:)
        real*8, allocatable :: un(:,:),u(:,:)
        real*8, allocatable :: En(:),Pn(:),Dn(:)
        real*8 t,tmp1,tmp2,dx,dt,xp,yp
        real*8 etmp1,etmp2
        real*8 ptmp1,ptmp2
        real*8 dtmp1,dtmp2,epsilon
        real*8 et,ett,ettt,etttt,pt,ptt,pttt,ptttt
        real*8 c,alphaP,b0,b1,pnec,peptc
        common tz

        ! speed of light
        c = 1.d0

        allocate(En(-ng:n+ng))
        allocate(Pn(-ng:n+ng))
        allocate(Dn(-ng:n+ng))

        En = un(-ng:n+ng,1)
        Pn = un(-ng:n+ng,2)
        Dn = un(-ng:n+ng,3)
        if (tz.eq.0) then ! initial boundary value problem
            ! left boundary: incomingWave
            call sechWave(etmp1,et,ett,ettt,etttt,ptmp1,pt,ptt,pttt,ptttt,dtmp1,xvec(0),t,epsilon,c,alphaP,b0,b1,pnec,peptc)
            En(0) = etmp1

            ! right boundary
            call sechWave(etmp2,et,ett,ettt,etttt,ptmp2,pt,ptt,pttt,ptttt,dtmp2,xvec(n),t,epsilon,c,alphaP,b0,b1,pnec,peptc)
            En(n) = etmp2

            if (iOrder.eq.2) then
                ! third order [3 -3 1]
                do ig=1,ng
                    En(-ig) = 3.*En(-ig+1)-3.*En(-ig+2)+1.*En(-ig+3)
                    En(n+ig) = 3.*En(n+ig-1)-3.*En(n+ig-2)+1.*En(n+ig-3)
                enddo
            else
                ! 5th order [ 5   -10    10    -5     1]
                do ig=1,ng
                    En(-ig) = 5.*En(-ig+1)-10.*En(-ig+2)+10.*En(-ig+3)-5.*En(-ig+4)+En(-ig+5)
                    En(n+ig) = 5.*En(n+ig-1)-10.*En(n+ig-2)+10.*En(n+ig-3)-5.*En(n+ig-4)+En(n+ig-5)

                    Pn(-ig) = 5.*Pn(-ig+1)-10.*Pn(-ig+2)+10.*Pn(-ig+3)-5.*Pn(-ig+4)+Pn(-ig+5)
                    Pn(n+ig) = 5.*Pn(n+ig-1)-10.*Pn(n+ig-2)+10.*Pn(n+ig-3)-5.*Pn(n+ig-4)+Pn(n+ig-5)

                    Dn(-ig) = 5.*Dn(-ig+1)-10.*Dn(-ig+2)+10.*Dn(-ig+3)-5.*Dn(-ig+4)+Dn(-ig+5)
                    Dn(n+ig) = 5.*Dn(n+ig-1)-10.*Dn(n+ig-2)+10.*Dn(n+ig-3)-5.*Dn(n+ig-4)+Dn(n+ig-5)
                enddo
            endif

        elseif (tz.eq.1) then! twilight zone
            call getExactSolution(tmp1,xvec(0),t)
            call getExactSolution(tmp2,xvec(n),t)
            En(0) = tmp1
            En(n) = tmp2

            if (iOrder.eq.2) then
                ! third order [3 -3 1]
                do ig=1,ng
                    En(-ig) = 3.*En(-ig+1)-3.*En(-ig+2)+1.*En(-ig+3)
                    En(n+ig) = 3.*En(n+ig-1)-3.*En(n+ig-2)+1.*En(n+ig-3)
                enddo
            else
                ! 5th order [ 5   -10    10    -5     1]
                do ig=1,ng
                    En(-ig) = 5.*En(-ig+1)-10.*En(-ig+2)+10.*En(-ig+3)-5.*En(-ig+4)+En(-ig+5)
                    En(n+ig) = 5.*En(n+ig-1)-10.*En(n+ig-2)+10.*En(n+ig-3)-5.*En(n+ig-4)+En(n+ig-5)

                    Pn(-ig) = 5.*Pn(-ig+1)-10.*Pn(-ig+2)+10.*Pn(-ig+3)-5.*Pn(-ig+4)+Pn(-ig+5)
                    Pn(n+ig) = 5.*Pn(n+ig-1)-10.*Pn(n+ig-2)+10.*Pn(n+ig-3)-5.*Pn(n+ig-4)+Pn(n+ig-5)

                    Dn(-ig) = 5.*Dn(-ig+1)-10.*Dn(-ig+2)+10.*Dn(-ig+3)-5.*Dn(-ig+4)+Dn(-ig+5)
                    Dn(n+ig) = 5.*Dn(n+ig-1)-10.*Dn(n+ig-2)+10.*Dn(n+ig-3)-5.*Dn(n+ig-4)+Dn(n+ig-5)
                enddo

            endif

        endif

        un(-ng:n+ng,1) = En
        un(-ng:n+ng,2) = Pn
        un(-ng:n+ng,3) = Dn

        deallocate(En,Pn,Dn)
        
    end subroutine setupBC

    subroutine interp1d(x,y,xp,yp) ! Lagrange interpolation
        implicit none
        real*8 xp,p,yp
        real*8, allocatable :: x(:),y(:)
        integer n,i,j

        n = size(x,1)

        ! print *, 'size of interpolation',n

        do i=1,n
            p=1.d0
            do j=1,n
                if (i.ne.j) then
                    p=p*(xp-x(j))/(x(i)-x(j))
                endif
            enddo
            yp=yp+p*y(i)
        enddo

    end subroutine interp1d

    subroutine getError(n,errE,errP,errN,t,xvec,ng,nP,nN,un) ! compute twilight zone errors
        implicit none
        real*8 t,errE,errP,errN,etmp,u
        integer ng,np,nN,i,k,n
        real*8, allocatable :: xvec(:)
        ! real*8, allocatable :: un(:,:)
        ! real*8 un(-1:n+1,3)
        real*8, allocatable :: un(:,:)

        errE = 0.
        errP = 0.
        errN = 0.
        do i=0,n
            call getExactSolution(u,xvec(i),t)
            etmp = abs(u-un(i,1))
            errE = max(errE,etmp)
            do k=1,nP
                etmp = abs(u-un(i,1+k))
                errP = max(errP,etmp)
            enddo
            do k=1,nN
                etmp = abs(u-un(i,1+nP+k))
                errN = max(errN,etmp)
            enddo
        enddo


    end subroutine getError

    ! subroutine getExactSolutionOld(u,x,t) ! exact degree 2 polynomials
    !     implicit none
    !     real*8 u,x,t
    !     real*8 a0, a1, a2, c0, c1, c2

    !     a0 = 1.d0
    !     a1 = 1.d0
    !     a2 = 1.d0
        
    !     c0 = 1.d0
    !     c1 = 1.d0
    !     c2 = 1.d0
        

    !     u = (a2*x**2+a1*x+a0)*(c2*t**2+c1*t+c0)

    ! end subroutine getExactSolutionOld

    subroutine getExactSolution(u,x,t) ! degree 4 polynomials
        implicit none
        real*8 u,x,t
        real*8 a0, a1, a2, a3, a4, c0, c1, c2,c3, c4

        a0 = 1.d0
        a1 = 1.d0
        a2 = 1.d0
        a3 = 1.d0
        a4 = 1.d0
        
        c0 = 1.d0
        c1 = 1.d0
        c2 = 1.d0
        c3 = 0.d0
        c4 = 0.d0
        

        u = (a4*x**4+a3*x**3+a2*x**2+a1*x+a0)*(c4*t**4+c3*t**3+c2*t**2+c1*t+c0)

    end subroutine getExactSolution

    subroutine timeStep1D(pnec,peptc,ng,n,alphaP,dx,dt,xvec,t,b0,b1,c,um,u,un,iOrder)
        implicit none
        real*8 dx,dt,t,c,Psum,alphaP,x,tmp1,tmp2
        integer iOrder,ng,n,i,j,k,tz,ig
        real*8, allocatable :: xvec(:)
        real*8, allocatable :: u(:,:),un(:,:)
        real*8, allocatable :: Em(:),E(:),En(:)
        real*8, allocatable :: Pm(:),P(:),Pn(:)
        real*8, allocatable :: Dm(:),D(:),Dn(:)
        real*8, allocatable :: fe(:),fet(:),fett(:),fexx(:)
        real*8, allocatable :: fp(:),fpt(:),fptt(:),fpxx(:)
        real*8, allocatable :: fn(:),fnt(:),fntt(:),fnttt(:)
        real*8 b0,b1,beta,pnec,peptc
        real*8 Pt,Ptt,Pttt,Ptttt,Et,Ett,Ettt,Nt,Dtt,Dttt,Dtttt ! Nt=Dt to avoid dt
        real*8 Pxx,Pxxn,Pxxm,Pxxtt,Dxx,Nx,Exx,Ex ! Nx=Dx to avoid dx
        real*8 elap4,elap2,elap2m,elap2n,elapsq
        real*8 um(-2:n+2,3)

        common tz

        allocate(Em(-ng:n+ng))
        allocate(E(-ng:n+ng))
        allocate(En(-ng:n+ng))
        allocate(Pm(-ng:n+ng))
        allocate(P(-ng:n+ng))
        allocate(Pn(-ng:n+ng))
        allocate(Dm(-ng:n+ng))
        allocate(D(-ng:n+ng))
        allocate(Dn(-ng:n+ng))

        beta = 1./(1.+0.5*b1*dt)

        allocate(fe(-ng:n+ng))
        allocate(fet(-ng:n+ng))
        allocate(fett(-ng:n+ng))
        allocate(fexx(-ng:n+ng))
        allocate(fp(-ng:n+ng))
        allocate(fpt(-ng:n+ng))
        allocate(fptt(-ng:n+ng))
        allocate(fpxx(-ng:n+ng))
        allocate(fn(-ng:n+ng))
        allocate(fnt(-ng:n+ng))
        allocate(fntt(-ng:n+ng))
        allocate(fnttt(-ng:n+ng))

        ! forcing functions
        if (tz==1) then
            ! call getForcing(b0,b1,n,ng,nP,nN,fe,fp,fn,fnt,xvec,t-dt,pnec,peptc,c,alphaP)
            call getForcing(b0,b1,n,ng,nP,nN,fe,fet,fett,fexx,fp,fpt,fptt,fpxx,fn,fnt,fntt,fnttt,xvec,t-dt,pnec,peptc,c,alphaP)
        else ! not twilight zone
            do i=0,n
                fe(i) = 0.
                fet(i) = 0.
                fett(i) = 0.
                fexx(i) = 0.
                fp(i) = 0.
                fpt(i) = 0.
                fptt(i) = 0.
                fpxx(i) = 0.
                fn(i) = 0.
                fnt(i) = 0.
                fntt(i) = 0.
                fnttt(i) = 0.
            enddo   
        endif

        Em = um(-ng:n+ng,1)
        E  = u(-ng:n+ng,1)
        Pm = um(-ng:n+ng,2)
        P  = u(-ng:n+ng,2)
        Dm = um(-ng:n+ng,3)
        D  = u(-ng:n+ng,3)

        !-------------------------------------
        ! second order update
        !-------------------------------------
        ig = 0
        ! if (iOrder.eq.2) then
        !     ig = 0
        ! else
        !     ig = 1
        ! endif

        do i=0-ig,n+ig ! loop over space
            x = xvec(i)
            ! second order update
            ! P
            Pn(i) = 2.*P(i)-Pm(i)+0.5*dt*b1*Pm(i)-dt**2*b0*P(i)+dt**2*fp(i)
            Pn(i) = Pn(i)+dt**2*pnec*D(i)*E(i)
            Pn(i) = beta*Pn(i) ! adjust
            Psum = Pn(i)-2.*P(i)+Pm(i)
            un(i,2) = Pn(i)
        
            ! E
            if (iOrder.eq.2)then
                Exx = (E(i+1)-2.*E(i)+E(i-1))/dx**2
            else
                Exx = (-E(i+2)/12.+4.*E(i+1)/3. -5.*E(i)/2.+4.*E(i-1)/3.-E(i-2)/12.)/(dx**2)
            endif
            En(i) = 2.*E(i)-Em(i) + dt**2*c**2*Exx-alphaP*Psum+dt**2*fe(i)
            un(i,1) = En(i)

            ! Nt
            Pt = (Pn(i)-Pm(i))/(2.*dt)
            Nt = peptc*E(i)*Pt + fn(i)

            ! Ntt
            Et = (En(i)-Em(i))/(2.*dt)
            Ptt = (Pn(i)-2.*P(i)+Pm(i))/(dt**2)
            Dtt =  peptc*Et*Pt + peptc*E(i)*Ptt + fnt(i)

            ! N
            Dn(i) = D(i)+dt*Nt+dt**2/2.*Dtt
            un(i,3) = Dn(i)
        enddo
        !-------------------------------------
        ! fourth order update
        !-------------------------------------
        if (iOrder.eq.4) then
            ! extrapolate

            ! 
            do i=0,n
                Pt = (Pn(i)-Pm(i))/(2.*dt)
                ! Nt = (Dn(i)-Dm(i))/(2.*dt)
                Nt = peptc*E(i)*Pt + fn(i)
                Et = (En(i)-Em(i))/(2.*dt)
                ! Ptt = (Pn(i)-2*P(i)+Pm(i))/(dt**2)
                Ptt = -b1*Pt - b0*P(i) + pnec*D(i)*E(i) + fp(i)
                ! Dtt = (Dn(i)-2*D(i)+Dm(i))/(dt**2)
                Dtt =  peptc*Et*Pt + peptc*E(i)*Ptt + fnt(i)
                Ett = (En(i)-2*E(i)+Em(i))/(dt**2)
                ! Exx = (E(i+1)-2.*E(i)+E(i-1))/dx**2
                ! Ett = c**2*Exx-alphaP*Ptt + fe(i)
                Pttt = -b1*Ptt -b0*Pt + pnec*(D(i)*Et+Nt*E(i)) + fpt(i)
                Ptttt = -b1*Pttt-b0*Ptt+pnec*(D(i)*Ett+2.*Nt*Et+Dtt*E(i)) + fptt(i)
                ! print *, 'Et,Ett',Et,Ett
                ! print * , 'Pt,Ptt',Pt,Ptt,Pttt,Ptttt
                ! print *, 'Dt,Dtt',Nt,Dtt
                ! print *, 'E',En(i),E(i),Em(i)
                ! print *, 'P',Pn(i),P(i),Pm(i)
                ! print *, 'D',Dn(i),D(i),Dm(i)

                ! P
                Pn(i) = 2.*P(i) - Pm(i) + 0.5*dt*b1*Pm(i) - dt**2*b0*P(i) + dt**2*fp(i) + dt**4/12.*Ptttt + dt**4/6.*b1*Pttt
                Pn(i) = Pn(i)+dt**2*pnec*D(i)*E(i)
                Pn(i) = beta*Pn(i) ! adjust
                Psum = Pn(i)-2.*P(i)+Pm(i)
                un(i,2) = Pn(i)

                ! Pxx = (P(i-1) - 2.*P(i) + P(i+1))/dx**2
                ! Exx = (E(i+1)-2.*E(i)+E(i-1))/dx**2
                ! Dxx = (D(i+1)-2.*D(i)+D(i-1))/dx**2
                ! Ex = (E(i+1)-E(i-1))/(2.*dx)
                ! Nx = (D(i+1)-D(i-1))/(2.*dx)
                ! Pxxm = (Pm(i-1) - 2.*Pm(i) + Pm(i+1))/dx**2
                Pxxm = (-Pm(i+2)/12.+4.*Pm(i+1)/3. -5.*Pm(i)/2.+4.*Pm(i-1)/3.-Pm(i-2)/12.)/(dx**2)
                Pxx = (-P(i+2)/12.+4.*P(i+1)/3. -5.*P(i)/2.+4.*P(i-1)/3.-P(i-2)/12.)/(dx**2)
                Exx = (-E(i+2)/12.+4.*E(i+1)/3. -5.*E(i)/2.+4.*E(i-1)/3.-E(i-2)/12.)/(dx**2)
                Dxx = (-D(i+2)/12.+4.*D(i+1)/3. -5.*D(i)/2.+4.*D(i-1)/3.-D(i-2)/12.)/(dx**2)
                Nx = (-D(i+2)+8.*D(i+1)-8.*D(i-1)+D(i-2))/(12*dx)
                Ex = (-E(i+2)+8.*E(i+1)-8.*E(i-1)+E(i-2))/(12*dx)
                Pxxn = 2.*Pxx-Pxxm+0.5*dt*b1*Pxxm-dt**2*b0*Pxx+dt**2*fpxx(i)
                Pxxn = Pxxn + pnec*dt**2*(Dxx*E(i)+2.*Nx*Ex+D(i)*Exx)
                Pxxn = Pxxn*beta
                ! print *, 'Pxxn',Pxxn-(Pn(i-1) - 2.*Pn(i) + Pn(i+1))/dx**2
                ! Pxxn = (Pn(i-1) - 2.*Pn(i) + Pn(i+1))/dx**2
                Pxxtt = (Pxxn-2.*Pxx+Pxxm)/dt**2
                ! print *, 'Pxxtt',Pxxtt

                ! E
                elap4 = (-E(i+2)/12.+4.*E(i+1)/3. -5.*E(i)/2.+4.*E(i-1)/3.-E(i-2)/12.)/(dx**2)
                elapsq = (E(i+2)-4.*E(i+1)+6.*E(i)-4.*E(i-1)+E(i-2))/(dx**4)
                En(i) = 2.*E(i)-Em(i)+dt**2*c**2*elap4-alphaP*Psum+dt**2*fe(i)&
                +dt**4/12.*(c**4*elapsq-c**2*alphaP*Pxxtt+c**2*fexx(i)+fett(i))

                ! E
                ! elap4 = (-E(i+2)/12.+4.*E(i+1)/3.-5.*E(i)/2.+4.*E(i-1)/3.-E(i-2)/12.)/dx**2
                ! ! print *, 'elap4',elap4,Exx,elap4-Exx
                ! elapsq = (E(i+2)-4.*E(i+1)+6.*E(i)-4.*E(i-1)+E(i-2))/dx**4
                ! ! print *, 'elapsq',elapsq
                ! En(i) = 2.*E(i)-Em(i) + dt**2*c**2*elap4 - alphaP*Psum + dt**2*fe(i) &
                !          + dt**4/12.*(c**4*elapsq - c**2*alphaP*Pxxtt+c**2*fexx(i)+fett(i))
                ! print *, 'Etttt1',(c**4*elapsq - c**2*alphaP*Pxxtt+c**2*fexx(i)+fett(i))
                ! !En(i) = 2.*E(i)-Em(i) + dt**2*c**2*elap4 - alphaP*Psum + dt**2*fe(i) &
                ! !        + dt**4/12.*(c**4*elapsq - c**2*alphaP*Pxxtt + c**2*fexx(i) + fett(i))

                ! elap2m = (Em(i+1)-2.*Em(i)+Em(i-1))/dx**2
                ! elap2 = (E(i+1)-2.*E(i)+E(i-1))/dx**2
                elap2m = (-Em(i+2)/12.+4.*Em(i+1)/3. -5.*Em(i)/2.+4.*Em(i-1)/3.-Em(i-2)/12.)/(dx**2)
                elap2 = (-E(i+2)/12.+4.*E(i+1)/3. -5.*E(i)/2.+4.*E(i-1)/3.-E(i-2)/12.)/(dx**2)
                elap2n = 2.*elap2-elap2m + dt**2*c**2*elapsq-alphaP*Pxxtt*dt**2+dt**2*fexx(i)
                Ettt = c**2*(elap2n-elap2m)/(2*dt)-alphaP*Pttt + fet(i)
                un(i,1) = En(i)

                ! D
                Pt = (Pn(i)-Pm(i))/(2.*dt) - dt**2/6.*Pttt
                Ptt = -b1*Pt - b0*P(i) + pnec*D(i)*E(i) + fp(i)
                Pttt = -b1*Ptt - b0*Pt + pnec*D(i)*Et + pnec*Nt*E(i) + fpt(i)
                Ptttt = -b1*Pttt - b0*Ptt + pnec*D(i)*Ett + 2.*pnec*Nt*Et + pnec*Dtt*E(i) + fptt(i)
                ! Ptt = (Pn(i)-2*P(i)+Pm(i))/(dt**2) - dt**2/12.*Ptttt
                Et = (En(i)-Em(i))/(2.*dt) - dt**2/6.*Ettt

                Nt = peptc*E(i)*Pt + fn(i)
                Dtt = peptc*(Et*Pt + E(i)*Ptt) + fnt(i)
                Dttt = peptc*(Ett*Pt + 2.*Et*Ptt+ E(i)*Pttt) + fntt(i)
                Dtttt = peptc*(Ettt*Pt + 3.*Ett*Ptt + 3.*Et*Pttt+ E(i)*Ptttt) + fnttt(i)

                ! print *, t,Nt,Dtt,Dttt,Dtttt

                ! N
                Dn(i) = D(i)+dt*Nt+dt**2/2.*Dtt+dt**3/6.*Dttt+dt**4/24.*Dtttt
                un(i,3) = Dn(i)


            enddo
        endif

        deallocate(Em,E,En,Pm,P,Pn,Dm,D,Dn,fe,fet,fett,fexx,fp,fpt,fptt,fpxx,fn,fnt,fntt,fnttt)

    end subroutine timeStep1D

    ! subroutine timeStep1DOld(pnec,peptc,ng,n,alphaP,dx,dt,xvec,t,b0,b1,c,um,u,un,iOrder)
    !     implicit none
    !     real*8 dx,dt,t,c,Psum,alphaP,x,tmp1,tmp2
    !     integer iOrder,ng,n,i,j,k,tz
    !     real*8, allocatable :: xvec(:)
    !     ! real*8, allocatable :: um(:,:),u(:,:)!,un(:,:)
    !     real*8, allocatable :: Em(:),E(:),En(:)
    !     real*8, allocatable :: Pm(:),P(:),Pn(:)
    !     real*8, allocatable :: Dm(:),D(:),Dn(:)
    !     real*8, allocatable :: fe(:)
    !     real*8, allocatable :: fp(:)
    !     real*8, allocatable :: fn(:),fnt(:)
    !     real*8 b0,b1,beta,pnec,peptc,Nt,Ntt
    !     ! real*8 un(-1:n+1,3)
    !     real*8 um(-2:n+2,3)!,u(-1:n+1,3)
    !     real*8, allocatable :: u(:,:)
    !     real*8, allocatable :: un(:,:)

    !     common tz

    !     allocate(Em(-ng:n+ng))
    !     allocate(E(-ng:n+ng))
    !     allocate(En(-ng:n+ng))
    !     allocate(Pm(-ng:n+ng))
    !     allocate(P(-ng:n+ng))
    !     allocate(Pn(-ng:n+ng))
    !     allocate(Dm(-ng:n+ng))
    !     allocate(D(-ng:n+ng))
    !     allocate(Dn(-ng:n+ng))

    !     beta = 1./(1.+0.5*b1*dt)

    !     allocate(fe(-ng:n+ng))
    !     allocate(fp(-ng:n+ng))
    !     allocate(fn(-ng:n+ng))
    !     allocate(fnt(-ng:n+ng))

    !     ! forcing functions
    !     if (tz==1) then
    !         call getForcing(b0,b1,n,ng,nP,nN,fe,fp,fn,fnt,xvec,t-dt,pnec,peptc,c,alphaP)
    !     else ! not twilight zone
    !         do i=0,n
    !             fe(i) = 0.
    !             fp(i) = 0.
    !             fn(i) = 0.
    !             fnt(i) = 0.
    !         enddo   
    !     endif

    !     Em = um(-ng:n+ng,1)
    !     E  = u(-ng:n+ng,1)
    !     Pm = um(-ng:n+ng,2)
    !     P  = u(-ng:n+ng,2)
    !     Dm = um(-ng:n+ng,3)
    !     D  = u(-ng:n+ng,3)

    !     do i=0,n ! loop over space
    !         x = xvec(i)
    !         ! second order update
    !         ! P
    !         Pn(i) = 2.*P(i)-Pm(i)+0.5*dt*b1*Pm(i)-dt**2*b0*P(i)+dt**2*fp(i)
    !         Pn(i) = Pn(i)+dt**2*pnec*D(i)*E(i)
    !         Pn(i) = beta*Pn(i) ! adjust
    !         Psum = Pn(i)-2.*P(i)+Pm(i)
    !         un(i,2) = Pn(i)
        
    !         ! E
    !         En(i) = 2.*E(i)-Em(i) + dt**2*c**2*(E(i+1)-2.*E(i)+E(i-1))/dx**2-alphaP*Psum+dt**2*fe(i)
    !         un(i,1) = En(i)

    !         ! Nt
    !         Nt = peptc*E(i)*(Pn(i)-Pm(i))/(2.*dt) + fn(i)

    !         ! Ntt
    !         Ntt =  peptc*(En(i)-Em(i))/(2.*dt)*(Pn(i)-Pm(i))/(2.*dt) &
    !              + peptc*E(i)*(Pn(i)-2.*P(i)+Pm(i))/(dt**2) &
    !              + fnt(i)

    !         ! N
    !         Dn(i) = D(i)+dt*Nt+dt**2/2.*Ntt
    !         un(i,3) = Dn(i)

    !     enddo

    !     ! fourth order update
    !     ! if (iOrder.eq.4) then

    !     ! endif

    !     deallocate(Em,E,En,Pm,P,Pn,Dm,D,Dn,fe,fp,fn,fnt)

    ! end subroutine timeStep1DOld

    subroutine getForcing(b0,b1,n,ng,nP,nN,fe,fet,fett,fexx,fp,fpt,fptt,fpxx,fn,fnt,fntt,fnttt,xvec,t,pnec,peptc,c,alphaP)
        implicit none

        integer nP,nN,ng,n,j,k,i
        real*8, allocatable :: xvec(:)
        real*8, allocatable :: fe(:),fet(:),fett(:),fexx(:)
        real*8, allocatable :: fp(:),fpt(:),fptt(:),fpxx(:)
        real*8, allocatable :: fn(:),fnt(:),fntt(:),fnttt(:)
        real*8 t,x,u,ux,uxx,uxxxx,ut,utt,uttt,utttt,utxx,uttxx,c,alphaP
        real*8 a0, a1, a2, a3, a4, c0, c1, c2,c3, c4
        real*8 b0,b1,beta,pnec,peptc

        a0 = 1.d0
        a1 = 1.d0
        a2 = 1.d0
        a3 = 1.d0
        a4 = 1.d0
        
        c0 = 1.d0
        c1 = 1.d0
        c2 = 1.d0
        c3 = 0.d0
        c4 = 0.d0
        

        ! u = (a4*x**4+a3*x**3+a2*x**2+a1*x+a0)*(c4*t**4+c3*t**3+c2*t**2+c1*t+c0)

        do i=-ng,n+ng
            x   = xvec(i)

            u = (a4*x**4+a3*x**3+a2*x**2+a1*x+a0)*(c4*t**4+c3*t**3+c2*t**2+c1*t+c0)

            ux  = (4*a4*x**3+3*a3*x**2+2*a2*x+a1)*(c4*t**4+c3*t**3+c2*t**2+c1*t+c0)
            uxx = (4*3*a4*x**2+3*2*a3*x+2*a2)*(c4*t**4+c3*t**3+c2*t**2+c1*t+c0)
            uxxxx = (4*3*2*a4)*(c4*t**4+c3*t**3+c2*t**2+c1*t+c0)

            ut    = (a4*x**4+a3*x**3+a2*x**2+a1*x+a0)*(4*c4*t**3+3*c3*t**2+2*c2*t+c1)
            ! utx    = (4*a4*x**3+3*a3*x**2+2*a2*x+a1)*(4*c4*t**3+3*c3*t**2+2*c2*t+c1)
            utxx    = (4*3*a4*x**2+3*2*a3*x+2*a2)*(4*c4*t**3+3*c3*t**2+2*c2*t+c1)
            uttxx    = (4*3*a4*x**2+3*2*a3*x+2*a2)*(4*3*c4*t**2+3*2*c3*t+2*c2)

            utt   = (a4*x**4+a3*x**3+a2*x**2+a1*x+a0)*(4*3*c4*t**2+3*2*c3*t+2*c2)
            uttt  = (a4*x**4+a3*x**3+a2*x**2+a1*x+a0)*(4*3*2*c4*t+3*2*c3)
            utttt = (a4*x**4+a3*x**3+a2*x**2+a1*x+a0)*(4*3*2*c4)

            ! fe
            fe(i) = utt-c**2*uxx+alphaP*utt
            fet(i)  = uttt-c**2*utxx+alphaP*uttt
            fett(i) = utttt-c**2*uttxx+alphaP*utttt
            fexx(i) = uttxx-c**2*uxxxx+alphaP*uttxx

            ! fp
            fp(i) = utt+b1*ut+b0*u-pnec*u*u
            fpt(i) = uttt+b1*utt+b0*ut-pnec*(u*ut+ut*u)
            fptt(i) = utttt+b1*uttt+b0*utt-pnec*(u*utt+2.*ut*ut+utt*u)
            fpxx(i) = uttxx+b1*utxx+b0*uxx-pnec*(u*uxx+2.*ux*ux+uxx*u)

            ! fn
            fn(i)    = ut-peptc*u*ut
            fnt(i)   = utt-peptc*(ut*ut+u*utt)
            fntt(i)  = uttt-peptc*(utt*ut+2.*ut*utt+u*uttt)
            fnttt(i) = utttt-peptc*(uttt*ut+3.*utt*utt+3.*ut*uttt+u*utttt)
        enddo

    end subroutine getForcing

    ! subroutine getForcingOld(b0,b1,n,ng,nP,nN,fe,fp,fn,fnt,xvec,t,pnec,peptc,c,alphaP)
    !     implicit none

    !     integer nP,nN,ng,n,j,k,i
    !     real*8, allocatable :: xvec(:)
    !     real*8, allocatable :: fe(:)
    !     real*8, allocatable :: fp(:)
    !     real*8, allocatable :: fn(:),fnt(:)
    !     real*8 t,x,u,ux,uxx,ut,utt,c,alphaP
    !     real*8 a0, a1, a2, c0, c1, c2
    !     real*8 b0,b1,beta,pnec,peptc

    !     a0 = 1.d0
    !     a1 = 1.d0
    !     a2 = 1.d0

    !     c0 = 1.d0
    !     c1 = 1.d0
    !     c2 = 1.d0

    !     do i=0,n
    !         x   = xvec(i)
    !         u   = (a2*x**2+a1*x+a0)*(c2*t**2+c1*t+c0)
    !         ux  = (2.*a2*x+a1)*(c2*t**2+c1*t+c0)
    !         uxx = (2.*(a2))*(c2*t**2+c1*t+c0)
    !         ut  = (a2*x**2+a1*x+a0)*(2.*c2*t+c1)
    !         utt = (a2*x**2+a1*x+a0)*2.*c2

    !         ! fe
    !         fe(i) = utt-c**2*uxx+alphaP*utt

    !         ! fp
    !         fp(i) = utt+b1*ut+b0*u-pnec*u*u

    !         ! fn
    !         fn(i) = ut-peptc*u*ut
    !         fnt(i) = utt-peptc*ut*ut-peptc*u*utt
    !     enddo

    ! end subroutine getForcingOld

end subroutine mbe1d













