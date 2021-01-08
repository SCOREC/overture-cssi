subroutine mbe1d(epsilon,dt,xa,xb,n,tfinal,uout)
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
    real*8 um(-1:n+1,3),uout(-1:n+1,3)
    real*8, save, allocatable :: un(:,:)
    real*8, save, allocatable :: u(:,:)
    real*8, save :: tcurrent
    
    real*8, allocatable :: xvec(:)
    ! real*8, allocatable :: um(:,:),u(:,:)!,un(:,:)
    real*8 b0,b1,pnec,peptc,etmp,ptmp,dtmp,et,ett,pt,ptt

    common tz

    if (.not.allocated(u)) allocate(u(-1:n+1,3))
    if (.not.allocated(un)) allocate(un(-1:n+1,3))

    ! output file or not
    writeFile = 0

    ! call getarg(1,arg)
    ! read(arg,*) n
    ! n = 2**n-1

    iOrder=2 ! currently only support second order version
    ng = iOrder/2

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
    print *, "time step 1", dt
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
    print *, "Fortran 1D code, tfinal diff 1", Nt*dt-tfinal
    if (abs(Nt*dt-tfinal)>1.d-10) Nt = floor(abs(tfinal)/dt)
    print *, "Fortran 1D code, tfinal diff 2", Nt*dt-tfinal
    print *, "time step 2", dt

    ! print *, dx, dt, Nt, tfinal
    print *, tfinal, dt,'Number of time steps is:', Nt
    
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
        call setupInitials(n,ng,nP,nN,xvec,u,um,dt,dx,epsilon)
        t = 0.
        un = u ! value at t=0.
        u = um

        uout = un
        do i=-ng,n+ng
            uout(i,2) = epsilon*un(i,2)
        enddo

        !
        tcurrent = t
        print *, 'Case 1',t

    elseif (Nt.eq.1 .and. dt<0.) then
        t = dt
        call setupInitials(n,ng,nP,nN,xvec,u,um,-dt,dx,epsilon)
        un = um ! value at t=0.

        uout = un
        do i=-ng,n+ng
            uout(i,2) = epsilon*un(i,2)
        enddo

        !
        tcurrent = t
        print *, 'Case 2',t

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

        print *, 'Case 3',tcurrent, tfinal

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

            call setupBC(dx,dt,t,ng,n,nP,nN,xvec,u,un,epsilon)

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

        print *, 'Case 4',tcurrent,tfinal

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
            call sechWave(etmp,et,ett,ptmp,pt,ptt,dtmp,xvec(i),t,epsilon)
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

    print *, 'Max errors of (E,P,D) in space and time at time:', t,Nt
    print *, eE,eP,eN

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

    subroutine setupInitials(n,ng,nP,nN,xvec,u,um,dt,dx,epsilon)
        implicit none
        integer i,j,ng,nP,nN,n,tz
        real*8, allocatable :: xvec(:)
        ! real*8, allocatable :: um(:,:),u(:,:)
        real*8 um(-1:n+1,3),u(-1:n+1,3)
        real*8 dt,tmp1,tmp2,dx,epsilon
        real*8 etmp1,ptmp1,dtmp1,etmp2,ptmp2,dtmp2,et,ett,pt,ptt
        common tz
           
        do i=-ng,n+ng
            
            if (tz.eq.0) then 
                call sechWave(etmp1,et,ett,ptmp1,pt,ptt,dtmp1,xvec(i),-dt,epsilon) ! filling initial conditions
                call sechWave(etmp2,et,ett,ptmp2,pt,ptt,dtmp2,xvec(i),0.d0,epsilon)
                um(i,1) = etmp1
                ! um(i,1) = etmp2-et*dt+dt**2/2.*ett
                u(i,1)  = etmp2
                um(i,2) = ptmp1
                ! um(i,2) = ptmp2-pt*dt+dt**2/2.*ptt
                u(i,2)  = ptmp2
                um(i,3) = dtmp1
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

    subroutine sechWave(e,et,ett,p,pt,ptt,d,x,t,epsilon) ! incoming wave
        implicit none
        real*8 x,t,eta,eps,U,x0,e,p,d,epsilon,et,ett,pt,ptt

        eta = 1.
        U = 1./2.
        x0 = 0.
        eps = epsilon/2.*sqrt(eta/(U-U*U))
        ! E
        e = 2.*sqrt(eta*U/(1.-U))/cosh(eps*(x-x0-U*t))*sin(x-t)
        et = -2*sqrt(eta*U/(1 - U))*eps*U/cosh(eps*(-U*t + x - x0))*tanh(eps*(-U*t + x - x0))*sin(-x + t) & 
            - 2*sqrt(eta*U/(1 - U))/cosh(eps*(-U*t + x - x0))*cos(-x + t)
        ett = -2*sqrt(eta*U/(1 - U))*eps**2*U**2/cosh(eps*(-U*t + x - x0))*(tanh(eps*(-U*t + x - x0)))**2*sin(-x + t) &
             + 2*sqrt(eta*U/(1 - U))*eps**2*U**2/cosh(eps*(-U*t + x - x0))*(1 - (tanh(eps*(-U*t + x - x0)))**2)*sin(-x + t) &
             - 4*sqrt(eta*U/(1 - U))*eps*U/cosh(eps*(-U*t + x - x0))*tanh(eps*(-U*t + x - x0))*cos(-x + t) &
             + 2*sqrt(eta*U/(1 - U))/cosh(eps*(-U*t + x - x0))*sin(-x + t)
        ! P
        p = 2.*tanh(eps*(x-x0-U*t))/cosh(eps*(x-x0-U*t))*cos(x-t)
        pt = -2*eps*U*(1 - (tanh(eps*(-U*t + x - x0)))**2)/cosh(eps*(-U*t + x - x0))*cos(-x + t) &
            + 2*(tanh(eps*(-U*t + x - x0)))**2*eps*U/cosh(eps*(-U*t + x - x0))*cos(-x + t) &
            - 2*tanh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))*sin(-x + t)
        ptt = -10*eps**2*U**2*tanh(eps*(-U*t + x - x0))*(1 - (tanh(eps*(-U*t + x - x0)))**2)/cosh(eps*(-U*t + x - x0))*cos(-x + t) &
             + 4*eps*U*(1 - (tanh(eps*(-U*t + x - x0)))**2)/cosh(eps*(-U*t + x - x0))*sin(-x + t) &
             + 2*(tanh(eps*(-U*t + x - x0)))**3*eps**2*U**2/cosh(eps*(-U*t + x - x0))*cos(-x + t) &
             - 4*(tanh(eps*(-U*t + x - x0)))**2*eps*U/cosh(eps*(-U*t + x - x0))*sin(-x + t) &
             - 2*tanh(eps*(-U*t + x - x0))/cosh(eps*(-U*t + x - x0))*cos(-x + t)
        ! D
        d = 1-2./(cosh(eps*(x-x0-U*t)))**2

    end subroutine sechWave

    subroutine setupBC(dx,dt,t,ng,n,nP,nN,xvec,u,un,epsilon) ! set up bc for E
        implicit none
        integer ng,n,nP,nN,ig,tz,c
        real*8, allocatable :: xvec(:)
        ! real*8, allocatable :: u(:,:)
        real*8, allocatable :: En(:)
        real*8 t,tmp1,tmp2,dx,dt,xp,yp
        real*8 etmp1,etmp2
        real*8 ptmp1,ptmp2
        real*8 dtmp1,dtmp2,epsilon
        real*8 et,ett,ptt
        common tz
        real*8, allocatable :: u(:,:)
        real*8, allocatable :: un(:,:)

        ! speed of light
        c = 1.d0

        allocate(En(-ng:n+ng))

        En = un(-ng:n+ng,1)
        if (tz.eq.0) then ! initial boundary value problem
            ! left boundary: incomingWave
            call sechWave(etmp1,et,ett,ptmp1,pt,ptt,dtmp1,xvec(0),t,epsilon)
            En(0) = etmp1

            ! right boundary
            call sechWave(etmp2,et,ett,ptmp2,pt,ptt,dtmp2,xvec(n),t,epsilon)
            En(n) = etmp2

            ! extrapolation
            ! third order [3 -3 1]
            do ig=1,ng
                En(-ig) = 3.*En(-ig+1)-3.*En(-ig+2)+1.*En(-ig+3)
                En(n+ig) = 3.*En(n+ig-1)-3.*En(n+ig-2)+1.*En(n+ig-3)
            enddo
            ! fourth order accurate extrap coef = [4 -6 4 -1]
            ! do ig=1,ng
            !     En(-ig) = 4.*En(-ig+1)-6.*En(-ig+2)+4.*En(-ig+3)-1.*En(-ig+4)
            !     En(n+ig) = 4.*En(n+ig-1)-6.*En(n+ig-2)+4.*En(n+ig-3)-1.*En(n+ig-4)
            ! enddo

        elseif (tz.eq.1) then! twilight zone
            call getExactSolution(tmp1,xvec(0),t)
            call getExactSolution(tmp2,xvec(n),t)
            En(0) = tmp1
            En(n) = tmp2

            ! cheat ghost points, need to be replaced by extrapolation
            ! do ig=1,ng
            !     call getExactSolution(tmp1,xvec(-ig),t)
            !     En(-ig) = tmp1
            !     call getExactSolution(tmp2,xvec(n+ig),t)
            !     En(n+ig) = tmp2
            ! enddo
            
            ! extrapolation
            ! third order [3 -3 1]
            do ig=1,ng
                En(-ig) = 3.*En(-ig+1)-3.*En(-ig+2)+1.*En(-ig+3)
                En(n+ig) = 3.*En(n+ig-1)-3.*En(n+ig-2)+1.*En(n+ig-3)
            enddo
            ! fourth order accurate extrap coef = [4 -6 4 -1]
            ! do ig=1,ng
            !     En(-ig) = 4.*En(-ig+1)-6.*En(-ig+2)+4.*En(-ig+3)-1.*En(-ig+4)
            !     En(n+ig) = 4.*En(n+ig-1)-6.*En(n+ig-2)+4.*En(n+ig-3)-1.*En(n+ig-4)
            ! enddo
        endif

        un(-ng:n+ng,1) = En

        deallocate(En)
        
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

    subroutine getExactSolution(u,x,t) ! exact degree 2 polynomials
        implicit none
        real*8 u,x,t
        real*8 a0, a1, a2, c0, c1, c2

        a0 = 1.d0
        a1 = 1.d0
        a2 = 1.d0
        
        c0 = 1.d0
        c1 = 1.d0
        c2 = 1.d0
        

        u = (a2*x**2+a1*x+a0)*(c2*t**2+c1*t+c0)

    end subroutine getExactSolution

    subroutine timeStep1D(pnec,peptc,ng,n,alphaP,dx,dt,xvec,t,b0,b1,c,um,u,un,iOrder)
        implicit none
        real*8 dx,dt,t,c,Psum,alphaP,x,tmp1,tmp2
        integer iOrder,ng,n,i,j,k,tz
        real*8, allocatable :: xvec(:)
        ! real*8, allocatable :: um(:,:),u(:,:)!,un(:,:)
        real*8, allocatable :: Em(:),E(:),En(:)
        real*8, allocatable :: Pm(:),P(:),Pn(:)
        real*8, allocatable :: Dm(:),D(:),Dn(:)
        real*8, allocatable :: fe(:)
        real*8, allocatable :: fp(:)
        real*8, allocatable :: fn(:),fnt(:)
        real*8 b0,b1,beta,pnec,peptc,Nt,Ntt
        ! real*8 un(-1:n+1,3)
        real*8 um(-1:n+1,3)!,u(-1:n+1,3)
        real*8, allocatable :: u(:,:)
        real*8, allocatable :: un(:,:)

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
        allocate(fp(-ng:n+ng))
        allocate(fn(-ng:n+ng))
        allocate(fnt(-ng:n+ng))

        ! forcing functions
        if (tz==1) then
            call getForcing(b0,b1,n,ng,nP,nN,fe,fp,fn,fnt,xvec,t-dt,pnec,peptc,c,alphaP)
        else ! not twilight zone
            do i=0,n
                fe(i) = 0.
                fp(i) = 0.
                fn(i) = 0.
                fnt(i) = 0.
            enddo   
        endif

        Em = um(-ng:n+ng,1)
        E  = u(-ng:n+ng,1)
        Pm = um(-ng:n+ng,2)
        P  = u(-ng:n+ng,2)
        Dm = um(-ng:n+ng,3)
        D  = u(-ng:n+ng,3)

        do i=0,n ! loop over space
            x = xvec(i)
            ! second order update
            ! P
            Pn(i) = 2.*P(i)-Pm(i)+0.5*dt*b1*Pm(i)-dt**2*b0*P(i)+dt**2*fp(i)
            Pn(i) = Pn(i)+dt**2*pnec*D(i)*E(i)
            Pn(i) = beta*Pn(i) ! adjust
            Psum = Pn(i)-2.*P(i)+Pm(i)
            un(i,2) = Pn(i)
        
            ! E
            En(i) = 2.*E(i)-Em(i) + dt**2*c**2*(E(i+1)-2.*E(i)+E(i-1))/dx**2-alphaP*Psum+dt**2*fe(i)
            un(i,1) = En(i)

            ! Nt
            Nt = peptc*E(i)*(Pn(i)-Pm(i))/(2.*dt) + fn(i)

            ! Ntt
            Ntt =  peptc*(En(i)-Em(i))/(2.*dt)*(Pn(i)-Pm(i))/(2.*dt) &
                 + peptc*E(i)*(Pn(i)-2.*P(i)+Pm(i))/(dt**2) &
                 + fnt(i)

            ! N
            Dn(i) = D(i)+dt*Nt+dt**2/2.*Ntt
            un(i,3) = Dn(i)

        enddo

        ! fourth order update
        ! if (iOrder.eq.4) then

        ! endif

        deallocate(Em,E,En,Pm,P,Pn,Dm,D,Dn,fe,fp,fn,fnt)

    end subroutine timeStep1D

    subroutine getForcing(b0,b1,n,ng,nP,nN,fe,fp,fn,fnt,xvec,t,pnec,peptc,c,alphaP)
        implicit none

        integer nP,nN,ng,n,j,k,i
        real*8, allocatable :: xvec(:)
        real*8, allocatable :: fe(:)
        real*8, allocatable :: fp(:)
        real*8, allocatable :: fn(:),fnt(:)
        real*8 t,x,u,ux,uxx,ut,utt,c,alphaP
        real*8 a0, a1, a2, c0, c1, c2
        real*8 b0,b1,beta,pnec,peptc

        a0 = 1.d0
        a1 = 1.d0
        a2 = 1.d0

        c0 = 1.d0
        c1 = 1.d0
        c2 = 1.d0

        do i=0,n
            x   = xvec(i)
            u   = (a2*x**2+a1*x+a0)*(c2*t**2+c1*t+c0)
            ux  = (2.*a2*x+a1)*(c2*t**2+c1*t+c0)
            uxx = (2.*(a2))*(c2*t**2+c1*t+c0)
            ut  = (a2*x**2+a1*x+a0)*(2.*c2*t+c1)
            utt = (a2*x**2+a1*x+a0)*2.*c2

            ! fe
            fe(i) = utt-c**2*uxx+alphaP*utt

            ! fp
            fp(i) = utt+b1*ut+b0*u-pnec*u*u

            ! fn
            fn(i) = ut-peptc*u*ut
            fnt(i) = utt-peptc*ut*ut-peptc*u*utt
        enddo

    end subroutine getForcing

end subroutine mbe1d













