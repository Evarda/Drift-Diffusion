! Drift-Diffusion Code
! Michelle King
! Written 2017, Reformat 2019
    
program driftDiffusion
    use globalConstants
    implicit none
    
    ! Bernoulli Function
    real (kind=16) :: BB

    ! Declare Values for Subroutines
        ! Declare for computeSurfaceDensities subroutine
        real :: ns
        real :: ps
        real :: psis
        
    ! Declare Fundamental Constants
        kbT=kb*T
        Vt = kbT / e   ! Thermal Voltage
        Vs=0           ! Applied Voltage
    
    ! Compute Surface Densities ns, ps and Calulate Surface Potential psis
        call computeSurfaceDensities(ns, ps, psis) 
    
    ! Generate Mesh
        call generateMesh()
        
        
        ! Normalization Constants
        ! Normalization for Diffusion Coefficient [cm^2]
        dif0=1
        ! Current Normalization [A/cm^2]
        J0=e*ni*dif0/Ldi
        ! Lifetime Normalization [s]
        tau_0=Ldi**2/dif0
    
    !------------------------------------------Gumel's Decoupled Relaxation Scheme/Method--------------------------------------------!
    
    ! Solve the poisson equation with a guess for the quasi-fermi levels (use applied voltage as initial guess
    ! The potential is used to update the bernouli functions
    ! The equations (notes) are solved to provide an update for the quasi-fermi levels, that enter into the poisson equation
    
    !! Solve the Poisson Equation with a Guess for the Quasi-Fermi Levels (Use Applied Voltage as Initial Guess)
    
    !! Make Doping Array
    
    ! Allocate memory size nmax for array
    allocate (dop(nmax))
    
    ! Normalize Doping Concentrations with ni
    ! Set first  half of dop(i) with doping levels Na/ni
    ! Set second half of dop(i) with doping levels Nd/ni
    do i = 1, nmax
        if(i.le.(nmax/2)) then
            dop(i)=-Na/ni
        else
            dop(i)=Nd/ni
        endif
    enddo
    
    ! Guess for Potential from Charge Neutrality
    ! [2: Eq 2.29 (Factored)]
    ! Put into Quasi-Fermi Level Relationships [2: Eq 3.72]
    ! QUESTION: How are minus signs found? - Try working it out with the appropriate +/- Na or Nd
    
    allocate(psi(nmax))
    
    do i=1, nmax
        y=.5*dop(i)
        if (y.gt.0.) then
            y2=(y**2)
            z= 1/y2
            psi(i)=log(y*(1+sqrt(1+z)))
        else
            y2=(y**2)
            z= 1/y2
            psi(i)=log(y*(1-sqrt(1+z)))
        endif
        !if (abs(psi(i)).gt.9999999999) print *, "psi(", i, ") is very large"
        !print *, "psi(", i, ")", psi(i)
    enddo
    
    ! Potential Boundary Conditions
    psip=psi(1)
    psin=psi(nmax)
    
    ! Concentrations of Particles Boundary Conditions
    y=Nd/(2*ni)
    nn=y+sqrt(y**2+1)
    y=-Na/(2*ni)
    np=y+sqrt(y**2+1)
    pn=1/nn
    pp=1/np

    
    !---------------------------------------------------------Iteration Loop---------------------------------------------------------!
    
    ! Allocate Memory for Arrays
    allocate (a(nmax), b(nmax), c(nmax), alpha(nmax), beta(nmax), f(nmax), g(nmax), x(nmax), oldpsi(nmax))
    
    dx2=dx**2               ! Store squared value for future calculations
    k_iter=0                ! Clear iteration counter
    flag_conv = .false.     ! Clear convergence flag
    crit_conv = 1d-6        ! Set calue of convergence
    
    do while(.not.flag_conv)

        k_iter=k_iter+1                         ! Increase iteration counter
        print *, "Loop Check: iter=", k_iter    ! Print iteration counter
        
        do i=1,nmax
            oldpsi(i)=psi(i)                    ! Store previous iteration's potential
        enddo
        
        ! Populate A matrix containing arrays a, b, c
        do i=2, nmax-1
            a(i)=(-2/dx2)-exp(psi(i))-exp(-psi(i))
            b(i)=1/(dx2)
            c(i)=1/(dx2)
            f(i)=exp(psi(i))-exp(-psi(i))-dop(i)-(psi(i)*(exp(psi(i))+exp(-psi(i))))
        enddo
        
        ! Ohmic Contacts -> Dirichlet Boundary Conditions
        a(1)=1
        b(1)=0
        c(1)=0
        f(1)=psip
    
        a(nmax)=1
        b(nmax)=0
        c(nmax)=0
        f(nmax)=psin ! a(n)=1 forces f(n) to be psi(n) at point n

        call LUDecomp()
        
        ! Save calculated potential
        do i=1,nmax
            psi(i)=x(i)
        enddo
    
        ! Clear Delta (Convergence Check Variable)
        delta=0
    
        ! Set Delta as Maximum Change in New and Old Arrays
        do i=1,nmax
            if (abs(oldpsi(i)-psi(i)).gt.delta) delta=abs(oldpsi(i)-psi(i))
        enddo
        
        print *, "Loop Check: delta=", delta
        
        ! Check if Delta has reached Critical Convergence, change flag to true if convergence is achieved
        if (delta.lt.crit_conv) then
            flag_conv = .true.
            iterno = k_iter
            goto 888
        endif
    enddo 
    
888 print *, "Convergence Acheived, Iteration Number=", iterno
    print *, ""
    
    !------------------------------------------------Output of Equilibrium Solution--------------------------------------------------!   
    
    ! Allocate arrays for electron concentration and hole concentration, electric field, density, and conduction band
    allocate (n(nmax), p(nmax), elfield(nmax), density(nmax), condband(nmax))
    
    ! Calculate electron and hole concentration
    do i=1, nmax
        n(i)=exp(psi(i))
        p(i)=exp(-psi(i))
    enddo
    
    ! Write Mesh Values
    open(unit=33, file='Data/mesh.txt', status='unknown')
    do i=1,nmax
        write(33, *) float(i-1)*dx*Ldi
    enddo
    close(33)
    
    ! Calculate electric field, deinsity, and conduction band
    do i=2, nmax-1
        elfield(i)=(psi(i)-psi(i+1))/dx
        density(i)=p(i)-n(i)-dop(i)
        condband(i)=0.5*Eg-kbT*psi(i)
    enddo
    
    ! Define Boundaries
    elfield(nmax)=elfield(nmax-1)
    elfield(1)=elfield(2)
    
    p(nmax)=p(nmax-1)
    p(1)=p(2)
    
    n(nmax)=n(nmax-1)
    n(1)=n(2)
    
    condband(nmax)=condband(nmax-1)
    condband(1)=condband(2)
    
    density(nmax)=density(nmax-1)
    density(1)=density(2)
    
    ! Write values for potiential, conduction band, electric field, electron and hold concentration, and density
    open(unit=33, file="Data/Equilibrium/potential_equil.txt",         status="unknown")
    open(unit=44, file="Data/Equilibrium/cond_band_equil.txt",         status="unknown")
    open(unit=55, file="Data/Equilibrium/electric_field_equil.txt",    status="unknown")
    open(unit=66, file="Data/Equilibrium/conc_p_equil.txt",            status="unknown")
    open(unit=77, file="Data/Equilibrium/conc_n_equil.txt",            status="unknown")
    open(unit=88, file="Data/Equilibrium/conc_den_equil.txt",          status="unknown")
    
    do i=1, nmax
        write(33,*) psi(i)
        write(44,*) condband(i)
        write(55,*) elfield(i)
        write(66,*) ni*p(i)
        write(77,*) ni*n(i)
        write(88,*) ni*density(i)
    enddo

    close(88)
    close(77)
    close(66)
    close(55)
    close(44)
    close(33)
    
    !---------------------------------------------------Non-Equilibrium Solution-----------------------------------------------------!
    
    ! Calculate Relaxation Times for SRH Concentration Dependent Lifetime Model
    TAUN1=TAUN0/(1+(Na/NSRHN))
    TAUN2=TAUN0/(1+(Nd/NSRHN))
    TAUP1=TAUP0/(1+(Na/NSRHP))
    TAUP2=TAUP0/(1+(Nd/NSRHP))
    
    ! Allocate memory for Electron and Hole Lifetimes
    allocate (TAUN(nmax), TAUP(nmax))
    
    ! Calculate electron and hole lifetimes
    do i=1, nmax
        if(i.le.(nmax/2)) then
            TAUN(i)=TAUN1/TAU_0
            TAUP(i)=TAUP1/TAU_0
        else
            TAUN(i)=TAUN2/TAU_0
            TAUP(i)=TAUP2/TAU_0
        endif
    enddo
    
    ! Calculate Mobilities with Arora Model for Low Field Mobility
    MUN1 = (MU1N*((T/300)**ALPHAN)+(MU2N*((T/300)**BETAN)/(1+(Na/(NCRITN*((T/300)**GAMMAN))))))
    MUN2 = (MU1N*((T/300)**ALPHAN)+(MU2N*((T/300)**BETAN)/(1+(Nd/(NCRITN*((T/300)**GAMMAN))))))
    MUP1 = (MU1P*((T/300)**ALPHAP)+(MU2P*((T/300)**BETAP)/(1+(Na/(NCRITP*((T/300)**GAMMAP))))))
    MUP2 = (MU1P*((T/300)**ALPHAP)+(MU2P*((T/300)**BETAP)/(1+(Nd/(NCRITP*((T/300)**GAMMAP))))))
    
    ! Calculate Saturation Velocity for Field Dependent Mobility [5]
    ! Also given in "Semi-Classical Transport Theory" Dragica Vasileska
    vso=2.4d7                               ! [cm/s]
    Ccon=0.8                                ! [None]
    Icon=600                                ! [K]
    vsatn = vso/(1+(Ccon*exp(T/Icon)))
    vsatp = vsatn                           ! Same as vsatn since vsat is only T dependent
    
    
    ! Write constant values to a file
    open(unit=33, file=trim("Data/constants.txt"), status="unknown")
    write(33, *) "MUN1 = ", MUN1
    write(33, *) "MUN2 = ", MUN2
    write(33, *) "MUP1 = ", MUP1
    write(33, *) "MUP2 = ", MUP2
    write(33, *) "vsatn = ", vsatn
    write(33, *) "kbT = ", kbT
    write(33, *) "k = ", kb
    write(33, *) "T = ", T
    write(33, *) "psin = ", psin
    write(33, *) "psip = ", psip
    write(33, *) "nn = ", nn
    write(33, *) "np = ", np
    write(33, *) "pn = ", pn
    write(33, *) "pp = ", pp
    write(33, *) "ni*nn = ", ni*nn
    write(33, *) "ni*np = ", ni*np
    write(33, *) "ni*pn = ", ni*pn
    write(33, *) "ni*pp = ", ni*pp
    close(33)
    
    !---------------------------------------------------------Apply Voltage----------------------------------------------------------!
    
    ! Allocate Memory for Diffusion constants, field dependent mobilities, voltages, and generation/recombination rate
    allocate (d_n(nmax), d_p(nmax), bb_n(nmax), bb_p(nmax), voltage(nvstep), &
    Rsrh(nmax), current(nmax), current_n(nvstep, nmax), current_p(nvstep, nmax))
    
    open(unit=11, file="Data/voltage.txt", status="unknown")
    
    ! Voltage Stepping Loop
    do k=1, nvstep
        k_iternoneq = 0                 ! Reset Iteration Number
        flag_conv = .false.             ! Reset Convergence Flag
        
        voltage(k)=delta_V*float(k)/kbT ! Calculate and set voltage for voltage(k) array
        print *, "Voltage Increase, voltage=", voltage(k)*kbT, "eV"
        
        write(11, *) voltage(k)*kbT
        
        Vs = voltage(k)                            ! Set Vs to currrent voltage array value (Used?)
        call computeSurfaceDensities(ns, ps, psis) ! Compute Surface Densities for the k iteration (Used?)

        
       
        
        ! Poisson Convergence for voltage(k)
        do while (.not.flag_conv)
            
            k_iternoneq=k_iternoneq+1   ! Increase iteration counter for nonequilibrium k
            print *, "Loop Check: iternoneq=", k_iternoneq
            
            ! Store old potential
            do i=1, nmax
                oldpsi(i)=psi(i)        ! Used in Calculation for Convergence
            enddo
            
            ! Poisson Solver for Potential
            do i=1, nmax-1
                a(i)=(-2/dx2)-n(i)-p(i)
                b(i)=1/(dx2)
                c(i)=1/(dx2)
                f(i)=n(i)-p(i)-dop(i)-psi(i)*(n(i)+p(i))
            enddo
            
            ! Boundary Conditions
            a(1)=1
            b(1)=0
            c(1)=0
            f(1)=psip+voltage(k)
    
            a(nmax)=1
            b(nmax)=0
            c(nmax)=0
            f(nmax)=psin
            
            call LUDecomp()
            
            ! Store calculated potential
            do i=1, nmax
                psi(i)=x(i)
            enddo
            
            ! Write Down Potential
            write (psiname, '("Data/Potential/psi",i4, "-",i4,".txt")') k, k_iternoneq
        
            open(unit=33, file=trim(psiname), status="unknown")

            do i=1, nmax
                write(33,*) psi(i)
            enddo

            close(33)
            
            ! Generation/Recombination Rate from SRH Concentration Dependent Lifetime Model
            do i=1, nmax
                Rsrh(i)=-(n(i)*p(i)-1)/(TAUP(i)*(n(i)+1)+TAUN(i)*(p(i)+1))
            enddo
            
            ! Calculate Field Dependent Mobility
            do i=2, nmax-1
                aa=(psi(i)-psi(i+1))/(dx*Ldi)
                
                if (i.lt.(nmax/2)) then
                    bb_n(i)=1+((MUN1*aa/vsatn)**fdmbetan)
                    bb_n(i)=MUN1/((bb_n(i))**(1/fdmbetan))
                    
                    bb_p(i)=1+((MUP1*aa/vsatp)**fdmbetap)
                    bb_p(i)=MUP1/((bb_p(i))**(1/fdmbetap))    

                else
                    bb_n(i)=1+((MUN2*aa/vsatn)**fdmbetan)
                    bb_n(i)=MUN2/((bb_n(i))**(1/fdmbetan))
                
                    bb_p(i)=1+((MUP2*aa/vsatp)**fdmbetap)
                    bb_p(i)=MUP2/((bb_p(i))**(1/fdmbetap))
                      
                endif
            enddo
        
            ! Calculate Diffusion Constants from Mobility         
            do i=2, nmax-1
                d_n(i)=(kbT)*bb_n(i)
                d_p(i)=(kbT)*bb_p(i)
            enddo
            
            ! Diffusion Constant Boundary Conditions
            d_n(1)=d_n(2)
            d_p(1)=d_p(2)
            
            d_n(nmax)=d_n(nmax-1)
            d_p(nmax)=d_p(nmax-1)

            ! Write Down Diffusion Constants for runs
            write (diffnname, '("Data/Diffusion/Diffn",i4, "-",i4,".txt")') k, k_iternoneq
            write (diffpname, '("Data/Diffusion/Diffp",i4, "-",i4,".txt")') k, k_iternoneq
            write (bbnname, '("Data/Mu/bbn",i4, "-",i4,".txt")') k, k_iternoneq
            write (bbpname, '("Data/Mu/bbp",i4, "-",i4,".txt")') k, k_iternoneq
            
            open(unit=33, file=trim(bbnname), status="unknown") 
            open(unit=44, file=trim(bbpname), status="unknown")
            open(unit=55, file=trim(diffnname), status="unknown")
            open(unit=66, file=trim(diffpname), status="unknown")

            do i=1, nmax
                write(33, *) bb_n(i)
                write(44, *) bb_p(i)
                write(55, *) d_n(i)
                write(66, *) d_p(i)
            enddo
        
            close(33) 
            close(44)
            close(55)
            close(66)
           
        ! Poisson Solver for Electron Concentration
        
            ! Define Arrays
            do i=2, nmax-1
                AA=BB((psi(i+1)-psi(i)))
                CC=BB((psi(i-1)-psi(i)))
                ZZ=BB((psi(i)-psi(i+1)))
                YY=BB((psi(i)-psi(i-1)))
                
                c(i)=0.5*(d_n(i+1)+d_n(i))*AA/dx2
                b(i)=0.5*(d_n(i)+d_n(i-1))*CC/dx2
                
                caux=0.5*(d_n(i+1)+d_n(i))*ZZ/dx2
                baux=0.5*(d_n(i)+d_n(i-1))*YY/dx2
                
                a(i)=-(caux+baux)
                
                f(i)=Rsrh(i)
            enddo
            
            ! Ohmic BCs
            c(1)=0
            c(nmax)=0
            
            b(1)=0
            b(nmax)=0
            
            a(1)=1
            a(nmax)=1
            
            f(1)=np
            f(nmax)=nn
            
            call LUDecomp()
            
            ! Store calculated values of electron concentration
            do i=1, nmax
                n(i)=x(i)
            enddo
            
            ! Write Down electron concentration
            write (nname, '("Data/n/n",i4, "-",i4,".txt")') k, k_iternoneq
        
            open(unit=33, file=trim(nname), status="unknown")

            do i=1, nmax
                write(33,*) n(i)
            enddo

            close(33)
            
        ! Hole Continuity
            ! Update Generation Rate
            do i=1, nmax
                Rsrh(i)=-(n(i)*p(i)-1)/(TAUP(i)*(n(i)+1)+TAUN(i)*(p(i)+1))
            enddo 
            
            ! Define Arrays
            do i=2, nmax-1            
                AA=BB((psi(i+1)-psi(i)))
                CC=BB((psi(i-1)-psi(i)))
                ZZ=BB((psi(i)-psi(i+1)))
                YY=BB((psi(i)-psi(i-1)))
                
                c(i)=0.5*(d_p(i+1)+d_p(i))*ZZ/dx2
                b(i)=0.5*(d_p(i)+d_p(i-1))*YY/dx2
                
                caux=0.5*(d_p(i+1)+d_p(i))*AA/dx2
                baux=0.5*(d_p(i)+d_p(i-1))*CC/dx2
                
                a(i)=-(caux+baux)
                
                f(i)=Rsrh(i)
            enddo
            
            ! Ohmic BCs
            c(1)=0
            c(nmax)=0
            
            b(1)=0
            b(nmax)=0
            
            a(1)=1
            a(nmax)=1
            
            f(1)=pp
            f(nmax)=pn
            
            call LUDecomp()
            
            ! Store calculated values of hole concentration
            do i=1, nmax
                p(i)=x(i)
            enddo
            
            ! Write Down hole concentration
            write (pname, '("Data/p/p",i4, "-",i4,".txt")') k, k_iternoneq
        
            open(unit=33, file=trim(pname), status="unknown")

            do i=1, nmax
                write(33,*) p(i)
            enddo

            close(33)
            
            ! Clear Delta (Convergence Check Variable)
            delta=0
    
            ! Set Delta as Maximum Change in New and Old Arrays
            do i=1, nmax
                if (abs(oldpsi(i)-psi(i)).gt.delta) then
                    delta=abs(oldpsi(i)-psi(i))
                endif
            enddo
            
            print *, "Loop Check: delta=", delta

            if (delta.lt.crit_conv) then
                flag_conv = .true.
                iterno = k_iternoneq
                goto 1888
            endif
            
        enddo ! Ends Poisson Convergence
        
1888    print *, "Convergence Acheived, Iteration Number=", iterno
        print *, "V=", voltage(k)*kbT, "eV"
        print *, ""
        
        do i=2, nmax-1
            elfield(i)=(psi(i)-psi(i+1))/dx
            density(i)=p(i)-n(i)-dop(i)
            condband(i)=0.5*Eg-kbT*psi(i)
        enddo
        
        !BCS for values calculated in main part of loop
        p(nmax)=p(nmax-1)
        p(1)=p(2)
        n(nmax)=n(nmax-1)
        n(1)=n(2)
        
        !BCS for values calculated after main part of loop
        elfield(nmax)=elfield(nmax-1)
        elfield(1)=elfield(2)
        density(nmax)=density(nmax-1)
        density(1)=density(2)
        condband(nmax)=condband(nmax-1)
        condband(1)=condband(2)
        
        ! Write Outpute for Non-equilibrium Solution
        
            write (nonequilibriumcondband,      '("Data/Non-Equilibrium/condband",i4, ".txt")') k
            write (nonequilibriumelectricfield, '("Data/Non-Equilibrium/electricfield",i4, ".txt")') k
            write (nonequilibriumconcp,         '("Data/Non-Equilibrium/concp",i4, ".txt")') k
            write (nonequilibriumconcn,         '("Data/Non-Equilibrium/concn",i4, ".txt")') k
            write (nonequilibriumconcden,       '("Data/Non-Equilibrium/concden",i4, ".txt")') k
            write (nonequilibriumquasifermin,    '("Data/Non-Equilibrium/quasifermin",i4, ".txt")') k
            write (nonequilibriumquasifermip,    '("Data/Non-Equilibrium/quasifermip",i4, ".txt")') k
        
            open(unit=44, file=trim(nonequilibriumcondband),         status="unknown")
            open(unit=55, file=trim(nonequilibriumelectricfield),    status="unknown")
            open(unit=66, file=trim(nonequilibriumconcp),            status="unknown")
            open(unit=77, file=trim(nonequilibriumconcn),            status="unknown")
            open(unit=88, file=trim(nonequilibriumconcden),          status="unknown")
            open(unit=99, file=trim(nonequilibriumquasifermin),      status="unknown")
            open(unit=101,file=trim(nonequilibriumquasifermip),      status="unknown")
        
        do i=1, nmax

            write(44, *) condband(i)
            write(55, *) elfield(i)
            write(66, *) ni*p(i)
            write(77, *) ni*n(i)
            write(88, *) ni*density(i)
            write(99, *)  kbT*(-psi(i)+log(n(i)))
            write(101, *) kbT*(-psi(i)-log(p(i)))
     
        enddo
        
            close(44)
            close(55)
            close(66)
            close(77)
            close(88)
            close(99)
        
            
        ! Current Reset
        current(k)=0
    
        ! Current by average current
        do j=2, nmax-1
            aa=n(j  )*BB(psi(j  )-psi(j-1))-n(j-1)*BB(psi(j-1)-psi(j  ))  ! Equiv to Jn(i-1/2)*dx
            cc=n(j+1)*BB(psi(j+1)-psi(j  ))-n(j  )*BB(psi(j  )-psi(j+1)) ! Equiv to Jn(i+1/2)*dx
            current_n(k, j)=(0.5*(d_n(j)+d_n(j-1))*aa+0.5*(d_n(j)+d_n(j+1))*cc)/(2*dx) ! Jn(i) is average of Jn(i-1/2) and Jn(i+1/2)
            
            aa=p(j  )*BB(psi(j-1)-psi(j  ))-p(j-1)*BB(psi(j  )-psi(j-1))  ! Equiv to Jp(i-1/2)*dx
            cc=p(j+1)*BB(psi(j  )-psi(j+1))-p(j  )*BB(psi(j+1)-psi(j  )) ! Equiv to Jp(i+1/2)*dx
            current_p(k, j)=(0.5*(d_p(j)+d_p(j-1))*aa+0.5*(d_p(j)+d_p(j+1))*cc)/(2*dx) ! Jp(i) is average of Jp(i-1/2) and Jp(i+1/2)
            
            current(k)=current(k)+current_n(k, j)+current_p(k, j)
        enddo
        ! Average Current
        current(k)=current(k)/float(nmax-2)       
        ! Current n and p BCs
        current_n(k,1)=current_n(k,2)
        current_n(k,nmax)=current_n(k,nmax-1)
        
        current_p(k,1)=current_p(k,2)
        current_p(k,nmax)=current_p(k,nmax-1)
        
     
        
    enddo ! Ends Voltage Stepping
    
    close(11)
    
    
        
    open(unit=44, file="Data/current.txt",  status="unknown")
        
    print *, "Enter Current Write Loop"
    do k=1, nvstep
    write(44, *) current(k)*J0
    print *, "current(k) written"
    
    write (namecurrentn, '("Data/Non-Equilibrium/currentn",i4, ".txt")') k
    write (namecurrentp, '("Data/Non-Equilibrium/currentp",i4, ".txt")') k
    
    open(unit=55, file=trim(namecurrentn), status="unknown")
    open(unit=66, file=trim(namecurrentp), status="unknown")
        print *, "Enter inner current loop"
        do j=1, nmax
            write(55, *) current_n(k,j)*J0
            write(66, *) current_p(k,j)*J0
        enddo
        print *, "Exit inner current loop"
    close(55)
    close(66)
    
    
    enddo
    
    print *, "Exit Current Write Loop"
    
    close(44)
    
    

    
    end program driftDiffusion

    
!************************************************************************************************************************************!
    
    subroutine computeSurfaceDensities(ns, ps, psis)
    use globalConstants
    
    real, intent(out) :: ns
    real, intent(out) :: ps
    real, intent(out) :: psis
    
    ! Boundary Physics - Ohmic Contacts
    ! Atlas User's Manual, Nov 1998, Page 3-23
    
    ! Equation 3-104
    ns=.5*((Nd-Na)+sqrt(((Nd-Na)**2)+(4*(ni**2))))
    
    ! Equation 3-105
    ps=(ni**2)/ns

    ! Equation 3-106ni^2/
    psis=Vs+(((k*T)/e)*log(ns/ni))

    end subroutine computeSurfaceDensities

!************************************************************************************************************************************!
    real (kind=16) function BB(x)
    implicit real*16 (a-h, o-z)
    x1_b=-79.0188
    x2_b=-9d-6
    x3_b=1d-6
    x4_b=79.0188
    x5_b=11433.4628
    !x5_b=9.487

    if(x.lt.x1_b) then
        BB=-x
    else if ((x.ge.x1_b).and.(x.lt.x2_b)) then
        BB=x/(exp(x)-1)
    else if ((x.ge.x2_b).and.(x.lt.x3_b)) then
        BB=1-x/2
    else if ((x.ge.x3_b).and.(x.lt.x4_b)) then
        BB=x*exp(-x)/(1-exp(-x))
    else if ((x.ge.x4_b).and.(x.lt.x5_b)) then
        BB=x*exp(-x)
    else
        BB=0
        print *, "CAUTION: BERNOULLI FUNCTION = 0"
        print *, "x=", x
    endif
    return
    end function

!************************************************************************************************************************************!
    
!************************************************************************************************************************************!
