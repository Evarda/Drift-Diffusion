subroutine generateMesh
    use globalConstants
    implicit none

    ! Uniform Mesh
        ! [1] Mesh Size Slide
        ! [0] Unit Analysis for Equation Unit Analysis
    
        ! QUESTION: Where does this come from? I found this equation for Acceptors/Donors, but not for intrinsic?
        ! Intrinsic Debye Length (Used in Normalization?)
    Ldi=sqrt((eps*kb*T)/(ni*e))
    
    !! Acceptor Extrinsic Debye Length
    Lda=sqrt((eps*kb*T)/(Na*e))
    
    !! Donor Extrinsic Debye Length
    Ldd=sqrt((eps*kb*T)/(Nd*e))
    
    !! dx as Minimum Extrinsic Debye Length
    ! meshScale < 1 acts as scale factor to refine mesh
    dx=meshScale*min(Lda, Ldd)

    !! [2] Page 222 of gives code in Computations Section
    
    !! Built in Voltage                     [2: Eq 5.10]
    ! Vbi = ((kb*T)/q) * log((Na*Nd)/(ni**2))
    ! Know: 1 eV = e J and J/C = V
    ! kT/q = (kb*T*[eV/K]*[K])/(e*[C]) = (kb*T*[eV])/(e*[C]) * (e [J])/[eV] = (kb*T*[J])/[C] = kb*T [V]
    Vbi=(kb*T)*log((Na*Nd)/(ni**2))
    
    !! Distance from 0 to n Equilibrium     [2: Eq 5.30a]
    ! ks*ep0=eps
    xn0=sqrt(((2*eps)/e)*(Na/(Nd*(Na+Nd)))*Vbi)
    
    !! Distance from 0 to p Equilibrium     [2: Eq 5.22]
    xp0=xn0*(Nd/Na)
    
    !! Distance from n and p Equilibrium    [2: Eq 5.31]
    w0=xn0+xp0
    
    !! Device Length
    ! Set xmax as 10* size of 0 to equilibrium distance
    xmax=deviceLength*max(xn0,xp0)
    
    ! Set nmax Points
    nmax=int(xmax/dx)
    
    ! Re-Normalize dx
    ! QUESTION: What purpose does this serve? Why do we need to renormalize?
    dx=dx/Ldi
end subroutine