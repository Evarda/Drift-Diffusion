module globalConstants

    ! General/Given Constants
    real (kind=16):: e=1.60217662d-19       ! Charge of Electron                    [C]
    real (kind=16):: kb=8.6173303d-5        ! Boltzmann Constant                    [eV/K]
    real (kind=16):: eps=1.05d-12           ! Permitivity of Silicon Diode          [F/cm]

    real (kind=16) :: T=300                 ! Temperature                           [K]
    real (kind=16) :: ni=1.5d10             ! Concentration of Intrinsic Carriers   [cm^-3]
    real (kind=16) :: Na=1.0d16             ! Concentration of Acceptors            [cm^-3]
    real (kind=16) :: Nd=1.0d17             ! Concentration of Donors               [cm^-3]
    real (kind=16) :: Eg=1.12               ! Gap nmax                              [eV]
    real (kind=16) :: Vs                    ! Applied Voltage                       [eV]
    real (kind=16) :: Vt                    ! Thermal Voltage                       [eV]

    real (kind=16) :: kbT

    ! Calculated Values
    real (kind=16) :: Vbi                   ! Built in Voltage                      [V]
    real (kind=16) :: xn0                   ! Distance from 0 to n Equilibrium      [cm]
    real (kind=16) :: xp0                   ! Distance from 0 to p Equilibrium      [cm]
    real (kind=16) :: w0                    ! Distance from n and p Equilibriums    [cm]

    ! Boundary Conditions for Potential
    real (kind=16) :: psip                  ! psi at end of p side                  []
    real (kind=16) :: psin                  ! psi at end of n side                  []

    ! Debye Lengths (for Uniform Mesh)
    real (kind=16) :: Ldi                   ! Debye Length for Intrinsic Carriers   [cm]
    real (kind=16) :: Lda                   ! Debye Length for Acceptors            [cm]
    real (kind=16) :: Ldd                   ! Debye Length for Donors               [cm]

    ! Uniform Mesh
    real :: meshScale = 0.05                ! Scale Factor for Mesh Size            [cm]
    real (kind=16) :: dx                    ! Mesh Spacing                          [m]
    real (kind=16) :: dx2                   ! Mesh Spacing Squared                  [m^2]
    real (kind=16) :: xmax                  ! Mesh Size, Device Length              [m]
    integer :: nmax                         ! Max Number of Points in Array         [None]

    real :: deviceLength = 10               ! Scale Factor for Device Length

    ! Boundary Conditions for Concentration Values
    real (kind=16) :: nn
    real (kind=16) :: np
    real (kind=16) :: pn
    real (kind=16) :: pp

    ! Convergence Check Values
    real (kind=16) :: delta
    real (kind=16) :: crit_conv
    logical :: flag_conv

    ! SRH Generation-Recombination Parameters
    ! From Atlas User's Manual, Nov 1998, Page 3-61
    real :: TAUN0=1d-7          ! Electron Lifetime                     [s]
    real :: TAUP0=1d-7          ! Hole     Lifetime                     [s]
    real :: NSRHN=5d16
    real :: NSRHP=5d16

    ! Carrier Lifetimes
    real (kind=16) :: TAUN1
    real (kind=16) :: TAUN2
    real (kind=16) :: TAUP1
    real (kind=16) :: TAUP2

    ! Variables as a Function of Position
    real (kind=16), dimension(:), allocatable :: V          ! Voltage at position in junction       [V]

    ! LU Decomposition
    ! From Numerical Solution of Poisson's Equation - nanoHUB.org
    ! Define allocatable arrays of future size nmax
    real (kind=16), dimension(:), allocatable :: a          ! Array for a diagonal values in A
    real (kind=16), dimension(:), allocatable :: b          ! Array for b diagonal values in A
    real (kind=16), dimension(:), allocatable :: c          ! Array for c diagonal values in A - same in U

    real (kind=16), dimension(:), allocatable :: alpha      ! Array for main  diagonal in U
    real (kind=16), dimension(:), allocatable :: beta       ! Array for lower diagonal in L
    real (kind=16), dimension(:), allocatable :: f          ! Array for values of Forcing Function
    real (kind=16), dimension(:), allocatable :: g          ! Array for values of Decomposition Function (Equivalent to y in slides)
    real (kind=16), dimension(:), allocatable :: x          ! Array for values of Position                                                      [cm]

    real (kind=16), dimension(:), allocatable :: psi        ! Array for values of Electrostatic Potential
    real (kind=16), dimension(:), allocatable :: n          ! Array for values of Concentration of Electrons/ni     n=exp( psi)
    real (kind=16), dimension(:), allocatable :: p          ! Array for values of Concentration of Holes/ni         p=exp(-psi)
    real (kind=16), dimension(:), allocatable :: oldpsi     ! Array for values of old Electrostatic Potential

    real (kind=16), dimension(:), allocatable :: voltage    ! Array for values of Applied Voltage
    real (kind=16), dimension(:), allocatable :: Rsrh       ! Shockley-Read-Hall Recombination/Generation Rate

    real (kind=16), dimension(:), allocatable :: dop        ! Array for values of Doping

    ! Arrays of Calculated Output Values
    real (kind=16), dimension(:), allocatable :: elfield    ! Array for values of Electric Field
    real (kind=16), dimension(:), allocatable :: density    ! Array for values of Charge Density
    real (kind=16), dimension(:), allocatable :: condband   ! Array for values of Electric Field
    real (kind=16), dimension(:), allocatable :: current    ! Array for values of Current

    real (kind=16), dimension(:,:), allocatable :: current_n   ! Array for values of Currentn
    real (kind=16), dimension(:,:), allocatable :: current_p   ! Array for values of Currentp

    real (kind=16), dimension(:), allocatable :: TAUN       ! Array for values of Electron Lifetime
    real (kind=16), dimension(:), allocatable :: TAUP       ! Array for values of Hole Lifetime

    ! Do Loops
    real (kind=16) :: y                                     ! Intermediate Values for Psi Calculation
    real (kind=16) :: y2
    real (kind=16) :: z
    integer :: i                                            ! Step Counter for Do Loops
    integer :: j                                            ! Step Counter for Current
    integer :: k                                            ! Step Counter for Voltage
    integer :: k_iter                                       ! Step Counter for iterations for Poisson Convergence (Equilibrium)
    integer :: k_iternoneq                                  ! Step Counter for iterations for Poisson Convergence (Non-Equilibrium)
    integer :: iterno                                       ! Number of iterations
    integer :: iternoneq                                       ! Number of iterations

    ! Normalization Constants
    real (kind=16) :: dif0                                  ! Normalization fo Diffusion Coefficient [cm^2]
    real (kind=16) :: J0                                    ! Current Normaralization [A/cm^2]
    real (kind=16) :: tau_0                                 ! Lifetime Normalization [s]

    ! Arora Model for Low Field Mobility
    ! Atlas User's Manual, Nov 1998, Page 3-35, 3-36
    ! Table 3-19
    real (kind=16) :: MU1N=88.0       ! [cm^2/(V*s)]
    real (kind=16) :: MU1P=54.3       ! [cm^2/(V*s)]
    real (kind=16) :: MU2N=1252.0     ! [cm^2/(V*s)]
    real (kind=16) :: MU2P=407.0      ! [cm^2/(V*s)]
    real (kind=16) :: ALPHAN=-0.57    ! [None]
    real (kind=16) :: ALPHAP=-0.57    ! [None]
    real (kind=16) :: BETAN=-2.33     ! [None]
    real (kind=16) :: BETAP=-2.23     ! [None]
    real (kind=16) :: GAMMAN=2.546    ! [None]
    real (kind=16) :: GAMMAP=2.546    ! [None]
    real (kind=16) :: NCRITN=1.432d17 ! [cm^-3]
    real (kind=16) :: NCRITP=2.67d17  ! [cm^-3]

    ! Values for Mobility Constants
    real (kind=16) :: MUN1
    real (kind=16) :: MUN2
    real (kind=16) :: MUP1
    real (kind=16) :: MUP2

    !Non-Equilibrium
    ! Temperature Dependence of the Saturation Electron Drift Velocity
    real (kind=16) :: vsatn                     ! Saturation Velocity on n side
    real (kind=16) :: vsatp                     ! Saturation Velocity on p side
    real :: vso
    real :: Ccon
    real :: Icon

    ! Voltage Parameters
    real :: Vmax = 1d0                          ! Voltage Max
    integer :: nvstep = 60                      ! Non-Eq Voltage Step
    real :: delta_V = 1d-2

    ! Field-Dependent Mobility (Dragica Slides are wrong on nanohub)
    real :: fdmbetan = 1
    real :: fdmbetap = 2

    ! Diffusion Constants
    real (kind=16), dimension(:), allocatable :: d_n
    real (kind=16), dimension(:), allocatable :: d_p

    ! Field Dependent Mobilities
    real (kind=16), dimension(:), allocatable :: bb_n
    real (kind=16), dimension(:), allocatable :: bb_p

    ! Intermediate Store Values of Bernoulli Functions and in
    real (kind=16) :: AA
    real (kind=16) :: CC
    real (kind=16) :: ZZ
    real (kind=16) :: YY
    real (kind=16) :: baux
    real (kind=16) :: caux

    ! Output File Names
    character(100) :: diffnname
    character(100) :: diffpname
    character(100) :: psiname
    character(100) :: nname
    character(100) :: pname
    character(100) :: Rsrhname
    character(100) :: bbnname
    character(100) :: bbpname

    character(100) :: nonequilibriumcondband
    character(100) :: nonequilibriumelectricfield
    character(100) :: nonequilibriumconcp
    character(100) :: nonequilibriumconcn
    character(100) :: nonequilibriumconcden
    character(100) :: nonequilibriumquasifermin
    character(100) :: nonequilibriumquasifermip

    character(100) :: namecurrent
    character(100) :: namecurrentn
    character(100) :: namecurrentp

  end module