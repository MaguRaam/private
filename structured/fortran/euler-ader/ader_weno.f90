module ader_weno
    
    implicit none
    public

    ! ---------------------------------------------------------------------------------------------------------
    ! Basic stuff
    ! ---------------------------------------------------------------------------------------------------------

    integer, parameter :: nVar   = 4        ! Number of variables in the PDE system
    integer, parameter :: nDim   = 2        ! Number of dimensions
    real, parameter    :: CFL = 0.9 ! CFL condition < 1.0
    real, parameter    :: m_pi = 4.0*atan(1.0) 
    
    ! ---------------------------------------------------------------------------------------------------------
    ! Parameters controlling the simulation
    ! ---------------------------------------------------------------------------------------------------------

    integer            :: N                   ! Degree of approximation
    integer            :: IMAX, JMAX          ! Number of cells in each space dimension
    integer            :: timestep            ! Number of the current time step
    integer            :: WriteInterval       ! Number of time steps after which to write output data
    integer            :: bL, bR, bB, bT      ! Boundary condtions
    real               :: xL(nDim), xR(nDim)  ! Computational domain
    real               :: dx(nDim), dt        ! Vector of mesh spacings in each dimension and the time step
    real               :: time, tend          ! Current time  and final time
    real, allocatable  :: x(:,:,:)            ! Coordinates of cell centers
    integer            :: ICType              ! The test case to run
    character(len=200) :: BaseFile            ! Basic filename to write the results

    ! Supported boundary conditions: I-inflow, T-transmissive, R-reflective, P-periodic

    ! ---------------------------------------------------------------------------------------------------------
    ! Stencil related data
    ! ---------------------------------------------------------------------------------------------------------

    ! Basic stuff

    integer :: nDOF(0:nDim)           ! Number of degrees of freedom in each direction (0 is for time)
    integer :: nStencils              ! Number of stencils for WENO reconstruction
    integer :: nGhostCells            ! Number of ghost cells on each side of the domain to apply boundary conditions

    ! Arrays holding the solution, grid, rhs and fluxes

    real, allocatable  :: uh(:,:,:)         ! Solution vector (nVar, nCells+2*nGhostCells) -> also contains ghost values
    real, allocatable  :: duh(:,:,:)        ! RHS for updating the solution
    real, allocatable  :: wh(:,:,:)         ! Primitive variables in each cell
    real, allocatable  :: qBnd(:,:,:,:,:)   ! Boundary-extrapolated data for the state vector Q in the element
    real, allocatable  :: FBnd(:,:,:,:,:)   ! Boundary-extrapolated data for the normal flux F in the element
    real, allocatable  :: Fu(:,:,:)         ! Upwind flux in the x-direction
    real, allocatable  :: Gu(:,:,:)         ! Upwind flux in the y-direction

    ! Basis functions

    real, allocatable  :: xiGPN(:), wGPN(:)       ! Legendre-Gauss quadrature points in interval [0,1]
    real, allocatable  :: FLcoeff(:), FRcoeff(:)  ! Basis function values on the left and side of unit cell

    ! WENO related data

    real, parameter   :: lambda_s  = 1.0     ! Weight for the sided stencils
    real, parameter   :: lambda    = 100.0   ! Weight for the centered stencils
    real, allocatable :: lin_wt(:)           ! Normalized linear weights of each stencil
    real, parameter   :: small_num = 1.0e-12 ! Small number (to avoid division by zero)
    real, allocatable :: OS_M(:,:)           ! Oscillation indicator matrix
    real, allocatable :: iML(:,:), iMCL(:,:) ! Inverse coefficient matrix for left and center left stencil
    real, allocatable :: iMR(:,:), iMCR(:,:) ! Inverse coefficient matrix for right and center right stencil

    ! ADER related data

    real, allocatable :: Kxi(:,:)            ! Element stiffness matrix
    real, allocatable :: K1(:,:), iK1(:,:)   ! F1 - Ktau and its inverse
    real, allocatable :: F0(:)               ! Time flux matrices

    ! Important info and parameters concerning the governing PDE system

    type tEquations
        real :: GAMMA
    end type tEquations

    type(tEquations)   :: EQN

end module ader_weno

