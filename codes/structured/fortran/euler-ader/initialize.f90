subroutine initialize
    use ader_weno
    implicit none
    ! Local variables
    integer :: i, j, p, m, ii, jj, alpha, q, center, L, iGP, VMAX(nDim)
    real, dimension(0:N) :: phi_p, phi_m, phi
    real, dimension(N+1,N+1) :: MS, FRm
    real, dimension(nVar) :: q0
    real, dimension(N+1) :: phi_xi, phi_i, phiL, phiR
    real :: xiL, xiR, shift, corner(nDim), xGP(nDim), total

    ! Set other parameters to default/calculated values

    timestep = 0
    time = 0.0
    VMAX = (/ IMAX, JMAX /)
    dx = (xR-xL)/VMAX
    ndoF = N+1

    if (mod(N, 2) .eq. 0) then
        nStencils = 3
    else
        nStencils = 4
    end if

    nGhostCells = 2*(N+1) -1

    if (ICType .eq. 1) then
        BaseFile = 'SodShockTube'
    else if (ICType .eq. 2) then
        BaseFile = 'LaxShockTube'
    else if (ICType .eq. 3) then
        BaseFile = 'ShuOsherShockTube'
    else if (ICType .eq. 4) then
        BaseFile = 'SmoothVortex'
    else
        BaseFile = 'Sol'
    end if

    ! Allocate memory

    allocate ( x(nDim, IMAX, JMAX) )
    allocate ( uh(nVar, 1-nGhostCells:IMAX+nGhostCells, 1-nGhostCells:JMAX+nGhostCells) )
    allocate ( duh(nVar, IMAX, JMAX) )
    allocate ( wh(nVar, IMAX, JMAX) )
    allocate ( qBnd(nVar, 4, ndoF(2), 0:IMAX+1, 0:JMAX+1) )
    allocate ( FBnd(nVar, 4, ndoF(2), 0:IMAX+1, 0:JMAX+1) )
    allocate ( Fu(nVar, IMAX+1 , JMAX ) )
    allocate ( Gu(nVar, IMAX , JMAX+1 ) )

    ! Find the cell center coordinates

    do j = 1, JMAX
        do i = 1, IMAX
            x(1,i,j) = xL(1) + (real(i-1) + 0.5)*dx(1)
            x(2,i,j) = xL(2) + (real(j-1) + 0.5)*dx(2)
        end do
    end do

    ! Basis functions

    xiL = 0.0
    xiR = 1.0

    allocate(  xiGPN(N+1) )
    allocate(  wGPN(N+1) )
    allocate(  FLcoeff(N+1) )
    allocate(  FRcoeff(N+1) )

    call leggauss(xiL, xiR, xiGPN, wGPN, N+1)

    ! Find the basis function values at the left and the right end of the unit cell

    do i = 1, N+1
        call lagrange_basis_values(i, xiL, phi)
        phiL(i) = phi(0)
        call lagrange_basis_values(i, xiR, phi)
        phiR(i) = phi(0)
    end do

    FLcoeff = phiL
    FRcoeff = phiR

    ! WENO Related data
    
    allocate (lin_wt(nStencils)) ! Calculate linear weights

    lin_wt(1) = lambda_s
    lin_wt(2) = lambda_s
    lin_wt(3) = lambda

    if (nStencils .eq. 4) then
        lin_wt(4) = lambda
    end if

    total = sum(lin_wt)

    lin_wt = lin_wt/total

    allocate( OS_M(N+1, N+1) ) ! Initialize the oscillation indicator matrix

    OS_M = 0.0

    do p = 1, N+1
        do m = 1, N+1
            do q = 1, N+1
                call lagrange_basis_values(p, xiGPN(q), phi_p)
                call lagrange_basis_values(m, xiGPN(q), phi_m)
                do alpha = 1, N
                    OS_M(p,m) = OS_M(p,m) + wGPN(q)*phi_p(alpha)*phi_m(alpha)
                end do
            end do
        end do
    end do

    allocate( iMCL(N+1, N+1) ) ! Initialize center left stencil coefficient matrix
    allocate( iML(N+1, N+1) )  ! Initialize left-sided coefficient matrix
    allocate( iMCR(N+1, N+1) ) ! Initialize center right-sided coefficient matrix
    allocate( iMR(N+1, N+1) )  ! Initialize right-sided coefficient matrix

    ! Centered left side stencil

    MS = 0.0
    L = floor(real(N/2)) + 1

    do i = 1, N+1
        center = i - 1 - L
        do j = 1, N+1
            do q = 1, N+1
                shift = xiGPN(q) + real(center)
                call lagrange_basis_values(j, shift, phi)
                MS(i,j) = MS(i,j) + wGPN(q)*phi(0)
            end do
        end do
    end do

    call matrix_inverse(N+1,MS,iMCL)

    ! Centered right side stencil (Use this stencil for odd-ordered schemes)

    MS = 0.0
    L = floor(real(N/2))

    do i = 1, N+1
        center = i - 1 - L
        do j = 1, N+1
            do q = 1, N+1
                shift = xiGPN(q) + real(center)
                call lagrange_basis_values(j, shift, phi)
                MS(i,j) = MS(i,j) + wGPN(q)*phi(0)
            end do
        end do
    end do

    call matrix_inverse(N+1,MS,iMCR)

    ! Left sided stencil

    MS = 0.0
    L = N

    do i = 1, N+1
        center = i - 1 - L
        do j = 1, N+1
            do q = 1, N+1
                shift = xiGPN(q) + real(center)
                call lagrange_basis_values(j, shift, phi)
                MS(i,j) = MS(i,j) + wGPN(q)*phi(0)
            end do
        end do
    end do

    call matrix_inverse(N+1,MS,iML)

    ! Right sided stencil

    MS = 0.0
    L = 0

    do i = 1, N+1
        center = i - 1 - L
        do j = 1, N+1
            do q = 1, N+1
                shift = xiGPN(q) + real(center)
                call lagrange_basis_values(j, shift, phi)
                MS(i,j) = MS(i,j) + wGPN(q)*phi(0)
            end do
        end do
    end do

    call matrix_inverse(N+1,MS,iMR)

    ! ADER related data

    allocate( Kxi(N+1, N+1) ) ! Element stiffness matrix
    allocate( K1(N+1, N+1) ) ! Element stiffness matrix
    allocate( iK1(N+1, N+1) ) ! Element stiffness matrix
    allocate( F0(N+1) ) ! Time flux matrix

    Kxi = 0.0

    do iGP = 1, N+1
        call BaseFunc1D(phi_i,phi_xi,xiGPN(iGP))
        do i = 1, N+1
            do j = 1, N+1
                Kxi(i,j) = Kxi(i,j) + wGPN(iGP)*phi_xi(i)*phi_i(j)
            end do
        end do
    end do

    do i = 1, N+1
        do j = 1, N+1
            FRm(i,j) = phiR(i)*phiR(j)   ! Left contribution to the right flux matrix   (m = left  of the interface)
        end do
    end do

    ! The time flux matrices for the ADER-DG predictor method are given by the principle of upwinding in time (causality principle)

    K1 = FRm - Kxi ! upwinding in time = information comes from smaller times
    F0 = phiL   ! upwinding in time = information comes from smaller times

    call matrix_inverse(N+1,K1,iK1)

    ! Initialize the solution

    do j = 1, JMAX
        do i = 1, IMAX
            uh(:,i,j) = 0.0
            corner(:) = x(:,i,j) - 0.5*dx(:) ! coordinate of the lower left corner of the cell
            do jj = 1, nDOF(2)
                do ii = 1, nDOF(1)
                    xGP = corner + (/ xiGPN(ii), xiGPN(jj) /)*dx(:)
                    call initial_field(xGP, q0)
                    uh(:,i,j) = uh(:,i,j) + wGPN(ii)*wGPN(jj)*q0
                end do
            end do
        end do
    end do

end subroutine initialize
