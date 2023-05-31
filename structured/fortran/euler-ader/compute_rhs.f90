subroutine os_indicator(ul, OS)
    use ader_weno, only : N, OS_M
    implicit none
    ! Argument list
    real, intent(in) :: ul(N+1)
    real, intent(out) :: OS
    ! Local variables
    integer :: p, m

    OS = 0.0

    do p = 1, N+1
        do m = 1, N+1
            OS = OS + OS_M(p,m)*ul(p)*ul(m)
        end do
    end do

end subroutine os_indicator

!-----------------------------------------------------------------------
! Given stencil in 1D, reconstruct the solution
!-----------------------------------------------------------------------

subroutine reconstruct_1d(lu, luh)
    use ader_weno
    implicit none
    ! Argument list
    real, intent(in) :: lu(-N:N)  ! Stencil containing average values
    real, intent(out) :: luh(N+1) ! Unknown coefficients to be found
    ! Local variables
    integer :: L, iDOF
    real :: ustencil(N+1)         ! Local stencil
    real, dimension(N+1) :: luh_l, luh_r, luh_cl, luh_cr
    real :: OS_L, OS_R, OS_CL, OS_CR
    real :: w(nStencils), total

    ! Left biased stencil

    L = -N;

    ustencil = lu(L:L+N)
    luh_l = matmul(iML, ustencil)
    call os_indicator(luh_l, OS_L)
    w(1) = lin_wt(1)/(OS_L + small_num)**4

    ! Right biased stencil

    L = 0

    ustencil = lu(L:L+N)
    luh_r = matmul(iMR, ustencil)
    call os_indicator(luh_r, OS_R)
    w(2) = lin_wt(2)/(OS_R + small_num)**4

    ! Center right stencil

    L = -floor(real(N/2))

    ustencil = lu(L:L+N)
    luh_cr = matmul(iMCR, ustencil)
    call os_indicator(luh_cr, OS_CR)
    w(3) = lin_wt(3)/(OS_CR + small_num)**4

    if (nStencils .eq. 4) then
        L = -(floor(real(N/2)) + 1)

        ustencil = lu(L:L+N)
        luh_cl = matmul(iMCL, ustencil)
        call os_indicator(luh_cl, OS_CL)
        w(4) = lin_wt(4)/(OS_CL + small_num)**4

    end if

    ! Normalize the weights

    total = sum(w)
    w = w/total

    do idoF = 1, N+1
        luh(iDOF) = w(1)*luh_l(iDOF) + w(2)*luh_r(iDOF) + w(3)*luh_cr(iDOF)
        if(nStencils .eq. 4) then
            luh(iDOF) = luh(iDOF) + w(4)*luh_cl(iDOF)
        end if
    end do

end subroutine reconstruct_1d


subroutine reconstruct_2d(lu, luh)
    use ader_weno
    implicit none
    ! Argument list
    real, intent(in) :: lu(-N:N,-N:N)  ! Stencil containing average values
    real, intent(out) :: luh(N+1,N+1) ! Unknown coefficients to be found
    ! Local variables
    real :: luhx(N+1, -N:N), ustencil(-N:N)
    integer :: i, j

    do j = -N,N
        do i = -N,N
            ustencil(i) = lu(i,j)
        end do
        call reconstruct_1d(ustencil, luhx(:,j))
    end do

    do i = 1,N+1
        ustencil = luhx(i,:)
        call reconstruct_1d(ustencil, luh(i,:))
    end do

end subroutine reconstruct_2d

subroutine compute_rhs

    use ader_weno
    implicit none

    real, dimension(nVar, ndoF(1), ndoF(2)):: uhi, qhi  ! Local solution in the cell
    real :: luh(N+1,N+1), lu(-N:N,-N:N)                          ! Local stencils for reconstruction
    integer :: i, j, k, l, iVar
    real :: nv(nDim)
    real :: Flux(nVar)

    do j = 0, JMAX+1
        do i = 0, IMAX+1
            do iVar = 1, nVar

                do l = -N, N
                    do k = -N, N
                        lu(k,l) = uh(iVar, i+k, j+l)
                    end do
                end do

                call reconstruct_2d(lu, luh)

                do l = 1, ndoF(2)
                    do k = 1, ndoF(1)
                        uhi(iVar,k,l) = luh(k,l)
                    end do
                end do

            end do

            call ader_space_time_predictor(qhi(:,:,:), qBnd(:,:,:,i,j), FBnd(:,:,:,i,j), uhi(:,:,:))

        end do
    end do

    ! Find upwind flux in x-direction

    nv(1) = 1.0; nv(2) = 0.0

    do j = 1, JMAX
        do i = 1, IMAX+1
            Fu(:,i,j) = 0.0
            do l = 1, nDOF(2)
                call PDELLFFLux(qBnd(:,2,l,i-1,j), FBnd(:,2,l,i-1,j), qBnd(:,1,l,i,j), FBnd(:,1,l,i,j), nv, Flux)
                Fu(:,i,j) = Fu(:,i,j) + wGPN(l)*Flux(:)
            end do
        end do
    end do

    ! Find upwind flux in y-direction

    nv(1) = 0.0; nv(2) = 1.0

    do j = 1, JMAX+1
        do i = 1, IMAX
            Gu(:,i,j) = 0.0
            do k = 1, nDOF(1)
                call PDELLFFLux(qBnd(:,4,k,i,j-1), FBnd(:,4,k,i,j-1), qBnd(:,3,k,i,j), FBnd(:,3,k,i,j), nv, Flux)
                Gu(:,i,j) = Gu(:,i,j) + wGPN(k)*Flux(:)
            end do
        end do
    end do

    ! Find the update coefficients

    do j = 1, JMAX
        do i = 1, IMAX
            do iVar = 1, nVar
                duh(iVar,i,j) = -(1.0/dx(1))*( Fu(iVar,i+1,j) - Fu(iVar,i,j) ) - (1.0/dx(2))*( Gu(iVar,i,j+1) - Gu(iVar,i,j) )
            end do
        end do
    end do

end subroutine compute_rhs

