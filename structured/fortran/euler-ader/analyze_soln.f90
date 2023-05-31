subroutine analyze_soln
    USE ader_weno
    implicit none
    ! Local variables
    integer :: i, j, ii, jj
    real :: q0(nVar)
    real :: corner(nDim), xGP(nDim), error(IMAX, JMAX), density

    do j = 1, JMAX
        do i = 1, IMAX

            density = 0.0
            corner(:) = x(:,i,j) - 0.5*dx(:) ! coordinate of the lower left corner of the cell
            do jj = 1, nDOF(2)
                do ii = 1, nDOF(1)
                    xGP = corner + (/ xiGPN(ii), xiGPN(jj) /)*dx(:)
                    call initial_field(xGP, q0)
                    density = density + wGPN(ii)*wGPN(jj)*q0(1)
                end do
            end do

            error(i,j) = abs(uh(1,i,j) - density)

        end do
    end do

    print *, maxval(error)

end subroutine analyze_soln
