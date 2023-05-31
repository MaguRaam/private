subroutine apply_boundary_conditions
    use ader_weno
    ! Local variables
    integer :: i, j, ilhs, irhs, iVar, oned_begin, oned_end

    oned_begin = 1
    oned_end = IMAX

    ! Boundary conditions on the left

    do j = 1-nGhostCells,JMAX+nGhostCells

        do i = 1, nGhostCells

            ilhs = oned_begin - i

            if (bL .eq. 1) then       ! Transmissive boundary
                irhs = oned_begin
                do iVar = 1, nVar
                    uh(iVar, ilhs, j) = uh(iVar, irhs, j)
                end do
            end if

            if (bL .eq. 2) then       ! Reflective boundary
                irhs = oned_begin - 1 + i
                do iVar = 1, nVar
                    uh(iVar, ilhs, j) = uh(iVar, irhs, j)
                end do
                uh(1, ilhs, j) = -uh(1, ilhs, j)
            end if

            if (bL .eq. 3) then       ! Periodic boundary
                irhs = oned_end + 1 - i
                do iVar = 1, nVar
                    uh(iVar, ilhs, j) = uh(iVar, irhs, j)
                end do
            end if
        end do
    end do

    ! Boundary conditions on the right

    do j = 1-nGhostCells,JMAX+nGhostCells

        do i = 1, nGhostCells

            ilhs = oned_end + i

            if (bR .eq. 1) then       ! Transmissive boundary
                irhs = oned_end
                do iVar = 1, nVar
                    uh(iVar, ilhs, j) = uh(iVar, irhs, j)
                end do
            end if

            if (bR .eq. 2) then       ! Reflective boundary
                irhs = oned_end - 1 + i
                do iVar = 1, nVar
                    uh(iVar, ilhs, j) = uh(iVar, irhs, j)
                end do
                uh(1, ilhs, j) = -uh(1, ilhs, j)
            end if

            if (bR .eq. 3) then       ! Periodic boundary
                irhs = oned_begin + i - 1
                do iVar = 1, nVar
                    uh(iVar, ilhs, j) = uh(iVar, irhs, j)
                end do
            end if
        end do
    end do

    oned_begin = 1
    oned_end = JMAX

    ! Boundary conditions on the bottom

    do j = 1, nGhostCells
        do i = 1-nGhostCells,IMAX+nGhostCells

            ilhs = oned_begin - j

            if (bB .eq. 1) then       ! Transmissive boundary
                irhs = oned_begin
                do iVar = 1, nVar
                    uh(iVar, i, ilhs) = uh(iVar, i, irhs)
                end do
            end if

            if (bB .eq. 2) then       ! Reflective boundary
                irhs = oned_begin - 1 + j
                do iVar = 1, nVar
                    uh(iVar, i, ilhs) = uh(iVar, i, irhs)
                end do
                uh(2, ilhs, j) = -uh(2, ilhs, j)
            end if

            if (bB .eq. 3) then       ! Periodic boundary
                irhs = oned_end + 1 - j
                do iVar = 1, nVar
                    uh(iVar, i, ilhs) = uh(iVar, i, irhs)
                end do
            end if
        end do
    end do

    ! Boundary conditions on the top

    do j = 1, nGhostCells

        do i = 1-nGhostCells,IMAX+nGhostCells

            ilhs = oned_end + j

            if (bR .eq. 1) then       ! Transmissive boundary
                irhs = oned_end
                do iVar = 1, nVar
                    uh(iVar, i, ilhs) = uh(iVar, i, irhs)
                end do
            end if

            if (bR .eq. 2) then       ! Reflective boundary
                irhs = oned_end - 1 + j
                do iVar = 1, nVar
                    uh(iVar, i, ilhs) = uh(iVar, i, irhs)
                end do
                uh(1, ilhs, j) = -uh(1, ilhs, j)
            end if

            if (bR .eq. 3) then       ! Periodic boundary
                irhs = oned_begin + j - 1
                do iVar = 1, nVar
                    uh(iVar, i, ilhs) = uh(iVar, i, irhs)
                end do
            end if
        end do
    end do


end subroutine apply_boundary_conditions
