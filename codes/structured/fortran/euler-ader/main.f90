program ADERWENO
    use ader_weno
    implicit none
    ! Local variables
    integer :: i,j,iVar

    call read_input
    call initialize

    ! Main loop in time

    do while (time < tend)

        if (mod(timestep,WriteInterval) .eq. 0) then
            call write_data
        end if

        call calc_time_step
        call apply_boundary_conditions
        call compute_rhs

        do j = 1, JMAX
            do i = 1, IMAX
                do iVar = 1, nVar
                    uh(iVar,i,j) = uh(iVar,i,j) + dt*duh(iVar,i,j)
                end do
            end do
        end do

        print *, timestep, 'time = ', time

        time = time + dt
        timestep = timestep + 1

    end do

    print *, timestep, 'time = ', time

    call analyze_soln
    call write_data

end program ADERWENO
