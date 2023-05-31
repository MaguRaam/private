subroutine calc_time_step 
    use ader_weno
    implicit none
    
    ! Local variables
    
    integer  :: iDim, i, j
    real :: eig(nVar)
    real :: nv(nDim,nDim), denom
    
    ! Normal vectors pointing into the two space dimensions
    
    nv = 0.0
    
    do i = 1, nDim
        nv(i,i) = 1.0
    end do
    
    dt = 1.0e20

    do j = 1, JMAX
        do i = 1, IMAX

            denom = 0.0
            do iDim = 1, nDim
                call PDEEigenvalues(uh(:,i,j),nv(:,iDim),eig)
                denom = denom + maxval(abs(eig))/dx(iDim)
            end do

            dt = min(dt, CFL/denom )

        end do
    end do
    
    if (time+dt>tend) then
        dt = tend - time
    end if
    
end subroutine calc_time_step
