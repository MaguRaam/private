! Lagrange polynomial basis

subroutine lagrange_basis_values(center, xi, phi)
    use ader_weno
    implicit none
    ! Argument list
    integer, intent(in) :: center
    real, intent(in)    :: xi
    real, intent(out)   :: phi(0:N)
    ! Local variables
    integer :: i, k
    real    :: lagrange_weight, tmp_lagrange_weight, v, k_faculty

    tmp_lagrange_weight = 1.0
    k_faculty = 1.0

    ! First find the denominator

    do i = 1, N+1
        IF (i .NE. center) THEN
            tmp_lagrange_weight = tmp_lagrange_weight*(xiGPN(center) - xiGPN(i))
        END IF
    end do

    lagrange_weight = 1./tmp_lagrange_weight

    phi(0) = 1.0d0

    do i = 1, N
        phi(i) = 0.0d0
    end do

    do i = 0, N
        if (i+1 .ne. center) then
            v = xi - xiGPN(i+1)
            do k = N, 1, -1
                phi(k) = (phi(k)*v + phi(k - 1));
            end do
            phi(0) = phi(0)*v;
        end if 
    end do

    do k = 0, N
        phi(k) = phi(k)*k_faculty*lagrange_weight
        k_faculty = k_faculty*real(k+1)
    end do

end subroutine lagrange_basis_values


subroutine BaseFunc1D(phi,phi_xi,xi)
   use ader_weno
   implicit none
   ! Argument list
   real, intent(in ) :: xi                              ! coordinate in [0,1] where to evaluate the basis
   real, intent(out) :: phi(N+1), phi_xi(N+1)           ! the basis and its derivative w.r.t. xi
   ! Local variables
   integer :: i,j,m
   real    :: tmp
   
   ! Initialize variables
   phi      = 1.0
   phi_xi   = 0.0
   
   ! Lagrange polynomial and its derivative
   do m = 1, N+1
      do j = 1, N+1
         if (j .eq. m) cycle
         phi(m) = phi(m)*(xi-xiGPN(j))/(xiGPN(m)-xiGPN(j))
      end do
      do i = 1, N+1
         if (i .eq. m) cycle
         tmp = 1.0
         do j = 1, N+1
            IF (j .eq. i) CYCLE
            IF (j .eq. m) CYCLE
            tmp = tmp*(xi-xiGPN(j))/(xiGPN(m)-xiGPN(j))
         end do
         phi_xi(m) = phi_xi(m) + tmp/(xiGPN(m)-xiGPN(i))
      end do
   end do
   
end subroutine BaseFunc1D
