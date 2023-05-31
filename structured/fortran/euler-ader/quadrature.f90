! quadrature.f90
! author: sunder

!-----------------------------------------------------------------------
! Various quadrature routines
!-----------------------------------------------------------------------

pure subroutine leggauss(x1, x2, x, w, n)
    use ader_weno, ONLY : m_pi
    implicit none
    ! Argument list
    integer, intent(in) ::  n
    real, intent(in)    :: x1, x2
    real, intent(out)   :: x(n), w(n)
    ! Local variables
    real :: EPS
    integer  :: i,j,m
    real :: p1, p2, p3, pp, xl, xm, z, z1

    parameter (EPS=1.0e-15)

    m  = (n+1)/2
    xm = 0.5*(x2+x1)
    xl = 0.5*(x2-x1)
    do i=1,m
        z = cos(m_pi*(i-0.25)/(n+0.5))
    1   continue 
        p1 = 1.0
        p2 = 0.0
        do j = 1,n
            p3 = p2
            p2 = p1
            p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/real(j)
        end do
        pp = real(n)*(z*p1-p2)/(z*z-1.0)
        z1 = z
        z  = z1-p1/pp
        if (abs(z-z1) .gt. EPS) goto 1
        x(i)    = xm-xl*z
        x(n+1-i)= xm+xl*z
        w(i)    = 2.0*xl/((1.-z*z)*pp*pp)
        w(n+1-i)= w(i)
    end do
    return
end subroutine leggauss

