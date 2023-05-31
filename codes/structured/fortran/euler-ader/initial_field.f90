subroutine initial_field(xGP, u0)
    use ader_weno
    implicit none
    ! Argument list
    real, intent(in ) :: xGP(nDim)    ! spatial position vector
    real, intent(out) :: u0(nVar)     ! initial data vector in terms of conserved variables

    ! Local variables
    real :: v0(nVar), rho0, prs0, vlx0, vly0, kappa, rr2, exp_rr2, tempaa, tempab, xx, yy

    xx = xGP(1)
    yy = xGP(2)

    if (ICType .eq. 1) then ! Sod Shock Tube
        if (xx .le. 0.5 ) then
            v0(1) = 1.0
            v0(2) = 0.0 
            v0(3) = 0.0 
            v0(4) = 1.0 
        else
            v0(1) = 0.125
            v0(2) = 0.0 
            v0(3) = 0.0 
            v0(4) = 0.1
        end if
    else if (ICType .eq. 2) then ! Lax Shock Tube
        if (xx .le. 0.5 ) then
            v0(1) = 0.445
            v0(2) = 0.698 
            v0(3) = 0.0 
            v0(4) = 3.528 
        else
            v0(1) = 0.5
            v0(2) = 0.0 
            v0(3) = 0.0 
            v0(4) = 0.571
        end if
    else if (ICType .eq. 3) then ! Shu-osher oscillatory shock tube
        if (xx .le. -4.0 ) then
            v0(1) = 3.857143
            v0(2) = 2.629369 
            v0(3) = 0.0 
            v0(4) = 10.33333
        else
            v0(1) = 1.0 + 0.2*sin(5.0*m_pi*xx)
            v0(2) = 0.0
            v0(3) = 0.0 
            v0(4) = 1.0
        end if
    else if (ICType .eq. 4) then ! Smooth Isentropic Vortex 
        rho0 = 1.0
        prs0 = 1.0
        vlx0 = 1.0
        vly0 = 1.0
        kappa = 5.0

        rr2 = xx**2 + yy**2
        exp_rr2 = exp ( 0.5 * ( 1.0 - rr2) )

        tempaa = - exp_rr2**2 * kappa**2 * ( EQN%GAMMA - 1.0)/ ( 8.0 * EQN%GAMMA * m_pi**2)
        tempab = tempaa + prs0 / rho0
        tempab = tempab * rho0**EQN%GAMMA / prs0

        v0(1) = tempab**( 1.0 / ( EQN%GAMMA - 1.0))
        v0(2) = vlx0 - yy * kappa * exp_rr2 / ( 2.0 * m_pi)
        v0(3) = vly0 + xx * kappa * exp_rr2 / ( 2.0 * m_pi)
        v0(4) = v0(1) * (tempaa + prs0 / rho0)
    else ! Degenerate problem 
        v0(1) = 1.0 + 0.2*sin(m_pi*(xx))
        v0(2) = 1.0
        v0(3) = 1.0
        v0(4) = 1.0
    end if
    
    if (xx .ge. 0.5 .and. yy .ge. 0.5) then
        v0(1) = 0.5197
        v0(2) = 0.1
        v0(3) = 0.1 
        v0(4) = 0.4
    else if (xx .le. 0.5 .and. yy .ge. 0.5) then 
        v0(1) = 1.0
        v0(2) = -0.6259
        v0(3) = 0.1
        v0(4) = 1.0
    else if (xx .LE. 0.5 .AND. yy .LE. 0.5) then
        v0(1) = 0.8
        v0(2) = 0.1
        v0(3) = 0.1 
        v0(4) = 1.0
    else
        v0(1) = 1.0
        v0(2) = 0.1
        v0(3) = -0.6259 
        v0(4) = 1.0
    end if

    call PDEPrim2Cons(v0, u0)

end subroutine initial_field
