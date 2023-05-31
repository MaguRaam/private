! Read the input data from the parameters file

subroutine read_input
    use ader_weno
    implicit none

    open(unit = 1, file = 'params.ini', status = 'unknown')

    read(1,*)EQN%GAMMA
    read(1,*)N
    read(1,*)xL(1)
    read(1,*)xL(2)
    read(1,*)xR(1)
    read(1,*)xR(2)
    read(1,*)IMAX
    read(1,*)JMAX
    read(1,*)bL
    read(1,*)bR
    read(1,*)bB
    read(1,*)bT
    read(1,*)tend
    read(1,*)WriteInterval
    read(1,*)ICType

    close(1)

end subroutine
