FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fdefault-real-8
SRC=ader_weno.f90 \
    analyze_soln.f90 \
    basis.f90 \
    bconds.f90 \
    calc_time_step.f90 \
    compute_rhs.f90 \
    initial_field.f90 \
    initialize.f90 \
    la.f90 \
    pde.f90 \
    quadrature.f90 \
    read_input.f90 \
    space_time_predictor.f90 \
    write_data.f90 \
    main.f90

OBJ=${SRC:.f90=.o}

%.o: %.f90
	${FC} ${FFLAGS} -o $@ -c $<
hype: ${OBJ}
	${FC} ${FFLAGS} -o $@ ${OBJ}
clean: 
	@rm -f *.mod *.o hype
distclean: 
	@rm -f *.mod *.o *.vtk hype
