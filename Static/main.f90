program main
    use read_data_module  ! Import the module containing the subroutine
    use hammer_module
    use shape_module
    ! use assembly_module
    ! use particle_generator_module
    use sparse
    use solver_module
    use writing_acadview_module
    use sigma_module


    implicit none

    ! Declare the variables to hold the data
    character(len=100) :: filename

    integer :: nnodes,nmats,nthick, nelem3, nelem6, nelem10, nelems,nforces,npressures,ndisplas, gaprox, npe
    real, allocatable ::  MATS(:,:) , THICKS(:), ELEMS(:,:), LOAD(:,:), VINC(:,:)
    real(8), allocatable :: NODES(:,:)


    character(len=100) :: ep

    
    type(sparse_matrix):: Hessian_sp
    type(sparse_matrix):: Q_sp




    real(8), allocatable :: U(:)
    real(8), allocatable :: d_Global(:)
    real(8) :: l_0, G_c



    integer :: nLoad_steps
    real(8) :: tolerance

    

    ! real(8),dimension(:),allocatable::U
! =====================================================================================================================================================!
! =====================================================================================================================================================!
! gfortran -O3 -ffast-math -fopenmp -ffree-line-length-none sparse_set.f90 main.f90 read_data.f90 hammer.f90 shape.f90 Boundary.f90 LinearMomentum.f90 Algebra.f90 sigma.f90 phase-field.f90 Convergence.f90 solver.f90 writing_acadview.f90 -o program.out -L. -lsuperlu -lblasslu -llapack -lblas
! gfortran -O0 -fcheck=all -finit-real=snan -ffpe-trap=invalid,zero,overflow -fbacktrace -Wall -fopenmp -ffree-line-length-none sparse_set.f90 main.f90 read_data.f90 hammer.f90 shape.f90 Boundary.f90 LinearMomentum.f90 Algebra.f90 sigma.f90 phase-field.f90 Convergence.f90 solver.f90 writing_acadview.f90 -o program.out -L. -lsuperlu -lblasslu -llapack -lblas

    ! ====================================================!
    ! ====================================================!
    filename = "PHASEFIELD.txt"
    ! EP = "PLANE_STRESS" or "PLANE_STRAIN"
    ep = "PLANE_STRESS"
    l_0 = 0.01d0
    G_c = 2.7d0
    ! G_c = 0.0003d0
    ! ====================================================!
    ! ====================================================!

    call read_input_file(filename,nnodes, nmats, nthick, nelem3, nelem6, nelem10, &
                        nelems, nforces, npressures, ndisplas, &
                        NODES, MATS, THICKS, ELEMS, LOAD, VINC, &
                        gaprox, npe)


    ! Print the values
    print *, "Data read successfully"



    nLoad_steps=20
    tolerance=1.0e-6

    allocate(U(2*nnodes))
    allocate(d_Global(nnodes))

    CALL prepare_to_use(Hessian_sp,2*nnodes,((2*npe*npe+npe)*nelems)) !
    CALL prepare_to_use( Q_sp , nnodes , ((npe**2+npe)/2)*nelems)


    call MainSolver (Hessian_sp , Q_sp , ELEMS, NODES, LOAD, VINC, gaprox, npe, nelems, nnodes, ep, nLoad_steps, tolerance,U , d_Global, l_0, G_c)

    

    ! call appending(ELEMS, NODES, U, nelems, nnodes, npe, gaprox)


end program main
