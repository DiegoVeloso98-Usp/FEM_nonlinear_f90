program main
    use read_data_module  ! Import the module containing the subroutine
    use hammer_module
    use shape_module
    ! use assembly_module
    ! use particle_generator_module
    use sparse
    use solver_module
    use writing_acadview_module


    implicit none

    ! Declare the variables to hold the data
    character(len=100) :: filename

    integer :: nnodes,nmats,nthick, nelem3, nelem6, nelem10, nelems,nforces,npressures,ndisplas, gaprox, npe
    real, allocatable ::  MATS(:,:) , THICKS(:), ELEMS(:,:), LOAD(:,:), VINC(:,:)
    real(8), allocatable :: NODES(:,:)


    character(len=100) :: ep


    real(8) :: hammer_points(3, 12)
    integer :: nph


    real(8), dimension(:,:), allocatable :: FFORMA, COORD




    type(sparse_matrix):: Hessian_sp
    type(sparse_matrix):: Q_sp




    real(8), allocatable :: U(:)
    real(8), allocatable :: d_Global(:)



    integer :: nLoad_steps
    real(8) :: tolerance
    real(8) :: analysis_time , damping , gamma , beta,tmax
    integer :: type
    

    ! real(8),dimension(:),allocatable::U
! =====================================================================================================================================================!
! gfortran -ffast-math -fopenmp -ffree-line-length-none sparse_set.f90 main.f90 read_data.f90 hammer.f90 shape.f90 solver.f90 writing_acadview.f90 -o program.out -L. -lsuperlu -lblasslu -llapack -lblas
! =====================================================================================================================================================!
! gfortran -O3 -ffast-math -fopenmp -ffree-line-length-none sparse_set.f90 main.f90 read_data.f90 hammer.f90 shape.f90 solver.f90 writing_acadview.f90 -o program.out -L. -lsuperlu -lblasslu -llapack -lblas
! =====================================================================================================================================================!
! =====================================================================================================================================================!
    ! Call the subroutine to read the input file and store values in the declared variables
    filename = "SIMPLES_10El.txt"
    call read_input_file(filename,nnodes, nmats, nthick, nelem3, nelem6, nelem10, &
                        nelems, nforces, npressures, ndisplas, &
                        NODES, MATS, THICKS, ELEMS, LOAD, VINC, &
                        gaprox, npe)



       
    ! Print the values
    print *, "Data read successfully"

    call hammer(nph, hammer_points)
    ! print *, hammer_points(1,:)


    ! Allocate FFORMA based on your needs (example: 3xsize_)
    allocate(FFORMA(npe,npe))
    allocate(COORD(npe,2))
    call shapeFunc(gaprox, npe, FFORMA, COORD)


    nLoad_steps=1600

    tolerance=1.0e-6

    allocate(U(2*nnodes))
    allocate(d_Global(nnodes))

    CALL prepare_to_use(Hessian_sp,2*nnodes,((2*npe*npe+npe)*nelems)) !
    CALL prepare_to_use( Q_sp , nnodes , npe**2+npe)


    ! EP = "PLANE_STRESS" or "PLANE_STRAIN"
    ep = "PLANE_STRESS"
    damping = 0.0d0
    analysis_time = 0.8d0
    gamma = 0.5d0
    beta = 0.25d0
    type = 3   ! TYPE OF LOADING (0 = INCREASE --> CONSTANT       )
                ! TYPE OF LOADING (1 = INCREASE --> DECREASE --> 0 )
                ! TYPE OF LOADING (2 = MAXIMUM  --> DECREASE --> 0 )
                ! TYPE OF LOADING (3 = CONSTANT                    )
                ! TYPE OF LOADING (4 = PULSE                       ) 
    tmax = 0.2d0 ! TIME WHERE THE LOAD IS IN ITS MAXIMUM VALUE
    call solverPosicional (analysis_time,damping,beta,gamma,tmax,type,Hessian_sp , Q_sp , ELEMS, NODES, LOAD, VINC, FFORMA, gaprox, npe, nelems, nnodes, nph, hammer_points, ep, nLoad_steps, tolerance,U , d_Global )


    ! call appending(ELEMS, NODES, U, nelems, nnodes, npe, gaprox)


end program main
