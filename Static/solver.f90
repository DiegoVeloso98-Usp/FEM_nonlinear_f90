module solver_module
    use shape_module
    use read_data_module
    use sparse
    use writing_acadview_module
    use LinearMomentum_module
    use sigma_module
    use phaseField_module
    use convergence_module
    implicit none


contains

! subroutine phase_Field
    

subroutine MainSolver(Hessian_sp , Q_sp , ELEMS, NODES, LOAD, VINC, gaprox, npe, nelems, nnodes, ep, nLoad_steps, tolerance,U , d_n, l_0, G_c)
    implicit none

    integer, intent(in) :: npe, gaprox, nelems, nnodes, nLoad_steps
    character(len=100), intent(in) :: ep
    real(8), intent(in) :: tolerance
    real, intent(in) :: ELEMS(:,:), LOAD(:,:), VINC(:,:)
    real(8) :: NODES(:,:),NODES_1(nnodes,2)
    real(8) :: delta_Y(2*nnodes),NODES_1_vec(2*nnodes),NODES_0_vec(2*nnodes)
    real(8), intent(inout) :: U(:)


    real(8), intent(in) :: l_0, G_c
    real(8), intent(inout) :: d_n(:)
    real(8) :: d_i(nnodes) , delta_d_i(nnodes) , q_v(nnodes)

    real(8) :: CAUCHY(nnodes,4), VonMises(nnodes), Cauchy_EigenValues(nnodes,2)


    ! Local variables
    integer ::  i, step,iter,iter_total, staggered

    real(8) ::  RES,tolerance_stagg
    
    real(8) ::  FGLOBAL_ext(2*nnodes), F_int_previous(2*nnodes)

    integer, allocatable :: indsort(:)
    real(8), allocatable:: PA(:)
    integer(4), allocatable :: IR(:), JC(:)

    type(sparse_matrix) :: Hessian_sp   !
    type(sparse_matrix) :: Q_sp   !

    !============================================!
    ! VARIABLES FOR WRITING NAMLIST AT EACH STEP !
    !============================================!
    ! character(len=8) :: fmt='(I5.5)' ! format descriptor for writing integers to string
    ! character(len=8) :: stepstr
    !============================================!

    integer :: max_iterations

    namelist /RESTART/ U
    
    ! call phaseField(Q_sp , q_v, d_i , ELEMS , NODES , NODES_1 , l_0 , G_c , gaprox , npe , nnodes , nelems ,  ep)
    ! n=nnodes
    ! call sparse_COO_to_CCS(Q_sp, n, PA,IR,)

    CAUCHY=0.0d0
    



    !==================================================!
    !=== FIRST SOLUTION TRY IS THE INITIAL POSITION ===!
    !==================================================!
    NODES_1=NODES
    U=0.0d0
    
    !==================================================!
    
    do i=1,nnodes
        NODES_0_vec(2*i-1)=NODES(i,1)
        NODES_0_vec(2*i)=NODES(i,2)
    end do

    max_iterations = 1000
    iter_total=0

    call initialize_file(nelems, nnodes, npe, gaprox, ELEMS, NODES)
    d_n=0.0d0
    ! write(*,*) "test1"
    call initialize_sparse_COO_to_CCS(Q_sp,nnodes, IR, JC, indsort, d_n , ELEMS , NODES , NODES_1 , l_0 , G_c , gaprox , npe , nnodes , nelems ,  ep)
    ! write(*,*) "test2"
    do step=1,nLoad_steps
        ! Apply global forces

        FGLOBAL_ext = 0.0d0
        call globalForce(LOAD, FGLOBAL_ext,step,nLoad_steps)
        call prescribedDisplacement(VINC, NODES_1,nLoad_steps)

        !========================================!

        ! call Plot_Displacement_At_Load(U,LOAD,VINC,step,nLoad_steps)

        delta_d_i=0.0d0
        RES=1.0d0
        d_i=0.0d0
        d_i=d_n+delta_d_i
        ! d_i=0.0d0
        staggered=0
        tolerance_stagg=1.0e-4
        DO WHILE (RES > tolerance_stagg)
            staggered=staggered+1



            call posicional_FEM(F_int_previous , Hessian_sp , d_i , ELEMS , NODES,  VINC,  npe, nelems, nnodes, ep, tolerance,gaprox, FGLOBAL_ext, NODES_1,NODES_1_vec, delta_Y,iter)
            !========================================================================================!
            !== AFTER SOLUTION FROM posicional_FEM , THE NEW NODES POSITION, NODES_1 IS PASSED TO ===!
            !==   phase_Field SUBROUTINE TO ASSEMBLE THE PHASE-FIELD MATRIX Q_sp and VECTOR q_v    ==!
            !========================================================================================!
            ! print*, "-----------------------------------------------------"
            ! print*, "|  Linear Momentum Solved, Assembling Phase Field   |"
            ! print*, "-----------------------------------------------------"
            
            call phaseField(Q_sp , q_v, d_n , ELEMS , NODES , NODES_1 , l_0 , G_c , gaprox , npe , nnodes , nelems ,  ep)


            ! write(*,*) "-----------------------------------------------------"
            ! write(*,*) "Sum q_v = ", sum(q_v)
            ! write(*,*) "q_v = ", q_v
            ! stop

            ! write(*,*) "test3"
            call PSOR_Solver(nnodes ,Q_sp, q_v, delta_d_i, IR , JC , indsort)
            ! write(*,*) "test4"


            ! print*, "PSOR Solved"

            call clear_data(Q_sp)

            ! ! print*, "Calling DirichletBoundaryCondition"

            d_i = d_n + delta_d_i

            call DirichletBoundaryCondition(d_n , d_i , delta_d_i )


            
            ! ! print*, "Dirichlet Boundary Condition Applied "
            
            
            
            ! print*, "Sum Hessian_sp = ", sum(Hessian_sp%val)
            ! print*, "Sum Q_sp = ", sum(Q_sp%val)
            ! print*, "Sum q_v = ", sum(q_v)
            ! print*, "sum of d_i = ", sum(d_i)
            ! ! print*, "Checking Convergence"
            
            call Convergence_check(RES , F_int_previous , d_i ,ELEMS, NODES,VINC, npe, nelems, nnodes, ep, gaprox, NODES_1)

            ! print*, "Convergence Checked"
            
            iter_total=iter_total+iter
            
            write(*,*) "-----------------------------------------------------"
            write(*,*) "Step:",step,"  |  Iteration: ", staggered, "  |  Residual: ", RES

            ! RES=0.0d0

        end do ! RESIDUAL

        !== AFTER CONVERGENCE OF d_i, THE SOLUTION d_n IS UPDATED ===! 
                d_n=d_i                                              !
        !============================================================!
        CAUCHY=0.0d0


        call cauchy_Stress(CAUCHY , VonMises , ELEMS , NODES , npe, nelems, nnodes, ep,gaprox,  NODES_1,d_n)
        call Stress_EigenValues(CAUCHY,Cauchy_EigenValues)
        call appending(CAUCHY , VonMises , Cauchy_EigenValues , U , d_n ,  nnodes)
        print*, "Cauchy Stress Calculated"
        ! if (mod(step, 20) == 0) then
        !     call appending(CAUCHY , VonMises , Cauchy_EigenValues , U , d_n ,  nnodes)
        ! end if

        ! call appending(CAUCHY , VonMises , Cauchy_EigenValues , U , nnodes)
         

        !_ DISPLACEMENT AT STEP _!
        U=NODES_1_vec-NODES_0_vec
        !========================!

        ! write (stepstr, fmt) step ! converting integer to string using a 'internal file'
        ! OPEN(99, ACTION = 'write', FILE = 'displacement_'//trim(stepstr)//'.nml')
        ! write(nml=RESTART, unit=99)
        ! CLOSE(99)

    end do ! LOAD STEPS

    call close_file()





    print*, "Solution converged in ", iter_total, " iterations"
    print*, "============================================"
end subroutine MainSolver



end module solver_module