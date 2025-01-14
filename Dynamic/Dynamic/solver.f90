module solver_module
    use shape_module
    use read_data_module
    use sparse
    use writing_acadview_module
    ! use phaseField_module
    implicit none


contains

! subroutine phase_Field
    

subroutine solverPosicional(analysis_time,damping,beta,gamma,tmax,type,Hessian_sp , Q_sp , ELEMS, NODES, LOAD, VINC, FFORMA, gaprox, npe, nelems, nnodes, nph, hammer, ep, nLoad_steps, tolerance,U , d_Global )
    implicit none

    integer, intent(in) :: npe, nph, gaprox, nelems, nnodes, nLoad_steps,type
    character(len=100), intent(in) :: ep
    real(8), intent(in) ::  hammer(3, nph), FFORMA(:,:), tolerance 
    real, intent(in) :: ELEMS(:,:), LOAD(:,:), VINC(:,:)
    real(8) :: NODES(:,:),NODES_1(nnodes,2)
    real(8) :: delta_Y(2*nnodes),NODES_1_vec(2*nnodes),NODES_0_vec(2*nnodes)
    real(8), intent(inout) :: U(:)
    real(8), intent(inout) :: d_Global(:)
    real(8),intent(in) :: analysis_time,beta , gamma, damping,tmax
    ! Local variables
    integer ::  ih,  i, step,iter,iter_total,ino,jno,idir,ine
    real(8) :: ksi, eta


    real(8) ::  dt,time
    real(8) ::  FGLOBAL_ext(2*nnodes) , ACEL(2*nnodes), VEL(2*nnodes), Q_s(2*nnodes), R_s(2*nnodes)
    real(8), allocatable :: MATRIX_result_fi(:,:), MATRIX_result_dfidksi(:,:), MATRIX_result_dfideta(:,:)

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
    max_iterations = 1000

    

    allocate(MATRIX_result_fi(nph, npe))
    allocate(MATRIX_result_dfidksi(nph, npe))
    allocate(MATRIX_result_dfideta(nph, npe))


    MATRIX_result_fi = 0.0d0
    MATRIX_result_dfidksi = 0.0d0
    MATRIX_result_dfideta = 0.0d0


    do ih = 1, nph
        ksi = hammer(1, ih)
        eta = hammer(2, ih)
        MATRIX_result_fi(ih, :) = fi(FFORMA, ksi, eta, gaprox)
        MATRIX_result_dfidksi(ih, :) = dfidksi(FFORMA, ksi, eta, gaprox)
        MATRIX_result_dfideta(ih, :) = dfideta(FFORMA, ksi, eta, gaprox)
        
    end do


    do i=1,nnodes
        NODES_0_vec(2*i-1)=NODES(i,1)
        NODES_0_vec(2*i)=NODES(i,2)
    end do

    NODES_1=NODES
    NODES_1_vec=NODES_0_vec

    iter_total=0


    ACEL=0.0d0
    VEL=0.0d0
    dt=analysis_time/real(nLoad_steps,kind=8)
    U=0.0d0
    call initialize_file(nelems, nnodes, npe, gaprox, ELEMS, NODES)
    do step=1,nLoad_steps

        time=real(step,kind=8)*dt
        Q_s = 0.0d0
        R_s = 0.0d0

        Q_s = NODES_1_vec/(beta*dt**2) + VEL/(beta*dt) + ACEL*(1/(2*beta)-1)
        R_s = VEL + ACEL*dt*(1-gamma)

        ! Apply global forces
        FGLOBAL_ext = 0.0d0
        call globalForce(LOAD, FGLOBAL_ext,time,tmax,type,step)
        ! call prescribedDisplacement(VINC, NODES_1,nLoad_steps)

        !========================================!
        print*,"calling posicional_FEM"
        call Plot_Displacement_At_Load(U,LOAD,time,analysis_time)

        CALL posicional_FEM(Q_s , R_s ,VEL,ACEL,beta, gamma, dt, Hessian_sp,ELEMS, NODES,  VINC,  npe, nelems, nnodes, nph, hammer, ep, Step, tolerance, FGLOBAL_ext, NODES_1,NODES_1_vec, MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta,delta_Y,iter,damping)
        
        ACEL=0.0d0
        VEL=0.0d0
        ACEL = NODES_1_vec/(beta*dt**2) - Q_s
        VEL = NODES_1_vec*gamma/(beta*dt) + R_s - gamma * dt * Q_s 
        ! print*,"NODES_1_vec after",NODES_1_vec
        ! stop

        !========================================================================================!
        !== AFTER SOLUTION FROM posicional_FEM , THE NEW NODES POSITION, NODES_1 IS PASSED TO ===!
        !== phase_Field SUBROUTINE TO CALCULATE THE PHASE-FIELD VECTOR d_Global               ===!
        !========================================================================================!
        
        ! CALL phase_Field(Q_sp , d_n , ELEMS , NODES , hammer, NODES_1 , l_0 , G_c , MATRIX_result_fi , MATRIX_result_dfidksi , MATRIX_result_dfideta , npe , nnodes , nelems , nph , ep)
        
        iter_total=iter_total+iter



        !_ DISPLACEMENT AT STEP _!
        U=NODES_1_vec-NODES_0_vec
        !========================!
        call appending(U, nnodes)

        ! write (stepstr, fmt) step ! converting integer to string using a 'internal file'
        ! OPEN(99, ACTION = 'write', FILE = 'displacement_'//trim(stepstr)//'.nml')
        ! write(nml=RESTART, unit=99)
        ! CLOSE(99)

    end do ! LOAD STEPS



    print*, "Solution converged in ", iter_total, " iterations"
    print*, "============================================"
    call close_file()
end subroutine solverPosicional

subroutine posicional_FEM(Q_s , R_s , VEL, ACEL,beta, gamma, dt,Hessian_sp,ELEMS, NODES,  VINC,  npe, nelems, nnodes, nph, hammer, ep, Step, tolerance, FGLOBAL_ext, NODES_1,NODES_1_vec, MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta,delta_Y,iter, damping)

    integer, intent(in) :: npe, nph,  nelems, nnodes, Step
    character(len=100), intent(in) :: ep
    real(8), intent(in) ::  hammer(:,:),  tolerance 
    real, intent(in) :: ELEMS(:,:),  VINC(:,:)
    real(8),intent(in) :: NODES(:,:)
    real(8) :: delta_Y(:)
    real(8), intent(inout) :: NODES_1(:,:),NODES_1_vec(:)

    real(8), intent(in) :: Q_s(:), R_s(:)
    real(8), intent(IN) :: VEL(:), ACEL(:)
    real(8), intent(in) :: beta, gamma, dt, damping

    ! Local variables
    integer :: iel, ih, ine, idir, jdir, i,j,jno,ino
    integer, intent(out):: iter
    real(8) :: ksi, eta, peso, error
    real(8) :: result_fi(npe), result_dfidksi(npe), result_dfideta(npe)
    real(8) :: J0
    real:: h, bx, by , density
    real(8) :: MFIFI(2*npe,2*npe), mfi(npe,npe), MASS_el(2*npe,2*npe) , C_damping(2*npe,2*npe) , H_local(2*npe,2*npe)
    real(8) :: NODES_1_vec_elem(2*npe), Q_s_elem(2*npe), R_s_elem(2*npe)
    integer :: idx(2 * npe)
    
    real(8) :: Hessian(2*nnodes, 2*nnodes), Hessian_local(2*npe, 2*npe)
    real(8) :: cx_0(npe), cy_0(npe), cx_1(npe), cy_1(npe)
    real(8) :: FGLOBAL_int(2*nnodes), FGLOBAL_ext(2*nnodes), g_residue(2*nnodes)
    real(8) :: A0(2,2), A1(2,2), A0_inv(2,2), A(2,2), C(2,2), E_green(2,2), S_Piola(2,2), Constitutive(4,4)
    real(8) :: E_gr_voigt(4), SigPiola_voigt(4), DE_ab(2,2), DS_gz(2,2), DE_1(2,2), DE_2(2,2), DS_gz_voigt(4)
    real(8) :: D2E(2,2), h_abgz, f_int_i,  DE_gz(2,2), D2E_1(2,2), D2E_2(2,2)
    real(8) :: DE_gz_voigt(4)
    real(8) :: DA1_ab(2,2), DA1_gz(2,2)
    real(8) :: Identity(2,2)
    real(8) :: eh1(2*npe,2*npe),eh2(2*npe,2*npe),FLOCAL_int(2*npe),F_int_dyn(2*npe)
    real(8) :: MATRIX_result_fi(:,:), MATRIX_result_dfidksi(:,:), MATRIX_result_dfideta(:,:)

    type(sparse_matrix) :: Hessian_sp   !
    integer::indexes_sparse(2*npe)   

    real(8) :: dnrm2
    

iter=0
error=1.0d0



do while (error > tolerance)



    Hessian = 0.0d0
    FGLOBAL_int = 0.0d0
    eh1=0.0d0
    eh2=0.0d0 
    
    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Q_s,R_s,gamma,beta,dt,NODES_1_vec,Hessian_sp, FGLOBAL_int, NODES, NODES_1, MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta, ELEMS, FGLOBAL_ext, npe, nnodes, nph, hammer, nelems, ep,damping)
    !$OMP DO
    do iel=1,nelems
        call materialProperties(ELEMS, iel, npe, ep, h, bx, by, density , Constitutive)
        h_abgz=0.0d0
        ! COORDINATES OF NODES OF THE ELEMENT
        do ine = 1, npe
        cx_0(ine) = NODES(int(ELEMS(iel, ine)), 1)
        cy_0(ine) = NODES(int(ELEMS(iel, ine)), 2)

        cx_1(ine) = NODES_1(int(ELEMS(iel, ine)), 1)
        cy_1(ine) = NODES_1(int(ELEMS(iel, ine)), 2)
        end do

        
        do i = 1, npe
            idx(2 * i - 1) = 2 * (int(ELEMS(iel, i)) - 1) + 1
            idx(2 * i)     = 2 * (int(ELEMS(iel, i)) - 1) + 2
        end do
        NODES_1_vec_elem = NODES_1_vec(idx)
        Q_s_elem = Q_s(idx)
        R_s_elem = R_s(idx)
        ! print*, "Q_s", Q_s
        ! print*, "NODES_1_vec_elem", NODES_1_vec_elem
        ! print*, "Q_s_elem", Q_s_elem
        ! print*, "R_s_elem", R_s_elem
        ! STOP

        
        Hessian_local=0.0d0    ! Local Hessian equal to zero for each element, but updated for every integration point
        FLOCAL_int=0.0d0    ! Local force equal to zero for each element, but updated for every integration point
        F_int_dyn=0.0d0
        MASS_el=0.0d0
        H_local=0.0d0
        do ih=1,nph
            
            ksi=hammer(1,ih)
            eta=hammer(2,ih)
            peso=hammer(3,ih)


            result_fi      = MATRIX_result_fi(ih, :)
            result_dfidksi = MATRIX_result_dfidksi(ih, :)
            result_dfideta = MATRIX_result_dfideta(ih, :)


            ! A0(eq. 5.33) -> J0=det(A0)
            A0=0.0
            A1=0.0
            do ino=1,npe
                A0(1,1)=A0(1,1)+cx_0(ino)*result_dfidksi(ino)
                A0(1,2)=A0(1,2)+cx_0(ino)*result_dfideta(ino)
                A0(2,1)=A0(2,1)+cy_0(ino)*result_dfidksi(ino)
                A0(2,2)=A0(2,2)+cy_0(ino)*result_dfideta(ino)
                

                ! A1 w/ Y1    (eq. 5.34)    
                A1(1,1)=A1(1,1)+cx_1(ino)*result_dfidksi(ino)
                A1(1,2)=A1(1,2)+cx_1(ino)*result_dfideta(ino)
                A1(2,1)=A1(2,1)+cy_1(ino)*result_dfidksi(ino)
                A1(2,2)=A1(2,2)+cy_1(ino)*result_dfideta(ino)
            end do



            J0=A0(1,1)*A0(2,2)-A0(1,2)*A0(2,1)                   

            ! A0_inv , A1  ->  A  ->  C  ->  E_green(eq. 5.44)  -->  SigPiola(eq. 5.66)
            A0_inv(1,1) =  A0(2,2) / J0
            A0_inv(1,2) = -A0(1,2) / J0
            A0_inv(2,1) = -A0(2,1) / J0
            A0_inv(2,2) =  A0(1,1) / J0
            A=0.0
            A=matmul(A1,A0_inv)
            C=0.0
            C=matmul(transpose(A),A)
            ! Identity matrix
            Identity=0.0
            do i=1,2
                Identity(i,i)=1
            end do

            E_green=0.5*(C-Identity)

            E_gr_voigt(1)=E_green(1,1)
            E_gr_voigt(2)=E_green(2,2)
            E_gr_voigt(3)=2*E_green(1,2)
            E_gr_voigt(4)=2*E_green(2,1)


            SigPiola_voigt=matmul(Constitutive,E_gr_voigt)
            S_Piola(1,1)=SigPiola_voigt(1)
            S_Piola(2,2)=SigPiola_voigt(2)
            S_Piola(1,2)=SigPiola_voigt(3)
            S_Piola(2,1)=SigPiola_voigt(4)

            DA1_ab=0.0
            DE_1=0.0
            DE_2=0.0
            DE_ab=0.0
            do ino=1,npe
                do idir=1,2 

                    ! DA1_alpha,beta  (eq. 5.72)
                    DA1_ab(1,1)=result_dfidksi(ino)*(2-idir)
                    DA1_ab(1,2)=result_dfideta(ino)*(2-idir)
                    DA1_ab(2,1)=result_dfidksi(ino)*(idir-1)
                    DA1_ab(2,2)=result_dfideta(ino)*(idir-1)

                    
                    ! DE_alpha,beta (eq. 5.73)
                    DE_1=matmul(transpose(A0_inv),matmul(transpose(DA1_ab),A))

                    
                    DE_2=matmul(transpose(A),matmul(DA1_ab,A0_inv))

                    
                    DE_ab=0.5*(DE_1+DE_2)


                    ! fl_alpha,beta (eq. 5.74)
                    f_int_i=0.0d0
                    do i=1,2
                        do j=1,2
                            f_int_i=f_int_i+h*(DE_ab(i,j)*S_Piola(i,j))
                        end do
                    end do

                    FLOCAL_int(2*(ino-1)+idir)=FLOCAL_int(2*(ino-1)+idir)+f_int_i*peso*J0

                    do jno=1,npe
                        do jdir=1,2

                            ! DA_gamma,zeta (eq. 5.72)
                            DA1_gz(1,1)=result_dfidksi(jno)*(2-jdir)
                            DA1_gz(1,2)=result_dfideta(jno)*(2-jdir)
                            DA1_gz(2,1)=result_dfidksi(jno)*(jdir-1)
                            DA1_gz(2,2)=result_dfideta(jno)*(jdir-1)
                            
                            
                            
                            ! DE_gamma,zeta (eq. 5.73)
                            DE_1=matmul(transpose(A0_inv),matmul(transpose(DA1_gz),A))

                            
                            DE_2=matmul(transpose(A),matmul(DA1_gz,A0_inv))

                            
                            DE_gz=0.5*(DE_1+DE_2)



                            ! DS_alpha,beta (eq. 5.85)
                            DE_gz_voigt(1)=DE_gz(1,1)
                            DE_gz_voigt(2)=DE_gz(2,2)
                            DE_gz_voigt(3)=2*DE_gz(1,2)
                            DE_gz_voigt(4)=2*DE_gz(2,1)

                            
                            DS_gz_voigt=matmul(Constitutive,DE_gz_voigt)
                            
                            DS_gz=0.0

                            DS_gz(1,1)=DS_gz_voigt(1)
                            DS_gz(2,2)=DS_gz_voigt(2)
                            DS_gz(1,2)=DS_gz_voigt(3)
                            DS_gz(2,1)=DS_gz_voigt(4)


                            D2E_1=0.0
                            D2E_2=0.0
                            D2E=0.0

                            D2E_1=matmul(transpose(A0_inv),matmul(transpose(DA1_ab),matmul(DA1_gz,A0_inv)))
                            D2E_2= matmul(transpose(A0_inv),matmul(transpose(DA1_gz),matmul(DA1_ab,A0_inv)))
                            D2E=0.5*(D2E_1+D2E_2)

                            ! h_alpha,beta,gamma,zeta (eq. 5.88)
                    
                            h_abgz=0.0d0
                            do i=1,2
                                do j=1,2
                                    h_abgz=h_abgz+h*(DE_ab(i,j)*DS_gz(i,j)+S_Piola(i,j)*D2E(i,j))   

                                end do
                            end do


                            Hessian_local(2*(ino-1)+idir,2*(jno-1)+jdir)=Hessian_local(2*(ino-1)+idir,2*(jno-1)+jdir)+h_abgz*peso*J0
                                    


                            !=====================!
                            !=====================!
                            
                        end do ! j DIRECTIONS
                    end do ! j ELEM NODES
                end do ! i DIRECTIONS
            end do ! i ELEM NODES

                !======================================!
                !===========  MASS MATRIX  ============!
                !======================================!

                !==  OUTER PRODUCT OF VECTOR FI  ===!
            mfi=0.0d0
            do ino=1,npe
                do jno=1,npe
                    mfi(ino,jno)=MATRIX_result_fi(ih,ino)*MATRIX_result_fi(ih,jno)
                end do
            end do

            MFIFI=0.0d0
            do ino=1,npe
                do jno=1,npe
                    do idir=1,2
                        MFIFI(2*(ino - 1)+idir,2*(jno-1)+idir)=mfi(ino,jno)
                    end do
                end do
            end do
            ! print*, "MFIFI"
            ! do i=1,2*npe
            !     print*, MFIFI(i,:)
            ! end do
            ! stop
            MASS_el=MASS_el + MFIFI * h * peso * J0 * density
            
            
            ! print*, "density", density
            ! print*, "C_damping"
            ! do i=1,2*npe
            !     print*, C_damping(i,:)
            ! end do
            ! print*, "MASS_el"
            ! do i=1,2*npe
            !     print*, MASS_el(i,:)
            ! end do
            ! stop
            

            



        end do ! HAMMER POINTS

        !$OMP CRITICAL
        C_damping = damping * MASS_el
        H_local=Hessian_local + MASS_el/(beta*dt**2) + C_damping*(gamma/(beta*dt))

        F_int_dyn =  (1/(beta*dt**2))*matmul(MASS_el,NODES_1_vec_elem)-matmul(MASS_el,Q_s_elem)+(gamma/(beta*dt))*matmul(C_damping,NODES_1_vec_elem)+matmul(C_damping,R_s_elem)-(gamma*dt)*matmul(C_damping,Q_s_elem)
        ! print*, "FLOCAL_int", FLOCAL_int

        FLOCAL_int = FLOCAL_int + F_int_dyn


        ! if (iter==3) then
        !     PRINT*, "h ", h, "peso ", peso, "  J0  ", J0,"  density  ", density


        !     print*, "MASS_el"
        !     do i=1,2*npe
        !         print*, MASS_el(i,:)
        !     end do

        !     print*, "C_damping"
        !     do i=1,2*npe
        !         print*, C_damping(i,:)
        !     end do

        !     Print*, "NODES_1_vec_elem", NODES_1_vec_elem
        !     Print*, "Q_s_elem", Q_s_elem
        !     Print*, "R_s_elem", R_s_elem
        !     STOP
        ! end if
  
        !!$OMP CRITICAL
        do ino=1,npe
            do idir=1,2
                j=2*(int(ELEMS(iel,ino))-1)+idir
                FGLOBAL_int(j)=FGLOBAL_int(j)+FLOCAL_int(2*(ino-1)+idir)
                indexes_sparse(2*(ino-1)+idir)=j

            end do
        end do
        !!$OMP END CRITICAL  


        !!$OMP CRITICAL
        call add_matrix(Hessian_sp, H_local, indexes_sparse, 2*npe)
        !$OMP END CRITICAL  

    end do ! ELEMENTS
    !$OMP END DO
    !$OMP END PARALLEL
    
    g_residue=0.0d0
    g_residue=-FGLOBAL_int+FGLOBAL_ext


    ! SOLVE SYSTEM
    call assemble_sparse_matrix(Hessian_sp,timeit=.false.)
    call boundaryCondition(VINC, Hessian_sp, g_residue)

    ! Solve the system for me for matrix Hessian and vector g
    delta_Y = 0.0d0
    delta_Y = g_residue

    print*, "Solving system"
    call solve_system_of_equation(Hessian_sp, g_residue, delta_Y )
    print*, "System solved"


    ! print*, "delta_Y", delta_Y
    ! stop
    do i = 1, nnodes
        do j=1,2
            NODES_1(i,j)=NODES_1(i,j)+delta_Y(2*(i-1)+j)
        end do
    end do

    ! CURRENT POSITION VECTOR
    do i=1,nnodes
        NODES_1_vec(2*i-1)=NODES_1(i,1)
        NODES_1_vec(2*i)=NODES_1(i,2)
    end do
    
    error = dnrm2(2*nnodes, delta_Y, 1) / dnrm2(2*nnodes, NODES_1_vec, 1)
    
    iter=iter+1
    print*, ""
    print*, "| Error: ", error, "|    Step: ",step, "|    Iteration: ", iter, " |"
    print*, "============================================"
    ! if (iter > 5) then
    !     print*, "Solution did not converge"
    !     print*,"FGLOBAL",FGLOBAL_ext
    !     print*,"g_residue",g_residue
    !     print*,"FGLOBAL_int",FGLOBAL_int
    !     PRINT*, "NODES_1_vec", NODES_1_vec
    !     PRINT*, "MASS_el"
    !     DO i=1,2*npe
    !         PRINT*, MASS_el(i,:)
    !     END DO
    !     PRINT*, "Q_s_elem", Q_s_elem
    !     PRINT*,"NODES_1_vec_elem",NODES_1_vec_elem
    !     stop
    ! end if

    call clear_data(Hessian_sp)

end do ! RESIDUAL


end subroutine posicional_FEM






subroutine globalForce(LOAD, FGLOBAL,t,tmax,type,step)
    implicit none
    real, intent(in) :: LOAD(:,:)
    real(8), intent(inout) :: FGLOBAL(:)
    integer :: i, ino, idir
    real(8), intent(in) :: t, tmax
    integer, intent(in) :: type,step
    real:: value
     

    ! print*, "Applying global forces"
    do i = 1, size(LOAD, 1)
        ino=int(LOAD(i, 1))
        idir=int(LOAD(i, 2))
        value=LOAD(i, 3)


        if (type==0) then                                       !  (P)
            if (t<tmax) then                                    !  ^
                FGLOBAL(2*(ino-1)+idir)=value*t/tmax            !  |        *******
            else                                                !  |      *
                FGLOBAL(2*(ino-1)+idir)=value                   !  |    *
            end if                                              !  |  *
                                                                !  |*____________>(t)


        else if (type==1) then                                  !  (P)
            if (t<tmax) then                                    !  ^
                FGLOBAL(2*(ino-1)+idir)=value*(t)/(tmax)        !  |        *
            else if (tmax<=t .and. t<=2*tmax) then                !  |      *   *
                FGLOBAL(2*(ino-1)+idir)=value*(2*tmax-t)/(tmax) !  |    *       *
            else if (t>2*tmax) then                             !  |  *           *
                FGLOBAL(2*(ino-1)+idir)=0.0d0                   !  |*              *******
                                                                !  |_______________________>(t)
            end if                                              

        else if (type==2) then                                      !  (P)
            if (t<tmax) then                                        !  ^
                FGLOBAL(2*(ino-1)+idir)=value                       !  |*******
            else if (t>tmax) then                                   !  |        *
                FGLOBAL(2*(ino-1)+idir)=value*(1-(t-tmax)/(tmax))    !  |          *
            else                                                    !  |            *
                FGLOBAL(2*(ino-1)+idir)=0.0d0                       !  |              *******
            end if                                                  !  |_______________________>(t)


        else if(type==3) then                                       !  (P)
                                                                    !  ^
            FGLOBAL(2*(ino-1)+idir)=value                           !  |********************
                                                                    !  |        
                                                                    !  |          
                                                                    !  |_______________________>(t)


        else if (type==4) then                                       !  (P)
            if (step==1) then                                        !  ^
                FGLOBAL(2*(ino-1)+idir)=value                        !  |***
            else                                                     !  |  *
                FGLOBAL(2*(ino-1)+idir)=0.0d0                        !  |  *
            end if                                                   !  |  ******************** 
                                                                     !  |_______________________>(t)
        end if

    end do
end subroutine globalForce

subroutine boundaryCondition(VINC, KGLOBAL, FGLOBAL)
    real, intent(in) :: VINC(:,:)
    real(8), intent(inout) ::  FGLOBAL(:)
    type(sparse_matrix), intent(inout):: KGLOBAL
    integer :: i, ino, idir, type_,ngl
    real(8):: value_ !, nlg
    

    Print*, "Applying boundary conditions"
    do i = 1, size(VINC, 1)
        ino=int(VINC(i, 1))
        idir=int(VINC(i, 2))
        value_=VINC(i, 3)
        type_=int(VINC(i, 4))
        
        !_________________!
        ngl=2*(ino-1)+idir
        !_________________!

        ! RIGID SUPPORT !
        if (type_==1) then
            FGLOBAL(ngl) = 0.d0
            call set_value_to_row(KGLOBAL,ngl,0.0d0)
            call set_value_to_col(KGLOBAL,ngl,0.0d0)
            
            call set_value_in_term(KGLOBAL,ngl,ngl,1.0d0)

        ! ELASTIC SUPPORT !
        else
            ! I dont think this is right, row and column probably wont be the same value "ngl"
            call sum_value_in_term(KGLOBAL,ngl,ngl,value_)

        end if
        
    end do
    Print*, "Boundary conditions applied ! "
end subroutine boundaryCondition

subroutine prescribedDisplacement(VINC, NODES_1,nLoad_steps)
    real, intent(in) :: VINC(:,:)
    real(8), intent(inout) ::  NODES_1(:,:)
    integer, intent(in) :: nLoad_steps
    integer :: i, ino, idir, type_
    real(8):: value_ !, nlg
    

    Print*, "Applying Prescribed Displacement"

    do i = 1, size(VINC, 1)
        ino=int(VINC(i, 1))
        idir=int(VINC(i, 2))
        value_=VINC(i, 3)
        type_=int(VINC(i, 4))
        

        ! RIGID SUPPORT !
        if (type_==1) then
            NODES_1(ino,idir) = NODES_1(ino,idir) + (1.0d0/real(nLoad_steps,kind=8))*value_

        end if
        
    end do
    Print*, "Prescribed Displacements Applied ! "
end subroutine prescribedDisplacement

end module solver_module