module LinearMomentum_module
    use shape_module
    use read_data_module
    use sparse
    use writing_acadview_module
    use Boundary_module
    use hammer_module
    USE phaseField_module
    implicit none



contains

subroutine posicional_FEM(FGLOBAL_int , Hessian_sp,d_iter,ELEMS, NODES,  VINC,  npe, nelems, nnodes, ep, tolerance,gaprox, FGLOBAL_ext, NODES_1,NODES_1_vec, delta_Y,iter)
    implicit none

    integer, intent(in) :: npe, nelems, nnodes, gaprox
    integer :: nph
    character(len=100), intent(in) :: ep
    real(8) ::  tolerance 
    real(8), allocatable :: hammer_points(:,:)
    real, intent(in) :: ELEMS(:,:),  VINC(:,:)
    real(8),intent(in) :: NODES(:,:)
    real(8) :: delta_Y(:)
    real(8), intent(inout) :: NODES_1(:,:),NODES_1_vec(:)

    ! Local variables
    integer :: iel, ih, ine, idir, jdir, i,j,jno,ino
    integer, intent(out):: iter
    real(8) :: peso, error
    real(8) :: result_fi(npe), result_dfidksi(npe), result_dfideta(npe)
    real(8) :: J0
    real:: h, bx, by

    
    real(8) ::  Hessian_local(2*npe, 2*npe)
    real(8) :: cx_0(npe), cy_0(npe), cx_1(npe), cy_1(npe)
    real(8) :: FGLOBAL_ext(2*nnodes), g_residue(2*nnodes)
    real(8),intent(out) :: FGLOBAL_int(2*nnodes)
    real(8) :: A0(2,2), A1(2,2), A0_inv(2,2), A(2,2), C(2,2), E_green(2,2), S_Piola(2,2), Constitutive(4,4)
    real(8) :: E_gr_voigt(4), SigPiola_voigt(4), DE_ab(2,2), DS_gz(2,2), DE_1(2,2), DE_2(2,2), DS_gz_voigt(4)
    real(8) :: D2E(2,2), h_abgz, f_int_i,  DE_gz(2,2), D2E_1(2,2), D2E_2(2,2)
    real(8) :: DE_gz_voigt(4)
    real(8) :: DA1_ab(2,2), DA1_gz(2,2)
    real(8) :: Identity(2,2)
    real(8) :: FLOCAL_int(2*npe)
    real(8),allocatable :: MATRIX_result_fi(:,:), MATRIX_result_dfidksi(:,:), MATRIX_result_dfideta(:,:)


    real(8), dimension(:), intent(in) :: d_iter
    real(8), dimension(npe) :: d_el
    real(8) :: d_i

    type(sparse_matrix) :: Hessian_sp   !
    integer::indexes_sparse(2*npe)   

    real(8) :: dnrm2

    ! INTEGER,allocatable:: IR(:), JC(:)
    ! REAL(8),allocatable :: PA(:)
    ! INTEGER :: n
    


call hammer(nph, hammer_points)

call shapeFunc_Matrix(npe,gaprox,MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta)

iter=0
error=1.0d0
do while (error > tolerance)

    do i=1,nnodes
        if (d_iter(i) >= 1.0d0) then
            print*, "ERROR: d_iter(i) >= 1.0d0"
        end if
    end do
    FGLOBAL_int = 0.0d0

    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Hessian_sp, FGLOBAL_int, NODES, NODES_1, MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta, ELEMS, FGLOBAL_ext, npe, nnodes, nph, hammer_points, nelems, ep,d_iter)
    !$OMP DO
    do iel=1,nelems
        call materialProperties(ELEMS, iel, npe, ep, h, bx, by, Constitutive)

        ! COORDINATES OF NODES OF THE ELEMENT
        do ine = 1, npe
        cx_0(ine) = NODES(int(ELEMS(iel, ine)), 1)
        cy_0(ine) = NODES(int(ELEMS(iel, ine)), 2)

        cx_1(ine) = NODES_1(int(ELEMS(iel, ine)), 1)
        cy_1(ine) = NODES_1(int(ELEMS(iel, ine)), 2)

        d_el(ine) = d_iter(int(ELEMS(iel, ine)))

        end do

        Hessian_local=0.0d0    ! Local Hessian equal to zero for each element, but updated for every integration point
        FLOCAL_int=0.0d0    ! Local force equal to zero for each element, but updated for every integration point
        h_abgz=0.0d0
        do ih=1,nph

            peso=hammer_points(3,ih)

            result_fi      = MATRIX_result_fi(ih, :)
            result_dfidksi = MATRIX_result_dfidksi(ih, :)
            result_dfideta = MATRIX_result_dfideta(ih, :)

            ! A0(eq. 5.33) -> J0=det(A0)
            A0=0.0
            A1=0.0
            d_i=0.0d0
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

                d_i = d_i + d_el(ino)*result_fi(ino)
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



            !===============================================================!
            !============ DEGRADATED CONSTITUTIVE RELATIONSHIP =============!
            !=================           AT2 MODEL         =================!
            !===============================================================!

            SigPiola_voigt=((1.0d0 - d_i)**2)*matmul(Constitutive,E_gr_voigt)

            !===============================================================!
            !===============================================================!

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

                    ! print*, "test12"
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

                            
                            DS_gz_voigt=((1.0d0 - d_i)**2)*matmul(Constitutive,DE_gz_voigt)
                            
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

                        end do ! j DIRECTIONS
                    end do ! j ELEM NODES
                end do ! i DIRECTIONS
            end do ! i ELEM NODES
        end do ! HAMMER POINTS
        
        !$OMP CRITICAL
        do ino=1,npe
            do idir=1,2
                j=2*(int(ELEMS(iel,ino))-1)+idir
                FGLOBAL_int(j)=FGLOBAL_int(j)+FLOCAL_int(2*(ino-1)+idir)
                indexes_sparse(2*(ino-1)+idir)=j

            end do
        end do
        !$OMP END CRITICAL  


        !$OMP CRITICAL
        call add_matrix(Hessian_sp, Hessian_local, indexes_sparse, 2*npe)
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

    ! print*, "Solving system"
    call solve_system_of_equation(Hessian_sp, g_residue, delta_Y )
    ! print*, "System solved"

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
    ! print*, ""
    ! print*, "| Error: ", error, "|    Step: ",step, "|    Iteration: ", iter, " |"
    ! print*, "============================================"

    ! error=0.0d0
    if (error /= error) then
        print*, "ERROR: 'error' is NaN"
        stop
    end if
    

    ! Check for divergence
    if (error > 10000) then
        print*, "ERROR: DIVERGENCE"
        stop
    end if

    ! n = size(FGLOBAL_ext)
    ! call sparse_COO_to_CCS(Hessian_sp, n, PA, IR, JC)

    call clear_data(Hessian_sp)

    ! write(*,*) "FGLOBAL_int: ", FGLOBAL_int(1100:1200)
    ! STOP

end do ! RESIDUAL

end subroutine posicional_FEM
    
end module LinearMomentum_module


