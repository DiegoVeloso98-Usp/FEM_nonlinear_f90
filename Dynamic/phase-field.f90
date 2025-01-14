module phaseField_module
    use shape_module
    use read_data_module
    use sparse
    implicit none


contains

subroutine phaseField(Q_sp , d_n , ELEMS , NODES , hammer, NODES_1 , l_0 , G_c , MATRIX_result_fi , MATRIX_result_dfidksi , MATRIX_result_dfideta , npe , nnodes , nelems , nph , ep)
    integer, intent(in) :: npe , nnodes , nelems , nph
    character(len=100) :: ep
    real, dimension(:,:), intent(in) :: ELEMS
    real(8), dimension(:,:), intent(in) :: hammer 
    real(8), dimension(:,:), intent(in) :: NODES, NODES_1
    real(8), dimension(:,:), intent(in) :: MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta


    !=====  LOCAL VARIABLES  ======!
    integer :: iel, ih, ino, i, j , ine
    real :: bx, by, h
    real(8) :: J0 , peso
    real(8),dimension(4,4) :: Constitutive
    real(8), dimension(2,2) :: A1 , A0 , A0_inv , A , C , E_green , SigPiola , Identity
    real(8), dimension(4) :: E_gr_voigt , SigPiola_voigt
    real(8), dimension(npe) :: cx_0, cy_0, cx_1, cy_1
    real(8), dimension(npe) :: result_fi, result_dfidksi, result_dfideta
    

    !==== PHASE FIELD VARIABLES =====!
    real(8) :: Energy_density
    type(sparse_matrix) :: Q_sp
    real(8), dimension(:), intent(in) :: d_n
    real(8), dimension(npe) :: d_el
    real(8), dimension(npe) :: q_v_el
    real(8), intent(in) :: l_0, G_c
    real(8), dimension(npe,npe) :: PSI_el , FI_el , Q_el
    real(8), dimension(nnodes) :: q_v
    real(8), dimension(2,npe) :: B_de
    real(8), dimension(1,npe) :: N_de
    integer, dimension(npe) :: indexes_sparse


    do iel=1,nelems
        call materialProperties(ELEMS, iel, npe, ep, h, bx, by, Constitutive)
        ! COORDINATES OF NODES OF THE ELEMENT
        do ine = 1, npe
        cx_0(ine) = NODES(int(ELEMS(iel, ine)), 1)
        cy_0(ine) = NODES(int(ELEMS(iel, ine)), 2)

        cx_1(ine) = NODES_1(int(ELEMS(iel, ine)), 1)
        cy_1(ine) = NODES_1(int(ELEMS(iel, ine)), 2)

        d_el(ine) = d_n(int(ELEMS(iel, ine)))
        end do
        
        ! Hessian_local=0.0d0    ! Local Hessian equal to zero for each element, but updated for every integration point
        ! FLOCAL_int=0.0d0    ! Local force equal to zero for each element, but updated for every integration point
        do ih=1,nph
            
            peso=hammer(3,ih)


            result_fi      = MATRIX_result_fi(ih, :)
            result_dfidksi = MATRIX_result_dfidksi(ih, :)
            result_dfideta = MATRIX_result_dfideta(ih, :)
            ! ASSEMBLY SHAPE MATRICES B_d,e and N_d,e !
            B_de(1,:) = result_dfidksi(:)
            B_de(2,:) = result_dfideta(:)
            N_de(1,:) = result_fi(:)


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

            SigPiola(1,1)=SigPiola_voigt(1)
            SigPiola(2,2)=SigPiola_voigt(2)
            SigPiola(1,2)=SigPiola_voigt(3)
            SigPiola(2,1)=SigPiola_voigt(4)

            do i=1,2
                do j=1,2
                    Energy_density=Energy_density+0.5*E_green(i,j)*SigPiola(i,j)
                end do
            end do

            !==================================================!
            !=  MATRIX PSI_el, FI_el, Q_el  (size npe x npe)  =!
            !=            VECTOR qv (size npe x 1)            =!
            !==================================================!

            PSI_el = PSI_el + 2*Energy_density*matmul(transpose(N_de),N_de)*J0*peso
            FI_el = FI_el + ((l_0**(-1))*matmul(transpose(N_de),N_de) + l_0*matmul(transpose(B_de),B_de))*J0*peso
            Q_el = Q_el + PSI_el + FI_el * G_c
            q_v_el = q_v_el + ( matmul(Q_el, d_el) - 2*Energy_density*transpose(N_de)*J0*peso )

            do ino=1,npe
                indexes_sparse(ino) = int(ELEMS(iel, ino))
            end do
        end do ! HAMMER POINTS

        call add_matrix(Q_sp, Q_el, indexes_sparse, 2*npe)

    end do ! ELEMENTS

    ! SOLVE SYSTEM
    call assemble_sparse_matrix(Q_sp,timeit=.false.)

    ! call boundaryCondition(VINC, Q_sp, q_v)

    
    ! call solve_system_of_equation(Q_sp, g_residue, delta_Y ) 

end subroutine phaseField



! subroutine PSOR_Solver ( Q_sp , q_v , d_Global)
!     ! PSOR SOLVER


!     n=size(q_v)
!     allocate(x_km1(n))
!     allocate(f(n))

!     x = 0.0d0
!     do iter=1, 100000
!         x_km1 = x

!         do jcol=1, n
!             Mx = 0.0d0
!             Ldx = 0.0d0
!             D_m1 = 0.0d0

!         end do 
            
!     end do

! end subroutine PSOR_Solver

end module phaseField_module
