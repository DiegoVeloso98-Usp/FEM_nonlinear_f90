module convergence_module
    use LinearMomentum_module
    use boundary_module
    implicit none
    
contains

    !==============================|
    !_ tolerance is set to 2.0d0   |
    !  so posicional_FEM iterates  |
    !  only once !!!               |
    !==============================|

subroutine Convergence_check(RES , F_int_previous , d_iter,ELEMS, NODES,VINC, npe, nelems, nnodes, ep, gaprox, NODES_1)
    implicit none
    real(8), intent(inout) :: F_int_previous(:)
    real(8), allocatable :: F_int_current(:)
    integer :: nnodes , npe , nelems , gaprox
    character(len=100), intent(in) :: ep
    real(8),intent(out) :: RES

    real, intent(in) :: ELEMS(:,:), VINC(:,:)
    real(8),intent(in) :: NODES(:,:)
    real(8), intent(inout) :: NODES_1(:,:)

    real(8), dimension(:), intent(in) :: d_iter

    real(8) :: dnrm2
    
    
    allocate(F_int_current(2*nnodes))
    
    call internalForce_Current_Iter(F_int_current , d_iter,ELEMS, NODES, npe, nelems, nnodes, ep, gaprox, NODES_1)

    call boundary_F_int(VINC, F_int_current)
    CALL boundary_F_int(VINC, F_int_previous)
    !===========================================================!
    !                                                           |
    !   RES = (F_int_current - F_int_previous) / F_int_current  |
    !                                                           |
    !===========================================================!

    ! RES =( dnrm2(2*nnodes, F_int_current, 1) - dnrm2(2*nnodes, F_int_previous, 1) ) !/ dnrm2(2*nnodes, F_int_previous, 1)
    RES = dnrm2(2*nnodes, F_int_current, 1) 

    ! write(*,*) "sum F_int_current", sum(F_int_current)
    ! write(*,*) "F_int_current", F_int_current
    ! stop

    RES = ABS(RES)

    deallocate(F_int_current)

end subroutine Convergence_check




subroutine internalForce_Current_Iter(FGLOBAL_int , d_iter,ELEMS, NODES, npe, nelems, nnodes, ep, gaprox, NODES_1)
        implicit none
    
        integer, intent(in) :: npe, nelems, nnodes, gaprox
        integer :: nph
        character(len=100), intent(in) :: ep
        real(8), allocatable :: hammer_points(:,:)
        real, intent(in) :: ELEMS(:,:)
        real(8),intent(in) :: NODES(:,:)
        real(8), intent(inout) :: NODES_1(:,:)
    
        ! Local variables
        integer :: iel, ih, ine, idir,  i,j,ino
        real(8) :: ksi, eta, peso
        real(8) :: result_fi(npe), result_dfidksi(npe), result_dfideta(npe)
        real(8) :: J0
        real:: h, bx, by
    
        
        
        real(8) :: cx_0(npe), cy_0(npe), cx_1(npe), cy_1(npe)
        real(8),intent(out) :: FGLOBAL_int(2*nnodes)
        real(8) :: A0(2,2), A1(2,2), A0_inv(2,2), A(2,2), C(2,2), E_green(2,2), S_Piola(2,2), Constitutive(4,4)
        real(8) :: E_gr_voigt(4), SigPiola_voigt(4), DE_ab(2,2), DE_1(2,2), DE_2(2,2)
        real(8) :: f_int_i
        real(8) :: DA1_ab(2,2)
        real(8) :: Identity(2,2)
        real(8) :: FLOCAL_int(2*npe)
        real(8),allocatable :: MATRIX_result_fi(:,:), MATRIX_result_dfidksi(:,:), MATRIX_result_dfideta(:,:)
    
    
        real(8), dimension(:), intent(in) :: d_iter
        real(8), dimension(npe) :: d_el
        real(8) :: d_i
    
    
        
        call hammer(nph, hammer_points)
        
        call shapeFunc_Matrix(npe,gaprox,MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta)
    


        ! print*, "test3"
        FGLOBAL_int = 0.0d0

        ! print*, "test4"
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED( FGLOBAL_int, NODES, NODES_1, MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta, ELEMS, npe, nnodes, nph, hammer_points, nelems, ep,d_iter)
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
    
            FLOCAL_int=0.0d0    ! Local force equal to zero for each element, but updated for every integration point
            do ih=1,nph
                ! print*, "test6.1"
                ksi=hammer_points(1,ih)
                eta=hammer_points(2,ih)
                peso=hammer_points(3,ih)
    
                ! print*, "test6.2"
                result_fi      = MATRIX_result_fi(ih, :)
                result_dfidksi = MATRIX_result_dfidksi(ih, :)
                result_dfideta = MATRIX_result_dfideta(ih, :)
    
                ! print*, "test7"
                ! A0(eq. 5.33) -> J0=det(A0)
                A0=0.0
                A1=0.0
                d_i=0.0d0
                do ino=1,npe
                    A0(1,1)=A0(1,1)+cx_0(ino)*result_dfidksi(ino)
                    A0(1,2)=A0(1,2)+cx_0(ino)*result_dfideta(ino)
                    A0(2,1)=A0(2,1)+cy_0(ino)*result_dfidksi(ino)
                    A0(2,2)=A0(2,2)+cy_0(ino)*result_dfideta(ino)
                    
                    ! print*, "test7.1"
                    ! A1 w/ Y1    (eq. 5.34)    
                    A1(1,1)=A1(1,1)+cx_1(ino)*result_dfidksi(ino)
                    A1(1,2)=A1(1,2)+cx_1(ino)*result_dfideta(ino)
                    A1(2,1)=A1(2,1)+cy_1(ino)*result_dfidksi(ino)
                    A1(2,2)=A1(2,2)+cy_1(ino)*result_dfideta(ino)
    
                    ! print*, "test7.1"
                    ! print*, "d_el(ino)", d_el(ino)
                    ! print*, "result_fi(ino)", result_fi(ino)
                    ! print*, "d_i", d_i
                    d_i = d_i + d_el(ino)*result_fi(ino)
                    ! print*, "test7.3"
                end do
    
                ! print*, "test8"
    
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
    
                ! print*, "test9"
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
    
                    end do ! i DIRECTIONS
                end do ! i ELEM NODES
            end do ! HAMMER POINTS
            
            !$OMP CRITICAL
            do ino=1,npe
                do idir=1,2
                    j=2*(int(ELEMS(iel,ino))-1)+idir
                    FGLOBAL_int(j)=FGLOBAL_int(j)+FLOCAL_int(2*(ino-1)+idir)    
                end do
            end do
            !$OMP END CRITICAL  
    
        end do ! ELEMENTS
        !$OMP END DO
        !$OMP END PARALLEL

    
    end subroutine internalForce_Current_Iter


subroutine boundary_F_int(VINC,  FGLOBAL)
    
        implicit none
    
        real, intent(in) :: VINC(:,:)
        real(8), intent(inout) ::  FGLOBAL(:)
        integer :: i, ino, idir, type_,ngl
        real(8):: value_ !, nlg
        
    
        ! Print*, "Applying boundary conditions"
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
                FGLOBAL(ngl) = 0.0d0
            end if
        end do
    
    end subroutine boundary_F_int
    
end module convergence_module