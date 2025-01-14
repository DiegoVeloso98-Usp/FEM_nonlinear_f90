module sigma_module
    use shape_module
    use read_data_module
    use Algebra_module

    implicit none

    contains

subroutine cauchy_Stress(CAUCHY ,VonMises, ELEMS, NODES, npe, nelems, nnodes, ep,gaprox,  NODES_1,d_Global)
    implicit none

    integer, intent(in) :: npe, nelems, nnodes, gaprox
    integer :: nph
    character(len=100), intent(in) :: ep
    real(8), allocatable :: hammer_points(:,:)
    real, intent(in) :: ELEMS(:,:)
    real(8),intent(in) :: NODES(:,:)
    real(8), intent(inout) :: NODES_1(:,:)

    ! Local variables
    integer :: iel, ih, ine,  i,j,ino
    real(8) :: ksi, eta, peso
    real(8) :: result_fi(npe), result_dfidksi(npe), result_dfideta(npe)
    real(8) :: J0 , J_sig
    real:: h, bx, by

    
    real(8) :: cx_0(npe), cy_0(npe), cx_1(npe), cy_1(npe)
    real(8) :: A0(2,2), A1(2,2), A0_inv(2,2), A(2,2), C(2,2), E_green(2,2), S_Piola(2,2), Constitutive(4,4)
    real(8) :: E_gr_voigt(4), SigPiola_voigt(4)
    real(8) :: Identity(2,2)
    real(8),allocatable :: MATRIX_result_fi(:,:), MATRIX_result_dfidksi(:,:), MATRIX_result_dfideta(:,:)

    real(8) :: Cauchy_i(2,2)
    real(8), dimension(:), intent(in) :: d_Global
    real(8), dimension(npe) :: d_el
    real(8) :: d_i

    real(8),allocatable :: L_il(:,:)
    real(8),allocatable :: SigX(:), SigY(:), SigXY(:)
    real(8),allocatable :: v_aux_X(:), v_aux_Y(:), v_aux_XY(:)
    real(8),allocatable :: cauchy_x(:), cauchy_y(:), cauchy_xy(:)
    real(8),intent(inout) :: CAUCHY(nnodes,4), VonMises(nnodes)
    real(8),allocatable :: M_il(:,:), M_inv(:,:)
    
    ! print*, "test1"


        
    call hammer(nph, hammer_points)
    call shapeFunc_Matrix(npe,gaprox,MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta)



    allocate(L_il(nph,npe))    
    allocate(M_il(npe,npe), M_inv(npe,npe))

    !==============================================!
    !==============================================!
    allocate(v_aux_X(npe),v_aux_Y(npe),v_aux_XY(npe))
    allocate(cauchy_x(npe),cauchy_y(npe),cauchy_xy(npe))
    allocate(SigX(nph), SigY(nph), SigXY(nph))


    call shapeFunc_Matrix(npe,gaprox,MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta)

    do i=1,nph
        result_fi = MATRIX_result_fi(i, :)
        do j=1,npe
            L_il(i,j)=result_fi(j)
        end do
    end do

    M_il = matmul(transpose(L_il),L_il)

    call matrix_inverse(M_il, M_inv)

! print*, "test1"



    cauchy_x=0.0d0
    cauchy_y=0.0d0
    cauchy_xy=0.0d0


    CAUCHY=0.0d0
    do iel=1,nelems
        call materialProperties(ELEMS, iel, npe, ep, h, bx, by, Constitutive)
        !  print*, "d_iter: ", d_iter
        
        ! print*, "test1"
        ! COORDINATES OF NODES OF THE ELEMENT
        do ine = 1, npe
        cx_0(ine) = NODES(int(ELEMS(iel, ine)), 1)
        cy_0(ine) = NODES(int(ELEMS(iel, ine)), 2)

        cx_1(ine) = NODES_1(int(ELEMS(iel, ine)), 1)
        cy_1(ine) = NODES_1(int(ELEMS(iel, ine)), 2)

        d_el(ine) = d_Global(int(ELEMS(iel, ine)))

        end do

        ! allocate(SigX(nph), SigY(nph), SigXY(nph))
        do ih=1,nph

            ksi=hammer_points(1,ih)
            eta=hammer_points(2,ih)
            peso=hammer_points(3,ih)

            result_fi      = MATRIX_result_fi(ih, :)
            result_dfidksi = MATRIX_result_dfidksi(ih, :)
            result_dfideta = MATRIX_result_dfideta(ih, :)

            ! A0(eq. 5.33) -> J0=det(A0)
            A0=0.0
            A1=0.0
            d_i=0.0d0
            ! print*, "test2"
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

            ! print*, "test3"

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


            ! print*, "test4"

            E_green=0.5*(C-Identity)

            E_gr_voigt(1)=E_green(1,1)
            E_gr_voigt(2)=E_green(2,2)
            E_gr_voigt(3)=2*E_green(1,2)
            E_gr_voigt(4)=2*E_green(2,1)



            !===============================================================!
            !============ DEGRADATED CONSTITUTIVE RELATIONSHIP =============!
            !=================           AT2 MODEL         =================!
            !===============================================================!
            ! print*, "test5"
            SigPiola_voigt=((1.0d0 - d_i)**2)*matmul(Constitutive,E_gr_voigt)

            !===============================================================!
            !===============================================================!

            S_Piola(1,1)=SigPiola_voigt(1)
            S_Piola(2,2)=SigPiola_voigt(2)
            S_Piola(1,2)=SigPiola_voigt(3)
            S_Piola(2,1)=SigPiola_voigt(4)

            ! print*, "test6"
            J_sig = A(1,1)*A(2,2)-A(1,2)*A(2,1)
            Cauchy_i=(1.0d0/J_sig)*matmul(A,matmul(S_Piola,transpose(A)))

            SigX(ih)=Cauchy_i(1,1)
            SigY(ih)=Cauchy_i(2,2)
            SigXY(ih)=Cauchy_i(1,2)
        end do ! HAMMER POINTS
! 
        ! print*, "test7"
        v_aux_X=matmul(transpose(L_il),SigX)
        v_aux_Y=matmul(transpose(L_il),SigY)
        v_aux_XY=matmul(transpose(L_il),SigXY)



        
        cauchy_x = matmul(M_inv , v_aux_X)
        
        cauchy_y = matmul(M_inv , v_aux_Y)
        
        cauchy_xy = matmul(M_inv , v_aux_XY)

        
        do ino=1,npe
            CAUCHY(int(ELEMS(iel,ino)),1) = CAUCHY(int(ELEMS(iel,ino)),1) + cauchy_x(ino)
            CAUCHY(int(ELEMS(iel,ino)),2) = CAUCHY(int(ELEMS(iel,ino)),2) + cauchy_y(ino)
            CAUCHY(int(ELEMS(iel,ino)),3) = CAUCHY(int(ELEMS(iel,ino)),3) + cauchy_xy(ino)
            CAUCHY(int(ELEMS(iel,ino)),4) = CAUCHY(int(ELEMS(iel,ino)),4) + 1
        end do

        ! deallocate(v_aux_X,v_aux_Y,v_aux_XY)
        ! deallocate(SigX,SigY,SigXY)
        ! deallocate(cauchy_x,cauchy_y,cauchy_xy)
    end do ! ELEMENTS

    deallocate(v_aux_X,v_aux_Y,v_aux_XY)

    deallocate(SigX,SigY,SigXY)

    deallocate(cauchy_x,cauchy_y,cauchy_xy)

    deallocate(L_il,M_il,M_inv)

    call VonMises_Stress(CAUCHY, VonMises)
end subroutine cauchy_Stress

subroutine VonMises_Stress(CAUCHY, VonMises)
    real(8), intent(in) :: CAUCHY(:,:)
    integer :: i,n
    real(8),intent(inout) :: VonMises(:)
    n=size(CAUCHY,1)
    VonMises=0.0d0
    do i=1,n
        VonMises(i)=sqrt(0.5d0*((CAUCHY(i,1)-CAUCHY(i,2))/CAUCHY(i,4))**2 + 3*(CAUCHY(i,3)/CAUCHY(i,4))**2)
    end do  
end subroutine VonMises_Stress

subroutine Stress_EigenValues(CAUCHY, Cauchy_EigenValues)
    real(8), intent(in) :: CAUCHY(:,:)
    real(8) :: cauchy_i(1,3)
    real(8),intent(inout) :: Cauchy_EigenValues(:,:)
    integer :: i,n
    real(8) :: eigen1, eigen2
    n=size(CAUCHY,1)
    Cauchy_EigenValues=0.0d0
    do i=1,n
        cauchy_i(1,:)=CAUCHY(i,3)/CAUCHY(i,4)
        !      DIRECT FORMULA FOR EIGENVALUES OF 2X2 MATRIX : lambda
        eigen1 = 0.5d0*(cauchy_i(1,1)+cauchy_i(1,2)) + sqrt((0.5d0*(cauchy_i(1,1)-cauchy_i(1,2)))**2 + (cauchy_i(1,3))**2) 
        eigen2 = 0.5d0*(cauchy_i(1,1)+cauchy_i(1,2)) - sqrt((0.5d0*(cauchy_i(1,1)-cauchy_i(1,2)))**2 + (cauchy_i(1,3))**2 )

        Cauchy_EigenValues(i,1) = eigen1
        Cauchy_EigenValues(i,2) = eigen2
        ! print*, "Cauchy_i: ",cauchy_i 
    end do

    ! stop
end subroutine Stress_EigenValues

end module sigma_module