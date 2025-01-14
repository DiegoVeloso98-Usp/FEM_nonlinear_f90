module phaseField_module
    use shape_module
    use read_data_module
    use sparse
    use hammer_module
    implicit none


contains

subroutine phaseField(Q_sp , q_v , d_n , ELEMS , NODES , NODES_1 , l_0 , G_c , gaprox , npe , nnodes , nelems ,  ep)
    implicit none
    
    integer :: nph
    integer, intent(in) :: npe , nnodes , nelems , gaprox
    character(len=100) :: ep
    real, dimension(:,:), intent(in) :: ELEMS
    real(8), dimension(:,:), allocatable :: hammer_points
    real(8), dimension(:,:), intent(in) :: NODES, NODES_1
    real(8), dimension(:,:), allocatable :: MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta


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
    real(8) :: d_i
    real(8), dimension(npe) :: q_v_el, q_v_el_aux
    real(8), intent(in) :: l_0, G_c
    real(8), dimension(npe,npe) :: PSI_el , FI_el , Q_el
    real(8), dimension(nnodes),intent(out) :: q_v
    real(8), dimension(2,npe) :: B_de_aux , B_de
    real(8), dimension(1,npe) :: N_de
    real(8) :: N_de_vec(npe) , B_de_de(npe,npe) , N_de_de(npe,npe)
    integer, dimension(npe) :: indexes_sparse

    call hammer(nph, hammer_points)
    call shapeFunc_Matrix(npe, gaprox, MATRIX_result_fi, MATRIX_result_dfidksi, MATRIX_result_dfideta)

    q_v=0.0d0
    do iel=1,nelems
        q_v_el=0.0d0
        q_v_el_aux=0.0d0
        Q_el=0.0d0
        call materialProperties(ELEMS, iel, npe, ep, h, bx, by, Constitutive)
        ! COORDINATES OF NODES OF THE ELEMENT
        do ine = 1, npe
        cx_0(ine) = NODES(int(ELEMS(iel, ine)), 1)
        cy_0(ine) = NODES(int(ELEMS(iel, ine)), 2)

        cx_1(ine) = NODES_1(int(ELEMS(iel, ine)), 1)
        cy_1(ine) = NODES_1(int(ELEMS(iel, ine)), 2)

        d_el(ine) = d_n(int(ELEMS(iel, ine)))
        end do
        

        do ih=1,nph
            
            peso=hammer_points(3,ih)


            result_fi      = MATRIX_result_fi(ih, :)
            result_dfidksi = MATRIX_result_dfidksi(ih, :)
            result_dfideta = MATRIX_result_dfideta(ih, :)
            ! ASSEMBLY SHAPE MATRICES B_d,e and N_d,e !
            B_de_aux(1,:) = result_dfidksi(:)
            B_de_aux(2,:) = result_dfideta(:)
            N_de(1,:) = result_fi(:)

            N_de_de = matmul(transpose(N_de),N_de)
            




            ! A0(eq. 5.33) -> J0=det(A0)
            A0=0.0d0
            A1=0.0d0
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


            ! print *, "d_i: ", d_i
            ! stop

            J0=A0(1,1)*A0(2,2)-A0(1,2)*A0(2,1)                   


            ! A0_inv , A1  ->  A  ->  C  ->  E_green(eq. 5.44)  -->  SigPiola(eq. 5.66)
            A0_inv(1,1) =  A0(2,2) / J0
            A0_inv(1,2) = -A0(1,2) / J0
            A0_inv(2,1) = -A0(2,1) / J0
            A0_inv(2,2) =  A0(1,1) / J0

            !===============================================================!
            !==  B_de_aux FROM ADIMENSIONAL SYSTEM ksi-eta TO x-y SYSTEM  ==!
            !================  B_de = A0_inv * B_de_aux  ===================!
            !===============================================================!

            B_de = matmul(A0_inv , B_de_aux)
            B_de_de = matmul(transpose(B_de),B_de)

            A=0.0
            A=matmul(A1,A0_inv)
            C=0.0
            C=matmul(transpose(A),A)
            ! Identity matrix
            Identity=0.0d0
            do i=1,2
                Identity(i,i)=1.0d0
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


            SigPiola_voigt=matmul(Constitutive,E_gr_voigt)


            !===============================================================!
            !===============================================================!

            SigPiola(1,1)=SigPiola_voigt(1)
            SigPiola(2,2)=SigPiola_voigt(2)
            SigPiola(1,2)=SigPiola_voigt(3)
            SigPiola(2,1)=SigPiola_voigt(4)


                    !== PIOLA-KIRCHHOFF TO CAUCHY ===!
            ! J_sig = A(1,1)*A(2,2)-A(1,2)*A(2,1)
            ! Cauchy_i=(1.0d0/J_sig)*matmul(A,matmul(SigPiola,transpose(A)))
            

            Energy_density=0.0
            do i=1,2
                do j=1,2
                    Energy_density=Energy_density+0.5*E_green(i,j)*SigPiola(i,j)
                end do
            end do

            !==================================================!
            !=  MATRIX PSI_el, FI_el, Q_el  (size npe x npe)  =!
            !=            VECTOR qv (size npe x 1)            =!
            !==================================================!
            !    Outer Product N_de_de = (N_de ⊗ N_de^T)      !
            

            N_de_vec = N_de(1,:)
            ! print *, "N_de_vec: ", N_de_vec
            ! stop


            PSI_el =  2*Energy_density*matmul(transpose(N_de),N_de) * J0 * peso * h
            FI_el = ((1.0d0 / l_0) * matmul(transpose(N_de), N_de) + l_0 * matmul(transpose(B_de), B_de)) * J0 * peso * h


            Q_el = Q_el + PSI_el + FI_el * G_c

            ! write(*,*) "peso",peso
            ! write(*,*) "J0",J0
            ! write(*,*) "h",h
            ! write(*,*) "N_de_vec",N_de_vec
            ! write(*,*) "Energy_density",Energy_density
            ! write(*,*) "q_v_el",q_v_el
            ! stop
            ! q_v_el = q_v_el + ( matmul(Q_el, d_el) - 2*Energy_density * N_de_vec * J0 * peso * h )
            q_v_el_aux = q_v_el_aux - 2*Energy_density * N_de_vec * J0 * peso * h 


            ! print *, "TEST 1"
            ! print*, "Shape of PSI_el: ", shape(PSI_el)
            ! print*, "Shape of FI_el: ", shape(FI_el)
            ! print*, "Shape of Q_el: ", shape(Q_el)
            ! STOP


        end do ! HAMMER POINTS

        q_v_el = matmul(Q_el, d_el) + q_v_el_aux
        do ino=1,npe
            indexes_sparse(ino) = int(ELEMS(iel, ino))
            q_v(int(ELEMS(iel, ino))) = q_v(int(ELEMS(iel, ino))) + q_v_el(ino)
        end do


        ! print *, ELEMS(iel, :)

        ! print*,"FI_el(i,:)"
        ! do i=1,npe
        !     print*, FI_el(i,:)
        ! end do
        ! if (iel == 1) then !1137
        !     write(*,*) "Q_el"
        !     do i=1,npe
        !         write(*,*) Q_el(i,:)
        !     end do

        !     write(*,*) "q_v_el",q_v_el

        !     stop
        ! end if



        ! print *, "B_de(i,:)"
        ! do i=1,2
        !     print *, B_de(i,:)
        ! end do
        ! stop

        ! write(*,*) "test1.1.1"
        call add_matrix(Q_sp, Q_el, indexes_sparse, npe)

    end do ! ELEMENTS
    
    ! SOLVE SYSTEM
    ! write(*,*) "test1.1.2"
    call assemble_sparse_matrix(Q_sp,timeit=.false.)

    ! call boundaryCondition(VINC, Q_sp, q_v)

    
    ! call solve_system_of_equation(Q_sp, g_residue, delta_Y ) 

end subroutine phaseField





subroutine PSOR_Solver(nnodes ,M, q, x, IR , JC , indsort)
implicit none

type(sparse_matrix) :: M
real(8), intent(in) :: q(:)
real(8), intent(out) :: x(:)
real(8), allocatable :: x_km1(:) 
real(8), allocatable :: f(:)
integer,intent(in) :: nnodes
integer :: iter, n
integer :: k, Jcol, Irow , i
real(8) :: D_m1, Mx, Ldx
real(8) :: err_x, err_ap, err_a0
logical :: check1, check2, check3

real(8), allocatable:: PA(:)
integer(4), intent(in) :: IR(:), JC(:)
integer, allocatable,intent(IN) :: indsort(:)



n = size(q)

! call sparse_COO_to_CCS(M, n, PA, IR, JC, indsort)
! Update the values for PA_ccs using the latest spmat%val
! write(*,*) 'update_PA_ccs'
call update_PA_ccs(nnodes,M, indsort, PA)
! write(*,*) 'update_PA_ccs done'
! Proceed with computations using PA_ccs, IR_ccs, JC_ccs
! OPEN(UNIT=1, FILE='PA.txt', STATUS='UNKNOWN')
! OPEN(UNIT=2, FILE='JC.txt', STATUS='UNKNOWN')
! OPEN(UNIT=3, FILE='IR.txt', STATUS='UNKNOWN')

! do i = 1, size(PA)
!     write(1, '(F12.6)') PA(i)
! end do

! do i = 1, size(JC)
!     write(2, '(I10)') JC(i)
! end do

! do i = 1, size(IR)
!     write(3, '(I10)') IR(i)
! end do
! stop

! write(*,*) 'sum of PA = ', sum(PA)
! STOP

allocate(x_km1(n))
allocate(f(n))

if(size(q) .ne. (size(JC)-1)) then
    write(*,*) ' projected SOR : ERROR '
    write(*,*) ' matrix and input vector dimensions do not match '
    write(*,*) ' size(q) = ', size(q), ' size(JC) = ', size(JC)
    read(*,*)
    stop
end if

if(size(x) .ne. (size(JC)-1)) then
    write(*,*) ' projected SOR : ERROR '
    write(*,*) ' matrix and output vector dimensions do not match '
    read(*,*)
    stop
end if

! write(fout,'(A11,A8,A11,A24,A15,A15,A15)') &
!   '|','*',':','sub-iter','err_d','err_a0','err_ap'

!  PROJECTED SUCCESSIVE OVER-RELAXATION

! write(*,*) 'IR', IR
! STOP

x = 0.d0               ! initialize solution 

! ************************************************************************ !
!  ITERATION : start

do iter = 1, 100000    

    x_km1 = x              ! update solution @ previous iteration x_km1

    ! -------------------------------------------------------------------- !
    !   SOLUTION COMPONENT (i.e. M - columns) LOOP 

    do jcol = 1, n                

        Mx    = 0.d0           ! initialize scalar product M*x_km1     
        Ldx   = 0.d0           ! initialize scalar product L*(x_k-x_km1)
        D_m1  = 0.d0           ! initialize scalar    M(j,j)^-1

    ! - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    !   M - ROWS LOOP
    ! print*, 'JC(1)', JC(1)
    ! print*, 'JC(last)', JC(size(JC))
    ! stop
        do k = JC(jcol),(JC(jcol + 1) - 1)  
        ! write(*,*) 'k', k
            irow      = IR(k)
            Mx        = Mx +(PA(k))*(x_km1(irow))  
    
            ! diagonal term
            if(irow .eq. jcol) then                    
                D_m1 =(PA(k))**(-1d0)     
            end if
    
            ! strict upper triangular term
            if(irow .lt. jcol) then   
                Ldx  = Ldx +(PA(k))*(x(irow) - x_km1(irow))
            end if

        end do   
    ! - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - !

        x(jcol) = max(0d0, x_km1(jcol) -(D_m1)*(Mx + q(jcol) + Ldx))
        
    end do 
    ! -------------------------------------------------------------------- !


    !  CONVERGENCE CRITERIA

    err_x   = 0.d0
    err_a0  = 0.d0
    err_ap  = 0.d0

    ! MAX NORM OF THE SOLUTION VARIATION             
    err_x   = maxval(abs(x-x_km1))
    
    ! MAX NORM OF THE ACTIVATION VIOLATION
    f= 0.0d0
    call dot_matrix_vector(M, x, f)
    f       = f + q 

    do jcol = 1, n
        if(x(jcol) .gt. 0d0) then          
            err_ap = max(err_ap, abs(f(jcol))) ! x > 0 --> act = 0 
        else
            err_a0 = min(err_a0, f(jcol))      ! x = 0 --> act > 0 
        endif 
    end do

    !check1 = err_x  .lt. 1d-4
    check1 = err_x  .lt. 1d-6
    check2 = err_a0 .lt. 1d-9
    check3 = err_ap .lt. 1d-9

    if(check1 .and. check2 .and. check3) exit

end do
    
    ! write(fout,'(A11,  A8, A11, I24, es15.2, es15.2, es15.2)') &
	! 		      '|', '*', ':', iter, err_x, err_a0, err_ap

	!  ITERATION : end
	! ************************************************************************ !

end subroutine PSOR_Solver




subroutine update_PA_ccs(nnodes,spmat, indsort, PA_ccs)
    implicit none

    type(sparse_matrix)::spmat
    integer,intent(in)::nnodes

    real(8), allocatable:: val(:),val_aux(:), lower_triangular_pa(:), val_full(:)
    integer(4), allocatable:: row(:),row_aux(:), col(:),col_aux(:), lower_triangular_ir(:), lower_triangular_jc(:), row_full(:), col_full(:)
    integer:: i,count,n1,n2,n_full,diagonal
    integer, allocatable :: indsort(:)
    integer, allocatable :: row_full_sorted(:), col_full_sorted(:)
    real(8), allocatable :: val_full_sorted(:)

    integer :: ind, ir, ir_prev, jc, jc_prev
    real(8) :: pa, pa_prev
    logical :: reordering

    real(8), allocatable,intent(out) :: PA_ccs(:)
    integer :: nterms

    ! write(*,*) 'X 1'

    val_aux=spmat%val
    col_aux=spmat%col
    row_aux=spmat%row
    nterms=spmat%nterms

    allocate(val(nterms),row(nterms),col(nterms))
    do i=1,(nterms)
        val(i)=val_aux(i)
        row(i)=row_aux(i)
        col(i)=col_aux(i)
    end do

    ! write(*,*) 'X 2'

    allocate(lower_triangular_ir(size(val)-nnodes),lower_triangular_jc(size(val)-nnodes),lower_triangular_pa(size(val)-nnodes))
    lower_triangular_ir=0
    lower_triangular_jc=0
    lower_triangular_pa=0.0d0

    ! Loop over the arrays
    count = 0
    diagonal=0
    do i = 1, size(val)
        if (row(i) /= col(i)) then  ! Not a diagonal element
            count = count + 1
            lower_triangular_ir(count) = col(i)
            lower_triangular_jc(count) = row(i)
            lower_triangular_pa(count) = val(i)
        else if (row(i) == col(i)) then  ! Diagonal element
            diagonal=diagonal+1
        end if
    end do

!   write(*,*) 'X 3'


    ! Determine the sizes
    n1 = size(val)                      ! Assuming 'val', 'row', 'col' are of the same size
    n2 = size(lower_triangular_pa)      ! Size of the lower triangular arrays
    n_full = n1 + n2

    ! Allocate concatenated arrays
    allocate(row_full(n_full))
    allocate(col_full(n_full))
    allocate(val_full(n_full))

    ! Concatenate 'row' and 'lower_triangular_ir' into 'row_full'
    row_full(1:n1) = row
    
    ! write(*,*) 'size indsort', size(indsort)
    ! stop

    row_full(n1+1:n_full) = lower_triangular_ir

! write(*,*) 'X 4'
    ! stop


    ! Concatenate 'col' and 'lower_triangular_jc' into 'col_full'
    col_full(1:n1) = col
    col_full(n1+1:n_full) = lower_triangular_jc

    ! Concatenate 'val' and 'lower_triangular_pa' into 'val_full'
    val_full(1:n1) = val
    val_full(n1+1:n_full) = lower_triangular_pa


    deallocate(row,col,val,lower_triangular_ir,lower_triangular_jc,lower_triangular_pa)



! write(*,*) 'X 5'

    allocate(row_full_sorted(n_full))
    allocate(col_full_sorted(n_full))
    allocate(val_full_sorted(n_full))
    
    row_full_sorted = row_full(indsort)
    col_full_sorted = col_full(indsort)
    val_full_sorted = val_full(indsort)

    
    row_full = row_full_sorted
    col_full = col_full_sorted
    val_full = val_full_sorted


    
    deallocate(row_full_sorted)
    deallocate(col_full_sorted)
    deallocate(val_full_sorted)

! write(*,*) 'X 5'
    
    ! Initialize the reordering flag

reordering = .true.
do while (reordering)
    reordering = .false.
    do ind = 2, size(row_full)
        ir = row_full(ind)
        ir_prev = row_full(ind - 1)
        jc = col_full(ind)
        jc_prev = col_full(ind - 1)
        pa = val_full(ind)
        pa_prev = val_full(ind - 1)
        if ((ir == ir_prev) .and. (jc_prev > jc)) then
            reordering = .true.
            ! Swap col_full values
            col_full(ind - 1) = jc
            col_full(ind) = jc_prev
            ! Swap val_full values
            val_full(ind - 1) = pa
            val_full(ind) = pa_prev
        end if
    end do
end do


! write(*,*) 'X 6'

allocate(PA_ccs(size(val_full)))
PA_ccs=val_full



end subroutine update_PA_ccs

subroutine initialize_sparse_COO_to_CCS(Q_sp,n, IR_ccs, JC_ccs, indsort, d_n , ELEMS , NODES , NODES_1 , l_0 , G_c , gaprox , npe , nnodes , nelems ,  ep)
    implicit none
    integer, intent(in) :: npe , nnodes , nelems , gaprox
    character(len=100) :: ep
    real, dimension(:,:), intent(in) :: ELEMS
    real(8), dimension(:,:), intent(in) :: NODES, NODES_1


    !=====  LOCAL VARIABLES  ======!
    

    !==== PHASE FIELD VARIABLES =====!
    type(sparse_matrix) :: Q_sp
    real(8), dimension(:), intent(in) :: d_n
    real(8), intent(in) :: l_0, G_c
    real(8), dimension(nnodes) :: q_v
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: IR_ccs(:), JC_ccs(:), indsort(:)

    ! Local variables
    integer :: nterms
    integer, allocatable :: row_full(:), col_full(:)
    integer :: n_full, n1, n2
    integer ::  count, ind_cum, count_cum
    integer :: diagonal
    integer :: ind
    logical :: reordering
    integer :: ir, ir_prev, jc, jc_prev , i


    real(8), allocatable:: val(:),val_aux(:), lower_triangular_pa(:), val_full(:)
    integer(4), allocatable:: row(:),row_aux(:), col(:),col_aux(:), lower_triangular_ir(:), lower_triangular_jc(:)
    
    integer, allocatable :: row_full_sorted(:), col_full_sorted(:)
    real(8), allocatable :: val_full_sorted(:)


    ! Note: PA related variables are omitted here
    ! write(*,*) "test 1.1"
    call phaseField(Q_sp , q_v , d_n , ELEMS , NODES , NODES_1 , l_0 , G_c , gaprox , npe , nnodes , nelems ,  ep)
    ! write(*,*) "test 1.2"
    ! Extract matrix data
    val_aux = Q_sp%val
    col_aux = Q_sp%col
    row_aux = Q_sp%row
    nterms = Q_sp%nterms

    allocate(val(nterms),row(nterms),col(nterms))
    do i=1,(nterms)
        val(i)=val_aux(i)
        row(i)=row_aux(i)
        col(i)=col_aux(i)
    end do



    allocate(lower_triangular_ir(size(val)-n),lower_triangular_jc(size(val)-n),lower_triangular_pa(size(val)-n))
    lower_triangular_ir=0
    lower_triangular_jc=0
    lower_triangular_pa=0.0d0

    ! Loop over the arrays
    count = 0
    diagonal=0
    do i = 1, size(val)
        if (row(i) /= col(i)) then  ! Not a diagonal element
            count = count + 1
            lower_triangular_ir(count) = col(i)
            lower_triangular_jc(count) = row(i)
            lower_triangular_pa(count) = val(i)
        else if (row(i) == col(i)) then  ! Diagonal element
            diagonal=diagonal+1
        end if
    end do



    ! Determine the sizes
    n1 = size(val)                      ! Assuming 'val', 'row', 'col' are of the same size
    n2 = size(lower_triangular_pa)      ! Size of the lower triangular arrays
    n_full = n1 + n2

    ! Allocate concatenated arrays
    allocate(row_full(n_full))
    allocate(col_full(n_full))
    allocate(val_full(n_full))

    ! Concatenate 'row' and 'lower_triangular_ir' into 'row_full'
    row_full(1:n1) = row



    row_full(n1+1:n_full) = lower_triangular_ir


    ! Concatenate 'col' and 'lower_triangular_jc' into 'col_full'
    col_full(1:n1) = col
    col_full(n1+1:n_full) = lower_triangular_jc

    ! Concatenate 'val' and 'lower_triangular_pa' into 'val_full'
    val_full(1:n1) = val
    val_full(n1+1:n_full) = lower_triangular_pa


    deallocate(row,col,val,lower_triangular_ir,lower_triangular_jc,lower_triangular_pa)

    allocate(indsort(n_full))
    indsort = [(i, i = 1, n_full)]  ! indsort = [1, 2, ..., n_full]
    
    call argsort(n_full, row_full, indsort)



    allocate(row_full_sorted(n_full))
    allocate(col_full_sorted(n_full))
    allocate(val_full_sorted(n_full))
    
    row_full_sorted = row_full(indsort)
    col_full_sorted = col_full(indsort)
    val_full_sorted = val_full(indsort)



    
    row_full = row_full_sorted
    col_full = col_full_sorted
    val_full = val_full_sorted



    
    deallocate(row_full_sorted)
    deallocate(col_full_sorted)


    
    ! Initialize the reordering flag

    reordering = .true.
    do while (reordering)
        reordering = .false.
        do ind = 2, size(row_full)
            ir = row_full(ind)
            ir_prev = row_full(ind - 1)
            jc = col_full(ind)
            jc_prev = col_full(ind - 1)

            if ((ir == ir_prev) .and. (jc_prev > jc)) then
                reordering = .true.
                ! Swap col_full values
                col_full(ind - 1) = jc
                col_full(ind) = jc_prev

            end if
        end do
    end do




    allocate(IR_ccs(size(val_full)),JC_ccs(n+1))
    IR_ccs=col_full

    count_cum=1
    ind_cum=1
    JC_ccs(1)=1   

    ! stop

    do ind=2,size(row_full)
        ir=row_full(ind)
        ir_prev=row_full(ind-1)
        if (ir==ir_prev) then
            count_cum=count_cum+1
        else
            
            ind_cum=ind_cum+1
            JC_ccs(ind_cum)=JC_ccs(ind_cum-1)+count_cum

            count_cum=1

        end if
    end do


    ind_cum=ind_cum+1

    JC_ccs(ind_cum)=JC_ccs(ind_cum-1)+count_cum

    call clear_data(Q_sp)

end subroutine initialize_sparse_COO_to_CCS






subroutine argsort(n, array, index)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: array(n)
    integer, intent(inout) :: index(n)
    integer :: i, j, temp_idx
    integer :: temp_value

    ! write(*,*)"___________________"
    ! write(*,*)"n",n
    ! stop
    do i = 2, n
        temp_idx = index(i)
        temp_value = array(temp_idx)
        j = i - 1
        ! write(*,*)"___________________"
        ! write(*,*)"j",j
        ! write(*,*)"index(j)",index(j)
        ! write(*,*)"size of index",size(index)
        ! write(*,*)"size of array",size(array)
        do while (j >= 1 .and. ((array(index(j)) > temp_value) .or. &
            (array(index(j)) == temp_value .and. index(j) > temp_idx)))
            index(j + 1) = index(j)
            j = j - 1
            if (j==0) then
                exit
            end if
        end do
        index(j + 1) = temp_idx
    end do
end subroutine argsort





end module phaseField_module
