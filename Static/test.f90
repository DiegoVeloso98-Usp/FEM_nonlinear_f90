

program test_sparse
    use sparse
    ! use test2
    implicit none

    type(sparse_matrix)::spmat

    real(8), allocatable:: val(:), lower_triangular_pa(:), val_full(:)
    integer(4), allocatable:: row(:), col(:), lower_triangular_ir(:), lower_triangular_jc(:), row_full(:), col_full(:)
    integer:: n,i,count,n1,n2,n_full
    integer, allocatable :: indsort(:)
    integer, allocatable :: row_full_sorted(:), col_full_sorted(:)
    real(8), allocatable :: val_full_sorted(:)

    integer :: ind, ir, ir_prev, jc, jc_prev
    real(8) :: pa, pa_prev
    logical :: reordering

    integer(4), allocatable :: JC_ccs(:), IR_ccs(:)
    real(8), allocatable :: PA_ccs(:)
    integer :: count_cum, ind_cum

    call prepare_to_use(spmat,6,13)
    

    ! Here we only want to test the method to transform the sparseSET matrix
    ! to be sent to superLU

    ! First row
    spmat%row(1) = 1; spmat%col(1) = 1; spmat%val(1) = 11.0d0
    spmat%row(2) = 1; spmat%col(2) = 4; spmat%val(2) = 14.0d0
    spmat%row(3) = 1; spmat%col(3) = 5; spmat%val(3) = 15.0d0
    spmat%row(4) = 1; spmat%col(4) = 6; spmat%val(4) = 16.0d0    
    ! Second row
    spmat%row(5) = 2; spmat%col(5) = 2; spmat%val(5) = 22.0d0
    spmat%row(6) = 2; spmat%col(6) = 4; spmat%val(6) = 24.0d0
    spmat%row(7) = 2; spmat%col(7) = 5; spmat%val(7) = 25.0d0
    ! Third row
    spmat%row(8) = 3; spmat%col(8) = 3; spmat%val(8) = 33.0d0
    spmat%row(9) = 3; spmat%col(9) = 5; spmat%val(9) = 35.0d0
    ! Fourth row
    spmat%row(10) = 4; spmat%col(10) = 4; spmat%val(10) = 44.0d0
    spmat%row(11) = 4; spmat%col(11) = 6; spmat%val(11) = 46.0d0    
    ! Fifth row
    spmat%row(12) = 5; spmat%col(12) = 5; spmat%val(12) = 55.0d0
    ! Sixth row
    spmat%row(13) = 6; spmat%col(13) = 6; spmat%val(13) = 66.0d0  

    spmat%cum(1) = 1
    spmat%cum(2) = 5    
    spmat%cum(3) = 8
    spmat%cum(4) = 10
    spmat%cum(5) = 12
    spmat%cum(6) = 13
    spmat%nterms = 13


    val=spmat%val
    col=spmat%col
    row=spmat%row
    ! n=spmat%nrows (NÃO USAR. RETORNA O TAMANHO ALOCADO PARA A MATRIZ, QUE PODE SER MAIOR QUE SUA DIMENSÃO)
    ! n=size(spmat%cum) (TAMBÉM RETORNA O TAMANHO ALOCADO PARA A MATRIZ)
    n=6
    ! write(*,*)"___________________"
    ! write(*,*)"val",val
    ! write(*,*)"row",row
    ! write(*,*)"col",col
    ! write(*,*)"n",n

    allocate(lower_triangular_ir(size(val)-n),lower_triangular_jc(size(val)-n),lower_triangular_pa(size(val)-n))
    lower_triangular_ir=0
    lower_triangular_jc=0
    lower_triangular_pa=0.0d0

    ! Loop over the arrays
    count = 0
    do i = 1, size(val)
        if (row(i) /= col(i)) then  ! Not a diagonal element
            count = count + 1
            lower_triangular_ir(count) = col(i)
            lower_triangular_jc(count) = row(i)
            lower_triangular_pa(count) = val(i)
        end if
    end do
    ! write(*,*)"___________________"
    ! write(*,*)"lower_triangular_ir",lower_triangular_ir
    ! write(*,*)"lower_triangular_jc",lower_triangular_jc
    ! write(*,*)"lower_triangular_pa",lower_triangular_pa



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

    ! write(*,*)"___________________"
    ! write(*,*)"row_full",row_full
    ! write(*,*)"col_full",col_full
    ! write(*,*)"val_full",val_full
    ! write(*,*)"size(row)",size(row_full)
    ! write(*,*)"size(col)",size(col_full)
    ! write(*,*)"size(val)",size(val_full)

    
    allocate(indsort(n_full))
    indsort = [(i, i = 1, n_full)]  ! indsort = [1, 2, ..., n_full]
    
    ! Step 2: Sort the indices based on 'row_full'
    call argsort(n_full, row_full, indsort)

    ! write(*,*)"___________________"
    ! write(*,*)"size(indsort)",size(indsort)
    ! write(*,*)"indsort",indsort
    ! stop
    allocate(row_full_sorted(n_full))
    allocate(col_full_sorted(n_full))
    allocate(val_full_sorted(n_full))
    
    row_full_sorted = row_full(indsort)
    col_full_sorted = col_full(indsort)
    val_full_sorted = val_full(indsort)
    
    row_full = row_full_sorted
    col_full = col_full_sorted
    val_full = val_full_sorted
    ! write(*,*)"___________________"
    ! write(*,*)"row_full sorted ",row_full
    ! write(*,*)"col_full sorted",col_full
    ! write(*,*)"val_full sorted",val_full
    ! stop

    
    deallocate(row_full_sorted)
    deallocate(col_full_sorted)
    deallocate(val_full_sorted)

    ! write(*,*)"___________________"
    ! write(*,*)"row_full",row_full
    ! write(*,*)"col_full",col_full
    ! write(*,*)"val_full",val_full
    
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

write(*,*)"___________________"
write(*,*)"row_full",row_full
write(*,*)"col_full",col_full
write(*,*)"val_full",val_full
allocate(PA_ccs(size(val_full)),IR_ccs(size(val_full)),JC_ccs(n+1))
PA_ccs=val_full
IR_ccs=col_full

count_cum=1
ind_cum=1
JC_ccs(1)=1   
do ind=2,size(row_full)
    ir=row_full(ind)
    ir_prev=row_full(ind-1)
    
    if (ir==ir_prev) then
        count_cum=count_cum+1
    else
        ind_cum=ind_cum+1
        JC_ccs(ind_cum)=JC_ccs(ind_cum-1)+count_cum
        IF (ind_cum==size(JC_ccs)) THEN
            write(*,*)"______________________________________"
            write(*,*)"Error, already last position of JC| ir=",ir," | ir_prev=",ir_prev
            write(*,*)"______________________________________"
        end if
        write(*,*)" ir=",ir," | ir_prev=",ir_prev
        count_cum=1
        write(*,*)"LOOPING row_full array : || Column index",ind_cum," | Elements in column ",count_cum
    end if
end do

write(*,*)"Done looping row_full array"

ind_cum=ind_cum+1

write(*,*)"Writing last position of JC || Column index",ind_cum," | Elements in column ",count_cum

JC_ccs(ind_cum)=JC_ccs(ind_cum-1)+count_cum

! write(*,*)"___________________"
! write(*,*)"PA_ccs",PA_ccs
! write(*,*)"IR_ccs",IR_ccs
! write(*,*)"JC_ccs",JC_ccs



end program test_sparse

subroutine argsort(n, array, index)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: array(n)
    integer, intent(inout) :: index(n)
    integer :: i, j, temp_idx
    integer :: temp_value

    do i = 2, n
        temp_idx = index(i)
        temp_value = array(temp_idx)
        j = i - 1
        do while (j >= 1 .and. ((array(index(j)) > temp_value) .or. &
            (array(index(j)) == temp_value .and. index(j) > temp_idx)))
            index(j + 1) = index(j)
            j = j - 1
        end do
        index(j + 1) = temp_idx
    end do
end subroutine argsort