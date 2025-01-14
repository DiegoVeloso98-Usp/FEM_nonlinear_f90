
module test2
    
    use sparse
    implicit none

    
contains




subroutine COO_upper_to_CCS_lower(val, row, col, nnz, n, PA, IR, JC)
    implicit none
    integer, intent(in) :: nnz
    integer, intent(in) :: n
    real(8), intent(in) :: val(nnz)
    integer, intent(in) :: row(nnz), col(nnz)
    real(8), allocatable, intent(out) :: PA(:)
    integer, allocatable, intent(out) :: IR(:), JC(:)
    integer :: i, count, nz_lower
    real(8), allocatable :: val_lower(:)
    integer, allocatable :: row_lower(:), col_lower(:)
    integer, allocatable :: perm(:)
    integer :: idx, col_idx
    integer, allocatable :: col_sorted(:)

    ! Step 1: Generate Lower Triangular COO Data
    allocate(val_lower(nnz))
    allocate(row_lower(nnz))
    allocate(col_lower(nnz))
    count = 0

    do i = 1, nnz
        if (row(i) == col(i)) then
            ! Diagonal element
            count = count + 1
            val_lower(count) = val(i)
            row_lower(count) = row(i)
            col_lower(count) = col(i)
        else if (row(i) < col(i)) then
            ! Off-diagonal upper triangular element
            count = count + 1
            val_lower(count) = val(i)
            row_lower(count) = col(i)   ! Swap indices
            col_lower(count) = row(i)
        end if
        ! Elements where row > col are ignored
    end do

    nz_lower = count

    ! Step 2: Sort the Lower Triangular COO Data by Column and Row
    allocate(perm(nz_lower))
    do i = 1, nz_lower
        perm(i) = i
    end do

    call sort_coo_by_col_row(nz_lower, col_lower, row_lower, val_lower, perm)

    ! Apply the sorted order
    allocate(PA(nz_lower))
    allocate(IR(nz_lower))
    allocate(col_sorted(nz_lower))

    do i = 1, nz_lower
        PA(i) = val_lower(perm(i))
        IR(i) = row(i)
        col_sorted(i) = col_lower(perm(i))
    end do

    ! Step 3: Construct JC
    allocate(JC(n + 1))
    JC = 0
    JC(1) = 1  ! Start index for the first column

    do i = 1, nz_lower
        col_idx = col_sorted(i)
        if (JC(col_idx) == 0) then
            JC(col_idx) = i
        end if
    end do

    ! Fill in JC for columns with no entries
    do i = 2, n + 1
        if (JC(i) == 0) then
            JC(i) = JC(i - 1)
        end if
    end do
    JC(n + 1) = nz_lower + 1  ! End pointer

    ! Deallocate temporary arrays
    deallocate(val_lower, row_lower, col_lower, perm, col_sorted)
end subroutine COO_upper_to_CCS_lower


subroutine sort_coo_by_col_row(nnz, col_coo, row_coo, val_coo, perm)
    implicit none
    integer, intent(in) :: nnz
    integer, intent(in) :: col_coo(:), row_coo(:)
    real(8), intent(in) :: val_coo(:)
    integer, intent(inout) :: perm(:)
    integer :: i, j, key_perm
    integer :: key_col, key_row

    do i = 2, nnz
        key_perm = perm(i)
        key_col = col_coo(key_perm)
        key_row = row_coo(key_perm)
        j = i - 1

        do while (j >= 1 .and. (col_coo(perm(j)) > key_col  .or. &
                                (col_coo(perm(j)) == key_col .and. row_coo(perm(j)) > key_row)))
            perm(j + 1) = perm(j)
            j = j - 1
        end do
        perm(j + 1) = key_perm
    end do
end subroutine sort_coo_by_col_row


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
    
end module test2
