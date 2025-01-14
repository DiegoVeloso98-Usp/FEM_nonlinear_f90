!-----------------------------------------------------------------------------!
!   sparseSET - Symmetric sparse matrix manipulation routines                 !
!   Authors: Dorival Piedade Neto and Rodrigo Ribeiro Paccola                 !
!   University of Sao Paulo - Sao Carlos Engineering School - Nov 2020        !
!   (This version is based in the version released in May 2012)               !
!   (Changed sparse solver - it was MA67, now using SuperLU)                  !
!   This will be considered the version 1.0 (previous ones were 0.x version)  !
!   This file is released according to the BSD license (the same used by the  !
!   files from superLU, which are also found in this repository)              !
!-----------------------------------------------------------------------------!

module sparse
    implicit none

    type sparse_matrix
        !private
        integer(4)::nterms = 0
        integer(4)::nrows
        integer(4),dimension(:),allocatable::cum
        integer(4),dimension(:),allocatable::row, col
        real(8),dimension(:),allocatable::val
    end type sparse_matrix    

    integer(4),dimension(:),allocatable::sp_row, sp_col
    real(8),dimension(:),allocatable::sp_val, sp_vector
    integer(4)::sp_nterms

    contains
    
    subroutine prepare_to_use(sp_matrix, nvar, nterm)
        ! Prepares the sparse matrix to be build by matrices contributions
        ! nvar - number of variables (unknowns) of the system of equation
        ! nterm -number of terms to build the sparse matrix * See A.1.
        implicit none
        type(sparse_matrix),intent(inout)::sp_matrix
        integer(4),intent(in)::nvar,nterm
        ! Deallocating any possibly allocated sparse matrix intern arrays
        if (allocated(sp_matrix%cum)) deallocate(sp_matrix%cum)
        if (allocated(sp_matrix%row)) deallocate(sp_matrix%row)
        if (allocated(sp_matrix%col)) deallocate(sp_matrix%col)
        if (allocated(sp_matrix%val)) deallocate(sp_matrix%val)
        ! Allocating sparse matrix inter arrays
        allocate(sp_matrix%cum(nvar))
        allocate(sp_matrix%row(nterm), sp_matrix%col(nterm), &
                 sp_matrix%val(nterm))
        sp_matrix%nrows=nvar
        sp_matrix%cum = 0
    end subroutine prepare_to_use

    ! A.1. The nterm is the number of terms to be added to the sparse matrix to
    !      build it. It does not need to be exact, but if the actual number of 
    !      terms added to the matrix is greater than nterm, an error will occur
    !      since the space reserved is smaller than needed . So if this number
    !      is not precisilly know, its better to over estimated it.
    !      Since simetry is considered to use less memory and speed up the 
    !      sparse maytrix assembly, in order to compute nterm, one must use the
    !      following rule to compute it:
    !      n = (ndof^2+ndof) / 2, in which n is the local matrix number of
    !      terms, and ndof is the number of degrees of freedom of the local 
    !      matrix (equals to its number of rows/columns). 

    subroutine clear_data(sp_matrix)
        ! 'Clears' all stored data so the sparse matrix can be assemble again
        ! (In fact just sets the nterms and cum to zero, so new inputs will
        ! over write old data!)
        implicit none
        type(sparse_matrix),intent(inout)::sp_matrix
        sp_matrix%nterms = 0
        sp_matrix%cum = 0
    end subroutine clear_data

    subroutine deallocate_sparse_matrix(sp_matrix)
        ! Deallocates the sparse matrix
        implicit none
        type(sparse_matrix),intent(inout)::sp_matrix
        sp_matrix%nterms = 0
        sp_matrix%nrows = 0
        if (allocated(sp_matrix%cum)) deallocate(sp_matrix%cum)
        if (allocated(sp_matrix%row)) deallocate(sp_matrix%row)
        if (allocated(sp_matrix%col)) deallocate(sp_matrix%col)
        if (allocated(sp_matrix%val)) deallocate(sp_matrix%val)
    end subroutine deallocate_sparse_matrix

    subroutine add_matrix(sp_matrix, matrix, indexes, n)
        ! Inserts the matrix terms in the sparse_matrix, using the indexes info.
        ! The sparse matrix must be already prepared to be used; matrix is an
        ! (n x n) real(8) array and indexes is an integer (n) array.
        implicit none
        type(sparse_matrix),intent(inout)::sp_matrix
        integer(4), intent(in)::n
        real(8), dimension(n,n), intent(in)::matrix
        integer(4), dimension(n), intent(in)::indexes
        integer(4)::i, j, num, ind, row_ind, col_ind, aux
        num = (n ** 2 + n) / 2
        ind = sp_matrix%nterms
        do i = 1, n
            do j = 1, i
                ind = ind + 1
                row_ind = indexes(j)
                col_ind = indexes(i)
                if (col_ind .lt. row_ind) then
                    aux = row_ind
                    row_ind = col_ind
                    col_ind = aux
                end if
                sp_matrix%row(ind) = row_ind
                sp_matrix%col(ind) = col_ind
                sp_matrix%val(ind) = matrix(j,i)
            enddo
        enddo
        sp_matrix%nterms = ind
    end subroutine add_matrix

    subroutine add_full_matrix(sp_matrix, matrix, indexes, n)
        ! Inserts the matrix terms in the sparse_matrix, using the indexes info.
        ! The sparse matrix must be already prepared to be used; matrix is an
        ! (n x n) real(8) array and indexes is an integer (n) array.
        implicit none
        type(sparse_matrix),intent(inout)::sp_matrix
        integer(4), intent(in)::n
        real(8), dimension(:,:), intent(in)::matrix
        integer(4), dimension(:), intent(in)::indexes
        integer(4)::i, j, num, ind, row_ind, col_ind
        num = (n ** 2 + n) / 2
        ind = sp_matrix%nterms
        do i = 1, n
            do j = 1, n
                row_ind = indexes(i)
                col_ind = indexes(j)
                if (col_ind .ge. row_ind) then
                    ind = ind + 1
                    sp_matrix%row(ind) = row_ind
                    sp_matrix%col(ind) = col_ind
                    sp_matrix%val(ind) = matrix(i,j)
                end if
            enddo
        enddo
        sp_matrix%nterms = ind
    end subroutine add_full_matrix

    subroutine add_sliced_matrix(sp_matrix, row, col, val, n)
        ! Inserts the sliced matrix terms in the sparse_matrix.
        ! The sparse matrix must be already prepared to be used; matrix is an
        ! (n) real(8) array and indexes is an integer (n) array.
        implicit none
        type(sparse_matrix),intent(inout)::sp_matrix
        integer(4), intent(in)::n
        real(8), dimension(n), intent(in)::val
        integer(4), dimension(n), intent(in)::row, col
        integer(4)::i, ind, row_ind, col_ind
        ind = sp_matrix%nterms
        do i = 1, n
            row_ind = row(i)
            col_ind = col(i)
            if (col_ind.ge.row_ind) then
                ind = ind + 1
                sp_matrix%row(ind) = row_ind
                sp_matrix%col(ind) = col_ind
                sp_matrix%val(ind) = val(i)
            end if
        enddo
        sp_matrix%nterms = ind
    end subroutine add_sliced_matrix

    subroutine shell_sort(n,a,b,c)
        ! Sorts arrays a and b (integer*4) and c (real*8) into ascending order
        ! of a by Shell's method. Adapted from the book "Numerical Recipies in
        ! in Fortran 77"
        implicit none
        integer(4),intent(in)::n
        integer(4),dimension(n),intent(inout)::a,b
        real(8),dimension(n),intent(inout)::c
        integer(4)::i,j,inc,va,vb
        real(8)::vc
        inc=1   ! Determine the starting increment.
1       inc=3*inc+1
        if(inc.le.n)goto 1
2       continue    ! Loop over the partial sorts.
            inc=inc/3
            do i=inc+1,n ! Outer loop of straight insertion.
                va=a(i)
                vb=b(i)
                vc=c(i)
                j=i
3               if(a(j-inc).gt.va)then   ! Inner loop of straight insertion.
                    a(j)=a(j-inc)
                    b(j)=b(j-inc)
                    c(j)=c(j-inc)
                    j=j-inc
                    if(j.le.inc)goto 4
                    goto 3
                endif
4               a(j)=va
                b(j)=vb
                c(j)=vc
            enddo
        if(inc.gt.1)goto 2
        return
    end subroutine shell_sort

    subroutine row_data(sp_matrix, row_number, index_begin, index_end, nterms)
        !    Returns the index of the begin and of the end of the row..
        implicit none
        type(sparse_matrix),intent(in)::sp_matrix
        integer(4),intent(in)::row_number
        integer(4),intent(out)::index_begin,index_end,nterms
        index_begin = sp_matrix%cum(row_number)
        if (row_number.lt.sp_matrix%nrows) then
            index_end = sp_matrix%cum(row_number+1) - 1
            nterms = sp_matrix%cum(row_number+1) - sp_matrix%cum(row_number)
        else
            index_end = sp_matrix%nterms
            nterms = sp_matrix%nterms - sp_matrix%cum(row_number) + 1
        endif
    end subroutine row_data

    subroutine assemble_sparse_matrix(sp_matrix, timeit)
        !   Assembles the sparse matrix, sorting its rows and columns. Also
        ! sum all terms in the same position and prepares it to be used by
        ! the direct sparse solver. After the matrix is assemble is possible
        ! to change its terms, impose a value in its rows and columns and
        ! multiply it with a vector, but if new terms no new terms can be
        ! included yet.
        implicit none
        type(sparse_matrix),intent(inout)::sp_matrix
        integer(4)::nterms,i,j,ind_i,ind_j,current_ind,current_col,counter
        integer(4)::accumulated, current_row_nterms, max_ind
        logical(4),optional::timeit
        logical(4)::eval_time
        real(8)::ti,tf
        if (present(timeit)) then
            eval_time = timeit
        else
            eval_time = .false.
        endif
        ! First, sort the terms by its row number
        call cpu_time(ti)
        nterms = sp_matrix%nterms
        call shell_sort(nterms,sp_matrix%row(1:nterms),sp_matrix%col(1:nterms)&
                       ,sp_matrix%val(1:nterms))
        call cpu_time(tf)
        if (eval_time) print *, 'Time to sort rows:',(tf-ti)
        ! Count number of terms per row (including repeated terms (same column)
        call cpu_time(ti)
        do i=1,nterms
            sp_matrix%cum(sp_matrix%row(i)) = &
            sp_matrix%cum(sp_matrix%row(i)) + 1
        enddo
        accumulated = 1
        do i=1,sp_matrix%nrows
              current_row_nterms = sp_matrix%cum(i) 
              sp_matrix%cum(i) = accumulated
              accumulated = accumulated + current_row_nterms
        enddo
        call cpu_time(tf)
        if (eval_time) print *, 'Time to count terms per rows:',(tf-ti)
        call cpu_time(ti)
        ! Sorting each row by its column number
        do i=1, sp_matrix%nrows
            call row_data(sp_matrix, i, ind_i, ind_j,nterms) 
            call shell_sort(nterms, sp_matrix%col(ind_i: ind_j), &
                 sp_matrix%row(ind_i: ind_j), sp_matrix%val(ind_i: ind_j))
            ! After sorting all columns, them sum all repeated terms in the row
            current_ind = ind_i
            current_col = sp_matrix%col(ind_i)
            do j=ind_i+1, ind_j
                if (sp_matrix%col(j).eq.current_col) then
                    sp_matrix%val(current_ind) = sp_matrix%val(current_ind) + &
                    sp_matrix%val(j)
                    sp_matrix%row(j)=0 !Term to be eliminated (repeated)
                else
                    current_col = sp_matrix%col(j)
                    current_ind = j
                endif
            enddo
        enddo
        call cpu_time(tf)
        if (eval_time) print *, 'Time to sort cols and sum terms:',(tf-ti)
        ! Finally, assemble the sparse matrix with no repeated terms and update
        ! the number of terms per row.
        call cpu_time(ti)
        sp_matrix%cum = 0
        counter = 0
        nterms = sp_matrix%nterms
        do i=1,nterms
            if (sp_matrix%row(i).ne.0) then
                counter = counter + 1
                sp_matrix%cum(sp_matrix%row(i)) = &
                sp_matrix%cum(sp_matrix%row(i)) + 1
                sp_matrix%row(counter) = sp_matrix%row(i)
                sp_matrix%col(counter) = sp_matrix%col(i)
                sp_matrix%val(counter) = sp_matrix%val(i)
            endif
        enddo
        ! Allocating the global arrays sp_row, sp_col and sp_val
        max_ind = maxval(sp_matrix%cum)  ! At this point, cum stores the number
        max_ind = 3 * max_ind            ! We have computed for rows only ...
        call allocate_global_arrays(max_ind)
        ! Finally, compute the accumulated value (sp_matrix%cum)
        sp_matrix%nterms = counter
        accumulated = 1
        do i=1,sp_matrix%nrows
            current_row_nterms = sp_matrix%cum(i) 
            sp_matrix%cum(i) = accumulated
            accumulated = accumulated + current_row_nterms
        enddo
        call cpu_time(tf)
        if (eval_time) print *, 'Time to rebuild shortened sparse matrix:',(tf-ti)
    end subroutine assemble_sparse_matrix

    subroutine allocate_global_arrays(nterms)
        ! Allocates the global arrays sp_row, sp_col, sp_val if nterms is
        ! greater than the current size of these arrays.
        integer(4),intent(in)::nterms
        if (allocated(sp_row)) then
            if (size(sp_row).gt.nterms) then
                return
            else
                deallocate(sp_row, sp_col, sp_val)
            endif
        endif
        allocate(sp_row(nterms), sp_col(nterms), sp_val(nterms))
    end subroutine allocate_global_arrays

    subroutine optimize_storage(sp_matrix)
        ! Optimizes storage (memory) reallocating the row, col and data arrays
        ! to the exact size needed to store the assembled sparse matrix.
        implicit none
        type(sparse_matrix),intent(inout)::sp_matrix
        integer(4),dimension(:),allocatable::row_col
        real(8),dimension(:),allocatable::val
        integer(4)::nterms
        nterms = sp_matrix%nterms
        allocate(row_col(nterms))
        ! Row
        row_col = sp_matrix%row(1:nterms)
        deallocate(sp_matrix%row)
        allocate(sp_matrix%row(nterms))
        sp_matrix%row = row_col
        ! Column
        row_col = sp_matrix%col(1:nterms)
        deallocate(sp_matrix%col)
        allocate(sp_matrix%col(nterms))
        sp_matrix%col = row_col
        ! Deallocating row_col and allocating val
        deallocate(row_col)
        allocate(val(nterms))
        ! Values
        val = sp_matrix%val(1:nterms)
        deallocate(sp_matrix%val)
        allocate(sp_matrix%val(nterms))
        sp_matrix%val = val
        deallocate(val)
    end subroutine optimize_storage

    subroutine print_sparse_matrix(sp_matrix)
        ! Prints the sparse matrix in the screen.
        implicit none
        type(sparse_matrix),intent(in)::sp_matrix
        integer(4)::i
        write(*,*) ''
        write(*,*) '     index      row      col          value'
        write(*,*) '-----------------------------------------------'
        do i=1, sp_matrix%nterms
            write(*,10) i, sp_matrix%row(i), sp_matrix%col(i), sp_matrix%val(i)
        enddo
        10 format(3i10,2x,es16.8)
    end subroutine print_sparse_matrix

    subroutine print_memory_storage_report(sp_matrix)
        ! Prints the sparse matrix memory storage report
        implicit none
        type(sparse_matrix),intent(in)::sp_matrix
        integer(4)::nterms,max_terms
        real(8)::current_mem, total_mem
        nterms = sp_matrix%nterms
        max_terms = size(sp_matrix%row)
        current_mem = 16.0d0 * real(nterms) / 1024.0d0 ** 2
        total_mem = 16.0d0 * real(max_terms) / 1024.0d0 ** 2
        write(*,*)''
        write(*,10) 'Sparse matrix containing ',sp_matrix%nrows, ' unknows.'
        write(*,20) 'Current use: ', nterms, ' terms (',current_mem,' MB)'
        write(*,20) 'Total size: ', max_terms, ' terms (',total_mem,' MB)'
        write(*,30) 'Using ',(100.0d0 * current_mem / total_mem),' % of the total allocated space'
        write(*,*)''
        10 format(a,i7,a)
        20 format(a,i10,a,f12.5,a)
        30 format(a,f6.2,a)
    end subroutine print_memory_storage_report

    subroutine dot_matrix_vector(sp_matrix, vector, solution)
        ! Performs the dot product of the sparse matrix (A) and a given vector
        ! x, resulting in a vector y (y = A x) . Result (y) is returned in the 
        ! solution vector given by the user (same size of vector).
        ! Note that since sp_matrix is symmetric, the present operation is 
        ! identical to the operation (y = x^T A).
        implicit none
        type(sparse_matrix),intent(in)::sp_matrix
        real(8),dimension(sp_matrix%nrows),intent(in)::vector
        real(8),dimension(sp_matrix%nrows),intent(out)::solution
        integer(4)::i
        solution = 0.0d0
        do i=1, sp_matrix%nterms
            solution(sp_matrix%row(i)) = solution(sp_matrix%row(i)) + &
            sp_matrix%val(i) * vector(sp_matrix%col(i))
            if (sp_matrix%row(i).ne.sp_matrix%col(i)) &
                solution(sp_matrix%col(i)) = solution(sp_matrix%col(i)) + &
                sp_matrix%val(i) * vector(sp_matrix%row(i))
        enddo
    end subroutine dot_matrix_vector

    subroutine solve_system_of_equation(sp_matrix, vector, solution)
        ! Solves the system of equation A x = b, given by the sp_matrix
        ! (A) and the vector (b), resulting in the solution (x), using the 
        ! direct solver SuperLU (sequential) version 5.2.x (April 8, 2016).
        implicit none
        type(sparse_matrix),intent(in)::sp_matrix
        real(8),dimension(sp_matrix%nrows),intent(in)::vector
        real(8),dimension(sp_matrix%nrows),intent(out)::solution
        call superlu_solver(sp_matrix%nrows, sp_matrix%nterms, &
             sp_matrix%row, sp_matrix%col, sp_matrix%val, &
             vector, solution)
    end subroutine solve_system_of_equation

    function get_term_index(sp_matrix, row, col) result(ind)
        ! Returns the sparse matrix storage index for the term in row and col.
        ! If term is not stored in sp_matrix (zero entry not stored, or is
        ! above the diagonal (not stored either), returns 0.
        implicit none
        type(sparse_matrix), intent(in)::sp_matrix
        integer(4),intent(in)::row, col
        integer(4)::ind, inf, sup, nterms
        call row_data(sp_matrix, row, inf, sup, nterms)
        do ind=inf,sup
            if (sp_matrix%col(ind).eq.col) then
                return
            elseif (sp_matrix%col(ind).gt.col) then
                exit
            endif
        end do
        ind = 0
        return 
    end function get_term_index

    function get_term_value(sp_matrix, row, col) result(val)
        ! Returns the sparse matrix storage value for the term in row and col.
        ! If term is not stored in sp_matrix (zero entry not stored, or is
        ! above the diagonal (not stored either), returns 0.0d0.
        implicit none
        type(sparse_matrix), intent(in)::sp_matrix
        integer(4),intent(in)::row, col
        integer(4)::ind
        real(8)::val
        ind = get_term_index(sp_matrix, row, col)
        if (ind.ne.0) then
            val = sp_matrix%val(ind)
        else
            val = 0.0d0
        endif
    end function get_term_value

    subroutine sum_value_in_term(sp_matrix,row,col,value)
        ! Sums the given value in the sparse matrix term, if it exists
        implicit none
        type(sparse_matrix),intent(inout)::sp_matrix
        integer(4),intent(in)::row, col
        real(8),intent(in)::value
        integer(4)::ind
        ind = get_term_index(sp_matrix, row, col)
        sp_matrix%val(ind) = sp_matrix%val(ind) + value
    end subroutine sum_value_in_term

    subroutine set_value_in_term(sp_matrix,row,col,value)
        ! Sets the given value in the sparse matrix term, if it exists
        implicit none
        type(sparse_matrix),intent(inout)::sp_matrix
        integer(4),intent(in)::row, col
        real(8),intent(in)::value
        integer(4)::ind
        ind = get_term_index(sp_matrix, row, col)
        sp_matrix%val(ind) = value
    end subroutine set_value_in_term

    subroutine penalize_term(sp_matrix,ind)
        implicit none
        type(sparse_matrix),intent(inout)::sp_matrix
        integer(4),intent(in)::ind   
        sp_matrix%val(sp_matrix%cum(ind)) = 10.0d12
    end subroutine penalize_term

    subroutine get_row(sp_matrix, row)
        ! Returns the row terms of the sparse matrix. The 'row' parameter is
        ! the row number. The row terms are returned in the global variables
        ! sp_row, sp_col and sp_data, from index 1 to the index 'sp_nterms'.
        ! WARNING: Only the terms above the diagonal are return! If it is
        ! needed the whole row, one must also get the terms under diagonal
        ! by getting the column of the same index of row and exchange row and
        ! col indexes.
        implicit none
        type(sparse_matrix), intent(in)::sp_matrix
        integer(4),intent(in)::row
        integer(4)::ind, inf, sup, nterms, counter
        call row_data(sp_matrix, row, inf, sup, nterms)
        counter = 0
        do ind = inf,sup
            counter = counter + 1
            sp_row(counter) = sp_matrix%row(ind)
            sp_col(counter) = sp_matrix%col(ind)
            sp_val(counter) = sp_matrix%val(ind)
        enddo
        sp_nterms = nterms
    end subroutine get_row

    subroutine get_col(sp_matrix, col)
        ! Returns the col terms of the sparse matrix. The 'col' parameter is
        ! the row number. The row terms are returned in the global variables
        ! sp_row, sp_col and sp_data, from index 1 to the index 'sp_nterms'.
        ! WARNING: Only the terms above the diagonal are return! If it is
        ! needed the whole column, one must also get the terms under diagonal
        ! by getting the row of the same index of col and exchange row and
        ! col indexes.
       implicit none
        type(sparse_matrix), intent(in)::sp_matrix
        integer(4),intent(in)::col
        integer(4)::ind, inf, sup, nterms, row_ind, counter
        counter = 0
        do row_ind = 1, col
            call row_data(sp_matrix, row_ind, inf, sup, nterms)
            if ((sp_matrix%col(inf).le.col).and.(sp_matrix%col(sup).ge.col)) then
                do ind=inf,sup
                    if (sp_matrix%col(ind).eq.col) then
                        counter = counter + 1
                        sp_row(counter) = sp_matrix%row(ind)
                        sp_col(counter) = sp_matrix%col(ind)
                        sp_val(counter) = sp_matrix%val(ind)
                        exit
                    elseif (sp_matrix%col(ind).gt.col) then
                        exit
                    endif
                end do
            endif
        enddo
        sp_nterms = counter
    end subroutine get_col

    subroutine set_value_to_row(sp_matrix, row, value)
        ! Sets the given value to all sparse matrix terms in the given row.
        implicit none
        type(sparse_matrix), intent(inout)::sp_matrix
        integer(4),intent(in)::row
        real(8),intent(in)::value
        integer(4)::ind, inf, sup, nterms
        call row_data(sp_matrix, row, inf, sup, nterms)
        do ind = inf,sup
            sp_matrix%val(ind) = value
        enddo
    end subroutine set_value_to_row

    subroutine set_value_to_col(sp_matrix, col,value)
        ! Sets the given value to all sparse matrix terms in the given column.
        implicit none
        type(sparse_matrix),intent(inout)::sp_matrix
        integer(4),intent(in)::col
        real(8),intent(in)::value
        integer(4)::ind, inf, sup, nterms, row_ind
        do row_ind = 1, col
            call row_data(sp_matrix, row_ind, inf, sup, nterms)
            if ((sp_matrix%col(inf).le.col).and.(sp_matrix%col(sup).ge.col)) then
                do ind=inf,sup
                    if (sp_matrix%col(ind).eq.col) then
                        sp_matrix%val(ind) = value
                        exit
                    elseif (sp_matrix%col(ind).gt.col) then
                        exit
                    endif
                end do
            endif
        enddo
    end subroutine set_value_to_col

    function get_diagonal(sp_matrix) result(diag)
        implicit none
        type(sparse_matrix), intent(in)::sp_matrix
        real(8), dimension(sp_matrix%nrows)::diag
        integer(4)::i,idiag
        diag = 0.0d0
        do i=1,sp_matrix%nrows
            idiag = sp_matrix%cum(i)
            if (sp_matrix%row(idiag).eq.sp_matrix%col(idiag)) then
                diag(i) = sp_matrix%val(idiag)
            endif
        enddo
    end function get_diagonal

    subroutine solve_system_of_equation_cg(A, b, x)
        ! Solve system of equation using conjugate gradient
        implicit none
        ! Dummy Arguments
        type(sparse_matrix),intent(in)::A
        real(8),dimension(:),intent(in)::b
        real(8),dimension(:),intent(out)::x
        ! Local variables
        integer(4)::k, maxit, t
        real(8),dimension(:),allocatable::r,z,p,Ap,d,Ax
        real(8)::alpha, beta, r_x_z, r_x_z_next
        real(8)::ti, tf
        real(8)::TOLERANCE=1.0d-9
        t = size(b)
        maxit = int(1.1*t)
        allocate(r(t),d(t),Ax(t))
        ! First steps
        d = get_diagonal(A)
        x = b / d
        call dot_matrix_vector(A, x, Ax)
        r = b - Ax
        if (maxval(dabs(r)).lt.TOLERANCE) then
            deallocate(r, d, Ax)
            return
        endif
        allocate(z(t),p(t),Ap(t))
        k = 0
        z = r / d
        p = z
        call cpu_time(ti)
        r_x_z = dot_product(r,z)
        do while (k < maxit)
            call dot_matrix_vector(A, p, Ap)
            alpha = r_x_z / dot_product(p,Ap)
            x = x + alpha * p
            r = r - alpha * Ap
            if (maxval(dabs(r)).lt.TOLERANCE) exit
            z = r / d
            r_x_z_next = dot_product(r,z)
            beta = r_x_z_next / r_x_z
            p = z + beta * p
            k = k + 1
            r_x_z = r_x_z_next
        enddo
        call cpu_time(tf)
        ! print *, "residual vector -> ",r
        deallocate(r,p,z,d,Ap,Ax)
        if (k < maxit) then
            print 10, k,(tf-ti)
        else
            print 20, k
        endif
        10 format("Conjugate gradient solver converged after ",i5," iterations in ",f9.3," seconds!")
        20 format("Conjugate gradient solver did not converge in ",i5," iterations!")
    end subroutine solve_system_of_equation_cg

end module sparse