module Algebra_module
    implicit none
    
contains
subroutine matrix_inverse(M, Minv)
    implicit none
    ! Input arguments
    real(8), intent(in) :: M(:, :) ! Input matrix to be inverted
    
    ! Output arguments
    real(8), intent(out) :: Minv(:, :) ! Inverted matrix
    integer :: info               ! Output status (0 for success)
    
    ! Local variables
    integer,allocatable :: ipiv(:)           ! Pivot indices from dgetrf
    real(8),allocatable :: work(:)  ! Work array for dgetri
    integer :: lwork ,n           ! Optimal size of the work array
    
    ! Copy the input matrix to Minv (since LAPACK modifies the input matrix)
    n=size(M,1)
    allocate(ipiv(n))
    Minv = M
    
    ! Perform LU factorization using dgetrf
    call dgetrf(n, n, Minv, n, ipiv, info)
    if (info /= 0) then
        print *, "Error in LU factorization (dgetrf): info =", info
        return
    end if

    ! Allocate minimal workspace for the query
    allocate(work(1))
    lwork = -1
    call dgetri(n, Minv, n, ipiv, work, lwork, info)
    if (info /= 0) then
        print *, "Error querying workspace size for dgetri: info =", info
        return
    end if
    
    ! Optimal workspace size
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    
    ! Compute the inverse using dgetri
    call dgetri(n, Minv, n, ipiv, work, lwork, info)
    if (info /= 0) then
        print *, "Error in matrix inversion (dgetri): info =", info
        return
    end if

    ! Deallocate work array
! Deallocate work arrays
    deallocate(work)
    deallocate(ipiv)

end subroutine matrix_inverse

    
end module Algebra_module