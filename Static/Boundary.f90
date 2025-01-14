module Boundary_module
    use sparse
    implicit none
    
contains
subroutine globalForce(LOAD, FGLOBAL,step,nLoad_steps)
    implicit none

    real, intent(in) :: LOAD(:,:)
    real(8), intent(inout) :: FGLOBAL(:)
    integer :: i, ino, idir
    integer, intent(in) :: step, nLoad_steps
    real:: value
    ! print*, "Applying global forces"
    do i = 1, size(LOAD, 1)
        ino=int(LOAD(i, 1))
        idir=int(LOAD(i, 2))
        value=LOAD(i, 3)

        FGLOBAL(2*(ino-1)+idir)=value*real(step,kind=8)/real(nLoad_steps,kind=8)
    end do
end subroutine globalForce

subroutine boundaryCondition(VINC, KGLOBAL, FGLOBAL)
    implicit none

    real, intent(in) :: VINC(:,:)
    real(8), intent(inout) ::  FGLOBAL(:)
    type(sparse_matrix), intent(inout):: KGLOBAL
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
    ! Print*, "Boundary conditions applied ! "
end subroutine boundaryCondition

subroutine prescribedDisplacement(VINC, NODES_1,nLoad_steps)
    implicit none
    
    real, intent(in) :: VINC(:,:)
    real(8), intent(inout) ::  NODES_1(:,:)
    integer, intent(in) :: nLoad_steps
    integer :: i, ino, idir, type_
    real(8):: value_ !, nlg
    

    ! Print*, "Applying Prescribed Displacement"

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
    ! Print*, "Prescribed Displacements Applied ! "
end subroutine prescribedDisplacement


subroutine DirichletBoundaryCondition(dn , di , delta_d_i)
    implicit none
    integer :: i
    real(8), intent(inout) :: di(:) 
    real(8), intent(inout) ::  delta_d_i(:)

    real(8), intent(in) :: dn(:)
    real(8) :: dmax
    real(8), parameter :: d_tol = 0.9999999999d0
    dmax=maxval(di)
    if (dmax>=d_tol) then
        do i=1,size(di)
            if (di(i)>=d_tol) then
                delta_d_i(i)= max(d_tol-dn(i),0.0d0)
                di(i) = d_tol
            end if
        end do
    end if
end subroutine DirichletBoundaryCondition
    
end module Boundary_module