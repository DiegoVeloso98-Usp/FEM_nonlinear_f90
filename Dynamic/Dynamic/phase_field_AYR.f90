module phase_field
    implicit none
    
contains
    interface CDphasefield
        module procedure CDphasefield_function
        module procedure CDphasefield_function_gradient
    end interface
    private CDphasefield_function, CDphasefield_function_gradient


    interface phase_dissipation_residual
        module procedure phase_dissipation_residual_point
        module procedure phase_dissipation_residual_gauss
    end interface
    private phase_dissipation_residual_point, &
            phase_dissipation_residual_gauss
     !"FORCE" VECTOR qv
    subroutine PSOR
!===============================================================!
!===============================================================!

    subroutine CDphasefield_function(d, dg)
        ! element nodal phase-field 
        ! d   = (d_1,...,d_a,...,d_{nen}) with  nen = number of element nodes
        ! a = a-th element node
        ! nodal phase-field 
        ! d_a = (d_{a,1},...,d_{a,i},...,d_{a,npf}) with npf = number of phase field
        ! i = i-th phase field
        !  gauss phase-field
        !       dg = (dg_1,...,dg_i,...,dg_{npf})
      
        implicit none
      
        real(8), intent(in) :: d(nen)                 ! nodal phase field
        real(8), intent(out) :: dg(ng)                    ! gauss phase field
        integer :: g, a
      
        dg = 0d0
        do g = 1, ng
          do a = 1, nen
            dg(g) = dg(g) + N(a,g)*d(a)
          end do
        end do
      
        return
      end subroutine CDphasefield_function
      
      subroutine CDphasefield_function_gradient(e, d, dg, graddg)
        ! element nodal phase-field 
        ! d = (d_1,...,d_a,...,d_{nen}) with nen = number of element nodes
        !   a = a-th element node
        ! nodal phase-field 
        ! d_a = (d_{a,1},...,d_{a,i},...,d_{a,npf})    
        !   with npf = number of phase field
        !   i = i-th phase field
        ! gauss phase-field
        !   dg = (dg_1,...,dg_i,...,dg_{npf})
        !   dg_i =(d_i, d_{i,x_1}, ..., d_{i,x_{idim}}, ..., d_{i,x_{ndim}})
      
        implicit none
      
        integer, intent(in) :: e                    ! element label
        real(8), intent(in) :: d(nen)               ! nodal phase field
        real(8), intent(out) :: dg(ng)              ! gauss phase field
        real(8), intent(out) :: graddg(ndim,ng)     ! gauss phase field gradient
        integer :: g, a, idim
      
        dg     = 0d0
        graddg = 0d0
        do g = 1, ng
          do a = 1, nen
            dg(g) = dg(g) + N(a,g)*d(a)
      
            do idim = 1, ndim
              graddg(idim,g) = graddg(idim,g) + DSH0(idim,a,g,e)*d(a)
            end do
      
          end do
        end do
      
        return
      end subroutine CDphasefield_function_gradient
    !===============================================================!
    !===============================================================!



    !================================================!
    call phase_dissipation_residual(dg, graddg, resdisf, resdisg)
    !================================================!

    !================================================!

    !================================================!


    !===============================================================!
    !===============================================================!
    subroutine phase_dissipation_residual_gauss(d, gradd, disf, disg)
        implicit none
        
        real(8), intent(in) :: d(:)       ! phase field
        real(8), intent(in) :: gradd(:,:)   ! phase field gradient
        real(8), intent(out) :: disf(:)   ! dissipation residual function
        real(8), intent(out) :: disg(:,:)   ! dissipation residual gradient
        integer :: g, ng
        
        ng = size(d,1)
        
        do g = 1, ng
            call phase_dissipation_residual_point(d(g), gradd(:,g), disf(g), disg(:,g))
        end do
        
        return
    end subroutine phase_dissipation_residual_gauss
    !===============================================================!
    !===============================================================!


    !===============================================================!
    !===============================================================!
    subroutine phase_dissipation_residual_point(d, gradd, disf, disg)
        ! -----------------------------------------------------
        !    OUTPUT (first derivarives of dissipation)
        !    function term Gc/l0/cw*d(w(di))/d(di)
        !    gradient term Gc*l0/cw*d(|grad(di)|^2)/d(di,xidim)
        ! -----------------------------------------------------
      
        implicit none
        
        real(8), intent(in) :: d              ! phase field
        real(8), intent(in) :: gradd(ndim)    ! phase field gradient
        real(8), intent(out) :: disf          ! dissipation residual function
        real(8), intent(out) :: disg(ndim)    ! dissipation residual gradient
        integer :: idim
      
        select case(AT)
            
        case(1)
      
          disf = Gc/l0/cw(AT)
          do idim = 1, ndim
            disg(idim) = 2d0*Gc*l0/cw(AT)*gradd(idim)
          end do
          
        case(2)
        
          disf = 2d0*Gc/l0/cw(AT)*d
          do idim = 1, ndim
            disg(idim) = 2d0*Gc*l0/cw(AT)*gradd(idim)
          end do
          
        end select
        
        return
      end subroutine phase_dissipation_residual_point
      !===============================================================!
      !===============================================================!

    
end module phase_field
    
   
