module hammer_module
    implicit none
    contains

    subroutine hammer(nph, hammer_points)
        implicit none
        
        integer, intent(out) :: nph
        real(8), allocatable, intent(out) :: hammer_points(:,:)

        ! nph = 12
        ! allocate(hammer_points(3, 12))

        ! Assign the number of points
        

        ! ! Assign values to the hammer array in row-major order (like Python)
        ! hammer_points(1,:) = [0.501426509658179, 0.249286745170910, 0.249286745170910, &
        !     0.873821971016996, 0.063089014491502, 0.063089014491502, &
        !     0.053145049844816, 0.310352451033785, 0.636502499121399, &
        !     0.310352451033785, 0.636502499121399, 0.053145049844816]

        ! hammer_points(2,:) = [0.249286745170910, 0.249286745170910, 0.501426509658179, &
        !     0.063089014491502, 0.063089014491502, 0.873821971016996, &
        !     0.310352451033785, 0.636502499121399, 0.053145049844816, &
        !     0.053145049844816, 0.310352451033785, 0.636502499121399]

        ! hammer_points(3,:) = [0.116786275726379/2.0, 0.116786275726379/2.0, 0.116786275726379/2.0, &
        !     0.050844906370207/2.0, 0.050844906370207/2.0, 0.050844906370207/2.0, &
        !     0.082851075618374/2.0, 0.082851075618374/2.0, 0.082851075618374/2.0, &
        !     0.082851075618374/2.0, 0.082851075618374/2.0, 0.082851075618374/2.0]


        ! nph = 7
        ! allocate(hammer_points(3, nph))
        
        !     hammer_points(1,1) = 1.d0/3.d0
        !     hammer_points(2,1) = 1.d0/3.d0
        !     ! hammer_points(3,3) = 1.d0/3.d0
        !     hammer_points(3,1) = 0.11250d0
    
        !     hammer_points(1,2) = 0.797426985353087d0
        !     hammer_points(2,2) = 0.101286507323456d0
        !     ! hammer_points(2,3) = 0.101286507323456d0
        !     hammer_points(3,2) = 0.125939180544827d0/2.d0
    
        !     hammer_points(1,3) = 0.101286507323456d0
        !     hammer_points(2,3) = 0.797426985353087d0
        !     ! hammer_points(3,3) = 0.101286507323456d0
        !     hammer_points(3,3) = 0.125939180544827d0/2.d0
    
        !     hammer_points(1,4) = 0.101286507323456d0
        !     hammer_points(2,4) = 0.101286507323456d0
        !     ! hammer_points(4,3) = 0.797426985353087d0
        !     hammer_points(3,4) = 0.125939180544827d0/2.d0
    
        !     hammer_points(1,5) = 0.470142064105115d0
        !     hammer_points(2,5) = 0.470142064105115d0
        !     ! hammer_points(5,3) = 0.059715871789770d0
        !     hammer_points(3,5) = 0.132394152788506d0/2.d0
    
        !     hammer_points(1,6) = 0.059715871789770d0
        !     hammer_points(2,6) = 0.470142064105115d0
        !     ! hammer_points(6,3) = 0.470142064105115d0
        !     hammer_points(3,6) = 0.132394152788506d0/2.d0
    
        !     hammer_points(1,7) = 0.470142064105115d0
        !     hammer_points(2,7) = 0.059715871789770d0
        !     ! hammer_points(7,3) = 0.470142064105115d0
        !     hammer_points(3,7) = 0.132394152788506d0/2.d0


    nph = 3
    allocate(hammer_points(3, nph))

    hammer_points(1,1) = 1.d0/6.d0
    hammer_points(2,1) = 1.d0/6.d0
    hammer_points(3,1) = 1.d0/6.d0

    hammer_points(1,2) = 2.d0/3.d0
    hammer_points(2,2) = 1.d0/6.d0
    hammer_points(3,2) = 1.d0/6.d0

    hammer_points(1,3) = 1.d0/6.d0
    hammer_points(2,3) = 2.d0/3.d0
    hammer_points(3,3) = 1.d0/6.d0
    

        
    end subroutine hammer


end module hammer_module
