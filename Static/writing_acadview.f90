module writing_acadview_module
    implicit none

    integer :: unit_number
    logical :: file_opened = .false.
    contains


    subroutine initialize_file(nelems, nnodes, npe, gaprox, ELEMS, NODES)
        implicit none

        ! Input parameters
        integer, intent(in) :: nelems, nnodes, npe, gaprox
        real, intent(in) :: ELEMS(:,:)
        real(8), intent(in) :: NODES(:,:)

        ! Local variables
        integer :: i, j, total_nodes, total_elements
        integer :: io_status
        character(len=100) :: path
        integer, allocatable :: list(:)

        total_nodes = nnodes
        total_elements = nelems

        ! Open the file only once
        if (.not. file_opened) then
            path = 'output.ogl'
            open(newunit=unit_number, file=path, status='replace', action='write', iostat=io_status)
            if (io_status /= 0) then
                print *, "Error opening file."
                stop
            end if
            file_opened = .true.

            ! Write header information
            write(unit_number, '(A)') 'arquivo de entrada'
            write(unit_number, '(A)') 'n.nos n.elems n.listas'
            write(unit_number, '(A)') '#'
            write(unit_number, '(I10, I10, I10)') total_nodes, total_elements, 5
            write(unit_number, '(A)') 'coordx coordy coordz deslx desly deslz'
            write(unit_number, '(A)') '#'

            ! Write NODES
            do i = 1, nnodes
                write(unit_number, '(F18.8 , F18.8 , F18.8 , F18.8 , F18.8 , F18.8)') &
                    NODES(i, 1), NODES(i, 2), 0.0, 0.0, 0.0, 0.0
            end do

            ! Write element information
            write(unit_number, '(A)') 'tpelem (1-barra/2-triang/3-quad) grauaprox nó1 nó2...nó_n group'
            write(unit_number, '(A)') '#'

            ! Write ELEMS
            do i = 1, nelems
                allocate(list(npe + 3))
                list(1) = 2
                list(2) = gaprox
                do j = 1, npe
                    list(j + 2) = int(ELEMS(i, j))
                end do
                list(npe + 3) = 0
                write(unit_number, '(13I10)') list
                deallocate(list)
            end do
        end if

    end subroutine initialize_file

    subroutine appending(CAUCHY,VonMises, PrincipalStress,  U , d_Global, nnodes)
        implicit none

        ! Input parameters
        integer, intent(in) :: nnodes
        real(8), intent(in) :: U(:),CAUCHY(:,:),VonMises(:),PrincipalStress(:,:),d_Global(:)

        ! Local variables
        integer :: i

        ! Write displacement in x-direction
        ! write(unit_number, '(A)') '#'
        ! write(unit_number, '(A)') 'desl.x'
        ! do i = 1, nnodes
        !     write(unit_number, '(4F18.8)') U(2*i-1), U(2*i), 0.0, U(2*i-1)
        ! end do

        ! ! Write displacement in y-direction
        ! write(unit_number, '(A)') '#'
        ! write(unit_number, '(A)') 'desl.y'
        ! do i = 1, nnodes
        !     write(unit_number, '(4F18.8)') U(2*i-1), U(2*i), 0.0, U(2*i)
        ! end do

        ! write(unit_number, '(A)') '#'
        ! write(unit_number, '(A)') 'sigma.x'
        ! do i = 1, nnodes
        !     write(unit_number, '(4F18.8)') U(2*i-1), U(2*i), 0.0, CAUCHY(i,1)/CAUCHY(i,4)
        ! end do

        write(unit_number, '(A)') '#'
        write(unit_number, '(A)') 'sigma.y'
        do i = 1, nnodes
            write(unit_number, '(4F18.8)') U(2*i-1), U(2*i), 0.0, CAUCHY(i,2)/CAUCHY(i,4)
        end do

        ! write(unit_number, '(A)') '#'
        ! write(unit_number, '(A)') 'sigma.xy'
        ! do i = 1, nnodes
        !     write(unit_number, '(4F18.8)') U(2*i-1), U(2*i), 0.0, CAUCHY(i,3)/CAUCHY(i,4)
        ! end do

        ! write(unit_number, '(A)') '#'
        ! write(unit_number, '(A)') 'VonMises'
        ! do i = 1, nnodes
        !     write(unit_number, '(4F18.8)') U(2*i-1), U(2*i), 0.0, VonMises(i)
        ! end do



        ! write(unit_number, '(A)') '#'
        ! write(unit_number, '(A)') 'Principal Stress 1'
        ! do i = 1, nnodes
        !     write(unit_number, '(4F18.8)') U(2*i-1), U(2*i), 0.0, PrincipalStress(i,1)
        ! end do

        ! write(unit_number, '(A)') '#'
        ! write(unit_number, '(A)') 'Principal Stress 2'
        ! do i = 1, nnodes
        !     write(unit_number, '(4F18.8)') U(2*i-1), U(2*i), 0.0, PrincipalStress(i,2)
        ! end do

        
        write(unit_number, '(A)') '#'
        write(unit_number, '(A)') 'Phase-Field'
        do i = 1, nnodes
            write(unit_number, '(4F18.8)') U(2*i-1), U(2*i), 0.0, d_Global(i)
        end do

    end subroutine appending

    subroutine close_file()
        implicit none

        ! Close the file if it's open
        if (file_opened) then
            close(unit_number)
            file_opened = .false.
            print *, "File closed successfully."
        end if

    end subroutine close_file


    subroutine Plot_Displacement_At_Load(U,LOAD,VINC,time,analysis_time)
        implicit none

        real(8), intent(in) :: U(:)
        real, intent(in) :: LOAD(:,:)
        real, intent(in) :: VINC(:,:)
        integer, intent(in) :: time, analysis_time

        real(8) :: Ux , Uy

        open(10,file='deslocamento.txt',status='unknown')
        
        ! CONSIDERING THAT JUST 1 NODE IS LOADED

        ! Ux = U(2*int( VINC(5,1) ) - 1)
        ! Uy = U(2*int( LOAD(2,1) ))

        Ux = (U(2*318 - 1)+U(2*317-1))/2.0d0
        Uy = (U(2*318)+U(2*317))/2.0d0

        write(10,*) time , Ux , Uy

        if (time == analysis_time) then
            close(10)
        end if
    
    end subroutine Plot_Displacement_At_Load

    
end module writing_acadview_module