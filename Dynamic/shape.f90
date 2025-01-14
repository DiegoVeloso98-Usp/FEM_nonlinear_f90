module shape_module
    implicit none
    contains

    function fi (FFORMA,KSI,ETA,gaprox) result(evaluated_result)
        implicit none
        real(8),dimension(:,:),intent(in)::FFORMA
        real(8),intent(in)::KSI,ETA
        real(8),dimension(((gaprox+1)*(gaprox+2)/2))::evaluated_result
        real(8),dimension(((gaprox+1)*(gaprox+2)/2))::evaluated_coef
        integer,intent(in) :: gaprox
        integer::il,ic,index

        index=0
        do il=1,(gaprox+1)
            do ic=1,(gaprox+1-il+1)
                index=index+1
                evaluated_coef(index)=ksi**(ic-1)*eta**(il-1)
            end do
        end do

        evaluated_result=matmul(FFORMA,evaluated_coef)
        ! print *,"fi inside shape_mod", evaluated_result
    end function fi


    function dfidksi (FFORMA,KSI,ETA,gaprox) result(evaluated_result)
        implicit none
        real(8),dimension(:,:),intent(in)::FFORMA
        real(8),intent(in)::KSI,ETA
        real(8),dimension(((gaprox+1)*(gaprox+2)/2))::evaluated_result
        real(8),dimension(((gaprox+1)*(gaprox+2)/2))::evaluated_coef
        integer,intent(in) :: gaprox
        integer::il,ic,index

        index=0
        do il=1,(gaprox+1)
            do ic=1,(gaprox+1-il+1)
                index=index+1
                if (ic>1) then
                    evaluated_coef(index)=(ic-1)*ksi**(ic-2)*eta**(il-1)
                else
                    evaluated_coef(index)=0.0
                end if
            end do
        end do
        evaluated_result=matmul(FFORMA,evaluated_coef)

    end function dfidksi

    PURE function dfideta (FFORMA,KSI,ETA,gaprox) result(evaluated_result)
        implicit none
        real(8),dimension(:,:),intent(in)::FFORMA
        real(8),intent(in)::ksi,eta
        real(8),dimension(((gaprox+1)*(gaprox+2)/2))::evaluated_result
        real(8),dimension(((gaprox+1)*(gaprox+2)/2))::evaluated_coef
        integer,intent(in) :: gaprox
        integer::il,ic,index

        index=0
        do il=1,(gaprox+1)
            do ic=1,(gaprox+1-il+1)
                index=index+1
                if (il>1) then
                    evaluated_coef(index)=ksi**(ic-1)*(il-1)*eta**(il-2)
                else
                    evaluated_coef(index)=0.0
                end if
            end do
        end do

        evaluated_result=matmul(FFORMA,evaluated_coef)
    end function dfideta


    subroutine shapeFunc(gaprox,npe,FFORMA,COORD)
        implicit none
        integer,intent(in)::gaprox,npe
        real(8),dimension(npe,npe), intent(out)::FFORMA

        real(8),dimension(npe,2),intent(out)::COORD
        real(8),dimension(npe,npe)::MAT,MATCOPY

        integer::il,ic,i,ino,index,info
        real(8),dimension(npe)::vec
        real(8)::termo1,termo2,ksi,eta

        integer, allocatable :: IPIV(:)


        index=0
        do il=1,(gaprox+1)
            do ic=1,(gaprox+1-il+1)
                index=index+1
                ksi=1.0d0/gaprox*(ic-1)
                eta=1.0d0/gaprox*(il-1)
                COORD(index,1)=ksi
                COORD(index,2)=eta
            end do
        end do

        index=0
        do il=1,(gaprox+1)
            do ic=1,(gaprox+1-il+1)
                index=index+1
                do ino=1,npe
                    ksi=COORD(ino,1)
                    eta=COORD(ino,2)
                    if (ic==1) then
                        termo1=1
                    else
                        termo1=ksi**(ic-1)
                    end if
                    if (il==1) then
                        termo2=1
                    else
                        termo2=eta**(il-1)
                    end if
                    MAT(ino,index)=termo1*termo2
                end do
            end do
        end do
        
        allocate(ipiv(npe))
        do i=1,npe
            vec=0.0d0
            vec(i)=1.0d0
            MATCOPY=MAT
            call dgesv(npe, 1, MATCOPY, npe, IPIV, vec, npe, INFO)
            FFORMA(i,:) = vec
        end do
    end subroutine shapeFunc

end module shape_module