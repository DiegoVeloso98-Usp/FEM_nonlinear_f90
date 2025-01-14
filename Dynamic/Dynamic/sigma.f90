module sigma_module
    use shape_module
    use read_data_module

    implicit none

    contains

    subroutine stressPiola(ELEMS, NODES,FFORMA, U, nnodes, nelems, npe, ep, gaprox,SIGMA)

    real, dimension(:,:), intent(in) :: ELEMS, NODES
    real(8), intent(inout) :: FFORMA(:,:)
    real(8), dimension(:), intent(in) :: U
    integer, intent(in) :: nnodes, nelems, npe, gaprox, ep

    real(8), dimension(:,:), allocatable :: SIGMA

    integer :: iel, ih, ine, jn, idir, jdir, in_,index
    real(8) :: ksi, eta, detjac
    real(8) ::  dxdksi, dxdeta, dydksi, dydeta, dudksi, dudeta, dvdksi, dvdeta, dksidx, dksidy, detadx, detady
    real(8) :: J(2,2), JINV(2,2), DX(2,4), DY(2,4), DFI(4, 2*npe) , MXX(2,2), MXY(2,2), MYX(2,2), MYY(2,2)
    real(8) :: du(npe), dv(npe), result_fi(npe), result_dfidksi(npe), result_dfideta(npe)
    real(8) :: work(2), COORD(npe,2)
    real(8) :: cx(npe), cy(npe), vfi(npe), KLOCAL(2*npe, 2*npe), FLOCAL(2*npe)
    real(8) :: d11, d12, d21, d22, d33, d44
    real(8) :: epsilon_x, epsilon_y, epsilon_xy, sigmax, sigmay, talxy
    real    :: h, bx, by
    integer :: ipiv(2), info, lwork


    do iel=1,nelems
        index=1
        do ine=1,npe
            cx(index)=NODES(int(ELEMS(iel,ine)-1),1)
            cy(index)=NODES(int(ELEMS(iel,ine)-1),2)
            du(index)=U(2*int(ELEMS(iel,ine)-1)+1)
            dv(index)=U(2*int(ELEMS(iel,ine)-1)+2)
        end do
        
        do ih=1,npe
            ksi=COORD(ih,1)
            eta=COORD(ih,2)


            result_fi=fi(FFORMA,ksi,eta,gaprox)
            result_dfidksi=dfidksi(FFORMA,ksi,eta,gaprox)
            result_dfideta=dfideta(FFORMA,ksi,eta,gaprox)

            dxdksi = dot_product(result_dfidksi, cx)
            dxdeta = dot_product(result_dfideta, cx)
            dydksi = dot_product(result_dfidksi, cy)
            dydeta = dot_product(result_dfideta, cy)

            dudksi = dot_product(result_dfidksi, du)
            dudeta = dot_product(result_dfideta, du)
            dvdksi = dot_product(result_dfidksi, dv)
            dvdeta = dot_product(result_dfideta, dv)

            J(1,1) = dxdksi
            J(1,2) = dxdeta
            J(2,1) = dydksi
            J(2,2) = dydeta

            call dgetrf(2, 2, J, 2, ipiv, info)
            call dgetri(2, J, 2, ipiv, work, lwork, info)
            JINV = J

            dksidx=JINV(1,1)
            dksidy=JINV(1,2)
            detadx=JINV(2,1)
            detady=JINV(2,2)
            
            call materialProperties(ELEMS, iel, npe, ep, d11, d12, d21, d22, d33, d44, h, bx, by)

            epsilon_x = dudksi*dksidx + dudeta*detadx
            epsilon_y = dvdksi*dksidy + dvdeta*detady
            epsilon_xy = 0.5*(dudksi*dksidy+dudeta*detady+dvdksi*dksidx+dvdeta*detadx)

            sigmax=d11*epsilon_x+d12*epsilon_y
            sigmay=d21*epsilon_x+d22*epsilon_y
            talxy=d33*epsilon_xy

            SIGMA(int(ELEMS(iel,ih)),:) = SIGMA(int(ELEMS(iel,ih)),:) + [sigmax, sigmay, talxy]
        end do
    end do
    
    end subroutine stressPiola

end module sigma_module