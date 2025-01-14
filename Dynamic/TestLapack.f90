program linearEquation
implicit none

real(8)::A(3,3),b(3)
integer::i,j,pivot(3),info

A(1,1)=3.1d0
A(1,2)=1.3d0
A(1,3)=-5.7d0
A(2,1)=1.0d0
A(2,2)=-6.9d0
A(2,3)=5.8d0
A(3,1)=3.4d0
A(3,2)=7.2d0
A(3,3)=-8.8d0

b(1)=-1.3d0
b(2)=-0.1d0
b(3)=1.8d0

call dgesv(3,1,A,3,pivot,b,3,info)

do i=1,3
    print *, b(i)
end do

end program linearEquation
