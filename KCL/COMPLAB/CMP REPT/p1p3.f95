program problemone
implicit none

real (8) :: x0,x,a,VARa
integer (4) :: i,j,k

a=3.5
VARa= (3.8-3.5)/2000.0

do k=1,2000
a=a+VARa

      do j=1,50

            call random_number (x0)
            do i=1,100
            x=a*x0*(1-x0)
            x0=x 
            end do

            write (100,*) a,x

      end do

end do

end program
