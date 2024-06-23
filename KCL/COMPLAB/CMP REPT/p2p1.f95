program montecarlointegration
implicit none
real (8) :: pi,r,s
real (8) :: x,y
integer (4) :: Nhit, Ntot, i, N

N=10000
Nhit=0
Ntot=0
pi= acos(dble(-1))


do i= 1,N

call random_number (r)
call random_number (s)

x=r*(pi/2)
y=s

if ((sin(x))**4>y) then

Nhit=Nhit+1

endif

Ntot=Ntot+1
end do

write(*,*) 'Ntot =', Ntot
write(*,*) 'Nhit =', Nhit
write(*,*) 'Area =', (dble(Nhit)/dble(Ntot))*(pi/2)

end program


