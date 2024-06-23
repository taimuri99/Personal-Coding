program montecarlosphere
implicit none
real (8) :: r,s,t,V
real (8) :: x,y,z
integer (4) :: Nhit, Ntot, i, N

N=10000
Nhit=0
Ntot=0


do i= 1,N

call random_number (r)
call random_number (s)
call random_number (t)

x= -1 + 2*r
y= -1 + 2*s
z= -1 + 2*t

V = x**2 + y**2 + z**2 

if (1**2>V .and. 0.75**2<V) then

 Nhit=Nhit+1 

endif

Ntot=Ntot+1

end do

write(*,*) 'Ntot =', Ntot
write(*,*) 'Nhit =', Nhit
write(*,*) 'Volume =', (dble(Nhit)/dble(Ntot))*(dble(8))

end program


