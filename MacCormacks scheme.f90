program Maccorrmack

implicit none
!program to compute quasi 1D isentropic Nozzle flow in a converging diverging nozzle
!this is the exact probelm from Anderson's textbook on CFD 
integer,parameter::n=31  !Number of grid points
integer::i,xmin,xmax,dx_,time,nsteps
real::gammar,k,s,C,dx,dt

real,dimension(1:n)::rho,T,v,A,as,rho_pr,T_pr,Tpr,v_pr,rho_cor,T_cor,v_cor,rhoav,Tav,vav,x
real,dimension(1:n)::deltat,M,rho_old,v_old,T_old,P


xmin=0.0
xmax=3.0

dx=real((xmax-xmin))/(n-1)
dx_=(xmax-xmin)/(n-1)
C=0.5 ! Courant number

gammar=1.40
k=1/gammar
s=gammar-1


!initial conditions for the Nozzle

do i=1,n

x(i)=(i-1)*dx
A(i)=1+2.2*((x(i)-1.5)**2)
rho(i)=1-0.3146*(x(i))
T(i)=1-0.2314*(x(i))
v(i)=(0.1+1.09*x(i))*(T(i)**0.5)
as(i)=T(i)**0.5 
deltat(i)=(C*dx)/(as(i)+v(i))

end do

dt=MINVAL(deltat)


nsteps=1400


open(1,file="Mac.dat",status='replace')
open(2,file="Mac.plt",status='replace')



do time=1,nsteps

rho_old=rho
v_old=v
T_old=T

!Prediction step

do i=2,n-1

rho_pr(i)=(-rho(i)*((v(i+1)-v(i))/dx))- (rho(i)*v(i))*((LOG(A(i+1))-LOG(A(i)))/(dx))-v(i)*((rho(i+1)-rho(i))/(dx))

v_pr(i)=(-v(i))*((v(i+1)-v(i))/(dx))- (k*(((T(i+1)-T(i))/(dx))+((T(i)/rho(i))*((rho(i+1)-rho(i))/dx))))

T_pr(i)=(-v(i))*((T(i+1)-T(i))/dx)-  s*((T(i))*(((v(i+1)-v(i))/dx)+((v(i)*((LOG(A(i+1))-LOG(A(i)))/dx)))))

rho(i)=rho(i)+(rho_pr(i)*dt)

v(i)=v(i)+(v_pr(i)*dt)

T(i)=T(i)+(T_pr(i)*dt)

end do

!Corrector step


do i=2,n-1

 rho_cor(i)=-rho(i)*((v(i)-v(i-1))/(dx))-  (rho(i)*v(i))*((LOG(A(i))-LOG(A(i-1)))/(dx))  -(v(i)*((rho(i)-rho(i-1))/(dx)))

 v_cor(i)=(-v(i))*((v(i)-v(i-1))/(dx))-  (k*(((T(i)-T(i-1))/(dx))+((T(i)/rho(i))*((rho(i)-rho(i-1))/dx))))

 T_cor(i)=(-v(i))*((T(i)-T(i-1))/(dx))-  s*((T(i))*(((v(i)-v(i-1))/dx)+  ((v(i)*((LOG(A(i))-LOG(A(i-1)))/dx)))))

end do

!calculating average slope from the above steps

do i=2,n-1

rhoav(i)=0.5*(rho_pr(i)+rho_cor(i))
vav(i)=0.5*(v_pr(i)+v_cor(i))
Tav(i)=0.5*(T_pr(i)+T_cor(i))

rho(i)=rho_old(i)+(rhoav(i)*dt)
v(i)=v_old(i)+(vav(i)*dt)
T(i)=T_old(i)+(Tav(i)*dt)

end do

!Applying boundary conditions

v(1)=2*v(2)-v(3)
rho(1)=1.0
T(1)=1.0

v(n)=2*v(n-1)-v(n-2)
rho(n)=2*rho(n-1)-rho(n-2)
T(n)=2*T(n-1)-T(n-2)

!Calculating Mach Number and pressure
M=v/sqrt(T)

P=rho*T

!Writing final results to file
if(time.EQ.nsteps)then

do i=1,n

write(1,*)x(i),A(i),rho(i),v(i),T(i),P(i),M(i)

end do

end if

end do

!writing plot files and calling gnuplot

write(2,*)"set xlabel 'Nozzle length'"
write(2,*)"set ylabel 'Parameter'"
write(2,*)"set title 'Quasi 1D Isentropic flow through a CD nozzle'"
write(2,*)"set grid"
write(2,*)"set size ratio 1.5"
write(2,*)'plot "Mac.dat" using  1:3  with line lt rgb "blue" title "Density",\'
write(2,*)'"Mac.dat" using  1:4  with line lt rgb "green" title "velocity",\'
write(2,*)'"Mac.dat" using  1:5  with line lt rgb "cyan" title "Temperature",\'
write(2,*)'"Mac.dat" using  1:6  with line lt rgb "orange" title "Pressure",\'
write(2,*)'"Mac.dat" using  1:7  with line lt rgb "red" title "Mach number"'



call system('gnuplot -p Mac.plt')

close(1)
close(2)

end program



