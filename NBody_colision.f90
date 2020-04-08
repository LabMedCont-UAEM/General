program NBODY
implicit none

integer, parameter :: NTDOM=200
real, dimension(:), allocatable :: X1,X2,Y1,Y2,U1,U2,V1,V2,FX,FY
integer :: Nplot,t,Tmax,i,n
real :: dt,xrand,yrand,thrand,PI,MAGU,LX,LY,R0,PSI,raction,dist,time

ALLOCATE(X1(1:NTDOM),X2(1:NTDOM),Y1(1:NTDOM),Y2(1:NTDOM),&
	U1(1:NTDOM),U2(1:NTDOM),V1(1:NTDOM),V2(1:NTDOM),&
	FX(1:NTDOM),FY(1:NTDOM))
LX=1.0
LY=1.0
R0=0.0
PSI=1e-4
raction=0.05

MAGU=0.05
PI=ACOS(-1.0)

dt=1.0e-4
Tmax=60.0/dt
Nplot=500

time=0.0
do n=1,NTDOM
time=time+dt
CALL RANDOM_NUMBER(xrand)
CALL RANDOM_NUMBER(yrand)
CALL RANDOM_NUMBER(thrand)
x1(n)=LX*XRAND
y1(n)=LY*YRAND
u1(n)=MAGU*cos(2.0*PI*thrand)
v1(n)=MAGU*sin(2.0*PI*thrand)
enddo

x2=x1 ; y2=y1
u2=u1 ; v2=v1

do t=1,Tmax
FX=0.0 ; FY=0.0


do i=1,NTDOM

do n=1,NTDOM
dist=sqrt( (x1(i)-x1(n))**2 + (y1(i)-y1(n))**2) 
IF((dist.le.raction).and.(i.ne.n)) THEN
FX(i)=FX(i)+PSI*( ((raction/dist)**(4))- ((raction/dist)**(2)))*(x1(i)-x1(n))/(dist*dist)
FY(i)=FY(i)+PSI*( ((raction/dist)**(4))- ((raction/dist)**(2)))*(Y1(i)-Y1(n))/(dist*dist)
END IF
enddo

X2(i)=X1(i)+dt*U2(i)+dt*(FX(i))
Y2(i)=Y1(i)+dt*V2(i)+dt*(FY(i))

IF((X2(i).le.R0).or.(X2(i).ge.(Lx-R0))) THEN
	U2(i)=-(1.0)*U2(i)
	V2(i)=(1.0-0.0)*V2(i)
ELSE IF((Y2(i).le.R0).or.(Y2(i).ge.(Ly-R0))) THEN
	U2(i)=(1.0-0.0)*U2(i)
	V2(i)=-(1.0)*V2(i)
END IF
enddo

1 format('Time=',F6.3x,'seg, U=',F14.9,' m/s',1x,', V=',F14.9,' m/s')
if(mod(t,Nplot).eq.0) then
!print 1, time,MAXVAL(U2),MAXVAL(V1)
call ANIMACION(X1,Y1,U1,V1,Nplot,t,NTDOM)
end if

U1=U2
V1=V2
X1=X2
Y1=Y2
enddo

DEALLOCATE(X1,X2,Y1,Y2,U1,U2,V1,V2)
print*, 'End program'
end program

subroutine ANIMACION(X1,Y1,U1,V1,Nplot,t,NTDOM)
implicit none

integer, intent(in) :: Nplot, t,NTDOM
real, dimension(1:NTDOM), intent(inout) :: X1,Y1,U1,V1
integer :: n

integer :: NCOUNT
character :: EXT*4, DESTINY*512, NCTEXT*4, OUTFILE*512, FNAME*16

EXT='.csv'
DESTINY='/home/garzon/Escritorio/anim/'
FNAME='Frame_'
NCOUNT=t/Nplot
WRITE(NCTEXT,'(I4.4)') NCOUNT 

OUTFILE=trim(DESTINY)//trim(FNAME)//trim(NCTEXT)//EXT

OPEN(10,file=OUTFILE)
100 format(I6,',',E12.4,',',E12.4,',',E12.4,',',E12.4)
write(10,*) 'NPART,XPOS,YPOS,UVEL,VVEL'
do n=1,NTDOM
write(10,100) n,X1(n), Y1(n), U1(n), V1(n)
enddo
close(10)
return
end subroutine
