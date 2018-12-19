program Laplace
include 'mpi.h'
implicit none

integer :: ierr, numtasks, taskid

call mpi_init(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)


integer ::  taskid, numtasks, ierr ,M
integer ::  No, N, ids,ide,jds,jde
integer :: ite,its,jte,jts,ims,ime, jms 
integer :: jme,isize, jsize,ib,jb
integer :: i,j 

ids = 1
ide =1000
jds = 1
jde = 1100
! who am I? what is my 



call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
! find the largest integer N such that N^2 <= numtasks


No = floor(sqrt(dble(numtasks)))
N = INT(No)
M = INT(N)

print *,'task ',taskid,' domain ',ids,':',ide,' ',jds,':',jde,' tasks ',numtask, &
  ' proc grid ',M,' by ',N

!Indices for each tile
jb = 1+taskid/M
! M=5 taskid = 0: => jb = 1+(0)/5 = 1 + 0 = 1
! M=5 taskid = 4: => jb = 1+(4)/5 = 1 + 0 = 1
! M=5 taskid = 5: => jb = 1+(5)/5 = 1 + 1 = 2
ib = 1+mod(taskid,M)


isize = divideup(ide-ids+1, M)
jsize = divideup(jde-jds+1, N)

! compute my ib,jb,its,ite,jts,jte
its = ids + isize*(ib-1)
jts = jds + jsize*(jb-1)

ite = min(ide,ids - 1 + isize*ib)
jte = min(jde,jds - 1 + jsize*jb)

! set ims,ime,jms,jme

ims = its-1
ime = ite+1
jms = jts-1
jme = jte+1



!Defining arrays to store boundary conditions
REAL, DIMENSION(100), TARGET :: tbc, bbc
real, dimension(ime), target :: leftbc, rbc

!Defining pointers to assign to each targer

REAL, DIMENSION(100), POINTER :: ptt, ptb
real, dimension(ime), pointer :: ptleft, ptr

ptt => tbc
ptb => bbc
ptleft => leftbc
ptr => rbc


!Defining  boundary conditions for the top of the domain
do i = 1,jme
tbc(i) = 100
end do

!Defining  boundary conditions for the bottom of the domain


do i = 1,jme
bbc(i) = 100
end do

!Defining  boundary conditions for the left side of the domain


do i = 1,ime
lbc(i) = 100
end do

!Defining  boundary conditions for the right side of the domain


do i = 1, ime
rbc(i) = 100
end do


! have ids:ide jds:jde = domain dimensions
! compute isize, jsize =
! tile size in directions i and j
! number of points in i direction
! divided by number of domains in i direction


print *,'task ',taskid,' block ',ib,jb, ' tile ', its,':',ite,' ',jts,':',jte


!possibly change to function
call work(taskid, ib, jb, its, ite, jts, jte, ims, ime, jms, jme, N, M, leftbc, rbc, tbc, bbc)

call MPI_FINALIZE(ierr)

end program Laplace

integer function divideup(m,n)
! return m/n rounded up
! divideup = ceil(real(m)/real(n))
divideup = (m+n-1)/n
! m =n : (m + n -1)/n = (2*n-1)/n = 1
! m =n+1 : (m + n -1)/n = (n+1 + n-1)/n = (2n)/n = 2
end function divideup







subroutine work(taskid, ib,jb,its,ite,jts,jte,ims,ime,jms,jme, N, M, leftbc, rbc, tbc, bbc)

integer:: i, j, taskid,ib, jb, its,ite,jts,jte    
! starting and ending indices of the submatrix I am updating
integer::  ims,ime,jms,jme, N, M, turn 
     
 ! starting and ending indices of the submatrix I have in memory - at least by 1 larger on each side

double precision:: omega,h, lbc, rbc tbc, bbc

turn = 0
omega = 0.25

h = 0.001


double precision, dimension(ime,jme), target :: u,f,v
double precision, dimension(ime,jme), pointer :: ptu,ptf,ptv

      ptu => u
      ptv => v
      ptf => f
 
! dynamic memory allocation

!  assign subdomains to processors in M by N grid

!                         jb  ->
!                   1               N 
!        1       |----------------------|
!                | 1    M+1             | 
!   ib |         | 2    M+2             | 
!      v         |                      |
!                | M    M+M  ....   M*N |
!        M       |-----------------------


print *,ib,jb
! I am subdomain ib,jb on processor k
! I am respondible for updating the array u with the solution for i = its:ite and j=jts:jte
! my arrays are dimensioned as ims:ime, jms:jme 
!
! initialize u,f
! do the following in iteration loop

! x: if jb>1 then receive u(its:ite,jts-1) from subdomain ib,jb-1 on processor k-M 
!    if jb<N then send u(its:ite,jte+1) to subdomain ib,jb+1 on processor k+M
! x: if jb=1 then assume that u(its:ite,jts-1) has boundary values - set them first!
! y: do not need to worry about corner because I have only 5 point difference scheme
! z: etc.

!                  jts                jte
!               y zzzzzzzzzzzzzzzzzzzzzzzz 
!         its   x |----------------------|----------------
! ib,jb-1       x |        ib,jb         |  ij, jb+1
!               x |                      |
!               x |                      |
!         ite   x |----------------------|----------------
!                 |                      |
!                 |       ib+1, jb



! main work loop - compute  v = u - omega * (f - A*u)
! need also values with i=its-1 and ite+1 and j=jts-1 and jte+1

if (turn .eq. 0) then
do i = ims, ime
do j = jms, jme
u(i,j) = 0
v(i,j) = 0
f(i,j) = 0
end do
end do
end if 
!Corner Cases
if (jb > 1) then

call MPI_RECV(u(its:ite, jms), ite, MPI_DOUBLE, taskid - M, 0, MPI_COMM_WORLD, ierr)
call MPI_SEND(u(its:ite, jts), ite, MPI_DOUBLE, taskid -M, 0, MPI_COMM_WORLD, ierr)

if (jb<N)

call MPI_RECV(u(its:ite, jme), ite, MPI_DOUBLE, taskid + M, 1, MPI_COMM_WORLD, ierr)
call MPI_SEND(u(its:ite, jte), ite, MPI_DOUBLE, taskid +M, 1, MPI_COMM_WORLD, ierr)

else if (jb==1) then

call MPI_RECV(u(its:ite, jme), ite, MPI_DOUBLE, taskid +M, 1, MPI_COMM_WORLD, ierr)
call MPI_SEND(u(its:ite, jte), ite, MPI_DOUBLE, taskid +M, 1, MPI_COMM_WORLD, ierr)

u(its:ime, jms) = lbc

else if (jb==N) then 

call MPI_RECV(u(its:ite, jms), ite, MPI_DOUBLE, taskid - M, 3, MPI_COMM_WORLD, ierr)
call MPI_SEND(u(its:ite, jts), ite, MPI_DOUBLE, taskid -M, 3, MPI_COMM_WORLD, ierr)

u(its:ime,jme) = rbc 

else if (ib>1) then

call MPI_RECV(u(ims, jts:jte), jte, MPI_DOUBLE, taskid - 1, 2, MPI_COMM_WORLD, ierr)
call MPI_SEND(u(its, jts:jte), jte, MPI_DOUBLE, taskid -1, 2, MPI_COMM_WORLD, ierr)
 
else if (ib<M) then

call MPI_RECV(u(ime, jts:jte), jte, MPI_DOUBLE, taskid + 1, 3, MPI_COMM_WORLD, ierr)
call MPI_SEND(u(ite, jts:jte), jte, MPI_DOUBLE, taskid +1, 3, MPI_COMM_WORLD, ierr)


else if (ib == 1) then

call MPI_RECV(u(ime, jts:jte), jte, MPI_DOUBLE, taskid - 1, 2, MPI_COMM_WORLD, ierr)
call MPI_SEND(u(ite, jts:jte), jte, MPI_DOUBLE, taskid -1, 2, MPI_COMM_WORLD, ierr)

u(ims,jms:jme) = tbc

else if (ib == M) then

call MPI_RECV(u(ims, jts:jte), jte, MPI_DOUBLE, taskid - 1, 2, MPI_COMM_WORLD, ierr)
call MPI_SEND(u(its, jts:jte), jte, MPI_DOUBLE, taskid -1, 2, MPI_COMM_WORLD, ierr)


u(ime, jms:jme) = bbc
end if

!Wait for all processes to reach this point
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
 
do j=jts,jte
      do i=its,ite
         v(i,j) = u(i,j) - omega * (f(i,j)*h*h - (4*u(i,j)-u(i+1,j) - u(i-1,j) - u(i,j+1) - u(i,j-1)))
      enddo
enddo

Print *, 'The tile with indices', ib, jb, 'has an array of:' v

! u = v
do j=jts,jte
    do i=its,ite
        u(i,j) = v(i,j)
    enddo
enddo

! compute also the error (residual, or the change in one iteration)
! add up the errors over subdomains and print the result

! end iteration loop here

! think how to get the results out to a file!

end subroutine work

