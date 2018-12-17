program Laplace
include 'mpi.h'

integer ierr numtasks taskid

call mpi_init(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)

function create_array(ids, ide,jds, jde)





if (taskid == 0)
end if

a

call MPI_SEND(a, size(a),
function create_array(ids, ide,jds, jde)

function create_array(ids, ide,jds, jde)
integer:: taskid, numtasks, ierr, M, N, ids,ide,jds,jde
double precision :: lbc, rbc, tbc, bbc

lbc = 100
rbc = 100
tbc = 100
bbc = 100

ids = 1
ide =1000
jds = 1
jde = 1100
! who am I? what is my 



 call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
! find the largest integer N such that N^2 <= numtasks
 N = floor(sqrt(double(numtasks)))
 M = N

print *,'task ',taskid,' domain ',ids,':',ide,' ',jds,':',jde,' tasks ',numtask, &
  ' proc grid ',M,' by ',N
 ib = 1+taskid/M
! M=5 taskid = 0: => jb = 1+(0)/5 = 1 + 0 = 1
! M=5 taskid = 4: => jb = 1+(4)/5 = 1 + 0 = 1
! M=5 taskid = 5: => jb = 1+(5)/5 = 1 + 1 = 2
jb = 1+mod(taskid,M)


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

! have ids:ide jds:jde = domain dimensions
! compute isize, jsize = tile size in directions i and j
! number of points in i direction divided by number of domains in i direction
print *,'task ',taskid,' block ',ib,jb, ' tile ', its,':',ite,' ',jts,':',jte


!possibly change to function
call work(ib, jb, its, ite, jts, jte, ims, ime, jms, jme, N, M, lbc, rbc, tbc, bbc)
end program project

subroutine work(ib,jb,its,ite,jts,jte,ims,ime,jms,jme, N, M, lbc, rbc, tbc, bbc)

integer:: its,ite,jts,jte     ! starting and ending indices of the submatrix I am updating
integer:: ims,ime,jms,jme     ! starting and ending indices of the submatrix I have in memory - at least by 1 larger on each side

double precision, dimension(ims:ime,jms:jme):: u,f,v ! dynamic memory allocation

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

!Corner Cases
if (jb > 1) then

call MPI_RECV(u(its:ite, jms), ite, MPI_DOUBLE, taskid - M, 0, MPI_COMM_WORLD, ierr)

else if (jb < N) then

call MPI_SEND(u(its:ite, jme), ite, MPI_DOUBLE, taskid +M, 0, MPI_COMM_WORLD, ierr)

else if (jb==1) then

u(ims:ime, ims) = lbc

else if (ib>1) then

call MPI_RECV(u(ims, jts:jte), jte, MPI_DOUBLE, taskid -1,1, MPI_COMM_WORLD, ierr)

else if (ib < N)
 
u(ims,jms:jme) = 100
u(ims:ime, jms) = 100
u(its:ite,jts:jte) =0
call MPI_RECV(u(ime, jts
call MPI_RECV(BUFFER1, size(BUFFER1), MPI_DOUBLE, M+1

else if (ib == M .AND. jb == 1) then

u(ime,jms:jme) = 100
u(ims:ime, jms) = 100
u(ims:ite,jts:jme) =0


else if (ib == 1 .AND. jb == N) then

u(ims,jms:jme) = 100
u(ims:ime, jme) = 100
u(its:ime,jms:jte) =0

else if (ib== M .AND. jb == N) then

u(ime,jms:jme) = 100
u(ims:ime, jme) = 100
u(ims:ite,jms:jte) =0

!EDGE CASES
else if (ib .NEQ. 1 .AND. ib .NEQ. M .AND. jb ==1) then

u(ims:ime, jms) = 100
u(ims:ime,jts:jme) =0

else if (ib .NEQ. 1 .AND. ib .NEQ. M .AND. jb == N) then

u(ims:ime, jme) = 100
u(ims:ime,jms:jte) =0


else if (ib == 1 .AND. jb .NEQ. 1 .AND. jb .NEQ. N) then

u(ims, jms:jme) = 100
u(its:ime,jms:jme) =0

else if (ib == M .AND. jb .NEQ. 1 .AND. jb .NEQ. N) then

u(ime, jms:jme) = 100
u(ims:ite,jms:jme) =0


!Interior Cases
else if (ib .NEQ.1 .AND. ib .NEQ. M .AND. jb .NEQ. 1 .AND. jb .NEQ. N) then

u(ims:ime, jms:jme) = 0

end if 
do j=jts,jte
    do i=its,ite
        v(i,j) = u(i,j) - omega * (f(i,j)*h*h - (4*u(u,j)-u(i+1,j) - u(i-1,j) - u(i,j+1) - u(i,j-1)))
    enddo
enddo

! u = v
do j=jts,jte
    do i=its,ite
        u(i,j) = u(i,j)
    enddo
enddo

! compute also the error (residual, or the change in one iteration)
! add up the errors over subdomains and print the result

! end iteration loop here

! think how to get the results out to a file!

end subroutine work
