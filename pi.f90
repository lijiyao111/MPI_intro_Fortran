program pi
use mpi
implicit none

integer(kind=4) :: npe_wrld,rnk_wrld,ierr
integer(kind=8) :: n=0,i
real (kind=8) :: del_x,x_left,pi_piece,pi_approx,time_start,time_end,x

! setup MPI
call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,npe_wrld,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rnk_wrld,ierr)

! read and broadcast total number of intervals
if(rnk_wrld==0) then
  print *,'Enter the total number of intervals :'
  read(*,*) n
endif

!n=100000

call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

time_start=MPI_WTIME() ! wall clock timer start

! integrate subinterval
del_x=1.0d0/DBLE(n); x_left=DBLE(rnk_wrld)/DBLE(npe_wrld);
pi_piece=0.0d0
do i=1,n/npe_wrld
  x=x_left+del_x*(DBLE(i)-0.5d0)
  pi_piece=pi_piece+del_x*(4.0d0/(1.0d0+x**2))
enddo

! gather the pieces of the pi
call MPI_REDUCE(pi_piece,pi_approx,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

time_end=MPI_WTIME() ! wall clock timer stop

! print the approximate value
if(rnk_wrld==0) then
  write(*,"(A12,F24.20)"),"pi_approx= ",pi_approx
endif

call MPI_BARRIER(MPI_COMM_WORLD,ierr)
print "(A12,F14.10,A12,I10)","time = ",time_end-time_start," in process ", rnk_wrld

call MPI_FINALIZE(ierr)


end program pi
