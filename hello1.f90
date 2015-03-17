program hello1
use mpi
integer npe_wrld, &! number of processes within the world communicator
	rnk_wrld, &! rank of process within the world communicator
	ierr

real(kind=8) :: timstart,timend

call MPI_INIT(ierr)  ! initialize MPI environment
call MPI_COMM_RANK(MPI_COMM_WORLD,rnk_wrld,ierr)  ! determine world size
call MPI_COMM_SIZE(MPI_COMM_WORLD,npe_wrld,ierr)   ! determine rank within world

timstart=MPI_WTIME()
print *, "Hello world! I'm process ", rnk_wrld," out of ",&
npe_wrld, " processes."

timend=MPI_WTIME()
write(*,*) 'elapsed CPU time (s): ',timend-timstart 

call MPI_FINALIZE(ierr)  ! terminate MPI environment


end program hello1
