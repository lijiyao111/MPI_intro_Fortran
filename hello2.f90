program hello2
use mpi
implicit none
integer :: npe_wrld, rnk_wrld
integer :: n,name_len,ierr
real (kind=8) :: wall_tick,time_start,time_end
character(len=128) :: proc_name



call MPI_INIT(ierr)

call MPI_COMM_SIZE(MPI_COMM_WORLD,npe_wrld,ierr)

call MPI_COMM_RANK(MPI_COMM_WORLD,rnk_wrld,ierr)

call MPI_GET_PROCESSOR_NAME (proc_name,name_len,ierr)
!print *,'proc_name ',proc_name,'name_len',name_len

wall_tick=mpi_wtick()

if(rnk_wrld==0) write(*,'(A13,F12.8)'), "wall_tick = ",wall_tick

time_start=MPI_WTIME()


do n=0,npe_wrld-1         ! do some useless loop, just to make the output of process rank in order
  if(rnk_wrld .eq. n) then
    time_end=mpi_wtime()

    write(*,*), "hello prom proc =", rnk_wrld," of ", npe_wrld, "running on ",trim(proc_name),&
" time = ",time_end-time_start

  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)  ! every process stop here to wairt for every one to reach this point, then move on
enddo

call MPI_FINALIZE(ierr)

end program hello2
