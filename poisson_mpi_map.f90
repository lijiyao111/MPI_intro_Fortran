program poisson
use mpi

implicit none
integer :: npe_wrld,rnk_wrld,ierr
integer, parameter :: &
  n=72, &
  iblk_max=3, &
  jblk_max=3, &
  i_max=n/iblk_max, &
  j_max=n/jblk_max
integer :: i,j,ib,jb,proc,iblk,jblk,nghbr_count,req,edge,iter,nproc
integer :: proc_map(0:iblk_max+1,0:jblk_max+1),nghbr_list(4)

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,npe_wrld,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rnk_wrld,ierr)

do nproc=0,npe_wrld-1
  if(rnk_wrld==nproc) then

    print *, rnk_wrld,npe_wrld

    proc_map(:,:)=-1
    proc=0

    do jb=1,jblk_max
        do ib=1,iblk_max
            proc_map(ib,jb)=proc
            proc=proc+1
        end do
    end do

    iblk=1+mod(rnk_wrld,iblk_max)
    jblk=1+(rnk_wrld-mod(rnk_wrld,iblk_max))/iblk_max

    print *, 'iblk',iblk,'jblk',jblk

    nghbr_list(:)=(/ proc_map(iblk+1,jblk),proc_map(iblk,jblk+1),&
                    proc_map(iblk-1,jblk),proc_map(iblk,jblk-1) /)

    print *, 'nghbr_list : '
    write(*,*) (nghbr_list(i),i=1,4)

    nghbr_count=count(nghbr_list(:) /= -1)

    print *, 'nghbr_count', nghbr_count
    print *,'******************************'

  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

enddo


call MPI_FINALIZE(ierr)



end program poisson
