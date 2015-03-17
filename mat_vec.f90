program mat_vec
use mpi
implicit none

integer,parameter :: rows=10,cols=100
integer :: npe_wrld,rnk_wrld,master,i,j,count_rows,sender,row_index,ierr
integer :: status(MPI_STATUS_SIZE)
real (kind=8) :: a(rows,cols),b(cols),c(rows),buffer(cols),ans,time_start,time_end
integer :: Msendcount=0,Mreceivcount=0 ! count number of messages sent and received in the Master
integer :: Ssendcount=0,Sreceivcount=0 ! count number of messages sent and received in the Slave

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,npe_wrld,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rnk_wrld,ierr)

time_start=MPI_WTIME()
master=0
!! master
if(rnk_wrld==master) then  
  do j=1,cols ! make an arbitrary marix a and vector b
    b(j)=1.0d0
    do i=1,rows
      a(i,j)=DBLE(i+j)
    enddo
  enddo
  call MPI_BCAST(b,cols,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  count_rows=0
  do i=1,npe_wrld-1   ! first, send the first npe_wrld-1 rows to the slave processes
    print *, "i",i
    if(count_rows<rows) then
      print *,'count_rows', count_rows
      do j=1,cols
        buffer(j)=a(i,j)
     enddo
     call MPI_SEND(buffer,cols,MPI_DOUBLE_PRECISION,i,i,MPI_COMM_WORLD,ierr)
    !  Msendcount=Msendcount+1
    !print *,"Msendcount ",Msendcount,"from ",rnk_wrld
     count_rows=count_rows+1
    else  ! in case the number of rows is less than number of processes
     call MPI_SEND(MPI_BOTTOM,0,MPI_DOUBLE_PRECISION,i,0,MPI_COMM_WORLD,ierr)
    endif
  enddo

  do i=1,rows  ! second, after one slave process finished its job, it send the result to the master and get assigned a new row, untill there is no more row to work on
    print *, "i",i
    print *,'count_rows', count_rows
    call MPI_RECV(ans,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
    !  Mreceivcount=Mreceivcount+1
    !  print *,"Mreceivcount",Mreceivcount,'from ',rnk_wrld
    sender=status(MPI_SOURCE) ! get the rank of the sender 
    row_index=status(MPI_TAG) ! tag value in status is the row index
    print *,"sender,row_index",sender,row_index
    !print *,"MPI_SOURCE,MPI_TAG",MPI_SOURCE,MPI_TAG
    c(row_index)=ans
    !print *, 'row_index,c',row_index,c(row_index)

    if(count_rows<rows) then ! more wok to be done, send another row
      do j=1,cols
        buffer(j)=a(count_rows+1,j)
      enddo
      call MPI_SEND(buffer,cols,MPI_DOUBLE_PRECISION,sender,count_rows+1,MPI_COMM_WORLD,ierr)
      !Msendcount=Msendcount+1
      !print *,"Msendcount",Msendcount,'from ',rnk_wrld

      count_rows=count_rows+1
      print *,'continue to send to slave'
    else ! tell sender that there is no more work
      call MPI_SEND(MPI_BOTTOM,0,MPI_DOUBLE_PRECISION,sender,0,MPI_COMM_WORLD,ierr)
      Msendcount=Msendcount+1
      !print *,"Msendcount",Msendcount,'from ',rnk_wrld
      !print *,'no more send to slave'
    endif
  enddo
  
  time_end=MPI_WTIME()
  write(*,*),"Master elapsed CPU time (s):",time_end-time_start

!! slave
else
  call MPI_BCAST(b,cols,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  do
    call MPI_RECV(buffer,cols,MPI_DOUBLE_PRECISION,master,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
    !  Sreceivcount=Sreceivcount+1
     ! print *,"Sreceivcount",Sreceivcount,'from ',rnk_wrld
    !print *, 'received TAG:',status(MPI_TAG)

    if(status(MPI_TAG)==0) EXIT ! there is no more work

    row_index=status(MPI_TAG) ! tag value status is the row index
    ans=0.0d0
    do i=1,cols
      ans=ans+buffer(i)*b(i)
    enddo
    call MPI_SEND(ans,1,MPI_DOUBLE_PRECISION,master,row_index,MPI_COMM_WORLD,ierr)
     ! Ssendcount=Ssendcount+1
     ! print *,"Ssendcount",Ssendcount,'from ',rnk_wrld
  enddo
  time_end=MPI_WTIME()
  write(*,*),"Slave elapsed CPU time (s):",time_end-time_start,rnk_wrld

endif

call MPI_FINALIZE(ierr)

end program mat_vec
