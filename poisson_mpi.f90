program poisson
use mpi

implicit none
integer :: npe_wrld,rnk_wrld,ierr,msg_tag
integer, parameter :: &
  n=72, & ! global number of grid points along an edge
  iblk_max=3, & ! number domain decomposition blocks in the i-direction
  jblk_max=3, & ! number domain decomposition blocks in the j-direction
  i_max=n/iblk_max, &! local number of grid-points in the  i-direction
  j_max=n/jblk_max ! local number of grid-points in the  j-direction
integer :: i,j,ib,jb,proc,iblk,jblk,nghbr_count,req,edge,iter,nproc
integer :: proc_map(0:iblk_max+1,0:jblk_max+1),nghbr_list(4)
real(kind=kind(1.0d0)),allocatable :: alph(:,:),tmpry(:,:),beta(:,:)
! define buf_node
type buf_node
    real(kind=kind(1.0d0)),pointer :: send(:),recv(:)
end type buf_node
type(buf_node) :: buf(4)

integer,allocatable :: send_req(:),recv_req(:)
integer,allocatable :: send_status(:,:),recv_status(:,:)

character*10 :: liblk,ljblk
integer :: filenum1,filenum2
character*50 :: filename1,filename2

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,npe_wrld,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rnk_wrld,ierr)

allocate(alph(0:i_max+1,0:j_max+1))
allocate(tmpry(1:i_max,1:j_max),beta(1:i_max,1:j_max))
alph=0.0d0
tmpry=0.0d0
beta=0.0d0
print *, size(alph)


print *, rnk_wrld,npe_wrld

! allocate memory for send and recv buffers
allocate(buf(1)%send(j_max),buf(1)%recv(j_max)) ! East
allocate(buf(2)%send(i_max),buf(2)%recv(i_max)) ! South
allocate(buf(3)%send(j_max),buf(3)%recv(j_max)) ! West
allocate(buf(4)%send(i_max),buf(4)%recv(i_max)) ! North

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! set proc_map
proc_map(:,:)=-1
proc=0

do jb=1,jblk_max
    do ib=1,iblk_max
        proc_map(ib,jb)=proc
        proc=proc+1
    end do
end do

! determine position of the local process on the proc_map
iblk=1+mod(rnk_wrld,iblk_max)
jblk=1+(rnk_wrld-mod(rnk_wrld,iblk_max))/iblk_max

print *, 'iblk',iblk,'jblk',jblk

if(iblk==2 .AND. jblk==2) then  ! Set the source, in the middle sub box
    !beta(i_max/2,j_max/2)=-1.0d0
    beta(i_max/2,:)=-1.0d0
endif

! count the number of neighboring blocks
nghbr_list(:)=(/ proc_map(iblk+1,jblk),proc_map(iblk,jblk+1),&
                proc_map(iblk-1,jblk),proc_map(iblk,jblk-1) /)

print *, 'nghbr_list : '
write(*,*) (nghbr_list(i),i=1,4)

nghbr_count=count(nghbr_list(:) /= -1)

print *, 'nghbr_count', nghbr_count

print *,'******************************'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(send_req(nghbr_count),recv_req(nghbr_count))
allocate(send_status(MPI_STATUS_SIZE,nghbr_count))
allocate(recv_status(MPI_STATUS_SIZE,nghbr_count))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! start iteration
do iter =1,5000

    !!!!!!!!! Send messages
    req=0
    send_req(:)=-999
    do edge=1,4
        if (nghbr_list(edge)/=-1) then
            if(edge==1) buf(edge)%send(:)=alph(i_max,1:j_max) ! East
            if(edge==2) buf(edge)%send(:)=alph(1:i_max,j_max) ! South
            if(edge==3) buf(edge)%send(:)=alph(1,1:j_max) ! West
            if(edge==4) buf(edge)%send(:)=alph(1:i_max,1) ! North

            msg_tag=(npe_wrld+1)*rnk_wrld+nghbr_list(edge)+1
            req=req+1

            call MPI_ISEND(buf(edge)%send,SIZE(buf(edge)%send(:)),&
                            MPI_DOUBLE_PRECISION,nghbr_list(edge),msg_tag,&
                            MPI_COMM_WORLD,send_req(req),ierr)
        endif
    enddo

    !!!!!!!!!!! Receive messages
    req=0
    recv_req(:)=-999
    do edge=1,4
        if(nghbr_list(edge)/=-1) then
            buf(edge)%recv(:)=0.0d0  ! Edge of the overall grid, in this case, the boundary condition is U=0

            msg_tag=(npe_wrld+1)*nghbr_list(edge)+rnk_wrld+1
            req=req+1

            call MPI_IRECV(buf(edge)%recv,SIZE(buf(edge)%recv(:)),&
                            MPI_DOUBLE_PRECISION,nghbr_list(edge),msg_tag,&
                            MPI_COMM_WORLD,recv_req(req),ierr)
        endif
    enddo
    !write(*,*) '-- request:', (recv_req(i),i=1,nghbr_count)

    !!!!!!!!!!!!!! wait for messages to complete
    send_status(:,:)=-999
    recv_status(:,:)=-999
    call MPI_WAITALL(nghbr_count,send_req,send_status,ierr)
    call MPI_WAITALL(nghbr_count,recv_req,recv_status,ierr)

    !write(*,*) '++ request:', (recv_req(i),i=1,nghbr_count)
    !write(*,*) '++ source:', (recv_status(MPI_SOURCE,i),i=1,nghbr_count)
    !write(*,*) '++ tag:', (recv_status(MPI_TAG,i),i=1,nghbr_count)


    !!!!!!!!!!!!!!! unpack buffers
    do edge=1,4
        if(nghbr_list(edge)/=-1) then
            !print *, 'buf size', size(buf(edge)%recv)
            !print *,'imax, jmax', i_max,j_max
            !print *, 'alph edge size',size(alph(i_max+1,1:j_max))
            if(edge==1) alph(i_max+1,1:j_max) = buf(edge)%recv(:) ! East
            if(edge==2) alph(1:i_max,j_max+1) = buf(edge)%recv(:) ! South
            if(edge==3) alph(0,1:j_max) = buf(edge)%recv(:) ! West
            if(edge==4) alph(1:i_max,0) = buf(edge)%recv(:) ! North
        endif
    enddo


    !!!!!!!!!!!!!!! jacobi iteration
    do i=1,i_max
        do j=1,j_max
            tmpry(i,j)=(alph(i-1,j)+alph(i+1,j)+alph(i,j-1)+alph(i,j+1)-beta(i,j))/4.0d0
        enddo
    enddo
    alph(1:i_max,1:j_max)=tmpry(1:i_max,1:j_max)
   !!!!!!!!!!!!!!! end iteration

end do ! end iteration


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  write sub grid into seperate files
write(liblk,'(I5)') iblk
write(ljblk,'(I5)') jblk
print *, 'liblk',liblk,'ljblk',ljblk
filename1='source_I'//trim(adjustl(liblk))//'_J'//trim(adjustl(ljblk))//'.txt'
filename2='field_I'//trim(adjustl(liblk))//'_J'//trim(adjustl(ljblk))//'.txt'
write(*,*) filename1
write(*,*) filename2

filenum1=21+rnk_wrld
filenum2=61+rnk_wrld
open(unit=filenum1,file=filename1,status='replace')
open(unit=filenum2,file=filename2,status='replace')

do i=1,i_max
    write(filenum1,'(1000f10.5)') (beta(i,j),j=1,j_max)
    write(filenum2,'(1000f10.5)') (alph(i,j),j=1,j_max)
enddo

close(filenum1)
close(filenum2)

!deallocate(alpha)
!deallocate(beta)

call MPI_FINALIZE(ierr)


!deallocate(send_req,recv_req)
!deallocate(send_status,recv_status)
! note MPI arrays can only be deallocated after MPI_FINALIZE
! arrays not related with MPI functions can be deallocated before MPI_FINALIZE

end program poisson
