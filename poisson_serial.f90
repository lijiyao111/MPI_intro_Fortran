program main
implicit none
integer :: n0,n,nx,ny
integer :: i,j,l,iblk_max,jblk_max
integer :: iter
real(kind=kind(1.0d0)) h
real(kind=kind(1.0d0)),allocatable :: alph(:,:),beta(:,:),tmpry(:,:)

n0=72
n=n0+2

iblk_max=3
jblk_max=3

allocate(alph(n,n),beta(n,n))
allocate(tmpry(n,n))

alph=0.0d0
beta=0.0d0

h=1.0d0/dble(n0+1)
    print *, h

! source
nx=n0/iblk_max
ny=n0/jblk_max

beta(n0/2+1,ny+1+1:2*ny+1)=-1.0d0;
!beta=beta*(h**2)*(-1)/(h**2)

!do i=1,n
! write(*,'(1000f10.5)') (beta(i,j),j=1,n)
!end do
!do i=1,n
! write(*,'(1000f10.5)') (alph(i,j),j=1,n)
!end do
write(*,*) 'before iteration!'


!iter=int(((n0+1)/pi)**2)*3
iter=5000

do l=1,iter
  
  do i=2,n-1
    do j=2,n-1
      tmpry(i,j)=(alph(i-1,j)+alph(i+1,j)+alph(i,j-1)+alph(i,j+1)-beta(i,j))/4.0d0
    end do
  end do
  alph(:,:)=tmpry(:,:)
end do

!do i=1,n
! write(*,'(1000f10.5)') (beta(i,j),j=1,n)
!end do
!do i=1,n
! write(*,'(1000f10.5)') (alph(i,j),j=1,n)
!end do


open(21,file='source.txt',status='replace')
open(22,file='field.txt',status='replace')

do i=2,n-1
   write(21,'(1000f10.5)') (beta(i,j),j=2,n-1)
   write(22,'(1000f10.5)') (alph(i,j),j=2,n-1)
end do

close(21)
close(22)


deallocate(alph,beta)

end program main
