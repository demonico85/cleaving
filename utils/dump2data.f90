program createbox

implicit none

integer :: ios,natoms1,nline,k,totatoms,indx,i,j,nargs,argint(3)
integer,pointer :: id1(:),id2(:),type1(:),type2(:)
real(kind=8) :: charge,a,b,c,xed,yed,zed
real(kind=8),pointer :: x(:),y(:),z(:)
character(len=50) :: file1,file2,output,arg
character(len=100) :: inpline


output="inputStep3.lmp"

    nargs = iargc()
    if(nargs .eq. 0) then
        stop "two input files"
    end if


call getarg( 1, arg )
read(arg,*) file1
call getarg( 2, arg )
read(arg,*) file2
call getarg( 3, arg )
read(arg,*) nline

call getarg( 4, arg )
read(arg,*) xed
call getarg( 5, arg )
read(arg,*) yed
call getarg( 6, arg )
read(arg,*) zed

open(unit=11,file=file1, status="old",iostat=ios)
if(ios .ne. 0)stop "File1 not found"
open(unit=10,file=file2, status="old",iostat=ios)
if(ios .ne. 0)stop "File2 not found"

output="newdata.data"
open(unit=20,file=output, status="replace")

write(*,*)
write(*,*)'*********************************************************'
write(*,*)
write(*,*)' First  file MUST be the dump'
write(*,*)' Second File MUST be the data'
write(*,*)


read(10,*)
read(10,*)
read(10,*)natoms1

rewind(10)

allocate(x(natoms1),y(natoms1),z(natoms1))
i=0
do 
    i=i+1
    read(11,*,iostat=ios)x(i),y(i),z(i)
    if(ios .ne. 0)exit
end do

do i=1,natoms1
    write(22,*)i,x(i),y(i),z(i)
end do

i=0
j=0
do 
    i=i+1
    if( i .le. nline+1) then ! +1 perche' c'e' la riga vuota dopo
        read(10,'(A)')inpline
        write(20,*)inpline

    elseif( i .le. nline+1+natoms1) then
        j=j+1
        read(10,*)(argint(k),k=1,3),charge,a,b,c
!if(i < 30) write(*,*)i,j,charge,x(j),y(j),z(j),a,b,c
        write(20,*)(argint(k),k=1,3),charge,x(j)*xed,y(j)*yed,z(j)*zed
    else
        read(10,'(A)',iostat=ios)inpline
        if(ios .ne. 0) exit
        write(20,*)inpline
    end if

end do



end program





