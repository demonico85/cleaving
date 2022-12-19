program createbox

implicit none

integer :: ios,natoms1,b,totatoms,indx,i,a
integer, pointer :: id1(:)
real(kind=8) :: xlo,xhi,ylo,yhi,zlo1,zhi1
real(kind=8) :: x1,y1,z1,vx1,vy1,vz1
character(len=50) :: file1,output


file1="Fstep3.1.data"
output="inputStep4.lmp"

open(unit=10,file=file1, status="old",iostat=ios)
if(ios .ne. 0)stop "File1 not found"

open(unit=20,file=output, status="replace")


read(10,*)
read(10,*)
read(10,*)natoms1
read(10,*)
read(10,*)
read(10,*)xlo,xhi
read(10,*)ylo,yhi
read(10,*)zlo1,zhi1
do i=1,14
    read(10,*)
end do

allocate(id1(natoms1))

id1 = 0



indx=0
do i=1,natoms1
    read(10,*) a,b,x1,y1,z1
    if(b .lt. 5) then
        indx=indx+1
    end if
end do


write(20,*)"LAMMPS data file "
write(20,*)
write(20,*)indx," atoms"
write(20,*)"4 atom types"
write(20,*)
write(20,99)xlo,xhi,"xlo xhi"
write(20,99)ylo,yhi,"ylo yhi"
write(20,99)zlo1,zhi1,"zlo zhi"
write(20,*)
write(20,*)"Atoms # atomic"
write(20,*)

rewind(10)

read(10,*)
read(10,*)
read(10,*)natoms1
read(10,*)
read(10,*)
read(10,*)xlo,xhi
read(10,*)ylo,yhi
read(10,*)zlo1,zhi1
do i=1,14
    read(10,*)
end do

id1=0


indx=1
do i=1,natoms1
    read(10,*) a,b,x1,y1,z1
    if(b .lt. 5) then
        write(20,100)indx,b,x1,y1,z1
        id1(a)=indx
        indx=indx+1
    else if(b .eq. 9)then
        write(20,100)indx,5,x1,y1,z1
        id1(i)=1
        indx=indx+1
    end if
end do


read(10,*)
read(10,*)
read(10,*)
write(20,*)
write(20,*)"Velocities"
write(20,*)

indx=1
do i=1,natoms1
    read(10,*) b,vx1,vy1,vz1
    if(id1(b) .gt. 0)then
        write(20,101)id1(b),vx1,vy1,vz1
    end if
end do



 close(21,status="delete")
99  format(2(F12.8,3X),3X,A)
100 format(I5,3X,I5,3(F10.6,3X))
101 format(I5,3X,3(F10.6,3X))
end program









































