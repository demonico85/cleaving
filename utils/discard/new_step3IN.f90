program createbox

implicit none

integer :: ios,natoms1,natoms2,b,totatoms,indx,i,bb,nargs
integer,pointer :: id1(:),id2(:),type1(:),type2(:)
real(kind=8) :: xlo,xhi,ylo,yhi,zlo1,zhi1,zlo2,zhi2,totzlo,totzhi,cleavwall,halfz,zlength,u,v,w
real(kind=8),pointer :: x1(:),y1(:),z1(:),x2(:),y2(:),z2(:),vx1(:),vy1(:),vz1(:),vx2(:),vy2(:),vz2(:)
character(len=50) :: file1,file2,output,arg


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
read(arg,*,iostat=ios) cleavwall
if(ios .ne. 0)stop "Error in the input no cleavwall parameter"

if(nargs .eq. 4) then
    call getarg( 4, arg )
    read(arg,*) output
end if


! cleavwall=16.810782908130228

open(unit=10,file=file1, status="old",iostat=ios)
if(ios .ne. 0)stop "File1 not found"
open(unit=11,file=file2, status="old",iostat=ios)
if(ios .ne. 0)stop "File2 not found"

write(*,*)
write(*,*)'*********************************************************'
write(*,*)
write(*,*)' First  file MUST be the crystal'
write(*,*)' Second File MUST be the liquid'
write(*,*)
open(unit=27,file="prova.xyz", status="replace")
open(unit=20,file=output, status="replace")
open(unit=21,file="vel.tmp", status="replace")

read(10,*)
read(10,*)
read(10,*)natoms1
read(10,*)
read(10,*)
read(10,*)xlo,xhi
read(10,*)ylo,yhi
read(10,*)zlo1,zhi1
read(10,*)
read(10,*)
read(10,*)
read(10,*)
read(10,*)
read(10,*)
read(10,*)

read(11,*)
read(11,*)
read(11,*)natoms2
read(11,*)
read(11,*)
read(11,*)xlo,xhi
read(11,*)ylo,yhi
read(11,*)zlo2,zhi2
read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*)

allocate(x1(natoms1),y1(natoms1),z1(natoms1))
allocate(x2(natoms1),y2(natoms2),z2(natoms2))

allocate(vx1(natoms1),vy1(natoms1),vz1(natoms1))
allocate(vx2(natoms1),vy2(natoms2),vz2(natoms2))

allocate(id1(natoms1))
allocate(id2(natoms2))

do i=1,natoms1
    read(10,*) id1(i),b,x1(i),y1(i),z1(i)
end do
read(10,*)
read(10,*)
read(10,*)


do i=1,natoms1
    read(10,*) b,vx1(i),vy1(i),vz1(i)
end do




do i=1,natoms2
    read(11,*) id2(i),b,x2(i),y2(i),z2(i)
end do
read(11,*)
read(11,*)
read(11,*)


do i=1,natoms2
    read(11,*) b,vx2(i),vy2(i),vz2(i)
end do




totatoms=natoms1+natoms1+natoms2+natoms2
totzlo=zlo1
halfz = cleavwall - zlo1
zlength = zhi1-zlo1
totzhi=zhi1+zlength

write(*,*)' CLeaving Wall Position1:     ', cleavwall
write(*,*)' CLeaving Wall Position2:     ', cleavwall+zlength
write(*,*)
write(*,*)'*********************************************************'

write(*,*)"Zlength",zlength
write(20,*)"LAMMPS data file "
write(20,*)
write(20,*)totatoms," atoms"
write(20,*)"4 atom types"
write(20,*)
write(20,99)xlo,xhi,"xlo xhi"
write(20,99)ylo,yhi,"ylo yhi"
write(20,99)totzlo,totzhi,"zlo zhi"
write(20,*)
write(20,*)"Atoms # atomic"
write(20,*)

indx=1
! A,B1

bb=1

do i=1,natoms1
    if(z1(i) .lt. cleavwall) then
        write(20,100)indx,1,x1(i),y1(i),z1(i)
!write(27,*)"O",x1(i),y1(i),z1(i)
    else
        write(20,100)indx,3,x1(i),y1(i),z1(i)
write(27,*)"N",x1(i),y1(i),z1(i)
    end if
    write(21,*)indx,vx1(i),vy1(i),vz1(i)
    indx=indx+1
    bb=bb+1
end do
write(*,*)"A, B1",indx-1,bb-1


! beta1, alpha
bb=1
do i=1,natoms2
    if(z2(i) .lt. cleavwall) then
        write(20,100)indx,4,x2(i),y2(i),z2(i) !+halfz
        write(21,*)indx,vx2(i),vy2(i),vz2(i)
        indx=indx+1
        bb=bb+1
write(27,*)"S",x2(i),y2(i),z2(i)
    else
        write(20,100)indx,2,x2(i),y2(i),z2(i)
        write(21,*)indx,vx2(i),vy2(i),vz1(i)
        indx=indx+1
        bb=bb+1
!write(27,*)"H",x2(i),y2(i),z2(i)
     end if
end do
write(*,*)"alpha",indx-1,bb-1

! beta1, alpha
bb=1
do i=1,natoms2
    if(z2(i) .lt. cleavwall) then
        write(20,100)indx,2,x2(i),y2(i),z2(i)+zlength
        write(21,*)indx,vx2(i),vy2(i),vz2(i)
        indx=indx+1
        bb=bb+1
!write(27,*)"C",x2(i),y2(i),z2(i)+zlength
    else
        write(20,100)indx,4,x2(i),y2(i),z2(i)+zlength
        write(21,*)indx,vx2(i),vy2(i),vz1(i)
        indx=indx+1
        bb=bb+1
!write(27,*)"He",x2(i),y2(i),z2(i)+zlength
     end if
end do

write(*,*)"beta",indx-1,bb-1

! A1,B


bb=1
do i=1,natoms1
    if(z1(i) .lt. cleavwall) then
        write(20,100)indx,3,x1(i),y1(i),z1(i)+zlength
!write(27,*)"Fe",x1(i),y1(i),z1(i)+zlength
    else
        write(20,100)indx,1,x1(i),y1(i),z1(i)+zlength
!write(27,*)"Hg",x1(i),y1(i),z1(i)+zlength
    end if
    write(21,*)indx,vx1(i),vy1(i),vz1(i)
    indx=indx+1
    bb=bb+1
end do
write(*,*)"A1,B",indx-1,bb-1


write(20,*)
write(20,*)"Velocities"
write(20,*)

indx=indx-1

 rewind(21)

do i=1,indx
    read(21,*)b,u,v,w
    write(20,101)b,u,v,w
end do

 close(21,status="delete")
99  format(2(F12.8),3X,A)
100 format(I5,3X,I5,3(F10.6,3X))
101 format(I5,3X,3(F10.6,3X))

  end program









































