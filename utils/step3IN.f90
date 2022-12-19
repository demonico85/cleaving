program createbox

implicit none

integer :: ios,natoms1,natoms2,b,totatoms,indx,i,bb,nargs,newtypelo,newtypehi,kk
integer,pointer :: id1(:),id2(:),type1(:),type2(:),mapping(:),mp2(:),id1v(:),id2v(:)
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
!open(unit=22,file="mapping.inp", status="replace")
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
allocate(id1v(natoms1))
allocate(id2v(natoms2))

allocate(type1(natoms1))
allocate(type2(natoms2))
do i=1,natoms1
    read(10,*) b,type1(i),x1(i),y1(i),z1(i)
    id1(i)=b
end do
read(10,*)
read(10,*)
read(10,*)


do i=1,natoms1
    read(10,*) b,vx1(b),vy1(b),vz1(b)
    id1v(b)=i
end do



do i=1,natoms2
    read(11,*) b,type2(i),x2(i),y2(i),z2(i)
    id2(i)=b
end do
read(11,*)
read(11,*)
read(11,*)


do i=1,natoms2
    read(11,*) b,vx2(b),vy2(b),vz2(b)
    id2v(b)=i
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
write(20,*)"8 atom types"
write(20,*)
write(20,99)xlo,xhi,"xlo xhi"
write(20,99)ylo,yhi,"ylo yhi"
write(20,99)totzlo,totzhi,"zlo zhi"
write(20,*)
write(20,*)"Atoms # atomic"
write(20,*)

allocate(mapping(totatoms))
mapping=0
allocate(mp2(totatoms))
indx=1

do i=1,natoms1
    kk=id1(i)
    if(z1(i) .lt. cleavwall) then
        if(type1(i) .eq. 1) then
                newtypelo=1
                newtypehi=6
        else
                newtypelo=9
                newtypehi=10
        end if
        write(20,100)indx,newtypelo,x1(i),y1(i),z1(i)
        write(20,100)indx+1,newtypehi,x1(i),y1(i),z1(i)+zlength
!        kk=idv1(i)
        write(21,101)indx,vx1(kk),vy1(kk),vz1(kk)
        write(21,101)indx+1,vx1(kk),vy1(kk),vz1(kk)
        mapping(indx)=indx+1
        indx=indx+2
    end if
end do

do i=1,natoms1
    kk=id1(i)
    if(z1(i) .gt. cleavwall) then
        if(type1(i) .eq. 1) then
                newtypelo=5
                newtypehi=4
        else
                newtypelo=10
                newtypehi=9
        end if
        write(20,100)indx,newtypehi,x1(i),y1(i),z1(i)+zlength
        write(20,100)indx+1,newtypelo,x1(i),y1(i),z1(i)
        write(21,101)indx,vx1(kk),vy1(kk),vz1(kk)
        write(21,101)indx+1,vx1(kk),vy1(kk),vz1(kk)
        mapping(indx)=indx+1
        indx=indx+2
    end if
end do


do i=1,natoms2
    kk=id2(i)
    if(z2(i) .gt. cleavwall) then
        write(20,100)indx,2,x2(i),y2(i),z2(i) 
        write(20,100)indx+1,8,x2(i),y2(i),z2(i)+zlength
        write(21,101)indx,vx2(kk),vy2(kk),vz2(kk)
        write(21,101)indx+1,vx2(kk),vy2(kk),vz2(kk)
        mapping(indx)=indx+1
        indx=indx+2
    end if
end do

do i=1,natoms2
    kk=id2(i)
    if(z2(i) .lt. cleavwall) then
        write(20,100)indx,3,x2(i),y2(i),z2(i) +zlength
        write(20,100)indx+1,7,x2(i),y2(i),z2(i)
        write(21,101)indx,vx2(kk),vy2(kk),vz2(kk)
        write(21,101)indx+1,vx2(kk),vy2(kk),vz2(kk)
        mapping(indx)=indx+1
        indx=indx+2
    end if
end do


write(20,*)
write(20,*)"Velocities"
write(20,*)

indx=indx-1

 rewind(21)

do i=1,indx
    read(21,*)b,u,v,w
    write(20,101)b,u,v,w
end do

!do i=1,indx
!    if(mapping(i) .gt. 0)write(22,*)i,mapping(i)
!end do


 close(21,status="delete")
99  format(2(F14.8,3X),3X,A)
100 format(I5,3X,I5,4X,3(F22.16,3X))
101 format(I5,3X,3(F22.16,3X))

  end program









































