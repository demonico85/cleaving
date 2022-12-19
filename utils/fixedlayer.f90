program CrystalBlock

!-----------------------------------------------------------------------
! Create a block of FCC crystal with orientation 'Orient', 
! size 'nUcells(1)-by-nUcells(2)-by-nUcells(3)' basic blocks,
! and density rhoc
! Return global parameters: aa, zmin, region, regionH, nAtom
! Uses subroutines: CrystalLayrs
!-----------------------------------------------------------------------
!  use Globals

  implicit none


integer, parameter :: IOFILE=110,lrt=0,NATMAX=10000,lmps=120
integer, dimension(3) :: nUcells
integer :: i,j,k,m,nwells,nwAlls(2),nargs
integer :: nfrom,natm,nto,changecell
real(8) :: rhoc,a,b,c,red,last,kk
real(8) :: aa,cx,cy,cz,factr,R(3,NATMAX),zmin,cleaving_plane,shift(3),first
real(8),dimension(3) :: region,regionH
character(len=50) :: file1,arg

    nargs = iargc()
    if(nargs .eq. 0) then
        stop 'Noinputs'
    end if


call getarg( 1, arg )
read(arg,*) file1

rhoc=-100

call getarg( 2, arg )
read(arg,*) rhoc


changecell=0

i=1
do 
  call getarg( i, arg )
  if(arg .eq. '-c') then

    if(nargs .ne. 6)then
        stop 'Not nought argumernts for cell'
    end if

    changecell=1

    call getarg( i+1, arg )
    read(arg,*) a
    call getarg( i+2, arg )
    read(arg,*) b
    call getarg( i+3, arg )
    read(arg,*) c
    exit
  end if
  i=i+1
  if(i .eq. nargs)exit
end do






if(rhoc .lt. 0)stop "ERROR: you need the density"

open(unit=100,file="fxdlyr.tmp",status='replace')


aa = (4.d0/rhoc)**(1.d0/3.d0)

if(file1 == "111")then

    nUcells(1)=11
    nUcells(2)=6
    nUcells(3)=12

    if(changecell .eq. 1)then
        nUcells(1)=a
        nUcells(2)=b
        nUcells(3)=c
    end if

    region(1) = dble(nUcells(1))*aa/sqrt(2.d0)
    region(2) = dble(nUcells(2))*aa*sqrt(1.5d0)
    region(3) = dble(nUcells(3))*aa*sqrt(3.d0)

    cx = region(1)/dble(nUcells(1));  cy = cx*sqrt(3.d0)/12.d0
    cz = cx/sqrt(6.d0)
! La distanza del primo layer per (111) 
! dble(2*k-1)*cz con k = 1
    first = cz

! La incremento un po' sulla quarta cifra per essere sicuro che l'if < distanza prenda tutti gli atomi
!    red = first*0.0001
    first = first + 0.5*first 

! Distanza dell'ultimo layer
    kk=3.0*dble(nUcells(3))
    last=(2.0*kk-1.0)*cz

! La  diminuisco un po' sulla quarta cifra
!    red = last*0.0001
    last = last - 0.5*first


else if(file1 == "110")then

    nUcells(1)=6
    nUcells(2)=9
    nUcells(3)=24


    if(changecell .eq. 1)then
        nUcells(1)=a
        nUcells(2)=b
        nUcells(3)=c
    end if

    region(1) = dble(nUcells(1))*aa*sqrt(2.d0)
    region(2) = dble(nUcells(2))*aa
    region(3) = dble(nUcells(3))*aa/sqrt(2.d0)

    cy = region(2)/dble(nUcells(2));  cx = cy*sqrt(2.d0)/8.d0
    cz = cy/sqrt(2.d0)
! La distanza del primo layer per (110) 
!(.25d0 + .5d0*dble(k-1))*cz con k = 1
    first = 0.25*cz

! La incremento un po' sulla quarta cifra per essere sicuro che l'if < distanza prenda tutti gli atomi
!    red = first*0.0001
    first = first + 0.5*first 

! Distanza dell'ultimo layer
    kk=2.0*dble(nUcells(3))-1.0
    last = (.25d0 + .5d0*kk)*cz

! La  diminuisco un po' sulla quarta cifra
!    red = last*0.0001
    last = last - 0.5*first

else if(file1 == "100")then

    nUcells(1)=9
    nUcells(2)=9
    nUcells(3)=20

    if(changecell .eq. 1)then
        nUcells(1)=a
        nUcells(2)=b
        nUcells(3)=c
    end if

    region(1) = dble(nUcells(1))*aa
    region(2) = dble(nUcells(2))*aa
    region(3) = dble(nUcells(3))*aa

    cx = region(1)/dble(nUcells(1));  cy = cx;  cz = cx
! La distanza del primo layer per (100) 
!  (.25d0 + .5d0*dble(k-1))*cz con k = 1
    first = .25d0*cz

! La incremento un po' sulla quarta cifra per essere sicuro che l'if < distanza prenda tutti gli atomi
!    red = first*0.0001
    first = first + 0.5*first

! Distanza dell'ultimo layer
    kk=2.0*dble(nUcells(3))-1.0
    last = (.25d0 + .5d0*kk)*cz

! La  diminuisco un po' sulla quarta cifra
    red = last*0.0001
    last = last - 0.5*first


end if 

write(100,*)"cx",cx
write(100,*)"cy",cy
write(100,*)"cz",cz

write(100,*)" type crystal layer ",file1
write(100,*)"xlo",0.0
write(100,*)"xhi",region(1)
write(100,*)"ylo",0.0
write(100,*)"yhi",region(2)
write(100,*)"zlo1",0.0
write(100,*)"zhi1",0.0+first
write(100,*)"zlo2",last
!write(100,*)"zlo2",region(3)-first
write(100,*)"zhi2",region(3)
write(100,*)"zlo1a",0.0+region(3)
write(100,*)"zhi1a",0.0+first+region(3)
!write(100,*)"zlo2a",region(3)-first+region(3)
write(100,*)"zlo2a",last+region(3)
write(100,*)"zhi2a",region(3)+region(3)
write(100,*)

 close(15)

end program
