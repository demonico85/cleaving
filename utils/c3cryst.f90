program calcw

implicit none

integer :: ios,ntypes,iexit,i,j,matrixdim,jj,b,nargs,kk,a,filen
real*8,pointer :: intermatrix(:,:)
real*8 :: totinteraction,lambda,delta,DD
character(len=30) :: namefile,outfile,arg,serout,namelam
character(len=8) :: namebit


ntypes=2
matrixdim=ntypes+1 ! ntypes+zero col
lambda=0.0


namefile="inters3.1.dat"
open(unit=10,file=namefile,status="old",iostat=ios)
if(ios .ne. 0)then
   namefile="inters2.1.dat"
   open(unit=10,file=namefile,status="old",iostat=ios)
   if(ios .ne. 0)then
      namefile="inters1.1.dat"
      open(unit=10,file=namefile,status="old",iostat=ios)
      if(ios .ne. 0)then
        write(*,*)"ERROR: unkown .dat files"
        write(*,*)"Exiting...."
        stop
      else
         filen=1
      end if
   else
      filen=2
   end if
else
   filen=3
end if 


select case (filen)
case (1)
   namebit="inters1."
case (2)
   namebit="inters2."
case (3)
   namebit="inters3."
end select


allocate(intermatrix(matrixdim,matrixdim))


do jj=1,500
    write(namefile,*)namebit,jj,".dat"
    call StripSpaces(namefile)
    open(unit=10,file=namefile,status="old",iostat=ios)
    if(ios .ne. 0)then
        write(namefile,*)namebit,jj-1,".dat"
        call StripSpaces(namefile)
        write(*,*)"Last read file: ",namefile
        exit
    end if
    write(namelam,*)"lambda.",jj,".dat"
    call StripSpaces(namelam)
    open(unit=30,file=namelam,status="old",iostat=ios)
    if(ios .eq. 0 )then 
    write(serout,*)"ave.F.",jj,".out"
    call StripSpaces(serout)
    open(unit=20,file=serout,status="replace")
    write(20,*)
    write(20,*)

    read(10,*)
    read(10,*)
    read(10,*)

    read(30,*)
    read(30,*)
    read(30,*)a,lambda

    iexit=1
    kk=0
do 
    read(10,*,iostat=ios)b
    if(ios .ne. 0)then
        iexit = 0
        exit
    end if
    do i=1,matrixdim
        read(10,*,iostat=ios)a,(intermatrix(i,j),j=1,matrixdim)
        if(ios .ne. 0)then
            iexit = 0
            exit
        end if
    end do


    if(iexit .eq. 0)exit

    DD = intermatrix(2,3) 

    kk=kk+1
   write(20,*)kk,kk,kk,DD,lambda

end do

end if
 close(10)
 close(20)
 close(30)
end do


contains

    subroutine StripSpaces(string)
    character(len=*) :: string
    integer :: stringLen 
    integer :: last, actual

    stringLen = len (string)
    last = 1
    actual = 1

    do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        endif
    end do

    end subroutine






end program


